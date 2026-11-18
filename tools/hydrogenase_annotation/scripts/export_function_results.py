"""Export completed run predictions and their existing evidence without relabeling."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows,fields):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        w.writerows({k:json.dumps(row.get(k),ensure_ascii=False) if isinstance(row.get(k),(dict,list)) else row.get(k,'') for k in fields} for row in rows)


def export(run,ledger):
    run=run.resolve();ledger=read(ledger);workflow=read(run/'workflow.json')
    assert workflow['status']=='complete' and not workflow['annotation_reused']
    assert len(workflow['stages'])==7 and all(s['exit_code']==0 for s in workflow['stages'])
    for name,digest in workflow['output_sha256'].items():assert sha(run/name)==digest,name
    manifest=read(run/'functions/function_rules/manifest.json')
    for path,digest in manifest['input_sha256'].items():assert sha(Path(path))==digest,path
    assert workflow['function_rules_version']==ledger['function_rules_version']
    inputs=manifest['input_sha256']
    for name,field in [('config/function_rules_34.json','function_rules_sha256'),('research_workflow/function_rules.py','function_engine_sha256')]:
        matches=[digest for p,digest in inputs.items() if p.endswith('/hydrogenase_annotation/'+name)]
        assert matches==[ledger[field]],name
    annotation=read(run/'annotation/run_manifest.json')
    assert annotation['input_files'][0]['sha256']==ledger['sha256']
    genes={r['id']:r for r in table(run/'annotation/genes.tsv')}
    proteins={};pid=None
    for line in (run/'annotation/proteins.faa').read_text().splitlines():
        if line.startswith('>'):pid=line[1:].split()[0];proteins[pid]=''
        else:proteins[pid]+=line.strip()
    assert set(genes)==set(proteins)
    profiles={pid:[] for pid in genes}
    for d in table(run/'domains/query_domains.tsv'):
        if d['passes'].lower()=='true':profiles[d['protein_id']].append(d)
    evidence=read(run/'functions/function_rules/evidence.json')
    predictions=table(run/'functions/functional_annotations.tsv')
    summary=[];members=[];sequence_roles={};selected={};donor_comparison=[]
    raw_fields='qseqid sseqid pident evalue bitscore qstart qend sstart send qlen slen qseq sseq'.split()
    donor_hits=[dict(zip(raw_fields,line.split('\t'))) for line in (run/'functions/function_rules/member_hits.tsv').read_text().splitlines() if '\tDONOR__' in line]
    accepted={'reference_architecture_supported','modular_function_prediction','candidate_only','external_reconstitution_potential','ambiguous_architecture'}
    for prediction in predictions:
        pid=prediction['candidate_id'];classes=json.loads(prediction['function_class_ids'])
        candidates=evidence[pid];chosen=[]
        for fid in classes:
            options=[r for r in candidates if r['function_id']==fid and r['status'] in accepted]
            assert options,(pid,fid)
            chosen.append(max(options,key=lambda r:r['score']))
        if not chosen and candidates:chosen=[max(candidates,key=lambda r:(r['status']=='partial_architecture',r['score']))]
        selected[pid]=chosen
        complete_candidate=next((r for r in chosen if not classes and r['status']=='partial_architecture' and not r['missing_required']),None)
        summary.append(dict(protein=pid.rsplit('_',1)[-1],full_protein_id=pid,start=genes[pid]['start'],end=genes[pid]['end'],
            strand=genes[pid]['strand'],length_aa=len(proteins[pid]),identity=prediction['identity_status'],
            classes=classes,prediction=prediction['function_class_prediction'],confidence=prediction['function_rule_confidence'],
            status=prediction['function_rule_status'],resolution=prediction['function_rule_resolution'],
            candidate_function_id=complete_candidate['function_id'] if complete_candidate else '',
            candidate_prediction='候选：'+complete_candidate['label'] if complete_candidate else '',
            candidate_confidence=complete_candidate['function_confidence'] if complete_candidate else '',
            candidate_scope='all_members_detected_association_unresolved' if complete_candidate else '',
            coarse_preferred=prediction.get('preferred_function',''),
            architectures=[dict(function=r['function_id'],architecture=r['architecture'],status=r['status'],
                                roles=len(r['members']),proteins=len({m['protein_id'] for m in r['members']}),
                                missing=r['missing_required'],blockers=r['blockers'],limits=r['unmeasured_qualifiers']+r['risks']) for r in chosen]))
        for architecture in chosen:
            for group in ['members','optional_members','maturation_candidates']:
                for member in architecture[group]:
                    qid=member['protein_id'];assert qid in proteins
                    if member['role']=='HdrA':
                        hits={h['sseqid']:h for h in donor_hits if h['qseqid']==qid}
                        h2=hits.get('DONOR__H2_linked_HdrA_reference');alternative=hits.get('DONOR__F420_Fd_HdrA2_core_reference')
                        if h2 and alternative:
                            donor_comparison.append(dict(anchor=pid.rsplit('_',1)[-1],HdrA=qid,H2_reference_identity=float(h2['pident']),
                                F420_Fd_reference_identity=float(alternative['pident']),identity_margin=round(float(h2['pident'])-float(alternative['pident']),3),
                                identity_margin_hint=5.0,positive_margin_required=False,
                                alternative_bitscore_ratio=round(float(alternative['bitscore'])/float(h2['bitscore']),4),
                                scope='Uncalibrated reference hint; near ties do not veto complete-module predictions; no substrate assay.'))
                    sequence_roles.setdefault(qid,set()).add(member['role'])
                    members.append(dict(anchor=pid.rsplit('_',1)[-1],function=architecture['function_id'],architecture=architecture['architecture'],
                        selected_class=architecture['function_id'] in classes,group=group,role=member['role'],protein=qid.rsplit('_',1)[-1],
                        full_protein_id=qid,start=genes[qid]['start'],end=genes[qid]['end'],length_aa=len(proteins[qid]),
                        query_interval=member['query_interval'],identity=None if member.get('identity_scope') else member['identity'],
                        identity_scope=member.get('identity_scope','reference_alignment'),reference_coverage=member['reference_coverage'],
                        support_route=member.get('support_route','strict_catalytic_reference_alignment'),profiles=profiles[qid],
                        coordination=member.get('coordination_evidence',{}).get(member['role'],[])))
    write(run/'functional_summary.tsv',summary,list(summary[0]) if summary else ['protein','classes','prediction'])
    write(run/'functional_members.tsv',members,list(members[0]) if members else ['anchor','function','role','protein'])
    (run/'functional_member_sequences.faa').write_text(''.join('>'+pid+' roles='+','.join(sorted(roles))+'\n'+proteins[pid]+'\n' for pid,roles in sequence_roles.items()))
    sites=table(run/'annotation/catalytic_sites.tsv')
    write(run/'functional_catalytic_sites.tsv',[s for s in sites if s['protein_id'] in sequence_roles],list(sites[0]) if sites else ['protein_id','status'])
    (run/'functional_target_evidence.json').write_text(json.dumps(selected,ensure_ascii=False,indent=2)+'\n')
    write(run/'HdrA_donor_comparison.tsv',donor_comparison,list(donor_comparison[0]) if donor_comparison else ['anchor','HdrA','identity_margin'])
    report=dict(status='complete',seven_stages_passed=True,workflow_output_hashes_verified=True,annotation_reused=False,
        original_input=ledger['original_source'],input_sha256=ledger['sha256'],genome_bp=ledger['length_bp'],contigs=ledger['contigs'],
        predicted_proteins=len(proteins),core_anchors=workflow['core_anchors'],broad_candidates=workflow['candidates'],
        function_rules_version=workflow['function_rules_version'],rule_and_engine_match_input_ledger=True,
        class_counts=workflow['function_class_predictions'],function_evaluated=len(summary),
        with_class_prediction=sum(bool(r['classes']) for r in summary),
        with_supported_function=sum(r['status'] in {'reference_architecture_supported','modular_function_prediction'} for r in summary),
        with_complete_members_unresolved_association=sum(bool(r['candidate_function_id']) for r in summary),
        results=summary,scope='Offline sequence functional predictions. Confidence is qualitative, not a measured accuracy probability.')
    (run/'genome_annotation_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    text=[run.parent.name+' 氢酶注释结果','输入：'+str(report['genome_bp'])+' bp，单条环状序列；预测蛋白 '+str(len(proteins))+' 条。',
          '原始FNA七阶段均完成；规则版本与输入清单中的哈希均已核实。',
          '催化锚点：'+str(report['core_anchors'])+'；具体功能预测：'+str(report['with_class_prediction'])+'/'+str(len(summary))+'。','']
    for r in summary:
        text.extend(['蛋白'+r['protein']+'：'+r['prediction'],'身份：'+r['identity']+'；功能置信度：'+r['confidence']+'；判定层级：'+r['resolution']])
        if r['candidate_function_id']:text.append(r['candidate_prediction']+'；候选置信度：'+r['candidate_confidence']+'；成员齐全，跨模块关联未决；未正式赋予该功能标签。')
        text.extend([json.dumps(r['architectures'],ensure_ascii=False),''])
    if donor_comparison:text.append('HdrA供体参考竞争：'+json.dumps(donor_comparison,ensure_ascii=False))
    text.append('电子分叉、直接底物、离子选择性和体内方向等限定按逐条证据记录；未进行本样本的实验测活。')
    (run/'annotation_summary.txt').write_text('\n'.join(text)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='results'},ensure_ascii=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--input-ledger',type=Path,required=True)
    args=p.parse_args();export(args.run,args.input_ledger)
