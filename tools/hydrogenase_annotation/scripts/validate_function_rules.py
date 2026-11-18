#!/usr/bin/env python3
"""Replay curated members through BLAST and the same rules used by the genome workflow."""
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'research_workflow'),str(ROOT/'scripts')]
from function_rules import annotate,load_rules,sha,table,write


def validate(out):
    out=Path(out).resolve()
    if out.exists():raise FileExistsError('Use fresh validation output')
    out.mkdir(parents=True)
    model,_,_=load_rules();refs=model['references'];genes={};proteins={};identities={};expected=[]
    for cls in model['classes']:
        for arch in cls['architectures']:
            sid=arch['system']
            for m in arch['members']:
                rid=m['reference'];pid=sid+'__'+rid;r=refs[rid]
                proteins[pid]=r['sequence']
                genes[pid]=dict(id=pid,contig=sid+'__'+r['accession'],start=r['start'],end=r['end'],strand=r['strand'],bin_id=sid,contig_length=0,topology='linear')
                if m['category']=='催化亚基':
                    identities[pid]=dict(candidate_id=pid,family=arch['center'],identity_status='sequence_supported_uncalibrated_core')
                    expected.append(dict(candidate_id=pid,function_id=cls['id'],architecture=arch['id']))
    target=next(a for c in model['classes'] if c['id']=='Fe08' for a in c['architectures'] if a['system']=='FeFe_Mthermo_HydABCDEF')
    controls=[]
    for case in ['baseline','without_maturation','without_complex_HydE','HydB_fragment','catalytic_ligand_change','cross_bin_HydE','maturase_replaces_complex_HydE','unconfirmed_identity','HydE_HydF_fusion']:
        scope='control_'+case;cat_id='';role_ids={}
        for m in target['members']:
            if case=='without_maturation' and m['requirement']=='maturation':continue
            if case=='without_complex_HydE' and m['role']=='HydE':continue
            if case=='HydE_HydF_fusion' and m['role']=='HydF':continue
            rid=m['reference'];r=refs[rid];pid=scope+'__'+rid;seq=r['sequence'];role_ids[m['role']]=pid
            if case=='HydB_fragment' and m['role']=='HydB':seq=seq[:100]
            if case=='catalytic_ligand_change' and m['category']=='催化亚基':
                positions=[p for p,res in r['metal_sites'] if res in {'C','U'}];assert positions
                chars=list(seq)
                for position in positions:chars[position-1]='S'
                seq=''.join(chars)
            if case=='maturase_replaces_complex_HydE' and m['role']=='HydE':
                actual=next(x for x in target['members'] if x['role']=='hydE' and x['requirement']=='maturation')
                seq=refs[actual['reference']]['sequence']
            if case=='HydE_HydF_fusion' and m['role']=='HydE':
                other=next(x for x in target['members'] if x['role']=='HydF')
                seq+='GGG'+refs[other['reference']]['sequence']
            proteins[pid]=seq
            genes[pid]=dict(id=pid,contig=scope+'__'+r['accession'],start=r['start'],end=r['end'],strand=r['strand'],bin_id=scope,contig_length=0,topology='linear')
            if case=='cross_bin_HydE' and m['role']=='HydE':genes[pid]['bin_id']='different_bin'
            if m['category']=='催化亚基':
                cat_id=pid;identities[pid]=dict(candidate_id=pid,family='FeFe',identity_status='evidence_insufficient' if case=='unconfirmed_identity' else 'sequence_supported_uncalibrated_core')
        expected_status='reference_architecture_supported' if case in {'baseline','without_maturation','HydE_HydF_fusion'} else 'excluded_identity_or_center' if case in {'catalytic_ligand_change','unconfirmed_identity'} else 'partial_architecture'
        controls.append(dict(case=case,candidate_id=cat_id,architecture=target['id'],expected_status=expected_status))
    fhl=next(a for c in model['classes'] if c['id']=='Ni15' for a in c['architectures'] if a['system']=='NiFe_FHL1_Hyc')
    for counter in model.get('role_counters',{}).get('Fdh_catalytic',[]):
        scope='control_'+counter;cat_id=''
        for m in fhl['members']:
            r=refs[m['reference']];pid=scope+'__'+m['reference']
            proteins[pid]=refs[counter]['sequence'] if m['role']=='FdhF' else r['sequence']
            genes[pid]=dict(id=pid,contig=scope+'__'+r['accession'],start=r['start'],end=r['end'],strand=r['strand'],bin_id=scope,contig_length=0,topology='linear')
            if m['category']=='催化亚基':
                cat_id=pid;identities[pid]=dict(candidate_id=pid,family='NiFe',identity_status='sequence_supported_uncalibrated_core')
        controls.append(dict(case=counter,candidate_id=cat_id,architecture=fhl['id'],expected_status='partial_architecture'))
    annotation=out/'reference_input';annotation.mkdir()
    (annotation/'proteins.faa').write_text(''.join('>'+pid+'\n'+seq+'\n' for pid,seq in proteins.items()))
    write(annotation/'genes.tsv',list(genes.values()),list(next(iter(genes.values()))))
    (annotation/'run_manifest.json').write_text(json.dumps(dict(input_mode='bin',normalized_proteins_sha256=sha(annotation/'proteins.faa'))))
    summaries,manifest=annotate(annotation,out/'replay',list(identities.values()),2)
    evidence=json.loads((out/'replay/evidence.json').read_text());results=[]
    for e in expected:
        r=next((x for x in evidence[e['candidate_id']] if x['architecture']==e['architecture']),{})
        results.append(dict(**e,status=r.get('status','not_evaluated'),missing=r.get('missing_required',[]),blockers=r.get('blockers',[]),selected=e['function_id'] in summaries[e['candidate_id']]['function_class_ids']))
    write(out/'architecture_replay.tsv',results,['candidate_id','function_id','architecture','status','missing','blockers','selected'])
    for c in controls:
        row=next(x for x in evidence[c['candidate_id']] if x['architecture']==c['architecture'])
        c['observed_status']=row['status'];c['passed']=row['status']==c['expected_status']
    write(out/'technical_controls.tsv',controls,['case','candidate_id','architecture','expected_status','observed_status','passed'])
    report=dict(status='completed',classes=len({r['function_id'] for r in results}),architectures=len(results),
        state_counts=dict(Counter(r['status'] for r in results)),selected_architecture_rows=sum(r['selected'] for r in results),
        missing_or_excluded=[r for r in results if r['status'] not in {'reference_architecture_supported','candidate_only','external_reconstitution_potential','ambiguous_architecture'}],
        technical_controls_passed=sum(c['passed'] for c in controls),technical_controls_total=len(controls),
        scope='Development replay of source reference sequences. Core identity supplied explicitly; this does not measure discovery recall or independent genome functional accuracy.',
        policy_sha256=sha(ROOT/'config/function_rules_34.json'),engine_sha256=sha(ROOT/'research_workflow/function_rules.py'))
    (out/'summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False),flush=True)
    assert report['classes']>=1 and report['architectures']>=report['classes']
    assert all(c['passed'] for c in controls),controls
    return report


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    validate(p.parse_args().out)
