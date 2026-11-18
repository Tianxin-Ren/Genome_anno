"""Execute sequence-bound functional architectures, independently of maturation presence."""
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from annotate import association_allowed, local_cluster
from reference_reactions import _coordinate_distance

FIELDS=['function_class_ids','function_class_prediction','function_rule_status','function_rule_confidence',
        'function_rule_structure_confidence','function_rule_alternatives','function_rule_evidence','function_rule_resolution']
SUPPORTED_IDENTITIES={'calibrated_high','sequence_supported_uncalibrated_core','sequence_supported_uncalibrated_noncanonical'}
HIT_FIELDS='qseqid sseqid pident evalue bitscore qstart qend sstart send qlen slen qseq sseq'.split()
EXCLUSIVE_GROUPS=[{'Fe03','Fe04','Fe05','Fe07','Fe08'},{'Ni12','Ni12R','Ni14','Ni24'}]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def write(path,rows,fields):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        w.writerows({k:json.dumps(r.get(k),ensure_ascii=False) if isinstance(r.get(k),(list,dict)) else r.get(k,'') for k in fields} for r in rows)


def load_rules():
    path=ROOT/'config/function_rules_34.json';model=json.loads(path.read_text())
    bank=ROOT/'references/processed/function_rules/members.faa'
    if sha(bank)!=model['fasta_sha256']:raise ValueError('Functional member sequence bank changed; rebuild rules')
    ids=[c['id'] for c in model['classes']]
    if not ids or len(ids)!=len(set(ids)):raise ValueError('Expected unique non-empty functional classes')
    for rid,r in model['references'].items():
        if hashlib.sha256(r['sequence'].encode()).hexdigest()!=r['sha256']:raise ValueError('Reference sequence changed: '+rid)
    return model,path,bank


def alignment_evidence(raw,reference,query,policy):
    """Count aligned residues, not the bounding span across deleted domains."""
    rp=int(raw['sstart']);qp=int(raw['qstart']);mapped={}
    for q,s in zip(raw['qseq'],raw['sseq']):
        if q!='-' and s!='-':mapped[rp]=qp
        if q!='-':qp+=1
        if s!='-':rp+=1
    coverage=len(mapped)/len(reference['sequence'])
    changes=[]
    for position,residue in reference['metal_sites']:
        observed=query[mapped[position]-1] if position in mapped else '-'
        compatible=observed==residue or {observed,residue}<={'C','U'}
        if not compatible:changes.append(dict(reference_position=position,expected=residue,observed=observed))
    domain_spans=[]
    for d in reference['domains']:
        positions={p for lo,hi in d['segments'] for p in range(lo,hi+1)}
        covered=len(positions & mapped.keys())/len(positions) if positions else 0
        domain_spans.append(dict(accession=d['accession'],aligned_reference_domain_coverage=round(covered,4)))
    failures=[]
    if coverage<policy['min_reference_coverage']:failures.append('required_member_incomplete')
    if any(d['aligned_reference_domain_coverage']<.7 for d in domain_spans):failures.append('reference_domain_span_incomplete')
    if changes:failures.append('incompatible_catalytic_ligand')
    coordination={role:[dict(reference_position=p,query_position=mapped.get(p),observed=query[mapped[p]-1] if p in mapped else '-') for p in positions]
                  for role,positions in reference.get('role_coordination_sites',{}).items()}
    return dict(protein_id=raw['qseqid'],reference_id=raw['sseqid'],identity=float(raw['pident']),
        reference_coverage=round(coverage,4),bitscore=float(raw['bitscore']),
        query_interval=[int(raw['qstart']),int(raw['qend'])],reference_interval=[int(raw['sstart']),int(raw['send'])],
        ligand_changes=changes,reference_domain_span_evidence=domain_spans,coordination_evidence=coordination,
        domain_evidence_scope='Aligned reference domains; query HMM evidence remains in the original domains stage.',
        admissible=not failures,failures=failures)


def load_query_profiles(directory,annotation,proteins):
    """Reuse the existing hash-bound HMM stage; a fold alone never assigns function."""
    if directory is None or not (Path(directory)/'query_domain_manifest.json').is_file():return {},{}
    from query_domain_organization import verify_completed
    directory=Path(directory).resolve();report=verify_completed(directory)
    if str(Path(annotation).resolve()) not in report['source_runs']:raise ValueError('Role domains belong to another annotation')
    profiles=defaultdict(list)
    for d in table(directory/'query_domains.tsv'):
        if d['protein_id'] not in proteins or d['query_id']!=hashlib.sha256(proteins[d['protein_id']].encode()).hexdigest():
            raise ValueError('Role profile sequence ownership differs')
        if d['passes'].lower()=='true':profiles[d['protein_id']].append(d)
    return profiles,{str(directory/name):sha(directory/name) for name in ['query_domain_manifest.json','query_domains.tsv']}


def member_support(hit,role,model,genes):
    policy=model.get('member_policies',{}).get(role,{})
    sites=hit.get('coordination_evidence',{}).get(role,[])
    fraction=sum(s['observed']=='C' for s in sites)/len(sites) if sites else None
    if fraction is not None and fraction<policy.get('min_conserved_ligands',1):
        return None,['required_redox_coordination_incompatible']
    if hit['admissible'] and hit['identity']>=model['search']['min_identity']:
        return dict(hit,support_route='full_reference_alignment',coordination_fraction=fraction),[]
    if 'min_identity' not in policy or hit['identity']<policy['min_identity']:
        return None,hit['failures'] or ['sequence_identity_below_role_rule']
    if policy.get('motif') and not sites:return None,['reviewed_coordination_mapping_unavailable']
    if any(x not in {'required_member_incomplete','reference_domain_span_incomplete'} for x in hit['failures']):return None,hit['failures']
    if str(genes[hit['protein_id']].get('partial','False')).lower()=='true':return None,['partial_protein_cannot_supply_complete_role']
    allowed=set(policy['profiles'])
    covered=[d for d in hit['reference_domain_span_evidence'] if d['accession'].split('.')[0] in allowed and d['aligned_reference_domain_coverage']>=policy['min_reference_domain_coverage']]
    proof=[];lo,hi=hit['query_interval']
    for d in model.get('_query_profiles',{}).get(hit['protein_id'],[]):
        a,b=int(d['query_start']),int(d['query_end'])
        if d['accession'].split('.')[0] in allowed and float(d['model_coverage'])>=policy['min_profile_coverage'] and max(0,min(hi,b)-max(lo,a)+1)/(b-a+1)>=.7:
            proof.append(d)
    if not covered or not proof:return None,['independent_complete_role_profile_missing']
    return dict(hit,admissible=True,support_route='aligned_domain_and_query_profile',legacy_failures=hit['failures'],failures=[],
                query_profile_evidence=proof,coordination_fraction=fraction),[]


def fused_scaffold_hits(role,reference,anchor,cat,model,genes):
    policy=model.get('member_policies',{}).get(role,{})
    if not policy.get('profile_only_catalytic_fusion') or str(genes[anchor].get('partial','False')).lower()=='true':return []
    result=[]
    for d in model.get('_query_profiles',{}).get(anchor,[]):
        if d['accession'].split('.')[0] not in policy['profiles'] or float(d['model_coverage'])<policy['min_profile_coverage']:continue
        interval=[int(d['query_start']),int(d['query_end'])]
        if interval[1]>=cat['query_interval'][0]:continue
        result.append(dict(protein_id=anchor,reference_id=reference,identity=0.0,identity_scope='profile_only_not_a_sequence_identity_measurement',
            reference_coverage=0.0,bitscore=float(d['domain_bitscore']),query_interval=interval,reference_interval=[],
            admissible=True,failures=[],support_route='nonoverlapping_scaffold_catalytic_fusion',query_profile_evidence=[d],coordination_fraction=None))
    return result


def donor_comparisons(hits,model):
    ids=model.get('HdrA_donor_comparison',{});classes={}
    margin_hint=model.get('modular_associations',{}).get('Ni13',{}).get('identity_margin_hint',5.0)
    for name,rid in ids.items():
        bypid={}
        for h in hits.get(rid,[]):
            if h['identity']<30 or h['reference_coverage']<.65:continue
            if h['protein_id'] not in bypid or h['bitscore']>bypid[h['protein_id']]['bitscore']:bypid[h['protein_id']]=h
        classes[name]=bypid
    result={}
    for pid in set().union(*(set(v) for v in classes.values())):
        positive=classes.get('H2',{}).get(pid);alternative=classes.get('F420_Fd_alternative',{}).get(pid)
        comparable=positive and alternative and overlaps(positive,[alternative])
        margin=positive['identity']-alternative['identity'] if comparable else None
        ratio=alternative['bitscore']/positive['bitscore'] if comparable and positive['bitscore']>0 else None
        result[pid]=dict(identity_margin=margin,identity_margin_hint=margin_hint,alternative_bitscore_ratio=ratio,
            donor_state='not_comparable' if margin is None else 'H2_favored' if margin>=margin_hint else 'alternative_favored' if margin<=-margin_hint else 'donor_unresolved',
            H2_alignment=positive,alternative_alignment=alternative,scope='Compared HdrA regions; identity margin is an uncalibrated reference hint, not a biochemical boundary or substrate assay.')
    return result


def check_association(cls,arch,anchor,selected,hits,genes,model,reference_groups,mode):
    span=model['search']['max_local_span_bp'];byrole={h['role']:h for h in selected};reasons=[]
    if cls['id']!='Ni13':
        for group in reference_groups:
            gs=[genes[byrole[r]['protein_id']] for r in group if r in byrole]
            if any(not local(a,b,span) for a in gs for b in gs):reasons.append('member_linkage_unresolved')
        return sorted(set(reasons)),{}
    for group in [['MvhA','MvhG'],['HdrB','HdrC']]:
        gs=[genes[byrole[r]['protein_id']] for r in group if r in byrole]
        if any(not local(a,b,span) for a in gs for b in gs):reasons.append('functional_module_linkage_unresolved')
    detail=dict(scope='Separate local hydrogenase and HdrBC modules within the permitted genome/bin; no physical interaction proof.')
    if 'MvhD' in byrole and not local(genes[anchor],genes[byrole['MvhD']['protein_id']],span):
        a=byrole.get('HdrA');d=byrole['MvhD']
        fused=bool(a and a['protein_id']==d['protein_id'] and not overlaps(a,[d]))
        policy=model.get('modular_associations',{}).get('Ni13',{})
        linked=bool(a and a['protein_id']!=d['protein_id'] and local_cluster(
            [genes[a['protein_id']],genes[d['protein_id']]],
            dict(local_max_span_bp=span,local_max_intergenic_gap_bp=policy.get('HdrA_MvhD_max_gap_bp',2500))))
        # Donor identity belongs to HdrA, including when MvhD is a separate protein.
        comparison=model.get('_donor_comparisons',{}).get(a['protein_id'],{}) if a else {}
        state=comparison.get('donor_state','not_comparable')
        strong_alternative=(state=='alternative_favored' and
            (comparison.get('alternative_bitscore_ratio') or 0)>=policy.get('alternative_bitscore_ratio',1.10))
        detail.update(remote_relay_module=True,remote_relay_fusion=fused,HdrA_MvhD_local_pair=linked,
            donor_comparison=comparison,donor_specificity='H2_reference_favored' if state=='H2_favored' else 'unresolved',
            support_basis='complete_redox_modules_with_fused_or_linked_HdrA_MvhD')
        if not (fused or linked):reasons.append('remote_MvhD_requires_linked_or_fused_HdrA')
        if state=='not_comparable':reasons.append('remote_HdrA_comparison_unavailable')
        if strong_alternative:reasons.append('strong_competing_F420_Fd_HdrA_reference')
        catalytic=next(m for m in arch['members'] if m['category']=='催化亚基')
        competing={h['protein_id'] for h in hits.get(catalytic['reference'],[]) if h['admissible'] and h['identity']>=model['search']['min_identity'] and association_allowed(genes[anchor],genes[h['protein_id']],mode)}
        if len(competing)>1:reasons.append('remote_donor_pairing_ambiguous_between_cores')
    return sorted(set(reasons)),detail


def local(a,b,span):
    if a['id']==b['id']:return True
    if a.get('contig')!=b.get('contig'):return False
    if any(int(g.get('start',0))<1 or int(g.get('strand',0)) not in {-1,1} for g in [a,b]):return False
    return _coordinate_distance(a,b)<=span


def overlaps(hit,selected):
    for other in selected:
        if hit['protein_id']!=other['protein_id']:continue
        a,b=hit['query_interval'];c,d=other['query_interval']
        if max(a,c)<=min(b,d):return True
    return False


def required_groups(arch,refs,span):
    """Reference proximity is an association template, never proof of transcription."""
    members=[m for m in arch['members'] if m['requirement']=='required']
    groups=[]
    for m in sorted(members,key=lambda m:(refs[m['reference']]['accession'],refs[m['reference']]['start'])):
        r=refs[m['reference']]
        if groups:
            last=refs[groups[-1][-1]['reference']]
            if last['accession']==r['accession'] and r['start']-last['end']<=span:
                groups[-1].append(m);continue
        groups.append([m])
    return [[m['role'] for m in g] for g in groups]


def evaluate_architecture(cls,arch,anchor,hits,genes,identity,model,best_catalytic_identity,mode,existing_context=None):
    refs=model['references'];policy=model['search'];members=arch['members']
    catalytic=next(m for m in members if m['category']=='催化亚基')
    cat_options=[h for h in hits.get(catalytic['reference'],[]) if h['protein_id']==anchor and h['identity']>=policy['min_identity']]
    if not cat_options:return None
    cat=max(cat_options,key=lambda h:h['bitscore'])
    if cat['identity']+policy['catalytic_identity_margin']<best_catalytic_identity:return None
    selected=[dict(cat,role=catalytic['role'])];missing=[];risks=[];optional=[];maturation=[];rejections={}
    groups=required_groups(arch,refs,policy['max_local_span_bp'])
    if arch['association']=='external':groups=[]
    alternatives=[m for m in members if m['requirement']=='alternative_external_partner']
    for member in sorted(members,key=lambda m:(m['requirement']!='required',m['role'])):
        if member is catalytic:continue
        kind=member['requirement'];options=[];rejected=set()
        candidates=[*hits.get(member['reference'],[]),*fused_scaffold_hits(member['role'],member['reference'],anchor,cat,model,genes)]
        for raw_hit in candidates:
            if raw_hit.get('support_route')=='nonoverlapping_scaffold_catalytic_fusion':h=raw_hit;why=[]
            else:h,why=member_support(raw_hit,member['role'],model,genes)
            if h is None:
                rejected.update(why);continue
            gene=genes[h['protein_id']]
            if re.fullmatch(r'FdhF\d*|FdhA|Fdh2_A',member['role'],re.I):
                counters=[c for rid in model.get('role_counters',{}).get('Fdh_catalytic',[]) for c in hits.get(rid,[]) if c['protein_id']==h['protein_id'] and c['admissible'] and c['bitscore']>h['bitscore']*1.05]
                if counters:
                    rejected.update('stronger_nonformate_partner_reference:'+c['reference_id'] for c in counters);continue
            if not association_allowed(genes[anchor],gene,mode):
                rejected.add('cross_bin_pairing' if mode=='bin' else 'outside_input_association_scope');continue
            linked=local(genes[anchor],gene,policy['max_local_span_bp'])
            # A closer competing hydrogenase must not donate its partners to this anchor.
            competitors=[other for other in hits.get(catalytic['reference'],[]) if other['admissible'] and other['protein_id']!=anchor]
            if kind=='required' and any(local(genes[other['protein_id']],gene,policy['max_local_span_bp']) and
                 genes[anchor].get('contig')!=gene.get('contig') for other in competitors):continue
            options.append((not linked,-h['bitscore'],h['protein_id'],h))
        choice=None
        for *_,h in sorted(options):
            if kind in {'required','alternative_external_partner'} and overlaps(h,selected):
                rejected.add('overlapping_roles_on_same_query_region');continue
            choice=dict(h,role=member['role'],category=member['category']);break
        if kind=='maturation':
            if choice:maturation.append(choice)
            continue
        if kind=='optional':
            if choice:optional.append(choice)
            continue
        if choice:selected.append(choice)
        elif kind=='required':
            missing.append(member['role']);rejections[member['role']]=sorted(rejected) or ['member_not_detected']
    if alternatives and not any(h['role'] in {m['role'] for m in alternatives} for h in selected):missing.append('one_external_NADP_partner')
    association_blockers,association_detail=check_association(cls,arch,anchor,selected,hits,genes,model,groups,mode)
    if len(groups)>1:risks.append('reference_members_are_distributed; genome association does not establish a physical complex')
    if arch['association']=='external':risks.append('external_reconstitution; native partnership is not inferred')
    blockers=[]
    if identity not in SUPPORTED_IDENTITIES:blockers.append('unconfirmed_catalytic_identity')
    if cat['failures']:blockers.extend(cat['failures'])
    if missing:blockers.append('required_member_incomplete')
    blockers.extend(association_blockers)
    context=existing_context or {}
    if cls['id']=='H01' and (str(context.get('preferred_family','')).lower()=='hmdii' or context.get('function_prediction')=='candidate_H2_THF_derivative_hydrogenation_potential'):
        blockers.append('HmdII_branch_cannot_inherit_H4MPT_function')
    status='reference_architecture_supported'
    if 'HmdII_branch_cannot_inherit_H4MPT_function' in blockers:status='excluded_incompatible_function'
    elif blockers:status='partial_architecture' if identity in SUPPORTED_IDENTITIES and not cat['failures'] else 'excluded_identity_or_center'
    elif arch['support_mode']=='candidate_only':status='candidate_only'
    elif arch['association']=='external':status='external_reconstitution_potential'
    variant=any(h.get('support_route','full_reference_alignment')!='full_reference_alignment' for h in selected) or association_detail.get('remote_relay_module',False)
    qualifiers=list(cls['unresolved_qualifiers'])
    if not blockers and association_detail.get('remote_relay_module') and association_detail.get('donor_specificity')=='unresolved':
        qualifiers.append('HdrA_donor_specificity_not_resolved_by_two_reference_sequences')
        risks.append('Near-tied HdrA donor references limit donor specificity; complete hydrogenase, relay and Hdr modules support the genome-level coupling prediction.')
    resolution='reference_architecture_prediction'
    label=cls['label']
    if not blockers and variant and status=='reference_architecture_supported':
        status='modular_function_prediction';resolution='genome_level_coupling_prediction' if cls['id']=='Ni13' else 'reaction_module_prediction'
        label=model.get('module_prediction_labels',{}).get(cls['id'],label)
        risks.append('Functional modules supported; exact polypeptide organization or physical coupling differs from the reference.')
        if cls['id'] in {'Ni19','Ni20'}:
            status='reference_context_only';risks.append('Methanogenic physiological context cannot be inherited from a divergent module match.')
    scores=[h['identity']/100*h['reference_coverage'] for h in selected if h['role'] in {m['role'] for m in members if m['requirement']=='required'}]
    required_count=sum(m['requirement']=='required' for m in members)
    score=.5*cat['identity']/100*cat['reference_coverage']+.5*sum(scores)/required_count
    return dict(candidate_id=anchor,function_id=cls['id'],label=label,architecture=arch['id'],system=arch['system'],prediction_scope=resolution,
        status=status,score=round(score,5),missing_required=missing,member_rejections=rejections,blockers=sorted(set(blockers)),members=selected,
        optional_members=optional,maturation_candidates=maturation,reference_association_groups=groups,modular_association=association_detail,
        maturation_scope='Genome/bin search candidates, possibly shared; absence is not a reaction exclusion.',
        unmeasured_qualifiers=qualifiers,risks=risks,
        citations=arch['citations'],reference_limitations=arch['limitations'],
        structure_confidence='high' if not blockers and min(h['identity'] for h in selected)>=80 else 'medium' if not blockers else 'low',
        function_confidence='medium' if status in {'reference_architecture_supported','modular_function_prediction'} and not any(h.get('coordination_fraction') is not None and h['coordination_fraction']<1 for h in selected) else 'low',
        confidence_scope='qualitative_sequence_prediction_not_calibrated_probability')


def decide(rows,margin=.08):
    if not rows:return dict(zip(FIELDS,[[],'功能类别未定','no_matching_architecture','unresolved','unresolved',[],'function_rules/evidence.json','not_assigned']))
    byclass={}
    for row in rows:
        old=byclass.get(row['function_id'])
        rank=lambda r:(r['status'] in {'reference_architecture_supported','modular_function_prediction'},r['status'] in {'candidate_only','external_reconstitution_potential'},not r['blockers'],r['score'])
        if old is None or rank(row)>rank(old):byclass[row['function_id']]=row
    for group in EXCLUSIVE_GROUPS:
        available=[r for r in rows if r['function_id'] in group]
        best=max((r['score'] for r in available),default=0)
        for fid,r in byclass.items():
            if fid in group and r['score']+margin<best and r['status'] in {'reference_architecture_supported','modular_function_prediction'}:
                r['status']='excluded_competing_architecture';r['blockers'].append('stronger_incompatible_reference_architecture')
        supported=[r for r in byclass.values() if r['function_id'] in group and r['status'] in {'reference_architecture_supported','modular_function_prediction'}]
        signature=lambda r: 'NAD_Fd' if r['function_id'] in {'Fe03','Fe04'} else 'NiFe_NADP' if r['function_id'] in {'Ni12','Ni12R'} else r['function_id']
        # One assayed reference can have multiple compatible reaction capabilities.
        if len({signature(r) for r in supported})>1 and len({r['system'] for r in supported})>1:
            for r in supported:
                r['status']='ambiguous_architecture';r['function_confidence']='low'
                r['risks'].append('mechanism_specificity_unresolved_between_reference_architectures')
    selected=sorted((r for r in byclass.values() if r['status'] in {'reference_architecture_supported','modular_function_prediction'}),key=lambda r:(-r['score'],r['function_id']))
    if not selected:
        selected=sorted((r for r in byclass.values() if r['status'] in {'candidate_only','external_reconstitution_potential','ambiguous_architecture'}),key=lambda r:(-r['score'],r['function_id']))[:1]
    alternatives=[dict(function_id=r['function_id'],architecture=r['architecture'],status=r['status'],score=r['score'],missing_required=r['missing_required'],blockers=r['blockers']) for r in sorted(byclass.values(),key=lambda r:-r['score'])]
    partial=any(r['status']=='partial_architecture' for r in rows)
    return dict(zip(FIELDS,[[r['function_id'] for r in selected],
        '；'.join(r['label']+('（候选）' if r['status'] not in {'reference_architecture_supported','modular_function_prediction'} else '') for r in selected) or ('催化核心保留，具体功能所需成员或关联未齐' if partial else '具体功能未赋值；身份、位点或分支限制见证据'),
        selected[0]['status'] if selected else 'partial_architecture' if partial else 'identity_or_center_unresolved',
        selected[0]['function_confidence'] if selected else 'unresolved',
        selected[0]['structure_confidence'] if selected else 'low',alternatives,'function_rules/evidence.json',selected[0].get('prediction_scope','reference_architecture_prediction') if selected else 'not_assigned']))


def annotate(annotation,out,identity_rows,threads=2,domain_evidence=None):
    annotation,out=Path(annotation).resolve(),Path(out).resolve()
    if out.exists():raise FileExistsError('Use a fresh function-rule output directory')
    if not isinstance(threads,int) or isinstance(threads,bool) or threads<1:raise ValueError('Positive threads required')
    model,policy_path,bank=load_rules()
    genes={r['id']:r for r in table(annotation/'genes.tsv')}
    proteins={r.id:str(r.seq).upper().rstrip('*') for r in SeqIO.parse(annotation/'proteins.faa','fasta')}
    if set(genes)!=set(proteins):raise ValueError('Functional rule genes/proteins ownership mismatch')
    for g in genes.values():
        for k in ['start','end','strand','contig_length']:g[k]=int(g.get(k) or 0)
    run=json.loads((annotation/'run_manifest.json').read_text());mode=run.get('input_mode','assembly')
    inputs={str(p):sha(p) for p in [annotation/'run_manifest.json',annotation/'genes.tsv',annotation/'proteins.faa',policy_path,bank,Path(__file__)]}
    domain_directory=Path(domain_evidence) if domain_evidence else out.parent.parent/'domains'
    model['_query_profiles'],domain_inputs=load_query_profiles(domain_directory,annotation,proteins)
    inputs.update(domain_inputs)
    expected=run.get('normalized_proteins_sha256')
    if expected and expected!=sha(annotation/'proteins.faa'):raise ValueError('Functional rule input protein hash mismatch')
    out.mkdir(parents=True)
    query=out/'search_U_to_X.faa';query.write_text(''.join('>'+pid+'\n'+seq.replace('U','X')+'\n' for pid,seq in proteins.items()))
    raw=out/'member_hits.tsv'
    command=['blastp','-query',str(query),'-subject',str(bank),'-out',str(raw),'-outfmt','6 '+' '.join(HIT_FIELDS),
             '-evalue',str(model['search']['max_evalue']),'-max_target_seqs','5000','-max_hsps','1','-seg','no','-comp_based_stats','0','-num_threads',str(threads)]
    with (out/'blast.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    hits=defaultdict(list)
    for line in raw.read_text().splitlines():
        fields=line.split('\t')
        if len(fields)!=len(HIT_FIELDS):raise ValueError('Malformed functional member alignment')
        h=dict(zip(HIT_FIELDS,fields));rid,pid=h['sseqid'],h['qseqid']
        if rid not in model['references'] or pid not in proteins:raise ValueError('Unknown sequence in member search')
        if any(not math.isfinite(float(h[k])) for k in ['pident','evalue','bitscore']):raise ValueError('Invalid member search score')
        if float(h['pident'])<min([model['search']['min_identity'],*[p.get('min_identity',100) for p in model.get('member_policies',{}).values()]]):continue
        hits[rid].append(alignment_evidence(h,model['references'][rid],proteins[pid],model['search']))
    model['_donor_comparisons']=donor_comparisons(hits,model)
    catalysts={m['reference'] for c in model['classes'] for a in c['architectures'] for m in a['members'] if m['category']=='催化亚基'}
    best=defaultdict(float)
    for rid in catalysts:
        for h in hits[rid]:
            if h['admissible'] and h['identity']>=model['search']['min_identity']:best[h['protein_id']]=max(best[h['protein_id']],h['identity'])
    results={};summaries={}
    for identity in identity_rows:
        qid=identity['candidate_id']
        if qid not in genes:raise ValueError('Function identity belongs to an unknown protein: '+qid)
        rows=[]
        for cls in model['classes']:
            for arch in cls['architectures']:
                r=evaluate_architecture(cls,arch,qid,hits,genes,identity['identity_status'],model,best[qid],mode,identity)
                if r:rows.append(r)
        summaries[qid]=decide(rows,model['search']['architecture_margin']);results[qid]=rows
    outputs=[out/'predictions.tsv',out/'evidence.json']
    write(outputs[0],[dict(candidate_id=q,**r) for q,r in summaries.items()],['candidate_id',*FIELDS])
    outputs[1].write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n')
    if any(sha(Path(p))!=s for p,s in inputs.items()):raise ValueError('Functional rule input changed during execution')
    manifest=dict(status='complete',version=model['version'],classes_configured=len(model['classes']),
        architectures_configured=sum(len(c['architectures']) for c in model['classes']),candidates=len(summaries),
        predicted_class_counts=dict(Counter(fid for r in summaries.values() for fid in r['function_class_ids'])),
        status_counts=dict(Counter(r['function_rule_status'] for r in summaries.values())),
        input_sha256=inputs,identity_rows_sha256=hashlib.sha256(json.dumps(identity_rows,sort_keys=True,ensure_ascii=False).encode()).hexdigest(),
        output_sha256={p.name:sha(p) for p in [*outputs,raw]},command=command,
        scope='Offline sequence functional predictions. Catalytic identity is supplied by the existing workflow. Maturation, exact reaction mechanism and physiological qualifiers remain separate.')
    (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    return summaries,manifest


def load_predictions(directory,annotation):
    directory=Path(directory);manifest=json.loads((directory/'manifest.json').read_text())
    if manifest.get('status')!='complete':raise ValueError('Functional rule execution incomplete')
    source=str(Path(annotation).resolve()/'run_manifest.json')
    if source not in manifest['input_sha256']:raise ValueError('Functional rules belong to another annotation')
    for path,digest in manifest['input_sha256'].items():
        if sha(path)!=digest:raise ValueError('Functional rule input changed: '+path)
    for name,digest in manifest['output_sha256'].items():
        if sha(directory/name)!=digest:raise ValueError('Functional rule output changed: '+name)
    return {r['candidate_id']:{k:(json.loads(r[k]) if k in {'function_class_ids','function_rule_alternatives'} else r[k]) for k in FIELDS} for r in table(directory/'predictions.tsv')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--annotation',type=Path,required=True);p.add_argument('--functions',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--threads',type=int,default=2);p.add_argument('--domain-evidence',type=Path)
    args=p.parse_args();_,manifest=annotate(args.annotation,args.out,table(args.functions/'functional_annotations.tsv'),args.threads,args.domain_evidence)
    print(json.dumps({k:manifest[k] for k in ['status','classes_configured','architectures_configured','candidates','predicted_class_counts']},ensure_ascii=False))
