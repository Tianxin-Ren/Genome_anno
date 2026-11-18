"""Reviewed electron-role homology independent of catalytic architecture.

A role is a sequence/domain observation; proximity does not assign a physical
partner, hydrogen turnover, electron donor/acceptor, or nucleotide specificity.
"""
import csv,hashlib,itertools,json,math,subprocess,tempfile
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT,load_blast,association_allowed,BLAST_FIELDS,local_cluster
from localization import sha,verify_source_file
from tool_paths import resolve_tool
from nife_core_projections import _metric
def load_pyridine_model(root=ROOT):
 return json.loads((Path(root)/'research_workflow/pyridine_function_model.json').read_text())
from query_domain_organization import load_models,verify_completed

VERSION='module-role-evidence-v1'
ROLE_PROFILES={'NuoF_flavin_FeS':('PF01512','PF10589'),
               'NuoE_Fer2':('PF01257',),
               'SHI_gamma_NAD_FAD':('PF00970','PF00175','PF10418'),
               'Fdh_molybdo':('PF04879','PF00384','PF01568')}
SCOPE='Reviewed electron-role homology only; association tiers are context, not H2 coupling or physical complex membership.'
def table(path):
 with Path(path).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def seqhash(seq):return hashlib.sha256(seq.encode()).hexdigest()
def load_model(root=ROOT):
 base=load_pyridine_model(root)
 path=Path(root)/'research_workflow/fdh_role_model.json'; fdh=json.loads(path.read_text())
 if fdh.get('version')!='fdh-role-evidence-v2' or fdh.get('role_class')!='Fdh_molybdo':raise ValueError('FDH role model scope mismatch')
 fdh['_source_sequence']=_fdh_source_sequence(root,dict(fdh_role_model=fdh))
 fdh['_counter_sequences']=_fdh_counter_sequences(root,dict(fdh_role_model=fdh))
 return dict(base,fdh_role_model=fdh)
def _references(root):
 rows=json.loads((root/'references/processed/references.json').read_text());meta={r['reference_id']:r for r in rows}
 if len(meta)!=len(rows):raise ValueError('Duplicate reference metadata')
 compatible={};owners={r:set() for r in ROLE_PROFILES}
 for rid,m in meta.items():
  role=m.get('reference_role_evidence',{}).get('role_class')
  reviewed=m.get('reference_role_review_status')=='supported_reference_family'
  compatible[rid]={role} if reviewed and role else set()
  if reviewed and role in ROLE_PROFILES:
   profiles=m['reference_role_evidence'].get('profile_evidence',[])
   e=[c for p in profiles if p['accession'].split('.')[0]=='PF01257' and p.get('supported') for c in p.get('components',[])]
   if role=='NuoF_flavin_FeS' and e:
    compatible[rid].add('NuoE_Fer2')
    # A profile boundary alone does not define a complete biochemical subunit.
   else:owners[role].add(rid)
 if not all(owners.values()):raise ValueError('Reviewed role prototype class is empty')
 return meta,compatible,owners

def _profile_proof(pid,rows,positions,seq,profiles,models,mincov):
 found=defaultdict(list)
 for d in rows:
  acc=d['accession'].split('.')[0]
  if d['protein_id']!=pid or acc not in profiles or str(d['passes']).lower()!='true':continue
  m=models[acc];ga=[float(m['sequence_ga']),float(m['domain_ga'])]
  if d['accession']!=m['accession'] or d['model_sha256']!=m['sha256']:raise ValueError('Role profile model mismatch')
  if d['query_id']!=seqhash(seq) or int(d['query_length'])!=len(seq):raise ValueError('Role profile query mismatch')
  if [float(d['sequence_ga']),float(d['domain_ga'])]!=ga:raise ValueError('Role profile GA mismatch')
  scores=[float(d[k]) for k in ('sequence_bitscore','domain_bitscore','model_coverage')]
  if not all(math.isfinite(x) for x in scores) or not 0<=scores[2]<=1:raise ValueError('Role profile score invalid')
  lo,hi=int(d['query_start']),int(d['query_end'])
  if not 1<=lo<=hi<=len(seq):raise ValueError('Role profile coordinates invalid')
  if scores[0]<ga[0] or scores[1]<ga[1] or scores[2]<mincov:continue
  if len(positions&set(range(lo,hi+1)))/(hi-lo+1)>=.7:found[acc].append(d)
 if not all(found[p] for p in profiles):return []
 for combination in itertools.product(*(found[p] for p in profiles)):
  if all(int(a['query_end'])<int(b['query_start']) for a,b in zip(combination,combination[1:])):return list(combination)
 return []

def _fdh_source_sequence(root,model):
 spec=model['fdh_role_model']['source_positive'];seq=''
 with (Path(root)/spec['source_genbank']).open() as handle:
  for record in SeqIO.parse(handle,'genbank'):
   for feature in record.features:
    if feature.type=='CDS' and feature.qualifiers.get('protein_id')==[spec['protein_id']]:
     seq=feature.qualifiers.get('translation',[''])[0].rstrip('*')
 if not seq or seqhash(seq)!=spec['sequence_sha256']:raise ValueError('FDH source member sequence mismatch')
 return seq

def _fdh_counter_sequences(root,model):
 out={}
 for spec in model['fdh_role_model']['source_counters']:
  with (Path(root)/spec['fasta']).open() as handle:records=list(SeqIO.parse(handle,'fasta'))
  if len(records)!=1:raise ValueError('FDH counter reference ownership mismatch')
  seq=str(records[0].seq).upper().rstrip('*')
  if seqhash(seq)!=spec['sequence_sha256']:raise ValueError('FDH counter reference sequence mismatch')
  out[spec['reference_id']]=(seq,spec['role_class'])
 return out

def _map_fdhf_sites(seq,hsp,site):
 rp=int(hsp.get('sstart',0)); qp=int(hsp.get('qstart',0)); mapped={}
 targets={'cys_or_sec':site['canonical_fdhf_site'],'his':site['canonical_fdhf_his'],'arg':site['canonical_fdhf_arg']}
 for q,r in zip(hsp.get('qseq',''),hsp.get('sseq','')):
  if r!='-' and q!='-':
   for name,target in targets.items():
    if rp==target:mapped[name]={'position':qp,'residue':seq[qp-1]}
  if r!='-':rp+=1
  if q!='-':qp+=1
 return mapped

def _fdh_chemistry(seq,canonical_hsp,model):
 site=model['fdh_role_model']['site_mapping']
 mapped=_map_fdhf_sites(seq,canonical_hsp or {},site)
 c=mapped.get('cys_or_sec',{}).get('residue','-'); h=mapped.get('his',{}).get('residue','-'); r=mapped.get('arg',{}).get('residue','-')
 compatible=c in {'C','U'}
 blockers=[] if compatible else ['FDH_Cys_Sec_site_not_supported']
 if h!='H':blockers.append('FDH_canonical_His_support_missing')
 if r!='R':blockers.append('FDH_canonical_Arg_support_missing')
 return dict(status='formate_class_compatible' if compatible else 'unresolved',site_evidence=dict(positions={k:v['position'] for k,v in mapped.items()},residues={'cys':c,'his':h,'arg':r},canonical_reference='NiFe_FHL1_Hyc__FdhF'),blockers=blockers,confidence='medium' if compatible and h=='H' and r=='R' else 'low')

def _tier(anchor,member,mode,linkage):
 if not association_allowed(anchor,member,mode):return 'excluded_input_scope'
 if anchor['id']==member['id']:return 'same_protein'
 if anchor['contig']!=member['contig']:return 'genome_or_bin_context_only'
 same_strand=anchor['strand']==member['strand']
 # Measure proximity in both orientations, then retain the association tier.
 if local_cluster([anchor,member if same_strand else dict(member,strand=anchor['strand'])],linkage):
  return 'same_contig_same_strand_local' if same_strand else 'same_contig_opposite_strand_local'
 return 'same_contig_distant'

def _coding_proof(pid, seq, gene, coding_completeness, source_inputs):
 """Accept only a validated, current-source coding proof for this member."""
 if not coding_completeness or pid not in coding_completeness:
  return None
 proof=coding_completeness[pid]
 artifact=proof.get('source_artifact',{}) if isinstance(proof,dict) else {}
 path=str(Path(str(artifact.get('path',''))).resolve()) if artifact.get('path') else ''
 allowed={str(Path(p).resolve()):h for p,h in (source_inputs or {}).items()}
 if not path or path not in allowed or artifact.get('sha256')!=allowed[path]:
  return None
 try:
  from coding_completeness import supports
  if not supports(proof,pid,expected_protein_sha256=seqhash(seq),gene=gene): return None
 except Exception:
  return None
 return {'status':proof['status'],'record_digest':proof['record_digest'],'source_artifact':dict(artifact),
         'protein_sha256':proof['protein_sha256'],'translation_match_mode':proof.get('translation_match_mode')}

def analyse(raw,domains,genes,seqs,refseqs,meta,compatible,owners,model,profile_models,candidates,mode,coding_completeness=None,source_inputs=None,fdh_aux_raw=None,fdh_aux_refseqs=None):
 """Normalize validated saved evidence; no searches or phenotype labels."""
 by={g['id']:dict(g,**{k:int(g[k] or 0) for k in ('start','end','strand')}) for g in genes}
 if len(by)!=len(genes) or set(by)!=set(seqs):raise ValueError('Role gene/protein ownership mismatch')
 if len(set(candidates))!=len(candidates) or set(candidates)-set(by):raise ValueError('Role context candidate ownership mismatch')
 owned_refs=set().union(*owners.values());byraw=defaultdict(list);bydom=defaultdict(list)
 for h in raw:
  if h['qseqid'] not in by or h['sseqid'] not in refseqs:raise ValueError('Role HSP ownership mismatch')
  byraw[h['qseqid']].append(h)
 for d in domains:
  if d['protein_id'] not in by:raise ValueError('Role domain ownership mismatch')
  bydom[d['protein_id']].append(d)
 normal=lambda s:s.upper().replace('U','X').replace('O','X')
 search_seqs={q:normal(s) for q,s in seqs.items()};search_refs={q:normal(s) for q,s in refseqs.items()}
 fdh_model=model.get('fdh_role_model',{});fdh_source=fdh_model.get('_source_sequence','');fdh_counters=fdh_model.get('_counter_sequences',{})
 aux_refs={k:normal(v) for k,v in (fdh_aux_refseqs or {}).items()} if 'normal' in locals() else {}
 gene_roles={};fusion=[]
 for pid,hits in byraw.items():
  if not any(h['sseqid'] in owned_refs and h['evalue']<=model['max_evalue'] for h in hits):continue
  qualified=[];validated=[]
  for h in hits:
   if h['evalue']>model['max_evalue']:continue
   metric=_metric(h,search_seqs,search_refs)
   validated.append(dict(hsp=h,metric=metric,positions=set(metric['query_positions'])))
   if metric['identity']>=model['min_identity'] and metric['reference_coverage']>=model['min_coverage']:
    qualified.append(dict(hsp=h,metric=metric,positions=set(metric['query_positions'])))
  roles=[]
  for role,prototype_ids in owners.items():
   options=[]
   for hit in validated:
    h=hit['hsp'];rid=h['sseqid']
    if rid not in prototype_ids or hit['metric']['reference_coverage']<model['min_coverage']:continue
    proof=_profile_proof(pid,bydom[pid],hit['positions'],seqs[pid],ROLE_PROFILES[role],profile_models,model['min_coverage'])
    options.append(dict(hit,profiles=proof))
   if not options:continue
   options.sort(key=lambda x:(bool(x['profiles']),x['metric']['identity'],x['metric']['reference_coverage'],x['hsp']['bitscore']),reverse=True)
   best=options[0];rivals=[];same=[]
   for other in qualified:
    rid=other['hsp']['sseqid']
    if role in compatible.get(rid,set()):same.append(rid);continue
    overlap=len(best['positions']&other['positions'])/len(best['positions']) if best['positions'] else 0
    if overlap>=.7 and other['metric']['identity']>=best['metric']['identity']*model['competition_identity_ratio']:
     rivals.append(dict(reference_id=rid,role_class=meta[rid].get('reference_role_evidence',{}).get('role_class','unreviewed'),paired_identity=other['metric']['identity'],owned_resource_overlap=overlap,raw_hsp=other['hsp']))
   blockers=[]
   if best['metric']['identity']<model['min_identity']:blockers.append('paired_identity_below_role_support_threshold')
   if not best['profiles']:blockers.append('complete_owned_ordered_GA_domains_not_supported')
   source_partial=by[pid].get('partial')
   partial_flag=str(source_partial).lower()=='true'
   coding_proof=_coding_proof(pid,seqs[pid].rstrip('*'),by[pid],coding_completeness,source_inputs) if partial_flag else None
   if partial_flag and coding_proof is None:blockers.append('partial_protein')
   if rivals:blockers.append('overlapping_reference_role_competition')
   role_hit=dict(role_class=role,status='supported_role_homology' if not blockers else 'unresolved_role_homology',blockers=blockers,
     best=dict(reference_id=best['hsp']['sseqid'],paired_identity=best['metric']['identity'],reference_coverage=best['metric']['reference_coverage'],whole_query_coverage=best['metric']['query_coverage'],query_interval=[min(best['positions']),max(best['positions'])],query_positions=sorted(best['positions']),profiles=best['profiles'],raw_hsp=best['hsp']),
     source_partial=source_partial,coding_proof_used=coding_proof,
     positive_reference_audit=[dict(reference_id=x['hsp']['sseqid'],paired_identity=x['metric']['identity'],reference_coverage=x['metric']['reference_coverage'],complete_profile_organization=bool(x['profiles'])) for x in options],
     compatible_reference_hits=sorted(set(same)),competition=rivals,resource_overlap=[],scope=SCOPE)
   if role=='Fdh_molybdo':
    canonical=next((x for x in options if x['hsp']['sseqid']=='NiFe_FHL1_Hyc__FdhF'),None)
    chemistry=_fdh_chemistry(seqs[pid],canonical['hsp'] if canonical else None,model)
    counter_hits=[];positive_aux=[]
    for haux in (fdh_aux_raw or []):
     if haux['qseqid']!=pid or haux['sseqid'] not in aux_refs:continue
     metric=_metric(haux,{pid:normal(seqs[pid])},aux_refs)
     record=dict(reference_id=haux['sseqid'],paired_identity=metric['identity'],reference_coverage=metric['reference_coverage'],query_coverage=metric['query_coverage'],raw_hsp=haux)
     spec=next((x for x in fdh_model.get('source_counters',[])+[fdh_model.get('source_positive',{})] if x.get('reference_id')==haux['sseqid']),{})
     if metric['identity']>=model['min_identity'] and metric['reference_coverage']>=model['min_coverage']:
      if haux['sseqid']==fdh_model.get('source_positive',{}).get('reference_id'):positive_aux.append(record)
      else:counter_hits.append(dict(record,role_class=spec.get('role_class','unreviewed')))
    best_positive=max([role_hit['best']['paired_identity']]+[x['paired_identity'] for x in positive_aux])
    stronger=[x for x in counter_hits if x['paired_identity']>=best_positive*model['competition_identity_ratio']]
    if role_hit['status']!='supported_role_homology' or stronger:
     chemistry['status']='unresolved';chemistry['confidence']='low'
     chemistry['blockers'].append('non_fdh_competition_or_role_status_unresolved' if stronger else 'fdh_role_status_unresolved')
    source_profiles=_profile_proof(pid,bydom[pid],set(range(1,len(seqs[pid])+1)),seqs[pid],fdh_model.get('profiles',[]),profile_models,model['min_coverage'])
    source_compat=[]
    legacy_source_match=seqhash(seqs[pid])==fdh_model.get('source_positive',{}).get('sequence_sha256')
    if role_hit['status']=='supported_role_homology' and source_profiles and chemistry['status']=='formate_class_compatible' and not stronger and (positive_aux or legacy_source_match) and (not partial_flag or coding_proof):
     p=positive_aux[0] if positive_aux else role_hit['best']
     source_compat=[dict(source_sequence_sha256=fdh_model['source_positive']['sequence_sha256'],reference_id=fdh_model['source_positive']['reference_id'],status='source_like_module_supported',paired_identity=p['paired_identity'],reference_coverage=p.get('reference_coverage'),query_coverage=p.get('query_coverage',p.get('whole_query_coverage')),source_profile_proof=source_profiles)]
    preferred_counter=max(counter_hits,key=lambda x:(x['paired_identity'],x['reference_coverage'],x['reference_id']),default=None)
    preferred_class=preferred_counter['role_class'] if preferred_counter and preferred_counter['paired_identity']>best_positive else 'Fdh_molybdo'
    role_hit.update(preferred_role_class=preferred_class,major_alternatives=counter_hits,fdh_chemistry=chemistry,source_member_compatibility=source_compat,
      fdh_competition_evidence=dict(query_id=pid,query_sha256=seqhash(seqs[pid]),query_length=len(seqs[pid]),fdh_max_identity=best_positive,
       strongest_counter_identity=max([float(x['paired_identity']) for x in counter_hits] or [0.0]),competition_identity_ratio=model['competition_identity_ratio'],
       preferred_counter_class=preferred_class if preferred_class!='Fdh_molybdo' else None,profile_accessions=sorted({p['accession'].split('.')[0] for p in (best.get('profiles') or [])})))
    if stronger:role_hit['blockers'].append('non_fdh_molybdo_competition_unresolved');role_hit['status']='unresolved_role_homology'
   roles.append(role_hit)
  for a,b in itertools.combinations(roles,2):
   overlap=set(a['best']['query_positions'])&set(b['best']['query_positions'])
   if overlap:a['resource_overlap'].append(b['role_class']);b['resource_overlap'].append(a['role_class'])
   elif a['status']==b['status']=='supported_role_homology':fusion.append(dict(protein_id=pid,role_classes=[a['role_class'],b['role_class']],status='disjoint_role_regions_observed',scope='sequence architecture, not a complete functional complex'))
  if roles:gene_roles[pid]=roles
 # Build the independent NuoE/NuoF-like input-module hypothesis from only
 # already supported role observations.  This is a role-combination context,
 # never a physical-complex or hydrogen-turnover claim.
 role_by_gene={pid:{r['role_class']:r for r in rs if r['status']=='supported_role_homology'} for pid,rs in gene_roles.items()}
 fused={pid for pid,rs in role_by_gene.items() if {'NuoE_Fer2','NuoF_flavin_FeS'}<=set(rs)}
 modules=[]
 for pid in sorted(fused):
  e,f=role_by_gene[pid]['NuoE_Fer2'],role_by_gene[pid]['NuoF_flavin_FeS']
  if set(e['best']['query_positions']) & set(f['best']['query_positions']): continue
  modules.append(dict(module_id='', status='supported_role_module_homology',
   members=[dict(protein_id=pid,role_class='NuoE_Fer2',association_tier='same_protein',contig=by[pid]['contig'],start=by[pid]['start'],end=by[pid]['end'],strand=by[pid]['strand']),dict(protein_id=pid,role_class='NuoF_flavin_FeS',association_tier='same_protein',contig=by[pid]['contig'],start=by[pid]['start'],end=by[pid]['end'],strand=by[pid]['strand'])],
   resources=[dict(protein_id=pid,role_class='NuoE_Fer2',query_interval=e['best']['query_interval']),dict(protein_id=pid,role_class='NuoF_flavin_FeS',query_interval=f['best']['query_interval'])],
   internal_linkage='same_protein_disjoint_role_regions', scope='complete role combination candidate; not physical assembly, H2 coupling, turnover or nucleotide specificity'))
 split=[]
 for epid,ers in role_by_gene.items():
  if 'NuoE_Fer2' not in ers or epid in fused: continue
  for fpid,frs in role_by_gene.items():
   if fpid==epid or 'NuoF_flavin_FeS' not in frs or fpid in fused: continue
   eg,fg=by[epid],by[fpid]
   if eg['contig']!=fg['contig'] or eg['strand']!=fg['strand']: continue
   tier=_tier(eg,fg,mode,model['linkage'])
   if tier!='same_contig_same_strand_local' or not association_allowed(eg,fg,mode): continue
   split.append((epid,fpid,tier))
 if split:
  e_degree=defaultdict(int);f_degree=defaultdict(int)
  for epid,fpid,_ in split:e_degree[epid]+=1;f_degree[fpid]+=1
  for epid,fpid,tier in sorted(split):
   modules.append(dict(module_id='', status='supported_role_module_homology' if e_degree[epid]==1 and f_degree[fpid]==1 else 'ambiguous_role_module_homology',
    members=[dict(protein_id=epid,role_class='NuoE_Fer2',association_tier=tier,contig=by[epid]['contig'],start=by[epid]['start'],end=by[epid]['end'],strand=by[epid]['strand']),dict(protein_id=fpid,role_class='NuoF_flavin_FeS',association_tier=tier,contig=by[fpid]['contig'],start=by[fpid]['start'],end=by[fpid]['end'],strand=by[fpid]['strand'])],
    resources=[dict(protein_id=epid,role_class='NuoE_Fer2',query_interval=role_by_gene[epid]['NuoE_Fer2']['best']['query_interval']),dict(protein_id=fpid,role_class='NuoF_flavin_FeS',query_interval=role_by_gene[fpid]['NuoF_flavin_FeS']['best']['query_interval'])],
    internal_linkage=tier, scope='complete role combination candidate; not physical assembly, H2 coupling, turnover or nucleotide specificity'))
 contexts={};bycontig=defaultdict(list);bybin=defaultdict(list)
 for module_index,module in enumerate(modules,1): module['module_id']=f'input_module_{module_index}'
 total_roles=sum(len(rs) for rs in gene_roles.values())
 for pid in gene_roles:
  bycontig[by[pid]['contig']].append(pid);bybin[by[pid].get('bin_id','')].append(pid)
 for cid in candidates:
  pool=(list(gene_roles) if mode=='genome' else bybin[by[cid]['bin_id']] if mode=='bin' and by[cid].get('bin_id') else [cid] if mode=='transcript' else bycontig[by[cid]['contig']])
  allowed=[pid for pid in pool if pid in gene_roles and association_allowed(by[cid],by[pid],mode)]
  links=[dict(protein_id=pid,role_class=r['role_class'],status=r['status'],resource=r['best']['query_interval'],association_tier=_tier(by[cid],by[pid],mode,model['linkage'])) for pid in allowed for r in gene_roles[pid]]
  allowed_ids=set(allowed)
  input_contexts=[]
  for module in modules:
   if {x['protein_id'] for x in module['members']}<=allowed_ids:
    context_module=dict(module)
    context_module['anchor_association']=[dict(protein_id=x['protein_id'],association_tier=_tier(by[cid],by[x['protein_id']],mode,model['linkage'])) for x in module['members']]
    input_contexts.append(context_module)
  contexts[cid]=dict(role_contexts=links,excluded_role_count=total_roles-sum(len(gene_roles[p]) for p in allowed),input_module_contexts=input_contexts,scope=SCOPE)
 return dict(version=VERSION,competition_identity_ratio=model['competition_identity_ratio'],gene_roles=gene_roles,candidate_contexts=contexts,observed_fusions=fusion,scope=SCOPE,
   input_modules=modules,
   fdh_source_member_classes={model['fdh_role_model']['source_positive']['sequence_sha256']:'Fdh_molybdo'},
   role_classes={k:dict(profiles=list(v),prototype_reference_ids=sorted(owners[k])) for k,v in ROLE_PROFILES.items()},
   limitations=['Unrepresented role classes are not negative evidence.','A complete profile organization supports role homology, not cofactor occupancy, hydrogenase membership or turnover.','Distant, opposite-strand and excluded-bin relations are reported, not converted into physical assembly.'])

def fdh_role_competition(inventory,protein_id,sequence):
 """Project only a current, resource-overlapping non-FDH winner for this query."""
 qsha=seqhash(sequence);base={'status':'not_observed','query_id':protein_id,'query_sha256':qsha,'evidence':[],'scope':'current_query_bound_FDH_role_competition_only'}
 if not isinstance(inventory,dict):return base
 rows=inventory.get('gene_roles',{}).get(protein_id,[])
 row=next((r for r in rows if r.get('role_class')=='Fdh_molybdo'),None)
 if not row:return base
 profiles=row.get('best',{}).get('profiles') or []
 if not profiles or {p.get('accession','').split('.')[0] for p in profiles} != {'PF04879','PF00384','PF01568'}:return base
 for p in profiles:
  if p.get('query_id')!=qsha or int(p.get('query_length',-1))!=len(sequence):
   raise ValueError('FDH role profile query binding mismatch')
 evidence=row.get('fdh_competition_evidence') or {}
 if evidence.get('query_id')!=protein_id or evidence.get('query_sha256')!=qsha or int(evidence.get('query_length',-1))!=len(sequence):return base
 try:fdh_identity=float(evidence['fdh_max_identity']);ratio=float(evidence['competition_identity_ratio'])
 except (KeyError,TypeError,ValueError):return base
 best=row.get('best',{})
 norm=sequence.upper().replace('U','X').replace('O','X');resources=set(best.get('query_positions',[]))
 valid=[];close=[]
 for counter in row.get('major_alternatives',[]) or []:
  h=counter.get('raw_hsp',{});qstart,qend=int(h.get('qstart',0)),int(h.get('qend',0))
  qseq=str(h.get('qseq','')).replace('-','').upper().replace('U','X').replace('O','X')
  if h.get('qseqid')!=protein_id or int(h.get('qlen',-1))!=len(sequence) or not (1<=qstart<=qend<=len(sequence)) or qseq!=norm[qstart-1:qend]:continue
  overlap=len(resources & set(range(qstart,qend+1)))/len(resources) if resources else 0
  if overlap<.7:continue
  evidence=dict(reference_id=counter.get('reference_id'),role_class=counter.get('role_class'),paired_identity=float(counter.get('paired_identity',0)),resource_overlap=overlap)
  if counter.get('role_class')==row.get('preferred_role_class') and fdh_identity < evidence['paired_identity']*ratio:valid.append(evidence)
  elif counter.get('role_class')==row.get('preferred_role_class'):close.append(evidence)
 if valid:
  winner=max(valid,key=lambda x:(x['paired_identity'],x['reference_id'] or ''))
  base.update(status='non_fdh_preferred',blocker='preferred_non_FDH_role_competition',evidence=[dict(fdh_max_identity=fdh_identity,strongest_counter=winner,competition_identity_ratio=ratio)])
 elif close:
  base.update(status='ambiguous',evidence=[dict(fdh_max_identity=fdh_identity,close_counters=close,competition_identity_ratio=ratio)])
 else:
  base.update(status='fdh_preferred',evidence=[dict(fdh_max_identity=fdh_identity,strongest_counter_identity=max([float(x.get('paired_identity',0)) for x in row.get('major_alternatives',[]) or []] or [0]),competition_identity_ratio=ratio)])
 return base

def prepare(annotation,domain_dir,candidates,output,root=ROOT,coding_completeness=None):
 annotation,domain_dir,output,root=map(lambda p:Path(p).resolve(),(annotation,domain_dir,output,root))
 if output.exists():raise FileExistsError('Use a fresh role-evidence directory')
 output.mkdir(parents=True,exist_ok=False)
 manifest=json.loads((annotation/'run_manifest.json').read_text())
 if manifest['status']!='complete':raise ValueError('Role annotation source is incomplete')
 report=verify_completed(domain_dir,root)
 if str(annotation) not in report['source_runs']:raise ValueError('Role domains belong to another annotation')
 model=load_model(root);_,profile_by_name,_,domain_bindings=load_models(root)
 profile_models={p['accession'].split('.')[0]:p for p in profile_by_name.values()}
 paths=[annotation/n for n in ('blast.tsv','proteins.faa','genes.tsv','run_manifest.json')]
 for n in ('proteins.faa','genes.tsv'):verify_source_file(annotation,manifest,n)
 if manifest.get('output_sha256',{}).get('blast.tsv'):verify_source_file(annotation,manifest,'blast.tsv')
 if manifest.get('reference_fasta_sha256')!=sha(root/'references/processed/references.faa'):raise ValueError('Role HSP reference library differs from source search')
 paths += [domain_dir/'query_domain_manifest.json',domain_dir/'query_domains.tsv',root/'references/processed/references.json',root/'references/processed/references.faa',root/'research_workflow/pyridine_function_model.json',Path(__file__),root/'research_workflow/nife_core_projections.py']
 fdh_paths=[root/'research_workflow/fdh_role_model.json']
 paths += fdh_paths
 source_inputs={str(Path(x['path']).resolve()):x['sha256'] for x in manifest.get('input_files',[]) if isinstance(x,dict) and x.get('path') and x.get('sha256')}
 if coding_completeness:
  paths.append(root/'research_workflow/coding_completeness.py')
  for proof in coding_completeness.values():
   if isinstance(proof,dict) and proof.get('source_artifact',{}).get('path'):
    path=Path(proof['source_artifact']['path']).resolve()
    if str(path) in source_inputs and path.is_file():paths.append(path)
  before={str(p):sha(p) for p in paths};before.update(domain_bindings)
 meta,compatible,owners=_references(root)
 with (annotation/'proteins.faa').open() as h:records=list(SeqIO.parse(h,'fasta'))
 seqs={r.id:str(r.seq).upper().rstrip('*') for r in records}
 if len(seqs)!=len(records):raise ValueError('Duplicate role query protein')
 with (root/'references/processed/references.faa').open() as h:records=list(SeqIO.parse(h,'fasta'))
 refseqs={r.id:str(r.seq).upper().rstrip('*') for r in records}
 if len(refseqs)!=len(records) or set(refseqs)!=set(meta):raise ValueError('Role reference ownership mismatch')
 for rid in set().union(*owners.values()):
  if meta[rid]['sequence_sha256']!=seqhash(refseqs[rid]):raise ValueError('Reviewed role reference sequence mismatch')
 ids=[c if isinstance(c,str) else c['protein_id'] for c in candidates]
 domains=[d for d in table(domain_dir/'query_domains.tsv') if d['source_run']==str(annotation)]
 raw=load_blast(annotation/'blast.tsv',dict(blast_evalue=model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
 aux_raw=[];aux_refs={}
 fdh_model=model.get('fdh_role_model',{})
 aux_refs[fdh_model['source_positive']['reference_id']]=fdh_model['_source_sequence']
 aux_refs.update({rid:seq for rid,(seq,_) in fdh_model.get('_counter_sequences',{}).items()})
 with tempfile.TemporaryDirectory(prefix='fdh_aux_') as td:
  td=Path(td);qfa=td/'query.faa';sfa=td/'fdh_subjects.faa';outfa=td/'fdh_auxiliary.tsv'
  with qfa.open('w') as h:
   for pid,seq in seqs.items():h.write(f'>{pid}\n{seq.upper().replace("U","X").replace("O","X")}\n')
  with sfa.open('w') as h:
   for rid,seq in aux_refs.items():h.write(f'>{rid}\n{seq.upper().replace("U","X").replace("O","X")}\n')
  binary=resolve_tool(root, 'blastp')
  cmd=[str(binary),'-query',str(qfa),'-subject',str(sfa),'-evalue',str(model['max_evalue']),'-max_target_seqs',str(len(aux_refs)),'-max_hsps','1','-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(outfa)]
  subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
  aux_raw=load_blast(outfa,dict(blast_evalue=model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
  (output/'fdh_auxiliary_alignments.tsv').write_text(outfa.read_text())
  (output/'fdh_auxiliary_references.faa').write_text(sfa.read_text())
  aux_manifest={'command':cmd,'reference_ids':list(aux_refs),'query_count':len(seqs),'alignment_sha256':sha(output/'fdh_auxiliary_alignments.tsv'),'reference_fasta_sha256':sha(output/'fdh_auxiliary_references.faa')}
  (output/'fdh_auxiliary_manifest.json').write_text(json.dumps(aux_manifest,indent=2)+'\n')
  inventory=analyse(raw,domains,table(annotation/'genes.tsv'),seqs,refseqs,meta,compatible,owners,model,profile_models,ids,manifest['input_mode'],coding_completeness=coding_completeness,source_inputs=source_inputs,fdh_aux_raw=aux_raw,fdh_aux_refseqs=aux_refs)
 if any(sha(p)!=h for p,h in before.items()):raise ValueError('Role source changed during analysis')
 inventory['input_sha256']=before;inventory['fdh_auxiliary']=aux_manifest
 (output/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
 with (output/'gene_roles.tsv').open('w') as h:
  w=csv.DictWriter(h,fieldnames=['protein_id','role_evidence_json'],delimiter='\t');w.writeheader()
  w.writerows(dict(protein_id=pid,role_evidence_json=json.dumps(rs)) for pid,rs in inventory['gene_roles'].items())
 return inventory,before
