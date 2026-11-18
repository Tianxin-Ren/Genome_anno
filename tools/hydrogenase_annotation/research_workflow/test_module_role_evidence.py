"""Actual cached role evidence and causal ownership/GA/scope checks."""
import copy,json,tempfile,unittest
from pathlib import Path
from Bio import SeqIO
from annotate import load_blast
from query_domain_organization import load_models,annotate_domains
from module_role_evidence import ROOT,prepare,analyse,_references,table,load_model,fdh_role_competition,seqhash
ROUND=ROOT/'results/round108_Hadar12_release_v2'
def source(name):
 case=ROUND/name;annotation=Path(json.loads((case/'workflow.json').read_text())['annotation'])
 return annotation,case/'domains'
def dataset(name):
 a,d=source(name);model=load_model(ROOT);meta,compatible,owners=_references(ROOT)
 with (a/'proteins.faa').open() as h:seqs={r.id:str(r.seq) for r in SeqIO.parse(h,'fasta')}
 with (ROOT/'references/processed/references.faa').open() as h:refs={r.id:str(r.seq) for r in SeqIO.parse(h,'fasta')}
 _,models,_,_=load_models(ROOT)
 return dict(raw=load_blast(a/'blast.tsv',dict(blast_evalue=model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0)),domains=table(d/'query_domains.tsv'),genes=table(a/'genes.tsv'),seqs=seqs,refseqs=refs,meta=meta,compatible=compatible,owners=owners,model=model,profile_models={p['accession'].split('.')[0]:p for p in models.values()},candidates=[],mode='genome')
def role(inv,name):return next((pid,r) for pid,rs in inv['gene_roles'].items() for r in rs if r['role_class']==name)
class RoleEvidenceTests(unittest.TestCase):
 def test_fdh_projection_keeps_close_auxiliary_and_disjoint_explanations(self):
  base=ROOT/'results/round149_before_Pa_FHL_native_2NAP_v1'
  inv=json.loads((base/'functions/independent_roles/inventory.json').read_text());pid='FHL92_011'
  with (base/'annotation/proteins.faa').open() as h:seq=next(str(r.seq) for r in SeqIO.parse(h,'fasta') if r.id==pid)
  row=next(r for r in inv['gene_roles'][pid] if r['role_class']=='Fdh_molybdo')
  row['fdh_competition_evidence']=dict(query_id=pid,query_sha256=seqhash(seq),query_length=len(seq),fdh_max_identity=38.67,competition_identity_ratio=.9)
  original=copy.deepcopy(inv)
  self.assertEqual(fdh_role_competition(inv,pid,seq)['status'],'non_fdh_preferred')
  self.assertEqual(inv,original)
  # The auxiliary source may be much closer than the four generic prototypes.
  row['fdh_competition_evidence']['fdh_max_identity']=99
  self.assertEqual(fdh_role_competition(inv,pid,seq)['status'],'ambiguous')
  row['fdh_competition_evidence']['fdh_max_identity']=100;row['preferred_role_class']='Fdh_molybdo'
  self.assertEqual(fdh_role_competition(inv,pid,seq)['status'],'fdh_preferred')
  fusion=copy.deepcopy(original);r=fusion['gene_roles'][pid][0];joined=seq+seq
  for p in r['best']['profiles']:p.update(query_id=seqhash(joined),query_length=len(joined))
  r['fdh_competition_evidence'].update(query_sha256=seqhash(joined),query_length=len(joined))
  counter=next(c for c in r['major_alternatives'] if c['paired_identity']==100)
  h=counter['raw_hsp'];h.update(qstart=h['qstart']+len(seq),qend=h['qend']+len(seq),qlen=len(joined))
  r['major_alternatives']=[counter]
  self.assertNotEqual(fdh_role_competition(fusion,pid,joined)['status'],'non_fdh_preferred')
 def test_fdh_role_competition_requires_current_bound_hsp(self):
  self.assertEqual(fdh_role_competition({},'q','M')['status'],'not_observed')
  inv=self.fdh_case('round145_KIST612_registered_old_consumer_v1','ELI_0994');row=next(r for r in inv['gene_roles']['ELI_0994'] if r['role_class']=='Fdh_molybdo')
  base=ROOT/'results/round145_KIST612_registered_old_consumer_v1';
  with (base/'annotation'/'proteins.faa').open() as h: seqs={r.id:str(r.seq) for r in SeqIO.parse(h,'fasta')}
  seq=seqs['ELI_0994']; row['preferred_role_class']='nitrate_reductase_counter'; h=copy.deepcopy(row['best']['raw_hsp']); h['qseqid']='ELI_0994'; row['major_alternatives']=[{'reference_id':'counter','role_class':'nitrate_reductase_counter','paired_identity':row['best']['paired_identity']*1.2,'raw_hsp':h}]; row['fdh_competition_evidence']={'query_id':'ELI_0994','query_sha256':seqhash(seq),'query_length':len(seq),'fdh_max_identity':row['best']['paired_identity'],'competition_identity_ratio':.9}; inv['gene_roles']['ELI_0994']=[row]
  result=fdh_role_competition(inv,'ELI_0994',seq);self.assertEqual(result['status'],'non_fdh_preferred');self.assertEqual(result['blocker'],'preferred_non_FDH_role_competition')
  stale=copy.deepcopy(inv);stale['gene_roles']['ELI_0994'][0]['best']['profiles'][0]['query_id']='stale'
  with self.assertRaises(ValueError):fdh_role_competition(stale,'ELI_0994',seq)
 def test_fdh_model_requires_bound_distinct_references(self):
  from unittest.mock import patch
  original=json.loads((ROOT/'research_workflow/fdh_role_model.json').read_text())
  for change in ('empty_sources','unbound_counter','duplicate_reference_id'):
   with self.subTest(change=change),tempfile.TemporaryDirectory() as t:
    root=Path(t);(root/'research_workflow').mkdir();model=copy.deepcopy(original)
    if change=='empty_sources':model['sources']={}
    elif change=='unbound_counter':model['source_counters'][0]['fasta']='unbound.faa'
    else:model['source_counters'][0]['reference_id']=model['source_positive']['reference_id']
    (root/'research_workflow/fdh_role_model.json').write_text(json.dumps(model))
    with patch('module_role_evidence.load_pyridine_model',return_value={}):
     with self.assertRaises(ValueError):load_model(root)

 def fdh_case(self,case,pid):
  base=ROOT/'results'/case;annotation=base/'annotation';domains=base/'domains';model=load_model(ROOT);meta,compatible,owners=_references(ROOT)
  with (annotation/'proteins.faa').open() as h:seqs={r.id:str(r.seq) for r in SeqIO.parse(h,'fasta')}
  with (ROOT/'references/processed/references.faa').open() as h:refs={r.id:str(r.seq) for r in SeqIO.parse(h,'fasta')}
  _,models,_,_=load_models(ROOT);raw=load_blast(annotation/'blast.tsv',dict(blast_evalue=model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
  return analyse(raw,table(domains/'query_domains.tsv'),table(annotation/'genes.tsv'),seqs,refs,meta,compatible,owners,model,{p['accession'].split('.')[0]:p for p in models.values()},[pid],'genome')

 def test_fdh_source_member_and_chemistry_are_exported(self):
  inv=self.fdh_case('round145_KIST612_registered_old_consumer_v1','ELI_0994')
  hit=next(r for r in inv['gene_roles']['ELI_0994'] if r['role_class']=='Fdh_molybdo')
  self.assertEqual(hit['status'],'supported_role_homology')
  self.assertEqual(hit['fdh_chemistry']['status'],'formate_class_compatible')
  self.assertEqual(hit['fdh_chemistry']['site_evidence']['residues'],{'cys':'C','his':'H','arg':'R'})
  self.assertEqual(hit['source_member_compatibility'][0]['source_sequence_sha256'],'e4fbacef321d8e9645dc090d8aaaf8df376f28f98910b665c71175dfd2ab4699')
  self.assertEqual(inv['fdh_source_member_classes']['e4fbacef321d8e9645dc090d8aaaf8df376f28f98910b665c71175dfd2ab4699'],'Fdh_molybdo')

 def test_fdh_cofactor_site_variant_keeps_family_but_withdraws_chemistry(self):
  base=ROOT/'results/round145_KIST612_A_02_FdhA_sequence_variant_after_v2'
  with tempfile.TemporaryDirectory() as t:
   domains=Path(t)/'domains';annotate_domains([base/'annotation'],domains)
   inv,_=prepare(base/'annotation',domains,['ELI_0994_variant'],Path(t)/'inventory',root=ROOT)
  hit=next(r for r in inv['gene_roles']['ELI_0994_variant'] if r['role_class']=='Fdh_molybdo')
  self.assertEqual(hit['status'],'supported_role_homology')
  self.assertEqual(hit['fdh_chemistry']['status'],'formate_class_compatible')
  self.assertTrue(hit['source_member_compatibility'])
  self.assertEqual(hit['source_member_compatibility'][0]['reference_id'],'FDH_source__ADO35985.1')

 def test_fdh_role_model_keeps_non_fdh_counters_as_competitors(self):
  model=load_model(ROOT);counters=model['fdh_role_model']['_counter_sequences']
  self.assertEqual({v[1] for v in counters.values()},{'DMSO_reductase_counter','nitrate_reductase_counter','polysulfide_reductase_counter'})
  self.assertEqual(len(model['fdh_role_model']['profiles']),4)
 def copy_gene(self,data,pid,new_id,contig=None):
  gene=copy.deepcopy(next(g for g in data['genes'] if g['id']==pid));gene['id']=new_id
  if contig is not None:gene['contig']=contig
  data['genes'].append(gene);data['seqs'][new_id]=data['seqs'][pid]
  data['raw'] += [dict(h,qseqid=new_id) for h in list(data['raw']) if h['qseqid']==pid]
  data['domains'] += [dict(d,protein_id=new_id) for d in list(data['domains']) if d['protein_id']==pid]

 def test_shared_local_partner_is_ambiguous_but_separate_pairs_are_independent(self):
  data=dataset('Amobile_A');base=analyse(**data)
  e=role(base,'NuoE_Fer2')[0];f=role(base,'NuoF_flavin_FeS')[0]
  shared=copy.deepcopy(data);self.copy_gene(shared,f,'second_F')
  modules=analyse(**shared)['input_modules']
  self.assertEqual(len(modules),2)
  self.assertTrue(all(m['status']=='ambiguous_role_module_homology' for m in modules))
  separate=copy.deepcopy(data)
  self.copy_gene(separate,e,'other_E','other_contig');self.copy_gene(separate,f,'other_F','other_contig')
  modules=analyse(**separate)['input_modules']
  self.assertEqual(len(modules),2)
  self.assertTrue(all(m['status']=='supported_role_module_homology' for m in modules))

 def run_case(self,name,candidates=None):
  a,d=source(name)
  with tempfile.TemporaryDirectory() as t:
   d=Path(t)/'domains';annotate_domains([a],d)
   return prepare(a,d,candidates or [],Path(t)/'inventory')[0]
 def test_native_amobile_complete_roles(self):
  inv=self.run_case('q01_Amobile_missing_relay')
  for name in ('NuoF_flavin_FeS','NuoE_Fer2'):self.assertEqual(role(inv,name)[1]['status'],'supported_role_homology')
  self.assertEqual(len(inv['input_modules']),1)
  self.assertEqual(inv['input_modules'][0]['status'],'supported_role_module_homology')
 def test_candidate_context_contains_only_allowed_input_module_members(self):
  inv=self.run_case('q01_Amobile_missing_relay',['q0005'])
  self.assertEqual(len(inv['candidate_contexts']['q0005']['input_module_contexts']),1)
  self.assertEqual({m['protein_id'] for m in inv['candidate_contexts']['q0005']['input_module_contexts'][0]['members']},{'q0005','q0004'})
  tiers={x['protein_id']:x['association_tier'] for x in inv['candidate_contexts']['q0005']['input_module_contexts'][0]['anchor_association']}
  self.assertEqual(tiers['q0005'],'same_protein');self.assertEqual(tiers['q0004'],'same_contig_same_strand_local')
  self.assertEqual({m['role_class'] for m in inv['input_modules'][0]['members']},{'NuoE_Fer2','NuoF_flavin_FeS'})

 def test_actual_edge_complete_coding_proof_overrides_only_partial_blocker(self):
  case=ROOT/'results/round119_complete_edge_first_v1/Syntrophus_complete_edge'
  annotation=case/'annotation'; pid='edge_Syntrophus_2'
  domains=Path(self.enterContext(tempfile.TemporaryDirectory()))/'domains'
  annotate_domains([annotation],domains)
  rows=json.loads((case/'functions/coding_completeness/coding_completeness.json').read_text())['rows']
  proof={r['protein_id']:r for r in rows}
  with tempfile.TemporaryDirectory() as t:
   inv,_=prepare(annotation,domains,[pid],Path(t)/'inventory',coding_completeness=proof)
  hit=inv['gene_roles'][pid][0]
  self.assertEqual(hit['status'],'supported_role_homology')
  self.assertTrue(hit['coding_proof_used'])
  self.assertEqual(hit['source_partial'],'True')
  self.assertEqual(hit['coding_proof_used']['record_digest'],proof[pid]['record_digest'])

  for bad in ({}, {pid:dict(proof[pid],protein_sha256='0'*64)},
              {pid:dict(proof[pid],source_artifact=dict(proof[pid]['source_artifact'],path=str(case/'annotation'/'proteins.faa')))},
              {pid:dict(proof[pid],end=3736)},
              {pid:dict(proof[pid],terminal_stop_present=False)}) :
   with tempfile.TemporaryDirectory() as t:
    inv,_=prepare(annotation,domains,[pid],Path(t)/'inventory',coding_completeness=bad)
   hit=inv['gene_roles'][pid][0]
   self.assertEqual(hit['status'],'unresolved_role_homology')
   self.assertIn('partial_protein',hit['blockers'])
 def test_long_fusion_retains_complete_disjoint_roles_without_hydrogenase(self):
  case=ROOT/'results/round110_fusion_probe_v1';annotation=Path(json.loads((case/'workflow.json').read_text())['annotation'])
  with tempfile.TemporaryDirectory() as t:
   d=Path(t)/'domains';annotate_domains([annotation],d)
   inv,_=prepare(annotation,d,[],Path(t)/'inventory')
  pid,f=role(inv,'NuoF_flavin_FeS');epid,e=role(inv,'NuoE_Fer2')
  self.assertEqual(pid,epid);self.assertEqual(f['status'],e['status']);self.assertEqual(e['status'],'supported_role_homology')
  self.assertLess(e['best']['whole_query_coverage'],.3);self.assertGreater(e['best']['reference_coverage'],.8)
  self.assertFalse(set(f['best']['query_positions'])&set(e['best']['query_positions']))
  self.assertTrue(inv['observed_fusions']);self.assertEqual(inv['candidate_contexts'],{})
  self.assertEqual(inv['input_modules'][0]['status'],'supported_role_module_homology')
 def test_actual_hoxf_regions_keep_low_identity_e_visible_without_claiming_complete_fusion(self):
  inv=self.run_case('BBS_Hox2_A_R84');pid,f=role(inv,'NuoF_flavin_FeS');eid,e=role(inv,'NuoE_Fer2')
  self.assertEqual(pid,eid);self.assertEqual(f['status'],'supported_role_homology')
  self.assertEqual(e['status'],'unresolved_role_homology');self.assertIn('paired_identity_below_role_support_threshold',e['blockers'])
  self.assertTrue(e['best']['profiles']);self.assertLess(e['best']['paired_identity'],30)
  self.assertFalse(set(f['best']['query_positions'])&set(e['best']['query_positions']))
  self.assertFalse(any(x['protein_id']==pid for x in inv['observed_fusions']))
 def test_actual_shi_gamma_organization(self):
  inv=self.run_case('u834437');rs=[r for p,values in inv['gene_roles'].items() if p.endswith('_632') for r in values if r['role_class']=='SHI_gamma_NAD_FAD']
  self.assertEqual(len(rs),1);self.assertEqual(rs[0]['status'],'supported_role_homology')
  self.assertEqual({d['accession'].split('.')[0] for d in rs[0]['best']['profiles']},{'PF00970','PF00175','PF10418'})
 def test_nonhydrogenase_nfn_has_roles_but_no_function_claim(self):
  inv=self.run_case('9N5U_NfnABC_only')
  self.assertEqual(role(inv,'NuoF_flavin_FeS')[1]['status'],'supported_role_homology');self.assertEqual(role(inv,'NuoE_Fer2')[1]['status'],'supported_role_homology')
  self.assertNotIn('function_prediction',inv);self.assertEqual(inv['candidate_contexts'],{})
 def test_missing_ga_or_partial_withholds_but_retains_role(self):
  data=dataset('Amobile_A');base=analyse(**data);pid=role(base,'NuoF_flavin_FeS')[0]
  missing=copy.deepcopy(data);missing['domains']=[d for d in missing['domains'] if not(d['protein_id']==pid and d['accession'].split('.')[0]=='PF10589')]
  r=role(analyse(**missing),'NuoF_flavin_FeS')[1];self.assertNotEqual(r['status'],'supported_role_homology');self.assertIn('complete_owned_ordered_GA_domains_not_supported',r['blockers'])
  self.assertFalse(analyse(**missing)['input_modules'])
  partial=copy.deepcopy(data);next(g for g in partial['genes'] if g['id']==pid)['partial']='True'
  r=role(analyse(**partial),'NuoF_flavin_FeS')[1];self.assertIn('partial_protein',r['blockers'])
  self.assertFalse(analyse(**partial)['input_modules'])
 def test_invalid_model_query_or_ga_binding_is_rejected(self):
  data=dataset('Amobile_A');pid=role(analyse(**data),'NuoF_flavin_FeS')[0]
  for key,value in [('model_sha256','wrong'),('query_id','wrong'),('sequence_ga','0')]:
   with self.subTest(key=key):
    changed=copy.deepcopy(data);d=next(d for d in changed['domains'] if d['protein_id']==pid and d['accession'].split('.')[0]=='PF01512');d[key]=value
    with self.assertRaisesRegex(ValueError,'mismatch'):analyse(**changed)
 def test_domain_scores_and_order_are_not_replaced_by_pass_flag(self):
  data=dataset('Amobile_A');pid=role(analyse(**data),'NuoF_flavin_FeS')[0]
  for mutation in ('score','order'):
   changed=copy.deepcopy(data);d=next(d for d in changed['domains'] if d['protein_id']==pid and d['accession'].split('.')[0]=='PF10589')
   if mutation=='score':d['domain_bitscore']='0'
   else:d.update(query_start='1',query_end='85')
   self.assertNotEqual(role(analyse(**changed),'NuoF_flavin_FeS')[1]['status'],'supported_role_homology')
 def test_a_competing_reference_withdrawn_from_role_review_cannot_be_silently_ignored(self):
  data=dataset('Amobile_A');rid='FeFe_HytABCDE__HytB'
  data['owners']['NuoF_flavin_FeS'].remove(rid);data['compatible'][rid]=set()
  r=role(analyse(**data),'NuoF_flavin_FeS')[1]
  self.assertIn('overlapping_reference_role_competition',r['blockers']);self.assertTrue(any(x['reference_id']==rid for x in r['competition']))
 def test_opposite_strand_distant_and_cross_bin_relations(self):
  data=dataset('Amobile_A');base=analyse(**data);f=role(base,'NuoF_flavin_FeS')[0];e=role(base,'NuoE_Fer2')[0]
  data['candidates']=[e];data['mode']='assembly';fg=next(g for g in data['genes'] if g['id']==f);eg=next(g for g in data['genes'] if g['id']==e);fg['strand']=str(-int(eg['strand']))
  out=analyse(**data);ctx=out['candidate_contexts'][e]
  self.assertTrue(any(r['protein_id']==f and r['association_tier']=='same_contig_opposite_strand_local' for r in ctx['role_contexts']))
  self.assertFalse(out['input_modules'])
  far=copy.deepcopy(data);g=next(g for g in far['genes'] if g['id']==f);g['start']=str(int(g['start'])+100000);g['end']=str(int(g['end'])+100000)
  self.assertTrue(any(r['protein_id']==f and r['association_tier']=='same_contig_distant' for r in analyse(**far)['candidate_contexts'][e]['role_contexts']))
  binned=copy.deepcopy(data);binned['mode']='bin'
  for g in binned['genes']:g['bin_id']='B' if g['id']==f else 'A'
  out=analyse(**binned);ctx=out['candidate_contexts'][e]
  self.assertTrue(out['gene_roles'][f]);self.assertFalse(any(r['protein_id']==f for r in ctx['role_contexts']));self.assertGreater(ctx['excluded_role_count'],0)
  self.assertFalse(out['input_modules'])
  split_contig=dataset('Amobile_A');next(g for g in split_contig['genes'] if g['id']==f)['contig']='remote_contig'
  self.assertFalse(analyse(**split_contig)['input_modules'])
