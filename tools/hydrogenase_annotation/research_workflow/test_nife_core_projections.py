import copy,csv,json,unittest,sys,tempfile,shutil
from unittest.mock import patch
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from Bio import SeqIO
import nife_core_projections as ncp
from nife_core_projections import ROOT,_metric,_counter_search,_shared_large_roles,prepare_fallbacks
from mvh_hdr_function import _ag_profile,load_model,prepare_ag_core
from annotate import load_blast

class NiFeCoreProjectionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.root=Path(ROOT);cls.model=load_model(ROOT)
 def _case(self):
  a=self.root/'results/round83_Hadar_FNA_v1/u1686407/annotation';d=self.root/'results/round87_HDCR_release_v1/u1686407/domains/query_domains.tsv'
  genes=list(csv.DictReader((a/'genes.tsv').open(),delimiter='\t')); domains=list(csv.DictReader(d.open(),delimiter='\t')); anchors=[r for r in csv.DictReader((a/'catalytic_anchors.tsv').open(),delimiter='\t') if r['catalytic_family']=='NiFe']; prog={r['candidate_id']:r for r in csv.DictReader((self.root/'results/round87_HDCR_release_v1/u1686407/functions/candidate_progression.tsv').open(),delimiter='\t')}; ids={x['protein_id']:prog.get(x['protein_id'],{}).get('identity_status') for x in anchors}; seqs={r.id:str(r.seq) for r in SeqIO.parse((a/'proteins.faa').open(),'fasta')}; refs={r.id:str(r.seq) for r in SeqIO.parse((self.root/'references/processed/references.faa').open(),'fasta')}; raw=load_blast(a/'blast.tsv',dict(blast_evalue=self.model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0)); return a,genes,domains,anchors,ids,seqs,refs,raw
 def test_real_hadar_three_projection_cases(self):
  a,g,d,anc,ids,seqs,refs,raw=self._case(); out,inputs=prepare_ag_core(a,anc,g,d,ids,self.model); self.assertTrue(out);self.assertTrue(inputs); qs=[q for q in out if any(q.endswith('_'+x) for x in ('114','1386','1566'))];self.assertEqual(len(qs),3)
  for q in qs:
   r=out[q];self.assertEqual(r['status'],'sequence_supported_candidate');self.assertEqual(len(r['evidence']['hypotheses']),1);small=r['evidence']['hypotheses'][0][1];self.assertEqual(small['role'],'MvhG');self.assertTrue(small['reference_id'].startswith('NiFe_'));self.assertNotIn('Tonn',small['reference_id']);self.assertNotIn('Coo',small['reference_id']);self.assertNotIn('activity',r)
 def _current_real_case(self,genome,anchor_suffixes):
  a=self.root/f'results/round83_Hadar_FNA_v1/{genome}/annotation';dpath=self.root/f'results/round103_Hadar12_release_v1/{genome}/domains/query_domains.tsv';cat=self.root/f'results/round103_Hadar12_release_v1/{genome}/candidate_function_catalog.tsv'
  with (a/'genes.tsv').open() as f: genes=list(csv.DictReader(f,delimiter='\t'))
  with dpath.open() as f: domains=list(csv.DictReader(f,delimiter='\t'))
  with (a/'catalytic_anchors.tsv').open() as f: anchors=[r for r in csv.DictReader(f,delimiter='\t') if r['catalytic_family']=='NiFe' and any(r['protein_id'].endswith('_'+s) for s in anchor_suffixes)]
  with cat.open() as f: catalog={r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
  ids={x['protein_id']:catalog[x['protein_id']]['identity_status'] for x in anchors};return a,genes,domains,anchors,ids
 def test_current_R104_real_Hadar_core_regressions_are_core_only(self):
  for genome,suffix in [('u1202559','594'),('u1605023','465')]:
   with self.subTest(genome=genome):
    a,g,d,anc,ids=self._current_real_case(genome,[suffix]);out,_=prepare_ag_core(a,anc,g,d,ids,self.model);r=out[anc[0]['protein_id']];self.assertEqual(r['status'],'sequence_supported_candidate');self.assertEqual(len(r['evidence']['hypotheses']),1);self.assertEqual({c['role'] for c in r['evidence']['hypotheses'][0]},{'MvhA','MvhG'});self.assertEqual(r['evidence']['assembly_status'],'unresolved');self.assertNotIn('function_prediction',r);self.assertNotIn('reaction_prediction',r);self.assertNotIn('query_activity',r)
 def test_current_R104_real_boundaries_withdraw_core_without_relaxing_rules(self):
  a,g,d,anc,ids=self._current_real_case('u1605023',['465']);base,_=prepare_ag_core(a,anc,g,d,ids,self.model);qid=anc[0]['protein_id'];pid=base[qid]['evidence']['hypotheses'][0][1]['protein_id'];
  missing=[x for x in d if not(x['protein_id']==pid and x['accession'].split('.')[0]==self.model['ag_core']['small_profile'])];out,_=prepare_ag_core(a,anc,g,missing,ids,self.model);self.assertEqual(out[qid]['status'],'unresolved')
  lowga=copy.deepcopy(d);x=next(x for x in lowga if x['protein_id']==pid and x['accession'].split('.')[0]==self.model['ag_core']['small_profile']);x['sequence_bitscore']='0';out,_=prepare_ag_core(a,anc,g,lowga,ids,self.model);self.assertEqual(out[qid]['status'],'unresolved')
  pg=copy.deepcopy(g);next(x for x in pg if x['id']==pid)['partial']='True';out,_=prepare_ag_core(a,anc,pg,d,ids,self.model);self.assertEqual(out[qid]['status'],'unresolved')
  cross=copy.deepcopy(g);next(x for x in cross if x['id']==pid)['contig']='other_contig';out,_=prepare_ag_core(a,anc,cross,d,ids,self.model);self.assertEqual(out[qid]['status'],'unresolved')
  self.assertEqual(base[qid]['evidence']['hypotheses'][0][0]['query_domain_end']-base[qid]['evidence']['hypotheses'][0][0]['query_domain_start']+1,base[qid]['evidence']['hypotheses'][0][0]['resources'][qid][-1]-base[qid]['evidence']['hypotheses'][0][0]['resources'][qid][0]+1)
 def test_shared_large_certificate_is_strict_and_actual_reference_backed(self):
  with (self.root/'references/processed/references.faa').open() as f: refs={r.id:str(r.seq) for r in SeqIO.parse(f,'fasta')}
  allowed,inputs=_shared_large_roles(self.model,refs);self.assertTrue(allowed);self.assertTrue(inputs);self.assertTrue(all(r.startswith('NiFe_') for r in allowed));badrefs=dict(refs);badrefs[next(iter(allowed))]='X';self.assertRaisesRegex(ValueError,'sequence certificate mismatch',_shared_large_roles,self.model,badrefs)
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);(root/'references/processed').mkdir(parents=True);shutil.copyfile(self.root/'references/processed/references.json',root/'references/processed/references.json');cert=json.loads((self.root/'experiments/round104_hadar11/large_core_functional_equivalence.json').read_text());cert['rows'][0]['direct_ligand_positions']=cert['rows'][0]['direct_ligand_positions'][:3];(root/'cert.json').write_text(json.dumps(cert));m=copy.deepcopy(self.model);m['ag_core']['reference_projection']['shared_large_role_review_path']='cert.json';m['sources']['cert.json']=ncp.sha(root/'cert.json')
   with patch.object(ncp,'ROOT',root): self.assertRaisesRegex(ValueError,'ligand/structure certificate incomplete',_shared_large_roles,m,refs)
 def test_partial_or_missing_profile_withdraws_real_case(self):
  a,g,d,anc,ids,seqs,refs,raw=self._case();base,_=prepare_ag_core(a,anc,g,d,ids,self.model);q=next(q for q in base if q.endswith('_114') or q.endswith('_1386') or q.endswith('_1566'));pid=base[q]['evidence']['hypotheses'][0][1]['protein_id'];missing=[x for x in d if not(x['protein_id']==pid and x['accession'].split('.')[0]==self.model['ag_core']['small_profile'])];out,_=prepare_ag_core(a,anc,g,missing,ids,self.model);self.assertEqual(out[q]['status'],'unresolved');pg=copy.deepcopy(g);next(x for x in pg if x['id']==pid)['partial']='True';out,_=prepare_ag_core(a,anc,pg,d,ids,self.model);self.assertEqual(out[q]['status'],'unresolved')
 def test_ga_and_query_ownership_tampering_raises(self):
  a,g,d,anc,ids,seqs,refs,raw=self._case();base,_=prepare_ag_core(a,anc,g,d,ids,self.model);q=next(q for q in base if q.endswith('_114') or q.endswith('_1386') or q.endswith('_1566'));c=base[q]['evidence']['hypotheses'][0][1];mut=copy.deepcopy(d);x=next(x for x in mut if x['protein_id']==c['protein_id'] and x['accession'].split('.')[0]==self.model['ag_core']['small_profile']);x['model_sha256']='wrong';
  with self.assertRaisesRegex(ValueError,'model mismatch'):_ag_profile(c,mut,seqs,self.model)
  mut=copy.deepcopy(d);x=next(x for x in mut if x['protein_id']==c['protein_id'] and x['accession'].split('.')[0]==self.model['ag_core']['small_profile']);x['query_id']='wrong'
  with self.assertRaisesRegex(ValueError,'query mismatch'):_ag_profile(c,mut,seqs,self.model)
 def test_metric_real_self_and_malformed(self):
  sid='NiFe_FHL1_Hyc__HycG';s= self._case()[6][sid];h=dict(qseqid='anonymous',sseqid=sid,pident=100.,length=len(s),mismatch=0,gapopen=0,qstart=1,qend=len(s),sstart=1,send=len(s),evalue=0.,bitscore=float(len(s)),qlen=len(s),slen=len(s),qseq=s,sseq=s);self.assertEqual(_metric(h,{'anonymous':s},self._case()[6])['identity'],100.);self.assertRaises(ValueError,_metric,dict(h,qlen=len(s)+1),{'anonymous':s},self._case()[6]);self.assertRaises(ValueError,_metric,dict(h,bitscore=float('nan')),{'anonymous':s},self._case()[6])
 def test_anonymous_counter_self_search_404_residues(self):
  bank=self.root/'research_workflow/nife_small_counter_bank.faa';c={r.id:str(r.seq) for r in SeqIO.parse(bank.open(),'fasta')};self.assertEqual(sum(map(len,c.values())),404);q={'anonymous_a':c['counter_Ec_NuoB'],'anonymous_b':c['counter_Mmazei_FpoB']};m={'max_evalue':1e-10,'ag_core':{'reference_projection':{'counter_bank':'research_workflow/nife_small_counter_bank.faa'}}};hits,_,audit=_counter_search(q,c,m);self.assertEqual(audit['reference_residues'],404);self.assertIn('counter_Ec_NuoB',{x['raw_hsp']['sseqid'] for x in hits['anonymous_a']});self.assertIn('counter_Mmazei_FpoB',{x['raw_hsp']['sseqid'] for x in hits['anonymous_b']});self.assertTrue(any(x['identity']==100. and x['raw_hsp']['sseqid']=='counter_Ec_NuoB' for x in hits['anonymous_a']));self.assertTrue(any(x['identity']==100. and x['raw_hsp']['sseqid']=='counter_Mmazei_FpoB' for x in hits['anonymous_b']))
if __name__=='__main__':unittest.main()
