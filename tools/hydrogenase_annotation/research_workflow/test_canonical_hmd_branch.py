"""Actual remote Hmd diagnostic, with separate transfer and branch claims."""
import copy,csv,json,shutil,tempfile,unittest
from pathlib import Path
import canonical_hmd
from annotate import load_blast
from Bio import SeqIO
from explain_candidates import sequence_status
ROOT=Path(__file__).resolve().parents[1]
ANNOTATION=ROOT/'results/round111_Hmd_function_challenge_v1/Hmd_direct_assay_diagnostic/annotation'
def inputs(annotation=ANNOTATION):
 features=json.loads((annotation/'confirmation/features.json').read_text())
 with (annotation/'genes.tsv').open() as h:genes=list(csv.DictReader(h,delimiter='\t'))
 return features,genes
class CanonicalBranchTests(unittest.TestCase):
 def test_remote_experimental_Hmd_gets_branch_potential_without_80_percent_transfer(self):
  features,genes=inputs();result,bound=canonical_hmd.prepare(ANNOTATION,features,genes,ROOT);r=result[features[0]['protein_id']]
  self.assertEqual(r['core_blockers'],['full_canonical_reference_match_not_supported'])
  self.assertEqual(r['reaction_status'],'unresolved')
  self.assertEqual(r['branch_core_status'],'sequence_supported_candidate')
  self.assertEqual(r['branch_reaction_status'],'compatible_reaction_candidate_uncalibrated')
  self.assertGreater(r['branch_comparison']['margin'],30)
  self.assertEqual(len(r['branch_comparison']['canonical_reference_scores']),2)
  self.assertTrue(r['branch_comparison']['alternative_reference_scores'])
  self.assertIn(str(ANNOTATION/'confirmation/feature_manifest.json'),bound)
  self.assertEqual(sequence_status({},features[0],{'Fe_reference_identity':{'canonical':r}},{}),'sequence_supported_uncalibrated_core')
 def test_a_withdrawn_canonical_reference_is_not_ignored_as_an_alternative(self):
  features,_=inputs();model=json.loads((ROOT/'research_workflow/canonical_hmd_model.json').read_text())
  records=json.loads((ROOT/'references/processed/catalytic_confirmation/prototypes.json').read_text())
  with (ANNOTATION/'proteins.faa').open() as h:queries={r.id:str(r.seq) for r in SeqIO.parse(h,'fasta')}
  hits=load_blast(ANNOTATION/'confirmation/prototype_blast.tsv',dict(blast_evalue=model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
  r=canonical_hmd.branch_comparison(features[0],hits,queries,records,['prototype_e94f14b526651b1c783a'],model)
  self.assertIn('canonical_branch_does_not_beat_all_alternative_prototypes',r['blockers'])
  self.assertLess(r['margin'],0)
 def test_partial_core_stays_unresolved_even_with_concordant_references(self):
  features,genes=inputs();genes=copy.deepcopy(genes);genes[0]['partial']='True'
  r=canonical_hmd.prepare(ANNOTATION,features,genes,ROOT)[0][features[0]['protein_id']]
  self.assertEqual(r['branch_core_status'],'unresolved');self.assertIn('partial_or_unknown_CDS',r['branch_core_blockers'])
 def test_missing_competitor_rows_cannot_reuse_saved_search_binding(self):
  features,genes=inputs()
  with tempfile.TemporaryDirectory() as t:
   dest=Path(t)/'annotation';shutil.copytree(ANNOTATION,dest)
   p=dest/'confirmation/prototype_blast.tsv';lines=p.read_text().splitlines();p.write_text('\n'.join(l for l in lines if 'prototype_bd2495c6c1d36b1b3866' not in l)+'\n')
   r=canonical_hmd.prepare(dest,features,genes,ROOT)[0][features[0]['protein_id']]
   self.assertEqual(r['branch_reaction_status'],'unresolved');self.assertIn('complete_current_prototype_search_not_verified',r['branch_core_blockers'])
 def test_unknown_bank_binding_cannot_establish_branch_competition(self):
  features,genes=inputs()
  with tempfile.TemporaryDirectory() as t:
   dest=Path(t)/'annotation';shutil.copytree(ANNOTATION,dest);p=dest/'confirmation/feature_manifest.json';fm=json.loads(p.read_text());fm['model_manifest_sha256']='not_current';p.write_text(json.dumps(fm))
   r=canonical_hmd.prepare(dest,features,genes,ROOT)[0][features[0]['protein_id']]
   self.assertEqual(r['branch_core_status'],'unresolved')
 def test_HmdII_comparators_do_not_gain_canonical_branch_potential(self):
  for name in ('canonical_4YT4','native_A'):
   with self.subTest(name=name):
    workflow=ROOT/'results/round75_HmdII_controls_v2'/name/'workflow.json';a=Path(json.loads(workflow.read_text())['annotation']);features,genes=inputs(a)
    r=next(iter(canonical_hmd.prepare(a,features,genes,ROOT)[0].values()))
    self.assertEqual(r['branch_reaction_status'],'unresolved')
    self.assertTrue(set(r['branch_comparison']['blockers'])&{'complete_concordant_canonical_reference_group_not_supported','canonical_branch_does_not_beat_all_alternative_prototypes'})
if __name__=='__main__':unittest.main()
