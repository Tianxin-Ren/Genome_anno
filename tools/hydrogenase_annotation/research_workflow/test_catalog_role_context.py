"""The broad catalog consumes bound role context without promoting identity."""
import csv,json,shutil,tempfile,unittest
from pathlib import Path
from candidate_screen import catalog,sha,ROOT

class CatalogRoleContextTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.out=Path(self.tmp.name)
  self.prior=ROOT/'results/round109_general_tool_pilot_v5/Chlamydomonas_NDUFS2_sequence'
  with (self.prior/'candidate_function_catalog.tsv').open() as h:self.original=next(csv.DictReader(h,delimiter='\t'))
  self.functions=self.out/'functions';self.functions.mkdir();(self.functions/'independent_roles').mkdir()
  shutil.copyfile(self.prior/'functions/functional_annotations.tsv',self.functions/'functional_annotations.tsv')
  self.manifest=json.loads((self.prior/'functions/manifest.json').read_text())
  source=json.loads((self.prior/'candidates/candidate_screen_manifest.json').read_text())
  key=str(Path(source['annotation'])/'run_manifest.json');self.source_key=key
  self.roles={'candidate_contexts':{self.original['candidate_id']:{'role_contexts':[],'excluded_role_count':0,'scope':'reviewed role context only'}},'scope':'role evidence only','input_sha256':{key:source['input_sha256'][key]}}
  self._save()
 def _save(self):
  p=self.functions/'independent_roles/inventory.json';p.write_text(json.dumps(self.roles));self.manifest['module_role_inventory_sha256']=sha(p);(self.functions/'manifest.json').write_text(json.dumps(self.manifest))
 def run_catalog(self):return catalog(self.prior/'candidates',self.functions,self.out/'catalog.tsv',self.prior/'candidate_confirmation')
 def test_unconfirmed_candidate_gets_context_without_promotion(self):
  r=self.run_catalog()[0]
  for key in ['identity_status','function_prediction','function_confidence','supported_attributes']:self.assertEqual(r[key],self.original[key])
  self.assertEqual(r['electron_role_context'],self.roles['candidate_contexts'][r['candidate_id']])
  self.assertIn('functions/independent_roles/inventory.json',r['evidence'])
 def test_tampered_inventory_fails_before_catalog_write(self):
  (self.functions/'independent_roles/inventory.json').write_text('{}')
  with self.assertRaisesRegex(ValueError,'changed or missing'):self.run_catalog()
  self.assertFalse((self.out/'catalog.tsv').exists())
 def test_same_candidate_id_from_another_annotation_is_rejected(self):
  self.roles['input_sha256'][self.source_key]='wrong';self._save()
  with self.assertRaisesRegex(ValueError,'another annotation'):self.run_catalog()
