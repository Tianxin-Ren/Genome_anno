"""Verify that the main catalog preserves the bound, scoped reference detail."""
import csv,json,tempfile,unittest
from pathlib import Path
from candidate_screen import catalog,sha
ROOT=Path(__file__).resolve().parents[1]
def write(path,rows):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
class CatalogReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();p=Path(self.tmp.name);self.s=p/'screen';self.f=p/'functions';self.s.mkdir();self.f.mkdir();self.out=p/'catalog.tsv';self.a=p/'annotation';self.a.mkdir();src=self.a/'run_manifest.json';src.write_text('{}')
        write(self.s/'candidate_screen.tsv',[dict(protein_id='q',source_sets='["existing_anchor"]',discovery_status='existing_anchor',family_candidates='["FeFe"]',rejected_reason='')])
        (self.s/'candidate_screen_manifest.json').write_text(json.dumps(dict(annotation=str(self.a),output_sha256={'candidate_screen.tsv':sha(self.s/'candidate_screen.tsv')},input_sha256={str(src):sha(src)})))
        contracts=json.loads((ROOT/'config/reaction_contracts.json').read_text())['contracts'];c=next(c for c in contracts if c['contract_id']=='Ref_Mthermo_HydABCDEF_NADP_oxidation');self.cap={**c,'scope':'reference_assay_conditions_only'}
        self.row=dict(candidate_id='q',identity_status='calibrated_high',function_prediction='reference_sequence_matched_reaction_potential',coarse_function_prediction='candidate_H2_pyridine_nucleotide_coupling',confidence='exact_reference_sequence_potential',supported_attributes='["Reference_assay_reaction_capabilities"]',reference_reaction_status='exact_reference_sequence_potential',reference_contract_ids=json.dumps([c['contract_id']]),reference_reaction_keys=json.dumps([c['reaction_key']]),reference_assay_reactions=json.dumps([c['molecular_reaction']]),reference_assay_capabilities=json.dumps([self.cap]),NAD_NADP_selectivity='unresolved',in_vivo_direction='unknown',query_activity='not_tested')
        (self.f/'reference_reaction_evidence.json').write_text(json.dumps({'q':dict(status='exact_reference_sequence_potential',capabilities=[self.cap])}));self.save()
    def tearDown(self):self.tmp.cleanup()
    def save(self):
        write(self.f/'functional_annotations.tsv',[self.row]);(self.f/'manifest.json').write_text(json.dumps(dict(functional_annotations_sha256=sha(self.f/'functional_annotations.tsv'),reference_reaction_output_sha256=sha(self.f/'reference_reaction_evidence.json'),input_sha256={str(self.a/'run_manifest.json'):sha(self.a/'run_manifest.json')})))
    def test_reference_detail_reaches_main_table_without_changing_query_scope(self):
        rows=catalog(self.s,self.f,self.out);r=rows[0];self.assertEqual(json.loads(r['reference_reaction_keys']),['H2_NADP_reduction_without_required_Fd']);self.assertEqual(r['coarse_function_prediction'],'candidate_H2_pyridine_nucleotide_coupling');self.assertEqual(r['query_activity'],'not_tested');self.assertEqual(r['NAD_NADP_selectivity'],'unresolved')
        with self.out.open() as handle:final=next(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(final['function_prediction'],'reference_sequence_matched_reaction_potential');self.assertEqual(json.loads(final['reference_assay_capabilities'])[0]['electron_acceptors'],['NADP+'])
    def test_changed_reference_evidence_is_rejected(self):
        (self.f/'reference_reaction_evidence.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'evidence'):catalog(self.s,self.f,self.out)

    def add_class_rules(self):
        from function_rules import FIELDS
        directory=self.f/'function_rules';directory.mkdir()
        values=dict(function_class_ids=['Fe08'],function_class_prediction='H2 / NADP coupling',function_rule_status='reference_architecture_supported',function_rule_confidence='medium',function_rule_structure_confidence='high',function_rule_alternatives=[],function_rule_evidence='function_rules/evidence.json',function_rule_resolution='reference_architecture_prediction')
        projected={k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in values.items()}
        write(directory/'predictions.tsv',[dict(candidate_id='q',**projected)])
        (directory/'manifest.json').write_text(json.dumps(dict(status='complete',input_sha256={str(self.a/'run_manifest.json'):sha(self.a/'run_manifest.json')},output_sha256={'predictions.tsv':sha(directory/'predictions.tsv')})))
        self.row.update(projected);self.save()
        manifest=json.loads((self.f/'manifest.json').read_text());manifest['function_rules_manifest_sha256']=sha(directory/'manifest.json')
        (self.f/'manifest.json').write_text(json.dumps(manifest))

    def test_specific_class_rules_reach_main_catalog(self):
        self.add_class_rules();r=catalog(self.s,self.f,self.out)[0]
        self.assertEqual(r['function_class_ids'],['Fe08'])
        self.assertEqual(r['function_rule_confidence'],'medium')
        self.assertEqual(r['query_activity'],'not_tested')

    def test_tampered_class_rule_projection_is_rejected(self):
        self.add_class_rules();(self.f/'function_rules/predictions.tsv').write_text('changed')
        with self.assertRaisesRegex(ValueError,'rule output changed'):catalog(self.s,self.f,self.out)
    def test_generic_label_cannot_overwrite_exact_reference_detail(self):
        self.row['function_prediction']='candidate_H2_pyridine_nucleotide_coupling';self.save()
        with self.assertRaisesRegex(ValueError,'overwritten'):catalog(self.s,self.f,self.out)

    def test_coding_locus_potential_uses_distinct_main_label(self):
        self.row['reference_reaction_status']='exact_reference_coding_locus_potential'
        self.row['function_prediction']='reference_coding_locus_matched_reaction_potential'
        self.row['confidence']='exact_reference_coding_locus_potential'
        (self.f/'reference_reaction_evidence.json').write_text(json.dumps({'q':dict(status='exact_reference_coding_locus_potential',capabilities=[self.cap])}))
        self.save(); write(self.f/'functional_annotations.tsv',[self.row])
        final=catalog(self.s,self.f,self.out)
        self.assertEqual(final[0]['function_prediction'],'reference_coding_locus_matched_reaction_potential')
        self.assertEqual(final[0]['reference_reaction_status'],'exact_reference_coding_locus_potential')

    def test_coding_locus_potential_rejects_exact_sequence_label(self):
        self.row['reference_reaction_status']='exact_reference_coding_locus_potential'
        self.row['function_prediction']='reference_sequence_matched_reaction_potential'
        (self.f/'reference_reaction_evidence.json').write_text(json.dumps({'q':dict(status='exact_reference_coding_locus_potential',capabilities=[self.cap])}))
        self.save(); write(self.f/'functional_annotations.tsv',[self.row])
        with self.assertRaisesRegex(ValueError,'main label'):
            catalog(self.s,self.f,self.out)
    def test_detail_from_another_contract_cannot_be_spliced_into_main_table(self):
        self.row['reference_reaction_keys']='["invented_reaction"]';self.save()
        with self.assertRaisesRegex(ValueError,'differs'):catalog(self.s,self.f,self.out)
if __name__=='__main__':unittest.main()
