import copy
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from hashlib import sha256
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
sys.path.insert(0,str(Path(__file__).resolve().parent))
from explain_candidates import identity_reasons,sequence_status,explain,table


class CandidateExplanationTests(unittest.TestCase):
    def setUp(self):
        self.parameters=dict(min_domain_coverage=.8,min_core_similarity=40,min_competitor_margin=5)
        self.feature=dict(protein_id='q',family='FeFe',eligible=False,eligibility_reasons=['reviewed_single_chain_center_or_span_not_supported'],
            core_domain_coverage=.96,positive_score=61,competitor_margin=10,positive_prototype='positive',competitor_prototype='negative',
            source_active_center=dict(active_center_status='ligands_substituted'))
        self.alt=dict(candidate_model_id='Aw_H_cubane_3C1Ser_literature',pattern_status='literature_ligand_pattern_candidate',query_partial='False')

    def test_noncanonical_sequence_support_does_not_change_calibrated_identity(self):
        original=copy.deepcopy(self.feature);anchor=dict(catalytic_identity_confidence='not_high_confidence')
        self.assertEqual(sequence_status(anchor,self.feature,self.alt,self.parameters),'sequence_supported_uncalibrated_noncanonical')
        self.assertEqual(self.feature,original)
        reasons=identity_reasons(self.feature,self.parameters)
        self.assertEqual(reasons[0]['type'],'observed_conflict')
        self.assertEqual(reasons[0]['claim'],'calibrated_high_identity')

    def test_unknown_partial_and_bare_substitution_do_not_gain_noncanonical_support(self):
        for alt in [{},dict(self.alt,query_partial='True'),dict(self.alt,query_partial=''),dict(self.alt,candidate_model_id='other')]:
            self.assertNotEqual(sequence_status({},self.feature,alt,self.parameters),'sequence_supported_uncalibrated_noncanonical')

    def test_joint_identity_reasons_are_not_reduced_to_one_gate(self):
        self.feature.update(core_domain_coverage=.2,positive_score=20,competitor_margin=-2)
        codes={r['code'] for r in identity_reasons(self.feature,self.parameters)}
        self.assertEqual(codes,{'reviewed_single_chain_center_or_span_not_supported','domain_coverage_below_rule','core_similarity_below_rule','nonpositive_competitor_margin'})

    def test_missing_competitor_detail_retains_margin_reason(self):
        self.feature.update(competitor_margin=-2,competitor_evidence=None)
        reason=next(r for r in identity_reasons(self.feature,self.parameters) if r['code']=='nonpositive_competitor_margin')
        self.assertEqual(reason['type'],'observed_conflict')
        self.assertEqual(reason['detail']['margin'],-2)

    def test_zero_margin_without_any_reference_is_missing_rule_not_conflict(self):
        feature=copy.deepcopy(self.feature)
        feature.update(positive_score=0,competitor_score=0,competitor_margin=0,
                       positive_prototype='',competitor_prototype='',
                       positive_evidence=None,competitor_evidence=None)
        before=copy.deepcopy(feature)
        anchor={}
        status=sequence_status(anchor,feature,{},self.parameters)
        reason=next(r for r in identity_reasons(feature,self.parameters)
                     if r['code']=='nonpositive_competitor_margin')
        self.assertEqual(reason['type'],'reference_rule_missing')
        self.assertIn('No comparable competitor reference',reason['detail']['scope'])
        self.assertEqual(status,'evidence_insufficient')
        self.assertEqual(feature,before)

    def test_indeterminate_negative_margin_is_truth_unavailable(self):
        self.feature.update(competitor_margin=-2,
            competitor_evidence={'label':'indeterminate','source':'prototype_review'})
        reason=next(r for r in identity_reasons(self.feature,self.parameters)
                     if r['code']=='nonpositive_competitor_margin')
        self.assertEqual(reason['type'],'truth_unavailable')
        self.assertEqual(reason['detail']['margin'],-2)
        self.assertEqual(reason['detail']['competitor_label'],'indeterminate')

    def test_domain_coverage_does_not_prove_sequence_fragmentation(self):
        self.feature['eligibility_reasons']=[];self.feature['core_domain_coverage']=.2
        self.assertEqual(identity_reasons(self.feature,self.parameters)[0]['type'],'outside_calibration_scope')
        self.feature['eligibility_reasons']=['partial_CDS']
        self.assertEqual(identity_reasons(self.feature,self.parameters)[0]['type'],'sequence_incomplete')

    def test_partial_profile_family_support_remains_uncalibrated(self):
        self.feature.update(core_domain_coverage=.33,eligibility_reasons=[],domain_evidence=[dict(accession='PF02906',profile='Fe_hyd_lg_C',model_coverage=.7,supported=True,query_core_overlap=.9,
            components=[dict(sequence_bitscore=60,sequence_ga=22,domain_bitscore=60,domain_ga=22)])])
        self.feature['source_active_center']['direct_ligand_status']='present'
        self.assertEqual(sequence_status({},self.feature,{},self.parameters),'sequence_supported_uncalibrated_core')
        self.assertEqual(sequence_status({'catalytic_identity_confidence':'calibrated_high'},self.feature,{},self.parameters),'calibrated_high')
        self.assertEqual(sequence_status({},self.feature,{},self.parameters,{'uncertain_positive_prototypes':[{'prototype_id':'positive'}]}),'hydrogenase_like_lineage_unconfirmed')

    def test_partial_attributes_survive_unresolved_center_and_reaction(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'source';source.mkdir();(source/'confirmation').mkdir()
            (source/'proteins.faa').write_text('>q\nMKK\n>partner\nMKK\n>unlinked\nMKK\n')
            (source/'genes.tsv').write_text(
                'id\tcontig\tstart\tend\tstrand\tpartial\n'
                'q\tq\t1\t9\t1\tFalse\n'
                'partner\tpartner\t1\t9\t1\tFalse\n'
                'unlinked\tunlinked\t1\t9\t1\tFalse\n')
            (source/'blast.tsv').write_text('')
            def save(name,rows):
                with (source/name).open('w') as handle:
                    writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
            save('catalytic_anchors.tsv',[dict(protein_id='q',catalytic_family='FeFe',catalytic_identity_confidence='not_high_confidence',active_center_status='ligands_substituted',
                missing_roles='no_reviewed_structure_rule',flags='family_candidate_only',structure_id='unresolved',structure_status='incomplete_or_ambiguous',structure_confidence='low',function_status='unassigned')])
            save('family_candidates.tsv',[dict(protein_id='q')])
            save('complex_members.tsv',[dict(anchor='q',module_category='electron_transfer_or_accessory',role='HydB',protein_id='partner',reference_protein='reference',query_domain_start='1',query_domain_end='400')])
            save('domain_hits.tsv',[dict(protein_id=q,accession='PF00175.27',profile='NAD_binding_1',passes='True',query_start='100',query_end='200',model_coverage='.9',independent_evalue='1e-20') for q in ['partner','unlinked']])
            (source/'confirmation/features.json').write_text(json.dumps([self.feature]))
            thresholds=json.loads((Path(__file__).resolve().parents[1]/'config/thresholds.json').read_text())
            (source/'run_manifest.json').write_text(json.dumps(dict(status='complete',input_mode='genome',thresholds=thresholds,input_files=[])))
            rows,attrs=explain(source,Path(directory)/'out',rules={'FeFe':self.parameters})
            self.assertIn('reference_rule_missing',rows[0]['missing_evidence_type'])
            self.assertNotIn('partner_not_detected',rows[0]['missing_evidence_type'])
            byname={r['attribute']:r for r in attrs}
            index=table(Path(directory)/'out'/'functional_annotations.tsv')[0]
            self.assertNotIn('candidate_complex_architecture',json.loads(index['supported_attributes']))
            self.assertEqual(index['confidence'],'unresolved')
            self.assertEqual(byname['electron_transfer_module_presence']['status'],'reference_like_modules_detected')
            self.assertEqual(byname['nucleotide_binding_domain_in_associated_module']['status'],'domain_profile_supported')
            self.assertEqual([r['protein_id'] for r in byname['nucleotide_binding_domain_in_associated_module']['value']],['partner'])
            for name in ['nucleotide_selectivity','ferredoxin_requirement','electron_coupling']:self.assertEqual(byname[name]['status'],'unresolved')
            self.assertEqual(table(source/'catalytic_anchors.tsv')[0]['function_status'],'unassigned')
            with self.assertRaises(FileExistsError):explain(source,Path(directory)/'out',rules={'FeFe':self.parameters})
            with self.assertRaises(ValueError):explain(source,Path(directory)/'absent',candidate_ids={'absent'},rules={'FeFe':self.parameters})

            manifest=json.loads((Path(directory)/'out/manifest.json').read_text())
            self.assertEqual(manifest['fdh_coding_hypotheses']['status'],'not_assessed_no_generated_Prodigal_GFF')
            (source/'prodigal.gff').write_text('# mock generated GFF; producer validation tested separately\n')
            def coding_evidence(annotation,output,model,root):
                output.mkdir()
                evidence=dict(status='DNA_candidates_evaluated',scope='conditional only',
                    candidates=[dict(status='coding_locus_supported_conditional')],rejected_proposals=[])
                (output/'coding_locus_evidence.json').write_text(json.dumps(evidence))
                (output/'manifest.json').write_text('{}')
                return evidence,{str(p):sha256(p.read_bytes()).hexdigest() for p in output.iterdir()}
            with patch('fdh_coding_locus.prepare',side_effect=coding_evidence) as producer:
                standalone=Path(directory)/'standalone'
                new_rows,new_attrs=explain(source,standalone,rules={'FeFe':self.parameters})
                producer.assert_called_once()
            self.assertEqual(new_rows,rows)
            self.assertEqual(new_attrs,attrs)
            manifest=json.loads((standalone/'manifest.json').read_text())
            summary=manifest['fdh_coding_hypotheses']
            self.assertEqual(summary['source'],'independent_coding_locus_search')
            self.assertEqual(summary['supported_conditional_count'],1)
            self.assertEqual(summary['sha256'],sha256((standalone/summary['path']).read_bytes()).hexdigest())
            self.assertIsNone(manifest['fhl_coding_locus_output_sha256'])
            self.assertIsNone(manifest['fhl_coding_association_output_sha256'])
            self.assertFalse((standalone/'fhl_function_evidence.json').exists())


if __name__=='__main__':unittest.main()
