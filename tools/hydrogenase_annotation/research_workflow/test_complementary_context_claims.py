"""Partial scaffold evidence cannot become identity, reaction or selected members."""
import copy
import unittest
from claim_assertions import build_assertion


class ComplementaryContextClaims(unittest.TestCase):
    def setUp(self):
        self.record = dict(candidate_id='a', family='FeFe', identity_status='sequence_supported_uncalibrated_core',
            coarse_function_prediction='H2_core_supported_partner_unresolved', preferred_function='hydrogenase_core_candidate',
            preferred_members=[dict(role='candidate_anchor',protein_id='a')])
        self.members = [dict(role='catalytic_core',protein_id='a'),dict(role='complementary_structural_context',protein_id='b')]
        self.attr = dict(candidate_id='a', attribute='FeFe_complementary_scaffold_context', status='conditional_partial_context',
            value=dict(anchor_id='a',status='conditional',candidates=[dict(members=self.members)],selected_members=self.members),
            evidence_source='complementary_core_context_evidence.json', evidence_dependencies=['exact protein HSPs','local genomic linkage'],
            limitations='Partial scaffold only; assembly and activity unknown')

    def test_partial_hypothesis_and_deletion_leave_identity_and_function_unchanged(self):
        baseline = build_assertion(self.record)
        result = build_assertion(self.record,[self.attr])
        self.assertEqual(result['association_status'],'conditional')
        self.assertEqual(result['member_evidence'],baseline['member_evidence'])
        self.assertEqual(result['current_supported_annotation'],baseline['current_supported_annotation'])
        self.assertEqual(result['preferred_function_hypothesis'],baseline['preferred_function_hypothesis'])
        for key in ('catalytic_identity','electron_partner','coupling_mechanism','reaction_direction_conditions'):
            self.assertEqual(result['dimensions'][key],baseline['dimensions'][key])
        association=result['dimensions']['component_association']
        self.assertEqual(association['hypothesis_members'],self.members)
        self.assertEqual(association['hypothesis_scope'],'conditional_partial_structural_context')
        self.assertTrue(all(r['evidence_role']=='constraint' and r['evidence_group']=='association_and_context' for r in association['module_evidence']))
        withdrawn=build_assertion(self.record,[])
        self.assertNotIn('hypothesis',withdrawn['dimensions']['component_association'])

    def test_ambiguous_context_has_no_selected_hypothesis_members(self):
        attr=copy.deepcopy(self.attr);attr['status']='ambiguous_partial_context';attr['value']['status']='ambiguous'
        attr['value']['selected_members']=[]
        result=build_assertion(self.record,[attr])
        self.assertEqual(result['association_status'],'ambiguous')
        self.assertNotIn('hypothesis_members',result['dimensions']['component_association'])
        self.assertEqual(result['member_evidence'],self.record['preferred_members'])

    def test_supported_architecture_is_not_downgraded_or_replaced(self):
        architecture=dict(candidate_id='a',attribute='candidate_complex_architecture',status='complete_supported',value={})
        result=build_assertion(self.record,[architecture,self.attr])
        self.assertEqual(result['association_status'],'supported')
        self.assertNotIn('hypothesis',result['dimensions']['component_association'])
        self.assertEqual(len(result['dimensions']['component_association']['module_evidence']),1)

    def test_missing_or_wrong_owner_and_unsupported_identity_do_not_create_hypothesis(self):
        for change in (dict(anchor_id='wrong'),dict(candidates=[])):
            attr=copy.deepcopy(self.attr);attr['value'].update(change)
            self.assertNotIn('hypothesis',build_assertion(self.record,[attr])['dimensions']['component_association'])
        result=build_assertion(dict(self.record,identity_status='evidence_insufficient'),[self.attr])
        self.assertNotIn('hypothesis',result['dimensions']['component_association'])


if __name__=='__main__':unittest.main()
