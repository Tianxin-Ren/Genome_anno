import copy
import csv
import unittest
from pathlib import Path

from mbh_function import _predict, load_model, prepare


class Group4FunctionalProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parents[1]
        cls.model = load_model(cls.root)
        cls.annotation = cls.root / 'results/round83_Hadar_FNA_v1/u1202559/annotation'
        cls.domains_path = cls.root / 'results/round88_NiFe_core_release_v1/u1202559/domains/query_domains.tsv'
        with (cls.annotation / 'catalytic_anchors.tsv').open() as h:
            cls.anchors = list(csv.DictReader(h, delimiter='\t'))
        with cls.domains_path.open() as h:
            cls.domains = list(csv.DictReader(h, delimiter='\t'))
        with (cls.root / 'results/round88_NiFe_core_release_v1/u1202559/functions/candidate_progression.tsv').open() as h:
            cls.identities = {r['candidate_id']: r['identity_status'] for r in csv.DictReader(h, delimiter='\t')}
        cls.results, cls.inputs = prepare(cls.annotation, cls.anchors, cls.domains, cls.identities, cls.model, root=cls.root)
        cls.supported_id = next(q for q, r in cls.results.items() if r['support_basis'] == 'joint_JKLN_functional_role_projection')

    def supported_evidence(self):
        evidence = copy.deepcopy(self.results[self.supported_id]['evidence'])
        # The functional projection intentionally collapses reviewed equivalent
        # role audits. Keep the direct evaluator fixture at that projected state.
        for audit in evidence.get('audits', {}).values():
            audit['functional_status'] = 'no_qualified_rival'
        return evidence

    def test_real_u1202559_801_retains_coarse_fd_compatible_projection(self):
        result = self.results[self.supported_id]
        self.assertEqual(result['core_coupling_status'], 'reference_like_H2_Fd_coupling_potential')
        self.assertEqual(result['support_basis'], 'joint_JKLN_functional_role_projection')
        self.assertEqual(result['evidence']['scaffold_route']['status'], 'joint_architecture_supported')
        self.assertEqual(result['evidence']['functional_projection']['scope'].split(';')[0], 'Complete JKLN-like functional projection')

    def test_real_nonprojectable_cases_remain_unresolved(self):
        unresolved = [q for q, r in self.results.items() if r['core_coupling_status'] == 'unresolved']
        self.assertTrue(unresolved)
        for q in unresolved:
            self.assertIn('unresolved', self.results[q]['support_basis'])

    def test_projection_does_not_assign_direction_or_obligatory_fd(self):
        result = self.results[self.supported_id]
        limitations = result['limitations']
        self.assertIn('reaction direction', limitations)
        self.assertIn('Fd compatibility does not establish obligatory Fd use', limitations)
        self.assertNotIn('NAD specificity', result['support_basis'])

    def test_missing_n_blocks_core_and_coupling(self):
        evidence = self.supported_evidence(); evidence['roles']['MbhN'] = []
        result = _predict(evidence, self.model)
        self.assertEqual(result['h2_core_status'], 'sequence_supported_candidate')
        self.assertEqual(result['core_coupling_status'], 'unresolved')
        self.assertIn('missing_or_ambiguous_MbhN', result['fd_blockers'])

    def test_partial_n_blocks_relay_without_assigning_activity(self):
        evidence = self.supported_evidence(); evidence['roles']['MbhN'][0]['partial'] = True
        evidence['partial_members'] = True; evidence['core_partial_members'] = False; evidence['relay_partial_members'] = True
        result = _predict(evidence, self.model)
        self.assertEqual(result['fd_relay_status'], 'unresolved')
        self.assertEqual(result['core_coupling_status'], 'unresolved')
        self.assertIn('partial_relay_member_sequence', result['fd_blockers'])

    def test_shifted_k_scaffold_profile_is_not_owned(self):
        evidence = self.supported_evidence(); evidence['roles']['MbhK'][0]['domain_supported'] = False
        result = _predict(evidence, self.model)
        self.assertEqual(result['h2_core_status'], 'unresolved')
        self.assertIn('required_domain_missing_MbhK', result['blockers'])

    def test_cys_j_or_n_pattern_break_blocks_core(self):
        for role in ('MbhJ', 'MbhN'):
            evidence = self.supported_evidence(); evidence['roles'][role][0]['pattern_evidence']['status'] = 'pattern_gap_or_substitution'
            result = _predict(evidence, self.model)
            self.assertEqual(result['h2_core_status'], 'unresolved' if role == 'MbhJ' else 'sequence_supported_candidate')
            field = 'blockers' if role == 'MbhJ' else 'fd_blockers'
            self.assertIn('reviewed_residue_pattern_not_supported_' + role, result[field])

    def test_non_nife_competition_audit_blocks_role(self):
        evidence = self.supported_evidence(); evidence['audits']['MbhJ']['functional_status'] = 'qualified_rival'
        result = _predict(evidence, self.model)
        self.assertEqual(result['h2_core_status'], 'unresolved')
        self.assertIn('alternate_family_competition_unresolved_MbhJ', result['blockers'])

    def test_independent_catalytic_identity_is_required(self):
        evidence = self.supported_evidence(); evidence['identity_status'] = 'not_a_catalytic_anchor'
        result = _predict(evidence, self.model)
        self.assertEqual(result['h2_core_status'], 'unresolved')
        self.assertIn('catalytic_family_sequence_support_insufficient', result['blockers'])

    def test_shared_core_resource_blocks_but_unrelated_resource_does_not(self):
        evidence = self.supported_evidence(); evidence['resource_conflict'] = True; evidence['resource_conflict_roles'] = ['MbhJ']
        result = _predict(evidence, self.model)
        self.assertEqual(result['h2_core_status'], 'unresolved')
        self.assertIn('cross_anchor_shared_query_resource', result['blockers'])
        evidence = self.supported_evidence(); evidence['resource_conflict'] = True; evidence['resource_conflict_roles'] = ['HdrA']
        result = _predict(evidence, self.model)
        self.assertEqual(result['h2_core_status'], 'sequence_supported_candidate')


if __name__ == '__main__':
    unittest.main()
