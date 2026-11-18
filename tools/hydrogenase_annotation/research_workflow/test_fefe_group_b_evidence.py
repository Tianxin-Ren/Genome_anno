import json
import unittest
from pathlib import Path

from fefe_group_b_evidence import MODEL, RULE, _site_positions, record_supported, sha


class AssayedLineageTests(unittest.TestCase):
    def test_gaps_and_substitutions_are_distinct_from_coordinate_support(self):
        hit = {'qstart': 10, 'sstart': 20, 'qseq': 'C-ACT', 'sseq': 'CCACT'}
        self.assertEqual(_site_positions(hit, [20, 21, 23, 24]), [10, None, 12, None])

    def test_both_assayed_references_and_negative_margin_are_required(self):
        path = Path(__file__).resolve().parents[1]/MODEL
        model = json.loads(path.read_text())
        names = list(model['assayed_references'])
        positions = [10, 60, 170, 174]
        evidence = {name: dict(query_ligand_positions=positions, core_coverage=.95,
                               similarity=48, evalue=1e-42) for name in names}
        proof = dict(rule=RULE, model_sha256=sha(path), references=evidence,
                     best_known_negative_similarity=25, negative_prototype_count=20)
        profile = dict(blockers=['low_profile_posterior'], assayed_lineage_corroboration=proof,
                       sites=[dict(pp=pp, query_position=pos, within_domain=True, residue_status='present')
                              for pp, pos in zip('9***', positions)])
        feature = dict(source_active_center=dict(direct_ligand_status='ambiguous', center_organization='single_polypeptide'))
        parameters = dict(min_domain_coverage=.8, min_core_similarity=40, min_competitor_margin=5)
        self.assertTrue(record_supported(feature, profile, parameters))
        proof['best_known_negative_similarity'] = 44
        self.assertFalse(record_supported(feature, profile, parameters))
        proof['best_known_negative_similarity'] = 25
        evidence.pop(names[0])
        self.assertFalse(record_supported(feature, profile, parameters))
