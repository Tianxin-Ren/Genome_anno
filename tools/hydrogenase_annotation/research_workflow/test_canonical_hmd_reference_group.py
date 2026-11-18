import csv
import json
import unittest
from pathlib import Path

import canonical_hmd

ROOT = Path(__file__).resolve().parents[1]


def _inputs(annotation):
    features = json.loads((annotation / 'confirmation/features.json').read_text())
    with (annotation / 'genes.tsv').open() as handle:
        genes = list(csv.DictReader(handle, delimiter='\t'))
    return features, genes


class CanonicalHmdReferenceGroupTests(unittest.TestCase):
    def test_new_m_kandleri_group_positive_is_accepted(self):
        annotation = ROOT / 'results/round80_anonymous955_v1/annotation'
        features, genes = _inputs(annotation)
        row = next(f for f in features if f['protein_id'] == 'anon955_0843')
        result, _ = canonical_hmd.prepare(annotation, [row], genes, ROOT)
        self.assertEqual(result[row['protein_id']]['core_status'], 'reference_supported_uncalibrated')
        self.assertEqual(result[row['protein_id']]['reaction_status'], 'compatible_reaction_candidate_uncalibrated')
        self.assertEqual(result[row['protein_id']]['selected_canonical_positive_prototype'], 'prototype_e94f14b526651b1c783a')

    def test_original_native_hmd_remains_supported(self):
        annotation = ROOT / 'results/round75_HmdII_controls_v2/canonical_3F47/annotation'
        features, genes = _inputs(annotation)
        result, _ = canonical_hmd.prepare(annotation, features, genes, ROOT)
        self.assertEqual(result[features[0]['protein_id']]['core_status'], 'reference_supported_uncalibrated')

    def test_hmdii_controls_do_not_gain_canonical_reaction(self):
        for name in ('canonical_4YT4', 'native_A'):
            workflow = ROOT / 'results/round75_HmdII_controls_v2' / name / 'workflow.json'
            annotation = Path(json.loads(workflow.read_text())['annotation'])
            features, genes = _inputs(annotation)
            result, _ = canonical_hmd.prepare(annotation, features, genes, ROOT)
            self.assertTrue(all(r['reaction_status'] == 'unresolved' for r in result.values()), name)


if __name__ == '__main__':
    unittest.main()
