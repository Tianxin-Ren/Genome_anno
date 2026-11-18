import unittest
import json
import hashlib
import tempfile
from pathlib import Path

from candidate_screen import (CATALOG_FIELDS, PREFERRED_CATALOG_FIELDS,
                              PREDICTION_FIELDS, _preferred_prediction,
                              _validate_preferred_manifest, ROOT, _preferred_policy)


class PreferredCatalogProjectionTests(unittest.TestCase):
    def test_preferred_columns_are_visible_before_legacy_labels(self):
        self.assertEqual(CATALOG_FIELDS[0], 'candidate_id')
        self.assertEqual(CATALOG_FIELDS[1:1 + len(PREFERRED_CATALOG_FIELDS)], PREFERRED_CATALOG_FIELDS)
        self.assertLess(CATALOG_FIELDS.index('preferred_function'), CATALOG_FIELDS.index('function_prediction'))

    def test_shared_redox_or_dna_evidence_does_not_force_hydrogenase_function(self):
        row = {
            'candidate_id': 'DNA_1', 'candidate_kind': 'dna_locus',
            'identity_status': 'unconfirmed', 'function_prediction': 'potential_homolog_function_unresolved',
            'screen_families': ['NiFe_or_Nuo_related'],
        }
        prediction = _preferred_prediction(row, {'screen_status': 'DNA_locus_unresolved'})
        self.assertIn(prediction['prediction_confidence'], {'low', 'medium', 'high'})
        self.assertNotEqual(prediction['preferred_function'], 'candidate_H2_pyridine_nucleotide_coupling')
        self.assertTrue(prediction['preferred_function'])

    def test_bound_sidecar_tamper_is_rejected(self):
        row = {field: ([] if field.endswith(('support', 'members', 'alternatives', 'risks'))
                       else {} if field == 'experimental_unknowns' else 'x')
               for field in PREDICTION_FIELDS}
        row['candidate_id'] = 'p1'
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            sidecar = directory / 'preferred_predictions.json'
            sidecar.write_text(json.dumps({'p1': row}, ensure_ascii=False))
            digest = hashlib.sha256(sidecar.read_bytes()).hexdigest()
            manifest = {
                'preferred_prediction_output_sha256': digest,
                'preferred_prediction_policy_version': _preferred_policy.VERSION,
                'preferred_prediction_policy_sha256': hashlib.sha256(
                    (ROOT / 'research_workflow/preferred_prediction.py').read_bytes()).hexdigest(),
            }
            self.assertTrue(_validate_preferred_manifest(directory, manifest, [row]))
            sidecar.write_text(json.dumps({'p1': dict(row, preferred_function='tampered')}))
            with self.assertRaisesRegex(ValueError, 'output changed'):
                _validate_preferred_manifest(directory, manifest, [row])


if __name__ == '__main__':
    unittest.main()
