import json, tempfile, unittest
from pathlib import Path
import fes_refinement


class TestFeSRefinement(unittest.TestCase):
    def test_model_is_pinned_and_core_unchanged(self):
        spec, sites = fes_refinement._check_model()
        self.assertEqual(len(sites), 36)
        self.assertEqual(spec['development_regression']['fixed24']['TP'], 15)
        self.assertEqual(sites[0:34], json.loads((fes_refinement.CORE_DIR/'sites.json').read_text()))

    def test_tampered_script_is_rejected(self):
        original = fes_refinement.SCRIPT
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / 'fes_sites.py'
            p.write_bytes(original.read_bytes() + b'\n# tamper\n')
            fes_refinement.SCRIPT = p
            with self.assertRaises(ValueError):
                fes_refinement._check_model()
        fes_refinement.SCRIPT = original

    def test_real_annotation_smoke(self):
        annotation = fes_refinement.ROOT/'results/hadar_pilot_round33_final/current/u1686407'
        with tempfile.TemporaryDirectory() as td:
            manifest = fes_refinement.refine(annotation, Path(td), threads=1)
            self.assertEqual(manifest['site_count'], 36)
            self.assertTrue((Path(td)/'fes_sites.tsv').exists())
            text = (Path(td)/'fes_sites.tsv').read_text()
            self.assertIn('unlinked_FeS_homolog_not_hydrogenase_identity', text)
            self.assertIn('unknown', text)


if __name__ == '__main__':
    unittest.main()
