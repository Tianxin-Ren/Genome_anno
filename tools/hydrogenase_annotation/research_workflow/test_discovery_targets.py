import tempfile
import unittest
import json
from pathlib import Path
from discovery_evidence import load_targets, search_discovery
from Bio import SeqIO


class DiscoveryTargetTests(unittest.TestCase):
    def test_full_sequence_hit_is_retained_without_domain_threshold(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / 'x.tbl'
            p.write_text('query001 - NiFe - 2e-40 120.5 0.1 2e-40 119.0 0.0 1 1 0 0 1 1 1 1 1 -\n')
            rows = load_targets(p, {'NiFe': 'NiFe'})
        self.assertEqual(rows[0]['protein_id'], 'query001')
        self.assertEqual(rows[0]['family'], 'NiFe')
        self.assertAlmostEqual(rows[0]['sequence_bitscore'], 120.5)

    def test_invalid_full_sequence_rows_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / 'x.tbl'
            p.write_text('query001 - NiFe - 2e-40 120.5 0.1\n')
            with self.assertRaises(ValueError):
                load_targets(p, {'NiFe': 'NiFe'})
            p.write_text('query001 - Unknown - 2e-40 120.5 0.1 2e-40 119.0 0.0 1 1 0 0 1 1 1 1 1 -\n')
            with self.assertRaises(ValueError):
                load_targets(p, {'NiFe': 'NiFe'})
            p.write_text('query001 - NiFe - nan 120.5 0.1 2e-40 119.0 0.0 1 1 0 0 1 1 1 1 1 -\n')
            with self.assertRaises(ValueError):
                load_targets(p, {'NiFe': 'NiFe'})

    def test_real_hmmsearch_writes_full_sequence_table(self):
        root = Path(__file__).resolve().parents[1]
        source = root / 'references/processed/NiFe.faa'
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            query = out / 'query.faa'
            with source.open() as handle:
                record = next(SeqIO.parse(handle, 'fasta'))
            SeqIO.write([record], query, 'fasta')
            thresholds = json.loads((root / 'config/thresholds.json').read_text())
            search_discovery(query, out, root, 1, thresholds)
            self.assertTrue((out / 'family.tbl').is_file())
            self.assertTrue((out / 'family.domtbl').is_file())
            self.assertTrue((out / 'family.tbl').read_text().startswith('#'))
            targets = load_targets(out / 'family.tbl', {r['profile']: r['family'] for r in json.loads((root / 'references/processed/profile_manifest.json').read_text())})
            self.assertIn(record.id, {r['protein_id'] for r in targets})


if __name__ == '__main__':
    unittest.main()
