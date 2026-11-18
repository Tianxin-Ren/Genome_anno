import csv
import copy
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import canonical_hmd
from coding_completeness import prepare as coding_evidence, supports as complete_cds_support

ROOT = Path(__file__).resolve().parents[1]
CONTROLS = ROOT / 'results/round75_HmdII_controls_v2'


def inputs(name):
    ann = Path(json.loads((CONTROLS / name / 'workflow.json').read_text())['annotation'])
    features = json.loads((ann / 'confirmation/features.json').read_text())
    with (ann / 'genes.tsv').open() as handle:
        genes = list(csv.DictReader(handle, delimiter='\t'))
    return ann, features, genes


class CanonicalHmdTests(unittest.TestCase):
    def test_stopless_complete_region_retains_potential_without_complete_cds(self):
        jobs = json.loads((ROOT / 'experiments/round124_prior_function_regression/run_status_v1.json').read_text())['jobs']
        for case, expected in [('Hmd_missing_stop', 'reference_supported_uncalibrated'),
                               ('Hmd_Cys172_to_Ala', 'unresolved')]:
            command = next(j for j in jobs if j['name'] == case)['commands'][0]['command']
            ann = Path(command[command.index('--annotation') + 1])
            features = json.loads((ann / 'confirmation/features.json').read_text())
            with (ann / 'genes.tsv').open() as handle:
                genes = list(csv.DictReader(handle, delimiter='\t'))
            with tempfile.TemporaryDirectory() as directory:
                coding, _ = coding_evidence(ann, Path(directory) / 'coding', root=ROOT)
                result, _ = canonical_hmd.prepare(ann, features, genes, ROOT, coding_completeness=coding)
            q, = result
            self.assertEqual(result[q]['core_status'], expected, case)
            if case == 'Hmd_missing_stop':
                self.assertEqual(result[q]['reaction_status'], 'compatible_reaction_candidate_uncalibrated')
                self.assertTrue(result[q]['coding_region_used'])
                self.assertFalse(result[q]['coding_completeness_used'])
                self.assertEqual(coding[q]['status'], 'incomplete')
                self.assertEqual(coding[q]['reasons'], ['terminal_stop_missing'])
                self.assertTrue(coding[q]['source_partial'])
                self.assertFalse(complete_cds_support(coding[q], q))
                without_source, _ = canonical_hmd.prepare(ann, features, genes, ROOT)
                self.assertEqual(without_source[q]['core_status'], 'unresolved')
                altered = copy.deepcopy(genes)
                altered[0]['start'] = '2'
                mismatched, _ = canonical_hmd.prepare(ann, features, altered, ROOT, coding_completeness=coding)
                self.assertEqual(mismatched[q]['core_status'], 'unresolved')
            else:
                self.assertEqual(result[q]['reaction_status'], 'unresolved')

    def test_real_3F47_native_core_and_reaction_are_distinct_claims(self):
        ann, features, genes = inputs('canonical_3F47')
        results, bound = canonical_hmd.prepare(ann, features, genes, ROOT)
        row = results[features[0]['protein_id']]
        self.assertEqual(row['core_status'], 'reference_supported_uncalibrated')
        self.assertEqual(row['reaction_status'], 'compatible_reaction_candidate_uncalibrated')
        self.assertTrue(row['remapped_sites'])
        self.assertGreater(len(bound), 10)
        self.assertEqual({a['attribute'] for a in canonical_hmd.attributes(row)}, {'canonical_Hmd_FeGP_core','FeGP_H4MPT_hydrogenation_potential','canonical_Hmd_branch_core','canonical_H4MPT_branch_reaction_potential'})

    def test_HmdII_controls_do_not_receive_canonical_reaction(self):
        for name in ['canonical_4YT4', 'native_A']:
            ann, features, genes = inputs(name)
            results, _ = canonical_hmd.prepare(ann, features, genes, ROOT)
            for row in results.values():
                self.assertEqual(row['reaction_status'], 'unresolved', name)

    def test_DTH_identity_remains_supported_by_its_own_rule(self):
        ann, features, genes = inputs('native_A')
        results, _ = canonical_hmd.identity_evidence(ann, features, genes, ROOT)
        row = results[features[0]['protein_id']]
        self.assertEqual(row['hmdii']['hmdii_status'], 'reference_supported_uncalibrated')
        self.assertEqual(row['canonical']['reaction_status'], 'unresolved')

    def test_unknown_and_partial_CDS_do_not_support_complete_core(self):
        ann, features, genes = inputs('canonical_3F47')
        for value in [None, 'True']:
            altered = copy.deepcopy(genes)
            altered[0]['partial'] = value
            result, _ = canonical_hmd.prepare(ann, features, altered, ROOT)
            self.assertEqual(result[features[0]['protein_id']]['core_status'], 'unresolved')

    def test_stale_positive_alignment_score_is_rejected(self):
        ann, features, genes = inputs('canonical_3F47')
        features[0]['positive_evidence']['bitscore'] -= 1
        with self.assertRaises(ValueError):
            canonical_hmd.prepare(ann, features, genes, ROOT)

    def test_changed_reviewed_source_is_rejected(self):
        with patch('canonical_hmd.sha', return_value='changed'):
            with self.assertRaises(ValueError):
                canonical_hmd.prepare(*inputs('canonical_3F47'), ROOT)

    def test_changed_raw_site_alignment_is_rejected(self):
        ann, features, genes = inputs('canonical_3F47')
        with tempfile.TemporaryDirectory() as directory:
            dest = Path(directory) / 'annotation'
            shutil.copytree(ann, dest)
            path = dest / 'active_center_blast.tsv'
            rows = path.read_text().splitlines()
            for i, row in enumerate(rows):
                cells = row.split('\t')
                if cells[1] == '3F47__entity_1':
                    cells[12] = 'X' + cells[12][1:]
                    rows[i] = '\t'.join(cells)
            path.write_text('\n'.join(rows) + '\n')
            with self.assertRaises(ValueError):
                canonical_hmd.prepare(dest, features, genes, ROOT)

    def test_missing_raw_site_alignment_cannot_reuse_saved_ligands(self):
        ann, features, genes = inputs('canonical_3F47')
        with tempfile.TemporaryDirectory() as directory:
            dest = Path(directory) / 'annotation'
            shutil.copytree(ann, dest)
            (dest / 'active_center_blast.tsv').unlink()
            result, _ = canonical_hmd.prepare(dest, features, genes, ROOT)
            row = result[features[0]['protein_id']]
            self.assertEqual(row['core_status'], 'unresolved')
            self.assertIn('raw_canonical_site_alignment_missing', row['core_blockers'])


if __name__ == '__main__':
    unittest.main()
