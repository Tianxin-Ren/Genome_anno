import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from research_workflow.maturation_observation import ROOT, audit


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MaturationObservationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.annotation = Path(self.temp.name) / 'annotation'
        self.annotation.mkdir()
        self.output = Path(self.temp.name) / 'out'
        self.model = ROOT / 'references/processed/maturation_domain_audit/maturation_domains.hmm'
        self.write('proteins.faa', '>q1\nMPEPTIDE\n>q2\nMSECOND\n')
        self.write('genes.tsv', 'id\tcontig\tbin_id\tpartial\tcoordinate_evidence\ttopology\nq1\tc1\tb1\tFalse\tgenomic_CDS\tcircular\nq2\tc2\tb2\tFalse\tgenomic_CDS\tcircular\n')
        self.write('domain_hits.tsv', 'protein_id\taccession\tpasses\n')
        self.write('maturation.tsv', 'protein_id\tmaturation_gene\n')
        self.write('maturation_domains/hmmsearch.log', f'# query HMM file: {self.model}\n# target sequence database: {self.annotation / "proteins.faa"}\n[ok]\n')
        self.pin()

    def write(self, name, text):
        path = self.annotation / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def pin(self, mode='genome'):
        manifest = dict(status='complete', input_mode=mode,
                        normalized_proteins_sha256=sha(self.annotation / 'proteins.faa'),
                        maturation_domain_policy_sha256=sha(ROOT / 'config/maturation_domain_policy.json'),
                        maturation_domain_model_sha256=sha(self.model),
                        output_sha256={name: sha(self.annotation / name) for name in ('genes.tsv', 'domain_hits.tsv', 'maturation.tsv')})
        self.write('run_manifest.json', json.dumps(manifest))

    def result(self):
        return audit(self.annotation, self.output)

    def test_search_scope_and_family_rule_are_separate(self):
        result = self.result()
        self.assertTrue(result['search_proof']['complete'])
        self.assertTrue(all(row['search_status'] == 'searched_not_detected' and row['searched_protein_count'] == 2 for row in result['rows']))
        self.assertEqual(result['family_rule_scope']['Fe_Hmd']['status'], 'not_applicable')

    def test_incomplete_log_cannot_prove_absence(self):
        path = self.annotation / 'maturation_domains/hmmsearch.log'
        path.write_text(path.read_text().replace('[ok]', ''))
        result = self.result()
        self.assertEqual(result['status'], 'not_evaluated')
        self.assertTrue(all(row['search_status'] == 'not_evaluated' for row in result['rows']))

    def test_missing_log_cannot_complete_observation(self):
        (self.annotation / 'maturation_domains/hmmsearch.log').unlink()
        result = self.result()
        self.assertEqual(result['status'], 'not_evaluated')
        self.assertFalse(result['search_proof']['complete'])

    def test_open_boundary_is_separate_from_completed_search(self):
        path = self.annotation / 'genes.tsv'
        path.write_text(path.read_text().replace('circular', 'linear'))
        self.pin()
        result = self.result()
        self.assertEqual(result['status'], 'complete')
        self.assertTrue(all(row['genomic_observability'] == 'open_contig_boundaries' for row in result['rows']))

    def test_bound_table_tamper_is_rejected(self):
        with (self.annotation / 'maturation.tsv').open('a') as handle:
            handle.write('q1\thypA\n')
        with self.assertRaisesRegex(ValueError, 'input changed: maturation.tsv'):
            self.result()

    def test_protein_tamper_is_rejected(self):
        self.write('proteins.faa', '>q1\nMOTHER\n')
        with self.assertRaisesRegex(ValueError, 'input changed: proteins.faa'):
            self.result()

    def test_bins_cannot_borrow_maturases(self):
        self.write('maturation.tsv', 'protein_id\tmaturation_gene\nq1\thypA\n')
        self.write('domain_hits.tsv', 'protein_id\taccession\tpasses\nq1\tPF01155.1\tTrue\n')
        self.pin(mode='bin')
        rows = [row for row in self.result()['rows'] if row['maturation_role'] == 'hypA']
        self.assertEqual([(row['scope_id'], row['search_status'], row['matched_protein_ids']) for row in rows],
                         [('bin:b1', 'detected', ['q1']), ('bin:b2', 'searched_not_detected', [])])
        self.pin(mode='genome')
        rows = [row for row in self.result()['rows'] if row['maturation_role'] == 'hypA']
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['searched_protein_count'], 2)

    def test_hypf_requires_the_policy_any_of_architecture(self):
        self.write('maturation.tsv', 'protein_id\tmaturation_gene\nq1\thypF\n')
        self.write('domain_hits.tsv', 'protein_id\taccession\tpasses\nq1\tPF07503\tTrue\n')
        self.pin()
        role = lambda: next(row for row in self.result()['rows'] if row['maturation_role'] == 'hypF')
        self.assertEqual(role()['evidence_status'], 'homology_only')
        with (self.annotation / 'domain_hits.tsv').open('a') as handle:
            handle.write('q1\tPF17788.1\tTrue\n')
        self.pin()
        self.assertEqual(role()['profile_architecture_protein_ids'], ['q1'])

    def test_coordinate_uncertainty_does_not_hide_a_completed_search(self):
        path = self.annotation / 'genes.tsv'
        path.write_text(path.read_text().replace('genomic_CDS', 'protein_only'))
        self.pin()
        self.assertTrue(all(row['genomic_observability'] == 'coordinates_unavailable' and row['search_status'] == 'searched_not_detected' for row in self.result()['rows']))

    def test_legacy_unpinned_table_is_not_an_absence(self):
        path = self.annotation / 'run_manifest.json'
        manifest = json.loads(path.read_text())
        manifest['output_sha256'].pop('domain_hits.tsv')
        path.write_text(json.dumps(manifest))
        result = self.result()
        self.assertEqual(result['status'], 'not_evaluated')
        self.assertTrue(all(row['search_status'] == 'not_evaluated' for row in result['rows']))

    def test_rehashed_unknown_protein_is_rejected(self):
        self.write('maturation.tsv', 'protein_id\tmaturation_gene\nwrong\thypA\n')
        self.pin()
        with self.assertRaisesRegex(ValueError, 'unknown protein'):
            self.result()


if __name__ == '__main__':
    unittest.main()
