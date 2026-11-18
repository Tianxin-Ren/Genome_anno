import copy
import csv
import hashlib
import json
import unittest
from pathlib import Path

from Bio import SeqIO

from mbh_function import _validate_hit
from mvh_hdr_function import (
    ROOT, _ag_profile, _ag_reference_audit, load_model, predict_ag_core,
    prepare_ag_core,
)
from annotate import load_blast


class MvhAgCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(ROOT)
        cls.model = load_model(ROOT)

    def _case(self, name):
        ann = self.root / 'results/round83_Hadar_FNA_v1' / name / 'annotation'
        domains_path = self.root / 'results/round87_HDCR_release_v1' / name / 'domains/query_domains.tsv'
        with (ann / 'genes.tsv').open() as h: genes = list(csv.DictReader(h, delimiter='\t'))
        with (ann / 'catalytic_anchors.tsv').open() as h:
            anchors = [r for r in csv.DictReader(h, delimiter='\t') if r['catalytic_family'] == 'NiFe']
        with domains_path.open() as h: domains = list(csv.DictReader(h, delimiter='\t'))
        progression = self.root / 'results/round87_HDCR_release_v1' / name / 'functions/candidate_progression.tsv'
        with progression.open() as h: progression_rows = {r['candidate_id']: r for r in csv.DictReader(h, delimiter='\t')}
        identities = {a['protein_id']: progression_rows.get(a['protein_id'], {}).get('identity_status') for a in anchors}
        return ann, genes, anchors, domains, identities

    def test_real_hadar_two_case_prepare_ag_core(self):
        observed = {}
        for name in ('u377048', 'u499975'):
            ann, genes, anchors, domains, identities = self._case(name)
            results, inputs = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
            self.assertTrue(results)
            self.assertTrue(inputs)
            observed[name] = [r['status'] for r in results.values()]
            self.assertTrue(any(status == 'sequence_supported_candidate' for status in observed[name]))
        self.assertEqual(observed['u377048'], ['sequence_supported_candidate'] * 3)
        # SHI/Frh catalytic-role equivalents resolve this coarse core, not its electron partner.
        self.assertEqual(observed['u499975'], ['sequence_supported_candidate'] * 2)

    def test_ag_core_does_not_require_mvh_d_or_hdr_partner(self):
        ann, genes, anchors, domains, identities = self._case('u377048')
        results, _ = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
        for result in results.values():
            self.assertEqual(result['status'], 'sequence_supported_candidate')
            self.assertNotIn('activity', result)
            self.assertIn('no standalone AG activity', result['limitations'])

    def test_missing_or_partial_mvhg_blocks_ag_core(self):
        ann, genes, anchors, domains, identities = self._case('u377048')
        baseline, _ = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
        target = next(q for q, r in baseline.items() if r['status'] == 'sequence_supported_candidate')
        g = baseline[target]['evidence']['hypotheses'][0][1]['protein_id']
        missing = [d for d in domains if d['protein_id'] != g]
        results, _ = prepare_ag_core(ann, anchors, genes, missing, identities, self.model)
        self.assertEqual(results[target]['status'], 'unresolved')
        self.assertIn('owned_complete_small_subunit_profile_missing', results[target]['blockers'])

        partial_genes = copy.deepcopy(genes)
        next(row for row in partial_genes if row['id'] == g)['partial'] = 'True'
        results, _ = prepare_ag_core(ann, anchors, partial_genes, domains, identities, self.model)
        self.assertEqual(results[target]['status'], 'unresolved')

    def test_small_profile_requires_real_model_ga_and_query_ownership(self):
        ann, genes, anchors, domains, identities = self._case('u377048')
        result, _ = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
        component = result[next(iter(result))]['evidence']['hypotheses'][0][1]
        with (ann / 'proteins.faa').open() as h: seqs = {r.id: str(r.seq) for r in SeqIO.parse(h, 'fasta')}
        bad_hash = copy.deepcopy(domains)
        d = next(d for d in bad_hash if d['protein_id'] == component['protein_id'] and d['accession'].split('.')[0] == self.model['ag_core']['small_profile'])
        d['model_sha256'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'model mismatch'):
            _ag_profile(component, bad_hash, seqs, self.model)
        bad_query = copy.deepcopy(domains); d = next(d for d in bad_query if d['protein_id'] == component['protein_id'] and d['accession'].split('.')[0] == self.model['ag_core']['small_profile']); d['query_id'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'query mismatch'):
            _ag_profile(component, bad_query, seqs, self.model)
        bad_ga = copy.deepcopy(domains); d = next(d for d in bad_ga if d['protein_id'] == component['protein_id'] and d['accession'].split('.')[0] == self.model['ag_core']['small_profile']); d['sequence_ga'] = '0'
        with self.assertRaisesRegex(ValueError, 'GA changed'):
            _ag_profile(component, bad_ga, seqs, self.model)
        bad_coverage = copy.deepcopy(domains); d = next(d for d in bad_coverage if d['protein_id'] == component['protein_id'] and d['accession'].split('.')[0] == self.model['ag_core']['small_profile']); d['model_coverage'] = 'nan'
        with self.assertRaises((ValueError, TypeError)):
            _ag_profile(component, bad_coverage, seqs, self.model)

    def test_small_profile_rejects_nan_ga_instead_of_accepting_it(self):
        ann, genes, anchors, domains, identities = self._case('u377048')
        result, _ = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
        component = result[next(iter(result))]['evidence']['hypotheses'][0][1]
        with (ann / 'proteins.faa').open() as h: seqs = {r.id: str(r.seq) for r in SeqIO.parse(h, 'fasta')}
        mutated = copy.deepcopy(domains); d = next(d for d in mutated if d['protein_id'] == component['protein_id'] and d['accession'].split('.')[0] == self.model['ag_core']['small_profile']); d['sequence_ga'] = 'nan'
        with self.assertRaises((ValueError, TypeError)):
            _ag_profile(component, mutated, seqs, self.model)

    def test_non_nife_reference_can_compete_without_being_called_mvh(self):
        ann, genes, anchors, domains, identities = self._case('u377048')
        result, _ = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
        target = next(r for r in result.values() if r['status'] == 'sequence_supported_candidate')
        component = target['evidence']['hypotheses'][0][0]
        with (ann / 'proteins.faa').open() as h: seqs = {r.id: str(r.seq) for r in SeqIO.parse(h, 'fasta')}
        with (self.root / 'references/processed/references.faa').open() as h: refs = {r.id: str(r.seq) for r in SeqIO.parse(h, 'fasta')}
        raw = load_blast(ann / 'blast.tsv', dict(blast_evalue=self.model['max_evalue'], min_identity=0, min_query_coverage=0, min_reference_coverage=0))
        # Engineering control: bind the real query sequence to an existing FeFe
        # reference ID through a valid synthetic HSP.  This exercises competition
        # without changing the reviewed competition ratio or model policy.
        rival_id = 'FeFe_TmHydABC__HydA'
        rival_len = min(len(seqs[component['protein_id']]), len(refs[rival_id]))
        refs = dict(refs); refs[rival_id] = seqs[component['protein_id']][:rival_len]
        engineered = dict(qseqid=component['protein_id'], sseqid=rival_id, pident=100.0,
                          length=rival_len, qstart=1, qend=rival_len, sstart=1, send=rival_len,
                          evalue=1e-40, bitscore=10000.0, qlen=len(seqs[component['protein_id']]),
                          slen=rival_len, qseq=seqs[component['protein_id']][:rival_len],
                          sseq=seqs[component['protein_id']][:rival_len])
        audit = _ag_reference_audit(component, raw + [engineered], seqs, refs,
                                    set(self.model['ag_core']['large_core_reference_ids']), self.model)
        self.assertEqual(audit['status'], 'unresolved')
        self.assertTrue(any(not x['reviewed_core_reference'] for x in audit['competing_other_function_references']))

    def test_malformed_real_blast_alignment_is_rejected_by_sequence_binding(self):
        ann, genes, anchors, domains, identities = self._case('u377048')
        result, _ = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
        target = next(r for r in result.values() if r['status'] == 'sequence_supported_candidate')
        component = target['evidence']['hypotheses'][0][0]
        with (ann / 'proteins.faa').open() as h: seqs = {r.id: str(r.seq) for r in SeqIO.parse(h, 'fasta')}
        with (self.root / 'references/processed/references.faa').open() as h: refs = {r.id: str(r.seq) for r in SeqIO.parse(h, 'fasta')}
        raw = load_blast(ann / 'blast.tsv', dict(blast_evalue=self.model['max_evalue'], min_identity=0, min_query_coverage=0, min_reference_coverage=0))
        hit = next(h for h in raw if h['qseqid'] == component['protein_id'] and h['sseqid'].endswith('__MvhA'))
        hit['qseq'] = ('A' if hit['qseq'][0] != 'A' else 'C') + hit['qseq'][1:]
        with self.assertRaisesRegex(ValueError, 'bound qseq'):
            _ag_reference_audit(component, raw, seqs, refs, set(self.model['ag_core']['large_core_reference_ids']), self.model)

    def test_core_resource_conflict_blocks_only_when_it_is_an_ag_overlap(self):
        ann, genes, anchors, domains, identities = self._case('u377048')
        result, _ = prepare_ag_core(ann, anchors, genes, domains, identities, self.model)
        target = next(r for r in result.values() if r['status'] == 'sequence_supported_candidate')
        blocked = copy.deepcopy(target['evidence']); blocked['resource_conflicts'] = [{'other_anchor': 'other', 'source': 'AG'}]
        self.assertEqual(predict_ag_core(blocked)['status'], 'unresolved')
        unrelated = copy.deepcopy(target['evidence']); unrelated['resource_conflicts'] = []
        unrelated['hdr_overlap'] = {'other_anchor': 'other', 'source': 'Hdr'}
        self.assertEqual(predict_ag_core(unrelated)['status'], 'sequence_supported_candidate')


if __name__ == '__main__':
    unittest.main()
