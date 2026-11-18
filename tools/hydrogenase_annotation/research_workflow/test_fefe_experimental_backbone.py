import copy
import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from fefe_experimental_backbone import ROOT, _digest, _evaluate, _load_model, prepare, record_supported, verify_saved
from identity_scope import FEFE_EXPERIMENTAL_BACKBONE_RULE, fefe_experimental_backbone_support, sequence_core_support


class FeFeExperimentalBackboneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(ROOT)
        cls.annotation = cls.root / 'results/round80_anonymous955_v1/annotation'
        cls.parameters = json.loads((cls.root / 'config/catalytic_confirmation_thresholds.json').read_text())['families']['FeFe']
        with (cls.annotation / 'confirmation/features.json').open() as h:
            cls.features = {f['protein_id']: f for f in json.load(h)}
        cls.expected = {'anon955_' + x for x in ('0091','0209','0350','0367','0552','0665','0692','0720','0748','0798')}

    def _prepare(self, features=None, annotation=None):
        tmp = tempfile.TemporaryDirectory()
        ann = annotation or self.annotation
        out = Path(tmp.name) / 'backbone'
        records, hashes = prepare(ann, features or self.features, out, self.parameters, root=self.root)
        return tmp, ann, out, records, hashes

    def test_model_and_reference_mapping_are_bound(self):
        model, refs, context, rules, mapper = _load_model(self.root)
        self.assertEqual(model['positive_reference'], 'TamHydS')
        self.assertEqual(model['negative_reference'], 'NARF_Q6DHP6')
        self.assertEqual(model['reference_positions'], [138, 193, 293, 297])
        self.assertEqual(model['template_positions'], [300, 355, 499, 503])
        self.assertEqual(model['profile_match_columns'], [77, 130, 272, 276])
        self.assertEqual(set(refs), {'TamHydS', 'NARF_Q6DHP6'})
        self.assertIn('min_window_coverage', context)
        self.assertEqual(rules, self.parameters)
        self.assertEqual(mapper['reference_id'], 'TamHydS')

    def test_real_anonymous955_recovers_expected_rescues(self):
        tmp, _, out, records, hashes = self._prepare()
        with tmp:
            supported = {q for q, r in records.items() if r['status'] == 'experimental_backbone_pattern_supported'}
            self.assertTrue(self.expected <= supported)
            self.assertEqual(len(supported), 24)
            self.assertEqual(len([q for q in supported if self.features[q]['source_active_center']['direct_ligand_status'] == 'not_covered']), 10)
            self.assertEqual(len([q for q in supported if self.features[q]['source_active_center']['direct_ligand_status'] == 'present']), 14)
            self.assertTrue(hashes)
            self.assertTrue((out / 'backbone_evidence.json').is_file())

    def test_prepare_does_not_mutate_features(self):
        snapshot = copy.deepcopy(self.features)
        tmp, _, _, _, _ = self._prepare()
        tmp.cleanup()
        self.assertEqual(self.features, snapshot)

    def test_empty_core_is_unresolved_but_malformed_nonempty_core_is_rejected(self):
        model = {'positive_reference': 'positive', 'negative_reference': 'negative', 'max_evalue': 1.0}
        hits = {'positive': {'qlen': 4, 'slen': 4, 'sseqid': 'positive', 'evalue': 0.0, 'bitscore': 1.0, 'pident': 100.0}}
        args = ('query', {'core_positions': []}, hits, 'AAAA', {}, model, {'positive': 'AAAA'}, {}, {}, {}, 'model')
        with patch('fefe_experimental_backbone._validate_hit'):
            self.assertIsNone(_evaluate(*args))
            with self.assertRaisesRegex(ValueError, 'query-core positions invalid'):
                _evaluate('query', {'core_positions': [0]}, hits, 'AAAA', {}, model, {'positive': 'AAAA'}, {}, {}, {}, 'model')
            for malformed in (None, False, ''):
                with self.assertRaisesRegex(ValueError, 'query-core positions invalid'):
                    _evaluate('query', {'core_positions': malformed}, hits, 'AAAA', {}, model, {'positive': 'AAAA'}, {}, {}, {}, 'model')

    def test_saved_record_verifies_against_actual_inputs(self):
        tmp, ann, out, records, _ = self._prepare()
        with tmp:
            checked, manifest = verify_saved(out, annotation=ann, features=self.features,
                                             parameters=self.parameters, root=self.root)
            self.assertEqual(checked, records)
            self.assertEqual(manifest['status'], 'complete')

    def test_raw_alignment_tamper_is_rejected_after_output_hash_repair(self):
        tmp, ann, out, _, _ = self._prepare()
        with tmp:
            blast = out / 'blast.tsv'; lines = blast.read_text().splitlines()
            fields = lines[0].split('\t'); fields[9] = str(float(fields[9]) + 1.0); lines[0] = '\t'.join(fields); blast.write_text('\n'.join(lines) + '\n')
            manifest_path = out / 'manifest.json'; manifest = json.loads(manifest_path.read_text())
            manifest['output_sha256']['blast.tsv'] = hashlib.sha256(blast.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
            with self.assertRaises(ValueError):
                verify_saved(out, annotation=ann, features=self.features, parameters=self.parameters, root=self.root)

    def test_feature_source_tamper_is_rejected(self):
        tmp, ann, out, _, _ = self._prepare()
        with tmp:
            tampered = copy.deepcopy(self.features)
            pid = next(iter(self.expected & set(tampered))); tampered[pid]['core_positions'] = list(tampered[pid].get('core_positions', [])) + [1]
            with self.assertRaisesRegex(ValueError, 'feature source'):
                prepare(ann, tampered, Path(tmp.name) / 'second', self.parameters, root=self.root)

    def test_partial_cds_is_excluded_using_real_copied_annotation(self):
        tmp = tempfile.TemporaryDirectory(); ann = Path(tmp.name) / 'annotation'; shutil.copytree(self.annotation, ann)
        genes = ann / 'genes.tsv'; rows = list(csv.DictReader(genes.open(), delimiter='\t'))
        pid = 'anon955_0209'; row = next(r for r in rows if r['id'] == pid); row['partial'] = 'True'
        with genes.open('w', newline='') as h:
            writer = csv.DictWriter(h, fieldnames=rows[0].keys(), delimiter='\t'); writer.writeheader(); writer.writerows(rows)
        out = Path(tmp.name) / 'out'; records, _ = prepare(ann, self.features, out, self.parameters, root=self.root)
        self.assertNotIn(pid, records)
        tmp.cleanup()

    def test_record_support_rejects_wrong_center_context_and_model_binding(self):
        tmp, _, _, records, _ = self._prepare()
        with tmp:
            pid = next(iter(self.expected & set(records))); record = records[pid]
            wrong = copy.deepcopy(record); wrong['backbone_span_status'] = 'unsupported'
            self.assertFalse(record_supported(wrong, self.features[pid], self.parameters, root=self.root))
            wrong = copy.deepcopy(record); wrong['model_sha256'] = 'wrong'
            self.assertFalse(record_supported(wrong, self.features[pid], self.parameters, root=self.root))
            wrong = copy.deepcopy(record); wrong.pop('positive_reference_sha256', None)
            self.assertFalse(record_supported(wrong, self.features[pid], self.parameters, root=self.root))

    def test_record_support_rejects_narf_competition_and_missing_sites(self):
        tmp, _, _, records, _ = self._prepare()
        with tmp:
            pid = next(iter(self.expected & set(records))); record = records[pid]
            negative = copy.deepcopy(record); negative['NARF_countercomparison']['similarity'] = negative['positive_comparison']['similarity'] + 10
            self.assertFalse(record_supported(negative, self.features[pid], self.parameters, root=self.root))
            missing = copy.deepcopy(record); missing['sites'] = missing['sites'][:-1]
            self.assertFalse(record_supported(missing, self.features[pid], self.parameters, root=self.root))

    def test_duplicate_fasta_id_is_rejected_on_eligible_path(self):
        tmp = tempfile.TemporaryDirectory(); ann = Path(tmp.name) / 'annotation'; shutil.copytree(self.annotation, ann)
        fasta = ann / 'proteins.faa'
        with open(fasta, 'a') as h: h.write('>anon955_0209\nAAAA\n')
        with self.assertRaisesRegex(ValueError, 'ownership'):
            prepare(ann, self.features, Path(tmp.name) / 'out', self.parameters, root=self.root)
        tmp.cleanup()

    def test_identity_scope_backbone_gate_respects_score_and_margin_boundaries(self):
        tmp, _, _, records, _ = self._prepare()
        with tmp:
            pid = next(iter(self.expected & set(records))); feature = copy.deepcopy(self.features[pid]); alt = {'FeFe_experimental_backbone': records[pid], 'query_partial': False}
            for field, value in [('positive_score', 39.99), ('competitor_margin', 4.99)]:
                changed = copy.deepcopy(feature); changed[field] = value
                rebound = copy.deepcopy(records[pid]); rebound['feature_sha256'] = _digest(changed)
                self.assertFalse(fefe_experimental_backbone_support(changed, dict(alt, FeFe_experimental_backbone=rebound), self.parameters))

    def test_identity_scope_backbone_gate_rejects_old_center_conflicts_and_partial_alt(self):
        tmp, _, _, records, _ = self._prepare()
        with tmp:
            pid = next(iter(self.expected & set(records))); base = copy.deepcopy(self.features[pid])
            for status in ['substituted', 'ambiguous']:
                changed = copy.deepcopy(base); changed['source_active_center']['direct_ligand_status'] = status
                rebound = copy.deepcopy(records[pid]); rebound['feature_sha256'] = _digest(changed)
                self.assertFalse(fefe_experimental_backbone_support(changed, {'FeFe_experimental_backbone': rebound, 'query_partial': False}, self.parameters))
            present = copy.deepcopy(base); present['source_active_center']['direct_ligand_status'] = 'present'
            rebound = copy.deepcopy(records[pid]); rebound['feature_sha256'] = _digest(present)
            self.assertTrue(fefe_experimental_backbone_support(present, {'FeFe_experimental_backbone': rebound, 'query_partial': False}, self.parameters))
            rebound = copy.deepcopy(records[pid]); rebound['feature_sha256'] = _digest(base)
            self.assertFalse(fefe_experimental_backbone_support(base, {'FeFe_experimental_backbone': rebound, 'query_partial': True}, self.parameters))
            self.assertFalse(sequence_core_support(base, {'FeFe_experimental_backbone': rebound, 'query_partial': False}, self.parameters, keep_center_context=True))

    def test_record_gate_rejects_nan_countercomparison_score(self):
        tmp, _, _, records, _ = self._prepare()
        with tmp:
            pid = next(iter(self.expected & set(records))); record = copy.deepcopy(records[pid]); record['NARF_countercomparison']['similarity'] = float('nan')
            self.assertFalse(record_supported(record, self.features[pid], self.parameters, root=self.root))

    def test_real_tamhyds_positive_annotation_is_supported_without_identity_upgrade(self):
        ann = self.root / 'results/round90_FeFe_backbone_release_v1/BB_q0008_TamHydS_native_positive/annotation'
        with (ann / 'confirmation/features.json').open() as h: feature = json.load(h)[0]
        tmp = tempfile.TemporaryDirectory(); records, _ = prepare(ann, {feature['protein_id']: feature}, Path(tmp.name) / 'out', self.parameters, root=self.root)
        record = records[feature['protein_id']]
        self.assertEqual(record['status'], 'experimental_backbone_pattern_supported')
        self.assertEqual(feature['source_active_center']['direct_ligand_status'], 'present')
        self.assertTrue(record_supported(record, feature, self.parameters, root=self.root))
        tmp.cleanup()


if __name__ == '__main__':
    unittest.main()
