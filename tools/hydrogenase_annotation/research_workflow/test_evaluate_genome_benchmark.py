import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from evaluate_genome_benchmark import (ROOT, PROTECTED, metric, score_metrics,
                                       score_unit, validate_manifest, sequence_hash,
                                       is_specific_function)


class GenomeBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.seq = "M" + "A" * 20
        self.unit = {
            "unit_id": "positive",
            "biological_source_id": "source-positive",
            "independence_stratum": "heldout",
            "input": {"type": "proteins", "mode": "transcript"},
            "truth": {"hydrogenase": True, "target_AA_sha256": sequence_hash(self.seq),
                      "function_evaluable": True, "compatible_function_predictions": ["candidate_H2_F420_reduction"]},
        }

    def _output(self, seq=None, label="candidate_H2_F420_reduction", identity="calibrated_high"):
        out = self.root / "out"; (out / "annotation").mkdir(parents=True)
        sequence = self.seq if seq is None else seq
        with (out / "annotation/proteins.faa").open("w") as h:
            h.write(">target\n" + sequence + "\n")
        with (out / "candidate_function_catalog.tsv").open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=["candidate_id", "identity_status", "function_prediction"], delimiter="\t")
            w.writeheader(); w.writerow({"candidate_id": "target", "identity_status": identity, "function_prediction": label})
        return out

    def test_positive_and_function_correct(self):
        row = score_unit(self.unit, self._output())
        self.assertTrue(row["retained"]); self.assertTrue(row["identity_supported"]); self.assertTrue(row["function_correct"])
        self.assertEqual(score_metrics([row])["function_evaluable_positive"]["recall"], 1.0)

    def test_missing_target_is_discovery_fn(self):
        row = score_unit(self.unit, self._output(seq="M" + "C" * 20))
        self.assertFalse(row["retained"]); self.assertEqual(row["diagnostic"], "target_AA_not_found_in_called_proteome")
        self.assertEqual(score_metrics([row])["discovery"]["FN"], 1)

    def test_wrong_function_is_not_correct(self):
        row = score_unit(self.unit, self._output(label="candidate_H2_quinone_coupling"))
        self.assertTrue(row["function_assigned"]); self.assertFalse(row["function_correct"])
        metrics = score_metrics([row])["function_evaluable_positive"]
        self.assertEqual(metrics["FP"], 1)
        self.assertEqual(metrics["FN"], 1)
        self.assertIsNotNone(metrics["precision_wilson95"])
        self.assertIsNotNone(metrics["recall_wilson95"])

    def test_negative_specific_function_is_false_positive(self):
        negative = dict(self.unit, unit_id="negative", biological_source_id="source-negative",
                        truth={"hydrogenase": False, "target_AA_sha256": sequence_hash(self.seq), "function_evaluable": False, "compatible_function_predictions": []})
        row = score_unit(negative, self._output())
        self.assertTrue(row["function_assigned"])
        self.assertEqual(score_metrics([row])["negative_function"]["FP"], 1)

    def test_no_zero_denominator_success(self):
        m = metric(0, 0, 0, 0)
        self.assertIsNone(m["precision"]); self.assertIsNone(m["recall"]); self.assertIsNone(m["precision_wilson95"])
        self.assertIsNone(score_metrics([])["function_end_to_end_recall"])

    def test_partial_truth_does_not_score_other_candidates(self):
        # The scorer uses only the reviewed target hash; an unlisted candidate
        # in the catalog cannot become a negative truth observation.
        out = self._output();
        with (out / "candidate_function_catalog.tsv").open("a") as h:
            h.write("other\tcalibrated_high\tcandidate_H2_F420_reduction\n")
        row = score_unit(self.unit, out)
        self.assertEqual(row["target_protein_id"], "target")
        self.assertEqual(len(score_metrics([row])["discovery"]), 9)
        self.assertEqual(row['unscored_catalog_candidates'], 1)

    def test_unretained_target_does_not_remove_an_unscored_candidate(self):
        out = self._output()
        with (out/'candidate_function_catalog.tsv').open('w') as h:
            h.write('candidate_id\tidentity_status\tfunction_prediction\nother\tunresolved\thydrogenase_candidate_identity_unresolved\n')
        row = score_unit(self.unit, out)
        self.assertFalse(row['retained'])
        self.assertEqual(row['unscored_catalog_candidates'], 1)

    def test_incomplete_function_truth_keeps_assignment_coverage_separate(self):
        unit = dict(self.unit, truth=dict(self.unit['truth'], function_evaluable=False, compatible_function_predictions=[]))
        row = score_unit(unit, self._output())
        metrics = score_metrics([row])
        self.assertEqual(metrics['automatic_assignment_coverage_all_positive'], 1.0)
        self.assertEqual(metrics['function_evaluable_positive_count'], 0)
        self.assertIsNone(metrics['function']['precision'])

    def test_membrane_and_recoding_potentials_are_not_specific_partner_assignments(self):
        from evaluate_genome_benchmark import is_specific_function
        self.assertFalse(is_specific_function('candidate_membrane_associated_H2_redox_system'))
        self.assertFalse(is_specific_function('candidate_H2_formate_pathway_conditional_on_recoding'))
        self.assertTrue(is_specific_function('candidate_H2_ferredoxin_compatible_core'))

    def test_native_locus_mapping_keeps_start_boundary_differences(self):
        out = self._output(seq='M'+'C'*20)
        dna=self.root/'native.fna';dna.write_text('>native_chr\n'+'ATG'*110+'\n')
        truth=dict(self.unit['truth'],target_locus=dict(contig_sequence_sha256=sequence_hash('ATG'*110),start=1,end=300,strand=1))
        unit=dict(self.unit,input=dict(type='genome',path=str(dna),mode='genome'),truth=truth)
        genes=out/'annotation/genes.tsv'
        genes.write_text('id\tcontig\tstart\tend\tstrand\ntarget\tnative_chr\t25\t300\t1\n')
        row=score_unit(unit,out)
        self.assertTrue(row['retained']);self.assertEqual(row['target_hashes_found'],0)
        self.assertEqual(row['target_mapping_method'],'native_locus_same_frame_containment')
        self.assertEqual(row['diagnostic'],'native_gene_boundary_differs')
        for contig,start,end,strand in [('other',25,300,1),('native_chr',26,300,1),('native_chr',25,300,-1),('native_chr',301,330,1)]:
            genes.write_text(f'id\tcontig\tstart\tend\tstrand\ntarget\t{contig}\t{start}\t{end}\t{strand}\n')
            self.assertFalse(score_unit(unit,out)['retained'])

    def test_native_locus_mapping_uses_reverse_reading_frame(self):
        out=self._output(seq='M'+'C'*20)
        dna=self.root/'native.fna';dna.write_text('>native_chr\n'+'ATG'*110+'\n')
        truth=dict(self.unit['truth'],target_locus=dict(contig_sequence_sha256=sequence_hash('ATG'*110),start=1,end=300,strand=-1))
        unit=dict(self.unit,input=dict(type='genome',path=str(dna),mode='genome'),truth=truth)
        genes=out/'annotation/genes.tsv';genes.write_text('id\tcontig\tstart\tend\tstrand\ntarget\tnative_chr\t1\t276\t-1\n')
        self.assertTrue(score_unit(unit,out)['retained'])
        genes.write_text('id\tcontig\tstart\tend\tstrand\ntarget\tnative_chr\t1\t277\t-1\n')
        self.assertFalse(score_unit(unit,out)['retained'])

    def test_missing_workflow_artifact_fails_instead_of_becoming_fn(self):
        out = self.root / "missing"
        with self.assertRaises(FileNotFoundError):
            score_unit(self.unit, out)

    def test_duplicate_protein_hashes_keep_all_ids(self):
        out = self._output()
        with (out / "annotation/proteins.faa").open("a") as h:
            h.write(">duplicate\n" + self.seq + "\n")
        row = score_unit(self.unit, out)
        self.assertEqual(row["target_hashes_found"], 2)
        self.assertTrue(row["retained"])

    def test_duplicate_calls_do_not_hide_a_wrong_function_behind_a_correct_one(self):
        out=self._output()
        with (out/'annotation/proteins.faa').open('a') as h:h.write('>duplicate\n'+self.seq+'\n')
        with (out/'candidate_function_catalog.tsv').open('a') as h:h.write('duplicate\tcalibrated_high\tcandidate_H2_quinone_coupling\n')
        row=score_unit(self.unit,out)
        self.assertTrue(row['multiple_target_calls']);self.assertTrue(row['function_assigned'])
        self.assertFalse(row['function_correct'])
        self.assertEqual(score_metrics([row])['function']['FP'],1)
        self.assertEqual(score_metrics([row])['function']['FN'],1)

    def test_function_label_allowlist(self):
        self.assertFalse(is_specific_function("candidate_function_rule_missing"))
        self.assertTrue(is_specific_function("candidate_H2_ferredoxin_compatible_core"))
        self.assertFalse(is_specific_function("candidate_unknown_new_rule"))

    def test_manifest_binding_tamper_fails(self):
        admission = self.root / "admission.json"; admission.write_text("{}\n")
        inp = self.root / "x.faa"; inp.write_text(">x\nMA\n")
        gff = self.root / "x.gff"; gff.write_text("##gff-version 3\n")
        snapshot = ROOT / "revisions/research_round108_checked/snapshot.json"
        guard = {"metadata_path": str(PROTECTED), "metadata_sha256": hashlib.sha256(PROTECTED.read_bytes()).hexdigest(), "protected_hash_count": 14, "sequence_content_read": False}
        manifest = {"state": "frozen_before_prediction", "model_snapshot": {"path": str(snapshot), "sha256": hashlib.sha256(snapshot.read_bytes()).hexdigest()}, "protected_guard": guard,
                    "units": [{"unit_id": "x", "biological_source_id": "s", "admission": {"path": str(admission), "sha256": hashlib.sha256(admission.read_bytes()).hexdigest()}, "input": {"type": "proteins", "path": str(inp), "sha256": "bad", "mode": "transcript", "gff": {"path": str(gff), "sha256": hashlib.sha256(gff.read_bytes()).hexdigest()}}, "truth": {"hydrogenase": False}}]}
        path = self.root / "manifest.json"; path.write_text(json.dumps(manifest))
        with self.assertRaises(ValueError): validate_manifest(path)


if __name__ == "__main__":
    unittest.main()
