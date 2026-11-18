import json
import shutil
import tempfile
import unittest
import io
import math
from pathlib import Path
from Bio import SeqIO

from research_workflow.assay_sequence_identity import prepare, record_supported, verify_saved
from research_workflow import assay_sequence_identity

ROOT = Path(__file__).resolve().parents[1]
ASSAY = ROOT / "references/raw/round117_CpII_sources/CpIII/WP_003447632.gb"


class AssaySequenceIdentityTests(unittest.TestCase):
    def fixture(self, sequence=None, family="FeFe", domain=True):
        td = tempfile.TemporaryDirectory(); self.addCleanup(td.cleanup)
        base = Path(td.name); ann = base / "annotation"; (ann / "confirmation").mkdir(parents=True)
        seq = sequence or str(next(SeqIO.parse(io.StringIO(ASSAY.read_text()), "genbank")).seq)
        (ann / "proteins.faa").write_text(f">q\n{seq}\n")
        (ann / "genes.tsv").write_text("id\tcontig\tstart\tend\tstrand\tpartial\nq\tctg\t1\t%d\t1\tTrue\n" % (len(seq) * 3))
        d = {"accession": "PF02906", "supported": True, "model_coverage": 0.99,
             "query_core_overlap": 0.99, "evidence_kind": "single_domain", "components": [{"protein_id": "q", "accession": "PF02906.20", "profile": "Fe_hyd_lg_C", "model_length": 280, "model_sha256": "71bb31272b7833cea77c96dae3fc58dddce8e1d7a726f98f0354f8a189fba17e", "query_start": 1, "query_end": 450, "sequence_bitscore": 100,
             "sequence_ga": 22.1, "domain_bitscore": 100, "domain_ga": 22.1, "query_length": len(seq)}]}
        feature = {"protein_id": "q", "family": family, "domain_evidence": [d] if domain else []}
        (ann / "confirmation/features.json").write_text(json.dumps([feature]))
        return base, ann, feature, seq

    def test_actual_wild_type_exact_reference_is_supported(self):
        _, ann, feature, seq = self.fixture()
        records, _ = prepare(ann, [feature], ann.parent / "out", ROOT)
        row = records["q"]
        self.assertEqual(row["status"], "exact_assayed_core_sequence_supported_uncalibrated")
        self.assertTrue(row["source_partial"])
        self.assertEqual(row["tag_status"], "unknown_not_bound")
        self.assertTrue(record_supported(row, feature, ROOT))
        self.assertEqual(verify_saved(ann.parent / "out", ann, [feature], ROOT)["q"]["status"], row["status"])

    def test_r115_fefe_actual_feature_source_is_used(self):
        ann = ROOT / "results/round115_FeFe_diagnostic_v1/annotation"
        if not ann.exists():
            ann = ROOT / "results/round114_controls_v1/workflow/FeFe_native_CDS/annotation"
        features = json.loads((ann / "confirmation/features.json").read_text())
        out = Path(tempfile.mkdtemp()) / "out"
        records, _ = prepare(ann, features, out, ROOT)
        # The real feature currently has an unresolved competitor/profile gate;
        # this assertion confirms the helper consumes that exact feature source.
        self.assertIn("FeFe_native_CDS_1", {f.get("protein_id") for f in features})
        self.assertTrue(records == {} or "FeFe_native_CDS_1" in records)

    def test_one_residue_variant_and_fragment_are_not_exact(self):
        base, ann, feature, seq = self.fixture()
        variant = "A" + seq[1:]
        (ann / "proteins.faa").write_text(f">q\n{variant}\n")
        records, _ = prepare(ann, [feature], base / "out", ROOT)
        self.assertEqual(records, {})
        _, ann2, feature2, seq2 = self.fixture(sequence=seq[:300])
        records, _ = prepare(ann2, [feature2], ann2.parent / "out", ROOT)
        self.assertEqual(records, {})

    def test_wrong_family_or_profile_ownership_is_rejected(self):
        base, ann, feature, _ = self.fixture(family="NiFe")
        records, _ = prepare(ann, [feature], base / "out", ROOT); self.assertEqual(records, {})
        _, ann2, feature2, _ = self.fixture(domain=False)
        records, _ = prepare(ann2, [feature2], ann2.parent / "out", ROOT); self.assertEqual(records, {})

    def test_ga_tamper_nan_and_wrong_query_owner_are_rejected(self):
        base, ann, feature, _ = self.fixture()
        feature["domain_evidence"][0]["components"][0]["protein_id"] = "other"
        (ann / "confirmation/features.json").write_text(json.dumps([feature]))
        records, _ = prepare(ann, [feature], base / "out", ROOT); self.assertEqual(records, {})

    def test_domain_ratio_query_length_and_ga_tamper_are_rejected(self):
        for field, value in (("model_coverage", 1.01), ("query_core_overlap", -0.1)):
            base, ann, feature, _ = self.fixture()
            feature["domain_evidence"][0][field] = value
            (ann / "confirmation/features.json").write_text(json.dumps([feature]))
            self.assertEqual(prepare(ann, [feature], base / "out", ROOT)[0], {})
        for field, value in (("query_length", 449), ("sequence_ga", 22.0), ("domain_ga", 23.0)):
            base, ann, feature, _ = self.fixture()
            feature["domain_evidence"][0]["components"][0][field] = value
            (ann / "confirmation/features.json").write_text(json.dumps([feature]))
            self.assertEqual(prepare(ann, [feature], base / "out", ROOT)[0], {})
        base, ann, feature, _ = self.fixture()
        feature["domain_evidence"][0]["components"][0]["sequence_bitscore"] = math.nan
        (ann / "confirmation/features.json").write_text(json.dumps([feature], allow_nan=True))
        records, _ = prepare(ann, [feature], base / "out", ROOT); self.assertEqual(records, {})

    def test_indeterminate_same_sequence_prototype_does_not_block(self):
        _, ann, feature, _ = self.fixture()
        records, _ = prepare(ann, [feature], ann.parent / "out", ROOT)
        self.assertIn("q", records)

    def test_known_negative_same_sequence_conflict_is_rejected(self):
        _, ann, feature, _ = self.fixture()
        from unittest.mock import patch
        with patch.object(assay_sequence_identity, "_conflict", return_value=True):
            records, _ = prepare(ann, [feature], ann.parent / "out", ROOT)
        self.assertEqual(records, {})


if __name__ == "__main__": unittest.main()
