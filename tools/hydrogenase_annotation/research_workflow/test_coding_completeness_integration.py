"""Small integration checks for coding-completeness as an identity gate."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from research_workflow.coding_completeness import prepare
from research_workflow.explain_candidates import sequence_status


ROOT = Path(__file__).resolve().parents[1]


class CodingCompletenessIntegrationTests(unittest.TestCase):
    def fixture(self, dna=True, partial="True", protein="MV", dna_override=None):
        td = tempfile.TemporaryDirectory(); self.addCleanup(td.cleanup)
        base = Path(td.name); ann = base / "annotation"; ann.mkdir()
        source = base / "source.fna"
        dna_sequence = dna_override or ("ATGGTTTAA" if protein == "MV" else "ATGCTGTAA")
        if dna:
            source.write_text(">ctg\n" + dna_sequence + "\n")
            inputs = [{"path": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}]
        else:
            inputs = []
        (ann / "proteins.faa").write_text(f">p1\n{protein}\n")
        (ann / "genes.tsv").write_text(
            "id\tcontig\tstart\tend\tstrand\tpartial\n"
            f"p1\tctg\t1\t9\t1\t{partial}\n")
        (ann / "prodigal.gff").write_text(
            "##gff-version 3\n# Sequence Data: ctg\n"
            "# Model Data: source;transl_table=11;\n"
            "ctg\tProdigal\tCDS\t1\t9\t.\t+\t0\tID=p1\n")
        (ann / "run_manifest.json").write_text(json.dumps({
            "status": "complete", "input_mode": "fna", "input_files": inputs}))
        rows, _ = prepare(ann, base / "coding", root=ROOT)
        return base, ann, rows["p1"]

    def feature(self, row, ligand="present"):
        return {
            "protein_id": "p1", "family": "NiFe", "protein_sha256": row["protein_sha256"],
            "eligibility_reasons": ["partial_CDS"], "positive_prototype": "P",
            "positive_score": 10.0, "core_domain_coverage": 1.0,
            "competitor_margin": 2.0, "competitor_prototype": "N",
            "source_active_center": {
                "center_organization": "single_polypeptide", "direct_ligand_status": ligand,
                "active_center_status": "template_ligand_set_supported"},
            "domain_evidence": [{
                "accession": "PF00374", "profile": "nife",
                "model_coverage": 1.0, "supported": True, "query_core_overlap": 1.0,
                "components": [{"sequence_bitscore": 10.0, "sequence_ga": 5.0,
                                "domain_bitscore": 10.0, "domain_ga": 5.0}]}]
        }

    def parameters(self):
        return {"min_domain_coverage": 0.7, "min_core_similarity": 5.0,
                "min_competitor_margin": 1.0}

    def test_bound_dna_proof_recovers_uncalibrated_status_only(self):
        _, ann, row = self.fixture(partial="True")
        feature = self.feature(row)
        status = sequence_status(
            {"catalytic_identity_confidence": "not_high_confidence"}, feature,
            {"query_partial": True}, self.parameters(), coding_completeness=row,
            gene={"contig": "ctg", "start": "1", "end": "9", "strand": "1"})
        self.assertEqual(status, "sequence_supported_uncalibrated_core")
        self.assertTrue(row["source_partial"])
        self.assertIn("partial_CDS", feature["eligibility_reasons"])

    def test_coding_proof_does_not_restore_missing_catalytic_cys(self):
        _, _, row = self.fixture(partial="True", protein="ML")
        feature = self.feature(row, ligand="substituted")
        status = sequence_status(
            {"catalytic_identity_confidence": "not_high_confidence"}, feature,
            {"query_partial": True}, self.parameters(), coding_completeness=row,
            gene={"contig": "ctg", "start": "1", "end": "9", "strand": "1"})
        self.assertNotEqual(status, "sequence_supported_uncalibrated_core")
        self.assertIn(status, {"evidence_insufficient", "sequence_support_center_model_unresolved"})
        self.assertEqual(row["status"], "complete_under_bound_translation_model")

    def test_no_dna_partial_stays_blocked(self):
        _, _, row = self.fixture(dna=False, partial="True")
        feature = self.feature(row)
        status = sequence_status(
            {"catalytic_identity_confidence": "not_high_confidence"}, feature,
            {"query_partial": True}, self.parameters(), coding_completeness=row,
            gene={"contig": "ctg", "start": "1", "end": "9", "strand": "1"})
        self.assertEqual(row["status"], "unknown")
        self.assertEqual(status, "evidence_insufficient")

    def test_missing_start_or_stop_cannot_recover_partial_identity(self):
        for broken in ("GCGGTTTAA", "ATGGTTAAA"):
            _, _, row = self.fixture(partial="True", dna_override=broken)
            self.assertNotEqual(row["status"], "complete_under_bound_translation_model")
            status = sequence_status(
                {"catalytic_identity_confidence": "not_high_confidence"}, self.feature(row),
                {"query_partial": True}, self.parameters(), coding_completeness=row,
                gene={"contig": "ctg", "start": "1", "end": "9", "strand": "1"})
            self.assertNotEqual(status, "sequence_supported_uncalibrated_core")

    def test_calibrated_high_is_never_added_by_coding_proof(self):
        _, _, row = self.fixture(partial="True")
        status = sequence_status(
            {"catalytic_identity_confidence": "calibrated_high"}, self.feature(row),
            {"query_partial": True}, self.parameters(), coding_completeness=row)
        self.assertEqual(status, "calibrated_high")
        self.assertNotEqual(status, "sequence_supported_uncalibrated_core")


if __name__ == "__main__":
    unittest.main()
