import csv
import json
import math
import unittest
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

from mbh_function import ROOT, _hyf_alias, load_model, predict, prepare


class HyfCoreAliasTests(unittest.TestCase):
    """Trust-boundary tests for the reviewed, unresolved EcHyfG alias."""

    @classmethod
    def setUpClass(cls):
        cls.root = Path(ROOT)
        cls.model = load_model(cls.root)
        cls.annotation = cls.root / "results/round83_Hadar_FNA_v1/u1859624/annotation"
        cls.domains_path = cls.root / "results/round90_FeFe_backbone_release_v1/u1859624/domains/query_domains.tsv"
        with (cls.annotation / "catalytic_anchors.tsv").open() as h:
            cls.anchors = list(csv.DictReader(h, delimiter="\t"))
        with cls.domains_path.open() as h:
            cls.domains = list(csv.DictReader(h, delimiter="\t"))
        with (cls.root / "results/round90_FeFe_backbone_release_v1/u1859624/functions/candidate_progression.tsv").open() as h:
            progression = {r["candidate_id"]: r for r in csv.DictReader(h, delimiter="\t")}
        cls.identities = {a["protein_id"]: progression.get(a["protein_id"], {}).get("identity_status") for a in cls.anchors}

    def _cert_context(self):
        review = json.loads((self.root / self.model["functional_projection_review"]).read_text())
        refs = {r["reference_id"]: r for r in json.loads((self.root / "references/processed/references.json").read_text())}
        cert = _hyf_alias(review, refs, self.model, self.root)
        return review, refs, cert

    def _real_results(self):
        return prepare(self.annotation, self.anchors, self.domains, self.identities, self.model, root=self.root)[0]

    def test_certificate_is_real_unresolved_mbh_l_role_with_four_cys_sites(self):
        review, refs, cert = self._cert_context()
        self.assertEqual(cert["reference_id"], "unresolved__ecoli_MG1655_hyfG")
        self.assertEqual(cert["role"], "MbhL")
        self.assertEqual(cert["reference_kind_retained"], "unresolved")
        self.assertEqual([s["source_position"] for s in cert["homologous_Cys_sites"]], [241, 244, 531, 534])
        self.assertGreaterEqual(cert["paired_identity"], self.model["min_identity"])
        self.assertEqual(refs[cert["reference_id"]]["kind"], "unresolved")

    def test_real_u1859624_722_gets_coarse_fd_potential_but_keeps_exact_rival(self):
        results = self._real_results()
        qid = next(pid for pid in results if pid.endswith("_722"))
        out = results[qid]
        self.assertEqual(out["core_coupling_status"], "reference_like_H2_Fd_coupling_potential")
        self.assertEqual(out["evidence"]["functional_projection"]["reference_groups"]["MbhL"].count("unresolved__ecoli_MG1655_hyfG"), 1)
        self.assertEqual(out["evidence"]["functional_projection"]["homology_alias_reviews"]["unresolved__ecoli_MG1655_hyfG"]["reference_kind_retained"], "unresolved")
        # The reviewed alias is only in the coarse projection; the exact Mbh
        # lineage evidence remains separately unresolved and retains blockers.
        self.assertTrue(out["evidence"]["exact_reference_evidence"]["core_blockers"])
        self.assertIn("MbhL", out["evidence"]["exact_reference_evidence"]["roles"])

    def test_missing_n_non_nife_rival_or_weak_identity_stays_unresolved(self):
        out = next(v for k, v in self._real_results().items() if k.endswith("_722"))
        base = deepcopy(out["evidence"])
        missing_n = deepcopy(base)
        missing_n["roles"]["MbhN"] = []
        self.assertEqual(predict(missing_n, self.model)["core_coupling_status"], "unresolved")
        rival = deepcopy(base)
        rival["audits"]["MbhL"]["functional_status"] = "qualified_non_NiFe_rival"
        self.assertEqual(predict(rival, self.model)["core_coupling_status"], "unresolved")
        weak = deepcopy(base)
        weak["identity_status"] = "evidence_insufficient"
        self.assertEqual(predict(weak, self.model)["core_coupling_status"], "unresolved")

    def test_alignment_metric_tamper_is_rejected_after_outer_source_binding(self):
        review, refs, cert = self._cert_context()
        altered = deepcopy(cert)
        altered["paired_identity"] = cert["paired_identity"] - 1.0
        cert_path = self.root / review["homology_aliases"][cert["reference_id"]]
        original_read_text = Path.read_text

        def read_tampered(path, *args, **kwargs):
            if path == cert_path:
                return json.dumps(altered)
            return original_read_text(path, *args, **kwargs)

        with patch.object(Path, "read_text", new=read_tampered):
            with self.assertRaisesRegex(ValueError, "homology changed"):
                _hyf_alias(review, refs, self.model, self.root)

    def test_cys_certificate_and_native_genbank_binding_are_independently_checked(self):
        review, refs, cert = self._cert_context()
        altered = deepcopy(cert)
        altered["homologous_Cys_sites"][0]["alias_residue"] = "A"
        cert_path = self.root / review["homology_aliases"][cert["reference_id"]]
        original_read_text = Path.read_text

        def read_tampered(path, *args, **kwargs):
            if path == cert_path:
                return json.dumps(altered)
            return original_read_text(path, *args, **kwargs)

        with patch.object(Path, "read_text", new=read_tampered):
            with self.assertRaisesRegex(ValueError, "Cys correspondence changed"):
                _hyf_alias(review, refs, self.model, self.root)

        from Bio import SeqIO
        original_parse = SeqIO.parse

        def broken_genbank(path, fmt):
            records = list(original_parse(path, fmt))
            if fmt == "genbank":
                for record in records:
                    for feature in record.features:
                        if feature.type == "CDS" and feature.qualifiers.get("locus_tag") == [cert["source_locus"]]:
                            feature.qualifiers["translation"] = ["X"]
            return iter(records)

        with patch("mbh_function.SeqIO.parse", side_effect=broken_genbank):
            with self.assertRaisesRegex(ValueError, "native CDS binding changed"):
                _hyf_alias(review, refs, self.model, self.root)

    def test_model_source_hash_binding_rejects_even_valid_certificate(self):
        review, refs, cert = self._cert_context()
        altered_model = deepcopy(self.model)
        source_name = cert["source_fasta"]
        altered_model["sources"][source_name] = "0" * 64
        with self.assertRaisesRegex(ValueError, "Unbound HyfG alias source"):
            _hyf_alias(review, refs, altered_model, self.root)

    def test_alias_model_source_hash_guard_cannot_be_bypassed_by_nonfinite_metrics(self):
        review, refs, cert = self._cert_context()
        from annotate import load_blast
        raw = load_blast(self.root / cert["alignment_path"], dict(blast_evalue=self.model["max_evalue"], min_identity=0, min_query_coverage=0, min_reference_coverage=0))
        altered = deepcopy(raw)
        altered[0]["bitscore"] = math.nan
        with patch("mbh_function.load_blast", return_value=altered):
            with self.assertRaisesRegex(ValueError, "alignment score invalid"):
                _hyf_alias(review, refs, self.model, self.root)


if __name__ == "__main__":
    unittest.main()
