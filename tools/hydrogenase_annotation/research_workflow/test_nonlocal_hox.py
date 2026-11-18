import copy
import csv
import json
import unittest
from pathlib import Path

from pyridine_function import guard_cross_anchor_coupling, promote_nonlocal_sequence_potential
from pyridine_function import load_model, recover_module_members
from annotate import load_blast
from preferred_prediction import predict as preferred_predict


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "experiments/round155_Mcaps_Bath/global_role_diagnostic_v1.json"


class NonlocalHoxPredictionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = json.loads(FIXTURE.read_text())

    def valid(self):
        result = copy.deepcopy(self.payload["result"])
        result["evidence"]["roles"]["HoxH"][0]["domain_profiles"] = ["PF27537"]
        with (ROOT / "results/round155_Bath_A_v1/domains/query_domains.tsv").open() as handle:
            result["evidence"]["roles"]["HoxH"][0]["domain_profile_evidence"] = [
                r for r in csv.DictReader(handle, delimiter="\t") if r["protein_id"] == "q0001" and r["accession"].startswith("PF27537.")]
        option = {"provenance": {"input_scope": "nonlocal_sequence", "missing_structure": [],
                                  "ambiguous": False, "truncated": False}}
        return result, option

    def test_complete_distant_hox_roles_get_low_sequence_potential(self):
        result, option = self.valid()
        # Runtime genes.tsv coordinates are strings, unlike the diagnostic JSON.
        for item in result["evidence"]["coordinates"]["roles"].values():
            item["gene"].update({key: str(item["gene"][key]) for key in ("start", "end", "strand")})
        self.assertTrue(promote_nonlocal_sequence_potential(result, option, "genome", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))
        self.assertEqual(result["coupling_potential_status"], "unresolved")
        self.assertEqual(result["nonlocal_sequence_potential_status"], "sequence_supported_nonlocal_potential")
        preferred = preferred_predict({"candidate_id": "q0001", "family": "NiFe", "identity_status": "sequence_supported_uncalibrated_core"},
                                      evidence={"claims": {"pyridine": result}})
        self.assertEqual(preferred["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(preferred["prediction_confidence"], "low")
        self.assertIn("nonlocal_member_association", preferred["prediction_risks"])

    def test_missing_or_ambiguous_or_competing_roles_are_rejected(self):
        for mutation in ("missing", "ambiguous", "competition"):
            result, option = self.valid()
            if mutation == "missing":
                del result["evidence"]["roles"]["HoxU"]
            elif mutation == "ambiguous":
                result["role_audit"]["HoxF"]["candidate_count"] = 2
            else:
                result["module_competition_blockers"] = ["multiple_role_partner_candidates_HoxF"]
            self.assertFalse(promote_nonlocal_sequence_potential(result, option, "genome", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}), mutation)

    def test_scope_and_local_physical_status_remain_strict(self):
        result, option = self.valid()
        self.assertFalse(promote_nonlocal_sequence_potential(result, option, "assembly", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))
        result["evidence"]["physical_association_supported"] = True
        self.assertFalse(promote_nonlocal_sequence_potential(result, option, "genome", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))

    def test_nonlocal_scope_requires_intact_fuy_cassette_and_declared_bin(self):
        result, option = self.valid()
        result["evidence"]["coordinates"]["roles"]["HoxY"]["gene"]["contig"] = "other"
        self.assertFalse(promote_nonlocal_sequence_potential(result, option, "genome", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))
        result, option = self.valid()
        for role in ("HoxF", "HoxU", "HoxY"):
            result["evidence"]["coordinates"]["roles"][role]["gene"]["bin_id"] = "b1"
        result["evidence"]["coordinates"]["anchor"]["bin_id"] = "b1"
        option["provenance"]["selected_bin_ids"] = {"q0001": "b1", "q0002": "b1", "q0003": "b1", "q0004": "b1"}
        self.assertTrue(promote_nonlocal_sequence_potential(result, option, "bin", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))
        result, option = self.valid()
        for role in ("HoxF", "HoxU", "HoxY"):
            result["evidence"]["coordinates"]["roles"][role]["gene"]["bin_id"] = "b1"
        result["evidence"]["coordinates"]["anchor"]["bin_id"] = "b1"
        result["evidence"]["coordinates"]["roles"]["HoxY"]["gene"]["bin_id"] = "b2"
        option["provenance"]["selected_bin_ids"] = {"q0001": "b1", "q0002": "b1", "q0003": "b1", "q0004": "b2"}
        self.assertFalse(promote_nonlocal_sequence_potential(result, option, "bin", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))
        result, option = self.valid()
        result["evidence"]["roles"]["HoxH"][0]["partial"] = True
        result["evidence"]["coding_completeness_proof_validated"] = True
        result["evidence"]["coding_completeness_proof_used"] = ["q0001"]
        self.assertTrue(promote_nonlocal_sequence_potential(result, option, "genome", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))
        result, option = self.valid()
        result["evidence"]["core_partial_members"] = True
        self.assertFalse(promote_nonlocal_sequence_potential(result, option, "genome", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))

    def test_cross_anchor_guard_withdraws_nonlocal_hypothesis(self):
        result, _ = self.valid()
        result["evidence"]["roles"]["HoxF"][0]["query_domain_start"] = 10
        result["evidence"]["roles"]["HoxF"][0]["query_domain_end"] = 100
        result["nonlocal_sequence_potential"] = True
        result["coupling_blockers"] = []
        rival = {"coupling_blockers": [], "evidence": {"roles": {"other": [{"protein_id": "q0002", "query_domain_start": 20, "query_domain_end": 90}]}}}
        guard_cross_anchor_coupling({"a": result, "b": rival})
        self.assertFalse(result["nonlocal_sequence_potential"])

    def test_generic_nife_or_wrong_core_profile_is_not_promoted(self):
        result, option = self.valid()
        result["evidence"]["roles"]["HoxH"][0]["domain_profiles"] = ["PF02906"]
        result["evidence"]["roles"]["HoxH"][0]["domain_profile_evidence"] = []
        self.assertFalse(promote_nonlocal_sequence_potential(result, option, "genome", {"local_max_span_bp": 30000, "local_max_intergenic_gap_bp": 2500}))

    def test_actual_bath_raw_hits_resolve_all_four_roles(self):
        run = ROOT / "results/round155_Bath_A_v1/annotation"
        refs = json.loads((ROOT / "references/processed/references.json").read_text())
        thresholds = json.loads((ROOT / "config/thresholds.json").read_text())
        with (run / "genes.tsv").open() as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        genes = [{**g, "start": int(g["start"]), "end": int(g["end"]), "strand": int(g["strand"])} for g in genes]
        raw = load_blast(run / "blast.tsv", thresholds, refs)
        anchor = {**next(g for g in genes if g["id"] == "q0001"),
                  "protein_id": "q0001", "catalytic_family": "NiFe", "input_mode": "genome"}
        options = recover_module_members(anchor, raw, genes, load_model(ROOT), ROOT,
                                          structure_ids=["NiFe_HoxFUYH"], input_scope="nonlocal_sequence")
        self.assertEqual(len(options), 1)
        self.assertEqual({m["role"] for m in options[0]["members"]}, {"HoxH", "HoxY", "HoxF", "HoxU"})


if __name__ == "__main__":
    unittest.main()
