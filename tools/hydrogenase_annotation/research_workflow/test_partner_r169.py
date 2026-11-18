import unittest

from preferred_prediction import _component_prediction


def _base():
    return {
        "preferred_function": "candidate_H2_pyridine_nucleotide_coupling",
        "preferred_family": "NiFe",
        "prediction_kind": "coarse_function",
        "prediction_confidence": "medium",
        "prediction_summary": "existing core",
        "prediction_support": [],
        "preferred_members": [{"role": "HydA", "protein_id": "anchor"}],
        "prediction_alternatives": [],
        "prediction_risks": [],
    }


def _hypothesis(partner, *, supported=True):
    return {
        "status": "source_equivalent_components_potential",
        "structure_id": "HydAB_source_equivalent",
        "matched_members": [
            {"role": "HydA", "protein_id": "anchor"},
            {"role": "HydB", "protein_id": partner},
        ],
        "linkage_blockers": ["whole_reference_tuple_local_linkage_missing"],
        "strict_tuple_blockers": ["strict_tuple_onlylinkage"],
        "capabilities": ([{
            "scope": "sequence_component_function_hypothesis",
            "reaction_key": "nucleotide_coupling",
            "electron_acceptors": ["pyridine_nucleotide"],
            "electron_donors": ["protons"],
        }] if supported else [{
            "scope": "sequence_component_function_hypothesis",
            "molecular_reaction": "redox activity not typed",
        }]),
    }


def _run(hypotheses):
    row = {"candidate_id": "anchor", "identity_status": "calibrated_high", "preferred_family": "NiFe"}
    return _component_prediction(row, [{"claim_name": "reference", "sequence_component_hypotheses": hypotheses}], _base())


class R169PartnerRankingTests(unittest.TestCase):
    def test_single_strong_local_and_remote_are_conditional_not_rejected(self):
        local = _run([_hypothesis("local")])
        remote = _run([_hypothesis("remote")])
        self.assertEqual(local["preferred_function"], remote["preferred_function"])
        self.assertEqual(local["preferred_members"][-1]["protein_id"], "local")
        self.assertEqual(remote["preferred_members"][-1]["protein_id"], "remote")

    def test_equal_supported_same_function_is_member_ambiguous(self):
        out = _run([_hypothesis("local"), _hypothesis("remote")])
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["preferred_members"], [{"role": "HydA", "protein_id": "anchor"}])
        self.assertIn("partner_assignment_ambiguous", out["prediction_risks"])
        alternatives = [{m["protein_id"] for m in a.get("preferred_members", [])}
                        for a in out["prediction_alternatives"]]
        self.assertTrue(any(ids & {"local", "remote"} for ids in alternatives))

    def test_weak_only_does_not_promote_and_unique_support_withdraws(self):
        weak = _run([_hypothesis("weak", supported=False)])
        self.assertIsNone(weak)
        mixed = _run([_hypothesis("strong"), _hypothesis("weak", supported=False)])
        self.assertIsNotNone(mixed)
        withdrawn = _run([_hypothesis("weak", supported=False)])
        self.assertIsNone(withdrawn)

    def test_rename_and_input_order_have_same_normalized_outcome(self):
        forward = _run([_hypothesis("local"), _hypothesis("remote")])
        reverse = _run([_hypothesis("remote"), _hypothesis("local")])
        self.assertEqual((forward["preferred_function"], forward["prediction_kind"],
                          sorted(forward["prediction_risks"])),
                         (reverse["preferred_function"], reverse["prediction_kind"],
                          sorted(reverse["prediction_risks"])))
        renamed = _run([_hypothesis("z_copy"), _hypothesis("a_copy")])
        self.assertEqual(renamed["preferred_function"], forward["preferred_function"])
        self.assertEqual(renamed["prediction_kind"], forward["prediction_kind"])
        self.assertIn("partner_assignment_ambiguous", renamed["prediction_risks"])
        # Alternative order is deterministic after normalizing the member
        # identifiers (the selected candidate itself is not biologically
        # meaningful in this tie).
        def normalized(out, rename):
            def members(item):
                return tuple(sorted((m.get("role"), rename.get(m.get("protein_id"), m.get("protein_id")))
                                   for m in item.get("preferred_members", [])))
            return sorted(members(a) for a in out["prediction_alternatives"] if a.get("preferred_members"))
        self.assertEqual(normalized(forward, {"local": "a", "remote": "b"}),
                         normalized(reverse, {"local": "a", "remote": "b"}))


if __name__ == "__main__":
    unittest.main()
