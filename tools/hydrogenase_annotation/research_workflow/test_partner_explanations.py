import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from pyridine_role_composition import compare_partner_explanations, _shareable_conflict


def hypothesis(anchor, role, protein, positions=(1, 2)):
    return {
        "protein_ids": [protein],
        "resources": {protein: list(positions)},
        "roles": [{"role": role, "protein_id": protein,
                   "resources": {protein: list(positions)}}],
    }


class PartnerExplanationTests(unittest.TestCase):
 def test_empty_partner_set_is_explicitly_unresolved(self):
    result = compare_partner_explanations("core", [])
    self.assertEqual(result["status"], "unresolved")
    self.assertIsNone(result["selected"])


 def test_two_complete_explanations_are_not_silently_ranked(self):
    result = compare_partner_explanations(
        "core", [hypothesis("core", "relay", "a"), hypothesis("core", "relay", "b")]
    )
    self.assertEqual(result["status"], "ambiguous")
    self.assertIsNone(result["selected"])
    self.assertEqual({x["protein_ids"][0] for x in result["explanations"]}, {"a", "b"})


 def test_remote_and_local_are_competing_whole_explanations(self):
    result = compare_partner_explanations(
        "core", [hypothesis("core", "relay", "local")],
        [hypothesis("core", "relay", "remote")],
    )
    self.assertEqual(result["status"], "ambiguous")
    self.assertEqual({x["kind"] for x in result["explanations"]}, {"local", "remote"})

 def test_compact_role_record_and_role_specific_sharing(self):
    option = hypothesis("core", "relay", "a")
    option["roles"] = [{"role": "relay", "protein_id": "a", "query_interval": [1, 2]}]
    result = compare_partner_explanations("core", [option])
    self.assertEqual(result["explanations"][0]["roles"], {"a": "relay"})
    self.assertTrue(_shareable_conflict("ferredoxin", "ferredoxin"))
    self.assertFalse(_shareable_conflict("ferredoxin", "structural_relay"))


if __name__ == "__main__":
    unittest.main()

