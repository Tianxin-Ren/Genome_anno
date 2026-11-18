import copy
import unittest

from pyridine_function import _strict_predict, load_model, predict


def evidence():
    roles = {r: [{"protein_id": r, "identity": 55, "query_coverage": .95,
                  "reference_coverage": .95, "partial": False,
                  "evalue_verified": True, "domain_organization_supported": True}]
             for r in ("SHI_alpha", "SHI_delta", "SHI_beta", "SHI_gamma")}
    return {"structure_id": "NiFe_HyhNADPH", "identity_status": "sequence_supported_uncalibrated_core",
            "roles": roles, "association_supported": True, "core_association_supported": True,
            "module_association_supported": True, "coupling_linkage_supported": False,
            "partial_members": False, "module_competition_audited": True,
            "module_competition_blockers": ["multiple_role_partner_candidates_SHI_beta"]}


class SelectedPyridineHypothesisTests(unittest.TestCase):
    def setUp(self):
        self.model = load_model()

    def test_strict_result_is_unchanged_and_selected_child_is_positive(self):
        ev = evidence()
        ev.update(physical_association_supported=True, core_physical_association_supported=True,
                  module_physical_association_supported=True, actual_resource_overlap=False,
                  selected_core_anchor_present=True, role_assignment_ambiguous=True)
        strict = _strict_predict(ev, self.model)
        out = predict(ev, self.model)
        self.assertEqual(strict["coupling_potential_status"], "unresolved")
        self.assertEqual({k:v for k,v in out.items() if k!="selected_hypothesis"},strict)
        self.assertEqual(out["selected_hypothesis"]["strict_child"]["coupling_potential_status"],
                         "local_coupling_potential_uncalibrated")

    def test_physical_false_or_real_overlap_blocks_selected_child(self):
        for key, value in (("physical_association_supported", False), ("actual_resource_overlap", True),
                           ("module_competition_blockers", ["incomplete_domain_organization_SHI_beta"])):
            ev = evidence()
            ev.update(physical_association_supported=True, core_physical_association_supported=True,
                      module_physical_association_supported=True, actual_resource_overlap=False,
                      selected_core_anchor_present=True, role_assignment_ambiguous=True)
            ev[key] = value
            self.assertNotIn("selected_hypothesis", predict(ev, self.model))

    def test_selected_child_still_requires_domain_coverage_partial_and_family_checks(self):
        for change in ("domain", "coverage", "partial", "family_competition", "unknown_overlap", "foreign_anchor"):
            ev=evidence()
            ev.update(physical_association_supported=True,core_physical_association_supported=True,
                      module_physical_association_supported=True,actual_resource_overlap=False,
                      selected_core_anchor_present=True,role_assignment_ambiguous=True)
            member=ev["roles"]["SHI_beta"][0]
            if change=="domain":member["domain_organization_supported"]=False
            if change=="coverage":member["reference_coverage"]=.1
            if change=="partial":member["partial"]=True
            if change=="family_competition":ev["module_competition_blockers"].append("adjacent_family_competition_unresolved_SHI_beta")
            if change=="unknown_overlap":ev.pop("actual_resource_overlap")
            if change=="foreign_anchor":ev["selected_core_anchor_present"]=False
            self.assertNotIn("selected_hypothesis",predict(ev,self.model),change)

    def test_no_ambiguity_has_no_selected_hypothesis(self):
        ev = evidence()
        ev.update(physical_association_supported=True, core_physical_association_supported=True,
                  module_physical_association_supported=True, actual_resource_overlap=False,
                  selected_core_anchor_present=True, role_assignment_ambiguous=False)
        self.assertNotIn("selected_hypothesis", predict(ev, self.model))


if __name__ == "__main__":
    unittest.main()
