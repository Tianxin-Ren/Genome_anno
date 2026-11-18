import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from claim_assertions import build_assertion, from_results


class ClaimAssertionTests(unittest.TestCase):
    def test_split_center_is_scoped_and_does_not_claim_partner_or_mechanism(self):
        attrs = [{"candidate_id": "split", "attribute": "NiFe_split_catalytic_core_potential",
                  "status": "sequence_supported_candidate",
                  "value": {"members": [{"role": "VhuA", "protein_id": "a"},
                                           {"role": "VhuU", "protein_id": "u"}],
                            "limitations": ["full_complex_unresolved"]}}]
        out = build_assertion({"candidate_id": "split",
                               "identity_status": "sequence_supported_uncalibrated_split_center",
                               "coarse_function_prediction": "NiFe_split_catalytic_center_supported_function_unresolved",
                               "preferred_function": "NiFe_split_catalytic_center_supported_function_unresolved"}, attrs)
        self.assertEqual(out["current_supported_annotation"],
                         "NiFe_split_catalytic_center_supported_function_unresolved")
        self.assertEqual(out["dimensions"]["catalytic_identity"]["status"], "supported")
        self.assertEqual(out["dimensions"]["catalytic_identity"]["scope"], "split_catalytic_center_only")
        self.assertEqual(out["member_evidence"], attrs[0]["value"]["members"])
        self.assertEqual(out["dimensions"]["electron_partner"]["status"], "not_assessed")
        self.assertEqual(out["dimensions"]["coupling_mechanism"]["status"], "not_assessed")

    def test_split_identity_without_member_proof_does_not_invent_members(self):
        out = build_assertion({"candidate_id": "split",
                               "identity_status": "sequence_supported_uncalibrated_split_center",
                               "coarse_function_prediction": "NiFe_split_catalytic_center_supported_function_unresolved",
                               "preferred_function": "NiFe_split_catalytic_center_supported_function_unresolved",
                               "preferred_members": [{"protein_id": "synthetic"}]})
        self.assertEqual(out["member_evidence"], [])
        self.assertEqual(out["preferred_function_hypothesis"], "possible_hydrogenase")
        self.assertEqual(out["current_supported_annotation"], "")
        self.assertEqual(out["dimensions"]["catalytic_identity"]["status"], "unresolved")

    def test_co_context_is_conditional_and_withdraws_on_route_change(self):
        row = {"candidate_id": "co1", "identity_status": "calibrated_high",
               "coarse_function_prediction": "H2_core_supported_partner_unresolved",
               "preferred_function": "candidate_CO_linked_H2_pathway"}
        attr = {"candidate_id": "co1", "attribute": "candidate_CO_linked_H2_pathway",
                "status": "candidate_CO_linked_H2_pathway",
                "value": {"entry_route": "hyc_core_local_CODH_context"},
                "evidence_source": "coo_function_evidence.json"}
        out = build_assertion(row, [attr])
        self.assertEqual(out["current_supported_annotation"], "H2_core_supported_partner_unresolved")
        self.assertEqual(out["association_status"], "conditional")
        self.assertEqual(out["dimensions"]["electron_partner"]["hypothesis"], "carbon_monoxide_associated")
        withdrawn = build_assertion(row, [{**attr, "value": {"entry_route": "legacy_CooH_L_U_X_K"}}])
        self.assertNotIn("hypothesis", withdrawn["dimensions"]["electron_partner"])
        self.assertNotEqual(withdrawn["association_status"], "conditional")
        alternate = build_assertion({**row, "preferred_function": "candidate_H2_formate_coupling"}, [attr])
        self.assertNotEqual(alternate["association_status"], "conditional")

    def test_identity_survives_unresolved_partner(self):
        row = {"candidate_id": "u1", "identity_status": "sequence_supported_uncalibrated_core",
               "function_prediction": "candidate_H2_pyridine_nucleotide_coupling",
               "preferred_function": "candidate_H2_pyridine_nucleotide_coupling",
               "query_activity": "not_tested"}
        out = build_assertion(row, [{"candidate_id": "u1", "attribute": "candidate_complex_architecture",
                                    "status": "unresolved", "value": "partner_unknown"}])
        self.assertEqual(out["dimensions"]["catalytic_identity"]["status"], "supported")
        self.assertEqual(out["association_status"], "unresolved")
        self.assertEqual(out["dimensions"]["reaction_direction_conditions"]["status"], "unresolved")

    def test_no_evidence_does_not_create_specific_mechanism(self):
        out = build_assertion({"candidate_id": "u2", "identity_status": "unconfirmed"})
        self.assertEqual(out["dimensions"]["coupling_mechanism"]["status"], "not_assessed")
        self.assertEqual(out["novelty_status"], "not_assessed")

    def test_evaluated_lineage_and_center_uncertainty_are_unresolved(self):
        for value in ("hydrogenase_like_lineage_unconfirmed", "sequence_support_center_model_unresolved"):
            with self.subTest(value=value):
                identity = build_assertion({"identity_status": value})["dimensions"]["catalytic_identity"]
                self.assertEqual(identity["status"], "unresolved")
                self.assertEqual(identity["value"], value)
        for row in ({}, {"identity_status": ""}, {"identity_status": "not_assessed"}):
            identity = build_assertion(row)["dimensions"]["catalytic_identity"]
            self.assertEqual(identity["status"], "not_assessed")

    def test_catalytic_region_support_is_local_to_the_region(self):
        for value in ("sequence_supported_uncalibrated_core_region", "sequence_supported_core_region"):
            with self.subTest(value=value):
                out = build_assertion({"identity_status": value})
                identity = out["dimensions"]["catalytic_identity"]
                self.assertEqual(identity["status"], "supported")
                self.assertEqual(identity["value"], value)
                self.assertEqual(identity["scope"], "catalytic_core_region_only")
                self.assertEqual(out["dimensions"]["electron_partner"]["status"], "not_assessed")
                self.assertEqual(out["association_status"], "not_assessed")

    def test_blocked_and_ambiguous_are_not_support(self):
        attrs = [
            {"candidate_id": "u3", "attribute": "electron_coupling", "status": "blocked", "value": "x"},
            {"candidate_id": "u3", "attribute": "candidate_complex_architecture", "status": "partner_ambiguous", "value": "x"},
        ]
        out = build_assertion({"candidate_id": "u3", "identity_status": "unconfirmed"}, attrs)
        self.assertEqual(out["dimensions"]["coupling_mechanism"]["status"], "unresolved")
        self.assertEqual(out["association_status"], "unresolved")

    def test_association_name_does_not_create_ion_or_coupling_claim(self):
        attrs = [{"candidate_id": "u4", "attribute": "candidate_complex_architecture",
                  "status": "unresolved", "value": "locality_unknown"}]
        out = build_assertion({"candidate_id": "u4", "identity_status": "unconfirmed"}, attrs)
        self.assertEqual(out["dimensions"]["coupling_mechanism"]["status"], "not_assessed")

    def test_membrane_membership_is_separate_from_ion_pumping_and_competition(self):
        row = {"preferred_function": "candidate_membrane_associated_H2_redox_system"}
        attrs = [{"attribute": "Group4_local_membrane_H2_redox_potential",
                  "status": "local_membrane_H2_redox_potential_uncalibrated"}]
        out = build_assertion(row, attrs)
        self.assertEqual(out["association_status"], "supported")
        self.assertEqual(out["dimensions"]["localization_topology"]["status"], "supported")
        self.assertIn("ion_translocation", out["dimensions"]["coupling_mechanism"]["unresolved_attributes"])
        row["prediction_risks"] = ["cross_anchor_partner_ownership_ambiguous"]
        self.assertEqual(build_assertion(row, attrs)["association_status"], "ambiguous")

    def test_partner_hypothesis_requires_matching_evidence(self):
        out = build_assertion({"preferred_function": "candidate_H2_quinone_coupling"}, [])
        self.assertEqual(out["dimensions"]["electron_partner"]["status"], "not_assessed")
        self.assertNotIn("hypothesis", out["dimensions"]["electron_partner"])

    def test_unresolved_matching_partner_is_not_promoted(self):
        attrs = [{"attribute": "huc_local_h2_quinone_coupling_potential",
                  "status": "unresolved"},
                 {"attribute": "Pyridine_nucleotide_electron_module", "status": "supported"}]
        out = build_assertion({"preferred_function": "candidate_H2_quinone_coupling"}, attrs)
        partner = out["dimensions"]["electron_partner"]
        self.assertEqual(partner["status"], "unresolved")
        self.assertEqual(partner["hypothesis"], "quinone")
        self.assertEqual(partner["hypothesis_scope"], "conditional_unresolved")
        self.assertEqual(partner["hypothesis_evidence"][0]["status"], "unresolved")

    def test_fusion_partner_scope_requires_exact_reference_assay(self):
        row = {"preferred_function": "candidate_FeFe_rubrerythrin_H2_H2O2_coupling"}
        attr = {"attribute": "FeFe_rubrerythrin_function", "status": "assay_supported_candidate",
                "value": {"reference_assay": {"electron_acceptor": "H2O2"}}}
        dimensions = build_assertion(row, [attr])["dimensions"]
        partner = dimensions["electron_partner"]
        self.assertEqual(partner["status"], "supported")
        self.assertEqual(partner["hypothesis"], "hydrogen_peroxide")
        self.assertEqual(partner["hypothesis_scope"], "exact_reference_assay_potential")
        self.assertEqual(dimensions["reaction_direction_conditions"]["status"], "unresolved")
        attr["value"] = {}
        self.assertNotIn("hypothesis", build_assertion(row, [attr])["dimensions"]["electron_partner"])

    def test_explicit_hya_hyb_and_hdcr_attributes_support_partner_hypotheses(self):
        for label, attribute, status, hypothesis in (
            ("quinone", "Hya_like_Hyd_local_H2_quinone_coupling_potential",
             "local_quinone_coupling_potential_uncalibrated", "quinone"),
            ("quinone", "Hyb_like_Hyd_local_H2_quinone_coupling_potential",
             "local_quinone_coupling_potential_uncalibrated", "quinone"),
            ("formate", "H2_formate_coupling_potential", "candidate_H2_formate_coupling", "formate_associated"),
        ):
            with self.subTest(attribute=attribute):
                row = {"preferred_function": f"candidate_H2_{label}_coupling"}
                attrs = [{"attribute": attribute, "status": status}]
                partner = build_assertion(row, attrs)["dimensions"]["electron_partner"]
                self.assertEqual(partner["hypothesis"], hypothesis)
                self.assertEqual(partner["hypothesis_scope"], "sequence_prediction")

    def test_hya_hyb_full_coupling_supports_association_but_partial_paths_do_not(self):
        for attribute in (
            "hya_like_hyd_local_h2_quinone_coupling_potential",
            "hyb_like_hyd_local_h2_quinone_coupling_potential",
        ):
            with self.subTest(attribute=attribute):
                attrs = [
                    {"candidate_id": "x", "attribute": attribute,
                     "status": "local_quinone_coupling_potential_uncalibrated"},
                    {"candidate_id": "x", "attribute": "candidate_complex_architecture",
                     "status": "unresolved"},
                ]
                out = build_assertion({"candidate_id": "x",
                                       "preferred_function": "candidate_H2_quinone_coupling"}, attrs)
                self.assertEqual(out["dimensions"]["component_association"]["status"], "supported")
                self.assertEqual(out["dimensions"]["coupling_mechanism"]["status"], "not_assessed")
                self.assertEqual(out["dimensions"]["reaction_direction_conditions"]["status"], "unresolved")
                self.assertIn("ion_translocation", out["dimensions"]["coupling_mechanism"]["unresolved_attributes"])

                shared = build_assertion(
                    {"candidate_id": "x", "preferred_function": "candidate_H2_quinone_coupling"}, attrs,
                    {"status": "shared_role_conditionally_allowed"},
                )
                self.assertEqual(shared["association_status"], "conditional")
                ambiguous = build_assertion(
                    {"candidate_id": "x", "preferred_function": "candidate_H2_quinone_coupling"}, attrs,
                    {"status": "ambiguous"},
                )
                self.assertEqual(ambiguous["association_status"], "ambiguous")

        for attribute in (
            "hya_like_hyd_quinone_h2_reaction_core",
            "hya_like_hyd_quinone_fes_relay",
            "hya_like_hyd_quinone_membrane_module",
        ):
            with self.subTest(attribute=attribute):
                out = build_assertion(
                    {"candidate_id": "x", "preferred_function": "candidate_H2_quinone_coupling"},
                    [{"candidate_id": "x", "attribute": attribute,
                      "status": "sequence_supported_candidate"},
                     {"candidate_id": "x", "attribute": "candidate_complex_architecture",
                      "status": "unresolved"}],
                )
                self.assertEqual(out["dimensions"]["component_association"]["status"], "unresolved")

    def test_predicted_partner_scope_tracks_conditional_component_association(self):
        row = {"candidate_id": "split", "identity_status": "sequence_supported_uncalibrated_core",
               "preferred_function": "candidate_H2_pyridine_nucleotide_coupling"}
        attrs = [{"candidate_id": "split", "attribute": "Role_composed_H2_pyridine_coupling_potential",
                  "status": "local_coupling_potential_uncalibrated"}]
        supported = build_assertion(row, attrs)
        self.assertEqual(supported["association_status"], "supported")
        self.assertEqual(supported["dimensions"]["electron_partner"]["status"], "supported")
        self.assertEqual(supported["dimensions"]["electron_partner"]["hypothesis_scope"], "sequence_prediction")

        conditional = build_assertion(row, attrs,
                                      {"status": "shared_role_conditionally_allowed"})
        partner = conditional["dimensions"]["electron_partner"]
        self.assertEqual(conditional["association_status"], "conditional")
        self.assertEqual(partner["status"], "supported")
        self.assertEqual(partner["hypothesis"], "pyridine_nucleotide_unspecified")
        self.assertEqual(partner["hypothesis_scope"], "conditional_unresolved")
        self.assertEqual(partner["hypothesis_evidence"], supported["dimensions"]["electron_partner"]["hypothesis_evidence"])
        self.assertEqual(conditional["dimensions"]["catalytic_identity"]["status"], "supported")

        ambiguous = build_assertion(row, attrs, {"status": "ambiguous"})
        self.assertEqual(ambiguous["association_status"], "ambiguous")
        self.assertEqual(ambiguous["dimensions"]["electron_partner"]["hypothesis_scope"], "conditional_unresolved")

    def test_exact_reference_partner_scope_is_not_downgraded_by_association(self):
        row = {"candidate_id": "ref", "preferred_function": "candidate_FeFe_rubrerythrin_H2_H2O2_coupling"}
        attrs = [{"candidate_id": "ref", "attribute": "FeFe_rubrerythrin_function",
                  "status": "assay_supported_candidate",
                  "value": {"reference_assay": {"electron_acceptor": "H2O2"}}}]
        out = build_assertion(row, attrs, {"status": "ambiguous"})
        partner = out["dimensions"]["electron_partner"]
        self.assertEqual(out["association_status"], "ambiguous")
        self.assertEqual(partner["hypothesis_scope"], "exact_reference_assay_potential")

    def test_pyridine_module_alone_does_not_link_an_electron_partner(self):
        module = {"attribute": "Pyridine_nucleotide_electron_module", "status": "sequence_supported_candidate"}
        link = {"attribute": "Role_composed_H2_pyridine_coupling_potential",
                "status": "local_coupling_potential_uncalibrated"}
        for preferred in ("candidate_H2_pyridine_nucleotide_coupling", "hydrogenase_core_candidate"):
            with self.subTest(preferred=preferred):
                row = {"preferred_function": preferred, "identity_status": "sequence_supported_uncalibrated_core"}
                linked = build_assertion(row, [module, link])["dimensions"]
                out = build_assertion(row, [module])
                partner = out["dimensions"]["electron_partner"]
                self.assertEqual(linked["electron_partner"]["status"], "supported")
                self.assertEqual(partner["status"], "unresolved")
                self.assertEqual(partner["evidence"], [])
                self.assertNotIn("hypothesis", partner)
                self.assertEqual(partner["module_evidence"], linked["electron_partner"]["module_evidence"])
                self.assertEqual(partner["evidence_ids"], [partner["module_evidence"][0]["evidence_id"]])
                self.assertEqual(out["dimensions"]["catalytic_identity"], linked["catalytic_identity"])
                self.assertEqual(out["preferred_function_hypothesis"], preferred)

    def test_membrane_label_requires_topology_evidence(self):
        row = {"preferred_function": "candidate_membrane_associated_H2_redox_system"}
        self.assertNotIn("hypothesis", build_assertion(row)["dimensions"]["localization_topology"])

    def test_top_level_preferred_function_remains_compatibility_field(self):
        out = build_assertion({"preferred_function": "candidate_H2_quinone_coupling"}, [])
        self.assertEqual(out["preferred_function_hypothesis"], "candidate_H2_quinone_coupling")

    def test_evidence_ids_bind_rows_and_withdraw_only_the_dependent_hypothesis(self):
        row = {"candidate_id": "core", "identity_status": "sequence_supported_uncalibrated_core",
               "preferred_function": "candidate_H2_quinone_coupling"}
        quinone = {"candidate_id": "core", "attribute": "huc_local_h2_quinone_coupling_potential",
                   "status": "supported", "value": '{"members":["A","B"]}', "evidence_source": "huc.json"}
        other = {"attribute": "Pyridine_nucleotide_electron_module", "status": "supported"}
        reference = {"attribute": "Reference_assay_reaction_capabilities", "status": "inherited_reference_potential"}
        constraint = {"attribute": "huc_local_h2_quinone_coupling_potential", "status": "unresolved",
                      "evidence_source": "alternative_architecture.json"}
        attrs = [quinone, other, reference, constraint, dict(quinone, value={"members": ["A", "B"]})]
        dimensions = build_assertion(row, attrs)["dimensions"]
        partner = dimensions["electron_partner"]
        self.assertEqual(partner["evidence_ids"], build_assertion(row, attrs[::-1])["dimensions"]["electron_partner"]["evidence_ids"])
        self.assertEqual(len(partner["hypothesis_evidence_ids"]), 1)
        self.assertIn(partner["hypothesis_evidence_ids"][0], partner["evidence_ids"])
        self.assertEqual(partner["constraint_evidence"][0]["evidence_role"], "constraint")
        self.assertEqual(dimensions["reaction_direction_conditions"]["evidence"][0]["evidence_role"], "reference_only")
        self.assertEqual(dimensions["reaction_direction_conditions"]["status"], "unresolved")
        revised = build_assertion(row, [other, reference, constraint])["dimensions"]
        self.assertEqual(revised["electron_partner"]["hypothesis_scope"], "conditional_unresolved")
        self.assertNotIn(partner["hypothesis_evidence_ids"][0], revised["electron_partner"]["evidence_ids"])
        self.assertEqual(revised["catalytic_identity"], dimensions["catalytic_identity"])

    def test_unresolved_partner_competition_is_exported_without_erasing_architecture(self):
        row = {"candidate_id": "core", "identity_status": "sequence_supported_uncalibrated_core",
               "preferred_function": "candidate_H2_pyridine_nucleotide_coupling",
               "coarse_function_prediction": "candidate_H2_pyridine_nucleotide_coupling"}
        attrs = [{"candidate_id": "core", "attribute": "candidate_complex_architecture",
                  "status": "supported", "value": "{\"members\":[\"core\"]}",
                  "evidence_source": "architecture"}]
        assignment = {"status": "unresolved", "selected": None,
                      "reason": "exclusive_resource_overlap_with_other_anchor",
                      "explanations": [{"kind": "local", "protein_ids": ["core"]}]}
        out = build_assertion(row, attrs, assignment)
        self.assertEqual(out["association_status"], "supported")
        self.assertEqual(out["dimensions"]["component_association"]["status"], "supported")
        self.assertEqual(out["dimensions"]["component_association"]["competition"]["status"], "unresolved")
        self.assertEqual(out["dimensions"]["component_association"]["competition"]["reason"], assignment["reason"])

    def test_claim_producer_manifest_binds_the_files_parsed(self):
        with tempfile.TemporaryDirectory() as tmp:
            functions = Path(tmp)
            tables = {
                "functional_annotations.tsv": "candidate_id\tidentity_status\na\tunresolved\n",
                "function_attributes.tsv": "candidate_id\tattribute\tstatus\na\tquery_activity\tnot_tested\n",
                "preferred_predictions.json": '{"a":{"preferred_function":"unresolved"}}',
            }
            for name, value in tables.items():
                (functions / name).write_text(value)
            manifest = {"status": "claim_specific_explanation_complete",
                        "functional_annotations_sha256": hashlib.sha256(tables["functional_annotations.tsv"].encode()).hexdigest(),
                        "preferred_prediction_output_sha256": hashlib.sha256(tables["preferred_predictions.json"].encode()).hexdigest(),
                        "pyridine_role_composition_output_sha256": None}
            path = functions / "manifest.json"
            path.write_text(json.dumps(manifest))
            self.assertEqual(json.loads(from_results(functions).read_text())["source_binding"]["status"],
                             "legacy_inputs_unpinned")
            manifest["function_attributes_sha256"] = hashlib.sha256(tables["function_attributes.tsv"].encode()).hexdigest()
            path.write_text(json.dumps(manifest))
            out = from_results(functions)
            self.assertEqual(json.loads(out.read_text())["source_binding"]["status"], "producer_manifest_verified")
            original = out.read_bytes()
            (functions / "function_attributes.tsv").write_text(tables["function_attributes.tsv"].replace("not_tested", "supported"))
            with self.assertRaisesRegex(ValueError, "Claim producer output changed"):
                from_results(functions)
            self.assertEqual(out.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
