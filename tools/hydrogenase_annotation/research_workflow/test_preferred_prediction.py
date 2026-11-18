import csv
import json
import unittest
from pathlib import Path

from .preferred_prediction import PREDICTION_FIELDS, predict


def _row(**kw):
    base = {"family": "FeFe", "identity_status": "sequence_supported_uncalibrated_core",
            "coarse_function_prediction": "H2_core_supported_partner_unresolved",
            "supported_attributes": ["Hydrogen_reaction_core"]}
    base.update(kw)
    return base


def _pyr(**kw):
    value = {"hydrogen_reaction_core_status": "sequence_supported_candidate",
             "pyridine_nucleotide_module_status": "sequence_supported_candidate",
             "association_supported": True,
             "roles": {"HydA": [{"protein_id": "hydA"}]},
             "role_audit": {}}
    value.update(kw)
    return {"claims": {"pyridine": value}}


class PreferredPredictionTests(unittest.TestCase):
    def test_split_center_projects_identity_only_with_exact_two_producer_members(self):
        row = {"candidate_id": "split", "family": "NiFe",
               "identity_status": "sequence_supported_uncalibrated_split_center",
               "coarse_function_prediction": "NiFe_split_catalytic_center_supported_function_unresolved"}
        evidence = {"attributes": [{"attribute": "NiFe_split_catalytic_core_potential",
                                     "status": "sequence_supported_candidate",
                                     "value": {"members": [{"role": "VhuA", "protein_id": "a"},
                                                              {"role": "VhuU", "protein_id": "u"}]}}]}
        out = predict(row, evidence=evidence)
        self.assertEqual(out["preferred_function"],
                         "NiFe_split_catalytic_center_supported_function_unresolved")
        self.assertEqual({m["protein_id"] for m in out["preferred_members"]}, {"a", "u"})
        self.assertIn("split_catalytic_center_only", out["prediction_risks"])
        self.assertIn("electron_partner_unresolved", out["prediction_risks"])

    def test_split_center_without_proof_has_no_synthetic_members(self):
        out = predict({"candidate_id": "split", "family": "NiFe",
                       "identity_status": "sequence_supported_uncalibrated_split_center",
                       "coarse_function_prediction": "NiFe_split_catalytic_center_supported_function_unresolved"})
        self.assertEqual(out["preferred_function"], "unresolved_candidate")
        self.assertEqual(out["preferred_members"], [])

    def test_split_only_does_not_override_independent_legacy_function(self):
        out = predict({"candidate_id": "split", "family": "NiFe",
                       "identity_status": "sequence_supported_uncalibrated_split_center",
                       "coarse_function_prediction": "candidate_H2_formate_coupling"})
        self.assertEqual(out["preferred_function"], "candidate_H2_formate_coupling")

    def test_cooperative_support_is_explicit_without_upgrading_strict_status(self):
        evidence=_pyr(confidence_hint='low',composition_compensation={'applied':True,
                      'roles':['SHI_delta','SHI_gamma','SHI_beta'],'strict_core_status':'unresolved'})
        out=predict(_row(family='NiFe'),evidence=evidence)
        self.assertEqual(out['preferred_function'],'candidate_H2_pyridine_nucleotide_coupling')
        self.assertEqual(out['prediction_confidence'],'low')
        self.assertIn('邻近亚基组合证据',out['prediction_summary'])
        item=next(s for s in out['prediction_support'] if isinstance(s,dict) and s.get('claim')=='cooperative_role_compensation')
        self.assertEqual(item['strict_core_status'],'unresolved')
        self.assertEqual(item['roles'],['SHI_delta','SHI_gamma','SHI_beta'])

    def test_hhy_context_projects_function_and_preserves_low_confidence(self):
        row = _row(candidate_id='L', family='NiFe', identity_status='calibrated_high')
        context = {'status': 'candidate_H2_oxidation_potential', 'structure_id': 'NiFe_Hhy1h',
                   'confidence': 'low', 'roles': {'HhyL': [{'protein_id': 'L'}], 'HhyS': [{'protein_id': 'S'}]},
                   'risks': ['HhyS_reference_ligand_pattern_conflict'], 'blockers': []}
        result = predict(row, evidence={'claims': {'hhy_context': context}})
        self.assertEqual(result['preferred_function'], 'candidate_H2_oxidation_potential')
        self.assertEqual(result['prediction_confidence'], 'low')
        self.assertIn('参考Cys模式存在差异', result['prediction_summary'])
        self.assertNotIn('partner_or_electron_module_unresolved', result['prediction_risks'])
        self.assertEqual({m['protein_id'] for m in result['preferred_members']}, {'L', 'S'})
        self.assertEqual(result['experimental_unknowns']['in_vivo_direction'], 'unknown')
        context['status'] = 'not_supported'
        self.assertEqual(predict(row, evidence={'claims': {'hhy_context': context}})['preferred_function'], 'hydrogenase_core_candidate')

    def _artificial_mediator_case(self, **changes):
        row = _row(identity_status="sequence_supported_uncalibrated_core",
                   reference_reaction_status="exact_reference_sequence_potential", candidate_id="hydA")
        member = {"role": "HydA", "protein_id": "hydA"}
        cap = {"reference_evidence_scope": "purified_protein_artificial_mediator_reference_potential",
               "contract_id": "artificial-h2", "electron_donors": ["methyl_viologen_reduced"],
               "electron_acceptors": ["protons"], "hydrogen_role": "evolution"}
        ref = {"claim_name": "reference", "status": "exact_reference_sequence_potential",
               "matched_structure_id": "ref", "capabilities": [cap], "tuple_matches": [member],
               "tuple_hypotheses": [{"structure_id": "ref", "matched_members": [member], "blockers": []}]}
        row.update(changes)
        return row, {"claims": {"reference": ref}}

    def test_artificial_mediator_reference_potential_is_explicit_and_typed(self):
        row, evidence = self._artificial_mediator_case()
        out = predict(row, evidence=evidence)
        self.assertEqual(out["preferred_function"], "reference_sequence_matched_reaction_potential")
        self.assertIn("人工介质条件下的产H2", out["prediction_summary"])
        self.assertIn("甲基紫精还原态供电子产H2", out["prediction_summary"])
        self.assertIn("artificial_mediator_not_native_coupling", out["prediction_risks"])
        self.assertNotIn("NAD", out["prediction_summary"])
        self.assertNotIn("ferredoxin", out["prediction_summary"])

    def test_artificial_mediator_reference_requires_scope_typed_participants_and_identity(self):
        for change in ("scope", "typed", "identity"):
            row, evidence = self._artificial_mediator_case()
            cap = evidence["claims"]["reference"]["capabilities"][0]
            if change == "scope":
                cap.pop("reference_evidence_scope")
            elif change == "typed":
                cap["electron_acceptors"] = []
            else:
                row["identity_status"] = "unconfirmed"
            out = predict(row, evidence=evidence)
            self.assertNotEqual(out["preferred_function"], "reference_sequence_matched_reaction_potential", change)

    def test_artificial_mediator_does_not_create_nucleotide_or_ferredoxin_function(self):
        row, evidence = self._artificial_mediator_case()
        out = predict(row, evidence=evidence)
        self.assertNotIn(out["preferred_function"], {"candidate_H2_pyridine_nucleotide_coupling", "candidate_H2_ferredoxin_compatible_core"})

    def test_profile_supported_relay_keeps_explicit_low_confidence_function(self):
        claim={'status':'local_coupling_potential_uncalibrated','entry_route':'profile_supported_low_identity_role_complement',
               'confidence_hint':'low','selection_risks':['reference_contact_difference:C299->H'],
               'evidence':{'core_status':'sequence_supported_candidate','module_status':'sequence_supported_candidate',
                           'hypotheses':[{'roles':[{'role':role,'protein_id':pid} for role,pid in
                               [('SHI_alpha','a'),('SHI_delta','d'),('SHI_gamma','g'),('SHI_beta','b')]]}]}}
        out=predict(_row(family='NiFe'),evidence={'claims':{'composition':claim}})
        self.assertEqual(out['preferred_function'],'candidate_H2_pyridine_nucleotide_coupling')
        self.assertEqual(out['prediction_confidence'],'low');self.assertIn('低同一性中继',out['prediction_summary'])
        self.assertIn('reference_contact_difference:C299->H',out['prediction_risks'])
        self.assertEqual({m['protein_id'] for m in out['preferred_members']},{'a','d','g','b'})
    def test_rejected_fdh_partner_explanation_reaches_main_prediction(self):
        base=Path(__file__).resolve().parents[1]/'results/round149_after_Pa_FHL_native_2NAP_v1'
        claim=json.loads((base/'functions/fhl_function_evidence.json').read_text())['FHL92_003']
        out=predict(_row(family='NiFe'),evidence={'claims':{'fhl':claim}})
        self.assertEqual(out['preferred_function'],'hydrogenase_core_candidate')
        risk=next(r for r in out['prediction_risks'] if isinstance(r,dict) and r.get('code')=='preferred_non_FDH_role_competition')
        self.assertEqual(risk['protein_id'],'FHL92_011')
        self.assertEqual(risk['preferred_role_class'],'nitrate_reductase_counter')
        self.assertEqual(risk['scope'],'query_partner_interpretation_only')
        self.assertIn('其他酶功能解释',out['prediction_summary'])
        self.assertTrue(any(a.get('kind')=='preferred_other_partner_role' for a in out['prediction_alternatives']))
        claim['evidence']['excluded_FdhF_candidates'][0]['independent_role_competition']['status']='ambiguous'
        out=predict(_row(family='NiFe'),evidence={'claims':{'fhl':claim}})
        self.assertFalse(any(isinstance(r,dict) and r.get('code')=='preferred_non_FDH_role_competition' for r in out['prediction_risks']))
    def test_duplicate_relay_is_preferred_with_alternative(self):
        ev = _pyr(coupling_potential_status="unresolved",
                  coupling_blockers=["distributed_relay_pairing_unresolved_HydF"],
                  roles={"HydA": [{"protein_id": "hydA"}], "HydF": [{"protein_id": "hydf1"}]},
                  role_audit={"HydF": {"candidate_count": 2, "selected_protein_ids": ["hydf1"]}})
        out = predict(_row(), evidence=ev)
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["prediction_confidence"], "medium")
        self.assertTrue(out["prediction_alternatives"])
        self.assertEqual(out["prediction_confidence_scope"], "qualitative_evidence_not_calibrated_probability")

    def test_missing_module_does_not_force_coupling(self):
        out = predict(_row(), evidence=_pyr(pyridine_nucleotide_module_status="missing",
                                             module_association_supported=False,
                                             blockers=["missing_module_domain"]))
        self.assertEqual(out["preferred_function"], "hydrogenase_core_candidate")

    def test_remote_or_replacement_negative(self):
        out = predict(_row(), evidence=_pyr(coupling_blockers=["cross_contig_unsupported_link"]))
        self.assertEqual(out["preferred_function"], "hydrogenase_core_candidate")
        self.assertNotEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")

    def test_cross_anchor_competition_retains_low_preferred_call(self):
        out = predict(_row(), evidence=_pyr(coupling_blockers=["cross_anchor_competition"],
                                             alternative_architectures=[{"structure_id": "other", "status": "competing"}]))
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["prediction_confidence"], "low")

    def test_cross_anchor_only_is_not_same_resource_overlap(self):
        out = predict(_row(), evidence=_pyr(cross_anchor_conflicts=["remote_anchor"],
                                             resource_conflict=True))
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")

    def test_same_protein_overlap_blocks_coupling(self):
        out = predict(_row(), evidence=_pyr(resource_conflict=True,
                                             roles={"HydA": [{"protein_id": "p", "resources": [[1, 20], [10, 30]]}]}))
        self.assertNotEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")

    def test_shared_redox_domain_is_not_hydrogenase(self):
        out = predict({"identity_status": "unconfirmed", "supported_attributes": ["shared_redox_domain"]})
        self.assertEqual(out["prediction_kind"], "electron_module_only")
        self.assertNotIn("hydrogenase", out["preferred_function"])

    def test_shared_domain_does_not_erase_trusted_core(self):
        out = predict({"family": "NiFe", "identity_status": "calibrated_high",
                       "supported_attributes": ["Hydrogen_reaction_core", "electron_transfer_module_presence"]})
        self.assertEqual(out["preferred_function"], "hydrogenase_core_candidate")

    def test_rule_missing_is_not_a_function(self):
        out = predict({"family": "NiFe", "identity_status": "sequence_supported_uncalibrated_core",
                       "coarse_function_prediction": "candidate_function_rule_missing"})
        self.assertEqual(out["preferred_function"], "hydrogenase_core_candidate")

    def test_p00346_like_family_alone_is_not_possible_hydrogenase(self):
        out = predict({"family": "Nuo_related", "identity_status": "sequence_supported_uncalibrated_core",
                       "supported_attributes": ["electron_transfer_module_presence"]})
        self.assertEqual(out["prediction_kind"], "electron_module_only")

    def test_catalytic_profile_screen_can_yield_low_possible_family(self):
        out = predict({"identity_status": "unconfirmed", "supported_attributes": []},
                      evidence={"screen": {"family_candidates": ["PF02906"]}})
        self.assertEqual(out["preferred_function"], "possible_hydrogenase")
        self.assertEqual(out["preferred_family"], "FeFe")

    def _q06_selected_hypothesis(self, **changes):
        from pyridine_function import predict as evaluate,load_model
        root=Path(__file__).resolve().parents[1]
        claim=json.loads((root/"results/round101_Fdh2_release_v3/q06_duplicate_C/functions/pyridine_function_evidence.json").read_text())["q0004"]
        ev=claim["evidence"]
        ev.update(physical_association_supported=True,core_physical_association_supported=True,
                  module_physical_association_supported=True,actual_resource_overlap=False,
                  selected_core_anchor_present=True,role_assignment_ambiguous=True)
        ev.update(changes)
        return evaluate(ev,load_model())

    def test_q06_complete_selected_module_is_preferred_with_competition(self):
        claim=self._q06_selected_hypothesis()
        out=predict({"family":"FeFe","identity_status":"calibrated_high"},evidence={"claims":{"pyridine":claim}})
        self.assertEqual(claim["pyridine_nucleotide_module_status"],"unresolved")
        self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["prediction_confidence"],"medium")
        self.assertIn("partner_assignment_ambiguous",out["prediction_risks"])
        self.assertEqual({x["protein_id"] for x in out["preferred_members"]},{"q0001","q0002","q0003","q0004"})

    def test_selected_module_withdraws_for_missing_domain_or_invalid_association(self):
        for changes in ({"physical_association_supported":False},{"module_physical_association_supported":False},
                        {"actual_resource_overlap":True},{"module_competition_blockers":["target_reference_not_supported_HndC"]}):
            claim=self._q06_selected_hypothesis(**changes)
            self.assertNotIn("selected_hypothesis",claim)
            out=predict(_row(),evidence={"pyridine_result":claim})
            self.assertNotEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")

    def test_selected_hypothesis_retains_other_supported_function_and_ownership_risk(self):
        claim=self._q06_selected_hypothesis()
        out=predict(_row(coarse_function_prediction="candidate_H2_quinone_coupling"),evidence={"pyridine_result":claim})
        self.assertEqual(out["preferred_function"],"candidate_H2_quinone_coupling")
        claim["evidence"]["cross_anchor_conflicts"]=[{"rival_anchor":"other"}]
        out=predict(_row(),evidence={"pyridine_result":claim})
        self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["prediction_confidence"],"low")

    def test_selected_child_selection_risk_survives_preferred_projection(self):
        claim = self._q06_selected_hypothesis()
        child = dict(claim)
        child.pop("selected_hypothesis", None)
        child["coupling_potential_status"] = "local_coupling_potential_uncalibrated"
        child["pyridine_nucleotide_module_status"] = "sequence_supported_candidate"
        child["coupling_blockers"] = []
        child["selection_risks"] = ["child_specific_contact_risk"]
        claim["selected_hypothesis"] = {
            "strict_child": child
        }
        out = predict(_row(), evidence={"claims": {"pyridine": claim}})
        self.assertEqual(out["preferred_function"],
                         "candidate_H2_pyridine_nucleotide_coupling")
        self.assertIn("child_specific_contact_risk", out["prediction_risks"])

    def test_hmdii_is_not_canonical_h2(self):
        out = predict(_row(family="HmdII", coarse_function_prediction="candidate_H2_H4MPT_hydrogenation_potential"))
        self.assertEqual(out["preferred_family"], "HmdII")
        self.assertNotEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")

    def test_canonical_hmd_is_not_reclassified_as_hmdii(self):
        out = predict(_row(family="Hmd", coarse_function_prediction="candidate_H2_H4MPT_hydrogenation_potential"))
        self.assertEqual(out["preferred_function"], "candidate_H2_H4MPT_hydrogenation_potential")

    def test_existing_frh_and_hdr_claims_survive_core_fallback(self):
        for label in ("candidate_H2_F420_reduction", "candidate_H2_Hdr_coupling"):
            out = predict(_row(coarse_function_prediction=label, function_prediction=label), evidence={})
            self.assertEqual(out["preferred_function"], label)

    def test_coo_members_use_only_the_supported_substrate_pair(self):
        claim={"status":"candidate_CO_linked_H2_pathway","evidence":{
            "core":{"roles":{"CooH":"core","CooL":"relay"}},
            "substrate_modules":[{"roles":{"CooS":"substrate","CooF":"carrier"}}],
            "excluded_substrate_modules":[{"roles":{"CooS":"excluded"}}],
            "competing_anchor_ids":["rival"]}}
        row=_row(family="NiFe",coarse_function_prediction="candidate_CO_linked_H2_pathway")
        out=predict(row,evidence={"claims":{"coo":claim}})
        self.assertEqual({m["protein_id"] for m in out["preferred_members"]},
                         {"core","relay","substrate","carrier"})
        claim["evidence"]["substrate_modules"].append({"roles":{"CooS":"another"}})
        out=predict(row,evidence={"claims":{"coo":claim}})
        self.assertFalse(out["preferred_members"])

    def test_strong_legacy_quinone_remains_primary_with_composition_alternative(self):
        root = Path(__file__).resolve().parents[1]
        path = root / "results/round124_hadar_function_v1_retry4/u1055931/functions/pyridine_role_composition_evidence.json"
        composition = next(iter(json.loads(path.read_text()).values()))
        out = predict({"family": "NiFe", "coarse_function_prediction": "candidate_H2_quinone_coupling"},
                      evidence={"claims":{"composition":composition}})
        self.assertEqual(out["preferred_function"], "candidate_H2_quinone_coupling")
        self.assertTrue(out["prediction_alternatives"])

    def test_reference_assay_capabilities_can_support_high(self):
        out = predict(_row(reference_reaction_status="exact_reference_sequence_potential",
                           reference_assay_capabilities=[{"contract_id": "assay-1", "reaction_key": "NADH_H2_oxidation",
                                                          "electron_coupling": "pyridine_nucleotide_electron_transfer",
                                                          "electron_donors":["NADH"],"electron_acceptors":["protons"]}],
                           coarse_function_prediction="candidate_H2_pyridine_nucleotide_coupling"), evidence=_pyr())
        self.assertEqual(out["prediction_confidence"], "high")

    def _reference_case(self, case):
        root=Path(__file__).resolve().parents[1]
        targets=json.loads((root/"experiments/round135_function_regression/target_scoring_v1.json").read_text())["targets"]
        q=next(t["target_id"] for t in targets if t["case"]==case)
        path=root/"results/round139_function_regression_v3/prior30"/case/"functions"
        with (path/"functional_annotations.tsv").open() as f:
            row=next(r for r in csv.DictReader(f,delimiter="\t") if r["candidate_id"]==q)
        claims={name:json.loads((path/file).read_text())[q] for name,file in
                (("reference","reference_reaction_evidence.json"),("pyridine","pyridine_function_evidence.json"))}
        return row,{"claims":claims}

    def test_actual_tlitoralis_A_B_prefer_measured_reference_partner(self):
        for case in ("Tlitoralis_A_assayed_AF039208","Tlitoralis_B_whole_genome"):
            row,ev=self._reference_case(case);out=predict(row,evidence=ev)
            self.assertEqual(out["preferred_function"],"reference_sequence_matched_reaction_potential")
            self.assertEqual(out["prediction_confidence"],"high")
            self.assertIn("铁氧还蛋白",out["prediction_summary"])
            self.assertEqual({m["role"] for m in out["preferred_members"]},{"HydA","HydB","HydG","HydD"})
            self.assertEqual({m["protein_id"] for m in out["preferred_members"]},
                             {m["protein_id"] for m in ev["claims"]["reference"]["tuple_matches"]})
            self.assertTrue(any(a.get("preferred_function")=="candidate_H2_pyridine_nucleotide_coupling" and a["prediction_confidence"]=="medium" for a in out["prediction_alternatives"]))
            self.assertEqual(out["experimental_unknowns"]["query_activity"],"not_tested")

    def test_non_h2_source_reaction_retains_typed_electron_direction(self):
        from preferred_prediction import _reference_reaction_summary
        summary = _reference_reaction_summary([{
            "hydrogen_role": "not_a_net_H2_reaction", "electron_donors": ["formate"],
            "electron_acceptors": ["NAD+", "ferredoxin_oxidized"],
            "limitations": "NADP specificity not established",
        }])
        self.assertEqual(summary, "formate供电子，还原NAD+、Fd氧化态")
        self.assertNotIn("NADP", summary)

    def test_reference_summary_preserves_typed_coupled_assay_scope(self):
        from preferred_prediction import _reference_reaction_summary
        direct = {"hydrogen_role": "evolution", "electron_donors": ["NADPH"],
                  "electron_acceptors": ["protons"], "direct_acceptor_status": "reference_reported",
                  "limitations": "A coupled assay elsewhere must not relabel this direct reaction."}
        coupled = {**direct, "electron_donors": ["NADH"],
                   "direct_acceptor_status": "reference_reported_coupled"}
        self.assertEqual(_reference_reaction_summary([direct, coupled]),
                         "NADPH供电子产H2；NADH供电子产H2（耦联测定）")

    def test_equivalent_reference_member_loci_keep_function_and_expose_alternatives(self):
        root = Path(__file__).resolve().parents[1]
        path = root / "results/round145_KIST612_A_after_v2/functions"
        with (path / "functional_annotations.tsv").open() as handle:
            row = next(r for r in csv.DictReader(handle, delimiter="\t") if r["candidate_id"] == "ELI_0847")
        reference = json.loads((path / "reference_reaction_evidence.json").read_text())["ELI_0847"]
        pyridine = json.loads((path / "pyridine_function_evidence.json").read_text())["ELI_0847"]
        reference["member_selection_ambiguities"] = [{"role":"FdhA", "selected_protein_id":"ELI_0994",
            "alternative_protein_ids":["ELI_0994_dup"], "query_member_assay_status":"not_tested",
            "selection_basis":"source_equivalent_complete_members_ranked_by_context_then_id",
            "literal_member_assignment_status":"not_unique_sequence_equivalent_loci"}]
        out = predict(row, evidence={"claims":{"reference":reference, "pyridine":pyridine}})
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["prediction_confidence"], "medium")
        self.assertIn("中等置信度", out["prediction_summary"])
        self.assertNotIn("高置信度", out["prediction_summary"])
        self.assertIn("formate供电子，还原NAD+", out["prediction_summary"])
        self.assertIn("同序列备选", out["prediction_summary"])
        self.assertEqual(next(m["protein_id"] for m in out["preferred_members"] if m["role"]=="FdhA"), "ELI_0994")
        alt = next(a for a in out["prediction_alternatives"] if a.get("role")=="FdhA")
        self.assertEqual(alt["candidates"], [{"protein_ids":["ELI_0994"],"selected":True},
                                              {"protein_ids":["ELI_0994_dup"],"selected":False}])
        self.assertIn("sequence_equivalent_partner_locus_ambiguous", out["prediction_risks"])
        self.assertEqual(out["experimental_unknowns"]["query_activity"], "not_tested")
        self.assertEqual(out["prediction_confidence_scope"], "qualitative_evidence_not_calibrated_probability")

    def test_equivalent_member_loci_reference_override_is_medium(self):
        row, ev = self._reference_case("Tlitoralis_A_assayed_AF039208")
        reference = ev["claims"]["reference"]
        member = next(m for m in reference["tuple_matches"] if m["protein_id"]!=row["candidate_id"])
        reference["member_selection_ambiguities"] = [{"role":member["role"], "selected_protein_id":member["protein_id"],
            "alternative_protein_ids":["equivalent_copy"], "selection_basis":"context_then_id",
            "query_member_assay_status":"not_tested", "literal_member_assignment_status":"not_unique_sequence_equivalent_loci"}]
        out = predict(row, evidence=ev)
        self.assertEqual(out["preferred_function"], "reference_sequence_matched_reaction_potential")
        self.assertEqual(out["prediction_confidence"], "medium")
        self.assertIn("中等置信度", out["prediction_summary"])
        self.assertTrue(any(a.get("role")==member["role"] for a in out["prediction_alternatives"]))

    def _component_case(self):
        root=Path(__file__).resolve().parents[1]
        broken=root/"results/round146_KIST612_A_core_group_broken_after_v1/functions"
        with (broken/"functional_annotations.tsv").open() as h:
            row=next(r for r in csv.DictReader(h,delimiter="\t") if r["candidate_id"]=="ELI_0847")
        reference=json.loads((broken/"reference_reaction_evidence.json").read_text())["ELI_0847"]
        original=json.loads((root/"results/round145_KIST612_A_after_v2/functions/reference_reaction_evidence.json").read_text())["ELI_0847"]
        members=original["tuple_matches"]
        for member in members:
            if member["role"] in {"HydB","HydC"}:member["protein_id"]+="_moved"
        caps=original["capabilities"]
        for cap in caps:cap["scope"]="sequence_component_function_hypothesis"
        reference["sequence_component_hypotheses"]=[{
            "status":"source_equivalent_components_potential", "structure_id":original["matched_structure_id"],
            "matched_members":members, "capabilities":caps,
            "linkage_blockers":["reference_member_group_local_linkage_missing_0"], "hard_blockers":[],
            "strict_tuple_blockers":["missing_or_ambiguous_exact_reference_role_HydB"],
            "member_selection_ambiguities":[], "query_activity":"not_tested"}]
        return row,{"claims":{"reference":reference}}

    def test_complete_components_predict_function_without_exact_source_attribution(self):
        row,evidence=self._component_case();original=json.loads(json.dumps(row))
        out=predict(row,evidence=evidence)
        self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["prediction_confidence"],"medium")
        self.assertIn("完整来源同序列组件",out["prediction_summary"])
        self.assertIn("formate供电子，还原NAD+",out["prediction_summary"])
        self.assertEqual({m["role"] for m in out["preferred_members"]},{"HydA","HydB","HydC","FdhA"})
        self.assertIn("ELI_0846_moved",{m["protein_id"] for m in out["preferred_members"]})
        self.assertIn("source_assay_tuple_local_linkage_unresolved",out["prediction_risks"])
        self.assertEqual(out["experimental_unknowns"]["query_activity"],"not_tested")
        self.assertEqual(row,original)
        self.assertEqual(row["reference_reaction_status"],"not_supported")
        self.assertEqual(json.loads(row["reference_contract_ids"]),[])

        component=evidence["claims"]["reference"]["sequence_component_hypotheses"][0]
        component["linkage_blockers"]=[]
        component["member_selection_ambiguities"]=[{"role":"HydB", "selected_protein_id":"ELI_0846_moved", "alternative_protein_ids":["equivalent_HydB_copy"]}]
        out=predict(row,evidence=evidence)
        self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
        self.assertIn("source_assay_tuple_member_assignment_unresolved",out["prediction_risks"])
        self.assertTrue(any(a.get("role")=="HydB" for a in out["prediction_alternatives"]))
        component["member_selection_ambiguities"][0]["equally_supported_loci"]=["equivalent_HydB_copy"]
        out=predict(row,evidence=evidence)
        self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
        self.assertIn("partner_assignment_ambiguous",out["prediction_risks"])
        self.assertEqual({m["role"] for m in out["preferred_members"]},{"HydA","HydC","FdhA"})
        self.assertEqual(next(a for a in out["prediction_alternatives"] if a.get("role")=="HydB")["selected_protein_ids"],[])

    def test_component_inference_does_not_bypass_identity_or_hard_proof(self):
        for change in ("untrusted_core","hard_block","missing_anchor","non_linkage_error","untyped_reaction","wrong_scope"):
            with self.subTest(change=change):
                row,evidence=self._component_case();h=evidence["claims"]["reference"]["sequence_component_hypotheses"][0]
                if change=="untrusted_core":row["identity_status"]="identity_insufficient"
                elif change=="hard_block":h["hard_blockers"]=["member_resource_reuse_HydB"]
                elif change=="missing_anchor":h["matched_members"]=[m for m in h["matched_members"] if m["protein_id"]!=row["candidate_id"]]
                elif change=="non_linkage_error":h["linkage_blockers"]=["partial_or_unknown_member_HydB"]
                elif change=="wrong_scope":h["capabilities"][0]["scope"]="reference_assay_conditions_only"
                else:h["capabilities"]=[{"scope":"sequence_component_function_hypothesis","molecular_reaction":"NAD activity not measured"}]
                out=predict(row,evidence=evidence)
                self.assertFalse(any(isinstance(s,dict) and s.get("claim")=="source_equivalent_component_function" for s in out["prediction_support"]))

    def _source_like_case(self):
        root=Path(__file__).resolve().parents[1]
        path=root/"results/round147_KIST612_A_02_FdhA_sequence_variant_after_v1/functions"
        with (path/"functional_annotations.tsv").open() as handle:
            row=next(r for r in csv.DictReader(handle,delimiter="\t") if r["candidate_id"]=="ELI_0847")
        reference=json.loads((path/"reference_reaction_evidence.json").read_text())[row["candidate_id"]]
        original=json.loads((root/"results/round145_KIST612_A_after_v2/functions/reference_reaction_evidence.json").read_text())[row["candidate_id"]]
        members=original["tuple_matches"];fdh=next(m for m in members if m["role"]=="FdhA")
        fdh.update(protein_id="ELI_0994_variant",query_sha256="c9f61abd6742167337ee8aeb0875e7c24f6816167eafef561798c5b04e2fabb6",
                   matched_by="supported_FDH_role_homology",source_equivalent=False,paired_identity=99.8882681564,reference_coverage=1.0)
        reference["source_like_component_hypotheses"]=[{
            "status":"source_like_formate_component_potential","scope":"source_like_FDH_role_combination",
            "structure_id":original["matched_structure_id"],"matched_members":members,"hard_blockers":[],"linkage_risks":[],
            "supported_functions":["candidate_H2_formate_coupling"],"confidence":"medium","homologous_roles":["FdhA"],
            "reference_reaction_details_transferred":False,"member_selection_ambiguities":[]}]
        return row,{"claims":{"reference":reference}}

    def test_nonexact_FDH_partner_adds_coarse_function_without_source_reactions(self):
        row,evidence=self._source_like_case();original=json.loads(json.dumps(row))
        out=predict(row,evidence=evidence)
        self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
        self.assertIn("甲酸耦联潜力",out["prediction_summary"])
        self.assertNotIn("formate供电子产H2",out["prediction_summary"])
        self.assertNotIn("NAD+、Fd氧化态",out["prediction_summary"])
        fdh=next(m for m in out["preferred_members"] if m["role"]=="FdhA")
        self.assertEqual(fdh["protein_id"],"ELI_0994_variant")
        self.assertIs(fdh["source_equivalent"],False)
        support=next(s for s in out["prediction_support"] if isinstance(s,dict) and s.get("claim")=="source_like_FDH_role_combination")
        self.assertEqual(support["supported_functions"],["candidate_H2_formate_coupling"])
        self.assertIs(support["reference_reaction_details_transferred"],False)
        self.assertEqual(row,original)
        self.assertEqual(json.loads(row["reference_contract_ids"]),[])
        row.update(coarse_function_prediction="H2_core_supported_partner_unresolved",function_prediction="H2_core_supported_partner_unresolved")
        self.assertEqual(predict(row,evidence=evidence)["preferred_function"],"candidate_H2_formate_coupling")

    def _source_like_module_case(self, prefer=True, source_label="candidate_H2_pyridine_nucleotide_coupling"):
        row = _row(candidate_id="query_core", coarse_function_prediction="candidate_H2_pyridine_nucleotide_coupling")
        participants = [{"electron_donors": ["NADH"], "electron_acceptors": ["protons"]}]
        if source_label == "candidate_H2_formate_coupling":
            participants = ["formate", "protons"]
        members = [{"role": "HydA", "protein_id": "query_core", "source_equivalent": True},
                   {"role": "FdhA", "protein_id": "near_fdh", "matched_by": "supported_local_module_homology",
                    "source_equivalent": False, "query_sha256": "q" * 64, "reference_sha256": "r" * 64,
                    "source_alignment": {"identity": 99.1, "query_coverage": 1.0, "reference_coverage": .98,
                                         "query_sha256": "q" * 64, "reference_sha256": "r" * 64,
                                         "method": "local_module_homology", "alignments": [{"query": [1, 10], "reference": [1, 10]}]},
                    "source_alignment_dominates_role_template": prefer,
                    "role_module_evidence": {"complete": True, "role": "FdhA"}}]
        hypothesis = {"status": "source_like_module_component_potential", "scope": "source_like_module_function_prior",
                      "hard_blockers": [], "reference_reaction_details_transferred": False,
                      "source_contracts_transferred": False, "capabilities": [], "prefer_source_function": prefer,
                      "matched_members": members, "source_participants": participants,
                      "source_sources": ["primary-source"]}
        return row, {"claims": {"reference": {"source_like_component_hypotheses": [hypothesis]}}}

    def test_source_like_module_can_prefer_near_source_function_and_retain_NAD_alternative(self):
        row, evidence = self._source_like_module_case(prefer=True, source_label="candidate_H2_formate_coupling")
        out = predict(row, evidence=evidence)
        self.assertEqual(out["preferred_function"], "candidate_H2_formate_coupling")
        self.assertEqual(out["prediction_confidence"], "medium")
        self.assertTrue(any(a.get("preferred_function") == "candidate_H2_pyridine_nucleotide_coupling"
                            for a in out["prediction_alternatives"]))
        self.assertEqual(next(m for m in out["preferred_members"] if m["role"] == "FdhA")["source_alignment"]["query_coverage"], 1.0)
        self.assertIn("role_module_evidence", next(s for s in out["prediction_support"] if isinstance(s, dict)))
        self.assertNotIn("NADP", out["prediction_summary"])

    def test_source_like_module_low_confidence_is_not_raised_by_source_preference(self):
        row, evidence = self._source_like_module_case(
            prefer=True, source_label="candidate_H2_formate_coupling")
        hypothesis = evidence["claims"]["reference"]["source_like_component_hypotheses"][0]
        hypothesis["confidence"] = "low"
        hypothesis["selection_risks"] = ["reference_contact_difference:C299->H"]
        out = predict(row, evidence=evidence)
        self.assertEqual(out["preferred_function"], "candidate_H2_formate_coupling")
        self.assertEqual(out["prediction_confidence"], "low")
        self.assertIn("reference_contact_difference:C299->H", out["prediction_risks"])
        self.assertTrue(any(a.get("preferred_function") ==
                            "candidate_H2_pyridine_nucleotide_coupling"
                            for a in out["prediction_alternatives"]))

    def test_source_like_module_dominance_false_keeps_current_and_adds_source_alternative(self):
        row, evidence = self._source_like_module_case(prefer=False, source_label="candidate_H2_formate_coupling")
        out = predict(row, evidence=evidence)
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertTrue(any(a.get("preferred_function") == "candidate_H2_formate_coupling"
                            for a in out["prediction_alternatives"]))
        self.assertEqual(out["prediction_confidence"], "medium")

    def test_source_like_module_requires_scope_core_and_fresh_nonidentical_homolog(self):
        for change in ("scope", "hard", "core", "same", "typed", "dominance", "parent_hash"):
            row, evidence = self._source_like_module_case()
            hypothesis = evidence["claims"]["reference"]["source_like_component_hypotheses"][0]
            homolog = hypothesis["matched_members"][1]
            if change == "scope": hypothesis["scope"] = "source_like_FDH_role_combination"
            elif change == "hard": hypothesis["hard_blockers"] = ["missing_module"]
            elif change == "core": row["identity_status"] = "evidence_insufficient"
            elif change == "same": homolog["source_alignment"]["query_sha256"] = homolog["source_alignment"]["reference_sha256"]
            elif change == "dominance": homolog["source_alignment_dominates_role_template"] = not hypothesis["prefer_source_function"]
            elif change == "parent_hash": homolog["source_alignment"]["query_sha256"] = "0" * 64
            else: hypothesis["source_participants"] = ["NADP was measured"]
            out = predict(row, evidence=evidence)
            self.assertNotIn("source_like_module_component_function",
                             [s.get("claim") for s in out["prediction_support"] if isinstance(s, dict)], change)

    def test_source_like_FDH_inference_requires_own_homolog_and_scope(self):
        for change in ("hard_block","exact_protein","missing_anchor","details_transferred","wrong_scope"):
            with self.subTest(change=change):
                row,evidence=self._source_like_case();h=evidence["claims"]["reference"]["source_like_component_hypotheses"][0]
                if change=="hard_block":h["hard_blockers"]=["non_fdh_competition"]
                elif change=="exact_protein":next(m for m in h["matched_members"] if m["role"]=="FdhA")["source_equivalent"]=True
                elif change=="missing_anchor":h["matched_members"]=[m for m in h["matched_members"] if m["protein_id"]!=row["candidate_id"]]
                elif change=="details_transferred":h["reference_reaction_details_transferred"]=True
                else:h["scope"]="unreviewed"
                out=predict(row,evidence=evidence)
                self.assertFalse(any(isinstance(s,dict) and s.get("claim")=="source_like_FDH_role_combination" for s in out["prediction_support"]))

    def test_actual_ralbus_exact_same_coarse_reaction_is_propagated(self):
        root = Path(__file__).resolve().parents[1]
        path = root / "results/round144_Ralbus_A_after_v1/functions"
        with (path / "functional_annotations.tsv").open() as handle:
            row = next(r for r in csv.DictReader(handle, delimiter="\t") if r["candidate_id"] == "Rumal_2964")
        q = row["candidate_id"]
        ev = {"claims": {
            "reference": json.loads((path / "reference_reaction_evidence.json").read_text())[q],
            "pyridine": json.loads((path / "pyridine_function_evidence.json").read_text())[q],
        }}
        out = predict(row, evidence=ev)
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertIn("参考条件支持的具体反应潜力", out["prediction_summary"])
        self.assertTrue(any(x.get("claim") == "reference_reaction" for x in out["prediction_support"] if isinstance(x, dict)))
        self.assertEqual({m["protein_id"] for m in out["preferred_members"]},
                         {m["protein_id"] for m in ev["claims"]["reference"]["tuple_matches"]})

        pyridine = ev["claims"]["pyridine"]
        roles = pyridine.get("roles") or pyridine.get("evidence", {}).get("roles")
        roles["HydB"][0]["protein_id"] = "sequence_competing_HydB"
        changed = predict(row, evidence=ev)
        self.assertEqual(changed["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertTrue(any(any(m.get("protein_id") == "sequence_competing_HydB" for m in a.get("preferred_members", []))
                            for a in changed["prediction_alternatives"] if isinstance(a, dict)))

    def test_same_coarse_exact_reference_negative_bindings_do_not_propagate(self):
        root = Path(__file__).resolve().parents[1]
        path = root / "results/round144_Ralbus_A_after_v1/functions"
        with (path / "functional_annotations.tsv").open() as handle:
            row = next(r for r in csv.DictReader(handle, delimiter="\t") if r["candidate_id"] == "Rumal_2964")
        reference = json.loads((path / "reference_reaction_evidence.json").read_text())["Rumal_2964"]
        pyridine = json.loads((path / "pyridine_function_evidence.json").read_text())["Rumal_2964"]
        baseline = predict(row, evidence={"claims": {"pyridine": pyridine}})
        for change in ("status", "tuple", "prose"):
            ref = json.loads(json.dumps(reference))
            if change == "status":
                ref["status"] = "reference_associated_potential"
            elif change == "tuple":
                ref["tuple_matches"] = []
            else:
                ref["capabilities"] = [{"contract_id": "prose-only", "molecular_reaction": "NAD was not measured"}]
            out = predict(row, evidence={"claims": {"reference": ref, "pyridine": pyridine}})
            self.assertEqual(out["preferred_function"], baseline["preferred_function"], change)
            self.assertEqual(out["prediction_summary"], baseline["prediction_summary"], change)
            self.assertEqual(out["preferred_members"], baseline["preferred_members"], change)
            self.assertFalse(any(isinstance(x, dict) and x.get("claim") == "reference_reaction"
                                 for x in out["prediction_support"]), change)

    def test_unrelated_reference_does_not_raise_sequence_confidence(self):
        for caps in ([{"contract_id":"untyped"}], [{"electron_donors":["H2"],"electron_acceptors":["ferredoxin_oxidized"],
                     "molecular_reaction":"NAD was not measured; pyridine nucleotide specificity unknown"}]):
            out=predict(_row(reference_reaction_status="exact_reference_sequence_potential",reference_assay_capabilities=caps),evidence=_pyr())
            self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
            self.assertEqual(out["prediction_confidence"],"medium")

    def test_reference_requires_own_complete_tuple_and_exact_scope(self):
        for change in ("inherited","truncated_tuple","wrong_anchor","claim_unresolved"):
            row,ev=self._reference_case("Tlitoralis_A_assayed_AF039208");ref=ev["claims"]["reference"]
            if change=="inherited":row["reference_reaction_status"]=ref["status"]="reference_associated_potential"
            elif change=="truncated_tuple":ref["tuple_matches"].pop()
            elif change=="wrong_anchor":row["candidate_id"]="other_anchor"
            else:ref["status"]="not_supported"
            self.assertEqual(predict(row,evidence=ev)["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling",change)

    def test_reference_coding_locus_preserves_its_distinct_scope(self):
        row,ev=self._reference_case("Tlitoralis_A_assayed_AF039208");ref=ev["claims"]["reference"]
        row["reference_reaction_status"]=ref["status"]="exact_reference_coding_locus_potential"
        ref["tuple_matches"][0].update(raw_partial=True,coding_locus_proof_used=True)
        out=predict(row,evidence=ev)
        self.assertEqual(out["preferred_function"],"reference_coding_locus_matched_reaction_potential")
        self.assertIn("编码位点",out["prediction_summary"])
        self.assertTrue(out["preferred_members"][0]["source_partial"])
        self.assertTrue(out["preferred_members"][0]["coding_locus_proof_used"])

    def test_actual_tlit_variant_deleted_and_remote_members_do_not_inherit_assay(self):
        for case in ("Tlit_C_one_AA_variant","Tlit_C_remove_HydB","Tlit_C_relocated_HydB"):
            row,ev=self._reference_case(case);out=predict(row,evidence=ev)
            self.assertNotIn(out["preferred_function"],{"reference_sequence_matched_reaction_potential","reference_coding_locus_matched_reaction_potential"})
            self.assertNotEqual(out["prediction_confidence"],"high")
            if case=="Tlit_C_one_AA_variant":self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")
            else:self.assertEqual(out["preferred_function"],"hydrogenase_core_candidate")

    def test_experimental_unknowns_do_not_block_coarse_prediction(self):
        out = predict(_row(NAD_NADP_selectivity="unresolved", query_activity="not_tested"), evidence=_pyr())
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertEqual(out["experimental_unknowns"]["query_activity"], "not_tested")

    def test_partial_catalytic_profile_retains_possible_family_and_risk(self):
        out=predict({"family":"NiFe","identity_status":"evidence_insufficient"},evidence={
            "screen":{"raw_signal_evidence":[{"accession":"PF00374.25","passes":False,"model_coverage":.19}]},
            "gene":{"partial":True}})
        self.assertEqual(out["preferred_function"],"possible_hydrogenase")
        self.assertEqual(out["preferred_family"],"NiFe")
        self.assertIn("partial_gene_call",out["prediction_risks"])

    def test_real_role_intervals_overlap_even_with_cross_anchor_competition(self):
        evidence=_pyr(resource_conflict=True,cross_anchor_conflicts=[{"rival_anchor":"rival"}],
            roles={"HydA":[{"protein_id":"p","query_domain_start":1,"query_domain_end":400}],
                   "HydB":[{"protein_id":"p","query_domain_start":300,"query_domain_end":600}]})
        out=predict(_row(),evidence=evidence)
        self.assertNotEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")

    def test_incomplete_primary_can_use_atomic_context_limited_alternative(self):
        root = Path(__file__).resolve().parents[1]
        for name in ("q01_missing_B", "q03_partial_B"):
            path = root / f"results/round101_Fdh2_release_v3/{name}/functions/pyridine_function_evidence.json"
            claim = next(iter(json.loads(path.read_text()).values()))
            out = predict({"family": "FeFe", "coarse_function_prediction": "H2_core_supported_partner_unresolved"},
                          evidence={"pyridine_result": claim})
            self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
            self.assertEqual(out["prediction_confidence"], "low")
            self.assertTrue(any("context" in str(x).lower() for x in out["prediction_risks"]))
            self.assertTrue(any("FeFe_Sfructo_HndABCD" in str(x) for x in out["prediction_support"]))
            expected={"q0001","q0003","q0004"}|({"q0002"} if name=="q03_partial_B" else set())
            self.assertEqual({m["protein_id"] for m in out["preferred_members"]},expected)
            self.assertTrue(all("role" in m for m in out["preferred_members"]))
            if name=="q03_partial_B":self.assertTrue(next(m for m in out["preferred_members"] if m["role"]=="HndB")["source_partial"])

    def test_nested_alternatives_are_traversed_atomically(self):
        child = {"structure_id": "FeFe_nested", "hydrogen_reaction_core_status": "sequence_supported_candidate",
                 "pyridine_nucleotide_module_status": "sequence_supported_candidate",
                 "association_supported": True, "coupling_linkage_supported": True,
                 "roles": {"HydA": [{"protein_id": "nested"}]}}
        evidence = {"pyridine_result": {"structure_id": "FeFe_parent", "hydrogen_reaction_core_status": "unresolved",
                    "alternative_architectures": [{"structure_id": "level1", "alternative_architectures": [child]}]}}
        out = predict({"family": "FeFe"}, evidence=evidence)
        self.assertEqual(out["preferred_function"], "candidate_H2_pyridine_nucleotide_coupling")
        self.assertIn("FeFe_nested", " ".join(map(str, out["prediction_support"])))
        self.assertEqual(out["preferred_members"][0]["protein_id"], "nested")

    def test_actual_composition_supplies_core_relay_and_fusion_members(self):
        root=Path(__file__).resolve().parents[1]/"results/round138_function_regression_v2/hadar13/u1055931/functions"
        composition=json.loads((root/"pyridine_role_composition_evidence.json").read_text())
        q=next(q for q,c in composition.items() if c["status"]=="local_coupling_potential_uncalibrated")
        pyridine=json.loads((root/"pyridine_function_evidence.json").read_text())[q]
        out=predict(_row(family="NiFe",coarse_function_prediction="candidate_H2_pyridine_nucleotide_coupling"),
                    evidence={"claims":{"pyridine":pyridine,"composition":composition[q]}})
        members=out["preferred_members"]
        self.assertEqual(len(members),6)
        self.assertEqual(len({m["protein_id"] for m in members}),5)
        self.assertIn(q,{m["protein_id"] for m in members})
        fused=[m for m in members if m["protein_id"].endswith("_1117")]
        self.assertEqual(len(fused),2)
        self.assertLess(fused[0]["query_interval"][1],fused[1]["query_interval"][0])

    def test_alternatives_never_supply_each_others_missing_components(self):
        claim={"hydrogen_reaction_core_status":"unresolved","pyridine_nucleotide_module_status":"unresolved",
               "alternative_architectures":[
                   {"hydrogen_reaction_core_status":"sequence_supported_candidate","pyridine_nucleotide_module_status":"unresolved","association_supported":True},
                   {"hydrogen_reaction_core_status":"unresolved","pyridine_nucleotide_module_status":"sequence_supported_candidate","association_supported":True}]}
        out=predict(_row(),evidence={"claims":{"pyridine":claim}})
        self.assertNotEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")

    def test_actual_missing_electron_partner_and_intrinsic_relay_withdraw_coupling(self):
        root=Path(__file__).resolve().parents[1]/"results/round101_Fdh2_release_v3"
        for name in ("q02_missing_C","q04_remote_C","q05_cross_bin_C"):
            claim=next(iter(json.loads((root/name/"functions/pyridine_function_evidence.json").read_text()).values()))
            self.assertNotEqual(predict(_row(),evidence={"pyridine_result":claim})["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling",name)
        evidence=_pyr(coupling_blockers=["internal_FeFe_relay_organization_not_supported"])
        self.assertNotEqual(predict(_row(),evidence=evidence)["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")

    def test_function_level_status_selects_hdr_and_quinone_members(self):
        out=predict(_row(coarse_function_prediction="candidate_H2_Hdr_coupling"),evidence={"claims":{"mvh":{
            "coupling_status":"reference_like_H2_Hdr_coupling_hypothesis","roles":{"HdrA":[{"protein_id":"hdr"}]}}}})
        self.assertEqual(out["preferred_members"],[{"role":"HdrA","protein_id":"hdr"}])
        out=predict(_row(coarse_function_prediction="candidate_H2_quinone_coupling"),evidence={"claims":{
            "huc":{"coupling_potential_status":"unresolved","roles":{"wrong":[{"protein_id":"wrong"}]}},
            "quinone":{"coupling_potential_status":"local_quinone_coupling_potential_uncalibrated","roles":{"membrane":[{"protein_id":"right","query_domain_start":1,"query_domain_end":50}]},"role_audit":{"membrane":{"selected_protein_ids":["right"]}}}}})
        self.assertEqual(out["preferred_members"],[{"role":"membrane","protein_id":"right","query_interval":[1,50]}])

    def test_actual_hdr_members_and_capacity_alternatives_are_explicit(self):
        root=Path(__file__).resolve().parents[1]/"results/round138_function_regression_v2/hadar13"
        for case in ("u1564395","u147935","u667600"):
            data=json.loads((root/case/"functions/mvh_hdr_function_evidence.json").read_text())
            for q,claim in data.items():
                if claim["coupling_status"]!="reference_like_H2_Hdr_coupling_hypothesis":continue
                out=predict(_row(candidate_id=q,coarse_function_prediction="candidate_H2_Hdr_coupling"),evidence={"claims":{"mvh":claim}})
                members=out["preferred_members"]
                self.assertEqual({m["role"] for m in members},{"MvhA","MvhG","MvhD","HdrA","HdrB","HdrC"})
                self.assertIn(q,{m["protein_id"] for m in members})
                if case=="u147935":
                    self.assertEqual(next(m for m in members if m["role"]=="HdrA")["selection_basis"],"capacity_witness")
                    self.assertIn("HdrA_capacity_witness_not_physical_pairing",out["prediction_risks"])
                    self.assertTrue(any(a.get("role")=="HdrA" and len(a["candidates"])==3 for a in out["prediction_alternatives"]))
                if case=="u667600":
                    self.assertIn("competing_claim_not_resolved:"+q,out["prediction_risks"])
                    claim["evidence"]["global_function_assignment"]["blockers"].append("competing_claim_not_resolved:other_anchor")
                    out=predict(_row(candidate_id=q,coarse_function_prediction="candidate_H2_Hdr_coupling"),evidence={"claims":{"mvh":claim}})
                    self.assertNotIn("competing_claim_not_resolved:other_anchor",out["prediction_risks"])

    def test_single_protein_core_names_its_anchor(self):
        out=predict(_row(candidate_id="hmd",family="Hmd",coarse_function_prediction="candidate_H2_H4MPT_hydrogenation_potential"))
        self.assertEqual(out["preferred_members"],[{"role":"candidate_anchor","protein_id":"hmd"}])

    def test_uncalibrated_statistics_do_not_block_coarse_prediction(self):
        out=predict(_row(),evidence=_pyr(coupling_blockers=["outside_calibration_scope"],independent_validation="not_completed"))
        self.assertEqual(out["preferred_function"],"candidate_H2_pyridine_nucleotide_coupling")

    def test_fields_are_stable(self):
        self.assertEqual(set(PREDICTION_FIELDS), set(predict({})))

    def test_public_narf_negative_competitor_prefers_non_catalytic_shared_module(self):
        root = Path(__file__).resolve().parents[1]
        features = json.loads((root / "results/round142_NARF_before_v2/candidate_confirmation/candidate_annotation/confirmation/features.json").read_text())
        for feature in features:
            out = predict({"candidate_id": feature["protein_id"], "family": "FeFe", "identity_status": "unconfirmed"},
                          evidence={"identity_feature": feature})
            self.assertEqual(out["preferred_function"], "shared_redox_module_unresolved")
            self.assertEqual(out["prediction_kind"], "electron_module_only")
            self.assertEqual(out["prediction_confidence"], "medium")
            self.assertIn(feature["competitor_prototype"], str(out["prediction_support"]))

    def test_indeterminate_or_weak_negative_identity_does_not_trigger(self):
        base = {"competitor_margin": -20, "competitor_prototype": "p",
                "competitor_evidence": {"label": "negative", "prototype_id": "p", "core_coverage": .9, "evalue": 0}}
        for feature in (dict(base, competitor_evidence=dict(base["competitor_evidence"], label="indeterminate")), dict(base, competitor_evidence=dict(base["competitor_evidence"], core_coverage=.2)),
                        dict(base, competitor_evidence=dict(base["competitor_evidence"], evalue=1e-2))):
            out = predict({"family": "FeFe", "identity_status": "unconfirmed", "supported_attributes": ["H_cluster_redox_core_potential"]}, evidence={"identity_feature": feature})
            self.assertEqual(out["preferred_function"], "possible_hydrogenase")
        out=predict({"family":"FeFe","identity_status":"unconfirmed"},evidence={"identity_feature":dict(base,label="indeterminate")})
        self.assertEqual(out["preferred_function"],"shared_redox_module_unresolved")

    def test_negative_competitor_does_not_override_trusted_core(self):
        feature = {"competitor_margin": -20, "competitor_prototype": "p",
                   "competitor_evidence": {"label": "negative", "prototype_id": "p", "core_coverage": .9, "evalue": 0}}
        out = predict({"family": "FeFe", "identity_status": "calibrated_high", "supported_attributes": ["Hydrogen_reaction_core"]},
                      evidence={"identity_feature": feature})
        self.assertEqual(out["preferred_function"], "hydrogenase_core_candidate")


class LocalPartnerPreferenceTests(unittest.TestCase):
    def case(self,name):
        root=Path(__file__).resolve().parents[1]
        folder=root/f'results/round164_{name}_before_v1/functions'
        with (folder/'functional_annotations.tsv').open() as handle:
            row=next(r for r in csv.DictReader(handle,delimiter='\t') if r['candidate_id']=='qA')
        claims={key:json.loads((folder/file).read_text())['qA'] for key,file in (
            ('pyridine','pyridine_function_evidence.json'),('composition','pyridine_role_composition_evidence.json'),
            ('reference','reference_reaction_evidence.json'))}
        return row,dict(claims=claims)

    def test_local_supported_members_survive_remote_exact_components(self):
        row,evidence=self.case('competition');out=predict(row,evidence=evidence)
        self.assertEqual(out['preferred_function'],'candidate_H2_pyridine_nucleotide_coupling')
        self.assertEqual({m['protein_id'] for m in out['preferred_members']},{'qA','qL'})
        self.assertTrue(any({m['protein_id'] for m in a.get('preferred_members',[])}=={'qA','qR'}
                            for a in out['prediction_alternatives']))

    def test_remote_only_still_has_an_explicit_preferred_function(self):
        row,evidence=self.case('remote_exact');out=predict(row,evidence=evidence)
        self.assertEqual(out['preferred_function'],'candidate_H2_pyridine_nucleotide_coupling')
        self.assertEqual({m['protein_id'] for m in out['preferred_members']},{'qA','qR'})
        self.assertIn('source_assay_tuple_local_linkage_unresolved',out['prediction_risks'])

    def test_local_only_is_preserved(self):
        row,evidence=self.case('local_homolog');out=predict(row,evidence=evidence)
        self.assertEqual({m['protein_id'] for m in out['preferred_members']},{'qA','qL'})

    def test_local_confidence_is_not_raised_by_an_alternative_component_tuple(self):
        row,evidence=self.case('competition');evidence['claims']['composition']['confidence_hint']='low'
        out=predict(row,evidence=evidence)
        self.assertEqual(out['prediction_confidence'],'low')
        self.assertEqual({m['protein_id'] for m in out['preferred_members']},{'qA','qL'})

    def test_same_typed_member_ids_do_not_erase_low_confidence(self):
        row,evidence=self.case('competition');evidence['claims']['composition']['confidence_hint']='low'
        hypothesis=evidence['claims']['reference']['sequence_component_hypotheses'][0]
        for member in hypothesis['matched_members']:
            if member['protein_id']=='qR':member['protein_id']='qL'
        out=predict(row,evidence=evidence)
        self.assertEqual(out['prediction_confidence'],'low')
        self.assertEqual({m['protein_id'] for m in out['preferred_members']},{'qA','qL'})

    def test_unproven_local_member_list_cannot_disable_component_prediction(self):
        from preferred_prediction import _component_prediction
        row,evidence=self.case('competition')
        reference={**evidence['claims']['reference'],'claim_name':'reference'}
        fake=dict(preferred_function='candidate_H2_pyridine_nucleotide_coupling',preferred_family='FeFe',
                  preferred_members=[dict(role='HydA',protein_id='qA'),dict(role='HydB',protein_id='qL')],
                  prediction_support=[],prediction_risks=[],prediction_alternatives=[],prediction_confidence='low')
        out=_component_prediction(row,[reference],fake)
        self.assertEqual({m['protein_id'] for m in out['preferred_members']},{'qA','qR'})


class MembraneContextPredictionTests(unittest.TestCase):
    def test_group4_core_module_is_low_coarse_with_explicit_alternative(self):
        claim=dict(status='local_coupling_potential_uncalibrated',entry_route='independent_Group4_core_pyridine_module_composition',
                   confidence_hint='low',blockers=[],selection_risks=['reference_contact_difference:C299->H'],
                   evidence=dict(core_status='sequence_supported_candidate',module_status='sequence_supported_candidate',strict_core_status='unresolved',
                     hypotheses=[dict(roles=[dict(role=role,protein_id=role) for role in ('NiFe_catalytic_core','Group4_FeS_relay','SHI_beta','SHI_gamma')])]))
        out=predict(dict(candidate_id='NiFe_catalytic_core',family='NiFe',identity_status='calibrated_high'),evidence=dict(claims=dict(composition=claim)))
        self.assertEqual(out['preferred_function'],'candidate_H2_pyridine_nucleotide_coupling');self.assertEqual(out['prediction_confidence'],'low')
        self.assertEqual(len(out['preferred_members']),4);self.assertIn('核心连接架构未定',out['prediction_summary'])
        self.assertTrue(any(a.get('preferred_function')=='hydrogenase_core_candidate' for a in out['prediction_alternatives']))
        self.assertIn('reference_contact_difference:C299->H',out['prediction_risks'])

    def test_coarse_context_uses_its_own_members_and_low_confidence(self):
        child=dict(membrane_redox_status='local_membrane_H2_redox_potential_uncalibrated',confidence_hint='low',
                   selection_risks=['exact_Group4_architecture_unresolved'],evidence=dict(roles={r:[dict(protein_id=r)]
                    for r in ('H2_core','proximal_FeS_small','FeS_relay','membrane_antiporter','membrane_anchor')}))
        row=dict(candidate_id='H2_core',family='NiFe',identity_status='sequence_supported_uncalibrated_core',
                 coarse_function_prediction='candidate_membrane_associated_H2_redox_system')
        result=predict(row,evidence=dict(claims=dict(ech=dict(membrane_redox_status='unresolved'),membrane_context=child)))
        self.assertEqual(result['prediction_confidence'],'low')
        self.assertEqual(len(result['preferred_members']),5)
        self.assertIn('typed_claim:membrane_context',result['prediction_support'])
        self.assertIn('exact_Group4_architecture_unresolved',result['prediction_risks'])


class SingleChainReferenceFunctionTests(unittest.TestCase):
    def case(self):
        from explain_candidates import preferred_role_options
        root=Path(__file__).resolve().parents[1]
        folder=root/'results/round165_fixed_regression_v1/functions/regression__prior30__CpIII_B_whole_genome'
        annotation=root/'results/round117_function_conversion_v4/CpIII_B_whole_genome/annotation'
        def rows(path):
            with path.open() as h:return list(csv.DictReader(h,delimiter='\t'))
        q='CP009268.1_21';row=next(r for r in rows(folder/'functional_annotations.tsv') if r['candidate_id']==q)
        ref=json.loads((folder/'reference_reaction_evidence.json').read_text())[q]
        evidence=dict(claims=dict(reference=ref),attributes=[a for a in rows(folder/'function_attributes.tsv') if a['candidate_id']==q],
                      module_graph_nodes=preferred_role_options(rows(annotation/'module_graph_nodes.tsv'),q,'FeFe_CpI_monomer'))
        return row,evidence

    def test_actual_native_CpI_gets_coarse_Fd_and_keeps_assay_scope(self):
        row,evidence=self.case();original=json.dumps(evidence,sort_keys=True);out=predict(row,evidence=evidence)
        self.assertEqual(out['preferred_function'],'candidate_H2_ferredoxin_compatible_core')
        self.assertEqual(out['prediction_confidence'],'medium')
        self.assertEqual(out['preferred_members'],[dict(role='CpI',protein_id=row['candidate_id'])])
        self.assertEqual(out['experimental_unknowns']['query_activity'],'not_tested')
        self.assertIn('native_ferredoxin_partner_unidentified',out['prediction_risks'])
        self.assertTrue(any(a['preferred_function']=='hydrogenase_core_candidate' for a in out['prediction_alternatives']))
        self.assertEqual(row['reference_reaction_status'],'inherited_reference_potential')
        self.assertEqual(json.dumps(evidence,sort_keys=True),original)

    def test_identity_and_assigned_typed_source_reaction_are_required(self):
        for bad in ('identity','status','assignment','prose_only','NAD_coupled'):
            row,e=self.case();ref=e['claims']['reference']
            if bad=='identity':row['identity_status']='unconfirmed'
            if bad=='status':ref['status']='not_supported'
            if bad=='assignment':ref['inherited_assignments']=[]
            if bad=='prose_only':ref['capabilities'][0]['electron_donors']=[]
            if bad=='NAD_coupled':ref['capabilities'][0]['electron_coupling']='coupled_NAD_Fd'
            self.assertNotEqual(predict(row,evidence=e)['preferred_function'],'candidate_H2_ferredoxin_compatible_core',bad)

    def test_complete_single_chain_member_evidence_is_required(self):
        for bad in ('architecture','structure_binding','missing_nodes','extra_role','wrong_member','fragmented','not_jointly_feasible'):
            row,e=self.case();nodes=e['module_graph_nodes']
            if bad in {'architecture','structure_binding'}:
                a=next(a for a in e['attributes'] if a['attribute']=='candidate_complex_architecture')
                a['status' if bad=='architecture' else 'value']='unresolved'
            if bad=='missing_nodes':e['module_graph_nodes']=[]
            if bad=='extra_role':nodes.append(dict(role='external_partner',selected_protein_ids=['other'],candidates=[]))
            if bad=='wrong_member':nodes[0]['selected_protein_ids']=['other']
            if bad=='fragmented':nodes[0]['candidates'][0]['fragmented']=True
            if bad=='not_jointly_feasible':nodes[0]['candidates'][0]['jointly_feasible']=False
            self.assertEqual(predict(row,evidence=e)['preferred_function'],'hydrogenase_core_candidate',bad)

    def test_unrelated_unresolved_template_does_not_block_current_source(self):
        row,e=self.case();e['claims']['experimental_backbone']=dict(status='unresolved',blockers=['other_template_coverage_low'])
        self.assertEqual(predict(row,evidence=e)['preferred_function'],'candidate_H2_ferredoxin_compatible_core')

    def test_existing_specific_function_is_retained(self):
        row,e=self.case();row['coarse_function_prediction']='candidate_H2_pyridine_nucleotide_coupling'
        self.assertEqual(predict(row,evidence=e)['preferred_function'],'candidate_H2_pyridine_nucleotide_coupling')

    def test_low_confidence_is_not_raised_by_inherited_source(self):
        from .preferred_prediction import _single_chain_reference_prediction,_claims
        row,e=self.case();base=predict(row,evidence={});base['prediction_confidence']='low'
        self.assertEqual(_single_chain_reference_prediction(row,_claims(e),base,e)['prediction_confidence'],'low')


if __name__ == "__main__":
    unittest.main()
