"""Contract tests for the source-pinned reference-reaction projection.

These tests combine the native complex-extension bank with selected public
``results`` fixtures.  They exercise the boundary between identity/pyridine
evidence and the typed reference-reaction registry; they are not tests of the
legacy catalytic classifier.
"""

import csv
import hashlib
import json
import shutil
import tempfile
import unittest
import unittest.mock
from copy import deepcopy
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from coding_completeness import _record_digest, load_evidence as coding_load_evidence, prepare as coding_prepare


ROOT = Path(__file__).resolve().parents[1]
from reference_reactions import (assay_bank, prepare, tuple_proof, _coding_proof_valid,
                                 coding_locus_candidate_ids,
                                 _validate_admitted_bank,
                                 _validate_model_scope, _validate_member_groups,
                                 _source_module_prior, VERSION)  # noqa: E402


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _tsv(path, rows):
    rows = list(rows)
    fields = sorted({key for row in rows for key in row}) or ["_empty"]
    with Path(path).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


class ReferenceReactionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refs = json.loads((ROOT / "references/processed/complex_extensions/references.json").read_text())
        cls.refs_by_structure = {}
        for ref in cls.refs:
            cls.refs_by_structure.setdefault(ref["structure_id"], []).append(ref)
        cls.contracts = json.loads((ROOT / "config/reaction_contracts.json").read_text())["contracts"]
        cls.contract_by_id = {row["contract_id"]: row for row in cls.contracts}

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.annotation = Path(self.tmp.name) / "annotation"
        self.annotation.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def _fixture(self, structure, *, input_mode="genome", unknown_length=False,
                 mutate=None, partial=(), duplicate_role=None, bad_length=False,
                 cross_bin=(), identity_status="sequence_supported_uncalibrated_core",
                 anchor_status="candidate", potential_ids=None, labels=None):
        refs = self.refs_by_structure[structure]
        ids = {r["role"]: f"q_{i:02d}" for i, r in enumerate(refs, 1)}
        labels = labels or {}
        sequences = {ids[r["role"]]: r["sequence"] for r in refs}
        if mutate:
            qid = ids[mutate]
            seq = sequences[qid]
            sequences[qid] = seq[:-1] + ("A" if seq[-1] != "A" else "G")

        genes = []
        members = []
        for ref in refs:
            qid, role = ids[ref["role"]], ref["role"]
            contig = ref["contig"]
            start, end = ref["start"], ref["end"]
            length = 0 if unknown_length else max(end + 100, 2_000_000)
            if bad_length and role == refs[0]["role"]:
                length = end - 1
            bin_id = "b2" if role in set(cross_bin) else "b1"
            genes.append({
                "id": qid, "contig": contig, "start": str(start), "end": str(end),
                "strand": str(ref["strand"]), "gene": labels.get(role, role),
                "aliases": labels.get(role, ""), "protein_accession": ref["protein_accession"],
                "partial": "True" if role in set(partial) else "False",
                "contig_length": str(length), "bin_id": bin_id,
                "coordinate_evidence": "source_reference_coordinates",
            })
            members.append({
                "anchor": ids[refs[0]["role"]], "structure_id": structure, "role": role,
                "protein_id": qid, "reference_protein": ref["protein_accession"],
                "required_for_structure": "True", "scope": "same_contig",
                "bin_id": bin_id, "contig": contig, "start": str(start), "end": str(end),
                "strand": str(ref["strand"]), "identity": "100", "query_coverage": "1",
                "reference_coverage": "1", "reference_cysteines_conserved": "True",
                "cross_anchor_pairing_status": "no_overlapping_fixed_member_claim_detected",
            })
        duplicate_qid = None
        if duplicate_role:
            row = next(x for x in members if x["role"] == duplicate_role)
            duplicate_qid = "q_duplicate"
            sequences[duplicate_qid] = sequences[row["protein_id"]]
            genes.append(dict(next(g for g in genes if g["id"] == row["protein_id"]), id=duplicate_qid,
                              gene=duplicate_role, protein_accession="duplicate"))
            members.append(dict(row, protein_id=duplicate_qid, required_for_structure="True"))

        anchor_qid = ids[refs[0]["role"]]
        anchor = {
            "protein_id": anchor_qid, "structure_id": structure, "input_mode": input_mode,
            "bin_id": "b1", "catalytic_identity_confidence": "calibrated_high",
            "catalytic_status": anchor_status, "structure_status": "eligible",
            "active_center_status": "supported", "active_center_span_status": "supported",
            "module_organization_status": "supported", "module_graph_search_status": "supported",
            "cross_anchor_pairing_status": "no_overlapping_fixed_member_claim_detected",
            "potential_function_ids": ";".join(potential_ids or []),
            "function_status": "reference_sequence_supported_potential" if potential_ids else "unassigned", "flags": "",
        }
        self.annotation.joinpath("genes.tsv").parent.mkdir(exist_ok=True)
        _tsv(self.annotation / "genes.tsv", genes)
        _tsv(self.annotation / "complex_members.tsv", members)
        _tsv(self.annotation / "catalytic_anchors.tsv", [anchor])
        with (self.annotation / "proteins.faa").open("w") as fh:
            for qid, seq in sequences.items():
                fh.write(f">{qid}\n{seq}\n")

        roles = {}
        for ref in refs:
            roles.setdefault(ref["role"], []).append({"protein_id": ids[ref["role"]]})
        if duplicate_qid:
            roles.setdefault(duplicate_role, []).append({"protein_id": duplicate_qid})
        pyridine = {anchor_qid: {
            "structure_id": structure,
            "coupling_potential_status": "local_coupling_potential_uncalibrated",
            "evidence": {"roles": roles},
        }}
        identities = {anchor_qid: identity_status}
        self._manifest(input_mode=input_mode)
        return [anchor], identities, pyridine, ids

    def _manifest(self, *, tamper=False, input_mode="genome"):
        files = {name: _sha(self.annotation / name) for name in (
            "proteins.faa", "genes.tsv", "complex_members.tsv", "catalytic_anchors.tsv")}
        for name in ("reaction_capabilities.tsv", "reaction_contract_candidates.tsv"):
            if (self.annotation / name).exists():
                files[name] = _sha(self.annotation / name)
        manifest = {
            "status": "complete", "input_mode": input_mode,
            "thresholds": json.loads((ROOT / "config/thresholds.json").read_text()),
            "output_sha256": dict(files),
            "files": dict(files),
            "normalized_proteins_sha256": files["proteins.faa"],
        }
        if tamper:
            manifest["output_sha256"]["genes.tsv"] = "0" * 64
        (self.annotation / "run_manifest.json").write_text(json.dumps(manifest))

    def _run(self, structure, **kwargs):
        anchors, identities, pyridine, _ = self._fixture(structure, **kwargs)
        return prepare(self.annotation, anchors, identities, pyridine, root=ROOT)

    def _result(self, output, qid):
        results = output[0] if isinstance(output, tuple) else output
        return results[qid]

    def _tlit_public_inputs(self, case, qid):
        """Build public prepare inputs from the frozen R118 outputs."""
        if case == "Tlit_C_one_AA_variant":
            base = ROOT / "results/round118_reference_function_v2" / case
            pyridine_path = (ROOT / "results/round151_regression_v2/regression/prior30"
                             / case / "functions/pyridine_function_evidence.json")
        else:
            base = ROOT / "results/round118_Tlitoralis_transfer_v1" / case
            pyridine_path = base / "functions/pyridine_function_evidence.json"
        annotation = base / "annotation"
        evidence = json.loads(pyridine_path.read_text())[qid]
        role_map = {"SHI_alpha": "HydA", "SHI_delta": "HydD",
                    "SHI_gamma": "HydG", "SHI_beta": "HydB"}
        roles = {role_map[name]: [dict(member) for member in members]
                 for name, members in evidence["evidence"]["roles"].items()}
        pyridine = {qid: {
            "structure_id": "NiFe_Tlitoralis_HydBGDA",
            "hydrogen_reaction_core_status": evidence["hydrogen_reaction_core_status"],
            "coupling_potential_status": evidence["coupling_potential_status"],
            "evidence": {"roles": roles, "component_roles": roles,
                         "association_supported": True},
        }}
        with (annotation / "catalytic_anchors.tsv").open() as handle:
            anchors = list(csv.DictReader(handle, delimiter="\t"))
        if qid not in {row["protein_id"] for row in anchors}:
            # The native genome anchor is supplied by the exact HydA locus,
            # while the legacy catalytic table's nearest anchor has another ID.
            anchor = next(row for row in anchors if row["protein_id"] == "CP006670.1_1688")
        else:
            anchor = next(row for row in anchors if row["protein_id"] == qid)
        anchor = dict(anchor, protein_id=qid, structure_id="NiFe_Tlitoralis_HydBGDA",
                      input_mode="genome" if "whole_genome" in case else "assembly")
        return annotation, [anchor], {qid: "sequence_supported_uncalibrated_core"}, pyridine

    def test_tlitoralis_assay_bank_promotes_actual_A_and_B_without_identity_mutation(self):
        for case, qid in (("Tlitoralis_A_assayed_AF039208", "a01"),
                          ("Tlitoralis_B_whole_genome", "CP006670.1_668")):
            with self.subTest(case=case):
                annotation, anchors, identities, pyridine = self._tlit_public_inputs(case, qid)
                before = deepcopy(identities)
                results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
                row = results[qid]
                self.assertEqual(row["status"], "exact_reference_sequence_potential")
                self.assertEqual({cap["contract_id"] for cap in row["capabilities"]}, {
                    "Ref_NiFe_Tlitoralis_HydBGDA_oxidation",
                    "Ref_NiFe_Tlitoralis_HydBGDA_evolution"})
                self.assertEqual(identities, before)

    def test_tlit_variant_gets_source_like_fd_prior_without_exact_contract_transfer(self):
        annotation, anchors, identities, pyridine = self._tlit_public_inputs(
            "Tlit_C_one_AA_variant", "a01")
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        row = results["a01"]
        self.assertNotIn(row["status"], {"exact_reference_sequence_potential",
                                          "exact_reference_coding_locus_potential"})
        hypotheses = [h for h in row["source_like_component_hypotheses"]
                      if h["status"] == "source_like_module_component_potential"]
        self.assertEqual(len(hypotheses), 1)
        h = hypotheses[0]
        self.assertTrue(h["prefer_source_function"])
        self.assertAlmostEqual(h["matched_members"][-1]["source_alignment"]["identity"],
                               99.7, delta=0.5)
        self.assertAlmostEqual(h["matched_members"][-1]["source_alignment"]["query_coverage"], 1.0)
        self.assertAlmostEqual(h["matched_members"][-1]["source_alignment"]["reference_coverage"], 1.0)
        self.assertEqual(h["capabilities"], [])
        self.assertFalse(h["source_contracts_transferred"])
        self.assertFalse(h["reference_reaction_details_transferred"])
        self.assertEqual(row["capabilities"], [])

    def test_measured_Frh_Fdh3B_components_require_all_four_members_and_input_ownership(self):
        sid='NiFe_Tonnurineus_FrhAGB_Fdh3B'
        cid='Ref_NiFe_Tonnurineus_FrhAGB_Fdh3B_NADP_reduction'
        for case in ['intact_A','delete_Fdh3B','replace_Fdh3B','crossbin_Fdh3B']:
            folder=ROOT/f'results/round154_before_{case}_v1'
            with (folder/'annotation/catalytic_anchors.tsv').open() as h:
                anchors=list(csv.DictReader(h,delimiter='\t'))
            with (folder/'candidate_function_catalog.tsv').open() as h:
                identities={r['candidate_id']:r['identity_status'] for r in csv.DictReader(h,delimiter='\t')}
            results,_=prepare(folder/'annotation',anchors,identities,{},root=ROOT)
            result=results['FrhA']
            hypotheses=[h for h in result.get('sequence_component_hypotheses',[]) if h['structure_id']==sid]
            self.assertFalse(any(c['contract_id']==cid for c in result.get('capabilities',[])))
            if case=='intact_A':
                self.assertEqual(len(hypotheses),1)
                hypothesis=hypotheses[0]
                self.assertEqual(hypothesis['status'],'source_equivalent_components_potential')
                self.assertEqual(hypothesis['hard_blockers'],[])
                self.assertEqual({m['protein_id'] for m in hypothesis['matched_members']},{'FrhA','FrhG','FrhB','Fdh3B'})
                self.assertEqual(hypothesis['linkage_blockers'],['whole_reference_tuple_local_linkage_missing'])
                self.assertEqual([c['contract_id'] for c in hypothesis['capabilities']],[cid])
                self.assertEqual(hypothesis['capabilities'][0]['scope'],'sequence_component_function_hypothesis')
            else:
                self.assertFalse(any(h['status']=='source_equivalent_components_potential' for h in hypotheses),case)
                self.assertTrue(all(not h['capabilities'] for h in hypotheses),case)

    def test_source_module_prior_requires_complete_nonmodule_backbone_and_clean_pyridine_claim(self):
        annotation, anchors, identities, pyridine = self._tlit_public_inputs(
            "Tlit_C_one_AA_variant", "a01")
        baseline, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        baseline_h = next(h for h in baseline["a01"]["source_like_component_hypotheses"]
                          if h["status"] == "source_like_module_component_potential")

        def invoke(result=None, genes_override=None, sequences_override=None, current=None):
            refs, systems, _ = assay_bank(json.loads(
                (ROOT / "research_workflow/reference_reactions_model.json").read_text()), ROOT)
            rs = [r for r in refs if r["structure_id"] == "NiFe_Tlitoralis_HydBGDA"]
            component_roles = {"HydA": [{"protein_id": "a01"}],
                               "HydD": [{"protein_id": "a02"}],
                               "HydG": [{"protein_id": "a03"}],
                               "HydB": [{"protein_id": "a04"}]}
            with (annotation / "genes.tsv").open() as handle:
                genes_rows = list(csv.DictReader(handle, delimiter="\t"))
            gs = {g["id"]: g for g in genes_rows}
            for g in gs.values():
                for k in ("start", "end", "strand", "contig_length"):
                    g[k] = int(g[k])
            with (annotation / "proteins.faa").open() as handle:
                seqs = {r.id: str(r.seq) for r in SeqIO.parse(handle, "fasta")}
            if genes_override:
                genes_override(gs)
            if sequences_override:
                sequences_override(seqs)
            rr = result or {"entry_route": "exact_full_protein_lookup",
                            "evidence": {"roles": component_roles,
                                         "resource_conflict": False}}
            current = current or pyridine["a01"]
            current_map = {"a01": deepcopy(current)}
            return _source_module_prior("a01", rr, rs, component_roles, gs, seqs,
                                        "assembly", {"min_identity": 40,
                                                      "min_query_coverage": .8,
                                                      "min_reference_coverage": .8},
                                        current_map,
                                        [c for c in self.contracts if c["structure_id"] == "NiFe_Tlitoralis_HydBGDA"])

        # A broken backbone, resource conflict, or cross-bin ownership cannot
        # be repaired by the near-source module prior.
        for mutate in (lambda r: r["evidence"].update(resource_conflict=True),
                       lambda r: r["evidence"].update(actual_resource_overlap=True),
                       lambda r: r["evidence"]["roles"].update(HydG=[{"protein_id": "a04"}])):
            result = {"entry_route": "exact_full_protein_lookup",
                      "evidence": {"roles": {"HydA": [{"protein_id": "a01"}],
                                                "HydD": [{"protein_id": "a02"}],
                                                "HydG": [{"protein_id": "a03"}],
                                                "HydB": [{"protein_id": "a04"}]},
                                  "resource_conflict": False}}
            mutate(result)
            self.assertEqual(invoke(result=result), [])
        self.assertEqual(invoke(current=dict(pyridine["a01"], coupling_blockers=["x"])), [])
        self.assertEqual(invoke(current=dict(pyridine["a01"], coupling_potential_status="unresolved")), [])
        self.assertEqual(invoke(sequences_override=lambda seqs: seqs.update(a02=seqs["a02"][:-1] + "A")), [])
        self.assertEqual(invoke(genes_override=lambda gs: gs["a02"].update(contig="other_contig")), [])
        low_current = deepcopy(pyridine["a01"])
        low_current["confidence_hint"] = "low"
        low_current["selection_risks"] = ["reference_contact_difference:C299->H"]
        low_results, _ = prepare(annotation, anchors, identities,
                                 {"a01": low_current}, root=ROOT)
        low_prior = next(h for h in low_results["a01"]["source_like_component_hypotheses"]
                         if h["status"] == "source_like_module_component_potential")
        self.assertTrue(low_prior["prefer_source_function"])
        self.assertEqual(low_prior["confidence"], "low")
        self.assertIn("reference_contact_difference:C299->H", low_prior["selection_risks"])
        weak = deepcopy(pyridine["a01"])
        weak["evidence"]["roles"]["HydB"][0]["identity"] = "100"
        weak["evidence"]["roles"]["HydB"][0]["query_coverage"] = "1"
        weak_results, _ = prepare(annotation, anchors, identities, {"a01": weak}, root=ROOT)
        low = next(h for h in weak_results["a01"]["source_like_component_hypotheses"]
                   if h["status"] == "source_like_module_component_potential")
        self.assertFalse(low["prefer_source_function"])
        self.assertEqual(low["confidence"], "low")

    def _assay_bank_fixture(self, mutate=None):
        root = Path(self.tmp.name) / "assay_root"
        source = root / "references/raw/round118_Thermococcus_hyh/AF039208.gb"
        review = root / "experiments/round118_reference_function/Tlitoralis_HydBGDA_source_review.json"
        source.parent.mkdir(parents=True, exist_ok=True)
        review.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "references/raw/round118_Thermococcus_hyh/AF039208.gb", source)
        shutil.copy2(ROOT / "experiments/round118_reference_function/Tlitoralis_HydBGDA_source_review.json", review)
        data = json.loads((ROOT / "research_workflow/reference_assay_sets.json").read_text())
        # This fixture intentionally exercises only the historical Tlitoralis
        # set; the production registry may contain additional source-bound sets.
        data["sets"] = [s for s in data["sets"] if s["structure_id"] == "NiFe_Tlitoralis_HydBGDA"]
        if mutate:
            mutate(data)
        bank = root / "research_workflow/reference_assay_sets.json"
        bank.parent.mkdir(parents=True, exist_ok=True)
        bank.write_text(json.dumps(data))
        model = {"assay_bank": "research_workflow/reference_assay_sets.json", "sources": {
            "research_workflow/reference_assay_sets.json": _sha(bank),
            "references/raw/round118_Thermococcus_hyh/AF039208.gb": _sha(source),
            "experiments/round118_reference_function/Tlitoralis_HydBGDA_source_review.json": _sha(review),
        }}
        return root, model

    def test_assay_bank_rejects_unbound_or_corrupt_members_and_contracts(self):
        root, model = self._assay_bank_fixture()
        del model["sources"]["references/raw/round118_Thermococcus_hyh/AF039208.gb"]
        with self.assertRaises(ValueError):
            assay_bank(model, root)

        def bad_sequence(data):
            data["sets"][0]["members"][0]["sequence_sha256"] = "0" * 64
        root, model = self._assay_bank_fixture(bad_sequence)
        with self.assertRaises(ValueError):
            assay_bank(model, root)

        def duplicate_role(data):
            data["sets"][0]["members"][1]["role"] = data["sets"][0]["members"][0]["role"]
        root, model = self._assay_bank_fixture(duplicate_role)
        with self.assertRaises(ValueError):
            assay_bank(model, root)

        def foreign_contract(data):
            data["sets"][0]["contracts"][0]["structure_id"] = "unknown_structure"
        root, model = self._assay_bank_fixture(foreign_contract)
        with self.assertRaises(ValueError):
            assay_bank(model, root)

    def test_model_scope_rejects_empty_or_bad_source_binding(self):
        base = {
            "version": "reference-protein-set-assay-potential-v1",
            "admitted_systems": {"known": ["contract"]},
            "sources": {"source.json": "0" * 64},
        }
        empty = deepcopy(base)
        empty["admitted_systems"] = {}
        with self.assertRaises(ValueError):
            _validate_model_scope(empty)
        bad = deepcopy(base)
        bad["sources"]["source.json"] = "not-a-sha"
        with self.assertRaises(ValueError):
            _validate_model_scope(bad)

    def test_model_scope_rejects_unknown_structure_or_contract(self):
        model = {
            "admitted_systems": {"missing": ["contract"]},
        }
        with self.assertRaises(ValueError):
            _validate_admitted_bank(model, [], [])
        model = {"admitted_systems": {"known": ["missing"]}}
        systems = [{"structure_id": "known", "reference_contract_ids": ["missing"]}]
        with self.assertRaises(ValueError):
            _validate_admitted_bank(model, systems, [])

    def test_local_member_group_partition_is_strict(self):
        selected = [{"role": "HydA", "kind": "catalytic"},
                    {"role": "HydB", "kind": "module"},
                    {"role": "FdhA", "kind": "module"}]
        self.assertEqual(_validate_member_groups([["HydA", "HydB"], ["FdhA"]], selected),
                         [["HydA", "HydB"], ["FdhA"]])
        for groups in (None, 1, ["HydA"], [[]], [["HydA", "HydB"], ["HydB", "FdhA"]], [["HydA", "unknown"], ["FdhA"]], [["HydB"], ["FdhA"]]):
            if groups is None:
                continue
            with self.subTest(groups=groups):
                with self.assertRaises(ValueError):
                    _validate_member_groups(groups, selected)
        with self.assertRaises(ValueError):
            tuple_proof("anchor", {"evidence": {"roles": {}}}, selected,
                        {}, {}, "genome", {"local_max_span_bp": 5000}, member_groups=[None])

    def test_partial_reference_member_requires_bound_coding_proof(self):
        ref = self.refs_by_structure["NiFe_Amobile_HydABCSL"][0]
        pid = "partial_member"
        seq = ref["sequence"]
        genes = {pid: {"id": pid, "contig": ref["contig"], "start": ref["start"],
                       "end": ref["end"], "strand": ref["strand"],
                       "contig_length": ref["end"] + 100, "partial": "True"},
                 "anchor": {"id": "anchor", "contig": ref["contig"], "start": ref["start"],
                            "end": ref["end"], "strand": ref["strand"], "contig_length": ref["end"] + 100}}
        result = {"entry_route": "exact_full_protein_lookup", "evidence": {"roles": {ref["role"]: [{"protein_id": pid}]}}}
        errors, matched = tuple_proof("anchor", result, [dict(ref, kind="module")], genes,
                                      {pid: seq, "anchor": seq}, "genome",
                                      {"local_max_span_bp": 5000})
        self.assertIn("partial_or_unknown_member_" + ref["role"], errors)
        self.assertTrue(matched[0]["raw_partial"])
        self.assertFalse(matched[0]["coding_proof_used"])

    def test_partial_duplicate_is_competitor_only_with_bound_proof(self):
        gene = {"contig": "ctg", "start": 10, "end": 40, "strand": 1}
        seq = "M" + "A" * 9
        proof = {"source_artifact": {"path": "/source.fna", "sha256": "a" * 64}}
        with unittest.mock.patch("reference_reactions.coding_proof_supports", return_value=True):
            self.assertTrue(_coding_proof_valid("partial", gene, seq, {"partial": proof},
                                                {"/source.fna": "a" * 64}))
        self.assertFalse(_coding_proof_valid("partial", gene, seq, {},
                                             {"/source.fna": "a" * 64}))

    def test_amobile_exact_native_unknown_length_is_assay_potential(self):
        anchors, identities, pyridine, ids = self._fixture(
            "NiFe_Amobile_HydABCSL", input_mode="genome", unknown_length=True,
            anchor_status="candidate")
        results, bound = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        row = results[anchors[0]["protein_id"]]
        self.assertEqual(row["status"], "exact_reference_sequence_potential")
        self.assertEqual(row["query_activity"], "not_tested")
        self.assertEqual(row["in_vivo_direction"], "unknown")
        self.assertEqual({x["contract_id"] for x in row["capabilities"]}, {
            "Ref_Amobile_HydABCSL_oxidation", "Ref_Amobile_HydABCSL_evolution"})
        self.assertTrue(any(x.get("contig_length") == "not_supplied" for x in row["tuple_matches"]))
        self.assertIn("proteins.faa", " ".join(bound))

    def test_mthermo_exact_contract_is_typed_and_nadp_scoped(self):
        anchors, identities, pyridine, _ = self._fixture("FeFe_Mthermo_HydABCDEF")
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        row = results[anchors[0]["protein_id"]]
        self.assertEqual(row["status"], "exact_reference_sequence_potential")
        cap = row["capabilities"][0]
        self.assertEqual(cap["contract_id"], "Ref_Mthermo_HydABCDEF_NADP_oxidation")
        self.assertEqual(cap["reaction_key"], "H2_NADP_reduction_without_required_Fd")
        self.assertEqual(cap["scope"], "reference_assay_conditions_only")
        self.assertEqual(row["query_activity"], "not_tested")

    def test_kist_grouped_exact_tuple_accepts_remote_fdha_member(self):
        base = ROOT / "results/round145_KIST612_registered_old_consumer_v1/annotation"
        with (base / "catalytic_anchors.tsv").open() as handle:
            anchor = next(r for r in csv.DictReader(handle, delimiter="\t") if r["protein_id"] == "ELI_0847")
        anchor = dict(anchor, input_mode="genome")
        pyridine = json.loads((base.parent / "functions/pyridine_function_evidence.json").read_text())
        results, _ = prepare(base, [anchor], {"ELI_0847": "calibrated_high"}, pyridine, root=ROOT)
        row = results["ELI_0847"]
        self.assertEqual(row["status"], "exact_reference_sequence_potential")
        self.assertEqual(row["declared_member_groups"], [["HydA", "HydB", "HydC"], ["FdhA"]])
        self.assertEqual({m["role"] for m in row["tuple_matches"]}, {"HydA", "HydB", "HydC", "FdhA"})

    def test_group_scope_does_not_leak_between_kist_and_old_source(self):
        kist = ROOT / "results/round145_KIST612_registered_old_consumer_v1/annotation"
        old = ROOT / "results/round117_function_conversion_v4/Syntrophus_A_given_pair_v2/annotation"
        combined = Path(self.tmp.name) / "combined" / "annotation"
        shutil.copytree(kist, combined)
        for name in ("genes.tsv", "complex_members.tsv", "catalytic_anchors.tsv"):
            with (combined / name).open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            with (old / name).open() as handle:
                rows.extend(csv.DictReader(handle, delimiter="\t"))
            _tsv(combined / name, rows)
        with (combined / "proteins.faa").open("a") as out, (old / "proteins.faa").open() as source:
            out.write("\n" + source.read())
        manifest = json.loads((combined / "run_manifest.json").read_text())
        hashes = {name: _sha(combined / name) for name in ("proteins.faa", "genes.tsv", "complex_members.tsv", "catalytic_anchors.tsv")}
        manifest["output_sha256"].update(hashes)
        manifest["normalized_proteins_sha256"] = hashes["proteins.faa"]
        (combined / "run_manifest.json").write_text(json.dumps(manifest))
        with (combined / "catalytic_anchors.tsv").open() as handle:
            anchor_rows = list(csv.DictReader(handle, delimiter="\t"))
        k_anchor = dict(next(r for r in anchor_rows if r["protein_id"] == "ELI_0847"), input_mode="genome")
        old_anchor = dict(next(r for r in anchor_rows if r["protein_id"] == "q02"), input_mode="genome")
        k_pyr = json.loads((kist.parent / "functions/pyridine_function_evidence.json").read_text())
        old_pyr = json.loads((old.parent / "functions/pyridine_function_evidence.json").read_text())["q02"]
        results, _ = prepare(combined, [k_anchor, old_anchor],
                             {"ELI_0847": "calibrated_high", "q02": "sequence_supported_uncalibrated_noncanonical"},
                             {"ELI_0847": k_pyr["ELI_0847"], "q02": old_pyr}, root=ROOT)
        self.assertEqual(results["ELI_0847"]["status"], "exact_reference_sequence_potential")
        self.assertNotIn("declared_member_groups", results["q02"])

    def test_kist_duplicate_exact_fdhA_selects_nearest_and_retains_alternative(self):
        base = ROOT / "results/round145_KIST612_A_03_duplicateFdhA_after_v2/annotation"
        with (base / "catalytic_anchors.tsv").open() as handle:
            anchor = next(r for r in csv.DictReader(handle, delimiter="\t") if r["protein_id"] == "ELI_0847")
        anchor = dict(anchor, input_mode="assembly")
        pyridine = json.loads((base.parent / "functions/pyridine_function_evidence.json").read_text())
        results, _ = prepare(base, [anchor], {"ELI_0847": "calibrated_high"}, pyridine, root=ROOT)
        row = results["ELI_0847"]
        self.assertEqual(row["status"], "exact_reference_sequence_potential")
        ambiguity = row["member_selection_ambiguities"]
        self.assertEqual(ambiguity[0]["selected_protein_id"], "ELI_0994")
        self.assertEqual(ambiguity[0]["alternative_protein_ids"], ["ELI_0994_dup"])
        self.assertEqual(next(m["protein_id"] for m in row["tuple_matches"] if m["role"] == "FdhA"), "ELI_0994")

    def test_kist_coding_locus_core_proof_survives_group_candidate_validation(self):
        base = ROOT / "results/round145_KIST612_FNA_FdhA_same_bin_after_v3"
        annotation = base / "annotation"
        with (annotation / "catalytic_anchors.tsv").open() as handle:
            anchor = next(r for r in csv.DictReader(handle, delimiter="\t") if r["protein_id"] == "CP002273.1_759")
        anchor = dict(anchor, input_mode="bin")
        coding = json.loads((base / "functions/coding_completeness/coding_completeness.json").read_text())
        proofs = {row["protein_id"]: row for row in coding["rows"]}
        pyridine = json.loads((base / "functions/pyridine_function_evidence.json").read_text())
        results, _ = prepare(annotation, [anchor], {"CP002273.1_759": "calibrated_high"}, pyridine,
                             root=ROOT, coding_completeness=proofs)
        row = results["CP002273.1_759"]
        self.assertEqual(row["status"], "exact_reference_coding_locus_potential")
        self.assertTrue(next(m for m in row["tuple_matches"] if m["role"] == "HydA")["coding_locus_proof_used"])
        self.assertTrue(next(m for m in row["tuple_matches"] if m["role"] == "FdhA")["coding_locus_proof_used"])

    def test_component_hypothesis_keeps_source_tuple_when_only_locality_breaks(self):
        annotation = ROOT / "results/round145_KIST612_A_core_group_broken_after_v3/annotation"
        evidence = ROOT / "results/round146_KIST612_A_core_group_broken_after_v1/functions/pyridine_function_evidence.json"
        with (annotation / "catalytic_anchors.tsv").open() as handle:
            anchor = next(r for r in csv.DictReader(handle, delimiter="\t") if r["protein_id"] == "ELI_0847")
        claims = json.loads(evidence.read_text())
        results, _ = prepare(annotation, [dict(anchor, input_mode="genome")],
                             {"ELI_0847": "calibrated_high"}, claims, root=ROOT)
        row = results["ELI_0847"]
        self.assertEqual(row["status"], "not_supported")
        self.assertEqual(row["capabilities"], [])
        component = row["sequence_component_hypotheses"][0]
        self.assertEqual(component["status"], "source_equivalent_components_potential")
        self.assertEqual({m["role"] for m in component["matched_members"]}, {"HydA", "HydB", "HydC", "FdhA"})
        self.assertTrue(component["linkage_blockers"])
        self.assertEqual(component["hard_blockers"], [])
        self.assertTrue(all(c["scope"] == "sequence_component_function_hypothesis" for c in component["capabilities"]))

    def test_source_like_fdh_role_keeps_coarse_formate_potential_without_exact_tuple(self):
        base = ROOT / "results/round145_KIST612_A_02_FdhA_sequence_variant_after_v2/annotation"
        model = json.loads((ROOT / "research_workflow/reference_reactions_model.json").read_text())
        refs, systems, _ = assay_bank(model, ROOT)
        system = next(s for s in systems if s["structure_id"] == "FeFe_KIST612_FDH_HydABC")
        fdh = next(r for r in refs if r["structure_id"] == system["structure_id"] and r["role"] == "FdhA")
        with (base.parent / "domains" / "query_domains.tsv").open() as handle:
            domain_rows = list(csv.DictReader(handle, delimiter="\t"))
        wanted = {"PF13510", "PF04879", "PF00384", "PF01568"}
        source_profile_proof = [r for r in domain_rows if r["protein_id"] == "ELI_0994_variant" and r["accession"].split(".")[0] in wanted and r["passes"] == "True"]
        source_profile_proof.sort(key=lambda r: int(r["query_start"]))
        inventory = {
            "gene_roles": {"ELI_0994_variant": [{
                "role_class": "Fdh_molybdo", "preferred_role_class": "Fdh_molybdo", "status": "supported_role_homology",
                "fdh_chemistry": {"status": "formate_class_compatible", "confidence": "medium"},
                "source_member_compatibility": [{
                    "source_sequence_sha256": fdh["sequence_sha256"],
                    "status": "source_like_module_supported", "paired_identity": 99.9,
                    "reference_coverage": 1.0, "query_coverage": 1.0,
                    "source_profile_proof": source_profile_proof
                }]
            }]},
            "fdh_source_member_classes": {fdh["sequence_sha256"]: "Fdh_molybdo"}
        }
        with (base / "catalytic_anchors.tsv").open() as handle:
            anchor = next(r for r in csv.DictReader(handle, delimiter="\t") if r["protein_id"] == "ELI_0847")
        pyridine = json.loads((base.parent / "functions/pyridine_function_evidence.json").read_text())
        self.assertEqual(pyridine["ELI_0847"]["hydrogen_reaction_core_status"], "unresolved")
        model["version"] = "reference-assay-potential-v7"
        with unittest.mock.patch("reference_reactions.load_model", return_value=model):
            results, _ = prepare(base, [dict(anchor, input_mode="assembly")], {"ELI_0847": "calibrated_high"}, pyridine,
                                 root=ROOT, module_role_inventory=inventory)
        row = results["ELI_0847"]
        self.assertEqual(row["status"], "not_supported")
        self.assertEqual(row["capabilities"], [])
        component = row["source_like_component_hypotheses"][0]
        self.assertEqual(component["status"], "source_like_formate_component_potential")
        self.assertEqual(component["supported_functions"], ["candidate_H2_formate_coupling"])
        fdh_member = next(m for m in component["matched_members"] if m["role"] == "FdhA")
        self.assertEqual(fdh_member["matched_by"], "supported_FDH_role_homology")
        self.assertFalse(fdh_member["source_equivalent"])
        self.assertEqual(component["reference_reaction_details_transferred"], False)

    def test_component_witness_cannot_share_members_with_another_anchor(self):
        source = ROOT / "results/round145_KIST612_A_core_group_broken_after_v3/annotation"
        annotation = Path(self.tmp.name) / "shared_components" / "annotation"
        shutil.copytree(source, annotation)
        with (annotation / "genes.tsv").open() as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        genes.append(dict(next(row for row in genes if row["id"] == "ELI_0847"), id="ELI_0847_copy"))
        _tsv(annotation / "genes.tsv", genes)
        with (annotation / "proteins.faa").open() as handle:
            proteins = list(SeqIO.parse(handle, "fasta"))
        copy_record = deepcopy(next(record for record in proteins if record.id == "ELI_0847"))
        copy_record.id = "ELI_0847_copy"; copy_record.description = ""; proteins.append(copy_record)
        SeqIO.write(proteins, annotation / "proteins.faa", "fasta")
        with (annotation / "catalytic_anchors.tsv").open() as handle:
            anchors = list(csv.DictReader(handle, delimiter="\t"))
        anchors.append(dict(next(row for row in anchors if row["protein_id"] == "ELI_0847"), protein_id="ELI_0847_copy"))
        _tsv(annotation / "catalytic_anchors.tsv", anchors)
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        for name in ("genes.tsv", "proteins.faa", "catalytic_anchors.tsv"):
            manifest["output_sha256"][name] = _sha(annotation / name)
        manifest["normalized_proteins_sha256"] = _sha(annotation / "proteins.faa")
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        evidence = json.loads((ROOT / "results/round146_KIST612_A_core_group_broken_after_v1/functions/pyridine_function_evidence.json").read_text())
        evidence["ELI_0847_copy"] = deepcopy(evidence["ELI_0847"])
        evidence["ELI_0847_copy"]["evidence"]["anchor"] = "ELI_0847_copy"
        anchor_rows = [dict(row, input_mode="genome") for row in anchors]
        results, _ = prepare(annotation, anchor_rows,
                             {"ELI_0847": "calibrated_high", "ELI_0847_copy": "calibrated_high"},
                             evidence, root=ROOT)
        for qid in ("ELI_0847", "ELI_0847_copy"):
            component = results[qid]["sequence_component_hypotheses"][0]
            self.assertEqual(component["status"], "not_supported")
            self.assertIn("component_resources_shared_with_other_anchor", component["hard_blockers"])
            self.assertEqual(component["capabilities"], [])

    def test_equivalent_member_choice_is_stable_and_keeps_tuple_limits(self):
        base = ROOT / "results/round145_KIST612_A_03_duplicateFdhA_after_v2/annotation"
        originals = {name:(base/name).read_bytes() for name in
                     ("genes.tsv","proteins.faa","catalytic_anchors.tsv","complex_members.tsv","reaction_capabilities.tsv","run_manifest.json")}
        pyridine = json.loads((base.parent/"functions/pyridine_function_evidence.json").read_text())
        model = json.loads((ROOT/"research_workflow/reference_reactions_model.json").read_text())
        bank = assay_bank(model, ROOT)
        for change in ("tie_reversed_order", "equally_distant_loci", "invalid_copy", "second_role_ambiguous", "resource_conflict", "no_groups"):
            with self.subTest(change=change):
                for name,data in originals.items():(self.annotation/name).write_bytes(data)
                with (self.annotation/"genes.tsv").open() as h:genes = list(csv.DictReader(h, delimiter="\t"))
                with (self.annotation/"proteins.faa").open() as h:proteins = list(SeqIO.parse(h, "fasta"))
                duplicate = next(g for g in genes if g["id"]=="ELI_0994_dup")
                if change=="tie_reversed_order":
                    original = next(g for g in genes if g["id"]=="ELI_0994")
                    duplicate.update({k:original[k] for k in ("contig","strand","start","end")})
                    genes.reverse();proteins.reverse()
                elif change=="equally_distant_loci":
                    original = next(g for g in genes if g["id"]=="ELI_0994")
                    anchor_gene = next(g for g in genes if g["id"]=="ELI_0847")
                    start = 2*int(anchor_gene["start"])-int(original["start"])
                    duplicate.update(start=str(start), end=str(start+int(original["end"])-int(original["start"])))
                elif change=="invalid_copy":duplicate["end"] = str(int(duplicate["end"])+1)
                elif change=="second_role_ambiguous":
                    original = next(g for g in genes if g["id"]=="ELI_0846")
                    genes.append(dict(original, id="HydB_copy"))
                    record = deepcopy(next(p for p in proteins if p.id=="ELI_0846"))
                    record.id="HydB_copy";record.description="";proteins.append(record)
                _tsv(self.annotation/"genes.tsv",genes)
                SeqIO.write(proteins,self.annotation/"proteins.faa","fasta")
                manifest=json.loads(originals["run_manifest.json"])
                for name in ("genes.tsv","proteins.faa"):manifest["output_sha256"][name]=_sha(self.annotation/name)
                manifest["normalized_proteins_sha256"]=_sha(self.annotation/"proteins.faa")
                (self.annotation/"run_manifest.json").write_text(json.dumps(manifest))
                with (self.annotation/"catalytic_anchors.tsv").open() as h:
                    anchor=next(r for r in csv.DictReader(h,delimiter="\t") if r["protein_id"]=="ELI_0847")
                anchor=dict(anchor,input_mode="assembly")
                claims=deepcopy(pyridine);identities={"ELI_0847":"calibrated_high"}
                if change=="resource_conflict":
                    claims["ELI_0847"]={"structure_id":"FeFe_KIST612_FDH_HydABC", "evidence":{"cross_anchor_conflicts":[{"rival_anchor":"other_core"}]}}
                    identities["other_core"]="calibrated_high"
                refs,systems,contracts=deepcopy(bank)
                if change=="no_groups":
                    next(s for s in systems if s["structure_id"]=="FeFe_KIST612_FDH_HydABC").pop("local_member_groups")
                with unittest.mock.patch("reference_reactions.assay_bank",return_value=(refs,systems,contracts)):
                    results,_=prepare(self.annotation,[anchor],identities,claims,root=ROOT)
                row=results["ELI_0847"]
                if change in {"equally_distant_loci","second_role_ambiguous","resource_conflict","no_groups"}:
                    self.assertEqual(row["capabilities"],[])
                    if change=="equally_distant_loci":
                        component=row["sequence_component_hypotheses"][0]
                        self.assertEqual(component["status"],"source_equivalent_components_potential")
                        ambiguity=next(a for a in component["member_selection_ambiguities"] if a["role"]=="FdhA")
                        self.assertEqual(ambiguity["equally_supported_loci"],["ELI_0994_dup"])
                    if change=="second_role_ambiguous":
                        component=row["sequence_component_hypotheses"][0]
                        self.assertEqual(component["status"],"source_equivalent_components_potential")
                        self.assertTrue(any(a["role"]=="HydB" and "HydB_copy" in a["alternative_protein_ids"]
                                            for a in component["member_selection_ambiguities"]))
                else:
                    self.assertEqual(len(row["capabilities"]),4)
                    self.assertEqual(next(m["protein_id"] for m in row["tuple_matches"] if m["role"]=="FdhA"),"ELI_0994")
                    if change=="invalid_copy":
                        self.assertFalse(row.get("member_selection_ambiguities"))
                        self.assertTrue(any(h.get("member_selection_rejections") for h in row["tuple_hypotheses"]))
                    else:self.assertEqual(row["member_selection_ambiguities"][0]["alternative_protein_ids"],["ELI_0994_dup"])

    def test_tkodakarensis_assay_bridge_has_three_typed_capabilities(self):
        model = json.loads((ROOT / "research_workflow/reference_reactions_model.json").read_text())
        refs, _systems, _contracts = assay_bank(model, ROOT)
        refs = [r for r in refs if r["structure_id"] == "NiFe_Tkodakarensis_HyhBGSL"]
        ids = {r["role"]: "tk_" + r["role"] for r in refs}
        genes = []
        for index, ref in enumerate(refs):
            start = 100 + index * 1500
            genes.append(dict(id=ids[ref["role"]], contig="tk", start=str(start), end=str(start + 3 * len(ref["sequence"]) - 1),
                              strand="1", partial="False", contig_length="100000", bin_id="b1", coordinate_evidence="source_CDS"))
        _tsv(self.annotation / "genes.tsv", genes)
        _tsv(self.annotation / "complex_members.tsv", [])
        _tsv(self.annotation / "catalytic_anchors.tsv", [dict(protein_id=ids["HyhL"], structure_id="NiFe_Tkodakarensis_HyhBGSL", input_mode="genome", potential_function_ids="")])
        with (self.annotation / "proteins.faa").open("w") as handle:
            for ref in refs:
                handle.write(f">{ids[ref['role']]}\n{ref['sequence']}\n")
        files = {name: _sha(self.annotation / name) for name in ("proteins.faa", "genes.tsv", "complex_members.tsv", "catalytic_anchors.tsv")}
        run = dict(status="complete", input_mode="genome", thresholds=json.loads((ROOT / "config/thresholds.json").read_text()),
                   output_sha256=files, normalized_proteins_sha256=files["proteins.faa"])
        (self.annotation / "run_manifest.json").write_text(json.dumps(run))
        anchors = [dict(protein_id=ids["HyhL"], structure_id="NiFe_Tkodakarensis_HyhBGSL", input_mode="genome")]
        results, _ = prepare(self.annotation, anchors, {ids["HyhL"]: "sequence_supported_uncalibrated_core"}, {}, root=ROOT)
        self.assertEqual({c["contract_id"] for c in results[ids["HyhL"]]["capabilities"]}, {
            "Ref_NiFe_Tkodakarensis_HyhBGSL_NADP_oxidation",
            "Ref_NiFe_Tkodakarensis_HyhBGSL_NADPH_evolution",
            "Ref_NiFe_Tkodakarensis_HyhBGSL_NAD_oxidation"})
        self.assertTrue(all(c["scope"] == "reference_assay_conditions_only" for c in results[ids["HyhL"]]["capabilities"]))
        protein_path = self.annotation / "proteins.faa"
        text = protein_path.read_text()
        protein_path.write_text(text.replace(refs[-1]["sequence"], refs[-1]["sequence"][:-1] + ("A" if refs[-1]["sequence"][-1] != "A" else "G"), 1))
        files["proteins.faa"] = _sha(protein_path)
        run["output_sha256"] = files; run["normalized_proteins_sha256"] = files["proteins.faa"]
        (self.annotation / "run_manifest.json").write_text(json.dumps(run))
        mutated, _ = prepare(self.annotation, anchors, {ids["HyhL"]: "sequence_supported_uncalibrated_core"}, {}, root=ROOT)
        self.assertEqual(mutated[ids["HyhL"]]["status"], "not_supported")
        self.assertEqual(mutated[ids["HyhL"]]["capabilities"], [])

    def test_one_aa_variant_does_not_promote_or_change_identity(self):
        anchors, identities, pyridine, _ = self._fixture(
            "NiFe_Amobile_HydABCSL", mutate="HydB")
        before = deepcopy(identities)
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results[anchors[0]["protein_id"]]["status"], "not_supported")
        self.assertEqual(identities, before)

    def test_partial_duplicate_bad_length_and_cross_bin_are_rejected(self):
        cases = (
            {"partial": ("HydB",)}, {"duplicate_role": "HydB"},
            {"bad_length": True}, {"cross_bin": ("HydB",)},
        )
        for changes in cases:
            with self.subTest(changes=changes):
                if "cross_bin" in changes:
                    changes["input_mode"] = "bin"
                anchors, identities, pyridine, _ = self._fixture("NiFe_Amobile_HydABCSL", **changes)
                results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
                self.assertEqual(results[anchors[0]["protein_id"]]["status"], "not_supported")

    def test_unsupported_identity_rejects_exact_tuple(self):
        anchors, identities, pyridine, _ = self._fixture(
            "NiFe_Amobile_HydABCSL", identity_status="evidence_insufficient")
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results[anchors[0]["protein_id"]]["status"], "not_supported")

    def test_inherited_mthermo_assigned_row_is_typed(self):
        contract = self.contract_by_id["Ref_Mthermo_HydABCDEF_NADP_oxidation"]
        anchors, identities, _, _ = self._fixture(
            "FeFe_Mthermo_HydABCDEF", potential_ids=[contract["contract_id"]])
        pyridine = {}
        qid = anchors[0]["protein_id"]
        _tsv(self.annotation / "reaction_capabilities.tsv", [{
            "complex_anchor": qid, "function_id": contract["contract_id"],
            "structure_id": contract["structure_id"], "status": "reference_sequence_supported_potential",
            "source": contract["primary_source"], "molecular_reaction": contract["molecular_reaction"],
            "reference_reaction": contract["reaction_key"], "reference_direction": contract["hydrogen_role"],
            "query_net_direction": "unresolved",
            "capability_scope": "exact_reference_sequence_molecular_reaction_potential",
            "reference_evidence_scope": contract["reference_evidence_scope"],
        }])
        self._manifest(input_mode="genome")
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        row = results[qid]
        self.assertEqual(row["status"], "exact_reference_sequence_potential")
        self.assertEqual(row["inherited_assignments"][0]["contract_id"], contract["contract_id"])
        self.assertEqual(row["capabilities"][0]["electron_acceptors"], ["NADP+"])
        self.assertEqual(row["capabilities"][0]["primary_source"], contract["primary_source"])

    def test_exact_native_tuple_does_not_require_pyridine_gate(self):
        anchors, identities, _, _ = self._fixture("NiFe_Amobile_HydABCSL")
        results, _ = prepare(self.annotation, anchors, identities, {}, root=ROOT)
        row = results[anchors[0]["protein_id"]]
        self.assertEqual(row["status"], "exact_reference_sequence_potential")
        self.assertEqual({c["contract_id"] for c in row["capabilities"]}, {
            "Ref_Amobile_HydABCSL_oxidation", "Ref_Amobile_HydABCSL_evolution"})

    def test_duplicate_exact_role_is_rejected_even_if_generic_evidence_picks_one(self):
        anchors, identities, pyridine, _ = self._fixture(
            "NiFe_Amobile_HydABCSL", duplicate_role="HydB")
        # Simulate the old generic selector retaining only one of two equally
        # exact local proteins. The independent exact lookup must still fail
        # closed on the sequence-level ambiguity.
        qid = anchors[0]["protein_id"]
        pyridine[qid]["evidence"]["roles"]["HydB"] = pyridine[qid]["evidence"]["roles"]["HydB"][:1]
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        row = results[qid]
        self.assertEqual(row["status"], "not_supported")
        self.assertTrue(any("ambiguous_exact_reference_role_HydB" in b for b in row["blockers"]))

    def test_tampered_or_unknown_assigned_row_fails_closed(self):
        anchors, identities, pyridine, _ = self._fixture("FeFe_Mthermo_HydABCDEF")
        qid = anchors[0]["protein_id"]
        _tsv(self.annotation / "reaction_capabilities.tsv", [{
            "complex_anchor": qid, "function_id": "Ref_Mthermo_HydABCDEF_NADP_oxidation",
            "status": "reference_sequence_supported_potential",
        }])
        self._manifest(tamper=True)
        with self.assertRaises(ValueError):
            prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        self._manifest(input_mode="genome")
        _tsv(self.annotation / "reaction_capabilities.tsv", [{
            "complex_anchor": qid, "function_id": "unknown_contract", "status": "reference_sequence_supported_potential",
        }])
        self._manifest()
        with self.assertRaises(ValueError):
            prepare(self.annotation, anchors, identities, pyridine, root=ROOT)

    def test_candidate_contract_table_is_not_positive_evidence(self):
        # Deliberately break one member: a row in the candidate-only table must
        # not turn an unsupported tuple into a positive assigned capability.
        anchors, identities, pyridine, _ = self._fixture("FeFe_Mthermo_HydABCDEF", mutate="HydB")
        _tsv(self.annotation / "reaction_contract_candidates.tsv", [{
            "complex_anchor": anchors[0]["protein_id"],
            "function_id": "Ref_Mthermo_HydABCDEF_NADP_oxidation",
            "status": "candidate",
        }])
        self._manifest(input_mode="genome")
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results[anchors[0]["protein_id"]]["status"], "not_supported")
        self.assertEqual(results[anchors[0]["protein_id"]]["capabilities"], [])

    def test_hnd_is_not_admitted_as_new_precise_reference(self):
        anchors, identities, pyridine, _ = self._fixture("FeFe_Sfructo_HndABCD")
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results[anchors[0]["protein_id"]]["status"], "not_supported")

    def test_source_gene_labels_do_not_affect_exact_sequence_tuple(self):
        labels = {r["role"]: "renamed_product_label" for r in self.refs_by_structure["NiFe_Amobile_HydABCSL"]}
        anchors, identities, pyridine, _ = self._fixture("NiFe_Amobile_HydABCSL", labels=labels)
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results[anchors[0]["protein_id"]]["status"], "exact_reference_sequence_potential")

    def test_original_anchor_table_can_contain_unselected_anchor_subset(self):
        anchors, identities, pyridine, _ = self._fixture("NiFe_Amobile_HydABCSL")
        original = []
        with (self.annotation / "catalytic_anchors.tsv").open() as fh:
            original = list(csv.DictReader(fh, delimiter="\t"))
        original.append(dict(original[0], protein_id="outside_anchor", structure_id="FeFe_Sfructo_HndABCD"))
        _tsv(self.annotation / "catalytic_anchors.tsv", original)
        self._manifest(input_mode="genome")
        results, _ = prepare(self.annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(set(results), {anchors[0]["protein_id"]})

    def _copy_syntrophus_annotation(self, case):
        """Copy the frozen public annotation and rebind its mutable manifest."""
        source = ROOT / "results/round117_function_conversion_v4" / case / "annotation"
        target = Path(self.tmp.name) / f"{case}_{len(list(Path(self.tmp.name).iterdir()))}" / "annotation"
        shutil.copytree(source, target)
        manifest_path = target / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        hashes = dict(manifest.get("output_sha256", {}))
        for name in ("proteins.faa", "genes.tsv", "complex_members.tsv", "catalytic_anchors.tsv"):
            hashes[name] = _sha(target / name)
        manifest["output_sha256"] = hashes
        manifest["normalized_proteins_sha256"] = _sha(target / "proteins.faa")
        manifest_path.write_text(json.dumps(manifest))
        return target

    def _copy_pfur_annotation(self, case):
        """Copy the round-134 source annotation named by its frozen metadata."""
        meta = json.loads((ROOT / "results/round134_Pfur_source_only_v1" /
                           case / "metadata.json").read_text())
        source = Path(meta["source_annotation"])
        target = Path(self.tmp.name) / f"{case}_real_{len(list(Path(self.tmp.name).iterdir()))}" / "annotation"
        shutil.copytree(source, target)
        return target

    def _pfur_inputs(self, annotation):
        with (annotation / "catalytic_anchors.tsv").open() as handle:
            anchor = next(csv.DictReader(handle, delimiter="\t"))
        qid = anchor["protein_id"]
        run = json.loads((annotation / "run_manifest.json").read_text())
        anchor = dict(anchor, input_mode=run["input_mode"])
        with (annotation / "complex_members.tsv").open() as handle:
            members = [row for row in csv.DictReader(handle, delimiter="\t") if row["anchor"] == qid]
        roles = {}
        for row in members:
            roles.setdefault(row["role"], []).append({"protein_id": row["protein_id"]})
        pyridine = {qid: {"structure_id": anchor["structure_id"],
                          "coupling_potential_status": "local_coupling_potential_uncalibrated",
                          "evidence": {"roles": roles}}}
        ids = sorted({row["protein_id"] for row in members})
        output = Path(self.tmp.name) / f"pfur_coding_{qid}_{len(list(Path(self.tmp.name).iterdir()))}"
        coding_prepare(annotation, output, candidate_ids=ids, root=ROOT)
        _manifest, proofs = coding_load_evidence(output, annotation=annotation, root=ROOT)
        return anchor, qid, pyridine, proofs

    def _saciditrophicus_probe(self, annotation, anchor_id):
        with (annotation / "catalytic_anchors.tsv").open() as handle:
            anchor = next(row for row in csv.DictReader(handle, delimiter="\t")
                          if row["protein_id"] == anchor_id)
        mode = json.loads((annotation / "run_manifest.json").read_text())["input_mode"]
        anchor = dict(anchor, input_mode=mode)
        return [anchor], {anchor_id: "sequence_supported_uncalibrated_noncanonical"}, {}

    def _boundary_proof(self, annotation, protein_id):
        output = Path(self.tmp.name) / f"coding_{protein_id}_{len(list(Path(self.tmp.name).iterdir()))}"
        coding_prepare(annotation, output, candidate_ids=[protein_id], root=ROOT)
        _manifest, rows = coding_load_evidence(output, annotation=annotation, root=ROOT)
        return rows[protein_id]

    def _rebind_boundary_source(self, annotation, coding_codon=None):
        manifest_path = annotation / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        item = next(item for item in manifest["input_files"] if item["path"].endswith(".fna"))
        source = Path(item["path"])
        rebound = Path(self.tmp.name) / "rebound_source.fna"
        shutil.copy2(source, rebound)
        if coding_codon is not None:
            with rebound.open() as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            record = next(record for record in records if record.id == "CP000252.1")
            # The first four encoded residues precede the called negative-strand
            # query.  Change coding codon 2 (GAG) through its reverse complement.
            self.assertEqual(str(record.seq)[1122813:1122816], "CTC")
            record.seq = Seq(str(record.seq)[:1122813] + coding_codon + str(record.seq)[1122816:])
            with rebound.open("w") as handle:
                SeqIO.write(records, handle, "fasta")
        item["path"] = str(rebound.resolve())
        item["sha256"] = _sha(rebound)
        manifest_path.write_text(json.dumps(manifest))
        return rebound

    def test_saciditrophicus_actual_A_gets_both_typed_reference_potentials(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_A_given_pair_v2")
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, "q02")
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        row = results["q02"]
        self.assertEqual(row["status"], "exact_reference_sequence_potential")
        self.assertEqual({c["contract_id"] for c in row["capabilities"]}, {
            "Ref_FeFe_Saciditrophicus_HydAB_NAD_oxidation",
            "Ref_FeFe_Saciditrophicus_HydAB_NADH_evolution",
        })
        self.assertEqual(identities["q02"], "sequence_supported_uncalibrated_noncanonical")

    def test_syntrophus_B_bound_coding_locus_gets_two_caps_without_identity_mutation(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
        qid = "CP000252.1_979"
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, qid)
        proof = self._boundary_proof(annotation, qid)
        before = deepcopy(identities)
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT,
                             coding_completeness={qid: proof})
        row = results[qid]
        self.assertEqual(row["status"], "exact_reference_coding_locus_potential")
        self.assertEqual({c["contract_id"] for c in row["capabilities"]}, {
            "Ref_FeFe_Saciditrophicus_HydAB_NAD_oxidation",
            "Ref_FeFe_Saciditrophicus_HydAB_NADH_evolution"})
        self.assertEqual(identities, before)
        locus = next(m for m in row["tuple_matches"] if m["protein_id"] == qid)
        self.assertTrue(locus.get("query_translation_boundary_different"))
        self.assertTrue(locus.get("coding_locus_proof_used"))
        self.assertEqual(locus.get("coding_proof_digest"), proof["record_digest"])

    def test_syntrophus_B_synonymous_upstream_boundary_is_accepted(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
        self._rebind_boundary_source(annotation, coding_codon="TTC")
        qid = "CP000252.1_979"
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, qid)
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT,
                             coding_completeness={qid: self._boundary_proof(annotation, qid)})
        self.assertEqual(results[qid]["status"], "exact_reference_coding_locus_potential")

    def test_syntrophus_B_upstream_aa_change_or_stop_is_rejected(self):
        for codon in ("GCC", "TTA"):
            with self.subTest(codon=codon):
                annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
                self._rebind_boundary_source(annotation, coding_codon=codon)
                qid = "CP000252.1_979"
                anchors, identities, pyridine = self._saciditrophicus_probe(annotation, qid)
                results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT,
                                     coding_completeness={qid: self._boundary_proof(annotation, qid)})
                self.assertNotEqual(results.get(qid, {}).get("status"), "exact_reference_coding_locus_potential")
                self.assertEqual(results.get(qid, {}).get("capabilities", []), [])

    def test_syntrophus_B_missing_source_or_wrong_coordinates_reject_locus(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
        qid = "CP000252.1_979"
        proof = self._boundary_proof(annotation, qid)
        proof["source_artifact"]["path"] = str(Path(self.tmp.name) / "missing.fna")
        proof["record_digest"] = _record_digest(proof)
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, qid)
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT,
                             coding_completeness={qid: proof})
        self.assertNotEqual(results.get(qid, {}).get("status"), "exact_reference_coding_locus_potential")
        annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
        genes_path = annotation / "genes.tsv"
        stale_proof = self._boundary_proof(annotation, qid)
        with genes_path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        next(row for row in rows if row["id"] == qid)["start"] = str(int(next(row for row in rows if row["id"] == qid)["start"]) + 3)
        _tsv(genes_path, rows)
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        manifest["output_sha256"]["genes.tsv"] = _sha(genes_path)
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, qid)
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT,
                             coding_completeness={qid: stale_proof})
        self.assertNotEqual(results.get(qid, {}).get("status"), "exact_reference_coding_locus_potential")

    def test_syntrophus_B_forged_locus_payload_without_local_proof_is_rejected(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
        qid = "CP000252.1_979"
        anchors, identities, _ = self._saciditrophicus_probe(annotation, qid)
        pyridine = {qid: {"structure_id": "FeFe_Saciditrophicus_HydAB",
            "entry_route": "exact_coding_locus_lookup",
            "coupling_potential_status": "local_coupling_potential_uncalibrated",
            "coding_locus": {"query_protein_sha256": "0" * 64},
            "evidence": {"roles": {"HydA": [{"protein_id": qid}], "HydB": [{"protein_id": "CP000252.1_980"}]}}}}
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        self.assertNotEqual(results.get(qid, {}).get("status"), "exact_reference_coding_locus_potential")
        self.assertEqual(results.get(qid, {}).get("capabilities", []), [])

    def test_saciditrophicus_one_aa_beta_variant_has_no_exact_caps(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_A_given_pair_v2")
        protein_path = annotation / "proteins.faa"
        with protein_path.open() as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        q01 = next(record for record in records if record.id == "q01")
        q01.seq = Seq(("A" if str(q01.seq)[0] != "A" else "G") + str(q01.seq)[1:])
        with protein_path.open("w") as handle:
            SeqIO.write(records, handle, "fasta")
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        manifest["output_sha256"]["proteins.faa"] = _sha(protein_path)
        manifest["normalized_proteins_sha256"] = _sha(protein_path)
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, "q02")
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results["q02"]["status"], "not_supported")
        self.assertEqual(results["q02"]["capabilities"], [])

    def test_saciditrophicus_beta_removal_fails_closed(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_A_given_pair_v2")
        protein_path = annotation / "proteins.faa"
        with protein_path.open() as handle:
            records = [record for record in SeqIO.parse(handle, "fasta") if record.id != "q01"]
        with protein_path.open("w") as handle:
            SeqIO.write(records, handle, "fasta")
        genes_path = annotation / "genes.tsv"
        with genes_path.open() as handle:
            genes = [row for row in csv.DictReader(handle, delimiter="\t") if row["id"] != "q01"]
        _tsv(genes_path, genes)
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        manifest["output_sha256"]["proteins.faa"] = _sha(protein_path)
        manifest["output_sha256"]["genes.tsv"] = _sha(genes_path)
        manifest["normalized_proteins_sha256"] = _sha(protein_path)
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, "q02")
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results["q02"]["status"], "not_supported")
        self.assertEqual(results["q02"]["capabilities"], [])

    def test_saciditrophicus_exact_beta_relocated_beyond_local_range_is_rejected(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_A_given_pair_v2")
        path = annotation / "genes.tsv"
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            if row["id"] == "q01":
                length = int(row["end"]) - int(row["start"])
                row["start"], row["end"] = "2800000", str(2800000 + length)
        _tsv(path, rows)
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        manifest["output_sha256"]["genes.tsv"] = _sha(path)
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, "q02")
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        self.assertEqual(results["q02"]["status"], "not_supported")
        self.assertEqual(results["q02"]["capabilities"], [])

    def test_syntrophus_B_601aa_core_cannot_inherit_605aa_reference_caps(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
        qid = "CP000252.1_979"
        with (annotation / "proteins.faa").open() as handle:
            proteins = {record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")}
        self.assertIn(qid, proteins)
        self.assertEqual(len(proteins[qid]), 601)
        anchors, identities, pyridine = self._saciditrophicus_probe(annotation, qid)
        results, _ = prepare(annotation, anchors, identities, pyridine, root=ROOT)
        self.assertNotIn("FeFe_Saciditrophicus_HydAB", {
            row.get("matched_structure_id") for row in results.values()})
        self.assertNotIn(qid, results)

    def test_coding_locus_candidate_ids_is_genome_only_and_boundary_shaped(self):
        annotation = self._copy_syntrophus_annotation("Syntrophus_B_whole_genome_v2")
        self.assertIn("CP000252.1_979", coding_locus_candidate_ids(annotation, root=ROOT))
        assembly = self._copy_syntrophus_annotation("Syntrophus_A_given_pair_v2")
        self.assertEqual(coding_locus_candidate_ids(assembly, root=ROOT), [])

    def test_pfur_real_source_annotations_cover_exact_and_two_coding_locus_cases(self):
        expected = {
            "A_PfurHII_given_members_v2": "exact_reference_sequence_potential",
            "B_PfurHII_native_locus_v2": "exact_reference_coding_locus_potential",
            "B_PfurHII_whole_genome_v2": "exact_reference_coding_locus_potential",
        }
        for case, status in expected.items():
            with self.subTest(case=case):
                annotation = self._copy_pfur_annotation(case)
                anchor, qid, pyridine, proofs = self._pfur_inputs(annotation)
                results, _ = prepare(annotation, [anchor],
                                     {qid: "sequence_supported_uncalibrated_noncanonical"},
                                     pyridine, root=ROOT, coding_completeness=proofs)
                self.assertEqual(results[qid]["status"], status)
                self.assertEqual({c["contract_id"] for c in results[qid]["capabilities"]}, {
                    "Ref_NiFe_Pfuriosus_HII_BGDA_NADP_oxidation",
                    "Ref_NiFe_Pfuriosus_HII_BGDA_NADPH_evolution",
                    "Ref_NiFe_Pfuriosus_HII_BGDA_NAD_oxidation",
                    "Ref_NiFe_Pfuriosus_HII_BGDA_NADH_evolution",
                })

    def test_pfur_real_coding_proof_rejects_member_aa_change_and_source_bound_change(self):
        annotation = self._copy_pfur_annotation("B_PfurHII_native_locus_v2")
        anchor, qid, pyridine, proofs = self._pfur_inputs(annotation)
        with (annotation / "proteins.faa").open() as handle:
            proteins = list(SeqIO.parse(handle, "fasta"))
        with (annotation / "complex_members.tsv").open() as handle:
            member = next(row["protein_id"] for row in csv.DictReader(handle, delimiter="\t")
                          if row["anchor"] == qid and row["role"] == "SHI_delta")
        record = next(r for r in proteins if r.id == member)
        record.seq = Seq(str(record.seq)[:80] + ("A" if str(record.seq)[80] != "A" else "G") + str(record.seq)[81:])
        with (annotation / "proteins.faa").open("w") as handle:
            SeqIO.write(proteins, handle, "fasta")
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        manifest["output_sha256"]["proteins.faa"] = _sha(annotation / "proteins.faa")
        manifest["normalized_proteins_sha256"] = _sha(annotation / "proteins.faa")
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        results, _ = prepare(annotation, [anchor],
                             {qid: "sequence_supported_uncalibrated_noncanonical"},
                             pyridine, root=ROOT, coding_completeness=proofs)
        self.assertNotEqual(results[qid]["status"], "exact_reference_coding_locus_potential")
        self.assertEqual(results[qid]["capabilities"], [])

    def test_pfur_native_locus_without_dna_proof_and_stale_coordinates_rejects(self):
        annotation = self._copy_pfur_annotation("B_PfurHII_native_locus_v2")
        anchor, qid, pyridine, proofs = self._pfur_inputs(annotation)
        no_proof, _ = prepare(annotation, [anchor],
                              {qid: "sequence_supported_uncalibrated_noncanonical"},
                              pyridine, root=ROOT)
        self.assertNotEqual(no_proof[qid]["status"], "exact_reference_coding_locus_potential")
        with (annotation / "genes.tsv").open() as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        next(row for row in genes if row["id"] == "PfurHII_native_locus_diagnostic_1")["start"] = "4"
        _tsv(annotation / "genes.tsv", genes)
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        manifest["output_sha256"]["genes.tsv"] = _sha(annotation / "genes.tsv")
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        stale, _ = prepare(annotation, [anchor],
                           {qid: "sequence_supported_uncalibrated_noncanonical"},
                           pyridine, root=ROOT, coding_completeness=proofs)
        self.assertNotEqual(stale[qid]["status"], "exact_reference_coding_locus_potential")
        self.assertEqual(stale[qid]["capabilities"], [])

    def test_pfur_longer_query_with_internal_alternate_start_keeps_four_caps(self):
        """An upstream ATG may precede the source TTG without changing its CDS."""
        annotation = self._copy_pfur_annotation("B_PfurHII_whole_genome_v2")
        anchor, qid, pyridine, _ = self._pfur_inputs(annotation)
        member = next(row["protein_id"] for row in csv.DictReader(
            (annotation / "complex_members.tsv").open(), delimiter="\t")
                      if row["anchor"] == qid and row["role"] == "SHI_beta")
        proteins = list(SeqIO.parse((annotation / "proteins.faa").open(), "fasta"))
        original = next(r for r in proteins if r.id == member)
        reference = str(original.seq)
        self.assertTrue(reference.startswith("MRY"))
        original.seq = Seq("ML" + reference[1:])
        with (annotation / "proteins.faa").open("w") as handle:
            SeqIO.write(proteins, handle, "fasta")
        genes = list(csv.DictReader((annotation / "genes.tsv").open(), delimiter="\t"))
        gene = next(row for row in genes if row["id"] == member)
        self.assertEqual(gene["start"], "1250021")
        gene["start"] = "1250018"
        _tsv(annotation / "genes.tsv", genes)
        manifest = json.loads((annotation / "run_manifest.json").read_text())
        source_item = next(item for item in manifest["input_files"] if item["path"].endswith(".fna"))
        source_copy = Path(self.tmp.name) / "pfur_longer_query_source.fna"
        with Path(source_item["path"]).open() as handle:
            source_records = list(SeqIO.parse(handle, "fasta"))
        source_record = next(record for record in source_records if record.id == "NC_003413.1")
        self.assertNotEqual(str(source_record.seq)[1250017:1250020], "ATG")
        source_record.seq = Seq(str(source_record.seq)[:1250017] + "ATG" + str(source_record.seq)[1250020:])
        with source_copy.open("w") as handle:
            SeqIO.write(source_records, handle, "fasta")
        source_item["path"] = str(source_copy.resolve())
        source_item["sha256"] = _sha(source_copy)
        manifest["output_sha256"]["proteins.faa"] = _sha(annotation / "proteins.faa")
        manifest["output_sha256"]["genes.tsv"] = _sha(annotation / "genes.tsv")
        manifest["normalized_proteins_sha256"] = _sha(annotation / "proteins.faa")
        (annotation / "run_manifest.json").write_text(json.dumps(manifest))
        output = Path(self.tmp.name) / "pfur_longer_query_coding"
        coding_prepare(annotation, output,
                       candidate_ids=[member, qid, "NC_003413.1_1364", "NC_003413.1_1365"],
                       root=ROOT)
        _manifest, proofs = coding_load_evidence(output, annotation=annotation, root=ROOT)
        self.assertEqual(proofs[member]["status"], "complete_under_bound_translation_model")
        self.assertEqual(proofs[member]["translation_match_mode"], "normalized_initiator_M")
        results, _ = prepare(annotation, [anchor],
                             {qid: "sequence_supported_uncalibrated_noncanonical"},
                             pyridine, root=ROOT, coding_completeness=proofs)
        row = results[qid]
        self.assertEqual(row["status"], "exact_reference_coding_locus_potential")
        self.assertEqual(len(row["capabilities"]), 4)
        beta = next(m for m in row["tuple_matches"] if m["role"] == "HydB")
        self.assertEqual(beta["missing_n_terminal_residues"], 0)
        self.assertEqual(beta["extra_n_terminal_residues"], 1)
        self.assertEqual(beta["boundary_translation_mode"], "query_has_extra_n_terminal_residues")
        self.assertEqual(beta["encoded_reference_start"], 1250021)
        self.assertEqual(beta["encoded_reference_end"], 1251025)
        self.assertNotEqual(beta["query_protein_sha256"], beta["reference_protein_sha256"])

    def test_pfur_bound_fna_is_valid_for_assembly_and_bin_modes(self):
        for mode in ("assembly", "bin"):
            with self.subTest(mode=mode):
                annotation = self._copy_pfur_annotation("B_PfurHII_native_locus_v2")
                manifest_path = annotation / "run_manifest.json"
                manifest = json.loads(manifest_path.read_text())
                manifest["input_mode"] = mode
                genes = list(csv.DictReader((annotation / "genes.tsv").open(), delimiter="\t"))
                for row in genes:
                    row["bin_id"] = "b1"
                _tsv(annotation / "genes.tsv", genes)
                manifest["output_sha256"]["genes.tsv"] = _sha(annotation / "genes.tsv")
                manifest_path.write_text(json.dumps(manifest))
                anchor, qid, pyridine, proofs = self._pfur_inputs(annotation)
                results, _ = prepare(annotation, [anchor],
                                     {qid: "sequence_supported_uncalibrated_noncanonical"},
                                     pyridine, root=ROOT, coding_completeness=proofs)
                self.assertEqual(results[qid]["status"], "exact_reference_coding_locus_potential")
                self.assertEqual(len(results[qid]["capabilities"]), 4)

    def test_pfur_bin_mode_mixed_required_member_does_not_complete_tuple(self):
        annotation = self._copy_pfur_annotation("B_PfurHII_native_locus_v2")
        manifest_path = annotation / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["input_mode"] = "bin"
        genes = list(csv.DictReader((annotation / "genes.tsv").open(), delimiter="\t"))
        for row in genes:
            row["bin_id"] = "b1"
        beta = next(row["protein_id"] for row in csv.DictReader(
            (annotation / "complex_members.tsv").open(), delimiter="\t")
                    if row["role"] == "SHI_beta")
        next(row for row in genes if row["id"] == beta)["bin_id"] = "b2"
        _tsv(annotation / "genes.tsv", genes)
        manifest["output_sha256"]["genes.tsv"] = _sha(annotation / "genes.tsv")
        manifest_path.write_text(json.dumps(manifest))
        anchor, qid, pyridine, proofs = self._pfur_inputs(annotation)
        results, _ = prepare(annotation, [anchor],
                             {qid: "sequence_supported_uncalibrated_noncanonical"},
                             pyridine, root=ROOT, coding_completeness=proofs)
        self.assertNotEqual(results[qid]["status"], "exact_reference_coding_locus_potential")
        self.assertEqual(results[qid]["capabilities"], [])


if __name__ == "__main__":
    unittest.main()
