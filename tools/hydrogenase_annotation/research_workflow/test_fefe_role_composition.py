"""Focused FeFe fused E/F role-composition tests on the R117 Syntrophus evidence."""
import copy
import csv
import hashlib
import json
import unittest
import tempfile
from pathlib import Path

from Bio import SeqIO

from pyridine_role_composition import _fefe_fused_roles, prepare
from query_domain_organization import annotate_domains

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "results/round117_Syntrophus_first_v2/A_given_pair_v2"


def _fixture():
    annotation = CASE / "annotation"
    with (annotation / "proteins.faa").open() as handle:
        aa = {r.id: str(r.seq).upper() for r in SeqIO.parse(handle, "fasta")}
    with (annotation / "genes.tsv").open() as handle:
        genes = {}
        for row in csv.DictReader(handle, delimiter="\t"):
            genes[row["id"]] = dict(row, start=int(row["start"]), end=int(row["end"]), strand=int(row["strand"]))
    with (annotation / "catalytic_anchors.tsv").open() as handle:
        anchor = next(csv.DictReader(handle, delimiter="\t"))
    anchor = dict(anchor, input_mode="assembly", catalytic_family="FeFe")
    pyridine = json.loads((CASE / "functions/pyridine_function_evidence.json").read_text())["q02"]
    inventory = json.loads((CASE / "functions/independent_roles/inventory.json").read_text())
    return annotation, aa, genes, anchor, pyridine, inventory, json.loads((ROOT / "research_workflow/pyridine_function_model.json").read_text())


class FeFeRoleCompositionTests(unittest.TestCase):
    def test_actual_partner_controls_preserve_core_and_relay(self):
        from explain_candidates import explain
        for name in ('Syntrophus_B_beta_deleted','Syntrophus_B_beta_remote'):
            case=ROOT/'experiments/round218_research_fixtures/annotations'/('syntrophus_beta_deleted' if name.endswith('deleted') else 'syntrophus_beta_remote')
            with (case/'annotation/genes.tsv').open() as handle:
                pid=next(g['id'] for g in csv.DictReader(handle,delimiter='\t') if g['start']=='1121002')
            with tempfile.TemporaryDirectory() as tmp:
                out=Path(tmp)/'functions'
                domains=Path(tmp)/'domains'
                annotate_domains([case/'annotation'],domains)
                explain(case/'annotation',out,candidate_ids={pid},domain_evidence=domains)
                result=json.loads((out/'pyridine_function_evidence.json').read_text())[pid]
                self.assertEqual(result['hydrogen_reaction_core_status'],'sequence_supported_candidate')
                self.assertTrue(result['evidence']['internal_catalytic_relay_supported'])
                composed=json.loads((out/'pyridine_role_composition_evidence.json').read_text())[pid]
                self.assertEqual(composed['status'],'unresolved')
                self.assertIn('module_role_inventory_no_unique_local_fusion',composed['blockers'])

    def test_actual_prepare_accepts_A_and_B_and_external_resource_conflict(self):
        for case_name in ("A_given_pair_v2", "B_whole_genome_v2"):
            case = ROOT / "results/round117_Syntrophus_first_v2" / case_name
            annotation = case / "annotation"
            with (annotation / "catalytic_anchors.tsv").open() as handle:
                anchors = list(csv.DictReader(handle, delimiter="\t"))
            with (annotation / "genes.tsv").open() as handle:
                genes = list(csv.DictReader(handle, delimiter="\t"))
            pyridine = json.loads((case / "functions/pyridine_function_evidence.json").read_text())
            inventory = json.loads((case / "functions/independent_roles/inventory.json").read_text())
            model = json.loads((ROOT / "research_workflow/pyridine_function_model.json").read_text())
            source = ROOT / "research_workflow/pyridine_role_composition.py"
            for source_name in list(model["sources"]):
                if source_name.endswith("research_workflow/pyridine_role_composition.py"):
                    model["sources"][source_name] = hashlib.sha256(source.read_bytes()).hexdigest()
            result, _ = prepare(annotation, anchors, genes, [], {}, {}, pyridine, model,
                                root=ROOT, module_role_inventory=inventory)
            fe = next(a["protein_id"] for a in anchors if a["catalytic_family"] == "FeFe")
            self.assertEqual(result[fe]["status"], "local_coupling_potential_uncalibrated", case_name)
            h = result[fe]["evidence"]["hypotheses"][0]
            conflict = {"anchor": "other", "resources": {h["module_protein_id"]: h["resources"][h["module_protein_id"]][:3]}}
            result, _ = prepare(annotation, anchors, genes, [], {}, {}, pyridine, model,
                                root=ROOT, other_members=(conflict,), module_role_inventory=inventory)
            self.assertEqual(result[fe]["status"], "unresolved", case_name)
            self.assertIn("main_relay_resources_shared_with_other_anchor", result[fe]["blockers"])

    def test_actual_syntrophus_fused_roles_compose(self):
        _, aa, genes, anchor, pyridine, inventory, model = _fixture()
        result = _fefe_fused_roles("q02", anchor, {"status": "sequence_supported_candidate"},
                                   pyridine, genes, aa, inventory, set(), "assembly", model)
        self.assertEqual(result["status"], "local_coupling_potential_uncalibrated")
        h = result["evidence"]["hypotheses"][0]
        self.assertEqual(h["module_protein_id"], "q01")
        self.assertEqual({r["role_class"] for r in h["roles"]}, {"FeFe_catalytic_core", "NuoF_flavin_FeS", "NuoE_Fer2"})
        self.assertEqual(h["composition_kind"], "FeFe_independent_electron_role_composition_v2")

    def test_missing_fold_withholds_composition(self):
        _, aa, genes, anchor, pyridine, inventory, model = _fixture()
        bad = copy.deepcopy(inventory)
        role = next(r for r in bad["gene_roles"]["q01"] if r["role_class"] == "NuoF_flavin_FeS")
        role["best"]["profiles"] = []
        out = _fefe_fused_roles("q02", anchor, {"status": "sequence_supported_candidate"}, pyridine,
                                 genes, aa, bad, set(), "assembly", model)
        self.assertIn("module_role_profile_or_competition_unresolved_NuoF_flavin_FeS", out["blockers"])

    def test_overlapping_resources_withholds_composition(self):
        _, aa, genes, anchor, pyridine, inventory, model = _fixture()
        bad = copy.deepcopy(inventory)
        e = next(r for r in bad["gene_roles"]["q01"] if r["role_class"] == "NuoE_Fer2")
        f = next(r for r in bad["gene_roles"]["q01"] if r["role_class"] == "NuoF_flavin_FeS")
        e["best"]["query_positions"] = f["best"]["query_positions"][:10]
        out = _fefe_fused_roles("q02", anchor, {"status": "sequence_supported_candidate"}, pyridine,
                                 genes, aa, bad, set(), "assembly", model)
        self.assertIn("fused_module_role_resources_overlap", out["blockers"])

    def test_remote_or_partial_or_wrong_aa_withholds_composition(self):
        _, aa, genes, anchor, pyridine, inventory, model = _fixture()
        for change in ("remote", "partial", "wrong_aa"):
            g = copy.deepcopy(genes); seq = dict(aa)
            if change == "remote":
                g["q01"]["contig"] = "other_contig"
            elif change == "partial":
                g["q01"]["partial"] = "True"
            else:
                seq["q01"] = "Q" + seq["q01"][1:]
            out = _fefe_fused_roles("q02", anchor, {"status": "sequence_supported_candidate"}, pyridine,
                                     g, seq, inventory, set(), "assembly", model)
            self.assertEqual(out["status"], "unresolved", change)

    def test_competing_role_is_not_silently_ignored(self):
        _, aa, genes, anchor, pyridine, inventory, model = _fixture()
        bad = copy.deepcopy(inventory)
        role = next(r for r in bad["gene_roles"]["q01"] if r["role_class"] == "NuoE_Fer2")
        role["competition"] = [{"reference_id": "Nfn_competitor", "paired_identity": 99.0}]
        out = _fefe_fused_roles("q02", anchor, {"status": "sequence_supported_candidate"}, pyridine,
                                 genes, aa, bad, set(), "assembly", model)
        self.assertIn("module_role_profile_or_competition_unresolved_NuoE_Fer2", out["blockers"])

    def test_kist_independent_fefe_core_and_nuoef_module_compose(self):
        base = ROOT / "results/round143_KIST612_genbank_baseline_v1"
        annotation = base / "annotation"
        with (annotation / "catalytic_anchors.tsv").open() as handle:
            anchor = next(r for r in csv.DictReader(handle, delimiter="\t") if r["protein_id"] == "ELI_0847")
        anchor = dict(anchor, input_mode="genome", catalytic_family="FeFe")
        with (annotation / "genes.tsv").open() as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        with (base / "domains/query_domains.tsv").open() as handle:
            domains = list(csv.DictReader(handle, delimiter="\t"))
        pyridine = json.loads((base / "functions/pyridine_function_evidence.json").read_text())
        inventory = json.loads((base / "functions/independent_roles/inventory.json").read_text())
        features = {row["protein_id"]: row for row in json.loads((annotation / "confirmation/features.json").read_text())}
        model = json.loads((ROOT / "research_workflow/pyridine_function_model.json").read_text())
        source = ROOT / "research_workflow/pyridine_role_composition.py"
        for source_name in list(model["sources"]):
            if source_name.endswith("research_workflow/pyridine_role_composition.py"):
                model["sources"][source_name] = hashlib.sha256(source.read_bytes()).hexdigest()
        result, _ = prepare(annotation, [anchor], genes, domains,
                            {"ELI_0847": "sequence_supported_uncalibrated_noncanonical"}, {}, pyridine,
                            model, root=ROOT, module_role_inventory=inventory, identity_features=features)
        self.assertEqual(result["ELI_0847"]["status"], "local_coupling_potential_uncalibrated")
        hypothesis = result["ELI_0847"]["evidence"]["hypotheses"][0]
        self.assertIn("ELI_0843", hypothesis["protein_ids"])
        self.assertIn("ELI_0846", hypothesis["protein_ids"])

        # The independent anchor proof is deliberately tied to all three
        # ordered FeFe profiles; losing the relay-domain profile must withdraw
        # the composition rather than fall back to the old HydA search.
        no_relay = [d for d in domains if not (d.get("protein_id") == hypothesis["core_protein_id"] and d.get("accession", "").split(".")[0] == "PF10588")]
        lost, _ = prepare(annotation, [anchor], genes, no_relay,
                          {"ELI_0847": "sequence_supported_uncalibrated_noncanonical"}, {}, pyridine,
                          model, root=ROOT, module_role_inventory=inventory, identity_features=features)
        self.assertEqual(lost["ELI_0847"]["status"], "unresolved")

        tampered = [dict(d) for d in domains]
        next(d for d in tampered if d.get("protein_id") == "ELI_0847" and d.get("accession", "").split(".")[0] == "PF13510")["model_sha256"] = "tampered"
        with self.assertRaises(ValueError):
            prepare(annotation, [anchor], genes, tampered,
                    {"ELI_0847": "sequence_supported_uncalibrated_noncanonical"}, {}, pyridine,
                    model, root=ROOT, module_role_inventory=inventory, identity_features=features)

        misplaced = copy.deepcopy(features)
        misplaced["ELI_0847"]["core_positions"] = list(range(1, 30))
        wrong_span, _ = prepare(annotation, [anchor], genes, domains,
                                {"ELI_0847": "sequence_supported_uncalibrated_noncanonical"}, {}, pyridine,
                                model, root=ROOT, module_role_inventory=inventory, identity_features=misplaced)
        self.assertEqual(wrong_span["ELI_0847"]["status"], "unresolved")


if __name__ == "__main__":
    unittest.main()
