"""Regression tests for recovery of an independent pyridine module.

The fixtures are the real u834437 round-33 genome, BLAST output, and gene
coordinates.  These tests deliberately do not manufacture a catalytic HSP:
the local 632/633 SHI-like beta/gamma pair is recovered while 629 remains a
competitive Group-4 catalytic candidate.
"""

import csv
import copy
import json
import unittest
from pathlib import Path

from pyridine_function import (PYRIDINE_ROLES, annotation_evidence, attributes, load_model,
                               predict, recover_module_members)
from annotate import ROOT, load_blast


HERE = Path(__file__).resolve().parents[1]
ANNOTATION = HERE / "results/hadar_pilot_round33_final/current/u834437"


class IndependentPyridineRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = load_model()
        cls.system = "NiFe_HyhNADPH"
        with (ANNOTATION / "genes.tsv").open() as handle:
            cls.genes = list(csv.DictReader(handle, delimiter="\t"))
        cls.thresholds = json.loads((HERE / "config/thresholds.json").read_text())
        cls.raw = load_blast(ANNOTATION / "blast.tsv", cls.thresholds)
        with (HERE / "results/round94_coding_locus_release_v1/u834437/domains/query_domains.tsv").open() as handle:
            cls.domains = list(csv.DictReader(handle, delimiter="\t"))
        cls.anchor_id = next(g["id"] for g in cls.genes if g["id"].endswith("_629"))
        cls.beta_id = next(g["id"] for g in cls.genes if g["id"].endswith("_633"))
        cls.gamma_id = next(g["id"] for g in cls.genes if g["id"].endswith("_632"))
        cls.delta_id = next(g["id"] for g in cls.genes if g["id"].endswith("_626"))

    def _anchor(self, mode="genome"):
        return {"protein_id": self.anchor_id, "catalytic_family": "NiFe",
                "input_mode": mode, "bin_id": ""}

    def _recover(self, genes=None, raw=None, mode="genome"):
        return recover_module_members(
            self._anchor(mode), self.raw if raw is None else raw,
            self.genes if genes is None else genes, self.model, ROOT,
            structure_ids=[self.system], role_sets={self.system: PYRIDINE_ROLES[self.system]})

    def _predictions(self, options, genes=None, raw=None, domains=None):
        """Run recovered hypotheses through the real evidence/predict stage."""
        out = []
        for option in options:
            anchor = dict(option["anchor"], catalytic_identity_confidence="calibrated_high")
            ev = annotation_evidence(anchor, option["members"],
                self.genes if genes is None else genes, "calibrated_high",
                domains=self.domains if domains is None else domains,
                model=self.model, nodes=option.get("nodes", []),
                raw_hits=self.raw if raw is None else raw)
            out.append(predict(ev, self.model))
        return out

    def test_real_632_633_pair_recovers_without_catalytic_hsp(self):
        options = self._recover()
        predictions = self._predictions(options)
        self.assertTrue(options, "real beta/gamma HSPs should seed the module-only resolver")
        members = options[0]["members"]
        selected = {m["role"]: m["protein_id"] for m in members}
        self.assertEqual(selected.get("SHI_beta"), self.beta_id)
        self.assertEqual(selected.get("SHI_gamma"), self.gamma_id)
        self.assertEqual(options[0]["provenance"]["method"],
                         "existing_module_resolver_independent_electron_module")
        self.assertTrue(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                            for p in predictions))

    def test_module_recovery_does_not_synthesize_core_or_coupling(self):
        options = self._recover()
        self.assertTrue(options)
        # The recovery object is deliberately a module hypothesis only.  The
        # complete Hyh predict() path still requires alpha/delta and linkage;
        # no catalytic HSP is inserted into the recovered members.
        roles = {m["role"] for m in options[0]["members"]}
        self.assertNotIn("SHI_alpha", roles)
        self.assertNotIn("SHI_delta", roles)
        for p in self._predictions(options):
            self.assertEqual(p["hydrogen_reaction_core_status"], "unresolved")
            self.assertEqual(p["coupling_potential_status"], "unresolved")

    def test_supported_alternative_module_is_visible_without_combining_claims(self):
        module = self._predictions(self._recover())[0]
        primary = copy.deepcopy(module)
        primary['pyridine_nucleotide_module_status'] = 'unresolved'
        primary['hydrogen_reaction_core_status'] = 'sequence_supported_candidate'
        primary['alternative_architectures'] = [module]
        rows = attributes(primary)
        self.assertTrue(any(r['attribute']=='Pyridine_nucleotide_electron_module' and
                            r['status']=='sequence_supported_candidate' and
                            r['value']['structure_id']==self.system for r in rows))
        self.assertFalse(any(r['attribute']=='Local_pyridine_hydrogenase_coupling_potential' and
                             r['status']=='local_coupling_potential_uncalibrated' for r in rows))
        self.assertEqual(primary['coupling_potential_status'], 'unresolved')
        primary['pyridine_nucleotide_module_status'] = 'sequence_supported_candidate'
        self.assertEqual(sum(r['attribute']=='Pyridine_nucleotide_electron_module' and r['status']=='sequence_supported_candidate'
                             for r in attributes(primary)), 1)

    def test_module_support_is_independent_of_catalytic_identity(self):
        option = self._recover()[0]
        ev = annotation_evidence(option['anchor'], option['members'], self.genes,
            'unresolved', self.domains, self.model, option['nodes'], self.raw)
        result = predict(ev, self.model)
        self.assertEqual(result['pyridine_nucleotide_module_status'], 'sequence_supported_candidate')
        self.assertEqual(result['hydrogen_reaction_core_status'], 'unresolved')
        self.assertEqual(result['coupling_potential_status'], 'unresolved')

    def test_missing_beta_does_not_promote_gamma_alone(self):
        raw = [h for h in self.raw if not (h["qseqid"] == self.beta_id and
                                            h["sseqid"] == "NiFe_HyhNADPH__SHI_beta")]
        predictions = self._predictions(self._recover(raw=raw), raw=raw)
        self.assertFalse(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                             for p in predictions))

    def test_incomplete_gamma_domain_hit_is_not_a_module_replacement(self):
        # Preserve the real HSP but remove one required GA domain from the
        # actual domain table; raw alignment quality alone is insufficient.
        domains = [d for d in self.domains if not
                   (d["protein_id"] == self.gamma_id and d["accession"].startswith("PF00175"))]
        predictions = self._predictions(self._recover(), domains=domains)
        self.assertFalse(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                             for p in predictions))

    def test_incompatible_full_role_rival_blocks_module_resolution(self):
        rival = dict(next(h for h in self.raw if h["qseqid"] == self.gamma_id and
                          h["sseqid"] == "NiFe_HyhNADPH__SHI_gamma"))
        rival.update(sseqid="FeFe_TmHydABC__HydB", pident=38.0,
                     passes_thresholds=True, bitscore=rival["bitscore"] + 1)
        options = self._recover(raw=[*self.raw, rival])
        predictions = self._predictions(options, raw=[*self.raw, rival])
        self.assertFalse(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                             for p in predictions))

    def test_distant_module_is_not_local_to_anchor(self):
        genes = [dict(g) for g in self.genes]
        for g in genes:
            if g["id"] in {self.beta_id, self.gamma_id}:
                g["start"] = str(int(g["start"]) + 100000)
                g["end"] = str(int(g["end"]) + 100000)
        predictions = self._predictions(self._recover(genes=genes), genes=genes)
        self.assertFalse(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                             for p in predictions))

    def test_cross_bin_module_is_not_promoted(self):
        genes = [dict(g) for g in self.genes]
        for g in genes:
            if g["id"] in {self.beta_id, self.gamma_id}:
                g["bin_id"] = "other-bin"
        predictions = self._predictions(self._recover(genes=genes, mode="bin"), genes=genes)
        self.assertFalse(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                             for p in predictions))

    def test_anchor_with_only_catalytic_hit_has_no_independent_module(self):
        raw = [h for h in self.raw if h["qseqid"] == self.anchor_id]
        predictions = self._predictions(self._recover(raw=raw), raw=raw)
        self.assertFalse(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                             for p in predictions))

    def test_partial_member_is_not_promoted(self):
        genes = [dict(g) for g in self.genes]
        for g in genes:
            if g["id"] == self.gamma_id:
                g["partial"] = "True"
        predictions = self._predictions(self._recover(genes=genes), genes=genes)
        self.assertFalse(any(p["pyridine_nucleotide_module_status"] == "sequence_supported_candidate"
                             for p in predictions))


if __name__ == "__main__":
    unittest.main()
