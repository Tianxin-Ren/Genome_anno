import hashlib
import inspect
import json
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fefe_profile_ligands import ROOT, _load_model, _query_domain, _site_map, _stockholm, prepare, verify_saved


class FeFeProfileLigandTests(unittest.TestCase):
    """Integration tests use the checked-in PF02906 HMM and the 8RU6 HydA source."""

    @classmethod
    def setUpClass(cls):
        cls.model, cls.hmm, cls.reference, cls.reference_file = _load_model(ROOT)
        source = ROOT / "experiments/round85/FeFe_center_experiment_v1/new_model.json"
        cls.native = json.loads(source.read_text())["sequence"]

    def _feature(self, protein_id, length, start=1, end=None):
        end = length if end is None else end
        component = {
            "accession": "PF02906.20", "profile": "Fe_hyd_lg_C",
            "model_start": 1, "model_end": 280, "model_length": 280,
            "query_start": start, "query_end": end, "query_length": length,
            "model_coverage": 1.0,
            "sequence_bitscore": 378.9, "domain_bitscore": 378.9,
            "sequence_ga": 22.1, "domain_ga": 22.1,
            "model_sha256": self.model["hmm_sha256"],
            "official_component_model_sha256": self.model["hmm_sha256"],
            "source_sha256": self.model["hmm_sha256"], "protein_id": protein_id,
        }
        return {"protein_id": protein_id, "family": "FeFe", "domain_evidence": [{
            "accession": "PF02906", "profile": "Fe_hyd_lg_C", "model_coverage": 1.0,
            "query_core_overlap": 1.0, "supported": True, "evidence_kind": "single_domain",
            "query_start": start, "query_end": end, "components": [component],
        }]}

    def _fixture(self, records):
        tmp = tempfile.TemporaryDirectory(); root = Path(tmp.name)
        annotation = root / "annotation"; annotation.mkdir()
        with (annotation / "proteins.faa").open("w") as fh:
            for pid, seq in records.items(): fh.write(f">{pid}\n{seq}\n")
        features = {pid: self._feature(pid, len(seq)) for pid, seq in records.items()}
        return tmp, annotation, root / "result", features

    def test_8ru6_native_profile_maps_four_cysteines(self):
        tmp, annotation, output, features = self._fixture({"8RU6_native": self.native})
        with tmp:
            mapping, manifest = prepare(annotation, features, output, {}, root=ROOT)
            item = mapping["8RU6_native"]
            self.assertEqual(item["profile_ligand_status"], "profile_mapped_ligand_pattern_supported")
            self.assertEqual([s["residue"] for s in item["sites"]], ["C"] * 4)
            self.assertTrue(all(s["residue_status"] == "present" for s in item["sites"]))
            self.assertTrue(all(s["posterior_status"] == "high" for s in item["sites"]))
            self.assertTrue(manifest)

    def test_controls_distinguish_substitution_ambiguity_and_deletion(self):
        c179 = self.native[:178] + "A" + self.native[179:]
        x234 = self.native[:233] + "U" + self.native[234:]
        deleted = self.native[:377] + self.native[378:]
        records = {"native": self.native, "C179A": c179, "U234": x234, "delete378": deleted}
        tmp, annotation, output, features = self._fixture(records)
        with tmp:
            mapping, _ = prepare(annotation, features, output, {}, root=ROOT)
            self.assertEqual(mapping["C179A"]["sites"][0]["residue_status"], "substituted")
            self.assertEqual(mapping["U234"]["sites"][1]["residue_status"], "ambiguous")
            self.assertEqual(mapping["delete378"]["sites"][2]["residue_status"], "not_covered")
            for pid in ("C179A", "U234", "delete378"):
                self.assertNotEqual(mapping[pid]["profile_ligand_status"], "profile_mapped_ligand_pattern_supported")

    def test_no_eligible_nife_does_not_invoke_alignment(self):
        tmp, annotation, output, _ = self._fixture({"not_hydrogenase": self.native})
        features = {"not_hydrogenase": {"protein_id": "not_hydrogenase", "family": "NiFe"}}
        with tmp, patch("fefe_profile_ligands.subprocess.run", side_effect=AssertionError("unexpected alignment")):
            mapping, manifest = prepare(annotation, features, output, {}, root=ROOT)
        self.assertEqual(mapping, {}); self.assertEqual(manifest, {})

    def test_duplicate_fasta_identifier_is_rejected(self):
        tmp = tempfile.TemporaryDirectory()
        with tmp:
            fasta = Path(tmp.name) / "proteins.faa"
            fasta.write_text(">dup\nAAAA\n>dup\nCCCC\n")
            annotation = Path(tmp.name) / "annotation"; annotation.mkdir()
            (annotation / "proteins.faa").write_text(fasta.read_text())
            with self.assertRaises(ValueError): prepare(annotation, {}, Path(tmp.name) / "out", {}, root=ROOT)

    def test_domain_acceptance_rejects_identity_and_boundary_mutations(self):
        valid = self._feature("p", len(self.native))["domain_evidence"][0]
        for key, value in (("model_sha256", "wrong"), ("sequence_ga", 100.0),
                           ("domain_ga", 100.0), ("protein_id", "other")):
            bad = json.loads(json.dumps({"domain_evidence": [valid]}))
            bad["domain_evidence"][0]["components"][0][key] = value
            self.assertEqual(_query_domain(bad, self.model, protein_id="p", query_length=len(self.native)), [])
        bad = json.loads(json.dumps({"domain_evidence": [valid]}))
        bad["domain_evidence"][0]["query_start"] = 0
        bad["domain_evidence"][0]["components"][0]["query_start"] = 0
        self.assertEqual(_query_domain(bad, self.model, protein_id="p", query_length=len(self.native)), [])

    def test_low_profile_posterior_cannot_be_supported(self):
        tmp, annotation, output, features = self._fixture({"native": self.native})
        with tmp:
            prepare(annotation, features, output, {}, root=ROOT)
            alignment = output / "PF02906.alignment.sto"
            text = alignment.read_text()
            _, _, rf = _stockholm(text)
            target_col = [i for i, mark in enumerate(rf, 1) if mark.lower() == "x"][76]
            chunks_seen = 0
            changed = []
            for line in text.splitlines():
                if line.startswith("#=GR q_00000000 PP "):
                    prefix, row = line.rsplit(" ", 1)
                    chars = list(row)
                    if chunks_seen < target_col <= chunks_seen + len(chars):
                        chars[target_col - chunks_seen - 1] = "9"
                    chunks_seen += len(chars)
                    line = prefix + " " + "".join(chars)
                changed.append(line)
            text = "\n".join(changed) + "\n"
            alignment.write_text(text)
            manifest_path = output / "PF02906.manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["output_sha256"]["PF02906.alignment.sto"] = hashlib.sha256(alignment.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
            with self.assertRaises(ValueError):
                verify_saved(output, annotation=annotation, features=features, parameters={}, root=ROOT)

    def test_stockholm_parser_concatenates_chunks_and_pp_rows(self):
        text = """# STOCKHOLM 1.0
q_1 AC-
#=GR q_1 PP **.
#=GC RF xx.
q_1 DE
#=GR q_1 PP 89
#=GC RF xx
//
"""
        seqs, pp, rf = _stockholm(text)
        self.assertEqual(seqs["q_1"], "AC-DE"); self.assertEqual(pp["q_1"], "**.89"); self.assertEqual(rf, "xx.xx")

    def test_saved_mapping_tamper_is_reconstructed_not_hash_only(self):
        tmp, annotation, output, features = self._fixture({"native": self.native})
        with tmp:
            prepare(annotation, features, output, {}, root=ROOT)
            mapping_path = output / "PF02906.query_mapping.json"
            mapping = json.loads(mapping_path.read_text()); mapping["native"]["sites"][0]["residue"] = "A"
            mapping_path.write_text(json.dumps(mapping, indent=2, sort_keys=True))
            manifest_path = output / "PF02906.manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["output_sha256"]["PF02906.query_mapping.json"] = hashlib.sha256(mapping_path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
            with self.assertRaisesRegex(ValueError, "contradicts actual alignment"):
                if "annotation" in inspect.signature(verify_saved).parameters:
                    verify_saved(output, annotation=annotation, features=features, parameters={}, root=ROOT)
                else:
                    verify_saved(output, root=ROOT)


if __name__ == "__main__": unittest.main()
