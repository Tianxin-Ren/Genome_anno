import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from .split_core_identity import SPLIT_IDENTITY, SPLIT_LABEL, prepare


ROOT = Path(__file__).resolve().parents[1]
ANNOTATION = ROOT / "experiments/round199_JJ_inputs/r198_genbank/annotation"


def _inputs():
    anchors = list(csv.DictReader((ANNOTATION / "catalytic_anchors.tsv").open(), delimiter="\t"))
    genes = list(csv.DictReader((ANNOTATION / "genes.tsv").open(), delimiter="\t"))
    sequences = {}
    current = None
    chunks = []
    for line in (ANNOTATION / "proteins.faa").open():
        if line.startswith(">"):
            if current:
                sequences[current] = "".join(chunks)
            current, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if current:
        sequences[current] = "".join(chunks)
    return anchors, genes, sequences


class SplitCoreIdentityTests(unittest.TestCase):
    def test_archived_jj_proof(self):
        anchors, genes, sequences = _inputs()
        proofs, hashes = prepare(ANNOTATION, anchors, genes, sequences, ROOT)
        proof = proofs["MMJJ_11350"]
        self.assertEqual(proof["identity"], SPLIT_IDENTITY)
        self.assertEqual(proof["label"], SPLIT_LABEL)
        self.assertEqual({m["protein_id"] for m in proof["members"]}, {"MMJJ_11350", "MMJJ_11360"})
        self.assertEqual([(m["role"], m["protein_id"]) for m in proof["members"]],
                         [("VhuA", "MMJJ_11350"), ("VhuU", "MMJJ_11360")])
        for member in proof["members"]:
            start, end = member["query_interval"]
            self.assertTrue(1 <= start <= end <= len(sequences[member["protein_id"]]))
        self.assertEqual(len(proof["actual_sites"]), 4)
        self.assertEqual(proof["coupling_status"], "unresolved")
        self.assertTrue(any(p.endswith("/catalytic_center_candidates.tsv") for p in hashes))

    def _copy_annotation(self):
        tmp = Path(tempfile.mkdtemp())
        for name in ("run_manifest.json", "catalytic_anchors.tsv", "catalytic_center_candidates.tsv", "catalytic_center_components.tsv", "active_center_summary.tsv", "genes.tsv", "proteins.faa"):
            shutil.copy2(ANNOTATION / name, tmp / name)
        return tmp

    def test_missing_or_ambiguous_abstains(self):
        anchors, genes, sequences = _inputs()
        tmp = self._copy_annotation()
        (tmp / "catalytic_center_candidates.tsv").unlink()
        self.assertEqual(prepare(tmp, anchors, genes, sequences, ROOT)[0], {})
        tmp = self._copy_annotation()
        with (tmp / "catalytic_center_candidates.tsv").open("a") as fh:
            row = next(csv.DictReader((ANNOTATION / "catalytic_center_candidates.tsv").open(), delimiter="\t"))
            fh.write("\t".join(row.values()) + "\n")
        manifest = json.loads((tmp / "run_manifest.json").read_text())
        digest = hashlib.sha256((tmp / "catalytic_center_candidates.tsv").read_bytes()).hexdigest()
        manifest["output_sha256"]["catalytic_center_candidates.tsv"] = digest
        (tmp / "run_manifest.json").write_text(json.dumps(manifest))
        self.assertEqual(prepare(tmp, anchors, genes, sequences, ROOT)[0], {})

    def test_altered_sha_is_rejected(self):
        anchors, genes, sequences = _inputs()
        tmp = self._copy_annotation()
        with (tmp / "genes.tsv").open("a") as fh:
            fh.write("\n")
        with self.assertRaises(ValueError):
            prepare(tmp, anchors, genes, sequences, ROOT)


if __name__ == "__main__":
    unittest.main()
