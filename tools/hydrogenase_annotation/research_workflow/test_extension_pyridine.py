"""Contract tests for the source-pinned HydABCSL/HydABCDEF extensions.

These tests intentionally build the same anonymous FAA/GFF and raw-search
artifacts produced by the extension search.  They are a boundary test for the
claim-specific resolver, rather than a second classifier implementation.
"""
import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "research_workflow")]
from annotate import BLAST_FIELDS  # noqa: E402
from pyridine_function import load_model, prepare_extension  # noqa: E402


digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


SYSTEMS = {
    "NiFe_Amobile_HydABCSL": {
        "roles": ["HydL", "HydS", "HydA_ET_not_H_cluster", "HydB", "HydC"],
        "core": ["HydL", "HydS"], "module": ["HydB", "HydC"],
        "relay": ["HydA_ET_not_H_cluster"], "context": [],
        "domains": {
            # Coordinates are fixed from the R98 Amobile_A domain observation.
            "HydS": [("PF01058", 14, 162)],
            "HydB": [("PF01512", 159, 331), ("PF10589", 441, 525)],
            "HydC": [("PF01257", 12, 156)],
            "HydA_ET_not_H_cluster": [("PF13510", 2, 78), ("PF10588", 85, 124), ("PF22117", 139, 208)],
        },
    },
    "FeFe_Mthermo_HydABCDEF": {
        "roles": ["HydA", "HydD", "HydE", "HydF", "HydB", "HydC"],
        "core": ["HydA"], "module": ["HydB", "HydC"],
        "relay": ["HydF"], "context": ["HydD", "HydE"],
        "domains": {
            # Coordinates are fixed from the R98 Mthermo_A domain observation.
            "HydB": [("PF01512", 154, 325), ("PF10589", 436, 520)],
            "HydC": [("PF01257", 6, 150)],
            "HydF": [("PF13510", 4, 81), ("PF10588", 85, 122), ("PF22117", 139, 206)],
        },
    },
}


class ExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = load_model()
        bank = json.loads((ROOT / "references/processed/complex_extensions/references.json").read_text())
        cls.refs = {sid: [r for r in bank if r["structure_id"] == sid] for sid in SYSTEMS}

    def fixture(self, path, sid, omit=(), partial=(), cross_bin=(), input_mode="genome", bad_aa=False, remote=(), low_coverage=False, bad_significance=False, missing_domain=False):
        spec = SYSTEMS[sid]
        ext = path / "structure_extensions"
        ext.mkdir()
        refs = {r["role"]: r for r in self.refs[sid]}
        proteins, genes, hits, domains = [], [], [], []
        for i, role in enumerate(spec["roles"], 1):
            if role in omit:
                continue
            ref, pid = refs[role], f"q{i:04d}"
            seq = ref["sequence"]
            proteins.append(f">{pid}\n{seq}\n")
            offset = 50000 if role in remote else 0
            genes.append(dict(id=pid, contig="ctg1", start=ref["start"] + offset, end=ref["end"] + offset,
                              strand=ref["strand"], bin_id=("b2" if role in cross_bin else "b1"),
                              partial=(role in partial)))
            qseq = ("A" + seq[1:]) if bad_aa and role == spec["roles"][0] else seq
            hits.append(dict(qseqid=pid, sseqid=ref["reference_id"], pident=100,
                             length=len(seq), qstart=1, qend=len(seq), sstart=1,
                             send=len(seq), evalue=0, bitscore=len(seq) * 2,
                             qlen=len(qseq), slen=len(seq), qseq=qseq, sseq=seq))
            for acc, lo, hi in spec["domains"].get(role, []):
                if missing_domain and role == "HydB" and acc == "PF10589":
                    continue
                domains.append(dict(protein_id=pid, accession=acc, query_start=lo,
                                    query_end=hi, passes="True"))
        aa = "".join(proteins)
        (path / "proteins.faa").write_text(aa)
        (ext / "queries.faa").write_text(aa)
        (path / "blast.tsv").write_text("")
        (ext / "blast.tsv").write_text("".join("\t".join(str(h[k]) for k in BLAST_FIELDS) + "\n" for h in hits))
        if low_coverage or bad_significance:
            target = next(h for h in hits if h["sseqid"] == refs[spec["core"][0]]["reference_id"])
            if low_coverage:
                n = max(1, int(len(target["sseq"]) * 0.40))
                target.update(length=n, qend=n, send=n, qseq=target["qseq"][:n], sseq=target["sseq"][:n], evalue=1e-50)
            if bad_significance:
                target["evalue"] = 1.0
            (ext / "blast.tsv").write_text("".join("\t".join(str(h[k]) for k in BLAST_FIELDS) + "\n" for h in hits))
        with (path / "genes.tsv").open("w") as f:
            writer = csv.DictWriter(f, fieldnames=list(genes[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(genes)
        thresholds = json.loads((ROOT / "config/thresholds.json").read_text())
        (path / "run_manifest.json").write_text(json.dumps(dict(thresholds=thresholds, input_mode=input_mode)))
        manifest = dict(status="post_identity_structure_extension_search_complete",
            model_manifest_sha256=digest(ROOT / "references/processed/complex_extensions/manifest.json"),
            reaction_registry_sha256=digest(ROOT / "config/reaction_contracts.json"),
            script_sha256=digest(ROOT / "scripts/complex_extensions.py"), mode=input_mode,
            eligible_query_proteins=len(genes), source_proteins=len(genes), hits=len(hits),
            query_fasta_sha256=digest(ext / "queries.faa"), raw_alignment_sha256=digest(ext / "blast.tsv"))
        (ext / "search_manifest.json").write_text(json.dumps(manifest))
        anchor_role = spec["core"][0]
        anchor_id = next(i for i, r in enumerate(spec["roles"], 1) if r == anchor_role)
        anchor = dict(protein_id=f"q{anchor_id:04d}", structure_id=sid, input_mode=input_mode,
                      catalytic_identity_confidence="calibrated_high")
        return anchor, domains

    def evaluate(self, sid, **kwargs):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t); anchor, domains = self.fixture(p, sid, **kwargs)
            return prepare_extension(p, anchor, domains, self.model, ROOT, sid)[0]

    def test_native_systems_support_local_coupling(self):
        for sid in SYSTEMS:
            with self.subTest(sid=sid):
                self.assertEqual(self.evaluate(sid)["coupling_potential_status"], "local_coupling_potential_uncalibrated")

    def test_missing_relay_keeps_core_and_module(self):
        for sid, spec in SYSTEMS.items():
            with self.subTest(sid=sid):
                r = self.evaluate(sid, omit=spec["relay"])
                self.assertEqual(r["hydrogen_reaction_core_status"], "sequence_supported_candidate")
                self.assertEqual(r["pyridine_nucleotide_module_status"], "sequence_supported_candidate")
                self.assertEqual(r["coupling_potential_status"], "unresolved")

    def test_partial_context_or_relay_only_blocks_coupling(self):
        for sid, spec in SYSTEMS.items():
            role = (spec["context"] or spec["relay"])[0]
            with self.subTest(sid=sid):
                r = self.evaluate(sid, partial=[role])
                self.assertEqual(r["hydrogen_reaction_core_status"], "sequence_supported_candidate")
                self.assertEqual(r["pyridine_nucleotide_module_status"], "sequence_supported_candidate")
                self.assertEqual(r["coupling_potential_status"], "unresolved")

    def test_cross_bin_or_remote_does_not_upgrade_coupling(self):
        for sid, spec in SYSTEMS.items():
            with self.subTest(sid=sid):
                r = self.evaluate(sid, cross_bin=spec["module"], input_mode="bin")
                self.assertEqual(r["hydrogen_reaction_core_status"], "sequence_supported_candidate")
                self.assertEqual(r["pyridine_nucleotide_module_status"], "unresolved")
                self.assertEqual(r["coupling_potential_status"], "unresolved")

    def test_remote_coordinates_keep_core_but_not_module_or_coupling(self):
        for sid, spec in SYSTEMS.items():
            with self.subTest(sid=sid):
                r = self.evaluate(sid, remote=spec["module"])
                self.assertEqual(r["hydrogen_reaction_core_status"], "sequence_supported_candidate")
                self.assertEqual(r["pyridine_nucleotide_module_status"], "unresolved")
                self.assertEqual(r["coupling_potential_status"], "unresolved")

    def test_relay_is_not_called_a_second_catalytic_hydrogenase(self):
        r = self.evaluate("NiFe_Amobile_HydABCSL")
        roles = r["evidence"]["roles"]
        self.assertIn("HydA_ET_not_H_cluster", roles)
        ref = next(x for x in self.refs["NiFe_Amobile_HydABCSL"] if x["role"] == "HydA_ET_not_H_cluster")
        self.assertEqual(ref["kind"], "module")
        self.assertEqual(ref["center"], "NiFe")
        self.assertEqual(r["structure_id"], "NiFe_Amobile_HydABCSL")

    def test_raw_alignment_hash_is_bound(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t); anchor, domains = self.fixture(p, "NiFe_Amobile_HydABCSL")
            raw = p / "structure_extensions/blast.tsv"; raw.write_text(raw.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "provenance"):
                prepare_extension(p, anchor, domains, self.model, ROOT, "NiFe_Amobile_HydABCSL")

    def test_bound_query_residue_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t); anchor, domains = self.fixture(p, "NiFe_Amobile_HydABCSL", bad_aa=True)
            raw = p / "structure_extensions/blast.tsv"
            manifest = p / "structure_extensions/search_manifest.json"
            m = json.loads(manifest.read_text()); m["raw_alignment_sha256"] = digest(raw); manifest.write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError, "sequence"):
                prepare_extension(p, anchor, domains, self.model, ROOT, "NiFe_Amobile_HydABCSL")

    def test_significant_low_coverage_catalytic_hsp_keeps_module_only(self):
        r = self.evaluate("NiFe_Amobile_HydABCSL", low_coverage=True)
        self.assertEqual(r["pyridine_nucleotide_module_status"], "sequence_supported_candidate")
        self.assertEqual(r["hydrogen_reaction_core_status"], "unresolved")
        self.assertEqual(r["coupling_potential_status"], "unresolved")

    def test_missing_module_ga_domain_does_not_support_module(self):
        r = self.evaluate("NiFe_Amobile_HydABCSL", missing_domain=True)
        self.assertEqual(r["pyridine_nucleotide_module_status"], "unresolved")
        self.assertEqual(r["coupling_potential_status"], "unresolved")

    def test_insignificant_catalytic_anchor_does_not_start_extension(self):
        r = self.evaluate("NiFe_Amobile_HydABCSL", bad_significance=True)
        self.assertIsNone(r)

    def test_real_hadar1118_extension_regression(self):
        annotation = ROOT / "experiments/round218_research_fixtures/annotations/hadar_u1055931/annotation"
        domain_path = ROOT / "results/round104_Hadar11_release_v2/u1055931/domains/query_domains.tsv"
        anchor_path = annotation / "catalytic_anchors.tsv"
        if not annotation.exists() or not domain_path.exists() or not anchor_path.exists():
            self.skipTest("Hadar1118 release fixture is not present")
        with anchor_path.open() as h:
            anchor = next(dict(r) for r in csv.DictReader(h, delimiter="\t") if r["protein_id"].endswith("_1118"))
        anchor.update(input_mode="genome", catalytic_identity_confidence="sequence_supported_uncalibrated_core")
        with domain_path.open() as h:
            domains = list(csv.DictReader(h, delimiter="\t"))
        for d in domains:
            for k in ("query_start", "query_end", "query_length", "model_length"):
                if d.get(k): d[k] = int(float(d[k]))
        result, _ = prepare_extension(annotation, anchor, domains, self.model, ROOT, "NiFe_Amobile_HydABCSL")
        self.assertIsNotNone(result)
        roles = result["evidence"]["roles"]
        self.assertEqual({role for role, rows in roles.items() if rows}, {"HydB", "HydC"})
        self.assertEqual(result["hydrogen_reaction_core_status"], "unresolved")
        self.assertEqual(result["coupling_potential_status"], "unresolved")


if __name__ == "__main__":
    unittest.main()
