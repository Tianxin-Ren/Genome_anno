import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from research_workflow.coding_completeness import _record_digest, load_evidence, prepare, supports, supports_translated_region
from research_workflow.prepare_annotated_input import prepare as prepare_annotated_input


class CodingCompletenessTests(unittest.TestCase):
    def make(self, dna="ATGGTTTAA", protein="MV", partial="False", source=True, coords=None):
        d = Path(tempfile.mkdtemp())
        ann = d / "annotation"; ann.mkdir()
        source_path = d / "input.fna"
        if source:
            source_path.write_text(">ctg\n" + dna + "\n")
        if coords is None: coords = (1, len(dna))
        (ann / "proteins.faa").write_text(f">p1\n{protein}\n")
        (ann / "genes.tsv").write_text("id\tcontig\tstart\tend\tstrand\tpartial\n" +
                                      f"p1\tctg\t{coords[0]}\t{coords[1]}\t1\t{partial}\n")
        (ann / "prodigal.gff").write_text(
            "##gff-version 3\n# Sequence Data: seqhdr\n"
            "# Model Data: something;transl_table=11;\n"
            "ctg\tProdigal\tCDS\t1\t%d\t.\t+\t0\tID=p1\n" % len(dna))
        inp = []
        if source:
            inp = [{"path": str(source_path), "sha256": hashlib.sha256(source_path.read_bytes()).hexdigest()}]
        (ann / "run_manifest.json").write_text(json.dumps({"status": "complete", "input_files": inp, "input_mode": "fna"}))
        return d, ann

    def test_genbank_origin_spanning_cds_exports_unwrapped_gff(self):
        d = Path(tempfile.mkdtemp()); source = d / "source.gb"
        record = SeqRecord(Seq("GTTTAACCCATG"), id="ctg", name="ctg",
                           description="circular test")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "circular"
        record.features = [SeqFeature(CompoundLocation([
            SimpleLocation(9, 12, strand=1), SimpleLocation(0, 6, strand=1)], operator="join"),
            type="CDS", qualifiers={"locus_tag": ["p1"], "protein_id": ["p1"],
                                      "translation": ["MV"], "transl_table": ["11"]})]
        SeqIO.write(record, source, "genbank")
        out = d / "prepared"
        manifest = prepare_annotated_input(source, out)
        gff = (out / "genes.gff").read_text()
        self.assertEqual(manifest["status"], "explicit_source_translation_prepared")
        self.assertIn("##topology ctg circular", gff)
        self.assertIn("ctg\texplicit_source_CDS\tCDS\t10\t18\t", gff)
        ann = d / "annotation"; ann.mkdir()
        (ann / "proteins.faa").write_text(">p1\nMV\n")
        (ann / "genes.tsv").write_text(
            "id\tcontig\tstart\tend\tstrand\tpartial\ttopology\tcontig_length\n"
            "p1\tctg\t10\t18\t1\tFalse\tcircular\t12\n")
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        (ann / "run_manifest.json").write_text(json.dumps({
            "status": "complete", "input_mode": "assembly",
            "input_files": [{"path": str(source), "sha256": source_hash}]}))
        rows, _ = prepare(ann, d / "coding-out")
        self.assertEqual(rows["p1"]["status"], "complete_under_bound_translation_model")
        self.assertEqual(rows["p1"]["dna_length"], 9)
        proof = rows["p1"]
        self.assertTrue(supports(proof, "p1", proof["protein_sha256"],
                                 feature={"contig": "ctg", "start": "10", "end": "18", "strand": "1"},
                                 gene={"contig": "ctg", "start": "10", "end": "18", "strand": "1",
                                       "topology": "circular", "contig_length": "12"}))
        (ann / "genes.tsv").write_text(
            "id\tcontig\tstart\tend\tstrand\tpartial\ttopology\tcontig_length\n"
            "p1\tctg\t0\t0\t0\tFalse\tlinear\t12\n")
        rows, _ = prepare(ann, d / "linear-out")
        self.assertNotEqual(rows["p1"]["status"], "complete_under_bound_translation_model")

        record.seq = Seq("ATGGTTTAACCC")
        record.features[0].location = SimpleLocation(0, 9, strand=1)
        SeqIO.write(record, source, "genbank")
        manifest = json.loads((ann / "run_manifest.json").read_text())
        manifest["input_files"][0]["sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        (ann / "run_manifest.json").write_text(json.dumps(manifest))
        (ann / "genes.tsv").write_text(
            "id\tcontig\tstart\tend\tstrand\tpartial\ttopology\tcontig_length\n"
            "p1\tctg\t1\t9\t1\tFalse\tcircular\t12\n")
        rows, _ = prepare(ann, d / "ordinary-circular-out")
        proof = rows["p1"]
        self.assertTrue(supports(proof, "p1", proof["protein_sha256"]))

    def run_case(self, **kwargs):
        d, ann = self.make(**kwargs)
        rows, _ = prepare(ann, d / "out")
        return d, rows["p1"]

    def test_explicit_source_sec_proof_is_retranslated_and_bound(self):
        for topology, dna, location, qualifier, gene_start, gene_end, strand, positions, annotated in (
            ("linear", "ATGTGAGTGTAA", SimpleLocation(0, 12, strand=1),
             "(pos:4..6,aa:OTHER)", 1, 12, 1, [4, 5, 6], "MXV"),
            ("circular", "GAGTTTAACCCATGT", CompoundLocation([
                SimpleLocation(11, 15, strand=1), SimpleLocation(0, 8, strand=1)]),
             "(pos:join(15,1..2),aa:Sec)", 12, 23, 1, [15, 1, 2], "MUV"),
            ("circular", str(Seq("GAGTTTAACCCATGT").reverse_complement()), CompoundLocation([
                SimpleLocation(0, 4, strand=-1), SimpleLocation(7, 15, strand=-1)]),
             "(pos:complement(join(14..15,1)),aa:Sec)", 8, 19, -1, [1, 15, 14], "MUV"),
        ):
            with self.subTest(topology=topology, strand=strand), tempfile.TemporaryDirectory() as tmp:
                d = Path(tmp)
                ann = d / "annotation"
                ann.mkdir()
                source = d / "source.gb"
                record = SeqRecord(Seq(dna), id="ctg", name="ctg", description="annotated Sec test",
                                   annotations={"molecule_type": "DNA", "topology": topology})
                q = {"locus_tag": ["p1"], "translation": [annotated], "transl_table": ["11"],
                     "transl_except": [qualifier], "note": ["selenocysteine at 4..6"]}
                feature = SeqFeature(location, type="CDS", qualifiers=q)
                record.features = [feature]
                (ann / "genes.tsv").write_text(
                    "id\tcontig\tstart\tend\tstrand\tpartial\ttopology\tcontig_length\n"
                    f"p1\tctg\t{gene_start}\t{gene_end}\t{strand}\tFalse\t{topology}\t{len(dna)}\n")
                (ann / "proteins.faa").write_text(">p1\nMUV\n")

                def read_proof(name):
                    SeqIO.write(record, source, "genbank")
                    (ann / "run_manifest.json").write_text(json.dumps({"status": "complete",
                        "input_files": [{"path": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}]}))
                    rows, _ = prepare(ann, d / name)
                    return rows["p1"]

                proof = read_proof("verified")
                self.assertTrue(supports(proof, "p1", proof["protein_sha256"]))
                self.assertEqual(proof["translation_match_mode"], "explicit_genbank_sec")
                self.assertEqual(proof["source_sec_sites"][0]["genomic_positions_1based"], positions)
                self.assertEqual(proof["raw_stop_positions"], [2, 4])
                self.assertEqual(proof["stop_positions"], [4])
                self.assertEqual(proof["coding_completeness_scope"], "protein_coding_completeness_only")
                self.assertFalse(supports(proof, "p1", gene=dict(
                    contig="ctg", start=gene_start, end=gene_end, strand=strand, partial=True,
                    topology=topology, contig_length=len(dna))))
                for changed in (
                    dict(proof, source_sec_sites=[]), dict(proof, raw_stop_positions=[4]),
                    dict(proof, translation_match_mode="normalized_initiator_M", source_sec_sites=[]),
                    dict(proof, raw_translated_protein_sha256="0" * 64),
                ):
                    changed["record_digest"] = _record_digest(changed)
                    self.assertFalse(supports(changed, "p1", proof["protein_sha256"]))
                feature.qualifiers["transl_except"] = [qualifier, qualifier]
                self.assertFalse(supports(read_proof("duplicate"), "p1"))
                feature.qualifiers.pop("transl_except")
                self.assertFalse(supports(read_proof("unannotated"), "p1"))
                feature.qualifiers["transl_except"] = [qualifier]
                feature.qualifiers["pseudo"] = [""]
                self.assertFalse(supports(read_proof("pseudo"), "p1"))
        _, plain_fna = self.run_case(dna="ATGTGAGTTTAA", protein="MUV")
        self.assertEqual(plain_fna["reasons"], ["internal_stop"])

    def test_canonical_and_alternative_starts_normalize_to_m(self):
        for codon in ("ATG", "GTG", "TTG", "ATT"):
            d, row = self.run_case(dna=codon + "GTTTAA", protein="MV")
            self.assertEqual(row["status"], "complete_under_bound_translation_model", d)
            self.assertEqual(row["source_partial"], False)
            self.assertEqual(row["start_codon"], codon)

    def test_partial_is_retained_but_complete_proof_is_separate(self):
        _, row = self.run_case(partial="True")
        self.assertTrue(row["source_partial"])
        self.assertEqual(row["status"], "complete_under_bound_translation_model")

    def test_negative_coding_cases(self):
        cases = [
            ("GTTTAA", "V", "canonical_start_missing"),
            ("ATGGTT", "MV", "terminal_stop_missing"),
            ("ATGT AATAA".replace(" ", ""), "MN", "internal_stop"),
            ("ATGGTTT", "MV", "frameshift_cds_length_not_multiple_of_3"),
        ]
        for dna, protein, reason in cases:
            _, row = self.run_case(dna=dna, protein=protein)
            self.assertEqual(row["status"], "incomplete", reason)
            self.assertIn(reason, row["reasons"])

    def test_translated_region_support_is_only_missing_terminal_stop(self):
        d, row = self.run_case(dna="ATGGTT", protein="MV")
        self.assertTrue(supports_translated_region(row, "p1", "MV",
                                                   gene={"contig": "ctg", "start": "1", "end": "6", "strand": "1"}))
        self.assertFalse(supports_translated_region(row, "p1", "ML"))
        self.assertFalse(supports_translated_region(dict(row, reasons=["internal_stop"]), "p1", "MV"))
        source = d / "input.fna"
        source.write_text(">ctg\nATGGTT\nchanged\n")
        self.assertFalse(supports_translated_region(row, "p1", "MV"))

    def test_translated_region_rechecks_binding_even_with_recomputed_record_digest(self):
        from research_workflow.coding_completeness import _record_digest
        _, row = self.run_case(dna="ATGGTT", protein="MV", partial="True")
        self.assertFalse(supports(row, "p1", row['protein_sha256']))
        for key, value in [('coding_completeness_scope', 'activity'),
                           ('input_dna_sha256', '0' * 64), ('protein_length', 3),
                           ('start_codon', 'GTG'), ('stop_positions', [2])]:
            changed = dict(row, **{key: value})
            changed['record_digest'] = _record_digest(changed)
            self.assertFalse(supports_translated_region(changed, 'p1', 'MV'), key)
        self.assertFalse(supports_translated_region(row, 'p1', 'MV', gene={'start': 2}))
        _, row = self.run_case(dna='GTGGTT', protein='VV', partial='True')
        self.assertTrue(supports_translated_region(row, 'p1', 'VV'))
        changed = dict(row, source_partial=False)
        changed['record_digest'] = _record_digest(changed)
        self.assertFalse(supports_translated_region(changed, 'p1', 'VV'))

    def test_translation_mismatch(self):
        _, row = self.run_case(protein="ML")
        self.assertEqual(row["status"], "incomplete")
        self.assertIn("protein_translation_mismatch", row["reasons"])

    def test_partial_legal_alternate_start_matches_raw_translation(self):
        for codon, raw_aa in (("GTG", "V"), ("TTG", "L"), ("ATT", "I")):
            _, row = self.run_case(dna=codon + "GTTTAA", protein=raw_aa + "V", partial="True")
            self.assertEqual(row["status"], "complete_under_bound_translation_model")
            self.assertEqual(row["translation_match_mode"], "raw_legal_alternate_initiator")
            self.assertEqual(row["translated_protein_sha256"], row["protein_sha256"])
            self.assertNotEqual(row["raw_translated_protein_sha256"], row["normalized_translated_protein_sha256"])

    def test_supports_rejects_tampered_translation_mode(self):
        from research_workflow.coding_completeness import supports
        _, row = self.run_case(dna="GTGGTTTAA", protein="VV", partial="True")
        self.assertTrue(supports(row, "p1", row["protein_sha256"]))
        tampered = dict(row, translation_match_mode="normalized_initiator_M")
        tampered["record_digest"] = __import__("research_workflow.coding_completeness", fromlist=["_record_digest"])._record_digest(tampered)
        self.assertFalse(supports(tampered, "p1", row["protein_sha256"]))

    def test_partial_legal_alternate_start_can_still_match_normalized_M(self):
        _, row = self.run_case(dna="GTGGTTTAA", protein="MV", partial="True")
        self.assertEqual(row["status"], "complete_under_bound_translation_model")
        self.assertEqual(row["translation_match_mode"], "normalized_initiator_M")

    def test_alternate_start_wrong_first_or_internal_residue_is_rejected(self):
        _, row = self.run_case(dna="GTGGTTTAA", protein="ML", partial="True")
        self.assertEqual(row["status"], "incomplete")
        _, row = self.run_case(dna="GTGGCCTAA", protein="MV", partial="True")
        self.assertEqual(row["status"], "incomplete")

    def test_noncanonical_start_is_rejected_even_if_raw_matches(self):
        _, row = self.run_case(dna="GCGGTTTAA", protein="AV", partial="True")
        self.assertEqual(row["status"], "incomplete")

    def test_faa_only_is_unknown(self):
        _, row = self.run_case(source=False)
        self.assertEqual(row["status"], "unknown")
        self.assertIn("source_nucleotide_unavailable", row["reasons"])

    def test_declared_source_hash_mismatch_fails_closed(self):
        d, ann = self.make()
        source = d / "input.fna"
        source.write_text(">ctg\nATGGTTTAA\nchanged\n")
        with self.assertRaises(ValueError):
            prepare(ann, d / "out")

    def test_load_rejects_rewritten_bound_source(self):
        d, row = self.run_case()
        source = d / "input.fna"
        source.write_text(">ctg\nATGGTTTAA\nchanged\n")
        with self.assertRaises(ValueError):
            load_evidence(d / "out", annotation=d / "annotation")

    def test_load_rejects_wrong_owner_missing_input_and_tamper(self):
        d, _ = self.run_case()
        with self.assertRaises(ValueError):
            load_evidence(d / "out", annotation=d / "other")
        source = d / "input.fna"; source.unlink()
        with self.assertRaises(ValueError):
            load_evidence(d / "out", annotation=d / "annotation")
        d, _ = self.run_case()
        manifest = json.loads((d / "out" / "coding_completeness_manifest.json").read_text())
        manifest["script_sha256"] = "0" * 64
        (d / "out" / "coding_completeness_manifest.json").write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            load_evidence(d / "out", annotation=d / "annotation")

    def test_supports_requires_digest_and_current_identity(self):
        d, row = self.run_case(partial="True")
        self.assertTrue(supports(row, "p1", row["protein_sha256"], gene={"contig": "ctg", "start": 1, "end": 9, "strand": 1}))
        for key, value in (("protein_id", "other"), ("protein_sha256", "0" * 64), ("translation_table", 4), ("start", 2)):
            tampered = dict(row); tampered[key] = value
            # Digest is deliberately stale; even a recomputed digest cannot bypass identity checks.
            self.assertFalse(supports(tampered, "p1", row["protein_sha256"]))
        no_digest = dict(row); no_digest.pop("record_digest")
        self.assertFalse(supports(no_digest, "p1", row["protein_sha256"]))
        incomplete = dict(row); incomplete["status"] = "incomplete"
        self.assertFalse(supports(incomplete, "p1", row["protein_sha256"]))

    def test_missing_source_is_unknown_but_scan_mutation_fails(self):
        d, ann = self.make()
        source = d / "input.fna"; source.unlink()
        (ann / "run_manifest.json").write_text(json.dumps({"status": "complete", "input_files": [{"path": str(source), "sha256": "0" * 64}]}))
        # A missing source is a reportable unknown, not a fabricated complete proof.
        rows, _ = prepare(ann, d / "out")
        self.assertEqual(rows["p1"]["status"], "unknown")
        d, ann = self.make()
        original = __import__("research_workflow.coding_completeness", fromlist=["_read_fasta"])._read_fasta
        def mutate(path):
            value = original(path)
            if Path(path) == d / "input.fna":
                Path(path).write_text(">ctg\nATGGTTTAA\nchanged\n")
            return value
        with patch("research_workflow.coding_completeness._read_fasta", side_effect=mutate):
            with self.assertRaises(RuntimeError):
                prepare(ann, d / "out")


if __name__ == "__main__":
    unittest.main()
