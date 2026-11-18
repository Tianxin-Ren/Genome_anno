import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from dna_rescue import (_circular_lengths, _fasta, _genbank_subject, _normalize_seam_row,
                        _oriented_translation, _parse_hsp, _representation, load_candidates, prepare,
                        _least_rotation_offset, _native_hsp_span)


class DNARescueTests(unittest.TestCase):
    def test_iupac_inputs_preserve_bases_length_and_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequence = 'ACGTryswkmbdhvn'
            source = root/'input.fna'
            source.write_text('>c\n'+sequence+'\n')
            self.assertEqual(_fasta(source), {'c': sequence.upper()})
            record = SeqRecord(Seq(sequence), id='c', annotations={'molecule_type':'DNA'})
            genbank = root/'input.gb'
            SeqIO.write([record], genbank, 'genbank')
            subject = root/'subject.fna'
            _genbank_subject(genbank, subject)
            self.assertEqual(_fasta(subject), {'c': sequence.upper()})
            self.assertEqual(source.read_text(), '>c\n'+sequence+'\n')

    def test_gzip_fasta_is_read_with_the_same_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.fna.gz"
            with gzip.open(source, "wt") as handle:
                handle.write(">c\nACGTN\n")
            self.assertEqual(_fasta(source), {"c": "ACGTN"})

    def test_non_dna_symbols_are_still_rejected_in_both_formats(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for invalid in ('U','X','-','?','.'):
                with self.subTest(symbol=invalid):
                    source = root/'input.fna'
                    source.write_text('>c\nATG'+invalid+'AA\n')
                    with self.assertRaisesRegex(ValueError, 'invalid alphabet'):
                        _fasta(source)
                    record = SeqRecord(Seq('ATG'+invalid+'AA'), id='c', annotations={'molecule_type':'DNA'})
                    genbank = root/'input.gb'
                    SeqIO.write([record], genbank, 'genbank')
                    with self.assertRaisesRegex(ValueError, 'invalid nucleotide alphabet'):
                        _genbank_subject(genbank, root/'subject.fna')

    def test_iupac_translation_preserves_resolvable_unknown_and_stop_codons(self):
        sequence = 'GGYAARATHTGYMRNWSNBDVNNNTRATAR'
        expected = 'GKICXXXX**'
        for dna, frame in ((sequence, 1), (str(Seq(sequence).reverse_complement()), -1)):
            aa, nt, span = _oriented_translation(dna, 1, len(dna), frame, 11)
            self.assertEqual((aa, nt, span), (expected, sequence, (1, len(sequence))))

    def test_translation_tables_and_reverse_frame(self):
        aa11, _, _ = _oriented_translation("ATGTGATAA", 1, 9, 1, 11)
        aa4, _, _ = _oriented_translation("ATGTGATAA", 1, 9, 1, 4)
        self.assertEqual(aa11, "M**")
        self.assertEqual(aa4, "MW*")
        forward, _, _ = _oriented_translation("TTACAT", 1, 6, -1, 11)
        self.assertEqual(forward, "M*")

    def test_parser_rejects_bad_owner_and_query_alignment(self):
        refs = {"r": "MKK"}
        contigs = {"c": "ATGAAAAAA"}
        line = "r\tc\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1"
        row = _parse_hsp(line, refs, contigs)
        self.assertEqual(row["sframe"], 1)
        with self.assertRaisesRegex(ValueError, "query alignment"):
            _parse_hsp(line.replace("MKK\tMKK\t1", "MKM\tMKK\t1"), refs, contigs)
        with self.assertRaisesRegex(ValueError, "query or subject"):
            _parse_hsp(line.replace("\tc\t", "\tunknown\t"), refs, contigs)
        with self.assertRaisesRegex(ValueError, "frame phase"):
            _parse_hsp(line[:-1] + "2", refs, contigs)

    def test_representation_uses_dna_offset_and_reverse_gene_orientation(self):
        genes = {
            "plus": {"id": "plus", "contig": "c", "start": "10", "end": "24", "strand": "1"},
            "minus": {"id": "minus", "contig": "c", "start": "30", "end": "44", "strand": "-1"},
        }
        proteins = {"plus": "MKKKK", "minus": "MKKKK"}
        row = {"sseqid": "c", "sstart": 13, "send": 18, "sframe": 1, "translation": "KK"}
        self.assertEqual(_representation(row, genes, proteins, {"plus"}), "plus")
        reverse = {"sseqid": "c", "sstart": 36, "send": 41, "sframe": -1, "translation": "KK"}
        self.assertEqual(_representation(reverse, genes, proteins, {"minus"}), "minus")
        wrong_frame = dict(row, sframe=-1)
        self.assertEqual(_representation(wrong_frame, genes, proteins, {"plus"}), "")
        wrong_offset = dict(row, sstart=12, send=17)
        self.assertEqual(_representation(wrong_offset, genes, proteins, {"plus"}), "")

    def test_start_codon_normalization_is_local_and_table_aware(self):
        genes = {
            "plus": {"id": "plus", "contig": "p", "start": "1", "end": "9", "strand": "1"},
            "minus": {"id": "minus", "contig": "m", "start": "1", "end": "9", "strand": "-1"},
            "internal": {"id": "internal", "contig": "i", "start": "1", "end": "12", "strand": "1"},
            "att": {"id": "att", "contig": "a", "start": "1", "end": "9", "strand": "1"},
        }
        proteins = {"plus": "MKK", "minus": "MFF", "internal": "MKVK", "att": "MKK"}
        plus = {"sseqid": "p", "sstart": 1, "send": 9, "sframe": 1, "translation": "VKK"}
        self.assertEqual(_representation(plus, genes, proteins, {"plus"}, {"p": "GTGAAAAAA"}, 11), "plus")
        self.assertTrue(plus["represented_start_codon_normalized"])
        minus = {"sseqid": "m", "sstart": 1, "send": 9, "sframe": -1, "translation": "MFF"}
        # TTG is a canonical start in table 11; the raw reverse translation is LFF.
        minus["translation"] = "LFF"
        self.assertEqual(_representation(minus, genes, proteins, {"minus"}, {"m": "AAAAAACAA"}, 11), "minus")
        self.assertTrue(minus["represented_start_codon_normalized"])
        internal = {"sseqid": "i", "sstart": 7, "send": 12, "sframe": 1, "translation": "VK"}
        self.assertEqual(_representation(internal, genes, proteins, {"internal"}, {"i": "GTGAAAGTTAAA"}, 11), "internal")
        self.assertFalse(internal["represented_start_codon_normalized"])
        att = {"sseqid": "a", "sstart": 1, "send": 9, "sframe": 1, "translation": "IKK"}
        self.assertEqual(_representation(att, genes, proteins, {"att"}, {"a": "ATTAAAAAA"}, 11), "att")
        self.assertTrue(att["represented_start_codon_normalized"])

    def test_circular_origin_hsps_map_to_unwrapped_called_cds_on_both_strands(self):
        genes = {
            "plus": {"contig": "c", "start": "10", "end": "15", "strand": "1",
                     "topology": "circular", "contig_length": "12"},
            "minus": {"contig": "c", "start": "10", "end": "15", "strand": "-1",
                      "topology": "circular", "contig_length": "12"},
        }
        proteins = {"plus": "MK", "minus": "MK"}
        contigs = {"c": "A" * 12}
        plus_edge = {"sseqid": "c", "sstart": 10, "send": 12, "sframe": 1, "translation": "M"}
        plus_origin = {"sseqid": "c", "sstart": 1, "send": 3, "sframe": 1, "translation": "K"}
        minus_edge = {"sseqid": "c", "sstart": 10, "send": 12, "sframe": -1, "translation": "K"}
        minus_origin = {"sseqid": "c", "sstart": 1, "send": 3, "sframe": -1, "translation": "M"}
        self.assertEqual(_representation(plus_edge, genes, proteins, {"plus"}, contigs, 11), "plus")
        self.assertEqual(_representation(plus_origin, genes, proteins, {"plus"}, contigs, 11), "plus")
        self.assertEqual(_representation(minus_edge, genes, proteins, {"minus"}, contigs, 11), "minus")
        self.assertEqual(_representation(minus_origin, genes, proteins, {"minus"}, contigs, 11), "minus")

    def test_circular_origin_requires_declared_length_and_rejects_invalid_span(self):
        genes = {"x": {"contig": "c", "start": "10", "end": "15", "strand": "1",
                        "topology": "circular", "contig_length": "11"}}
        row = {"sseqid": "c", "sstart": 1, "send": 3, "sframe": 1, "translation": "K"}
        self.assertEqual(_representation(row, genes, {"x": "MK"}, {"x"}, {"c": "A" * 12}, 11), "")
        genes["x"]["contig_length"] = "12"
        genes["x"]["end"] = "23"
        self.assertEqual(_representation(row, genes, {"x": "MK"}, {"x"}, {"c": "A" * 12}, 11), "")

    def test_seam_hsp_keeps_raw_coordinates_and_adds_native_one_turn_span(self):
        row = {"sstart": 10, "send": 16, "slen": 24}
        _normalize_seam_row(row, 12)
        self.assertEqual((row["native_start"], row["native_end"]), (10, 16))
        self.assertTrue(row["circular_seam"])
        self.assertEqual((row["seam_sstart"], row["seam_send"], row["seam_slen"]), (10, 16, 24))

    def test_seam_hsp_skips_second_copy_and_non_crossing_hits(self):
        duplicate = {"sstart": 14, "send": 20, "slen": 24}
        self.assertIsNone(_normalize_seam_row(duplicate, 12))
        non_crossing = {"sstart": 2, "send": 8, "slen": 24}
        self.assertIsNone(_normalize_seam_row(non_crossing, 12))
        too_long = {"sstart": 10, "send": 23, "slen": 24}
        self.assertIsNone(_normalize_seam_row(too_long, 12))
        malformed = {"sstart": 0, "send": 13, "slen": 24}
        with self.assertRaisesRegex(ValueError, "Invalid circular seam"):
            _normalize_seam_row(malformed, 12)

    def test_circular_search_cut_is_deterministic_and_zero_based(self):
        self.assertEqual(_least_rotation_offset("BCABCA"), 2)
        self.assertEqual(_least_rotation_offset("AAAA"), 0)
        self.assertEqual(_least_rotation_offset("CAGT"), 1)

    def test_native_hsp_span_preserves_strand_and_recomputes_phase(self):
        plus = {"sstart": 10, "send": 16, "sframe": 1, "slen": 24}
        minus = {"sstart": 16, "send": 10, "sframe": -3, "slen": 24}
        self.assertEqual(_native_hsp_span(plus), (10, 16, 1))
        self.assertEqual(_native_hsp_span(minus), (10, 16, -3))
        self.assertEqual(_native_hsp_span(plus, offset=1), (11, 17, 2))
        self.assertEqual(_native_hsp_span(minus, offset=16), (2, 8, -2))

    def test_canonical_search_binds_native_circular_hsps_and_candidates(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)
        source.write_text(">c\nCTTATGAAAAAA\n>uncalled\nATGAAAAAA\n")
        run_path = annotation / "run_manifest.json"
        run = json.loads(run_path.read_text())
        run["input_files"][0]["sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        run_path.write_text(json.dumps(run))
        (annotation / "genes.tsv").write_text(
            "id\tcontig\tstart\tend\tstrand\ttopology\tcontig_length\n"
            "prot\tc\t4\t12\t1\tcircular\t12\n")

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            if subject.name == "nucleotide_search_subject.fna":
                self.assertEqual(_fasta(subject), {"c": "AAAAAACTTATG", "uncalled": "ATGAAAAAA"})
                out.write_text("r\tc\t20\t3\t1\t3\t1\t9\t1e-10\t20\t3\t12\tMKK\tKKL\t1\n"
                               "r\tc\t20\t3\t1\t3\t9\t1\t1e-10\t20\t3\t12\tMKK\tKFF\t-1\n")
            else:
                self.assertEqual(_fasta(subject), {"c": "AAAAAACTTATG" * 2})
                out.write_text("r\tc\t100\t3\t1\t3\t10\t18\t1e-10\t20\t3\t24\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log")
            log.write_text("fake\n")
            return ["tblastn", "-subject", str(subject), "-db_gencode", str(table)], log

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            manifest = prepare(annotation, screen, out, root=root)
        self.assertEqual(manifest["circular_search_offsets"], {"c": 6})
        self.assertEqual(manifest["represented_hsp_count"], 2)
        _, candidates = load_candidates(out, root=root)
        self.assertEqual(len(candidates), 2)
        self.assertEqual(manifest["distinct_native_span_count"], 2)
        self.assertTrue(all(r["translation_tables"] == [4, 11] for r in candidates))
        self.assertEqual({(r["start"], r["end"], r["frame"]) for r in candidates}, {(7, 15, 1), (7, 15, -1)})

        path = out / "dna_rescue_manifest.json"
        tampered = dict(manifest, circular_search_offsets={"c": 7})
        path.write_text(json.dumps(tampered))
        with self.assertRaisesRegex(ValueError, "offsets differ"):
            load_candidates(out, root=root)

        subject = out / "nucleotide_search_subject.fna"
        subject.write_text(subject.read_text().replace("AAAAAACTTATG", "TAAAAACTTATG"))
        manifest["output_sha256"][subject.name] = hashlib.sha256(subject.read_bytes()).hexdigest()
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "subject differs"):
            load_candidates(out, root=root)

    def test_circular_lengths_requires_explicit_matching_native_length(self):
        genes = {"p": {"contig": "c", "topology": "circular", "contig_length": "12"}}
        self.assertEqual(_circular_lengths(genes, {"c": "A" * 12}), {"c": 12})
        with self.assertRaisesRegex(ValueError, "differs"):
            _circular_lengths(genes, {"c": "A" * 11})

    def _fixture(self, source_suffix=".fna"):
        tmp = tempfile.TemporaryDirectory()
        root = Path(tmp.name)
        (root / "references/processed").mkdir(parents=True)
        (root / "references/processed/references.json").write_text(json.dumps([{"reference_id": "r", "kind": "catalytic", "center": "NiFe"}]))
        (root / "references/processed/references.faa").write_text(">r\nMKK\n")
        (root / "references/processed/catalytic_families.hmm").write_text("fake model\n")
        (root / "references/processed/profile_manifest.json").write_text(json.dumps([{"profile": "fake", "family": "FeFe"}]))
        (root / "config").mkdir()
        (root / "config/discovery_evidence.json").write_text(json.dumps({"evalue_reference_sequences": 1, "max_extra_query_gap_model_fraction": 0.5}))
        (root / "config/thresholds.json").write_text(json.dumps({"hmm_domain_evalue": 1e-5, "hmm_coverage": 0.55}))
        (root / "scripts").mkdir()
        (root / "scripts/discovery_evidence.py").write_text("# test\n")
        binary = root / ".env/bin"; binary.mkdir(parents=True)
        tblastn = binary / "tblastn"; tblastn.write_text("#!/bin/sh\necho 'tblastn 2.17.0 test'\n"); tblastn.chmod(0o755)
        hmmsearch = binary / "hmmsearch"
        hmmsearch.write_text("#!/bin/sh\ndom=; tbl=; while [ $# -gt 0 ]; do if [ \"$1\" = \"--domtblout\" ]; then dom=$2; shift 2; elif [ \"$1\" = \"--tblout\" ]; then tbl=$2; shift 2; else shift; fi; done; : > \"$dom\"; : > \"$tbl\"\n")
        hmmsearch.chmod(0o755)
        annotation = root / "annotation"; annotation.mkdir()
        source = root / ("input" + source_suffix)
        source.write_text(">c\nATGAAAAAA\n>uncalled\nATGAAAAAA\n")
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        (annotation / "prodigal.gff").write_text("# Model Data: version=Prodigal;transl_table=11;\n")
        (annotation / "run_manifest.json").write_text(json.dumps({"status": "complete", "input_mode": "genome", "input_files": [{"path": str(source), "sha256": digest}]}))
        (annotation / "genes.tsv").write_text("id\tcontig\tstart\tend\tstrand\nprot\tc\t1\t9\t1\n")
        (annotation / "proteins.faa").write_text(">prot\nMKK\n")
        screen = root / "screen"; screen.mkdir()
        (screen / "candidate_screen_manifest.json").write_text(json.dumps({"status": "complete", "annotation": str(annotation)}))
        (screen / "candidate_screen.tsv").write_text("protein_id\nprot\n")
        return tmp, root, annotation, screen, source

    def test_prepare_records_represented_and_unpredicted_dna(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tc\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n" +
                           "r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            manifest = prepare(annotation, screen, out, root=root)
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual(manifest["represented_hsp_count"], 2)
        self.assertEqual(manifest["candidate_count"], 1)
        self.assertEqual(manifest["distinct_native_span_count"], 1)
        m, candidates = load_candidates(out, root=root)
        self.assertEqual(m["status"], "complete")
        self.assertEqual(candidates[0]["translation_tables"], [4, 11])
        self.assertEqual(len(candidates[0]["source_hsp_ids"]), 2)
        self.assertTrue(all(c["identity_status"] == "unconfirmed" for c in candidates))
        self.assertTrue(all(c["function_prediction"] == "potential_homolog_function_unresolved" for c in candidates))
        self.assertTrue(all(c["scope"] == "translated_DNA_interval_not_reconstructed_CDS" for c in candidates))

        candidates[0]["translation_tables"] = [4]
        candidate_path = out / "dna_rescue_candidates.json"
        candidate_path.write_text(json.dumps(candidates))
        m["output_sha256"][candidate_path.name] = hashlib.sha256(candidate_path.read_bytes()).hexdigest()
        (out / "dna_rescue_manifest.json").write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError, "table provenance"):
            load_candidates(out, root=root)

    def test_same_span_keeps_distinct_table_peptides(self):
        tmp, root, _, _, source = self._fixture()
        self.addCleanup(tmp.cleanup)
        source.write_text(">uncalled\nATGTGAAAA\n")

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            peptide = "MWK" if table == 4 else "M*K"
            out.write_text(f"r\tuncalled\t66.7\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\t{peptide}\t1\n")
            log = out.with_suffix(".log")
            log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            prepare(out=out, genome=source, root=root)
        manifest, candidates = load_candidates(out, root=root)
        self.assertEqual(manifest["candidate_count"], 2)
        self.assertEqual(manifest["distinct_native_span_count"], 1)
        self.assertEqual({tuple(row["translation_tables"]) for row in candidates}, {(4,), (11,)})
        self.assertEqual(len({row["translation_sha256"] for row in candidates}), 2)
        self.assertEqual(len({row["DNA_sha256"] for row in candidates}), 1)

    def test_prepare_not_applicable_for_protein_input(self):
        tmp, root, annotation, screen, source = self._fixture(".faa")
        self.addCleanup(tmp.cleanup)
        out = root / "rescue"
        manifest = prepare(annotation, screen, out, root=root)
        self.assertEqual(manifest["status"], "not_applicable")
        self.assertIn("not_nucleotide", manifest["reason"])

    def test_missing_source_is_not_evaluated(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)
        source.unlink()
        out = root / "rescue"
        manifest = prepare(annotation, screen, out, root=root)
        self.assertEqual(manifest["status"], "not_evaluated")

    def test_source_mutation_during_scan_fails_closed(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def mutate(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            subject.write_text(subject.read_text() + ">tampered\nATG\n")
            return ["tblastn", "-db_gencode", str(table)], log

        with patch("dna_rescue._run_tblastn", side_effect=mutate), self.assertRaisesRegex(ValueError, "changed"):
            prepare(annotation, screen, root / "rescue", root=root)
        failed = json.loads((root / "rescue/dna_rescue_manifest.json").read_text())
        self.assertEqual(failed["status"], "failed")

    def test_raw_nucleotide_only_mode_has_no_protein_context(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        out = root / "raw_rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            manifest = prepare(out=out, genome=source, root=root)
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual(manifest["source_mode"], "raw_nucleotide_only")
        self.assertIsNone(manifest["annotation"])
        self.assertIsNone(manifest["screen"])
        checked, candidates = load_candidates(out, genome=source, root=root)
        self.assertTrue(candidates)
        self.assertTrue(all(c["identity_status"] == "unconfirmed" for c in candidates))

    def test_shared_verifier_rejects_bound_input_or_required_artifact_mutation(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            prepare(annotation, screen, out, root=root)
        source.write_text(source.read_text() + ">changed\nATG\n")
        with self.assertRaisesRegex(ValueError, "bound input"):
            load_candidates(out, annotation, screen, root=root)

        # Restore source and remove a required raw artifact; the verifier must
        # fail closed instead of silently returning the candidate JSON.
        source.write_text(">c\nATGAAAAAA\n>uncalled\nATGAAAAAA\n")
        (out / "tblastn_table11.log").unlink()
        with self.assertRaisesRegex(ValueError, "changed or missing"):
            load_candidates(out, annotation, screen, root=root)

    def test_shared_verifier_rejects_profile_input_mutation(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            prepare(out=out, genome=source, root=root)
        (root / "references/processed/catalytic_families.hmm").write_text("tampered\n")
        with self.assertRaisesRegex(ValueError, "profile scan input binding"):
            load_candidates(out, genome=source, root=root)

    def test_shared_verifier_rejects_profile_artifact_drop(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            prepare(out=out, genome=source, root=root)
        (out / "profile_scan/family.tbl").unlink()
        with self.assertRaisesRegex(ValueError, "profile scan output changed or missing"):
            load_candidates(out, genome=source, root=root)

    def test_raw_target_profile_evidence_uses_bound_hsp_length(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def fake_blast(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        def fake_targets(path, families):
            hsp = next(line[1:].split()[0] for line in (root / "rescue/dna_rescue_translations.faa").read_text().splitlines() if line.startswith(">"))
            return [{"protein_id": hsp, "family": "FeFe", "profile": "fake", "sequence_evalue": 1e-20,
                     "sequence_bitscore": 30.0, "sequence_bias": 0.0}]

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake_blast), patch("dna_rescue.load_targets", side_effect=fake_targets):
            prepare(out=out, genome=source, root=root)
        _manifest, candidates = load_candidates(out, genome=source, root=root)
        evidence = [e for c in candidates for e in c["profile_evidence"] if e["evidence_kind"] == "raw_target"]
        self.assertEqual(evidence[0]["query_length"], 3)
        self.assertEqual(evidence[0]["query_length_evidence"], "bound_HSP_translation")

    def test_domain_profile_length_mismatch_is_rejected(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)

        def fake_blast(root_arg, query, subject, table, out, threads, max_targets):
            out.write_text("r\tuncalled\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        def fake_domains(path, families):
            hsp = next(line[1:].split()[0] for line in (root / "rescue/dna_rescue_translations.faa").read_text().splitlines() if line.startswith(">"))
            return [{"protein_id": hsp, "family": "FeFe", "profile": "fake", "query_length": 4}]

        out = root / "rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake_blast), patch("dna_rescue.load_domains", side_effect=fake_domains):
            prepare(out=out, genome=source, root=root)
        with self.assertRaisesRegex(ValueError, "query length"):
            load_candidates(out, genome=source, root=root)

    def test_raw_genbank_creates_bound_subject_and_scans_dna(self):
        tmp, root, annotation, screen, source = self._fixture()
        self.addCleanup(tmp.cleanup)
        gb = root / "input.gb"
        record = SeqRecord(Seq("ATGAAAAAA"), id="gbc", name="gbc", description="test")
        record.annotations["molecule_type"] = "DNA"
        record.features.append(SeqFeature(FeatureLocation(0, 9), type="CDS", qualifiers={"transl_table": ["11"]}))
        SeqIO.write(record, gb, "genbank")

        def fake(root_arg, query, subject, table, out, threads, max_targets):
            self.assertEqual(subject.name, "nucleotide_subject.fna")
            out.write_text("r\tgbc\t100\t3\t1\t3\t1\t9\t1e-10\t20\t3\t9\tMKK\tMKK\t1\n")
            log = out.with_suffix(".log"); log.write_text("fake\n")
            return ["tblastn", "-db_gencode", str(table)], log

        out = root / "gb_rescue"
        with patch("dna_rescue._run_tblastn", side_effect=fake):
            manifest = prepare(out=out, genome=gb, root=root)
        self.assertEqual(manifest["source_format"], "genbank")
        self.assertTrue((out / "nucleotide_subject.fna").exists())
        load_candidates(out, genome=gb, root=root)


if __name__ == "__main__":
    unittest.main()
