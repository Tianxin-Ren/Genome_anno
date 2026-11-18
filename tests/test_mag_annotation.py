#!/usr/bin/env python3
"""
Unit and integration test suite for mag_annotation module scripts.
Verifies:
1. run_pyrodigal.py: GFF3, FAA, protein_metadata.tsv.gz (coordinate sorting & IDs)
2. dedup_proteins.py: Collision-checked dedup + run_metadata.json (dual statistical space)
3. shard_fasta.py: Residue-balanced sharding + shards_manifest.tsv (length stats)
4. merge_kofam_results.py: Primary artifact (protein_KO_all_strict) + Secondary (genome_KO_counts.long)
5. validate_invariance.py: Level 1, 2, 3 verification
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
import unittest
import pandas as pd
import numpy as np

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "modules", "mag_annotation", "scripts")
TEST_TMP = os.path.join(os.path.dirname(__file__), "..", "scratch", "test_mag_annotation_run")

class TestMagAnnotation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.makedirs(TEST_TMP, exist_ok=True)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(TEST_TMP):
            shutil.rmtree(TEST_TMP)

    def test_01_run_pyrodigal(self):
        """Test gene calling with coordinate-based sorting and metadata."""
        mock_fna = os.path.join(TEST_TMP, "mock_genome.fna")
        with open(mock_fna, "w") as f:
            # Synthetic DNA with two ORFs on contig1 and one on contig2
            f.write(">contig1\n")
            f.write("ATG" + "GCT"*35 + "TAA" + "NNNNNN" + "ATG" + "AAA"*40 + "TGA\n")
            f.write(">contig2\n")
            f.write("ATG" + "CCG"*30 + "TAG\n")

        out_prefix = os.path.join(TEST_TMP, "mock_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "run_pyrodigal.py"),
            "-i", mock_fna,
            "-g", "GENOME_001",
            "-o", out_prefix,
            "-m", "meta"
        ]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        faa_file = f"{out_prefix}.faa.gz"
        gff_file = f"{out_prefix}.gff.gz"
        meta_file = f"{out_prefix}_protein_metadata.tsv.gz"

        self.assertTrue(os.path.exists(faa_file))
        self.assertTrue(os.path.exists(gff_file))
        self.assertTrue(os.path.exists(meta_file))

        # Check metadata content
        with gzip.open(meta_file, "rt") as f:
            lines = [l.strip() for l in f if l.strip()]
            header = lines[0].split("\t")
            self.assertIn("gene_id", header)
            self.assertIn("protein_id", header)
            self.assertIn("translation_table", header)
            # Check deterministic gene_id format
            row1 = lines[1].split("\t")
            gene_id = row1[header.index("gene_id")]
            self.assertIn("GENOME_001___contig1___gene00001", gene_id)

        # Test auto gzip detection even when extension is .fa
        mock_fa_gz = os.path.join(TEST_TMP, "mock_compressed.fa")
        with gzip.open(mock_fa_gz, "wt") as gz_f:
            gz_f.write(">contig1\nATG" + "AAA"*35 + "TAA\n")
        out_prefix_gz = os.path.join(TEST_TMP, "mock_gz_out")
        cmd_gz = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "run_pyrodigal.py"),
            "-i", mock_fa_gz,
            "-g", "GENOME_002",
            "-o", out_prefix_gz
        ]
        subprocess.run(cmd_gz, check=True, capture_output=True, text=True)
        self.assertTrue(os.path.exists(f"{out_prefix_gz}.faa.gz"))

        # Test direct explicit destination flags (--out_faa, --out_gff, --out_meta)
        direct_faa = os.path.join(TEST_TMP, "direct", "proteins.faa.gz")
        direct_gff = os.path.join(TEST_TMP, "direct", "genes.gff.gz")
        direct_meta = os.path.join(TEST_TMP, "direct", "metadata.tsv.gz")
        cmd_direct = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "run_pyrodigal.py"),
            "-i", mock_fna,
            "-g", "GENOME_DIRECT",
            "--out_faa", direct_faa,
            "--out_gff", direct_gff,
            "--out_meta", direct_meta,
            "-m", "meta"
        ]
        subprocess.run(cmd_direct, check=True, capture_output=True, text=True)
        self.assertTrue(os.path.exists(direct_faa))
        self.assertTrue(os.path.exists(direct_gff))
        self.assertTrue(os.path.exists(direct_meta))

        # Test edge cases: tiny contigs (<60bp), degenerate chars, all-N contigs
        mock_dirty_fna = os.path.join(TEST_TMP, "mock_dirty.fna")
        with open(mock_dirty_fna, "w") as f:
            f.write(">c_tiny\nATGCA\n")
            f.write(">c_all_n\n" + "N"*80 + "\n")
            f.write(">c_weird_chars\nATG-*-*???123ATGCATGC\n")
            f.write(">c_normal\n" + "ATG" + "GCT"*35 + "TAA\n")
        dirty_out = os.path.join(TEST_TMP, "dirty_out")
        cmd_dirty = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "run_pyrodigal.py"),
            "-i", mock_dirty_fna,
            "-g", "GENOME_DIRTY",
            "-o", dirty_out,
            "-m", "meta"
        ]
        res_dirty = subprocess.run(cmd_dirty, capture_output=True, text=True)
        self.assertEqual(res_dirty.returncode, 0)
        self.assertTrue(os.path.exists(f"{dirty_out}.faa.gz"))
        self.assertTrue(os.path.exists(f"{dirty_out}.gff.gz"))
        self.assertTrue(os.path.exists(f"{dirty_out}_protein_metadata.tsv.gz"))

    def test_01c_auto_closed_genome_mode(self):
        """Auto mode selects single only for explicit closed-genome markers."""
        from modules.mag_annotation.scripts.run_pyrodigal import infer_auto_gene_calling_mode

        mode, reason = infer_auto_gene_calling_mode([
            ("u999856ctg_len-2079472_circular-yes_depth-1", "A" * 100),
        ])
        self.assertEqual(mode, "single")
        self.assertIn("explicit closed/complete markers", reason)

        mode, reason = infer_auto_gene_calling_mode([
            ("draft_contig_1", "A" * 100),
        ])
        self.assertEqual(mode, "meta")
        self.assertIn("lack closed/complete markers", reason)

        mode, _ = infer_auto_gene_calling_mode([
            ("circular-yes_contig_1", "A" * 100),
            ("unmarked_contig_2", "A" * 100),
        ])
        self.assertEqual(mode, "meta")

    def test_01b_pyrodigal_ipc_large_scale(self):
        """Test pyrodigal worker IPC with large payloads (>64KB) to verify no pipe deadlocks."""
        import types
        import pickle
        import multiprocessing as mp
        from modules.mag_annotation.scripts.run_pyrodigal import _pyrodigal_worker

        # Create mock pyrodigal module generating 2,500 genes (~300 KB pickled payload)
        mock_pyro = types.ModuleType("pyrodigal")
        mock_pyro.__version__ = "3.7.1"
        class MockGene:
            def __init__(self, i):
                self.begin = i * 300
                self.end = self.begin + 270
                self.strand = 1
                self.partial_begin = False
                self.partial_end = False
            def translate(self):
                return "M" + "A" * 89

        class MockFinder:
            def __init__(self, meta=True, **kwargs): pass
            def find_genes(self, b):
                return [MockGene(i) for i in range(2500)]

        mock_pyro.GeneFinder = MockFinder
        orig_pyro = sys.modules.get("pyrodigal")
        sys.modules["pyrodigal"] = mock_pyro

        try:
            # Test with fork context
            try:
                ctx = mp.get_context("fork")
            except ValueError:
                ctx = mp.get_context()

            parent_conn, child_conn = ctx.Pipe()
            clean_records = [("c1", "A" * 1000)]
            p = ctx.Process(target=_pyrodigal_worker, args=(clean_records, "meta", 11, "TEST_IPC", parent_conn, child_conn))
            p.start()
            child_conn.close()

            # Must poll within 5 seconds and not deadlock on 64KB OS pipe buffer
            has_data = parent_conn.poll(5)
            self.assertTrue(has_data, "Worker IPC timed out or deadlocked on pipe buffer!")

            status, path, ver = parent_conn.recv()
            p.join(timeout=5)
            self.assertEqual(status, "OK")
            self.assertEqual(ver, "3.7.1")
            self.assertTrue(os.path.exists(path))

            with open(path, "rb") as f:
                genes = pickle.load(f)
            os.remove(path)
            parent_conn.close()

            self.assertEqual(len(genes), 2500)
            self.assertEqual(genes[0]["translation"], "M" + "A" * 89)
        finally:
            if orig_pyro is not None:
                sys.modules["pyrodigal"] = orig_pyro
            elif "pyrodigal" in sys.modules:
                del sys.modules["pyrodigal"]

    def test_02_dedup_proteins(self):
        """Test 100% exact deduplication and dual-space provenance."""
        prot1 = os.path.join(TEST_TMP, "p1.faa")
        prot2 = os.path.join(TEST_TMP, "p2.faa")
        
        # Seq A appears in both p1 and p2 (duplicate)
        # Seq B appears only in p1
        with open(prot1, "w") as f:
            f.write(">GENOME_001___c1___prot00001\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE\n")
            f.write(">GENOME_001___c1___prot00002\nMKTIIALSYIFCLVFA\n")
        with open(prot2, "w") as f:
            f.write(">GENOME_002___c1___prot00001\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE*\n") # identical to seq A, with stop

        list_file = os.path.join(TEST_TMP, "protein_files.txt")
        with open(list_file, "w") as f:
            f.write(f"{prot1}\n{prot2}\n")

        dedup_dir = os.path.join(TEST_TMP, "dedup_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "dedup_proteins.py"),
            "-i", list_file,
            "-o", dedup_dir,
            "-c", "TestProteome"
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        meta_json = os.path.join(dedup_dir, "run_metadata.json")
        self.assertTrue(os.path.exists(meta_json))
        with open(meta_json) as jf:
            meta = json.load(jf)
            self.assertEqual(meta["raw_sequences_count"], 3)
            self.assertEqual(meta["unique_sequences_count"], 2)
            self.assertEqual(meta["hmmsearch_database_size_Z"], 2)
            self.assertAlmostEqual(meta["dedup_ratio"], 1/3, places=3)

    def test_03_shard_fasta(self):
        """Test residue-balanced sharding with length stats."""
        test_faa = os.path.join(TEST_TMP, "multi.faa")
        with open(test_faa, "w") as f:
            for i in range(20):
                f.write(f">seq_{i+1}\n" + "M" * 100 + "\n")

        shards_dir = os.path.join(TEST_TMP, "shards_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "shard_fasta.py"),
            "-i", test_faa,
            "-o", shards_dir,
            "--target_aa", "500", # 500 aa per shard -> 5 sequences per shard -> 4 shards
            "--target_proteins", "100"
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        manifest = os.path.join(shards_dir, "shards_manifest.tsv")
        self.assertTrue(os.path.exists(manifest))
        with open(manifest) as f:
            lines = f.readlines()
            # header + 4 shards
            self.assertEqual(len(lines), 5)
            header = lines[0].strip().split("\t")
            self.assertIn("max_length", header)
            self.assertIn("median_length", header)

    def test_04_merge_results(self):
        """Test merge_kofam_results producing primary asset and secondary sparse matrix with 80% threshold."""
        # Create mock shard tsv
        res_dir = os.path.join(TEST_TMP, "pykofam_shards")
        os.makedirs(res_dir, exist_ok=True)
        shard1_tsv = os.path.join(res_dir, "shard_0001.tsv")
        with open(shard1_tsv, "w") as f:
            f.write("query\ttarget\tscore\tthreshold\tscore_type\tevalue\tdescription\n")
            # UPROT_000000001 matches:
            # 1. K02111: strict pass (score 450.5 >= threshold 300.0) -> is_strict=True, rank 1, pass_type=score_pass
            f.write("UPROT_000000001\tK02111\t450.5\t300.0\tfull\t1e-90\tATP synthase subunit alpha\n")
            # 2. K02112: strict pass (score 320.0 >= threshold 250.0) -> is_strict=True, rank 2, pass_type=score_pass
            f.write("UPROT_000000001\tK02112\t320.0\t250.0\tfull\t1e-60\tATP synthase subunit beta\n")
            # 3. K02113: relaxed pass (score 255.0 >= 0.8 * 300.0 = 240.0) -> is_strict=False, rank 3, pass_type=score_pass
            f.write("UPROT_000000001\tK02113\t255.0\t300.0\tfull\t1e-45\tATP synthase subunit gamma\n")
            # 4. K02114: rejected (score 195.0 < 0.8 * 300.0 = 240.0) -> filtered out because UPROT_000000001 has score-passing hits!
            f.write("UPROT_000000001\tK02114\t195.0\t300.0\tfull\t1e-25\tATP synthase subunit delta\n")

            # UPROT_000000002 has NO score-passing hits -> Tier 2 E-value fallback:
            # 1. K00001: sub-threshold (score 180.0 < 0.8 * 300.0 = 240.0), evalue 1e-45 (BEST E-value!)
            f.write("UPROT_000000002\tK00001\t180.0\t300.0\tfull\t1e-45\talcohol dehydrogenase\n")
            # 2. K00002: sub-threshold (score 120.0 < 0.8 * 200.0 = 160.0), evalue 1e-20
            f.write("UPROT_000000002\tK00002\t120.0\t200.0\tfull\t1e-20\talcohol dehydrogenase (NADP+)\n")
            # 3. K00003: rejected by fallback cutoff (evalue 0.05 > 1e-3)
            f.write("UPROT_000000002\tK00003\t30.0\t200.0\tfull\t0.05\thomoserine dehydrogenase\n")

        # Create mock mapping
        mapping_file = os.path.join(TEST_TMP, "mock_mapping.tsv.gz")
        with gzip.open(mapping_file, "wt") as f:
            f.write("protein_id\tunique_id\taa_length\tgenome_id\n")
            f.write("GCA_001___c1___prot00001\tUPROT_000000001\t500\tGCA_001\n")
            f.write("GCA_001___c1___prot00002\tUPROT_000000002\t400\tGCA_001\n")

        merged_dir = os.path.join(TEST_TMP, "merged_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "merge_kofam_results.py"),
            "-s", res_dir,
            "-m", mapping_file,
            "-o", merged_dir,
            "--threshold_factor", "0.80",
            "--evalue_fallback_cutoff", "1e-3",
            "--export_dense"
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        strict_out = os.path.join(merged_dir, "protein_KO_all_strict.tsv")
        best_out = os.path.join(merged_dir, "protein_KO_best.tsv")
        long_out = os.path.join(merged_dir, "genome_KO_counts.long.tsv")

        self.assertTrue(os.path.exists(strict_out))
        self.assertTrue(os.path.exists(best_out))
        self.assertTrue(os.path.exists(long_out))

        # Check primary artifact multi-hit retention, fallback, and pass_type annotations
        with open(strict_out, "r", encoding="utf-8") as f:
            lines = [l.strip() for l in f if l.strip()]
            header = lines[0].split("\t")
            self.assertIn("is_strict", header)
            self.assertIn("score_margin", header)
            self.assertIn("pass_type", header)
            # header + 3 passing hits for prot1 (K02111, K02112, K02113) + 1 fallback hit for prot2 (K00001) = 5 lines
            self.assertEqual(len(lines), 5)

            # Hit 1 (prot1 best score)
            row1 = lines[1].split("\t")
            self.assertEqual(row1[0], "GCA_001___c1___prot00001")
            self.assertEqual(row1[3], "K02111")
            self.assertEqual(row1[9], "1")            # ko_rank
            self.assertEqual(row1[11], "True")        # is_strict
            self.assertEqual(row1[12], "score_pass")  # pass_type

            # Hit 2 (prot1 second strict)
            row2 = lines[2].split("\t")
            self.assertEqual(row2[3], "K02112")
            self.assertEqual(row2[9], "2")            # ko_rank
            self.assertEqual(row2[11], "True")        # is_strict
            self.assertEqual(row2[12], "score_pass")

            # Hit 3 (prot1 relaxed 85% hit, score 255 vs threshold 300 -> margin -45.0)
            row3 = lines[3].split("\t")
            self.assertEqual(row3[3], "K02113")
            self.assertEqual(row3[9], "3")            # ko_rank
            self.assertEqual(row3[10], "-45.0")       # score_margin
            self.assertEqual(row3[11], "False")       # is_strict = False (relaxed)
            self.assertEqual(row3[12], "score_pass")

            # Hit 4 (prot2 rescued via Tier 2 E-value fallback: K00001 with min evalue 1e-45)
            row4 = lines[4].split("\t")
            self.assertEqual(row4[0], "GCA_001___c1___prot00002")
            self.assertEqual(row4[3], "K00001")
            self.assertEqual(row4[9], "1")                 # ko_rank
            self.assertEqual(row4[10], "-120.0")           # score_margin = 180 - 300
            self.assertEqual(row4[11], "False")            # is_strict
            self.assertEqual(row4[12], "evalue_fallback")  # pass_type

        # Check derivative best hit
        with open(best_out, "r", encoding="utf-8") as f:
            lines = [l.strip() for l in f if l.strip()]
            self.assertEqual(len(lines), 3) # header + 1 best hit for prot1 + 1 fallback for prot2
            row_b1 = lines[1].split("\t")
            self.assertEqual(row_b1[0], "GCA_001___c1___prot00001")
            self.assertEqual(row_b1[3], "K02111")
            self.assertEqual(row_b1[10], "True")
            self.assertEqual(row_b1[11], "score_pass")

            row_b2 = lines[2].split("\t")
            self.assertEqual(row_b2[0], "GCA_001___c1___prot00002")
            self.assertEqual(row_b2[3], "K00001")
            self.assertEqual(row_b2[10], "False")
            self.assertEqual(row_b2[11], "evalue_fallback")

    def test_04b_merge_native_pykofamsearch_format(self):
        """Test merge_kofam_results with real pykofamsearch native output header and columns."""
        res_dir = os.path.join(TEST_TMP, "pykofam_native_shards")
        os.makedirs(res_dir, exist_ok=True)
        shard1_tsv = os.path.join(res_dir, "shard_0001.tsv")
        with open(shard1_tsv, "w") as f:
            f.write("id_protein\tid_ko\tthreshold\tscore\te-value\tdefinition\tenzyme_commission\n")
            f.write("UPROT_000162789\tK00001\t341.1\t333.992\t4.20643e-98\talcohol dehydrogenase [EC:1.1.1.1]\t{'1.1.1.1'}\n")

        mapping_file = os.path.join(TEST_TMP, "mock_mapping_native.tsv.gz")
        with gzip.open(mapping_file, "wt") as f:
            f.write("protein_id\tunique_id\taa_length\tgenome_id\n")
            f.write("GCA_001___c1___prot00001\tUPROT_000162789\t400\tGCA_001\n")

        merged_dir = os.path.join(TEST_TMP, "merged_native_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "merge_kofam_results.py"),
            "-s", res_dir,
            "-m", mapping_file,
            "-o", merged_dir,
            "--threshold_factor", "0.80"
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        strict_out = os.path.join(merged_dir, "protein_KO_all_strict.tsv")
        self.assertTrue(os.path.exists(strict_out))
        with open(strict_out, "r", encoding="utf-8") as f:
            lines = [l.strip() for l in f if l.strip()]
            self.assertEqual(len(lines), 2)
            row = lines[1].split("\t")
            # protein_id, gene_id, genome_id, ko, ko_definition, bitscore, curated_threshold, score_type, evalue, ko_rank, score_margin, is_strict, pass_type
            self.assertEqual(row[3], "K00001")
            self.assertEqual(row[4], "alcohol dehydrogenase [EC:1.1.1.1]")
            self.assertEqual(float(row[5]), 333.992)
            self.assertEqual(float(row[6]), 341.1)
            self.assertEqual(row[7], "full")
            self.assertAlmostEqual(float(row[8]), 4.20643e-98)
            self.assertEqual(row[11], "False") # 333.992 < 341.1
            self.assertEqual(row[12], "score_pass")



    def test_05_validate_invariance(self):
        """Test validate_invariance script across Level 1, 2, 3."""
        out_report = os.path.join(TEST_TMP, "val_report.json")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "validate_invariance.py"),
            "-o", out_report
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        self.assertTrue(os.path.exists(out_report))
        with open(out_report) as jf:
            rep = json.load(jf)
    def test_06_sample_resolution_and_defensive_formatting(self):
        """Test robust multi-sample resolution and defensive first_prot logic."""
        def resolve_active(cfg):
            cli_val = cfg.get("sample") or cfg.get("sample_name")
            if isinstance(cli_val, str) and cli_val.strip() and cli_val.strip().lower() not in ["all", "*"]:
                return [cli_val.strip()]
            for s_candidate in [cfg.get("samples"), cfg.get("sample")]:
                if isinstance(s_candidate, dict) and s_candidate:
                    return list(s_candidate.keys())
                elif isinstance(s_candidate, list) and s_candidate:
                    return [str(x).strip() for x in s_candidate if str(x).strip()]
                elif isinstance(s_candidate, str) and s_candidate.strip() and s_candidate.strip().lower() not in ["all", "*"]:
                    return [s_candidate.strip()]
            return ["default"]

        # Test case A: multiple samples configured, no single sample filter -> returns all samples
        config_multi = {
            "samples": {
                "hadar_Ytt": {"genome_dir": "/path/to/hadar"},
                "gtdb_r232_archaea_all_representatives": {"genome_dir": "/path/to/gtdb"}
            }
        }
        self.assertEqual(resolve_active(config_multi), ["hadar_Ytt", "gtdb_r232_archaea_all_representatives"])

        # Test case B: CLI override with sample=gtdb_r232_archaea_all_representatives
        config_multi_override = {
            "sample": "gtdb_r232_archaea_all_representatives",
            "samples": {
                "hadar_Ytt": {"genome_dir": "/path/to/hadar"},
                "gtdb_r232_archaea_all_representatives": {"genome_dir": "/path/to/gtdb"}
            }
        }
        self.assertEqual(resolve_active(config_multi_override), ["gtdb_r232_archaea_all_representatives"])

        # Test case C: defensive first_prot lambda with empty list
        prots_empty = []
        first_prot_fn = lambda prots: str(prots[0]) if prots else ""
        self.assertEqual(first_prot_fn(prots_empty), "")

        # Test case D: defensive first_prot lambda with non-empty list
        prots_mock = ["/path/to/prot.faa"]
        self.assertEqual(first_prot_fn(prots_mock), "/path/to/prot.faa")

    def test_07_ko_score_type_resolution(self):
        """Test accurate score_type resolution ('domain' vs 'full') using ko_list and database.pkl."""
        # Create mock ko_list TSV
        ko_list_file = os.path.join(TEST_TMP, "mock_ko_list.tsv")
        with open(ko_list_file, "w") as f:
            f.write("koid\tthreshold\tscore_type\tprofile_type\tF-measure\tn_seq\tn_seq_o\tscore_margin\tdefinition\n")
            f.write("K00001\t341.1\tdomain\tEC 1.1.1.1\t0.85\t100\t100\t0.0\talcohol dehydrogenase\n")
            f.write("K00002\t200.0\tfull\tEC 1.1.1.2\t0.90\t80\t80\t0.0\talcohol dehydrogenase (NADP+)\n")

        # Create mock shard output with native pykofamsearch header (no score_type column)
        res_dir = os.path.join(TEST_TMP, "pykofam_score_type_shards")
        os.makedirs(res_dir, exist_ok=True)
        shard_tsv = os.path.join(res_dir, "shard_0001.tsv")
        with open(shard_tsv, "w") as f:
            f.write("id_protein\tid_ko\tthreshold\tscore\te-value\tdefinition\tenzyme_commission\n")
            # Hit 1: K00001 (curated score_type is domain)
            f.write("UPROT_001\tK00001\t341.1\t350.0\t1e-100\talcohol dehydrogenase\t{'1.1.1.1'}\n")
            # Hit 2: K00002 (curated score_type is full)
            f.write("UPROT_002\tK00002\t200.0\t210.0\t1e-60\talcohol dehydrogenase (NADP+)\t{'1.1.1.2'}\n")
            # Hit 3: K99999 (not in ko_list -> should default to full)
            f.write("UPROT_003\tK99999\t150.0\t160.0\t1e-40\thypothetical protein\t{}\n")

        mapping_file = os.path.join(TEST_TMP, "mock_mapping_score_type.tsv.gz")
        with gzip.open(mapping_file, "wt") as f:
            f.write("protein_id\tunique_id\taa_length\tgenome_id\n")
            f.write("GCA_001___c1___prot00001\tUPROT_001\t400\tGCA_001\n")
            f.write("GCA_001___c1___prot00002\tUPROT_002\t400\tGCA_001\n")
            f.write("GCA_001___c1___prot00003\tUPROT_003\t400\tGCA_001\n")

        merged_dir = os.path.join(TEST_TMP, "merged_score_type_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "merge_kofam_results.py"),
            "-s", res_dir,
            "-m", mapping_file,
            "-o", merged_dir,
            "--ko_list", ko_list_file
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        strict_out = os.path.join(merged_dir, "protein_KO_all_strict.tsv")
        best_out = os.path.join(merged_dir, "protein_KO_best.tsv")
        self.assertTrue(os.path.exists(strict_out))
        self.assertTrue(os.path.exists(best_out))

        with open(strict_out, "r", encoding="utf-8") as f:
            lines = [l.strip() for l in f if l.strip()]
            header = lines[0].split("\t")
            st_idx = header.index("score_type")
            ko_idx = header.index("ko")

            rows = {line.split("\t")[ko_idx]: line.split("\t")[st_idx] for line in lines[1:]}
            self.assertEqual(rows.get("K00001"), "domain")
            self.assertEqual(rows.get("K00002"), "full")
            self.assertEqual(rows.get("K99999"), "full")

        with open(best_out, "r", encoding="utf-8") as f:
            lines = [l.strip() for l in f if l.strip()]
            header = lines[0].split("\t")
            st_idx = header.index("score_type")
            ko_idx = header.index("ko")

            rows = {line.split("\t")[ko_idx]: line.split("\t")[st_idx] for line in lines[1:]}
            self.assertEqual(rows.get("K00001"), "domain")
            self.assertEqual(rows.get("K00002"), "full")
            self.assertEqual(rows.get("K99999"), "full")

    def test_08_combine_mag_annotations(self):
        """Test combine_mag_annotations.py for merging two non-overlapping batches."""
        batch1_dir = os.path.join(TEST_TMP, "batch1")
        batch2_dir = os.path.join(TEST_TMP, "batch2")
        combined_dir = os.path.join(TEST_TMP, "batch_combined")

        # Create mock batch 1
        os.makedirs(os.path.join(batch1_dir, "merged"), exist_ok=True)
        os.makedirs(os.path.join(batch1_dir, "proteins", "all"), exist_ok=True)
        os.makedirs(os.path.join(batch1_dir, "proteins", "unique"), exist_ok=True)
        os.makedirs(os.path.join(batch1_dir, "genes"), exist_ok=True)
        os.makedirs(os.path.join(batch1_dir, "metadata"), exist_ok=True)

        with open(os.path.join(batch1_dir, "merged", "protein_KO_all_strict.tsv"), "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tko_rank\tscore_margin\tis_strict\tpass_type\n")
            f.write("GCA_001___c1___prot1\tGCA_001___c1___gene1\tGCA_001\tK00001\tadh\t350.0\t340.0\tdomain\t1e-90\t1\t10.0\tTrue\tscore_pass\n")
        with open(os.path.join(batch1_dir, "merged", "protein_KO_best.tsv"), "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tscore_margin\tis_strict\tpass_type\n")
            f.write("GCA_001___c1___prot1\tGCA_001___c1___gene1\tGCA_001\tK00001\tadh\t350.0\t340.0\tdomain\t1e-90\t10.0\tTrue\tscore_pass\n")
        with open(os.path.join(batch1_dir, "merged", "genome_KO_counts.long.tsv"), "w") as f:
            f.write("genome_id\tko\tcount\n")
            f.write("GCA_001\tK00001\t1\n")
        with gzip.open(os.path.join(batch1_dir, "proteins", "all", "GCA_001.faa.gz"), "wt") as f:
            f.write(">GCA_001___c1___prot1\nMKWVTFISLLLLFSSAYSRG\n")

        # Create mock batch 2 (non-overlapping)
        os.makedirs(os.path.join(batch2_dir, "merged"), exist_ok=True)
        os.makedirs(os.path.join(batch2_dir, "proteins", "all"), exist_ok=True)
        os.makedirs(os.path.join(batch2_dir, "proteins", "unique"), exist_ok=True)
        os.makedirs(os.path.join(batch2_dir, "genes"), exist_ok=True)
        os.makedirs(os.path.join(batch2_dir, "metadata"), exist_ok=True)

        with open(os.path.join(batch2_dir, "merged", "protein_KO_all_strict.tsv"), "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tko_rank\tscore_margin\tis_strict\tpass_type\n")
            f.write("GCA_002___c1___prot1\tGCA_002___c1___gene1\tGCA_002\tK00002\tadh2\t250.0\t200.0\tfull\t1e-60\t1\t50.0\tTrue\tscore_pass\n")
        with open(os.path.join(batch2_dir, "merged", "protein_KO_best.tsv"), "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tscore_margin\tis_strict\tpass_type\n")
            f.write("GCA_002___c1___prot1\tGCA_002___c1___gene1\tGCA_002\tK00002\tadh2\t250.0\t200.0\tfull\t1e-60\t50.0\tTrue\tscore_pass\n")
        with open(os.path.join(batch2_dir, "merged", "genome_KO_counts.long.tsv"), "w") as f:
            f.write("genome_id\tko\tcount\n")
            f.write("GCA_002\tK00002\t1\n")
        with gzip.open(os.path.join(batch2_dir, "proteins", "all", "GCA_002.faa.gz"), "wt") as f:
            f.write(">GCA_002___c1___prot1\nMKTIIALSYIFCLVFA\n")

        # Run combine script
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "combine_mag_annotations.py"),
            "-i", batch1_dir, batch2_dir,
            "-o", combined_dir
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        comb_strict = os.path.join(combined_dir, "merged", "protein_KO_all_strict.tsv")
        comb_best = os.path.join(combined_dir, "merged", "protein_KO_best.tsv")
        comb_long = os.path.join(combined_dir, "merged", "genome_KO_counts.long.tsv")
        self.assertTrue(os.path.exists(comb_strict))
        self.assertTrue(os.path.exists(comb_best))
        self.assertTrue(os.path.exists(comb_long))

        # Check combined content
        with open(comb_strict) as f:
            lines = [l.strip() for l in f if l.strip()]
            self.assertEqual(len(lines), 3) # header + GCA_001 + GCA_002
            self.assertIn("GCA_001", lines[1])
            self.assertIn("GCA_002", lines[2])

        with open(comb_long) as f:
            lines = [l.strip() for l in f if l.strip()]
            self.assertEqual(len(lines), 3) # header + GCA_001 + GCA_002

        # Check symlinks
        symlink_1 = os.path.join(combined_dir, "proteins", "all", "GCA_001.faa.gz")
        symlink_2 = os.path.join(combined_dir, "proteins", "all", "GCA_002.faa.gz")
        self.assertTrue(os.path.islink(symlink_1))
        self.assertTrue(os.path.islink(symlink_2))

    def test_09_postprocess_atp_synthase(self):
        """Test postprocess_atp_synthase.py across all biological edge cases and gating rules."""
        test_dir = os.path.join(TEST_TMP, "postprocess_test")
        os.makedirs(test_dir, exist_ok=True)

        meta_path = os.path.join(test_dir, "test_metadata.tsv")
        with open(meta_path, "w") as f:
            f.write("accession\tcheckm2_completeness\tcheckm2_contamination\tcheckm_completeness\tcheckm_contamination\tcontig_count\tgc_percentage\tgenome_size\tgtdb_representative\tgtdb_taxonomy\tgtdb_type_designation_ncbi_taxa\tgtdb_type_designation_ncbi_taxa_sources\tgtdb_type_species_of_genus\tl50_contigs\tmimag_high_quality\tmimag_low_quality\tmimag_medium_quality\tn50_contigs\tncbi_assembly_level\tncbi_genome_category\tncbi_isolation_source\tncbi_trna_count\ttrna_aa_count\ttrna_count\n")
            # 1. G_VA_STRICT: High quality, complete V/A
            f.write("GB_GCA_001.1\t98.0\t1.0\t98.0\t1.0\t10\t45.0\t2000000\tt\td__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei\t.\t.\t.\t2\tt\tf\tf\t500000\tComplete\tderived\thydrothermal\t45\t20\t45\n")
            # 2. G_VA_RESCUED: High quality, 8 strong + 1 operon-rescued borderline
            f.write("RS_GCF_002.1\t96.0\t2.0\t96.0\t2.0\t12\t42.0\t1800000\tt\td__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina barkeri\t.\t.\t.\t3\tt\tf\tf\t400000\tComplete\tderived\tsediment\t42\t20\t42\n")
            # 3. G_VA_NEAR: High quality, 8 strong + 1 unrescued borderline on diff contig
            f.write("GCA_003.1\t95.0\t1.5\t95.0\t1.5\t15\t40.0\t1900000\tf\td__Archaea;p__Halobacteriota;c__Methanococci;o__Methanococcales;f__Methanocaldococcaceae;g__Methanocaldococcus;s__Methanocaldococcus jannaschii\t.\t.\t.\t4\tt\tf\tf\t350000\tComplete\tderived\thot spring\t40\t20\t40\n")
            # 4. G_MIXED: Complete V/A + Partial F-like (Must be Mixed_evidence_ambiguous!)
            f.write("GCA_004.1\t97.0\t1.0\t97.0\t1.0\t8\t50.0\t2500000\tt\td__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanothrix;s__Methanothrix soehngenii\t.\t.\t.\t2\tt\tf\tf\t600000\tComplete\tderived\tanarobic sludge\t48\t20\t48\n")
            # 5. G_DUAL: Complete V/A + Complete F-like
            f.write("GCA_005.1\t99.0\t0.5\t99.0\t0.5\t5\t48.0\t2800000\tt\td__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina acetivorans\t.\t.\t.\t1\tt\tf\tf\t800000\tComplete\tderived\tmarine\t50\t20\t50\n")
            # 6. G_F_ONLY: Complete F-like + 0 V/A
            f.write("GCA_006.1\t95.0\t1.0\t95.0\t1.0\t10\t55.0\t2200000\tt\td__Archaea;p__Thermoproteota;c__Thermoprotei;o__Sulfolobales;f__Sulfolobaceae;g__Sulfolobus;s__Sulfolobus islandicus\t.\t.\t.\t3\tt\tf\tf\t450000\tComplete\tderived\tsolfatara\t44\t20\t44\n")
            # 7. G_PARTIAL_MAG: 5/9 units of V/A (Must be Ambiguous, NOT ATP-!)
            f.write("GCA_007.1\t85.0\t3.0\t85.0\t3.0\t40\t38.0\t1500000\tf\td__Archaea;p__Thermoproteota;c__Thermoprotei;o__Sulfolobales;f__Sulfolobaceae;g__Sulfolobus;s__Sulfolobus acidocaldarius\t.\t.\t.\t10\tf\tt\tf\t150000\tScaffold\tderived\tsoil\t35\t18\t35\n")
            # 8. G_HQ_NONE: 0 ATPase units, HQ genome (Valid ATP-)
            f.write("GCA_008.1\t95.0\t0.5\t95.0\t0.5\t6\t32.0\t800000\tt\td__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Nanoarchaeales;f__Nanoarchaeaceae;g__Nanoarchaeum;s__Nanoarchaeum equitans\t.\t.\t.\t1\tt\tf\tf\t500000\tComplete\tderived\tsymbiont\t30\t15\t30\n")
            # 9. G_LQ_NONE: 0 ATPase units, low quality genome (Must be Ambiguous, NOT ATP-!)
            f.write("GCA_009.1\t55.0\t2.0\t55.0\t2.0\t80\t30.0\t600000\tf\td__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Nanoarchaeales;f__Nanoarchaeaceae;g__Nanoarchaeum;s__Nanoarchaeum sp.\t.\t.\t.\t25\tf\tf\tt\t40000\tContig\tderived\tsediment\t20\t10\t20\n")

        # Create mock protein_KO_all_strict.tsv
        ko_strict_path = os.path.join(test_dir, "test_protein_KO_all_strict.tsv")
        with open(ko_strict_path, "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tko_rank\tscore_margin\tis_strict\tpass_type\n")
            # G_VA_STRICT: All 9 V/A units strong (> threshold 200)
            va_kos = ["K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02107", "K02123", "K02124"]
            for idx, ko in enumerate(va_kos):
                order = idx + 1
                f.write(f"GB_GCA_001.1___c1___prot{order:05d}\tGB_GCA_001.1___c1___gene{order:05d}\tGB_GCA_001.1\t{ko}\tV_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")

            # G_VA_RESCUED: 8 strong on c1, K02124 borderline (190.0/200.0, m=-0.05) also on c1 order 9 (operon rescued)
            for idx, ko in enumerate(va_kos[:-1]):
                order = idx + 1
                f.write(f"RS_GCF_002.1___c1___prot{order:05d}\tRS_GCF_002.1___c1___gene{order:05d}\tRS_GCF_002.1\t{ko}\tV_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")
            f.write("RS_GCF_002.1___c1___prot00009\tRS_GCF_002.1___c1___gene00009\tRS_GCF_002.1\tK02124\tV_ATPase_K\t190.0\t200.0\tfull\t1e-40\t1\t-10.0\tFalse\tscore_pass\n")

            # G_VA_NEAR: 8 strong on c1, K02124 borderline on c99 (unrescued borderline)
            for idx, ko in enumerate(va_kos[:-1]):
                order = idx + 1
                f.write(f"GCA_003.1___c1___prot{order:05d}\tGCA_003.1___c1___gene{order:05d}\tGCA_003.1\t{ko}\tV_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")
            f.write("GCA_003.1___c99___prot00001\tGCA_003.1___c99___gene00001\tGCA_003.1\tK02124\tV_ATPase_K\t190.0\t200.0\tfull\t1e-40\t1\t-10.0\tFalse\tscore_pass\n")

            # G_MIXED: All 9 V/A units + 5 F-like units (K02111, K02112, K02113, K02114, K02115)
            for idx, ko in enumerate(va_kos):
                order = idx + 1
                f.write(f"GCA_004.1___c1___prot{order:05d}\tGCA_004.1___c1___gene{order:05d}\tGCA_004.1\t{ko}\tV_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")
            f_kos = ["K02111", "K02112", "K02113", "K02114", "K02115", "K02108", "K02109", "K02110"]
            for idx, ko in enumerate(f_kos[:5]):
                order = idx + 1
                f.write(f"GCA_004.1___c2___prot{order:05d}\tGCA_004.1___c2___gene{order:05d}\tGCA_004.1\t{ko}\tF_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")

            # G_DUAL: All 9 V/A units + All 8 F-like units
            for idx, ko in enumerate(va_kos):
                order = idx + 1
                f.write(f"GCA_005.1___c1___prot{order:05d}\tGCA_005.1___c1___gene{order:05d}\tGCA_005.1\t{ko}\tV_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")
            for idx, ko in enumerate(f_kos):
                order = idx + 1
                f.write(f"GCA_005.1___c2___prot{order:05d}\tGCA_005.1___c2___gene{order:05d}\tGCA_005.1\t{ko}\tF_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")

            # G_F_ONLY: All 8 F-like units, 0 V/A
            for idx, ko in enumerate(f_kos):
                order = idx + 1
                f.write(f"GCA_006.1___c1___prot{order:05d}\tGCA_006.1___c1___gene{order:05d}\tGCA_006.1\t{ko}\tF_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")

            # G_PARTIAL_MAG: 5/9 V/A units on c1 (order 1 to 5)
            for idx, ko in enumerate(va_kos[:5]):
                order = idx + 1
                f.write(f"GCA_007.1___c1___prot{order:05d}\tGCA_007.1___c1___gene{order:05d}\tGCA_007.1\t{ko}\tV_ATPase_{ko}\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")

        # Run postprocess script
        out_dir = os.path.join(test_dir, "out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "postprocess_atp_synthase.py"),
            "-k", ko_strict_path,
            "-m", meta_path,
            "-o", out_dir,
            "--include_all_metadata_genomes"
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        master_tsv = os.path.join(out_dir, "genome_atp_synthase_master_table.tsv")
        long_tsv = os.path.join(out_dir, "genome_atp_synthase_modules.long.tsv")
        subunits_tsv = os.path.join(out_dir, "genome_atp_synthase_subunits.tsv")
        pruned_meta_tsv = os.path.join(out_dir, "pruned_metadata.tsv")

        self.assertTrue(os.path.exists(master_tsv))
        self.assertTrue(os.path.exists(long_tsv))
        self.assertTrue(os.path.exists(subunits_tsv))
        self.assertTrue(os.path.exists(pruned_meta_tsv))

        # Validate master table contents
        df_master = pd.read_csv(master_tsv, sep="\t")
        df_master.set_index("accession", inplace=True)

        # 1. G_VA_STRICT: Complete_strict, VA_only, ATP_POSITIVE, VA_ONLY
        self.assertEqual(df_master.loc["GB_GCA_001.1", "M00159_state"], "Complete_strict")
        self.assertEqual(df_master.loc["GB_GCA_001.1", "atp_synthase_genome_phenotype"], "VA_only")
        self.assertEqual(df_master.loc["GB_GCA_001.1", "atp_contrast1_cohort"], "ATP_POSITIVE")
        self.assertEqual(df_master.loc["GB_GCA_001.1", "atp_contrast2_cohort"], "VA_ONLY")

        # 2. G_VA_RESCUED: Complete_supported (operon context rescue), VA_only
        self.assertEqual(df_master.loc["RS_GCF_002.1", "M00159_state"], "Complete_supported")
        self.assertEqual(df_master.loc["RS_GCF_002.1", "atp_synthase_genome_phenotype"], "VA_only")

        # 3. G_VA_NEAR: Near_complete, VA_only
        self.assertEqual(df_master.loc["GCA_003.1", "M00159_state"], "Near_complete")
        self.assertEqual(df_master.loc["GCA_003.1", "atp_synthase_genome_phenotype"], "VA_only")

        # 4. G_MIXED: Mixed_evidence_ambiguous, EXCLUDED from Contrast 2
        self.assertEqual(df_master.loc["GCA_004.1", "atp_synthase_genome_phenotype"], "Mixed_evidence_ambiguous")
        self.assertEqual(df_master.loc["GCA_004.1", "atp_contrast2_cohort"], "EXCLUDED")

        # 5. G_DUAL: VA_plus_F_like, EXCLUDED from Contrast 2
        self.assertEqual(df_master.loc["GCA_005.1", "atp_synthase_genome_phenotype"], "VA_plus_F_like")
        self.assertEqual(df_master.loc["GCA_005.1", "atp_contrast2_cohort"], "EXCLUDED")

        # 6. G_F_ONLY: F_like_only, ATP_POSITIVE, F_LIKE_ONLY
        self.assertEqual(df_master.loc["GCA_006.1", "atp_synthase_genome_phenotype"], "F_like_only")
        self.assertEqual(df_master.loc["GCA_006.1", "atp_contrast1_cohort"], "ATP_POSITIVE")
        self.assertEqual(df_master.loc["GCA_006.1", "atp_contrast2_cohort"], "F_LIKE_ONLY")

        # 7. G_PARTIAL_MAG: Partial, Ambiguous, EXCLUDED from ATP_NEGATIVE!
        self.assertIn("Partial", df_master.loc["GCA_007.1", "M00159_state"])
        self.assertEqual(df_master.loc["GCA_007.1", "atp_contrast1_cohort"], "EXCLUDED")

        # 8. G_HQ_NONE: None_confident, ATP_NEGATIVE
        self.assertEqual(df_master.loc["GCA_008.1", "atp_synthase_genome_phenotype"], "None_confident")
        self.assertEqual(df_master.loc["GCA_008.1", "atp_contrast1_cohort"], "ATP_NEGATIVE")

        # 9. G_LQ_NONE: Low_quality_not_detected, Ambiguous, EXCLUDED from ATP_NEGATIVE!
        self.assertEqual(df_master.loc["GCA_009.1", "M00159_state"], "Low_quality_not_detected")
        self.assertEqual(df_master.loc["GCA_009.1", "atp_contrast1_cohort"], "EXCLUDED")

        # 10. Verify that all 4 ATP synthase modules (M00159, M00157, M00158, M00160) have identical 18 fields
        expected_suffixes = [
            "state", "confidence", "strict_comp", "weighted_comp",
            "n_strong", "n_supported_borderline", "n_borderline", "n_missing",
            "catalytic", "stalk", "membrane", "gene_clustered",
            "cluster_size_units", "cluster_span_bp", "max_gene_gap",
            "strand_consistency", "mean_margin", "missing_KOs"
        ]
        for mod in ["M00159", "M00157", "M00158", "M00160"]:
            for suffix in expected_suffixes:
                col = f"{mod}_{suffix}"
                self.assertIn(col, df_master.columns, f"Column {col} missing from master table!")
        
        # Verify M00158 and M00160 fields contain valid typed values
        for mod in ["M00158", "M00160"]:
            self.assertEqual(df_master.loc["GB_GCA_001.1", f"{mod}_state"], "Not_detected")
            self.assertFalse(bool(df_master.loc["GB_GCA_001.1", f"{mod}_catalytic"]))
            self.assertFalse(bool(df_master.loc["GB_GCA_001.1", f"{mod}_membrane"]))
            self.assertEqual(df_master.loc["GB_GCA_001.1", f"{mod}_n_strong"], 0)
            self.assertEqual(df_master.loc["GB_GCA_001.1", f"{mod}_confidence"], "High")

        print("✓ test_09_postprocess_atp_synthase passed all boundary gating validations and 4-module standardized field checks!")

    def test_10_analyze_atp_correlation(self):
        """Test analyze_atp_correlation.py for positive control isolation, Firth fallback, and CMH."""
        test_dir = os.path.join(TEST_TMP, "correlation_test")
        os.makedirs(test_dir, exist_ok=True)

        # Create master table with 12 genomes
        master_tsv = os.path.join(test_dir, "test_master.tsv")
        with open(master_tsv, "w") as f:
            f.write("genome_id\taccession\tcheckm2_completeness\tcheckm2_contamination\tgenome_size\tgtdb_taxonomy\tgtdb_representative\tatp_contrast1_cohort\tatp_contrast2_cohort\n")
            # 6 ATP_POSITIVE (3 in OrderA, 3 in OrderB)
            # 3 are VA_ONLY, 3 are F_LIKE_ONLY
            f.write("G_POS_1\tG_POS_1\t95.0\t1.0\t2000000\td__Archaea;o__OrderA\tt\tATP_POSITIVE\tVA_ONLY\n")
            f.write("G_POS_2\tG_POS_2\t95.0\t1.0\t2000000\td__Archaea;o__OrderA\tt\tATP_POSITIVE\tVA_ONLY\n")
            f.write("G_POS_3\tG_POS_3\t95.0\t1.0\t2000000\td__Archaea;o__OrderA\tt\tATP_POSITIVE\tVA_ONLY\n")
            f.write("G_POS_4\tG_POS_4\t95.0\t1.0\t2000000\td__Archaea;o__OrderB\tt\tATP_POSITIVE\tF_LIKE_ONLY\n")
            f.write("G_POS_5\tG_POS_5\t95.0\t1.0\t2000000\td__Archaea;o__OrderB\tt\tATP_POSITIVE\tF_LIKE_ONLY\n")
            f.write("G_POS_6\tG_POS_6\t95.0\t1.0\t2000000\td__Archaea;o__OrderB\tt\tATP_POSITIVE\tF_LIKE_ONLY\n")
            # 6 ATP_NEGATIVE (3 in OrderA, 3 in OrderB)
            f.write("G_NEG_1\tG_NEG_1\t95.0\t1.0\t1800000\td__Archaea;o__OrderA\tt\tATP_NEGATIVE\tEXCLUDED\n")
            f.write("G_NEG_2\tG_NEG_2\t95.0\t1.0\t1800000\td__Archaea;o__OrderA\tt\tATP_NEGATIVE\tEXCLUDED\n")
            f.write("G_NEG_3\tG_NEG_3\t95.0\t1.0\t1800000\td__Archaea;o__OrderA\tt\tATP_NEGATIVE\tEXCLUDED\n")
            f.write("G_NEG_4\tG_NEG_4\t95.0\t1.0\t1800000\td__Archaea;o__OrderB\tt\tATP_NEGATIVE\tEXCLUDED\n")
            f.write("G_NEG_5\tG_NEG_5\t95.0\t1.0\t1800000\td__Archaea;o__OrderB\tt\tATP_NEGATIVE\tEXCLUDED\n")
            f.write("G_NEG_6\tG_NEG_6\t95.0\t1.0\t1800000\td__Archaea;o__OrderB\tt\tATP_NEGATIVE\tEXCLUDED\n")

        # Create genome KO counts
        ko_counts_tsv = os.path.join(test_dir, "test_ko_counts.tsv")
        with open(ko_counts_tsv, "w") as f:
            f.write("genome_id\tko\tko_definition\n")
            # 1. K02117: Target definition KO (should be positive control)
            for g in ["G_POS_1", "G_POS_2", "G_POS_3"]:
                f.write(f"{g}\tK02117\tV_ATPase_A\n")
            # 2. K99001: Perfect candidate for ATP+ (Present in all 6 POS, absent in all 6 NEG) -> separation test!
            for g in ["G_POS_1", "G_POS_2", "G_POS_3", "G_POS_4", "G_POS_5", "G_POS_6"]:
                f.write(f"{g}\tK99001\tCore_Respiratory_Candidate\n")
            # 3. K99002: Present equally in all genomes (not significant)
            for g in ["G_POS_1", "G_POS_2", "G_POS_3", "G_POS_4", "G_POS_5", "G_POS_6",
                      "G_NEG_1", "G_NEG_2", "G_NEG_3", "G_NEG_4", "G_NEG_5", "G_NEG_6"]:
                f.write(f"{g}\tK99002\tHousekeeping_Polymerase\n")
            # 4. K99003: Specific candidate for VA_ONLY (Present in G_POS_1,2,3, absent in G_POS_4,5,6)
            for g in ["G_POS_1", "G_POS_2", "G_POS_3"]:
                f.write(f"{g}\tK99003\tVA_Specific_Coupling_Factor\n")

        try:
            import scipy
        except ImportError:
            self.skipTest("scipy is not installed; skipping test_10_analyze_atp_correlation")

        out_dir = os.path.join(test_dir, "out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "analyze_atp_correlation.py"),
            "-m", master_tsv,
            "-k", ko_counts_tsv,
            "-o", out_dir,
            "--min_ko_prevalence_count", "3",
            "--min_stratum_group_size", "2"
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        all_tsv = os.path.join(out_dir, "atp_gene_correlations_all_kos.tsv")
        pos_tsv = os.path.join(out_dir, "atp_gene_correlations_positive_controls.tsv")
        cand_tsv = os.path.join(out_dir, "atp_gene_correlations_candidates.tsv")
        top_tsv = os.path.join(out_dir, "atp_gene_correlations_top_associated.tsv")

        self.assertTrue(os.path.exists(all_tsv))
        self.assertTrue(os.path.exists(pos_tsv))
        self.assertTrue(os.path.exists(cand_tsv))
        self.assertTrue(os.path.exists(top_tsv))

        df_cand = pd.read_csv(cand_tsv, sep="\t")
        df_pos = pd.read_csv(pos_tsv, sep="\t")

        # Verify K02117 is isolated to positive controls and NEVER in candidates
        self.assertIn("K02117", df_pos["ko"].values)
        self.assertNotIn("K02117", df_cand["ko"].values)

        # Verify K99001 in Contrast 1 has complete separation handled via Firth fallback
        k99001_row = df_cand[(df_cand["ko"] == "K99001") & (df_cand["contrast"] == "ATP_any_presence")]
        self.assertEqual(len(k99001_row), 1)
        k99001_data = k99001_row.iloc[0]
        self.assertEqual(k99001_data["a"], 6)
        self.assertEqual(k99001_data["b"], 0)
        self.assertEqual(k99001_data["delta_prevalence"], 1.0)
        self.assertLess(k99001_data["fisher_p"], 0.01)
        self.assertIn(k99001_data["regression_method"], ["logistic", "firth_logistic"])
        self.assertTrue(k99001_data["regression_converged"])
        self.assertIn(k99001_data["candidate_tier"], ["Tier_1_robust", "Tier_2_raw_discovery"])

        print("✓ test_10_analyze_atp_correlation passed positive-control isolation, Firth separation, and CMH tests!")

    def test_11_postprocess_prevents_metadata_explosion(self):
        """Test that external metadata databases (e.g. 1000+ genomes) do not pollute small cohorts."""
        test_dir = os.path.join(TEST_TMP, "explosion_test")
        os.makedirs(test_dir, exist_ok=True)

        meta_path = os.path.join(test_dir, "large_mock_metadata.tsv")
        with open(meta_path, "w") as f:
            f.write("accession\tcheckm2_completeness\tcheckm2_contamination\tcheckm_completeness\tcheckm_contamination\tcontig_count\tgc_percentage\tgenome_size\tgtdb_representative\tgtdb_taxonomy\tgtdb_type_designation_ncbi_taxa\tgtdb_type_designation_ncbi_taxa_sources\tgtdb_type_species_of_genus\tl50_contigs\tmimag_high_quality\tmimag_low_quality\tmimag_medium_quality\tn50_contigs\tncbi_assembly_level\tncbi_genome_category\tncbi_isolation_source\tncbi_trna_count\ttrna_aa_count\ttrna_count\n")
            # Create 1000 genomes in metadata
            for i in range(1, 1001):
                f.write(f"GB_GCA_{i:07d}.1\t95.0\t1.0\t95.0\t1.0\t10\t45.0\t2000000\tt\td__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina mazei\t.\t.\t.\t2\tt\tf\tf\t500000\tComplete\tderived\thydrothermal\t45\t20\t45\n")

        # Create KO hits for ONLY 3 genomes
        ko_strict_path = os.path.join(test_dir, "small_cohort_protein_KO_all_strict.tsv")
        with open(ko_strict_path, "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tko_rank\tscore_margin\tis_strict\tpass_type\n")
            for g_id in ["GB_GCA_0000001.1", "GB_GCA_0000002.1", "GB_GCA_0000003.1"]:
                f.write(f"{g_id}___c1___prot00001\t{g_id}___c1___gene00001\t{g_id}\tK02117\tV_ATPase_A\t300.0\t200.0\tfull\t1e-60\t1\t100.0\tTrue\tscore_pass\n")

        out_dir = os.path.join(test_dir, "out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "postprocess_atp_synthase.py"),
            "-k", ko_strict_path,
            "-m", meta_path,
            "-o", out_dir
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        master_tsv = os.path.join(out_dir, "genome_atp_synthase_master_table.tsv")
        long_tsv = os.path.join(out_dir, "genome_atp_synthase_modules.long.tsv")
        subunits_tsv = os.path.join(out_dir, "genome_atp_synthase_subunits.tsv")
        pruned_meta_tsv = os.path.join(out_dir, "pruned_metadata.tsv")

        df_master = pd.read_csv(master_tsv, sep="\t")
        df_long = pd.read_csv(long_tsv, sep="\t")
        df_subunits = pd.read_csv(subunits_tsv, sep="\t")
        df_pruned = pd.read_csv(pruned_meta_tsv, sep="\t")

        # Crucial assertions: EXACT cohort sizes, 0 pollution from 1000 metadata genomes
        self.assertEqual(len(df_master), 3)
        self.assertEqual(len(df_long), 3 * 4)
        self.assertEqual(len(df_subunits), 3 * 41)
        self.assertEqual(len(df_pruned), 3)

        expected_genomes = {"GB_GCA_0000001.1", "GB_GCA_0000002.1", "GB_GCA_0000003.1"}
        self.assertEqual(set(df_master["genome_id"]), expected_genomes)
        self.assertEqual(set(df_long["genome_id"]), expected_genomes)
        self.assertEqual(set(df_subunits["genome_id"]), expected_genomes)
        self.assertEqual(set(df_pruned["accession"]), expected_genomes)

        print("✓ test_11_postprocess_prevents_metadata_explosion passed! 1000 metadata genomes did not pollute 3-genome cohort.")

if __name__ == "__main__":
    unittest.main()

