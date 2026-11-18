#!/usr/bin/env python3
"""
test_cog_annotation.py
Comprehensive test suite for COG2024-Compatible Footprint-Aware Annotation Engine.
"""

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "modules", "cog_annotation", "scripts")

sys.path.insert(0, SCRIPTS_DIR)
import assign_cog_footprints
import audit_database_provenance

class TestCOGAnnotation(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_cog_anno_")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_01_gate_0_database_provenance_audit(self):
        """Test Gate 0: Database provenance, checksums, and Cog_LE vs COG2024 intersection."""
        cddid_path = os.path.join(self.test_dir, "cddid.tbl")
        with open(cddid_path, "w") as f:
            f.write("# PSSM-Id\tAccession\tShort Name\tDescription\tPSSM-Length\n")
            f.write("223080\tCOG0001\tGluRS\tGlutamyl-tRNA synthetase\t450\n")
            f.write("223081\tCOG0002\tN-acetyl\tN-acetylglutamate synthase\t380\n")
            f.write("223082\tCOG0003\tOEC\tOEC enzyme\t250\n")
            f.write("223999\tCOG9999\tExtraCog\tLegacy CDD model\t300\n") # in Cog_LE only

        cog_def_path = os.path.join(self.test_dir, "cog-24.def.tab")
        with open(cog_def_path, "w") as f:
            f.write("COG0001\tJ\tGlutamyl-tRNA synthetase\tgltX\tTranslation\t12345\t1XYZ\n")
            f.write("COG0002\tE\tN-acetylglutamate synthase\targA\tArginine biosyn\t12346\t2XYZ\n")
            f.write("COG0003\tCH\tOEC enzyme\toecA\tEnergy\t12347\t3XYZ\n")
            f.write("COG8888\tT\tNewly added COG2024\tnewA\tSignaling\t12348\t4XYZ\n") # in COG2024 only

        cog_fun_path = os.path.join(self.test_dir, "cog-24.fun.tab")
        with open(cog_fun_path, "w") as f:
            f.write("INFORMATION STORAGE AND PROCESSING\n")
            f.write("[J]\t#FF0000\tTranslation, ribosomal structure and biogenesis\n")
            f.write("METABOLISM\n")
            f.write("[E]\t#00FF00\tAmino acid transport and metabolism\n")
            f.write("[C]\t#0000FF\tEnergy production and conversion\n")
            f.write("[H]\t#FFFF00\tCoenzyme transport and metabolism\n")
            f.write("CELLULAR PROCESSES AND SIGNALING\n")
            f.write("[T]\t#00FFFF\tSignal transduction mechanisms\n")

        bitscore_path = os.path.join(self.test_dir, "bitscore_specific.txt")
        with open(bitscore_path, "w") as f:
            f.write("# PSSM-Id / COG\tBitscore\n")
            f.write("223080\t85.5\n")
            f.write("COG0002\t110.0\n")

        audit_out = os.path.join(self.test_dir, "audit_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "audit_database_provenance.py"),
            "--cddid", cddid_path,
            "--cog_def", cog_def_path,
            "--cog_fun", cog_fun_path,
            "--bitscore_specific", bitscore_path,
            "-o", audit_out
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Gate 0 failed: {res.stderr}")

        report_file = os.path.join(audit_out, "cog_database_provenance_report.json")
        self.assertTrue(os.path.isfile(report_file))
        with open(report_file) as f:
            rep = json.load(f)

        self.assertEqual(rep["method_designation"], "CDD Cog_LE-based COG annotation mapped to COG2024 definitions")
        self.assertEqual(rep["n_cog2024"], 4)
        self.assertEqual(rep["n_cogle_cogs"], 4)
        self.assertEqual(rep["n_intersection"], 3)
        self.assertEqual(rep["n_only_cog2024"], 1) # COG8888
        self.assertEqual(rep["n_only_cogle"], 1)   # COG9999
        self.assertAlmostEqual(rep["coverage_fraction"], 0.75, places=4)
        self.assertEqual(rep["n_specific_bitscore_thresholds"], 2)
        self.assertIsNotNone(rep["database_checksums"]["cog_fun_sha256"])
        self.assertIsNotNone(rep["database_checksums"]["bitscore_specific_sha256"])

        # Check official category definition table
        cat_file = os.path.join(audit_out, "cog_category_definition.tsv")
        self.assertTrue(os.path.isfile(cat_file))
        with open(cat_file) as f:
            c_lines = [l.strip().split("\t") for l in f if l.strip()]
        self.assertEqual(c_lines[0], ["category_code", "category_name", "functional_class"])
        cat_codes = {r[0]: r[1] for r in c_lines[1:]}
        self.assertIn("J", cat_codes)
        self.assertIn("E", cat_codes)

        # Check mapping table
        mapping_file = os.path.join(audit_out, "cdd_pssm_mapping.tsv")
        self.assertTrue(os.path.isfile(mapping_file))
        with open(mapping_file, "r", encoding="utf-8") as f:
            lines = [l.strip().split("\t") for l in f]
        self.assertEqual(lines[0][0], "pssm_id")
        self.assertIn("specific_bitscore_threshold", lines[0])
        spec_col = lines[0].index("specific_bitscore_threshold")
        self.assertEqual(len(lines), 5) # header + 4 models
        pssm_row_map = {r[0]: r for r in lines[1:]}
        self.assertEqual(pssm_row_map["223080"][spec_col], "85.5")
        self.assertEqual(pssm_row_map["223081"][spec_col], "110.0")
        self.assertEqual(pssm_row_map["223082"][spec_col], "NA")

    def test_02_hsp_aggregation_and_match_codes(self):
        """Test collinear HSP-to-footprint aggregation and COG2024 75% rule (Match Code 0-3)."""
        mapping_path = os.path.join(self.test_dir, "cdd_pssm_mapping.tsv.gz")
        with gzip.open(mapping_path, "wt") as f:
            f.write("pssm_id\tcog_id\tpssm_length\tcdd_name\tcdd_desc\tcog_category\tcog_name\tspecific_bitscore_threshold\n")
            f.write("223080\tCOG0001\t300\tModel1\tDesc1\tJ\tName1\t100.0\n")

        # Two collinear HSPs on protein_1 (400 aa) matching PSSM 223080 (300 aa)
        # HSP1: query 20-150, pssm 15-145 (span 131 aa)
        # HSP2: query 165-340, pssm 155-300 (span 176 aa)
        # Combined footprint: query 20-340 (321 aa, query_cov = 321/400 = 80.25% >= 75%)
        # Combined profile: pssm 15-300 (286 aa, profile_cov = 286/300 = 95.3% >= 75%)
        # Expected: match_code = 0 (full_match)
        raw_hits_path = os.path.join(self.test_dir, "raw_hits.tsv")
        with open(raw_hits_path, "w") as f:
            f.write("protein_1\tgnl|CDD|223080\t45.0\t131\t400\t300\t20\t150\t15\t145\t1e-25\t95.0\n")
            f.write("protein_1\tgnl|CDD|223080\t42.0\t176\t400\t300\t165\t340\t155\t300\t1e-35\t120.0\n")

        assigned_out = os.path.join(self.test_dir, "assigned.tsv")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "assign_cog_footprints.py"),
            "-i", raw_hits_path,
            "-m", mapping_path,
            "-o", assigned_out,
            "--evalue", "1e-4"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Assignment failed: {res.stderr}")

        with open(assigned_out) as f:
            header = f.readline().strip().split("\t")
            rows = [l.strip().split("\t") for l in f if l.strip()]

        col = {h: i for i, h in enumerate(header)}
        self.assertEqual(len(rows), 1, "Expected exactly 1 aggregated footprint")
        r = rows[0]
        self.assertEqual(r[col["cog_id"]], "COG0001")
        self.assertEqual(r[col["qstart"]], "20")
        self.assertEqual(r[col["qend"]], "340")
        self.assertEqual(r[col["cog_match_code"]], "0")
        self.assertEqual(r[col["cog_match_class"]], "full_match")
        self.assertEqual(r[col["specific_bitscore_threshold"]], "100.0")
        self.assertEqual(r[col["passes_specific_threshold"]], "True")
        self.assertEqual(r[col["profile_coverage_class"]], "profile_complete")
        self.assertEqual(r[col["query_coverage_class"]], "query_complete")

    def test_03_competing_hits_and_ambiguity(self):
        """Test overlapping competing profiles (overlap >= 50%) and ambiguity score margin."""
        mapping_path = os.path.join(self.test_dir, "cdd_pssm_mapping.tsv.gz")
        with gzip.open(mapping_path, "wt") as f:
            f.write("pssm_id\tcog_id\tpssm_length\tcdd_name\tcdd_desc\tcog_category\tcog_name\n")
            f.write("223080\tCOG0001\t250\tModel1\tDesc1\tJ\tName1\n")
            f.write("223081\tCOG0002\t250\tModel2\tDesc2\tE\tName2\n")

        # Two competing hits on aa 50-280 (overlap > 80%)
        # Hit 1 (COG0001): bitscore 150.0
        # Hit 2 (COG0002): bitscore 140.0 -> delta_bitscore = 10.0, R_delta = 10/150 = 0.0667 < 0.15
        raw_hits_path = os.path.join(self.test_dir, "raw_competing.tsv")
        with open(raw_hits_path, "w") as f:
            f.write("protein_comp\tgnl|CDD|223080\t50.0\t230\t500\t250\t50\t280\t10\t240\t1e-50\t150.0\n")
            f.write("protein_comp\tgnl|CDD|223081\t48.0\t225\t500\t250\t55\t280\t15\t240\t1e-45\t140.0\n")

        assigned_out = os.path.join(self.test_dir, "assigned_comp.tsv")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "assign_cog_footprints.py"),
            "-i", raw_hits_path,
            "-m", mapping_path,
            "-o", assigned_out,
            "--ambiguity_ratio_cutoff", "0.15",
            "--ambiguity_delta_cutoff", "20.0"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Assignment failed: {res.stderr}")

        with open(assigned_out) as f:
            header = f.readline().strip().split("\t")
            rows = [l.strip().split("\t") for l in f if l.strip()]

        col = {h: i for i, h in enumerate(header)}
        self.assertEqual(len(rows), 1, "Competing hit should yield single primary footprint")
        r = rows[0]
        self.assertEqual(r[col["cog_id"]], "COG0001")
        self.assertEqual(r[col["region_status"]], "ambiguous")
        self.assertEqual(r[col["top2_cog"]], "COG0002")
        self.assertEqual(r[col["has_ambiguous_region"]], "True")
        self.assertEqual(r[col["has_competing_profiles"]], "True")
        self.assertEqual(r[col["architecture_class"]], "single_cog")

    def test_04_multidomain_architecture_preservation(self):
        """Test non-overlapping hits (overlap < 20%) retained as multi_cog architecture."""
        mapping_path = os.path.join(self.test_dir, "cdd_pssm_mapping.tsv.gz")
        with gzip.open(mapping_path, "wt") as f:
            f.write("pssm_id\tcog_id\tpssm_length\tcdd_name\tcdd_desc\tcog_category\tcog_name\n")
            f.write("223080\tCOG0001\t200\tModel1\tDesc1\tJ\tName1\n")
            f.write("223082\tCOG0003\t200\tModel3\tDesc3\tC\tName3\n")

        # Protein 800 aa:
        # Domain 1: aa 30-220 (COG0001)
        # Domain 2: aa 500-690 (COG0003)
        raw_hits_path = os.path.join(self.test_dir, "raw_multi.tsv")
        with open(raw_hits_path, "w") as f:
            f.write("protein_multi\tgnl|CDD|223080\t55.0\t190\t800\t200\t30\t220\t5\t195\t1e-60\t180.0\n")
            f.write("protein_multi\tgnl|CDD|223082\t52.0\t190\t800\t200\t500\t690\t10\t200\t1e-55\t170.0\n")

        assigned_out = os.path.join(self.test_dir, "assigned_multi.tsv")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "assign_cog_footprints.py"),
            "-i", raw_hits_path,
            "-m", mapping_path,
            "-o", assigned_out
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Assignment failed: {res.stderr}")

        with open(assigned_out) as f:
            header = f.readline().strip().split("\t")
            rows = [l.strip().split("\t") for l in f if l.strip()]

        col = {h: i for i, h in enumerate(header)}
        self.assertEqual(len(rows), 2, "Both non-overlapping footprints must be retained!")
        self.assertEqual(rows[0][col["cog_id"]], "COG0001")
        self.assertEqual(rows[1][col["cog_id"]], "COG0003")
        self.assertEqual(rows[0][col["architecture_class"]], "multi_cog")
        self.assertEqual(rows[0][col["is_multi_cog"]], "True")
        self.assertEqual(rows[0][col["cog_footprint_architecture"]], "COG0001[30-220];COG0003[500-690]")

    def test_05_strict_integer_length_validation(self):
        """Test strict integer slen check raises DatabaseSynchronizationError on mismatch."""
        mapping_path = os.path.join(self.test_dir, "cdd_pssm_mapping.tsv.gz")
        with gzip.open(mapping_path, "wt") as f:
            f.write("pssm_id\tcog_id\tpssm_length\tcdd_name\tcdd_desc\tcog_category\tcog_name\n")
            f.write("223080\tCOG0001\t300\tModel1\tDesc1\tJ\tName1\n") # Length is 300

        # Raw hit has slen=305 (mismatch)
        raw_hits_path = os.path.join(self.test_dir, "raw_mismatch.tsv")
        with open(raw_hits_path, "w") as f:
            f.write("protein_err\tgnl|CDD|223080\t50.0\t200\t500\t305\t50\t250\t10\t210\t1e-40\t130.0\n")

        assigned_out = os.path.join(self.test_dir, "assigned_mismatch.tsv")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "assign_cog_footprints.py"),
            "-i", raw_hits_path,
            "-m", mapping_path,
            "-o", assigned_out
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertNotEqual(res.returncode, 0, "Should fail when slen does not match cddid length!")
        self.assertIn("DatabaseSynchronizationError", res.stderr + res.stdout)

    def test_06_expansion_to_master_and_summary(self):
        """Test streaming expansion to protein_COG_master.tsv.gz and protein_COG_summary.tsv.gz."""
        # 1. Mock assigned shard
        shard_dir = os.path.join(self.test_dir, "assigned_shards")
        os.makedirs(shard_dir, exist_ok=True)
        assigned_file = os.path.join(shard_dir, "shard_0001.assigned.tsv.gz")

        cols = [
            "protein_id", "cog_id", "pssm_id", "cog_category", "cog_name",
            "qstart", "qend", "qlen", "sstart", "send", "slen", "aln_len",
            "query_cov", "profile_cov", "evalue", "bitscore",
            "specific_bitscore_threshold", "passes_specific_threshold",
            "cog_match_code", "cog_match_class", "profile_coverage_class", "query_coverage_class",
            "footprint_rank", "region_id", "region_status", "architecture_class",
            "has_competing_profiles", "has_ambiguous_region", "is_multi_cog", "is_repeated_cog",
            "top2_cog", "delta_bitscore", "bitscore_ratio", "delta_bitscore_fraction",
            "likely_truncated_gene", "truncation_evidence_count", "prodigal_partial", "near_contig_edge", "short_vs_cog",
            "cog_footprint_architecture"
        ]
        with gzip.open(assigned_file, "wt") as f:
            f.write("\t".join(cols) + "\n")
            # unique_1 has 2 footprints (multi_cog)
            f.write("unique_1\tCOG0001\t223080\tJ\tName1\t20\t200\t600\t10\t190\t200\t181\t0.3017\t0.905\t1e-50\t150.0\t85.5\tTrue\t1\tprofile_complete_query_partial\tprofile_complete\tquery_partial\t1\t20-200\tunique\tmulti_cog\tFalse\tFalse\tTrue\tFalse\tNA\tNA\tNA\tNA\tFalse\t1\tFalse\tFalse\tTrue\tCOG0001[20-200];COG0003[350-550]\n")
            f.write("unique_1\tCOG0003\t223082\tCH\tName3\t350\t550\t600\t10\t210\t210\t201\t0.335\t0.957\t1e-60\t175.0\tNA\tNA\t1\tprofile_complete_query_partial\tprofile_complete\tquery_partial\t2\t350-550\tunique\tmulti_cog\tFalse\tFalse\tTrue\tFalse\tNA\tNA\tNA\tNA\tFalse\t1\tFalse\tFalse\tTrue\tCOG0001[20-200];COG0003[350-550]\n")
            # unique_2 has 1 footprint
            f.write("unique_2\tCOG0002\t223081\tE\tName2\t10\t300\t320\t5\t295\t300\t291\t0.9094\t0.97\t1e-80\t220.0\t110.0\tTrue\t0\tfull_match\tprofile_complete\tquery_complete\t1\t10-300\tunique\tsingle_cog\tFalse\tFalse\tFalse\tFalse\tNA\tNA\tNA\tNA\tFalse\t0\tFalse\tFalse\tFalse\tCOG0002[10-300]\n")

        # 2. Mock protein_to_unique.tsv.gz (4 original proteins across 2 genomes, including 1 unannotated)
        mapping_file = os.path.join(self.test_dir, "protein_to_unique.tsv.gz")
        with gzip.open(mapping_file, "wt") as f:
            f.write("protein_id\tunique_id\taa_length\tgenome_id\n")
            f.write("prot_g1_1\tunique_1\t600\tgenome_1\n")
            f.write("prot_g1_2\tunique_2\t320\tgenome_1\n")
            f.write("prot_g2_1\tunique_1\t600\tgenome_2\n") # shared identical protein
            f.write("prot_g2_unannot\tunique_unannot\t150\tgenome_2\n") # unannotated

        master_out = os.path.join(self.test_dir, "master_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "expand_dedup_cog.py"),
            "-s", shard_dir,
            "-m", mapping_file,
            "-o", master_out
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Expansion failed: {res.stderr}")

        master_tsv = os.path.join(master_out, "protein_COG_master.tsv")
        summary_tsv = os.path.join(master_out, "protein_COG_summary.tsv")
        self.assertTrue(os.path.isfile(master_tsv))
        self.assertTrue(os.path.isfile(summary_tsv))

        # Master rows: prot_g1_1 (2 fps) + prot_g1_2 (1 fp) + prot_g2_1 (2 fps) = 5 rows
        with open(master_tsv, "r", encoding="utf-8") as f:
            m_lines = [l.strip().split("\t") for l in f if l.strip()]
        self.assertEqual(len(m_lines) - 1, 5)
        m_header = m_lines[0]
        self.assertIn("specific_bitscore_threshold", m_header)
        self.assertIn("passes_specific_threshold", m_header)
        m_spec_idx = m_header.index("specific_bitscore_threshold")
        m_pass_idx = m_header.index("passes_specific_threshold")
        # Row 1 (prot_g1_1 COG0001) should have 85.5 and True
        self.assertEqual(m_lines[1][m_spec_idx], "85.5")
        self.assertEqual(m_lines[1][m_pass_idx], "True")

        # Summary rows: exactly 4 proteins (100% accounting)
        with open(summary_tsv, "r", encoding="utf-8") as f:
            s_lines = [l.strip().split("\t") for l in f if l.strip()]
        self.assertEqual(len(s_lines) - 1, 4)
        unannot_row = [r for r in s_lines[1:] if r[1] == "prot_g2_unannot"][0]
        self.assertEqual(unannot_row[4], "0")
        self.assertEqual(unannot_row[5], "unassigned")

    def test_07_matrix_generation_and_category_conservation(self):
        """Test comparative matrices generation and exact protein-normalized conservation."""
        # Setup master TSV from previous expansion
        self.test_06_expansion_to_master_and_summary()
        master_tsv = os.path.join(self.test_dir, "master_out", "protein_COG_master.tsv")
        mat_out = os.path.join(self.test_dir, "matrices_out")

        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "generate_cog_matrices.py"),
            "-m", master_tsv,
            "-o", mat_out,
            "--export_dense"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Matrix generation failed: {res.stderr}")

        # Check category protein-normalized conservation:
        # In genome_1: 2 annotated proteins (prot_g1_1 and prot_g1_2) -> sum must be 2.0000
        # In genome_2: 1 annotated protein (prot_g2_1) -> sum must be 1.0000
        prot_norm_file = os.path.join(mat_out, "genome_COG_category_protein_normalized.tsv")
        with open(prot_norm_file) as f:
            header = f.readline().strip().split("\t")
            for line in f:
                parts = line.strip().split("\t")
                gid = parts[0]
                tot = float(parts[-1])
                if gid == "genome_1":
                    self.assertAlmostEqual(tot, 2.0, places=3)
                elif gid == "genome_2":
                    self.assertAlmostEqual(tot, 1.0, places=3)

        # Check category footprint-normalized conservation:
        # In genome_1: prot_g1_1 (2 fps) + prot_g1_2 (1 fp) = 3 fps -> sum must be 3.0000
        # In genome_2: prot_g2_1 (2 fps) = 2 fps -> sum must be 2.0000
        fp_norm_file = os.path.join(mat_out, "genome_COG_category_footprint_normalized.tsv")
        with open(fp_norm_file) as f:
            header = f.readline().strip().split("\t")
            for line in f:
                parts = line.strip().split("\t")
                gid = parts[0]
                tot = float(parts[-1])
                if gid == "genome_1":
                    self.assertAlmostEqual(tot, 3.0, places=3)
                elif gid == "genome_2":
                    self.assertAlmostEqual(tot, 2.0, places=3)

    def test_08_six_gate_qc_validator(self):
        """Test complete 6-gate QC validator."""
        self.test_01_gate_0_database_provenance_audit()
        self.test_07_matrix_generation_and_category_conservation()

        master_tsv = os.path.join(self.test_dir, "master_out", "protein_COG_master.tsv")
        summary_tsv = os.path.join(self.test_dir, "master_out", "protein_COG_summary.tsv")
        mapping_file = os.path.join(self.test_dir, "protein_to_unique.tsv.gz")
        prov_report = os.path.join(self.test_dir, "audit_out", "cog_database_provenance_report.json")
        mat_dir = os.path.join(self.test_dir, "matrices_out")
        qc_out = os.path.join(self.test_dir, "qc_out")

        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "validate_cog_qc.py"),
            "-m", master_tsv,
            "-s", summary_tsv,
            "-p", mapping_file,
            "-r", prov_report,
            "-d", mat_dir,
            "-o", qc_out
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"QC validation failed:\n{res.stderr}\n{res.stdout}")

        with open(os.path.join(qc_out, "cog_qc_report.json")) as f:
            rep = json.load(f)
        self.assertEqual(rep["overall_qc_status"], "PASSED")
        self.assertEqual(rep["gate_0_provenance"]["status"], "PASS")
        self.assertEqual(rep["gate_1_completeness"]["status"], "PASS")
        self.assertEqual(rep["gate_2_expansion_invariance"]["status"], "PASS")
        self.assertEqual(rep["gate_3_strict_length_qc"]["status"], "PASS")
        self.assertEqual(rep["gate_4_footprint_separation"]["status"], "PASS")
        self.assertEqual(rep["gate_5_category_conservation"]["status"], "PASS")

    def test_09_ground_truth_recovery_benchmark(self):
        """Test Tier 1 Official COG Member Recovery Benchmark with IoU and coordinate delta analysis."""
        self.test_06_expansion_to_master_and_summary()
        master_tsv = os.path.join(self.test_dir, "master_out", "protein_COG_master.tsv")

        # Mock ground truth for prot_g1_1
        # Official footprint: COG0001 at 25-205, match_code = 1
        gt_file = os.path.join(self.test_dir, "official_ground_truth.tsv")
        with open(gt_file, "w") as f:
            f.write("protein_id\tcog_id\tqstart\tqend\tmatch_code\n")
            f.write("prot_g1_1\tCOG0001\t25\t205\t1\n")

        bench_out = os.path.join(self.test_dir, "bench_out")
        cmd = [
            sys.executable,
            os.path.join(SCRIPTS_DIR, "benchmark_cog_recovery.py"),
            "-p", master_tsv,
            "-g", gt_file,
            "-o", bench_out
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Benchmark failed: {res.stderr}")

        with open(os.path.join(bench_out, "cog_recovery_benchmark_report.json")) as f:
            b_rep = json.load(f)
        self.assertEqual(b_rep["n_evaluated_proteins"], 1)
        self.assertEqual(b_rep["exact_cog_recovery_rate"], 1.0)
        self.assertGreater(b_rep["mean_footprint_iou"], 0.80)
        self.assertEqual(b_rep["match_code_concordance_rate"], 1.0)

    def test_10_multi_platform_database_resolution(self):
        """Test Gate 0: Platform configuration database overrides module config across gses0 and icdli1."""
        sys.path.insert(0, REPO_ROOT)
        from common.scripts.runtime import apply_runtime_config

        # 1. Test gses0 platform resolution
        gses0_config = {
            "databases": {
                "cdd": "/dssg/home/acct-schhui/gses0/Database/CDD",
                "cogle_db": "/dssg/home/acct-schhui/gses0/Database/CDD/Cog_LE",
                "cddid": "/dssg/home/acct-schhui/gses0/Database/CDD/cddid.tbl",
                "cog_def": "/dssg/home/acct-schhui/gses0/Database/COG/cog-24.def.tab",
                "cog_fun": "/dssg/home/acct-schhui/gses0/Database/COG/cog-24.fun.tab",
                "bitscore_specific": "/dssg/home/acct-schhui/gses0/Database/CDD/bitscore_specific.txt",
            },
            "cog_annotation": {
                "databases": {
                    "cddid": "",
                    "cog_def": "",
                    "cog_fun": "",
                    "cogle_db": "",
                    "bitscore_specific": "",
                }
            }
        }
        apply_runtime_config(gses0_config)
        cog_dbs = gses0_config["cog_annotation"]["databases"]
        self.assertEqual(cog_dbs["cddid"], "/dssg/home/acct-schhui/gses0/Database/CDD/cddid.tbl")
        self.assertEqual(cog_dbs["cogle_db"], "/dssg/home/acct-schhui/gses0/Database/CDD/Cog_LE")
        self.assertEqual(cog_dbs["cog_def"], "/dssg/home/acct-schhui/gses0/Database/COG/cog-24.def.tab")
        self.assertEqual(cog_dbs["cog_fun"], "/dssg/home/acct-schhui/gses0/Database/COG/cog-24.fun.tab")
        self.assertEqual(cog_dbs["bitscore_specific"], "/dssg/home/acct-schhui/gses0/Database/CDD/bitscore_specific.txt")

        # 2. Test platform overrides stale module paths (cross-cluster isolation)
        gses0_with_stale_module = {
            "databases": {
                "cddid": "/dssg/home/acct-schhui/gses0/Database/CDD/cddid.tbl",
                "cog_def": "/dssg/home/acct-schhui/gses0/Database/COG/cog-24.def.tab",
                "cogle_db": "/dssg/home/acct-schhui/gses0/Database/CDD/Cog_LE",
            },
            "cog_annotation": {
                "databases": {
                    "cddid": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/CDD/cddid.tbl",
                    "cog_def": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/cog-24.def.tab",
                    "cogle_db": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/CDD/Cog",
                }
            }
        }
        apply_runtime_config(gses0_with_stale_module)
        stale_dbs = gses0_with_stale_module["cog_annotation"]["databases"]
        self.assertEqual(stale_dbs["cddid"], "/dssg/home/acct-schhui/gses0/Database/CDD/cddid.tbl")
        self.assertEqual(stale_dbs["cogle_db"], "/dssg/home/acct-schhui/gses0/Database/CDD/Cog_LE")
        self.assertEqual(stale_dbs["cog_def"], "/dssg/home/acct-schhui/gses0/Database/COG/cog-24.def.tab")

if __name__ == "__main__":
    unittest.main()
