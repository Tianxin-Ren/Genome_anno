#!/usr/bin/env python3
"""
test_combine_kegg_cog.py
Unit and integration test suite for combine_kegg_cog_annotations.py.
Verifies:
1. 100% Proteome completeness in joint protein_annotation_master.tsv
2. Correct 4-way annotation status categorization (both_annotated, ko_only, cog_only, unassigned)
3. Translation of category letters to full descriptions & functional classes
4. Strict conservation and accuracy of genome-level macro statistics
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
COMBINE_SCRIPT = os.path.join(REPO_ROOT, "modules", "mag_annotation", "scripts", "combine_kegg_cog_annotations.py")


class TestCombineKeggCog(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_combine_kegg_cog_")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_joint_kegg_cog_integration(self):
        # 1. Create mock category definitions
        cat_def_path = os.path.join(self.test_dir, "cog_category_definition.tsv")
        with open(cat_def_path, "w", encoding="utf-8") as f:
            f.write("category_code\tcategory_name\tfunctional_class\n")
            f.write("J\tTranslation, ribosomal structure\tINFORMATION STORAGE AND PROCESSING\n")
            f.write("E\tAmino acid transport and metabolism\tMETABOLISM\n")
            f.write("C\tEnergy production and conversion\tMETABOLISM\n")
            f.write("H\tCoenzyme transport and metabolism\tMETABOLISM\n")

        # 2. Create mock CDD PSSM mapping
        pssm_mapping_path = os.path.join(self.test_dir, "cdd_pssm_mapping.tsv")
        with open(pssm_mapping_path, "w", encoding="utf-8") as f:
            f.write("pssm_id\tcog_id\tpssm_length\tcdd_name\tcdd_desc\tcog_category\tcog_name\tspecific_bitscore_threshold\n")
            f.write("223080\tCOG0001\t450\tGluRS\tGlutamyl-tRNA syn\tJ\tGlutamyl-tRNA synthetase\t85.5\n")
            f.write("223081\tCOG0002\t380\tN-acetyl\tN-acetyl synth\tE\tN-acetylglutamate synthase\t110.0\n")
            f.write("223082\tCOG0003\t250\tOEC\tOEC enzyme\tCH\tOEC enzyme complex\tNA\n")

        # 3. Create mock protein mapping (5 proteins across 2 genomes)
        # genome_A: prot_A_1, prot_A_2, prot_A_3
        # genome_B: prot_B_1, prot_B_2
        map_path = os.path.join(self.test_dir, "protein_to_unique.tsv.gz")
        with gzip.open(map_path, "wt", encoding="utf-8") as f:
            f.write("protein_id\tunique_id\taa_length\tgenome_id\n")
            f.write("genome_A___prot_1\tuniq_1\t450\tgenome_A\n")
            f.write("genome_A___prot_2\tuniq_2\t300\tgenome_A\n")
            f.write("genome_A___prot_3\tuniq_3\t250\tgenome_A\n")
            f.write("genome_B___prot_1\tuniq_4\t280\tgenome_B\n")
            f.write("genome_B___prot_2\tuniq_5\t120\tgenome_B\n")

        # 4. Create mock protein_KO_best.tsv
        # Hits: prot_A_1 (K01885), prot_A_2 (K00001), prot_B_1 (K00002)
        # Missing from KO: prot_A_3, prot_B_2
        ko_path = os.path.join(self.test_dir, "protein_KO_best.tsv")
        with open(ko_path, "w", encoding="utf-8") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tscore_margin\tis_strict\tpass_type\n")
            f.write("genome_A___prot_1\tgenome_A_1\tgenome_A\tK01885\tglutamyl-tRNA synthetase\t250.5\t200.0\tfull\t1e-60\t50.5\tTrue\tscore_pass\n")
            f.write("genome_A___prot_2\tgenome_A_2\tgenome_A\tK00001\talcohol dehydrogenase\t180.2\t150.0\tfull\t1e-45\t30.2\tTrue\tscore_pass\n")
            f.write("genome_B___prot_1\tgenome_B_1\tgenome_B\tK00002\tAKR1A1 dehydrogenase\t195.0\t170.0\tfull\t1e-50\t25.0\tTrue\tscore_pass\n")

        # 5. Create mock protein_COG_summary.tsv (1 row per protein)
        # Hits:
        # prot_A_1 -> COG0001 (has KO & COG -> both_annotated)
        # prot_A_2 -> unassigned (has KO, no COG -> ko_only)
        # prot_A_3 -> COG0002 (no KO, has COG -> cog_only)
        # prot_B_1 -> COG0003 (has KO & COG -> both_annotated)
        # prot_B_2 -> unassigned (no KO, no COG -> unassigned)
        cog_path = os.path.join(self.test_dir, "protein_COG_summary.tsv")
        cog_header = [
            "genome_id", "protein_id", "unique_protein_id", "aa_length",
            "n_footprints", "primary_cog", "cog_footprint_architecture", "architecture_class",
            "cog_match_codes", "has_competing_profiles", "has_ambiguous_region",
            "is_multi_cog", "is_repeated_cog", "likely_truncated_gene", "method_designation"
        ]
        with open(cog_path, "w", encoding="utf-8") as f:
            f.write("\t".join(cog_header) + "\n")
            f.write("genome_A\tgenome_A___prot_1\tuniq_1\t450\t1\tCOG0001\tCOG0001[20-440]\tsingle_cog\t0\tFalse\tFalse\tFalse\tFalse\tFalse\tCDD\n")
            f.write("genome_A\tgenome_A___prot_2\tuniq_2\t300\t0\tunassigned\tnone\tunassigned\tnone\tFalse\tFalse\tFalse\tFalse\tFalse\tCDD\n")
            f.write("genome_A\tgenome_A___prot_3\tuniq_3\t250\t1\tCOG0002\tCOG0002[10-240]\tsingle_cog\t1\tFalse\tFalse\tFalse\tFalse\tFalse\tCDD\n")
            f.write("genome_B\tgenome_B___prot_1\tuniq_4\t280\t1\tCOG0003\tCOG0003[15-265]\tsingle_cog\t0\tFalse\tFalse\tFalse\tFalse\tFalse\tCDD\n")
            f.write("genome_B\tgenome_B___prot_2\tuniq_5\t120\t0\tunassigned\tnone\tunassigned\tnone\tFalse\tFalse\tFalse\tFalse\tFalse\tCDD\n")

        # 6. Execute integration
        out_dir = os.path.join(self.test_dir, "combined_out")
        cmd = [
            sys.executable,
            COMBINE_SCRIPT,
            "-p", map_path,
            "-k", ko_path,
            "-c", cog_path,
            "--cat_def", cat_def_path,
            "-m", pssm_mapping_path,
            "-o", out_dir
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, f"Combine script failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

        # 7. Check Deliverable 1: protein_annotation_master.tsv
        master_file = os.path.join(out_dir, "protein_annotation_master.tsv")
        self.assertTrue(os.path.isfile(master_file), "Master table must exist and be uncompressed .tsv")

        with open(master_file, "r", encoding="utf-8") as f:
            lines = [l.strip().split("\t") for l in f if l.strip()]

        header = lines[0]
        self.assertIn("annotation_status", header)
        self.assertIn("ko_id", header)
        self.assertIn("cog_id", header)
        self.assertIn("cog_category_desc", header)
        self.assertIn("cog_functional_class", header)

        # Total rows must strictly equal 5 proteins (100% proteome completeness)
        self.assertEqual(len(lines) - 1, 5, "Every protein in proteome must have exactly 1 row")

        status_map = {r[header.index("protein_id")]: r[header.index("annotation_status")] for r in lines[1:]}
        self.assertEqual(status_map["genome_A___prot_1"], "both_annotated")
        self.assertEqual(status_map["genome_A___prot_2"], "ko_only")
        self.assertEqual(status_map["genome_A___prot_3"], "cog_only")
        self.assertEqual(status_map["genome_B___prot_1"], "both_annotated")
        self.assertEqual(status_map["genome_B___prot_2"], "unassigned")

        # Check multi-category parsing for COG0003 (categories 'CH')
        row_b1 = [r for r in lines[1:] if r[header.index("protein_id")] == "genome_B___prot_1"][0]
        self.assertEqual(row_b1[header.index("cog_category")], "CH")
        desc_col = row_b1[header.index("cog_category_desc")]
        self.assertIn("Energy production", desc_col)
        self.assertIn("Coenzyme transport", desc_col)

        # 8. Check Deliverable 2: genome_annotation_summary.tsv
        summary_file = os.path.join(out_dir, "genome_annotation_summary.tsv")
        self.assertTrue(os.path.isfile(summary_file), "Summary table must exist and be uncompressed .tsv")

        with open(summary_file, "r", encoding="utf-8") as f:
            s_lines = [l.strip().split("\t") for l in f if l.strip()]

        s_header = s_lines[0]
        self.assertEqual(len(s_lines) - 1, 2, "Must summarize both genomes")

        s_dict = {r[0]: {s_header[i]: r[i] for i in range(len(s_header))} for r in s_lines[1:]}

        # Validate genome_A macro stats
        gA = s_dict["genome_A"]
        self.assertEqual(int(gA["total_proteins"]), 3)
        self.assertEqual(int(gA["annotated_proteins"]), 3)
        self.assertEqual(int(gA["dual_annotated_proteins"]), 1)
        self.assertEqual(int(gA["ko_only_proteins"]), 1)
        self.assertEqual(int(gA["cog_only_proteins"]), 1)
        self.assertEqual(int(gA["unannotated_proteins"]), 0)
        self.assertAlmostEqual(float(gA["overall_annotation_rate"]), 1.0000, places=3)
        self.assertAlmostEqual(float(gA["ko_annotation_rate"]), 2.0 / 3.0, places=3)
        self.assertAlmostEqual(float(gA["cog_annotation_rate"]), 2.0 / 3.0, places=3)

        # Validate genome_B macro stats
        gB = s_dict["genome_B"]
        self.assertEqual(int(gB["total_proteins"]), 2)
        self.assertEqual(int(gB["annotated_proteins"]), 1)
        self.assertEqual(int(gB["dual_annotated_proteins"]), 1)
        self.assertEqual(int(gB["ko_only_proteins"]), 0)
        self.assertEqual(int(gB["cog_only_proteins"]), 0)
        self.assertEqual(int(gB["unannotated_proteins"]), 1)
        self.assertAlmostEqual(float(gB["overall_annotation_rate"]), 0.5000, places=3)
        self.assertAlmostEqual(float(gB["ko_annotation_rate"]), 0.5000, places=3)
        self.assertAlmostEqual(float(gB["cog_annotation_rate"]), 0.5000, places=3)


if __name__ == "__main__":
    unittest.main()
