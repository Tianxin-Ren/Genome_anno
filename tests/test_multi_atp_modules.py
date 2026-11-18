#!/usr/bin/env python3
import unittest
import tempfile
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from scripts.compile_atp_synteny_matrix import determine_overall_atp_state, compile_matrices

class TestMultiATPModules(unittest.TestCase):
    def test_determine_overall_atp_state(self):
        # Case 1: Canonical V/A-type retained, others absent
        m1 = {
            "M00159_inferred_loss": "INFERRED_LOSS_NO",
            "M00157_inferred_loss": "ANCESTRAL_ABSENT",
            "M00158_inferred_loss": "ANCESTRAL_ABSENT",
            "M00160_inferred_loss": "ANCESTRAL_ABSENT"
        }
        self.assertEqual(determine_overall_atp_state(m1), "CANONICAL_V_TYPE_RETAINED")

        # Case 2: V/A-type lost, F-type substituted
        m2 = {
            "M00159_inferred_loss": "INFERRED_LOSS_HIGH",
            "M00157_inferred_loss": "INFERRED_LOSS_NO",
            "M00158_inferred_loss": "ANCESTRAL_ABSENT",
            "M00160_inferred_loss": "ANCESTRAL_ABSENT"
        }
        self.assertEqual(determine_overall_atp_state(m2), "F_TYPE_SUBSTITUTED")

        # Case 3: Both V-type and F-type present
        m3 = {
            "M00159_inferred_loss": "INFERRED_LOSS_NO",
            "M00157_inferred_loss": "INFERRED_LOSS_NO",
            "M00158_inferred_loss": "ANCESTRAL_ABSENT",
            "M00160_inferred_loss": "ANCESTRAL_ABSENT"
        }
        self.assertEqual(determine_overall_atp_state(m3), "DUAL_V_AND_F_TYPE_PRESENT")

        # Case 4: Total loss (V-type lost, others absent)
        m4 = {
            "M00159_inferred_loss": "INFERRED_LOSS_HIGH",
            "M00157_inferred_loss": "ANCESTRAL_ABSENT",
            "M00158_inferred_loss": "ANCESTRAL_ABSENT",
            "M00160_inferred_loss": "ANCESTRAL_ABSENT"
        }
        self.assertEqual(determine_overall_atp_state(m4), "TOTAL_ATP_SYNTHASE_LOSS")

        # Case 5: All ancestral absent
        m5 = {
            "M00159_inferred_loss": "ANCESTRAL_ABSENT",
            "M00157_inferred_loss": "ANCESTRAL_ABSENT",
            "M00158_inferred_loss": "ANCESTRAL_ABSENT",
            "M00160_inferred_loss": "ANCESTRAL_ABSENT"
        }
        self.assertEqual(determine_overall_atp_state(m5), "COMPLETE_ATP_SYNTHASE_LOSS")

    def test_compile_matrices(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_id = "GCA_964582085.1"
            
            # Create mock Stage 4 TSVs for M00159 and M00157
            s4_m159 = os.path.join(tmpdir, "s4_m159.tsv")
            df159 = pd.DataFrame([{
                "target_id": sample_id,
                "module_id": "M00159",
                "module_name": "V/A-type ATPase, prokaryotes",
                "GTDB_taxonomy": "d__Archaea;p__Thermoproteota",
                "completeness": 95.0,
                "contamination": 1.2,
                "evolutionary_state": "INFERRED_LOSS_HIGH",
                "confidence": "HIGH",
                "cluster_state": "FULLY_ABSENT",
                "structural_event_type": "CLEAN_DELETION",
                "nearest_positive_ref": "GCA_041479225.1",
                "target_contig_id": "CAZOPV010000010.1",
                "target_contig": "CAZOPV010000010.1"
            }])
            df159.to_csv(s4_m159, sep="\t", index=False)

            s4_m157 = os.path.join(tmpdir, "s4_m157.tsv")
            df157 = pd.DataFrame([{
                "target_id": sample_id,
                "module_id": "M00157",
                "module_name": "F-type ATPase, prokaryotes and chloroplasts",
                "GTDB_taxonomy": "d__Archaea;p__Thermoproteota",
                "completeness": 95.0,
                "contamination": 1.2,
                "evolutionary_state": "ANCESTRAL_ABSENT",
                "confidence": "HIGH",
                "cluster_state": "ANCESTRAL_ABSENT",
                "structural_event_type": "ANCESTRAL_ABSENT",
                "nearest_positive_ref": "NONE",
                "target_contig_id": "CAZOPV010000010.1",
                "target_contig": "CAZOPV010000010.1"
            }])
            df157.to_csv(s4_m157, sep="\t", index=False)

            out_stacked = os.path.join(tmpdir, "stacked.tsv")
            out_matrix = os.path.join(tmpdir, "matrix.tsv")

            compile_matrices([s4_m159, s4_m157], out_stacked, out_matrix)

            self.assertTrue(os.path.exists(out_stacked))
            self.assertTrue(os.path.exists(out_matrix))

            # Verify stacked
            st_df = pd.read_csv(out_stacked, sep="\t")
            self.assertEqual(len(st_df), 2)
            self.assertIn("M00159", st_df["module_id"].values)
            self.assertIn("M00157", st_df["module_id"].values)

            # Verify matrix
            mat_df = pd.read_csv(out_matrix, sep="\t")
            self.assertEqual(len(mat_df), 1)
            row = mat_df.iloc[0]
            self.assertEqual(row["target_id"], sample_id)
            self.assertEqual(row["M00159_inferred_loss"], "INFERRED_LOSS_HIGH")
            self.assertEqual(row["M00157_inferred_loss"], "ANCESTRAL_ABSENT")
            self.assertEqual(row["overall_atp_synthase_state"], "TOTAL_ATP_SYNTHASE_LOSS")

    def test_rescue_graceful_on_empty_operon(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            from scripts.rescue_and_classify_gene_state import parse_args, main
            # Create empty query faa
            empty_faa = os.path.join(tmpdir, "empty_query.faa")
            with open(empty_faa, "w") as f:
                pass
            target_fa = os.path.join(tmpdir, "target.fa")
            with open(target_fa, "w") as f:
                f.write(">ctg1\nATGCATGCATGC\n")
            target_faa = os.path.join(tmpdir, "target.faa")
            with open(target_faa, "w") as f:
                f.write(">prot1\nMHAC\n")

            out_gene = os.path.join(tmpdir, "gene.tsv")
            out_cluster = os.path.join(tmpdir, "cluster.tsv")

            test_args = [
                "rescue_and_classify_gene_state.py",
                "--target_name", "TEST_TARGET",
                "--target_locus_fa", target_fa,
                "--target_faa", target_faa,
                "--query_operon_faa", empty_faa,
                "--out_gene_tsv", out_gene,
                "--out_cluster_tsv", out_cluster
            ]
            old_argv = sys.argv
            sys.argv = test_args
            try:
                main()
            finally:
                sys.argv = old_argv

            self.assertTrue(os.path.exists(out_cluster))
            c_df = pd.read_csv(out_cluster, sep="\t")
            self.assertEqual(c_df.iloc[0]["cluster_state"], "UNRESOLVED_NO_QUERY")

if __name__ == "__main__":
    unittest.main()
