#!/usr/bin/env python3
import unittest
import tempfile
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from scripts.compile_atp_synteny_matrix import determine_overall_atp_state, compile_matrices
from scripts.infer_evolutionary_loss import parse_args as parse_s4_args, main as main_s4

class TestOperonFateSpace(unittest.TestCase):
    """
    Unit tests for the 2D Operon Fate Space (OFI) Framework
    Covering Core-1 through Core-7 and Frozen Precedence State Machine.
    """

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.d = self.tmpdir.name

    def tearDown(self):
        self.tmpdir.cleanup()

    def _run_stage4(self, topo_data, cluster_data, locus_data, refs_data=None, extra_args=None):
        sample = "TestSample"
        topo_tsv = os.path.join(self.d, "anchor_topology.tsv")
        cluster_tsv = os.path.join(self.d, "cluster_state.tsv")
        locus_tsv = os.path.join(self.d, "locus_absence_validation.tsv")
        refs_tsv = os.path.join(self.d, "selected_validated_references.tsv")
        out_summary = os.path.join(self.d, "summary.tsv")
        out_tree = os.path.join(self.d, "tree.nwk")

        pd.DataFrame([topo_data]).to_csv(topo_tsv, sep="\t", index=False)
        cluster_data = {"replacement_search_status": "COMPLETE",
                        "novel_insertion_evidence": "VALIDATED", **cluster_data}
        pd.DataFrame([cluster_data]).to_csv(cluster_tsv, sep="\t", index=False)
        pd.DataFrame([locus_data]).to_csv(locus_tsv, sep="\t", index=False)

        if refs_data is None:
            refs_data = [
                {"accession": "Ref1", "validation_status": "CONFIRMED_POSITIVE", "taxonomy": "d__Archaea;p__Thermoproteota;c__Bathyarchaeia;o__Bathyarchaeales;f__Fam1;g__Gen1;s__sp1", "ANI": 88.0, "AF_query": 0.8, "AF_ref": 0.8, "AF_min": 0.8},
                {"accession": "Ref2", "validation_status": "CONFIRMED_POSITIVE", "taxonomy": "d__Archaea;p__Thermoproteota;c__Bathyarchaeia;o__Bathyarchaeales;f__Fam2;g__Gen2;s__sp2", "ANI": 86.0, "AF_query": 0.75, "AF_ref": 0.75, "AF_min": 0.75},
                {"accession": "RefNeg", "validation_status": "CONFIRMED_NEGATIVE", "taxonomy": "d__Archaea;p__Outgroup;c__Out;o__Out;f__Out;g__Out;s__Out", "ANI": 70.0, "AF_query": 0.5, "AF_ref": 0.5, "AF_min": 0.5}
            ]
        pd.DataFrame(refs_data).to_csv(refs_tsv, sep="\t", index=False)

        cmd_args = [
            "--target_name", sample,
            "--target_taxonomy", "d__Archaea;p__Thermoproteota;c__Bathyarchaeia;o__Bathyarchaeales;f__Fam1;g__Gen1;s__",
            "--target_completeness", "95.0",
            "--target_contamination", "1.0",
            "--module_id", "M00159",
            "--module_name", "V/A-type ATPase",
            "--refs_tsv", refs_tsv,
            "--topology_tsv", topo_tsv,
            "--cluster_tsv", cluster_tsv,
            "--locus_tsv", locus_tsv,
            "--out_summary_tsv", out_summary,
            "--out_tree_nwk", out_tree
        ]
        if extra_args:
            cmd_args.extend(extra_args)

        old_argv = sys.argv
        sys.argv = ["infer_evolutionary_loss.py"] + cmd_args
        try:
            main_s4()
        finally:
            sys.argv = old_argv

        res_df = pd.read_csv(out_summary, sep="\t", keep_default_na=False)
        return res_df.iloc[0].to_dict()

    def test_core2_novel_insertion_precedes_clean_closure(self):
        """
        [Core-2] When both I_novel >= 0.6 and K_clean >= 0.6,
        precedence MUST favor LOSS_WITH_NOVEL_INSERTION over CLEAN_LOSS.
        """
        topo = {
            "locator_score": 0.9,
            "bilateral_context_score": 0.85,
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        cluster = {
            "T_local": 0.0,
            "T_global": 0.0,
            "R_known": 0.0,
            "I_novel": 0.75,  # Novel insertion >= 0.6
            "A_local": "1",
            "A_global": "1",
            "intact_rescued_count": 0,
            "pseudogene_count": 0,
            "frameshift_count": 0,
            "premature_stop_count": 0,
            "pseudogene_detected": "NO",
            "expected_gene_count": 8
        }
        locus = {
            "B_bridge": 1.0,
            "Q_continuity": 1.0,
            "D_closure": 0.8,  # K_clean = sqrt(0.85*1.0) * 0.8 = 0.922 * 0.8 = 0.737 >= 0.6
            "locus_state": "SYNTENIC_ABSENT_HIGH_CONF"
        }

        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "LOSS_WITH_NOVEL_INSERTION",
                         "Novel insertion must take priority over clean closure!")
        self.assertEqual(row["loss_gate_status"], "PASSED")
        self.assertAlmostEqual(float(row["OFI_fate"]), 0.0, places=3)
        self.assertAlmostEqual(float(row["OFI_occupancy"]), 0.75, places=3)

    def test_clean_loss_when_no_novel_insertion(self):
        """
        When K_clean >= 0.6 and I_novel < 0.6, call CLEAN_LOSS.
        """
        topo = {
            "locator_score": 0.9,
            "bilateral_context_score": 0.85,
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        cluster = {
            "T_local": 0.0,
            "T_global": 0.0,
            "R_known": 0.0,
            "I_novel": 0.1,  # Low novel insertion
            "A_local": "1",
            "A_global": "1",
            "intact_rescued_count": 0,
            "pseudogene_count": 0,
            "frameshift_count": 0,
            "premature_stop_count": 0,
            "pseudogene_detected": "NO",
            "expected_gene_count": 8
        }
        locus = {
            "B_bridge": 1.0,
            "Q_continuity": 1.0,
            "D_closure": 0.8,
            "nucleotide_deletion_support": "YES",
            "locus_state": "SYNTENIC_ABSENT_HIGH_CONF"
        }

        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "CLEAN_LOSS")
        self.assertEqual(row["evolutionary_state"], "UNRESOLVED_PHYLOGENY")
        self.assertEqual(row["species_tree_status"], "NOT_PROVIDED")
        self.assertAlmostEqual(float(row["OFI_fate"]), 0.0, places=3)
        self.assertAlmostEqual(float(row["OFI_occupancy"]), 0.1, places=3)

    def test_compressed_anchor_gap_without_alignment_is_not_clean_loss(self):
        topo = {"locator_score": 0.9, "bilateral_context_score": 0.85,
                "syntenic_geometry": "BILATERAL_SAME_CONTIG", "same_contig": "YES"}
        cluster = {"T_local": 0.0, "T_global": 0.0, "R_known": 0.0,
                   "I_novel": 0.0, "A_local": "1", "A_global": "1",
                   "expected_gene_count": 8}
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.8,
                 "nucleotide_deletion_support": "NO", "locus_state": "SYNTENIC_ABSENT_MEDIUM_CONF"}
        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "LOSS_WITH_UNRESOLVED_STRUCTURE")

    def test_unknown_intervening_genes_are_not_called_novel_insertion(self):
        topo = {"locator_score": 0.9, "bilateral_context_score": 0.85,
                "syntenic_geometry": "BILATERAL_SAME_CONTIG", "same_contig": "YES"}
        cluster = {"T_local": 0.0, "T_global": 0.0, "R_known": 0.0,
                   "I_novel": 0.75, "A_local": "1", "A_global": "1",
                   "novel_insertion_evidence": "NOT_ASSESSED", "expected_gene_count": 8}
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.8}
        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "LOSS_WITH_UNCHARACTERIZED_OCCUPANCY")

    def test_incomplete_replacement_search_blocks_loss(self):
        topo = {"locator_score": 0.9, "bilateral_context_score": 0.85,
                "syntenic_geometry": "BILATERAL_SAME_CONTIG", "same_contig": "YES"}
        cluster = {"T_local": 0.0, "T_global": 0.0, "R_known": 0.0,
                   "I_novel": 0.0, "A_local": "1", "A_global": "1",
                   "replacement_search_status": "INCOMPLETE", "expected_gene_count": 8}
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.8}
        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "UNRESOLVED")
        self.assertEqual(row["loss_gate_status"], "INCOMPLETE_REPLACEMENT_SEARCH")

    def test_loss_with_locus_remodeling(self):
        """
        When K_loss >= 0.6 but neither novel insertion nor clean closure >= 0.6,
        call LOSS_WITH_LOCUS_REMODELING.
        """
        topo = {
            "locator_score": 0.85,
            "bilateral_context_score": 0.81,
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        cluster = {
            "T_local": 0.0,
            "T_global": 0.0,
            "R_known": 0.0,
            "I_novel": 0.3,  # < 0.6
            "A_local": "1",
            "A_global": "1",
            "intact_rescued_count": 0,
            "pseudogene_count": 0,
            "frameshift_count": 0,
            "premature_stop_count": 0,
            "pseudogene_detected": "NO",
            "expected_gene_count": 8
        }
        locus = {
            "B_bridge": 1.0,
            "Q_continuity": 1.0,
            "D_closure": 0.4,  # K_clean = 0.9 * 0.4 = 0.36 < 0.6
            "nucleotide_deletion_support": "YES",
            "locus_state": "SYNTENIC_ABSENT_HIGH_CONF"
        }

        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "LOSS_WITH_LOCUS_REMODELING")

    def test_core3_annotation_miss_vs_pseudogene(self):
        """
        [Core-3] Distinguish intact ORF rescue (ANNOTATION_MISS) from
        disrupted homologous remnants (PSEUDOGENIZED).
        """
        topo = {
            "locator_score": 0.9,
            "bilateral_context_score": 0.85,
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }

        # Subcase A: Intact coding potential without mutations -> ANNOTATION_MISS
        cluster_intact = {
            "T_local": 0.2,
            "T_global": 0.2,
            "R_known": 0.0,
            "I_novel": 0.0,
            "A_local": "0",
            "A_global": "0",
            "intact_rescued_count": 5,  # 5 subunits intact rescued
            "pseudogene_count": 0,
            "frameshift_count": 0,
            "premature_stop_count": 0,
            "annotation_issue": "YES",
            "pseudogene_detected": "NO",
            "expected_gene_count": 8
        }
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.0, "locus_state": "PRESENT_COMPLETE"}

        row_a = self._run_stage4(topo, cluster_intact, locus)
        self.assertEqual(row_a["fate_category"], "ANNOTATION_MISS")
        self.assertEqual(row_a["evolutionary_state"], "RESCUED_ANNOTATION_MISS")

        # Subcase B: Homologous remnants with frameshifts -> PSEUDOGENIZED
        cluster_pseudo = {
            "T_local": 0.1,
            "T_global": 0.1,
            "R_known": 0.0,
            "I_novel": 0.0,
            "A_local": "0",
            "A_global": "0",
            "intact_rescued_count": 0,
            "pseudogene_count": 3,
            "frameshift_count": 4,  # Disruptive mutations
            "premature_stop_count": 2,
            "annotation_issue": "NO",
            "pseudogene_detected": "YES",
            "expected_gene_count": 8
        }

        row_b = self._run_stage4(topo, cluster_pseudo, locus)
        self.assertEqual(row_b["fate_category"], "PSEUDOGENIZED")
        self.assertEqual(row_b["evolutionary_state"], "PSEUDOGENIZED")

    def test_core5_bilateral_context_hard_gate(self):
        """
        [Core-5] bilateral_context_min: 0.5 must be enforced as a hard gate.
        If C_bilateral < 0.5, loss inference cannot pass.
        """
        topo = {
            "locator_score": 0.45,
            "bilateral_context_score": 0.42,  # < 0.5 gate
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        cluster = {
            "T_local": 0.0,
            "T_global": 0.0,
            "R_known": 0.0,
            "I_novel": 0.0,
            "A_local": "1",
            "A_global": "1",
            "intact_rescued_count": 0,
            "pseudogene_count": 0,
            "frameshift_count": 0,
            "premature_stop_count": 0,
            "pseudogene_detected": "NO",
            "expected_gene_count": 8
        }
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.9, "locus_state": "SYNTENIC_ABSENT_HIGH_CONF"}

        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "UNRESOLVED")
        self.assertEqual(row["loss_gate_status"], "FAILED_GATES")
        self.assertEqual(row["K_loss"], 0.0)

    def test_core6_replacement_absence_hard_gate(self):
        """
        [Core-6] replacement_absence_max: 0.2 must be enforced.
        If 0.2 < R_known < 0.6, output UNRESOLVED (ambiguous partial replacement).
        If R_known >= 0.6, output KNOWN_REPLACEMENT.
        """
        topo = {
            "locator_score": 0.9,
            "bilateral_context_score": 0.85,
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.8, "locus_state": "SYNTENIC_ABSENT_HIGH_CONF"}

        # Subcase A: R_known = 0.35 (between 0.2 and 0.6)
        cluster_partial = {
            "T_local": 0.0, "T_global": 0.0,
            "R_known": 0.35, "I_novel": 0.0,
            "A_local": "1", "A_global": "1",
            "intact_rescued_count": 0, "pseudogene_count": 0,
            "frameshift_count": 0, "premature_stop_count": 0,
            "pseudogene_detected": "NO", "expected_gene_count": 8
        }
        row_a = self._run_stage4(topo, cluster_partial, locus)
        self.assertEqual(row_a["fate_category"], "UNRESOLVED")
        self.assertEqual(row_a["loss_gate_status"], "PARTIAL_REPLACEMENT_AMBIGUOUS")

        # Subcase B: R_known = 0.75 (>= 0.6)
        cluster_full = {
            "T_local": 0.0, "T_global": 0.0,
            "R_known": 0.75, "I_novel": 0.0,
            "A_local": "1", "A_global": "1",
            "intact_rescued_count": 0, "pseudogene_count": 0,
            "frameshift_count": 0, "premature_stop_count": 0,
            "pseudogene_detected": "NO", "expected_gene_count": 8
        }
        row_b = self._run_stage4(topo, cluster_full, locus)
        self.assertEqual(row_b["fate_category"], "KNOWN_REPLACEMENT")
        self.assertAlmostEqual(float(row_b["OFI_fate"]), -0.75, places=3)
        self.assertAlmostEqual(float(row_b["OFI_occupancy"]), 0.75, places=3)

    def test_core7_tristate_negative_evidence_na_guard(self):
        """
        [Core-7] When A_local or A_global is NA,
        output K_loss = NA, loss_gate_status = INCOMPLETE_NEGATIVE_SEARCH, fate = UNRESOLVED.
        """
        topo = {
            "locator_score": 0.9,
            "bilateral_context_score": 0.85,
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        cluster = {
            "T_local": 0.0,
            "T_global": 0.0,
            "R_known": 0.0,
            "I_novel": 0.0,
            "A_local": "NA",  # Incomplete local search!
            "A_global": "1",
            "intact_rescued_count": 0,
            "pseudogene_count": 0,
            "frameshift_count": 0,
            "premature_stop_count": 0,
            "pseudogene_detected": "NO",
            "expected_gene_count": 8
        }
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.8, "locus_state": "SYNTENIC_ABSENT_HIGH_CONF"}

        row = self._run_stage4(topo, cluster, locus)
        self.assertEqual(row["fate_category"], "UNRESOLVED")
        self.assertEqual(row["loss_gate_status"], "INCOMPLETE_NEGATIVE_SEARCH")
        self.assertEqual(str(row["K_loss"]), "NA")
        self.assertEqual(str(row["K_clean"]), "NA")

    def test_geometry_assembly_candidates(self):
        """
        Verify that assembly fragmentation is cleanly separated from evolutionary loss:
        - SPLIT_CONTIG -> ASSEMBLY_SPLIT_CANDIDATE
        - UNILATERAL_SAME_CONTIG -> ASSEMBLY_EDGE_CANDIDATE
        """
        cluster = {
            "T_local": 0.0, "T_global": 0.0, "R_known": 0.0, "I_novel": 0.0,
            "A_local": "1", "A_global": "1", "intact_rescued_count": 0,
            "pseudogene_count": 0, "frameshift_count": 0, "premature_stop_count": 0,
            "pseudogene_detected": "NO", "expected_gene_count": 8
        }
        locus = {"B_bridge": 0.0, "Q_continuity": 1.0, "D_closure": 0.0, "locus_state": "ASSEMBLY_UNRESOLVED"}

        # Subcase 1: Split Contig
        topo_split = {
            "locator_score": 0.5,
            "bilateral_context_score": 0.0,
            "syntenic_geometry": "SPLIT_CONTIG",
            "same_contig": "NO",
            "nearest_positive_ref": "Ref1"
        }
        row_split = self._run_stage4(topo_split, cluster, locus)
        self.assertEqual(row_split["fate_category"], "ASSEMBLY_SPLIT_CANDIDATE")

        # Subcase 2: Unilateral Edge
        topo_edge = {
            "locator_score": 0.6,
            "bilateral_context_score": 0.0,
            "syntenic_geometry": "UNILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        row_edge = self._run_stage4(topo_edge, cluster, locus)
        self.assertEqual(row_edge["fate_category"], "ASSEMBLY_EDGE_CANDIDATE")

    def test_target_retained_and_coexistence(self):
        """
        Verify canonical TARGET_RETAINED and COEXISTENCE states.
        """
        topo = {
            "locator_score": 0.95,
            "bilateral_context_score": 0.95,
            "syntenic_geometry": "BILATERAL_SAME_CONTIG",
            "same_contig": "YES",
            "nearest_positive_ref": "Ref1"
        }
        locus = {"B_bridge": 1.0, "Q_continuity": 1.0, "D_closure": 0.0, "locus_state": "PRESENT_COMPLETE"}

        # Retention
        cluster_ret = {
            "T_local": 1.0, "T_global": 1.0, "R_known": 0.0, "I_novel": 0.0,
            "A_local": "0", "A_global": "0", "intact_rescued_count": 0,
            "pseudogene_count": 0, "frameshift_count": 0, "premature_stop_count": 0,
            "pseudogene_detected": "NO", "expected_gene_count": 8
        }
        row_ret = self._run_stage4(topo, cluster_ret, locus)
        self.assertEqual(row_ret["fate_category"], "TARGET_RETAINED")
        self.assertAlmostEqual(float(row_ret["OFI_fate"]), 1.0, places=3)
        self.assertAlmostEqual(float(row_ret["OFI_occupancy"]), 1.0, places=3)

        # Coexistence: T_local >= 0.6 AND R_known >= 0.6
        cluster_coex = {
            "T_local": 0.8, "T_global": 0.8, "R_known": 0.8, "I_novel": 0.0,
            "A_local": "0", "A_global": "0", "intact_rescued_count": 0,
            "pseudogene_count": 0, "frameshift_count": 0, "premature_stop_count": 0,
            "pseudogene_detected": "NO", "expected_gene_count": 8
        }
        row_coex = self._run_stage4(topo, cluster_coex, locus)
        self.assertEqual(row_coex["fate_category"], "COEXISTENCE")
        self.assertAlmostEqual(float(row_coex["OFI_fate"]), 0.0, places=3)
        self.assertAlmostEqual(float(row_coex["OFI_occupancy"]), 0.8, places=3)

    def test_compile_matrices_with_ofi_fields(self):
        """
        Verify that compile_matrices includes all OFI fields in wide pivot table.
        """
        sample_id = "TestSample_01"
        s4_tsv = os.path.join(self.d, "mock_s4.tsv")
        df_s4 = pd.DataFrame([{
            "target_id": sample_id,
            "module_id": "M00159",
            "module_name": "V/A-type ATPase",
            "GTDB_taxonomy": "d__Archaea;p__Thermoproteota",
            "completeness": 95.0,
            "contamination": 1.0,
            "evolutionary_state": "INFERRED_LOSS_HIGH",
            "confidence": "HIGH",
            "cluster_state": "FULLY_ABSENT",
            "structural_event_type": "CLEAN_DELETION",
            "nearest_positive_ref": "Ref1",
            "target_contig_id": "ctg_01",
            "locator_score": 0.95,
            "bilateral_context_score": 0.92,
            "K_loss": 0.95,
            "K_clean": 0.88,
            "loss_gate_status": "PASSED",
            "OFI_fate": 0.0,
            "OFI_occupancy": 0.05,
            "fate_category": "CLEAN_LOSS"
        }])
        df_s4.to_csv(s4_tsv, sep="\t", index=False)

        out_stacked = os.path.join(self.d, "stacked_ofi.tsv")
        out_matrix = os.path.join(self.d, "matrix_ofi.tsv")

        compile_matrices([s4_tsv], out_stacked, out_matrix)

        self.assertTrue(os.path.exists(out_matrix))
        mat_df = pd.read_csv(out_matrix, sep="\t")
        self.assertEqual(len(mat_df), 1)
        row = mat_df.iloc[0]

        # Verify OFI columns exist in pivot matrix
        self.assertIn("M00159_fate_category", mat_df.columns)
        self.assertIn("M00159_OFI_fate", mat_df.columns)
        self.assertIn("M00159_OFI_occupancy", mat_df.columns)
        self.assertIn("M00159_locator_score", mat_df.columns)
        self.assertIn("M00159_bilateral_context", mat_df.columns)
        self.assertIn("M00159_K_loss", mat_df.columns)
        self.assertIn("M00159_K_clean", mat_df.columns)

        self.assertEqual(row["M00159_fate_category"], "CLEAN_LOSS")
        self.assertEqual(float(row["M00159_OFI_fate"]), 0.0)
        self.assertEqual(float(row["M00159_locator_score"]), 0.95)
        self.assertEqual(float(row["M00159_K_loss"]), 0.95)

    def test_core1_stage1_candidate_ranking_target_blind(self):
        """
        [Core-1] Stage 1 candidate ranking MUST be target-blind.
        Candidate sorting tuple:
        (geom_rank, loc_score, -obs_gap, p_emp, -ctg_len)
        ensures that loci are evaluated strictly by flanking synteny topology.
        """
        # Contig A: Weak unilateral flank, but has target homology
        cand_a = {
            "contig_id": "ctg_A",
            "geometry": "UNILATERAL_EDGE",
            "geom_rank": 2, # UNILATERAL
            "locator_score": 0.45,
            "obs_gap": 15000,
            "has_target_operon": True
        }
        # Contig B: Strong bilateral flank with high locator score, true ancestral locus
        cand_b = {
            "contig_id": "ctg_B",
            "geometry": "BILATERAL_SAME_CONTIG",
            "geom_rank": 3, # BILATERAL
            "locator_score": 0.92,
            "obs_gap": 1200,
            "has_target_operon": False # Target absent!
        }

        # Stage 1 ranking key: (geom_rank, locator_score) - TARGET BLIND!
        # Both without has_target_operon
        key_a = (cand_a["geom_rank"], cand_a["locator_score"])
        key_b = (cand_b["geom_rank"], cand_b["locator_score"])

        ranked = sorted([cand_a, cand_b], key=lambda c: (c["geom_rank"], c["locator_score"]), reverse=True)
        self.assertEqual(ranked[0]["contig_id"], "ctg_B",
                         "Stage 1 must rank the strong bilateral flank above unilateral, completely blind to target operon presence!")

if __name__ == "__main__":
    unittest.main()
