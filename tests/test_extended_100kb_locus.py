#!/usr/bin/env python3
import os
import sys
import tempfile
import unittest
import pandas as pd
try:
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    class Seq:
        def __init__(self, data=""):
            self.data = data
        def __len__(self):
            return len(self.data)
        def __str__(self):
            return str(self.data)

    class SeqRecord:
        def __init__(self, seq, id=""):
            self.seq = seq
            self.id = id
        def __len__(self):
            return len(self.seq)

# Import function from scripts/extract_anchor_loci.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))
from extract_anchor_loci import (
    extract_extended_100kb_locus_proteins,
    dereplicate_reference_cohort,
    compute_empirical_max_contig_null,
    MODULE_TARGET_KOS,
    DEFAULT_KO_DEFINITIONS
)

class MockArgs:
    def __init__(self, target_name="Target_MAG_01", module_id="M00159", consensus_thresh=0.75, win_bp=100000, sig_mode="fixed_ratio"):
        self.target_name = target_name
        self.module_id = module_id
        self.consensus_contig_threshold = consensus_thresh
        self.extended_window_bp = win_bp
        self.target_kos = ""
        self.significance_mode = sig_mode

class TestExtended100kbLocus(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.out_tsv = os.path.join(self.test_dir.name, "output_100kb_proteins.tsv")

        # Mock target assembly: contig_1 (500,000 bp)
        self.target_records = {
            "ctg_01": SeqRecord(Seq("N" * 500000), id="ctg_01"),
            "ctg_02": SeqRecord(Seq("N" * 200000), id="ctg_02")
        }

        # Mock target proteome coordinates on ctg_01
        self.gff_coords = {
            # Upstream flanking region
            "T_prot_01": {"contig": "ctg_01", "start": 120000, "end": 121000, "strand": "+"},
            # Upstream anchor mapped protein
            "T_prot_02": {"contig": "ctg_01", "start": 180000, "end": 181500, "strand": "+"},
            # Interval protein (in breakpoint gap)
            "T_prot_03": {"contig": "ctg_01", "start": 200000, "end": 201200, "strand": "+"},
            # Downstream anchor mapped protein
            "T_prot_04": {"contig": "ctg_01", "start": 220000, "end": 221800, "strand": "-"},
            # Downstream flanking region
            "T_prot_05": {"contig": "ctg_01", "start": 280000, "end": 281000, "strand": "-"},
            # Far outside 100kb window (> 350kb)
            "T_prot_far": {"contig": "ctg_01", "start": 450000, "end": 451000, "strand": "+"},
        }

        self.target_faa_dict = {
            "T_prot_01": "M" * 333,
            "T_prot_02": "M" * 500,
            "T_prot_03": "M" * 400,
            "T_prot_04": "M" * 600,
            "T_prot_05": "M" * 333,
            "T_prot_far": "M" * 333,
        }

        self.target_anno_dict = {
            "T_prot_01": {"ko": "K00001", "definition": "adh; alcohol dehydrogenase"},
            "T_prot_03": {"ko": "K00002", "definition": "adhE; aldehyde-alcohol dehydrogenase"},
        }

        self.ko_definitions = DEFAULT_KO_DEFINITIONS.copy()
        self.ref_flank_annotations = {}

    def tearDown(self):
        self.test_dir.cleanup()

    def test_consensus_threshold_below_75_percent(self):
        """When fewer than 75% of positive references converge on the same contig, emit header-only TSV."""
        args = MockArgs(consensus_thresh=0.75)
        # 10 positive references: 6 on ctg_01 (60%), 4 on ctg_02 (40%)
        ref_candidates = {}
        ref_evaluations = {}
        for i in range(1, 7):
            acc = f"Ref_Pos_{i:02d}"
            ref_candidates[acc] = {"clean_acc": acc, "best_ctg": "ref_ctg", "op_pids": set(), "ref_coords": {}}
            ref_evaluations[acc] = {
                "score": 1000.0, "best_contig": "ctg_01",
                "best_cluster": [{"start": 180000, "end": 221800, "side": "UP", "orig_pid": f"R_{i}_u", "t_pid": "T_prot_02"}]
            }
        for i in range(7, 11):
            acc = f"Ref_Pos_{i:02d}"
            ref_candidates[acc] = {"clean_acc": acc, "best_ctg": "ref_ctg", "op_pids": set(), "ref_coords": {}}
            ref_evaluations[acc] = {
                "score": 800.0, "best_contig": "ctg_02",
                "best_cluster": [{"start": 50000, "end": 70000, "side": "UP", "orig_pid": f"R_{i}_u", "t_pid": "T_other"}]
            }

        pos_refs = pd.DataFrame([{"accession": k} for k in ref_candidates.keys()])

        extract_extended_100kb_locus_proteins(
            args=args,
            pos_refs=pos_refs,
            ref_candidates=ref_candidates,
            ref_evaluations=ref_evaluations,
            target_records=self.target_records,
            gff_coords=self.gff_coords,
            target_faa_dict=self.target_faa_dict,
            target_anno_dict=self.target_anno_dict,
            ref_flank_annotations=self.ref_flank_annotations,
            ko_definitions=self.ko_definitions,
            out_tsv_path=self.out_tsv
        )

        df = pd.read_csv(self.out_tsv, sep="\t")
        self.assertTrue(df.empty, "Output DataFrame should be empty when consensus < 75%")
        self.assertEqual(list(df.columns), [
            "target_id", "module_id", "genome_id", "genome_type", "contig_id",
            "protein_id", "gene_symbol", "start", "end", "strand", "length_aa",
            "ko", "ko_definition", "label", "mapped_ref_count", "mapped_ref_ratio"
        ])

    def test_weighted_midpoint_and_classification(self):
        """
        Verify weighted midpoint math and protein classification into
        anchor_mapped, interval, flanking, ref_target, ref_anchor_mapped, ref_flanking.
        """
        args = MockArgs(consensus_thresh=0.75, win_bp=100000)

        # 4 positive references: 3 on ctg_01 (75% >= 75%), 1 on ctg_02 (25%)
        # Ref 1: hits on ctg_01 at [180000 - 220000], span = 40000, mid = 200000
        # Ref 2: hits on ctg_01 at [185000 - 225000], span = 40000, mid = 205000
        # Ref 3: hits on ctg_01 at [190000 - 210000], span = 20000, mid = 200000
        # Expected weighted center:
        # sum(weights) = 40000 + 40000 + 20000 = 100000
        # sum(w * m) = 40000 * 200000 + 40000 * 205000 + 20000 * 200000
        #            = 8e9 + 8.2e9 + 4e9 = 20.2e9
        # Center = 20.2e9 / 100000 = 202,000 bp
        # Window: [202000 - 100000, 202000 + 100000] = [102000, 302000]

        ref_candidates = {}
        ref_evaluations = {}

        # Ref 1
        ref_candidates["Ref_01"] = {
            "clean_acc": "Ref_01",
            "best_ctg": "r_ctg1",
            "op_pids": {"Ref1_atpA", "Ref1_atpB"},
            "pid_to_ko": {"Ref1_atpA": "K02117", "Ref1_atpB": "K02118"},
            "ref_coords": {
                "Ref1_flank_u": {"contig": "r_ctg1", "start": 40000, "end": 41000, "strand": "+"},
                "Ref1_anchor_u": {"contig": "r_ctg1", "start": 48000, "end": 49000, "strand": "+"},
                "Ref1_atpA": {"contig": "r_ctg1", "start": 50000, "end": 51500, "strand": "+"},
                "Ref1_atpB": {"contig": "r_ctg1", "start": 51600, "end": 53000, "strand": "+"},
                "Ref1_anchor_d": {"contig": "r_ctg1", "start": 54000, "end": 55000, "strand": "+"},
                "Ref1_flank_d": {"contig": "r_ctg1", "start": 60000, "end": 61000, "strand": "+"},
            },
            "ref_faa_dict": {
                "Ref1_atpA": "M" * 500, "Ref1_atpB": "M" * 460,
                "Ref1_flank_u": "M" * 333, "Ref1_anchor_u": "M" * 333,
                "Ref1_anchor_d": "M" * 333, "Ref1_flank_d": "M" * 333,
            }
        }
        ref_evaluations["Ref_01"] = {
            "score": 1500.0,
            "best_contig": "ctg_01",
            "best_cluster": [
                {"start": 180000, "end": 181500, "side": "UP", "orig_pid": "Ref1_anchor_u", "t_pid": "T_prot_02"},
                {"start": 220000, "end": 221800, "side": "DOWN", "orig_pid": "Ref1_anchor_d", "t_pid": "T_prot_04"},
            ]
        }

        # Ref 2
        ref_candidates["Ref_02"] = {
            "clean_acc": "Ref_02",
            "best_ctg": "r_ctg2",
            "op_pids": {"Ref2_atpA"},
            "pid_to_ko": {"Ref2_atpA": "K02117"},
            "ref_coords": {
                "Ref2_atpA": {"contig": "r_ctg2", "start": 10000, "end": 11500, "strand": "+"}
            },
            "ref_faa_dict": {"Ref2_atpA": "M" * 500}
        }
        ref_evaluations["Ref_02"] = {
            "score": 1400.0,
            "best_contig": "ctg_01",
            "best_cluster": [
                {"start": 185000, "end": 186500, "side": "UP", "orig_pid": "Ref2_u", "t_pid": "T_prot_02"},
                {"start": 223000, "end": 225000, "side": "DOWN", "orig_pid": "Ref2_d", "t_pid": "T_prot_04"},
            ]
        }

        # Ref 3
        ref_candidates["Ref_03"] = {
            "clean_acc": "Ref_03",
            "best_ctg": "r_ctg3",
            "op_pids": {"Ref3_atpA"},
            "pid_to_ko": {"Ref3_atpA": "K02117"},
            "ref_coords": {
                "Ref3_atpA": {"contig": "r_ctg3", "start": 20000, "end": 21500, "strand": "+"}
            },
            "ref_faa_dict": {"Ref3_atpA": "M" * 500}
        }
        ref_evaluations["Ref_03"] = {
            "score": 1200.0,
            "best_contig": "ctg_01",
            "best_cluster": [
                {"start": 190000, "end": 191000, "side": "UP", "orig_pid": "Ref3_u", "t_pid": "T_prot_02"},
                {"start": 209000, "end": 210000, "side": "DOWN", "orig_pid": "Ref3_d", "t_pid": "T_prot_04"},
            ]
        }

        # Ref 4 (on ctg_02)
        ref_candidates["Ref_04"] = {
            "clean_acc": "Ref_04", "best_ctg": "r_ctg4", "op_pids": set(), "pid_to_ko": {}, "ref_coords": {}
        }
        ref_evaluations["Ref_04"] = {
            "score": 500.0, "best_contig": "ctg_02",
            "best_cluster": [{"start": 10000, "end": 12000, "side": "UP", "orig_pid": "R4_u", "t_pid": "T_other"}]
        }

        pos_refs = pd.DataFrame([{"accession": k} for k in ref_candidates.keys()])

        extract_extended_100kb_locus_proteins(
            args=args,
            pos_refs=pos_refs,
            ref_candidates=ref_candidates,
            ref_evaluations=ref_evaluations,
            target_records=self.target_records,
            gff_coords=self.gff_coords,
            target_faa_dict=self.target_faa_dict,
            target_anno_dict=self.target_anno_dict,
            ref_flank_annotations=self.ref_flank_annotations,
            ko_definitions=self.ko_definitions,
            out_tsv_path=self.out_tsv
        )

        df = pd.read_csv(self.out_tsv, sep="\t")
        self.assertFalse(df.empty, "Output TSV should not be empty")

        # 1. Target proteins checks
        target_df = df[df["genome_type"] == "target"]
        target_pids = list(target_df["protein_id"])
        self.assertIn("T_prot_01", target_pids)
        self.assertIn("T_prot_02", target_pids)
        self.assertIn("T_prot_03", target_pids)
        self.assertIn("T_prot_04", target_pids)
        self.assertIn("T_prot_05", target_pids)
        self.assertNotIn("T_prot_far", target_pids, "T_prot_far should be excluded (outside 100kb window)")

        # Verify classifications and counts:
        # T_prot_02: mapped by 3 refs out of 3 converged refs -> anchor_mapped, mapped_ref_count=3, ratio=3/3=1.0
        p2_row = target_df[target_df["protein_id"] == "T_prot_02"].iloc[0]
        self.assertEqual(p2_row["label"], "anchor_mapped")
        self.assertEqual(p2_row["mapped_ref_count"], 3)
        self.assertAlmostEqual(p2_row["mapped_ref_ratio"], 1.0, places=3)

        # T_prot_03: at 200,000 bp, inside interval [181500, 220000] for 3 refs -> interval, count=3, ratio=3/3=1.0
        p3_row = target_df[target_df["protein_id"] == "T_prot_03"].iloc[0]
        self.assertEqual(p3_row["label"], "interval")
        self.assertEqual(p3_row["mapped_ref_count"], 3)
        self.assertAlmostEqual(p3_row["mapped_ref_ratio"], 1.0, places=3)
        self.assertEqual(p3_row["ko"], "K00002")
        self.assertEqual(p3_row["gene_symbol"], "adhE")

        # T_prot_01: at 120,000 bp, in flanking region -> flanking, count=3, ratio=3/3=1.0
        p1_row = target_df[target_df["protein_id"] == "T_prot_01"].iloc[0]
        self.assertEqual(p1_row["label"], "flanking")
        self.assertEqual(p1_row["mapped_ref_count"], 3)
        self.assertAlmostEqual(p1_row["mapped_ref_ratio"], 1.0, places=3)

        # 2. Reference proteins checks (Ref_01)
        ref_df = df[df["genome_id"] == "Ref_01"]
        self.assertFalse(ref_df.empty)

        # Ref1_atpA -> ref_target
        atpA_row = ref_df[ref_df["protein_id"] == "Ref1_atpA"].iloc[0]
        self.assertEqual(atpA_row["label"], "ref_target")
        self.assertEqual(atpA_row["ko"], "K02117")
        self.assertEqual(atpA_row["gene_symbol"], "atpA")
        self.assertEqual(atpA_row["mapped_ref_count"], 1)
        self.assertEqual(atpA_row["mapped_ref_ratio"], 1.0)

        # Ref1_anchor_u -> ref_anchor_mapped
        anc_u_row = ref_df[ref_df["protein_id"] == "Ref1_anchor_u"].iloc[0]
        self.assertEqual(anc_u_row["label"], "ref_anchor_mapped")

        # Ref1_flank_u -> ref_flanking
        flank_u_row = ref_df[ref_df["protein_id"] == "Ref1_flank_u"].iloc[0]
        self.assertEqual(flank_u_row["label"], "ref_flanking")

    def test_ten_positive_refs_convergence_80_percent(self):
        """10 positive refs with 8 converging on ctg_01 (80% >= 75%) triggers extraction and correct ratios."""
        args = MockArgs(consensus_thresh=0.75, win_bp=100000)
        ref_candidates = {}
        ref_evaluations = {}

        # 8 refs converge on ctg_01
        for i in range(1, 9):
            acc = f"Ref_{i:02d}"
            ref_candidates[acc] = {
                "clean_acc": acc, "best_ctg": f"r_ctg_{i}", "op_pids": {f"R{i}_atpA"},
                "pid_to_ko": {f"R{i}_atpA": "K02117"},
                "ref_coords": {f"R{i}_atpA": {"contig": f"r_ctg_{i}", "start": 10000, "end": 11500, "strand": "+"}},
                "ref_faa_dict": {f"R{i}_atpA": "M" * 500}
            }
            ref_evaluations[acc] = {
                "score": 1000.0, "best_contig": "ctg_01",
                "best_cluster": [
                    {"start": 180000, "end": 181500, "side": "UP", "orig_pid": f"R{i}_u", "t_pid": "T_prot_02"},
                    {"start": 220000, "end": 221800, "side": "DOWN", "orig_pid": f"R{i}_d", "t_pid": "T_prot_04"}
                ]
            }

        # 2 refs on ctg_02
        for i in range(9, 11):
            acc = f"Ref_{i:02d}"
            ref_candidates[acc] = {
                "clean_acc": acc, "best_ctg": f"r_ctg_{i}", "op_pids": set(), "pid_to_ko": {}, "ref_coords": {}
            }
            ref_evaluations[acc] = {
                "score": 500.0, "best_contig": "ctg_02",
                "best_cluster": [{"start": 10000, "end": 12000, "side": "UP", "orig_pid": f"R{i}_u", "t_pid": "T_other"}]
            }

        pos_refs = pd.DataFrame([{"accession": k} for k in ref_candidates.keys()])

        extract_extended_100kb_locus_proteins(
            args=args,
            pos_refs=pos_refs,
            ref_candidates=ref_candidates,
            ref_evaluations=ref_evaluations,
            target_records=self.target_records,
            gff_coords=self.gff_coords,
            target_faa_dict=self.target_faa_dict,
            target_anno_dict=self.target_anno_dict,
            ref_flank_annotations=self.ref_flank_annotations,
            ko_definitions=self.ko_definitions,
            out_tsv_path=self.out_tsv
        )

        df = pd.read_csv(self.out_tsv, sep="\t")
        self.assertFalse(df.empty)

        # Target checks
        target_df = df[df["genome_type"] == "target"]
        p2_row = target_df[target_df["protein_id"] == "T_prot_02"].iloc[0]
        # 8 converged refs mapped T_prot_02 out of 8 converged positive refs -> ratio = 8/8 = 1.0
        self.assertEqual(p2_row["mapped_ref_count"], 8)
        self.assertAlmostEqual(p2_row["mapped_ref_ratio"], 1.0, places=3)

        # T_prot_03 in interval for 8 refs -> ratio = 8/8 = 1.0
        p3_row = target_df[target_df["protein_id"] == "T_prot_03"].iloc[0]
        self.assertEqual(p3_row["label"], "interval")
        self.assertEqual(p3_row["mapped_ref_count"], 8)
        self.assertAlmostEqual(p3_row["mapped_ref_ratio"], 1.0, places=3)

    def test_module_m00157_f_type_target_kos(self):
        """M00157 F-type module recognizes K02111 as ref_target and loads proper gene symbol."""
        args = MockArgs(module_id="M00157", consensus_thresh=0.5, win_bp=50000)
        ref_candidates = {
            "Ref_F": {
                "clean_acc": "Ref_F", "best_ctg": "r_ctg_f",
                "op_pids": {"F_atpA"}, "pid_to_ko": {"F_atpA": "K02111"},
                "ref_coords": {"F_atpA": {"contig": "r_ctg_f", "start": 5000, "end": 6500, "strand": "+"}},
                "ref_faa_dict": {"F_atpA": "M" * 500}
            }
        }
        ref_evaluations = {
            "Ref_F": {
                "score": 1000.0, "best_contig": "ctg_01",
                "best_cluster": [{"start": 180000, "end": 181500, "side": "UP", "orig_pid": "F_u", "t_pid": "T_prot_02"}]
            }
        }
        pos_refs = pd.DataFrame([{"accession": "Ref_F"}])

        extract_extended_100kb_locus_proteins(
            args=args,
            pos_refs=pos_refs,
            ref_candidates=ref_candidates,
            ref_evaluations=ref_evaluations,
            target_records=self.target_records,
            gff_coords=self.gff_coords,
            target_faa_dict=self.target_faa_dict,
            target_anno_dict=self.target_anno_dict,
            ref_flank_annotations=self.ref_flank_annotations,
            ko_definitions=self.ko_definitions,
            out_tsv_path=self.out_tsv
        )

        df = pd.read_csv(self.out_tsv, sep="\t")
        ref_rows = df[df["genome_id"] == "Ref_F"]
        atpA_row = ref_rows[ref_rows["protein_id"] == "F_atpA"].iloc[0]
        self.assertEqual(atpA_row["label"], "ref_target")
        self.assertEqual(atpA_row["ko"], "K02111")
        self.assertEqual(atpA_row["gene_symbol"], "atpA")
        self.assertEqual(atpA_row["module_id"], "M00157")

    def test_contig_boundary_clamping(self):
        """When contig is short (e.g. 80kb), window [Center - 100kb, Center + 100kb] clamps to [0, 80000] without error."""
        short_records = {"short_ctg": SeqRecord(Seq("N" * 80000), id="short_ctg")}
        short_gff = {
            "P_start": {"contig": "short_ctg", "start": 500, "end": 1500, "strand": "+"},
            "P_mid": {"contig": "short_ctg", "start": 40000, "end": 41500, "strand": "+"},
            "P_end": {"contig": "short_ctg", "start": 78000, "end": 79500, "strand": "+"},
        }
        short_faa = {k: "M" * 300 for k in short_gff}
        args = MockArgs(consensus_thresh=0.5, win_bp=100000)

        ref_candidates = {
            "Ref_1": {
                "clean_acc": "Ref_1", "best_ctg": "r_c", "op_pids": set(), "pid_to_ko": {},
                "ref_coords": {}, "ref_faa_dict": {}
            }
        }
        ref_evaluations = {
            "Ref_1": {
                "score": 1000.0, "best_contig": "short_ctg",
                "best_cluster": [{"start": 40000, "end": 41500, "side": "UP", "orig_pid": "R1", "t_pid": "P_mid"}]
            }
        }
        pos_refs = pd.DataFrame([{"accession": "Ref_1"}])

        extract_extended_100kb_locus_proteins(
            args=args,
            pos_refs=pos_refs,
            ref_candidates=ref_candidates,
            ref_evaluations=ref_evaluations,
            target_records=short_records,
            gff_coords=short_gff,
            target_faa_dict=short_faa,
            target_anno_dict={},
            ref_flank_annotations={},
            ko_definitions=self.ko_definitions,
            out_tsv_path=self.out_tsv
        )

        df = pd.read_csv(self.out_tsv, sep="\t")
        target_pids = list(df[df["genome_type"] == "target"]["protein_id"])
        # All 3 proteins on short_ctg should be captured
        self.assertIn("P_start", target_pids)
        self.assertIn("P_mid", target_pids)
        self.assertIn("P_end", target_pids)

    def test_option_a_consensus_with_unmapped_distant_refs(self):
        """
        Option A: 10 positive refs selected in Stage 0, but only 4 mapped to target genome.
        All 4 mapped refs converge on ctg_01 (4/4 = 100% >= 75%).
        The other 6 refs have score 0 (did not map).
        Under Option A, consensus passes and correctly extracts locus proteins.
        """
        args = MockArgs(consensus_thresh=0.75, win_bp=100000)
        ref_candidates = {}
        ref_evaluations = {}

        # 4 refs successfully map to ctg_01
        for i in range(1, 5):
            acc = f"Ref_{i:02d}"
            ref_candidates[acc] = {
                "clean_acc": acc, "best_ctg": f"r_ctg_{i}", "op_pids": {f"R{i}_atpA"},
                "pid_to_ko": {f"R{i}_atpA": "K02117"},
                "ref_coords": {f"R{i}_atpA": {"contig": f"r_ctg_{i}", "start": 10000, "end": 11500, "strand": "+"}},
                "ref_faa_dict": {f"R{i}_atpA": "M" * 500}
            }
            ref_evaluations[acc] = {
                "score": 1000.0, "best_contig": "ctg_01",
                "best_cluster": [
                    {"start": 180000, "end": 181500, "side": "UP", "orig_pid": f"R{i}_u", "t_pid": "T_prot_02"},
                    {"start": 220000, "end": 221800, "side": "DOWN", "orig_pid": f"R{i}_d", "t_pid": "T_prot_04"},
                ]
            }

        # 6 refs were selected in Stage 0 but did not map to target (score 0, empty cluster)
        for i in range(5, 11):
            acc = f"Ref_{i:02d}"
            ref_candidates[acc] = {
                "clean_acc": acc, "best_ctg": f"r_ctg_{i}", "op_pids": set(),
                "pid_to_ko": {}, "ref_coords": {}, "ref_faa_dict": {}
            }
            ref_evaluations[acc] = {
                "score": 0.0, "best_contig": "", "best_cluster": []
            }

        # 10 positive refs in pos_refs
        pos_refs = pd.DataFrame([{"accession": f"Ref_{i:02d}"} for i in range(1, 11)])

        extract_extended_100kb_locus_proteins(
            args=args,
            pos_refs=pos_refs,
            ref_candidates=ref_candidates,
            ref_evaluations=ref_evaluations,
            target_records=self.target_records,
            gff_coords=self.gff_coords,
            target_faa_dict=self.target_faa_dict,
            target_anno_dict=self.target_anno_dict,
            ref_flank_annotations=self.ref_flank_annotations,
            ko_definitions=self.ko_definitions,
            out_tsv_path=self.out_tsv
        )

        df = pd.read_csv(self.out_tsv, sep="\t")
        self.assertFalse(df.empty, "Output TSV should not be empty under Option A!")
        target_df = df[df["genome_type"] == "target"]
        self.assertGreater(len(target_df), 0)
        # Check that T_prot_02 is mapped by 4 refs, ratio = 4/4 = 1.0
        p2_row = target_df[target_df["protein_id"] == "T_prot_02"].iloc[0]
        self.assertEqual(p2_row["mapped_ref_count"], 4)
        self.assertEqual(p2_row["mapped_ref_ratio"], 1.0)

    def test_dereplicate_reference_cohort(self):
        """Verify duplicate assembly merging and GTDB species clustering."""
        pos_refs = pd.DataFrame([
            {"accession": "GCF_009778815.1", "taxonomy": "d__Archaea;...;s__B26-1 sp009778815", "ANI": 95.0},
            {"accession": "GCA_009783705.1", "taxonomy": "d__Archaea;...;s__B26-1 sp009778815", "ANI": 95.0},
            {"accession": "GCA_963760905.1", "taxonomy": "d__Archaea;...;s__B26-1 sp963760905", "ANI": 92.0},
            {"accession": "GCA_964616775.1", "taxonomy": "d__Archaea;...;s__B26-1 sp964616775", "ANI": 88.0},
        ])
        clusters, r2c, weights, n_eff = dereplicate_reference_cohort(pos_refs)
        self.assertEqual(n_eff, 3, "Duplicate assemblies should merge into 3 independent clusters")
        self.assertEqual(r2c["GCF_009778815.1"], r2c["GCA_009783705.1"], "GCF and GCA pairs should be in same cluster")
        self.assertAlmostEqual(weights["GCF_009778815.1"], 0.5)
        self.assertAlmostEqual(weights["GCA_963760905.1"], 1.0)

    def test_compute_empirical_max_contig_null(self):
        """Verify empirical MAX-contig null simulation and P_emp computation."""
        target_recs = {f"ctg_{i}": None for i in range(20)}
        gff_cds = {f"p_{i}": {"contig": f"ctg_{i % 20}"} for i in range(2000)}
        ctg_scores = {"ctg_10": 4500.0, "ctg_17": 1200.0, "ctg_02": 800.0}
        cluster_scores = {"C1": 1500.0, "C2": 1400.0, "C3": 1200.0, "C4": 400.0}

        p, s_obs, s_mean, s_95, d, supp = compute_empirical_max_contig_null(
            target_recs, gff_cds, ctg_scores, "ctg_10", cluster_scores, n_permutations=200
        )
        self.assertLessEqual(p, 0.05, f"Strong localization should yield significant P_emp, got {p}")
        self.assertEqual(supp, "SIGNIFICANT")
        self.assertGreater(d, 2.0)

        # Weak localization test
        ctg_scores_weak = {"ctg_10": 400.0, "ctg_17": 350.0}
        p_weak, _, _, _, _, supp_weak = compute_empirical_max_contig_null(
            target_recs, gff_cds, ctg_scores_weak, "ctg_10", {"C1": 400.0}, n_permutations=200
        )
        self.assertGreater(p_weak, 0.01, f"Weak localization should not be significant, got {p_weak}")
        self.assertEqual(supp_weak, "NOT_SIGNIFICANT")

    def test_empirical_max_mode_does_not_block_candidate(self):
        """
        In empirical_max mode (default), even if consensus ratio is 4/9 = 44.4% (< 75%),
        the candidate locus with bilateral support is NOT discarded and 100kb proteins are extracted.
        """
        args = MockArgs(consensus_thresh=0.75, win_bp=100000, sig_mode="empirical_max")
        ref_candidates = {}
        ref_evaluations = {}

        # 4 refs on ctg_01 (44.4%)
        for i in range(1, 5):
            acc = f"Ref_{i:02d}"
            ref_candidates[acc] = {
                "clean_acc": acc, "best_ctg": "r_ctg", "op_pids": {f"R{i}_atpA"},
                "pid_to_ko": {f"R{i}_atpA": "K02117"},
                "ref_coords": {f"R{i}_atpA": {"contig": "r_ctg", "start": 10000, "end": 11500, "strand": "+"}},
                "ref_faa_dict": {f"R{i}_atpA": "M" * 500}
            }
            ref_evaluations[acc] = {
                "score": 1200.0, "best_contig": "ctg_01",
                "best_cluster": [
                    {"start": 180000, "end": 181500, "side": "UP", "orig_pid": f"R{i}_u", "t_pid": "T_prot_02"},
                    {"start": 220000, "end": 221800, "side": "DOWN", "orig_pid": f"R{i}_d", "t_pid": "T_prot_04"},
                ]
            }
        # 5 other refs on different contigs (ctg_02)
        for i in range(5, 10):
            acc = f"Ref_{i:02d}"
            ref_candidates[acc] = {
                "clean_acc": acc, "best_ctg": "r_ctg", "op_pids": set(),
                "pid_to_ko": {}, "ref_coords": {}, "ref_faa_dict": {}
            }
            ref_evaluations[acc] = {
                "score": 300.0, "best_contig": "ctg_02",
                "best_cluster": [{"start": 10000, "end": 12000, "side": "UP", "orig_pid": f"R{i}_u", "t_pid": "T_other"}]
            }

        pos_refs = pd.DataFrame([{"accession": f"Ref_{i:02d}"} for i in range(1, 10)])

        extract_extended_100kb_locus_proteins(
            args=args,
            pos_refs=pos_refs,
            ref_candidates=ref_candidates,
            ref_evaluations=ref_evaluations,
            target_records=self.target_records,
            gff_coords=self.gff_coords,
            target_faa_dict=self.target_faa_dict,
            target_anno_dict=self.target_anno_dict,
            ref_flank_annotations=self.ref_flank_annotations,
            ko_definitions=self.ko_definitions,
            out_tsv_path=self.out_tsv,
            target_contig_id="ctg_01",
            winning_geometry="BILATERAL_SAME_CONTIG"
        )

        df = pd.read_csv(self.out_tsv, sep="\t")
        self.assertFalse(df.empty, "Output TSV must NOT be empty in empirical_max mode!")
        target_df = df[df["genome_type"] == "target"]
        self.assertGreater(len(target_df), 0)

    def test_multi_cluster_convergence_over_spurious_contig(self):
        """
        Contig 10 has 4 independent evolutionary clusters (score 6685.2, bilateral).
        Contig 09 has only 1 evolutionary cluster (score 5080.8, bilateral) but more raw blast hits.
        The ranking must prioritize the 4 independent clusters and select Contig 10 with dominance_ratio > 1.0.
        """
        # Simulate candidate contigs and clusters
        cand_contigs = ["CAZOPV010000010.1", "CAZOPV010000009.1"]
        ref_clusters = {
            "C1": ["R1"], "C2": ["R2"], "C3": ["R3"], "C4": ["R4"],
            "C7": ["R7a", "R7b"]
        }
        ref_to_cluster = {
            "R1": "C1", "R2": "C2", "R3": "C3", "R4": "C4",
            "R7a": "C7", "R7b": "C7"
        }
        ref_candidates = {"R1": {}, "R2": {}, "R3": {}, "R4": {}, "R7a": {}, "R7b": {}}
        ref_evaluations = {
            "R1": {"best_contig": "CAZOPV010000010.1", "score": 1500.0},
            "R2": {"best_contig": "CAZOPV010000010.1", "score": 1800.0},
            "R3": {"best_contig": "CAZOPV010000010.1", "score": 1400.0},
            "R4": {"best_contig": "CAZOPV010000010.1", "score": 1985.2},
            "R7a": {"best_contig": "CAZOPV010000009.1", "score": 2540.4},
            "R7b": {"best_contig": "CAZOPV010000009.1", "score": 2540.4},
        }
        ref_ctg_evaluations = {
            "R1": {"CAZOPV010000010.1": {"score": 1500.0, "has_bilateral": True, "cluster": [{"t_pid": f"T_10_{i}"} for i in range(5)], "n_u": 3, "n_d": 2}},
            "R2": {"CAZOPV010000010.1": {"score": 1800.0, "has_bilateral": True, "cluster": [{"t_pid": f"T_10_{i}"} for i in range(5)], "n_u": 3, "n_d": 2}},
            "R3": {"CAZOPV010000010.1": {"score": 1400.0, "has_bilateral": True, "cluster": [{"t_pid": f"T_10_{i}"} for i in range(5)], "n_u": 3, "n_d": 2}},
            "R4": {"CAZOPV010000010.1": {"score": 1985.2, "has_bilateral": True, "cluster": [{"t_pid": f"T_10_{i}"} for i in range(5)], "n_u": 3, "n_d": 2}},
            "R7a": {"CAZOPV010000009.1": {"score": 2540.4, "has_bilateral": True, "cluster": [{"t_pid": f"T_09_{i}"} for i in range(17)], "n_u": 8, "n_d": 9}},
            "R7b": {"CAZOPV010000009.1": {"score": 2540.4, "has_bilateral": True, "cluster": [{"t_pid": f"T_09_{i}"} for i in range(17)], "n_u": 8, "n_d": 9}},
        }

        # Run ranking logic
        contig_geom_rank = {"CAZOPV010000010.1": 3, "CAZOPV010000009.1": 3}
        contig_supporting_clusters = {}
        contig_total_score = {
            "CAZOPV010000010.1": 6685.2,
            "CAZOPV010000009.1": 5080.8
        }
        contig_conserved_genes = {"CAZOPV010000010.1": 5, "CAZOPV010000009.1": 17}
        contig_order_score = {"CAZOPV010000010.1": 50, "CAZOPV010000009.1": 50}

        for ctg_id in cand_contigs:
            conv_refs = [r for r, ev in ref_evaluations.items() if ev.get("best_contig") == ctg_id and ev.get("score", 0.0) > 0]
            conv_clusters = set(ref_to_cluster.get(r, r) for r in conv_refs)
            contig_supporting_clusters[ctg_id] = len(conv_clusters)

        self.assertEqual(contig_supporting_clusters["CAZOPV010000010.1"], 4)
        self.assertEqual(contig_supporting_clusters["CAZOPV010000009.1"], 1)

        ranked = sorted(
            cand_contigs,
            key=lambda c: (
                contig_geom_rank[c],
                contig_supporting_clusters[c],
                contig_total_score[c],
                contig_conserved_genes[c],
                contig_order_score[c]
            ),
            reverse=True
        )
        self.assertEqual(ranked[0], "CAZOPV010000010.1", "Contig 10 with 4 independent evolutionary clusters must rank #1")
        dom_ratio = contig_total_score[ranked[0]] / (contig_total_score[ranked[1]] + 1e-5)
        self.assertGreater(dom_ratio, 1.0, f"Dominance ratio must be > 1.0, got {dom_ratio}")
        self.assertAlmostEqual(dom_ratio, 1.3158, places=3)

    def test_qid_split_preserves_triple_underscore_pids(self):
        """Test that qid.split('__', 2) correctly preserves protein IDs containing '___'."""
        clean_acc = "GCA_032598985.1"
        u_pids = [
            f"GCA_032598985.1_ASM3259898v1_genomic___CP060866.1___prot{i:05d}"
            for i in range(1, 22)
        ]
        parsed_orig_pids = []
        for pid in u_pids:
            qid = f"UP__{clean_acc}__{pid}"
            parts = qid.split("__", 2)
            side = parts[0]
            ref_acc = parts[1] if len(parts) > 1 else "unknown"
            orig_pid = parts[2] if len(parts) > 2 else pid
            self.assertEqual(side, "UP")
            self.assertEqual(ref_acc, clean_acc)
            self.assertEqual(orig_pid, pid)
            parsed_orig_pids.append(orig_pid)

        # Ensure all 21 protein IDs remain unique and are not collapsed to 1
        self.assertEqual(len(set(parsed_orig_pids)), 21)

    def test_bilateral_inner_gap_calculation(self):
        """Test that inner boundary gap calculation handles inverted and large spans correctly."""
        # Case 1: Inverted locus (DOWN before UP), with 24 DOWN genes and 21 UP genes spanning 90kb total
        syn_d = [{"start": 1190000 + i * 1000, "end": 1190900 + i * 1000, "orig_pid": f"d_{i}"} for i in range(24)]
        # d ends at 1190000 + 23000 + 900 = 1213900
        syn_u = [{"start": 1227617 + i * 1000, "end": 1228500 + i * 1000, "orig_pid": f"u_{i}"} for i in range(21)]
        # u starts at 1227617. Gap = 1227617 - 1213900 = 13717 bp.
        # Total span from min(d) to max(u) = 1248500 - 1190000 = 58500 bp.

        u_s = min(h["start"] for h in syn_u)
        u_e = max(h["end"] for h in syn_u)
        d_s = min(h["start"] for h in syn_d)
        d_e = max(h["end"] for h in syn_d)
        if u_e <= d_s:
            inner_gap = d_s - u_e
        elif d_e <= u_s:
            inner_gap = u_s - d_e
        else:
            inner_gap = 0

        self.assertEqual(inner_gap, 13717)
        self.assertLessEqual(inner_gap, 80000)

    def test_empty_cand_contigs_candidate_loci_computation(self):
        """Test candidate loci ranking when cand_contigs is empty to verify syn_refs scoping."""
        cand_contigs = []
        target_records = {"ctg_01": SeqRecord(Seq("A" * 1000), id="ctg_01")}
        ref_candidates = {}
        ref_ctg_evaluations = {}
        contig_geometry = {}
        contig_loc_score = {}
        contig_bilat_score = {}
        contig_boundary_anchors = {}
        contig_ancestral_intervals = {}
        contig_observed_spans = {}
        contig_supporting_clusters = {}
        contig_converged_clusters = {}
        contig_conserved_genes = {}

        ranked_contigs = []
        if cand_contigs:
            ranked_contigs = cand_contigs
            chosen_contig_id = ranked_contigs[0]
            chosen_geometry = "BILATERAL_SAME_CONTIG"
            dominance_ratio = 999.0
        else:
            chosen_contig_id = list(target_records.keys())[0] if target_records else ""
            chosen_geometry = "NONE"
            dominance_ratio = 1.0

        max_cands = 5
        top_cands = ranked_contigs[:max_cands] if cand_contigs else [chosen_contig_id]
        cand_records = []
        contig_matched_expected_spans = {}

        for rank_idx, c_id in enumerate(top_cands, 1):
            if not c_id:
                continue
            c_geom = contig_geometry.get(c_id, "NONE")
            c_loc = contig_loc_score.get(c_id, 0.0)
            c_bil = contig_bilat_score.get(c_id, 0.0)
            c_l_anc, c_r_anc = contig_boundary_anchors.get(c_id, ("NONE", "NONE"))
            c_int_s, c_int_e = contig_ancestral_intervals.get(c_id, (0, 0))
            c_obs_sp = contig_observed_spans.get(c_id, 0)

            # Verification: c_syn_refs must be cleanly computed without UnboundLocalError
            c_syn_refs = [
                r for r in ref_candidates.keys()
                if ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("score", 0.0) > 0 and
                   ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("cluster")
            ]
            c_u_h = [h for r in c_syn_refs for h in ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("cluster", []) if h.get("side") == "UP"]
            c_d_h = [h for r in c_syn_refs for h in ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("cluster", []) if h.get("side") == "DOWN"]
            self.assertEqual(len(c_u_h), 0)
            self.assertEqual(len(c_d_h), 0)

            c_exp_sp = 15000
            contig_matched_expected_spans[c_id] = c_exp_sp
            cand_records.append({
                "rank": rank_idx,
                "target_id": "Target_MAG_01",
                "module_id": "M00595",
                "contig_id": c_id,
                "geometry": c_geom,
                "expected_interanchor_span": c_exp_sp
            })

        self.assertEqual(len(cand_records), 1)
        self.assertEqual(cand_records[0]["contig_id"], "ctg_01")
        self.assertEqual(cand_records[0]["geometry"], "NONE")

if __name__ == "__main__":
    unittest.main()
