#!/usr/bin/env python3
"""
Unit and integration tests for Circular Contigs Illumina Short-Read Orthogonal
Abundance and Verification Pipeline (Tracks A, B, and C).
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path
import subprocess

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.validate_rotated_junction import (
    rotate_sequence,
    build_rotated_assembly,
    load_fasta_dict,
)
from scripts.compile_circular_illumina_qc import (
    parse_coverm_genome_tsv,
    parse_coverm_contig_tsv,
    parse_rotated_junction_qc_tsv,
    compile_summary,
)
from scripts.audit_circular_assembly import (
    parse_illumina_qc_row,
    evaluate_four_gates,
    evaluate_five_gates,
)
from scripts.plot_circular_audit import (
    load_audit_metrics,
    generate_verification_plot,
)


class TestCircularIllumina(unittest.TestCase):

    def test_rotate_sequence(self):
        """Test rotation of circular sequences at fraction L/2."""
        # 10 bp sequence: indices 0..9. L=10, M=5
        seq = "ABCDE12345"
        rot, m, j_rot = rotate_sequence(seq, fraction=0.5)
        self.assertEqual(m, 5)
        self.assertEqual(j_rot, 5)
        # Expected: seq[5:] + seq[:5] = "12345ABCDE"
        self.assertEqual(rot, "12345ABCDE")

        # Odd length: L=9, M=4. seq[:4] = "ABCD", seq[4:] = "E1234"
        seq_odd = "ABCDE1234"
        rot_odd, m_odd, j_rot_odd = rotate_sequence(seq_odd, fraction=0.5)
        self.assertEqual(m_odd, 4)
        self.assertEqual(rot_odd, "E1234ABCD")

    def test_build_rotated_assembly(self):
        """Test rotated assembly FASTA and coordinates TSV generation."""
        with tempfile.TemporaryDirectory() as td:
            asm_fa = Path(td) / "assembly.fasta"
            asm_fa.write_text(">ctg1\nACGTACGT\n>ctg2\n12345678\n>linear1\nGGGGTTTT\n")

            bins_dir = Path(td) / "circular_bins"
            bins_dir.mkdir()
            (bins_dir / "ctg1.fa").write_text(">ctg1\nACGTACGT\n")
            (bins_dir / "ctg2.fa").write_text(">ctg2\n12345678\n")

            out_asm = Path(td) / "rotated.fasta"
            out_coords = Path(td) / "coords.tsv"

            rotated_cnt = build_rotated_assembly(str(asm_fa), str(bins_dir), str(out_asm), str(out_coords), fraction=0.5)
            self.assertEqual(rotated_cnt, 2)
            self.assertTrue(out_asm.exists())
            self.assertTrue(out_coords.exists())

            rot_seqs = load_fasta_dict(str(out_asm))
            self.assertIn("ctg1", rot_seqs)
            self.assertIn("ctg2", rot_seqs)
            self.assertIn("linear1", rot_seqs)
            # ctg2: "12345678" -> "56781234"
            self.assertEqual(rot_seqs["ctg2"], "56781234")
            # linear1 remains unchanged
            self.assertEqual(rot_seqs["linear1"], "GGGGTTTT")

            coords_lines = out_coords.read_text().splitlines()
            self.assertEqual(len(coords_lines), 3)  # header + 2 contigs
            self.assertIn("rotated_junction_coord", coords_lines[0])

    def test_compile_circular_illumina_qc(self):
        """Test parsing CoverM A/B/C tracks, ratio calculation, flags and overall status."""
        with tempfile.TemporaryDirectory() as td:
            # Track A: genome abundance
            track_a = Path(td) / "track_a.tsv"
            track_a.write_text(
                "Genome\tRelative Abundance (%)\tMean\tTrimmed Mean\tCovered Fraction\tRead Count\n"
                "unmapped\t20.0\t0.0\t0.0\t0.0\t2000\n"
                "ctg1\t80.0\t50.0\t48.0\t0.995\t8000\n"
            )

            # Track B: contig coverage standard (95/75)
            track_b = Path(td) / "track_b.tsv"
            track_b.write_text(
                "Contig\tMean\tTrimmed Mean\tCovered Fraction\tVariance\tRead Count\n"
                "ctg1\t49.0\t47.0\t0.995\t12.0\t7900\n"
                "linear1\t10.0\t9.5\t0.850\t5.0\t1000\n"
            )

            # Track B: contig coverage strict (98/90)
            track_b_strict = Path(td) / "track_b_strict.tsv"
            track_b_strict.write_text(
                "Contig\tTrimmed Mean\tCovered Fraction\tRead Count\n"
                "ctg1\t46.0\t0.990\t7700\n"
            )

            # Track C: rotated junction QC
            track_c = Path(td) / "track_c.tsv"
            track_c.write_text(
                "contig_id\tcontig_length\trotated_junction_coord\tjunction_median_depth_100bp\t"
                "junction_median_depth_300bp\tjunction_median_depth_1kb\tlocal_flank_median_depth\t"
                "genome_median_depth\tjunction_local_ratio\tjunction_genome_ratio\t"
                "zero_coverage_run_bp\tspanning_reads_unique\tspanning_reads_unique_mapq20\t"
                "bridging_pairs_unique\tbridging_pairs_unique_mapq20\n"
                "ctg1\t10000\t5000\t45.0\t46.0\t46.5\t47.0\t48.0\t0.957\t0.938\t0\t15\t14\t22\t20\n"
            )

            out_summary = Path(td) / "illumina_summary.tsv"
            compile_summary(
                sample="test_sample",
                genome_abundance_tsv=str(track_a),
                contig_coverage_tsv=str(track_b),
                strict_coverage_tsv=str(track_b_strict),
                junction_qc_tsv=str(track_c),
                out_tsv=str(out_summary)
            )

            self.assertTrue(out_summary.exists())
            lines = out_summary.read_text().splitlines()
            self.assertEqual(len(lines), 2)
            headers = lines[0].split("\t")
            row = dict(zip(headers, lines[1].split("\t")))

            self.assertEqual(row["sample"], "test_sample")
            self.assertEqual(row["circular_genome"], "ctg1")
            self.assertAlmostEqual(float(row["relative_abundance_percent"]), 80.0)
            self.assertAlmostEqual(float(row["catalogue_unmapped_percent"]), 20.0)
            self.assertAlmostEqual(float(row["catalogue_mapped_percent"]), 80.0)
            # C = contig_trimmed / genome_trimmed = 47.0 / 48.0 = 0.979
            self.assertAlmostEqual(float(row["competitive_depth_ratio"]), 0.979, places=2)
            # S = strict_trimmed / standard_trimmed = 46.0 / 47.0 = 0.979
            self.assertAlmostEqual(float(row["strict_depth_ratio"]), 0.979, places=2)
            self.assertEqual(row["junction_spanning_reads_unique"], "15")
            self.assertEqual(row["flag_junction_gap"], "FALSE")
            self.assertEqual(row["flag_possible_homologous_recruitment"], "FALSE")
            self.assertEqual(row["overall_status"], "PASS_HIGH_CONFIDENCE")

    def test_missing_track_c_handling(self):
        """Test that missing Track C data outputs NA metrics and UNVERIFIED_NO_ROTATED_DATA without false conflict."""
        with tempfile.TemporaryDirectory() as td:
            track_a = Path(td) / "track_a.tsv"
            track_a.write_text(
                "Genome\tRelative Abundance (%)\tMean\tTrimmed Mean\tCovered Fraction\tRead Count\n"
                "u499\t50.0\t30.0\t29.0\t0.99\t5000\n"
            )
            track_b = Path(td) / "track_b.tsv"
            track_b.write_text(
                "Contig\tMean\tTrimmed Mean\tCovered Fraction\tVariance\tRead Count\n"
                "u499\t29.0\t28.0\t0.99\t5.0\t4900\n"
            )
            track_b_strict = Path(td) / "track_b_strict.tsv"
            track_b_strict.write_text(
                "Contig\tTrimmed Mean\tCovered Fraction\tRead Count\n"
                "u499\t27.0\t0.99\t4800\n"
            )
            # Empty Track C (or another contig only)
            track_c = Path(td) / "track_c.tsv"
            track_c.write_text(
                "contig\tlength\trotation_offset\trotated_junction_coord\t"
                "junction_median_depth_100bp\tjunction_median_depth_300bp\tjunction_median_depth_1kb\t"
                "local_flank_median_depth\twhole_contig_median_depth\tjunction_local_ratio\t"
                "junction_genome_ratio\tjunction_zero_coverage_run_bp\tjunction_spanning_reads\t"
                "junction_spanning_reads_unique\tjunction_bridging_pairs\tjunction_bridging_pairs_unique\n"
                "other_contig\t1000\t500\t500\t20.0\t20.0\t20.0\t20.0\t20.0\t1.0\t1.0\t0\t5\t5\t5\t5\n"
            )

            out_summary = Path(td) / "summary.tsv"
            compile_summary(
                sample="sample_test",
                genome_abundance_tsv=str(track_a),
                contig_coverage_tsv=str(track_b),
                strict_coverage_tsv=str(track_b_strict),
                junction_qc_tsv=str(track_c),
                out_tsv=str(out_summary)
            )

            lines = out_summary.read_text().splitlines()
            headers = lines[0].split("\t")
            row_map = {r.split("\t")[1]: dict(zip(headers, r.split("\t"))) for r in lines[1:]}
            self.assertIn("u499", row_map)
            row = row_map["u499"]
            self.assertEqual(row["circular_genome"], "u499")
            self.assertEqual(row["junction_median_depth_100bp"], "NA")
            self.assertEqual(row["junction_local_ratio"], "NA")
            self.assertEqual(row["junction_zero_coverage_run_bp"], "NA")
            self.assertEqual(row["junction_spanning_reads_unique"], "NA")
            self.assertEqual(row["flag_junction_depth_drop"], "FALSE")
            self.assertEqual(row["flag_no_junction_span"], "FALSE")
            self.assertEqual(row["overall_status"], "UNVERIFIED_NO_ROTATED_DATA")

            # Verify Gate 5 handling in audit_circular_assembly
            ill_metrics = parse_illumina_qc_row(str(out_summary), "u499")
            self.assertEqual(ill_metrics["Gate5_Illumina_Orthogonal_Status"], "NOT_EVALUATED")
            self.assertIn("No_Rotated_Junction_Data", ill_metrics["Gate5_Notes"])

            # Verify evaluation does not downgrade to SUSPECT_CIRCULAR_JUNCTION
            *_, grade = evaluate_five_gates(
                contig_len=1000000,
                spanning_unique_reads=25,
                support_fraction=0.85,
                clip_fraction=0.05,
                unresolved_zero_bp=0,
                has_high_risk_zero_block=False,
                clipping_sites_count=0,
                partner_reads=0,
                partner_fraction=0.0,
                partner_name="None",
                anchor_status="UNIQUE",
                primary_mean_depth=30.0,
                primary_breadth_1x=99.8,
                repeat_breadth_pct=0.0,
                longest_repeat_bp=0,
                comp=98.0,
                cont=1.0,
                r5s=1,
                r16s=1,
                r23s=1,
                trna_aa_count=20,
                replicon_state="CHROMOSOME_SUPPORTED",
                gate5_status=ill_metrics["Gate5_Illumina_Orthogonal_Status"]
            )
            self.assertNotEqual(grade, "SUSPECT_CIRCULAR_JUNCTION (Illumina Orthogonal Conflict)")
            self.assertIn("CIRCULAR_CHROMOSOME", grade)

    def test_parse_illumina_qc_row_and_gate5_evaluation(self):
        """Test parse_illumina_qc_row and Gate 5 status / grade arbitration."""
        with tempfile.TemporaryDirectory() as td:
            summary_tsv = Path(td) / "summary.tsv"
            summary_tsv.write_text(
                "sample\tcircular_genome\tlength\trelative_abundance_percent\t"
                "catalogue_mapped_percent\tcatalogue_unmapped_percent\tgenome_mean_depth\t"
                "genome_trimmed_mean_depth\tgenome_trimmed_to_mean_ratio\tgenome_covered_fraction\t"
                "genome_aligned_read_count\tgenome_anir_percent\tcontig_trimmed_mean_depth\t"
                "contig_covered_fraction\tcontig_anir_percent\tcompetitive_depth_ratio\t"
                "strict_trimmed_mean_depth\tstrict_covered_fraction\tstrict_depth_ratio\t"
                "junction_median_depth_100bp\tlocal_flank_median_depth\tjunction_local_ratio\t"
                "junction_zero_coverage_run_bp\tjunction_spanning_reads_unique\t"
                "junction_bridging_pairs_unique\tflag_low_relative_abundance\tflag_low_depth\t"
                "flag_low_breadth\tflag_possible_homologous_recruitment\tflag_identity_sensitive\t"
                "flag_junction_depth_drop\tflag_junction_gap\tflag_no_junction_span\toverall_status\n"
                "Bathy\tu100\t1500000\t5.2\t90.0\t10.0\t35.0\t34.5\t0.98\t0.998\t50000\t99.2\t"
                "34.0\t0.998\t99.2\t0.985\t33.0\t0.995\t0.971\t32.0\t34.0\t0.941\t0\t25\t30\t"
                "FALSE\tFALSE\tFALSE\tFALSE\tFALSE\tFALSE\tFALSE\tFALSE\tPASS_HIGH_CONFIDENCE\n"
                "Bathy\tu200\t800000\t0.05\t90.0\t10.0\t2.0\t1.8\t0.90\t0.80\t2000\t98.0\t"
                "0.5\t0.50\t98.0\t0.278\t0.4\t0.45\t0.80\t0.0\t2.0\t0.0\t150\t0\t0\t"
                "TRUE\tTRUE\tTRUE\tTRUE\tFALSE\tTRUE\tTRUE\tTRUE\tJUNCTION_SUSPECT\n"
            )

            # 1. Contig u100: PASS_HIGH_CONFIDENCE
            metrics_u100 = parse_illumina_qc_row(str(summary_tsv), "u100")
            self.assertEqual(metrics_u100["Gate5_Illumina_Orthogonal_Status"], "PASS")
            self.assertIn("Orthogonally_Supported", metrics_u100["Gate5_Notes"])
            self.assertEqual(metrics_u100["Illumina_Relative_Abundance_Pct"], "5.2")

            # Evaluate gates for u100: should upgrade to ORTHOGONALLY_VERIFIED
            *_, grade_u100 = evaluate_five_gates(
                contig_len=1500000,
                spanning_unique_reads=20,
                support_fraction=0.85,
                clip_fraction=0.05,
                unresolved_zero_bp=0,
                has_high_risk_zero_block=False,
                clipping_sites_count=0,
                partner_reads=0,
                partner_fraction=0.0,
                partner_name="None",
                anchor_status="UNIQUE",
                primary_mean_depth=30.0,
                primary_breadth_1x=99.8,
                repeat_breadth_pct=0.0,
                longest_repeat_bp=0,
                comp=98.0,
                cont=1.0,
                r5s=1,
                r16s=1,
                r23s=1,
                trna_aa_count=20,
                replicon_state="CHROMOSOME_SUPPORTED",
                gate5_status=metrics_u100["Gate5_Illumina_Orthogonal_Status"]
            )
            self.assertEqual(grade_u100, "ORTHOGONALLY_VERIFIED_HIGH_QUALITY_CIRCULAR_CHROMOSOME")

            # 2. Contig u200: JUNCTION_SUSPECT -> JUNCTION_CONFLICT
            metrics_u200 = parse_illumina_qc_row(str(summary_tsv), "u200")
            self.assertEqual(metrics_u200["Gate5_Illumina_Orthogonal_Status"], "JUNCTION_CONFLICT")
            self.assertIn("Junction_Zero_Gap(150bp)", metrics_u200["Gate5_Notes"])

            # Evaluate gates for u200: even if Gate 1-4 pass, Gate 5 conflict overrides!
            *_, grade_u200 = evaluate_five_gates(
                contig_len=800000,
                spanning_unique_reads=20,
                support_fraction=0.85,
                clip_fraction=0.05,
                unresolved_zero_bp=0,
                has_high_risk_zero_block=False,
                clipping_sites_count=0,
                partner_reads=0,
                partner_fraction=0.0,
                partner_name="None",
                anchor_status="UNIQUE",
                primary_mean_depth=30.0,
                primary_breadth_1x=99.8,
                repeat_breadth_pct=0.0,
                longest_repeat_bp=0,
                comp=98.0,
                cont=1.0,
                r5s=1,
                r16s=1,
                r23s=1,
                trna_aa_count=20,
                replicon_state="CHROMOSOME_SUPPORTED",
                gate5_status=metrics_u200["Gate5_Illumina_Orthogonal_Status"]
            )
            self.assertEqual(grade_u200, "SUSPECT_CIRCULAR_JUNCTION (Illumina Orthogonal Conflict)")

            # 3. Missing / Nonexistent TSV -> NOT_APPLICABLE
            metrics_na = parse_illumina_qc_row(None, "u300")
            self.assertEqual(metrics_na["Gate5_Illumina_Orthogonal_Status"], "NOT_APPLICABLE")
            self.assertEqual(metrics_na["Gate5_Notes"], "No_Short_Reads")
            self.assertEqual(metrics_na["Illumina_Relative_Abundance_Pct"], "NA")

            *_, grade_na = evaluate_five_gates(
                contig_len=1500000,
                spanning_unique_reads=20,
                support_fraction=0.85,
                clip_fraction=0.05,
                unresolved_zero_bp=0,
                has_high_risk_zero_block=False,
                clipping_sites_count=0,
                partner_reads=0,
                partner_fraction=0.0,
                partner_name="None",
                anchor_status="UNIQUE",
                primary_mean_depth=30.0,
                primary_breadth_1x=99.8,
                repeat_breadth_pct=0.0,
                longest_repeat_bp=0,
                comp=98.0,
                cont=1.0,
                r5s=1,
                r16s=1,
                r23s=1,
                trna_aa_count=20,
                replicon_state="CHROMOSOME_SUPPORTED",
                gate5_status=metrics_na["Gate5_Illumina_Orthogonal_Status"]
            )
            self.assertEqual(grade_na, "READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME")

    def test_audit_script_cli_with_illumina_qc(self):
        """End-to-end CLI execution of audit_circular_assembly.py with --illumina-qc-tsv."""
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "ctg_test.fa"
            fa.write_text(">ctg_test\n" + "ACGT" * 500 + "\n")  # 2000 bp
            out_tsv = Path(td) / "audit_out.tsv"

            # Create mock Illumina summary TSV
            ill_tsv = Path(td) / "illumina_summary.tsv"
            ill_tsv.write_text(
                "sample\tcircular_genome\tlength\trelative_abundance_percent\t"
                "catalogue_mapped_percent\tcatalogue_unmapped_percent\tgenome_mean_depth\t"
                "genome_trimmed_mean_depth\tgenome_trimmed_to_mean_ratio\tgenome_covered_fraction\t"
                "genome_aligned_read_count\tgenome_anir_percent\tcontig_trimmed_mean_depth\t"
                "contig_covered_fraction\tcontig_anir_percent\tcompetitive_depth_ratio\t"
                "strict_trimmed_mean_depth\tstrict_covered_fraction\tstrict_depth_ratio\t"
                "junction_median_depth_100bp\tlocal_flank_median_depth\tjunction_local_ratio\t"
                "junction_zero_coverage_run_bp\tjunction_spanning_reads_unique\t"
                "junction_bridging_pairs_unique\tflag_low_relative_abundance\tflag_low_depth\t"
                "flag_low_breadth\tflag_possible_homologous_recruitment\tflag_identity_sensitive\t"
                "flag_junction_depth_drop\tflag_junction_gap\tflag_no_junction_span\toverall_status\n"
                "sample1\tctg_test\t2000\t12.5\t95.0\t5.0\t45.0\t44.0\t0.978\t0.999\t15000\t99.5\t"
                "43.5\t0.999\t99.5\t0.989\t42.0\t0.998\t0.966\t41.0\t43.0\t0.953\t0\t18\t24\t"
                "FALSE\tFALSE\tFALSE\tFALSE\tFALSE\tFALSE\tFALSE\tFALSE\tPASS_HIGH_CONFIDENCE\n"
            )

            cmd = [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "audit_circular_assembly.py"),
                "--bin", "ctg_test",
                "--contig-fa", str(fa),
                "--whole-bam", str(Path(td) / "nonexistent.bam"),
                "--illumina-qc-tsv", str(ill_tsv),
                "--out-tsv", str(out_tsv)
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, 0, f"Script failed: {res.stderr}")
            self.assertTrue(out_tsv.exists())

            lines = out_tsv.read_text().splitlines()
            self.assertEqual(len(lines), 2)
            headers = lines[0].split("\t")
            values = lines[1].split("\t")
            self.assertEqual(len(headers), len(values))
            self.assertEqual(len(headers), 82)
            self.assertIn("Assembly_Usability_Grade", headers)
            self.assertIn("Circularity_Confidence_Tier", headers)
            self.assertIn("Biological_Complexity_Tier", headers)
            self.assertIn("Graph_Complexity_Tier", headers)

            row = dict(zip(headers, values))
            self.assertEqual(row["Contig_ID"], "ctg_test")
            self.assertEqual(row["Illumina_Relative_Abundance_Pct"], "12.5")
            self.assertEqual(row["Illumina_Competitive_Depth_Ratio_C"], "0.989")
            self.assertEqual(row["Gate5_Illumina_Orthogonal_Status"], "PASS")

    def test_plot_circular_audit_with_illumina_qc(self):
        """Test verification plot generation with Illumina QC metrics."""
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "ctg_plot.fa"
            fa.write_text(">ctg_plot\n" + "ACGT" * 1000 + "\n")  # 4000 bp

            audit_tsv = Path(td) / "audit.tsv"
            audit_tsv.write_text(
                "Contig_ID\tLength_bp\tPrimary_Mean_Depth\tJunction_Spanning_Unique_Reads\t"
                "Gate5_Illumina_Orthogonal_Status\tIllumina_Relative_Abundance_Pct\t"
                "Illumina_Genome_Trimmed_Depth\tIllumina_Competitive_Depth_Ratio_C\t"
                "Illumina_Rotated_Spanning_Reads_Unique\tIllumina_Rotated_Bridging_Pairs_Unique\t"
                "Final_Confidence_Grade\n"
                "ctg_plot\t4000\t35.0\t15\tPASS\t8.5\t42.0\t0.96\t18\t22\t"
                "ORTHOGONALLY_VERIFIED_HIGH_QUALITY_CIRCULAR_CHROMOSOME\n"
            )

            metrics = load_audit_metrics(str(audit_tsv), "ctg_plot")
            self.assertEqual(metrics["Gate5_Illumina_Orthogonal_Status"], "PASS")
            self.assertEqual(metrics["Illumina_Relative_Abundance_Pct"], 8.5)

            out_pdf = Path(td) / "plot.pdf"
            out_png = Path(td) / "plot.png"
            generate_verification_plot(
                bin_id="ctg_plot",
                contig_fa=str(fa),
                whole_bam=None,
                junction_bam=None,
                audit_tsv=str(audit_tsv),
                bakta_tsv=None,
                out_pdf=str(out_pdf),
                out_png=str(out_png),
                flank_len=500,
                max_pileup_reads=5
            )
            self.assertTrue(out_pdf.exists())
            self.assertTrue(out_png.exists())
            self.assertGreater(out_pdf.stat().st_size, 0)
    def test_contig_header_prefix_matching(self):
        """Test header prefix matching where assembly BAM has ctg_len suffix but circular bin has clean ID."""
        with tempfile.TemporaryDirectory() as td:
            # 1. Rotated assembly matching
            asm_fa = Path(td) / "assembly.fasta"
            asm_fa.write_text(">u377048ctg_len-1543357_circular-yes_depth-27-27-27\nACGTACGT\n>other_ctg\nGGGGTTTT\n")

            bins_dir = Path(td) / "circular_bins"
            bins_dir.mkdir()
            (bins_dir / "u377048.fa").write_text(">u377048\nACGTACGT\n")

            out_asm = Path(td) / "rotated.fasta"
            out_coords = Path(td) / "coords.tsv"

            rotated_cnt = build_rotated_assembly(str(asm_fa), str(bins_dir), str(out_asm), str(out_coords), fraction=0.5)
            self.assertEqual(rotated_cnt, 1)

            # 2. Compile summary matching
            track_a = Path(td) / "track_a.tsv"
            track_a.write_text(
                "Genome\tRelative Abundance (%)\tMean\tTrimmed Mean\tCovered Fraction\tRead Count\n"
                "unmapped\t50.0\t0.0\t0.0\t0.0\t5000\n"
                "u377048\t0.458\t27.2\t27.17\t0.9997\t292303\n"
            )

            track_b = Path(td) / "track_b.tsv"
            track_b.write_text(
                "Contig\tMean\tTrimmed Mean\tCovered Fraction\tVariance\tRead Count\n"
                "u377048ctg_len-1543357_circular-yes_depth-27-27-27\t27.0\t26.8\t0.9995\t10.0\t280000\n"
            )

            track_b_strict = Path(td) / "track_b_strict.tsv"
            track_b_strict.write_text(
                "Contig\tTrimmed Mean\tCovered Fraction\tRead Count\n"
                "u377048ctg_len-1543357_circular-yes_depth-27-27-27\t26.0\t0.9990\t270000\n"
            )

            track_c = Path(td) / "track_c.tsv"
            track_c.write_text(
                "contig_id\tcontig_length\trotated_junction_coord\tjunction_median_depth_100bp\t"
                "junction_median_depth_300bp\tjunction_median_depth_1kb\tlocal_flank_median_depth\t"
                "genome_median_depth\tjunction_local_ratio\tjunction_genome_ratio\t"
                "zero_coverage_run_bp\tspanning_reads_unique\tspanning_reads_unique_mapq20\t"
                "bridging_pairs_unique\tbridging_pairs_unique_mapq20\n"
                "u377048\t1543357\t771678\t25.0\t26.0\t26.5\t27.0\t27.0\t0.926\t0.926\t0\t28\t25\t35\t30\n"
            )

            out_summary = Path(td) / "summary.tsv"
            compile_summary(
                sample="Bathy",
                genome_abundance_tsv=str(track_a),
                contig_coverage_tsv=str(track_b),
                strict_coverage_tsv=str(track_b_strict),
                junction_qc_tsv=str(track_c),
                out_tsv=str(out_summary)
            )

            self.assertTrue(out_summary.exists())
            lines = out_summary.read_text().splitlines()
            self.assertEqual(len(lines), 2)
            row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            self.assertEqual(row["circular_genome"], "u377048")
            self.assertAlmostEqual(float(row["contig_trimmed_mean_depth"]), 26.8)
            self.assertAlmostEqual(float(row["competitive_depth_ratio"]), 0.986, places=2)
            self.assertEqual(row["overall_status"], "PASS_HIGH_CONFIDENCE")


if __name__ == "__main__":
    unittest.main()
