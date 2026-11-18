#!/usr/bin/env python3
"""
Unit and integration tests for scripts/audit_circular_assembly.py (Version 2.1).
Tests coverage decoupling, CIGAR-aware junction spanning, BLASTN repeat filtering,
canonical 20-AA tRNA normalization, geNomad replicon arbitration,
Methanothrix premature circularization detection, and Four-Gate classifications.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch
from pathlib import Path
import subprocess

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.audit_circular_assembly import (
    wilson_score_interval,
    compute_primary_coherence_metrics,
    find_zero_coverage_blocks,
    parse_blastn_repeats,
    parse_bakta_biological_markers,
    parse_checkm2,
    parse_gtdb,
    parse_genomad,
    resolve_bam_reference,
    evaluate_four_gates,
    analyze_targeted_synthetic_junction,
    find_systematic_clipping_sites,
    detect_methanothrix_partner_clustering,
    check_terminal_anchor_uniqueness,
    calculate_coverage_arrays,
    CANONICAL_20_AA,
)


class TestCircularAudit(unittest.TestCase):

    def test_wilson_score_interval(self):
        # Zero observations
        lower, upper = wilson_score_interval(0, 0)
        self.assertEqual((lower, upper), (0.0, 0.0))

        # 10 out of 10
        lower, upper = wilson_score_interval(10, 10)
        self.assertGreater(lower, 0.65)
        self.assertLessEqual(upper, 1.0)

        # 5 out of 10
        lower, upper = wilson_score_interval(5, 10)
        self.assertTrue(0.20 < lower < 0.50)
        self.assertTrue(0.50 < upper < 0.85)

    def test_compute_primary_coherence_metrics(self):
        # 100 bp contig with uniform 20x depth
        cov = [20] * 100
        mean_d, median_d, mad_d, cv_d, b1, b3 = compute_primary_coherence_metrics(cov, 100)
        self.assertAlmostEqual(mean_d, 20.0)
        self.assertAlmostEqual(median_d, 20.0)
        self.assertAlmostEqual(mad_d, 0.0)
        self.assertAlmostEqual(cv_d, 0.0)
        self.assertAlmostEqual(b1, 100.0)
        self.assertAlmostEqual(b3, 100.0)

        # Variable coverage with zero depth regions
        cov_var = [0] * 10 + [5] * 40 + [25] * 50
        mean_d, median_d, mad_d, cv_d, b1, b3 = compute_primary_coherence_metrics(cov_var, 100)
        self.assertAlmostEqual(b1, 90.0)
        self.assertAlmostEqual(b3, 90.0)
        self.assertGreater(cv_d, 0.0)

    def test_find_zero_coverage_blocks(self):
        # 1000 bp contig with a 50 bp gap and a 600 bp gap
        cov = [10] * 1000
        # gap 1: 100-150 (50 bp)
        for i in range(100, 150):
            cov[i] = 0
        # gap 2: 300-900 (600 bp, high-risk >500 bp)
        for i in range(300, 900):
            cov[i] = 0

        blocks, total_zero, max_zero, high_risk = find_zero_coverage_blocks(cov, 1000)
        self.assertEqual(len(blocks), 2)
        self.assertEqual(total_zero, 650)
        self.assertEqual(max_zero, 600)
        self.assertEqual(len(high_risk), 1)
        self.assertEqual(high_risk[0]["len"], 600)

    def test_parse_blastn_repeats(self):
        with tempfile.TemporaryDirectory() as td:
            blast_tsv = Path(td) / "blast.tsv"
            # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            blast_tsv.write_text(
                # Exact diagonal (should be excluded)
                "ctg1\tctg1\t100.0\t1000\t0\t0\t1\t1000\t1\t1000\t0.0\t1800\n"
                # Low identity <80% (should be filtered)
                "ctg1\tctg1\t75.0\t400\t50\t0\t200\t600\t1500\t1900\t1e-5\t300\n"
                # Short length <=200 bp (should be filtered)
                "ctg1\tctg1\t99.0\t150\t1\t0\t100\t250\t2000\t2150\t1e-10\t250\n"
                # Valid repeat 1: 300 bp, 95% id, positions 1000-1300
                "ctg1\tctg1\t95.0\t300\t15\t0\t1000\t1300\t5000\t5300\t1e-50\t550\n"
                # Valid repeat 2: overlapping with repeat 1 (1200-1600, 400 bp) -> merged 1000-1600 = 600 bp
                "ctg1\tctg1\t90.0\t400\t20\t0\t1200\t1600\t7000\t7400\t1e-60\t700\n"
            )

            breadth_pct, longest_rep, hit_count = parse_blastn_repeats(str(blast_tsv), contig_len=10000)
            self.assertEqual(hit_count, 2)
            self.assertEqual(longest_rep, 400)
            # Merged interval: 1000 to 1600 = 600 bp out of 10000 = 6.0%
            self.assertAlmostEqual(breadth_pct, 6.0)

    def test_parse_bakta_biological_markers(self):
        with tempfile.TemporaryDirectory() as td:
            bakta_tsv = Path(td) / "bakta.tsv"
            lines = [
                "# Sequence Id\tType\tProduct\tGene\tStart\tStop",
                "ctg1\trRNA\t5S ribosomal RNA\trrf\t100\t215",
                "ctg1\trRNA\t16S ribosomal RNA\trrs\t300\t1800",
                "ctg1\trRNA\t23S ribosomal RNA\trrl\t2000\t5000",
                # Canonical tRNAs
                "ctg1\ttRNA\ttRNA-Ala(tgc)\ttraA\t5100\t5175",
                "ctg1\ttRNA\ttRNA-Arg(ccg)\ttraR\t5200\t5275",
                "ctg1\ttRNA\ttRNA-Asn(gtt)\ttraN\t5300\t5375",
                "ctg1\ttRNA\ttRNA-Asp(gtc)\ttraD\t5400\t5475",
                "ctg1\ttRNA\ttRNA-Cys(gca)\ttraC\t5500\t5575",
                "ctg1\ttRNA\ttRNA-Gln(ttg)\ttraQ\t5600\t5675",
                "ctg1\ttRNA\ttRNA-Glu(ttc)\ttraE\t5700\t5775",
                "ctg1\ttRNA\ttRNA-Gly(gcc)\ttraG\t5800\t5875",
                "ctg1\ttRNA\ttRNA-His(gtg)\ttraH\t5900\t5975",
                "ctg1\ttRNA\ttRNA-Ile(gat)\ttraI\t6000\t6075",
                "ctg1\ttRNA\ttRNA-Ile2(cat)\ttraI2\t6100\t6175",  # Normalize Ile2 -> Ile
                "ctg1\ttRNA\ttRNA-Leu(caa)\ttraL\t6200\t6275",
                "ctg1\ttRNA\ttRNA-Lys(ttt)\ttraK\t6300\t6375",
                "ctg1\ttRNA\ttRNA-fMet(cat)\ttraM\t6400\t6475",  # Normalize fMet -> Met
                "ctg1\ttRNA\ttRNA-Phe(gaa)\ttraF\t6500\t6575",
                "ctg1\ttRNA\ttRNA-Pro(tgg)\ttraP\t6600\t6675",
                "ctg1\ttRNA\ttRNA-Ser(tga)\ttraS\t6700\t6775",
                "ctg1\ttRNA\ttRNA-Thr(ggt)\ttraT\t6800\t6875",
                "ctg1\ttRNA\ttRNA-Trp(cca)\ttraW\t6900\t6975",
                "ctg1\ttRNA\ttRNA-Tyr(gta)\ttraY\t7000\t7075",
                "ctg1\ttRNA\ttRNA-Val(tac)\ttraV\t7100\t7175",
                # Non-canonical / special (should be discarded)
                "ctg1\ttRNA\ttRNA-Sec(tca)\ttraSec\t7200\t7280",
                "ctg1\ttRNA\ttRNA-Pyl(cta)\ttraPyl\t7300\t7380",
            ]
            bakta_tsv.write_text("\n".join(lines) + "\n")

            r5s, r16s, r23s, trna_count = parse_bakta_biological_markers(str(bakta_tsv))
            self.assertEqual(r5s, 1)
            self.assertEqual(r16s, 1)
            self.assertEqual(r23s, 1)
            # All 20 canonical amino acids should be accounted for (20/20)
            self.assertEqual(trna_count, 20)

    def test_parse_bakta_official_9_columns(self):
        with tempfile.TemporaryDirectory() as td:
            bakta_tsv = Path(td) / "bakta_official.tsv"
            # Official Bakta output:
            # #Sequence Id \t Type \t Start \t Stop \t Strand \t Locus Tag \t Gene \t Product \t DbXrefs
            lines = [
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs",
                "u999856\trRNA\t100\t215\t+\tu999856_0001\trrf\t5S ribosomal RNA\tRFAM:RF00001",
                "u999856\trRNA\t300\t1800\t+\tu999856_0002\trrs\t16S ribosomal RNA\tRFAM:RF00177",
                "u999856\trRNA\t2000\t5000\t+\tu999856_0003\trrl\t23S ribosomal RNA\tRFAM:RF02541",
                "u999856\ttRNA\t5100\t5175\t+\tu999856_0004\ttraA\ttRNA-Ala(tgc)\tSO:0000254",
                "u999856\ttRNA\t5200\t5275\t+\tu999856_0005\ttraM\ttRNA-fMet(cat)\tSO:0000254",
                "u999856\tCDS\t6000\t7000\t+\tu999856_0006\tdnaA\tchromosomal replication initiator protein DnaA\tUniRef:UniRef50_P03004",
            ]
            bakta_tsv.write_text("\n".join(lines) + "\n")

            r5s, r16s, r23s, trna_count = parse_bakta_biological_markers(str(bakta_tsv))
            self.assertEqual(r5s, 1)
            self.assertEqual(r16s, 1)
            self.assertEqual(r23s, 1)
            self.assertEqual(trna_count, 2)  # Ala and Met (from fMet)

    def test_resolve_bam_reference(self):
        refs = [
            "u999856_len-2079472_circular-yes_depth-35.2_duplicated-0",
            "u499975_len-1500000_circular-yes",
            "contig_3",
            "bin.4.fa length=50000"
        ]
        # Prefix match
        self.assertEqual(
            resolve_bam_reference(refs, "u999856"),
            "u999856_len-2079472_circular-yes_depth-35.2_duplicated-0"
        )
        self.assertEqual(
            resolve_bam_reference(refs, "u499975"),
            "u499975_len-1500000_circular-yes"
        )
        # Exact match
        self.assertEqual(resolve_bam_reference(refs, "contig_3"), "contig_3")
        # Token match
        self.assertEqual(resolve_bam_reference(refs, "bin.4.fa"), "bin.4.fa length=50000")
        # Unmatched
        self.assertEqual(resolve_bam_reference(refs, "nonexistent"), "nonexistent")

    def test_parse_checkm2_and_gtdb_fuzzy_names(self):
        with tempfile.TemporaryDirectory() as td:
            checkm2_tsv = Path(td) / "checkm2.tsv"
            checkm2_tsv.write_text(
                "Name\tCompleteness\tContamination\n"
                "u999856.fa\t98.43\t1.13\n"
            )
            gtdb_tsv = Path(td) / "gtdb.tsv"
            gtdb_tsv.write_text(
                "user_genome\tclassification\n"
                "u999856.fa\td__Archaea;p__Thermoproteota;c__Bathyarchaeia;o__Bathyarchaeales;f__Bathyarchaeaceae;g__Bathyarchaeum;s__Bathyarchaeum sp.\n"
            )

            comp, cont = parse_checkm2(str(checkm2_tsv), "u999856")
            self.assertEqual(comp, 98.43)
            self.assertEqual(cont, 1.13)

            lineage = parse_gtdb(str(gtdb_tsv), "u999856")
            self.assertTrue("Bathyarchaeia" in lineage)

    def test_parse_genomad(self):
        with tempfile.TemporaryDirectory() as td:
            genomad_tsv = Path(td) / "genomad.tsv"
            genomad_tsv.write_text(
                "seq_name\tchromosome_score\tplasmid_score\tvirus_score\n"
                "ctg_chrom\t0.92\t0.05\t0.03\n"
                "ctg_plasmid\t0.08\t0.88\t0.04\n"
                "ctg_virus\t0.02\t0.08\t0.90\n"
                "ctg_ambig\t0.45\t0.42\t0.13\n"
            )

            # Chromosome supported
            raw, state, margin = parse_genomad(str(genomad_tsv), "ctg_chrom")
            self.assertEqual(state, "CHROMOSOME_SUPPORTED")
            self.assertAlmostEqual(margin, 0.92 - 0.05)

            # Plasmid supported
            raw, state, margin = parse_genomad(str(genomad_tsv), "ctg_plasmid")
            self.assertEqual(state, "PLASMID_SUPPORTED")

            # Virus supported
            raw, state, margin = parse_genomad(str(genomad_tsv), "ctg_virus")
            self.assertEqual(state, "VIRUS_SUPPORTED")

            # Ambiguous
            raw, state, margin = parse_genomad(str(genomad_tsv), "ctg_ambig")
            self.assertEqual(state, "AMBIGUOUS")

    def test_four_gates_classification(self):
        # 1. Top Tier: READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=2500000,
            spanning_unique_reads=15,
            support_fraction=0.85,
            clip_fraction=0.05,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=35.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=4.2,
            longest_repeat_bp=1200,
            comp=98.5,
            cont=0.8,
            r5s=2,
            r16s=2,
            r23s=2,
            trna_aa_count=20,
            replicon_state="CHROMOSOME_SUPPORTED",
        )
        self.assertEqual(g1, "PASS")
        self.assertEqual(g3, "PASS")
        self.assertEqual(grade, "READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME")

        # 2. Viral Replicon: READ_SUPPORTED_CIRCULAR_VIRAL_REPLICON
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=45000,
            spanning_unique_reads=12,
            support_fraction=0.90,
            clip_fraction=0.02,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=50.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=2.0,
            longest_repeat_bp=300,
            comp=None,
            cont=None,
            r5s=0,
            r16s=0,
            r23s=0,
            trna_aa_count=1,
            replicon_state="VIRUS_SUPPORTED",
        )
        self.assertEqual(g1, "PASS")
        self.assertEqual(grade, "READ_SUPPORTED_CIRCULAR_VIRAL_REPLICON")

        # 3. High Confidence Circular Plasmid
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=85000,
            spanning_unique_reads=8,
            support_fraction=0.75,
            clip_fraction=0.08,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=25.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=1.5,
            longest_repeat_bp=250,
            comp=None,
            cont=None,
            r5s=0,
            r16s=0,
            r23s=0,
            trna_aa_count=0,
            replicon_state="PLASMID_SUPPORTED",
        )
        self.assertEqual(g1, "PASS")
        self.assertEqual(grade, "HIGH_CONFIDENCE_CIRCULAR_PLASMID")

        # 4. Premature Circularization Suspect (Methanothrix phenomenon)
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=1800000,
            spanning_unique_reads=6,
            support_fraction=0.55,
            clip_fraction=0.10,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=14,
            partner_fraction=0.72,
            partner_name="contig_linear_002",
            anchor_status="UNIQUE",
            primary_mean_depth=30.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=5.0,
            longest_repeat_bp=800,
            comp=85.0,
            cont=1.0,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=16,
            replicon_state="CHROMOSOME_SUPPORTED",
        )
        self.assertEqual(g1, "FAIL")
        self.assertIn("Premature_Circularization_Risk", g1_note)
        self.assertEqual(grade, "PREMATURE_CIRCULARIZATION_SUSPECT (Consistent Partner Alignment)")

        # 5. Potential Misassembly (Internal clipping / zero-coverage gap)
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=2000000,
            spanning_unique_reads=10,
            support_fraction=0.80,
            clip_fraction=0.05,
            unresolved_zero_bp=850,
            has_high_risk_zero_block=True,
            clipping_sites_count=2,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=20.0,
            primary_breadth_1x=98.5,
            repeat_breadth_pct=3.0,
            longest_repeat_bp=400,
            comp=95.0,
            cont=2.0,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=19,
            replicon_state="CHROMOSOME_SUPPORTED",
        )
        self.assertEqual(g1, "FAIL")
        self.assertEqual(grade, "POTENTIAL_MISASSEMBLY (Physical Breakpoint / Unresolved Gap)")

        # 6. Insufficient Physical Evidence (Zero spanning reads)
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=1500000,
            spanning_unique_reads=0,
            support_fraction=0.0,
            clip_fraction=0.0,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=1.5,
            primary_breadth_1x=40.0,
            repeat_breadth_pct=0.0,
            longest_repeat_bp=0,
            comp=60.0,
            cont=0.0,
            r5s=0,
            r16s=0,
            r23s=0,
            trna_aa_count=5,
            replicon_state="AMBIGUOUS",
        )
        self.assertEqual(g1, "FAIL")
        self.assertEqual(grade, "INSUFFICIENT_PHYSICAL_EVIDENCE (Low Depth / Non-Informative)")

    def test_sv_polymorphism_gate1_pass(self):
        # Scenario: Contig has systematic clipping at 1,485,346 bp, but it is quantified as an authentic
        # in-vivo deletion polymorphism (AF_del = 20.3%, 294 deletion reads vs 1154 insertion reads).
        # Gate 1 should PASS_WITH_SV_HETEROGENEITY and overall grade should be READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME_WITH_POPULATION_HETEROGENEITY.
        sv_detail = [{
            "pos": 1485346,
            "clipped_reads": 294,
            "depth": 1448,
            "clipping_ratio": 0.203,
            "is_sv_polymorphism": True,
            "n_ins": 1154,
            "n_del": 294,
            "del_freq": 0.203,
            "del_partner_pos": 1489770,
            "del_span_bp": 4424
        }]
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=2079472,
            spanning_unique_reads=34,
            support_fraction=0.97,
            clip_fraction=0.03,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=1,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=2483.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=0.5,
            longest_repeat_bp=300,
            comp=99.85,
            cont=0.77,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=20,
            replicon_state="CHROMOSOME_SUPPORTED",
            gate5_status="NOT_APPLICABLE",
            clipping_sites_detail=sv_detail
        )
        self.assertEqual(g1, "PASS_WITH_SV_HETEROGENEITY")
        self.assertIn("Population_SV_Polymorphism", g1_note)
        self.assertEqual(grade, "READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME_WITH_POPULATION_HETEROGENEITY")

        # Test with Gate 5 Short-Read Recruitment Caution
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=2079472,
            spanning_unique_reads=34,
            support_fraction=0.97,
            clip_fraction=0.03,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=1,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=2483.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=0.5,
            longest_repeat_bp=300,
            comp=99.85,
            cont=0.77,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=20,
            replicon_state="CHROMOSOME_SUPPORTED",
            gate5_status="RECRUITMENT_CAUTION",
            clipping_sites_detail=sv_detail
        )
        self.assertEqual(g1, "PASS_WITH_SV_HETEROGENEITY")
        self.assertEqual(grade, "READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME_WITH_POPULATION_HETEROGENEITY (Short-Read Recruitment Caution)")

    def test_low_abundance_partner_tolerance(self):
        # Scenario: Contig has 10 reads clipped near terminus pointing to a 36 kb contig,
        # but the partner contig is an in-silico assembly chimera with only 0.15% relative depth (RDR=0.0015).
        # It should NOT trigger Premature_Circularization_Risk failure on the main chromosome.
        g1, g1_note, g2, g3, g3_note, g4, grade = evaluate_four_gates(
            contig_len=2079472,
            spanning_unique_reads=34,
            support_fraction=0.97,
            clip_fraction=0.03,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=10,
            partner_fraction=0.65,
            partner_name="u1804999ctg",
            anchor_status="UNIQUE",
            primary_mean_depth=2483.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=0.5,
            longest_repeat_bp=300,
            comp=99.85,
            cont=0.77,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=20,
            replicon_state="CHROMOSOME_SUPPORTED",
            is_low_abundance_partner=True,
            partner_rel_depth_ratio=0.0015
        )
        self.assertEqual(g1, "PASS")
        self.assertIn("Low_Abundance_Partner_Artifact", g1_note)
        self.assertNotIn("Premature_Circularization_Risk", g1_note)
        self.assertEqual(grade, "READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME")

    def test_terminal_anchor_low_abundance_partner_exemption(self):
        # When spanning reads have secondary hits pointing to a low-abundance partner artifact,
        # check_terminal_anchor_uniqueness should exempt them and report "UNIQUE" rather than AMBIGUOUS.
        spanning_qnames = ["read1", "read2"]
        mock_hit1 = MagicMock()
        mock_hit1.query_name = "read1"
        mock_hit1.mapping_quality = 60
        mock_hit1.has_tag.return_value = False

        mock_hit2 = MagicMock()
        mock_hit2.query_name = "read2"
        mock_hit2.mapping_quality = 0  # low MAPQ due to low-abundance partner duplicate
        mock_hit2.has_tag.side_effect = lambda t: t in ("SA", "XA")
        mock_hit2.get_tag.side_effect = lambda t: "u1804999ctg,100,+,50M,60;" if t == "SA" else "u1804999ctg,+100,50M,0;"

        mock_bam = MagicMock()
        mock_bam.references = ["u999856"]
        mock_bam.fetch.return_value = [mock_hit1, mock_hit2]

        with patch("scripts.audit_circular_assembly.pysam") as mock_pysam:
            mock_pysam.AlignmentFile.return_value.__enter__.return_value = mock_bam
            with patch("os.path.exists", return_value=True):
                status, ambig = check_terminal_anchor_uniqueness(
                    "mock.bam", "u999856", 2000000, spanning_qnames, 5000,
                    low_abundance_partners={"u1804999ctg"}
                )
                self.assertEqual(status, "UNIQUE")
                self.assertEqual(ambig, 0)

    def test_small_plasmid_dynamic_flank(self):
        # Flank logic: effective_flank = min(config_flank, L // 2)
        config_flank = 5000

        # Normal chromosome (2 Mb) -> 5000 bp
        chrom_len = 2000000
        self.assertEqual(min(config_flank, chrom_len // 2), 5000)

        # 6 kb plasmid -> 3000 bp
        p6k_len = 6000
        self.assertEqual(min(config_flank, p6k_len // 2), 3000)

        # 800 bp small plasmid -> 400 bp
        p800_len = 800
        self.assertEqual(min(config_flank, p800_len // 2), 400)

    def test_mock_targeted_synthetic_junction(self):
        # Test CIGAR-aware junction logic with mock pysam reads
        flank_len = 5000
        min_flank = 500
        J = flank_len  # 5000
        left_target = J - min_flank  # 4500
        right_target = J + min_flank  # 5500

        # Create mock read that spans continuously: start 4000, length 2000 -> covers 4000 to 6000
        mock_spanning_read = MagicMock()
        mock_spanning_read.is_unmapped = False
        mock_spanning_read.query_name = "read_span_01"
        mock_spanning_read.reference_start = 4000
        mock_spanning_read.reference_end = 6000
        mock_spanning_read.cigartuples = [(0, 2000)]  # 2000M

        # Create mock read that clips at J: start 4000, end 5000, soft clipped at right
        mock_clipped_read = MagicMock()
        mock_clipped_read.is_unmapped = False
        mock_clipped_read.query_name = "read_clip_02"
        mock_clipped_read.reference_start = 4000
        mock_clipped_read.reference_end = 5000
        mock_clipped_read.cigartuples = [(0, 1000), (4, 500)]  # 1000M 500S

        # Create mock ONT read with typical indels spanning junction: start 3500, end 6500
        mock_ont_read = MagicMock()
        mock_ont_read.is_unmapped = False
        mock_ont_read.query_name = "read_ont_indels_03"
        mock_ont_read.reference_start = 3500
        mock_ont_read.reference_end = 6500
        mock_ont_read.cigartuples = [(0, 400), (1, 2), (0, 300), (2, 1), (0, 500), (1, 1), (0, 800)]

        mock_bam = MagicMock()
        mock_bam.fetch.return_value = [mock_spanning_read, mock_clipped_read, mock_ont_read]

        with patch("scripts.audit_circular_assembly.pysam") as mock_pysam:
            mock_pysam.AlignmentFile.return_value.__enter__.return_value = mock_bam
            with patch("os.path.exists", return_value=True):
                (n_span, n_clip, n_inf,
                 s_frac, c_frac,
                 ci_l, ci_u, span_qnames) = analyze_targeted_synthetic_junction("mock.bam", flank_len, min_flank)

                self.assertEqual(n_span, 2)
                self.assertEqual(n_clip, 1)
                self.assertEqual(n_inf, 3)
                self.assertAlmostEqual(s_frac, 2.0 / 3.0)
                self.assertAlmostEqual(c_frac, 1.0 / 3.0)
                self.assertIn("read_span_01", span_qnames)
                self.assertIn("read_ont_indels_03", span_qnames)

    def test_mock_methanothrix_partner_clustering(self):
        # Test partner contig tracing with mock supplementary alignments
        contig_len = 1000000
        flank_len = 5000

        mock_reads = []
        for i in range(5):
            r = MagicMock()
            r.is_unmapped = False
            r.query_name = f"partner_read_{i}"
            r.reference_start = 100
            r.reference_end = 1500
            r.cigartuples = [(4, 200), (0, 1200)]  # Soft clip at start
            r.has_tag.side_effect = lambda tag: tag == "SA"
            r.get_tag.side_effect = lambda tag: "contig_linear_003,50000,+,100,10,0;" if tag == "SA" else None
            mock_reads.append(r)

        mock_bam = MagicMock()
        mock_bam.references = ["ctg1"]
        mock_bam.fetch.return_value = mock_reads

        with patch("scripts.audit_circular_assembly.pysam") as mock_pysam:
            mock_pysam.AlignmentFile.return_value.__enter__.return_value = mock_bam
            with patch("os.path.exists", return_value=True):
                (p_reads, p_name, p_frac,
                 p_mad, total_clipped,
                 p_depth, rdr, is_artifact) = detect_methanothrix_partner_clustering("mock.bam", "ctg1", contig_len, flank_len)
                self.assertEqual(p_name, "contig_linear_003")
                self.assertEqual(p_reads, 5)

    def test_checkm2_and_gtdb_parsing(self):
        with tempfile.TemporaryDirectory() as td:
            checkm2_tsv = Path(td) / "quality_report.tsv"
            checkm2_tsv.write_text(
                "Name\tCompleteness\tContamination\n"
                "bin_001.fa\t97.5\t1.2\n"
                "bin_002\t88.0\t6.5\n"
            )
            gtdb_tsv = Path(td) / "gtdbtk.summary.tsv"
            gtdb_tsv.write_text(
                "user_genome\tclassification\n"
                "bin_001\td__Bacteria;p__Bacteroidota;c__Bacteroidia\n"
                "bin_002.fa\td__Archaea;p__Methanobacteriota\n"
            )

            c1, k1 = parse_checkm2(str(checkm2_tsv), "bin_001")
            self.assertAlmostEqual(c1, 97.5)
            self.assertAlmostEqual(k1, 1.2)

            c2, k2 = parse_checkm2(str(checkm2_tsv), "bin_002")
            self.assertAlmostEqual(c2, 88.0)
            self.assertAlmostEqual(k2, 6.5)

            lin1 = parse_gtdb(str(gtdb_tsv), "bin_001")
            self.assertIn("p__Bacteroidota", lin1)

            lin2 = parse_gtdb(str(gtdb_tsv), "bin_002")
            self.assertIn("p__Methanobacteriota", lin2)

    def test_cli_end_to_end(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "ctg.fa"
            fa.write_text(">ctg1\n" + "ACGT" * 500 + "\n")  # 2000 bp
            bakta = Path(td) / "bakta.tsv"
            bakta.write_text("ctg1\trRNA\t5S ribosomal RNA\n")
            out_tsv = Path(td) / "audit_out.tsv"

            cmd = [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "audit_circular_assembly.py"),
                "--bin", "ctg1",
                "--contig-fa", str(fa),
                "--whole-bam", str(Path(td) / "nonexistent.bam"),
                "--bakta-tsv", str(bakta),
                "--out-tsv", str(out_tsv)
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, 0, f"Script failed: {res.stderr}")
            self.assertTrue(out_tsv.exists())
            content = out_tsv.read_text().splitlines()
            self.assertEqual(len(content), 2)
            headers = content[0].split("\t")
            values = content[1].split("\t")
            self.assertIn("Contig_ID", headers)
            self.assertIn("Final_Confidence_Grade", headers)
            row = dict(zip(headers, values))
            self.assertEqual(row["Contig_ID"], "ctg1")
            self.assertEqual(row["Length_bp"], "2000")

    def test_extract_terminal_qnames_script(self):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "ctg_6k.fa"
            # 6000 bp sequence: 3000 'A's + 3000 'C's
            fa.write_text(">ctg_6k\n" + "A" * 3000 + "C" * 3000 + "\n")
            out_junc = Path(td) / "junction_ref.fa"
            out_qnames = Path(td) / "qnames.txt"

            cmd = [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "extract_terminal_qnames.py"),
                "--contig", str(fa),
                "--bam", str(Path(td) / "nonexistent.bam"),
                "--bin", "ctg_6k",
                "--flank", "5000",
                "--out-junction", str(out_junc),
                "--out-qnames", str(out_qnames)
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, 0, f"Script failed: {res.stderr}")
            self.assertTrue(out_junc.exists())
            self.assertTrue(out_qnames.exists())

            junc_content = out_junc.read_text().splitlines()
            # Header should reflect effective flank 3000
            self.assertIn("eff_3000", junc_content[0])
            self.assertIn("junction_len_6000", junc_content[0])
            # Sequence should be last 3000 ('C'*3000) + first 3000 ('A'*3000)
            seq = "".join(junc_content[1:])
            self.assertEqual(len(seq), 6000)
            self.assertEqual(seq, "C" * 3000 + "A" * 3000)

    def test_recover_terminal_fastq_script(self):
        import subprocess
        import gzip
        with tempfile.TemporaryDirectory() as td:
            raw_fq_gz = Path(td) / "raw.fq.gz"
            records = [
                "@read_1 molecule 1\nACGTACGT\n+\nIIIIIIII\n",
                "@read_2 molecule 2\nTTTTGGGG\n+\nIIIIIIII\n",
                "@read_3 molecule 3\nCCCCAAAA\n+\nIIIIIIII\n",
                "@read_4 molecule 4\nGGGGCCCC\n+\nIIIIIIII\n",
            ]
            with gzip.open(raw_fq_gz, "wt") as gz:
                gz.write("".join(records))

            qnames_file = Path(td) / "qnames.txt"
            qnames_file.write_text("read_2\nread_4\n")

            out_reads_gz = Path(td) / "recovered.fq.gz"

            cmd = [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "recover_terminal_fastq.py"),
                "--qnames", str(qnames_file),
                "--raw-reads", str(raw_fq_gz),
                "--out-reads", str(out_reads_gz)
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, 0, f"Script failed: {res.stderr}")
            self.assertTrue(out_reads_gz.exists())

            with gzip.open(out_reads_gz, "rt") as gz:
                recovered_lines = gz.read().splitlines()
            # 2 reads = 8 lines
            self.assertEqual(len(recovered_lines), 8)
            recovered_headers = [line for line in recovered_lines if line.startswith("@")]
            self.assertEqual(recovered_headers, ["@read_2 molecule 2", "@read_4 molecule 4"])

    def test_plot_circular_audit_helpers(self):
        from scripts.plot_circular_audit import (
            load_fasta,
            load_audit_metrics,
            compute_gc_and_depth_bins,
            parse_bakta_markers_coordinates,
            generate_verification_plot
        )
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "ctg.fa"
            fa.write_text(">ctg1\nACGTACGTACGTACGT\n" * 100)
            l, s = load_fasta(str(fa))
            self.assertEqual(l, 1600)

            tsv = Path(td) / "audit.tsv"
            tsv.write_text(
                "Contig_ID\tLength_bp\tPrimary_Mean_Depth\tJunction_Spanning_Unique_Reads\tFinal_Confidence_Grade\n"
                "ctg1\t1600\t50.5\t12\tREAD_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME\n"
            )
            m = load_audit_metrics(str(tsv), "ctg1")
            self.assertEqual(float(m["Primary_Mean_Depth"]), 50.5)

            centers, gc_p, skew, cov = compute_gc_and_depth_bins(s, None, "ctg1", bin_size=200)
            self.assertEqual(len(centers), 8)
            self.assertAlmostEqual(gc_p[0], 50.0)

            bakta = Path(td) / "bakta.tsv"
            bakta.write_text(
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\n"
                "ctg1\trRNA\t100\t200\t+\tL1\trrs\t16S ribosomal RNA\n"
                "ctg1\ttRNA\t300\t370\t-\tL2\ttrnE\ttRNA-Glu(ttc)\n"
            )
            rrna_c, trna_c = parse_bakta_markers_coordinates(str(bakta))
            self.assertEqual(len(rrna_c), 1)
            self.assertEqual(len(trna_c), 1)

            # Test generation of plot files (PDF + PNG)
            out_pdf = Path(td) / "test_plot.pdf"
            out_png = Path(td) / "test_plot.png"
            generate_verification_plot(
                bin_id="ctg1",
                contig_fa=str(fa),
                whole_bam=None,
                junction_bam=None,
                audit_tsv=str(tsv),
                bakta_tsv=str(bakta),
                out_pdf=str(out_pdf),
                out_png=str(out_png),
                flank_len=500,
                max_pileup_reads=10
            )
            self.assertTrue(out_pdf.exists())
            self.assertTrue(out_png.exists())
            self.assertGreater(out_pdf.stat().st_size, 0)
            self.assertGreater(out_png.stat().st_size, 0)

    def test_audit_breakpoint_and_gap_coordinate_reporting(self):
        """Verify that evaluate_four_gates records specific coordinates in Gate1_Notes."""
        clipping_sites = [{"pos": 50000, "clipped_reads": 20, "depth": 22, "clipping_ratio": 0.909}]
        zero_gaps = [{"start": 100000, "end": 101000, "len": 1000}]

        gate1_status, gate1_notes, _, _, _, _, final_grade = evaluate_four_gates(
            contig_len=200000,
            spanning_unique_reads=25,
            support_fraction=0.85,
            clip_fraction=0.05,
            unresolved_zero_bp=1000,
            has_high_risk_zero_block=True,
            clipping_sites_count=1,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=30.0,
            primary_breadth_1x=99.5,
            repeat_breadth_pct=0.0,
            longest_repeat_bp=0,
            comp=98.0,
            cont=1.0,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=20,
            replicon_state="CHROMOSOME_SUPPORTED",
            clipping_sites_detail=clipping_sites,
            zero_gaps_detail=zero_gaps
        )

        self.assertEqual(gate1_status, "FAIL")
        self.assertEqual(final_grade, "POTENTIAL_MISASSEMBLY (Physical Breakpoint / Unresolved Gap)")
        self.assertIn("pos=50000", gate1_notes)
        self.assertIn("100000-101000", gate1_notes)

    def test_plot_circular_audit_event_extraction_and_misassembly_plot(self):
        """Verify that plot script extracts coordinates and renders misassembly markers."""
        from scripts.plot_circular_audit import (
            extract_audit_events_from_metrics,
            load_audit_metrics,
            generate_verification_plot,
        )
        with tempfile.TemporaryDirectory() as td:
            tsv = Path(td) / "misassembly_audit.tsv"
            tsv.write_text(
                "Contig_ID\tLength_bp\tPrimary_Mean_Depth\tJunction_Spanning_Unique_Reads\t"
                "Final_Confidence_Grade\tGate1_Notes\tInternal_Clipping_Sites_Detail\tZero_Coverage_Gaps_Detail\n"
                "u377048\t1543357\t27.0\t28\tPOTENTIAL_MISASSEMBLY (Physical Breakpoint / Unresolved Gap)\t"
                "Systematic_Internal_Clipping_Sites(1:pos=820500(24/27:89%));High_Risk_Zero_Coverage_Gap(801200-801900(700bp))\t"
                "820500:24/27:88.9%\t801200-801900:700bp\n"
            )
            metrics = load_audit_metrics(str(tsv), "u377048")
            self.assertEqual(metrics["Internal_Clipping_Sites_Detail"], "820500:24/27:88.9%")
            self.assertEqual(metrics["Zero_Coverage_Gaps_Detail"], "801200-801900:700bp")

            clip_bps, zero_gaps = extract_audit_events_from_metrics(metrics)
            self.assertEqual(len(clip_bps), 1)
            self.assertEqual(clip_bps[0]["pos"], 820500)
            self.assertEqual(clip_bps[0]["count"], 24)
            self.assertAlmostEqual(clip_bps[0]["ratio"], 0.889)

            self.assertEqual(len(zero_gaps), 1)
            self.assertEqual(zero_gaps[0]["start"], 801200)
            self.assertEqual(zero_gaps[0]["end"], 801900)
            self.assertEqual(zero_gaps[0]["len"], 700)

            # Test generation of plot files with these alerts
            fa = Path(td) / "u377048.fa"
            fa.write_text(">u377048\n" + "ACGT" * 25000 + "\n")  # 100 kb sequence
            out_pdf = Path(td) / "u377048_verification_plot.pdf"
            out_png = Path(td) / "u377048_verification_plot.png"
            generate_verification_plot(
                bin_id="u377048",
                contig_fa=str(fa),
                whole_bam=None,
                junction_bam=None,
                audit_tsv=str(tsv),
                bakta_tsv="",
                out_pdf=str(out_pdf),
                out_png=str(out_png),
                flank_len=500,
                max_pileup_reads=10
            )
            self.assertTrue(out_pdf.exists())
            self.assertTrue(out_png.exists())
            self.assertGreater(out_pdf.stat().st_size, 0)
            self.assertGreater(out_png.stat().st_size, 0)

    def test_calculate_coverage_arrays_pysam_attributes(self):
        """Verify that calculate_coverage_arrays correctly handles AlignedSegment without is_primary attribute."""
        mock_primary = MagicMock()
        mock_primary.is_unmapped = False
        mock_primary.is_secondary = False
        mock_primary.is_supplementary = False
        mock_primary.reference_start = 0
        mock_primary.cigartuples = [(0, 50)]
        del mock_primary.is_primary  # Ensure accessing is_primary would raise AttributeError

        mock_secondary = MagicMock()
        mock_secondary.is_unmapped = False
        mock_secondary.is_secondary = True
        mock_secondary.is_supplementary = False
        mock_secondary.reference_start = 25
        mock_secondary.cigartuples = [(0, 50)]
        del mock_secondary.is_primary

        mock_bam = MagicMock()
        mock_bam.references = ["ctg1"]
        mock_bam.fetch.return_value = [mock_primary, mock_secondary]

        import scripts.audit_circular_assembly as aca
        mock_pysam = MagicMock()
        mock_pysam.AlignmentFile.return_value.__enter__.return_value = mock_bam

        with patch.object(aca, "pysam", mock_pysam):
            with tempfile.NamedTemporaryFile(suffix=".bam") as tf:
                p_cov, a_cov = calculate_coverage_arrays(tf.name, "ctg1", 100)

                # Primary read covers [0, 50)
                self.assertEqual(int(p_cov[10]), 1)
                self.assertEqual(int(p_cov[60]), 0)

                # All cov includes secondary read [25, 75)
                self.assertEqual(int(a_cov[10]), 1)
                self.assertEqual(int(a_cov[30]), 2)
                self.assertEqual(int(a_cov[60]), 1)
                self.assertEqual(int(a_cov[80]), 0)

    def test_audit_script_zero_blocks_reporting(self):
        """Verify that audit_circular_assembly CLI handles gaps_to_report without NameError when blocks >= 100 bp."""
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "ctg_gap.fa"
            fa.write_text(">ctg_gap\n" + "ACGT" * 500 + "\n")  # 2000 bp
            out_tsv = Path(td) / "out_gap.tsv"

            # Run CLI with no reads mapped -> 2000 bp zero coverage block (high risk >500 bp)
            cmd = [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "audit_circular_assembly.py"),
                "--bin", "ctg_gap",
                "--contig-fa", str(fa),
                "--whole-bam", str(Path(td) / "nonexistent.bam"),
                "--out-tsv", str(out_tsv)
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, 0, f"Script failed: {res.stderr}")
            self.assertTrue(out_tsv.exists())
            content = out_tsv.read_text().splitlines()
            self.assertEqual(len(content), 2)
            row = dict(zip(content[0].split("\t"), content[1].split("\t")))
            self.assertIn("0-2000", row["Zero_Coverage_Gaps_Detail"])


class TestDecoupledCircularAuditArchitecture(unittest.TestCase):
    """
    Comprehensive tests for the Two-Dimensional Decoupled Audit Architecture
    (Circularity Confidence vs Population / Graph Complexity).
    """

    def test_calibration_u667600_textbook_positive(self):
        """Case 1: u667600 textbook positive (VERIFIED_CIRCULAR, CLONAL_SIMPLE)."""
        res = evaluate_four_gates(
            contig_len=3200000,
            spanning_unique_reads=63,
            support_fraction=0.692,
            clip_fraction=0.04,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=72.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=2.1,
            longest_repeat_bp=450,
            comp=99.2,
            cont=0.5,
            r5s=2,
            r16s=2,
            r23s=2,
            trna_aa_count=20,
            replicon_state="CHROMOSOME_SUPPORTED",
            gate5_status="PASS",
            gate5_metrics={"Illumina_Rotated_Spanning_Reads_Unique": 73, "Illumina_Rotated_Bridging_Pairs_Unique": 99}
        )
        self.assertEqual(res.circularity.tier, "VERIFIED_CIRCULAR")
        self.assertEqual(res.circularity.evidence_type, "CONCORDANT_DUAL_OMICS")
        self.assertEqual(res.complexity.bio_tier, "CLONAL_SIMPLE")
        self.assertEqual(res.complexity.graph_tier, "GRAPH_SIMPLE")
        self.assertEqual(res.assembly_grade, "HIGH_CONFIDENCE_CIRCULAR")
        self.assertEqual(res.final_grade, "ORTHOGONALLY_VERIFIED_HIGH_QUALITY_CIRCULAR_CHROMOSOME")

    def test_calibration_u999856_multi_complex_benchmark(self):
        """Case 2: u999856 multi-complex benchmark (VERIFIED_CIRCULAR, MULTI_COMPLEX)."""
        sv_detail = [{
            "pos": 1485346,
            "clipped_reads": 294,
            "depth": 1448,
            "clipping_ratio": 0.203,
            "is_sv_polymorphism": True,
            "n_ins": 1154,
            "n_del": 294,
            "del_freq": 0.203,
            "del_span_bp": 4424
        }]
        res = evaluate_four_gates(
            contig_len=2079472,
            spanning_unique_reads=34,
            support_fraction=0.97,
            clip_fraction=0.03,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=1,
            partner_reads=10,
            partner_fraction=0.65,
            partner_name="u1804999",
            anchor_status="UNIQUE",
            primary_mean_depth=2483.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=0.5,
            longest_repeat_bp=300,
            comp=99.85,
            cont=0.77,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=20,
            replicon_state="CHROMOSOME_SUPPORTED",
            gate5_status="RECRUITMENT_CAUTION",
            clipping_sites_detail=sv_detail,
            is_low_abundance_partner=True,
            partner_rel_depth_ratio=0.0015
        )
        self.assertEqual(res.circularity.tier, "VERIFIED_CIRCULAR")
        self.assertEqual(res.circularity.evidence_type, "LONG_READ_DEFINITIVE")
        self.assertEqual(res.complexity.bio_tier, "MULTI_COMPLEX")
        self.assertEqual(res.complexity.graph_tier, "LOW_FREQUENCY_ALTERNATIVE_PATH")
        self.assertEqual(res.complexity.rdr_signal, "RARE_GRAPH_SIGNAL")
        self.assertEqual(res.assembly_grade, "HIGH_CONFIDENCE_CIRCULAR")
        self.assertIn("WITH_POPULATION_COMPLEXITY", res.final_grade)
        self.assertIn("Short-Read Recruitment Caution", res.final_grade)

    def test_calibration_u377048_internal_clipping_with_verified_circularity(self):
        """Case 3: u377048 internal clipping does NOT break circularity (VERIFIED_CIRCULAR, STRUCTURALLY_HETEROGENEOUS)."""
        clipping_sites = [{
            "pos": 820500,
            "clipped_reads": 24,
            "depth": 27,
            "clipping_ratio": 0.889,
            "is_sv_polymorphism": False
        }]
        res = evaluate_four_gates(
            contig_len=1543357,
            spanning_unique_reads=28,
            support_fraction=0.82,
            clip_fraction=0.06,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=1,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=27.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=1.0,
            longest_repeat_bp=200,
            comp=97.5,
            cont=1.2,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=19,
            replicon_state="CHROMOSOME_SUPPORTED",
            gate5_status="PASS",
            clipping_sites_detail=clipping_sites
        )
        self.assertEqual(res.circularity.tier, "VERIFIED_CIRCULAR")
        self.assertEqual(res.circularity.evidence_type, "CONCORDANT_DUAL_OMICS")
        self.assertEqual(res.complexity.bio_tier, "STRUCTURALLY_HETEROGENEOUS")
        self.assertEqual(res.assembly_grade, "HIGH_CONFIDENCE_CIRCULAR")
        self.assertEqual(res.final_grade, "ORTHOGONALLY_VERIFIED_HIGH_QUALITY_CIRCULAR_CHROMOSOME_WITH_STRUCTURAL_HETEROGENEITY")

    def test_calibration_u333097_orthogonal_reinforcement(self):
        """Case 4: u333097 marginal long-read spanning reinforced by short reads (SUPPORTED_CIRCULAR, ORTHOGONAL_REINFORCEMENT)."""
        res = evaluate_four_gates(
            contig_len=450000,
            spanning_unique_reads=2,
            support_fraction=0.35,
            clip_fraction=0.12,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=0,
            partner_fraction=0.0,
            partner_name="None",
            anchor_status="UNIQUE",
            primary_mean_depth=8.5,
            primary_breadth_1x=99.8,
            repeat_breadth_pct=0.5,
            longest_repeat_bp=150,
            comp=None,
            cont=None,
            r5s=0,
            r16s=0,
            r23s=0,
            trna_aa_count=2,
            replicon_state="PLASMID_SUPPORTED",
            gate5_status="PASS",
            gate5_metrics={"Illumina_Rotated_Spanning_Reads_Unique": 15, "Illumina_Rotated_Bridging_Pairs_Unique": 22}
        )
        self.assertEqual(res.circularity.tier, "SUPPORTED_CIRCULAR")
        self.assertEqual(res.circularity.evidence_type, "ORTHOGONAL_REINFORCEMENT")
        self.assertEqual(res.complexity.bio_tier, "CLONAL_SIMPLE")
        self.assertEqual(res.assembly_grade, "HIGH_CONFIDENCE_CIRCULAR")
        self.assertEqual(res.final_grade, "ORTHOGONALLY_VERIFIED_CIRCULAR_PLASMID")

    def test_calibration_false_circular_negative_control(self):
        """Case 5: False circular negative control strictly rejected (REJECTED_LINEAR_OR_BROKEN, REJECTED_NON_CIRCULAR)."""
        res = evaluate_four_gates(
            contig_len=1800000,
            spanning_unique_reads=0,
            support_fraction=0.0,
            clip_fraction=0.0,
            unresolved_zero_bp=1200,
            has_high_risk_zero_block=True,
            clipping_sites_count=3,
            partner_reads=25,
            partner_fraction=0.85,
            partner_name="contig_linear_break",
            anchor_status="UNIQUE",
            primary_mean_depth=3.0,
            primary_breadth_1x=65.0,
            repeat_breadth_pct=0.0,
            longest_repeat_bp=0,
            comp=40.0,
            cont=5.0,
            r5s=0,
            r16s=0,
            r23s=0,
            trna_aa_count=4,
            replicon_state="AMBIGUOUS",
            gate5_status="JUNCTION_CONFLICT",
            partner_rel_depth_ratio=0.95
        )
        self.assertEqual(res.circularity.tier, "REJECTED_LINEAR_OR_BROKEN")
        self.assertEqual(res.assembly_grade, "REJECTED_NON_CIRCULAR")
        self.assertEqual(res.final_grade, "POTENTIAL_MISASSEMBLY (Physical Breakpoint / Unresolved Gap)")

    def test_calibration_ambiguous_middle_case(self):
        """Case 6: Ambiguous middle case with high-frequency competing partner (AMBIGUOUS_CIRCULAR, PUTATIVE_CIRCULAR)."""
        res = evaluate_four_gates(
            contig_len=1200000,
            spanning_unique_reads=12,
            support_fraction=0.65,
            clip_fraction=0.10,
            unresolved_zero_bp=0,
            has_high_risk_zero_block=False,
            clipping_sites_count=0,
            partner_reads=15,
            partner_fraction=0.55,
            partner_name="contig_competitor",
            anchor_status="UNIQUE",
            primary_mean_depth=25.0,
            primary_breadth_1x=100.0,
            repeat_breadth_pct=1.0,
            longest_repeat_bp=200,
            comp=92.0,
            cont=1.0,
            r5s=1,
            r16s=1,
            r23s=1,
            trna_aa_count=18,
            replicon_state="CHROMOSOME_SUPPORTED",
            gate5_status="NOT_APPLICABLE",
            partner_rel_depth_ratio=0.85
        )
        self.assertEqual(res.circularity.tier, "AMBIGUOUS_CIRCULAR")
        self.assertEqual(res.complexity.graph_tier, "COMPETING_GRAPH_BRANCH")
        self.assertEqual(res.complexity.rdr_signal, "POPULATION_COMPETITOR_SIGNAL")
        self.assertEqual(res.assembly_grade, "PUTATIVE_CIRCULAR")
        self.assertNotEqual(res.circularity.tier, "VERIFIED_CIRCULAR")

    def test_plot_circular_audit_u999856_evidence_chain(self):
        """Verify that u999856 generates a 4-panel evidence-chain figure without false breakpoints."""
        from scripts.plot_circular_audit import (
            extract_comprehensive_audit_events,
            extract_audit_events_from_metrics,
            load_audit_metrics,
            generate_verification_plot,
        )
        with tempfile.TemporaryDirectory() as td:
            tsv = Path(td) / "u999856_audit.tsv"
            tsv.write_text(
                "Contig_ID\tLength_bp\tEffective_Junction_Flank_bp\t"
                "Primary_Mean_Depth\tPrimary_Median_Depth\tPrimary_Depth_CV\t"
                "Primary_Coverage_Breadth_1x_Pct\tPrimary_Coverage_Breadth_3x_Pct\tTotal_Zero_Coverage_bp\t"
                "Junction_Spanning_Unique_Reads\tJunction_Informative_Unique_Reads\t"
                "Junction_Support_Fraction\tJunction_Support_CI_95_Lower\tJunction_Support_CI_95_Upper\t"
                "Junction_Clip_Fraction\tTerminal_Anchor_Uniqueness\t"
                "Systematic_Internal_Clipping_Sites\tInternal_Clipping_Sites_Detail\t"
                "Population_SV_Polymorphisms_Count\tPopulation_SV_Polymorphisms_Detail\t"
                "CheckM2_Completeness\tCheckM2_Contamination\t"
                "rRNA_5S_Count\trRNA_16S_Count\trRNA_23S_Count\ttRNA_Canonical_20_AA_Count\t"
                "GTDB_Lineage\tGate1_Physical_Structural_Support\tGate1_Notes\t"
                "Gate5_Illumina_Orthogonal_Status\tGate5_Notes\t"
                "Illumina_Relative_Abundance_Pct\tIllumina_Genome_Trimmed_Depth\tIllumina_Competitive_Depth_Ratio_C\t"
                "Illumina_Rotated_Junction_Depth_100bp\tIllumina_Rotated_Junction_Local_Ratio\t"
                "Illumina_Rotated_Spanning_Reads_Unique\tIllumina_Rotated_Bridging_Pairs_Unique\t"
                "Final_Confidence_Grade\n"
                "u999856\t2079472\t5000\t"
                "2483.0\t2650.0\t0.21\t"
                "100.0\t100.0\t0\t"
                "1895\t2281\t"
                "0.831\t0.815\t0.846\t"
                "0.026\tUNIQUE\t"
                "0\t1485361:858/102:841.2%\t"
                "1\t1485361:AF_del=28.1%(del_bp=4409)\t"
                "98.43\t1.13\t"
                "1\t1\t1\t18\t"
                "d__Archaea;p__Thermoproteota;c__Bathyarchaeia;o__Bathyarchaeales;f__Bathyarchaeaceae;g__Bathyarchaeum;s__Bathyarchaeum sp026014565\t"
                "PASS_WITH_SV_HETEROGENEITY\tPopulation_SV_Polymorphism(1_loci:pos=1485361(AF_del=28.1%))\t"
                "RECRUITMENT_CAUTION\tRecruitment_Drop(C=0.336);No_Junction_Span(Span=0,Bridge=0)\t"
                "72.366\t2863.0\t0.336\t"
                "161\t0.904\t"
                "0\t0\t"
                "HIGH_CONFIDENCE_CIRCULAR_WITH_POPULATION_COMPLEXITY\n"
            )
            metrics = load_audit_metrics(str(tsv), "u999856")
            clip_bps, sv_sites, zero_gaps = extract_comprehensive_audit_events(metrics)

            # Assert 0 unresolved breakpoints, exactly 1 population SV, and 0 zero gaps
            self.assertEqual(len(clip_bps), 0)
            self.assertEqual(len(sv_sites), 1)
            self.assertEqual(sv_sites[0]["pos"], 1485361)
            self.assertAlmostEqual(sv_sites[0]["af_del"], 0.281, places=2)
            self.assertEqual(sv_sites[0]["del_bp"], 4409)
            self.assertEqual(len(zero_gaps), 0)

            # Test generation of 4-panel figure
            fa = Path(td) / "u999856.fa"
            fa.write_text(">u999856\n" + "ACGT" * 50000 + "\n")  # 200 kb mock sequence
            out_pdf = Path(td) / "u999856_evidence_chain.pdf"
            out_png = Path(td) / "u999856_evidence_chain.png"
            generate_verification_plot(
                bin_id="u999856",
                contig_fa=str(fa),
                whole_bam=None,
                junction_bam=None,
                audit_tsv=str(tsv),
                bakta_tsv="",
                out_pdf=str(out_pdf),
                out_png=str(out_png),
                flank_len=5000,
                max_pileup_reads=30
            )
            self.assertTrue(out_pdf.exists())
            self.assertTrue(out_png.exists())
            self.assertGreater(out_pdf.stat().st_size, 1000)
            self.assertGreater(out_png.stat().st_size, 1000)

    def test_plot_circular_audit_clean_circular_genome(self):
        """Verify that a clean circular genome (0 SV, 0 breakpoints) generates an intact evidence-chain figure."""
        from scripts.plot_circular_audit import (
            extract_comprehensive_audit_events,
            load_audit_metrics,
            generate_verification_plot,
        )
        with tempfile.TemporaryDirectory() as td:
            tsv = Path(td) / "clean_sample_audit.tsv"
            tsv.write_text(
                "Contig_ID\tLength_bp\tEffective_Junction_Flank_bp\t"
                "Primary_Mean_Depth\tPrimary_Median_Depth\tPrimary_Depth_CV\t"
                "Primary_Coverage_Breadth_1x_Pct\tPrimary_Coverage_Breadth_3x_Pct\tTotal_Zero_Coverage_bp\t"
                "Junction_Spanning_Unique_Reads\tJunction_Informative_Unique_Reads\t"
                "Junction_Support_Fraction\tJunction_Support_CI_95_Lower\tJunction_Support_CI_95_Upper\t"
                "Junction_Clip_Fraction\tTerminal_Anchor_Uniqueness\t"
                "Systematic_Internal_Clipping_Sites\tInternal_Clipping_Sites_Detail\t"
                "Population_SV_Polymorphisms_Count\tPopulation_SV_Polymorphisms_Detail\t"
                "CheckM2_Completeness\tCheckM2_Contamination\t"
                "rRNA_5S_Count\trRNA_16S_Count\trRNA_23S_Count\ttRNA_Canonical_20_AA_Count\t"
                "GTDB_Lineage\tGate1_Physical_Structural_Support\tGate1_Notes\t"
                "Final_Confidence_Grade\n"
                "u123456\t1400000\t5000\t"
                "518.0\t520.0\t0.28\t"
                "100.0\t100.0\t0\t"
                "125\t180\t"
                "0.694\t0.623\t0.758\t"
                "0.015\tUNIQUE\t"
                "0\tNone\t"
                "0\tNone\t"
                "99.0\t0.5\t"
                "1\t1\t1\t20\t"
                "d__Bacteria\tPASS\tNone\t"
                "VERIFIED_CIRCULAR\n"
            )
            metrics = load_audit_metrics(str(tsv), "u123456")
            clip_bps, sv_sites, zero_gaps = extract_comprehensive_audit_events(metrics)

            # Assert 0 unresolved breakpoints, 0 population SVs, and 0 zero gaps
            self.assertEqual(len(clip_bps), 0)
            self.assertEqual(len(sv_sites), 0)
            self.assertEqual(len(zero_gaps), 0)

            # Generate 3-panel figure for 1.40 Mb genome
            fa = Path(td) / "u123456.fa"
            fa.write_text(">u123456\n" + "ACGT" * 35000 + "\n")
            out_pdf = Path(td) / "u123456_evidence_chain.pdf"
            out_png = Path(td) / "u123456_evidence_chain.png"
            generate_verification_plot(
                bin_id="u123456",
                contig_fa=str(fa),
                whole_bam=None,
                junction_bam=None,
                audit_tsv=str(tsv),
                bakta_tsv="",
                out_pdf=str(out_pdf),
                out_png=str(out_png),
                flank_len=5000,
                max_pileup_reads=30
            )
            self.assertTrue(out_pdf.exists())
            self.assertTrue(out_png.exists())
            self.assertGreater(out_pdf.stat().st_size, 1000)
            self.assertGreater(out_png.stat().st_size, 1000)

    def test_plot_circular_audit_bundle_workflow(self):
        """
        Verify the client-server decoupled workflow:
        1. Extract and serialize compact plot bundle (<100 KB JSON)
        2. Render publication-grade figure directly from bundle without BAMs
        3. Verify CLI commands (--dump-bundle and --bundle)
        """
        from scripts.plot_circular_audit import (
            extract_plot_bundle,
            save_plot_bundle,
            load_plot_bundle,
            render_plot_from_bundle,
        )

        with tempfile.TemporaryDirectory() as td:
            tsv = Path(td) / "audit.tsv"
            tsv.write_text(
                "Contig_ID\tLength_bp\tEffective_Junction_Flank_bp\t"
                "Primary_Mean_Depth\tPrimary_Median_Depth\tPrimary_Depth_CV\t"
                "Primary_Coverage_Breadth_1x_Pct\tPrimary_Coverage_Breadth_3x_Pct\tTotal_Zero_Coverage_bp\t"
                "Junction_Spanning_Unique_Reads\tJunction_Informative_Unique_Reads\t"
                "Junction_Support_Fraction\tJunction_Support_CI_95_Lower\tJunction_Support_CI_95_Upper\t"
                "Junction_Clip_Fraction\tTerminal_Anchor_Uniqueness\t"
                "Systematic_Internal_Clipping_Sites\tInternal_Clipping_Sites_Detail\t"
                "Population_SV_Polymorphisms_Count\tPopulation_SV_Polymorphisms_Detail\t"
                "CheckM2_Completeness\tCheckM2_Contamination\t"
                "rRNA_5S_Count\trRNA_16S_Count\trRNA_23S_Count\ttRNA_Canonical_20_AA_Count\t"
                "GTDB_Lineage\tGate1_Physical_Structural_Support\tGate1_Notes\t"
                "Final_Confidence_Grade\n"
                "u999856\t2079472\t5000\t"
                "153.2\t150.0\t0.21\t"
                "100.0\t100.0\t0\t"
                "1895\t2281\t"
                "0.831\t0.815\t0.846\t"
                "0.012\tUNIQUE\t"
                "0\tNone\t"
                "1\tpos=1485361:del_bp=4409:af_del=0.281:ref_reads=1230:del_reads=480\t"
                "98.43\t1.13\t"
                "1\t1\t1\t18\t"
                "d__Archaea\tPASS_WITH_SV_HETEROGENEITY\tNone\t"
                "VERIFIED_CIRCULAR\n"
            )
            fa = Path(td) / "u999856.fa"
            fa.write_text(">u999856\n" + "ATCG" * 50000 + "\n")

            # 1. Test extraction & serialization
            bundle = extract_plot_bundle(
                bin_id="u999856",
                contig_fa=str(fa),
                whole_bam=None,
                junction_bam=None,
                audit_tsv=str(tsv),
                bakta_tsv=None,
                flank_len=5000,
                max_pileup_reads=36
            )
            self.assertEqual(bundle["bin_id"], "u999856")
            self.assertEqual(bundle["contig_len"], 200000)
            self.assertEqual(len(bundle["population_svs"]), 1)
            self.assertAlmostEqual(bundle["population_svs"][0]["pos"], 1485361)

            bundle_json = Path(td) / "u999856_bundle.json"
            save_plot_bundle(bundle, str(bundle_json))
            self.assertTrue(bundle_json.exists())
            sz_kb = bundle_json.stat().st_size / 1024.0
            self.assertLess(sz_kb, 300.0)  # Must be lightweight <300 KB

            # 2. Test deserialization and rendering without BAMs
            loaded = load_plot_bundle(str(bundle_json))
            self.assertEqual(loaded["bin_id"], "u999856")
            out_pdf = Path(td) / "u999856_bundle_render.pdf"
            out_png = Path(td) / "u999856_bundle_render.png"
            render_plot_from_bundle(loaded, out_pdf=str(out_pdf), out_png=str(out_png))
            self.assertTrue(out_pdf.exists())
            self.assertTrue(out_png.exists())
            self.assertGreater(out_pdf.stat().st_size, 1000)
            self.assertGreater(out_png.stat().st_size, 1000)

            # 3. Test CLI --dump-bundle
            cli_bundle_json = Path(td) / "cli_dumped_bundle.json"
            script_path = str(PROJECT_ROOT / "scripts" / "plot_circular_audit.py")
            cmd_dump = [
                sys.executable, script_path,
                "--bin", "u999856",
                "--contig-fa", str(fa),
                "--audit-tsv", str(tsv),
                "--dump-bundle", str(cli_bundle_json)
            ]
            res_dump = subprocess.run(cmd_dump, capture_output=True, text=True)
            self.assertEqual(res_dump.returncode, 0, f"Error: {res_dump.stderr}")
            self.assertTrue(cli_bundle_json.exists())

            # 4. Test CLI --bundle rendering
            cli_pdf = Path(td) / "cli_from_bundle.pdf"
            cli_png = Path(td) / "cli_from_bundle.png"
            cmd_render = [
                sys.executable, script_path,
                "--bundle", str(cli_bundle_json),
                "--out-pdf", str(cli_pdf),
                "--out-png", str(cli_png)
            ]
            res_render = subprocess.run(cmd_render, capture_output=True, text=True)
            self.assertEqual(res_render.returncode, 0, f"Error: {res_render.stderr}")
            self.assertTrue(cli_pdf.exists())
            self.assertTrue(cli_png.exists())


if __name__ == "__main__":
    unittest.main()


