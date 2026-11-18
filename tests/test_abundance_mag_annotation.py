#!/usr/bin/env python3
"""
tests/test_abundance_mag_annotation.py

Unit tests for MAG abundance pipeline integrated with mag_annotation:
1. scripts/combine_mag_reference.py: multi-MAG FASTA + GFF3 synchronization.
2. scripts/calculate_tpm.py: MAG extraction and TPM calculation.
3. scripts/calculate_tpm_ratio.py: RNA/DNA ratio calculation.
4. scripts/merge_abundance.py: decimal sample names support (e.g. C1_0-0.5).
5. scripts/summarize_abundance_by_gene.py: master gene x abundance x annotation table.
"""

import os
import shutil
import sys
import tempfile
import unittest
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

from scripts import combine_mag_reference
from scripts import calculate_tpm
from scripts import calculate_tpm_ratio
from scripts import merge_abundance
from scripts import summarize_abundance_by_gene


class TestAbundanceMagAnnotation(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_abun_mag_")

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_01_combine_mag_reference(self):
        """Test building synchronized combined FASTA and GFF3."""
        ref_dir = os.path.join(self.test_dir, "mags")
        genes_dir = os.path.join(self.test_dir, "genes")
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(genes_dir, exist_ok=True)

        # Create mock genomes
        mag1_fa = os.path.join(ref_dir, "MAG_01.fa")
        with open(mag1_fa, "w") as f:
            f.write(">contig_A\nATGCGATCGATCGATCGA\n>contig_B\nATGAAACCCGGGTTTAAA\n")

        mag2_fa = os.path.join(ref_dir, "MAG_02.fa")
        with open(mag2_fa, "w") as f:
            f.write(">contig_A\nATGTTTCCCGGGAAATTT\n>contig_C\nATGGCGGCGGCGGCGGCG\n")

        # Create mock GFF for MAG_01 and MAG_02 (from mag_annotation)
        mag1_gff = os.path.join(genes_dir, "MAG_01.gff")
        with open(mag1_gff, "w") as f:
            f.write("##gff-version 3\n")
            f.write("contig_A\tpyrodigal\tCDS\t1\t18\t.\t+\t0\tID=MAG_01___contig_A___gene00001;protein_id=MAG_01___contig_A___prot00001\n")
            f.write("contig_B\tpyrodigal\tCDS\t1\t18\t.\t+\t0\tID=MAG_01___contig_B___gene00001;protein_id=MAG_01___contig_B___prot00001\n")

        mag2_gff = os.path.join(genes_dir, "MAG_02.gff")
        with open(mag2_gff, "w") as f:
            f.write("##gff-version 3\n")
            f.write("contig_A\tpyrodigal\tCDS\t1\t18\t.\t+\t0\tID=MAG_02___contig_A___gene00001;protein_id=MAG_02___contig_A___prot00001\n")
            f.write("contig_C\tpyrodigal\tCDS\t1\t18\t.\t+\t0\tID=MAG_02___contig_C___gene00001;protein_id=MAG_02___contig_C___prot00001\n")

        out_fa = os.path.join(self.test_dir, "combined.fa")
        out_gff = os.path.join(self.test_dir, "combined.gff")

        sys.argv = [
            "combine_mag_reference.py",
            "--ref-dir", ref_dir,
            "--genes-dir", genes_dir,
            "--out-fa", out_fa,
            "--out-gff", out_gff,
            "--mag-ext", "fa"
        ]
        combine_mag_reference.main()

        self.assertTrue(os.path.exists(out_fa))
        self.assertTrue(os.path.exists(out_gff))

        # Check FASTA headers: duplicate 'contig_A' was disambiguated with prefix
        headers = []
        with open(out_fa) as f:
            for line in f:
                if line.startswith(">"):
                    headers.append(line.strip()[1:])

        self.assertIn("MAG_01___contig_A", headers)
        self.assertIn("MAG_02___contig_A", headers)

        # Check GFF lines match FASTA headers
        gff_seqnames = []
        gff_gene_ids = []
        with open(out_gff) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                gff_seqnames.append(parts[0])
                self.assertIn("ID=", parts[8])
                self.assertIn("locus_tag=", parts[8])
                for item in parts[8].split(";"):
                    if item.startswith("ID="):
                        gff_gene_ids.append(item.split("=")[1])

        self.assertIn("MAG_01___contig_A", gff_seqnames)
        self.assertIn("MAG_01___contig_A___gene00001", gff_gene_ids)
        self.assertIn("MAG_02___contig_A", gff_seqnames)
        self.assertIn("MAG_02___contig_A___gene00001", gff_gene_ids)

    def test_combine_mag_reference_missing_gff_error(self):
        """Test that missing mag_annotation GFF strictly raises FileNotFoundError (no re-prediction)."""
        ref_dir = os.path.join(self.test_dir, "missing_ref")
        genes_dir = os.path.join(self.test_dir, "empty_genes")
        os.makedirs(ref_dir, exist_ok=True)
        os.makedirs(genes_dir, exist_ok=True)

        mag_fa = os.path.join(ref_dir, "MAG_MISSING.fa")
        with open(mag_fa, "w") as f:
            f.write(">contig_1\nATGCGATCGATCGATCGAAATTTCCCGGG\n")

        out_fa = os.path.join(self.test_dir, "fail_combined.fa")
        out_gff = os.path.join(self.test_dir, "fail_combined.gff")

        sys.argv = [
            "combine_mag_reference.py",
            "--ref-dir", ref_dir,
            "--genes-dir", genes_dir,
            "--out-fa", out_fa,
            "--out-gff", out_gff,
            "--mag-ext", "fa"
        ]
        with self.assertRaises(FileNotFoundError):
            combine_mag_reference.main()

    def test_02_calculate_tpm_with_mag_extraction(self):
        """Test featureCounts parsing and MAG extraction in calculate_tpm.py."""
        counts_tsv = os.path.join(self.test_dir, "sample.featureCounts.tsv")
        summary_tsv = os.path.join(self.test_dir, "sample.featureCounts.tsv.summary")
        final_tsv = os.path.join(self.test_dir, "sample.final.tsv")

        with open(summary_tsv, "w") as f:
            f.write("Status\tsample.bam\nAssigned\t1000\nUnassigned_Unmapped\t500\n")

        with open(counts_tsv, "w") as f:
            f.write("# Program:featureCounts\n")
            f.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam\n")
            f.write("MAG_01___contig_A___gene00001\tcontig_A\t1\t1000\t+\t1000\t100\n")
            f.write("MAG_02___contig_B___gene00001\tcontig_B\t1\t500\t-\t500\t50\n")

        sys.argv = [
            "calculate_tpm.py",
            "--counts", counts_tsv,
            "--summary", summary_tsv,
            "--output", final_tsv
        ]
        calculate_tpm.main()

        self.assertTrue(os.path.exists(final_tsv))
        df = pd.read_csv(final_tsv, sep="\t")
        self.assertIn("MAG", df.columns)
        self.assertEqual(df.loc[0, "MAG"], "MAG_01")
        self.assertEqual(df.loc[1, "MAG"], "MAG_02")
        self.assertIn("TPM_Lib", df.columns)
        self.assertGreater(df.loc[0, "TPM_Lib"], 0)

    def test_03_merge_abundance_decimal_samples(self):
        """Test merge_abundance.py preserves sample names with decimals (e.g. C1_0-0.5)."""
        s1_file = os.path.join(self.test_dir, "C1_0-0.5.final.tsv")
        s2_file = os.path.join(self.test_dir, "C2_3-3.5.final.tsv")

        df1 = pd.DataFrame({
            "MAG": ["MAG1", "MAG2"],
            "LocalGene": ["G1", "G2"],
            "Chr": ["chr1", "chr2"],
            "Contig": ["chr1", "chr2"],
            "Start": [1, 1],
            "End": [100, 200],
            "Strand": ["+", "-"],
            "Length": [100, 200],
            "Counts": [10, 20],
            "TPM_Lib": [50.0, 100.0]
        })
        df2 = pd.DataFrame({
            "MAG": ["MAG1", "MAG2"],
            "LocalGene": ["G1", "G2"],
            "Chr": ["chr1", "chr2"],
            "Contig": ["chr1", "chr2"],
            "Start": [1, 1],
            "End": [100, 200],
            "Strand": ["+", "-"],
            "Length": [100, 200],
            "Counts": [15, 25],
            "TPM_Lib": [75.0, 125.0]
        })
        df1.to_csv(s1_file, sep="\t", index=False)
        df2.to_csv(s2_file, sep="\t", index=False)

        merged_out = os.path.join(self.test_dir, "merged_dna_tpm.tsv")
        sys.argv = [
            "merge_abundance.py",
            "--type", "dna",
            "--inputs", s1_file, s2_file,
            "--output", merged_out
        ]
        merge_abundance.main()

        self.assertTrue(os.path.exists(merged_out))
        merged_df = pd.read_csv(merged_out, sep="\t")
        # Ensure sample column headers were NOT truncated to C1_0-0 or C2_3-3
        self.assertIn("C1_0-0.5_Counts", merged_df.columns)
        self.assertIn("C2_3-3.5_Counts", merged_df.columns)
        self.assertIn("C1_0-0.5_TPM_Lib", merged_df.columns)
        self.assertIn("C2_3-3.5_TPM_Lib", merged_df.columns)

    def test_04_summarize_abundance_master_table(self):
        """Test master table generation: each gene is one row, annotations and abundances are columns."""
        # Create mock ratio file
        ratio_file = os.path.join(self.test_dir, "C1_0-0.5.ratio.tsv")
        ratio_df = pd.DataFrame({
            "MAG": ["MAG_01", "MAG_01"],
            "LocalGene": ["MAG_01___ctg1___gene00001", "MAG_01___ctg1___gene00002"],
            "Chr": ["ctg1", "ctg1"],
            "Start": [1, 501],
            "End": [500, 1000],
            "Strand": ["+", "+"],
            "Length": [500, 500],
            "DNA_Counts": [100, 200],
            "DNA_TPM_Lib": [25.0, 50.0],
            "RNA_Counts": [50, 150],
            "RNA_TPM_Lib": [12.5, 37.5],
            "Ratio": [0.5, 0.75]
        })
        ratio_df.to_csv(ratio_file, sep="\t", index=False)

        # Create mock annotation file (protein_KO_best.tsv)
        anno_file = os.path.join(self.test_dir, "protein_KO_best.tsv")
        anno_df = pd.DataFrame({
            "protein_id": ["MAG_01___ctg1___prot00001", "MAG_01___ctg1___prot00002"],
            "gene_id": ["MAG_01___ctg1___gene00001", "MAG_01___ctg1___gene00002"],
            "genome_id": ["MAG_01", "MAG_01"],
            "ko": ["K03531", "K03590"],
            "ko_definition": ["cell division protein FtsZ", "cell division protein FtsA"],
            "bitscore": [450.0, 510.0],
            "evalue": [1e-100, 1e-120]
        })
        anno_df.to_csv(anno_file, sep="\t", index=False)

        master_out = os.path.join(self.test_dir, "gene_abundance_and_annotation.tsv")
        sys.argv = [
            "summarize_abundance_by_gene.py",
            "--merged-ratio", ratio_file,
            "--annotation-tsv", anno_file,
            "--output", master_out
        ]
        summarize_abundance_by_gene.main()

        self.assertTrue(os.path.exists(master_out))
        master_df = pd.read_csv(master_out, sep="\t")

        # Verification of format: Each gene one row
        self.assertEqual(len(master_df), 2)
        self.assertIn("gene_id", master_df.columns)
        self.assertIn("KO", master_df.columns)
        self.assertIn("KO_Definition", master_df.columns)
        self.assertIn("C1_0-0.5_DNA_Counts", master_df.columns)
        self.assertIn("C1_0-0.5_DNA_TPM", master_df.columns)
        self.assertIn("C1_0-0.5_RNA_Counts", master_df.columns)
        self.assertIn("C1_0-0.5_RNA_TPM", master_df.columns)
        self.assertIn("C1_0-0.5_Ratio", master_df.columns)
        self.assertIn("Mean_Ratio", master_df.columns)

        self.assertEqual(master_df.loc[0, "KO"], "K03531")
        self.assertEqual(master_df.loc[0, "KO_Definition"], "cell division protein FtsZ")
        self.assertEqual(master_df.loc[0, "C1_0-0.5_Ratio"], 0.5)

    def test_05_summarize_abundance_from_merged_file(self):
        """Test master table generation when input is a multi-sample merged_ratio.tsv file."""
        merged_file = os.path.join(self.test_dir, "merged_ratio.tsv")
        merged_df = pd.DataFrame({
            "MAG": ["MAG_01", "MAG_01"],
            "LocalGene": ["MAG_01___ctg1___gene00001", "MAG_01___ctg1___gene00002"],
            "Chr": ["ctg1", "ctg1"],
            "Contig": ["ctg1", "ctg1"],
            "Start": [1, 501],
            "End": [500, 1000],
            "Strand": ["+", "+"],
            "Length": [500, 500],
            "C1_0-0.5_Ratio": [0.5, 0.75],
            "C1_0-0.5_DNA_Counts": [100, 200],
            "C1_0-0.5_DNA_TPM_Lib": [25.0, 50.0],
            "C1_0-0.5_RNA_Counts": [50, 150],
            "C1_0-0.5_RNA_TPM_Lib": [12.5, 37.5],
            "C2_3-3.5_Ratio": [1.2, 0.8],
            "C2_3-3.5_DNA_Counts": [120, 220],
            "C2_3-3.5_DNA_TPM_Lib": [30.0, 55.0],
            "C2_3-3.5_RNA_Counts": [60, 160],
            "C2_3-3.5_RNA_TPM_Lib": [36.0, 44.0],
        })
        merged_df.to_csv(merged_file, sep="\t", index=False)

        anno_file = os.path.join(self.test_dir, "protein_KO_best.tsv")
        anno_df = pd.DataFrame({
            "protein_id": ["MAG_01___ctg1___prot00001", "MAG_01___ctg1___prot00002"],
            "gene_id": ["MAG_01___ctg1___gene00001", "MAG_01___ctg1___gene00002"],
            "genome_id": ["MAG_01", "MAG_01"],
            "ko": ["K03531", "K03590"],
            "ko_definition": ["cell division protein FtsZ", "cell division protein FtsA"],
            "bitscore": [450.0, 510.0],
            "evalue": [1e-100, 1e-120]
        })
        anno_df.to_csv(anno_file, sep="\t", index=False)

        master_out = os.path.join(self.test_dir, "master_from_merged.tsv")
        sys.argv = [
            "summarize_abundance_by_gene.py",
            "--merged-ratio", merged_file,
            "--annotation-tsv", anno_file,
            "--output", master_out
        ]
        summarize_abundance_by_gene.main()

        self.assertTrue(os.path.exists(master_out))
        res_df = pd.read_csv(master_out, sep="\t")

        # Crucial assertions: Samples are separate columns, NOT collapsed into merged_
        self.assertNotIn("merged_DNA_TPM", res_df.columns)
        self.assertIn("C1_0-0.5_DNA_TPM", res_df.columns)
        self.assertIn("C2_3-3.5_DNA_TPM", res_df.columns)
        self.assertIn("C1_0-0.5_RNA_TPM", res_df.columns)
        self.assertIn("C2_3-3.5_RNA_TPM", res_df.columns)
        self.assertIn("C1_0-0.5_Ratio", res_df.columns)
        self.assertIn("C2_3-3.5_Ratio", res_df.columns)

        # Values are accurate non-zeros
        self.assertEqual(res_df.loc[0, "C1_0-0.5_DNA_TPM"], 25.0)
        self.assertEqual(res_df.loc[0, "C2_3-3.5_DNA_TPM"], 30.0)
        self.assertEqual(res_df.loc[0, "C1_0-0.5_Ratio"], 0.5)
        self.assertEqual(res_df.loc[0, "C2_3-3.5_Ratio"], 1.2)


if __name__ == "__main__":
    unittest.main()
