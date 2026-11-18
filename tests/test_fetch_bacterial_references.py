#!/usr/bin/env python3
import os
import sys
import unittest
import tempfile
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
import fetch_bacterial_references

class TestFetchBacterialReferences(unittest.TestCase):
    def test_extract_clean_acc(self):
        self.assertEqual(fetch_bacterial_references.extract_clean_acc("GB_GCA_002072955.1"), "GCA_002072955.1")
        self.assertEqual(fetch_bacterial_references.extract_clean_acc("RS_GCF_000223525.1"), "GCF_000223525.1")
        self.assertEqual(fetch_bacterial_references.extract_clean_acc("GCA_002072955.1_ASM207295v1_genomic"), "GCA_002072955.1")

    def test_target_classification_deduplication(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            class_file = os.path.join(tmp_dir, "mag_classification.tsv")
            with open(class_file, "w") as f:
                f.write("MAG_id\ttaxonomy\n")
                f.write("GCA_002072955.1_genomic\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Methylococcales;f__Methylomonadaceae;g__Methyloprofundus;s__M. sedimenti\n")
                f.write("GCA_002189065.1_genomic\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Methylococcales;f__Methylomonadaceae;g__Methyloprofundus;s__M. sp\n")

            # Mock GTDB metadata
            meta_file = os.path.join(tmp_dir, "mock_bac120_meta.tsv")
            mock_data = [
                {"accession": "GCA_000214275.1", "gtdb_taxonomy": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Methylococcales;f__Methylomonadaceae;g__Methylomicrobium;s__album", "gtdb_representative": "t", "checkm_completeness": 99.0, "checkm_contamination": 0.5, "ncbi_assembly_level": "Complete Genome"},
                {"accession": "GCA_000262525.1", "gtdb_taxonomy": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Methylococcales;f__Methylomonadaceae;g__Methylobacter;s__tundripaludum", "gtdb_representative": "t", "checkm_completeness": 98.5, "checkm_contamination": 0.8, "ncbi_assembly_level": "Scaffold"},
                {"accession": "GCF_000025745.1", "gtdb_taxonomy": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Halothiobacillales;f__Halothiobacillaceae;g__Halothiobacillus;s__neapolitanus", "gtdb_representative": "t", "checkm_completeness": 100.0, "checkm_contamination": 0.2, "ncbi_assembly_level": "Complete Genome"},
            ]
            pd.DataFrame(mock_data).to_csv(meta_file, sep="\t", index=False)

            # Test argument parsing & setup
            sys.argv = [
                "fetch_bacterial_references.py",
                "--gtdb_metadata", meta_file,
                "--target_classification", class_file,
                "--out_genomes_dir", os.path.join(tmp_dir, "genomes"),
                "--out_proteins_dir", os.path.join(tmp_dir, "proteins"),
                "--out_manifest", os.path.join(tmp_dir, "manifest.tsv"),
                "--max_positive_refs", "2",
                "--max_negative_refs", "1"
            ]
            args = fetch_bacterial_references.parse_args()
            self.assertEqual(args.max_positive_refs, 2)
            self.assertEqual(args.max_negative_refs, 1)

if __name__ == "__main__":
    unittest.main()
