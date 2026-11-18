#!/usr/bin/env python3
import unittest
import os
import sys
from unittest.mock import MagicMock

for mod in ["Bio", "Bio.Seq", "Bio.SeqIO", "Bio.SeqFeature", "Bio.SeqRecord", "pandas"]:
    if mod not in sys.modules:
        sys.modules[mod] = MagicMock()

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scripts.extract_anchor_loci import register_genome_alias_keys, lookup_ref_flanking_anno

class TestRefAnnotationBridging(unittest.TestCase):
    def test_bridging_prodigal_to_gtdb(self):
        # Simulated hit info from GTDB protein_KO_all_strict.tsv
        hit_info = {
            "ko": "K15888",
            "definition": "putative membrane protein",
            "score": "329.5",
            "threshold": "202.9",
            "evalue": "9.2e-97",
            "rank": "1"
        }

        # The raw protein ID in GTDB KO table
        raw_gtdb_pid = "GCA_963760905.1_SRP231805_co4_129_genomic___CAZOPV010000010.1___prot00018"
        genome_submap = {}

        register_genome_alias_keys(raw_gtdb_pid, hit_info, genome_submap)

        anno_map = {
            "GCA_963760905.1": genome_submap,
            "__by_genome__": {
                "GCA_963760905.1": genome_submap
            }
        }

        # Query 1: Prodigal style ID from ref_faa (CAZOPV010000010.1_18)
        hit1 = lookup_ref_flanking_anno("CAZOPV010000010.1_18", anno_map, genome_acc="GCA_963760905.1")
        self.assertIsNotNone(hit1)
        self.assertEqual(hit1["ko"], "K15888")

        # Query 2: Simple suffix (prot00018)
        hit2 = lookup_ref_flanking_anno("prot00018", anno_map, genome_acc="GCA_963760905.1")
        self.assertIsNotNone(hit2)
        self.assertEqual(hit2["ko"], "K15888")

        # Query 3: Short integer (_18 or 18)
        hit3 = lookup_ref_flanking_anno("_18", anno_map, genome_acc="GCA_963760905.1")
        self.assertIsNotNone(hit3)
        self.assertEqual(hit3["ko"], "K15888")

        # Query 4: Full GTDB ID
        hit4 = lookup_ref_flanking_anno(raw_gtdb_pid, anno_map, genome_acc="GCA_963760905.1")
        self.assertIsNotNone(hit4)
        self.assertEqual(hit4["ko"], "K15888")

if __name__ == "__main__":
    unittest.main()
