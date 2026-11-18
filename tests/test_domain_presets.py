#!/usr/bin/env python3
import os
import sys
import unittest
try:
    import yaml
except ImportError:
    yaml = None

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CONFIG_PATH = os.path.join(REPO_ROOT, "projects", "config_synteny", "synteny_validation.yaml")

class TestDomainPresets(unittest.TestCase):
    def setUp(self):
        if yaml is None:
            self.skipTest("PyYAML is not installed in current python environment")

    def test_config_yaml_has_domain_presets(self):
        self.assertTrue(os.path.exists(CONFIG_PATH), f"Config file not found: {CONFIG_PATH}")
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        
        syn_cfg = cfg.get("synteny_validation", {})
        self.assertIn("domain", syn_cfg, "domain key missing in synteny_validation")
        self.assertIn("domain_presets", syn_cfg, "domain_presets key missing in synteny_validation")
        
        presets = syn_cfg["domain_presets"]
        self.assertIn("archaea", presets, "archaea domain preset missing")
        self.assertIn("bacteria", presets, "bacteria domain preset missing")

    def test_archaea_preset_configuration(self):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        
        arch = cfg["synteny_validation"]["domain_presets"]["archaea"]
        self.assertIn("ar53_metadata_r232.tsv", arch.get("gtdb_metadata", ""))
        self.assertIn("archaea_all/genomes", arch.get("gtdb_genomes_dir", ""))
        self.assertIn("archaea_all/proteins", arch.get("gtdb_proteins_dir", ""))
        self.assertIn("protein_KO_all_strict.tsv", arch.get("gtdb_ko_table", ""))
        self.assertFalse(arch.get("only_representative", True))
        self.assertEqual(arch.get("default_preset"), "atp_synthase")
        self.assertEqual(arch.get("max_positive_refs"), 5)
        self.assertEqual(arch.get("max_negative_refs"), 2)
        self.assertFalse(arch.get("auto_fetch", True))

    def test_bacteria_preset_configuration(self):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        
        bac = cfg["synteny_validation"]["domain_presets"]["bacteria"]
        self.assertIn("bac120_metadata_r232.tsv", bac.get("gtdb_metadata", ""))
        self.assertIn("bacteria_all/genomes", bac.get("gtdb_genomes_dir", ""))
        self.assertIn("bacteria_all/proteins", bac.get("gtdb_proteins_dir", ""))
        self.assertTrue(bac.get("only_representative", False))
        self.assertEqual(bac.get("default_preset"), "sox_system")
        self.assertEqual(bac.get("max_positive_refs"), 10)
        self.assertEqual(bac.get("max_negative_refs"), 1)
        self.assertTrue(bac.get("auto_fetch", False))
        self.assertIn("MOB_abundance/data/MAGs", bac.get("target_genomes_dir", ""))
        self.assertIn("NS_MOB/proteins/all", bac.get("target_proteins_dir", ""))
        self.assertIn("NS_MOB/genes", bac.get("target_genes_dir", ""))
        self.assertIn("mag_classification", bac.get("target_classification", ""))

    def test_select_and_validate_refs_domain_arguments(self):
        sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
        import select_and_validate_refs
        
        # Test CLI parser accepts --domain and --auto_fetch / --no_auto_fetch
        parser = select_and_validate_refs.parse_args
        # Simulate CLI arguments
        test_args = [
            "--target_name", "test_sample",
            "--target_fa", "test.fa",
            "--target_taxonomy", "d__Archaea;p__Thermoproteota",
            "--gtdb_metadata", "dummy_meta.tsv",
            "--gtdb_genomes_dir", "dummy_dir",
            "--out_refs_tsv", "out.tsv",
            "--out_anchors_tsv", "anchors.tsv",
            "--out_dir", "out_dir",
            "--domain", "archaea",
            "--no_auto_fetch"
        ]
        import unittest.mock
        with unittest.mock.patch("sys.argv", ["select_and_validate_refs.py"] + test_args):
            args = select_and_validate_refs.parse_args()
            self.assertEqual(args.domain, "archaea")
            self.assertFalse(args.auto_fetch)

        test_args_bac = [
            "--target_name", "test_sample_bac",
            "--target_fa", "test.fa",
            "--target_taxonomy", "d__Bacteria;p__Pseudomonadota",
            "--gtdb_metadata", "dummy_meta.tsv",
            "--gtdb_genomes_dir", "dummy_dir",
            "--out_refs_tsv", "out.tsv",
            "--out_anchors_tsv", "anchors.tsv",
            "--out_dir", "out_dir",
            "--domain", "bacteria",
            "--auto_fetch"
        ]
        with unittest.mock.patch("sys.argv", ["select_and_validate_refs.py"] + test_args_bac):
            args = select_and_validate_refs.parse_args()
            self.assertEqual(args.domain, "bacteria")
            self.assertTrue(args.auto_fetch)

if __name__ == "__main__":
    unittest.main()
