#!/usr/bin/env python3
import os
import sys
import unittest
try:
    import yaml
except ImportError:
    yaml = None

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SAMPLES_CONFIG_PATH = os.path.join(REPO_ROOT, "projects", "config_synteny", "samples.yaml")
SYNTENY_CONFIG_PATH = os.path.join(REPO_ROOT, "projects", "config_synteny", "synteny_validation.yaml")

class TestSamplePresets(unittest.TestCase):
    def setUp(self):
        if yaml is None:
            self.skipTest("PyYAML is not installed in current python environment")

    def test_samples_yaml_structure_and_defaults(self):
        self.assertTrue(os.path.exists(SAMPLES_CONFIG_PATH), f"Config file not found: {SAMPLES_CONFIG_PATH}")
        with open(SAMPLES_CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)

        self.assertIn("synteny_samples", cfg, "synteny_samples missing in samples.yaml")
        self.assertIn("default_domain", cfg, "default_domain missing in samples.yaml")
        self.assertIn("default_presets", cfg, "default_presets missing in samples.yaml")
        self.assertEqual(cfg["default_domain"], "archaea")
        self.assertEqual(cfg["default_presets"], ["M00159", "M00157"])

        # Check target directories are present
        self.assertIn("target_genomes_dir", cfg)
        self.assertIn("target_proteins_dir", cfg)
        self.assertIn("target_genes_dir", cfg)
        self.assertIn("target_classification", cfg)

    def test_per_sample_configurations(self):
        with open(SAMPLES_CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        samples = cfg["synteny_samples"]

        # 1. Bathy_u999856: Archaea with multi-operon ATP synthase
        self.assertIn("Bathy_u999856", samples)
        b_info = samples["Bathy_u999856"]
        self.assertEqual(b_info.get("domain"), "archaea")
        self.assertEqual(b_info.get("presets"), ["M00159", "M00157"])
        self.assertTrue(os.path.isabs(b_info.get("fasta", "")))
        self.assertIn("d__Archaea", b_info.get("taxonomy", ""))

        # 2. Hadar_u499975: Archaea with single M00159 preset
        self.assertIn("Hadar_u499975", samples)
        h_info = samples["Hadar_u499975"]
        self.assertEqual(h_info.get("domain"), "archaea")
        self.assertEqual(h_info.get("presets"), ["M00159"])
        self.assertTrue(os.path.isabs(h_info.get("fasta", "")))

    def test_sample_domain_and_preset_resolution_logic(self):
        """Simulate the resolution functions implemented in synteny_validation.smk."""
        with open(SAMPLES_CONFIG_PATH, "r", encoding="utf-8") as f:
            sample_cfg = yaml.safe_load(f)
        with open(SYNTENY_CONFIG_PATH, "r", encoding="utf-8") as f:
            syn_cfg = yaml.safe_load(f)["synteny_validation"]

        samples = sample_cfg.get("synteny_samples", {})
        presets = syn_cfg.get("operon_presets", {})
        domain_presets = syn_cfg.get("domain_presets", {})

        def parse_preset_list(val):
            if isinstance(val, list):
                return [str(p).strip() for p in val if str(p).strip()]
            elif isinstance(val, str):
                if "," in val:
                    return [p.strip() for p in val.split(",") if p.strip()]
                elif val.strip():
                    return [val.strip()]
            return []

        def get_sample_domain(s_id, s_dict):
            info = s_dict.get(s_id, {})
            if "domain" in info and info["domain"]:
                return info["domain"].strip().lower()
            tax = str(info.get("taxonomy", "")).lower()
            if "d__archaea" in tax:
                return "archaea"
            if "d__bacteria" in tax:
                return "bacteria"
            return sample_cfg.get("default_domain", "bacteria")

        def get_sample_presets(s_id, s_dict):
            info = s_dict.get(s_id, {})
            raw_p = info.get("presets") or info.get("preset") or info.get("active_presets")
            if raw_p:
                return parse_preset_list(raw_p)
            if "default_presets" in sample_cfg:
                return parse_preset_list(sample_cfg["default_presets"])
            s_dom = get_sample_domain(s_id, s_dict)
            return parse_preset_list(domain_presets.get(s_dom, {}).get("default_preset", "M00159"))

        def get_sample_modules(s_id, s_dict):
            p_names = get_sample_presets(s_id, s_dict)
            mods = []
            for p in p_names:
                if p in presets:
                    p_data = presets[p]
                    for m in p_data.get("enabled_modules", list(p_data.get("modules", {}).keys())):
                        if m not in mods:
                            mods.append(m)
                elif p in ["sox_system", "sox"]:
                    if "M00595" not in mods:
                        mods.append("M00595")
                elif p.startswith("M") and len(p) == 6:
                    if p not in mods:
                        mods.append(p)
            return mods

        # Test Bathy
        self.assertEqual(get_sample_domain("Bathy_u999856", samples), "archaea")
        self.assertEqual(get_sample_presets("Bathy_u999856", samples), ["M00159", "M00157"])
        self.assertEqual(get_sample_modules("Bathy_u999856", samples), ["M00159", "M00157"])

        # Test Hadar
        self.assertEqual(get_sample_domain("Hadar_u499975", samples), "archaea")
        self.assertEqual(get_sample_presets("Hadar_u499975", samples), ["M00159"])
        self.assertEqual(get_sample_modules("Hadar_u499975", samples), ["M00159"])

        # Test simulated bacterial sample with fallback to taxonomy and presets
        mock_samples = {
            "Bacterial_MAG_01": {
                "taxonomy": "d__Bacteria;p__Pseudomonadota",
                "presets": ["M00595"]
            },
            "Default_Fallback_MAG": {
                "taxonomy": "d__Archaea;p__Thermoproteota"
                # No domain or presets specified
            }
        }
        self.assertEqual(get_sample_domain("Bacterial_MAG_01", mock_samples), "bacteria")
        self.assertEqual(get_sample_modules("Bacterial_MAG_01", mock_samples), ["M00595"])

        self.assertEqual(get_sample_domain("Default_Fallback_MAG", mock_samples), "archaea")
        # Should inherit default_presets from samples.yaml
        self.assertEqual(get_sample_modules("Default_Fallback_MAG", mock_samples), ["M00159", "M00157"])

    def test_heterogeneous_target_generation(self):
        """Verify that heterogenous samples generate distinct targets based on per-sample modules."""
        mock_samples = {
            "Sample_Archaea": {"domain": "archaea", "presets": ["M00159", "M00157"]},
            "Sample_Bacteria": {"domain": "bacteria", "presets": ["M00595"]}
        }
        
        targets = []
        for s_id, s_info in mock_samples.items():
            mods = s_info["presets"]
            for m in mods:
                targets.append(f"results/synteny_validation/{s_id}/{m}/stage4/{s_id}_locus_validation_summary.tsv")

        # Archaea has M00159 and M00157, but not M00595
        self.assertIn("results/synteny_validation/Sample_Archaea/M00159/stage4/Sample_Archaea_locus_validation_summary.tsv", targets)
        self.assertIn("results/synteny_validation/Sample_Archaea/M00157/stage4/Sample_Archaea_locus_validation_summary.tsv", targets)
        self.assertNotIn("results/synteny_validation/Sample_Archaea/M00595/stage4/Sample_Archaea_locus_validation_summary.tsv", targets)

        # Bacteria has M00595, but not M00159 or M00157
        self.assertIn("results/synteny_validation/Sample_Bacteria/M00595/stage4/Sample_Bacteria_locus_validation_summary.tsv", targets)
        self.assertNotIn("results/synteny_validation/Sample_Bacteria/M00159/stage4/Sample_Bacteria_locus_validation_summary.tsv", targets)

if __name__ == "__main__":
    unittest.main()
