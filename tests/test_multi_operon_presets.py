#!/usr/bin/env python3
import os
import sys
import unittest
try:
    import yaml
except ImportError:
    yaml = None
import pandas as pd
from io import StringIO

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CONFIG_PATH = os.path.join(REPO_ROOT, "projects", "config_synteny", "synteny_validation.yaml")
M00595_YAML = os.path.join(REPO_ROOT, "knowledge", "kegg", "modules", "M00595.yaml")

class TestMultiOperonPresets(unittest.TestCase):
    def setUp(self):
        if yaml is None:
            self.skipTest("PyYAML is not installed in current python environment")

    def test_config_yaml_has_presets(self):
        self.assertTrue(os.path.exists(CONFIG_PATH), f"Config file not found: {CONFIG_PATH}")
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        
        syn_cfg = cfg.get("synteny_validation", {})
        self.assertIn("operon_presets", syn_cfg, "operon_presets missing in synteny_validation.yaml")
        presets = syn_cfg["operon_presets"]
        
        # Verify both presets exist
        self.assertIn("sox_system", presets, "sox_system preset missing")
        self.assertIn("atp_synthase", presets, "atp_synthase preset missing")
        
        # Verify active_preset is set
        active_preset = syn_cfg.get("active_preset")
        self.assertIn(active_preset, presets, f"active_preset {active_preset} not in operon_presets")
        
    def test_sox_preset_specification(self):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        sox = cfg["synteny_validation"]["operon_presets"]["sox_system"]
        self.assertIn("M00595", sox["enabled_modules"])
        m_def = sox["modules"]["M00595"]
        expected_kos = {"K17224", "K17222", "K17223", "K17227", "K17226", "K17225", "K22622"}
        self.assertEqual(set(m_def["target_kos"]), expected_kos)
        self.assertGreater(m_def.get("expected_locus_length", 0), 0)

    def test_atp_synthase_preset_specification(self):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        atp = cfg["synteny_validation"]["operon_presets"]["atp_synthase"]
        self.assertIn("M00159", atp["enabled_modules"])
        for mod in ["M00159", "M00157", "M00158", "M00160"]:
            self.assertIn(mod, atp["modules"])
            self.assertTrue(len(atp["modules"][mod]["target_kos"]) > 0)

    def test_m00595_kegg_yaml(self):
        self.assertTrue(os.path.exists(M00595_YAML), f"M00595.yaml not found: {M00595_YAML}")
        with open(M00595_YAML, "r", encoding="utf-8") as f:
            kegg_mod = yaml.safe_load(f)
        self.assertEqual(kegg_mod.get("module_id"), "M00595")
        expected_kos = {"K17224", "K17222", "K17223", "K17227", "K17226", "K17225", "K22622"}
        self.assertEqual(set(kegg_mod.get("key_kos", [])), expected_kos)

    def test_scripts_registered_m00595(self):
        sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
        import extract_anchor_loci
        import rescue_and_classify_gene_state
        
        self.assertIn("M00595", extract_anchor_loci.MODULE_TARGET_KOS)
        self.assertIn("M00595", rescue_and_classify_gene_state.MODULE_TARGET_KOS)
        
        # Verify Sox KOs in extract_anchor_loci
        sox_kos = extract_anchor_loci.MODULE_TARGET_KOS["M00595"]
        self.assertIn("K17224", sox_kos)
        self.assertIn("K17222", sox_kos)
        self.assertIn("K22622", sox_kos)

    def test_compile_matrix_with_sox(self):
        sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
        import compile_atp_synteny_matrix
        
        # Test determine_overall_atp_state for a SOX module
        dummy_sox_dict = {
            "M00595": {
                "loss": "INFERRED_LOSS_HIGH",
                "cluster": "FULLY_ABSENT",
                "conf": "HIGH",
                "fate": "CLEAN_LOSS"
            }
        }
        res = compile_atp_synteny_matrix.determine_overall_atp_state(dummy_sox_dict)
        self.assertEqual(res, "CLEAN_LOSS")

    def test_individual_module_presets(self):
        """Verify each ATP synthase and SOX module exists as an independent preset named by module ID."""
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        presets = cfg["synteny_validation"]["operon_presets"]
        
        expected_modules = ["M00159", "M00157", "M00158", "M00160", "M00595"]
        for mod in expected_modules:
            self.assertIn(mod, presets, f"Individual preset {mod} missing in operon_presets")
            p_data = presets[mod]
            self.assertIn(mod, p_data.get("enabled_modules", []))
            self.assertIn(mod, p_data.get("modules", {}))
            m_info = p_data["modules"][mod]
            self.assertGreater(len(m_info.get("target_kos", [])), 0)
            self.assertGreater(m_info.get("expected_locus_length", 0), 0)

    def test_multi_preset_activation_resolution(self):
        """Verify that multi-preset activation aggregates multiple modules and definitions."""
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        presets = cfg["synteny_validation"]["operon_presets"]

        # Case 1: Multiple presets as a list
        active_list = ["M00159", "M00157"]
        mod_defs = {}
        enabled = []
        for p in active_list:
            if p in presets:
                p_data = presets[p]
                mod_defs.update(p_data.get("modules", {}))
                for m in p_data.get("enabled_modules", []):
                    if m not in enabled:
                        enabled.append(m)
        self.assertEqual(enabled, ["M00159", "M00157"])
        self.assertIn("M00159", mod_defs)
        self.assertIn("M00157", mod_defs)

        # Case 2: Comma-separated string
        raw_str = "M00159, M00595"
        parsed = [p.strip() for p in raw_str.split(",") if p.strip()]
        enabled_str = []
        for p in parsed:
            if p in presets:
                for m in presets[p].get("enabled_modules", []):
                    if m not in enabled_str:
                        enabled_str.append(m)
        self.assertEqual(enabled_str, ["M00159", "M00595"])

    def test_file_index_cache_in_select_refs(self):
        """Verify that index_directory writes to and reads from file_index_cache JSON in O(1) time."""
        import tempfile
        import json
        from scripts.select_and_validate_refs import index_directory

        with tempfile.TemporaryDirectory() as tmp_dir:
            # Create a mock directory with dummy files
            mock_genome_dir = os.path.join(tmp_dir, "genomes")
            os.makedirs(mock_genome_dir, exist_ok=True)
            mock_fna = os.path.join(mock_genome_dir, "GCA_000123456.1_genomic.fna.gz")
            with open(mock_fna, "w") as f:
                f.write(">ctg1\nATCG\n")

            cache_file = os.path.join(tmp_dir, "cache", "test_file_index.json")

            # 1. First run: scans directory and generates cache
            idx1 = index_directory(mock_genome_dir, (".fna", ".fna.gz"), cache_file=cache_file)
            self.assertIn("GCA_000123456.1", idx1)
            self.assertTrue(os.path.exists(cache_file))

            # 2. Second run: reads directly from cache in O(1) time
            # Delete mock genome file to prove it loads from cache without scanning!
            os.remove(mock_fna)
            idx2 = index_directory(mock_genome_dir, (".fna", ".fna.gz"), cache_file=cache_file)
            self.assertEqual(idx1, idx2)
            self.assertIn("GCA_000123456.1", idx2)

    def test_alternative_kos_args_in_rescue(self):
        """Verify that rescue_and_classify_gene_state accepts --alternative_kos and --alternative_module."""
        from scripts.rescue_and_classify_gene_state import parse_args
        import sys

        old_argv = sys.argv
        sys.argv = [
            "rescue_and_classify_gene_state.py",
            "--target_name", "TestSample",
            "--target_locus_fa", "/dev/null",
            "--target_faa", "/dev/null",
            "--query_operon_faa", "/dev/null",
            "--alternative_module", "M00157",
            "--alternative_kos", "K02111,K02112,K02113",
            "--out_gene_tsv", "/dev/null",
            "--out_cluster_tsv", "/dev/null"
        ]
        try:
            args = parse_args()
            self.assertEqual(args.alternative_module, "M00157")
            self.assertEqual(args.alternative_kos, "K02111,K02112,K02113")
        finally:
            sys.argv = old_argv

    def test_plot_operon_fate_space(self):
        """Verify that compile_atp_synteny_matrix generates 2D Operon Fate Space PDF and PNG."""
        import tempfile
        from scripts.compile_atp_synteny_matrix import compile_matrices

        with tempfile.TemporaryDirectory() as tmp_dir:
            s4_tsv = os.path.join(tmp_dir, "s4.tsv")
            df = pd.DataFrame([{
                "target_id": "SampleA",
                "module_id": "M00159",
                "module_name": "V/A-type",
                "GTDB_taxonomy": "d__Archaea;p__Thermoproteota",
                "completeness": 98.0,
                "contamination": 0.5,
                "evolutionary_state": "TARGET_RETAINED",
                "confidence": "HIGH",
                "cluster_state": "COMPLETE",
                "structural_event_type": "NONE",
                "nearest_positive_ref": "Ref1",
                "target_contig_id": "ctg1",
                "locator_score": 0.95,
                "bilateral_context_score": 0.95,
                "K_loss": 0.0,
                "K_clean": 0.0,
                "loss_gate_status": "NOT_APPLICABLE_RETAINED",
                "OFI_fate": 1.0,
                "OFI_occupancy": 1.0,
                "fate_category": "TARGET_RETAINED"
            }])
            df.to_csv(s4_tsv, sep="\t", index=False)

            out_stacked = os.path.join(tmp_dir, "stacked.tsv")
            out_matrix = os.path.join(tmp_dir, "matrix.tsv")
            out_pdf = os.path.join(tmp_dir, "test_ofi.pdf")
            out_png = os.path.join(tmp_dir, "test_ofi.png")

            compile_matrices([s4_tsv], out_stacked, out_matrix, out_plot_pdf=out_pdf, out_plot_png=out_png)

            self.assertTrue(os.path.exists(out_pdf), "OFI PDF plot should be generated")
            self.assertTrue(os.path.exists(out_png), "OFI PNG plot should be generated")
            self.assertGreater(os.path.getsize(out_pdf), 100)
            self.assertGreater(os.path.getsize(out_png), 100)

    def test_vectorized_empirical_permutation_speed(self):
        """Verify that vectorized compute_empirical_max_contig_null returns valid p-values and dominance ratios."""
        from scripts.extract_anchor_loci import compute_empirical_max_contig_null

        target_records = {"ctg1": "NNN", "ctg2": "NNN", "ctg3": "NNN"}
        gff_coords = {
            "p1": {"contig": "ctg1"}, "p2": {"contig": "ctg1"},
            "p3": {"contig": "ctg2"}, "p4": {"contig": "ctg3"}
        }
        ctg_scores = {"ctg1": 150.0, "ctg2": 20.0, "ctg3": 10.0}
        cluster_scores = {"cl1": 80.0, "cl2": 70.0}

        p_emp, s_obs, s_null_mean, s_null_95, dominance_ratio, stat_support = compute_empirical_max_contig_null(
            target_records,
            gff_coords,
            ctg_scores,
            winning_ctg_id="ctg1",
            cluster_scores=cluster_scores,
            n_permutations=200,
            alpha=0.01,
            random_seed=42
        )

        self.assertLessEqual(p_emp, 1.0)
        self.assertGreaterEqual(p_emp, 0.0)
        self.assertEqual(s_obs, 150.0)
        self.assertGreater(dominance_ratio, 1.0)
        self.assertIn(stat_support, ["SIGNIFICANT", "NOT_SIGNIFICANT"])

    def test_flanking_ko_cache_persistence(self):
        """Verify that load_ref_flanking_annotations_from_table reads from and writes to flanking_ko_cache."""
        import tempfile
        from scripts.extract_anchor_loci import load_ref_flanking_annotations_from_table

        with tempfile.TemporaryDirectory() as tmp_dir:
            cache_file = os.path.join(tmp_dir, "flanking_ko_cache.tsv")
            with open(cache_file, "w") as f:
                f.write("protein_id\tko\tgenome_id\tko_definition\tbitscore\tevalue\tko_rank\n")
                f.write("GCA_001_p1\tK02111\tGCA_001\tatpA\t450.0\t1e-90\t1\n")

            annos = load_ref_flanking_annotations_from_table(cache_file, ["GCA_001"], {"GCA_001_p1"})
            self.assertIn("GCA_001_p1", annos)
            self.assertEqual(annos["GCA_001_p1"]["ko"], "K02111")

if __name__ == "__main__":
    unittest.main()
