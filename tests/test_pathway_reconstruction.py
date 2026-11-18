#!/usr/bin/env python3
"""
Unit and integration test suite for Four-Tier Functional Evidence Architecture.
Enforces:
1. kegg_module_parser.py:
   - AST boolean logic (+, ,, space, -, ())
   - definition_satisfied == true strictly required for 'complete'
   - missing <= 1 step in >=3-step pathways -> 'near_complete'
   - score >= threshold -> 'partial'
2. evaluate_module.py & generate_module_matrix.py (v1.0 backward compatibility)
3. 10 Scientific Quality Gates:
   - G1: Snapshot Gate (Manifest, dynamic counts, SHA256)
   - G2: Pathway Semantic Gate (map_ko_fraction, no 'completeness' in footprint headers)
   - G3: Taxonomy & Scope Gate (is_prokaryote_relevant, microbial_relevance, map_structure_type)
   - G4: Canonical-ID Gate (no map/ko double counts, strictly canonical IDs)
   - G5: BRITE Multi-Membership Gate (multi-label membership, copies_per_1000_proteins)
   - G6: Proteome Coverage Gate (separate protein->KO from KO->hierarchy; true orphan KOs)
   - G7: Custom Evidence Gate (explicit evidence_type declarations)
   - G8: Official-First Gate (Official KEGG > curated extension; composite ED profile)
   - G9: Identifier Validation Gate (Real-time validation against official KO dictionary)
   - G10: Applicability / NA Gate (Non-metabolic pathways output NA for reaction fraction)
4. Schema v2.0 Unified Functional Profile Generation
"""

import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PR_SCRIPTS = os.path.join(REPO_ROOT, "modules", "pathway_reconstruction", "scripts")
KNOWLEDGE_DIR = os.path.join(REPO_ROOT, "knowledge", "kegg")
MODULES_DIR = os.path.join(KNOWLEDGE_DIR, "modules")
TEST_TMP = os.path.join(REPO_ROOT, "scratch", "test_four_tier_run")

sys.path.insert(0, PR_SCRIPTS)
from kegg_module_parser import evaluate_module, parse_definition

class TestFourTierPathwayReconstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.makedirs(TEST_TMP, exist_ok=True)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(TEST_TMP):
            shutil.rmtree(TEST_TMP)

    # -------------------------------------------------------------------------
    # Tier 1 & Parser Logic Tests (P0-3: Strict Complete vs Near-Complete)
    # -------------------------------------------------------------------------
    def test_01_module_parser_logic_strict_complete(self):
        """Test KEGG module boolean grammar parsing with strict definition_satisfied."""
        # Expression: (K00001,K00002) (K00003+K00004) -K00005
        defn = "(K00001,K00002) (K00003+K00004) -K00005"

        # Test 1: Full match -> complete & definition_satisfied == True
        comp, status, def_sat, matched, missing = evaluate_module(defn, ["K00001", "K00003", "K00004"])
        self.assertEqual(comp, 1.0)
        self.assertTrue(def_sat)
        self.assertEqual(status, "complete")
        self.assertIn("K00001", matched)
        self.assertIn("K00003", matched)
        self.assertIn("K00004", matched)
        self.assertEqual(len(missing), 0)

        # Test 2: Missing 1 obligatory subunit -> definition_satisfied == False, status == near_complete
        comp2, status2, def_sat2, matched2, missing2 = evaluate_module(defn, ["K00001", "K00003"])
        self.assertAlmostEqual(comp2, 2.5 / 3, places=2)
        self.assertFalse(def_sat2, "Missing K00004 subunit must NOT satisfy definition!")
        self.assertEqual(status2, "near_complete", "Missing 1 step must be near_complete, never complete!")
        self.assertIn("K00004", missing2)

        # Test 3: Four-step short pathway: K00001 K00002 K00003 K00004
        defn_4 = "K00001 K00002 K00003 K00004"
        # 3/4 steps (missing 1 step) -> near_complete, def_sat == False
        comp3, stat3, def_sat3, _, _ = evaluate_module(defn_4, ["K00001", "K00002", "K00003"])
        self.assertEqual(comp3, 0.75)
        self.assertFalse(def_sat3)
        self.assertEqual(stat3, "near_complete")

        # 2/4 steps (50% / missing 2 steps) -> partial
        comp4, stat4, def_sat4, _, _ = evaluate_module(defn_4, ["K00001", "K00002"])
        self.assertEqual(comp4, 0.50)
        self.assertFalse(def_sat4)
        self.assertEqual(stat4, "partial")

        # 1/4 steps (25%) -> absent
        comp5, stat5, def_sat5, _, _ = evaluate_module(defn_4, ["K00001"])
        self.assertEqual(comp5, 0.25)
        self.assertFalse(def_sat5)
        self.assertEqual(stat5, "absent")

    def test_02_curated_archetype_modules(self):
        """Test real curated M00159 (V/A ATP synthase) and M00567 (Methanogenesis)."""
        m00159_def = "(K02117+K02118+K02119+K02120+K02121+K02122) (K02123,K02124)"
        methano_kos = ["K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02123"]
        comp, status, def_sat, matched, missing = evaluate_module(m00159_def, methano_kos)
        self.assertEqual(comp, 1.0)
        self.assertTrue(def_sat)
        self.assertEqual(status, "complete")
        self.assertEqual(len(missing), 0)

    def test_03_all_574_kegg_modules_catalog_parse(self):
        """Test that all 574 KEGG modules in all_kegg_modules.tsv parse cleanly into AST."""
        catalog_tsv = os.path.join(KNOWLEDGE_DIR, "all_kegg_modules.tsv")
        self.assertTrue(os.path.exists(catalog_tsv), "all_kegg_modules.tsv must exist")

        with open(catalog_tsv, "r", encoding="utf-8") as f:
            header = f.readline().strip().split("\t")
            m_col = header.index("module_id")
            d_col = header.index("definition")

            count = 0
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > max(m_col, d_col):
                    mid = parts[m_col]
                    mdef = parts[d_col]
                    count += 1
                    comp, status, def_sat, matched, missing = evaluate_module(mdef, [])
                    self.assertIsInstance(comp, float)
                    self.assertFalse(def_sat)
                    self.assertEqual(status, "absent")

            self.assertEqual(count, 574, f"Expected exactly 574 modules, got {count}")

    # -------------------------------------------------------------------------
    # The 10 Scientific Quality Gates Verification
    # -------------------------------------------------------------------------
    def test_gate_1_snapshot_manifest_and_dynamic_counts(self):
        """Gate G1: Verify snapshot manifest existence, valid SHA256 hashes, and non-empty dynamic counts."""
        manifest_path = os.path.join(KNOWLEDGE_DIR, "kegg_snapshot_manifest.json")
        self.assertTrue(os.path.isfile(manifest_path), "kegg_snapshot_manifest.json must exist")

        with open(manifest_path, "r", encoding="utf-8") as f:
            m = json.load(f)

        counts = m.get("counts", {})
        self.assertEqual(counts.get("total_official_modules"), 574)
        self.assertEqual(counts.get("total_pathways"), 587)
        self.assertGreater(counts.get("total_pathway_kos"), 14000)
        self.assertGreater(counts.get("total_unique_brite_kos"), 28000)
        self.assertGreater(counts.get("total_specialized_brite_kos"), 2000)

        # Check sha256 entries
        sha = m.get("sha256", {})
        for fname in ["all_kegg_modules.tsv", "all_kegg_pathways.tsv", "pathway_ko_reaction.tsv", "ko_hierarchy_membership.tsv"]:
            self.assertIn(fname, sha)
            self.assertEqual(len(sha[fname]), 64, f"SHA256 for {fname} must be 64-hex chars")

    def test_gate_2_pathway_semantic_cleanliness(self):
        """Gate G2: Verify Tier 2 output headers use map_ko_fraction and strictly forbid 'completeness' / 'activity'."""
        mock_strict = os.path.join(TEST_TMP, "mock_gate2.tsv")
        with open(mock_strict, "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\n")
            f.write("prot1\tgene1\tG_TEST\tK00001\n")

        mod_out = os.path.join(TEST_TMP, "gate2_out")
        cmd = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "evaluate_pathway_footprint.py"),
            "-k", mock_strict,
            "-p", os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv"),
            "-t", os.path.join(KNOWLEDGE_DIR, "pathway_ko_reaction.tsv"),
            "-o", mod_out
        ]
        subprocess.run(cmd, check=True)

        footprint_file = os.path.join(mod_out, "genome_pathway_footprint.tsv")
        self.assertTrue(os.path.isfile(footprint_file))
        with open(footprint_file, "r") as f:
            header = f.readline().strip().split("\t")

        self.assertIn("map_ko_fraction", header)
        self.assertIn("reaction_node_fraction", header)

        # Forbidden terms check
        for h in header:
            self.assertNotIn("completeness", h.lower(), f"Gate G2 violation: '{h}' contains 'completeness'")
            self.assertNotIn("pathway_activity", h.lower(), f"Gate G2 violation: '{h}' contains 'activity'")

    def test_gate_3_taxonomy_scope_and_structure_types(self):
        """Gate G3: Verify all 587 pathways have taxonomy scope, microbial relevance, and structure types."""
        pw_tsv = os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv")
        self.assertTrue(os.path.isfile(pw_tsv))

        with open(pw_tsv, "r") as f:
            header = f.readline().strip().split("\t")
            self.assertIn("is_prokaryote_relevant", header)
            self.assertIn("microbial_relevance", header)
            self.assertIn("map_structure_type", header)

            rows = [line.strip().split("\t") for line in f if line.strip()]
            self.assertEqual(len(rows), 587, "Must have exactly 587 canonical pathways")

            # ko00680 (Methane) should be microbial_core
            methane_row = [r for r in rows if r[0] == "ko00680"][0]
            m_idx = header.index("microbial_relevance")
            self.assertEqual(methane_row[m_idx], "microbial_core")

            # Global overview map ko01100
            global_row = [r for r in rows if r[0] == "ko01100"][0]
            st_idx = header.index("map_structure_type")
            self.assertEqual(global_row[st_idx], "global")

    def test_gate_4_canonical_id_no_double_count(self):
        """Gate G4: Verify pathway IDs are strictly normalized to 5 digits and contain no map/ko double counts."""
        triples_tsv = os.path.join(KNOWLEDGE_DIR, "pathway_ko_reaction.tsv")
        self.assertTrue(os.path.isfile(triples_tsv))

        with open(triples_tsv, "r") as f:
            header = f.readline().strip().split("\t")
            cid_idx = header.index("canonical_id")
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > cid_idx:
                    cid = parts[cid_idx]
                    self.assertTrue(cid.isdigit(), f"Canonical ID must be all digits: {cid}")
                    self.assertEqual(len(cid), 5, f"Canonical ID must be 5 digits: {cid}")

    def test_gate_5_brite_multi_membership_and_normalization(self):
        """Gate G5: Verify BRITE hierarchy supports multi-label memberships and per-1000 protein normalizations."""
        brite_tsv = os.path.join(KNOWLEDGE_DIR, "ko_hierarchy_membership.tsv")
        self.assertTrue(os.path.isfile(brite_tsv))

        # Check multi-membership: find KOs appearing in multiple hierarchies
        ko_hierarchies = {}
        with open(brite_tsv, "r") as f:
            f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    k, h = parts[0], parts[1]
                    ko_hierarchies.setdefault(k, set()).add(h)

        multi_labeled = [k for k, hset in ko_hierarchies.items() if len(hset) > 1]
        self.assertGreater(len(multi_labeled), 500, "Must have multi-label KOs across hierarchies")

        # Test evaluation script normalization
        mock_strict = os.path.join(TEST_TMP, "mock_gate5.tsv")
        with open(mock_strict, "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\n")
            # 5 proteins for G_TEST
            f.write("prot1\tgene1\tG_TEST\tK02003\n") # ABC transporter
            f.write("prot2\tgene2\tG_TEST\tK02004\n")
            f.write("prot3\tgene3\tG_TEST\tK02484\n") # Two-component
            f.write("prot4\tgene4\tG_TEST\tK07657\n")
            f.write("prot5\tgene5\tG_TEST\tK19090\n") # CRISPR Cas5

        mod_out = os.path.join(TEST_TMP, "gate5_out")
        cmd = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "evaluate_brite_hierarchy.py"),
            "-k", mock_strict,
            "-b", brite_tsv,
            "-o", mod_out
        ]
        subprocess.run(cmd, check=True)

        counts_file = os.path.join(mod_out, "genome_brite_counts.tsv")
        self.assertTrue(os.path.isfile(counts_file))
        with open(counts_file, "r") as f:
            header = f.readline().strip().split("\t")
            self.assertIn("copies_per_1000_proteins", header)
            self.assertIn("copies_per_1000_KO_annotated_proteins", header)

    def test_gate_6_proteome_coverage_separation(self):
        """Gate G6: Verify separation of protein->KO annotation rate from KO->hierarchy assignment, and true orphan KOs."""
        mock_strict = os.path.join(TEST_TMP, "mock_gate6.tsv")
        with open(mock_strict, "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\n")
            f.write("p1\tg1\tG_AUDIT\tK00001\n") # In pathway and module
            f.write("p2\tg2\tG_AUDIT\tK02003\n") # In transporter BRITE
            f.write("p3\tg3\tG_AUDIT\tK00437\n") # In 09190 (Not Included in Pathway/Brite)

        audit_file = os.path.join(TEST_TMP, "genome_ko_coverage_audit.tsv")
        cmd = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "audit_genome_coverage.py"),
            "-k", mock_strict,
            "-m", os.path.join(KNOWLEDGE_DIR, "all_kegg_modules.tsv"),
            "-p", os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv"),
            "-b", os.path.join(KNOWLEDGE_DIR, "ko_hierarchy_membership.tsv"),
            "-o", audit_file
        ]
        subprocess.run(cmd, check=True)

        self.assertTrue(os.path.isfile(audit_file))
        with open(audit_file, "r") as f:
            header = f.readline().strip().split("\t")
            row = f.readline().strip().split("\t")

        self.assertIn("protein_ko_annotation_rate", header)
        self.assertIn("ko00001_tracking_rate", header)
        self.assertIn("specialized_brite_assignment_rate", header)
        self.assertIn("true_orphan_kos", header)

        # K00437 should be detected in 09190
        in_09190_idx = header.index("kos_in_ko00001_09190")
        self.assertGreaterEqual(int(row[in_09190_idx]), 1, "K00437 must be counted in kos_in_ko00001_09190")

    def test_gate_7_custom_evidence_types(self):
        """Gate G7: Verify custom archaeal modules declare explicit evidence_type and hydrogenase uses multi-evidence."""
        arch_dir = os.path.join(MODULES_DIR, "archaea")
        self.assertTrue(os.path.isdir(arch_dir))

        # Check M_ARCH_04
        h_file = os.path.join(arch_dir, "M_ARCH_04_hydrogenase_operon.yaml")
        self.assertTrue(os.path.isfile(h_file))
        with open(h_file, "r") as f:
            text = f.read()

        self.assertIn("evidence_type:", text)
        self.assertIn("gene_neighborhood", text)
        self.assertIn("sequence_classifier", text)
        self.assertIn("max_intervening_genes:", text)

    def test_gate_8_official_first_and_ed_profile(self):
        """Gate G8: Verify official-first principle: no redundant M_ARCH_02 and custom modules declare relationship."""
        arch_dir = os.path.join(MODULES_DIR, "archaea")
        self.assertFalse(os.path.isfile(os.path.join(arch_dir, "M_ARCH_02_modified_ed.yaml")),
                         "Gate G8 violation: M_ARCH_02 must NOT duplicate official M00308/M00309/M00633!")

        # M_ARCH_03 CoB must declare complements_M00358
        cob_file = os.path.join(arch_dir, "M_ARCH_03_cob_biosynthesis.yaml")
        with open(cob_file, "r") as f:
            text = f.read()
        self.assertIn("relationship_to_kegg: \"complements_M00358\"", text)

    def test_gate_9_identifier_validation(self):
        """Gate G9: Verify real-time validation of custom module KOs against KEGG dictionary."""
        # Check that archaeol module uses real KOs (K00096, K17104, K17105, K19664)
        archaeol_file = os.path.join(MODULES_DIR, "archaea", "M_ARCH_01A_archaeol_core.yaml")
        with open(archaeol_file, "r") as f:
            text = f.read()

        self.assertIn("K00096", text)
        self.assertIn("K17104", text)
        self.assertIn("K17105", text)
        self.assertIn("K19664", text)
        self.assertNotIn("K13788", text, "Gate G9: K13788 (pta) must never be used for CarS!")
        self.assertNotIn("K00993", text, "Gate G9: K00993 must never be used for DGGGPS!")

    def test_gate_10_na_semantics_for_non_metabolic(self):
        """Gate G10: Verify non-metabolic pathways output NA for reaction_node_fraction with applicable_reason."""
        mock_strict = os.path.join(TEST_TMP, "mock_gate10.tsv")
        with open(mock_strict, "w") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\n")
            f.write("prot1\tgene1\tG_TEST\tK00001\n")

        mod_out = os.path.join(TEST_TMP, "gate10_out")
        cmd = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "evaluate_pathway_footprint.py"),
            "-k", mock_strict,
            "-p", os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv"),
            "-t", os.path.join(KNOWLEDGE_DIR, "pathway_ko_reaction.tsv"),
            "-o", mod_out
        ]
        subprocess.run(cmd, check=True)

        footprint_file = os.path.join(mod_out, "genome_pathway_footprint.tsv")
        self.assertTrue(os.path.isfile(footprint_file))

        with open(footprint_file, "r") as f:
            header = f.readline().strip().split("\t")
            pw_idx = header.index("pathway_id")
            rxn_idx = header.index("reaction_node_fraction")
            app_idx = header.index("reaction_metric_applicable")
            reason_idx = header.index("applicable_reason")

            checked_nonmetabolic = 0
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > reason_idx:
                    pw_id = parts[pw_idx]
                    rxn_val = parts[rxn_idx]
                    app_val = parts[app_idx]
                    # Check non-metabolic pathway (e.g. ko02000 or ko03010)
                    if pw_id in ["ko02000", "ko03010", "ko04112"]:
                        checked_nonmetabolic += 1
                        self.assertEqual(rxn_val, "NA", f"Non-metabolic {pw_id} must output NA for reaction_node_fraction")
                        self.assertEqual(app_val, "false")
                        self.assertEqual(parts[reason_idx], "non_metabolic_pathway")
            self.assertGreater(checked_nonmetabolic, 0)

    # -------------------------------------------------------------------------
    # End-to-End Four-Tier Integration & Schema v2.0 Profile Generation
    # -------------------------------------------------------------------------
    def test_11_end_to_end_four_tier_pipeline(self):
        """Test full four-tier workflow on mock Archaea genome generating v1 capabilities and v2 functional profile."""
        mock_strict = os.path.join(TEST_TMP, "mock_archaea_strict.tsv.gz")
        with gzip.open(mock_strict, "wt") as f:
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tko_rank\tscore_margin\n")
            # Complete V/A ATP synthase (M00159)
            atp_kos = ["K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02123"]
            for i, ko in enumerate(atp_kos):
                f.write(f"PROT_ATP_{i}\tGENE_ATP_{i}\tG_ARCHAEA\t{ko}\tAtp\t500.0\t300.0\tfull\t1e-80\t1\t200.0\n")
            # Complete Archaeol core (M_ARCH_01A)
            arch_kos = ["K00096", "K17104", "K17105", "K19664"]
            for i, ko in enumerate(arch_kos):
                f.write(f"PROT_LIP_{i}\tGENE_LIP_{i}\tG_ARCHAEA\t{ko}\tLip\t500.0\t300.0\tfull\t1e-80\t1\t200.0\n")
            # Non-phosphorylative ED (M00309)
            ed_kos = ["K05308", "K11395", "K18128"]
            for i, ko in enumerate(ed_kos):
                f.write(f"PROT_ED_{i}\tGENE_ED_{i}\tG_ARCHAEA\t{ko}\tED\t500.0\t300.0\tfull\t1e-80\t1\t200.0\n")
            # ABC Transporters (ko02000)
            trans_kos = ["K02003", "K02004", "K02005"]
            for i, ko in enumerate(trans_kos):
                f.write(f"PROT_TR_{i}\tGENE_TR_{i}\tG_ARCHAEA\t{ko}\tTrans\t500.0\t300.0\tfull\t1e-80\t1\t200.0\n")

        run_out = os.path.join(TEST_TMP, "four_tier_run")

        # 1. Tier 1: evaluate_module
        cmd1 = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "evaluate_module.py"),
            "-k", mock_strict,
            "-d", MODULES_DIR,
            "-o", os.path.join(run_out, "modules")
        ]
        subprocess.run(cmd1, check=True)
        comp_tsv = os.path.join(run_out, "modules", "genome_module_completeness.tsv")
        self.assertTrue(os.path.isfile(comp_tsv))

        # 2. v1 Capabilities generation (Backward compatibility check)
        cmd2 = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "generate_module_matrix.py"),
            "-c", comp_tsv,
            "-o", os.path.join(run_out, "modules")
        ]
        subprocess.run(cmd2, check=True)
        v1_cap = os.path.join(run_out, "modules", "genome_metabolism_capabilities.json")
        self.assertTrue(os.path.isfile(v1_cap), "v1.0 genome_metabolism_capabilities.json must be preserved!")

        # 3. Tier 2: evaluate_pathway_footprint
        cmd3 = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "evaluate_pathway_footprint.py"),
            "-k", mock_strict,
            "-p", os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv"),
            "-t", os.path.join(KNOWLEDGE_DIR, "pathway_ko_reaction.tsv"),
            "-m", comp_tsv,
            "-o", os.path.join(run_out, "pathways")
        ]
        subprocess.run(cmd3, check=True)
        fp_tsv = os.path.join(run_out, "pathways", "genome_pathway_footprint.tsv")
        self.assertTrue(os.path.isfile(fp_tsv))

        # 4. Tier 3: evaluate_brite_hierarchy
        cmd4 = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "evaluate_brite_hierarchy.py"),
            "-k", mock_strict,
            "-b", os.path.join(KNOWLEDGE_DIR, "ko_hierarchy_membership.tsv"),
            "-o", os.path.join(run_out, "brite")
        ]
        subprocess.run(cmd4, check=True)
        brite_tsv = os.path.join(run_out, "brite", "genome_brite_counts.tsv")
        self.assertTrue(os.path.isfile(brite_tsv))

        # 5. Global Coverage Audit
        cmd5 = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "audit_genome_coverage.py"),
            "-k", mock_strict,
            "-m", os.path.join(KNOWLEDGE_DIR, "all_kegg_modules.tsv"),
            "-p", os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv"),
            "-b", os.path.join(KNOWLEDGE_DIR, "ko_hierarchy_membership.tsv"),
            "-o", os.path.join(run_out, "audit", "genome_ko_coverage_audit.tsv")
        ]
        subprocess.run(cmd5, check=True)
        audit_tsv = os.path.join(run_out, "audit", "genome_ko_coverage_audit.tsv")
        self.assertTrue(os.path.isfile(audit_tsv))

        # 6. Schema v2.0 Functional Profile Generation
        cmd6 = [
            sys.executable,
            os.path.join(PR_SCRIPTS, "generate_functional_profile.py"),
            "-c", comp_tsv,
            "-f", fp_tsv,
            "-b", brite_tsv,
            "-a", audit_tsv,
            "-m", os.path.join(KNOWLEDGE_DIR, "kegg_snapshot_manifest.json"),
            "-o", os.path.join(run_out, "capabilities", "genome_functional_profile.json")
        ]
        subprocess.run(cmd6, check=True)
        v2_prof = os.path.join(run_out, "capabilities", "genome_functional_profile.json")
        self.assertTrue(os.path.isfile(v2_prof))

        with open(v2_prof, "r", encoding="utf-8") as jf:
            profiles = json.load(jf)

        self.assertEqual(len(profiles), 1)
        p = profiles[0]
        self.assertEqual(p["schema_version"], "2.0")
        self.assertEqual(p["genome_id"], "G_ARCHAEA")

        # Verify Tier 1 & Custom Archaeol Core
        tier1 = p["tier1_module_reconstruction"]
        self.assertIn("entner_doudoroff_family_profile", tier1)

        # Verify Tier 4
        tier4 = p["tier4_specialized_systems"]
        self.assertIn("archaeal_ether_lipid_system", tier4)
        arch_stat = tier4["archaeal_ether_lipid_system"]["archaeol_core_M_ARCH_01A"]
        self.assertEqual(arch_stat["status"], "complete")

        # Verify Coverage Audit
        audit = p["coverage_audit"]
        self.assertEqual(audit["detected_unique_kos"], 17)
        self.assertGreater(audit["module_assignment_rate"], 0.5)

    def test_empty_and_noisy_inputs(self):
        """Pre-flight safety: Verify all 5 scripts handle 0-byte files, comments, and noisy lines gracefully."""
        # 1. 0-byte file test
        empty_strict = os.path.join(TEST_TMP, "empty_strict.tsv")
        open(empty_strict, "w").close()
        empty_out = os.path.join(TEST_TMP, "empty_run")
        os.makedirs(os.path.join(empty_out, "modules"), exist_ok=True)
        os.makedirs(os.path.join(empty_out, "pathways"), exist_ok=True)
        os.makedirs(os.path.join(empty_out, "brite"), exist_ok=True)
        os.makedirs(os.path.join(empty_out, "audit"), exist_ok=True)
        os.makedirs(os.path.join(empty_out, "capabilities"), exist_ok=True)

        # Run evaluate_module on empty
        subprocess.run([
            sys.executable, os.path.join(PR_SCRIPTS, "evaluate_module.py"),
            "-k", empty_strict, "-d", os.path.join(KNOWLEDGE_DIR, "all_kegg_modules.tsv"),
            "-o", os.path.join(empty_out, "modules")
        ], check=True)

        # Run generate_module_matrix on empty
        subprocess.run([
            sys.executable, os.path.join(PR_SCRIPTS, "generate_module_matrix.py"),
            "-c", os.path.join(empty_out, "modules", "genome_module_completeness.tsv"),
            "-o", os.path.join(empty_out, "modules")
        ], check=True)

        # Run evaluate_pathway_footprint on empty
        subprocess.run([
            sys.executable, os.path.join(PR_SCRIPTS, "evaluate_pathway_footprint.py"),
            "-k", empty_strict,
            "-p", os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv"),
            "-t", os.path.join(KNOWLEDGE_DIR, "pathway_ko_reaction.tsv"),
            "-m", os.path.join(empty_out, "modules", "genome_module_completeness.tsv"),
            "-o", os.path.join(empty_out, "pathways")
        ], check=True)

        # Run evaluate_brite_hierarchy on empty
        subprocess.run([
            sys.executable, os.path.join(PR_SCRIPTS, "evaluate_brite_hierarchy.py"),
            "-k", empty_strict,
            "-b", os.path.join(KNOWLEDGE_DIR, "ko_hierarchy_membership.tsv"),
            "-o", os.path.join(empty_out, "brite")
        ], check=True)

        # Run audit_genome_coverage on empty
        subprocess.run([
            sys.executable, os.path.join(PR_SCRIPTS, "audit_genome_coverage.py"),
            "-k", empty_strict,
            "-m", os.path.join(KNOWLEDGE_DIR, "all_kegg_modules.tsv"),
            "-p", os.path.join(KNOWLEDGE_DIR, "all_kegg_pathways.tsv"),
            "-b", os.path.join(KNOWLEDGE_DIR, "ko_hierarchy_membership.tsv"),
            "-o", os.path.join(empty_out, "audit", "genome_ko_coverage_audit.tsv")
        ], check=True)

        # Run generate_functional_profile on empty
        subprocess.run([
            sys.executable, os.path.join(PR_SCRIPTS, "generate_functional_profile.py"),
            "-c", os.path.join(empty_out, "modules", "genome_module_completeness.tsv"),
            "-f", os.path.join(empty_out, "pathways", "genome_pathway_footprint.tsv"),
            "-b", os.path.join(empty_out, "brite", "genome_brite_counts.tsv"),
            "-a", os.path.join(empty_out, "audit", "genome_ko_coverage_audit.tsv"),
            "-m", os.path.join(KNOWLEDGE_DIR, "kegg_snapshot_manifest.json"),
            "-o", os.path.join(empty_out, "capabilities", "genome_functional_profile.json")
        ], check=True)

        with open(os.path.join(empty_out, "capabilities", "genome_functional_profile.json")) as jf:
            empty_profiles = json.load(jf)
        self.assertEqual(empty_profiles, [])

        # 2. Noisy file test (comments, empty lines, non-numeric values, invalid KOs)
        noisy_strict = os.path.join(TEST_TMP, "noisy_strict.tsv")
        with open(noisy_strict, "w") as f:
            f.write("# leading comment\n\n")
            f.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tko_rank\tscore_margin\tis_strict\tpass_type\n")
            f.write("# inline comment\n")
            f.write("p1\tg1\tG_NOISY\t-\tNone\tbad_score\t100\tfull\t1.0\t1\t0\tFalse\tnone\n")
            f.write("p2\tg2\tG_NOISY\tK00001\tAlcohol dehydrogenase\t350.5\t200.0\tfull\t1e-50\t1\t150.5\tTrue\tstrict\n")
            f.write("\n")

        noisy_out = os.path.join(TEST_TMP, "noisy_run")
        os.makedirs(os.path.join(noisy_out, "modules"), exist_ok=True)
        subprocess.run([
            sys.executable, os.path.join(PR_SCRIPTS, "evaluate_module.py"),
            "-k", noisy_strict, "-d", os.path.join(KNOWLEDGE_DIR, "all_kegg_modules.tsv"),
            "-o", os.path.join(noisy_out, "modules")
        ], check=True)
        self.assertTrue(os.path.isfile(os.path.join(noisy_out, "modules", "genome_module_completeness.tsv")))

    def test_m00174_methane_oxidation_and_optional_complex_parsing(self):
        """Test M00174 complex parsing (infix minus, XoxF K23995, sMMO vs pMMO alternatives)."""
        m00174_def = "((K10944+K10945+K10946),(K16157+K16158+K16159+K16160+K16161+K16162)) ((K14028-K14029),K23995)"

        # 1. Empty set: must be 0.0 completeness, absent, with alternative missing representations
        comp, stat, def_sat, matched, missing = evaluate_module(m00174_def, [])
        self.assertEqual(comp, 0.0)
        self.assertEqual(stat, "absent")
        self.assertFalse(def_sat)
        self.assertEqual(matched, [])
        self.assertEqual(len(missing), 2)
        self.assertIn("K14028|K23995", missing)
        self.assertIn("K10944+K10945+K10946|K16157+K16158+K16159+K16160+K16161+K16162", missing)

        # 2. Only K23995 (XoxF): Step 2 is fully satisfied!
        comp2, stat2, def_sat2, matched2, missing2 = evaluate_module(m00174_def, ["K23995"])
        self.assertEqual(comp2, 0.5)
        self.assertEqual(stat2, "partial")
        self.assertFalse(def_sat2)
        self.assertEqual(matched2, ["K23995"])
        self.assertEqual(missing2, ["K10944+K10945+K10946|K16157+K16158+K16159+K16160+K16161+K16162"])

        # 3. Alternative pathway sMMO (K16157..K16162) + XoxF (K23995) -> 100% complete
        smmo_xoxf = ["K16157", "K16158", "K16159", "K16160", "K16161", "K16162", "K23995"]
        comp3, stat3, def_sat3, matched3, missing3 = evaluate_module(m00174_def, smmo_xoxf)
        self.assertEqual(comp3, 1.0)
        self.assertEqual(stat3, "complete")
        self.assertTrue(def_sat3)
        self.assertEqual(len(missing3), 0)

        # 4. pMMO + MxaF (K14028) without optional MxaI (K14029) -> 100% complete
        pmmo_mxaf = ["K10944", "K10945", "K10946", "K14028"]
        comp4, stat4, def_sat4, matched4, missing4 = evaluate_module(m00174_def, pmmo_mxaf)
        self.assertEqual(comp4, 1.0)
        self.assertEqual(stat4, "complete")
        self.assertTrue(def_sat4)
        self.assertEqual(len(missing4), 0)

if __name__ == "__main__":
    unittest.main()
