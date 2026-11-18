#!/usr/bin/env python3
"""
generate_functional_profile.py
Builds Schema v2.0 genome_functional_profile.json for the Four-Tier Evidence Architecture.
Preserves existing genome_metabolism_capabilities.json (v1.0) 100% untouched.

Includes:
- Tier 1: Deterministic Module Reconstruction (complete/near_complete/partial)
  and composite ED_pathway_family_profile (M00008, M00308, M00309, M00633).
- Tier 2: Pathway Footprint Highlights (map_ko_fraction & reaction_node_fraction).
- Tier 3: BRITE Molecular-System Profiling (Transporters, Two-component, Defense/CRISPR).
- Tier 4: Specialized Systems (Archaeol core, GDGT/Tes potential, CoB, Hydrogenase operons).
"""

import argparse
import gzip
import json
import os
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def load_manifest(manifest_path):
    if manifest_path and os.path.isfile(manifest_path):
        with open(manifest_path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {"snapshot_date": "unknown", "kegg_release": "KEGG Release 115.0+"}

def main():
    parser = argparse.ArgumentParser(description="Generate Schema v2.0 genome_functional_profile.json.")
    parser.add_argument("-c", "--module_completeness", required=True, help="Path to genome_module_completeness.tsv")
    parser.add_argument("-f", "--pathway_footprint", required=True, help="Path to genome_pathway_footprint.tsv")
    parser.add_argument("-b", "--brite_counts", required=True, help="Path to genome_brite_counts.tsv")
    parser.add_argument("-a", "--coverage_audit", required=True, help="Path to genome_ko_coverage_audit.tsv")
    parser.add_argument("-m", "--manifest", default="knowledge/kegg/kegg_snapshot_manifest.json", help="Path to snapshot manifest")
    parser.add_argument("-o", "--out_json", required=True, help="Output path for genome_functional_profile.json")

    args = parser.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)

    manifest_data = load_manifest(args.manifest)

    # 1. Load Coverage Audits
    # genome_id -> audit dict
    audits = {}
    if os.path.isfile(args.coverage_audit):
        with open_file(args.coverage_audit, "rt") as f:
            first_line = f.readline()
            if first_line:
                header = first_line.strip().split("\t")
                if "genome_id" in header:
                    col_map = {c: i for i, c in enumerate(header)}
                    for line in f:
                        p = line.strip().split("\t")
                        if len(p) >= len(header):
                            gid = p[col_map["genome_id"]]
                            audits[gid] = {
                                "total_predicted_proteins": int(p[col_map["total_predicted_proteins"]]),
                                "proteins_with_ko": int(p[col_map["proteins_with_ko"]]),
                                "protein_ko_annotation_rate": float(p[col_map["protein_ko_annotation_rate"]]),
                                "detected_unique_kos": int(p[col_map["detected_unique_kos"]]),
                                "kos_in_modules": int(p[col_map["kos_in_modules"]]),
                                "module_assignment_rate": float(p[col_map["module_assignment_rate"]]),
                                "kos_in_pathways": int(p[col_map["kos_in_pathways"]]),
                                "pathway_assignment_rate": float(p[col_map["pathway_assignment_rate"]]),
                                "kos_in_specialized_brite": int(p[col_map["kos_in_specialized_brite"]]),
                                "specialized_brite_assignment_rate": float(p[col_map["specialized_brite_assignment_rate"]]),
                                "kos_in_ko00001_09190": int(p[col_map["kos_in_ko00001_09190"]]),
                                "ko00001_tracking_rate": float(p[col_map["ko00001_tracking_rate"]]),
                                "true_orphan_kos": int(p[col_map["true_orphan_kos"]]),
                                "true_orphan_rate": float(p[col_map["true_orphan_rate"]])
                            }

    # 2. Load Module Results
    # genome_id -> {module_id -> (completeness, status)}
    modules_data = {}
    if os.path.isfile(args.module_completeness):
        with open_file(args.module_completeness, "rt") as f:
            first_line = f.readline()
            if first_line:
                header = first_line.strip().split("\t")
                if "genome_id" in header and "module_id" in header and "completeness" in header and "status" in header:
                    g_idx = header.index("genome_id")
                    m_idx = header.index("module_id")
                    c_idx = header.index("completeness")
                    s_idx = header.index("status")
                    for line in f:
                        p = line.strip().split("\t")
                        if len(p) > max(g_idx, m_idx, c_idx, s_idx):
                            gid = p[g_idx]
                            mid = p[m_idx]
                            comp = float(p[c_idx])
                            stat = p[s_idx]
                            modules_data.setdefault(gid, {})[mid] = {"completeness": comp, "status": stat}

    # 3. Load Pathway Footprints
    # genome_id -> {pathway_id -> {map_ko_fraction, reaction_node_fraction, scope}}
    pathways_data = {}
    if os.path.isfile(args.pathway_footprint):
        with open_file(args.pathway_footprint, "rt") as f:
            first_line = f.readline()
            if first_line:
                header = first_line.strip().split("\t")
                if "genome_id" in header and "pathway_id" in header:
                    col_map = {c: i for i, c in enumerate(header)}
                    for line in f:
                        p = line.strip().split("\t")
                        if len(p) >= len(header):
                            gid = p[col_map["genome_id"]]
                            pw_id = p[col_map["pathway_id"]]
                            pathways_data.setdefault(gid, {})[pw_id] = {
                                "name": p[col_map["pathway_name"]],
                                "map_ko_fraction": float(p[col_map["map_ko_fraction"]]),
                                "reaction_node_fraction": p[col_map["reaction_node_fraction"]],
                                "microbial_relevance": p[col_map["microbial_relevance"]]
                            }

    # 4. Load BRITE System Profiles
    # genome_id -> {system_name -> {copies_per_1000_proteins, unique_KO_count}}
    brite_data = {}
    if os.path.isfile(args.brite_counts):
        with open_file(args.brite_counts, "rt") as f:
            first_line = f.readline()
            if first_line:
                header = first_line.strip().split("\t")
                if "genome_id" in header and "hierarchy_id" in header and "system_name" in header:
                    col_map = {c: i for i, c in enumerate(header)}
                    for line in f:
                        p = line.strip().split("\t")
                        if len(p) >= len(header):
                            gid = p[col_map["genome_id"]]
                            hier_id = p[col_map["hierarchy_id"]]
                            sys_name = p[col_map["system_name"]]
                            if hier_id in ["ko02000", "ko02022", "ko02048"]:
                                sys_key = f"{hier_id}___{sys_name}"
                                brite_data.setdefault(gid, {})[sys_key] = {
                                    "hierarchy_id": hier_id,
                                    "system_name": sys_name,
                                    "unique_KO_count": int(p[col_map["unique_KO_count"]]),
                                    "copies_per_1000_proteins": float(p[col_map["copies_per_1000_proteins"]]),
                                    "presence": int(p[col_map["presence"]])
                                }

    # 5. Assemble Profiles
    all_genomes = sorted(set(audits.keys()) | set(modules_data.keys()) | set(pathways_data.keys()) | set(brite_data.keys()))
    profiles = []

    for gid in all_genomes:
        g_mods = modules_data.get(gid, {})
        g_pws = pathways_data.get(gid, {})
        g_brite = brite_data.get(gid, {})
        g_audit = audits.get(gid, {})

        # Composite ED Pathway Family Profile (Official-First)
        ed_profile = {
            "classical_ED_M00008": g_mods.get("M00008", {"status": "absent", "completeness": 0.0}),
            "semi_phosphorylative_ED_euryarchaeota_M00308": g_mods.get("M00308", {"status": "absent", "completeness": 0.0}),
            "non_phosphorylative_ED_archaea_M00309": g_mods.get("M00309", {"status": "absent", "completeness": 0.0}),
            "crenarchaeal_branched_semi_ED_M00633": g_mods.get("M00633", {"status": "absent", "completeness": 0.0}),
            "dominant_variant": "none"
        }
        # Determine dominant variant
        best_stat = "absent"
        best_variant = "none"
        for vname, mid in [("non_phosphorylative", "M00309"), ("semi_phosphorylative", "M00308"), ("crenarchaeal_branched", "M00633"), ("classical", "M00008")]:
            stat = g_mods.get(mid, {}).get("status", "absent")
            if stat == "complete":
                best_variant = vname
                break
            elif stat in ["near_complete", "partial"] and best_variant == "none":
                best_variant = vname
        ed_profile["dominant_variant"] = best_variant

        # Tier 4 Specialized Systems
        archaeol_stat = g_mods.get("M_ARCH_01A", {"status": "absent", "completeness": 0.0})
        cob_stat = g_mods.get("M_ARCH_03", {"status": "absent", "completeness": 0.0})
        hyd_stat = g_mods.get("M_ARCH_04", {"status": "absent", "completeness": 0.0})

        specialized_systems = {
            "archaeal_ether_lipid_system": {
                "archaeol_core_M_ARCH_01A": archaeol_stat,
                "gdgt_tetraether_potential": "evaluated_via_tes" if archaeol_stat["status"] in ["complete", "near_complete"] else "precursor_absent"
            },
            "methanogenesis_cofactor_systems": {
                "coenzyme_B_M_ARCH_03": cob_stat,
                "coenzyme_M_M00358": g_mods.get("M00358", {"status": "absent", "completeness": 0.0}),
                "coenzyme_F420_M00378": g_mods.get("M00378", {"status": "absent", "completeness": 0.0})
            },
            "hydrogenase_operon_system": hyd_stat
        }

        # Filter top microbial pathways
        top_pathways = {
            pw: info for pw, info in g_pws.items()
            if info["microbial_relevance"] == "microbial_core" and info["map_ko_fraction"] > 0.0
        }

        profile_entry = {
            "schema_version": "2.0",
            "genome_id": gid,
            "provenance": {
                "snapshot_date": manifest_data.get("snapshot_date"),
                "kegg_release": manifest_data.get("kegg_release")
            },
            "coverage_audit": g_audit,
            "tier1_module_reconstruction": {
                "complete_modules_count": sum(1 for m in g_mods.values() if m["status"] == "complete"),
                "near_complete_modules_count": sum(1 for m in g_mods.values() if m["status"] == "near_complete"),
                "partial_modules_count": sum(1 for m in g_mods.values() if m["status"] == "partial"),
                "entner_doudoroff_family_profile": ed_profile
            },
            "tier2_pathway_footprints": top_pathways,
            "tier3_brite_profiles": g_brite,
            "tier4_specialized_systems": specialized_systems
        }
        profiles.append(profile_entry)

    with open(args.out_json, "w", encoding="utf-8") as out:
        json.dump(profiles, out, indent=2)

    print(f"Schema v2.0 functional profiles generated for {len(profiles)} genomes -> {args.out_json}")

if __name__ == "__main__":
    main()
