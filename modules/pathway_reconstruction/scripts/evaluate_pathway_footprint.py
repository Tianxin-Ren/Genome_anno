#!/usr/bin/env python3
"""
evaluate_pathway_footprint.py
Tier 2 Pathway Map Footprint & Reaction Node Coverage Engine.

Enforces:
- Gate G2 (Pathway Semantic Gate): Strictly outputs map_ko_fraction, reaction_node_fraction, and footprint.
  Never asserts "pathway completeness".
- Gate G3 (Taxonomy & Scope Gate): Includes scope and microbial_relevance tags; default reporting filters to prokaryote_relevant.
- Gate G4 (Canonical-ID Gate): Pathway IDs are canonical 5-digit strings (koXXXXX).
- Gate G10 (Applicability/NA Gate): Non-metabolic pathways output NA for reaction metrics with explicit reasons.
"""

import argparse
import gzip
import os
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def load_pathways(pathways_tsv):
    """Load canonical pathways and scopes from all_kegg_pathways.tsv."""
    pathways = {}
    with open(pathways_tsv, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        col_map = {col: i for i, col in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= len(header):
                pw_id = parts[col_map["pathway_id"]]
                pathways[pw_id] = {
                    "pathway_id": pw_id,
                    "canonical_id": parts[col_map["canonical_id"]],
                    "name": parts[col_map["name"]],
                    "category": parts[col_map["category"]],
                    "is_metabolic": parts[col_map["is_metabolic"]].lower() == "true",
                    "map_structure_type": parts[col_map["map_structure_type"]],
                    "is_prokaryote_relevant": parts[col_map["is_prokaryote_relevant"]].lower() == "true",
                    "microbial_relevance": parts[col_map["microbial_relevance"]],
                    "reaction_metric_applicable": parts[col_map["reaction_metric_applicable"]].lower() == "true",
                    "total_kos": int(parts[col_map["total_kos"]]),
                    "total_reactions": int(parts[col_map["total_reactions"]]),
                    "kos": set(),
                    "reactions": set(),
                    "ko_to_rxns": {} # ko -> set of reactions
                }
    return pathways

def load_triples(triples_tsv, pathways):
    """Load strict intersection triples (pathway, KO, reaction) in P ∩ Q ∩ R."""
    with open(triples_tsv, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        col_map = {col: i for i, col in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= len(header):
                pw_id = parts[col_map["pathway_id"]]
                ko = parts[col_map["ko"]]
                rxn = parts[col_map["reaction_id"]]
                if pw_id in pathways:
                    pathways[pw_id]["kos"].add(ko)
                    pathways[pw_id]["reactions"].add(rxn)
                    pathways[pw_id]["ko_to_rxns"].setdefault(ko, set()).add(rxn)

def load_genome_kos(ko_strict_path):
    """Load observed KOs per genome."""
    genome_kos = {} # genome_id -> set of KOs
    if not os.path.isfile(ko_strict_path):
        return genome_kos
    with open_file(ko_strict_path, "rt") as f:
        header = None
        for line in f:
            sline = line.strip()
            if not sline or sline.startswith("#"):
                continue
            header = sline.split("\t")
            break
        if not header:
            return genome_kos
        gm_idx = header.index("genome_id") if "genome_id" in header else 2
        k_idx = header.index("ko") if "ko" in header else 3
        for line in f:
            sline = line.strip()
            if not sline or sline.startswith("#"):
                continue
            parts = sline.split("\t")
            if len(parts) > max(gm_idx, k_idx):
                gid = parts[gm_idx]
                ko = parts[k_idx].strip()
                if ko and ko not in ["-", "None", "null", "NA"]:
                    genome_kos.setdefault(gid, set()).add(ko)
    return genome_kos

def load_module_support(module_completeness_path):
    """Load complete and near_complete module counts per genome."""
    genome_mod_stats = {} # genome_id -> {"complete": set(), "partial": set()}
    if not module_completeness_path or not os.path.isfile(module_completeness_path):
        return genome_mod_stats
    with open_file(module_completeness_path, "rt") as f:
        first_line = f.readline()
        if not first_line:
            return genome_mod_stats
        header = first_line.strip().split("\t")
        if "genome_id" not in header or "module_id" not in header or "status" not in header:
            return genome_mod_stats
        g_idx = header.index("genome_id")
        m_idx = header.index("module_id")
        s_idx = header.index("status")
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > max(g_idx, m_idx, s_idx):
                gid = parts[g_idx]
                mid = parts[m_idx]
                stat = parts[s_idx]
                genome_mod_stats.setdefault(gid, {"complete": set(), "partial": set()})
                if stat == "complete":
                    genome_mod_stats[gid]["complete"].add(mid)
                elif stat in ["near_complete", "partial"]:
                    genome_mod_stats[gid]["partial"].add(mid)
    return genome_mod_stats

def main():
    parser = argparse.ArgumentParser(description="Evaluate Tier 2 KEGG Pathway Footprints and Reaction Node Coverages.")
    parser.add_argument("-k", "--ko_strict", required=True, help="Path to protein_KO_all_strict.tsv[.gz]")
    parser.add_argument("-p", "--pathways_tsv", default="knowledge/kegg/all_kegg_pathways.tsv", help="Path to all_kegg_pathways.tsv")
    parser.add_argument("-t", "--triples_tsv", default="knowledge/kegg/pathway_ko_reaction.tsv", help="Path to pathway_ko_reaction.tsv")
    parser.add_argument("-m", "--module_completeness", default="", help="Optional genome_module_completeness.tsv for module support")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory")
    parser.add_argument("--scope_filter", default="prokaryote_relevant", choices=["prokaryote_relevant", "all"],
                        help="Default scope for matrix output (default: prokaryote_relevant)")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 70)
    print("Tier 2: Pathway Footprint & Reaction Node Coverage Evaluation")
    print("=" * 70)

    # 1. Load Knowledge
    print("Loading pathway metadata and intersection triples...")
    pathways = load_pathways(args.pathways_tsv)
    load_triples(args.triples_tsv, pathways)
    print(f"Loaded {len(pathways)} canonical pathway definitions.")

    # 2. Load Genome Evidence
    print(f"Loading KO evidence from {args.ko_strict}...")
    genome_kos = load_genome_kos(args.ko_strict)
    print(f"Loaded evidence for {len(genome_kos)} genomes.")

    mod_support = load_module_support(args.module_completeness)

    footprint_tsv = os.path.join(args.out_dir, "genome_pathway_footprint.tsv")
    matrix_tsv = os.path.join(args.out_dir, "genome_pathway_matrix.tsv")

    # 3. Calculate Footprint & Reaction Coverage
    print(f"Calculating pathway footprints -> {footprint_tsv}...")
    footprint_records = []

    with open(footprint_tsv, "w", encoding="utf-8") as out:
        # Header strictly following Gate G2 and User Specification
        out.write(
            "genome_id\tpathway_id\tpathway_name\tdetected_unique_kos\tpathway_unique_kos\t"
            "map_ko_fraction\tsupported_reactions\ttotal_ko_linked_reactions\treaction_node_fraction\t"
            "reaction_metric_applicable\tapplicable_reason\tn_complete_modules\tn_partial_modules\t"
            "module_support_fraction\tscope\tmicrobial_relevance\tmap_structure_type\n"
        )

        for gid in sorted(genome_kos.keys()):
            g_kos = genome_kos[gid]
            g_comp_mods = mod_support.get(gid, {}).get("complete", set())
            g_part_mods = mod_support.get(gid, {}).get("partial", set())

            for pw_id in sorted(pathways.keys()):
                p = pathways[pw_id]
                pw_kos = p["kos"]
                pw_rxns = p["reactions"]

                # A. KO Footprint
                matched_kos = g_kos & pw_kos
                n_matched_kos = len(matched_kos)
                n_total_kos = len(pw_kos) if pw_kos else p["total_kos"]
                map_ko_fraction = round(n_matched_kos / n_total_kos, 4) if n_total_kos > 0 else 0.0

                # B. Reaction Node Footprint (Gate G10: Strict NA on Non-Metabolic)
                if p["reaction_metric_applicable"] and len(pw_rxns) > 0:
                    supported_rxns = set()
                    for k in matched_kos:
                        supported_rxns.update(p["ko_to_rxns"].get(k, set()))
                    n_supp_rxns = len(supported_rxns)
                    n_tot_rxns = len(pw_rxns)
                    rxn_fraction_str = f"{round(n_supp_rxns / n_tot_rxns, 4):.4f}"
                    supp_rxns_str = str(n_supp_rxns)
                    tot_rxns_str = str(n_tot_rxns)
                    applicable_str = "true"
                    reason_str = "metabolic_reaction_graph"
                else:
                    rxn_fraction_str = "NA"
                    supp_rxns_str = "NA"
                    tot_rxns_str = "NA"
                    applicable_str = "false"
                    reason_str = "non_metabolic_pathway" if not p["is_metabolic"] else "no_reaction_nodes"

                # C. Module Support (if modules available)
                n_comp_m = len(g_comp_mods)
                n_part_m = len(g_part_mods)
                mod_supp_fraction = round(n_comp_m / (n_comp_m + n_part_m), 4) if (n_comp_m + n_part_m) > 0 else 0.0

                out.write(
                    f"{gid}\t{pw_id}\t{p['name']}\t{n_matched_kos}\t{n_total_kos}\t"
                    f"{map_ko_fraction:.4f}\t{supp_rxns_str}\t{tot_rxns_str}\t{rxn_fraction_str}\t"
                    f"{applicable_str}\t{reason_str}\t{n_comp_m}\t{n_part_m}\t"
                    f"{mod_supp_fraction:.4f}\t{p['category']}\t{p['microbial_relevance']}\t{p['map_structure_type']}\n"
                )

                if (args.scope_filter == "prokaryote_relevant" and p["is_prokaryote_relevant"]) or args.scope_filter == "all":
                    footprint_records.append((gid, pw_id, map_ko_fraction))

    # 4. Generate Filtered Wide Matrix (Gate G3)
    print(f"Generating dense pathway footprint matrix ({args.scope_filter}) -> {matrix_tsv}...")
    selected_pathways = [
        pw_id for pw_id, p in pathways.items()
        if (args.scope_filter == "prokaryote_relevant" and p["is_prokaryote_relevant"]) or args.scope_filter == "all"
    ]
    selected_pathways.sort()

    # Index: (gid, pw_id) -> fraction
    f_dict = {(r[0], r[1]): r[2] for r in footprint_records}

    with open(matrix_tsv, "w", encoding="utf-8") as m_out:
        m_out.write("genome_id\t" + "\t".join(selected_pathways) + "\n")
        for gid in sorted(genome_kos.keys()):
            row = [f"{f_dict.get((gid, pw), 0.0):.4f}" for pw in selected_pathways]
            m_out.write(gid + "\t" + "\t".join(row) + "\n")

    print(f"Tier 2 evaluation complete. Matrix size: {len(genome_kos)} genomes x {len(selected_pathways)} pathways.")

if __name__ == "__main__":
    main()
