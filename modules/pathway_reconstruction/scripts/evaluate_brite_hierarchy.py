#!/usr/bin/env python3
"""
evaluate_brite_hierarchy.py
Tier 3 BRITE Functional Hierarchy & Molecular-System Profiling Engine.

Enforces:
- Gate G5 (BRITE Multi-Membership Gate): Handles multi-label hierarchy placements without double-counting sums.
- Normalization: Reports unique_KO_count, protein_count, presence, copies_per_1000_proteins,
  and copies_per_1000_KO_annotated_proteins to prevent genome-size bias.
- Dedicated System Profiling: Extracts Transporters (ko02000), Two-component systems (ko02022),
  and Prokaryotic Defense / CRISPR-Cas (ko02048).
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

def load_brite_memberships(brite_tsv):
    """
    Load variable-depth BRITE memberships.
    ko -> list of {hierarchy_id, path_id, depth, node_name, full_path}
    """
    ko_memberships = {}
    hier_systems = {} # (hierarchy_id, system_name) -> set of KOs

    with open(brite_tsv, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        col_map = {col: i for i, col in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= len(header):
                ko = parts[col_map["ko"]]
                hier_id = parts[col_map["hierarchy_id"]]
                depth = int(parts[col_map["depth"]])
                node_name = parts[col_map["node_name"]]
                full_path_str = parts[col_map["full_path"]]
                
                try:
                    path_list = json.loads(full_path_str)
                except Exception:
                    path_list = [node_name]

                # Identify system name (typically depth 1 or 2 in specialized, or depth 2 in ko00001)
                system_name = path_list[1] if len(path_list) > 1 else path_list[0]
                subsystem_name = path_list[2] if len(path_list) > 2 else system_name

                rec = {
                    "hierarchy_id": hier_id,
                    "depth": depth,
                    "system_name": system_name,
                    "subsystem_name": subsystem_name,
                    "full_path": path_list
                }
                ko_memberships.setdefault(ko, []).append(rec)
                hier_systems.setdefault((hier_id, system_name, subsystem_name), set()).add(ko)

    return ko_memberships, hier_systems

def load_genome_evidence(ko_strict_path):
    """
    Load KO hits and protein evidence per genome.
    Returns:
      genome_ko_counts: genome_id -> {ko -> count}
      genome_proteins: genome_id -> set of protein_ids
    """
    genome_ko_counts = {}
    genome_proteins = {} # genome_id -> set of protein_ids
    if not os.path.isfile(ko_strict_path):
        return genome_ko_counts, genome_proteins

    with open_file(ko_strict_path, "rt") as f:
        header = None
        for line in f:
            sline = line.strip()
            if not sline or sline.startswith("#"):
                continue
            header = sline.split("\t")
            break
        if not header:
            return genome_ko_counts, genome_proteins

        p_idx = header.index("protein_id") if "protein_id" in header else 0
        gm_idx = header.index("genome_id") if "genome_id" in header else 2
        k_idx = header.index("ko") if "ko" in header else 3

        for line in f:
            sline = line.strip()
            if not sline or sline.startswith("#"):
                continue
            parts = sline.split("\t")
            if len(parts) > max(p_idx, gm_idx, k_idx):
                prot_id = parts[p_idx]
                gid = parts[gm_idx]
                ko = parts[k_idx].strip()
                genome_proteins.setdefault(gid, set()).add(prot_id)
                if ko and ko not in ["-", "None", "null", "NA"]:
                    genome_ko_counts.setdefault(gid, {}).setdefault(ko, 0)
                    genome_ko_counts[gid][ko] += 1

    return genome_ko_counts, genome_proteins

def main():
    parser = argparse.ArgumentParser(description="Evaluate Tier 3 BRITE Functional Hierarchy and Molecular Systems.")
    parser.add_argument("-k", "--ko_strict", required=True, help="Path to protein_KO_all_strict.tsv[.gz]")
    parser.add_argument("-b", "--brite_tsv", default="knowledge/kegg/ko_hierarchy_membership.tsv", help="Path to ko_hierarchy_membership.tsv")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 70)
    print("Tier 3: BRITE Hierarchy & Molecular-System Profiling Engine")
    print("=" * 70)

    # 1. Load BRITE Memberships
    print(f"Loading BRITE multi-label hierarchy from {args.brite_tsv}...")
    ko_memberships, hier_systems = load_brite_memberships(args.brite_tsv)
    print(f"Loaded {len(ko_memberships)} KO membership entries across {len(hier_systems)} system classes.")

    # 2. Load Genome Evidence
    print(f"Loading genome evidence from {args.ko_strict}...")
    genome_ko_counts, genome_proteins = load_genome_evidence(args.ko_strict)
    print(f"Loaded {len(genome_ko_counts)} genomes.")

    counts_tsv = os.path.join(args.out_dir, "genome_brite_counts.tsv")
    matrix_tsv = os.path.join(args.out_dir, "genome_brite_matrix.tsv")

    # 3. Calculate Normalized Counts per System
    print(f"Profiling molecular systems and calculating normalized abundances -> {counts_tsv}...")
    matrix_records = {} # (gid, system_key) -> copies_per_1000_proteins
    all_system_keys = set()

    with open(counts_tsv, "w", encoding="utf-8") as out:
        out.write(
            "genome_id\thierarchy_id\tsystem_name\tsubsystem_name\t"
            "unique_KO_count\ttotal_system_KOs\tprotein_count\tpresence\t"
            "copies_per_1000_proteins\tcopies_per_1000_KO_annotated_proteins\tobserved_kos\n"
        )

        for gid in sorted(genome_ko_counts.keys()):
            g_kos = genome_ko_counts[gid]
            n_ko_proteins = sum(g_kos.values())
            # Estimate total proteins (at least distinct KO proteins, baseline normalization)
            n_total_proteins = max(len(genome_proteins.get(gid, set())), n_ko_proteins)

            for (hier_id, sys_name, subsys_name), sys_kos in sorted(hier_systems.items()):
                # Intersect
                matched_kos = set(g_kos.keys()) & sys_kos
                u_ko_count = len(matched_kos)
                tot_sys_kos = len(sys_kos)
                presence = 1 if u_ko_count > 0 else 0

                prot_count = sum(g_kos[k] for k in matched_kos)
                
                # Gate G5: Normalized counts per 1000 proteins
                copies_per_1000_p = round((prot_count * 1000.0) / n_total_proteins, 4) if n_total_proteins > 0 else 0.0
                copies_per_1000_ko = round((prot_count * 1000.0) / n_ko_proteins, 4) if n_ko_proteins > 0 else 0.0

                ko_str = ",".join(sorted(matched_kos))

                out.write(
                    f"{gid}\t{hier_id}\t{sys_name}\t{subsys_name}\t"
                    f"{u_ko_count}\t{tot_sys_kos}\t{prot_count}\t{presence}\t"
                    f"{copies_per_1000_p:.4f}\t{copies_per_1000_ko:.4f}\t{ko_str}\n"
                )

                # Collect for wide matrix (focus on specialized systems ko02000, ko02022, ko02048)
                if hier_id in ["ko02000", "ko02022", "ko02048"]:
                    sys_key = f"{hier_id}___{sys_name}"
                    all_system_keys.add(sys_key)
                    matrix_records[(gid, sys_key)] = copies_per_1000_p

    # 4. Generate Specialized System Wide Matrix
    sorted_sys_keys = sorted(all_system_keys)
    print(f"Generating specialized system matrix ({len(sorted_sys_keys)} systems) -> {matrix_tsv}...")
    with open(matrix_tsv, "w", encoding="utf-8") as m_out:
        m_out.write("genome_id\t" + "\t".join(sorted_sys_keys) + "\n")
        for gid in sorted(genome_ko_counts.keys()):
            row = [f"{matrix_records.get((gid, k), 0.0):.4f}" for k in sorted_sys_keys]
            m_out.write(gid + "\t" + "\t".join(row) + "\n")

    print(f"Tier 3 profiling complete. Outputs saved to {args.out_dir}")

if __name__ == "__main__":
    main()
