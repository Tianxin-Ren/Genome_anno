#!/usr/bin/env python3
"""
audit_genome_coverage.py
Multi-Tier Functional Coverage Audit Engine.

Enforces:
- Gate G6 (Proteome Coverage Gate): Strictly separates:
  1) protein_ko_annotation_rate = proteins_with_KO / total_predicted_proteins
  2) ko00001_tracking_rate = detected_KOs_in_ko00001 / detected_unique_KOs
- True Orphan KO Definition:
  Identifies true orphan KOs as KOs without module, pathway, or specialized BRITE context,
  strictly isolating unclassified KOs from umbrella 09190 tracking.
"""

import argparse
import gzip
import json
import os
import re
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def load_reference_sets(modules_tsv, pathways_tsv, brite_tsv):
    """Load reference KO universes across tiers."""
    module_kos = set()
    if os.path.isfile(modules_tsv):
        with open_file(modules_tsv, "rt") as f:
            for line in f:
                module_kos.update(re.findall(r'K\d{5}', line))

    pathway_kos = set()
    if os.path.isfile(pathways_tsv):
        # We can extract from all_kegg_pathways.tsv or pathway_ko_reaction.tsv
        with open_file(pathways_tsv, "rt") as f:
            header = f.readline().strip().split("\t")
            # If all_kegg_pathways.tsv, check if pathway_ko_reaction exists in same dir
            triples_path = os.path.join(os.path.dirname(pathways_tsv), "pathway_ko_reaction.tsv")
            if os.path.isfile(triples_path):
                with open_file(triples_path, "rt") as tf:
                    tf.readline()
                    for tline in tf:
                        parts = tline.strip().split("\t")
                        if len(parts) >= 3:
                            pathway_kos.add(parts[2])

    ko00001_kos = set()
    specialized_brite_kos = set()
    kos_in_09190 = set()

    if os.path.isfile(brite_tsv):
        with open_file(brite_tsv, "rt") as f:
            header = f.readline().strip().split("\t")
            k_idx = header.index("ko") if "ko" in header else 0
            h_idx = header.index("hierarchy_id") if "hierarchy_id" in header else 1
            p_idx = header.index("full_path") if "full_path" in header else -1

            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    ko = parts[k_idx]
                    hid = parts[h_idx]
                    if hid == "ko00001":
                        ko00001_kos.add(ko)
                        if p_idx != -1 and p_idx < len(parts) and "09190 Not Included in Pathway or Brite" in parts[p_idx]:
                            kos_in_09190.add(ko)
                    else:
                        specialized_brite_kos.add(ko)

    return {
        "module_kos": module_kos,
        "pathway_kos": pathway_kos,
        "ko00001_kos": ko00001_kos,
        "specialized_brite_kos": specialized_brite_kos,
        "kos_in_09190": kos_in_09190
    }

def load_genome_data(ko_strict_tsv):
    """
    Parse protein_KO_all_strict.tsv.
    Returns:
      genome_proteins: genome_id -> set of protein_ids
      genome_kos: genome_id -> set of unique KOs
    """
    genome_proteins = {}
    genome_kos = {}
    if not os.path.isfile(ko_strict_tsv):
        return genome_proteins, genome_kos

    with open_file(ko_strict_tsv, "rt") as f:
        header = None
        for line in f:
            sline = line.strip()
            if not sline or sline.startswith("#"):
                continue
            header = sline.split("\t")
            break
        if not header:
            return genome_proteins, genome_kos

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
                    genome_kos.setdefault(gid, set()).add(ko)

    return genome_proteins, genome_kos

def load_total_predicted_proteins(metadata_tsv):
    """Load total predicted proteins from protein_metadata or dedup_stats if available."""
    total_counts = {}
    if metadata_tsv and os.path.isfile(metadata_tsv):
        with open_file(metadata_tsv, "rt") as f:
            header = f.readline().strip().split("\t")
            gm_idx = header.index("genome_id") if "genome_id" in header else -1
            if gm_idx != -1:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) > gm_idx:
                        gid = parts[gm_idx]
                        total_counts[gid] = total_counts.get(gid, 0) + 1
    return total_counts

def main():
    parser = argparse.ArgumentParser(description="Audit proteome-to-KO and KO-to-hierarchy functional coverage.")
    parser.add_argument("-k", "--ko_strict", required=True, help="Path to protein_KO_all_strict.tsv[.gz]")
    parser.add_argument("-m", "--modules_tsv", default="knowledge/kegg/all_kegg_modules.tsv", help="Path to all_kegg_modules.tsv")
    parser.add_argument("-p", "--pathways_tsv", default="knowledge/kegg/all_kegg_pathways.tsv", help="Path to all_kegg_pathways.tsv")
    parser.add_argument("-b", "--brite_tsv", default="knowledge/kegg/ko_hierarchy_membership.tsv", help="Path to ko_hierarchy_membership.tsv")
    parser.add_argument("--protein_metadata", default="", help="Optional protein_metadata.tsv.gz to get total proteome size")
    parser.add_argument("-o", "--out_tsv", required=True, help="Output path for genome_ko_coverage_audit.tsv")

    args = parser.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)

    print("=" * 70)
    print("Multi-Tier Functional Coverage Audit (Gate G6 Enforcement)")
    print("=" * 70)

    print("Loading reference knowledge sets...")
    ref = load_reference_sets(args.modules_tsv, args.pathways_tsv, args.brite_tsv)
    print(f"Reference Module KOs: {len(ref['module_kos'])}")
    print(f"Reference Pathway KOs: {len(ref['pathway_kos'])}")
    print(f"Reference ko00001 KOs: {len(ref['ko00001_kos'])}")
    print(f"Reference Specialized BRITE KOs: {len(ref['specialized_brite_kos'])}")
    print(f"Reference KOs in 09190 (Not Included in Pathway/Brite): {len(ref['kos_in_09190'])}")

    print(f"Loading genome KO evidence from {args.ko_strict}...")
    genome_prots, genome_kos = load_genome_data(args.ko_strict)
    total_predicted = load_total_predicted_proteins(args.protein_metadata)

    with open(args.out_tsv, "w", encoding="utf-8") as out:
        out.write(
            "genome_id\t"
            "total_predicted_proteins\tproteins_with_ko\tprotein_ko_annotation_rate\t"
            "detected_unique_kos\t"
            "kos_in_modules\tmodule_assignment_rate\t"
            "kos_in_pathways\tpathway_assignment_rate\t"
            "kos_in_specialized_brite\tspecialized_brite_assignment_rate\t"
            "kos_in_ko00001_09190\tko00001_tracking_rate\t"
            "true_orphan_kos\ttrue_orphan_rate\t"
            "orphan_ko_list\n"
        )

        for gid in sorted(genome_kos.keys()):
            g_kos = genome_kos[gid]
            n_det_kos = len(g_kos)

            prots_with_ko = len(genome_prots.get(gid, set()))
            tot_pred = total_predicted.get(gid, prots_with_ko)
            prot_anno_rate = round(prots_with_ko / tot_pred, 4) if tot_pred > 0 else 0.0

            # Intersections
            in_mod = g_kos & ref["module_kos"]
            in_pw = g_kos & ref["pathway_kos"]
            in_spec = g_kos & ref["specialized_brite_kos"]
            in_09190 = g_kos & ref["kos_in_09190"]
            in_ko00001 = g_kos & ref["ko00001_kos"]

            # True orphan KOs: KOs without module, without pathway, and without specialized BRITE
            true_orphans = g_kos - (ref["module_kos"] | ref["pathway_kos"] | ref["specialized_brite_kos"])

            # Rates
            mod_rate = round(len(in_mod) / n_det_kos, 4) if n_det_kos > 0 else 0.0
            pw_rate = round(len(in_pw) / n_det_kos, 4) if n_det_kos > 0 else 0.0
            spec_rate = round(len(in_spec) / n_det_kos, 4) if n_det_kos > 0 else 0.0
            ko00001_rate = round(len(in_ko00001) / n_det_kos, 4) if n_det_kos > 0 else 0.0
            orphan_rate = round(len(true_orphans) / n_det_kos, 4) if n_det_kos > 0 else 0.0

            orphan_list_str = ",".join(sorted(true_orphans))

            out.write(
                f"{gid}\t"
                f"{tot_pred}\t{prots_with_ko}\t{prot_anno_rate:.4f}\t"
                f"{n_det_kos}\t"
                f"{len(in_mod)}\t{mod_rate:.4f}\t"
                f"{len(in_pw)}\t{pw_rate:.4f}\t"
                f"{len(in_spec)}\t{spec_rate:.4f}\t"
                f"{len(in_09190)}\t{ko00001_rate:.4f}\t"
                f"{len(true_orphans)}\t{orphan_rate:.4f}\t"
                f"{orphan_list_str}\n"
            )

    print(f"Coverage audit complete across {len(genome_kos)} genomes -> {args.out_tsv}")

if __name__ == "__main__":
    main()
