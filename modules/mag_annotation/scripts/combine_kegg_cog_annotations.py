#!/usr/bin/env python3
"""
combine_kegg_cog_annotations.py
High-throughput streaming integration of KEGG (PyKofamSearch) and COG (COG2024 / Cog_LE) annotations.

Inputs:
  --protein_mapping (-p): protein_to_unique.tsv.gz (complete proteome catalog)
  --ko_tsv (-k): protein_KO_best.tsv (single best KO per protein from KOfam)
  --cog_tsv (-c): protein_COG_summary.tsv (collapsed 1-row-per-protein summary from COG engine)
  --cat_def: cog_category_definition.tsv (category code -> name & functional class)
  --pssm_mapping (-m): cdd_pssm_mapping.tsv (optional, cog_id -> name & category)
  --out_dir (-o): destination directory for combined deliverables

Outputs:
  1. protein_annotation_master.tsv:
     100% complete proteome master table with single-row-per-protein dual annotations:
     genome_id | protein_id | gene_id | aa_length | ko_id | ko_definition | ko_bitscore | ko_evalue | ko_is_strict |
     cog_id | cog_name | cog_category | cog_category_desc | cog_functional_class | cog_match_codes | cog_architecture |
     annotation_status
     where annotation_status in {"both_annotated", "ko_only", "cog_only", "unassigned"}.

  2. genome_annotation_summary.tsv:
     Per-genome macro statistics:
     genome_id | total_proteins | annotated_proteins | dual_annotated_proteins | ko_only_proteins |
     cog_only_proteins | unannotated_proteins | overall_annotation_rate | ko_annotation_rate |
     cog_annotation_rate | top_cog_categories
"""

import argparse
import gzip
import os
import sys
import time
from collections import Counter, defaultdict


def open_file(path, mode="rt"):
    """Transparently opens plain or gzipped files."""
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)


def load_category_definitions(cat_def_path):
    """
    Loads category_code -> (category_name, functional_class) from cog_category_definition.tsv.
    """
    cat_map = {}
    if not cat_def_path or not os.path.isfile(cat_def_path):
        return cat_map

    with open_file(cat_def_path, "rt") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                code = parts[0].strip().upper()
                name = parts[1].strip()
                fclass = parts[2].strip()
                cat_map[code] = (name, fclass)
    return cat_map


def load_pssm_mapping(pssm_mapping_path):
    """
    Loads cog_id -> (cog_name, cog_category) from cdd_pssm_mapping.tsv.
    """
    cog_info = {}
    if not pssm_mapping_path or not os.path.isfile(pssm_mapping_path):
        return cog_info

    with open_file(pssm_mapping_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            cid = parts[col_idx["cog_id"]].strip().upper() if "cog_id" in col_idx else ""
            if not cid:
                continue
            cname = parts[col_idx["cog_name"]].strip() if "cog_name" in col_idx else ""
            cat = parts[col_idx["cog_category"]].strip() if "cog_category" in col_idx else "S"
            if cid not in cog_info or (cname and not cog_info[cid][0]):
                cog_info[cid] = (cname, cat)
    return cog_info


def load_kegg_annotations(ko_tsv_path):
    """
    Loads protein_KO_best.tsv into a compact dictionary.
    Returns: dict {protein_id: (gene_id, ko, ko_definition, bitscore, evalue, is_strict)}
    """
    ko_map = {}
    if not ko_tsv_path or not os.path.isfile(ko_tsv_path):
        print(f"WARNING: KEGG best KO file not found at: {ko_tsv_path}. Proceeding with COG only.")
        return ko_map

    print(f"INFO: Loading KEGG annotations from {ko_tsv_path}...")
    t0 = time.time()
    count = 0
    with open_file(ko_tsv_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        p_idx = col_idx.get("protein_id", 0)
        g_idx = col_idx.get("gene_id", 1)
        k_idx = col_idx.get("ko", 3)
        d_idx = col_idx.get("ko_definition", 4)
        b_idx = col_idx.get("bitscore", 5)
        e_idx = col_idx.get("evalue", 8)
        s_idx = col_idx.get("is_strict", 10)

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            pid = parts[p_idx].strip()
            gid = parts[g_idx].strip() if len(parts) > g_idx else pid
            ko = parts[k_idx].strip() if len(parts) > k_idx else "unassigned"
            kdef = parts[d_idx].strip() if len(parts) > d_idx else "none"
            bitscore = parts[b_idx].strip() if len(parts) > b_idx else "NA"
            evalue = parts[e_idx].strip() if len(parts) > e_idx else "NA"
            is_strict = parts[s_idx].strip() if len(parts) > s_idx else "False"

            ko_map[pid] = (gid, ko, kdef, bitscore, evalue, is_strict)
            count += 1

    print(f"INFO: Loaded {count:,} KEGG KO annotations in {time.time() - t0:.2f}s.")
    return ko_map


def infer_gene_id(protein_id):
    """Infers standard gene_id from structured protein_id."""
    if "___prot_" in protein_id:
        # e.g. GCA_000001___prot_1 -> GCA_000001_1
        prefix, suffix = protein_id.split("___prot_", 1)
        return f"{prefix}_{suffix}"
    if "___" in protein_id:
        return protein_id.replace("___", "_")
    return protein_id


def main():
    parser = argparse.ArgumentParser(
        description="Integrate KEGG and COG annotations into unified master and genome summary deliverables."
    )
    parser.add_argument("-p", "--protein_mapping", required=True,
                        help="protein_to_unique.tsv.gz (complete proteome catalog)")
    parser.add_argument("-k", "--ko_tsv", required=False, default="",
                        help="protein_KO_best.tsv from KOfam merge")
    parser.add_argument("-c", "--cog_tsv", required=True,
                        help="protein_COG_summary.tsv from COG expansion")
    parser.add_argument("--cat_def", required=False, default="",
                        help="cog_category_definition.tsv from Gate 0 audit")
    parser.add_argument("-m", "--pssm_mapping", required=False, default="",
                        help="cdd_pssm_mapping.tsv from Gate 0 audit")
    parser.add_argument("-o", "--out_dir", required=True,
                        help="Destination directory for combined deliverables")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    t_start = time.time()
    print("==============================================================================")
    print("KEGG + COG Annotation Integration Engine")
    print("==============================================================================")

    # 1. Load Reference Dictionaries
    cat_map = load_category_definitions(args.cat_def)
    cog_info = load_pssm_mapping(args.pssm_mapping)
    ko_map = load_kegg_annotations(args.ko_tsv)

    # 2. Output Paths
    master_path = os.path.join(args.out_dir, "protein_annotation_master.tsv")
    summary_path = os.path.join(args.out_dir, "genome_annotation_summary.tsv")

    master_header = [
        "genome_id", "protein_id", "gene_id", "aa_length",
        "ko_id", "ko_definition", "ko_bitscore", "ko_evalue", "ko_is_strict",
        "cog_id", "cog_name", "cog_category", "cog_category_desc", "cog_functional_class",
        "cog_match_codes", "cog_architecture", "annotation_status"
    ]

    summary_header = [
        "genome_id", "total_proteins", "annotated_proteins", "dual_annotated_proteins",
        "ko_only_proteins", "cog_only_proteins", "unannotated_proteins",
        "overall_annotation_rate", "ko_annotation_rate", "cog_annotation_rate",
        "top_cog_categories"
    ]

    # Per-genome metrics accumulators
    genome_totals = Counter()
    genome_dual = Counter()
    genome_ko_only = Counter()
    genome_cog_only = Counter()
    genome_unannot = Counter()
    genome_cog_cats = defaultdict(Counter)

    BUFFER_SIZE = 50000
    out_buffer = []

    print(f"INFO: Streaming COG summary from {args.cog_tsv}...")
    total_processed = 0
    annotated_count = 0

    with open(master_path, "w", encoding="utf-8") as master_out:
        master_out.write("\t".join(master_header) + "\n")

        with open_file(args.cog_tsv, "rt") as cog_f:
            c_header = cog_f.readline().strip().split("\t")
            c_idx = {name: i for i, name in enumerate(c_header)}

            cg_gid = c_idx.get("genome_id", 0)
            cg_pid = c_idx.get("protein_id", 1)
            cg_len = c_idx.get("aa_length", 3)
            cg_pcog = c_idx.get("primary_cog", 5)
            cg_arch = c_idx.get("cog_footprint_architecture", 6)
            cg_code = c_idx.get("cog_match_codes", 8)

            for line in cog_f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                gid = parts[cg_gid].strip()
                pid = parts[cg_pid].strip()
                aa_len = parts[cg_len].strip() if len(parts) > cg_len else "0"
                primary_cog = parts[cg_pcog].strip() if len(parts) > cg_pcog else "unassigned"
                cog_arch = parts[cg_arch].strip() if len(parts) > cg_arch else "none"
                cog_codes = parts[cg_code].strip() if len(parts) > cg_code else "none"

                # 1. Resolve KEGG hit
                if pid in ko_map:
                    k_gene, ko_id, ko_def, ko_bit, ko_eval, ko_strict = ko_map[pid]
                else:
                    k_gene = infer_gene_id(pid)
                    ko_id = "unassigned"
                    ko_def = "none"
                    ko_bit = "NA"
                    ko_eval = "NA"
                    ko_strict = "False"

                # 2. Resolve COG metadata
                has_cog = (primary_cog != "unassigned" and primary_cog != "none" and primary_cog.startswith("COG"))
                has_ko = (ko_id != "unassigned" and ko_id != "none" and ko_id.startswith("K"))

                if has_cog:
                    cname, ccat = cog_info.get(primary_cog, ("none", "S"))
                    cat_letters = [c for c in ccat if c.isalpha()] if ccat else ["S"]
                    descs = [cat_map[c][0] for c in cat_letters if c in cat_map]
                    fclasses = [cat_map[c][1] for c in cat_letters if c in cat_map]

                    cog_name = cname if cname else "none"
                    cog_category = ccat if ccat else "S"
                    cog_cat_desc = "; ".join(descs) if descs else "none"
                    cog_func_class = "; ".join(dict.fromkeys(fclasses)) if fclasses else "none"
                else:
                    cog_name = "none"
                    cog_category = "none"
                    cog_cat_desc = "none"
                    cog_func_class = "none"
                    cog_arch = "none"
                    cog_codes = "none"

                # 3. Determine joint annotation status
                if has_ko and has_cog:
                    status = "both_annotated"
                    genome_dual[gid] += 1
                    annotated_count += 1
                elif has_ko and not has_cog:
                    status = "ko_only"
                    genome_ko_only[gid] += 1
                    annotated_count += 1
                elif not has_ko and has_cog:
                    status = "cog_only"
                    genome_cog_only[gid] += 1
                    annotated_count += 1
                else:
                    status = "unassigned"
                    genome_unannot[gid] += 1

                genome_totals[gid] += 1
                if has_cog and cog_category != "none":
                    for c in cog_category:
                        if c.isalpha():
                            genome_cog_cats[gid][c] += 1

                # Row construction
                row = [
                    gid, pid, k_gene, aa_len,
                    ko_id, ko_def, ko_bit, ko_eval, ko_strict,
                    primary_cog, cog_name, cog_category, cog_cat_desc, cog_func_class,
                    cog_codes, cog_arch, status
                ]
                out_buffer.append("\t".join(str(x) for x in row) + "\n")
                total_processed += 1

                if len(out_buffer) >= BUFFER_SIZE:
                    master_out.write("".join(out_buffer))
                    out_buffer.clear()

        if out_buffer:
            master_out.write("".join(out_buffer))
            out_buffer.clear()

    print(f"INFO: Wrote {total_processed:,} proteins to {master_path}.")

    # 4. Generate genome_annotation_summary.tsv
    print(f"INFO: Generating genome macro summary for {len(genome_totals)} genomes...")
    with open(summary_path, "w", encoding="utf-8") as summary_out:
        summary_out.write("\t".join(summary_header) + "\n")

        for gid in sorted(genome_totals.keys()):
            tot = genome_totals[gid]
            dual = genome_dual[gid]
            ko_only = genome_ko_only[gid]
            cog_only = genome_cog_only[gid]
            unannot = genome_unannot[gid]
            annot = dual + ko_only + cog_only

            overall_rate = annot / max(1, tot)
            ko_rate = (dual + ko_only) / max(1, tot)
            cog_rate = (dual + cog_only) / max(1, tot)

            # Top categories (e.g. "J:150,E:120,C:95")
            top_cats = genome_cog_cats[gid].most_common(5)
            top_cats_str = ",".join(f"{cat}:{cnt}" for cat, cnt in top_cats) if top_cats else "none"

            s_row = [
                gid, tot, annot, dual, ko_only, cog_only, unannot,
                f"{overall_rate:.4f}", f"{ko_rate:.4f}", f"{cog_rate:.4f}",
                top_cats_str
            ]
            summary_out.write("\t".join(str(x) for x in s_row) + "\n")

    elapsed = time.time() - t_start
    print("==============================================================================")
    print(f"SUCCESS: Integrated {total_processed:,} proteins across {len(genome_totals):,} genomes in {elapsed:.2f}s.")
    print(f"  - Annotated proteins: {annotated_count:,} ({annotated_count/max(1, total_processed)*100:.2f}%)")
    print(f"  - Master table: {master_path}")
    print(f"  - Genome summary: {summary_path}")
    print("==============================================================================")


if __name__ == "__main__":
    main()
