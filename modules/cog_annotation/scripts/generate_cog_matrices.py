#!/usr/bin/env python3
"""
generate_cog_matrices.py
Constructs comparative genomics matrices from protein_COG_master.tsv.gz:

1. genome_COG_protein_counts.long.tsv.gz (Protein Incidence: each protein x COG <= 1)
2. genome_COG_footprint_counts.long.tsv.gz (Footprint Frequency: counts repeated domains)
3. genome_COG_category_protein_normalized.tsv (Sum = N_annotated_proteins per genome)
4. genome_COG_category_footprint_normalized.tsv (Sum = N_retained_footprints per genome)
5. genome_COG_category_multilabel.tsv (Multi-label full counts)
6. cog_category_definition.tsv (NCBI 25 functional categories reference table)
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

# Standard NCBI 25 COG Functional Categories
NCBI_COG_CATEGORIES = [
    # INFORMATION STORAGE AND PROCESSING
    ("J", "Translation, ribosomal structure and biogenesis", "INFORMATION STORAGE AND PROCESSING"),
    ("A", "RNA processing and modification", "INFORMATION STORAGE AND PROCESSING"),
    ("K", "Transcription", "INFORMATION STORAGE AND PROCESSING"),
    ("L", "Replication, recombination and repair", "INFORMATION STORAGE AND PROCESSING"),
    ("B", "Chromatin structure and dynamics", "INFORMATION STORAGE AND PROCESSING"),
    # CELLULAR PROCESSES AND SIGNALING
    ("D", "Cell cycle control, cell division, chromosome partitioning", "CELLULAR PROCESSES AND SIGNALING"),
    ("Y", "Nuclear structure", "CELLULAR PROCESSES AND SIGNALING"),
    ("V", "Defense mechanisms", "CELLULAR PROCESSES AND SIGNALING"),
    ("T", "Signal transduction mechanisms", "CELLULAR PROCESSES AND SIGNALING"),
    ("M", "Cell wall/membrane/envelope biogenesis", "CELLULAR PROCESSES AND SIGNALING"),
    ("N", "Cell motility", "CELLULAR PROCESSES AND SIGNALING"),
    ("Z", "Cytoskeleton", "CELLULAR PROCESSES AND SIGNALING"),
    ("W", "Extracellular structures", "CELLULAR PROCESSES AND SIGNALING"),
    ("U", "Intracellular trafficking, secretion, and vesicular transport", "CELLULAR PROCESSES AND SIGNALING"),
    ("O", "Posttranslational modification, protein turnover, chaperones", "CELLULAR PROCESSES AND SIGNALING"),
    ("X", "Mobilome: prophages, transposons", "CELLULAR PROCESSES AND SIGNALING"),
    # METABOLISM
    ("C", "Energy production and conversion", "METABOLISM"),
    ("G", "Carbohydrate transport and metabolism", "METABOLISM"),
    ("E", "Amino acid transport and metabolism", "METABOLISM"),
    ("F", "Nucleotide transport and metabolism", "METABOLISM"),
    ("H", "Coenzyme transport and metabolism", "METABOLISM"),
    ("I", "Lipid transport and metabolism", "METABOLISM"),
    ("P", "Inorganic ion transport and metabolism", "METABOLISM"),
    ("Q", "Secondary metabolites biosynthesis, transport and catabolism", "METABOLISM"),
    # POORLY CHARACTERIZED
    ("R", "General function prediction only", "POORLY CHARACTERIZED"),
    ("S", "Function unknown", "POORLY CHARACTERIZED"),
]

ALL_CAT_LETTERS = [c[0] for c in NCBI_COG_CATEGORIES]

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def export_category_definitions(out_path):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("category_code\tcategory_name\tfunctional_class\n")
        for code, name, fclass in NCBI_COG_CATEGORIES:
            f.write(f"{code}\t{name}\t{fclass}\n")

def parse_categories(cat_str):
    """Clean and extract valid single-letter categories."""
    if not cat_str or cat_str == "-":
        return ["S"]
    cats = [c for c in cat_str.strip() if c in ALL_CAT_LETTERS]
    return cats if cats else ["S"]

def main():
    parser = argparse.ArgumentParser(
        description="Generate comparative genomics matrices from master COG table."
    )
    parser.add_argument("-m", "--master_tsv", required=True, help="protein_COG_master.tsv.gz")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for matrices")
    parser.add_argument("--export_dense", action="store_true", help="Also export dense genome x COG matrix")
    parser.add_argument("--category_def", default="", help="Optional path to official cog_category_definition.tsv from Gate 0")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Export reference table
    def_path = os.path.join(args.out_dir, "cog_category_definition.tsv")
    if args.category_def and os.path.isfile(args.category_def):
        import shutil
        shutil.copyfile(args.category_def, def_path)
    else:
        export_category_definitions(def_path)

    # Accumulators
    # Genome x COG protein incidence: (gid, cog_id) -> set of protein_ids
    cog_protein_incidence = defaultdict(set)
    # Genome x COG footprint count: (gid, cog_id) -> count of footprints
    cog_footprint_counts = defaultdict(int)
    # COG metadata: cog_id -> (category, name)
    cog_meta = {}

    # Category accumulators: gid -> dict of category -> float count
    cat_protein_normalized = defaultdict(lambda: defaultdict(float))
    cat_footprint_normalized = defaultdict(lambda: defaultdict(float))
    cat_multilabel = defaultdict(lambda: defaultdict(float))

    # All genomes observed
    all_genomes = set()

    # In-memory aggregation of footprints per protein to compute w_footprint
    # Since master table is sorted by genome/protein, we stream protein by protein
    print(f"INFO: Aggregating master table {args.master_tsv}...")

    current_prot_key = None # (gid, pid)
    current_prot_fps = []

    def flush_protein(gid, pid, fps):
        if not fps:
            return
        all_genomes.add(gid)
        n_fps = len(fps)
        w_footprint = 1.0 / max(1, n_fps)

        for fp in fps:
            cid = fp["cog_id"]
            cat_str = fp["cog_category"]
            cname = fp["cog_name"]

            if cid not in cog_meta:
                cog_meta[cid] = (cat_str, cname)

            cog_protein_incidence[(gid, cid)].add(pid)
            cog_footprint_counts[(gid, cid)] += 1

            cats = parse_categories(cat_str)
            k_cats = len(cats)

            # Multilabel count: +1 to each category
            for c in cats:
                cat_multilabel[gid][c] += 1.0

            # Footprint-normalized: 1.0 / k_cats to each category
            w_fp_cat = 1.0 / k_cats
            for c in cats:
                cat_footprint_normalized[gid][c] += w_fp_cat

            # Protein-normalized: w_footprint / k_cats to each category
            w_prot_cat = w_footprint / k_cats
            for c in cats:
                cat_protein_normalized[gid][c] += w_prot_cat

    with open_file(args.master_tsv, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            gid = parts[col_idx["genome_id"]].strip()
            pid = parts[col_idx["protein_id"]].strip()
            cid = parts[col_idx["cog_id"]].strip()
            cat = parts[col_idx["cog_category"]].strip() if "cog_category" in col_idx else "S"
            cname = parts[col_idx["cog_name"]].strip() if "cog_name" in col_idx else ""

            key = (gid, pid)
            if key != current_prot_key:
                if current_prot_key:
                    flush_protein(current_prot_key[0], current_prot_key[1], current_prot_fps)
                current_prot_key = key
                current_prot_fps = []

            current_prot_fps.append({
                "cog_id": cid,
                "cog_category": cat,
                "cog_name": cname
            })

        if current_prot_key:
            flush_protein(current_prot_key[0], current_prot_key[1], current_prot_fps)

    sorted_genomes = sorted(all_genomes)
    print(f"INFO: Aggregated metrics across {len(sorted_genomes)} genomes.")

    # 2. Write genome_COG_protein_counts.long.tsv (Protein Incidence)
    prot_long_path = os.path.join(args.out_dir, "genome_COG_protein_counts.long.tsv")
    with open_file(prot_long_path, "wt") as out:
        out.write("genome_id\tcog_id\tprotein_count\tfunctional_category\tcog_name\n")
        for (gid, cid), p_set in sorted(cog_protein_incidence.items()):
            cat, cname = cog_meta.get(cid, ("S", ""))
            out.write(f"{gid}\t{cid}\t{len(p_set)}\t{cat}\t{cname}\n")

    # 3. Write genome_COG_footprint_counts.long.tsv (Footprint Frequency)
    fp_long_path = os.path.join(args.out_dir, "genome_COG_footprint_counts.long.tsv")
    with open_file(fp_long_path, "wt") as out:
        out.write("genome_id\tcog_id\tfootprint_count\tfunctional_category\tcog_name\n")
        for (gid, cid), count in sorted(cog_footprint_counts.items()):
            cat, cname = cog_meta.get(cid, ("S", ""))
            out.write(f"{gid}\t{cid}\t{count}\t{cat}\t{cname}\n")

    # 4. Write Category Matrices
    def write_category_matrix(out_path, cat_dict, is_float=True):
        with open(out_path, "w", encoding="utf-8") as out:
            out.write("genome_id\t" + "\t".join(ALL_CAT_LETTERS) + "\ttotal_sum\n")
            for gid in sorted_genomes:
                row_vals = []
                row_sum = 0.0
                for c in ALL_CAT_LETTERS:
                    val = cat_dict[gid].get(c, 0.0)
                    row_sum += val
                    if is_float:
                        row_vals.append(f"{val:.4f}")
                    else:
                        row_vals.append(str(int(round(val))))
                sum_str = f"{row_sum:.4f}" if is_float else str(int(round(row_sum)))
                out.write(f"{gid}\t" + "\t".join(row_vals) + f"\t{sum_str}\n")

    # 4a. Protein-Normalized Fractional Matrix
    prot_norm_path = os.path.join(args.out_dir, "genome_COG_category_protein_normalized.tsv")
    write_category_matrix(prot_norm_path, cat_protein_normalized, is_float=True)

    # 4b. Footprint-Normalized Fractional Matrix
    fp_norm_path = os.path.join(args.out_dir, "genome_COG_category_footprint_normalized.tsv")
    write_category_matrix(fp_norm_path, cat_footprint_normalized, is_float=True)

    # 4c. Multi-label Matrix
    multi_path = os.path.join(args.out_dir, "genome_COG_category_multilabel.tsv")
    write_category_matrix(multi_path, cat_multilabel, is_float=False)

    # 5. Optional Dense Genome x COG Matrix
    if args.export_dense:
        dense_path = os.path.join(args.out_dir, "genome_COG_protein_counts.dense.tsv")
        print("INFO: Exporting dense genome x COG matrix...")
        all_cogs = sorted(cog_meta.keys())
        with open_file(dense_path, "wt") as out:
            out.write("genome_id\t" + "\t".join(all_cogs) + "\n")
            for gid in sorted_genomes:
                vals = [str(len(cog_protein_incidence.get((gid, cid), set()))) for cid in all_cogs]
                out.write(f"{gid}\t" + "\t".join(vals) + "\n")

    print(f"SUCCESS: Generated comparative matrices in {args.out_dir}")

if __name__ == "__main__":
    main()
