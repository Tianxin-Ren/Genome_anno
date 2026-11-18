#!/usr/bin/env python3
"""
expand_dedup_cog.py
High-throughput streaming expansion from deduplicated unique proteins to
original genome-level proteins.

Outputs:
1. protein_COG_master.tsv.gz:
   Comprehensive per-footprint table across all MAGs with orthogonal dimensions.
2. protein_COG_summary.tsv.gz:
   Collapsed per-protein table (1 row per protein) guaranteeing 100% proteome accounting
   (unannotated proteins are recorded as unassigned).
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

def load_assigned_footprints(shard_files_or_dir):
    """
    Load assigned footprints from per-shard TSVs.
    Returns: dict {unique_id: [footprint_dict, ...]}
    """
    files = []
    if isinstance(shard_files_or_dir, list):
        for item in shard_files_or_dir:
            if os.path.isfile(item):
                files.append(item)
            elif os.path.isdir(item):
                for f in os.listdir(item):
                    if f.endswith(".assigned.tsv") or f.endswith(".assigned.tsv.gz"):
                        files.append(os.path.join(item, f))
    elif os.path.isdir(shard_files_or_dir):
        for f in os.listdir(shard_files_or_dir):
            if f.endswith(".assigned.tsv") or f.endswith(".assigned.tsv.gz"):
                files.append(os.path.join(shard_files_or_dir, f))
    elif os.path.isfile(shard_files_or_dir):
        files.append(shard_files_or_dir)

    files.sort()
    print(f"INFO: Loading assigned footprints from {len(files)} shard file(s)...")

    unique_fps = {}
    total_fps = 0

    for fpath in files:
        with open_file(fpath, "rt") as f:
            header = f.readline().strip().split("\t")
            col_idx = {name: i for i, name in enumerate(header)}
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                uid = parts[col_idx["protein_id"]].strip()
                rec = {name: parts[idx].strip() for name, idx in col_idx.items()}
                unique_fps.setdefault(uid, []).append(rec)
                total_fps += 1

    print(f"INFO: Loaded {total_fps} footprints across {len(unique_fps)} unique proteins.")
    return unique_fps

def main():
    parser = argparse.ArgumentParser(
        description="Expand unique protein COG assignments to all MAG proteins."
    )
    parser.add_argument("-s", "--shards_input", required=True, nargs="+", help="Directory or list of shard assigned TSVs")
    parser.add_argument("-m", "--mapping", required=True, help="protein_to_unique.tsv.gz from deduplication")
    parser.add_argument("-p", "--provenance_report", default="", help="Gate 0 provenance report JSON")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for master tables")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Load method designation and provenance
    method_desig = "CDD Cog_LE-based COG annotation mapped to COG2024 definitions"
    if args.provenance_report and os.path.isfile(args.provenance_report):
        try:
            with open(args.provenance_report, "r", encoding="utf-8") as pf:
                p_data = json.load(pf)
                method_desig = p_data.get("method_designation", method_desig)
        except Exception:
            pass

    # Load unique footprints
    unique_fps = load_assigned_footprints(args.shards_input)

    master_tsv_path = os.path.join(args.out_dir, "protein_COG_master.tsv")
    summary_tsv_path = os.path.join(args.out_dir, "protein_COG_summary.tsv")

    master_cols = [
        "genome_id", "protein_id", "unique_protein_id", "cog_id", "pssm_id", "cog_category", "cog_name",
        "qstart", "qend", "qlen", "sstart", "send", "slen", "aln_len",
        "query_cov", "profile_cov", "evalue", "bitscore",
        "specific_bitscore_threshold", "passes_specific_threshold",
        "cog_match_code", "cog_match_class", "profile_coverage_class", "query_coverage_class",
        "footprint_rank", "region_id", "region_status", "architecture_class",
        "has_competing_profiles", "has_ambiguous_region", "is_multi_cog", "is_repeated_cog",
        "top2_cog", "delta_bitscore", "bitscore_ratio", "delta_bitscore_fraction",
        "likely_truncated_gene", "truncation_evidence_count", "prodigal_partial", "near_contig_edge", "short_vs_cog",
        "cog_footprint_architecture", "method_designation"
    ]

    summary_cols = [
        "genome_id", "protein_id", "unique_protein_id", "aa_length",
        "n_footprints", "primary_cog", "cog_footprint_architecture", "architecture_class",
        "cog_match_codes", "has_competing_profiles", "has_ambiguous_region",
        "is_multi_cog", "is_repeated_cog", "likely_truncated_gene", "method_designation"
    ]

    print(f"INFO: Streaming expansion using mapping {args.mapping}...")
    total_proteins = 0
    annotated_proteins = 0
    total_expanded_footprints = 0

    with open_file(args.mapping, "rt") as map_f, \
         open_file(master_tsv_path, "wt") as master_out, \
         open_file(summary_tsv_path, "wt") as summary_out:

        master_out.write("\t".join(master_cols) + "\n")
        summary_out.write("\t".join(summary_cols) + "\n")

        # Skip header: protein_id \t unique_id \t aa_length \t genome_id
        map_header = map_f.readline().strip().split("\t")
        col_idx = {c: i for i, c in enumerate(map_header)}

        for line in map_f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            orig_pid = parts[col_idx["protein_id"]].strip()
            uid = parts[col_idx["unique_id"]].strip()
            aa_len = parts[col_idx["aa_length"]].strip() if "aa_length" in col_idx else "0"
            gid = parts[col_idx["genome_id"]].strip() if "genome_id" in col_idx else "unknown"

            total_proteins += 1
            fps = unique_fps.get(uid)

            if fps:
                annotated_proteins += 1
                total_expanded_footprints += len(fps)

                # Write master rows for each footprint
                for fp in fps:
                    row = [
                        gid, orig_pid, uid, fp["cog_id"], fp["pssm_id"], fp.get("cog_category", "S"), fp.get("cog_name", ""),
                        fp["qstart"], fp["qend"], fp["qlen"], fp["sstart"], fp["send"], fp["slen"], fp["aln_len"],
                        fp["query_cov"], fp["profile_cov"], fp["evalue"], fp["bitscore"],
                        fp.get("specific_bitscore_threshold", "NA"), fp.get("passes_specific_threshold", "NA"),
                        fp["cog_match_code"], fp["cog_match_class"], fp["profile_coverage_class"], fp["query_coverage_class"],
                        fp["footprint_rank"], fp["region_id"], fp["region_status"], fp["architecture_class"],
                        fp["has_competing_profiles"], fp["has_ambiguous_region"], fp["is_multi_cog"], fp["is_repeated_cog"],
                        fp["top2_cog"], fp["delta_bitscore"], fp["bitscore_ratio"], fp["delta_bitscore_fraction"],
                        fp["likely_truncated_gene"], fp["truncation_evidence_count"], fp["prodigal_partial"], fp["near_contig_edge"], fp["short_vs_cog"],
                        fp["cog_footprint_architecture"], method_desig
                    ]
                    master_out.write("\t".join(str(x) for x in row) + "\n")

                # Write summary row
                primary_cog = fps[0]["cog_id"]
                arch_str = fps[0]["cog_footprint_architecture"]
                arch_cls = fps[0]["architecture_class"]
                match_codes = ",".join(str(fp["cog_match_code"]) for fp in fps)
                has_comp = fps[0]["has_competing_profiles"]
                has_ambig = fps[0]["has_ambiguous_region"]
                is_multi = fps[0]["is_multi_cog"]
                is_rep = fps[0]["is_repeated_cog"]
                likely_trunc = any(str(fp["likely_truncated_gene"]).lower() == "true" for fp in fps)

                sum_row = [
                    gid, orig_pid, uid, aa_len,
                    len(fps), primary_cog, arch_str, arch_cls,
                    match_codes, has_comp, has_ambig,
                    is_multi, is_rep, likely_trunc, method_desig
                ]
                summary_out.write("\t".join(str(x) for x in sum_row) + "\n")

            else:
                # Unassigned protein accounting
                sum_row = [
                    gid, orig_pid, uid, aa_len,
                    0, "unassigned", "none", "unassigned",
                    "none", "False", "False",
                    "False", "False", "False", method_desig
                ]
                summary_out.write("\t".join(str(x) for x in sum_row) + "\n")

    print(
        f"SUCCESS: Expanded {total_proteins} total proteins ({annotated_proteins} annotated, "
        f"{total_expanded_footprints} total footprints). "
        f"Output written to:\n  - {master_tsv_path}\n  - {summary_tsv_path}"
    )

if __name__ == "__main__":
    main()
