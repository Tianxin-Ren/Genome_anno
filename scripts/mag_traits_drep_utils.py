#!/usr/bin/env python3
"""
scripts/mag_traits_drep_utils.py

Helper utilities for dRep block-wise processing in mag_traits module:
1. pool-manifest: Scans dereplicated genome directories from all blocks and generates
   a combined TSV manifest (mag_id \\t abs_path).
2. merge-genome-info: Merges multiple genomeInfo.csv tables into one, deduplicating headers.
3. filter-genome-info: Filters genomeInfo.csv to retain only round 1 winner genomes for round 2.
"""

import argparse
import csv
import glob
import os
import re
import sys


def cmd_pool_manifest(args):
    """
    Scans block dereplicated_genomes directories and creates a pooled manifest TSV.
    """
    out_manifest = os.path.abspath(args.out_manifest)
    out_dir = os.path.dirname(out_manifest)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    seen_mags = set()
    rows = []
    skipped_dups = 0

    for b_dir in args.block_dirs:
        if not b_dir or not os.path.isdir(b_dir):
            continue
        for ext in ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz"):
            for fa_path in sorted(glob.glob(os.path.join(b_dir, ext))):
                stem = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", os.path.basename(fa_path))
                if not stem:
                    continue
                if stem in seen_mags:
                    skipped_dups += 1
                    continue
                seen_mags.add(stem)
                rows.append((stem, os.path.abspath(fa_path)))

    with open(out_manifest, "w", encoding="utf-8") as f_out:
        for stem, abs_path in rows:
            f_out.write(f"{stem}\t{abs_path}\n")

    print(
        f"[pool-manifest] Staged {len(rows)} genomes into manifest: {out_manifest} "
        f"(skipped duplicate IDs: {skipped_dups})",
        file=sys.stderr,
    )


def cmd_merge_genome_info(args):
    """
    Merges multiple genomeInfo.csv files into one destination CSV.
    """
    out_path = os.path.abspath(args.out)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    header_written = False
    total_rows = 0

    with open(out_path, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.writer(out_f)
        for src in args.inputs:
            if not src or not os.path.exists(src) or os.path.getsize(src) == 0:
                continue
            with open(src, "r", encoding="utf-8", errors="ignore") as in_f:
                reader = csv.reader(in_f)
                header = next(reader, None)
                if header and not header_written:
                    writer.writerow(header)
                    header_written = True
                for row in reader:
                    if row:
                        writer.writerow(row)
                        total_rows += 1

    print(
        f"[merge-genome-info] Merged {len(args.inputs)} tables ({total_rows} rows) into {out_path}",
        file=sys.stderr,
    )


def cmd_filter_genome_info(args):
    """
    Filters round 1 genomeInfo.csv to include only genomes that survived dereplication.
    """
    dst_path = os.path.abspath(args.dst)
    dst_dir = os.path.dirname(dst_path)
    if dst_dir:
        os.makedirs(dst_dir, exist_ok=True)

    derep_dir = args.derep_dir
    src_path = args.src

    if not src_path or not os.path.exists(src_path) or os.path.getsize(src_path) == 0:
        # Create empty destination file
        open(dst_path, "w", encoding="utf-8").close()
        print(
            f"[filter-genome-info] Source genomeInfo missing or empty. Wrote empty '{dst_path}'.",
            file=sys.stderr,
        )
        return

    winners = set()
    winner_stems = set()
    if os.path.isdir(derep_dir):
        for w in os.listdir(derep_dir):
            if w.startswith("."):
                continue
            winners.add(w)
            stem = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", w)
            winner_stems.add(stem)

    rows = []
    header = None
    with open(src_path, "r", encoding="utf-8", errors="ignore") as in_f:
        reader = csv.reader(in_f)
        header = next(reader, None)
        for r in reader:
            if not r:
                continue
            g_raw = r[0].strip()
            g_stem = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", g_raw)
            if g_raw in winners or g_stem in winner_stems:
                fa_name = f"{g_stem}.fa" if f"{g_stem}.fa" in winners else g_raw
                rows.append([fa_name] + r[1:])

    with open(dst_path, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.writer(out_f)
        if header:
            writer.writerow(header)
        writer.writerows(rows)

    print(
        f"[filter-genome-info] Filtered {len(rows)} winner entries into {dst_path} "
        f"(out of {len(winners)} dereplicated genomes)",
        file=sys.stderr,
    )


def cmd_union_blocks(args):
    """
    Unions dereplicated outputs from all disjoint taxonomy blocks:
    1. Symlinks/copies all genomes from each block's dereplicated_genomes/ into out_derep/
    2. Merges Cdb.csv, prefixing cluster IDs with block_id to prevent collision
    3. Merges Widb.csv
    4. Merges genomeInfo.csv
    """
    import shutil

    out_derep = os.path.abspath(args.out_derep)
    out_cdb = os.path.abspath(args.out_cdb)
    out_widb = os.path.abspath(args.out_widb)
    out_ginfo = os.path.abspath(args.out_ginfo)

    os.makedirs(out_derep, exist_ok=True)
    os.makedirs(os.path.dirname(out_cdb), exist_ok=True)
    os.makedirs(os.path.dirname(out_widb), exist_ok=True)
    os.makedirs(os.path.dirname(out_ginfo), exist_ok=True)

    block_ids = []
    if os.path.isfile(args.block_list):
        with open(args.block_list, "r", encoding="utf-8") as f:
            block_ids = [line.strip() for line in f if line.strip()]

    total_genomes = 0
    cdb_header_written = False
    widb_header_written = False
    ginfo_header_written = False

    with open(out_cdb, "w", newline="", encoding="utf-8") as cdb_out, \
         open(out_widb, "w", newline="", encoding="utf-8") as widb_out, \
         open(out_ginfo, "w", newline="", encoding="utf-8") as ginfo_out:

        cdb_writer = csv.writer(cdb_out)
        widb_writer = csv.writer(widb_out)
        ginfo_writer = csv.writer(ginfo_out)

        for b_id in block_ids:
            b_dir = os.path.join(args.blocks_dir, b_id)
            derep_src = os.path.join(b_dir, "dereplicated_genomes")
            if os.path.isdir(derep_src):
                for fname in os.listdir(derep_src):
                    if fname.startswith("."):
                        continue
                    src_fa = os.path.join(derep_src, fname)
                    dst_fa = os.path.join(out_derep, fname)
                    if os.path.islink(dst_fa) or os.path.exists(dst_fa):
                        try:
                            os.unlink(dst_fa)
                        except OSError:
                            pass
                    try:
                        os.symlink(os.path.abspath(src_fa), dst_fa)
                    except OSError:
                        shutil.copyfile(src_fa, dst_fa)
                    total_genomes += 1

            # Merge Cdb
            cdb_path = os.path.join(b_dir, "data_tables", "Cdb.csv")
            if os.path.isfile(cdb_path) and os.path.getsize(cdb_path) > 0:
                with open(cdb_path, "r", encoding="utf-8", errors="ignore") as in_f:
                    reader = csv.reader(in_f)
                    h = next(reader, None)
                    if h and not cdb_header_written:
                        cdb_writer.writerow(h)
                        cdb_header_written = True
                    cluster_cols = []
                    if h:
                        for idx, col in enumerate(h):
                            if "cluster" in col.lower():
                                cluster_cols.append(idx)
                    if not cluster_cols and h and len(h) > 1:
                        cluster_cols = [1]
                    for row in reader:
                        if not row:
                            continue
                        for c_idx in cluster_cols:
                            if len(row) > c_idx and row[c_idx]:
                                row[c_idx] = f"{b_id}_{row[c_idx]}"
                        cdb_writer.writerow(row)

            # Merge Widb
            widb_path = os.path.join(b_dir, "data_tables", "Widb.csv")
            if not os.path.isfile(widb_path):
                widb_path = os.path.join(b_dir, "data_tables", "Wdb.csv")
            if os.path.isfile(widb_path) and os.path.getsize(widb_path) > 0:
                with open(widb_path, "r", encoding="utf-8", errors="ignore") as in_f:
                    reader = csv.reader(in_f)
                    h = next(reader, None)
                    if h and not widb_header_written:
                        widb_writer.writerow(h)
                        widb_header_written = True
                    cluster_cols = []
                    if h:
                        for idx, col in enumerate(h):
                            if "cluster" in col.lower():
                                cluster_cols.append(idx)
                    for row in reader:
                        if not row:
                            continue
                        for c_idx in cluster_cols:
                            if len(row) > c_idx and row[c_idx]:
                                row[c_idx] = f"{b_id}_{row[c_idx]}"
                        widb_writer.writerow(row)

            # Merge genomeInfo
            ginfo_path = os.path.join(b_dir, "genomeInfo.csv")
            if os.path.isfile(ginfo_path) and os.path.getsize(ginfo_path) > 0:
                with open(ginfo_path, "r", encoding="utf-8", errors="ignore") as in_f:
                    reader = csv.reader(in_f)
                    h = next(reader, None)
                    if h and not ginfo_header_written:
                        ginfo_writer.writerow(h)
                        ginfo_header_written = True
                    for row in reader:
                        if row:
                            ginfo_writer.writerow(row)

    print(
        f"[union-blocks] Merged {len(block_ids)} disjoint blocks: {total_genomes} winner genomes into {out_derep}",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Helper utility script for dRep operations in mag_traits module."
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # Subcommand: pool-manifest
    p_pool = subparsers.add_parser(
        "pool-manifest",
        help="Scan dereplication output directories and build a pooled manifest TSV.",
    )
    p_pool.add_argument(
        "--out-manifest",
        required=True,
        help="Path to output manifest TSV (mag_id \\t path)",
    )
    p_pool.add_argument(
        "--block-dirs",
        nargs="+",
        required=True,
        help="List of dereplicated_genomes directories from blocks",
    )
    p_pool.set_defaults(func=cmd_pool_manifest)

    # Subcommand: merge-genome-info
    p_merge = subparsers.add_parser(
        "merge-genome-info",
        help="Merge multiple genomeInfo CSV files into one combined file.",
    )
    p_merge.add_argument(
        "--out",
        required=True,
        help="Destination path for merged genomeInfo CSV",
    )
    p_merge.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="List of source genomeInfo CSV files",
    )
    p_merge.set_defaults(func=cmd_merge_genome_info)

    # Subcommand: filter-genome-info
    p_filter = subparsers.add_parser(
        "filter-genome-info",
        help="Filter genomeInfo CSV to keep only dereplicated winner genomes.",
    )
    p_filter.add_argument(
        "--derep-dir",
        required=True,
        help="Path to dereplicated_genomes directory containing winner FASTA files",
    )
    p_filter.add_argument(
        "--src",
        required=True,
        help="Path to source genomeInfo CSV",
    )
    p_filter.add_argument(
        "--dst",
        required=True,
        help="Path to destination filtered genomeInfo CSV",
    )
    p_filter.set_defaults(func=cmd_filter_genome_info)

    # Subcommand: union-blocks
    p_union = subparsers.add_parser(
        "union-blocks",
        help="Directly union dereplicated genomes and metadata tables from disjoint taxonomy blocks.",
    )
    p_union.add_argument(
        "--block-list",
        required=True,
        help="Path to block_list.txt",
    )
    p_union.add_argument(
        "--blocks-dir",
        required=True,
        help="Parent directory containing block output subdirectories",
    )
    p_union.add_argument(
        "--out-derep",
        required=True,
        help="Destination directory for unioned dereplicated FASTA files",
    )
    p_union.add_argument(
        "--out-cdb",
        required=True,
        help="Destination path for combined Cdb.csv",
    )
    p_union.add_argument(
        "--out-widb",
        required=True,
        help="Destination path for combined Widb.csv",
    )
    p_union.add_argument(
        "--out-ginfo",
        required=True,
        help="Destination path for combined genomeInfo.csv",
    )
    p_union.set_defaults(func=cmd_union_blocks)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

