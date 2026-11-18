#!/usr/bin/env python3
"""
scripts/pool_mag_genes_and_proteins.py

Gathers predicted genes (.ffn, .gff) and proteins (.faa) from partitioned
batch subdirectories into consolidated top-level directories:
- proteins/all/ : all predicted proteins (<mag>.faa)
- genes/        : all predicted nucleotide genes (<mag>.ffn) and annotations (<mag>.gff)

Uses zero-copy absolute symlinks by default for instant execution and zero disk duplication.
"""

import argparse
import os
import shutil
import sys
import time


def parse_args():
    parser = argparse.ArgumentParser(
        description="Pool predicted gene and protein files across batch folders."
    )
    parser.add_argument(
        "--batches_dir",
        required=True,
        help="Path to batches directory (e.g. results/mag_traits/batches)",
    )
    parser.add_argument(
        "--proteins_dir",
        required=True,
        help="Destination directory for proteins (e.g. results/mag_traits/proteins/all)",
    )
    parser.add_argument(
        "--genes_dir",
        required=True,
        help="Destination directory for genes (e.g. results/mag_traits/genes)",
    )
    parser.add_argument(
        "--batch_list",
        default="",
        help="Optional path to batch_list.txt. If omitted, scans batches_dir for batch_* folders.",
    )
    parser.add_argument(
        "--done_file",
        default="",
        help="Optional flag file to touch upon completion",
    )
    parser.add_argument(
        "--copy",
        action="store_true",
        default=False,
        help="Copy files instead of symlinking (default: False)",
    )
    return parser.parse_args()


def link_or_copy(src, dst, use_copy=False):
    src_abs = os.path.realpath(src)
    if os.path.islink(dst):
        try:
            if os.path.realpath(dst) == src_abs:
                return  # already correctly linked
            os.unlink(dst)
        except OSError:
            pass
    elif os.path.exists(dst):
        if use_copy:
            return  # already exists
        try:
            os.unlink(dst)
        except OSError:
            pass

    if use_copy:
        shutil.copyfile(src_abs, dst)
    else:
        try:
            os.symlink(src_abs, dst)
        except OSError:
            # Fallback to copy if symlinks not supported
            shutil.copyfile(src_abs, dst)


def main():
    args = parse_args()
    start_time = time.time()

    batches_dir = os.path.abspath(args.batches_dir)
    proteins_dir = os.path.abspath(args.proteins_dir)
    genes_dir = os.path.abspath(args.genes_dir)

    os.makedirs(proteins_dir, exist_ok=True)
    os.makedirs(genes_dir, exist_ok=True)

    # Determine batches to process
    batch_dirs = []
    if args.batch_list and os.path.isfile(args.batch_list):
        with open(args.batch_list, "r", encoding="utf-8") as f:
            for line in f:
                b = line.strip()
                if b:
                    b_path = os.path.join(batches_dir, b)
                    if os.path.isdir(b_path):
                        batch_dirs.append(b_path)
    else:
        if os.path.isdir(batches_dir):
            for name in sorted(os.listdir(batches_dir)):
                b_path = os.path.join(batches_dir, name)
                if os.path.isdir(b_path) and (name.startswith("batch_") or os.path.isdir(os.path.join(b_path, "faa"))):
                    batch_dirs.append(b_path)

    print(f"INFO: Found {len(batch_dirs)} batch directories to pool.", file=sys.stderr)

    total_faa = 0
    total_ffn = 0
    total_gff = 0

    for b_path in batch_dirs:
        # 1. Pool FAA -> proteins_dir
        faa_dir = os.path.join(b_path, "faa")
        if os.path.isdir(faa_dir):
            for fname in os.listdir(faa_dir):
                if fname.endswith((".faa", ".faa.gz")):
                    src = os.path.join(faa_dir, fname)
                    dst = os.path.join(proteins_dir, fname)
                    link_or_copy(src, dst, use_copy=args.copy)
                    total_faa += 1

        # 2. Pool FFN -> genes_dir
        ffn_dir = os.path.join(b_path, "ffn")
        if os.path.isdir(ffn_dir):
            for fname in os.listdir(ffn_dir):
                if fname.endswith((".ffn", ".ffn.gz")):
                    src = os.path.join(ffn_dir, fname)
                    dst = os.path.join(genes_dir, fname)
                    link_or_copy(src, dst, use_copy=args.copy)
                    total_ffn += 1

        # 3. Pool GFF -> genes_dir
        gff_dir = os.path.join(b_path, "gff")
        if os.path.isdir(gff_dir):
            for fname in os.listdir(gff_dir):
                if fname.endswith((".gff", ".gff.gz")):
                    src = os.path.join(gff_dir, fname)
                    dst = os.path.join(genes_dir, fname)
                    link_or_copy(src, dst, use_copy=args.copy)
                    total_gff += 1

    elapsed = time.time() - start_time
    summary_msg = (
        f"INFO: Successfully pooled across {len(batch_dirs)} batches in {elapsed:.2f}s:\n"
        f"  - Proteins (.faa): {total_faa} -> {proteins_dir}\n"
        f"  - Genes (.ffn):    {total_ffn} -> {genes_dir}\n"
        f"  - Genes (.gff):    {total_gff} -> {genes_dir}\n"
    )
    print(summary_msg, file=sys.stderr)

    if args.done_file:
        done_dir = os.path.dirname(os.path.abspath(args.done_file))
        if done_dir:
            os.makedirs(done_dir, exist_ok=True)
        with open(args.done_file, "w", encoding="utf-8") as f:
            f.write(summary_msg)


if __name__ == "__main__":
    main()
