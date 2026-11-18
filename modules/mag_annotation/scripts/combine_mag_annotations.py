#!/usr/bin/env python3
"""
combine_mag_annotations.py
Combines two or more non-overlapping mag_annotation result directories into a unified result directory.

Key Operations:
1. Merges primary and secondary tabular assets in `merged/`:
   - `protein_KO_all_strict.tsv` (streaming merge, validates column compatibility)
   - `protein_KO_best.tsv`
   - `genome_KO_counts.long.tsv`
   - Optional: regenerates dense matrices if present in inputs
2. Zero-copy symlinking for upstream per-genome files:
   - `proteins/all/{genome}.faa.gz`
   - `genes/{genome}.gff.gz`
   - `metadata/{genome}_protein_metadata.tsv.gz`
3. Verifies zero genome-ID collisions across batches.
"""

import argparse
import gzip
import os
import shutil
import sys
import time

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def find_asset(dir_path, sub_dir, filename):
    """Finds uncompressed or gzipped version of an asset."""
    p1 = os.path.join(dir_path, sub_dir, filename)
    if os.path.exists(p1):
        return p1
    p2 = p1 + ".gz"
    if os.path.exists(p2):
        return p2
    return None

def merge_tsv_files(input_files, output_file, name_desc="TSV"):
    """
    Streams and concatenates multiple TSV files, preserving the header from the first file.
    Validates header consistency.
    """
    if not input_files:
        return 0

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    first_header = None
    total_records = 0
    buffer = []
    BUFFER_SIZE = 100000

    out_mode = "wt" if output_file.endswith(".gz") else "w"
    open_out = gzip.open if output_file.endswith(".gz") else open

    with open_out(output_file, out_mode, encoding="utf-8") as out_f:
        for idx, in_file in enumerate(input_files):
            if not os.path.exists(in_file):
                print(f"Warning: {in_file} not found, skipping.", file=sys.stderr)
                continue

            with open_file(in_file, "rt") as in_f:
                header = in_f.readline()
                if not header:
                    continue
                header_strip = header.strip()
                if first_header is None:
                    first_header = header_strip
                    out_f.write(header)
                else:
                    if header_strip != first_header:
                        print(f"Warning: Header mismatch in {in_file} for {name_desc}!", file=sys.stderr)

                for line in in_f:
                    if not line.strip():
                        continue
                    buffer.append(line)
                    total_records += 1
                    if len(buffer) >= BUFFER_SIZE:
                        out_f.write("".join(buffer))
                        buffer.clear()

        if buffer:
            out_f.write("".join(buffer))
            buffer.clear()

    return total_records

def link_genome_files(src_dir, dst_dir, file_exts, desc="files"):
    """
    Creates relative/absolute symlinks for per-genome files from src_dir into dst_dir.
    Detects and reports genome collisions.
    """
    if not os.path.exists(src_dir):
        return 0

    os.makedirs(dst_dir, exist_ok=True)
    linked_count = 0

    for root, _, files in os.walk(src_dir):
        for f in files:
            for ext in file_exts:
                if f.endswith(ext):
                    src_f = os.path.join(root, f)
                    dst_f = os.path.join(dst_dir, f)
                    if os.path.exists(dst_f) or os.path.islink(dst_f):
                        if os.path.islink(dst_f) and os.path.realpath(dst_f) == os.path.realpath(src_f):
                            pass # already linked
                        else:
                            raise ValueError(f"Genome collision detected! File '{f}' already exists in destination {dst_dir} from another batch.")
                    else:
                        os.symlink(os.path.abspath(src_f), dst_f)
                    linked_count += 1
                    break

    return linked_count

def main():
    parser = argparse.ArgumentParser(description="Combine multiple non-overlapping mag_annotation result directories.")
    parser.add_argument("-i", "--inputs", nargs="+", required=True,
                        help="List of input mag_annotation result directories to combine")
    parser.add_argument("-o", "--out_dir", required=True,
                        help="Destination combined mag_annotation result directory")
    parser.add_argument("--link_upstream", action="store_true", default=True,
                        help="Create fast symlinks for proteins/all, genes, metadata (default: True)")
    parser.add_argument("--no_link_upstream", dest="link_upstream", action="store_false",
                        help="Do not create symlinks for upstream per-genome files")

    args = parser.parse_args()

    input_dirs = [os.path.abspath(d) for d in args.inputs]
    out_dir = os.path.abspath(args.out_dir)

    print("==============================================================================")
    print("MAG Annotation Batch Combiner")
    print("==============================================================================")
    print(f"Input directories ({len(input_dirs)}):")
    for d in input_dirs:
        if not os.path.exists(d):
            print(f"ERROR: Input directory does not exist: {d}", file=sys.stderr)
            sys.exit(1)
        print(f"  - {d}")
    print(f"Destination directory: {out_dir}")
    print("==============================================================================")

    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    # 1. Merge Primary & Secondary Derivative Tables in merged/
    print("\n[Step 1/3] Merging canonical tables in 'merged/'...")
    merged_out = os.path.join(out_dir, "merged")
    os.makedirs(merged_out, exist_ok=True)

    assets_to_merge = [
        ("protein_KO_all_strict.tsv", "Strict multi-hit annotations"),
        ("protein_KO_best.tsv", "Single best KO annotations"),
        ("genome_KO_counts.long.tsv", "Genome KO count sparse long table"),
    ]

    has_dense = False
    for fname, desc in assets_to_merge:
        in_files = []
        is_gz = False
        for d in input_dirs:
            found = find_asset(d, "merged", fname)
            if found:
                in_files.append(found)
                if found.endswith(".gz"):
                    is_gz = True

        out_fname = fname + (".gz" if is_gz else "")
        target_path = os.path.join(merged_out, out_fname)
        
        print(f"  -> Merging {desc} ({fname})...")
        cnt = merge_tsv_files(in_files, target_path, name_desc=desc)
        print(f"     Merged {cnt:,} rows -> {target_path}")

    # Check for dense matrices in inputs
    dense_count_files = [find_asset(d, "merged", "genome_KO_counts.dense.tsv") for d in input_dirs]
    dense_count_files = [f for f in dense_count_files if f]
    if len(dense_count_files) == len(input_dirs) and dense_count_files:
        print("  -> Regenerating combined dense genome x KO count and presence matrices from long table...")
        long_combined = os.path.join(merged_out, "genome_KO_counts.long.tsv")
        if not os.path.exists(long_combined) and os.path.exists(long_combined + ".gz"):
            long_combined = long_combined + ".gz"
        
        dense_c_out = os.path.join(merged_out, "genome_KO_counts.dense.tsv")
        dense_p_out = os.path.join(merged_out, "genome_KO_presence.dense.tsv")
        
        genomes_set = set()
        kos_set = set()
        genome_ko_counts = {}

        with open_file(long_combined, "rt") as lf:
            _ = lf.readline()
            for line in lf:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    g, k, c = parts[0], parts[1], parts[2]
                    genomes_set.add(g)
                    kos_set.add(k)
                    genome_ko_counts[(g, k)] = int(c)

        sorted_kos = sorted(kos_set)
        sorted_genomes = sorted(genomes_set)

        with open(dense_c_out, "w", encoding="utf-8") as c_f, \
             open(dense_p_out, "w", encoding="utf-8") as p_f:
            c_f.write("genome_id\t" + "\t".join(sorted_kos) + "\n")
            p_f.write("genome_id\t" + "\t".join(sorted_kos) + "\n")
            for g in sorted_genomes:
                cnts = [str(genome_ko_counts.get((g, k), 0)) for k in sorted_kos]
                pres = ["1" if genome_ko_counts.get((g, k), 0) > 0 else "0" for k in sorted_kos]
                c_f.write(f"{g}\t" + "\t".join(cnts) + "\n")
                p_f.write(f"{g}\t" + "\t".join(pres) + "\n")
        print(f"     Exported dense matrices across {len(sorted_genomes):,} genomes and {len(sorted_kos):,} KOs.")

    # 2. Merge unique mapping table in proteins/unique/
    print("\n[Step 2/3] Merging mapping tables in 'proteins/unique/'...")
    map_files = [find_asset(d, os.path.join("proteins", "unique"), "protein_to_unique.tsv.gz") for d in input_dirs]
    map_files = [f for f in map_files if f]
    if map_files:
        unique_out = os.path.join(out_dir, "proteins", "unique")
        os.makedirs(unique_out, exist_ok=True)
        target_map = os.path.join(unique_out, "protein_to_unique.tsv.gz")
        cnt = merge_tsv_files(map_files, target_map, name_desc="protein_to_unique mapping")
        print(f"  -> Combined protein mapping table ({cnt:,} entries) -> {target_map}")

    # 3. Create Zero-Copy Symlinks for per-genome files
    if args.link_upstream:
        print("\n[Step 3/3] Creating zero-copy symlinks for per-genome assets...")
        per_genome_tasks = [
            (os.path.join("proteins", "all"), [".faa.gz", ".faa"], "Proteins (.faa.gz)"),
            ("genes", [".gff.gz", ".gff"], "Genes (.gff.gz)"),
            ("metadata", ["_protein_metadata.tsv.gz", "_protein_metadata.tsv"], "Protein metadata (.tsv.gz)"),
        ]

        for subpath, exts, desc in per_genome_tasks:
            target_sub = os.path.join(out_dir, subpath)
            total_linked = 0
            for d in input_dirs:
                src_sub = os.path.join(d, subpath)
                n = link_genome_files(src_sub, target_sub, exts, desc=desc)
                total_linked += n
            print(f"  -> Symlinked {total_linked:,} {desc} -> {target_sub}")
    else:
        print("\n[Step 3/3] Skipping upstream per-genome symlinks (--no_link_upstream).")

    elapsed = time.time() - t0
    print("\n==============================================================================")
    print(f"✅ Successfully combined {len(input_dirs)} datasets into:")
    print(f"   {out_dir}")
    print(f"Elapsed time: {elapsed:.2f} seconds.")
    print("==============================================================================")

if __name__ == "__main__":
    main()
