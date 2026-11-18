#!/usr/bin/env python3
"""
scripts/prepare_drep_inputs.py

High-performance input preparation script for dRep.
Handles 70,000+ genomes in seconds without incurring bash fork overhead
or hitting Linux command line length limits (ARG_MAX / E2BIG).

Key Features:
1. Fast staging: Creates symlinks for uncompressed FASTA files and streams
   decompression for .gz files into a staging directory.
2. Formats all staged genome basenames uniformly as `<mag_id>.fa`.
3. Writes `genomes_list.txt` containing absolute paths, suitable for `dRep -g <file>`.
4. Intelligently parses and aligns `genomeInfo.csv` to ensure 100% of staged
   genomes have corresponding completeness & contamination entries (preventing
   dRep runtime fatal crashes when quality tables miss genomes).
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare input genomes and genomeInfo table for dRep at scale."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to raw_mags_manifest.tsv (tab-separated: mag_id \\t path)",
    )
    parser.add_argument(
        "--stage-dir",
        required=True,
        help="Directory to stage symlinked or decompressed .fa files",
    )
    parser.add_argument(
        "--out-list",
        required=True,
        help="Path to output text file listing staged genome paths (one per line)",
    )
    parser.add_argument(
        "--genome-info-src",
        default="",
        help="Path to user-provided genomeInfo file (CSV or TSV, optional)",
    )
    parser.add_argument(
        "--genome-info-dst",
        default="",
        help="Path to output standardized genomeInfo.csv (optional)",
    )
    return parser.parse_args()


def stage_genomes(manifest_path, stage_dir, out_list_path):
    os.makedirs(stage_dir, exist_ok=True)
    out_list_dir = os.path.dirname(os.path.abspath(out_list_path))
    if out_list_dir:
        os.makedirs(out_list_dir, exist_ok=True)

    staged_paths = []
    staged_mags = set()
    skipped_count = 0

    with open(manifest_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            mag_id, src_path = parts[0].strip(), parts[1].strip()
            if not mag_id or not src_path:
                continue

            target_fa = os.path.join(stage_dir, f"{mag_id}.fa")
            src_abs = os.path.abspath(src_path)

            if src_path.endswith(".gz"):
                if not os.path.exists(src_path):
                    print(
                        f"WARNING: File not found for MAG '{mag_id}': {src_path}",
                        file=sys.stderr,
                    )
                    skipped_count += 1
                    continue
                try:
                    if os.path.islink(target_fa) or os.path.exists(target_fa):
                        os.unlink(target_fa)
                except OSError:
                    pass
                with gzip.open(src_path, "rb") as f_in, open(target_fa, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out, length=1024 * 1024)
            else:
                try:
                    os.symlink(src_abs, target_fa)
                except FileExistsError:
                    try:
                        os.unlink(target_fa)
                        os.symlink(src_abs, target_fa)
                    except OSError:
                        pass
                except OSError:
                    if os.path.exists(src_path):
                        shutil.copyfile(src_path, target_fa)
                    else:
                        print(
                            f"WARNING: File not found for MAG '{mag_id}': {src_path}",
                            file=sys.stderr,
                        )
                        skipped_count += 1
                        continue

            abs_target = os.path.abspath(target_fa)
            staged_paths.append(abs_target)
            staged_mags.add(mag_id)

    if not staged_paths:
        sys.exit(
            f"ERROR: No valid genomes could be staged from manifest: {manifest_path}"
        )

    with open(out_list_path, "w", encoding="utf-8") as f_out:
        for p in staged_paths:
            f_out.write(p + "\n")

    print(
        f"INFO: Successfully staged {len(staged_paths)} genomes in '{stage_dir}' (skipped: {skipped_count}).",
        file=sys.stderr,
    )
    print(f"INFO: Genomes list written to '{out_list_path}'.", file=sys.stderr)
    return staged_mags


def process_genome_info(src_path, dst_path, staged_mags):
    if not dst_path:
        return

    dst_dir = os.path.dirname(os.path.abspath(dst_path))
    if dst_dir:
        os.makedirs(dst_dir, exist_ok=True)

    if not src_path or not os.path.isfile(src_path) or os.path.getsize(src_path) == 0:
        # Create empty destination file
        open(dst_path, "w", encoding="utf-8").close()
        print(
            f"INFO: No valid source genomeInfo provided. Wrote empty '{dst_path}'.",
            file=sys.stderr,
        )
        return

    # Detect delimiter
    with open(src_path, "r", encoding="utf-8", errors="ignore") as f:
        sample = f.read(4096)
        delim = "\t" if "\t" in sample and "," not in sample else ","

    quality_map = {}
    with open(src_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader, None)
        if not header:
            open(dst_path, "w", encoding="utf-8").close()
            return

        header_lower = [h.strip().lower() for h in header]
        g_col = next(
            (
                i
                for i, h in enumerate(header_lower)
                if any(k in h for k in ["genome", "name", "id", "mag", "bin"])
            ),
            0,
        )
        c_col = next((i for i, h in enumerate(header_lower) if "comp" in h), 1)
        t_col = next((i for i, h in enumerate(header_lower) if "cont" in h), 2)

        # Detect if header is actually a data row (all digits or missing keywords)
        first_row_is_header = any(
            re.search(r"[a-zA-Z]", str(header[c]))
            for c in [c_col, t_col]
            if c < len(header)
        )

        def record_row(r):
            if len(r) > max(g_col, c_col, t_col):
                g_raw = r[g_col].strip()
                g_stem = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", g_raw)
                comp = r[c_col].strip()
                cont = r[t_col].strip()
                quality_map[g_raw] = (comp, cont)
                quality_map[g_stem] = (comp, cont)

        if not first_row_is_header:
            record_row(header)

        for r in reader:
            if r:
                record_row(r)

    rows = []
    missing_count = 0
    # Sort for reproducibility
    for mag_id in sorted(staged_mags):
        fa_name = f"{mag_id}.fa"
        if mag_id in quality_map:
            comp, cont = quality_map[mag_id]
        elif fa_name in quality_map:
            comp, cont = quality_map[fa_name]
        else:
            missing_count += 1
            # Fallback quality values to prevent dRep from failing
            comp, cont = "0", "100"
        rows.append([fa_name, comp, cont])

    if missing_count > 0:
        print(
            f"WARNING: {missing_count}/{len(staged_mags)} staged genomes were missing from '{src_path}'. "
            f"Filled default quality (completeness=0, contamination=100) to ensure dRep stability.",
            file=sys.stderr,
        )

    with open(dst_path, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["genome", "completeness", "contamination"])
        writer.writerows(rows)

    print(
        f"INFO: Aligned quality table with {len(rows)} entries written to '{dst_path}'.",
        file=sys.stderr,
    )


def main():
    args = parse_args()
    staged_mags = stage_genomes(args.manifest, args.stage_dir, args.out_list)
    if args.genome_info_dst:
        process_genome_info(args.genome_info_src, args.genome_info_dst, staged_mags)


if __name__ == "__main__":
    main()
