#!/usr/bin/env python3
"""
Convert sample TSV table (RUN_id, Store_address) into standard Genome_anno samples.yaml.
Usage:
    python3 scripts/tsv_to_samples_yaml.py \
        --tsv samples.tsv \
        --project-root /dssg/home/acct-schhui/gses0/liliuyang/osma/test \
        --output projects/config_read/samples.yaml
"""

import argparse
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert TSV sample sheet to Genome_anno samples.yaml format."
    )
    parser.add_argument(
        "--tsv",
        "-i",
        required=True,
        help="Input TSV file path containing columns RUN_id and Store_address",
    )
    parser.add_argument(
        "--project-root",
        "-p",
        default="/dssg/home/acct-schhui/gses0/liliuyang/osma/test",
        help="Project root directory path",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="projects/config_read/samples.yaml",
        help="Target output samples.yaml file path",
    )
    parser.add_argument(
        "--subfolder",
        default="00.QC",
        help="QC subdirectory name under Store_address (default: 00.QC)",
    )
    parser.add_argument(
        "--prefix-replace",
        default="",
        help="Optional prefix mapping, e.g. '/dssg/home:/union/home'",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if not os.path.exists(args.tsv):
        print(f"Error: Input TSV file not found: {args.tsv}", file=sys.stderr)
        sys.exit(1)

    with open(args.tsv, "r", encoding="utf-8") as f:
        lines = [line.strip().split("\t") for line in f if line.strip() and not line.startswith("#")]

    if not lines:
        print(f"Error: Input TSV file is empty: {args.tsv}", file=sys.stderr)
        sys.exit(1)

    header = [h.strip() for h in lines[0]]
    run_col = 0
    store_col = 1
    has_header = False

    for idx, col in enumerate(header):
        col_lower = col.lower()
        if col_lower in ("run_id", "sample", "sample_id", "id", "sampleid"):
            run_col = idx
            has_header = True
        elif col_lower in ("store_address", "store_addr", "address", "path", "dir", "directory"):
            store_col = idx
            has_header = True

    data_rows = lines[1:] if has_header else lines
    sample_records = []

    src_prefix = ""
    dst_prefix = ""
    if args.prefix_replace and ":" in args.prefix_replace:
        src_prefix, dst_prefix = args.prefix_replace.split(":", 1)

    for row in data_rows:
        if len(row) > max(run_col, store_col):
            run_id = row[run_col].strip()
            store_addr = row[store_col].strip()
            if run_id and store_addr:
                if src_prefix and src_prefix in store_addr:
                    store_addr = store_addr.replace(src_prefix, dst_prefix, 1)
                r1_path = os.path.join(store_addr, args.subfolder, f"{run_id}.clean.1.fq.gz")
                r2_path = os.path.join(store_addr, args.subfolder, f"{run_id}.clean.2.fq.gz")

                # Auto-fallback check if file doesn't exist
                if not os.path.exists(r1_path):
                    if "/dssg/home/" in r1_path:
                        alt_r1 = r1_path.replace("/dssg/home/", "/union/home/", 1)
                        if os.path.exists(alt_r1):
                            r1_path = alt_r1
                            r2_path = r2_path.replace("/dssg/home/", "/union/home/", 1)
                    elif "/union/home/" in r1_path:
                        alt_r1 = r1_path.replace("/union/home/", "/dssg/home/", 1)
                        if os.path.exists(alt_r1):
                            r1_path = alt_r1
                            r2_path = r2_path.replace("/union/home/", "/dssg/home/", 1)

                sample_records.append((run_id, r1_path, r2_path))

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    with open(args.output, "w", encoding="utf-8") as out:
        out.write("# ==============================================================================\n")
        out.write("# Samples Configuration (config_read)\n")
        out.write("# Auto-generated from TSV input via scripts/tsv_to_samples_yaml.py\n")
        out.write("# ==============================================================================\n\n")
        out.write(f'project_root: "{args.project_root}"\n\n')
        out.write("data:\n")
        out.write('  samples: "./samples.yaml"\n\n')
        out.write("samples:\n")
        for sample_id, r1, r2 in sample_records:
            out.write(f'  "{sample_id}":\n')
            out.write("    short_read:\n")
            out.write(f'      r1: "{r1}"\n')
            out.write(f'      r2: "{r2}"\n\n')

    print(f"Successfully generated {args.output} with {len(sample_records)} samples.")


if __name__ == "__main__":
    main()
