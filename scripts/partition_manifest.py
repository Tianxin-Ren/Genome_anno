#!/usr/bin/env python3
"""
scripts/partition_manifest.py

Partition a large MAG manifest (e.g. 70,000 genomes) into multiple smaller
block manifests for parallel dRep dereplication:
1. Taxonomy-guided partition (recommended): Groups MAGs by GTDB-Tk Phylum (p__)
   and Class (c__) so each block is a biologically disjoint taxon. Disjoint taxa
   can be dereplicated independently in parallel without requiring an explosive
   all-vs-all pool consolidation.
2. Random partition (fallback): Evenly splits MAGs into fixed-size chunks (e.g. 10,000 each).
"""

import argparse
from collections import defaultdict
import csv
import os
import random
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Partition a MAG manifest into block manifests for parallel dRep."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to input raw_mags_manifest.tsv (format: mag_id <TAB> fasta_path).",
    )
    parser.add_argument(
        "--taxonomy",
        default="",
        help="Optional path to GTDB-Tk classification file (TSV/CSV with 'bin' and 'classification' columns). "
             "If provided and non-empty, partitions by Phylum/Class instead of random blocks.",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=10000,
        help="Maximum number of genomes per block (default: 10000).",
    )
    parser.add_argument(
        "--min-phylum-size",
        type=int,
        default=500,
        help="Minimum number of genomes for a phylum to form its own block (default: 500). "
             "Phyla with fewer genomes are aggregated into combined 'tax_minor_phyla' blocks "
             "to prevent small taxa from fragmenting jobs and consuming heavy cluster resources.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory to store block TSV files.",
    )
    parser.add_argument(
        "--out-list",
        required=True,
        help="Output text file listing generated block IDs (one per line).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for deterministic shuffling (default: 42).",
    )
    parser.add_argument(
        "--shuffle",
        action="store_true",
        default=True,
        help="Deterministically shuffle MAGs within random blocks (default: True).",
    )
    parser.add_argument(
        "--no-shuffle",
        dest="shuffle",
        action="store_false",
        help="Preserve manifest order without shuffling.",
    )
    return parser.parse_args()


def load_taxonomy(taxonomy_path):
    """
    Parses GTDB-Tk TSV/CSV and returns a dict mapping bin_id / stem -> (phylum, class).
    """
    if not taxonomy_path or not os.path.isfile(taxonomy_path) or os.path.getsize(taxonomy_path) == 0:
        return {}

    with open(taxonomy_path, "r", encoding="utf-8", errors="ignore") as f:
        sample = f.read(4096)
        delim = "\t" if "\t" in sample and "," not in sample else ","

    tax_map = {}
    with open(taxonomy_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader, None)
        if not header:
            return {}

        header_lower = [h.strip().lower() for h in header]
        bin_col = next(
            (i for i, h in enumerate(header_lower) if any(k in h for k in ["bin", "genome", "name", "id", "mag"])),
            0
        )
        tax_col = next(
            (i for i, h in enumerate(header_lower) if any(k in h for k in ["classification", "tax", "lineage"])),
            len(header) - 1
        )

        for row in reader:
            if not row or len(row) <= max(bin_col, tax_col):
                continue
            bin_id = row[bin_col].strip()
            tax_str = row[tax_col].strip()
            if not bin_id:
                continue

            p_m = re.search(r"p__([^;]+)", tax_str)
            c_m = re.search(r"c__([^;]+)", tax_str)
            phylum = p_m.group(1).strip() if p_m else ""
            cls = c_m.group(1).strip() if c_m else ""
            if phylum.lower() in ("unclassified", "unknown", "none", "na", ""):
                phylum = ""
            if cls.lower() in ("unclassified", "unknown", "none", "na", ""):
                cls = ""

            stem = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", bin_id)
            if phylum:
                tax_map[bin_id] = (phylum, cls)
                tax_map[stem] = (phylum, cls)

    return tax_map


def partition_by_taxonomy(entries, tax_map, block_size, min_phylum_size, out_dir, out_list):
    """
    Partitions manifest entries by GTDB Phylum and Class.
    If a phylum has >= min_phylum_size genomes:
      - If <= block_size, it forms a single block (tax_p__<Phylum>).
      - If > block_size, it is subdivided by class.
    Phyla with < min_phylum_size genomes are aggregated into combined 'tax_minor_phyla' block(s).
    """
    phylum_groups = defaultdict(lambda: defaultdict(list))
    unclassified_list = []

    for line in entries:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        mag_id = parts[0].strip()
        stem = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", mag_id)

        if mag_id in tax_map:
            p, c = tax_map[mag_id]
            phylum_groups[p][c].append(line)
        elif stem in tax_map:
            p, c = tax_map[stem]
            phylum_groups[p][c].append(line)
        else:
            unclassified_list.append(line)

    block_ids = []
    minor_phyla_lines = []

    # Separate major phyla from minor phyla
    major_phyla = {}
    for phylum, cls_dict in phylum_groups.items():
        total_p_lines = sum(len(lines) for lines in cls_dict.values())
        if min_phylum_size > 0 and total_p_lines < min_phylum_size:
            for c_lines in cls_dict.values():
                minor_phyla_lines.extend(c_lines)
        else:
            major_phyla[phylum] = cls_dict

    # Process major phyla
    for phylum in sorted(major_phyla.keys()):
        cls_dict = major_phyla[phylum]
        total_p_lines = sum(len(lines) for lines in cls_dict.values())
        safe_p = re.sub(r"[^a-zA-Z0-9_-]", "_", phylum)

        if total_p_lines <= block_size:
            # Whole phylum fits in one block
            bid = f"tax_p__{safe_p}"
            block_ids.append(bid)
            block_file = os.path.join(out_dir, f"{bid}.tsv")
            with open(block_file, "w", encoding="utf-8") as f:
                for c_lines in cls_dict.values():
                    f.writelines(c_lines)
            print(f"[partition_manifest]   -> {bid}: {total_p_lines} genomes ({block_file})", file=sys.stderr)
        else:
            # Phylum exceeds block_size, subdivide by class
            major_classes = {}
            minor_class_lines = []

            for cls, c_lines in cls_dict.items():
                if min_phylum_size > 0 and len(c_lines) < min_phylum_size:
                    minor_class_lines.extend(c_lines)
                else:
                    major_classes[cls] = c_lines

            # Process major classes (each has >= min_phylum_size genomes)
            for cls in sorted(major_classes.keys()):
                c_lines = major_classes[cls]
                safe_c = re.sub(r"[^a-zA-Z0-9_-]", "_", cls) if cls else "unclassified_class"
                if len(c_lines) <= block_size:
                    bid = f"tax_p__{safe_p}_c__{safe_c}"
                    block_ids.append(bid)
                    block_file = os.path.join(out_dir, f"{bid}.tsv")
                    with open(block_file, "w", encoding="utf-8") as f:
                        f.writelines(c_lines)
                    print(f"[partition_manifest]   -> {bid}: {len(c_lines)} genomes ({block_file})", file=sys.stderr)
                else:
                    # Class exceeds block_size, chunk into multiple parts
                    num_sub = (len(c_lines) + block_size - 1) // block_size
                    for part_idx in range(num_sub):
                        bid = f"tax_p__{safe_p}_c__{safe_c}_part{part_idx:03d}"
                        chunk = c_lines[part_idx * block_size : (part_idx + 1) * block_size]
                        block_ids.append(bid)
                        block_file = os.path.join(out_dir, f"{bid}.tsv")
                        with open(block_file, "w", encoding="utf-8") as f:
                            f.writelines(chunk)
                        print(f"[partition_manifest]   -> {bid}: {len(chunk)} genomes ({block_file})", file=sys.stderr)

            # Process minor classes within this phylum
            if minor_class_lines:
                if len(minor_class_lines) >= min_phylum_size:
                    bid = f"tax_p__{safe_p}_c__other"
                    block_ids.append(bid)
                    block_file = os.path.join(out_dir, f"{bid}.tsv")
                    with open(block_file, "w", encoding="utf-8") as f:
                        f.writelines(minor_class_lines)
                    print(f"[partition_manifest]   -> {bid}: {len(minor_class_lines)} genomes from small classes in {phylum} ({block_file})", file=sys.stderr)
                else:
                    minor_phyla_lines.extend(minor_class_lines)
                    print(f"[partition_manifest] Aggregated {len(minor_class_lines)} genomes from tiny classes in {phylum} into tax_minor_phyla.", file=sys.stderr)

    # Process aggregated minor phyla (< min_phylum_size)
    if minor_phyla_lines:
        print(f"[partition_manifest] Aggregated {len(minor_phyla_lines)} genomes from small phyla (< {min_phylum_size} genomes) into combined block(s).", file=sys.stderr)
        if len(minor_phyla_lines) <= block_size:
            bid = "tax_minor_phyla"
            block_ids.append(bid)
            block_file = os.path.join(out_dir, f"{bid}.tsv")
            with open(block_file, "w", encoding="utf-8") as f:
                f.writelines(minor_phyla_lines)
            print(f"[partition_manifest]   -> {bid}: {len(minor_phyla_lines)} genomes ({block_file})", file=sys.stderr)
        else:
            num_sub = (len(minor_phyla_lines) + block_size - 1) // block_size
            for part_idx in range(num_sub):
                bid = f"tax_minor_phyla_part{part_idx:03d}"
                chunk = minor_phyla_lines[part_idx * block_size : (part_idx + 1) * block_size]
                block_ids.append(bid)
                block_file = os.path.join(out_dir, f"{bid}.tsv")
                with open(block_file, "w", encoding="utf-8") as f:
                    f.writelines(chunk)
                print(f"[partition_manifest]   -> {bid}: {len(chunk)} genomes ({block_file})", file=sys.stderr)

    # Process unclassified genomes if any
    if unclassified_list:
        print(f"[partition_manifest] Found {len(unclassified_list)} genomes without taxonomy annotation.", file=sys.stderr)
        if len(unclassified_list) <= block_size:
            bid = "tax_unclassified"
            block_ids.append(bid)
            block_file = os.path.join(out_dir, f"{bid}.tsv")
            with open(block_file, "w", encoding="utf-8") as f:
                f.writelines(unclassified_list)
            print(f"[partition_manifest]   -> {bid}: {len(unclassified_list)} genomes ({block_file})", file=sys.stderr)
        else:
            num_sub = (len(unclassified_list) + block_size - 1) // block_size
            for part_idx in range(num_sub):
                bid = f"tax_unclassified_part{part_idx:03d}"
                chunk = unclassified_list[part_idx * block_size : (part_idx + 1) * block_size]
                block_ids.append(bid)
                block_file = os.path.join(out_dir, f"{bid}.tsv")
                with open(block_file, "w", encoding="utf-8") as f:
                    f.writelines(chunk)
                print(f"[partition_manifest]   -> {bid}: {len(chunk)} genomes ({block_file})", file=sys.stderr)

    # Write block IDs list
    with open(out_list, "w", encoding="utf-8") as f:
        for bid in block_ids:
            f.write(f"{bid}\n")

    print(f"[partition_manifest] Successfully generated {len(block_ids)} taxonomy blocks. List written to {out_list}", file=sys.stderr)


def partition_random(entries, block_size, out_dir, out_list, shuffle, seed):
    """
    Evenly partitions entries into fixed-size random chunks.
    """
    total_mags = len(entries)
    if shuffle and total_mags > 1:
        rng = random.Random(seed)
        rng.shuffle(entries)
        print(f"[partition_manifest] Shuffled {total_mags} entries with seed {seed}", file=sys.stderr)

    num_blocks = (total_mags + block_size - 1) // block_size
    pad_len = max(3, len(str(num_blocks)))

    block_ids = []
    for i in range(num_blocks):
        chunk = entries[i * block_size : (i + 1) * block_size]
        if not chunk:
            break
        block_id = f"block_{i:0{pad_len}d}"
        block_ids.append(block_id)
        block_file = os.path.join(out_dir, f"{block_id}.tsv")
        with open(block_file, "w", encoding="utf-8") as f:
            f.writelines(chunk)
        print(f"[partition_manifest]   -> {block_id}: {len(chunk)} genomes ({block_file})", file=sys.stderr)

    with open(out_list, "w", encoding="utf-8") as f:
        for bid in block_ids:
            f.write(f"{bid}\n")

    print(f"[partition_manifest] Successfully generated {len(block_ids)} random blocks. List written to {out_list}", file=sys.stderr)


def main():
    args = parse_args()

    if not os.path.isfile(args.manifest):
        print(f"ERROR: Input manifest not found: {args.manifest}", file=sys.stderr)
        sys.exit(1)

    if args.block_size <= 0:
        print(f"ERROR: --block-size must be a positive integer, got {args.block_size}", file=sys.stderr)
        sys.exit(1)

    # Read manifest entries
    entries = []
    with open(args.manifest, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            entries.append(line if line.endswith("\n") else line + "\n")

    total_mags = len(entries)
    print(f"[partition_manifest] Read {total_mags} genomes from {args.manifest}", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_list)), exist_ok=True)

    if total_mags == 0:
        block_name = "block_000"
        block_file = os.path.join(args.out_dir, f"{block_name}.tsv")
        with open(block_file, "w", encoding="utf-8") as f:
            pass
        with open(args.out_list, "w", encoding="utf-8") as f:
            f.write(f"{block_name}\n")
        print(f"[partition_manifest] Created empty fallback block: {block_name}", file=sys.stderr)
        return

    # Try taxonomy-based partition if provided
    tax_map = load_taxonomy(args.taxonomy) if args.taxonomy else {}
    if tax_map:
        print(f"[partition_manifest] Loaded taxonomy annotations for {len(tax_map)} MAGs from {args.taxonomy}", file=sys.stderr)
        partition_by_taxonomy(entries, tax_map, args.block_size, args.min_phylum_size, args.out_dir, args.out_list)
    else:
        if args.taxonomy:
            print(f"[partition_manifest] WARNING: Taxonomy file {args.taxonomy} was empty or unparseable. Falling back to random partition.", file=sys.stderr)
        partition_random(entries, args.block_size, args.out_dir, args.out_list, args.shuffle, args.seed)


if __name__ == "__main__":
    main()
