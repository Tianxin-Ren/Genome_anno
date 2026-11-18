#!/usr/bin/env python3
import os
import sys
import argparse
import re
import glob

def find_mags(genomes_dir):
    dirs = [genomes_dir] if isinstance(genomes_dir, str) else genomes_dir
    exts = ["*.fna", "*.fasta", "*.fa", "*.fna.gz", "*.fasta.gz", "*.fa.gz"]
    mag_files = {}
    for d in dirs:
        if not os.path.exists(d):
            print(f"Genomes dir does not exist: {d}")
            continue
        for ext in exts:
            for f_path in glob.glob(os.path.join(d, ext)):
                basename = os.path.basename(f_path)
                if basename.endswith((".faa", ".faa.gz")):
                    continue
                name = re.sub(r'\.(fna|fasta|fa)(\.gz)?$', '', basename)
                if name not in mag_files:
                    mag_files[name] = f_path
    return mag_files

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genomes-dir", required=True, action="append")
    p.add_argument("--combined-tbl", required=True)
    args = p.parse_args()

    print("Scanning genomes dir...")
    mag_files = find_mags(args.genomes_dir)
    mag_names = sorted(mag_files.keys())
    mag_names_set = set(mag_names)
    print(f"Found {len(mag_names)} MAGs. Example names:")
    for name in mag_names[:5]:
        print(f"  - '{name}'")
    
    test_name = "DP8480007411BR_L01_100_bin.1"
    print(f"Is '{test_name}' in scanned MAGs? {test_name in mag_names_set}")

    print("\nReading combined HMM table...")
    if not os.path.exists(args.combined_tbl):
        print(f"ERROR: combined table does not exist: {args.combined_tbl}")
        return

    match_count = 0
    with open(args.combined_tbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            prot_id = parts[0]
            if "|" not in prot_id:
                continue
            mag_name = prot_id.split("|")[0]
            if test_name in prot_id:
                match_count += 1
                if match_count <= 5:
                    print(f"Match {match_count}:")
                    print(f"  Line: {line.strip()}")
                    print(f"  prot_id: '{prot_id}'")
                    print(f"  parsed mag_name: '{mag_name}'")
                    print(f"  parsed mag_name in set? {mag_name in mag_names_set}")

    print(f"\nTotal matches containing '{test_name}': {match_count}")

if __name__ == "__main__":
    main()
