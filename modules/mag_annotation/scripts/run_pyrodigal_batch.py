#!/usr/bin/env python3
"""Run the existing Pyrodigal caller over one Snakemake genome batch."""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def genome_id(path):
    name = os.path.basename(path)
    for suffix in (".fna.gz", ".fasta.gz", ".fa.gz", ".fna", ".fasta", ".fa"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return re.sub(r"\.(?:gz)$", "", name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--out-faa", required=True)
    parser.add_argument("--out-gff", required=True)
    parser.add_argument("--out-meta", required=True)
    parser.add_argument("--marker", required=True)
    parser.add_argument("--mode", choices=("auto", "meta", "single"), default="auto")
    parser.add_argument("--trans-table", type=int, default=11)
    parser.add_argument("--auto-complete-mode", choices=("meta", "single"), default="single")
    parser.add_argument("--auto-fallback-mode", choices=("meta", "single"), default="meta")
    args = parser.parse_args()

    ids = [genome_id(path) for path in args.genomes]
    if len(ids) != len(set(ids)):
        raise SystemExit("Duplicate genome IDs in batch would overwrite gene-calling outputs")

    caller = Path(__file__).resolve().parents[3] / "modules/mag_annotation/scripts/run_pyrodigal.py"
    for gid, path in zip(ids, args.genomes):
        cmd = [
            sys.executable, str(caller), "-i", path, "-g", gid,
            "--out_faa", os.path.join(args.out_faa, f"{gid}.faa.gz"),
            "--out_gff", os.path.join(args.out_gff, f"{gid}.gff.gz"),
            "--out_meta", os.path.join(args.out_meta, f"{gid}_protein_metadata.tsv.gz"),
            "-m", args.mode, "-t", str(args.trans_table),
            "--auto-complete-mode", args.auto_complete_mode,
            "--auto-fallback-mode", args.auto_fallback_mode,
        ]
        subprocess.run(cmd, check=True)

    marker = Path(args.marker)
    marker.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=marker.parent, delete=False) as handle:
        json.dump({"genome_count": len(ids), "genomes": ids}, handle)
        handle.write("\n")
        temp_path = handle.name
    os.replace(temp_path, marker)
    print(f"Completed gene calling for {len(ids)} genomes")


if __name__ == "__main__":
    main()
