#!/usr/bin/env python3
"""Run MMseqs2 linclust clustering on 100%-unique proteins."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess


def run(input_faa, output_dir, mmseqs, identity, coverage, coverage_mode, threads):
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, "cluster99")
    tmp_dir = os.path.join(output_dir, "tmp_cluster99")
    subprocess.run([
        mmseqs, "easy-linclust", input_faa, prefix, tmp_dir,
        "--min-seq-id", str(identity), "-c", str(coverage),
        "--cov-mode", str(coverage_mode), "--alignment-mode", "3",
        "--threads", str(threads),
    ], check=True)
    rep = f"{prefix}_rep_seq.fasta"
    clusters = f"{prefix}_cluster.tsv"
    if not os.path.exists(rep) or not os.path.exists(clusters):
        raise FileNotFoundError("MMseqs2 did not produce expected representative/cluster files")
    shutil.copyfile(rep, os.path.join(output_dir, "representatives.faa"))
    shutil.copyfile(clusters, os.path.join(output_dir, "clusters.tsv"))
    with open(os.path.join(output_dir, "parameters.tsv"), "w") as handle:
        handle.write("parameter\tvalue\nidentity\t%s\ncoverage\t%s\ncoverage_mode\t%s\n" % (identity, coverage, coverage_mode))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-faa", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--mmseqs", default="mmseqs")
    parser.add_argument("--identity", type=float, default=0.99)
    parser.add_argument("--coverage", type=float, default=0.95)
    parser.add_argument("--coverage-mode", type=int, default=0)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    run(args.input_faa, args.output_dir, args.mmseqs, args.identity, args.coverage, args.coverage_mode, args.threads)


if __name__ == "__main__":
    main()
