#!/usr/bin/env python3
"""
recover_terminal_fastq.py
--------------------------------------------------------------------------------
Extract reads matching target QNAMEs from raw FASTQ (.fq/.fastq or .fq.gz/.fastq.gz)
and write to gzipped FASTQ. Used as fallback when seqkit is unavailable.
"""

import os
import sys
import gzip
import argparse


def main():
    parser = argparse.ArgumentParser(description="Recover terminal reads matching QNAMEs from raw FASTQ.")
    parser.add_argument("--qnames", required=True, help="File containing target QNAMEs, one per line")
    parser.add_argument("--raw-reads", required=True, help="Path to raw FASTQ (plain or .gz)")
    parser.add_argument("--out-reads", required=True, help="Output .fq.gz path")

    args = parser.parse_args()

    target_ids = set()
    if os.path.exists(args.qnames):
        with open(args.qnames) as f:
            for line in f:
                item = line.strip()
                if item:
                    target_ids.add(item)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_reads)), exist_ok=True)

    if not target_ids or not os.path.exists(args.raw_reads):
        # Empty output
        with gzip.open(args.out_reads, "wt") as fout:
            pass
        return

    is_gz = args.raw_reads.endswith(".gz")
    opener = gzip.open if is_gz else open

    count = 0
    with opener(args.raw_reads, "rt") as fin, gzip.open(args.out_reads, "wt") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline()
            p = fin.readline()
            q = fin.readline()
            if not h.startswith("@"):
                continue
            qname = h[1:].split()[0]
            if qname in target_ids:
                fout.write(h + s + p + q)
                count += 1

    print(f"[recover_terminal_fastq] Recovered {count} reads matching {len(target_ids)} target QNAMEs")


if __name__ == "__main__":
    main()
