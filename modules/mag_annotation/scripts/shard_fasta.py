#!/usr/bin/env python3
"""
shard_fasta.py
Sequence-complexity-aware residue-balanced FASTA sharding.
Balances HMMER computational load primarily by total amino-acid residues (~100-150M aa),
with a secondary upper guardrail on protein count (default 500k).
Outputs:
  - shards/shard_{shard_id}.faa.gz
  - shards/shards_manifest.tsv (with mean, median, and max length statistics to guard against tail skew)
  - shards/shard_list.txt (plain list of shard IDs for Snakemake dynamic checkpoints)
"""

import argparse
import gzip
import hashlib
import os
import statistics
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def stream_fasta(file_path):
    curr_id = None
    curr_seq = []
    with open_file(file_path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if curr_id is not None:
                    yield curr_id, "".join(curr_seq)
                curr_id = line[1:].split()[0]
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_id is not None:
            yield curr_id, "".join(curr_seq)

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description="Residue-balanced FASTA sharding.")
    parser.add_argument("-i", "--input_fasta", required=True, help="Input FASTA (.faa or .faa.gz)")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for shards")
    parser.add_argument("--target_aa", type=int, default=120000000, help="Target amino acid residues per shard (default: 120M)")
    parser.add_argument("--target_proteins", type=int, default=500000, help="Maximum proteins per shard (default: 500,000)")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    manifest_path = os.path.join(args.out_dir, "shards_manifest.tsv")
    shard_list_path = os.path.join(args.out_dir, "shard_list.txt")

    current_shard_idx = 1
    current_proteins = [] # list of (header, seq, len)
    current_residues = 0

    shard_stats = []

    def flush_shard(shard_idx, proteins):
        if not proteins:
            return None
        shard_id = f"shard_{shard_idx:04d}"
        shard_file = os.path.join(args.out_dir, f"{shard_id}.faa.gz")
        
        lengths = [p[2] for p in proteins]
        p_count = len(proteins)
        aa_count = sum(lengths)
        mean_len = round(statistics.mean(lengths), 2)
        median_len = round(statistics.median(lengths), 2)
        max_len = max(lengths)

        with gzip.open(shard_file, "wt", encoding="utf-8") as out:
            for p_id, seq, _ in proteins:
                out.write(f">{p_id}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")

        compressed_size = os.path.getsize(shard_file)
        sha = compute_sha256(shard_file)

        return {
            "shard_id": shard_id,
            "shard_file": shard_file,
            "protein_count": p_count,
            "aa_count": aa_count,
            "mean_length": mean_len,
            "median_length": median_len,
            "max_length": max_len,
            "compressed_size": compressed_size,
            "sha256": sha
        }

    print(f"Streaming and sharding {args.input_fasta} (Target: {args.target_aa/1e6:.1f}M aa, Max: {args.target_proteins} prot/shard)...")

    for p_id, seq in stream_fasta(args.input_fasta):
        seq_len = len(seq)
        current_proteins.append((p_id, seq, seq_len))
        current_residues += seq_len

        if current_residues >= args.target_aa or len(current_proteins) >= args.target_proteins:
            st = flush_shard(current_shard_idx, current_proteins)
            if st:
                shard_stats.append(st)
                print(f"  Written {st['shard_id']}: {st['protein_count']} proteins, {st['aa_count']/1e6:.2f}M aa (max_len={st['max_length']})")
                current_shard_idx += 1
            current_proteins = []
            current_residues = 0

    # Flush remainder
    if current_proteins:
        st = flush_shard(current_shard_idx, current_proteins)
        if st:
            shard_stats.append(st)
            print(f"  Written final {st['shard_id']}: {st['protein_count']} proteins, {st['aa_count']/1e6:.2f}M aa (max_len={st['max_length']})")

    # Write manifest and shard_list
    with open(manifest_path, "w", encoding="utf-8") as mf, \
         open(shard_list_path, "w", encoding="utf-8") as lf:
        
        mf.write("shard_id\tprotein_count\taa_count\tmean_length\tmedian_length\tmax_length\tcompressed_size\tsha256\n")
        for s in shard_stats:
            mf.write(f"{s['shard_id']}\t{s['protein_count']}\t{s['aa_count']}\t{s['mean_length']}\t{s['median_length']}\t{s['max_length']}\t{s['compressed_size']}\t{s['sha256']}\n")
            lf.write(f"{s['shard_id']}\n")

    print(f"Sharding complete. Total shards: {len(shard_stats)}. Manifest written to {manifest_path}")

if __name__ == "__main__":
    main()
