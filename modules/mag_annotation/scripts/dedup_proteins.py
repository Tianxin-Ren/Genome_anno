#!/usr/bin/env python3
"""
dedup_proteins.py
Memory-bounded, collision-checked exact sequence deduplication for 50M+ proteins.
Key features:
- Canonical sequence normalization (uppercase, trailing stop '*' removal).
- Two-level collision check: SHA256 digest + length -> actual sequence equality verification.
- Comprehensive provenance tracking recording both raw proteome size and HMMER search universe Z.
- Outputs:
  1. unique_proteins.faa.gz
  2. protein_to_unique.tsv.gz
  3. run_metadata.json (statistical space provenance)
  4. dedup_stats.tsv
"""

import argparse
import gzip
import hashlib
import json
import os
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def stream_fasta(file_path):
    """Yield (header, sequence) pairs from FASTA file."""
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

def parse_genome_id(protein_id):
    if "___" in protein_id:
        return protein_id.split("___")[0]
    return "unknown"

def main():
    parser = argparse.ArgumentParser(description="Collision-checked exact protein deduplication.")
    parser.add_argument("-i", "--input_list", required=True, help="File with list of protein FASTA paths (or directory)")
    parser.add_argument("-o", "--out_dir", required=True, help="Directory for deduplicated outputs")
    parser.add_argument("-c", "--context", default="GTDB_r232_archaea", help="Proteome context name for provenance")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    input_files = []
    if os.path.isfile(args.input_list):
        with open(args.input_list, "r", encoding="utf-8") as f:
            for line in f:
                p = line.strip()
                if p and os.path.exists(p):
                    if os.path.isdir(p):
                        for root, _, files in os.walk(p):
                            for f in files:
                                if f.endswith((".faa", ".faa.gz", ".fasta", ".fa", ".fa.gz")):
                                    input_files.append(os.path.join(root, f))
                    else:
                        input_files.append(p)
    elif os.path.isdir(args.input_list):
        for root, _, files in os.walk(args.input_list):
            for f in files:
                if f.endswith((".faa", ".faa.gz", ".fasta", ".fa", ".fa.gz")):
                    input_files.append(os.path.join(root, f))
    input_files.sort()

    if not input_files:
        print(f"ERROR: No input FASTA files found from {args.input_list}", file=sys.stderr)
        sys.exit(1)

    unique_faa_path = os.path.join(args.out_dir, "unique_proteins.faa.gz")
    mapping_path = os.path.join(args.out_dir, "protein_to_unique.tsv.gz")
    metadata_json_path = os.path.join(args.out_dir, "run_metadata.json")
    stats_tsv_path = os.path.join(args.out_dir, "dedup_stats.tsv")

    # Primary hash index: hash_digest -> list of (seq_str, unique_id) to verify collisions
    hash_to_uniques = {}
    
    raw_count = 0
    unique_count = 0
    total_raw_residues = 0
    total_unique_residues = 0

    print(f"Starting collision-checked deduplication across {len(input_files)} FASTA files...")

    with gzip.open(unique_faa_path, "wt", encoding="utf-8") as u_out, \
         gzip.open(mapping_path, "wt", encoding="utf-8") as m_out:
        
        m_out.write("protein_id\tunique_id\taa_length\tgenome_id\n")

        for fidx, fpath in enumerate(input_files):
            for prot_id, raw_seq in stream_fasta(fpath):
                raw_count += 1
                seq_clean = raw_seq.upper().rstrip("*")
                seq_len = len(seq_clean)
                total_raw_residues += seq_len
                genome_id = parse_genome_id(prot_id)

                # SHA256 digest
                digest = hashlib.sha256(seq_clean.encode("ascii")).hexdigest()

                # Collision-check lookup
                assigned_unique_id = None
                candidates = hash_to_uniques.get(digest)
                if candidates is not None:
                    for cand_seq, cand_uid in candidates:
                        if len(cand_seq) == seq_len and cand_seq == seq_clean:
                            # Exact match confirmed
                            assigned_unique_id = cand_uid
                            break

                if assigned_unique_id is None:
                    # New unique sequence
                    unique_count += 1
                    total_unique_residues += seq_len
                    assigned_unique_id = f"UPROT_{unique_count:09d}"
                    if candidates is None:
                        hash_to_uniques[digest] = [(seq_clean, assigned_unique_id)]
                    else:
                        candidates.append((seq_clean, assigned_unique_id))
                    
                    # Write to unique FASTA
                    u_out.write(f">{assigned_unique_id} length={seq_len}\n")
                    for i in range(0, seq_len, 60):
                        u_out.write(seq_clean[i:i+60] + "\n")

                # Record mapping
                m_out.write(f"{prot_id}\t{assigned_unique_id}\t{seq_len}\t{genome_id}\n")

            if (fidx + 1) % 100 == 0 or (fidx + 1) == len(input_files):
                print(f"  Processed {fidx + 1}/{len(input_files)} files: {raw_count} raw proteins -> {unique_count} unique proteins.")

    dedup_ratio = 1.0 - (unique_count / raw_count) if raw_count > 0 else 0.0

    metadata = {
        "original_proteome_context": args.context,
        "raw_sequences_count": raw_count,
        "unique_sequences_count": unique_count,
        "raw_amino_acid_residues": total_raw_residues,
        "unique_amino_acid_residues": total_unique_residues,
        "dedup_ratio": round(dedup_ratio, 6),
        "hmmsearch_database_size_Z": unique_count,
        "hash_algorithm": "SHA256",
        "collision_check": "exact_character_verification"
    }

    with open(metadata_json_path, "w", encoding="utf-8") as jf:
        json.dump(metadata, jf, indent=2)

    with open(stats_tsv_path, "w", encoding="utf-8") as sf:
        sf.write("metric\tvalue\n")
        for k, v in metadata.items():
            sf.write(f"{k}\t{v}\n")

    print(f"Deduplication complete: {raw_count} raw proteins -> {unique_count} unique ({dedup_ratio*100:.2f}% reduction).")
    print(f"Provenance saved to {metadata_json_path}")

if __name__ == "__main__":
    main()
