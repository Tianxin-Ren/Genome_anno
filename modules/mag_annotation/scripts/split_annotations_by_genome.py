#!/usr/bin/env python3
"""Write the all-hit KEGG table separately for each genome."""

import argparse
import json
import os
import re
import tempfile
from collections import OrderedDict


def _genome_ids_from_protein_dir(proteins_dir):
    ids = set()
    for root, _, files in os.walk(proteins_dir):
        for name in files:
            for suffix in (".faa.gz", ".faa", ".fasta.gz", ".fasta", ".fa.gz", ".fa"):
                if name.endswith(suffix):
                    ids.add(name[:-len(suffix)])
                    break
    return ids


def split_table(source, destination_dir, label, genome_column="genome_id", genome_ids=None):
    os.makedirs(destination_dir, exist_ok=True)
    # Keep the descriptor count bounded; a 60k-genome run cannot keep one
    # output handle open per genome on a typical cluster node.
    handles = OrderedDict()
    max_open_handles = 256
    genomes = set(genome_ids or [])
    written = set()
    try:
        with open(source, "r", encoding="utf-8") as src:
            header = src.readline()
            if not header:
                raise ValueError(f"Input table is empty: {source}")
            columns = header.rstrip("\r\n").split("\t")
            if genome_column not in columns:
                raise ValueError(f"Missing {genome_column!r} column in {source}")
            genome_index = columns.index(genome_column)
            for line in src:
                if not line.strip():
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if genome_index >= len(fields) or not fields[genome_index]:
                    continue
                genome = fields[genome_index]
                if genome not in handles:
                    if len(handles) >= max_open_handles:
                        _, old_handle = handles.popitem(last=False)
                        old_handle.close()
                    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", genome)
                    path = os.path.join(destination_dir, f"{safe_name}.{label}.tsv")
                    # Append after an LRU eviction so a genome can be resumed
                    # without losing its header or previously written rows.
                    first_seen = genome not in written
                    handles[genome] = open(path, "w" if first_seen else "a", encoding="utf-8")
                    if first_seen:
                        handles[genome].write(header)
                    written.add(genome)
                else:
                    handles.move_to_end(genome)
                handles[genome].write(line)
    finally:
        for handle in handles.values():
            handle.close()
    # Include genomes with no KO rows as header-only files.
    for genome in sorted(genomes - written):
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", genome)
        path = os.path.join(destination_dir, f"{safe_name}.{label}.tsv")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(header)
    return genomes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kegg", required=True, help="All-hit KEGG table (protein_KO_all_strict.tsv)")
    parser.add_argument("--proteins_dir", default="", help="Directory of per-genome FAA files, used to include zero-hit genomes")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--marker", required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    genome_ids = _genome_ids_from_protein_dir(args.proteins_dir) if args.proteins_dir else None
    kegg_genomes = split_table(args.kegg, os.path.join(args.out_dir, "kegg"), "kegg", genome_ids=genome_ids)

    os.makedirs(os.path.dirname(os.path.abspath(args.marker)), exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=os.path.dirname(os.path.abspath(args.marker)), delete=False) as tmp:
        json.dump({"genome_count": len(kegg_genomes), "genomes": sorted(kegg_genomes)}, tmp, indent=2)
        tmp.write("\n")
        temp_path = tmp.name
    os.replace(temp_path, args.marker)
    print(f"Wrote KEGG per-genome tables for {len(kegg_genomes)} genomes under {args.out_dir}")


if __name__ == "__main__":
    main()
