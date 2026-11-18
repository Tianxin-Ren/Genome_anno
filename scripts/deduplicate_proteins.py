#!/usr/bin/env python3
"""Build a disk-backed, lossless 100%-identity protein catalog.

The combined FAA can contain hundreds of millions of records, so this script
does not retain all sequences or members in Python memory. SQLite stores one
representative per normalized SHA256 sequence and one row for every original
protein; the final FAA and mapping are streamed from that catalog.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import sqlite3
from functools import lru_cache


def fasta_records(path):
    header = None
    sequence = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:].strip()
                sequence = []
            elif header is not None:
                sequence.append(line.strip())
        if header is not None:
            yield header, "".join(sequence)


def normalize_sequence(sequence):
    return sequence.replace("*", "").replace(" ", "").upper()


def genome_from_protein(protein_id):
    return protein_id.split("|", 1)[0]


@lru_cache(maxsize=8)
def load_partial_status_for_genome(gff_root, genome):
    """Read partial flags for one genome at a time, not for the whole batch."""
    statuses = {}
    gff = os.path.join(gff_root, genome, f"{genome}.gff3")
    if not os.path.exists(gff):
        return statuses
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = dict(
                item.split("=", 1)
                for item in fields[8].split(";")
                if "=" in item
            )
            raw_id = attrs.get("ID", "")
            if raw_id.startswith("cds-"):
                raw_id = raw_id[4:]
            partial = attrs.get("partial", "")
            if partial in {"00", "0", "false", "False"}:
                status = "complete"
            elif partial in {"10", "01", "11", "1", "true", "True"}:
                status = "partial"
            else:
                status = "unknown"
            if raw_id:
                statuses[raw_id] = status
    return statuses


def init_database(path):
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.executescript(
        """
        CREATE TABLE representatives (
            protein_hash TEXT PRIMARY KEY,
            representative_id TEXT NOT NULL,
            representative_header TEXT NOT NULL,
            sequence TEXT NOT NULL,
            length INTEGER NOT NULL,
            partial_status TEXT NOT NULL
        );
        CREATE TABLE members (
            protein_id TEXT PRIMARY KEY,
            protein_hash TEXT NOT NULL,
            genome_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            length INTEGER NOT NULL,
            partial_status TEXT NOT NULL
        );
        CREATE INDEX members_hash_idx ON members(protein_hash);
        CREATE INDEX representatives_id_idx ON representatives(representative_id);
        """
    )
    return connection


def candidate_key(status, length, protein_id):
    # Complete > unknown > partial; then prefer longer and deterministic ID.
    status_rank = {"complete": 0, "unknown": 1, "partial": 2}.get(status, 1)
    return status_rank, -length, protein_id


def build_catalog(input_faa, output_faa, mapping_tsv, stats_tsv, database, gff_root=""):
    for output in (output_faa, mapping_tsv, stats_tsv, database):
        os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    if os.path.exists(database):
        os.remove(database)
    connection = init_database(database)
    total = 0
    skipped_empty = 0
    try:
        with connection:
            for header, raw_sequence in fasta_records(input_faa):
                protein_id = header.split()[0]
                sequence = normalize_sequence(raw_sequence)
                if not sequence:
                    skipped_empty += 1
                    continue
                genome = genome_from_protein(protein_id)
                gene_id = protein_id.split("|", 1)[1] if "|" in protein_id else protein_id
                genome_status = load_partial_status_for_genome(gff_root, genome)
                status = genome_status.get(gene_id, "unknown")
                protein_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
                existing = connection.execute(
                    "SELECT representative_id, representative_header, length, partial_status "
                    "FROM representatives WHERE protein_hash = ?",
                    (protein_hash,),
                ).fetchone()
                candidate = (protein_id, header, len(sequence), status)
                if existing is None or candidate_key(status, len(sequence), protein_id) < candidate_key(existing[3], existing[2], existing[0]):
                    connection.execute(
                        "INSERT OR REPLACE INTO representatives VALUES (?, ?, ?, ?, ?, ?)",
                        (protein_hash, protein_id, header, sequence, len(sequence), status),
                    )
                connection.execute(
                    "INSERT OR REPLACE INTO members VALUES (?, ?, ?, ?, ?, ?)",
                    (protein_id, protein_hash, genome, gene_id, len(sequence), status),
                )
                total += 1

        with open(output_faa, "w") as faa:
            for protein_hash, rep_id, header, sequence, _, _ in connection.execute(
                "SELECT protein_hash, representative_id, representative_header, sequence, length, partial_status "
                "FROM representatives ORDER BY protein_hash"
            ):
                faa.write(f">{header}\n")
                for start in range(0, len(sequence), 80):
                    faa.write(sequence[start:start + 80] + "\n")

        with open(mapping_tsv, "w", newline="") as mapping:
            writer = csv.writer(mapping, delimiter="\t")
            writer.writerow([
                "protein_id", "representative_id", "protein_hash", "genome_id",
                "gene_id", "length", "partial_status",
            ])
            writer.writerows(connection.execute(
                "SELECT m.protein_id, r.representative_id, m.protein_hash, m.genome_id, "
                "m.gene_id, m.length, m.partial_status "
                "FROM members m JOIN representatives r USING (protein_hash) "
                "ORDER BY m.protein_id"
            ))

        unique_count = connection.execute("SELECT COUNT(*) FROM representatives").fetchone()[0]
        with open(stats_tsv, "w") as stats:
            complete_count = connection.execute(
                "SELECT COUNT(*) FROM representatives WHERE partial_status = 'complete'"
            ).fetchone()[0]
            stats.write("metric\tvalue\n")
            stats.write(f"input_proteins\t{total}\n")
            stats.write(f"empty_proteins_skipped\t{skipped_empty}\n")
            stats.write(f"unique_proteins\t{unique_count}\n")
            stats.write(f"deduplication_ratio\t{(1 - unique_count / total) if total else 0:.8f}\n")
            stats.write(f"complete_representatives\t{complete_count}\n")
    finally:
        connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-faa", required=True)
    parser.add_argument("--output-faa", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--gff-root", default="")
    args = parser.parse_args()
    build_catalog(
        args.input_faa, args.output_faa, args.mapping, args.stats,
        args.database, args.gff_root,
    )


if __name__ == "__main__":
    main()
