#!/usr/bin/env python3
"""Expand HMM hits from representatives using a streaming SQLite catalog."""

from __future__ import annotations

import argparse
import csv
import os
import sqlite3
from collections import defaultdict


def load_members(mapping_path):
    # Compatibility fallback for older catalogs without catalog.sqlite.
    members = defaultdict(list)
    with open(mapping_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            members[row["representative_id"]].append(row["protein_id"])
    return members


def load_cluster_members(mapping_path):
    """Load MMseqs2 representative -> member IDs from cluster TSV."""
    members = defaultdict(list)
    with open(mapping_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                members[fields[0]].append(fields[1])
    return members


def expand(input_tbl, output_tbl, database=None, mapping_tsv=None, cluster_mapping=None):
    connection = sqlite3.connect(database) if database and os.path.exists(database) else None
    fallback = None if connection else load_members(mapping_tsv)
    cluster_members = load_cluster_members(cluster_mapping) if cluster_mapping else None
    os.makedirs(os.path.dirname(os.path.abspath(output_tbl)), exist_ok=True)
    expanded = 0
    try:
        with open(input_tbl) as source, open(output_tbl, "w") as target:
            for line in source:
                if line.startswith("#") or not line.strip():
                    target.write(line)
                    continue
                fields = line.split()
                if len(fields) < 3:
                    target.write(line)
                    continue
                representative = fields[2]
                unique_targets = cluster_members.get(representative, [representative]) if cluster_members else [representative]
                targets = []
                for unique_target in unique_targets:
                    if connection:
                        targets.extend(row[0] for row in connection.execute(
                            "SELECT m.protein_id FROM members m "
                            "JOIN representatives r USING (protein_hash) "
                            "WHERE r.representative_id = ? ORDER BY m.protein_id",
                            (unique_target,),
                        ))
                    else:
                        targets.extend(fallback.get(unique_target, [unique_target]))
                for protein_id in targets or [representative]:
                    out_fields = list(fields)
                    out_fields[2] = protein_id
                    target.write("\t".join(out_fields) + "\n")
                    expanded += 1
    finally:
        if connection:
            connection.close()
    return expanded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--database")
    parser.add_argument("--mapping")
    parser.add_argument("--cluster-mapping")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if not args.database and not args.mapping:
        parser.error("one of --database or --mapping is required")
    count = expand(args.input, args.output, args.database, args.mapping, args.cluster_mapping)
    print(f"Expanded {count:,} HMM rows from representative hits.")


if __name__ == "__main__":
    main()
