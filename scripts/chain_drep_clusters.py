#!/usr/bin/env python3
"""
Chain block-level dRep clusters with consolidated pool-level dRep clusters.

Provenance mapping:
  Original MAG -> Block Secondary Cluster -> Block Winner -> Pool Cluster -> Final Winner
"""

import os
import sys
import csv
import argparse
import re
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Chain block-level dRep clusters with pool dRep clusters into a global Cdb.csv and Widb.csv."
    )
    parser.add_argument(
        "--block-cdbs",
        nargs="+",
        required=True,
        help="List of paths to block Cdb.csv files.",
    )
    parser.add_argument(
        "--block-widbs",
        nargs="+",
        default=[],
        help="List of paths to block Widb.csv files (matching order with block-cdbs, or auto-paired).",
    )
    parser.add_argument(
        "--pool-cdb",
        default=None,
        help="Path to consolidated pool Cdb.csv (optional if single block).",
    )
    parser.add_argument(
        "--pool-widb",
        default=None,
        help="Path to consolidated pool Widb.csv (optional if single block).",
    )
    parser.add_argument(
        "--out-cdb",
        required=True,
        help="Path to write final consolidated master Cdb.csv.",
    )
    parser.add_argument(
        "--out-widb",
        required=True,
        help="Path to write final consolidated master Widb.csv.",
    )
    parser.add_argument(
        "--out-provenance",
        default=None,
        help="Optional path to write full provenance TSV mapping table.",
    )
    return parser.parse_args()


def read_csv_rows(path):
    """Read CSV file returning (headers, list of row dicts)."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return [], []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            return [], []
        headers = [h.strip() for h in header]
        rows = []
        for r in reader:
            if not r:
                continue
            row_dict = {headers[i]: r[i].strip() for i in range(min(len(headers), len(r)))}
            rows.append(row_dict)
    return headers, rows


def find_column(headers, candidate_names):
    """Find the first matching column name in headers (case-insensitive)."""
    h_lower = {h.lower(): h for h in headers}
    for c in candidate_names:
        if c.lower() in h_lower:
            return h_lower[c.lower()]
    return None


def norm_genome_key(g):
    """Normalize genome name for matching across directories/stems."""
    if not g:
        return ""
    base = os.path.basename(g.strip())
    return re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", base)


def main():
    args = parse_args()

    # Create destination parent directories
    os.makedirs(os.path.dirname(os.path.abspath(args.out_cdb)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_widb)), exist_ok=True)
    if args.out_provenance:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_provenance)), exist_ok=True)

    block_cdbs = [p for p in args.block_cdbs if os.path.isfile(p)]
    if not block_cdbs:
        print("ERROR: No valid block Cdb.csv files found.", file=sys.stderr)
        sys.exit(1)

    # 1. Handle Single Block fast-path
    if len(block_cdbs) == 1 and (not args.pool_cdb or not os.path.isfile(args.pool_cdb)):
        print(f"[chain_drep_clusters] Single block detected ({block_cdbs[0]}), directly copying.")
        import shutil
        shutil.copy2(block_cdbs[0], args.out_cdb)
        if args.block_widbs and os.path.isfile(args.block_widbs[0]):
            shutil.copy2(args.block_widbs[0], args.out_widb)
        elif args.pool_widb and os.path.isfile(args.pool_widb):
            shutil.copy2(args.pool_widb, args.out_widb)
        print("[chain_drep_clusters] Single block pass-through complete.")
        return

    # 2. Parse Pool dRep outputs
    # pool_genome_to_cluster: normalized_name -> (secondary_cluster, primary_cluster, raw_genome_name)
    pool_genome_to_cluster = {}
    pool_cdb_headers, pool_cdb_rows = read_csv_rows(args.pool_cdb) if args.pool_cdb else ([], [])
    if pool_cdb_rows:
        g_col = find_column(pool_cdb_headers, ["genome", "Genome"]) or "genome"
        sec_col = find_column(pool_cdb_headers, ["secondary_cluster", "Secondary_cluster", "cluster"]) or "secondary_cluster"
        prim_col = find_column(pool_cdb_headers, ["primary_cluster", "Primary_cluster"]) or "primary_cluster"

        for row in pool_cdb_rows:
            raw_g = row.get(g_col, "")
            k = norm_genome_key(raw_g)
            sec = row.get(sec_col, "")
            prim = row.get(prim_col, "")
            pool_genome_to_cluster[k] = (sec, prim, raw_g)

    # pool_cluster_to_winner: secondary_cluster -> (final_winner, full_row_dict)
    pool_cluster_to_winner = {}
    pool_widb_headers, pool_widb_rows = read_csv_rows(args.pool_widb) if args.pool_widb else ([], [])
    if pool_widb_rows:
        w_g_col = find_column(pool_widb_headers, ["genome", "Genome"]) or "genome"
        w_c_col = find_column(pool_widb_headers, ["cluster", "secondary_cluster", "Secondary_cluster"]) or "cluster"
        for row in pool_widb_rows:
            winner = row.get(w_g_col, "")
            clust = row.get(w_c_col, "")
            if clust:
                pool_cluster_to_winner[clust] = (winner, row)
    elif args.pool_widb:
        # Fallback to Wdb.csv if Widb.csv was empty (e.g. --ignoreGenomeQuality)
        wdb_fallback = os.path.join(os.path.dirname(args.pool_widb), "Wdb.csv")
        w_hdrs, w_rows = read_csv_rows(wdb_fallback)
        if w_rows:
            w_g_col = find_column(w_hdrs, ["genome", "Genome"]) or "genome"
            w_c_col = find_column(w_hdrs, ["cluster", "secondary_cluster", "Secondary_cluster"]) or "cluster"
            for row in w_rows:
                winner = row.get(w_g_col, "")
                clust = row.get(w_c_col, "")
                if clust:
                    pool_cluster_to_winner[clust] = (winner, row)

    # 3. Parse each Block's Cdb.csv and Widb.csv
    # Map block_widbs by parent directory or index
    widb_by_block = {}
    for w in args.block_widbs:
        if os.path.isfile(w):
            # Try to associate by block dirname
            bname = os.path.basename(os.path.dirname(os.path.dirname(w)))
            if not bname or bname in ["drep", "round1", "round1_blocks"]:
                bname = os.path.basename(os.path.dirname(w))
            widb_by_block[bname] = w

    provenance_records = []
    final_cdb_rows = []

    for idx, b_cdb in enumerate(block_cdbs):
        b_name = os.path.basename(os.path.dirname(os.path.dirname(b_cdb)))
        if not b_name or b_name in ["drep", "round1", "round1_blocks"]:
            b_name = f"block_{idx:03d}"

        # Find corresponding Widb
        b_widb = widb_by_block.get(b_name)
        if not b_widb and idx < len(args.block_widbs) and os.path.isfile(args.block_widbs[idx]):
            b_widb = args.block_widbs[idx]

        # Parse block Widb to map block cluster -> block winner
        block_cluster_to_winner = {}
        if b_widb and os.path.isfile(b_widb):
            w_hdrs, w_rows = read_csv_rows(b_widb)
            g_col = find_column(w_hdrs, ["genome", "Genome"]) or "genome"
            c_col = find_column(w_hdrs, ["cluster", "secondary_cluster"]) or "cluster"
            for r in w_rows:
                g = r.get(g_col, "")
                c = r.get(c_col, "")
                if c:
                    block_cluster_to_winner[c] = g

        # Fallback to block Wdb.csv if Widb had no entries
        if not block_cluster_to_winner and b_widb:
            wdb_fallback = os.path.join(os.path.dirname(b_widb), "Wdb.csv")
            if os.path.isfile(wdb_fallback):
                w_hdrs, w_rows = read_csv_rows(wdb_fallback)
                g_col = find_column(w_hdrs, ["genome", "Genome"]) or "genome"
                c_col = find_column(w_hdrs, ["cluster", "secondary_cluster"]) or "cluster"
                for r in w_rows:
                    g = r.get(g_col, "")
                    c = r.get(c_col, "")
                    if c:
                        block_cluster_to_winner[c] = g

        # Fallback to dereplicated_genomes directory inspection
        derep_dir = os.path.join(os.path.dirname(os.path.dirname(b_cdb)), "dereplicated_genomes")
        winner_files = set()
        if os.path.isdir(derep_dir):
            winner_files = {norm_genome_key(f) for f in os.listdir(derep_dir)}

        # Parse block Cdb
        c_hdrs, c_rows = read_csv_rows(b_cdb)
        if not c_rows:
            continue
        g_col = find_column(c_hdrs, ["genome", "Genome"]) or "genome"
        sec_col = find_column(c_hdrs, ["secondary_cluster", "cluster"]) or "secondary_cluster"
        prim_col = find_column(c_hdrs, ["primary_cluster"]) or "primary_cluster"

        # Complete block_cluster_to_winner using winner_files if needed
        if winner_files:
            for r in c_rows:
                raw_mag = r.get(g_col, "")
                norm_mag = norm_genome_key(raw_mag)
                b_sec = r.get(sec_col, "")
                if b_sec and b_sec not in block_cluster_to_winner:
                    if norm_mag in winner_files:
                        block_cluster_to_winner[b_sec] = raw_mag

        for r in c_rows:
            raw_mag = r.get(g_col, "")
            norm_mag = norm_genome_key(raw_mag)
            b_sec = r.get(sec_col, "")
            b_prim = r.get(prim_col, "")

            # Block winner for this cluster
            b_winner = block_cluster_to_winner.get(b_sec, raw_mag)
            norm_b_winner = norm_genome_key(b_winner)

            # Chain to Pool cluster
            if norm_b_winner in pool_genome_to_cluster:
                pool_sec, pool_prim, _ = pool_genome_to_cluster[norm_b_winner]
                final_winner, _ = pool_cluster_to_winner.get(pool_sec, (b_winner, {}))
            else:
                # Fallback: keep block cluster with prefix to ensure global uniqueness
                pool_sec = f"{b_name}_{b_sec}"
                pool_prim = f"{b_name}_{b_prim}" if b_prim else pool_sec
                final_winner = b_winner

            final_cdb_rows.append({
                "genome": raw_mag,
                "secondary_cluster": pool_sec,
                "primary_cluster": pool_prim,
            })

            provenance_records.append({
                "genome": raw_mag,
                "block_id": b_name,
                "block_cluster": b_sec,
                "block_winner": b_winner,
                "final_secondary_cluster": pool_sec,
                "final_primary_cluster": pool_prim,
                "final_winner": final_winner,
            })

    # 4. Write master Cdb.csv
    print(f"[chain_drep_clusters] Writing {len(final_cdb_rows)} entries to master Cdb: {args.out_cdb}")
    with open(args.out_cdb, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["genome", "secondary_cluster", "primary_cluster"])
        for r in final_cdb_rows:
            writer.writerow([r["genome"], r["secondary_cluster"], r["primary_cluster"]])

    # 5. Write master Widb.csv (from pool_widb, or consolidated)
    if args.pool_widb and os.path.isfile(args.pool_widb) and os.path.getsize(args.pool_widb) > 0:
        import shutil
        shutil.copy2(args.pool_widb, args.out_widb)
        print(f"[chain_drep_clusters] Copied pool Widb to: {args.out_widb}")
    else:
        # Fallback: create Widb from winners
        with open(args.out_widb, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["genome", "cluster"])
            seen_clusters = set()
            for rec in provenance_records:
                c = rec["final_secondary_cluster"]
                if c not in seen_clusters:
                    seen_clusters.add(c)
                    writer.writerow([rec["final_winner"], c])
        print(f"[chain_drep_clusters] Generated synthesized Widb to: {args.out_widb}")

    # 6. Write provenance TSV if requested
    if args.out_provenance:
        print(f"[chain_drep_clusters] Writing provenance table: {args.out_provenance}")
        with open(args.out_provenance, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "genome",
                "block_id",
                "block_cluster",
                "block_winner",
                "final_secondary_cluster",
                "final_primary_cluster",
                "final_winner",
            ])
            for rec in provenance_records:
                writer.writerow([
                    rec["genome"],
                    rec["block_id"],
                    rec["block_cluster"],
                    rec["block_winner"],
                    rec["final_secondary_cluster"],
                    rec["final_primary_cluster"],
                    rec["final_winner"],
                ])

    print("[chain_drep_clusters] Successfully finished cluster chaining.")


if __name__ == "__main__":
    main()
