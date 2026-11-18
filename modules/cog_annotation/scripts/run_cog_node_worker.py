#!/usr/bin/env python3
"""
run_cog_node_worker.py
High-performance node-level bundle manager for RPS-BLAST Cog_LE searches on SLURM 64c512g nodes.

Key features:
1. Node-local scratch staging:
   Copies Cog_LE.* database files to $SLURM_TMPDIR (or /tmp/cog_db) ONCE per node allocation.
2. Parallel worker execution:
   Launches concurrent worker processes across all 64 node cores (e.g. 8 workers x 8 threads or 5 workers x 12 threads).
3. Two-phase shard processing:
   - Phase A: High-throughput RPS-BLAST against staged local database.
   - Phase B: Footprint-aware assignment and arbitration (assign_cog_footprints.py).
"""

import argparse
import glob
import gzip
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import tempfile
import time

def stage_cog_database(shared_db_prefix, scratch_dir):
    """
    Stage Cog_LE.* files to node-local scratch directory once per node allocation.
    Returns: local path prefix to staged Cog_LE database.
    """
    src_dir = os.path.dirname(os.path.abspath(shared_db_prefix))
    base_prefix = os.path.basename(shared_db_prefix)
    src_files = glob.glob(os.path.join(src_dir, f"{base_prefix}.*"))

    if not src_files:
        raise FileNotFoundError(f"ERROR: No Cog_LE database files found matching '{shared_db_prefix}.*'")

    local_db_dir = os.path.join(scratch_dir, "cog_db")
    os.makedirs(local_db_dir, exist_ok=True)
    local_prefix = os.path.join(local_db_dir, base_prefix)

    # Check if already staged
    all_staged = True
    for sf in src_files:
        bf = os.path.basename(sf)
        dest = os.path.join(local_db_dir, bf)
        if not (os.path.isfile(dest) and os.path.getsize(dest) == os.path.getsize(sf)):
            all_staged = False
            break

    if all_staged:
        print(f"INFO: Cog_LE database already staged at {local_prefix}")
        return local_prefix

    print(f"INFO: Staging Cog_LE database from {src_dir} to {local_db_dir}...")
    t0 = time.time()
    for sf in src_files:
        bf = os.path.basename(sf)
        dest = os.path.join(local_db_dir, bf)
        tmp_dest = dest + f".tmp.{os.getpid()}"
        shutil.copyfile(sf, tmp_dest)
        os.rename(tmp_dest, dest)

    dur = time.time() - t0
    total_mb = sum(os.path.getsize(sf) for sf in src_files) / (1024 * 1024)
    print(f"INFO: Staged {len(src_files)} files ({total_mb:.1f} MB in {dur:.2f} s).")
    return local_prefix

def execute_shard_task(task_args):
    """Worker task executed per shard."""
    shard_fasta, raw_hits_tsv, assigned_tsv, db_prefix, threads, evalue, pssm_mapping, assign_script, blast_env = task_args
    shard_name = os.path.basename(shard_fasta).replace(".faa.gz", "").replace(".faa", "")
    print(f"[{os.getpid()}] Starting shard {shard_name} on {threads} threads...")

    t0 = time.time()
    os.makedirs(os.path.dirname(os.path.abspath(raw_hits_tsv)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(assigned_tsv)), exist_ok=True)

    # 1. RPS-BLAST
    query_path = shard_fasta
    cleanup_query = False
    if shard_fasta.endswith(".gz"):
        local_scratch = os.environ.get("SLURM_TMPDIR") or tempfile.gettempdir()
        tmp_fd, query_path = tempfile.mkstemp(prefix=f"{shard_name}_", suffix=".faa", dir=local_scratch)
        os.close(tmp_fd)
        with gzip.open(shard_fasta, "rb") as f_in, open(query_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        cleanup_query = True

    try:
        blast_cmd = ["rpsblast"]
        if blast_env:
            if os.path.isdir(blast_env):
                blast_cmd = ["conda", "run", "-p", blast_env, "rpsblast"]
            else:
                blast_cmd = ["conda", "run", "-n", blast_env, "rpsblast"]

        blast_args = blast_cmd + [
            "-db", db_prefix,
            "-query", query_path,
            "-out", raw_hits_tsv,
            "-evalue", str(evalue),
            "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
            "-num_threads", str(threads),
            "-seg", "no",
            "-comp_based_stats", "1"
        ]

        res = subprocess.run(blast_args, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"RPS-BLAST failed for {shard_name}:\n{res.stderr}")
    finally:
        if cleanup_query and os.path.exists(query_path):
            try:
                os.remove(query_path)
            except OSError:
                pass

    t_blast = time.time() - t0

    # 2. Footprint Assignment & Arbitration
    t1 = time.time()
    assign_cmd = [
        sys.executable, assign_script,
        "-i", raw_hits_tsv,
        "-m", pssm_mapping,
        "-o", assigned_tsv,
        "--evalue", str(evalue)
    ]

    res2 = subprocess.run(assign_cmd, capture_output=True, text=True)
    if res2.returncode != 0:
        raise RuntimeError(f"Assignment failed for {shard_name}:\n{res2.stderr}")

    t_assign = time.time() - t1
    print(f"[{os.getpid()}] Finished {shard_name}: BLAST {t_blast:.1f}s, Assign {t_assign:.1f}s")
    return shard_name

def main():
    parser = argparse.ArgumentParser(description="Node-level worker manager for COG annotation.")
    parser.add_argument("--shard_list", required=True, help="File listing shard FASTA paths to process on this node")
    parser.add_argument("--db_prefix", required=True, help="Shared filesystem Cog_LE database prefix")
    parser.add_argument("--pssm_mapping", required=True, help="cdd_pssm_mapping.tsv.gz from Gate 0")
    parser.add_argument("--raw_out_dir", required=True, help="Directory for raw RPS-BLAST hits")
    parser.add_argument("--assigned_out_dir", required=True, help="Directory for assigned footprint TSVs")
    parser.add_argument("--assign_script", required=True, help="Path to assign_cog_footprints.py")
    parser.add_argument("--scratch_dir", default="", help="Scratch directory (defaults to $SLURM_TMPDIR or /tmp)")
    parser.add_argument("--workers", type=int, default=8, help="Number of concurrent shard workers (default: 8)")
    parser.add_argument("--threads_per_worker", type=int, default=8, help="Threads per worker (default: 8)")
    parser.add_argument("--evalue", type=float, default=1e-4, help="E-value threshold (default: 1e-4)")
    parser.add_argument("--blast_env", default="", help="Conda environment for rpsblast")

    args = parser.parse_args()

    # Determine scratch directory
    scratch = args.scratch_dir or os.environ.get("SLURM_TMPDIR") or "/tmp"
    local_db_prefix = stage_cog_database(args.db_prefix, scratch)

    # Read shards to process
    shards = []
    with open(args.shard_list, "r", encoding="utf-8") as f:
        for line in f:
            p = line.strip()
            if p and os.path.exists(p):
                shards.append(p)

    if not shards:
        print(f"WARNING: No valid shard FASTAs found in {args.shard_list}")
        sys.exit(0)

    print(f"INFO: Node worker running {len(shards)} shards with {args.workers} workers x {args.threads_per_worker} threads.")

    tasks = []
    for sf in shards:
        bname = os.path.basename(sf).replace(".faa.gz", "").replace(".faa", "")
        raw_tsv = os.path.join(args.raw_out_dir, f"{bname}.hits.tsv")
        assign_tsv = os.path.join(args.assigned_out_dir, f"{bname}.assigned.tsv.gz")
        tasks.append((
            sf, raw_tsv, assign_tsv, local_db_prefix, args.threads_per_worker,
            args.evalue, args.pssm_mapping, args.assign_script, args.blast_env
        ))

    # Run worker pool
    with mp.Pool(processes=args.workers) as pool:
        results = pool.map(execute_shard_task, tasks)

    print(f"SUCCESS: Node worker finished processing {len(results)} shards.")

if __name__ == "__main__":
    main()
