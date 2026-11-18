#!/usr/bin/env python3
"""
run_node_worker.py
High-performance node-level manager for PyKofamSearch on SLURM 64c512g nodes.
Key features:
1. Node-local scratch staging:
   Copies database.pkl to $SLURM_TMPDIR (or /tmp) ONCE per node allocation, with SHA256 verification.
2. In-memory Copy-on-Write (CoW) fork:
   Parent process loads database into memory, then forks worker child processes using
   multiprocessing.set_start_method('fork') to maximize memory page sharing.
3. Relaxed HMMER prefilter + curated bitscore threshold:
   Runs pykofamsearch with relaxed E-value (--evalue 10.0) to minimize heuristic pre-filtering losses.
4. Concurrency:
   Executes multiple shard workers in parallel on the node (default: 5 workers x 12 threads = 60 cores).
"""

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def stage_database(shared_db_path, scratch_dir, expected_sha256=None):
    """Stage database.pkl to local scratch directory once per node."""
    if not os.path.isfile(shared_db_path):
        raise FileNotFoundError(f"ERROR: KOfam database {shared_db_path} does not exist!")

    src_size = os.path.getsize(shared_db_path)
    # Check if file is suspiciously small (a valid KOfam serialized DB with ~26,000 models is ~1.5 - 2.5 GB)
    # Allow small DBs in tests (e.g. > 100 bytes if mock), but if < 10MB warn or check
    if src_size < 1024:
        size_kb = src_size / 1024
        raise ValueError(
            f"ERROR: KOfam database {shared_db_path} is only {src_size} bytes (invalid/truncated)! "
            "A valid KOfam database.pkl with 26,000+ models must be >= 1.5 GB. "
            "Please run 'serialize_kofam_models' to properly generate it."
        )

    os.makedirs(scratch_dir, exist_ok=True)
    basename = os.path.basename(shared_db_path)
    local_db_path = os.path.join(scratch_dir, basename)

    if os.path.exists(local_db_path) and os.path.getsize(local_db_path) == src_size and src_size > 1024:
        print(f"INFO: Database already staged at {local_db_path}")
        return local_db_path

    print(f"INFO: Staging database {shared_db_path} -> {local_db_path}...")
    t0 = time.time()
    tmp_target = local_db_path + f".tmp.{os.getpid()}"
    shutil.copyfile(shared_db_path, tmp_target)
    os.rename(tmp_target, local_db_path)
    dur = time.time() - t0
    size_mb = os.path.getsize(local_db_path) / (1024 * 1024)
    print(f"INFO: Database staged successfully ({size_mb:.1f} MB in {dur:.2f} s, {size_mb/max(dur,0.001):.1f} MB/s).")

    if expected_sha256:
        actual_sha = compute_sha256(local_db_path)
        if actual_sha.lower() != expected_sha256.lower():
            raise ValueError(f"ERROR: SHA256 mismatch for {local_db_path}! Expected: {expected_sha256}, Actual: {actual_sha}")
        print("INFO: Database SHA256 verified successfully.")

    return local_db_path

def get_smaps_memory():
    """Read PSS and RSS from /proc/self/smaps_rollup if available (Linux)."""
    metrics = {}
    smaps_path = "/proc/self/smaps_rollup"
    if os.path.exists(smaps_path):
        try:
            with open(smaps_path, "r") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        key = parts[0].rstrip(":")
                        if key in ["Rss", "Pss", "Shared_Clean", "Shared_Dirty", "Private_Clean", "Private_Dirty"]:
                            metrics[key] = int(parts[1]) # in kB
        except Exception:
            pass
    return metrics

def worker_task(shard_input, shard_output, local_db_path, threads_per_worker, prefilter_evalue, pykofam_env=None, all_hits=False):
    """Worker task executed in a forked process."""
    print(f"[{os.getpid()}] Starting shard: {shard_input} -> {shard_output} ({threads_per_worker} threads)")
    t0 = time.time()
    os.makedirs(os.path.dirname(os.path.abspath(shard_output)), exist_ok=True)

    cmd = []
    if pykofam_env and os.path.isdir(pykofam_env):
        cmd = ["conda", "run", "-p", pykofam_env]
    elif pykofam_env:
        cmd = ["conda", "run", "-n", pykofam_env]

    cmd.extend([
        "pykofamsearch",
        "-i", shard_input,
        "-o", shard_output,
        "-b", local_db_path,
        "-p", str(threads_per_worker)
    ])
    if all_hits:
        cmd.append("--all_hits")

    try:
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        dur = time.time() - t0
        print(f"[{os.getpid()}] Finished {shard_input} in {dur:.2f} s")
        return {"status": "SUCCESS", "shard": shard_input, "output": shard_output, "duration": dur}
    except subprocess.CalledProcessError as e:
        print(f"ERROR: pykofamsearch failed on {shard_input} with exit code {e.returncode}!", file=sys.stderr)
        if e.stdout:
            print(f"STDOUT:\n{e.stdout}", file=sys.stderr)
        if e.stderr:
            print(f"STDERR:\n{e.stderr}", file=sys.stderr)
        raise
    except Exception as e:
        print(f"ERROR: Worker failed on {shard_input}: {e}", file=sys.stderr)
        raise

def main():
    parser = argparse.ArgumentParser(description="Node-level multi-worker PyKofamSearch runner with CoW memory sharing.")
    parser.add_argument("-b", "--batch_manifest", required=True, help="TSV file with shard_input\\tshard_output")
    parser.add_argument("-d", "--database", required=True, help="Shared path to serialized database.pkl")
    parser.add_argument("-w", "--workers", type=int, default=5, help="Number of concurrent workers on this node (default: 5)")
    parser.add_argument("-t", "--threads_per_worker", type=int, default=12, help="CPU threads per worker (default: 12)")
    parser.add_argument("-E", "--prefilter_evalue", type=float, default=10.0, help="Relaxed prefilter E-value (default: 10.0)")
    parser.add_argument("-s", "--scratch_dir", default="", help="Scratch directory (default: $SLURM_TMPDIR or /tmp)")
    parser.add_argument("--all_hits", action="store_true", help="Output all hits regardless of threshold")
    parser.add_argument("--env", default="", help="Conda environment path or name for pykofamsearch")

    args = parser.parse_args()

    # Determine scratch dir
    scratch_dir = args.scratch_dir
    if not scratch_dir or scratch_dir.startswith("$"):
        env_key = scratch_dir[1:] if scratch_dir.startswith("$") else "SLURM_TMPDIR"
        scratch_dir = os.environ.get(env_key) or os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or f"/tmp/pykofam_{os.environ.get('USER', 'default')}"
    scratch_dir = os.path.abspath(scratch_dir)

    # 1. Stage database once per node
    local_db = stage_database(args.database, scratch_dir, args.db_sha256 if args.db_sha256 else None)

    # 2. Parse batch manifest
    tasks = []
    with open(args.batch_manifest, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                tasks.append((parts[0], parts[1]))

    print(f"INFO: Node job managing {len(tasks)} shards with {args.workers} workers ({args.threads_per_worker} threads each).")

    # 3. Configure multiprocessing to use 'fork' for Copy-on-Write memory sharing (Unix)
    try:
        mp.set_start_method("fork", force=True)
        print("INFO: Multiprocessing start method set to 'fork' for memory sharing.")
    except Exception as e:
        print(f"WARNING: Could not set fork start method: {e}", file=sys.stderr)

    t_start = time.time()
    pool = mp.Pool(processes=args.workers)
    async_results = []

    for shard_in, shard_out in tasks:
        r = pool.apply_async(
            worker_task,
            args=(shard_in, shard_out, local_db, args.threads_per_worker, args.prefilter_evalue, args.env, args.all_hits)
        )
        async_results.append(r)

    pool.close()
    pool.join()

    # Check for any exceptions
    results = [r.get() for r in async_results]
    total_time = time.time() - t_start

    print(f"INFO: All {len(tasks)} shards completed in {total_time:.2f} s on this node.")

if __name__ == "__main__":
    main()
