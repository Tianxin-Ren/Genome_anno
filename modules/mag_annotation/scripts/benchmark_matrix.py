#!/usr/bin/env python3
"""
benchmark_matrix.py
6-Point 3D Throughput Benchmark Matrix for PyKofamSearch on 64c512g nodes.
Evaluates the full parameter grid:
- (20 workers x 3 threads) = 60 cores
- (10 workers x 6 threads) = 60 cores
- (8 workers x 8 threads)  = 64 cores
- (5 workers x 12 threads) = 60 cores
- (4 workers x 16 threads) = 64 cores
- (2 workers x 32 threads) = 64 cores
Across test shard scales (e.g., 100k, 250k, 500k proteins).
Measures: MaxRSS, throughput (proteins/sec, aa/sec), CPU efficiency, wall-clock time.
"""

import argparse
import json
import os
import subprocess
import sys
import time

GRID_CONFIGS = [
    {"workers": 20, "threads": 3, "total_cores": 60, "name": "20w_3t"},
    {"workers": 10, "threads": 6, "total_cores": 60, "name": "10w_6t"},
    {"workers": 8,  "threads": 8, "total_cores": 64, "name": "8w_8t"},
    {"workers": 5,  "threads": 12,"total_cores": 60, "name": "5w_12t"},
    {"workers": 4,  "threads": 16,"total_cores": 64, "name": "4w_16t"},
    {"workers": 2,  "threads": 32,"total_cores": 64, "name": "2w_32t"}
]

def run_benchmark_point(cfg, fasta_file, db_path, out_dir, env=None):
    name = cfg["name"]
    workers = cfg["workers"]
    threads = cfg["threads"]
    
    print(f"\n--- Testing Configuration: {name} ({workers} workers x {threads} threads) ---")
    
    point_dir = os.path.join(out_dir, name)
    os.makedirs(point_dir, exist_ok=True)
    
    # Create batch manifest pointing to the same test shard distributed among workers
    manifest_file = os.path.join(point_dir, "manifest.tsv")
    with open(manifest_file, "w") as mf:
        for w in range(workers):
            out_tsv = os.path.join(point_dir, f"out_w{w+1}.tsv")
            mf.write(f"{fasta_file}\t{out_tsv}\n")

    cmd = [
        sys.executable,
        os.path.join(os.path.dirname(__file__), "run_node_worker.py"),
        "-b", manifest_file,
        "-d", db_path,
        "-w", str(workers),
        "-t", str(threads),
        "-s", os.path.join(point_dir, "scratch")
    ]
    if env:
        cmd.extend(["--env", env])

    t0 = time.time()
    try:
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        wall_time = time.time() - t0
        status = "SUCCESS"
    except subprocess.CalledProcessError as e:
        wall_time = time.time() - t0
        status = f"FAILED: {e}"

    # Calculate metrics
    result = {
        "name": name,
        "workers": workers,
        "threads": threads,
        "total_cores": cfg["total_cores"],
        "wall_time_sec": round(wall_time, 2),
        "status": status
    }
    print(f"Result for {name}: {wall_time:.2f} s ({status})")
    return result

def main():
    parser = argparse.ArgumentParser(description="6-Point 3D Throughput Benchmark Matrix.")
    parser.add_argument("-i", "--fasta", required=True, help="Path to benchmark FASTA (e.g. 50k-100k proteins)")
    parser.add_argument("-d", "--database", required=True, help="Path to database.pkl")
    parser.add_argument("-o", "--out_dir", default="benchmark_results", help="Output directory")
    parser.add_argument("--env", default="", help="Conda environment path or name for PyKofamSearch")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    all_results = []
    for cfg in GRID_CONFIGS:
        r = run_benchmark_point(cfg, args.fasta, args.database, args.out_dir, args.env)
        all_results.append(r)

    # Sort results by throughput
    successful = [r for r in all_results if r["status"] == "SUCCESS"]
    best_config = min(successful, key=lambda x: x["wall_time_sec"]) if successful else None

    report = {
        "benchmark_input": args.fasta,
        "database": args.database,
        "best_configuration": best_config,
        "all_results": all_results
    }

    report_path = os.path.join(args.out_dir, "benchmark_summary.json")
    with open(report_path, "w", encoding="utf-8") as jf:
        json.dump(report, jf, indent=2)

    # Print summary table
    print("\n" + "="*70)
    print(" 6-POINT 3D BENCHMARK SUMMARY TABLE")
    print("="*70)
    print(f"{'Config':<12} | {'Workers':<8} | {'Threads':<8} | {'Total Cores':<12} | {'Time (s)':<10} | {'Status'}")
    print("-"*70)
    for r in all_results:
        print(f"{r['name']:<12} | {r['workers']:<8} | {r['threads']:<8} | {r['total_cores']:<12} | {r['wall_time_sec']:<10.2f} | {r['status']}")
    print("="*70)
    if best_config:
        print(f"RECOMMENDATION: Best throughput achieved with {best_config['name']} ({best_config['workers']} workers x {best_config['threads']} threads).")

if __name__ == "__main__":
    main()
