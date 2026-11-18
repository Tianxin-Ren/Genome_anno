#!/usr/bin/env python3
"""
benchmark_cog_recovery.py
Tier 1: Official NCBI COG Member Recovery Benchmark.
Compares predicted COG footprints against NCBI official ground-truth footprints.

Calculates:
1. Exact COG recovery rate and Top-1 recovery rate.
2. Footprint coordinate Intersection-over-Union (IoU).
3. Coordinate start and end delta distributions (|start_pred - start_off|, |end_pred - end_off|).
4. Match code (0-3) concordance matrix.
"""

import argparse
import gzip
import json
import os
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def calculate_interval_iou(s1, e1, s2, e2):
    """Calculate 1D Intersection over Union for coordinate spans."""
    intersection = max(0, min(e1, e2) - max(s1, s2) + 1)
    union = max(e1, e2) - min(s1, s2) + 1
    return intersection / max(1, union)

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark predicted COG footprints against official ground truth."
    )
    parser.add_argument("-p", "--predicted_master", required=True, help="protein_COG_master.tsv.gz")
    parser.add_argument("-g", "--ground_truth", required=True, help="Official ground truth TSV (.tsv or .tsv.gz)")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for benchmark report")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Load Ground Truth
    # Columns expected: protein_id, cog_id, qstart, qend, match_code (optional)
    print(f"INFO: Loading ground truth from {args.ground_truth}...")
    ground_truth = {}
    with open_file(args.ground_truth, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name.lower(): i for i, name in enumerate(header)}
        pid_k = "protein_id" if "protein_id" in col_idx else header[0]
        cid_k = "cog_id" if "cog_id" in col_idx else header[1]

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            pid = parts[col_idx.get("protein_id", 0)].strip()
            cid = parts[col_idx.get("cog_id", 1)].strip()
            s = int(parts[col_idx["qstart"]]) if "qstart" in col_idx else None
            e = int(parts[col_idx["qend"]]) if "qend" in col_idx else None
            m_code = int(parts[col_idx["match_code"]]) if "match_code" in col_idx else None

            ground_truth.setdefault(pid, []).append({
                "cog_id": cid,
                "qstart": s,
                "qend": e,
                "match_code": m_code
            })

    # 2. Load Predictions
    print(f"INFO: Loading predictions from {args.predicted_master}...")
    predictions = {}
    with open_file(args.predicted_master, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            pid = parts[col_idx["protein_id"]].strip()
            cid = parts[col_idx["cog_id"]].strip()
            s = int(parts[col_idx["qstart"]].strip())
            e = int(parts[col_idx["qend"]].strip())
            m_code = int(parts[col_idx["cog_match_code"]].strip()) if "cog_match_code" in col_idx else None
            rank = int(parts[col_idx["footprint_rank"]].strip()) if "footprint_rank" in col_idx else 1

            predictions.setdefault(pid, []).append({
                "cog_id": cid,
                "qstart": s,
                "qend": e,
                "match_code": m_code,
                "rank": rank
            })

    # 3. Compute Metrics on Common Test Set
    eval_pids = set(ground_truth.keys()) & set(predictions.keys())
    print(f"INFO: Evaluating {len(eval_pids)} proteins common to ground truth and predictions...")

    exact_matches = 0
    top1_matches = 0
    iou_scores = []
    start_deltas = []
    end_deltas = []
    concordant_match_codes = 0
    total_match_code_evals = 0

    for pid in eval_pids:
        gt_fps = ground_truth[pid]
        pred_fps = predictions[pid]

        gt_cogs = {fp["cog_id"] for fp in gt_fps}
        pred_cogs = {fp["cog_id"] for fp in pred_fps}

        if gt_cogs & pred_cogs:
            exact_matches += 1

        top1_pred = pred_fps[0]["cog_id"]
        if top1_pred in gt_cogs:
            top1_matches += 1

        # Match footprints for IoU and coordinate deltas
        for gt_fp in gt_fps:
            if gt_fp["qstart"] is None or gt_fp["qend"] is None:
                continue
            best_iou = 0.0
            best_s_delta = 999
            best_e_delta = 999
            for p_fp in pred_fps:
                iou = calculate_interval_iou(gt_fp["qstart"], gt_fp["qend"], p_fp["qstart"], p_fp["qend"])
                if iou > best_iou:
                    best_iou = iou
                    best_s_delta = abs(p_fp["qstart"] - gt_fp["qstart"])
                    best_e_delta = abs(p_fp["qend"] - gt_fp["qend"])
                    if gt_fp["match_code"] is not None and p_fp["match_code"] is not None:
                        total_match_code_evals += 1
                        if gt_fp["match_code"] == p_fp["match_code"]:
                            concordant_match_codes += 1

            if best_iou > 0:
                iou_scores.append(best_iou)
                start_deltas.append(best_s_delta)
                end_deltas.append(best_e_delta)

    n_eval = max(1, len(eval_pids))
    exact_rate = round(exact_matches / n_eval, 4)
    top1_rate = round(top1_matches / n_eval, 4)
    mean_iou = round(sum(iou_scores) / max(1, len(iou_scores)), 4) if iou_scores else 0.0
    mean_start_delta = round(sum(start_deltas) / max(1, len(start_deltas)), 2) if start_deltas else 0.0
    mean_end_delta = round(sum(end_deltas) / max(1, len(end_deltas)), 2) if end_deltas else 0.0
    code_concordance = round(concordant_match_codes / max(1, total_match_code_evals), 4) if total_match_code_evals > 0 else 1.0

    report = {
        "n_evaluated_proteins": len(eval_pids),
        "exact_cog_recovery_rate": exact_rate,
        "top1_recovery_rate": top1_rate,
        "mean_footprint_iou": mean_iou,
        "mean_start_coordinate_delta_aa": mean_start_delta,
        "mean_end_coordinate_delta_aa": mean_end_delta,
        "match_code_concordance_rate": code_concordance,
        "n_footprints_evaluated": len(iou_scores)
    }

    report_json_path = os.path.join(args.out_dir, "cog_recovery_benchmark_report.json")
    with open(report_json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    summary_tsv_path = os.path.join(args.out_dir, "cog_recovery_benchmark_summary.tsv")
    with open(summary_tsv_path, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        for k, v in report.items():
            f.write(f"{k}\t{v}\n")

    print(f"SUCCESS: Benchmark complete. Exact recovery = {exact_rate*100:.2f}%, Mean IoU = {mean_iou:.4f}")

if __name__ == "__main__":
    main()
