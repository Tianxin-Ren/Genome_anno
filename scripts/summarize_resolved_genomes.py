#!/usr/bin/env python3
"""Create one genome-level row from resolved HMM annotations and MAG metadata."""

import argparse
import csv
from multiprocessing import Pool
import re
import time
from pathlib import Path

import pandas as pd


EVIDENCE_LEVELS = ("high", "medium", "candidate", "rejected")
_TARGET_LOOKUP = {}
_ALLOW_LEGACY = False


def read_target_map(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"target_id", "target_name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"target map missing columns: {', '.join(sorted(missing))}")
    df["target_id"] = df["target_id"].str.strip()
    df["target_name"] = df["target_name"].str.strip()
    df = df[df["target_id"] != ""].copy()
    return df.drop_duplicates("target_id")


def read_mag_info(path, genome_column):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if genome_column not in df.columns:
        alternatives = [c for c in ("bin", "genome", "Genome_Name", "user_genome") if c in df.columns]
        if not alternatives:
            raise ValueError(f"MAG info table has no genome column '{genome_column}'")
        genome_column = alternatives[0]
    df[genome_column] = df[genome_column].astype(str).str.strip()
    df = df[df[genome_column] != ""].copy()
    return df, genome_column


def infer_genome_id(path):
    """Support both MAG/MAG_resolved.tsv and flat MAG_resolved.tsv layouts."""
    name = path.name
    suffix = "_advanced_hmm_resolved.tsv"
    if name.endswith(suffix):
        stem_id = name[: -len(suffix)]
        if path.parent.name == stem_id:
            return stem_id
        return stem_id
    return path.stem


def get_target_ids(row, include_legacy=False):
    values = []
    columns = ("Target_KO", "Target_Annotation")
    if include_legacy:
        columns += ("Score_Annotation", "Evalue_Annotation")
    for column in columns:
        value = str(row.get(column, "")).strip()
        if value and value != "-":
            values.append(value)
            values.extend(re.split(r"[|,;\s]+", value))
    return {v for v in values if v and v != "-"}


def quality_category(row):
    try:
        completeness = float(row.get("completeness", ""))
        contamination = float(row.get("contamination", ""))
    except (TypeError, ValueError):
        return "unknown"
    if completeness >= 90 and contamination <= 5:
        return "high_quality"
    if completeness >= 70 and contamination <= 10:
        return "medium_quality"
    return "low_quality"


def init_worker(target_lookup):
    global _TARGET_LOOKUP
    _TARGET_LOOKUP = target_lookup


def init_worker_with_options(target_lookup, allow_legacy):
    global _TARGET_LOOKUP, _ALLOW_LEGACY
    _TARGET_LOOKUP = target_lookup
    _ALLOW_LEGACY = allow_legacy


def summarize_one_file(path, target_lookup=None, allow_legacy=None):
    """Stream one TSV; avoid constructing a pandas DataFrame per genome."""
    started = time.monotonic()
    target_lookup = target_lookup if target_lookup is not None else _TARGET_LOOKUP
    if allow_legacy is None:
        allow_legacy = _ALLOW_LEGACY
    target_ids = set(target_lookup)
    matched = set()
    rows_read = 0
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        legacy_file = "Evidence_Level" not in (reader.fieldnames or [])
        if not reader.fieldnames or (legacy_file and not allow_legacy):
            raise ValueError(f"{path} lacks required column Evidence_Level")
        for row in reader:
            rows_read += 1
            hit_targets = sorted(get_target_ids(row, include_legacy=legacy_file) & target_ids)
            if not hit_targets:
                continue
            evidence = "candidate" if legacy_file else str(row.get("Evidence_Level", "candidate")).strip().lower()
            if evidence not in EVIDENCE_LEVELS:
                evidence = "candidate"
            protein_id = str(row.get("Protein_ID", row.get("Protein", ""))).strip()
            for target_id in hit_targets:
                matched.add((protein_id, target_id, evidence))

    by_level = {level: 0 for level in EVIDENCE_LEVELS}
    by_target = {}
    for _, target_id, evidence in matched:
        by_level[evidence] += 1
        by_target.setdefault(target_id, {level: 0 for level in EVIDENCE_LEVELS})
        by_target[target_id][evidence] += 1
    accepted = by_level["high"] + by_level["medium"] + by_level["candidate"]
    return infer_genome_id(path), {
            "resolved_file": str(path),
            "target_high_count": by_level["high"],
            "target_medium_count": by_level["medium"],
            "target_candidate_count": by_level["candidate"],
            "target_rejected_count": by_level["rejected"],
            "target_total_count": accepted,
            "target_all_records_count": len(matched),
            "target_summary": ", ".join(
                f"{target_lookup[t]}:{sum(by_target.get(t, {}).values())}" for t in sorted(target_ids)
                if sum(by_target.get(t, {}).values())
            ) or "None",
            "_debug_rows": rows_read,
            "_debug_seconds": time.monotonic() - started,
            "_debug_legacy": legacy_file,
        }


def summarize_resolved(resolved_dir, target_map, workers=1, progress_every=1000, debug=False, allow_legacy=False):
    target_lookup = dict(zip(target_map["target_id"], target_map["target_name"]))
    files = list(Path(resolved_dir).rglob("*_advanced_hmm_resolved.tsv"))
    print(f"Resolved files found: {len(files)}; workers: {workers}", flush=True)
    summaries = {}
    processed = 0
    matched_total = 0

    def consume(pair):
        nonlocal processed, matched_total
        genome, values = pair
        summaries[genome] = values
        processed += 1
        matched_total += values["target_all_records_count"]
        if debug and (processed <= 3 or values["_debug_seconds"] >= 10):
            print(
                f"[debug] {genome}: rows={values['_debug_rows']}, "
                f"seconds={values['_debug_seconds']:.2f}",
                flush=True,
            )
        if progress_every > 0 and (processed % progress_every == 0 or processed == len(files)):
            print(
                f"Processed {processed}/{len(files)} genomes; "
                f"target records: {matched_total}",
                flush=True,
            )

    if workers <= 1:
        for path in files:
            consume(summarize_one_file(path, target_lookup, allow_legacy))
    else:
        # TSV parsing is CPU-bound enough that threads are limited by the GIL.
        # Use processes, matching count_hydrogenases.py's fast path.
        with Pool(
            processes=workers,
            initializer=init_worker_with_options,
            initargs=(target_lookup, allow_legacy),
        ) as pool:
            chunksize = 10 if debug else 100
            for pair in pool.imap_unordered(summarize_one_file, files, chunksize=chunksize):
                consume(pair)
    return summaries


def main():
    parser = argparse.ArgumentParser(description="Summarize resolved HMM targets per genome.")
    parser.add_argument("--resolved-dir", required=True, help="Directory containing per-genome resolved TSV files")
    parser.add_argument("--mag-info-tsv", required=True, help="Genome metadata, quality and taxonomy TSV")
    parser.add_argument("--target-map", required=True, help="TSV with target_id and target_name columns")
    parser.add_argument("--out-genome", required=True, help="One-row-per-genome output TSV")
    parser.add_argument("--genome-column", default="bin", help="Genome ID column in mag-info-tsv")
    parser.add_argument("--workers", type=int, default=1,
                        help="Concurrent resolved-file readers; start with 8-16 on shared storage")
    parser.add_argument("--progress-every", type=int, default=1000,
                        help="Print progress every N genomes; 0 disables progress output")
    parser.add_argument("--debug", action="store_true",
                        help="Print phase timing, sample file timing and slow-file diagnostics")
    parser.add_argument("--allow-legacy", action="store_true",
                        help="Allow old resolved TSVs without Evidence_Level; count matching hits as candidate")
    args = parser.parse_args()

    started = time.monotonic()
    if args.debug:
        print("[debug] loading target map...", flush=True)
    target_map = read_target_map(args.target_map)
    if args.debug:
        print(f"[debug] target map loaded: {len(target_map)} targets", flush=True)
        print("[debug] loading MAG metadata...", flush=True)
    mag_info, genome_column = read_mag_info(args.mag_info_tsv, args.genome_column)
    if args.debug:
        print(f"[debug] MAG metadata loaded: {len(mag_info)} rows", flush=True)
        print("[debug] scanning resolved files...", flush=True)
    summaries = summarize_resolved(
        args.resolved_dir,
        target_map,
        max(1, args.workers),
        args.progress_every,
        args.debug,
        args.allow_legacy,
    )
    summary_columns = [
        "resolved_file", "target_high_count", "target_medium_count",
        "target_candidate_count", "target_rejected_count", "target_total_count",
        "target_all_records_count", "target_summary",
    ]
    summary_df = pd.DataFrame(
        [{genome_column: genome, **values} for genome, values in summaries.items()],
        columns=[genome_column] + summary_columns,
    )
    output = mag_info.merge(summary_df, on=genome_column, how="left")
    for col in (
        "target_high_count", "target_medium_count", "target_candidate_count",
        "target_rejected_count", "target_total_count", "target_all_records_count",
    ):
        output[col] = pd.to_numeric(output[col], errors="coerce").fillna(0).astype(int)
    output["target_summary"] = output["target_summary"].fillna("None")
    output["resolved_file"] = output["resolved_file"].fillna("")
    output["quality_category"] = output.apply(quality_category, axis=1)

    preferred = [
        genome_column, "completeness", "contamination", "quality_category", "GC",
        "classification", "warnings", "target_high_count", "target_medium_count",
        "target_candidate_count", "target_rejected_count", "target_total_count",
        "target_all_records_count", "target_summary", "resolved_file",
    ]
    ordered = [c for c in preferred if c in output.columns]
    ordered.extend(c for c in output.columns if c not in ordered)
    output = output[ordered]

    out_path = Path(args.out_genome)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(output)} genome rows to {out_path}")
    print(f"Resolved files summarized: {len(summaries)}; target profiles: {len(target_map)}")
    if args.debug:
        print(f"[debug] total seconds: {time.monotonic() - started:.2f}", flush=True)


if __name__ == "__main__":
    main()
