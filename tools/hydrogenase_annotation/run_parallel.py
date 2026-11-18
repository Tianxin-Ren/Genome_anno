#!/usr/bin/env python3
"""Parallel server batch runner and cohort aggregator for hydrogenase_annotation."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
SAMPLE_RE = re.compile(r"[^A-Za-z0-9_.-]+")

ACTIVE_PROCS: set[subprocess.Popen] = set()
ACTIVE_LOCK = threading.Lock()
STOP_EVENT = threading.Event()


def sanitize_sample_name(name: str) -> str:
    for suf in FASTA_SUFFIXES:
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    cleaned = SAMPLE_RE.sub("_", name).strip("._-")
    return cleaned or "sample"


def discover_inputs(
    input_dir: Path | None,
    manifest: Path | None,
    glob_pattern: str | None,
) -> list[dict[str, str]]:
    items: list[dict[str, str]] = []
    seen_samples: set[str] = set()

    if input_dir is not None:
        if not input_dir.is_dir():
            raise ValueError(f"--input-dir does not exist or is not a directory: {input_dir}")
        candidates = sorted(input_dir.glob(glob_pattern)) if glob_pattern else sorted(input_dir.iterdir())
        for path in candidates:
            if not path.is_file():
                continue
            if not glob_pattern and not any(path.name.endswith(s) for s in FASTA_SUFFIXES):
                continue
            sample = sanitize_sample_name(path.name)
            if sample in seen_samples:
                raise ValueError(f"Duplicate sample name {sample!r} from {path}")
            seen_samples.add(sample)
            items.append({"sample": sample, "path": str(path.resolve())})
    elif manifest is not None:
        if not manifest.is_file():
            raise ValueError(f"--manifest does not exist: {manifest}")
        if manifest.suffix.lower() == ".json":
            raw = json.loads(manifest.read_text(encoding="utf-8"))
            if not isinstance(raw, list):
                raise ValueError("JSON manifest must be a list of {sample, path} objects")
            for entry in raw:
                sample = str(entry["sample"]).strip()
                path = Path(entry["path"]).expanduser().resolve()
                if not sample or sample in seen_samples:
                    raise ValueError(f"Invalid or duplicate sample name in manifest: {sample!r}")
                if not path.is_file():
                    raise ValueError(f"Missing genome file for {sample}: {path}")
                seen_samples.add(sample)
                items.append({"sample": sample, "path": str(path)})
        else:
            with manifest.open("r", encoding="utf-8-sig", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    sample = (row.get("sample") or "").strip()
                    raw_path = (row.get("path") or "").strip()
                    raw_prot = (row.get("proteins") or "").strip()
                    raw_gff = (row.get("gff") or "").strip()
                    if not sample:
                        raise ValueError("TSV manifest requires 'sample' column")
                    if sample in seen_samples:
                        raise ValueError(f"Duplicate sample name in manifest: {sample!r}")
                    seen_samples.add(sample)
                    if raw_prot and raw_gff:
                        p_prot = Path(raw_prot).expanduser().resolve()
                        p_gff = Path(raw_gff).expanduser().resolve()
                        if not p_prot.is_file() or not p_gff.is_file():
                            raise ValueError(f"Missing proteins/gff for {sample}: {p_prot}, {p_gff}")
                        entry = {"sample": sample, "proteins": str(p_prot), "gff": str(p_gff)}
                        if raw_path:
                            p_dna = Path(raw_path).expanduser().resolve()
                            if not p_dna.is_file():
                                raise ValueError(f"Missing genome file for {sample}: {p_dna}")
                            entry["path"] = str(p_dna)
                        items.append(entry)
                    elif raw_path:
                        path = Path(raw_path).expanduser().resolve()
                        if not path.is_file():
                            raise ValueError(f"Missing genome file for {sample}: {path}")
                        items.append({"sample": sample, "path": str(path)})
                    else:
                        raise ValueError("TSV manifest requires either 'path' or ('proteins' and 'gff') columns")
    else:
        raise ValueError("Provide either --input-dir or --manifest")

    if not items:
        raise ValueError("No input genome files found")
    return items


def is_sample_complete(sample_dir: Path) -> bool:
    wf = sample_dir / "workflow.json"
    ann = sample_dir / "functions" / "functional_annotations.tsv"
    if not wf.is_file() or not ann.is_file():
        return False
    try:
        data = json.loads(wf.read_text(encoding="utf-8"))
        return data.get("status") == "complete"
    except (OSError, json.JSONDecodeError):
        return False


def resolve_python_and_env() -> tuple[str, dict[str, str]]:
    local_py = ROOT / ".env/bin/python"
    python_exe = str(local_py) if local_py.is_file() else sys.executable
    env = dict(os.environ)
    bin_dir = str(Path(python_exe).parent)
    env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    # Prevent BLAS/OpenMP thread oversubscription when running many workers concurrently.
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env.setdefault(var, "1")
    return python_exe, env


def prepare_input_file(source: Path, tmp_dir: Path, sample: str, suffix: str) -> tuple[Path, bool]:
    if not source.name.endswith(".gz"):
        return source, False
    tmp_dir.mkdir(parents=True, exist_ok=True)
    dest = tmp_dir / f"{sample}{suffix}"
    with gzip.open(source, "rb") as src_fh, dest.open("wb") as dst_fh:
        shutil.copyfileobj(src_fh, dst_fh)
    return dest, True


def run_one_sample(
    item: dict[str, str],
    out_root: Path,
    args: argparse.Namespace,
    python_exe: str,
    base_env: dict[str, str],
) -> dict[str, object]:
    sample = item["sample"]
    samples_dir = out_root / "samples"
    sample_dir = samples_dir / sample
    log_path = out_root / "logs" / f"{sample}.log"
    tmp_dir = out_root / "_tmp_inputs"

    if STOP_EVENT.is_set():
        return {"sample": sample, "status": "cancelled", "elapsed_s": 0.0, "error": "interrupted", "run_dir": str(sample_dir)}

    if args.resume and is_sample_complete(sample_dir):
        return {"sample": sample, "status": "complete", "elapsed_s": 0.0, "error": "", "run_dir": str(sample_dir), "skipped": True}

    if sample_dir.exists():
        if args.resume or args.overwrite_failed:
            ts = time.strftime("%Y%m%d_%H%M%S")
            archived = samples_dir / f"{sample}.incomplete_{ts}"
            shutil.move(str(sample_dir), str(archived))
        else:
            return {
                "sample": sample,
                "status": "failed",
                "elapsed_s": 0.0,
                "error": f"output directory exists (use --resume): {sample_dir}",
                "run_dir": str(sample_dir),
            }

    t0 = time.monotonic()
    tmp_files: list[Path] = []
    if "proteins" in item and "gff" in item:
        faa_path, is_tmp1 = prepare_input_file(Path(item["proteins"]), tmp_dir, sample, ".faa")
        gff_path, is_tmp2 = prepare_input_file(Path(item["gff"]), tmp_dir, sample, ".gff")
        if is_tmp1:
            tmp_files.append(faa_path)
        if is_tmp2:
            tmp_files.append(gff_path)
        input_args = ["--proteins", str(faa_path), "--gff", str(gff_path)]
        if "path" in item:
            fasta_path, is_tmp3 = prepare_input_file(Path(item["path"]), tmp_dir, sample, ".fna")
            if is_tmp3:
                tmp_files.append(fasta_path)
            input_args += ["--genome", str(fasta_path)]
    else:
        fasta_path, is_tmp = prepare_input_file(Path(item["path"]), tmp_dir, sample, ".fna")
        if is_tmp:
            tmp_files.append(fasta_path)
        input_args = ["--genome", str(fasta_path)]

    cmd = [
        python_exe,
        "-u",
        str(ROOT / "run_hydrogenase.py"),
        *input_args,
        "--mode",
        args.mode,
        "--out",
        str(sample_dir),
        "--threads",
        str(args.threads_per_worker),
        "--translation-table",
        str(args.translation_table),
        "--topology",
        args.topology,
    ]
    log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = None
    try:
        with log_path.open("w", encoding="utf-8") as log_fh:
            proc = subprocess.Popen(
                cmd,
                cwd=ROOT,
                env=base_env,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                start_new_session=True,
                text=True,
            )
            with ACTIVE_LOCK:
                ACTIVE_PROCS.add(proc)
            rc = proc.wait()
    except Exception as exc:
        elapsed = round(time.monotonic() - t0, 2)
        return {"sample": sample, "status": "failed", "elapsed_s": elapsed, "error": repr(exc), "run_dir": str(sample_dir)}
    finally:
        if proc is not None:
            with ACTIVE_LOCK:
                ACTIVE_PROCS.discard(proc)
        for tf in tmp_files:
            if tf.exists():
                tf.unlink(missing_ok=True)

    elapsed = round(time.monotonic() - t0, 2)
    if rc == 0 and is_sample_complete(sample_dir):
        return {"sample": sample, "status": "complete", "elapsed_s": elapsed, "error": "", "run_dir": str(sample_dir), "skipped": False}
    return {"sample": sample, "status": "failed", "elapsed_s": elapsed, "error": f"exit_code={rc}; see {log_path}", "run_dir": str(sample_dir)}


def load_rule_class_ids() -> list[str]:
    rules_path = ROOT / "config/function_rules_34.json"
    if not rules_path.is_file():
        return []
    data = json.loads(rules_path.read_text(encoding="utf-8"))
    return [c["id"] for c in data.get("classes", []) if c.get("id")]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv_atomic(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp_name, path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def aggregate_cohort(out_root: Path, run_records: dict[str, dict[str, object]] | None = None) -> dict[str, object]:
    samples_dir = out_root / "samples"
    class_ids = load_rule_class_ids()
    summary_rows: list[dict[str, object]] = []
    merged_annotations: list[dict[str, object]] = []
    matrix_rows: list[dict[str, object]] = []

    sample_dirs = sorted(
        p for p in samples_dir.iterdir()
        if p.is_dir() and ".incomplete_" not in p.name
    ) if samples_dir.is_dir() else []

    rule_versions: set[str] = set()
    for sdir in sample_dirs:
        sample = sdir.name
        rec = (run_records or {}).get(sample, {})
        complete = is_sample_complete(sdir)
        wf = {}
        if (sdir / "workflow.json").is_file():
            try:
                wf = json.loads((sdir / "workflow.json").read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                wf = {}

        if not complete:
            summary_rows.append({
                "sample": sample,
                "status": rec.get("status", wf.get("status", "incomplete")),
                "elapsed_s": rec.get("elapsed_s", ""),
                "total_candidates": wf.get("candidates", ""),
                "core_anchors": wf.get("core_anchors", ""),
                "predicted_classes": "",
                "class_counts_json": "{}",
                "unresolved_core_count": "",
                "error": rec.get("error", wf.get("error", "")),
                "run_dir": str(sdir),
            })
            continue

        if wf.get("function_rules_version"):
            rule_versions.add(str(wf["function_rules_version"]))

        genes = {r["gene_id"]: r for r in read_tsv(sdir / "annotation/genes.tsv") if r.get("gene_id")}
        rule_preds = {r["candidate_id"]: r for r in read_tsv(sdir / "functions/function_rules/predictions.tsv") if r.get("candidate_id")}
        ann_rows = read_tsv(sdir / "functions/functional_annotations.tsv")

        class_counter: Counter[str] = Counter()
        unresolved_core = 0

        for r in ann_rows:
            cid = r.get("candidate_id", "")
            raw_fids = r.get("function_class_ids") or "[]"
            try:
                fids = json.loads(raw_fids) if isinstance(raw_fids, str) else list(raw_fids)
            except json.JSONDecodeError:
                fids = []
            id_status = r.get("identity_status", "")
            for fid in fids:
                class_counter[fid] += 1
            if not fids and id_status in ("calibrated_high", "sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical"):
                unresolved_core += 1

            # Include non-trivial candidates in cohort_functional_annotations.tsv
            if fids or id_status not in ("unresolved", "potential_homolog", ""):
                g = genes.get(cid, {})
                rp = rule_preds.get(cid, {})
                merged_annotations.append({
                    "sample": sample,
                    "candidate_id": cid,
                    "contig": g.get("contig", ""),
                    "start": g.get("start", ""),
                    "end": g.get("end", ""),
                    "strand": g.get("strand", ""),
                    "identity_status": id_status,
                    "function_class_ids": json.dumps(fids, ensure_ascii=False),
                    "function_class_prediction": r.get("function_class_prediction", ""),
                    "function_rule_confidence": r.get("function_rule_confidence", ""),
                    "function_rule_status": r.get("function_rule_status", ""),
                    "preferred_function": r.get("preferred_function", ""),
                    "blockers": rp.get("blockers", ""),
                })

        summary_rows.append({
            "sample": sample,
            "status": "complete",
            "elapsed_s": rec.get("elapsed_s", ""),
            "total_candidates": wf.get("candidates", len(ann_rows)),
            "core_anchors": wf.get("core_anchors", ""),
            "predicted_classes": ";".join(sorted(class_counter)),
            "class_counts_json": json.dumps(dict(sorted(class_counter.items())), ensure_ascii=False),
            "unresolved_core_count": unresolved_core,
            "error": "",
            "run_dir": str(sdir),
        })

        mrow: dict[str, object] = {"sample": sample}
        all_cols = class_ids or sorted(class_counter)
        for cid_cls in all_cols:
            mrow[cid_cls] = class_counter.get(cid_cls, 0)
        for extra in sorted(set(class_counter) - set(all_cols)):
            mrow[extra] = class_counter[extra]
        matrix_rows.append(mrow)

    write_tsv_atomic(
        out_root / "cohort_summary.tsv",
        ["sample", "status", "elapsed_s", "total_candidates", "core_anchors", "predicted_classes", "class_counts_json", "unresolved_core_count", "error", "run_dir"],
        summary_rows,
    )
    write_tsv_atomic(
        out_root / "cohort_functional_annotations.tsv",
        ["sample", "candidate_id", "contig", "start", "end", "strand", "identity_status", "function_class_ids", "function_class_prediction", "function_rule_confidence", "function_rule_status", "preferred_function", "blockers"],
        merged_annotations,
    )
    matrix_cols = ["sample"] + (class_ids or sorted({k for m in matrix_rows for k in m if k != "sample"}))
    extra_cols = sorted({k for m in matrix_rows for k in m} - set(matrix_cols))
    write_tsv_atomic(out_root / "cohort_class_matrix.tsv", matrix_cols + extra_cols, matrix_rows)

    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "total_samples": len(summary_rows),
        "complete_samples": sum(1 for r in summary_rows if r["status"] == "complete"),
        "failed_samples": sum(1 for r in summary_rows if r["status"] != "complete"),
        "merged_annotation_rows": len(merged_annotations),
        "function_rules_versions": sorted(rule_versions),
    }
    (out_root / "cohort_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def handle_signal(signum, _frame) -> None:
    STOP_EVENT.set()
    with ACTIVE_LOCK:
        procs = list(ACTIVE_PROCS)
    for proc in procs:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    src = parser.add_mutually_exclusive_group()
    src.add_argument("--input-dir", type=Path, help="Directory of genome FASTA files (.fna/.fa/.fasta/.gz)")
    src.add_argument("--manifest", type=Path, help="JSON or TSV manifest with sample and path columns")
    parser.add_argument("--glob", type=str, default=None, help="Optional glob filter inside --input-dir (e.g. '*.fna')")
    parser.add_argument("--out", type=Path, required=True, help="Output directory for parallel runs and cohort summary")
    parser.add_argument("--workers", type=int, default=1, help="Number of genomes to run concurrently (e.g. 40)")
    parser.add_argument("--threads-per-worker", type=int, default=1, help="Threads passed to each genome run")
    parser.add_argument("--mode", choices=["genome", "assembly", "bin", "transcript"], default="genome")
    parser.add_argument("--topology", choices=["auto", "linear", "circular"], default="auto")
    parser.add_argument("--translation-table", type=int, default=11)
    parser.add_argument("--resume", action="store_true", help="Skip already-complete genomes and retry incomplete ones")
    parser.add_argument("--overwrite-failed", action="store_true", help="Archive existing incomplete sample dirs before running")
    parser.add_argument("--aggregate-only", action="store_true", help="Only aggregate existing results under --out")
    args = parser.parse_args(argv)

    if args.workers < 1 or args.threads_per_worker < 1:
        parser.error("--workers and --threads-per-worker must be >= 1")

    out_root = args.out.resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    if args.aggregate_only:
        manifest = aggregate_cohort(out_root)
        print(json.dumps(manifest, ensure_ascii=False))
        return 0

    items = discover_inputs(args.input_dir, args.manifest, args.glob)
    python_exe, base_env = resolve_python_and_env()

    old_sigint = signal.signal(signal.SIGINT, handle_signal)
    old_sigterm = signal.signal(signal.SIGTERM, handle_signal)
    run_records: dict[str, dict[str, object]] = {}
    total = len(items)
    completed = 0

    print(
        f"Starting parallel run: samples={total}, workers={args.workers}, "
        f"threads_per_worker={args.threads_per_worker}, mode={args.mode}",
        flush=True,
    )
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            future_map = {
                pool.submit(run_one_sample, item, out_root, args, python_exe, base_env): item["sample"]
                for item in items
            }
            for future in concurrent.futures.as_completed(future_map):
                sample = future_map[future]
                res = future.result()
                run_records[sample] = res
                completed += 1
                tag = "SKIP" if res.get("skipped") else res["status"].upper()
                print(
                    f"[{completed}/{total}] {tag} {sample} ({res.get('elapsed_s', 0)}s)"
                    + (f" -> {res['error']}" if res.get("error") else ""),
                    flush=True,
                )
    finally:
        signal.signal(signal.SIGINT, old_sigint)
        signal.signal(signal.SIGTERM, old_sigterm)

    manifest = aggregate_cohort(out_root, run_records=run_records)
    print("COHORT SUMMARY:", json.dumps(manifest, ensure_ascii=False), flush=True)
    return 0 if manifest["failed_samples"] == 0 and not STOP_EVENT.is_set() else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        raise SystemExit(2)
