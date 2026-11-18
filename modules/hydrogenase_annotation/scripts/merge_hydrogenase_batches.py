#!/usr/bin/env python3
"""Merge per-batch hydrogenase_annotation outputs and optionally cross-link with mag_annotation."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        return [], []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def merge_batches(batch_dirs: list[Path], out_dir: Path) -> dict[str, object]:
    summary_fields: list[str] = []
    summary_rows: list[dict[str, str]] = []
    ann_fields: list[str] = []
    ann_rows: list[dict[str, str]] = []
    matrix_fields: list[str] = []
    matrix_rows: list[dict[str, str]] = []
    rule_versions: set[str] = set()

    for bdir in batch_dirs:
        sf, sr = read_tsv(bdir / "cohort_summary.tsv")
        if sf and not summary_fields:
            summary_fields = sf
        summary_rows.extend(sr)

        af, ar = read_tsv(bdir / "cohort_functional_annotations.tsv")
        if af and not ann_fields:
            ann_fields = af
        ann_rows.extend(ar)

        mf, mr = read_tsv(bdir / "cohort_class_matrix.tsv")
        if mf and not matrix_fields:
            matrix_fields = mf
        else:
            for col in mf:
                if col not in matrix_fields:
                    matrix_fields.append(col)
        matrix_rows.extend(mr)

        cman = bdir / "cohort_manifest.json"
        if cman.is_file():
            try:
                data = json.loads(cman.read_text(encoding="utf-8"))
                for v in data.get("function_rules_versions", []):
                    rule_versions.add(str(v))
            except (OSError, json.JSONDecodeError):
                pass

    summary_rows.sort(key=lambda r: r.get("sample", ""))
    ann_rows.sort(key=lambda r: (r.get("sample", ""), r.get("candidate_id", "")))
    matrix_rows.sort(key=lambda r: r.get("sample", ""))

    if not summary_fields:
        summary_fields = [
            "sample", "status", "elapsed_s", "total_candidates", "core_anchors",
            "predicted_classes", "class_counts_json", "unresolved_core_count", "error", "run_dir",
        ]
    if not ann_fields:
        ann_fields = [
            "sample", "candidate_id", "contig", "start", "end", "strand",
            "identity_status", "function_class_ids", "function_class_prediction",
            "function_rule_confidence", "function_rule_status", "preferred_function", "blockers",
        ]
    if not matrix_fields:
        matrix_fields = ["sample"]

    write_tsv(out_dir / "cohort_summary.tsv", summary_fields, summary_rows)
    write_tsv(out_dir / "cohort_functional_annotations.tsv", ann_fields, ann_rows)
    write_tsv(out_dir / "cohort_class_matrix.tsv", matrix_fields, matrix_rows)

    manifest = {
        "batch_count": len(batch_dirs),
        "total_genomes": len(summary_rows),
        "complete_genomes": sum(1 for r in summary_rows if r.get("status") == "complete"),
        "failed_genomes": sum(1 for r in summary_rows if r.get("status") != "complete"),
        "annotation_rows": len(ann_rows),
        "function_rules_versions": sorted(rule_versions),
    }
    (out_dir / "cohort_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def cross_link_mag_annotation(out_dir: Path, mag_sample_dir: Path) -> int:
    ann_fields, ann_rows = read_tsv(out_dir / "cohort_functional_annotations.tsv")
    if not ann_rows or not mag_sample_dir.is_dir():
        return 0

    ko_best_path = mag_sample_dir / "merged" / "protein_KO_best.tsv"
    if not ko_best_path.is_file():
        ko_best_path = mag_sample_dir / "merged" / "protein_KO_best.tsv.gz"
    if not ko_best_path.is_file():
        return 0

    # Index KEGG best hits by protein_id
    ko_by_prot: dict[str, dict[str, str]] = {}
    with open_text(ko_best_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = row.get("protein_id") or row.get("query") or ""
            if pid:
                ko_by_prot[pid] = {
                    "kegg_ko": row.get("ko") or row.get("KO") or "",
                    "kegg_score": row.get("score") or row.get("bitscore") or "",
                    "kegg_evalue": row.get("evalue") or "",
                }

    # Index mag_annotation protein metadata by (genome_id, strand) for coordinate matching
    meta_by_genome: dict[str, list[dict[str, object]]] = defaultdict(list)
    meta_dir = mag_sample_dir / "metadata"
    if meta_dir.is_dir():
        target_genomes = {r.get("sample", "") for r in ann_rows}
        for mf in meta_dir.glob("*_protein_metadata.tsv*"):
            with open_text(mf) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    gid = row.get("genome_id") or ""
                    if gid not in target_genomes:
                        continue
                    try:
                        s = int(row["start"])
                        e = int(row["end"])
                    except (KeyError, ValueError, TypeError):
                        continue
                    meta_by_genome[gid].append({
                        "protein_id": row.get("protein_id", ""),
                        "contig": row.get("contig_id") or row.get("contig") or "",
                        "start": s,
                        "end": e,
                        "strand": row.get("strand", ""),
                    })

    integrated: list[dict[str, object]] = []
    for row in ann_rows:
        gid = row.get("sample", "")
        cid = row.get("candidate_id", "")
        out_row = dict(row)
        out_row.update({"mag_protein_id": "", "kegg_ko": "", "kegg_score": "", "kegg_evalue": ""})
        if cid and cid in ko_by_prot:
            out_row["mag_protein_id"] = cid
            out_row.update(ko_by_prot[cid])
            integrated.append(out_row)
            continue

        try:
            r_start = int(row.get("start") or 0)
            r_end = int(row.get("end") or 0)
        except ValueError:
            r_start, r_end = 0, 0

        best_ov = 0
        best_pid = ""
        for m in meta_by_genome.get(gid, []):
            if m["strand"] and row.get("strand") and m["strand"] != row["strand"]:
                continue
            ov = min(r_end, int(m["end"])) - max(r_start, int(m["start"])) + 1
            if ov > best_ov and ov >= 90:
                best_ov = ov
                best_pid = str(m["protein_id"])

        if best_pid:
            out_row["mag_protein_id"] = best_pid
            if best_pid in ko_by_prot:
                out_row.update(ko_by_prot[best_pid])
        integrated.append(out_row)

    out_fields = ann_fields + ["mag_protein_id", "kegg_ko", "kegg_score", "kegg_evalue"]
    write_tsv(out_dir / "hydrogenase_kegg_integrated.tsv", out_fields, integrated)
    return len(integrated)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", nargs="+", type=Path, required=True, help="Batch output directories")
    parser.add_argument("--out", type=Path, required=True, help="Merged sample output directory")
    parser.add_argument("--mag-sample-dir", type=Path, default=None, help="Optional mag_annotation/{sample} dir for KEGG cross-link")
    args = parser.parse_args()

    out_dir = args.out.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = merge_batches([b.resolve() for b in args.batches], out_dir)
    if args.mag_sample_dir is not None and args.mag_sample_dir.is_dir():
        linked = cross_link_mag_annotation(out_dir, args.mag_sample_dir.resolve())
        manifest["kegg_integrated_rows"] = linked
        (out_dir / "cohort_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
