#!/usr/bin/env python3
"""Plan and apply read-profiling migration to union cold storage."""

import argparse
import csv
import hashlib
import os
import re
import shutil
import tempfile
from pathlib import Path

SOURCE_ROOT = Path("/dssg/home")
DEST_ROOT = Path("/union/home")
FIELDS = ["sample", "category", "action", "logical_path", "physical_path", "union_destination", "bytes", "mtime_ns"]


def actual_root(root):
    return Path(os.path.realpath(root))


def union_destination(physical):
    return DEST_ROOT / Path(physical).relative_to(actual_root(SOURCE_ROOT))


def storage_alias(path, source, destination):
    path = str(path)
    return destination + path[len(source):] if path == source or path.startswith(source + "/") else path


def classify(path, category, sample):
    logical = Path(path)
    physical = Path(os.path.realpath(logical))
    exists = physical.is_file()
    if category == "delete_candidate":
        action = "missing"
        if logical.is_symlink():
            action = "review_symlink_intermediate"
        elif exists:
            action = "delete_intermediate" if (
                within(physical, actual_root(SOURCE_ROOT))
                or within(physical, actual_root(DEST_ROOT))
            ) else "review_other_storage"
    elif not exists:
        action = "missing"
    elif within(physical, actual_root(SOURCE_ROOT)):
        action = "move_to_union_and_link"
    elif within(physical, actual_root(DEST_ROOT)):
        action = "already_on_union"
    else:
        action = "review_other_storage"
    return {
        "sample": sample,
        "category": category,
        "action": action,
        "logical_path": str(logical),
        "physical_path": str(physical),
        "union_destination": str(union_destination(physical)) if action == "move_to_union_and_link" else "",
        "bytes": physical.stat().st_size if exists else "",
        "mtime_ns": physical.stat().st_mtime_ns if exists else "",
    }


def raw_reads(sample, store):
    """Match the sample's raw pair at the store root or one directory below."""
    found = {}
    norm_sample = re.sub(r"[-_]", "", sample).lower()
    for root in (Path(store), Path(storage_alias(store, str(SOURCE_ROOT), str(DEST_ROOT)))):
        if not root.is_dir():
            continue
        directories = [root]
        directories.extend(
            child for child in root.iterdir()
            if child.is_dir() and child.name not in {"00.QC", "05.Kraken2", "logs", "results"}
        )
        for directory in directories:
            for path in directory.iterdir():
                match = re.match(r"^(.*?)[._-]R?([12])(?:[._-].*)?\.(?:fq|fastq)\.gz$", path.name, re.I)
                if match and re.sub(r"[-_]", "", match.group(1)).lower() == norm_sample and path.is_file():
                    found.setdefault(os.path.realpath(path), path)
    return sorted(found.values())


def inventory(tsv, results_dir):
    rows = []
    seen = set()
    with open(tsv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"RUN_id", "Store_address"}.issubset(reader.fieldnames):
            raise ValueError("TSV must contain RUN_id and Store_address columns")
        for record in reader:
            sample, store = record["RUN_id"].strip(), record["Store_address"].strip()
            if not sample or not store:
                raise ValueError(f"Empty RUN_id or Store_address in {tsv}")
            store = storage_alias(store, str(DEST_ROOT), str(SOURCE_ROOT))
            raw = raw_reads(sample, store)
            for path in raw:
                rows.append(classify(path, "raw_reads", sample))
            if not raw:
                rows.append(classify(Path(store) / f"{sample}.raw_pair_NOT_FOUND", "raw_reads", sample))

            qc = Path(store) / "00.QC"
            union_qc = Path(storage_alias(store, str(SOURCE_ROOT), str(DEST_ROOT))) / "00.QC"

            def add_intermediate(name):
                candidates = {os.path.realpath(path): path for path in (qc / name, union_qc / name) if path.exists() or path.is_symlink()}
                for path in candidates.values() if candidates else (qc / name,):
                    rows.append(classify(path, "delete_candidate", sample))

            for mate in (1, 2):
                clean = qc / f"{sample}.clean.{mate}.fq.gz"
                if not clean.is_file() and (union_qc / clean.name).is_file():
                    clean = union_qc / clean.name
                rows.append(classify(clean, "clean_reads", sample))
                add_intermediate(f"{sample}.unmerged.{mate}.fq.gz")
            add_intermediate(f"{sample}.merged.fq.gz")
            for suffix in ("m", "p"):
                for extension in ("json", "html"):
                    add_intermediate(f"{sample}.fastp_{suffix}.{extension}")
            for extension in ("json", "html"):
                add_intermediate(f"{sample}.fastp.{extension}")

            sample_result = Path(results_dir) / sample
            for path in sorted(sample_result.glob("*.sylph.tsv")) + sorted(sample_result.glob("*.sylphmpa")):
                rows.append(classify(path, "sylph_result", sample))

    for pattern in ("merged_*.tsv", "profiling.*.tsv", "*.sylphmpa"):
        for path in sorted(Path(results_dir).glob(pattern)):
            if path not in seen:
                rows.append(classify(path, "sylph_result", "ALL"))
                seen.add(path)
    return rows


def within(path, root):
    path = Path(os.path.abspath(path))
    root = Path(os.path.abspath(root))
    return path == root or root in path.parents


def digest(path):
    checksum = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.digest()


def delete_name_is_expected(sample, name):
    return name in {
        f"{sample}.merged.fq.gz",
        f"{sample}.unmerged.1.fq.gz",
        f"{sample}.unmerged.2.fq.gz",
        *(f"{sample}.fastp_{part}.{ext}" for part in ("m", "p") for ext in ("json", "html")),
        *(f"{sample}.fastp.{ext}" for ext in ("json", "html")),
    }


def raw_mate(name):
    match = re.search(r"[._-]R?([12])(?:[._-].*)?\.(?:fq|fastq)\.gz$", name, re.I)
    return match.group(1) if match else None


def check_snapshot(path, row):
    stat = path.stat()
    if stat.st_size != int(row["bytes"]) or stat.st_mtime_ns != int(row["mtime_ns"]):
        raise RuntimeError(f"File changed since inventory: {path}")
    return stat


def preflight(rows):
    moves, deletions, retained = {}, {}, set()
    present = {}
    for row in rows:
        sample, category, action = row["sample"], row["category"], row["action"]
        source = Path(row["physical_path"])
        logical = Path(row["logical_path"])
        if action in ("move_to_union_and_link", "already_on_union"):
            if category not in ("raw_reads", "clean_reads", "sylph_result"):
                raise ValueError(f"Unexpected retained category: {row}")
            retained.add(source)
            present.setdefault(sample, {}).setdefault(category, 0)
            present[sample][category] += 1
            if category == "raw_reads":
                present[sample].setdefault("raw_mates", set()).add(raw_mate(logical.name))
        if action == "already_on_union":
            if not within(source, actual_root(DEST_ROOT)) or not source.is_file() or os.path.realpath(logical) != str(source):
                raise RuntimeError(f"Archived file is missing or changed: {logical}")
            check_snapshot(source, row)
        if action == "move_to_union_and_link":
            dest = Path(row["union_destination"])
            if not within(source, actual_root(SOURCE_ROOT)) or dest != union_destination(source):
                raise ValueError(f"Invalid migration paths: {source} -> {dest}")
            if dest.is_symlink():
                raise RuntimeError(f"Destination is a symlink: {dest}")
            if source.is_symlink() and dest.is_file() and os.path.samefile(source, dest):
                if not logical.is_file() or not os.path.samefile(logical, dest):
                    raise RuntimeError(f"Original path no longer reaches the archived file: {logical}")
                check_snapshot(dest, row)
                continue  # This row was completed by an earlier interrupted run.
            if not source.is_file() or os.path.realpath(logical) != str(source):
                raise RuntimeError(f"Source is missing or changed: {logical}")
            if check_snapshot(source, row).st_nlink != 1:
                raise RuntimeError(f"Source has hard links: {source}")
            if dest.exists():
                if not dest.is_file() or os.path.samefile(source, dest) or digest(source) != digest(dest):
                    raise RuntimeError(f"Conflicting destination: {dest}")
            moves[source] = (dest, row)
        elif action == "delete_intermediate":
            if category != "delete_candidate" or not (within(source, actual_root(SOURCE_ROOT)) or within(source, actual_root(DEST_ROOT))):
                raise ValueError(f"Invalid deletion path: {source}")
            if not delete_name_is_expected(sample, logical.name) or logical.is_symlink():
                raise ValueError(f"Unexpected intermediate: {logical}")
            if not logical.exists():
                continue  # Already removed by an earlier interrupted run.
            if not source.is_file() or os.path.realpath(logical) != str(source):
                raise RuntimeError(f"Intermediate is missing or changed: {logical}")
            check_snapshot(source, row)
            deletions[source] = (logical, row)
    if retained.intersection(deletions):
        raise RuntimeError("A deletion candidate is also a retained file")
    return moves, deletions, present


def move_and_link(source, dest, row):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if within(Path(os.path.realpath(dest.parent)), actual_root(SOURCE_ROOT)):
        raise RuntimeError(f"Union destination points back into dssg: {dest.parent}")
    if not dest.exists():
        with tempfile.NamedTemporaryFile(dir=dest.parent, prefix=f".{dest.name}.copy-", delete=False) as temp:
            temp_path = Path(temp.name)
            try:
                with source.open("rb") as original:
                    shutil.copyfileobj(original, temp, 8 * 1024 * 1024)
                temp.flush()
                os.fsync(temp.fileno())
            except BaseException:
                temp_path.unlink(missing_ok=True)
                raise
        try:
            shutil.copystat(source, temp_path)
            if digest(source) != digest(temp_path):
                raise RuntimeError(f"Copy verification failed: {source}")
            os.link(temp_path, dest)  # Publish without overwriting another file.
        finally:
            temp_path.unlink(missing_ok=True)
    if dest.is_symlink() or digest(source) != digest(dest):
        raise RuntimeError(f"Destination verification failed: {dest}")
    check_snapshot(source, row)
    link_tmp = source.with_name(f".{source.name}.archive-link-{os.getpid()}")
    try:
        os.symlink(dest, link_tmp)
        os.replace(link_tmp, source)
    finally:
        link_tmp.unlink(missing_ok=True)
    if not source.is_symlink() or not os.path.samefile(source, dest):
        raise RuntimeError(f"Source link verification failed: {source}")


def apply_plan(plan, report):
    if plan.resolve() == report.resolve():
        raise ValueError("Execution report must not overwrite the archive plan")
    with plan.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(FIELDS).issubset(reader.fieldnames):
            raise ValueError("Archive plan is missing required columns; regenerate it")
        rows = list(reader)
    moves, deletions, present = preflight(rows)
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["action", "sample", "path", "destination"])
        for source, (dest, row) in moves.items():
            print(f"Archiving {source} -> {dest}", flush=True)
            move_and_link(source, dest, row)
            writer.writerow(["moved_and_linked", row["sample"], source, dest])
            handle.flush()
        for source, (logical, row) in deletions.items():
            counts = present.get(row["sample"], {})
            if counts.get("raw_mates") != {"1", "2"} or counts.get("clean_reads", 0) < 2 or counts.get("sylph_result", 0) < 1:
                writer.writerow(["kept_incomplete_sample", row["sample"], source, ""])
                continue
            check_snapshot(source, row)
            print(f"Deleting intermediate {logical}", flush=True)
            logical.unlink()
            writer.writerow(["deleted_intermediate", row["sample"], logical, ""])
            handle.flush()
    print(f"Archive execution report: {report}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--apply-plan", type=Path, help="Move verified files, link them back, then delete approved intermediates")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.apply_plan:
        if args.tsv or args.results_dir or args.output:
            parser.error("--apply-plan cannot be combined with inventory options")
        apply_plan(args.apply_plan, args.report or args.apply_plan.with_name("archive_apply_report.tsv"))
        return
    if not args.tsv or not args.results_dir or not args.output:
        parser.error("--tsv, --results-dir and --output are required to make a plan")
    rows = inventory(args.tsv, args.results_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} inventory rows to {args.output}; no files were moved or deleted.")


if __name__ == "__main__":
    main()
