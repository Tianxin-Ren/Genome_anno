#!/usr/bin/env python3
"""Package hydrogenase_annotation into a clean, server-ready tarball."""

from __future__ import annotations

import argparse
import hashlib
import tarfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INCLUDE_FILES = [
    "run_hydrogenase.py",
    "run_parallel.py",
    "run_batch.py",
    "setup.sh",
    "SERVER_README.md",
]

INCLUDE_DIRS = [
    "config",
    "scripts",
    "research_workflow",
    "tests",
    "references/processed",
    "references/curation",
]

EXCLUDE_PARTS = {"__pycache__", ".git", ".DS_Store"}
EXCLUDE_SUFFIXES = {".pyc", ".pyo"}


def should_include(path: Path) -> bool:
    if not path.is_file():
        return False
    if any(part in EXCLUDE_PARTS for part in path.parts):
        return False
    if path.suffix in EXCLUDE_SUFFIXES:
        return False
    return True


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_package(out_path: Path, prefix: str = "hydrogenase_annotation") -> dict[str, object]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    added: list[str] = []
    with tarfile.open(out_path, "w:gz") as tar:
        for rel in INCLUDE_FILES:
            src = ROOT / rel
            if src.is_file():
                tar.add(src, arcname=f"{prefix}/{rel}")
                added.append(rel)
        for rel_dir in INCLUDE_DIRS:
            d = ROOT / rel_dir
            if not d.is_dir():
                continue
            for src in sorted(d.rglob("*")):
                if should_include(src):
                    rel = src.relative_to(ROOT).as_posix()
                    tar.add(src, arcname=f"{prefix}/{rel}")
                    added.append(rel)
    size_mb = round(out_path.stat().st_size / (1024 * 1024), 2)
    digest = sha256_file(out_path)
    return {"output": str(out_path), "files": len(added), "size_mb": size_mb, "sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default_out = ROOT.parent / f".nas/releases/hydrogenase_annotation_server_{time.strftime('%Y%m%d')}.tar.gz"
    parser.add_argument("--out", type=Path, default=default_out, help="Output .tar.gz path")
    args = parser.parse_args()
    info = build_package(args.out.resolve())
    print(f"Packaged {info['files']} files -> {info['output']} ({info['size_mb']} MB, sha256={info['sha256']})")


if __name__ == "__main__":
    main()
