#!/usr/bin/env python3
"""Build and atomically install the repository's pinned Prodigal source."""
from __future__ import annotations

import hashlib, json, os, shutil, subprocess, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tools" / "prodigal" / "source"
PIN = ROOT / "tools" / "prodigal" / "manifest.json"
INSTALLER = Path(__file__).resolve()
MARKER = "Prodigal V2.6.3-r161-boundsfix"

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def tree(pin: dict) -> str:
    actual = {p.name for p in SOURCE.iterdir()}
    expected = set(pin["files"])
    if actual != expected:
        raise RuntimeError(f"Prodigal source file set mismatch: actual={sorted(actual)} expected={sorted(expected)}")
    values = {}
    for name in sorted(expected):
        p = SOURCE / name
        if p.is_symlink() or not p.is_file() or sha(p) != pin["files"][name]:
            raise RuntimeError(f"Prodigal source tamper or missing file: {name}")
        values[name] = pin["files"][name]
    return hashlib.sha256(json.dumps(values, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

def install() -> bool:
    pin = json.loads(PIN.read_text())
    if pin.get('custom_version_marker') != MARKER:
        raise RuntimeError('Prodigal source version marker mismatch')
    source_sha = tree(pin)
    compiler = shutil.which("cc")
    make = shutil.which("make")
    if not compiler or not make:
        raise RuntimeError("existing cc and make are required")
    env = ROOT / ".env"
    target = env / "bin" / "prodigal"
    record = env / "prodigal_build.json"
    identity = {"source_tree_sha256": source_sha, "installer_sha256": sha(INSTALLER),
                "custom_version_marker": pin["custom_version_marker"],
                "compiler": compiler, "compiler_version": subprocess.check_output([compiler, '--version'], text=True).splitlines()[0],
                "build_command": [make, "-B", f"CC={compiler}", 'CFLAGS=-pedantic -Wall -O3', 'LDFLAGS=']}
    if target.is_file() and record.is_file():
        try:
            old = json.loads(record.read_text())
            if all(old.get(k) == v for k, v in identity.items()) and old.get("binary_sha256") == sha(target):
                return False
        except (OSError, ValueError):
            pass
    env.mkdir(parents=True, exist_ok=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="prodigal-build-", dir=env) as td:
        build = Path(td) / "source"
        shutil.copytree(SOURCE, build)
        subprocess.run(identity["build_command"], cwd=build, check=True)
        binary = build / "prodigal"
        if not binary.is_file():
            raise RuntimeError("Prodigal build produced no binary")
        probe = subprocess.run([str(binary), "-v"], capture_output=True, text=True, check=False)
        if probe.returncode != 0 or MARKER not in (probe.stdout + probe.stderr):
            raise RuntimeError(f"unexpected Prodigal -v output: {probe.stdout}{probe.stderr}")
        staged = env / ".prodigal.install.tmp"
        shutil.copy2(binary, staged)
        os.chmod(staged, 0o755)
        if target.is_file() and sha(target) != sha(staged):
            backup = env / f"prodigal.backup.{sha(target)}"
            if not backup.exists():
                shutil.copy2(target, backup)
            if sha(backup) != sha(target):
                raise RuntimeError('Existing Prodigal backup does not match its SHA')
        os.replace(staged, target)
    result = dict(identity, binary_sha256=sha(target), source_files=pin["files"])
    temp_record = record.with_suffix(".json.tmp")
    temp_record.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    os.replace(temp_record, record)
    return True

if __name__ == "__main__":
    print("installed" if install() else "already current")
