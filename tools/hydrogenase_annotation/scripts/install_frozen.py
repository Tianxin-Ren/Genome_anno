#!/usr/bin/env python3
"""Install an exported, hash-checked runtime into the workspace-local .env."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import stat
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _runtime_rows(manifest: dict) -> tuple[dict, list[dict], list[dict], list[dict]]:
    binary = manifest.get("binary") if isinstance(manifest.get("binary"), dict) else {"path": manifest.get("binary_path"), "sha256": manifest.get("binary_sha256")}
    sources = manifest.get("sources", manifest.get("source_files", []))
    if isinstance(sources, dict):
        sources = [{"path": p, "sha256": h} for p, h in sources.items()]
    hosts = manifest.get("host_libraries", manifest.get("host_libs", []))
    if isinstance(hosts, dict):
        hosts = [{"path": p, "sha256": h} for p, h in hosts.items()]
    bundled = manifest.get("bundled_libraries", [])
    if isinstance(bundled, dict):
        bundled = [{"path": p, "sha256": h} for p, h in bundled.items()]
    return binary, list(sources), list(hosts), list(bundled)


def _package_file(root: Path, value: str) -> Path:
    if not isinstance(value, str) or not value or Path(value).is_absolute():
        raise RuntimeError("pfsearch package path must be relative: " + repr(value))
    rel = Path(value)
    if ".." in rel.parts or rel.as_posix() != value:
        raise RuntimeError("unsafe pfsearch package path: " + repr(value))
    path = root / rel
    current = root
    for part in rel.parts:
        current = current / part
        if current.is_symlink():
            raise RuntimeError("symlinked pfsearch package path: " + value)
    if root.resolve() not in path.resolve().parents or not path.is_file():
        raise RuntimeError("missing pfsearch package file: " + value)
    return path


def verify_native_runtime(root: Path) -> dict:
    """Verify the pinned non-conda pfsearch runtime and host dylibs."""
    manifest_path = root / "tools/pfsearch_runtime/manifest.json"
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise RuntimeError("missing pfsearch runtime manifest")
    manifest = json.loads(manifest_path.read_text())
    binary, sources, hosts, bundled = _runtime_rows(manifest)
    if not binary.get("path") or not binary.get("sha256") or not sources or (not hosts and not bundled):
        raise RuntimeError("incomplete pfsearch runtime manifest: binary, sources, and host_libraries or bundled_libraries are required")
    rows = [binary] + sources
    for row in rows:
        path = _package_file(root, row.get("path"))
        if path.is_symlink() or sha256(path) != row.get("sha256"):
            raise RuntimeError("pfsearch packaged file missing or hash mismatch: " + str(row.get("path")))
    bundled_names = set()
    bundled_hashes = {}
    for row in bundled:
        path = _package_file(root, row.get("path"))
        name = path.name
        if name in bundled_names:
            raise RuntimeError("duplicate bundled pfsearch library basename: " + name)
        bundled_names.add(name)
        if sha256(path) != row.get("sha256"):
            raise RuntimeError("pfsearch bundled library hash mismatch: " + str(row.get("path")))
        bundled_hashes[name] = row["sha256"]
    host_hashes = {}
    for row in hosts:
        path = Path(row["path"])
        if not path.is_absolute() or path.is_symlink() or not path.is_file():
            raise RuntimeError("missing pfsearch host library: " + str(path))
        actual = sha256(path)
        if actual != row.get("sha256"):
            raise RuntimeError("pfsearch host library hash mismatch: " + str(path))
        host_hashes[str(path)] = actual
    installed = root / ".env/bin/pfsearch2.3"
    if installed.exists() and (installed.is_symlink() or sha256(installed) != binary["sha256"]):
        raise RuntimeError("installed pfsearch2.3 is missing or hash-mismatched")
    installed_libs = root / ".env/lib/pfsearch_runtime"
    if installed.exists() and bundled:
        for name, digest in bundled_hashes.items():
            path = installed_libs / name
            try:
                _package_file(root, str(path.relative_to(root)))
            except RuntimeError:
                raise RuntimeError("installed pfsearch bundled library missing or hash-mismatched: " + name)
    return {"manifest_sha256": sha256(manifest_path), "binary_path": binary["path"], "binary_sha256": binary["sha256"], "host_libraries": host_hashes, "bundled_libraries": bundled_hashes}


def validate_package(root: Path) -> dict:
    manifest_path = root / "package_manifest.json"
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise RuntimeError("missing package_manifest.json")
    manifest = json.loads(manifest_path.read_text())
    rows = manifest.get("files")
    if manifest.get("status") != "exported" or not isinstance(rows, list) or not rows:
        raise RuntimeError("invalid exported package manifest")
    seen = set()
    for row in rows:
        rel = row.get("path")
        if not isinstance(rel, str) or not rel or Path(rel).is_absolute():
            raise RuntimeError(f"unsafe package path: {rel!r}")
        candidate = Path(rel)
        if ".." in candidate.parts or candidate.as_posix() != rel:
            raise RuntimeError(f"unsafe package path: {rel!r}")
        if rel in seen or not row.get("sha256"):
            raise RuntimeError(f"duplicate package path: {rel}")
        seen.add(rel)
        path = root / candidate
        resolved = path.resolve()
        current = root
        components_symlinked = False
        for part in candidate.parts:
            current = current / part
            components_symlinked |= current.is_symlink()
        if root.resolve() not in resolved.parents or components_symlinked or path.is_symlink() or not path.is_file():
            raise RuntimeError(f"missing or symlinked package file: {rel}")
        if sha256(path) != row.get("sha256"):
            raise RuntimeError(f"package hash mismatch: {rel}")
    required = {"config/conda-osx-arm64.explicit-md5.txt", "run_hydrogenase.py", "scripts/install_frozen.py", "scripts/install_prodigal.py", "tools/prodigal/manifest.json", "tools/pfsearch_runtime/manifest.json"}
    listed = seen
    if not required.issubset(listed):
        raise RuntimeError("package manifest omits required installer/lock/entrypoint files")
    prodigal = json.loads((root / "tools/prodigal/manifest.json").read_text())
    source_files = {f"tools/prodigal/source/{name}" for name in prodigal.get("files", {})}
    if not source_files.issubset(listed):
        raise RuntimeError("package manifest omits pinned Prodigal source files")
    runtime = json.loads((root / "tools/pfsearch_runtime/manifest.json").read_text())
    binary, sources, _, bundled = _runtime_rows(runtime)
    runtime_files = {binary.get("path")} | {row.get("path") for row in sources} | {row.get("path") for row in bundled}
    if None in runtime_files or not runtime_files.issubset(listed):
        raise RuntimeError("package manifest omits pinned pfsearch binary/source files")
    return manifest


def preflight(root: Path, verify_only: bool) -> tuple[Path | None, Path | None]:
    if not verify_only and (platform.system() != "Darwin" or platform.machine() != "arm64"):
        raise RuntimeError("installer requires macOS arm64")
    conda = Path(shutil.which("conda") or "")
    cc = Path(shutil.which("cc") or "")
    make = Path(shutil.which("make") or "")
    if not verify_only and (not conda.is_file() or not cc.is_file() or not make.is_file()):
        raise RuntimeError("conda, cc, and make are required")
    lock = root / "config/conda-osx-arm64.explicit-md5.txt"
    if not lock.is_file() or lock.is_symlink():
        raise RuntimeError("missing exact explicit lock")
    return (conda if conda.is_file() else None), lock


def install(root: Path, manifest: dict, conda: Path, lock: Path) -> dict:
    env = root / ".env"
    if env.exists() or env.is_symlink():
        raise RuntimeError("refusing to overwrite existing .env")
    cache = root / ".conda-pkgs"
    cache.mkdir(exist_ok=True)
    env_vars = dict(os.environ, CONDA_NO_PLUGINS="true", CONDA_PKGS_DIRS=str(cache))
    failure = {"status": "failed", "package_manifest_sha256": sha256(root / "package_manifest.json"), "lock_sha256": sha256(lock), "installer_sha256": sha256(Path(__file__))}
    try:
        # Conda prefixes contain absolute shebangs and metadata; never rename a
        # successfully-created prefix from a staging path.
        subprocess.run([str(conda), "create", "-y", "-p", str(env), "--file", str(lock)], cwd=root, env=env_vars, check=True)
        build = subprocess.run([sys.executable, str(root / "scripts/install_prodigal.py")], cwd=root, env=env_vars, check=False, capture_output=True, text=True)
        if build.returncode:
            failure["prodigal_returncode"] = build.returncode
            (env / "install_frozen_failure.json").write_text(json.dumps(failure, indent=2) + "\n")
            raise RuntimeError("pinned Prodigal build failed; partial .env retained with failure receipt")
        runtime = json.loads((root / "tools/pfsearch_runtime/manifest.json").read_text())
        binary, _, _, bundled = _runtime_rows(runtime)
        target = env / "bin/pfsearch2.3"
        if target.exists():
            raise RuntimeError("refusing to overwrite existing pfsearch2.3")
        shutil.copyfile(root / binary["path"], target)
        os.chmod(target, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
        if sha256(target) != binary["sha256"]:
            raise RuntimeError("copied pfsearch2.3 hash mismatch")
        if bundled:
            bundled_target = env / "lib/pfsearch_runtime"
            if bundled_target.is_symlink() or (bundled_target.exists() and not bundled_target.is_dir()):
                raise RuntimeError("refusing symlinked bundled pfsearch library directory")
            bundled_target.mkdir(parents=True, exist_ok=True)
            if bundled_target.parent.is_symlink():
                raise RuntimeError("refusing symlinked bundled pfsearch library parent")
            for row in bundled:
                source = _package_file(root, row["path"])
                target = bundled_target / source.name
                if target.exists() or target.is_symlink():
                    raise RuntimeError("refusing to overwrite bundled pfsearch library: " + source.name)
                shutil.copyfile(source, target)
                os.chmod(target, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
                if sha256(target) != row["sha256"]:
                    raise RuntimeError("copied bundled pfsearch library hash mismatch: " + source.name)
        native = verify_native_runtime(root)
        receipt = {"status": "installed", "package_manifest_sha256": sha256(root / "package_manifest.json"), "lock_sha256": sha256(lock), "installer_sha256": sha256(Path(__file__)), "prodigal_receipt_sha256": sha256(env / "prodigal_build.json"), "native_runtime": native}
        (env / "install_frozen_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        return receipt
    except Exception as exc:
        failure["error"] = str(exc)
        if env.is_dir() and not (env / "install_frozen_failure.json").exists():
            (env / "install_frozen_failure.json").write_text(json.dumps(failure, indent=2) + "\n")
        raise


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    manifest = validate_package(root)
    conda, lock = preflight(root, args.verify_only)
    native = verify_native_runtime(root)
    if args.verify_only:
        print(json.dumps({"status": "package_verified", "file_count": len(manifest["files"]), "lock_sha256": sha256(lock), "native_runtime": native}))
        return 0
    print(json.dumps(install(root, manifest, conda, lock)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
