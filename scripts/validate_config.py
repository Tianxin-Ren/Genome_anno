#!/usr/bin/env python3
"""Validate layered project/platform configuration without running Snakemake."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any

import yaml

WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW_ROOT))

from common.scripts.runtime import apply_runtime_config, deep_merge, unresolved_variables


def read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise TypeError(f"Expected a YAML mapping: {path}")
    return value


def resolve_ref(value: str, project_dir: Path) -> Path:
    expanded = Path(os.path.expanduser(os.path.expandvars(value)))
    if expanded.is_absolute():
        return expanded
    project_candidate = project_dir / expanded
    if project_candidate.exists():
        return project_candidate
    return WORKFLOW_ROOT / expanded


def load_config(project_dir: Path, platform_arg: str | None) -> dict[str, Any]:
    project = read_yaml(project_dir / "pipeline.yaml")
    platform_ref = platform_arg or (project.get("runtime", {}) or {}).get("platform")
    platform = read_yaml(resolve_ref(platform_ref, project_dir)) if platform_ref else {}
    merged = deep_merge(platform, project)

    sample_ref = (merged.get("data", {}) or {}).get("samples", "samples.yaml")
    sample_path = resolve_ref(sample_ref, project_dir)
    if sample_path.exists():
        sample_cfg = read_yaml(sample_path)
        merged = deep_merge(merged, sample_cfg)
        if sample_cfg.get("project_root"):
            merged.setdefault("paths", {})["project_root"] = sample_cfg["project_root"]

    for module_name, module_info in (merged.get("modules", {}) or {}).items():
        if not isinstance(module_info, dict) or not module_info.get("enable"):
            continue
        module_ref = module_info.get("path")
        if not module_ref:
            continue
        module_path = resolve_ref(module_ref, project_dir)
        if not module_path.exists():
            raise FileNotFoundError(f"Enabled module '{module_name}' config not found: {module_path}")
        merged = deep_merge(merged, read_yaml(module_path))

    return dict(apply_runtime_config(merged))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, help="Project config directory")
    parser.add_argument("--platform", help="Override platform YAML")
    parser.add_argument(
        "--check-paths",
        action="store_true",
        help="Also require the Conda root and configured database paths to exist",
    )
    args = parser.parse_args()

    project_dir = Path(args.project).expanduser().resolve()
    try:
        config = load_config(project_dir, args.platform)
    except (FileNotFoundError, TypeError, ValueError, yaml.YAMLError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    unresolved = unresolved_variables(config)
    if unresolved:
        print("ERROR: unresolved environment variables:", file=sys.stderr)
        for location in unresolved:
            print(f"  - {location}", file=sys.stderr)
        return 2

    enabled = [
        name
        for name, value in (config.get("modules", {}) or {}).items()
        if isinstance(value, dict) and value.get("enable")
    ]
    paths = config.get("paths", {}) or {}
    print(f"project_root: {config.get('project_root', '')}")
    print(f"results:      {paths.get('results', '')}")
    print(f"logs:         {paths.get('logs', '')}")
    print(f"work:         {paths.get('work', '')}")
    print(f"conda_envs:   {config.get('conda_env_dir', '')}")
    print(f"modules:      {', '.join(enabled) if enabled else '(none)'}")

    if args.check_paths:
        missing: list[str] = []
        env_root = config.get("conda_env_dir")
        if env_root and not Path(env_root).exists():
            missing.append(f"paths.conda_envs -> {env_root}")
        for name, path in (config.get("databases", {}) or {}).items():
            if path and not Path(str(path)).exists():
                missing.append(f"databases.{name} -> {path}")
        if missing:
            print("ERROR: configured paths do not exist:", file=sys.stderr)
            for item in missing:
                print(f"  - {item}", file=sys.stderr)
            return 2

    print("configuration: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
