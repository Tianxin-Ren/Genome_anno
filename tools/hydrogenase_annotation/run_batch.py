#!/usr/bin/env python3
"""Run the existing hydrogenase workflow once per manifest sample."""
import argparse
import csv
import hashlib
import json
import os
import signal
import subprocess
import sys
import tempfile
import re
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - the supported runner is POSIX
    fcntl = None

from scripts.install_frozen import verify_native_runtime

ROOT = Path(__file__).resolve().parent
PIN_ROOTS = (ROOT / "scripts", ROOT / "research_workflow", ROOT / "config", ROOT / "references" / "processed", ROOT / "references" / "model_snapshots", ROOT / "references" / "curation")
CURRENT = None
INTERRUPTED = False
FULL_STAGES = ["annotation", "domains", "candidates", "dna_rescue", "candidate_confirmation", "functions", "fes_refinement"]
CLAIM_DIMENSIONS = ("catalytic_identity", "electron_partner", "component_association", "coupling_mechanism", "localization_topology", "reaction_direction_conditions")
CLAIM_ASSERTIONS = "functions/claim_assertions.json"
REQUIRED_OUTPUTS = ["candidate_function_catalog.tsv", "dna_rescue/dna_rescue_manifest.json", "candidates/candidate_screen_manifest.json", "candidate_confirmation/candidate_confirmation_manifest.json", "functions/functional_annotations.tsv", "functions/function_attributes.tsv", "functions/manifest.json", "domains/query_domain_manifest.json", "fes_refinement/manifest.json", CLAIM_ASSERTIONS, "maturation_observation.json"]
REQUIRED_ANNOTATION_OUTPUTS = ["genes.tsv", "complexes.tsv", "catalytic_anchors.tsv", "domain_hits.tsv"]
SAMPLE_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_claim_assertions(out, workflow):
    payload = json.loads((out / CLAIM_ASSERTIONS).read_text())
    if not isinstance(payload, dict) or not isinstance(payload.get("assertions"), dict):
        raise ValueError("claim assertions: expected an assertions object")
    if not workflow.get("claim_assertions_schema_version") or payload.get("schema_version") != workflow["claim_assertions_schema_version"]:
        raise ValueError("claim assertions: schema version mismatch")
    assertions = payload["assertions"]
    if payload.get("candidate_count") != len(assertions) or workflow.get("function_evaluated_candidates") != len(assertions):
        raise ValueError("claim assertions: candidate count mismatch")
    for candidate_id, assertion in assertions.items():
        if not isinstance(assertion, dict) or assertion.get("candidate_id") != candidate_id:
            raise ValueError(f"claim assertions: inconsistent candidate ID {candidate_id}")
        dimensions = assertion.get("dimensions")
        if not isinstance(dimensions, dict) or set(dimensions) != set(CLAIM_DIMENSIONS):
            raise ValueError(f"claim assertions: {candidate_id} missing six dimensions")
        for name, dimension in dimensions.items():
            if not isinstance(dimension, dict) or not isinstance(dimension.get("status"), str) or not dimension["status"]:
                raise ValueError(f"claim assertions: {candidate_id} missing {name} status")


def runtime_fingerprint():
    h = hashlib.sha256()
    files = []
    for root in PIN_ROOTS:
        if root.exists():
            files.extend(p for p in root.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
    declared = []
    for root in (ROOT / "config", ROOT / "research_workflow", ROOT / "references" / "model_snapshots"):
        if not root.exists():
            continue
        for source in root.rglob("*.json"):
            try:
                value = json.loads(source.read_text())
            except (OSError, ValueError):
                continue
            stack = [value]
            while stack:
                current = stack.pop()
                if isinstance(current, dict): stack.extend(current.values())
                elif isinstance(current, list): stack.extend(current)
                elif isinstance(current, str):
                    # JSON model/config files contain prose as well as paths.  Only
                    # follow unambiguous, repository-local path forms; otherwise a
                    # sentence mentioning ``../`` can become an invalid filename.
                    if not current.startswith(("./", "../", "scripts/", "research_workflow/", "config/", "references/", "experiments/", ".env/")):
                        continue
                    try:
                        candidate = Path(current)
                        if not candidate.is_absolute(): candidate = ROOT / current
                        if len(str(candidate)) <= 4096 and candidate.exists(): declared.append(candidate)
                    except (OSError, ValueError):
                        continue
    for candidate in declared:
        if candidate.is_file(): files.append(candidate)
        elif candidate.is_dir(): files.extend(p for p in candidate.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
    package_manifest = ROOT / "package_manifest.json"
    if package_manifest.exists():
        try:
            package = json.loads(package_manifest.read_text())
            package_files = package.get("files")
            if not isinstance(package_files, list):
                raise ValueError("package_manifest.json files must be a list")
            files.append(package_manifest)
            for entry in package_files:
                relative = entry.get("path") if isinstance(entry, dict) else entry
                if not isinstance(relative, str):
                    raise ValueError("package manifest path must be a string")
                path = Path(relative)
                if path.is_absolute() or ".." in path.parts:
                    raise ValueError(f"package manifest path escapes repository: {relative!r}")
                path = ROOT / path
                if not path.is_file():
                    raise ValueError(f"package manifest file is missing: {relative}")
                files.append(path)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            raise ValueError(f"invalid package_manifest.json: {error}") from error
    runtime_manifest = ROOT / "tools/pfsearch_runtime/manifest.json"
    native_runtime = None
    if runtime_manifest.exists():
        try:
            native_runtime = verify_native_runtime(ROOT)
            files.append(runtime_manifest)
            binary = Path(native_runtime["binary_path"])
            if binary.is_absolute() or ".." in binary.parts:
                raise ValueError("pfsearch binary path is unsafe")
            files.append(ROOT / binary)
        except Exception as error:
            raise ValueError(f"invalid pfsearch runtime manifest: {error}") from error
    files.extend((ROOT / ".env/conda-meta").glob("*.json") if (ROOT / ".env/conda-meta").exists() else [])
    for path in sorted(set(files)):
        rel = path.relative_to(ROOT).as_posix().encode()
        h.update(rel + b"\0" + digest(path).encode() + b"\n")
    if native_runtime:
        for path, sha in sorted(native_runtime.get("host_libraries", {}).items()):
            h.update(b"pfsearch-host\0" + path.encode() + b"\0" + sha.encode() + b"\n")
    lock = ROOT / "config/conda-osx-arm64.lock.txt"
    h.update(b"conda-lock\0" + digest(lock).encode())
    for relative in ["run_batch.py", "config/environment.yml"]:
        h.update(relative.encode() + b"\0" + digest(ROOT / relative).encode())
    for binary in ["python", "blastp", "hmmsearch", "prodigal", "mafft", "iqtree3", "pfsearch2.3"]:
        path = ROOT / ".env/bin" / binary
        if path.exists() and path.is_file():
            h.update(binary.encode() + b"\0" + digest(path).encode())
    return h.hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def load_manifest(path):
    data = json.loads(path.read_text())
    if not isinstance(data, list):
        raise ValueError("manifest must be a JSON list")
    seen_sample, seen_path = set(), set()
    result = []
    for item in data:
        if not isinstance(item, dict) or not isinstance(item.get("sample"), str) or not isinstance(item.get("path"), str):
            raise ValueError("each manifest item needs string sample and path")
        sample = item["sample"]
        if not SAMPLE_RE.fullmatch(sample) or sample in {".", ".."}:
            raise ValueError(f"invalid sample name: {sample!r}")
        source = Path(item["path"]).expanduser().resolve()
        if not sample or sample in seen_sample:
            raise ValueError(f"duplicate/empty sample: {sample!r}")
        if source in seen_path:
            raise ValueError(f"duplicate input path: {source}")
        if not source.is_file():
            raise ValueError(f"input is not a file: {source}")
        expected = item.get("sha256")
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
            raise ValueError(f"sha256 must be a 64-hex string for {sample}")
        actual = digest(source)
        if expected.lower() != actual:
            raise ValueError(f"input sha256 mismatch for {sample}: {actual}")
        seen_sample.add(sample); seen_path.add(source)
        result.append({"sample": sample, "path": str(source), "sha256": actual})
    if not result:
        raise ValueError("manifest must not be empty")
    return result


def workflow_ok(out, expected_pin, item):
    workflow = out / "workflow.json"
    if not workflow.is_file():
        return False, "missing workflow.json"
    try:
        data = json.loads(workflow.read_text())
        if data.get("status") != "complete":
            return False, f"workflow status={data.get('status')!r}"
        if data.get("script_sha256") != digest(ROOT / "run_hydrogenase.py"):
            return False, "run_hydrogenase.py hash mismatch"
        if data.get("annotation_reused") is not False:
            return False, "annotation was reused"
        stages = data.get("stages", [])
        if [s.get("stage") for s in stages] != FULL_STAGES or any(s.get("exit_code") != 0 for s in stages):
            return False, "full seven-stage exit failure/missing"
        output_hashes = data.get("output_sha256", {})
        missing = [relative for relative in REQUIRED_OUTPUTS if relative not in output_hashes]
        if missing:
            return False, "missing required output hashes: " + ",".join(missing)
        if data.get('function_rules_version'):
            required_rules=['functions/function_rules/'+name for name in ['manifest.json','predictions.tsv','evidence.json','member_hits.tsv']]
            if any(name not in output_hashes for name in required_rules):return False,'missing functional rule output hashes'
        for relative, expected in output_hashes.items():
            path = out / relative
            if not path.is_file() or digest(path) != expected:
                return False, f"output hash mismatch: {relative}"
        verify_claim_assertions(out, data)
        annotation_manifest = out / "annotation/run_manifest.json"
        if not annotation_manifest.is_file():
            return False, "missing annotation input manifest"
        inputs = json.loads(annotation_manifest.read_text()).get("input_files", [])
        if not any(x.get("path") == item["path"] and x.get("sha256") == item["sha256"] for x in inputs):
            return False, "annotation input DNA hash mismatch"
        annotation_hashes = json.loads(annotation_manifest.read_text()).get("output_sha256", {})
        if not annotation_hashes:
            return False, "annotation output hashes are empty"
        missing_annotation = [relative for relative in REQUIRED_ANNOTATION_OUTPUTS if relative not in annotation_hashes]
        if missing_annotation:
            return False, "missing essential annotation output hashes: " + ",".join(missing_annotation)
        for relative, expected in annotation_hashes.items():
            path = out / "annotation" / relative
            if not path.is_file() or digest(path) != expected:
                return False, f"annotation output hash mismatch: {relative}"
        if data.get("runtime_fingerprint") and data["runtime_fingerprint"] != expected_pin:
            return False, "workflow runtime fingerprint mismatch"
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return False, f"invalid workflow: {error}"
    return True, "verified"


def attempt_number(sample_dir):
    numbers = []
    for path in sample_dir.glob("attempt_*"):
        try:
            numbers.append(int(path.name.split("_", 1)[1]))
        except (IndexError, ValueError):
            continue
    return max(numbers, default=0) + 1


def command_for(item, out, args):
    return [sys.executable, str(ROOT / "run_hydrogenase.py"), "--genome", item["path"], "--out", str(out), "--threads", str(args.threads), "--mode", args.mode, "--topology", args.topology]


def run_sample(item, out, args):
    global CURRENT
    out.parent.mkdir(parents=True, exist_ok=True)
    command = command_for(item, out, args)
    log = (out.parent / (out.name + ".runner.log")).open("w")
    process = subprocess.Popen(command, cwd=ROOT, start_new_session=True, stdout=log, stderr=subprocess.STDOUT, text=True)
    CURRENT = process
    if INTERRUPTED:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        code = process.wait()
    finally:
        log.close()
        CURRENT = None
    return code, command


def summary_row(item, status, output="", error="", workflow=None):
    workflow = workflow or {}
    predictions = workflow.get("function_predictions", {})
    preferred = workflow.get("preferred_function_predictions", {})
    confidence = workflow.get("preferred_prediction_confidence", {})
    class_predictions=json.dumps(workflow.get('function_class_predictions',{}),sort_keys=True)
    claims = str(Path(output) / CLAIM_ASSERTIONS) if status == "complete" else ""
    claim_status, claim_count = ("valid", workflow["function_evaluated_candidates"]) if status == "complete" else ("", "")
    return {"sample": item["sample"], "status": status, "output": output, "candidates": workflow.get("candidates", ""), "core": workflow.get("core_anchors", ""), "function_count": workflow.get("function_evaluated_candidates", ""), "function_predictions": json.dumps(predictions, sort_keys=True), "preferred_function_predictions": json.dumps(preferred, sort_keys=True), "preferred_prediction_confidence": json.dumps(confidence, sort_keys=True), "function_class_predictions":class_predictions, "error": error, "claim_assertions": claims, "claim_status": claim_status, "claim_count": claim_count}


def write_summary(out, rows):
    fields = ["sample", "status", "output", "candidates", "core", "function_count", "function_predictions", "preferred_function_predictions", "preferred_prediction_confidence", "function_class_predictions", "error", "claim_assertions", "claim_status", "claim_count"]
    fd, name = tempfile.mkstemp(prefix="summary.tsv.", dir=out)
    try:
        with os.fdopen(fd, "w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
            writer.writeheader(); writer.writerows(rows); stream.flush(); os.fsync(stream.fileno())
        os.replace(name, out / "summary.tsv")
    finally:
        if os.path.exists(name): os.unlink(name)
    atomic_json(out / "summary.json", {"samples": rows, "count": len(rows)})


def handle_sigint(signum, frame):
    global INTERRUPTED
    INTERRUPTED = True
    if CURRENT is not None:
        try:
            os.killpg(CURRENT.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass


def main(argv=None):
    global INTERRUPTED
    INTERRUPTED = False
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--mode", choices=["genome", "assembly"], default="assembly")
    parser.add_argument("--topology", choices=["auto", "linear", "circular"], default="auto")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")
    manifest_path = args.manifest.resolve(); out = args.out.resolve()
    items = load_manifest(manifest_path)
    pin = runtime_fingerprint()
    args_fingerprint = {"threads": args.threads, "mode": args.mode, "topology": args.topology, "entrypoint": digest(ROOT / "run_hydrogenase.py"), "runtime": pin}
    lock_path = out / ".batch.lock"
    out.mkdir(parents=True, exist_ok=True)
    lock = lock_path.open("a+")
    if fcntl is None:
        raise RuntimeError("POSIX fcntl flock is required")
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError(f"batch output is locked: {out}")
    state_path = out / "batch_state.json"
    if args.resume:
        if not state_path.is_file():
            raise ValueError("--resume requires batch_state.json")
        state = json.loads(state_path.read_text())
        if state.get("manifest_path") != str(manifest_path) or state.get("manifest_sha256") != digest(manifest_path) or state.get("items") != items or state.get("args") != args_fingerprint:
            raise ValueError("manifest, input SHA, or runtime/arguments fingerprint changed")
        rows = state.get("rows", [])
    else:
        if state_path.exists() or any(p.name != ".batch.lock" for p in out.iterdir()):
            raise ValueError("output exists; use --resume for an existing batch")
        rows = []
        state = {"status": "running", "manifest_path": str(manifest_path), "manifest_sha256": digest(manifest_path), "items": items, "args": args_fingerprint, "rows": rows}
        atomic_json(state_path, state)
    by_sample = {row["sample"]: row for row in rows}
    old_signal = signal.signal(signal.SIGINT, handle_sigint)
    try:
        for item in items:
            row = by_sample.get(item["sample"])
            if INTERRUPTED:
                state["status"] = "interrupted"; state["current_sample"] = None; atomic_json(state_path, state); write_summary(out, state["rows"]); return 130
            if row and row.get("status") == "complete":
                ok, reason = workflow_ok(Path(row["output"]), pin, item)
                if ok:
                    continue
                row["status"] = "pending"; row["error"] = "stale complete: " + reason
            sample_dir = out / item["sample"]; attempt = sample_dir / f"attempt_{attempt_number(sample_dir):03d}"
            state["current_sample"] = item["sample"]
            atomic_json(state_path, state)
            try:
                code, command = run_sample(item, attempt, args)
                workflow = json.loads((attempt / "workflow.json").read_text()) if (attempt / "workflow.json").is_file() else {}
                verified, reason = workflow_ok(attempt, pin, item) if code == 0 else (False, f"exit_code={code}")
                status = "complete" if verified else "failed"
                row = summary_row(item, status, str(attempt), "" if status == "complete" else reason, workflow)
            except Exception as error:
                row = summary_row(item, "failed", str(attempt), repr(error))
            by_sample[item["sample"]] = row
            state["rows"] = [by_sample.get(x["sample"], summary_row(x, "pending")) for x in items]
            state["current_sample"] = None
            state["status"] = "running"
            atomic_json(state_path, state); write_summary(out, state["rows"])
            if INTERRUPTED:
                state["status"] = "interrupted"; atomic_json(state_path, state); write_summary(out, state["rows"]); return 130
        state["status"] = "complete" if all(r.get("status") == "complete" for r in state["rows"]) else "failed"
        atomic_json(state_path, state); write_summary(out, state["rows"])
        return 0 if state["status"] == "complete" else 1
    finally:
        signal.signal(signal.SIGINT, old_signal)
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN); lock.close()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(error, file=sys.stderr); raise SystemExit(2)
