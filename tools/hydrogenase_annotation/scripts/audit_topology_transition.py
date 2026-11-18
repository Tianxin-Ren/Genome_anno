"""Verify the R169 circular-input topology transition.

This bridge binds the R112 reporting state to the reviewed annotation and
graph topology changes while preserving the original calibration parameters.
"""
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lock_sha(lock):
    body = dict(lock)
    body.pop("lock_sha256", None)
    return hashlib.sha256((json.dumps(body, sort_keys=True, indent=2) + "\n").encode()).hexdigest()


def verify_topology_transition(root, path, current_stage, rules):
    root, path = Path(root), Path(path)
    if not path.is_absolute():
        path = root / path
    lock = json.loads(path.read_text())
    if lock.get("status") != "circular_topology_scope_preserving_transition_verified":
        return False
    if lock_sha(lock) != lock.get("lock_sha256"):
        return False
    baseline_sha = hashlib.sha256(json.dumps(rules.get("feature_stage_fingerprint"), sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if lock.get("baseline_stage_fingerprint_sha256") != baseline_sha:
        return False
    expected = hashlib.sha256(json.dumps(current_stage, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if lock.get("current_stage_fingerprint_sha256") != expected:
        return False
    parent = root / lock["parent_reporting_lock"]
    if sha(parent) != lock["parent_reporting_lock_sha256"]:
        return False
    parent_lock = json.loads(parent.read_text())
    if parent_lock.get("status") != "reporting_only_transition_verified":
        return False
    old_rules_path = root / parent_lock['baseline_thresholds']
    if sha(old_rules_path) != parent_lock['baseline_thresholds_sha256']:
        return False
    if json.loads(old_rules_path.read_text()) != {key:value for key,value in rules.items()
            if key not in {'reporting_transition','topology_transition'}}:
        return False
    for name, digest in lock["protected_sha256"].items():
        if sha(root / name) != digest:
            return False
    if sha(root / "config/catalytic_confirmation.json") != rules["policy_sha256"]:
        return False
    if sha(root / "references/processed/catalytic_confirmation/manifest.json") != rules["model_manifest_sha256"]:
        return False
    for name, digest in lock["parent_current_scripts_sha256"].items():
        if name != "scripts/annotate.py" and sha(root / name) != digest:
            return False
    if sha(root / "scripts/annotate.py") != lock["current_annotate_sha256"]:
        return False
    if sha(root / 'scripts/module_graph.py') != lock['current_module_graph_sha256']:
        return False
    if lock.get("baseline_annotate_sha256") != parent_lock.get("current_scripts_sha256", {}).get("scripts/annotate.py"):
        return False
    if lock.get("linear_replay_sha256") != sha(root / lock["linear_replay"]):
        return False
    if lock.get("circular_test_sha256") != sha(root / lock["circular_test"]):
        return False
    if lock.get("raw_fna_rescue_sha256") != sha(root / lock["raw_fna_rescue"]):
        return False
    opaque = lock['opaque955_calibration_control']
    required = {'proteins.faa','confirmation/features.json','catalytic_confirmation.tsv',
                'active_center_summary.tsv','catalytic_sites.tsv','catalytic_anchors.tsv'}
    if set(opaque['identical_artifacts_sha256']) != required:
        return False
    for name, digest in opaque['identical_artifacts_sha256'].items():
        if any(sha(root / opaque[key] / name) != digest for key in ('frozen_annotation','fresh_annotation')):
            return False
    return True


if __name__ == "__main__":
    import argparse
    from catalytic_confirmation import stage_fingerprint
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--lock", type=Path, required=True)
    args = parser.parse_args()
    rules = json.loads((args.root / "config/catalytic_confirmation_thresholds.json").read_text())
    if not verify_topology_transition(args.root, args.lock, stage_fingerprint(args.root), rules):
        raise SystemExit("circular-topology transition rejected")
    print(json.dumps({"status": "verified", "lock": str(args.lock)}, indent=2))
