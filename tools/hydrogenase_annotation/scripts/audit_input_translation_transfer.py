"""Verify the narrowly scoped translation-input migration.

This module is deliberately independent of catalytic_confirmation.py so it can
be called from the stale-calibration guard without an import cycle.  It does
not recalibrate, change thresholds, or approve a new biological function.
"""
import ast
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lock_sha(lock):
    payload = dict(lock)
    payload.pop("lock_sha256", None)
    return hashlib.sha256((json.dumps(payload, sort_keys=True, indent=2) + "\n").encode()).hexdigest()


def _subscript_translation():
    return ast.Subscript(
        value=ast.Subscript(
            value=ast.Name(id="q", ctx=ast.Load()),
            slice=ast.Constant(value="translation"), ctx=ast.Load()),
        slice=ast.Constant(value=0),
        ctx=ast.Load(),
    )


def _translation_normalized_tree(source):
    """Remove only the three registered translation-input deltas.

    The transformation is structural, not a text wildcard.  Any additional
    change remains visible in the resulting AST and therefore fails equality
    against the frozen pre-translation source.
    """
    tree = ast.parse(source)
    tree.body = [
        node for node in tree.body
        if not (isinstance(node, ast.ImportFrom) and node.module == "input_translation")
    ]
    normalized_count = sum(
        1 for node in ast.walk(tree)
        if isinstance(node, ast.Assign) and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Tuple)
        and isinstance(node.value, ast.Call)
        and isinstance(node.value.func, ast.Name)
        and node.value.func.id == "normalized_translation"
    )
    evidence_count = sum(1 for node in ast.walk(tree) if isinstance(node, ast.Call)
                         and isinstance(node.func, ast.Name) and node.func.id == "translation_evidence")
    found_manifest = set()
    class _DropTranslationAssignment(ast.NodeTransformer):
        def visit_Assign(self, node):
            if (len(node.targets) == 1 and isinstance(node.targets[0], ast.Tuple)
                    and isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Name)
                    and node.value.func.id == "normalized_translation"):
                return None
            if (len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == "translation"
                    and isinstance(node.value, ast.Call)
                    and isinstance(node.value.func, ast.Name)
                    and node.value.func.id == "translation_evidence"):
                return None
            return self.generic_visit(node)

    _DropTranslationAssignment().visit(tree)
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        if node.name == "extract_genbank":
            pass
        for child in ast.walk(node):
            if not isinstance(child, ast.Call):
                continue
            # The manifest is the only registered caller of translation_evidence.
            if isinstance(child.func, ast.Name) and child.func.id == "translation_evidence":
                continue
            for keyword in child.keywords:
                if keyword.arg in {
                    "translation_table_requested", "translation_table_legacy_semantics",
                    "translation_table_used", "translation_table_evidence",
                    "translation_exception_evidence",
                }:
                    found_manifest.add(keyword.arg)
            child.keywords = [k for k in child.keywords if k.arg not in {
                "translation_table_requested", "translation_table_legacy_semantics",
                "translation_table_used", "translation_table_evidence",
                "translation_exception_evidence",
            }]
            # Current extract_genbank uses the normalized local name in the
            # genes record; restore the frozen source expression exactly.
            if isinstance(child.func, ast.Name) and child.func.id == "dict":
                for keyword in child.keywords:
                    if keyword.arg == "sequence" and isinstance(keyword.value, ast.Name) and keyword.value.id == "sequence":
                        keyword.value = _subscript_translation()
    if (normalized_count, evidence_count) not in {(0, 0), (1, 1)}:
        raise ValueError(f"Unexpected translation delta count: normalized={normalized_count}, evidence={evidence_count}")
    if found_manifest not in (set(), {
        "translation_table_requested", "translation_table_legacy_semantics",
        "translation_table_used", "translation_table_evidence", "translation_exception_evidence",
    }):
        raise ValueError("Manifest translation evidence fields do not match the registered scope")
    return tree


def normalized_annotate(source):
    return ast.dump(_translation_normalized_tree(source), include_attributes=False)


def normalized_annotate_source(source):
    """Return a parseable source after the exact registered delta removal."""
    return ast.unparse(_translation_normalized_tree(source))


def verify_registered_transition(root, lock_path, current_stage_fingerprint, old_stage_fingerprint=None, rules=None):
    """Validate a registered migration and return true only for exact scope."""
    root = Path(root)
    lock_path = Path(lock_path)
    if not lock_path.is_absolute():
        lock_path = root / lock_path
    lock = json.loads(lock_path.read_text())
    if lock.get("status") != "input_translation_scope_preserving_transition_verified":
        return False
    if lock_sha(lock) != lock.get("lock_sha256"):
        return False
    current = root / "scripts/annotate.py"
    baseline = root / lock["baseline_annotate"]
    if sha(current) != lock["current_annotate_sha256"] or sha(baseline) != lock["baseline_annotate_sha256"]:
        return False
    if normalized_annotate(current.read_text()) != normalized_annotate(baseline.read_text()):
        return False
    current_stage_sha = hashlib.sha256(json.dumps(current_stage_fingerprint, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if lock["current_stage_fingerprint_sha256"] != current_stage_sha:
        return False
    if old_stage_fingerprint is not None:
        if old_stage_fingerprint.get("discovery") != current_stage_fingerprint.get("discovery"):
            return False
        changed = sorted(name for name in old_stage_fingerprint.get("files", {})
                         if old_stage_fingerprint["files"].get(name) != current_stage_fingerprint["files"].get(name))
        if changed != lock.get("allowed_stage_files"):
            return False
    if rules is not None:
        if lock.get("families") != rules.get("families"):
            return False
        if lock.get("policy_sha256") != rules.get("policy_sha256"):
            return False
        if lock.get("model_manifest_sha256") != rules.get("model_manifest_sha256"):
            return False
        if lock.get("calibration_summary_sha256") != rules.get("calibration_summary_sha256"):
            return False
    if lock["translation_script_sha256"] != sha(root / "scripts/input_translation.py"):
        return False
    for name, digest in lock["protected_sha256"].items():
        if sha(root / name) != digest:
            return False
    return True


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--lock", type=Path, required=True)
    args = parser.parse_args()
    lock = json.loads(args.lock.read_text())
    from catalytic_confirmation import stage_fingerprint
    if not verify_registered_transition(args.root, args.lock, stage_fingerprint(args.root)):
        raise SystemExit("input-translation transition rejected")
    print(json.dumps({"status": "verified", "lock": str(args.lock)}, indent=2))


if __name__ == "__main__":
    main()
