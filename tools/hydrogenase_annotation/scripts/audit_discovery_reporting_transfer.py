"""Strict, reporting-only bridge for the R112 discovery output additions."""
import ast
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lock_sha(lock):
    body = dict(lock)
    body.pop("lock_sha256", None)
    return hashlib.sha256((json.dumps(body, sort_keys=True, indent=2) + "\n").encode()).hexdigest()


class _ReportingDelta(ast.NodeTransformer):
    """Remove only the registered report-producing syntax from a source AST."""
    def __init__(self, name):
        self.name = name

    def visit_ImportFrom(self, node):
        if self.name == "discovery_evidence" and node.module == "pathlib":
            return None
        return node

    def visit_FunctionDef(self, node):
        if self.name == "discovery_evidence" and node.name == "load_targets":
            return None
        return self.generic_visit(node)

    def visit_If(self, node):
        direct_reporting = (self.name == "catalytic_confirmation" and
            any(isinstance(x, ast.ImportFrom) and x.module == "audit_discovery_reporting_transfer" for x in node.body) and
            any(isinstance(x, ast.Call) and isinstance(x.func, ast.Attribute) and
                isinstance(x.func.value, ast.Name) and x.func.value.id == "rules" and x.func.attr == "get" and
                x.args and isinstance(x.args[0], ast.Constant) and x.args[0].value == "reporting_transition"
                for x in ast.walk(node.test)))
        if direct_reporting:
            return None
        return self.generic_visit(node)

    def visit_Assign(self, node):
        if self.name == "discovery_evidence" and any(isinstance(t, ast.Name) and t.id == "tbl" for t in node.targets):
            return None
        if self.name == "audit_maturation_domains" and any(isinstance(t, ast.Name) and t.id == "target" for t in node.targets):
            return None
        return self.generic_visit(node)

    def visit_List(self, node):
        kept = []
        for item in node.elts:
            if isinstance(item, ast.Constant) and item.value == "--tblout":
                continue
            if isinstance(item, ast.Call) and isinstance(item.func, ast.Name) and item.func.id == "str":
                if len(item.args) == 1 and isinstance(item.args[0], ast.Name) and item.args[0].id in {"tbl", "target"}:
                    continue
            if self.name == "annotate" and isinstance(item, ast.Constant) and item.value in {"family.domtbl", "family.tbl", "blast.tsv"}:
                continue
            kept.append(self.visit(item))
        node.elts = kept
        return node


def normalized_reporting_source(name, text):
    """Return an AST dump after removing only the registered report deltas."""
    tree = _ReportingDelta(name).visit(ast.parse(text))
    return ast.dump(ast.fix_missing_locations(tree), include_attributes=False)


def _normalized_reporting_text(name, text):
    return ast.unparse(ast.fix_missing_locations(_ReportingDelta(name).visit(ast.parse(text))))


def _check_panel(root, lock, rules):
    report_path = Path(lock["feature_equivalence_report"])
    if not report_path.is_absolute():
        report_path = root / report_path
    if sha(report_path) != lock["feature_equivalence_report_sha256"]:
        return False
    report = json.loads(report_path.read_text())
    if report.get("status") != "passed" or report.get("panel_units") != 8 or report.get("feature_count") != 30:
        return False
    if len(report.get("panel", [])) != report["panel_units"]:
        return False
    for unit in report["panel"]:
        old = Path(unit["old_features_path"]); new = Path(unit["new_features_path"])
        if sha(old) != unit["old_features_sha256"] or sha(new) != unit["new_features_sha256"]:
            return False
        old_rows, new_rows = json.loads(old.read_text()), json.loads(new.read_text())
        if old_rows != new_rows or len(new_rows) != unit["feature_count"]:
            return False
        from catalytic_confirmation import qualifies
        decisions = {}
        for r in new_rows:
            p = rules.get("families", {}).get(r.get("family"), {})
            decisions[r["protein_id"]] = qualifies(r, p)
        if decisions != unit.get("new_decisions") or decisions != unit.get("old_decisions"):
            return False
        if not unit.get("features_identical") or not unit.get("decision_vectors_identical"):
            return False
    return True


def verify_reporting_transition(root, path, currentstage, rules):
    """Verify a registered reporting transition without changing calibration."""
    root, path = Path(root), Path(path)
    if not path.is_absolute():
        path = root / path
    lock = json.loads(path.read_text())
    if lock.get("status") != "reporting_only_transition_verified" or lock_sha(lock) != lock.get("lock_sha256"):
        return False
    scope = {"scripts/annotate.py", "scripts/discovery_evidence.py",
             "scripts/audit_maturation_domains.py", "scripts/catalytic_confirmation.py"}
    if (set(lock.get("current_scripts_sha256", {})) != scope or
            set(lock.get("baseline_scripts_sha256", {})) != scope or
            set(lock.get("baseline_scripts", {})) != scope or
            set(lock.get("translation_projection", {})) != scope):
        return False
    if lock["translation_projection"] != lock["baseline_scripts_sha256"]:
        return False
    current_sha = hashlib.sha256(json.dumps(currentstage, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if lock.get("current_stage_fingerprint_sha256") != current_sha:
        return False
    if lock.get("baseline_stage_fingerprint") != rules.get("feature_stage_fingerprint"):
        return False
    old_rules = root / lock["baseline_thresholds"]
    if sha(old_rules) != lock["baseline_thresholds_sha256"] or sha(old_rules) != lock["protected_rules_sha256"]:
        return False
    old = json.loads(old_rules.read_text())
    unchanged_rules = dict(rules); unchanged_rules.pop("reporting_transition", None); unchanged_rules.pop("topology_transition", None)
    if old != unchanged_rules:
        return False
    for name, digest in lock["current_scripts_sha256"].items():
        if sha(root / name) != digest:
            return False
    for name, digest in lock["baseline_scripts_sha256"].items():
        if sha(root / lock["baseline_scripts"][name]) != digest:
            return False
        if normalized_reporting_source(name.rsplit("/", 1)[-1][:-3], (root / name).read_text()) != normalized_reporting_source(name.rsplit("/", 1)[-1][:-3], (root / lock["baseline_scripts"][name]).read_text()):
            return False
    if sha(root / "scripts/audit_discovery_reporting_transfer.py") != lock["verifier_sha256"]:
        return False
    translation = root / lock["translation_lock"]
    if sha(translation) != lock["translation_lock_sha256"]:
        return False
    tl = json.loads(translation.read_text())
    if lock_sha(tl) != tl.get("lock_sha256") or tl.get("baseline_annotate_sha256") != sha(root / tl["baseline_annotate"]):
        return False
    for name, digest in tl.get("protected_sha256", {}).items():
        protected = root / lock.get("baseline_scripts", {}).get(name, name)
        if sha(protected) != digest:
            return False
    if tl.get("policy_sha256") != rules.get("policy_sha256") or tl.get("model_manifest_sha256") != rules.get("model_manifest_sha256") or tl.get("calibration_summary_sha256") != rules.get("calibration_summary_sha256"):
        return False
    from audit_input_translation_transfer import normalized_annotate
    if normalized_annotate(_normalized_reporting_text("annotate", (root / "scripts/annotate.py").read_text())) != normalized_annotate((root / tl["baseline_annotate"]).read_text()):
        return False
    projection = json.loads(json.dumps(currentstage))
    for name, digest in lock.get("translation_projection", {}).items():
        if name in projection.get("files", {}):
            projection["files"][name] = digest
        elif name in projection.get("discovery", {}).get("files", {}):
            projection["discovery"]["files"][name] = digest
        else:
            return False
    if hashlib.sha256(json.dumps(projection, sort_keys=True, separators=(",", ":")).encode()).hexdigest() != tl.get("current_stage_fingerprint_sha256"):
        return False
    return _check_panel(root, lock, rules)
