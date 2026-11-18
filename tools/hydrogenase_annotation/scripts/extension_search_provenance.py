"""Accept an older extension search only when its search source is unchanged."""
import ast
import hashlib
import json

from audit_topology_transition import lock_sha


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def search_source(path):
    tree = ast.parse(path.read_text())
    # The excluded function runs after the search manifest and HSPs are written.
    tree.body = [node for node in tree.body if not (
        isinstance(node, ast.FunctionDef) and node.name == "exact_reference_patches"
        or isinstance(node, ast.ImportFrom) and node.module == "circular_coordinates"
        and [(alias.name, alias.asname) for alias in node.names] == [("valid_cds_interval", None)]
    )]
    return ast.dump(tree)


def compatible_extension_search(root, source_digest):
    name = "scripts/complex_extensions.py"
    current = sha(root / name)
    if source_digest == current:
        return True
    rules = json.loads((root / "config/catalytic_confirmation_thresholds.json").read_text())
    lock_path = rules.get("topology_transition")
    if not lock_path:
        return False
    lock = json.loads((root / lock_path).read_text())
    review = lock.get("structure_extension_search_compatibility", {})
    archived = review.get("archived_source")
    if (not archived or source_digest != review.get("archived_sha256")
            or current != review.get("current_sha256")
            or lock_sha(lock) != lock.get("lock_sha256")):
        return False
    return sha(root / archived) == source_digest and search_source(root / archived) == search_source(root / name)
