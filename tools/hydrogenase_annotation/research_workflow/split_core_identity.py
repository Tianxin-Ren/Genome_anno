"""Conservative reader for sequence-supported split NiFe centers.

This module is deliberately a reader, not a predictor.  A proof is emitted only
when the producer tables and their model manifest are bound by ``run_manifest``.
Missing (legacy) tables simply produce no proof.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

SPLIT_IDENTITY = "sequence_supported_uncalibrated_split_center"
SPLIT_LABEL = "NiFe_split_catalytic_center_supported_function_unresolved"
RULE_VERSION = "split_core_identity_v1"


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _json(value: Any, default: Any = None) -> Any:
    if isinstance(value, (list, dict)):
        return value
    try:
        return json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return default


def _index(value: Any, key: str = "protein_id") -> dict[str, dict[str, Any]]:
    if isinstance(value, dict):
        return {str(k): (v if isinstance(v, dict) else {"sequence": v}) for k, v in value.items()}
    return {str(row.get(key) or row.get("id")): row for row in (value or []) if isinstance(row, dict) and (row.get(key) or row.get("id"))}


def _bound_files(annotation: Path, manifest: dict[str, Any], root: Path) -> dict[str, str] | None:
    """Verify the producer outputs named by run_manifest; return absolute-path hashes."""
    expected = manifest.get("output_sha256") or {}
    names = (
        "catalytic_anchors.tsv",
        "catalytic_center_candidates.tsv",
        "catalytic_center_components.tsv",
        "active_center_summary.tsv",
        "genes.tsv",
    )
    bound: dict[str, str] = {}
    manifest_path = annotation / "run_manifest.json"
    for name in names:
        path = annotation / name
        if not path.exists():
            return None
        digest = _sha(path)
        recorded = expected.get(name) or expected.get(f"annotation/{name}")
        if recorded and recorded != digest:
            raise ValueError(f"run_manifest hash mismatch: {name}")
        if recorded is None:
            return None
        bound[str(path.resolve())] = digest
    protein_path = annotation / "proteins.faa"
    if not protein_path.exists() or _sha(protein_path) != manifest.get("normalized_proteins_sha256"):
        raise ValueError("proteins.faa is not bound to run_manifest")
    bound[str(protein_path.resolve())] = _sha(protein_path)
    model_path = root / "references" / "processed" / "split_centers" / "manifest.json"
    recorded = manifest.get("split_center_model_manifest_sha256")
    if not model_path.exists() or not recorded:
        raise ValueError("split-center model manifest is not bound")
    model_digest = _sha(model_path)
    if recorded != model_digest:
        raise ValueError("split-center model manifest hash mismatch")
    bound[str(model_path.resolve())] = model_digest
    bound[str(manifest_path.resolve())] = _sha(manifest_path)
    reader_path = root / "research_workflow" / "split_core_identity.py"
    if reader_path.exists():
        bound[str(reader_path.resolve())] = _sha(reader_path)
    return bound


def _complete_gene(g: dict[str, Any]) -> bool:
    if str(g.get("partial", "")).lower() in {"true", "1", "yes"}:
        return False
    if str(g.get("coding_completeness", g.get("status", ""))).lower() in {"partial", "incomplete"}:
        return False
    return str(g.get("coordinate_evidence", "")).lower() in {"genomic_cds", "genbank_cds", "bound_genomic_cds"}


def prepare(annotation, anchors, genes, sequences, root):
    """Return ``(proofs_by_anchor_id, input_hashes)`` for valid split centers.

    ``proofs_by_anchor_id`` is intentionally empty for absent, legacy, ambiguous,
    or incomplete evidence.  It never asserts coupling, activity, substrate, or
    physical assembly.
    """
    annotation = Path(annotation)
    if annotation.is_file():
        annotation = annotation.parent
    root = Path(root)
    manifest_path = annotation / "run_manifest.json"
    required = ["catalytic_center_candidates.tsv", "catalytic_center_components.tsv", "active_center_summary.tsv", "catalytic_anchors.tsv", "run_manifest.json"]
    if not all((annotation / name).exists() for name in required):
        return {}, {}
    manifest = json.loads(manifest_path.read_text())
    hashes = _bound_files(annotation, manifest, root)
    if hashes is None:
        return {}, {}

    anchor_rows = _index(anchors)
    gene_rows = _index(genes)
    seq_rows = _index(sequences, "protein_id")
    anchor_table = _rows(annotation / "catalytic_anchors.tsv")
    candidate_rows = _rows(annotation / "catalytic_center_candidates.tsv")
    component_rows = _rows(annotation / "catalytic_center_components.tsv")
    summary_rows = _rows(annotation / "active_center_summary.tsv")
    all_candidates = []
    for row in candidate_rows:
        proteins = _json(row.get("component_proteins"), [])
        sites = _json(row.get("sites"), [])
        by_role = {}
        for site in sites:
            by_role.setdefault(site.get("center_component_id"), []).append(site)
        if (row.get("center_status") == "split_template_ligand_set_supported" and
                row.get("pairing_status") == "competitive" and
                row.get("organization") == "multiple_polypeptides" and
                len(proteins) == 2 and len(set(proteins)) == 2 and
                len(sites) == 4 and len(by_role) == 2 and all(len(v) == 2 for v in by_role.values()) and
                all(s.get("status") == "present" and s.get("residue_role") == "direct_metal_ligand" and s.get("query_position") for s in sites)):
            all_candidates.append((row, proteins, sites))

    proofs: dict[str, dict[str, Any]] = {}
    for anchor in anchor_table:
        anchor_id = anchor["protein_id"]
        if anchor_id not in anchor_rows or anchor.get("catalytic_family") != "NiFe":
            continue
        candidates = [(r, p, s) for r, p, s in all_candidates if r.get("catalytic_anchor") == anchor_id]
        if len(candidates) != 1:
            continue
        candidate, proteins, sites = candidates[0]
        # A selected split component may not be silently shared by another anchor.
        if any(r.get("catalytic_anchor") != anchor_id and set(proteins) & set(p) for r, p, _ in all_candidates):
            continue
        summary = next((r for r in summary_rows if r.get("protein_id") == anchor_id), None)
        if not summary or summary.get("active_center_status") != "multi_polypeptide_ligand_set_supported":
            continue
        if summary.get("center_organization") != "multiple_polypeptides" or _json(summary.get("center_component_proteins"), []) != proteins:
            continue
        comps = [r for r in component_rows if r.get("catalytic_anchor") == anchor_id and r.get("center_id") == candidate.get("center_id")]
        if len(comps) != 2 or {r.get("protein_id") for r in comps} != set(proteins) or any(r.get("component_status") != "template_ligand_set_supported" or r.get("template_span_status") != "supported" for r in comps):
            continue
        comp_by_protein = {r.get("protein_id"): r for r in comps}
        selected_genes = []
        for protein_id in proteins:
            g = gene_rows.get(protein_id)
            if g is not None and not g.get("protein_id"):
                g = dict(g)
                g["protein_id"] = protein_id
            if not g or protein_id not in seq_rows or not _complete_gene(g):
                break
            selected_genes.append(g)
        if len(selected_genes) != 2:
            continue
        if len({g.get("contig") for g in selected_genes}) != 1 or len({str(g.get("strand")) for g in selected_genes}) != 1:
            continue
        if any(not g.get("start") or not g.get("end") or not g.get("contig_length") or str(g.get("strand")) not in {"1", "-1"} or str(g.get("topology", "")).lower() not in {"linear", "circular"} for g in selected_genes):
            continue
        if anchor.get("bin_id") and len({g.get("bin_id", anchor.get("bin_id")) for g in selected_genes}) != 1:
            continue
        # Every component owns exactly two direct sites, and each observed residue
        # must agree with the caller-supplied exact protein sequence.
        if any(s.get("protein_id") not in comp_by_protein for s in sites):
            continue
        members = []
        coordinates = {}
        for gene in selected_genes:
            pid = gene["protein_id"]
            comp = comp_by_protein[pid]
            role = comp["component_id"]
            seq = seq_rows[pid].get("sequence", "")
            mapped = [s for s in sites if s["protein_id"] == pid]
            positions = {int(s["query_position"]) for s in mapped}
            if len(mapped) != 2 or len(positions) != 2 or any(
                s.get("center_component_id") != role or s.get("reference_id") != comp["reference_id"] or
                not 1 <= int(s["query_position"]) <= len(seq) or
                seq[int(s["query_position"]) - 1] != s.get("observed_residue")
                for s in mapped
            ):
                break
            spans = _json(comp.get("alignment_components"), [])
            if not spans:
                break
            start = min(int(h["query_start"]) for h in spans)
            end = max(int(h["query_end"]) for h in spans)
            if not 1 <= start <= min(positions) <= max(positions) <= end <= len(seq):
                break
            members.append(dict(role=role, protein_id=pid, query_interval=[start, end]))
            coordinates[pid] = {k: gene.get(k) for k in ("contig", "start", "end", "strand", "coordinate_evidence")}
        if len(members) != 2:
            continue
        proofs[anchor_id] = {
            "rule_version": RULE_VERSION,
            "identity": SPLIT_IDENTITY,
            "label": SPLIT_LABEL,
            "anchor": anchor_id,
            "members": members,
            "coordinates": coordinates,
            "actual_sites": sites,
            "template": {"center_id": candidate.get("center_id"), "source": candidate.get("source"), "source_sha256": candidate.get("source_sha256")},
            "coupling_status": "unresolved",
            "limitations": "sequence-supported split center only; function, activity, assembly, and coupling remain unresolved",
        }
    return proofs, hashes
