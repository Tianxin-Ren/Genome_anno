"""Build a conservative, multi-dimensional view of function predictions.

The legacy TSV labels remain the compatibility surface.  This sidecar keeps
identity, association, mechanism, topology and reaction direction separate so
that an unresolved partner does not erase a supported catalytic identity (or
promote a specific mechanism without evidence).
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
from pathlib import Path

VERSION = "claim-assertions-v4"
_ATTRIBUTE_TABLE = "function_attributes.tsv"
_FUNCTION_TABLE = "functional_annotations.tsv"
_PREFERRED_TABLE = "preferred_predictions.json"
_COMPOSITION_TABLE = "pyridine_role_composition_evidence.json"
_PRODUCER_PINS = {
    _FUNCTION_TABLE: "functional_annotations_sha256",
    _ATTRIBUTE_TABLE: "function_attributes_sha256",
    _PREFERRED_TABLE: "preferred_prediction_output_sha256",
    _COMPOSITION_TABLE: "pyridine_role_composition_output_sha256",
}


def _json(value):
    if isinstance(value, str):
        try:
            return json.loads(value)
        except (TypeError, ValueError):
            return value
    return value


def _attrs(rows, candidate_id):
    result = {}
    for row in rows or []:
        if not isinstance(row, dict):
            continue
        name = row.get("attribute", "")
        if name:
            owner = row.get("candidate_id", candidate_id)
            if candidate_id and owner != candidate_id:
                raise ValueError("Function attribute belongs to a different candidate")
            result.setdefault(name, []).append(dict(row, candidate_id=owner))
    return result


def _evidence_group(row):
    source = str(row.get("evidence_source") or "").lower()
    name = str(row.get("attribute") or "").lower()
    if name in {"candidate_co_linked_h2_pathway", "fefe_complementary_scaffold_context"}:
        return "association_and_context"
    if any(token in source or token in name for token in ("reaction", "nucleotide", "electron", "coupling", "activity")):
        return "reaction_and_mechanism"
    if any(token in source or token in name for token in ("complex", "partner", "association", "role", "hypoth")):
        return "association_and_context"
    if any(token in source or token in name for token in ("domain", "profile", "identity", "core", "sequence", "confirmation")):
        return "catalytic_sequence"
    return "other"


def _values(rows):
    result = []
    for row in rows:
        value = {"attribute": row.get("attribute"), "status": row.get("status"),
                 "value": _json(row.get("value")), "evidence_source": row.get("evidence_source"),
                 "evidence_dependencies": _json(row.get("evidence_dependencies")),
                 "limitations": row.get("limitations")}
        key = dict(value, candidate_id=row["candidate_id"])
        evidence_id = hashlib.sha256(json.dumps(
            key, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()
        role = "support" if _positive_status(value["status"]) else "constraint"
        if value["attribute"] == "Reference_assay_reaction_capabilities" or (
            value["attribute"] == "FeFe_rubrerythrin_function" and value["status"] == "assay_supported_candidate"
        ):
            role = "reference_only"
        result.append(dict(value, evidence_id=evidence_id, evidence_table=_ATTRIBUTE_TABLE,
                           evidence_role=role, evidence_group=_evidence_group(row)))
    return result


_NEGATIVE_STATUS_WORDS = ("unresolved", "unsupported", "incomplete", "missing", "blocked",
                          "ambiguous", "not_assessed", "not_tested", "conflict", "unknown")


def _positive_status(status):
    value = str(status or "").lower()
    if not value or any(word in value for word in _NEGATIVE_STATUS_WORDS):
        return False
    return value in {
        "supported", "complete", "complete_supported", "domain_profile_supported",
        "profile_in_candidate_architecture", "supported_role_module_homology",
        "reference_like_modules_detected", "sequence_supported_candidate",
        "sequence_supported_core_region", "composition_supported_reorganized",
        "coding_locus_supported_conditional", "local_coupling_potential_uncalibrated",
        "compatible_reaction_candidate_uncalibrated", "reference_supported_uncalibrated",
        "local_membrane_h2_redox_potential_uncalibrated", "local_quinone_coupling_potential_uncalibrated",
        "reference_like_h2_fd_coupling_potential", "candidate_h2_formate_coupling",
        "reference_like_module_uncalibrated", "reference_like_fes_relay_candidate",
        "candidate_co_linked_h2_pathway",
        "assay_supported_candidate", "sequence_supported_fusion",
    }


def _dimension(attrs, names):
    rows = [row for name, values in attrs.items() if name.lower() in names for row in values]
    supported = [row for row in rows if _positive_status(row.get("status"))]
    if supported:
        status = "supported"
    elif rows:
        status = "unresolved"
    else:
        status = "not_assessed"
    result = {"status": status, "evidence": _values(supported or rows)}
    if supported:
        constraints = _values(row for row in rows if not _positive_status(row.get("status")))
        if constraints:
            result["constraint_evidence"] = constraints
    return result


ELECTRON_PARTNERS = {
    "candidate_H2_pyridine_nucleotide_coupling": "pyridine_nucleotide_unspecified",
    "candidate_H2_ferredoxin_compatible_core": "ferredoxin_compatible",
    "candidate_H2_quinone_coupling": "quinone",
    "candidate_H2_F420_reduction": "F420",
    "candidate_H2_H4MPT_hydrogenation_potential": "H4MPT_derivative",
    "candidate_H2_THF_derivative_hydrogenation_potential": "THF_derivative",
    "candidate_H2_formate_coupling": "formate_associated",
    "candidate_FeFe_rubrerythrin_H2_H2O2_coupling": "hydrogen_peroxide",
}

_PARTNER_ATTRIBUTES = {
    "pyridine_nucleotide_unspecified": {"local_pyridine_hydrogenase_coupling_potential",
                                        "role_composed_h2_pyridine_coupling_potential"},
    "ferredoxin_compatible": {"mbh_jkln_fd_redox_core", "ech_h2_fd_coupling_potential"},
    "quinone": {"huc_local_h2_quinone_coupling_potential", "hya_like_hyd_local_h2_quinone_coupling_potential",
                "hyb_like_hyd_local_h2_quinone_coupling_potential"},
    "F420": {"frh_f420_reduction_potential"},
    "H4MPT_derivative": {"canonical_h4mpt_branch_reaction_potential", "fegp_h4mpt_hydrogenation_potential"},
    "THF_derivative": {"fegp_thf_derivative_hydrogenation_potential"},
    "formate_associated": {"nife_fhl_formate_h2_pathway", "tonn_formate_h2_coupling_potential",
                           "h2_formate_coupling_potential"},
}
_PYRIDINE_MODULE = {"pyridine_nucleotide_electron_module"}
_ALL_PARTNER_ATTRIBUTES = set().union(*_PARTNER_ATTRIBUTES.values())
_REGION_IDENTITY = {"sequence_supported_core_region", "sequence_supported_uncalibrated_core_region"}
_SPLIT_IDENTITY = "sequence_supported_uncalibrated_split_center"
_SPLIT_ATTRIBUTE = "NiFe_split_catalytic_core_potential"
_SPLIT_FUNCTION = "NiFe_split_catalytic_center_supported_function_unresolved"


def _split_proof(rows):
    for row in rows:
        if row.get("status") != "sequence_supported_candidate":
            continue
        value = _json(row.get("value"))
        members = value.get("members") if isinstance(value, dict) else None
        if not isinstance(members, list) or len(members) != 2:
            continue
        ids = [member.get("protein_id") for member in members
               if isinstance(member, dict) and member.get("protein_id")]
        if len(ids) == 2 and len(set(ids)) == 2:
            return value, members
    return None, []


def build_assertion(record, attribute_rows=(), partner_assignment=None):
    """Return one serialisable candidate assertion without inventing evidence."""
    record = dict(record or {})
    candidate_id = record.get("candidate_id", "")
    attrs = _attrs(attribute_rows, candidate_id)
    identity = record.get("identity_status") or "not_assessed"
    _, split_members = _split_proof(attrs.get(_SPLIT_ATTRIBUTE, []))
    split_valid = identity == _SPLIT_IDENTITY and bool(split_members)
    if identity in _REGION_IDENTITY or split_valid or identity in {
        "calibrated_high", "sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical",
    }:
        identity_status = "supported"
    elif identity == _SPLIT_IDENTITY or identity in {"unresolved", "unconfirmed", "evidence_insufficient",
                      "hydrogenase_like_lineage_unconfirmed", "sequence_support_center_model_unresolved"}:
        identity_status = "unresolved"
    else:
        identity_status = "not_assessed"
    catalytic_rows = [row for rows in attrs.values() for row in rows
                      if any(token in str(row.get("attribute", "")).lower()
                             for token in ("core", "identity"))]
    dimensions = {
        "catalytic_identity": {"status": identity_status, "value": identity, "family": record.get("family"),
                               "evidence": _values(catalytic_rows)},
        "electron_partner": _dimension(attrs, _ALL_PARTNER_ATTRIBUTES),
        "component_association": _dimension(attrs, {"candidate_complex_architecture", "role_composed_h2_pyridine_coupling_potential",
                                                    "mbh_jkln_fd_redox_core", "group4_local_membrane_h2_redox_potential",
                                                    "nife_large_small_core_potential", "mvh_hdr_coupling_hypothesis",
                                                    "fefe_rubrerythrin_function",
                                                    "hya_like_hyd_local_h2_quinone_coupling_potential",
                                                    "hyb_like_hyd_local_h2_quinone_coupling_potential"}),
        "coupling_mechanism": _dimension(attrs, {"electron_coupling", "bifurcation", "ion_translocation", "ion_species"}),
        "localization_topology": _dimension(attrs, {"group4_local_membrane_h2_redox_potential", "ech_membrane_organization",
                                                     "membrane_topology", "localization_topology"}),
        "reaction_direction_conditions": {
            "status": "unresolved",
            "value": {"direction": record.get("in_vivo_direction", "unknown"),
                      "query_activity": record.get("query_activity", "not_tested"),
                      "selectivity": record.get("NAD_NADP_selectivity", "unresolved")},
            "reference_assays": list(_json(record.get("reference_assay_capabilities")) or []),
            "evidence": _values(row for name in ("H2_reaction_role", "nucleotide_selectivity", "query_activity",
                                                 "Reference_assay_reaction_capabilities")
                                for row in attrs.get(name, [])),
        },
    }
    if identity in _REGION_IDENTITY:
        dimensions["catalytic_identity"]["scope"] = "catalytic_core_region_only"
    elif split_valid:
        dimensions["catalytic_identity"]["scope"] = "split_catalytic_center_only"
    limits = sorted(set(_json(record.get("unresolved_attributes")) or []))
    limits.extend(sorted({row.get("attribute") for rows in attrs.values() for row in rows
                          if row.get("status") in {"unresolved", "not_tested", "not_assessed_or_not_detected"} and row.get("attribute")}))
    alternatives = _json(record.get("prediction_alternatives", []))
    if not isinstance(alternatives, list):
        alternatives = [alternatives] if alternatives else []
    preferred_members = _json(record.get("preferred_members", []))
    if not isinstance(preferred_members, list):
        preferred_members = [preferred_members] if preferred_members else []
    risks = _json(record.get("prediction_risks", []))
    if not isinstance(risks, list):
        risks = [risks] if risks else []
    unknowns = _json(record.get("experimental_unknowns")) or {}
    coarse = record.get("coarse_function_prediction") or record.get("function_prediction", "")
    preferred = record.get("preferred_function", "")
    association_conditions = _json(record.get("prediction_support", []))
    if split_valid and not coarse:
        coarse = _SPLIT_FUNCTION
    if coarse == _SPLIT_FUNCTION and not split_valid:
        coarse = ""
    if preferred == _SPLIT_FUNCTION or (not preferred and identity == _SPLIT_IDENTITY):
        preferred_members = split_members if split_valid else []
    if preferred == _SPLIT_FUNCTION and not split_valid:
        preferred = "possible_hydrogenase"
        association_conditions = []
    partner = dimensions["electron_partner"]
    module_rows = _values(row for name, rows in attrs.items() if name.lower() in _PYRIDINE_MODULE
                          for row in rows)
    if module_rows:
        partner["module_evidence"] = module_rows
        if partner["status"] == "not_assessed":
            partner["status"] = "unresolved"
    hypothesis_electron = ELECTRON_PARTNERS.get(preferred)
    if hypothesis_electron in _PARTNER_ATTRIBUTES:
        matching = _dimension(attrs, _PARTNER_ATTRIBUTES[hypothesis_electron])
        if matching["evidence"]:
            scope = "sequence_prediction" if matching["status"] == "supported" else "conditional_unresolved"
            partner.update(
                hypothesis=hypothesis_electron, hypothesis_scope=scope,
                hypothesis_evidence=matching["evidence"],
                hypothesis_evidence_ids=sorted({row["evidence_id"] for row in matching["evidence"]}))
    for row in attrs.get('FeFe_rubrerythrin_function', []):
        fusion = _json(row.get('value'))
        if row.get('status') == 'assay_supported_candidate' and isinstance(fusion, dict) and fusion.get('reference_assay'):
            reference_evidence = _values([row])
            dimensions['reaction_direction_conditions']['reference_assays'].append(fusion['reference_assay'])
            dimensions['reaction_direction_conditions']['evidence'].extend(reference_evidence)
            dimensions['electron_partner'].update(status='supported', hypothesis='hydrogen_peroxide',
                                                   hypothesis_scope='exact_reference_assay_potential',
                                                   evidence=reference_evidence,
                                                   hypothesis_evidence_ids=[reference_evidence[0]["evidence_id"]])
    if (preferred == "candidate_membrane_associated_H2_redox_system"
            and dimensions["localization_topology"]["status"] == "supported"):
        dimensions["localization_topology"].update(hypothesis="membrane_associated", catalytic_face="unresolved")
    dimensions["coupling_mechanism"]["unresolved_attributes"] = ["electron_bifurcation", "ion_translocation", "ion_species"]
    # Preserve an explicitly supplied unresolved competition result as an
    # auditable constraint.  Independent architecture evidence must remain
    # supported; the unresolved partner state must not be silently dropped.
    if partner_assignment:
        dimensions["component_association"]["competition"] = partner_assignment
    association = dimensions["component_association"]["status"]
    association_reason = "claim_specific_architecture_evidence"
    co_context = []
    if preferred == "candidate_CO_linked_H2_pathway":
        for row in attrs.get(preferred, []):
            value = _json(row.get("value"))
            if (_positive_status(row.get("status")) and isinstance(value, dict)
                    and value.get("entry_route") == "hyc_core_local_CODH_context"):
                co_context.append(row)
    if co_context:
        co_evidence = _values(co_context)
        association, association_reason = "conditional", "conditional_CO_context_positive_source_evidence"
        dimensions["component_association"]["hypothesis"] = "CO_linked_H2_context"
        dimensions["component_association"]["hypothesis_scope"] = "conditional_genome_context"
        dimensions["component_association"]["hypothesis_evidence"] = co_evidence
        dimensions["electron_partner"].update(
            hypothesis="carbon_monoxide_associated",
            hypothesis_scope="conditional_genome_context",
            hypothesis_evidence=co_evidence,
            hypothesis_evidence_ids=sorted({row["evidence_id"] for row in co_evidence}))
    assignment_state = (partner_assignment or {}).get("status")
    if assignment_state == "ambiguous" or any(r in risks for r in
            ("cross_anchor_partner_ownership_ambiguous", "partner_assignment_ambiguous")):
        association, association_reason = "ambiguous", "competing_member_explanations"
    elif assignment_state == "shared_role_conditionally_allowed":
        association, association_reason = "conditional", "candidate_shared_role"
    elif "sequence_equivalent_partner_locus_ambiguous" in risks:
        association, association_reason = "ambiguous", "sequence_equivalent_partner_loci"
    if "Group4_core_connection_architecture_unresolved" in risks:
        association, association_reason = "conditional", "nonstandard_core_connection"
        novelty = "nonstandard_module_combination"
        coarse = "H2_core_supported_partner_unresolved"
    else:
        novelty = "reference_matched" if record.get("reference_reaction_status") == "exact_reference_sequence_potential" else "not_assessed"
    context_rows = []
    if identity_status == 'supported' and record.get('family') == 'FeFe':
        for row in attrs.get('FeFe_complementary_scaffold_context', []):
            value = _json(row.get('value'))
            if (row.get('status') in {'conditional_partial_context', 'ambiguous_partial_context', 'unresolved'}
                    and isinstance(value, dict) and value.get('anchor_id') == candidate_id and value.get('candidates')):
                context_rows.append(row)
    if context_rows:
        context_evidence = _values(context_rows)
        dimension = dimensions['component_association']
        dimension.setdefault('module_evidence', []).extend(context_evidence)
        # Partial scaffold homology is a separate hypothesis, never a selected functional partner.
        if association in {'not_assessed', 'unresolved'}:
            states = {_json(row['value']).get('status') for row in context_rows}
            if states & {'conditional', 'ambiguous'}:
                association = 'ambiguous' if 'ambiguous' in states else 'conditional'
                association_reason = 'partial_complementary_scaffold_context'
                dimension.update(hypothesis='complementary_catalytic_scaffold',
                    hypothesis_scope='conditional_partial_structural_context', hypothesis_evidence=context_evidence)
                if association == 'conditional' and len(context_rows) == 1:
                    dimension['hypothesis_members'] = _json(context_rows[0]['value']).get('selected_members', [])
    dimensions["component_association"]["status"] = association
    # A predicted partner/module match is not an unconditional coupling claim
    # when the member association is only conditional or has competing owners.
    # Keep the module evidence and supported identity intact; only narrow the
    # scope of the dependent prediction.  Exact reference-assay scopes are
    # already explicitly bounded and must not be downgraded here.
    if (association in {"conditional", "ambiguous"}
            and partner.get("hypothesis_scope") == "sequence_prediction"):
        partner["hypothesis_scope"] = "conditional_unresolved"
    for dimension in dimensions.values():
        dimension["evidence_ids"] = sorted({row["evidence_id"]
            for key in ("evidence", "constraint_evidence", "hypothesis_evidence", "module_evidence")
            for row in dimension.get(key, ())})
    for key, fields in (
        ("catalytic_identity", ("identity_status", "family")),
        ("reaction_direction_conditions", ("in_vivo_direction", "query_activity", "NAD_NADP_selectivity",
                                           "reference_assay_capabilities")),
    ):
        dimensions[key]["record_inputs"] = [{
            "file": _FUNCTION_TABLE, "candidate_id": candidate_id, "fields": [field for field in fields if field in record]}]
    if "preferred_function" in record:
        for key in ("electron_partner", "localization_topology"):
            dimensions[key]["record_inputs"] = [{
                "file": _PREFERRED_TABLE, "json_path": [candidate_id, "preferred_function"]}]
    dimensions["component_association"]["record_inputs"] = [
        {"file": _PREFERRED_TABLE, "json_path": [candidate_id, field]}
        for field in ("preferred_members", "prediction_risks", "prediction_support") if field in record]
    if partner_assignment:
        dimensions["component_association"]["record_inputs"].append({
            "file": _COMPOSITION_TABLE, "json_path": [candidate_id, "evidence", "partner_assignment"]})
    dependencies = {
        "catalytic_identity": ["sequence/core evidence", "identity_status"],
        "electron_partner": ["partner/module evidence", "electron-transfer compatibility"],
        "component_association": ["member roles", "locality or architecture evidence"],
        "coupling_mechanism": ["reaction/mechanism evidence"],
        "localization_topology": ["membrane/topology evidence"],
        "reaction_direction_conditions": ["query activity or reference reaction evidence"],
    }
    return {
        "candidate_id": candidate_id,
        "schema_version": VERSION,
        "current_supported_annotation": coarse,
        "preferred_function_hypothesis": preferred,
        "association_status": association,
        "association_reason": association_reason,
        "member_evidence": preferred_members,
        "association_conditions": association_conditions,
        "dimensions": dimensions,
        "claim_dependencies": dependencies,
        "limiting_evidence": limits,
        "alternative_explanations": alternatives,
        "prediction_risks": risks,
        "experimental_unknowns": unknowns,
        "novelty_status": record.get("novelty_status", novelty),
        "prediction_confidence": record.get("prediction_confidence", record.get("confidence", "unresolved")),
    }


def build_assertions(records, attributes=(), composition=None):
    by_id = {}
    for row in attributes or ():
        by_id.setdefault(row.get("candidate_id"), []).append(row)
    return {row["candidate_id"]: build_assertion(row, by_id.get(row["candidate_id"], []),
            (composition or {}).get(row["candidate_id"], {}).get("evidence", {}).get("partner_assignment"))
            for row in records}


def write_assertions(path: Path, records, attributes=(), composition=None, source_binding=None):
    path = Path(path)
    payload = {"schema_version": VERSION, "candidate_count": len(records or []),
               "source_binding": source_binding or {"status": "in_memory_unbound"},
               "assertions": build_assertions(records, attributes, composition)}
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    return path


def from_results(functions, output=None):
    functions = Path(functions)
    manifest_path = functions / "manifest.json"
    blobs = {name: (functions / name).read_bytes() for name in _PRODUCER_PINS
             if name != _COMPOSITION_TABLE or (functions / name).is_file()}
    hashes = {name: hashlib.sha256(data).hexdigest() for name, data in blobs.items()}
    manifest_data = manifest_path.read_bytes() if manifest_path.is_file() else None
    manifest = json.loads(manifest_data) if manifest_data is not None else {}
    for name, key in _PRODUCER_PINS.items():
        if key in manifest and manifest[key] != hashes.get(name):
            raise ValueError(f"Claim producer output changed: {name}")
    fully_bound = manifest.get("status") == "claim_specific_explanation_complete" and all(
        key in manifest for key in _PRODUCER_PINS.values())
    source_binding = {
        "status": "producer_manifest_verified" if fully_bound else "legacy_inputs_unpinned",
        "files_sha256": hashes,
        "producer_manifest_sha256": hashlib.sha256(manifest_data).hexdigest() if manifest_data is not None else None,
    }
    # Function evidence can embed a full candidate architecture search ledger.
    csv.field_size_limit(sys.maxsize)
    records = list(csv.DictReader(io.StringIO(blobs[_FUNCTION_TABLE].decode("utf-8"), newline=""), delimiter="\t"))
    attributes = list(csv.DictReader(io.StringIO(blobs[_ATTRIBUTE_TABLE].decode("utf-8"), newline=""), delimiter="\t"))
    preferred = json.loads(blobs[_PREFERRED_TABLE])
    ids = {row["candidate_id"] for row in records}
    if len(ids) != len(records) or ids != preferred.keys() or any(row["candidate_id"] not in ids for row in attributes):
        raise ValueError("Claim producer outputs have inconsistent candidate IDs")
    for record in records:
        record.update(preferred[record["candidate_id"]])
    composition = json.loads(blobs.get(_COMPOSITION_TABLE, b"{}"))
    return write_assertions(output or functions / "claim_assertions.json", records, attributes, composition, source_binding)
