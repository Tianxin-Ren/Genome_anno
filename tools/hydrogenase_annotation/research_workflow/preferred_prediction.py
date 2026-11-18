"""Qualitative preferred-function prediction layered on the legacy catalog.

This module intentionally does not recalibrate identity thresholds.  It turns
the evidence already collected by the workflow into a *best explanation* and
keeps unresolved experimental questions separate from sequence evidence.
"""
import json

PREDICTION_FIELDS = (
    "preferred_function", "preferred_family", "prediction_kind",
    "prediction_confidence", "prediction_summary", "prediction_support",
    "preferred_members", "prediction_alternatives", "prediction_risks",
    "experimental_unknowns", "prediction_confidence_scope",
)
VERSION = "preferred-prediction-v27"
_SCOPE = "qualitative_evidence_not_calibrated_probability"
_ARTIFICIAL_MEDIATOR_SCOPE = "purified_protein_artificial_mediator_reference_potential"


def _json(value):
    if isinstance(value, str):
        try:
            return json.loads(value)
        except (TypeError, ValueError):
            return value
    return value


def _text(value):
    value = _json(value)
    if isinstance(value, (list, tuple, set)):
        return " ".join(map(str, value))
    return str(value or "")


def _claims(evidence):
    """Return self-contained claim dictionaries, without reading any files."""
    evidence = evidence if isinstance(evidence, dict) else {}
    out = []
    if any(k in evidence for k in ("hydrogen_reaction_core_status", "pyridine_nucleotide_module_status", "core_status")):
        item = dict(evidence)
        item.setdefault("claim_name", "evidence")
        out.append(item)
    claims = _json(evidence.get("claims", {}))
    if isinstance(claims, dict):
        for name, claim in claims.items():
            if isinstance(claim, dict):
                item = dict(claim)
                nested = item.get("evidence")
                if isinstance(nested, dict):
                    merged = dict(nested)
                    merged.update(item)
                    item = merged
                if "preferred_role_options" not in item and evidence.get("preferred_role_options") is not None:
                    item["preferred_role_options"] = evidence.get("preferred_role_options")
                if name == "composition":
                    item.setdefault("coupling_potential_status", item.get("status", ""))
                    if item.get("status") == "local_coupling_potential_uncalibrated":
                        item.update(association_supported=True, core_association_supported=True,
                                    module_association_supported=True, coupling_linkage_supported=True)
                item.setdefault("claim_name", name)
                out.append(item)
    for key in ("pyridine_result", "pyridine_function_evidence", "pyridine"):
        value = _json(evidence.get(key))
        if isinstance(value, dict) and not any(value is c for c in out):
            item = dict(value)
            nested = item.get("evidence")
            if isinstance(nested, dict):
                merged = dict(nested)
                merged.update(item)
                item = merged
            if "preferred_role_options" not in item and evidence.get("preferred_role_options") is not None:
                item["preferred_role_options"] = evidence.get("preferred_role_options")
            item.setdefault("claim_name", key)
            out.append(item)
    for index, claim in enumerate(_json(evidence.get("alternate_hypotheses", [])) or []):
        if isinstance(claim, dict):
            item = dict(claim)
            item.setdefault("claim_name", f"alternate_hypothesis_{index}")
            item["is_alternate_hypothesis"] = True
            out.append(item)
    # Keep every architecture atomic.  An alternative's roles, statuses and
    # selected members must never be merged into its parent claim.
    index_parent = 0
    while index_parent < len(out):
        parent = out[index_parent]
        index_parent += 1
        alternatives = _json(parent.get("alternative_architectures", []))
        for index, alternative in enumerate(alternatives if isinstance(alternatives, list) else []):
            if not isinstance(alternative, dict):
                continue
            item = dict(alternative)
            nested = item.get("evidence")
            if isinstance(nested, dict):
                merged = dict(nested)
                merged.update(item)
                item = merged
            item["claim_name"] = f"{parent.get('claim_name', 'evidence')}/alternative_{index}"
            item["parent_claim_name"] = parent.get("claim_name", "evidence")
            item["is_alternate_hypothesis"] = True
            out.append(item)
    for claim in list(out):
        selected=claim.get("selected_hypothesis", {})
        child=selected.get("strict_child", {}) if isinstance(selected,dict) else {}
        if child.get("coupling_potential_status")!="local_coupling_potential_uncalibrated":continue
        item={**child.get("evidence", {}),**child}
        item.update(claim_name=claim.get("claim_name", "pyridine")+"/selected_hypothesis",
                    conditional_on_chosen_members=True,
                    preferred_role_options=claim.get("preferred_role_options", []),
                    selection_risks=[*child.get('selection_risks', []), "partner_assignment_ambiguous"],
                    cross_anchor_conflicts=claim.get("cross_anchor_conflicts", []))
        out.append(item)
    return out


def _status(claim, *keys):
    for key in keys:
        value = claim.get(key)
        if isinstance(value, dict):
            value = value.get("status")
        if value not in (None, "", False):
            return str(value)
    return ""


def _positive(status):
    s = str(status).lower()
    return any(x in s for x in ("supported", "complete", "positive", "potential")) and not any(
        x in s for x in ("unresolved", "unsupported", "missing", "partial", "blocked", "ambiguous"))


def _flatten_strings(value):
    value = _json(value)
    if isinstance(value, dict):
        result = []
        for k, v in value.items():
            result.append(str(k))
            result.extend(_flatten_strings(v))
        return result
    if isinstance(value, (list, tuple, set)):
        result = []
        for v in value:
            result.extend(_flatten_strings(v))
        return result
    return [str(value)] if value not in (None, "") else []


def _blockers(claim):
    keys = ("blockers", "coupling_blockers", "module_competition_blockers",
            "relay_blockers", "unresolved_attributes", "identity_blockers")
    return [x for key in keys for x in _flatten_strings(claim.get(key))]


def _hard(blockers):
    words = ("missing", "partial", "incomplete", "cross_contig", "cross-contig",
             "cross_bin", "cross-bin", "ligand", "resource_overlap", "resource_reuse",
             "role_resource", "target_not_supported", "not_supported",
             "association_missing", "sequence_support_insufficient", "remote_replacement")
    result = []
    for b in blockers:
        low = b.lower()
        # Distributed relay/pairing uncertainty is deliberately a soft
        # coupling blocker: it must not erase a coherent core+module call.
        if ("not_supported" in low or "not-supported" in low) and "internal" not in low and "intrinsic" not in low and any(x in low for x in ("relay", "distributed", "pairing")):
            continue
        if any(w in low for w in words):
            result.append(b)
    return result


def _soft(blockers):
    words = ("competition", "ambiguous", "pairing", "relay", "distributed", "alternative", "cross_anchor")
    return [b for b in blockers if any(w in b.lower() for w in words)]


def _family(row):
    value = _text(row.get("family") or row.get("preferred_family"))
    for name in ("NiFe", "FeFe", "HmdII", "Hmd", "Fdh", "FHL"):
        if name.lower() in value.lower():
            return name
    label = _text(row.get("function_prediction") or row.get("coarse_function_prediction"))
    for name in ("NiFe", "FeFe", "HmdII", "Hmd"):
        if name.lower() in label.lower():
            return name
    return value or "unresolved"


def _screen_family(row, evidence):
    """Map only catalytic-specific discovery profiles to a family."""
    values = []
    screen = evidence.get("screen", {}) if isinstance(evidence, dict) else {}
    dna = evidence.get("dna_rescue", {}) if isinstance(evidence, dict) else {}
    for obj in (row, screen, dna):
        for key in ("family_candidates", "families", "GA_profiles", "profiles", "profile_evidence", "domain_evidence", "raw_signal_evidence"):
            value = obj.get(key, []) if isinstance(obj, dict) else []
            values.extend(_flatten_strings(value))
    text = " ".join(map(str, values)).lower()
    if "pf00374" in text or "nife_related" in text: return "NiFe"
    if "pf02906" in text or "fefe_related" in text: return "FeFe"
    if "pf03201" in text or "pf22616" in text or "hmd_related" in text: return "Hmd"
    return "unresolved"


def _same_resource_overlap(claim):
    if claim.get("same_protein_resource_overlap") is True or claim.get("resource_overlap") is True:
        return True
    roles = claim.get("roles", {})
    if not isinstance(roles, dict): return False
    intervals = {}
    for rows in roles.values():
        for member in rows if isinstance(rows, list) else []:
            if not isinstance(member, dict): continue
            pid = member.get("protein_id")
            resources = member.get("resources") or member.get("query_resource_intervals") or []
            if member.get("query_domain_start") is not None and member.get("query_domain_end") is not None:
                resources = [[member["query_domain_start"],member["query_domain_end"]]]
            if isinstance(resources, dict): resources = list(resources.values())
            for interval in resources if isinstance(resources, list) else []:
                if isinstance(interval, dict): interval = interval.get("interval") or interval.get("query_interval")
                if isinstance(interval, (list, tuple)) and len(interval) >= 2:
                    try: intervals.setdefault(pid, []).append((int(interval[0]), int(interval[1])))
                    except (TypeError, ValueError): pass
    for spans in intervals.values():
        for i, (a, b) in enumerate(spans):
            for c, d in spans[i + 1:]:
                if max(a, c) <= min(b, d): return True
    return False


def _hdr_roles(claim):
    if isinstance(claim.get("evidence"),dict):claim={**claim["evidence"],**claim}
    if claim.get("coupling_status")!="reference_like_H2_Hdr_coupling_hypothesis":return {},[]
    candidates=claim.get("role_candidates",{})
    local=claim.get("local_relay_HdrBC_hypotheses",[])
    selected=list(local[0]) if len(local)==1 else [rows[0] for role,rows in candidates.items() if role!="HdrA" and len(rows)==1]
    anchor=next((r["protein_id"] for r in selected if r.get("role")=="MvhA"),None)
    assignment=claim.get("global_function_assignment",{})
    eligible=assignment.get("admissible_HdrA_edges",{}).get(anchor)
    pool=[r for r in candidates.get("HdrA",[]) if eligible is None or r["protein_id"] in eligible]
    witness=assignment.get("witness",{}).get(anchor)
    chosen=next((r for r in pool if r["protein_id"]==witness),None) if witness else pool[0] if len(pool)==1 else None
    if chosen:selected.append({**chosen,"selection_basis":"capacity_witness" if witness else "unique_sequence_candidate"})
    roles={}
    for row in selected:roles.setdefault(row["role"],[]).append(row)
    return roles,pool


def _members(claim):
    if isinstance(claim.get("evidence"),dict):claim={**claim["evidence"],**claim}
    roles = claim.get("roles", {})
    if not roles:roles,_=_hdr_roles(claim)
    if not roles and claim.get("status")=="candidate_CO_linked_H2_pathway":
        pairs=claim.get("substrate_modules",[])
        if len(pairs)==1:
            selected={**(claim.get("core") or {}).get("roles",{}),**pairs[0].get("roles",{})}
            roles={role:[{"protein_id":pid}] for role,pid in selected.items()}
    if not roles and isinstance(claim.get("hypotheses"), list):
        candidates = [x for x in claim["hypotheses"] if isinstance(x, dict) and isinstance(x.get("roles"), list)]
        if len(candidates)==1:
            roles = {}
            for member in candidates[0]["roles"]:
                if isinstance(member, dict) and member.get("role"):
                    roles.setdefault(member["role"], []).append(member)
    if not roles and claim.get("claim_name") == "mechanistic" and isinstance(claim.get("role_evidence"), dict):
        result_mech = []
        for role, rows in claim["role_evidence"].items():
            for entry in (rows if isinstance(rows, list) else [rows]):
                if isinstance(entry, dict) and entry.get("protein_id"):
                    result_mech.append({"role": role, "protein_id": entry["protein_id"]})
        unique_mech = {}
        for member in result_mech:
            unique_mech.setdefault((member["role"], member["protein_id"]), member)
        return list(unique_mech.values())
    if not roles and isinstance(claim.get("evidence"), dict):
        roles = claim["evidence"].get("roles", {})
    result = []
    if isinstance(roles, dict):
        for role, rows in roles.items():
            rows = rows if isinstance(rows, list) else [rows]
            if len(rows) == 1 and isinstance(rows[0], dict):
                pid = rows[0].get("protein_id") or rows[0].get("id")
                if pid:
                    item = {"role": role, "protein_id": pid}
                    if rows[0].get("query_domain_start") is not None and rows[0].get("query_domain_end") is not None:
                        item["query_interval"] = [rows[0]["query_domain_start"], rows[0]["query_domain_end"]]
                    elif isinstance(rows[0].get("resources"),dict):
                        positions=rows[0]["resources"].get(pid,[])
                        if positions:item["query_interval"]=[min(positions),max(positions)]
                    if rows[0].get("partial") is True:item["source_partial"]=True
                    if rows[0].get("selection_basis"):item["selection_basis"]=rows[0]["selection_basis"]
                    result.append(item)
    audit = claim.get("role_audit", {})
    if not audit and isinstance(claim.get("evidence"), dict):
        audit = claim["evidence"].get("role_audit", {})
    if isinstance(audit, dict):
        for role, info in audit.items():
            if not isinstance(info, dict):
                continue
            ids = info.get("selected_protein_ids") or []
            if isinstance(ids, str): ids = [ids]
            for pid in ids:
                item = {"role": role, "protein_id": pid}
                if item not in result: result.append(item)
    options = claim.get("preferred_role_options", [])
    for option in options if isinstance(options, list) else []:
        if not isinstance(option, dict): continue
        role = option.get("role")
        ids = option.get("selected_protein_ids") or []
        if isinstance(ids, str): ids = [ids]
        for pid in ids:
            item = {"role": role, "protein_id": pid}
            if role and item not in result: result.append(item)
    nodes = claim.get("module_graph_nodes") or claim.get("complexmembers") or []
    if isinstance(nodes, dict): nodes = list(nodes.values())
    for node in nodes if isinstance(nodes, list) else []:
        if not isinstance(node, dict): continue
        pid = node.get("protein_id") or node.get("candidate_id") or node.get("id")
        role = node.get("role") or node.get("member_role")
        if pid and role and {"role": role, "protein_id": pid} not in result:
            result.append({"role": role, "protein_id": pid})
    unique={}
    for member in result:unique.setdefault((member["role"],member["protein_id"]),member)
    return list(unique.values())


def _alternatives(claim, evidence):
    out = []
    hdr_roles,pool=_hdr_roles(claim)
    if len(pool)>1:
        chosen=[m["protein_id"] for m in hdr_roles.get("HdrA",[])]
        out.append({"role":"HdrA","selected_protein_ids":chosen,
                    "candidates":[{"protein_ids":[r["protein_id"]],"selected":r["protein_id"] in chosen} for r in pool],
                    "scope":"sequence_compatible_choices; capacity witness is not physical pairing"})
    for key in ("alternative_architectures", "alternatives", "competitive_alternatives", "preferred_role_options"):
        value = claim.get(key, evidence.get(key))
        value = _json(value)
        if isinstance(value, list): out.extend(value)
        elif value: out.append(value)
    audit = claim.get("role_audit", {})
    if not audit and isinstance(claim.get("evidence"), dict):
        audit = claim["evidence"].get("role_audit", {})
    if isinstance(audit, dict):
        for role, info in audit.items():
            if not isinstance(info, dict): continue
            ids = (info.get("competitive_alternatives") or info.get("alternatives") or
                   info.get("alternative_protein_ids") or info.get("candidate_protein_ids") or
                   info.get("candidate_ids"))
            if ids: out.append({"role": role, "alternatives": ids})
            elif info.get("candidate_count", 0) not in (0, 1):
                out.append({"role": role, "candidate_count": info.get("candidate_count")})
    compact = []
    for value in out:
        if isinstance(value, dict):
            # Alternative architecture records can contain full alignments;
            # expose only the auditable comparison identity and statuses.
            keep = {k: value[k] for k in ("structure_id", "claim_name", "status",
                    "hydrogen_reaction_core_status", "pyridine_nucleotide_module_status",
                    "coupling_potential_status", "role", "protein_id", "alternatives",
                    "candidate_count", "selected_protein_ids", "candidates", "scope") if k in value}
            members = _members(value)
            if members: keep["members"] = members
            compact.append(keep or value)
        else:
            compact.append(value)
    return compact


def _unknowns(row):
    return {
        "NAD_NADP_selectivity": row.get("NAD_NADP_selectivity") or "unknown",
        "bifurcation": "unknown", "electron_acceptor_or_donor_specificity": "unknown",
        "expression": "unknown", "cofactor_occupancy": "unknown",
        "physical_assembly": "unknown", "query_activity": row.get("query_activity") or "not_tested",
        "in_vivo_direction": row.get("in_vivo_direction") or "unknown",
    }


def _positive_hmdii(claims, row):
    if _text(row.get("family")).strip().lower() in {"hmdii", "hmd2", "hmd-ii"}:
        return True
    label = _text(row.get("coarse_function_prediction") or row.get("function_prediction"))
    if "candidate_H2_THF_derivative_hydrogenation_potential" in label:
        return True
    for claim in claims:
        status = _status(claim, "hmdii_status", "HmdII_status")
        if _positive(status):
            return True
    return False


def _positive_legacy(label):
    return label in {
        "candidate_H2_F420_reduction", "candidate_H2_Hdr_coupling",
        "candidate_H2_ferredoxin_compatible_core", "candidate_H2_quinone_coupling", "candidate_H2_cytochrome_c3_coupling",
        "candidate_H2_formate_coupling", "candidate_H2_pyridine_nucleotide_coupling",
        "candidate_H2_H4MPT_hydrogenation_potential", "candidate_H2_THF_derivative_hydrogenation_potential",
        "candidate_CO_linked_H2_pathway", "candidate_H2_sensor_regulation", "H2_regulatory_signal",
        "candidate_membrane_associated_H2_redox_system", "candidate_genome_H2_formate_pathway",
        "candidate_H2_formate_pathway_conditional_on_recoding",
        "reference_sequence_matched_reaction_potential", "reference_coding_locus_matched_reaction_potential",
        "candidate_reference_reaction_potential",
    }


def _legacy_summary(label, confidence):
    names = {
        "candidate_H2_F420_reduction": "F420还原耦联氢酶",
        "candidate_H2_Hdr_coupling": "Hdr耦联氢酶",
        "candidate_H2_ferredoxin_compatible_core": "铁氧还蛋白相容氢酶",
        "candidate_H2_quinone_coupling": "醌耦联氢酶",
        "candidate_H2_cytochrome_c3_coupling": "细胞色素c3耦联氢酶",
        "candidate_H2_formate_coupling": "甲酸耦联氢酶",
        "candidate_H2_pyridine_nucleotide_coupling": "核苷酸耦联氢酶",
        "candidate_H2_H4MPT_hydrogenation_potential": "H4MPT加氢氢酶",
        "candidate_H2_THF_derivative_hydrogenation_potential": "四氢叶酸衍生物加氢氢酶",
        "candidate_membrane_associated_H2_redox_system": "膜关联氢酶氧化还原系统",
        "candidate_CO_linked_H2_pathway": "CO关联氢代谢途径",
        "candidate_H2_sensor_regulation": "调控/感测型氢酶",
        "H2_regulatory_signal": "调控/感测型氢酶",
        "candidate_genome_H2_formate_pathway": "基因组支持的氢/甲酸代谢途径",
        "candidate_H2_formate_pathway_conditional_on_recoding": "依赖重编码的氢/甲酸代谢途径",
    }
    return f"预测为{names.get(label, label)}，{confidence}置信度；具体受体选择性和实验活性仍未确定。"


def _matching_claim(label, claims):
    branches = {
        "candidate_H2_pyridine_nucleotide_coupling": ("pyridine", "composition"),
        "candidate_H2_F420_reduction": ("frh",), "candidate_H2_Hdr_coupling": ("mvh",),
        "candidate_H2_quinone_coupling": ("huc", "quinone"),
        "candidate_H2_formate_coupling": ("hdcr", "fhl", "tonn"),
        "candidate_H2_ferredoxin_compatible_core": ("mbh", "ech"),
        "candidate_membrane_associated_H2_redox_system": ("ech", "membrane_context"),
        "candidate_CO_linked_H2_pathway": ("coo",),
        "candidate_genome_H2_formate_pathway": ("fhl",),
        "candidate_H2_formate_pathway_conditional_on_recoding": ("fhl",),
        "candidate_H2_H4MPT_hydrogenation_potential": ("canonical_hmd",),
    }
    for branch in branches.get(label, ()):
        for claim in claims:
            if claim.get("claim_name") == branch and _claim_supports(label, claim): return claim
    for claim in claims:
        if claim.get("claim_name") == "mechanistic" and claim.get("function_prediction") == label:
            return claim
    return None


def _claim_supports(label, claim):
    fields={
        "candidate_H2_pyridine_nucleotide_coupling": ("coupling_potential_status","status"),
        "candidate_H2_quinone_coupling": ("coupling_potential_status",),
        "candidate_H2_formate_coupling": ("coupling_status",),
        "candidate_H2_F420_reduction": ("F420_status",),
        "candidate_H2_Hdr_coupling": ("coupling_status",),
        "candidate_H2_ferredoxin_compatible_core": ("core_coupling_status","coupling_potential_status"),
        "candidate_membrane_associated_H2_redox_system": ("membrane_redox_status",),
        "candidate_CO_linked_H2_pathway": ("status",),
        "candidate_genome_H2_formate_pathway": ("coupling_status",),
        "candidate_H2_formate_pathway_conditional_on_recoding": ("coding_pathway_status",),
        "candidate_H2_H4MPT_hydrogenation_potential": ("reaction_status","branch_reaction_status"),
    }
    expected={
        "candidate_H2_pyridine_nucleotide_coupling":"local_coupling_potential_uncalibrated",
        "candidate_H2_quinone_coupling":"local_quinone_coupling_potential_uncalibrated",
        "candidate_H2_F420_reduction":"compatible_reaction_candidate_uncalibrated",
        "candidate_H2_Hdr_coupling":"reference_like_H2_Hdr_coupling_hypothesis",
        "candidate_H2_ferredoxin_compatible_core":"reference_like_H2_Fd_coupling_potential",
        "candidate_membrane_associated_H2_redox_system":"local_membrane_H2_redox_potential_uncalibrated",
        "candidate_H2_H4MPT_hydrogenation_potential":"compatible_reaction_candidate_uncalibrated",
    }.get(label,label)
    return any(claim.get(field)==expected for field in fields.get(label,()))


def _claim_record(claim):
    if not isinstance(claim, dict): return claim
    record = {"structure_id": claim.get("structure_id") or claim.get("claim_name", "unknown"),
              "claim_name": claim.get("claim_name", "unknown"),
              "core_status": _status(claim, "hydrogen_reaction_core_status", "core_status"),
              "module_status": _status(claim, "pyridine_nucleotide_module_status", "module_status"),
              "coupling_status": _status(claim, "coupling_potential_status", "coupling_status", "status")}
    members = _members(claim)
    if members: record["members"] = members
    return record


_REACTION_PARTNERS = {
    "candidate_H2_pyridine_nucleotide_coupling": ("核苷酸", ("nad", "pyridine_nucleotide")),
    "candidate_H2_ferredoxin_compatible_core": ("铁氧还蛋白", ("ferredoxin", "fdred", "fdox")),
    "candidate_H2_quinone_coupling": ("醌", ("quinone", "menaquin", "ubiquin")),
    "candidate_H2_F420_reduction": ("F420", ("f420",)),
    "candidate_H2_Hdr_coupling": ("异二硫键", ("heterodisulfide", "com_s_s_cob")),
    "candidate_H2_formate_coupling": ("甲酸", ("formate",)),
    "candidate_H2_H4MPT_hydrogenation_potential": ("H4MPT", ("h4mpt", "methenyl_h4mpt", "methylene_h4mpt")),
}


def _capability_matches(label, capabilities):
    words = _REACTION_PARTNERS.get(label, ("", ()))[1]
    # Only typed reaction participants count; prose can mention untested or
    # excluded partners and must not supply positive evidence.
    return any(str(p).lower().startswith(words)
               for c in capabilities if isinstance(c, dict)
               for key in ("electron_donors", "electron_acceptors")
               for p in c.get(key, [])) if isinstance(capabilities, list) else False


def _reference_reaction_summary(capabilities):
    """Describe only typed participants; never infer from assay prose."""
    names = {"ferredoxin_oxidized": "Fd氧化态", "ferredoxin_reduced": "Fd还原态",
             "protons": "H+", "NAD+": "NAD+", "NAD_oxidized": "NAD+",
             "NADP_oxidized": "NADP+", "NADP+": "NADP+", "NADH": "NADH",
             "NADPH": "NADPH", "H2": "H2",
             "FMN_oxidized": "FMN", "FAD_oxidized": "FAD",
             "methyl_viologen_reduced": "甲基紫精还原态",
             "methyl_viologen_oxidized": "甲基紫精氧化态"}
    pieces = []
    for capability in capabilities if isinstance(capabilities, list) else []:
        if not isinstance(capability, dict):
            continue
        raw_donors = [str(x) for x in capability.get("electron_donors", []) if x]
        raw_acceptors = [str(x) for x in capability.get("electron_acceptors", []) if x]
        donors = [names.get(x, x) for x in raw_donors]
        acceptors = [names.get(x, x) for x in raw_acceptors]
        role = capability.get("hydrogen_role")
        if not donors and not acceptors:
            continue
        donor_text = "、".join(donors)
        acceptor_text = "、".join(acceptors)
        if role == "oxidation":
            text = f"{donor_text}供电子"
            if acceptor_text: text += f"，还原{acceptor_text}"
        elif role == "evolution":
            text = f"{donor_text}供电子产H2"
            if raw_acceptors and raw_acceptors != ["protons"]: text += f"，受体为{acceptor_text}"
        else:
            text = f"{donor_text}供电子，还原{acceptor_text}"
        if capability.get("direct_acceptor_status") == "reference_reported_coupled":
            text += "（耦联测定）"
        pieces.append(text)
    return "；".join(pieces)


def _artificial_mediator_capability(capability):
    """Accept only explicit typed artificial-mediator H2 chemistry."""
    if not isinstance(capability, dict) or capability.get("reference_evidence_scope") != _ARTIFICIAL_MEDIATOR_SCOPE:
        return False
    donors = {str(x).lower() for x in capability.get("electron_donors", []) if x}
    acceptors = {str(x).lower() for x in capability.get("electron_acceptors", []) if x}
    role = capability.get("hydrogen_role")
    if not donors or not acceptors:
        return False
    return ((role == "evolution" and "protons" in acceptors)
            or (role == "oxidation" and "h2" in donors))


def _reference_override(row, claims, prediction):
    reference = next((c for c in claims if c.get("claim_name")=="reference"), {})
    status = reference.get("status")
    if status not in {"exact_reference_sequence_potential", "exact_reference_coding_locus_potential"}:
        return None
    if status != row.get("reference_reaction_status") or reference.get("blockers"):
        return None
    caps = reference.get("capabilities", [])
    supported = [label for label in _REACTION_PARTNERS if _capability_matches(label, caps)]
    artificial_caps = [c for c in caps if _artificial_mediator_capability(c)]
    if not supported and not artificial_caps:
        return None
    artificial_mediator = bool(artificial_caps) and not supported
    if artificial_mediator and row.get("identity_status") not in {
            "calibrated_high", "sequence_supported_uncalibrated_core",
            "sequence_supported_uncalibrated_noncanonical"}:
        return None
    members = reference.get("tuple_matches", [])
    anchor = row.get("candidate_id")
    if not anchor or not members or anchor not in {m.get("protein_id") for m in members}:
        return None
    selected = [h for h in reference.get("tuple_hypotheses", [])
                if h.get("structure_id")==reference.get("matched_structure_id") and not h.get("blockers")]
    binding = lambda ms: sorted((m.get("role", ""), m.get("protein_id", "")) for m in ms)
    if len(selected)!=1 or binding(selected[0].get("matched_members", []))!=binding(members):
        return None
    ref_label = ("reference_coding_locus_matched_reaction_potential" if status=="exact_reference_coding_locus_potential"
                 else "reference_sequence_matched_reaction_potential")
    source = "完整参考编码位点" if status=="exact_reference_coding_locus_potential" else "完整参考成员序列"
    partners = "、".join(_REACTION_PARTNERS[label][0] for label in supported)
    chosen = [{"role":m["role"], "protein_id":m["protein_id"],
               **({"source_partial":True} if m.get("raw_partial") else {}),
               **({"coding_locus_proof_used":True} if m.get("coding_locus_proof_used") else {})} for m in members]
    alternate = {k:prediction[k] for k in ("preferred_function", "prediction_confidence", "preferred_members", "prediction_support")}
    alternate["scope"] = "sequence_hypothesis; not an assay-negative claim for this partner"
    reference_support = {"claim":"reference_reaction", "status":status,
                         "structure_id":reference.get("matched_structure_id"),
                         "reaction_keys":[c.get("reaction_key") for c in caps if c.get("reaction_key")],
                         "contract_ids":[c.get("contract_id") for c in caps if c.get("contract_id")],
                         "scope":reference.get("claim_scope", reference.get("scope")),
                         "sources":list(dict.fromkeys(c.get("primary_source", "") for c in caps if c.get("primary_source")))}
    ambiguities = reference.get("member_selection_ambiguities", [])
    equivalent_alternatives = [{"role":a["role"], "selected_protein_ids":[a["selected_protein_id"]],
                               "candidates":[{"protein_ids":[pid], "selected":pid==a["selected_protein_id"]}
                                             for pid in [a["selected_protein_id"], *a["alternative_protein_ids"]]],
                               "scope":"sequence_equivalent_member_loci; context ranking does not prove physical pairing"}
                              for a in ambiguities]
    ambiguity_risks = ["sequence_equivalent_partner_locus_ambiguous"] if ambiguities else []
    ambiguity_summary = "伙伴位点存在同序列备选，首选按上下文排序，实际装配未定。" if ambiguities else ""
    if ambiguities:reference_support["member_selection_ambiguities"] = ambiguities
    if prediction["preferred_function"] in supported:
        exact_members = chosen
        previous_binding = {m.get("protein_id") for m in prediction.get("preferred_members", [])}
        exact_binding = {m.get("protein_id") for m in exact_members}
        alternatives = list(prediction.get("prediction_alternatives", []))
        if previous_binding != exact_binding:
            alternatives.insert(0, alternate)
        alternatives.extend(equivalent_alternatives)
        reaction = _reference_reaction_summary(caps)
        summary = prediction["prediction_summary"]
        confidence = "medium" if ambiguities else prediction["prediction_confidence"]
        if reaction:
            lead = summary.split('；', 1)[0]
            if ambiguities:lead = lead.replace("高置信度", "中等置信度").replace("低置信度", "中等置信度")
            summary = f"{lead}；参考条件支持的具体反应潜力：{reaction}。{ambiguity_summary}查询样本实测活性、天然伙伴专一性和体内方向仍未知。"
        return _result(prediction["preferred_function"], prediction["preferred_family"], prediction["prediction_kind"], confidence,
                       summary, [*prediction.get("prediction_support", []), reference_support], exact_members,
                       alternatives, [*prediction.get("prediction_risks", []), *ambiguity_risks, "reference_assay_conditions_only", "query_activity_not_tested"], row)
    confidence = "medium" if ambiguities else "high"
    if artificial_mediator:
        role = artificial_caps[0].get("hydrogen_role")
        direction = "产H2" if role == "evolution" else "H2氧化"
        reaction = _reference_reaction_summary(artificial_caps)
        return _result(ref_label, prediction["preferred_family"], "reference_function", confidence,
                       f"首选为人工介质条件下的{direction}参考反应潜力，{'中等' if ambiguities else '高'}置信度；{source}与来源测定支持。具体反应：{reaction}。{ambiguity_summary}天然伙伴、查询样本活性及体内方向仍未知。",
                       [reference_support], chosen, [alternate, *prediction["prediction_alternatives"], *equivalent_alternatives],
                       [*ambiguity_risks, "artificial_mediator_not_native_coupling", "reference_assay_conditions_only",
                        "query_activity_not_tested", "native_partner_specificity_unknown"], row)
    return _result(ref_label, prediction["preferred_family"], "reference_function", confidence,
                   f"首选预测为{partners}电子转移反应潜力，{'中等' if ambiguities else '高'}置信度；{source}与来源测定支持。{ambiguity_summary}天然伙伴、专一性及查询样本实测活性仍未知。",
                   [reference_support],
                   chosen, [alternate, *prediction["prediction_alternatives"], *equivalent_alternatives],
                   [*ambiguity_risks, "reference_assay_conditions_only", "query_activity_not_tested", "native_partner_specificity_unknown"], row)


def _single_chain_reference_prediction(row, claims, prediction, evidence):
    """Consume an assigned single-chain Fd pathway without claiming a query assay."""
    if prediction['preferred_function'] not in {'hydrogenase_core_candidate','possible_hydrogenase','unresolved_candidate'}:
        return None
    if row.get('identity_status') not in {'calibrated_high','sequence_supported_uncalibrated_core'}:return None
    reference=next((c for c in claims if c.get('claim_name')=='reference'),{})
    if reference.get('status')!='inherited_reference_potential' or row.get('reference_reaction_status')!=reference['status']:return None
    architectures={a.get('value') for a in evidence.get('attributes',[]) if a.get('attribute')=='candidate_complex_architecture' and a.get('status')=='complete_supported'}
    if len(architectures)!=1:return None
    sid=next(iter(architectures));nodes=evidence.get('module_graph_nodes',[]);anchor=row.get('candidate_id')
    if not anchor or not isinstance(nodes,list) or len(nodes)!=1 or nodes[0].get('selected_protein_ids')!=[anchor]:return None
    selected=[c for c in nodes[0].get('candidates',[]) if c.get('selected') is True]
    if (len(selected)!=1 or selected[0].get('protein_ids')!=[anchor] or selected[0].get('fragmented') is not False
            or selected[0].get('jointly_feasible') is not True):return None
    assigned={a.get('contract_id') for a in reference.get('inherited_assignments',[]) if a.get('status') in {'homology_supported_potential','reference_sequence_supported_potential'}}
    caps=[c for c in reference.get('capabilities',[]) if c.get('structure_id')==sid and c.get('contract_id') in assigned]
    if not caps or any(c.get('scope')!='reference_assay_conditions_only' or c.get('ferredoxin_requirement')!='direct_electron_partner'
                       or c.get('electron_coupling')!='single_reported_electron_pathway' for c in caps):return None
    for c in caps:
        if not ((c.get('hydrogen_role')=='evolution' and 'ferredoxin_reduced' in c.get('electron_donors',[]) and 'protons' in c.get('electron_acceptors',[]))
                or (c.get('hydrogen_role')=='oxidation' and 'H2' in c.get('electron_donors',[]) and 'ferredoxin_oxidized' in c.get('electron_acceptors',[]))):return None
    if any(_capability_matches(label,caps) for label in _REACTION_PARTNERS if label!='candidate_H2_ferredoxin_compatible_core'):return None
    confidence='low' if prediction['prediction_confidence']=='low' else 'medium'
    support=dict(claim='single_chain_reference_function',structure_id=sid,contract_ids=[c['contract_id'] for c in caps],
                 sources=list(dict.fromkeys(c.get('primary_source','') for c in caps)),scope='complete_single_chain_reference_homology_not_query_assay')
    alternative={k:prediction[k] for k in ('preferred_function','prediction_confidence','preferred_members','prediction_support')}
    risks=[r for r in prediction['prediction_risks'] if r!='partner_or_electron_module_unresolved']
    return _result('candidate_H2_ferredoxin_compatible_core',prediction['preferred_family'],'coarse_function',confidence,
                   f"预测为铁氧还蛋白相容氢酶，{'低' if confidence=='low' else '中等'}置信度；独立催化身份、完整单链参考架构及已赋值Fd反应潜力支持。实际Fd伙伴、选择性、查询活性及体内方向未定。",
                   [*prediction['prediction_support'],support],[dict(role=nodes[0]['role'],protein_id=anchor)],
                   [alternative,*prediction['prediction_alternatives']],
                   [*risks,'reference_homology_function_uncalibrated','native_ferredoxin_partner_unidentified','query_activity_not_tested'],row)


def _component_prediction(row, claims, prediction):
    """Locality limits source-tuple attribution, not complete chemical components."""
    if row.get("identity_status") not in {"calibrated_high", "sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical"}:
        return None
    reference = next((c for c in claims if c.get("claim_name")=="reference"), {})
    options = []
    for hypothesis in reference.get("sequence_component_hypotheses", []):
        linkage = hypothesis.get("linkage_blockers", [])
        if (hypothesis.get("status")!="source_equivalent_components_potential" or hypothesis.get("hard_blockers")
                or not hypothesis.get("strict_tuple_blockers")
                or any(b!="whole_reference_tuple_local_linkage_missing" and not b.startswith("reference_member_group_local_linkage_missing_") for b in linkage)):
            continue
        members = hypothesis.get("matched_members", [])
        if row.get("candidate_id") not in {m.get("protein_id") for m in members}:
            continue
        caps = hypothesis.get("capabilities", [])
        if not caps or any(c.get("scope")!="sequence_component_function_hypothesis" for c in caps):
            continue
        supported = [label for label in _REACTION_PARTNERS if _capability_matches(label, caps)]
        if not supported:
            continue
        label = prediction["preferred_function"] if prediction["preferred_function"] in supported else supported[0]
        chosen = [{"role":m["role"], "protein_id":m["protein_id"],
                   **({"source_partial":True} if m.get("raw_partial") else {}),
                   **({"coding_locus_proof_used":True} if m.get("coding_locus_proof_used") else {})} for m in members]
        ambiguities = hypothesis.get("member_selection_ambiguities", [])
        tied_roles = {a["role"] for a in ambiguities if a.get("equally_supported_loci")}
        chosen = [m for m in chosen if m["role"] not in tied_roles]
        alternatives = [{"role":a["role"], "selected_protein_ids":[] if a["role"] in tied_roles else [a["selected_protein_id"]],
                         "candidates":[{"protein_ids":[pid], "selected":a["role"] not in tied_roles and pid==a["selected_protein_id"]}
                                       for pid in [a["selected_protein_id"], *a["alternative_protein_ids"]]],
                         "scope":"sequence_equivalent_member_loci; context ranking does not prove physical pairing"} for a in ambiguities]
        support = {"claim":"source_equivalent_component_function", "status":hypothesis["status"],
                   "structure_id":hypothesis.get("structure_id"), "linkage_blockers":linkage,
                   "strict_tuple_blockers":hypothesis["strict_tuple_blockers"],
                   "reaction_keys":[c.get("reaction_key") for c in caps if c.get("reaction_key")],
                   "contract_ids":[c.get("contract_id") for c in caps if c.get("contract_id")],
                   "scope":"source_equivalent_components_function_hypothesis",
                   "sources":list(dict.fromkeys(c.get("primary_source", "") for c in caps if c.get("primary_source")))}
        association_risk = "source_assay_tuple_local_linkage_unresolved" if linkage else "source_assay_tuple_member_assignment_unresolved"
        association_summary = "参考基因组织未闭合" if linkage else "原参考成员关联仍受限"
        risks = [association_risk, "physical_assembly_unknown", "query_activity_not_tested"]
        if ambiguities:risks.append("sequence_equivalent_partner_locus_ambiguous")
        if tied_roles:risks.append("partner_assignment_ambiguous")
        result = _result(label, prediction["preferred_family"], "coarse_function", "medium",
                         f"{_legacy_summary(label, '中等').split('；', 1)[0]}；完整来源同序列组件支持，{association_summary}；组件证据支持的反应潜力：{_reference_reaction_summary(caps)}。实际配对、表达、辅因子占据与体内净方向未定。",
                         [support], chosen, alternatives, risks, row)
        options.append((label!=prediction["preferred_function"], len(linkage), len(ambiguities), hypothesis.get("structure_id", ""), result))
    if not options:
        return None
    # A complete source-equivalent explanation is not a ranking tie-breaker
    # for another complete explanation.  In particular, two hypotheses with
    # the same structure, linkage evidence and ambiguity burden used to be
    # resolved by list insertion order (R164 Syntrophus probe).  Keep the
    # ordering deterministic for diagnostics, but defer the member assignment
    # when the evidence score is genuinely tied.
    options.sort(key=lambda option: (
        *option[:4],
        tuple(sorted((m.get("role", ""), m.get("protein_id", ""))
                     for m in option[-1].get("preferred_members", []))),
    ))
    result=options[0][-1]
    tied = [option for option in options if option[:3] == options[0][:3]]
    bindings = {
        tuple(sorted((m.get("role", ""), m.get("protein_id", ""))
                     for m in option[-1].get("preferred_members", [])))
        for option in tied
    }
    if len(tied) > 1 and len(bindings) > 1:
        # Preserve only member assignments that every tied explanation agrees
        # on.  The coarse function remains valid because all retained options
        # passed the same capability checks; only the physical partner choice
        # is unresolved.
        member_sets = [
            {(m.get("role", ""), m.get("protein_id", "")): m
             for m in option[-1].get("preferred_members", [])}
            for option in tied
        ]
        common = sorted(set.intersection(*(set(s) for s in member_sets)))
        result = dict(result)
        result["preferred_members"] = [member_sets[0][key] for key in common]
        result["prediction_summary"] += " 同等证据的伙伴解释并列，具体成员配对未定。"
        result["prediction_risks"] = [
            *result.get("prediction_risks", []), "partner_assignment_ambiguous"
        ]
        # Keep every tied explanation available for audit, in deterministic
        # order, rather than silently selecting the first input row.
        result["prediction_alternatives"] = [
            *result.get("prediction_alternatives", []),
            *[option[-1] for option in tied[1:]],
        ]
    matching=_matching_claim(prediction['preferred_function'],claims)
    typed_ids={m['protein_id'] for m in _members(matching)} if matching else set()
    current_ids={m['protein_id'] for m in prediction['preferred_members']}
    alternative_ids={m['protein_id'] for m in result['preferred_members']}
    if (prediction['preferred_function']==result['preferred_function'] and len(typed_ids)>1
            and current_ids==typed_ids):
        support=dict(claim='supported_partner_combination_retained',evidence_source='typed_claim:'+matching['claim_name'],
                     alternative_member_ids=sorted(alternative_ids),
                     scope='Same coarse function; source-equivalent alternatives do not replace independently supported members')
        return {**prediction,
                'prediction_summary':prediction['prediction_summary']+' 已有角色组合及关联证据优先，来源组件的条件性解释保留为备选。',
                'prediction_support':[*prediction['prediction_support'],support],
                'prediction_alternatives':[*prediction['prediction_alternatives'],result,*[o[-1] for o in options[1:]]]}
    alternatives = [{k:other[k] for k in ("preferred_function","prediction_confidence","preferred_members","prediction_support","prediction_risks")}
                    for other in [prediction, *[option[-1] for option in options[1:]]]]
    if prediction["preferred_function"] in _REACTION_PARTNERS and prediction["preferred_function"]!=result["preferred_function"]:
        return {**prediction, "prediction_alternatives":[*prediction["prediction_alternatives"], result]}
    result["prediction_alternatives"].extend([*alternatives, *prediction["prediction_alternatives"]])
    return result


def _source_like_module_prediction(row, reference, prediction):
    """Consume typed source-like module evidence without transferring assay details."""
    if row.get("identity_status") not in {"calibrated_high", "sequence_supported_uncalibrated_core",
                                          "sequence_supported_uncalibrated_noncanonical"}:
        return None
    if any(isinstance(s, dict) and s.get("claim") in {"reference_reaction", "source_equivalent_component_function"}
           for s in prediction["prediction_support"]):
        return None
    records = []
    anchor = row.get("candidate_id")
    for hypothesis in reference.get("source_like_component_hypotheses", []):
        if (hypothesis.get("status") != "source_like_module_component_potential"
                or hypothesis.get("scope") != "source_like_module_function_prior"
                or hypothesis.get("hard_blockers")
                or hypothesis.get("reference_reaction_details_transferred") is not False
                or hypothesis.get("source_contracts_transferred") is not False
                or hypothesis.get("capabilities") != []
                or not isinstance(hypothesis.get("prefer_source_function"), bool)):
            continue
        members = hypothesis.get("matched_members", [])
        if not anchor or anchor not in {m.get("protein_id") for m in members}:
            continue
        homologs = [m for m in members if m.get("matched_by") == "supported_local_module_homology"]
        if len(homologs) != 1:
            continue
        homolog = homologs[0]
        alignment = homolog.get("source_alignment")
        dominates = homolog.get("source_alignment_dominates_role_template")
        if (homolog.get("source_equivalent") is not False or not isinstance(alignment, dict)
                or not isinstance(dominates, bool) or dominates != hypothesis["prefer_source_function"]
                or not alignment.get("query_sha256") or not alignment.get("reference_sha256")
                or alignment.get("query_sha256") != homolog.get("query_sha256")
                or alignment.get("reference_sha256") != homolog.get("reference_sha256")
                or alignment.get("query_sha256") == alignment.get("reference_sha256")
                or not all(alignment.get(key) is not None for key in ("identity", "query_coverage", "reference_coverage"))
                or any(m is not homolog and m.get("source_equivalent") is not True for m in members)):
            continue
        participants = hypothesis.get("source_participants", [])
        if isinstance(participants, list) and all(isinstance(item, str) for item in participants):
            typed = {"H2", "formate", "protons", "NADH", "NADPH", "NAD+", "NADP+",
                     "NAD_oxidized", "NADP_oxidized", "ferredoxin_oxidized", "ferredoxin_reduced",
                     "methyl_viologen_reduced", "methyl_viologen_oxidized"}
            participants = [{"electron_donors": [item for item in participants if item in typed]}]
        else:
            continue
        labels = [label for label in _REACTION_PARTNERS if _capability_matches(label, participants)]
        if not labels:
            continue
        chosen = [{"role": m.get("role"), "protein_id": m.get("protein_id"),
                   **({"source_alignment": alignment} if m is homolog else {})} for m in members]
        support = {"claim": "source_like_module_component_function", "status": hypothesis["status"],
                   "scope": hypothesis["scope"], "source_functions": labels,
                   "prefer_source_function": hypothesis["prefer_source_function"],
                   "source_alignment": alignment,
                   "role_module_evidence": homolog.get("role_module_evidence", {}),
                   "source_sources": list(hypothesis.get("source_sources", []))}
        confidence = "medium" if hypothesis["prefer_source_function"] and hypothesis.get('confidence') != 'low' else "low"
        text = (_legacy_summary(labels[0], '中等' if confidence == 'medium' else '低').split('；', 1)[0]
                + "；来源近似成员组合和当前电子模块支持；非同序列伙伴功能转移未经独立验证，"
                  "天然伙伴、选择性、查询活性、方向未知。")
        risks = ["source_module_function_not_independently_validated", "source_reaction_details_not_transferred",
                 "query_activity_not_tested", "native_partner_specificity_unknown", *hypothesis.get('selection_risks', [])]
        records.append({"labels": labels, "prefer": hypothesis["prefer_source_function"],
                        "result": _result(labels[0], prediction["preferred_family"], "coarse_function", confidence,
                                           text, [support], chosen, [], risks, row)})
    if not records:
        return None
    source_labels = {label for record in records for label in record["labels"]}
    current = prediction["preferred_function"]
    current_supported = current in source_labels
    prefer = any(record["prefer"] for record in records)
    if current_supported or not (prefer and len(source_labels) == 1):
        result = prediction
        for record in records:
            result["prediction_support"].extend(record["result"]["prediction_support"])
        result["prediction_risks"].extend(risk for record in records for risk in record["result"]["prediction_risks"])
        result["prediction_alternatives"].extend(record["result"] for record in records if current not in record["labels"])
        return result
    chosen = next(record["result"] for record in records if record["prefer"] and len(record["labels"]) == 1)
    old = {key: prediction[key] for key in ("preferred_function", "prediction_confidence", "preferred_members", "prediction_support")}
    chosen["prediction_alternatives"].insert(0, old)
    chosen["prediction_alternatives"].extend(record["result"] for record in records if record["result"] is not chosen)
    return chosen


def _source_like_prediction(row, claims, prediction):
    """Add a coarse partner function without inheriting the source reaction details."""
    if row.get("identity_status") not in {"calibrated_high", "sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical"}:
        return None
    if any(isinstance(s,dict) and s.get("claim") in {"reference_reaction","source_equivalent_component_function"} for s in prediction["prediction_support"]):
        return None
    reference=next((c for c in claims if c.get("claim_name")=="reference"), {})
    module_result = _source_like_module_prediction(row, reference, prediction)
    if module_result is not None:
        return module_result
    options=[]
    for h in reference.get("source_like_component_hypotheses", []):
        if (h.get("status")!="source_like_formate_component_potential" or h.get("hard_blockers")
                or h.get("scope")!="source_like_FDH_role_combination" or h.get("reference_reaction_details_transferred") is not False
                or "candidate_H2_formate_coupling" not in h.get("supported_functions", [])):
            continue
        members=h.get("matched_members", []);homologs=[m for m in members if m.get("matched_by")=="supported_FDH_role_homology"]
        if row.get("candidate_id") not in {m.get("protein_id") for m in members} or len(homologs)!=1:
            continue
        fdh=homologs[0]
        if fdh.get("source_equivalent") is not False or fdh.get("query_sha256")==fdh.get("reference_sha256"):
            continue
        confidence="low" if h.get("confidence")=="low" or h.get("linkage_risks") else "medium"
        support={"claim":"source_like_FDH_role_combination", "status":h["status"], "structure_id":h.get("structure_id"),
                 "supported_functions":h["supported_functions"], "confidence":confidence,
                 "homologous_roles":h.get("homologous_roles", []), "reference_reaction_details_transferred":False,
                 "paired_identity":fdh.get("paired_identity"), "reference_coverage":fdh.get("reference_coverage"),
                 "scope":h["scope"]}
        chosen=[{"role":m["role"],"protein_id":m["protein_id"],
                 **({"source_partial":True} if m.get("raw_partial") else {}),
                 **({"coding_locus_proof_used":True} if m.get("coding_locus_proof_used") else {}),
                 **({"matched_by":m["matched_by"],"source_equivalent":False} if m.get("matched_by")=="supported_FDH_role_homology" else {})} for m in members]
        alternatives=[{"role":a["role"],"selected_protein_ids":[a["selected_protein_id"]],
                       "candidates":[{"protein_ids":[pid],"selected":pid==a["selected_protein_id"]}
                                     for pid in [a["selected_protein_id"],*a["alternative_protein_ids"]]],
                       "scope":"source_like_role_alternatives; pairing not measured"} for a in h.get("member_selection_ambiguities", [])]
        risks=["FDH_partner_inferred_from_homology","source_reaction_details_not_transferred","query_activity_not_tested",*h.get("linkage_risks", [])]
        text=f"预测具有甲酸耦联潜力，{'低' if confidence=='low' else '中等'}置信度；完整来源核心与电子模块、FDH样伙伴的域和保守位点及竞争证据支持。具体反应方向、联合机制与实际装配未定。"
        result=_result("candidate_H2_formate_coupling",prediction["preferred_family"],"coarse_function",confidence,text,[support],chosen,alternatives,risks,row)
        options.append((confidence!="medium",-float(fdh.get("paired_identity",0)),h.get("structure_id", ""),result))
    if not options:
        return None
    options.sort(key=lambda option:option[:3]);result=options[0][-1]
    old={k:prediction[k] for k in ("preferred_function","prediction_confidence","preferred_members","prediction_support")}
    if prediction["preferred_function"]=="candidate_H2_pyridine_nucleotide_coupling":
        return {**prediction,"prediction_summary":prediction["prediction_summary"].split('；',1)[0]+"；同时"+result["prediction_summary"],
                "prediction_support":[*prediction["prediction_support"],*result["prediction_support"]],"preferred_members":result["preferred_members"],
                "prediction_alternatives":[*prediction["prediction_alternatives"],old,result],
                "prediction_risks":[*prediction["prediction_risks"],*result["prediction_risks"]]}
    if prediction["preferred_function"] in _REACTION_PARTNERS and prediction["preferred_function"]!=result["preferred_function"]:
        return {**prediction,"prediction_alternatives":[*prediction["prediction_alternatives"],result]}
    result["prediction_alternatives"].extend([old,*prediction["prediction_alternatives"],*[option[-1] for option in options[1:]]])
    return result


def _negative_identity_competitor(evidence):
    feature = evidence.get("identity_feature") if isinstance(evidence, dict) else None
    if not isinstance(feature, dict):
        decision = evidence.get("decision", {}) if isinstance(evidence, dict) else {}
        feature = decision.get("features") if isinstance(decision, dict) else None
    if not isinstance(feature, dict):
        return None
    competitor = feature.get("competitor_evidence")
    try:
        margin = float(feature.get("competitor_margin"))
        coverage = float(competitor.get("core_coverage"))
        evalue = float(competitor.get("evalue"))
    except (AttributeError, TypeError, ValueError):
        return None
    prototype = feature.get("competitor_prototype")
    if (margin > 0 or not isinstance(competitor, dict) or str(competitor.get("label", "")).lower() != "negative"
            or not prototype or competitor.get("prototype_id") != prototype or coverage < 0.8 or evalue > 1e-10):
        return None
    return {"prototype_id": prototype, "competitor_margin": margin,
            "core_coverage": coverage, "evalue": evalue,
            "comparison_scope": competitor.get("score_scope", "hash_bound_identity_feature")}


def _result(function, family, kind, confidence, summary, support=(), members=(), alternatives=(), risks=(), row=None):
    if not members and row and row.get("candidate_id") and (kind in {"core_function","family_function"} or function=="candidate_H2_H4MPT_hydrogenation_potential"):
        members=[{"role":"candidate_anchor","protein_id":row["candidate_id"]}]
    return {
        "preferred_function": function, "preferred_family": family, "prediction_kind": kind,
        "prediction_confidence": confidence, "prediction_summary": summary,
        "prediction_support": list(support), "preferred_members": list(members),
        "prediction_alternatives": list(alternatives), "prediction_risks": list(risks),
        "experimental_unknowns": _unknowns(row or {}), "prediction_confidence_scope": _SCOPE,
    }


_SPLIT_IDENTITY = "sequence_supported_uncalibrated_split_center"
_SPLIT_ATTRIBUTE = "NiFe_split_catalytic_core_potential"
_SPLIT_FUNCTION = "NiFe_split_catalytic_center_supported_function_unresolved"


def _split_center_proof(evidence):
    """Return the producer-owned split member proof, if it is well formed."""
    for item in (evidence or {}).get("attributes", ()):
        if not isinstance(item, dict) or item.get("attribute") != _SPLIT_ATTRIBUTE:
            continue
        if item.get("status") != "sequence_supported_candidate":
            continue
        value = _json(item.get("value"))
        if not isinstance(value, dict) or not isinstance(value.get("members"), list):
            continue
        members = value["members"]
        if len(members) != 2:
            continue
        ids = {member.get("protein_id") for member in members
               if isinstance(member, dict) and member.get("protein_id")}
        if len(ids) != 2:
            continue
        return value, members
    return None, []


def _predict(row, *, evidence=None):
    """Return a preferred qualitative prediction while preserving legacy input."""
    row = row if isinstance(row, dict) else {}
    evidence = evidence if isinstance(evidence, dict) else {}
    family = _family(row)
    if family in {"unresolved", "Fe"}: family = _screen_family(row, evidence) if _screen_family(row, evidence) != "unresolved" else family
    claims = _claims(evidence)
    legacy_label = _text(row.get("coarse_function_prediction") or row.get("function_prediction"))
    split_value, split_members = _split_center_proof(evidence)
    split_label = legacy_label == _SPLIT_FUNCTION
    split_identity = row.get("identity_status") == _SPLIT_IDENTITY
    if split_identity and split_members and (split_label or not legacy_label):
        risks = ["split_catalytic_center_only", "full_complex_unresolved",
                 "electron_partner_unresolved", "reaction_mechanism_unresolved"]
        limitations = split_value.get("limitations") or []
        risks.extend(limitations if isinstance(limitations, list) else [limitations])
        support = [{"attribute": _SPLIT_ATTRIBUTE, "status": "sequence_supported_candidate",
                    "member_count": len(split_members), "scope": "split_catalytic_center_only"}]
        return _result(_SPLIT_FUNCTION, "NiFe", "coarse_function", "low",
                       "支持NiFe催化中心由两个多肽成员共同构成；完整复合体、Hdr/其他电子伙伴、反应机制和实际活性仍未确定。",
                       support=support, members=split_members, risks=risks, row=row)
    # A split label without its producer proof is not a supported claim.  Do
    # not let it shadow an independent legacy function; an empty/synthetic
    # split label falls through to the ordinary conservative prediction.
    if split_label:
        legacy_label = ""
    if legacy_label == "candidate_FeFe_rubrerythrin_H2_H2O2_coupling" and family == "FeFe":
        fusion = next((c for c in claims if c.get('claim_name') == 'fefe_rubrerythrin'
                       and c.get('status') == 'assay_supported_candidate' and c.get('exact_reference') is True), None)
        if fusion:
            members = [dict(role=role, protein_id=row['candidate_id'], query_interval=fusion['domain_spans'][key])
                       for role, key in [('FeFe_H_cluster_domain','fefe'),('rubrerythrin_diiron_domain','rbr')]]
            return _result(legacy_label, family, 'reference_exact_potential', 'high',
                '与已测参考序列相同，支持体外 H₂ 氧化耦联 H₂O₂ 还原的反应潜力。查询活性、天然成熟与生理净方向未测。',
                [dict(claim='single_chain_assayed_reaction', assay=fusion['reference_assay'],
                      evidence_source='fefe_rubrerythrin_function_evidence.json')], members,
                risks=['reference_assay_scope_only','native_maturation_unknown'], row=row)
    claims.sort(key=lambda c: (_claim_supports("candidate_H2_pyridine_nucleotide_coupling",c),
                              not c.get("conditional_on_chosen_members",False),
                              not c.get("is_alternate_hypothesis",False)),reverse=True)
    # A claim is an architecture: core and module must come from this same
    # object.  Never combine positive fields from separate claim branches.
    for claim in claims:
        if (_positive_legacy(legacy_label) and legacy_label != "candidate_H2_pyridine_nucleotide_coupling"):
            continue
        core = _status(claim, "hydrogen_reaction_core_status", "core_status", "Hydrogen_reaction_core")
        module = _status(claim, "pyridine_nucleotide_module_status", "module_status", "Pyridine_nucleotide_electron_module")
        blockers = _blockers(claim)
        hard, soft = _hard(blockers), _soft(blockers)
        association = claim.get("association_supported")
        if association is None:
            association = claim.get("coupling_linkage_supported") is True
        if claim.get("module_association_supported") is False: association = False
        cross_anchor_only = bool(claim.get("cross_anchor_conflicts")) and not _same_resource_overlap(claim)
        if (claim.get("resource_conflict") is True and not cross_anchor_only) or _same_resource_overlap(claim):
            hard.append("resource_conflict")
        if _positive(core) and _positive(module) and association is True and not hard:
            alternatives = _alternatives(claim, evidence)
            members=_members(claim)
            for rival in claims:
                if rival is claim or not _claim_supports("candidate_H2_pyridine_nucleotide_coupling",rival):continue
                if rival.get("structure_id")==claim.get("structure_id") and _members(rival)==members:continue
                alternatives.append({**_claim_record(rival),"function_prediction":"candidate_H2_pyridine_nucleotide_coupling"})
            risks = list(dict.fromkeys(soft + claim.get("selection_risks", []) + [str(x) for x in claim.get("unresolved_attributes", [])]))
            if claim.get("cross_anchor_conflicts"):risks.append("cross_anchor_partner_ownership_ambiguous")
            confidence = "medium"
            ref = _text(row.get("reference_reaction_status"))
            caps = _json(row.get("reference_assay_capabilities"))
            if ref.startswith("exact_reference_") and not risks and _capability_matches("candidate_H2_pyridine_nucleotide_coupling", caps):
                confidence = "high"
            elif any("cross" in x.lower() or "alternative_function" in x.lower() for x in risks):
                confidence = "low"
            if claim.get('confidence_hint')=='low':confidence='low'
            context_limited = any("reference_context" in str(x).lower() or "context_unresolved" in str(x).lower() for x in blockers)
            if context_limited:
                confidence = "low"
                risks.extend(["reference_architecture_context_unresolved", "physical_assembly_unknown"])
                for role, context_members in claim.get("context_roles", {}).items():
                    if not context_members:risks.append("context_role_missing:"+role)
                    elif any(m.get("partial") is True for m in context_members):risks.append("context_role_partial:"+role)
            profile_note='低同一性中继由长覆盖同源、Fe–S域和参考配位模式共同支持；' if claim.get('entry_route')=='profile_supported_low_identity_role_complement' else ''
            compensation=claim.get('composition_compensation',{})
            support=[f"Hydrogen_reaction_core:{core}", f"Pyridine_nucleotide_electron_module:{module}", "association_supported", "architecture:"+(claim.get("structure_id") or claim.get("claim_name", "unknown")), "hypothesis:"+claim.get("claim_name", "evidence")]
            if claim.get('entry_route')=='independent_Group4_core_pyridine_module_composition':
                profile_note+='Group4核心与局部核苷酸模块共同支持，核心连接架构未定；'
                support.append(dict(claim='independent_Group4_core_module_composition',strict_core_status=claim.get('strict_core_status'),
                                    evidence_source='pyridine_role_composition_evidence.json',scope='coarse_contextual_coupling_not_canonical_SHI_architecture'))
                alternatives.append(dict(preferred_function='hydrogenase_core_candidate',
                    reason='核苷酸模块也可能通向其他氧化还原支路，核心到该模块的物理连接尚未验证。'))
            if compensation.get('applied') is True:
                profile_note+='邻近亚基组合证据补偿了弱同源或角色区间限制；'
                support.append(dict(claim='cooperative_role_compensation',roles=compensation.get('roles',[]),
                                    strict_core_status=compensation.get('strict_core_status'),
                                    evidence_source='pyridine_role_composition_evidence.json',
                                    scope='joint_local_roles_not_strict_reference_upgrade'))
            return _result("candidate_H2_pyridine_nucleotide_coupling", family, "coarse_function", confidence,
                           f"预测为{family if family != 'unresolved' else ''}氢酶的核苷酸耦联功能，{dict(high='高',medium='中等',low='低')[confidence]}置信度；{'完整参考架构尚未闭合；' if context_limited else ''}{profile_note}NAD/NADP选择性、分叉机制和实验活性仍未确定。",
                           support,
                           members, alternatives, risks, row)
    for claim in claims:
        if _positive_legacy(legacy_label) and legacy_label != "candidate_H2_pyridine_nucleotide_coupling":
            continue
        if (claim.get("nonlocal_sequence_potential") is True and
                claim.get("nonlocal_sequence_potential_status") == "sequence_supported_nonlocal_potential" and
                claim.get("nonlocal_sequence_function_status") == "candidate_H2_pyridine_nucleotide_coupling"):
            members = _members(claim)
            risks = list(dict.fromkeys(["nonlocal_member_association", "physical_assembly_unknown",
                                        "local_coupling_linkage_not_supported", "NAD_NADP_selectivity_unresolved",
                                        *claim.get("nonlocal_sequence_blockers", [])]))
            return _result("candidate_H2_pyridine_nucleotide_coupling", family, "coarse_function", "low",
                           "预测为核苷酸耦联氢酶序列功能潜力，低置信度；四角色同源和模块证据完整，但成员为非局部组合，物理组装及实际耦联未确定。",
                           ["complete_NiFe_HoxFUYH_role_homology", "pyridine_module_sequence_support", "distinct_selected_members", "nonlocal_sequence_hypothesis"],
                           members, risks=risks, row=row)
    # HmdII is a real competing interpretation, not canonical H2 function.
    if _positive_hmdii(claims, row):
        label = _text(row.get("coarse_function_prediction") or row.get("function_prediction"))
        function = label if label == "candidate_H2_THF_derivative_hydrogenation_potential" else "hmdII_hydrogenase_like_candidate"
        return _result(function, "HmdII", "family_function", "medium" if label else "low",
                       ("预测为四氢叶酸衍生物加氢潜力，中等置信度；HmdII特异性和实验活性仍未确定。" if label else "预测为HmdII样氢酶候选，不判定为典型H2氢酶功能。"),
                       ["HmdII-like evidence"], risks=["canonical_H2_function_not_supported"], row=row)
    label = _text(row.get("coarse_function_prediction") or row.get("function_prediction"))
    # The local HycE/G and FeS relay + CODH route is deliberately a low-confidence
    # preferred explanation.  It must not upgrade the legacy current
    # annotation to a reaction contract or infer CO direction/activity.
    co_context = next((c for c in claims
                       if c.get("status") == "candidate_CO_linked_H2_pathway"
                       and c.get("entry_route") == "hyc_core_local_CODH_context"), None)
    if co_context and legacy_label in {"hydrogenase_family_supported_function_unresolved",
                                      "H2_core_supported_partner_unresolved"}:
        members = _members(co_context)
        roles = {m.get("role") for m in members}
        relay_role = "HycF" if "HycF" in roles else "H2_FeS_relay"
        required = {"HycE", "HycG", relay_role, "CooS", "CODH_FeS_relay"}
        if roles == required and len(members) == len({m['protein_id'] for m in members}) == 5:
            support = [
                {"claim": "hyc_core_local_CODH_context", "status": co_context.get("status"),
                 "entry_route": co_context.get("entry_route"),
                 "evidence_source": "coo_function_evidence.json"},
                {"claim": "core_roles", "roles": ["HycE", "HycG", relay_role],
                 "evidence_source": ("coo_function_evidence.json" if relay_role == "H2_FeS_relay"
                                     else "fhl_function_evidence.json")},
                {"claim": "CODH_roles", "roles": ["CooS", "CODH_FeS_relay"],
                 "evidence_source": "coo_function_evidence.json"},
            ]
            risks = ["CO_sequence_support_limited", "local_context_association_conditional",
                     "CO_reaction_direction_unresolved", "query_activity_not_tested",
                     "mechanism_unresolved"]
            core = "局部HycE/G核心及氢酶侧Fe-S中继" if relay_role == "H2_FeS_relay" else "局部HycE/G/F核心"
            return _result("candidate_CO_linked_H2_pathway", family, "coarse_function", "low",
                           f"{core}与CODH相关CooS/Fe-S中继共同支持条件性CO关联氢代谢路径；序列支持和局部上下文有限，CO反应方向、查询活性及机制未确定。",
                           support, members,
                           alternatives=[{"preferred_function":"hydrogenase_core_candidate",
                                          "reason":"核心身份成立，但CO关联仍依赖条件性局部上下文。"}],
                           risks=risks, row=row)

    # Existing specialized claims are more informative than generic core rows.
    # This preserves Frh/Hdr/MBH/etc. instead of collapsing them to H2 core.
    if _positive_legacy(label):
        matching = _matching_claim(label, claims)
        capabilities = _json(row.get("reference_assay_capabilities"))
        confidence = "high" if _text(row.get("reference_reaction_status")).startswith("exact_reference_") and _capability_matches(label, capabilities) else "medium"
        if label in {"candidate_genome_H2_formate_pathway", "candidate_H2_formate_pathway_conditional_on_recoding", "candidate_CO_linked_H2_pathway"}:confidence = "low"
        risks=["partner_or_mechanism_detail_unresolved"]
        membrane_context = matching and matching.get("claim_name") == "membrane_context"
        if membrane_context:
            confidence = "low"
            risks.extend(matching.get("selection_risks", []))
        if label=="candidate_H2_Hdr_coupling" and matching and matching.get("HdrA_pairing_status")!="unique_sequence_candidate":
            risks.extend(["HdrA_partner_ambiguous","HdrA_capacity_witness_not_physical_pairing"])
        if label=="candidate_H2_Hdr_coupling" and matching:
            assignment=matching.get("global_function_assignment",{})
            if assignment.get("status") in {"unresolved_conflict_or_evidence","no_feasible_distinct_HdrA_assignment"}:
                risks.append("Hdr_global_assignment_"+assignment["status"])
            risks.extend(b for b in assignment.get("blockers",[]) if ":" not in b or b.partition(":")[2]==row.get("candidate_id"))
        alternatives=_alternatives(matching,evidence) if matching else []
        for claim in claims:
            if claim is not matching and _claim_supports("candidate_H2_pyridine_nucleotide_coupling",claim):
                alternatives.append({**_claim_record(claim),"function_prediction":"candidate_H2_pyridine_nucleotide_coupling"})
        return _result(label, family, "coarse_function", confidence,
                       ("预测为膜关联氢酶氧化还原系统，低置信度；局部催化核心、Fe–S与膜模块共同支持；具体电子伙伴、离子耦联和实验活性仍未确定。"
                        if membrane_context else _legacy_summary(label, {"high":"高","medium":"中等","low":"低"}[confidence])),
                       ["legacy_positive_coarse_claim"] + (["typed_claim:" + matching.get("claim_name", "unknown")] if matching else []),
                       _members(matching) if matching else (),
                       alternatives,
                       risks, row)
    # Core-only evidence remains useful even without a partner/module.
    for claim in claims:
        core = _status(claim, "hydrogen_reaction_core_status", "core_status", "Hydrogen_reaction_core")
        blockers = _blockers(claim)
        if _positive(core) and not any("not_supported" in x.lower() for x in blockers):
            return _result("hydrogenase_core_candidate", family, "core_function", "medium",
                           "预测为氢酶催化核心候选，中等置信度；伙伴或电子模块尚未充分解析。",
                           [f"Hydrogen_reaction_core:{core}"], _members(claim), _alternatives(claim, evidence), blockers, row)
    # A trusted catalytic identity with no specialized claim is still a core
    # prediction; this is distinct from a weak/shared redox-domain signal.
    attrs = _flatten_strings(row.get("supported_attributes"))
    noncatalytic_family = any(x in family.lower() for x in ("nuo", "narf", "redox_module"))
    if not noncatalytic_family and (_text(row.get("identity_status")) in {"calibrated_high", "sequence_supported_uncalibrated_core", "sequence_supported_uncalibrated_noncanonical"} or "Hydrogen_reaction_core" in attrs):
        return _result("hydrogenase_core_candidate", family, "core_function", "medium",
                       "预测为氢酶催化核心候选，中等置信度；伙伴或电子模块尚未充分解析。",
                       ["trusted_catalytic_identity"], risks=["partner_or_electron_module_unresolved"], row=row)
    shared = any(x.lower() in {"electron_transfer_module_presence", "nuo/fmn_redox_domain", "shared_redox_domain"} for x in attrs)
    negative = _negative_identity_competitor(evidence)
    if negative:
        return _result("shared_redox_module_unresolved", "unresolved", "electron_module_only", "medium",
                       "更符合已审定的非催化同源蛋白，具体功能未定；该比较证据不等同于查询样本的实测H2阴性。",
                       ["negative_competitor_prototype:" + negative["prototype_id"], negative],
                       risks=["catalytic_identity_competitor_preferred"], row=row)
    if (shared or "PF00346" in _text(evidence.get("screen", {}))) and _screen_family(row, evidence) == "unresolved":
        return _result("shared_redox_module_unresolved", "unresolved", "electron_module_only", "low",
                       "检测到共享氧化还原模块，但证据不足以判定为氢酶。", ["shared_redox_module_only"], risks=["hydrogenase_identity_not_supported"], row=row)
    catalytic_family = any(x in family.lower() for x in ("nife", "fefe", "hmd", "hydrogenase"))
    catalytic_signal = any(x.lower() in {"hydrogen_reaction_core", "h_cluster_redox_core_potential", "catalytic_core_region"} for x in attrs)
    screen_signal = _screen_family(row, evidence) != "unresolved"
    if catalytic_family and (catalytic_signal or screen_signal or any(x in _text(row.get("identity_status")).lower() for x in ("fragment", "distant", "weak", "potential"))):
        return _result("possible_hydrogenase", family, "possible_hydrogenase", "low",
                       f"可能为{family if family != 'unresolved' else ''}氢酶，低置信度；催化身份或关键证据不足。",
                       ["weak_or_distant_catalytic_signal"], risks=["identity_false_positive_or_fragment"], row=row)
    return _result("unresolved_candidate", family, "unresolved", "low", "当前证据不足以形成可靠的氢酶功能首选判断。", risks=["catalytic_identity_not_supported"], row=row)



def _context_function_prediction(row, claims, prediction):
    claim = next((c for c in claims if c.get('claim_name') == 'hhy_context'), {})
    if claim.get('status') != 'candidate_H2_oxidation_potential' or claim.get('blockers'):
        return None
    if row.get('identity_status') not in {'calibrated_high', 'sequence_supported_uncalibrated_core'}:
        return None
    if prediction['preferred_function'] not in {'hydrogenase_core_candidate', 'possible_hydrogenase', 'unresolved_candidate'}:
        return None
    members = _members(claim)
    if {m['role'] for m in members} != {'HhyL', 'HhyS'} or len(members) != 2:
        return None
    if next(m['protein_id'] for m in members if m['role'] == 'HhyL') != row.get('candidate_id'):
        return None
    confidence = 'low' if claim.get('confidence') == 'low' else 'medium'
    risk_text = '电子亚基参考Cys模式存在差异。' if 'HhyS_reference_ligand_pattern_conflict' in claim.get('risks', []) else ''
    summary = (f"预测为Hhy样NiFe氢酶的H2氧化/摄氢潜力，{'低' if confidence == 'low' else '中等'}置信度；"
               f"完整HhyL/HhyS、局部关联和参考全细胞遗传证据支持。{risk_text}直接电子受体、亲和力、查询活性及体内净方向未定。")
    support = {'claim': 'Hhy_context_function', 'status': claim['status'], 'structure_id': claim['structure_id'],
               'capability_scope': claim.get('capability_scope'), 'reference_evidence_scope': claim.get('reference_evidence_scope'),
               'primary_source': claim.get('primary_source'), 'scope': 'coarse_sequence_potential_not_query_activity'}
    alternatives = list(prediction.get('prediction_alternatives', []))
    if confidence == 'low':
        alternatives.append({'preferred_function': 'hydrogenase_core_candidate',
                             'reason': 'Electron relay function may be impaired; source-reference details remain conditional.'})
    risks = [r for r in prediction.get('prediction_risks', []) if r != 'partner_or_electron_module_unresolved']
    return _result('candidate_H2_oxidation_potential', 'NiFe', 'coarse_function', confidence, summary,
                   [*prediction.get('prediction_support', []), support], members, alternatives,
                   [*risks, *claim.get('risks', []), 'reference_whole_cell_context_only', 'query_activity_not_tested'], row)


def predict(row, *, evidence=None):
    """Attach auditable support and claim-specific risks to the best explanation."""
    evidence=evidence or {}
    prediction=_predict(row,evidence=evidence)
    claims=_claims(evidence)
    prediction=_reference_override(row, claims, prediction) or _component_prediction(row, claims, prediction) or prediction
    prediction=_single_chain_reference_prediction(row,claims,prediction,evidence) or prediction
    prediction=_source_like_prediction(row,claims,prediction) or prediction
    prediction=_context_function_prediction(row,claims,prediction) or prediction
    excluded_partners={}
    for claim in claims:
        if claim.get('claim_name') not in {'fhl','hdcr','tonn'} or claim.get('fdh_module_status')!='unresolved':
            continue
        roles=claim.get('roles',{})
        members=[*claim.get('excluded_FdhF_candidates',[]),*(m for role in ('FdhF','FdhF2','Fdh2_A') for m in roles.get(role,[]))]
        for member in members:
            proof=member.get('independent_role_competition',{})
            if proof.get('status')!='non_fdh_preferred' or proof.get('query_id')!=member.get('protein_id'):
                continue
            winner=next((x['strongest_counter'] for x in proof.get('evidence',[]) if x.get('strongest_counter')),None)
            if winner:
                excluded_partners[member['protein_id']]={
                    'claim':'FDH_partner_role','code':'preferred_non_FDH_role_competition',
                    'protein_id':member['protein_id'],'preferred_role_class':winner['role_class'],
                    'reference_id':winner['reference_id'],'paired_identity':winner['paired_identity'],
                    'scope':'query_partner_interpretation_only','source':claim['claim_name']+'_function_evidence.json'}
    for risk in excluded_partners.values():
        prediction['prediction_risks'].append(risk)
        prediction['prediction_alternatives'].append(dict(risk,kind='preferred_other_partner_role'))
    if excluded_partners and prediction['preferred_function']=='hydrogenase_core_candidate':
        prediction['prediction_summary']+=' FDH样候选有更强的其他酶功能解释；甲酸耦联伙伴仍未确定。'
    supported=_json(row.get("supported_attributes", []))
    for attribute in evidence.get("attributes", []):
        if attribute.get("attribute") in supported:
            prediction["prediction_support"].append({k:attribute[k] for k in ("attribute","status","evidence_source") if k in attribute})
    for reason in evidence.get("progression", {}).get("all_blocking_reasons", []):
        prediction["prediction_risks"].append({k:reason[k] for k in ("claim","code","type","source") if k in reason})
    if prediction["prediction_kind"] in {"possible_hydrogenase","unresolved"} and evidence.get("gene", {}).get("partial") in {True,"True"}:
        prediction["prediction_risks"].append("partial_gene_call")
    return prediction
