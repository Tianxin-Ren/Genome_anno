"""Unified mechanistic constraint engine for hydrogenase function inference.

Replaces organism-specific function templates with generalizable biochemical
constraints: confirmed catalytic core + local operon/fusion electron-transfer
role primitives + configurable residue-level gates.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RULE_VERSION = "mechanistic_rules_v1"


def _yes(val):
    return val is True or str(val).strip().lower() in {"true", "1", "yes"}


def load_rules(root=ROOT):
    path = Path(root) / "config/mechanistic_rules.json"
    return json.loads(path.read_text())


def _base_acc(acc):
    return str(acc or "").split(".")[0]


def _build_role_lookup(rules):
    lookup = {}
    for role, specs in rules["domain_roles"].items():
        if isinstance(specs, dict):
            for subrole, accs in specs.items():
                for a in accs:
                    lookup.setdefault(a, []).append((role, subrole))
        else:
            for a in specs:
                lookup.setdefault(a, []).append((role, None))
    return lookup


def _gene_distance(g1, g2, allow_circular=False):
    """Return non-negative intergenic bp distance (0 if overlapping) or None."""
    if not g1 or not g2 or g1.get("contig") != g2.get("contig"):
        return None
    s1, e1 = int(g1["start"]), int(g1["end"])
    s2, e2 = int(g2["start"]), int(g2["end"])
    if max(s1, s2) <= min(e1, e2):
        return 0
    linear = max(s1, s2) - min(e1, e2) - 1
    if allow_circular and g1.get("contig_length"):
        clen = int(g1["contig_length"])
        if clen > 0:
            wrap = min((s2 - e1 - 1) % clen, (s1 - e2 - 1) % clen)
            return min(linear, wrap)
    return linear


def _local_cluster_genes(anchor_id, genes_by_id, linkage, input_mode):
    """Collect anchor + same-contig/strand neighbors within max_distance_bp."""
    anchor_gene = genes_by_id.get(anchor_id)
    if not anchor_gene:
        return [anchor_id] if anchor_id else []
    if input_mode not in set(linkage.get("eligible_input_modes", ["genome", "assembly", "bin"])):
        return [anchor_id]
    max_dist = int(linkage.get("max_distance_bp", 5000))
    same_strand = bool(linkage.get("require_same_strand", True))
    allow_circ = bool(linkage.get("allow_circular", True)) and input_mode == "genome"
    cluster = []
    for gid, g in genes_by_id.items():
        if gid == anchor_id:
            cluster.append((0, gid))
            continue
        if same_strand and str(g.get("strand")) != str(anchor_gene.get("strand")):
            continue
        dist = _gene_distance(anchor_gene, g, allow_circular=allow_circ)
        if dist is not None and dist <= max_dist:
            cluster.append((dist, gid))
    cluster.sort(key=lambda x: (x[0], x[1]))
    return [gid for _, gid in cluster]


def _role_accessions(found_roles, role_name):
    return {_base_acc(ev.get("accession")) for ev in found_roles.get(role_name, [])}


def _infer_subtype(primary_label, family, found_roles):
    """Layer-2 structural/organizational subtype within an electron-transfer class."""
    if primary_label == "candidate_H2_quinone_coupling":
        quin_accs = _role_accessions(found_roles, "MEMBRANE_QUINONE_CYT")
        mat_accs = _role_accessions(found_roles, "MATURATION")
        if "PF07849" in quin_accs:
            return "NiFe_HucM_scaffold"
        if mat_accs & {"PF07449", "PF04809"}:
            return "NiFe_O2_tolerant_cytb"
        if "PF03916" in quin_accs:
            return "NiFe_HybA_polyFeS_quinone"
        if "PF00173" in quin_accs and "PF01292" not in quin_accs:
            return "NiFe_cytb5_quinone"
        return "NiFe_standard_cytb"
    if primary_label == "candidate_H2_cytochrome_c3_coupling":
        return "NiFe_periplasmic_cytochrome_c3"
    if primary_label == "candidate_membrane_associated_H2_redox_system":
        pump_accs = _role_accessions(found_roles, "MEMBRANE_PROTON_PUMP")
        if pump_accs & {"PF04039", "PF00420", "PF00662", "PF01899", "PF03334", "PF04066"}:
            return "NiFe_Group4_Mbh_Mrp"
        return "NiFe_Group4_Ech"
    if primary_label in {"candidate_genome_H2_formate_pathway", "candidate_H2_formate_coupling"}:
        if family == "FeFe":
            return "FeFe_HDCR_Fdh"
        if "MEMBRANE_PROTON_PUMP" in found_roles:
            return "NiFe_Group4_FHL_Fdh_pump"
        return "NiFe_soluble_Fdh_coupled"
    if primary_label == "candidate_CO_linked_H2_pathway":
        return "NiFe_Group4_CODH_Coo"
    if primary_label == "candidate_H2_Hdr_coupling":
        return "NiFe_Mvh_Hdr"
    if primary_label == "candidate_H2_F420_reduction":
        return "NiFe_Frh_F420"
    if primary_label == "candidate_H2_pyridine_nucleotide_coupling":
        pyr_accs = _role_accessions(found_roles, "PYRIDINE_FLAVIN_DIAPHORASE")
        fes_accs = _role_accessions(found_roles, "FES_RELAY")
        if family == "FeFe":
            if (pyr_accs & {"PF01512", "PF10588"}) or (fes_accs & {"PF10589", "PF13510"}):
                return "FeFe_NuoEF_diaphorase"
            if pyr_accs & {"PF07992", "PF10418", "PF00175", "PF00970"}:
                return "FeFe_FAD_diaphorase"
            return "FeFe_NuoEF_diaphorase"
        if pyr_accs & {"PF01512", "PF10588"}:
            return "NiFe_FMN_NAD_HoxFUYH"
        if pyr_accs & {"PF10418", "PF00175", "PF00970", "PF07992"}:
            return "NiFe_FAD_NADP_SH"
        return "NiFe_pyridine_other"
    if primary_label == "candidate_H2_ferredoxin_compatible_core":
        if "FES_RELAY" in found_roles:
            return "FeFe_F_domain_ferredoxin"
        return "FeFe_monomeric_ferredoxin"
    if primary_label == "H2_regulatory_signal":
        return "NiFe_H2_sensor_kinase"
    if primary_label == "candidate_H2_H4MPT_hydrogenation_potential":
        return "Fe_only_canonical_Hmd"
    if primary_label == "candidate_H2_THF_derivative_hydrogenation_potential":
        return "Fe_only_HmdII_paralog"
    return primary_label


def _index_domains(domains):
    by_protein = {}
    for d in domains:
        pid = d.get("protein_id")
        if not pid:
            continue
        if "passes" in d and not _yes(d["passes"]):
            continue
        by_protein.setdefault(pid, []).append(d)
    return by_protein


def infer_candidate_function(
    qid,
    family,
    identity_status,
    input_mode,
    genes_by_id,
    domains_by_protein,
    rules,
    role_lookup,
    site_gates=None,
    genome_acceptors=None,
):
    """Infer function and mechanistic evidence for one candidate protein."""
    gates = rules["gates"]
    confirmed_set = set(gates["confirmed_identity_statuses"])
    core_confirmed = identity_status in confirmed_set

    cluster_ids = _local_cluster_genes(qid, genes_by_id, rules["locus_linkage"], input_mode)
    found_roles = {}  # role -> list of evidence dicts
    partner_subtypes = {}  # subrole -> list of evidence dicts

    for gid in cluster_ids:
        is_fusion = gid == qid
        for d in domains_by_protein.get(gid, []):
            acc = _base_acc(d.get("accession"))
            for role, subrole in role_lookup.get(acc, []):
                if role == "CATALYTIC_CORE":
                    continue
                ev = {
                    "protein_id": gid,
                    "accession": d.get("accession"),
                    "profile": d.get("profile", ""),
                    "fusion": is_fusion,
                }
                found_roles.setdefault(role, []).append(ev)
                if subrole:
                    partner_subtypes.setdefault(subrole, []).append(ev)

    if gates.get("require_catalytic_core_confirmation", True) and not core_confirmed:
        return {
            "function_prediction": "hydrogenase_candidate_identity_unresolved",
            "coarse_function_prediction": "hydrogenase_candidate_identity_unresolved",
            "confidence": "unresolved",
            "supported_roles": sorted(found_roles),
            "partner_subtypes": sorted(partner_subtypes),
            "cluster_members": cluster_ids,
            "role_evidence": found_roles,
            "gate_notes": ["catalytic_core_identity_unconfirmed"],
            "secondary_modules": [],
            "NAD_NADP_selectivity": rules["selectivity_and_bifurcation"]["nad_nadp_selectivity_default"],
            "bifurcation": rules["selectivity_and_bifurcation"]["bifurcation_default"],
        }

    # Check residue-level gates
    gate_notes = []
    site_gates = site_gates or {}
    fdh_ok = True
    if "Fdh" in partner_subtypes and gates.get("fdhf_catalytic_site_gate", True):
        if site_gates.get("fdh_site_status") == "excluded":
            fdh_ok = False
            gate_notes.append("FdhF_catalytic_Cys_Sec_excluded")
        elif site_gates.get("fdh_site_status") != "supported":
            gate_notes.append("FdhF_catalytic_site_unverified")

    hmd_mode = None
    if family == "Fe" and gates.get("hmd_his14_gate", True):
        hmd_mode = site_gates.get("hmd_gate_status")
        if hmd_mode == "hmdii_paralog":
            gate_notes.append("HmdII_paralog_signature")
        elif hmd_mode == "canonical_hmd":
            gate_notes.append("canonical_Hmd_His14_supported")

    # Arbitrate according to configured priority_order
    matched = []
    has_pump = "MEMBRANE_PROTON_PUMP" in found_roles
    has_quin = "MEMBRANE_QUINONE_CYT" in found_roles
    has_cytc3 = "CYTOCHROME_C3" in found_roles
    has_pyr = "PYRIDINE_FLAVIN_DIAPHORASE" in found_roles
    has_fes = "PROXIMAL_FES" in found_roles or "FES_RELAY" in found_roles

    for rule_name in rules["arbitration"]["priority_order"]:
        if rule_name == "REGULATORY_SENSOR" and "REGULATORY_SENSOR" in found_roles:
            sensor_accs = _role_accessions(found_roles, "REGULATORY_SENSOR")
            if sensor_accs & {"PF00512", "PF02518", "PF07568", "PF07730"}:
                matched.append(("H2_regulatory_signal", "regulatory_sensor_domain_colocalized"))
        elif rule_name == "METABOLIC_COMPOSITE_FHL" and rules["arbitration"].get("allow_composite_fhl", True):
            if "Fdh" in partner_subtypes and fdh_ok and has_pump:
                matched.append(("candidate_genome_H2_formate_pathway", "Fdh_plus_membrane_proton_pump"))
        elif rule_name == "METABOLIC_PARTNER":
            if "Fdh" in partner_subtypes and fdh_ok:
                matched.append(("candidate_H2_formate_coupling", "Fdh_partner_colocalized"))
            if "CODH" in partner_subtypes:
                matched.append(("candidate_CO_linked_H2_pathway", "CODH_partner_colocalized"))
            if "Hdr" in partner_subtypes:
                matched.append(("candidate_H2_Hdr_coupling", "Hdr_partner_colocalized"))
            if "F420" in partner_subtypes:
                matched.append(("candidate_H2_F420_reduction", "F420_partner_colocalized"))
        elif rule_name == "MEMBRANE_PROTON_PUMP" and has_pump:
            matched.append(("candidate_membrane_associated_H2_redox_system", "membrane_proton_pump_colocalized"))
        elif rule_name == "MEMBRANE_QUINONE_CYT" and has_quin:
            matched.append(("candidate_H2_quinone_coupling", "membrane_quinone_or_cyt_colocalized"))
        elif rule_name == "CYTOCHROME_C3" and has_cytc3:
            matched.append(("candidate_H2_cytochrome_c3_coupling", "cytochrome_c3_colocalized"))
        elif rule_name == "PYRIDINE_FLAVIN_DIAPHORASE" and has_pyr:
            matched.append(("candidate_H2_pyridine_nucleotide_coupling", "pyridine_flavin_diaphorase_colocalized"))
        elif rule_name == "FE_ONLY_HMD" and family == "Fe":
            if hmd_mode == "hmdii_paralog":
                matched.append(("candidate_H2_THF_derivative_hydrogenation_potential", "Fe_hydrogenase_HmdII_branch"))
            else:
                matched.append(("candidate_H2_H4MPT_hydrogenation_potential", "Fe_hydrogenase_Hmd_branch"))

    if matched:
        primary_label, primary_reason = matched[0]
        secondary = [lbl for lbl, _ in matched[1:]] if rules["arbitration"].get("record_secondary_modules", True) else []
        confidence = "uncalibrated_function_potential"
    elif family == "FeFe" and rules["arbitration"].get("fefe_ferredoxin_default", True):
        primary_label = "candidate_H2_ferredoxin_compatible_core"
        primary_reason = "FeFe_ferredoxin_core_with_relay" if has_fes else "FeFe_monomeric_ferredoxin_core"
        secondary = []
        confidence = "uncalibrated_function_potential"
    elif (
        family == "NiFe"
        and has_fes
        and input_mode == "genome"
        and rules["locus_linkage"].get("allow_genome_unlinked_soluble_acceptors", False)
        and genome_acceptors
    ):
        if genome_acceptors.get("CYTOCHROME_C3"):
            primary_label = "candidate_H2_cytochrome_c3_coupling"
            primary_reason = "genome_level_unlinked_cytochrome_c3"
            gate_notes.append("genome_level_unlinked_cytochrome_c3")
            secondary = []
            confidence = "uncalibrated_function_potential"
        elif genome_acceptors.get("Hdr"):
            primary_label = "candidate_H2_Hdr_coupling"
            primary_reason = "genome_level_unlinked_Hdr"
            gate_notes.append("genome_level_unlinked_Hdr")
            secondary = []
            confidence = "uncalibrated_function_potential"
        else:
            primary_label = "H2_core_supported_partner_unresolved"
            primary_reason = "FeS_relay_present_without_terminal_acceptor_module"
            secondary = []
            confidence = "unresolved"
    elif has_fes:
        primary_label = "H2_core_supported_partner_unresolved"
        primary_reason = "FeS_relay_present_without_terminal_acceptor_module"
        secondary = []
        confidence = "unresolved"
    else:
        primary_label = "hydrogenase_family_supported_function_unresolved"
        primary_reason = "catalytic_core_without_colocalized_electron_modules"
        secondary = []
        confidence = "unresolved"

    return {
        "function_prediction": primary_label,
        "coarse_function_prediction": primary_label,
        "mechanistic_subtype": _infer_subtype(primary_label, family, found_roles),
        "confidence": confidence,
        "primary_reason": primary_reason,
        "supported_roles": sorted(found_roles),
        "partner_subtypes": sorted(partner_subtypes),
        "cluster_members": cluster_ids,
        "role_evidence": found_roles,
        "gate_notes": gate_notes,
        "secondary_modules": secondary,
        "NAD_NADP_selectivity": rules["selectivity_and_bifurcation"]["nad_nadp_selectivity_default"],
        "bifurcation": rules["selectivity_and_bifurcation"]["bifurcation_default"],
    }


def evaluate_all(anchors, progress_rows, genes_by_id, domains, root=ROOT, site_gates_by_id=None):
    """Run unified mechanistic inference across all candidate rows."""
    rules = load_rules(root)
    role_lookup = _build_role_lookup(rules)
    domains_by_protein = _index_domains(domains)
    anchor_mode = {a["protein_id"]: a.get("input_mode", "genome") for a in anchors}
    site_gates_by_id = site_gates_by_id or {}

    genome_acceptors = {}
    for pid, dlist in domains_by_protein.items():
        for d in dlist:
            for role, subrole in role_lookup.get(_base_acc(d.get("accession")), []):
                if role == "CYTOCHROME_C3":
                    genome_acceptors.setdefault("CYTOCHROME_C3", set()).add(pid)
                elif role == "METABOLIC_PARTNER" and subrole == "Hdr":
                    genome_acceptors.setdefault("Hdr", set()).add(pid)

    results = {}
    for row in progress_rows:
        qid = row["candidate_id"]
        results[qid] = infer_candidate_function(
            qid=qid,
            family=row["family"],
            identity_status=row["identity_status"],
            input_mode=anchor_mode.get(qid, "genome"),
            genes_by_id=genes_by_id,
            domains_by_protein=domains_by_protein,
            rules=rules,
            role_lookup=role_lookup,
            site_gates=site_gates_by_id.get(qid),
            genome_acceptors=genome_acceptors,
        )
    return results
