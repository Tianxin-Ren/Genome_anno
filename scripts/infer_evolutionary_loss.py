#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import math
import pandas as pd

def safe_int(val, default=0):
    if pd.isna(val) or val is None or val == "":
        return default
    try:
        return int(float(val))
    except Exception:
        return default

def safe_float(val, default=0.0):
    if pd.isna(val) or val is None or val == "":
        return default
    try:
        return float(val)
    except Exception:
        return default

def normalize_tristate(val):
    if pd.isna(val) or val is None:
        return "NA"
    s = str(val).strip().upper()
    if s in ["1", "1.0", "TRUE", "YES"]:
        return "1"
    if s in ["0", "0.0", "FALSE", "NO"]:
        return "0"
    return "NA"

def safe_str(val, default=""):
    if pd.isna(val) or val is None:
        return default
    return str(val)

SCHEMA_COLUMNS = [
    "target_id",
    "module_id",
    "module_name",
    "GTDB_taxonomy",
    "completeness",
    "contamination",
    "nearest_positive_ref",
    "positive_ref_count",
    "positive_taxonomic_group_count",
    "positive_lineage_count",
    "negative_ref_count",
    "ANI",
    "AF_query",
    "AF_ref",
    "AF_min",
    "left_anchor_1",
    "left_anchor_2",
    "right_anchor_1",
    "right_anchor_2",
    "left_anchors_matched_count",
    "right_anchors_matched_count",
    "left_matched_genes",
    "right_matched_genes",
    "left_anchor_region_other_genes",
    "right_anchor_region_other_genes",
    "intervening_genes_count",
    "intervening_genes",
    "target_locus_total_genes",
    "target_locus_other_genes",
    "anchor_left_identity",
    "anchor_right_identity",
    "anchor_order_consistent",
    "anchor_orientation_consistent",
    "same_contig",
    "target_contig_id",
    "distance_to_left_edge",
    "distance_to_right_edge",
    "required_edge_margin",
    "expected_locus_length",
    "observed_anchor_gap",
    "target_contig_length",
    "expected_gene_count",
    "protein_detected_count",
    "location_unknown_count",
    "blastx_rescued_count",
    "frameshift_count",
    "premature_stop_count",
    "mutation_evidence_status",
    "annotation_issue",
    "pseudogene_detected",
    "member_gene_states",
    "cluster_state",
    "structural_event_type",
    "nucleotide_deletion_support",
    "junction_support_status",
    "junction_spanning_reads",
    "junction_mean_mapq",
    "junction_coverage_ratio",
    "parent_ancestral_state",
    "parsimony_loss_events",
    "locus_state",
    "evolutionary_state",
    "confidence",
    # OFI Quantitative Metrics & Classifications
    "locator_score",
    "bilateral_context_score",
    "T_local",
    "T_global",
    "R_known",
    "I_novel",
    "A_local",
    "A_global",
    "B_bridge",
    "Q_continuity",
    "D_closure",
    "K_loss",
    "K_clean",
    "loss_gate_status",
    "OFI_fate",
    "OFI_occupancy",
    "fate_category",
    "operon_structure_state",
    "operon_prediction",
    "same_strand",
    "order_consistent",
    "max_intergenic_bp",
    "intervening_non_target_count",
    "species_tree_status",
    "target_parent_tree_support",
    "novel_insertion_evidence",
    "replacement_search_status",
    "uncharacterized_intervening_count"
]

def parse_args():
    parser = argparse.ArgumentParser(description="Stage 4: Fitch Ancestral-State Reconstruction, 2D Operon Fate Space (OFI) Classification and final report compilation.")
    parser.add_argument("--target_name", required=True, help="Target genome sample ID")
    parser.add_argument("--target_taxonomy", required=True, help="Target GTDB taxonomy string")
    parser.add_argument("--target_completeness", type=float, default=100.0, help="CheckM completeness")
    parser.add_argument("--target_contamination", type=float, default=0.0, help="CheckM contamination")
    parser.add_argument("--refs_tsv", required=True, help="Stage 0 TSV of selected and validated references")
    parser.add_argument("--topology_tsv", required=True, help="Stage 1 anchor topology TSV")
    parser.add_argument("--cluster_tsv", required=True, help="Stage 2 cluster rescue summary TSV")
    parser.add_argument("--locus_tsv", required=True, help="Stage 3 locus absence validation TSV")
    parser.add_argument("--tree_tool_env", default="", help="Conda environment for iqtree / fasttree")
    parser.add_argument("--species_tree_nwk", default="", help="Independently inferred, rooted species tree containing target and references")
    parser.add_argument("--species_tree_outgroup", default="", help="Sequence-based outgroup tip used to root an unrooted Newick tree")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads")
    parser.add_argument("--module_id", default="M00159", help="ATP synthase module ID (e.g. M00159, M00157)")
    parser.add_argument("--module_name", default="", help="ATP synthase module human-readable name")
    # OFI Threshold Parameters
    parser.add_argument("--localizable_locator_min", type=float, default=0.4, help="Minimum locator_score to localize candidate")
    parser.add_argument("--bilateral_context_min", type=float, default=0.5, help="Minimum C_bilateral for loss hard gate")
    parser.add_argument("--target_retained_min", type=float, default=0.6, help="Minimum T_local for target retention")
    parser.add_argument("--known_replacement_min", type=float, default=0.6, help="Minimum R_known for known replacement")
    parser.add_argument("--replacement_absence_max", type=float, default=0.2, help="Maximum R_known permitted for clean/true loss")
    parser.add_argument("--novel_insertion_min", type=float, default=0.6, help="Minimum I_novel for novel insertion loss")
    parser.add_argument("--loss_confidence_min", type=float, default=0.6, help="Minimum K_loss for definitive loss call")
    parser.add_argument("--clean_closure_min", type=float, default=0.6, help="Minimum K_clean for clean loss call")
    parser.add_argument("--out_summary_tsv", required=True, help="Output final summary table TSV")
    parser.add_argument("--out_tree_nwk", required=True, help="Output local phylogenetic tree in Newick format")
    return parser.parse_args()

class TreeNode:
    def __init__(self, name="", left=None, right=None, state=None):
        self.name = name
        self.left = left
        self.right = right
        self.state_set = set() if state is None else {state}
        self.final_state = None

def run_fitch_parsimony(root):
    """
    Executes Fitch Parsimony Ancestral State Reconstruction for binary state (0=absent, 1=present).
    Returns (parent_of_target_state, total_loss_events, total_gain_events).
    """
    def postorder(node):
        if node is None:
            return
        if node.left is None and node.right is None:
            return
        postorder(node.left)
        postorder(node.right)

        intersect = node.left.state_set.intersection(node.right.state_set)
        if intersect:
            node.state_set = intersect
        else:
            node.state_set = node.left.state_set.union(node.right.state_set)

    postorder(root)

    def resolve(root_state):
        loss_events = 0
        gain_events = 0
        root.final_state = root_state
        def preorder(node):
            nonlocal loss_events, gain_events
            for child in (node.left, node.right):
                if child is None:
                    continue
                child.final_state = node.final_state if node.final_state in child.state_set else next(iter(child.state_set))
                if node.final_state == "1" and child.final_state == "0":
                    loss_events += 1
                elif node.final_state == "0" and child.final_state == "1":
                    gain_events += 1
                preorder(child)
        preorder(root)
        parent_state = "?"
        def find_parent(node):
            nonlocal parent_state
            if node is None:
                return
            if (node.left and node.left.name == "Target") or (node.right and node.right.name == "Target"):
                parent_state = node.final_state
                return
            find_parent(node.left)
            find_parent(node.right)
        find_parent(root)
        return parent_state, loss_events, gain_events

    scenarios = [resolve(state) for state in sorted(root.state_set)]
    if len({scenario[0] for scenario in scenarios}) > 1:
        return "?", 0, 0
    return min(scenarios, key=lambda scenario: scenario[1] + scenario[2])

def build_local_topology_tree(target_name, pos_refs, neg_refs, out_nwk):
    """
    Legacy synthetic topology helper retained for older unit tests; not used by the workflow.
    Returns (root_node, positive_lineage_count).
    """
    pos_nodes = [TreeNode(name=p["accession"], state="1") for p in pos_refs]
    neg_nodes = [TreeNode(name=n["accession"], state="0") for n in neg_refs]
    target_node = TreeNode(name="Target", state="0")

    if not pos_nodes:
        root = TreeNode(name="Root", left=target_node, right=TreeNode(name="Outgroup", state="0"))
        with open(out_nwk, "w") as f: f.write("(Target:0.1,Outgroup:0.2);\n")
        return root, 0

    lineage_clusters = {}
    for p in pos_refs:
        parts = p.get("taxonomy", "").split(";")
        g = [x for x in parts if "g__" in x]
        l_name = g[0] if g else p["accession"]
        if l_name not in lineage_clusters:
            lineage_clusters[l_name] = []
        lineage_clusters[l_name].append(p["accession"])

    positive_lineage_count = len(lineage_clusters)

    if len(pos_nodes) >= 2:
        sister_pos = pos_nodes[0]
        other_pos = pos_nodes[1]
        internal_target_sister = TreeNode(name="Parent_Target_Sister", left=target_node, right=sister_pos)
        pos_clade = TreeNode(name="Pos_Lineage_Clade", left=internal_target_sister, right=other_pos)
        outgroup = neg_nodes[0] if neg_nodes else TreeNode(name="Outgroup", state="0")
        root = TreeNode(name="Root", left=pos_clade, right=outgroup)
        nwk_str = f"(({target_name}:0.05,{sister_pos.name}:0.05)Parent_Target:0.04,{other_pos.name}:0.08)Clade:0.06,{outgroup.name}:0.15);\n"
    else:
        root = TreeNode(name="Root", left=TreeNode(name="Parent_Target", left=target_node, right=pos_nodes[0]), right=TreeNode(name="Outgroup", state="0"))
        nwk_str = f"(({target_name}:0.05,{pos_nodes[0].name}:0.05)Parent_Target:0.05,Outgroup:0.15);\n"

    with open(out_nwk, "w") as f:
        f.write(nwk_str)

    return root, positive_lineage_count

def load_independent_species_tree(path, target_name, pos_refs, neg_refs, out_nwk, outgroup=""):
    """Load observed species topology; do not construct topology from phenotype labels."""
    if not path:
        with open(out_nwk, "w") as handle:
            handle.write("# No independently inferred species tree supplied\n")
        return None, 0, 0.0
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Configured species tree does not exist: {path}")
    from Bio import Phylo
    tree = Phylo.read(path, "newick")
    required = {target_name} | {str(r["accession"]) for r in pos_refs + neg_refs}
    present = {terminal.name for terminal in tree.get_terminals()}
    if not required.issubset(present):
        raise ValueError("Species tree must contain target plus every validated reference")
    if outgroup:
        if outgroup not in present or outgroup == target_name:
            raise ValueError("Specified outgroup is absent from species tree or is the target")
        tree.root_with_outgroup(outgroup)
    elif not tree.rooted:
        raise ValueError("Supply a rooted species tree or --species_tree_outgroup")
    for terminal in list(tree.get_terminals()):
        if terminal.name not in required:
            tree.prune(terminal)
    states = {target_name: "0"}
    states.update({str(r["accession"]): "1" for r in pos_refs})
    states.update({str(r["accession"]): "0" for r in neg_refs})
    def convert(clade):
        if clade.is_terminal():
            return TreeNode("Target" if clade.name == target_name else clade.name, state=states[clade.name])
        if len(clade.clades) != 2:
            raise ValueError("Polytomies must be resolved or left as unresolved phylogeny")
        return TreeNode(clade.name or "Internal", convert(clade.clades[0]), convert(clade.clades[1]))
    root = convert(tree.root)
    parent_map = {child: clade for clade in tree.find_clades() for child in clade.clades}
    target_tip = next(t for t in tree.get_terminals() if t.name == target_name)
    parent = parent_map.get(target_tip)
    support = float(parent.confidence) if parent and parent.confidence is not None else 0.0
    if 0.0 < support <= 1.0:
        support *= 100.0
    Phylo.write(tree, out_nwk, "newick")
    genera = {next((p for p in str(r.get("taxonomy", "")).split(";") if p.startswith("g__") and len(p) > 3), r["accession"]) for r in pos_refs}
    return root, len(genera), support

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out_summary_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_tree_nwk), exist_ok=True)

    # 1. Load upstream data
    refs_df = pd.read_csv(args.refs_tsv, sep="\t") if os.path.exists(args.refs_tsv) and os.path.getsize(args.refs_tsv) > 0 else pd.DataFrame()
    topo_df = pd.read_csv(args.topology_tsv, sep="\t") if os.path.exists(args.topology_tsv) and os.path.getsize(args.topology_tsv) > 0 else pd.DataFrame()
    cluster_df = pd.read_csv(args.cluster_tsv, sep="\t") if os.path.exists(args.cluster_tsv) and os.path.getsize(args.cluster_tsv) > 0 else pd.DataFrame()
    locus_df = pd.read_csv(args.locus_tsv, sep="\t") if os.path.exists(args.locus_tsv) and os.path.getsize(args.locus_tsv) > 0 else pd.DataFrame()

    topo = topo_df.iloc[0].to_dict() if not topo_df.empty else {}
    cluster = cluster_df.iloc[0].to_dict() if not cluster_df.empty else {}
    locus = locus_df.iloc[0].to_dict() if not locus_df.empty else {}

    pos_refs = refs_df[refs_df["validation_status"] == "CONFIRMED_POSITIVE"].to_dict("records") if not refs_df.empty and "validation_status" in refs_df.columns else []
    neg_refs = refs_df[refs_df["validation_status"] == "CONFIRMED_NEGATIVE"].to_dict("records") if not refs_df.empty and "validation_status" in refs_df.columns else []

    nearest_pos_acc = topo.get("nearest_positive_ref")
    if not nearest_pos_acc or nearest_pos_acc == "NONE" or (not refs_df.empty and nearest_pos_acc not in refs_df["accession"].values):
        nearest_pos_acc = pos_refs[0]["accession"] if pos_refs else "NONE"

    matched_ref = refs_df[refs_df["accession"] == nearest_pos_acc] if not refs_df.empty and "accession" in refs_df.columns else pd.DataFrame()
    if not matched_ref.empty:
        r_row = matched_ref.iloc[0]
        nearest_ani = r_row.get("ANI", 0.0)
        nearest_af_q = r_row.get("AF_query", 0.0)
        nearest_af_ref = r_row.get("AF_ref", 0.0)
        nearest_af_min = r_row.get("AF_min", 0.0)
    else:
        nearest_ani = pos_refs[0].get("ANI", 0.0) if pos_refs else 0.0
        nearest_af_q = pos_refs[0].get("AF_query", 0.0) if pos_refs else 0.0
        nearest_af_ref = pos_refs[0].get("AF_ref", 0.0) if pos_refs else 0.0
        nearest_af_min = pos_refs[0].get("AF_min", 0.0) if pos_refs else 0.0
    pos_tax_count = refs_df["positive_taxonomic_group_count"].iloc[0] if (not refs_df.empty and "positive_taxonomic_group_count" in refs_df.columns) else len(pos_refs)

    # 2. Build local phylogenetic tree & calculate positive_lineage_count
    root_node, tree_positive_lineages, target_parent_support = load_independent_species_tree(
        args.species_tree_nwk, args.target_name, pos_refs, neg_refs, args.out_tree_nwk,
        args.species_tree_outgroup)
    parent_state, loss_events, gain_events = run_fitch_parsimony(root_node) if root_node else ("?", 0, 0)

    # --------------------------------------------------------------------------
    # 4. Operon Fate Space (OFI) Quantitative Classification Engine
    # --------------------------------------------------------------------------
    locator_score = safe_float(topo.get("locator_score", 0.0))
    bilateral_context_score = safe_float(topo.get("bilateral_context_score", 0.0))
    syntenic_geometry = safe_str(topo.get("syntenic_geometry", "BILATERAL_SAME_CONTIG"))
    same_contig = safe_str(topo.get("same_contig", "YES"))

    T_local = safe_float(cluster.get("T_local", 0.0))
    T_global = safe_float(cluster.get("T_global", 0.0))
    R_known = safe_float(cluster.get("R_known", 0.0))
    I_novel = safe_float(cluster.get("I_novel", 0.0))
    novel_insertion_evidence = safe_str(cluster.get("novel_insertion_evidence", "NOT_ASSESSED"))
    replacement_search_status = safe_str(cluster.get("replacement_search_status", "INCOMPLETE"))
    A_local = normalize_tristate(cluster.get("A_local", "NA"))
    A_global = normalize_tristate(cluster.get("A_global", "NA"))

    intact_rescued_count = safe_int(cluster.get("intact_rescued_count", 0))
    pseudogene_count = safe_int(cluster.get("pseudogene_count", 0))
    frameshift_count = safe_int(cluster.get("frameshift_count", 0))
    premature_stop_count = safe_int(cluster.get("premature_stop_count", 0))
    pseudogene_detected = safe_str(cluster.get("pseudogene_detected", "NO"))
    expected_gene_count = safe_int(cluster.get("expected_gene_count", 8))

    B_bridge = safe_float(locus.get("B_bridge", 1.0 if same_contig == "YES" else 0.0))
    Q_continuity = safe_float(locus.get("Q_continuity", 1.0))
    D_closure = safe_float(locus.get("D_closure", 0.0))
    l_state = safe_str(locus.get("locus_state", "UNRESOLVED"))

    # Compute continuous OFI coordinates
    OFI_fate = round(T_local - R_known, 4)
    OFI_occupancy = round(max(T_local, R_known, I_novel), 4)

    # Execute Frozen State Machine (Precedence Order)
    loss_gate_status = "NOT_EVALUATED"
    K_loss_val = 0.0
    K_clean_val = 0.0

    if locator_score < args.localizable_locator_min:
        fate_category = "UNRESOLVED"
        loss_gate_status = "FAILED_LOCATOR_SCORE"
    elif syntenic_geometry == "SPLIT_CONTIG" or same_contig == "NO":
        fate_category = "ASSEMBLY_SPLIT_CANDIDATE"
        loss_gate_status = "FAILED_SPLIT_CONTIG"
    elif syntenic_geometry in ["UNILATERAL_SAME_CONTIG", "UNILATERAL_EDGE"] or bilateral_context_score == 0:
        fate_category = "ASSEMBLY_EDGE_CANDIDATE"
        loss_gate_status = "FAILED_UNILATERAL"
    # Bilateral ancestral locus evaluations:
    elif T_local >= args.target_retained_min and R_known >= args.known_replacement_min:
        fate_category = "COEXISTENCE"
        loss_gate_status = "NOT_APPLICABLE_COEXISTENCE"
    # [Core-3] Restore ANNOTATION_MISS for intact coding potential without disruptive mutations
    elif intact_rescued_count > 0 and (T_local >= args.target_retained_min or (T_local + (intact_rescued_count / max(1, expected_gene_count))) >= args.target_retained_min):
        fate_category = "ANNOTATION_MISS"
        loss_gate_status = "NOT_APPLICABLE_ANNOTATION_MISS"
    # Disrupted homologous remnant (pseudogene)
    elif pseudogene_detected == "YES" or (frameshift_count > 0 or premature_stop_count > 0 or pseudogene_count > 0):
        fate_category = "PSEUDOGENIZED"
        loss_gate_status = "NOT_APPLICABLE_PSEUDOGENIZED"
    elif T_local >= args.target_retained_min:
        fate_category = "TARGET_RETAINED"
        loss_gate_status = "NOT_APPLICABLE_RETAINED"
    elif T_local < 0.3 and T_global >= 0.7:
        fate_category = "RELOCATED"
        loss_gate_status = "NOT_APPLICABLE_RELOCATED"
    elif R_known >= args.known_replacement_min:
        fate_category = "KNOWN_REPLACEMENT"
        loss_gate_status = "NOT_APPLICABLE_REPLACEMENT"
    # Partial replacement / ambiguous occupancy
    elif R_known > args.replacement_absence_max:
        fate_category = "UNRESOLVED"
        loss_gate_status = "PARTIAL_REPLACEMENT_AMBIGUOUS"
    elif replacement_search_status != "COMPLETE":
        fate_category = "UNRESOLVED"
        loss_gate_status = "INCOMPLETE_REPLACEMENT_SEARCH"
    # [Core-7] Tri-state Negative Evidence Guard
    elif A_local == "NA" or A_global == "NA":
        fate_category = "UNRESOLVED"
        loss_gate_status = "INCOMPLETE_NEGATIVE_SEARCH"
        K_loss_val = "NA"
        K_clean_val = "NA"
    else:
        # Loss Hard Gates Pass Check:
        # B_bridge == 1, C_bilateral >= 0.5, R_known <= 0.2, A_local == "1", A_global == "1"
        gates_pass = (
            B_bridge == 1.0 and
            bilateral_context_score >= args.bilateral_context_min and
            R_known <= args.replacement_absence_max and
            A_local == "1" and
            A_global == "1"
        )
        if gates_pass:
            loss_gate_status = "PASSED"
            K_loss_val = round(math.sqrt(bilateral_context_score * Q_continuity), 4)
            K_clean_val = round(K_loss_val * D_closure, 4)

            if K_loss_val >= args.loss_confidence_min:
                # [Core-2] State Machine Priority: Novel Insertion strictly precedes Clean Closure!
                if I_novel >= args.novel_insertion_min:
                    fate_category = "LOSS_WITH_NOVEL_INSERTION" if novel_insertion_evidence == "VALIDATED" else "LOSS_WITH_UNCHARACTERIZED_OCCUPANCY"
                elif K_clean_val >= args.clean_closure_min and safe_str(locus.get("nucleotide_deletion_support", "NO")) == "YES":
                    fate_category = "CLEAN_LOSS"
                elif safe_str(locus.get("nucleotide_deletion_support", "NO")) != "YES":
                    fate_category = "LOSS_WITH_UNRESOLVED_STRUCTURE"
                else:
                    fate_category = "LOSS_WITH_LOCUS_REMODELING"
            else:
                fate_category = "UNRESOLVED"
        else:
            loss_gate_status = "FAILED_GATES"
            K_loss_val = 0.0
            K_clean_val = 0.0
            fate_category = "UNRESOLVED"

    # Evolutionary state & confidence determination
    confidence = "HIGH"
    if not pos_refs or (len(pos_refs) == 1 and pos_refs[0].get("accession") in ["NONE", "None", ""]):
        evolutionary_state = "UNRESOLVED_PHYLOGENY"
        confidence = "LOW"
    elif fate_category in ["CLEAN_LOSS", "LOSS_WITH_NOVEL_INSERTION", "LOSS_WITH_UNCHARACTERIZED_OCCUPANCY", "LOSS_WITH_LOCUS_REMODELING", "LOSS_WITH_UNRESOLVED_STRUCTURE"]:
        if root_node is None:
            evolutionary_state = "UNRESOLVED_PHYLOGENY"
            confidence = "LOW"
        elif parent_state == "1" and tree_positive_lineages >= 2 and target_parent_support >= 95:
            evolutionary_state = "INFERRED_LOSS_HIGH"
            confidence = "HIGH"
        elif parent_state == "1" and tree_positive_lineages == 1:
            evolutionary_state = "INFERRED_LOSS_MEDIUM"
            confidence = "MEDIUM"
        elif parent_state == "0":
            evolutionary_state = "ANCESTRAL_ABSENCE_POSSIBLE"
            confidence = "MEDIUM"
        elif parent_state == "?":
            evolutionary_state = "UNRESOLVED_ANCESTRAL_STATE"
            confidence = "LOW"
        else:
            evolutionary_state = "INFERRED_LOSS_MEDIUM"
            confidence = "MEDIUM"
    elif fate_category == "TARGET_RETAINED":
        evolutionary_state = "NOT_APPLICABLE"
        confidence = "HIGH"
    elif fate_category == "ANNOTATION_MISS":
        evolutionary_state = "RESCUED_ANNOTATION_MISS"
        confidence = "HIGH"
    elif fate_category == "PSEUDOGENIZED":
        evolutionary_state = "PSEUDOGENIZED"
        confidence = "HIGH"
    elif fate_category == "COEXISTENCE":
        evolutionary_state = "COEXISTENCE"
        confidence = "HIGH"
    elif fate_category == "RELOCATED":
        evolutionary_state = "RELOCATED"
        confidence = "MEDIUM"
    elif fate_category == "KNOWN_REPLACEMENT":
        evolutionary_state = "KNOWN_REPLACEMENT"
        confidence = "HIGH"
    elif fate_category in ["ASSEMBLY_EDGE_CANDIDATE", "ASSEMBLY_SPLIT_CANDIDATE"]:
        evolutionary_state = "UNRESOLVED"
        confidence = "LOW"
    else:
        evolutionary_state = "UNRESOLVED"
        confidence = "LOW"

    # 5. Compile final programmatic record
    summary_records = []
    seen_summary_refs = set()
    for _, t_row in topo_df.iterrows():
        ref_acc = str(t_row.get("nearest_positive_ref", ""))
        clean_ref_acc = ref_acc.replace("GB_", "").replace("RS_", "")
        if clean_ref_acc in seen_summary_refs:
            continue
        seen_summary_refs.add(clean_ref_acc)
        matched_ref = refs_df[refs_df["accession"] == ref_acc] if not refs_df.empty and "accession" in refs_df.columns else pd.DataFrame()
        if matched_ref.empty and not refs_df.empty and "accession" in refs_df.columns:
            matched_ref = refs_df[refs_df["accession"].str.replace("GB_", "").str.replace("RS_", "") == clean_ref_acc]

        if not matched_ref.empty:
            r_row = matched_ref.iloc[0]
            ref_ani = r_row.get("ANI", 0.0)
            ref_af_q = r_row.get("AF_query", 0.0)
            ref_af_ref = r_row.get("AF_ref", 0.0)
            ref_af_min = r_row.get("AF_min", 0.0)
        else:
            ref_ani = nearest_ani
            ref_af_q = nearest_af_q
            ref_af_ref = nearest_af_ref
            ref_af_min = nearest_af_min

        row_rec = {
            "target_id": args.target_name,
            "module_id": args.module_id,
            "module_name": args.module_name or args.module_id,
            "GTDB_taxonomy": args.target_taxonomy,
            "completeness": args.target_completeness,
            "contamination": args.target_contamination,
            "nearest_positive_ref": ref_acc or nearest_pos_acc,
            "positive_ref_count": len(pos_refs),
            "positive_taxonomic_group_count": pos_tax_count,
            "positive_lineage_count": tree_positive_lineages,
            "negative_ref_count": len(neg_refs),
            "ANI": ref_ani,
            "AF_query": ref_af_q,
            "AF_ref": ref_af_ref,
            "AF_min": ref_af_min,
            "left_anchor_1": t_row.get("left_anchor_1", topo.get("left_anchor_1", "YES")),
            "left_anchor_2": t_row.get("left_anchor_2", topo.get("left_anchor_2", "YES")),
            "right_anchor_1": t_row.get("right_anchor_1", topo.get("right_anchor_1", "YES")),
            "right_anchor_2": t_row.get("right_anchor_2", topo.get("right_anchor_2", "YES")),
            "left_anchors_matched_count": safe_int(t_row.get("left_anchors_matched_count", topo.get("left_anchors_matched_count", 0))),
            "right_anchors_matched_count": safe_int(t_row.get("right_anchors_matched_count", topo.get("right_anchors_matched_count", 0))),
            "left_matched_genes": safe_str(t_row.get("left_matched_genes", topo.get("left_matched_genes", ""))),
            "right_matched_genes": safe_str(t_row.get("right_matched_genes", topo.get("right_matched_genes", ""))),
            "left_anchor_region_other_genes": safe_int(t_row.get("left_anchor_region_other_genes", topo.get("left_anchor_region_other_genes", 0))),
            "right_anchor_region_other_genes": safe_int(t_row.get("right_anchor_region_other_genes", topo.get("right_anchor_region_other_genes", 0))),
            "intervening_genes_count": safe_int(t_row.get("intervening_genes_count", topo.get("intervening_genes_count", 0))),
            "intervening_genes": safe_str(t_row.get("intervening_genes", topo.get("intervening_genes", ""))),
            "target_locus_total_genes": safe_int(t_row.get("target_locus_total_genes", topo.get("target_locus_total_genes", 0))),
            "target_locus_other_genes": safe_int(t_row.get("target_locus_other_genes", topo.get("target_locus_other_genes", 0))),
            "anchor_left_identity": t_row.get("anchor_left_identity", topo.get("anchor_left_identity", 0.0)),
            "anchor_right_identity": t_row.get("anchor_right_identity", topo.get("anchor_right_identity", 0.0)),
            "anchor_order_consistent": t_row.get("anchor_order_consistent", topo.get("anchor_order_consistent", "YES")),
            "anchor_orientation_consistent": t_row.get("anchor_orientation_consistent", topo.get("anchor_orientation_consistent", "YES")),
            "same_contig": t_row.get("same_contig", topo.get("same_contig", "YES")),
            "target_contig_id": safe_str(t_row.get("target_contig", topo.get("target_contig", ""))),
            "distance_to_left_edge": t_row.get("distance_to_left_edge", topo.get("distance_to_left_edge", 0)),
            "distance_to_right_edge": t_row.get("distance_to_right_edge", topo.get("distance_to_right_edge", 0)),
            "required_edge_margin": t_row.get("required_edge_margin", topo.get("required_edge_margin", 15000)),
            "expected_locus_length": t_row.get("expected_locus_length", topo.get("expected_locus_length", 15000)),
            "observed_anchor_gap": t_row.get("observed_anchor_gap", topo.get("observed_anchor_gap", 0)),
            "target_contig_length": t_row.get("target_contig_length", topo.get("target_contig_length", 0)),
            "expected_gene_count": expected_gene_count,
            "protein_detected_count": safe_int(cluster.get("protein_detected_count", 0)),
            "location_unknown_count": safe_int(cluster.get("location_unknown_count", 0)),
            "blastx_rescued_count": safe_int(cluster.get("blastx_rescued_count", 0)),
            "frameshift_count": frameshift_count,
            "premature_stop_count": premature_stop_count,
            "mutation_evidence_status": cluster.get("mutation_evidence_status", "NOT_ASSESSED"),
            "annotation_issue": cluster.get("annotation_issue", "NO"),
            "pseudogene_detected": pseudogene_detected,
            "member_gene_states": cluster.get("member_gene_states", ""),
            "cluster_state": cluster.get("cluster_state", "FULLY_ABSENT"),
            "structural_event_type": locus.get("structural_event_type", "UNRESOLVED"),
            "nucleotide_deletion_support": locus.get("nucleotide_deletion_support", "NO"),
            "junction_support_status": locus.get("junction_support_status", "NO_READS"),
            "junction_spanning_reads": locus.get("junction_spanning_reads", 0),
            "junction_mean_mapq": locus.get("junction_mean_mapq", 0.0),
            "junction_coverage_ratio": locus.get("junction_coverage_ratio", 0.0),
            "parent_ancestral_state": parent_state,
            "parsimony_loss_events": loss_events,
            "locus_state": l_state,
            "evolutionary_state": evolutionary_state,
            "confidence": confidence,
            # OFI Metrics
            "locator_score": locator_score,
            "bilateral_context_score": bilateral_context_score,
            "T_local": T_local,
            "T_global": T_global,
            "R_known": R_known,
            "I_novel": I_novel,
            "A_local": A_local,
            "A_global": A_global,
            "B_bridge": B_bridge,
            "Q_continuity": Q_continuity,
            "D_closure": D_closure,
            "K_loss": K_loss_val,
            "K_clean": K_clean_val,
            "loss_gate_status": loss_gate_status,
            "OFI_fate": OFI_fate,
            "OFI_occupancy": OFI_occupancy,
            "fate_category": fate_category
        }
        summary_records.append(row_rec)

    if not summary_records:
        summary_records.append({
            "target_id": args.target_name,
            "module_id": args.module_id,
            "module_name": args.module_name or args.module_id,
            "GTDB_taxonomy": args.target_taxonomy,
            "completeness": args.target_completeness,
            "contamination": args.target_contamination,
            "nearest_positive_ref": nearest_pos_acc,
            "positive_ref_count": len(pos_refs),
            "positive_taxonomic_group_count": pos_tax_count,
            "positive_lineage_count": tree_positive_lineages,
            "negative_ref_count": len(neg_refs),
            "ANI": nearest_ani,
            "AF_query": nearest_af_q,
            "AF_ref": nearest_af_ref,
            "AF_min": nearest_af_min,
            "left_anchor_1": topo.get("left_anchor_1", "YES"),
            "left_anchor_2": topo.get("left_anchor_2", "YES"),
            "right_anchor_1": topo.get("right_anchor_1", "YES"),
            "right_anchor_2": topo.get("right_anchor_2", "YES"),
            "left_anchors_matched_count": safe_int(topo.get("left_anchors_matched_count", 0)),
            "right_anchors_matched_count": safe_int(topo.get("right_anchors_matched_count", 0)),
            "left_matched_genes": safe_str(topo.get("left_matched_genes", "")),
            "right_matched_genes": safe_str(topo.get("right_matched_genes", "")),
            "left_anchor_region_other_genes": safe_int(topo.get("left_anchor_region_other_genes", 0)),
            "right_anchor_region_other_genes": safe_int(topo.get("right_anchor_region_other_genes", 0)),
            "intervening_genes_count": safe_int(topo.get("intervening_genes_count", 0)),
            "intervening_genes": safe_str(topo.get("intervening_genes", "")),
            "target_locus_total_genes": safe_int(topo.get("target_locus_total_genes", 0)),
            "target_locus_other_genes": safe_int(topo.get("target_locus_other_genes", 0)),
            "anchor_left_identity": topo.get("anchor_left_identity", 0.0),
            "anchor_right_identity": topo.get("anchor_right_identity", 0.0),
            "anchor_order_consistent": topo.get("anchor_order_consistent", "YES"),
            "anchor_orientation_consistent": topo.get("anchor_orientation_consistent", "YES"),
            "same_contig": topo.get("same_contig", "YES"),
            "target_contig_id": safe_str(topo.get("target_contig", "")),
            "distance_to_left_edge": topo.get("distance_to_left_edge", 0),
            "distance_to_right_edge": topo.get("distance_to_right_edge", 0),
            "required_edge_margin": topo.get("required_edge_margin", 15000),
            "expected_locus_length": topo.get("expected_locus_length", 15000),
            "observed_anchor_gap": topo.get("observed_anchor_gap", 0),
            "target_contig_length": topo.get("target_contig_length", 0),
            "expected_gene_count": expected_gene_count,
            "protein_detected_count": safe_int(cluster.get("protein_detected_count", 0)),
            "location_unknown_count": safe_int(cluster.get("location_unknown_count", 0)),
            "blastx_rescued_count": safe_int(cluster.get("blastx_rescued_count", 0)),
            "frameshift_count": frameshift_count,
            "premature_stop_count": premature_stop_count,
            "mutation_evidence_status": cluster.get("mutation_evidence_status", "NOT_ASSESSED"),
            "annotation_issue": cluster.get("annotation_issue", "NO"),
            "pseudogene_detected": pseudogene_detected,
            "member_gene_states": cluster.get("member_gene_states", ""),
            "cluster_state": cluster.get("cluster_state", "FULLY_ABSENT"),
            "structural_event_type": locus.get("structural_event_type", "UNRESOLVED"),
            "nucleotide_deletion_support": locus.get("nucleotide_deletion_support", "NO"),
            "junction_support_status": locus.get("junction_support_status", "NO_READS"),
            "junction_spanning_reads": locus.get("junction_spanning_reads", 0),
            "junction_mean_mapq": locus.get("junction_mean_mapq", 0.0),
            "junction_coverage_ratio": locus.get("junction_coverage_ratio", 0.0),
            "parent_ancestral_state": parent_state,
            "parsimony_loss_events": loss_events,
            "locus_state": l_state,
            "evolutionary_state": evolutionary_state,
            "confidence": confidence,
            # OFI Metrics
            "locator_score": locator_score,
            "bilateral_context_score": bilateral_context_score,
            "T_local": T_local,
            "T_global": T_global,
            "R_known": R_known,
            "I_novel": I_novel,
            "A_local": A_local,
            "A_global": A_global,
            "B_bridge": B_bridge,
            "Q_continuity": Q_continuity,
            "D_closure": D_closure,
            "K_loss": K_loss_val,
            "K_clean": K_clean_val,
            "loss_gate_status": loss_gate_status,
            "OFI_fate": OFI_fate,
            "OFI_occupancy": OFI_occupancy,
            "fate_category": fate_category
        })

    for record in summary_records:
        for field in ("operon_structure_state", "operon_prediction", "same_strand",
                      "order_consistent", "max_intergenic_bp", "intervening_non_target_count"):
            record[field] = cluster.get(field, "UNRESOLVED" if field.startswith("operon_") else "NA")
        record["species_tree_status"] = "INDEPENDENT_TREE" if root_node else "NOT_PROVIDED"
        record["target_parent_tree_support"] = target_parent_support if root_node else "NA"
        record["novel_insertion_evidence"] = novel_insertion_evidence
        record["replacement_search_status"] = replacement_search_status
        record["uncharacterized_intervening_count"] = cluster.get("uncharacterized_intervening_count", "NA")
    final_df = pd.DataFrame(summary_records)[SCHEMA_COLUMNS]
    final_df.to_csv(args.out_summary_tsv, sep="\t", index=False)
    print(f"Stage 4 OFI Inference Completed: fate_category={fate_category}, OFI_fate={OFI_fate}, OFI_occupancy={OFI_occupancy}, K_loss={K_loss_val}, K_clean={K_clean_val}, loss_gate={loss_gate_status}", flush=True)

if __name__ == "__main__":
    main()
