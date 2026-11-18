#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import pandas as pd

MODULES_ORDER = ["M00159", "M00157", "M00158", "M00160"]

MODULE_NAMES = {
    "M00159": "V/A-type (Prokaryotic)",
    "M00157": "F-type (Prokaryotic)",
    "M00158": "F-type (Eukaryotic)",
    "M00160": "V-type (Eukaryotic)"
}

def parse_args():
    parser = argparse.ArgumentParser(description="Aggregates multi-module synteny validation summaries into a global stacked TSV and an ATP synthase pivot matrix.")
    parser.add_argument("--summary_tsvs", "--input_tsvs", nargs="*", default=[], help="List of individual Stage 4 summary TSVs")
    parser.add_argument("--search_dir", default="", help="Directory to recursively search for Stage 4 summary TSVs")
    parser.add_argument("--out_all_summary_tsv", "--out_stacked_tsv", required=True, help="Output concatenated summary TSV across all samples and modules")
    parser.add_argument("--out_matrix_tsv", required=True, help="Output pivot matrix TSV (one row per genome with columns for each module)")
    parser.add_argument("--out_plot_pdf", default="", help="Optional output PDF plot for 2D Operon Fate Space (OFI)")
    parser.add_argument("--out_plot_png", default="", help="Optional output PNG plot for 2D Operon Fate Space (OFI)")
    return parser.parse_args()

def determine_overall_atp_state(mod_dict):
    """
    Determines overall biological ATP synthase status across the 4 modules.
    Accepts either nested dict: { 'M00159': {'loss': ..., 'cluster': ..., 'fate': ...}, ... }
    or flat dict: { 'M00159_inferred_loss': ..., 'M00157_inferred_loss': ... }
    """
    # If neither M00159 nor M00157 is in mod_dict, this is a non-ATP operon or custom operon
    has_v = ("M00159_inferred_loss" in mod_dict) or ("M00159" in mod_dict)
    has_f = ("M00157_inferred_loss" in mod_dict) or ("M00157" in mod_dict)
    if not has_v and not has_f:
        if mod_dict:
            first_val = next(iter(mod_dict.values()))
            if isinstance(first_val, dict):
                return str(first_val.get("fate", first_val.get("loss", "OPERON_EVALUATED")))
            return "OPERON_EVALUATED"
        return "NOT_APPLICABLE"

    if "M00159_inferred_loss" in mod_dict:
        v_loss = str(mod_dict.get("M00159_inferred_loss", "UNRESOLVED"))
        v_clust = str(mod_dict.get("M00159_cluster_state", "UNRESOLVED"))
        v_fate = str(mod_dict.get("M00159_fate_category", v_loss))
        f_loss = str(mod_dict.get("M00157_inferred_loss", "UNRESOLVED"))
        f_clust = str(mod_dict.get("M00157_cluster_state", "UNRESOLVED"))
        f_fate = str(mod_dict.get("M00157_fate_category", f_loss))
    else:
        v_info = mod_dict.get("M00159", {})
        f_info = mod_dict.get("M00157", {})
        v_loss = v_info.get("loss", "UNRESOLVED")
        v_clust = v_info.get("cluster", "UNRESOLVED")
        v_fate = v_info.get("fate", v_loss)
        f_loss = f_info.get("loss", "UNRESOLVED")
        f_clust = f_info.get("cluster", "UNRESOLVED")
        f_fate = f_info.get("fate", f_loss)

    is_v_present = (v_clust in ["COMPLETE", "COMPLETE_AFTER_RESCUE", "PARTIAL_INTACT"] or 
                    v_loss in ["NOT_APPLICABLE", "INFERRED_LOSS_NO", "PRESENT", "RETAINED", "RESCUED_ANNOTATION_MISS"] or
                    v_fate in ["TARGET_RETAINED", "ANNOTATION_MISS", "COEXISTENCE"])
    is_f_present = (f_clust in ["COMPLETE", "COMPLETE_AFTER_RESCUE", "PARTIAL_INTACT"] or 
                    f_loss in ["NOT_APPLICABLE", "INFERRED_LOSS_NO", "PRESENT", "RETAINED", "RESCUED_ANNOTATION_MISS"] or
                    f_fate in ["TARGET_RETAINED", "ANNOTATION_MISS", "COEXISTENCE", "KNOWN_REPLACEMENT"])

    if is_v_present and is_f_present:
        return "DUAL_V_AND_F_TYPE_PRESENT"
    if is_v_present:
        return "CANONICAL_V_TYPE_RETAINED"
    if is_f_present:
        return "F_TYPE_SUBSTITUTED"

    # Both main prokaryotic types are missing or lost
    if "LOSS_HIGH" in v_loss or "LOSS_MEDIUM" in v_loss:
        if f_loss in ["ANCESTRAL_ABSENT", "NOT_APPLICABLE", "UNRESOLVED", "FULLY_ABSENT", "NOT_ANALYZED"] or f_clust in ["FULLY_ABSENT", "ANCESTRAL_ABSENT"]:
            return "TOTAL_ATP_SYNTHASE_LOSS"
        return "V_TYPE_LOST_WITHOUT_F_TYPE"

    if v_loss == "ANCESTRAL_ABSENT" and f_loss == "ANCESTRAL_ABSENT":
        return "COMPLETE_ATP_SYNTHASE_LOSS"

    if "PARTIAL_LOSS" in v_loss or "PARTIAL_LOSS" in f_loss or v_fate == "PSEUDOGENIZED":
        return "PARTIAL_DEGRADATION"

    if v_clust in ["FULLY_ABSENT", "ANCESTRAL_ABSENT"] and f_clust in ["FULLY_ABSENT", "ANCESTRAL_ABSENT"]:
        return "COMPLETE_ATP_SYNTHASE_LOSS"

    return "UNRESOLVED"

def compile_matrices(files, out_all_summary_tsv, out_matrix_tsv, out_plot_pdf="", out_plot_png=""):
    """
    Concatenates individual summary TSVs and generates wide pivot matrix TSV.
    """
    if not files:
        print("Notice: No summary TSVs provided or found. Writing empty outputs.", file=sys.stderr)
        os.makedirs(os.path.dirname(os.path.abspath(out_all_summary_tsv)), exist_ok=True)
        os.makedirs(os.path.dirname(os.path.abspath(out_matrix_tsv)), exist_ok=True)
        pd.DataFrame().to_csv(out_all_summary_tsv, sep="\t", index=False)
        pd.DataFrame().to_csv(out_matrix_tsv, sep="\t", index=False)
        return

    all_dfs = []
    for f in sorted(files):
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            continue
        try:
            df = pd.read_csv(f, sep="\t")
            if not df.empty:
                # Infer module_id if missing from column
                if "module_id" not in df.columns:
                    m_id = "M00159"
                    for cand_m in MODULES_ORDER:
                        if cand_m in f:
                            m_id = cand_m
                            break
                    df.insert(1, "module_id", m_id)
                    df.insert(2, "module_name", MODULE_NAMES.get(m_id, m_id))
                all_dfs.append(df)
        except Exception as e:
            print(f"Warning reading {f}: {e}", file=sys.stderr)

    if not all_dfs:
        print("Notice: All summary TSVs were empty.", file=sys.stderr)
        os.makedirs(os.path.dirname(os.path.abspath(out_all_summary_tsv)), exist_ok=True)
        os.makedirs(os.path.dirname(os.path.abspath(out_matrix_tsv)), exist_ok=True)
        pd.DataFrame().to_csv(out_all_summary_tsv, sep="\t", index=False)
        pd.DataFrame().to_csv(out_matrix_tsv, sep="\t", index=False)
        return

    combined_df = pd.concat(all_dfs, ignore_index=True)
    dedup_cols = [c for c in ["target_id", "module_id", "nearest_positive_ref"] if c in combined_df.columns]
    if dedup_cols:
        combined_df = combined_df.drop_duplicates(subset=dedup_cols)

    os.makedirs(os.path.dirname(os.path.abspath(out_all_summary_tsv)), exist_ok=True)
    combined_df.to_csv(out_all_summary_tsv, sep="\t", index=False)
    print(f"Wrote global stacked summary TSV: {out_all_summary_tsv} ({len(combined_df)} records)", flush=True)

    # Compile Pivot Matrix
    matrix_rows = []
    by_target = combined_df.groupby("target_id", sort=False)

    for target_id, group in by_target:
        first_row = group.iloc[0]
        rec = {
            "target_id": target_id,
            "GTDB_taxonomy": first_row.get("GTDB_taxonomy", ""),
            "completeness": first_row.get("completeness", 100.0),
            "contamination": first_row.get("contamination", 0.0),
        }

        mod_dict = {}
        active_modules = [m for m in MODULES_ORDER if m in group["module_id"].values] + [m for m in group["module_id"].unique() if m not in MODULES_ORDER]
        if not active_modules:
            active_modules = MODULES_ORDER
        for mod in active_modules:
            mod_rows = group[group["module_id"] == mod]
            if not mod_rows.empty:
                m_row = mod_rows.iloc[0]
                loss_val = m_row.get("evolutionary_state", m_row.get("inferred_evolutionary_loss", "UNRESOLVED"))
                conf_val = m_row.get("confidence", "LOW")
                clust_val = m_row.get("cluster_state", "UNRESOLVED")
                struct_val = m_row.get("structural_event_type", "UNRESOLVED")
                ref_val = m_row.get("nearest_positive_ref", "NONE")
                ctg_val = m_row.get("target_contig_id", m_row.get("target_contig", ""))
                fate_val = m_row.get("fate_category", loss_val)
                ofi_fate = m_row.get("OFI_fate", 0.0)
                ofi_occ = m_row.get("OFI_occupancy", 0.0)
                loc_score = m_row.get("locator_score", 0.0)
                bilat_score = m_row.get("bilateral_context_score", 0.0)
                k_loss = m_row.get("K_loss", 0.0)
                k_clean = m_row.get("K_clean", 0.0)

                rec[f"{mod}_inferred_loss"] = loss_val
                rec[f"{mod}_confidence"] = conf_val
                rec[f"{mod}_cluster_state"] = clust_val
                rec[f"{mod}_structural_event"] = struct_val
                rec[f"{mod}_nearest_ref"] = ref_val
                rec[f"{mod}_target_contig"] = ctg_val
                rec[f"{mod}_fate_category"] = fate_val
                rec[f"{mod}_operon_prediction"] = m_row.get("operon_prediction", "UNRESOLVED")
                rec[f"{mod}_operon_structure_state"] = m_row.get("operon_structure_state", "UNRESOLVED")
                rec[f"{mod}_OFI_fate"] = ofi_fate
                rec[f"{mod}_OFI_occupancy"] = ofi_occ
                rec[f"{mod}_locator_score"] = loc_score
                rec[f"{mod}_bilateral_context"] = bilat_score
                rec[f"{mod}_K_loss"] = k_loss
                rec[f"{mod}_K_clean"] = k_clean

                mod_dict[mod] = {
                    "loss": str(loss_val),
                    "cluster": str(clust_val),
                    "conf": str(conf_val),
                    "fate": str(fate_val)
                }
            else:
                rec[f"{mod}_inferred_loss"] = "NOT_ANALYZED"
                rec[f"{mod}_confidence"] = "NONE"
                rec[f"{mod}_cluster_state"] = "NOT_ANALYZED"
                rec[f"{mod}_structural_event"] = "NOT_ANALYZED"
                rec[f"{mod}_nearest_ref"] = "NONE"
                rec[f"{mod}_target_contig"] = ""
                rec[f"{mod}_fate_category"] = "NOT_ANALYZED"
                rec[f"{mod}_operon_prediction"] = "NOT_ANALYZED"
                rec[f"{mod}_operon_structure_state"] = "NOT_ANALYZED"
                rec[f"{mod}_OFI_fate"] = 0.0
                rec[f"{mod}_OFI_occupancy"] = 0.0
                rec[f"{mod}_locator_score"] = 0.0
                rec[f"{mod}_bilateral_context"] = 0.0
                rec[f"{mod}_K_loss"] = 0.0
                rec[f"{mod}_K_clean"] = 0.0

                mod_dict[mod] = {
                    "loss": "NOT_ANALYZED",
                    "cluster": "NOT_ANALYZED",
                    "conf": "NONE",
                    "fate": "NOT_ANALYZED"
                }

        rec["overall_operon_state"] = determine_overall_atp_state(mod_dict)
        rec["overall_atp_synthase_state"] = rec["overall_operon_state"]
        matrix_rows.append(rec)

    matrix_df = pd.DataFrame(matrix_rows)
    os.makedirs(os.path.dirname(os.path.abspath(out_matrix_tsv)), exist_ok=True)
    matrix_df.to_csv(out_matrix_tsv, sep="\t", index=False)
    print(f"Wrote ATP synthase synteny matrix: {out_matrix_tsv} ({len(matrix_df)} genomes)", flush=True)

    # Automatically generate 2D Operon Fate Space (OFI) plot if outputs provided
    if out_plot_pdf or out_plot_png:
        plot_operon_fate_space(combined_df, out_plot_pdf, out_plot_png)

def plot_operon_fate_space(combined_df, out_pdf="", out_png=""):
    """
    Renders 2D Operon Fate Space (OFI) landscape scatter plot.
    X-axis: OFI_fate = T_local - R_known ([-1.0, 1.0])
    Y-axis: OFI_occupancy = max(T_local, R_known, I_novel) ([0.0, 1.0])
    """
    if combined_df.empty or "OFI_fate" not in combined_df.columns:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        plot_df = combined_df.dropna(subset=["OFI_fate", "OFI_occupancy"]).copy()
        if plot_df.empty:
            return

        fig, ax = plt.subplots(figsize=(9, 7), dpi=300)

        # Standard category colors
        color_map = {
            "TARGET_RETAINED": "#2ca02c",             # Green
            "COEXISTENCE": "#17becf",                 # Cyan
            "ANNOTATION_MISS": "#8c564b",             # Brown
            "PSEUDOGENIZED": "#ff7f0e",               # Orange
            "KNOWN_REPLACEMENT": "#9467bd",           # Purple
            "LOSS_WITH_NOVEL_INSERTION": "#d62728",   # Red
            "CLEAN_LOSS": "#8c1515",                  # Dark red
            "LOSS_WITH_LOCUS_REMODELING": "#e377c2",  # Pink
            "RELOCATED": "#1f77b4",                   # Blue
            "ASSEMBLY_SPLIT_CANDIDATE": "#7f7f7f",   # Gray
            "ASSEMBLY_EDGE_CANDIDATE": "#bcbd22",    # Olive
            "UNRESOLVED": "#d3d3d3",                 # Light gray
        }

        # Quadrant background guides
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.axhline(0.6, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
        ax.axvline(0.6, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
        ax.axvline(-0.6, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)

        # Region quadrant labels
        ax.text(0.85, 0.96, "Target Retained\n(T_local ≥ 0.6)", ha="center", va="top", fontsize=8, color="#2ca02c", fontweight="bold", alpha=0.85)
        ax.text(-0.85, 0.96, "Known Replacement\n(R_known ≥ 0.6)", ha="center", va="top", fontsize=8, color="#9467bd", fontweight="bold", alpha=0.85)
        ax.text(0.0, 0.96, "Novel Insertion\n(I_novel ≥ 0.6)", ha="center", va="top", fontsize=8, color="#d62728", fontweight="bold", alpha=0.85)
        ax.text(0.0, 0.05, "Clean Loss\n(K_clean ≥ 0.6)", ha="center", va="bottom", fontsize=8, color="#8c1515", fontweight="bold", alpha=0.85)

        # Scatter points by category
        for cat, group in plot_df.groupby("fate_category"):
            c = color_map.get(cat, "#333333")
            # Slight jitter for visualization if multiple exact overlaps
            jitter_x = (np.random.rand(len(group)) - 0.5) * 0.02
            jitter_y = (np.random.rand(len(group)) - 0.5) * 0.02
            ax.scatter(group["OFI_fate"].astype(float) + jitter_x, 
                       group["OFI_occupancy"].astype(float) + jitter_y, 
                       color=c, label=f"{cat} (n={len(group)})",
                       s=45, alpha=0.75, edgecolors="none")

        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-0.05, 1.08)
        ax.set_xlabel("Operon Fate Index (OFI_fate = T_local - R_known)", fontsize=11, fontweight="bold")
        ax.set_ylabel("Operon Occupancy Index (OFI_occupancy)", fontsize=11, fontweight="bold")
        ax.set_title("2D Operon Fate Space (OFI) Multi-Genome Landscape", fontsize=13, fontweight="bold", pad=12)

        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, fontsize=8)
        plt.tight_layout()

        if out_pdf:
            os.makedirs(os.path.dirname(os.path.abspath(out_pdf)), exist_ok=True)
            plt.savefig(out_pdf, bbox_inches="tight")
            print(f"Saved Operon Fate Space PDF plot: {out_pdf}", flush=True)
        if out_png:
            os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
            plt.savefig(out_png, bbox_inches="tight", dpi=300)
            print(f"Saved Operon Fate Space PNG plot: {out_png}", flush=True)
        plt.close()
    except Exception as e:
        print(f"Notice: Operon Fate Space plotting skipped: {e}", file=sys.stderr)

def main():
    args = parse_args()

    files = list(args.summary_tsvs)
    if args.search_dir and os.path.isdir(args.search_dir):
        pattern = os.path.join(args.search_dir, "**", "*locus_validation_summary.tsv")
        found = glob.glob(pattern, recursive=True)
        for f in found:
            if f not in files and not f.endswith("all_modules_summary.tsv") and not f.endswith("atp_synthase_synteny_matrix.tsv"):
                files.append(f)

    compile_matrices(files, args.out_all_summary_tsv, args.out_matrix_tsv, args.out_plot_pdf, args.out_plot_png)

if __name__ == "__main__":
    main()
