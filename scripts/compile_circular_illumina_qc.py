#!/usr/bin/env python3
"""
compile_circular_illumina_qc.py
===============================
Merges results from:
1. Track A: coverm genome output (circular_genome_abundance.tsv)
2. Track B: coverm contig standard & strict outputs (circular_contigs.coverage.tsv & strict.tsv)
3. Track C: rotated junction QC output (rotated_junction_qc.tsv)

Computes diagnostic ratios:
- Competitive depth retention ratio: C = contig_trimmed / genome_trimmed
- Sensitivity ratio: S = strict_trimmed / standard_trimmed
- Robustness ratio: trimmed_to_mean_ratio = genome_trimmed / genome_mean

Evaluates 8 orthogonal boolean flags and assigns overall_status.
Produces: circular_illumina_qc_summary.tsv
"""

import os
import sys
import re
import argparse


def clean_val(v, default=0.0):
    if v is None:
        return default
    v_str = str(v).strip().rstrip("%")
    if v_str == "" or v_str.upper() in ("NA", "NAN", "NULL", "NONE"):
        return default
    try:
        return float(v_str)
    except ValueError:
        return default


def parse_coverm_tsv(tsv_path):
    """Parses CoverM TSV into dict {entity_id: {col_name: value}}."""
    data = {}
    if not tsv_path or not os.path.exists(tsv_path):
        return data

    with open(tsv_path) as f:
        header_line = f.readline()
        if not header_line:
            return data
        headers = [h.strip() for h in header_line.split("\t")]
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split("\t")]
            row = dict(zip(headers, parts))
            # Entity name is first column
            entity = parts[0]
            # Strip file extension if genome mode outputs filename
            for ext in (".fa", ".fasta", ".fna"):
                if entity.endswith(ext):
                    entity = entity[:-len(ext)]
                    break
            data[entity] = row
    return data


def parse_coverm_genome_tsv(genome_tsv):
    """
    Parses Track A coverm genome output.
    Finds sample-level unmapped / mapped percent, and per-genome metrics.
    """
    raw_data = parse_coverm_tsv(genome_tsv)
    unmapped_pct = 0.0
    mapped_pct = 100.0

    # Look for 'unmapped' row
    for k, row in raw_data.items():
        if k.lower() in ("unmapped", "unmapped_reads"):
            # find relative abundance column
            for col, val in row.items():
                if "relative" in col.lower() or "abundance" in col.lower():
                    unmapped_pct = clean_val(val, 0.0)
                    mapped_pct = max(0.0, 100.0 - unmapped_pct)
                    break

    genome_metrics = {}
    for gid, row in raw_data.items():
        if gid.lower() in ("unmapped", "unmapped_reads"):
            continue

        ra = 0.0
        mean_d = 0.0
        trim_d = 0.0
        cov_frac = 0.0
        reads_cnt = 0
        anir = 0.0

        for col, val in row.items():
            cl = col.lower()
            if "relative" in cl or "abundance" in cl:
                ra = clean_val(val, 0.0)
            elif "trimmed" in cl:
                trim_d = clean_val(val, 0.0)
            elif "mean" in cl and "trimmed" not in cl:
                mean_d = clean_val(val, 0.0)
            elif "covered" in cl or "fraction" in cl:
                cov_frac = clean_val(val, 0.0)
            elif "count" in cl or "reads" in cl:
                reads_cnt = int(clean_val(val, 0))
            elif "anir" in cl:
                anir_raw = clean_val(val, 0.0)
                anir = anir_raw * 100.0 if anir_raw <= 1.0 else anir_raw

        genome_metrics[gid] = {
            "relative_abundance_percent": round(ra, 3),
            "genome_mean_depth": round(mean_d, 2),
            "genome_trimmed_mean_depth": round(trim_d, 2),
            "genome_covered_fraction": round(cov_frac, 4),
            "genome_aligned_read_count": reads_cnt,
            "genome_anir_percent": round(anir, 2)
        }

    return genome_metrics, round(mapped_pct, 2), round(unmapped_pct, 2)


def parse_coverm_contig_tsv(contig_tsv):
    """Parses Track B coverm contig standard output."""
    raw_data = parse_coverm_tsv(contig_tsv)
    contig_metrics = {}
    for cid, row in raw_data.items():
        mean_d = 0.0
        trim_d = 0.0
        cov_frac = 0.0
        var_d = 0.0
        anir = 0.0

        for col, val in row.items():
            cl = col.lower()
            if "trimmed" in cl:
                trim_d = clean_val(val, 0.0)
            elif "mean" in cl and "trimmed" not in cl:
                mean_d = clean_val(val, 0.0)
            elif "covered" in cl or "fraction" in cl:
                cov_frac = clean_val(val, 0.0)
            elif "variance" in cl:
                var_d = clean_val(val, 0.0)
            elif "anir" in cl:
                anir_raw = clean_val(val, 0.0)
                anir = anir_raw * 100.0 if anir_raw <= 1.0 else anir_raw

        contig_metrics[cid] = {
            "contig_mean_depth": round(mean_d, 2),
            "contig_trimmed_mean_depth": round(trim_d, 2),
            "contig_covered_fraction": round(cov_frac, 4),
            "contig_coverage_variance": round(var_d, 2),
            "contig_anir_percent": round(anir, 2)
        }
    return contig_metrics


def parse_coverm_strict_tsv(strict_tsv):
    """Parses Track B coverm contig strict (98/90) output."""
    raw_data = parse_coverm_tsv(strict_tsv)
    strict_metrics = {}
    for cid, row in raw_data.items():
        trim_d = 0.0
        cov_frac = 0.0
        for col, val in row.items():
            cl = col.lower()
            if "trimmed" in cl:
                trim_d = clean_val(val, 0.0)
            elif "covered" in cl or "fraction" in cl:
                cov_frac = clean_val(val, 0.0)

        strict_metrics[cid] = {
            "strict_trimmed_mean_depth": round(trim_d, 2),
            "strict_covered_fraction": round(cov_frac, 4)
        }
    return strict_metrics


def parse_junction_qc_tsv(junc_tsv):
    """Parses Track C rotated_junction_qc.tsv."""
    data = {}
    if not junc_tsv or not os.path.exists(junc_tsv):
        return data

    with open(junc_tsv) as f:
        header = [h.strip() for h in f.readline().split("\t")]
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split("\t")]
            row = dict(zip(header, parts))
            cid = row.get("contig") or row.get("contig_id") or row.get("Contig") or parts[0]
            data[cid] = {
                "length": int(clean_val(row.get("length") or row.get("contig_length"), 0)),
                "junction_median_depth_100bp": clean_val(row.get("junction_median_depth_100bp"), 0.0),
                "junction_median_depth_300bp": clean_val(row.get("junction_median_depth_300bp"), 0.0),
                "local_flank_median_depth": clean_val(row.get("local_flank_median_depth"), 0.0),
                "junction_local_ratio": clean_val(row.get("junction_local_ratio"), 0.0),
                "junction_zero_coverage_run_bp": int(clean_val(row.get("junction_zero_coverage_run_bp") or row.get("zero_coverage_run_bp"), 0)),
                "junction_spanning_reads_unique": int(clean_val(row.get("junction_spanning_reads_unique") or row.get("spanning_reads_unique"), 0)),
                "junction_bridging_pairs_unique": int(clean_val(row.get("junction_bridging_pairs_unique") or row.get("bridging_pairs_unique"), 0)),
            }
    return data


parse_rotated_junction_qc_tsv = parse_junction_qc_tsv


def find_dict_entry(data_dict, cid):
    """Find dictionary entry matching cid directly or by header prefix/pattern."""
    if not data_dict:
        return {}
    if cid in data_dict:
        return data_dict[cid]
    for k, v in data_dict.items():
        if k.startswith(f"{cid}_") or k.startswith(f"{cid}ctg_") or k.startswith(f"{cid}.") or k == cid:
            return v
    for k, v in data_dict.items():
        tok = k.split()[0].split("_")[0]
        if tok == cid:
            return v
    for k, v in data_dict.items():
        if re.search(rf"(^|[^0-9a-zA-Z]){re.escape(cid)}([^0-9a-zA-Z]|$)", k):
            return v
    return {}


def compile_summary(
    sample,
    genome_abundance_tsv,
    contig_coverage_tsv,
    strict_coverage_tsv,
    junction_qc_tsv,
    out_tsv,
    low_ra_threshold=0.1,
    low_depth_threshold=5.0,
    low_breadth_threshold=0.95,
    competitive_warning_threshold=0.50,
    identity_sensitive_threshold=0.50,
    junction_warning_threshold=0.70,
    junction_fail_threshold=0.30
):
    """Compiles Track A, Track B, Track C into final summary table."""
    t_a_data, mapped_pct, unmapped_pct = parse_coverm_genome_tsv(genome_abundance_tsv)
    t_b_data = parse_coverm_contig_tsv(contig_coverage_tsv)
    t_b_strict = parse_coverm_strict_tsv(strict_coverage_tsv)
    t_c_data = parse_junction_qc_tsv(junction_qc_tsv)

    # Collect circular contig IDs from Track A (circular catalogue) and/or Track C (rotated junctions)
    all_cids = sorted(set(list(t_a_data.keys()) + list(t_c_data.keys())))
    if not all_cids:
        all_cids = sorted(set(t_b_data.keys()))

    rows = []
    for cid in all_cids:
        ga = find_dict_entry(t_a_data, cid)
        cb = find_dict_entry(t_b_data, cid)
        cbs = find_dict_entry(t_b_strict, cid)
        jc = find_dict_entry(t_c_data, cid)

        # Length
        contig_len = jc.get("length") or cb.get("length") or ga.get("length") or 0
        if not contig_len:
            m_len = re.search(r"len-(\d+)", cid)
            if m_len:
                contig_len = int(m_len.group(1))

        # Track A metrics
        ra_pct = ga.get("relative_abundance_percent", 0.0)
        g_mean = ga.get("genome_mean_depth", 0.0)
        g_trimmed = ga.get("genome_trimmed_mean_depth", 0.0)
        g_cov_frac = ga.get("genome_covered_fraction", 0.0)
        g_reads = ga.get("genome_aligned_read_count", 0)
        g_anir = ga.get("genome_anir_percent", 0.0)
        trimmed_to_mean = round(g_trimmed / max(0.001, g_mean), 3) if g_mean > 0 else 1.0

        # Track B metrics
        c_trimmed = cb.get("contig_trimmed_mean_depth", 0.0)
        c_cov_frac = cb.get("contig_covered_fraction", 0.0)
        c_anir = cb.get("contig_anir_percent", 0.0)

        # Track B strict
        s_trimmed = cbs.get("strict_trimmed_mean_depth", 0.0)
        s_cov_frac = cbs.get("strict_covered_fraction", 0.0)

        # Diagnostic ratios
        # C = Track B contig trimmed / Track A genome trimmed
        comp_ratio = round(c_trimmed / max(0.001, g_trimmed), 3) if g_trimmed > 0 else 1.0
        # S = strict trimmed / standard contig trimmed
        strict_ratio = round(s_trimmed / max(0.001, c_trimmed), 3) if c_trimmed > 0 else 1.0

        # Track C junction metrics
        has_junction_data = bool(jc) and (jc.get("junction_median_depth_100bp") is not None)
        if has_junction_data:
            j_depth_100 = jc.get("junction_median_depth_100bp", 0.0)
            d_local = jc.get("local_flank_median_depth", 0.0)
            j_local_r = jc.get("junction_local_ratio", 0.0)
            j_zero_run = jc.get("junction_zero_coverage_run_bp", 0)
            j_span_uniq = jc.get("junction_spanning_reads_unique", 0)
            j_bridge_uniq = jc.get("junction_bridging_pairs_unique", 0)

            # Boolean Flags
            flag_junc_depth_drop = (j_local_r < junction_warning_threshold)
            flag_junc_gap = (j_zero_run > 50)
            flag_no_span = (j_span_uniq == 0 and j_bridge_uniq == 0)
        else:
            j_depth_100 = "NA"
            d_local = "NA"
            j_local_r = "NA"
            j_zero_run = "NA"
            j_span_uniq = "NA"
            j_bridge_uniq = "NA"

            flag_junc_depth_drop = False
            flag_junc_gap = False
            flag_no_span = False

        # Boolean Flags
        flag_low_ra = (ra_pct < low_ra_threshold)
        flag_low_depth = (g_trimmed < low_depth_threshold)
        flag_low_breadth = (g_cov_frac < low_breadth_threshold)
        flag_possible_recruitment = (comp_ratio < competitive_warning_threshold)
        flag_identity_sensitive = (strict_ratio < identity_sensitive_threshold)

        # Overall Status Arbitration
        if g_trimmed < low_depth_threshold:
            status = "LOW_INFORMATION"
        elif has_junction_data and (flag_junc_gap or (j_local_r < junction_fail_threshold and flag_no_span)):
            status = "JUNCTION_SUSPECT"
        elif g_cov_frac < 0.85:
            status = "GLOBAL_MAPPING_SUSPECT"
        elif flag_possible_recruitment:
            status = "RECRUITMENT_CAUTION"
        elif flag_identity_sensitive:
            status = "IDENTITY_SENSITIVE"
        elif not has_junction_data:
            status = "UNVERIFIED_NO_ROTATED_DATA"
        elif not any([flag_junc_depth_drop, flag_junc_gap, flag_no_span, flag_low_depth, flag_low_breadth]):
            if g_trimmed >= 10.0 and j_local_r >= 0.70 and j_span_uniq >= 3 and j_bridge_uniq >= 5:
                status = "PASS_HIGH_CONFIDENCE"
            else:
                status = "PASS_MODERATE"
        else:
            status = "PASS_MODERATE"

        rows.append({
            "sample": sample,
            "circular_genome": cid,
            "length": contig_len,
            # Track A
            "relative_abundance_percent": ra_pct,
            "catalogue_mapped_percent": mapped_pct,
            "catalogue_unmapped_percent": unmapped_pct,
            "genome_mean_depth": g_mean,
            "genome_trimmed_mean_depth": g_trimmed,
            "genome_trimmed_to_mean_ratio": trimmed_to_mean,
            "genome_covered_fraction": g_cov_frac,
            "genome_aligned_read_count": g_reads,
            "genome_anir_percent": g_anir,
            # Track B
            "contig_trimmed_mean_depth": c_trimmed,
            "contig_covered_fraction": c_cov_frac,
            "contig_anir_percent": c_anir,
            "competitive_depth_ratio": comp_ratio,
            # Track B strict
            "strict_trimmed_mean_depth": s_trimmed,
            "strict_covered_fraction": s_cov_frac,
            "strict_depth_ratio": strict_ratio,
            # Track C
            "junction_median_depth_100bp": j_depth_100,
            "local_flank_median_depth": d_local,
            "junction_local_ratio": j_local_r,
            "junction_zero_coverage_run_bp": j_zero_run,
            "junction_spanning_reads_unique": j_span_uniq,
            "junction_bridging_pairs_unique": j_bridge_uniq,
            # Flags
            "flag_low_relative_abundance": "TRUE" if flag_low_ra else "FALSE",
            "flag_low_depth": "TRUE" if flag_low_depth else "FALSE",
            "flag_low_breadth": "TRUE" if flag_low_breadth else "FALSE",
            "flag_possible_homologous_recruitment": "TRUE" if flag_possible_recruitment else "FALSE",
            "flag_identity_sensitive": "TRUE" if flag_identity_sensitive else "FALSE",
            "flag_junction_depth_drop": "TRUE" if flag_junc_depth_drop else "FALSE",
            "flag_junction_gap": "TRUE" if flag_junc_gap else "FALSE",
            "flag_no_junction_span": "TRUE" if flag_no_span else "FALSE",
            # Status
            "overall_status": status
        })

    fields = [
        "sample", "circular_genome", "length",
        "relative_abundance_percent", "catalogue_mapped_percent", "catalogue_unmapped_percent",
        "genome_mean_depth", "genome_trimmed_mean_depth", "genome_trimmed_to_mean_ratio",
        "genome_covered_fraction", "genome_aligned_read_count", "genome_anir_percent",
        "contig_trimmed_mean_depth", "contig_covered_fraction", "contig_anir_percent",
        "competitive_depth_ratio",
        "strict_trimmed_mean_depth", "strict_covered_fraction", "strict_depth_ratio",
        "junction_median_depth_100bp", "local_flank_median_depth", "junction_local_ratio",
        "junction_zero_coverage_run_bp",
        "junction_spanning_reads_unique", "junction_bridging_pairs_unique",
        "flag_low_relative_abundance", "flag_low_depth", "flag_low_breadth",
        "flag_possible_homologous_recruitment", "flag_identity_sensitive",
        "flag_junction_depth_drop", "flag_junction_gap", "flag_no_junction_span",
        "overall_status"
    ]

    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    with open(out_tsv, "w") as out_f:
        out_f.write("\t".join(fields) + "\n")
        for r in rows:
            out_f.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")

    print(f"[compile_qc] Written {len(rows)} circular genomes summary to {out_tsv}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Compile circular Illumina abundance and junction QC summary.")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--genome-abundance-tsv", required=True, help="Track A coverm genome output")
    parser.add_argument("--contig-coverage-tsv", required=True, help="Track B coverm contig output")
    parser.add_argument("--strict-coverage-tsv", required=True, help="Track B coverm contig strict output")
    parser.add_argument("--junction-qc-tsv", required=True, help="Track C rotated junction QC output")
    parser.add_argument("--out-tsv", required=True, help="Output circular_illumina_qc_summary.tsv")

    # Thresholds
    parser.add_argument("--low-ra-threshold", type=float, default=0.1)
    parser.add_argument("--low-depth-threshold", type=float, default=5.0)
    parser.add_argument("--low-breadth-threshold", type=float, default=0.95)
    parser.add_argument("--competitive-warning-threshold", type=float, default=0.50)
    parser.add_argument("--identity-sensitive-threshold", type=float, default=0.50)
    parser.add_argument("--junction-warning-threshold", type=float, default=0.70)
    parser.add_argument("--junction-fail-threshold", type=float, default=0.30)

    args = parser.parse_args()
    compile_summary(
        args.sample,
        args.genome_abundance_tsv,
        args.contig_coverage_tsv,
        args.strict_coverage_tsv,
        args.junction_qc_tsv,
        args.out_tsv,
        low_ra_threshold=args.low_ra_threshold,
        low_depth_threshold=args.low_depth_threshold,
        low_breadth_threshold=args.low_breadth_threshold,
        competitive_warning_threshold=args.competitive_warning_threshold,
        identity_sensitive_threshold=args.identity_sensitive_threshold,
        junction_warning_threshold=args.junction_warning_threshold,
        junction_fail_threshold=args.junction_fail_threshold
    )


if __name__ == "__main__":
    main()
