#!/usr/bin/env python3
"""
assign_cog_footprints.py
High-throughput domain/footprint-aware COG assignment engine.

Key algorithms:
1. PSSM/CDD translation and strict integer length validation (slen_blast == PSSM_Length_cddid).
2. Collinear HSP-to-Footprint aggregation for multi-HSP alignments on the same query x PSSM.
3. Footprint-first greedy interval selection for highest-scoring non-overlapping footprints.
4. Three orthogonal classification dimensions:
   - Dimension 1: COG Match Quality (Match code 0-3 based on COG2024 75% rule on footprints).
   - Dimension 2: Region Evidence (unique, competing, ambiguous with continuous score margins).
   - Dimension 3: Protein Architecture (single_cog, multi_cog, repeated_cog, mixed).
5. Decoupled Truncation Evidence Vector (0-3 score).
"""

import argparse
import gzip
import os
import re
import sys

class DatabaseSynchronizationError(Exception):
    """Raised when BLAST slen does not match cddid.tbl PSSM-Length."""
    pass

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def load_pssm_mapping(mapping_path):
    """
    Load cdd_pssm_mapping.tsv.gz created by Gate 0.
    Columns: pssm_id, cog_id, pssm_length, cdd_name, cdd_desc, cog_category, cog_name.
    """
    pssm_map = {}
    with open_file(mapping_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            pssm_id = parts[col_idx["pssm_id"]].strip()
            cog_id = parts[col_idx["cog_id"]].strip()
            length = int(parts[col_idx["pssm_length"]].strip()) if "pssm_length" in col_idx else 0
            category = parts[col_idx["cog_category"]].strip() if "cog_category" in col_idx else "S"
            name = parts[col_idx["cog_name"]].strip() if "cog_name" in col_idx else ""
            spec_thresh = parts[col_idx["specific_bitscore_threshold"]].strip() if "specific_bitscore_threshold" in col_idx else "NA"
            pssm_map[pssm_id] = {
                "cog_id": cog_id,
                "length": length,
                "category": category,
                "name": name,
                "specific_bitscore_threshold": spec_thresh
            }
            # Also register cog_id directly in case BLAST output contains COG accession as sseqid
            if cog_id not in pssm_map:
                pssm_map[cog_id] = pssm_map[pssm_id]
    return pssm_map

def extract_pssm_id(sseqid):
    """
    Extract clean PSSM identifier or COG ID from BLAST sseqid string.
    Examples: 'gnl|CDD|223080' -> '223080', 'CDD:223080' -> '223080', 'COG0001' -> 'COG0001'.
    """
    if "|" in sseqid:
        parts = sseqid.split("|")
        return parts[-1].strip()
    if ":" in sseqid:
        return sseqid.split(":")[-1].strip()
    return sseqid.strip()

def aggregate_collinear_hsps(hsps):
    """
    Aggregate collinear HSPs belonging to the same (query_id, pssm_id).
    hsps: list of raw hit dicts for the same query x pssm, sorted by qstart.
    Returns: list of aggregated candidate footprint dicts.
    """
    if not hsps:
        return []
    if len(hsps) == 1:
        h = dict(hsps[0])
        h["n_hsps"] = 1
        return [h]

    # Sort strictly by qstart
    sorted_hsps = sorted(hsps, key=lambda x: (x["qstart"], x["sstart"]))
    footprints = []
    
    current_fp = dict(sorted_hsps[0])
    current_fp["n_hsps"] = 1
    current_fp["constituent_hsps"] = [(sorted_hsps[0]["qstart"], sorted_hsps[0]["qend"])]

    for next_h in sorted_hsps[1:]:
        # Collinearity criteria:
        # 1. Monotonic on both query and subject
        same_orientation = (next_h["sstart"] >= current_fp["sstart"]) and (next_h["send"] >= current_fp["send"])
        # 2. Minimal query overlap (<= 15%)
        q_overlap = max(0, current_fp["qend"] - next_h["qstart"] + 1)
        q_overlap_ratio = q_overlap / max(1, min(current_fp["qend"] - current_fp["qstart"] + 1, next_h["qend"] - next_h["qstart"] + 1))
        # 3. Minimal subject overlap (<= 15%)
        s_overlap = max(0, current_fp["send"] - next_h["sstart"] + 1)
        s_overlap_ratio = s_overlap / max(1, min(current_fp["send"] - current_fp["sstart"] + 1, next_h["send"] - next_h["sstart"] + 1))
        # 4. Inter-HSP distance <= 150 aa
        inter_gap = max(0, next_h["qstart"] - current_fp["qend"] - 1)

        if same_orientation and q_overlap_ratio <= 0.15 and s_overlap_ratio <= 0.15 and inter_gap <= 150:
            # Merge into current footprint
            current_fp["qstart"] = min(current_fp["qstart"], next_h["qstart"])
            current_fp["qend"] = max(current_fp["qend"], next_h["qend"])
            current_fp["sstart"] = min(current_fp["sstart"], next_h["sstart"])
            current_fp["send"] = max(current_fp["send"], next_h["send"])
            # Length = total span
            current_fp["length"] = current_fp["qend"] - current_fp["qstart"] + 1
            # Additive bitscore with small transition penalty
            current_fp["bitscore"] = round(current_fp["bitscore"] + next_h["bitscore"] - 2.0, 1)
            # Minimum evalue
            current_fp["evalue"] = min(current_fp["evalue"], next_h["evalue"])
            current_fp["n_hsps"] += 1
            current_fp["constituent_hsps"].append((next_h["qstart"], next_h["qend"]))
        else:
            footprints.append(current_fp)
            current_fp = dict(next_h)
            current_fp["n_hsps"] = 1
            current_fp["constituent_hsps"] = [(next_h["qstart"], next_h["qend"])]

    footprints.append(current_fp)
    return footprints

def calculate_overlap_ratio(fp1, fp2):
    """Compute spatial overlap ratio: overlap / min(len1, len2)."""
    overlap = max(0, min(fp1["qend"], fp2["qend"]) - max(fp1["qstart"], fp2["qstart"]) + 1)
    len1 = fp1["qend"] - fp1["qstart"] + 1
    len2 = fp2["qend"] - fp2["qstart"] + 1
    return overlap / max(1, min(len1, len2))

def evaluate_75_percent_match_code(profile_cov, query_cov):
    """
    Evaluate COG2024 official match code (0-3) based on 75% coverage rule.
    Returns: (match_code, match_class, profile_cov_class, query_cov_class)
    """
    prof_ge_75 = (profile_cov >= 0.75)
    query_ge_75 = (query_cov >= 0.75)

    if prof_ge_75 and query_ge_75:
        return 0, "full_match", "profile_complete", "query_complete"
    elif prof_ge_75 and not query_ge_75:
        return 1, "profile_complete_query_partial", "profile_complete", "query_partial"
    elif not prof_ge_75 and query_ge_75:
        return 2, "query_complete_profile_partial", "profile_partial", "query_complete"
    else:
        return 3, "partial_both", "profile_partial", "query_partial"

def process_query_hits(qseqid, q_hits, pssm_map, evalue_cutoff, min_aln_len, min_prof_cov,
                       ambiguity_ratio_cutoff, ambiguity_delta_cutoff, gff_metadata=None):
    """
    Process all raw hits for a single query protein.
    Returns list of assigned footprint records (dict).
    """
    if not q_hits:
        return []

    # 1. Group by pssm_id and aggregate collinear HSPs
    by_pssm = {}
    for h in q_hits:
        pssm_id = h["pssm_id"]
        by_pssm.setdefault(pssm_id, []).append(h)

    candidate_footprints = []
    for pssm_id, hsp_list in by_pssm.items():
        aggr_fps = aggregate_collinear_hsps(hsp_list)
        for fp in aggr_fps:
            slen = fp["slen"]
            qlen = fp["qlen"]
            aln_len = fp["qend"] - fp["qstart"] + 1
            query_cov = round(aln_len / max(1, qlen), 4)
            s_span = abs(fp["send"] - fp["sstart"]) + 1
            profile_cov = round(s_span / max(1, slen), 4)

            # Pre-filter candidate footprints
            if fp["evalue"] > evalue_cutoff:
                continue
            if aln_len < min_aln_len and profile_cov < min_prof_cov:
                continue

            fp["aln_len"] = aln_len
            fp["query_cov"] = query_cov
            fp["profile_cov"] = profile_cov
            candidate_footprints.append(fp)

    if not candidate_footprints:
        return []

    # 2. Footprint-First Selection: Greedy interval selection by bitscore desc, evalue asc
    candidate_footprints.sort(key=lambda x: (-x["bitscore"], x["evalue"], -(x["qend"] - x["qstart"])))

    retained_footprints = []
    for cand in candidate_footprints:
        is_overlapping = False
        for ret in retained_footprints:
            overlap_ratio = calculate_overlap_ratio(cand, ret)
            if overlap_ratio >= 0.50:
                is_overlapping = True
                break
        if not is_overlapping:
            retained_footprints.append(cand)

    # Sort retained footprints by qstart (linear protein order)
    retained_footprints.sort(key=lambda x: x["qstart"])

    # 3. Protein Architecture Dimension
    n_total_fp = len(retained_footprints)
    distinct_cogs = {fp["cog_id"] for fp in retained_footprints}
    n_distinct_cog = len(distinct_cogs)

    if n_total_fp == 1:
        architecture_class = "single_cog"
    elif n_distinct_cog == n_total_fp:
        architecture_class = "multi_cog"
    elif n_distinct_cog == 1:
        architecture_class = "repeated_cog"
    else:
        architecture_class = "mixed_multi_repeated"

    is_multi_cog = (n_distinct_cog > 1)
    is_repeated_cog = (n_total_fp > n_distinct_cog)
    cog_arch_str = ";".join([f"{fp['cog_id']}[{fp['qstart']}-{fp['qend']}]" for fp in retained_footprints])

    # 4. Region-Level Evidence & Match Quality per retained footprint
    assigned_records = []
    has_competing = False
    has_ambiguous = False

    for rank, fp in enumerate(retained_footprints, 1):
        region_id = f"{fp['qstart']}-{fp['qend']}"
        m_code, m_class, p_cov_cls, q_cov_cls = evaluate_75_percent_match_code(fp["profile_cov"], fp["query_cov"])

        # Find all candidate footprints competing for this region (overlap >= 0.50)
        competitors = []
        for cand in candidate_footprints:
            if cand is fp:
                continue
            ov = calculate_overlap_ratio(fp, cand)
            if ov >= 0.50:
                competitors.append(cand)

        if not competitors:
            region_status = "unique"
            top2_cog = "NA"
            delta_bitscore = "NA"
            bitscore_ratio = "NA"
            delta_bitscore_fraction = "NA"
        else:
            has_competing = True
            # Sort competitors by bitscore desc
            competitors.sort(key=lambda x: -x["bitscore"])
            top2 = competitors[0]
            top2_cog = top2["cog_id"]
            d_bit = round(fp["bitscore"] - top2["bitscore"], 1)
            b_ratio = round(top2["bitscore"] / max(0.1, fp["bitscore"]), 4)
            d_frac = round(d_bit / max(0.1, fp["bitscore"]), 4)

            delta_bitscore = str(d_bit)
            bitscore_ratio = str(b_ratio)
            delta_bitscore_fraction = str(d_frac)

            if d_frac < ambiguity_ratio_cutoff or d_bit < ambiguity_delta_cutoff:
                region_status = "ambiguous"
                has_ambiguous = True
            else:
                region_status = "competing"

        # Truncation evidence vector
        truncation_score = 0
        prodigal_partial = False
        near_contig_edge = False
        short_vs_cog = (fp["query_cov"] < 0.60)
        if short_vs_cog:
            truncation_score += 1

        if gff_metadata and qseqid in gff_metadata:
            meta = gff_metadata[qseqid]
            if meta.get("partial", "00") != "00":
                prodigal_partial = True
                truncation_score += 1
            if meta.get("contig_edge_dist", 999) < 50:
                near_contig_edge = True
                truncation_score += 1

        likely_truncated_gene = (truncation_score >= 2)

        p_info = pssm_map.get(fp["pssm_id"], {})
        spec_str = p_info.get("specific_bitscore_threshold", "NA")
        passes_spec = "NA"
        if spec_str != "NA":
            try:
                thresh_val = float(spec_str)
                passes_spec = str(fp["bitscore"] >= thresh_val)
            except ValueError:
                passes_spec = "NA"

        assigned_records.append({
            "protein_id": qseqid,
            "cog_id": fp["cog_id"],
            "pssm_id": fp["pssm_id"],
            "qstart": fp["qstart"],
            "qend": fp["qend"],
            "qlen": fp["qlen"],
            "sstart": fp["sstart"],
            "send": fp["send"],
            "slen": fp["slen"],
            "aln_len": fp["aln_len"],
            "query_cov": fp["query_cov"],
            "profile_cov": fp["profile_cov"],
            "evalue": fp["evalue"],
            "bitscore": fp["bitscore"],
            "specific_bitscore_threshold": spec_str,
            "passes_specific_threshold": passes_spec,
            "cog_match_code": m_code,
            "cog_match_class": m_class,
            "profile_coverage_class": p_cov_cls,
            "query_coverage_class": q_cov_cls,
            "footprint_rank": rank,
            "region_id": region_id,
            "region_status": region_status,
            "architecture_class": architecture_class,
            "has_competing_profiles": has_competing,
            "has_ambiguous_region": has_ambiguous,
            "is_multi_cog": is_multi_cog,
            "is_repeated_cog": is_repeated_cog,
            "top2_cog": top2_cog,
            "delta_bitscore": delta_bitscore,
            "bitscore_ratio": bitscore_ratio,
            "delta_bitscore_fraction": delta_bitscore_fraction,
            "likely_truncated_gene": likely_truncated_gene,
            "truncation_evidence_count": truncation_score,
            "prodigal_partial": prodigal_partial,
            "near_contig_edge": near_contig_edge,
            "short_vs_cog": short_vs_cog,
            "cog_footprint_architecture": cog_arch_str,
            "cog_category": fp.get("cog_category", "S"),
            "cog_name": fp.get("cog_name", "")
        })

    # Update has_competing_profiles and has_ambiguous_region across all records of this protein
    for rec in assigned_records:
        rec["has_competing_profiles"] = has_competing
        rec["has_ambiguous_region"] = has_ambiguous

    return assigned_records

def main():
    parser = argparse.ArgumentParser(
        description="Assign COG footprints from raw RPS-BLAST output with 3-dimensional classification."
    )
    parser.add_argument("-i", "--raw_hits", required=True, help="Raw RPS-BLAST TSV output (.tsv or .tsv.gz)")
    parser.add_argument("-m", "--pssm_mapping", required=True, help="cdd_pssm_mapping.tsv.gz from Gate 0")
    parser.add_argument("-o", "--out_tsv", required=True, help="Output assigned footprints TSV (.tsv or .tsv.gz)")
    parser.add_argument("--gff_meta", default="", help="Optional protein metadata TSV with prodigal partial flags")
    parser.add_argument("--evalue", type=float, default=1e-4, help="E-value cutoff (default: 1e-4)")
    parser.add_argument("--min_aln_len", type=int, default=40, help="Minimum alignment length (default: 40)")
    parser.add_argument("--min_prof_cov", type=float, default=0.30, help="Minimum profile coverage (default: 0.30)")
    parser.add_argument("--ambiguity_ratio_cutoff", type=float, default=0.15, help="R_delta ratio cutoff (default: 0.15)")
    parser.add_argument("--ambiguity_delta_cutoff", type=float, default=20.0, help="Delta bitscore cutoff (default: 20.0)")

    args = parser.parse_args()

    # Load mapping
    pssm_map = load_pssm_mapping(args.pssm_mapping)

    # Load optional GFF metadata
    gff_metadata = {}
    if args.gff_meta and os.path.isfile(args.gff_meta):
        with open_file(args.gff_meta, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    pid = parts[0].strip()
                    gff_metadata[pid] = {
                        "partial": parts[1].strip() if len(parts) > 1 else "00",
                        "contig_edge_dist": int(parts[2].strip()) if len(parts) > 2 and parts[2].strip().isdigit() else 999
                    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)

    # Stream raw RPS-BLAST hits
    # Format: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
    out_opener = gzip.open if args.out_tsv.endswith(".gz") else open
    with open_file(args.raw_hits, "rt") as in_f, out_opener(args.out_tsv, "wt", encoding="utf-8") as out_f:
        out_fields = [
            "protein_id", "cog_id", "pssm_id", "cog_category", "cog_name",
            "qstart", "qend", "qlen", "sstart", "send", "slen", "aln_len",
            "query_cov", "profile_cov", "evalue", "bitscore",
            "specific_bitscore_threshold", "passes_specific_threshold",
            "cog_match_code", "cog_match_class", "profile_coverage_class", "query_coverage_class",
            "footprint_rank", "region_id", "region_status", "architecture_class",
            "has_competing_profiles", "has_ambiguous_region", "is_multi_cog", "is_repeated_cog",
            "top2_cog", "delta_bitscore", "bitscore_ratio", "delta_bitscore_fraction",
            "likely_truncated_gene", "truncation_evidence_count", "prodigal_partial", "near_contig_edge", "short_vs_cog",
            "cog_footprint_architecture"
        ]
        out_f.write("\t".join(out_fields) + "\n")

        current_qid = None
        current_hits = []
        sync_checked = set()

        for line_idx, line in enumerate(in_f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue

            qseqid = parts[0].strip()
            sseqid = parts[1].strip()
            pident = float(parts[2].strip())
            length = int(parts[3].strip())
            qlen = int(parts[4].strip())
            slen = int(parts[5].strip())
            qstart = int(parts[6].strip())
            qend = int(parts[7].strip())
            sstart = int(parts[8].strip())
            send = int(parts[9].strip())
            evalue = float(parts[10].strip())
            bitscore = float(parts[11].strip())

            pssm_id = extract_pssm_id(sseqid)
            p_info = pssm_map.get(pssm_id, {})
            cog_id = p_info.get("cog_id", pssm_id)
            expected_slen = p_info.get("length", 0)

            # Strict integer length validation against cddid.tbl
            if expected_slen > 0 and pssm_id not in sync_checked:
                if slen != expected_slen:
                    raise DatabaseSynchronizationError(
                        f"Database Sync Error at line {line_idx}: PSSM {pssm_id} ({cog_id}) "
                        f"has RPS-BLAST slen={slen} but cddid.tbl PSSM_Length={expected_slen}! "
                        "The Cog_LE database and cddid.tbl come from inconsistent builds."
                    )
                sync_checked.add(pssm_id)

            hit_dict = {
                "qseqid": qseqid,
                "pssm_id": pssm_id,
                "cog_id": cog_id,
                "cog_category": p_info.get("category", "S"),
                "cog_name": p_info.get("name", ""),
                "qlen": qlen,
                "slen": slen,
                "qstart": qstart,
                "qend": qend,
                "sstart": sstart,
                "send": send,
                "evalue": evalue,
                "bitscore": bitscore
            }

            if qseqid != current_qid:
                if current_hits:
                    records = process_query_hits(
                        current_qid, current_hits, pssm_map,
                        args.evalue, args.min_aln_len, args.min_prof_cov,
                        args.ambiguity_ratio_cutoff, args.ambiguity_delta_cutoff, gff_metadata
                    )
                    for r in records:
                        out_f.write("\t".join(str(r[k]) for k in out_fields) + "\n")
                current_qid = qseqid
                current_hits = [hit_dict]
            else:
                current_hits.append(hit_dict)

        if current_hits:
            records = process_query_hits(
                current_qid, current_hits, pssm_map,
                args.evalue, args.min_aln_len, args.min_prof_cov,
                args.ambiguity_ratio_cutoff, args.ambiguity_delta_cutoff, gff_metadata
            )
            for r in records:
                out_f.write("\t".join(str(r[k]) for k in out_fields) + "\n")

    print(f"SUCCESS: Assigned COG footprints written to {args.out_tsv}")

if __name__ == "__main__":
    main()
