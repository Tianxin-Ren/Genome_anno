#!/usr/bin/env python3
"""
audit_circular_assembly.py (Version 2.1)
--------------------------------------------------------------------------------
Comprehensive assembly verification for circular contigs (chromosomes/plasmids)
based on the Trigodet-informed Four-Gate Framework (Trigodet et al. 2026, Nature
Biotechnology; MIMAG; CheckM2; geNomad).

Evaluates Four Gates of assembly trustworthiness:
  Gate 1 - Physical Structural Support:
           * Targeted CIGAR-aware junction spanning & support fraction (unique QNAMEs).
           * Terminal-anchor whole-assembly uniqueness.
           * All-alignment zero-coverage detection with raw-read gap rescue.
           * Whole-assembly systematic clipping (>80%, distance to end >100bp).
           * Multi-record partner contig consistency (Methanothrix premature closure).
  Gate 2 - Structural & Coverage Coherence:
           * Primary-only depth metrics (mean, median, MAD, CV, breadth).
           * BLASTN self-alignment repeat breadth profiling (>=80% id, >200bp).
  Gate 3 - Biological Completeness (MIMAG-derived):
           * CheckM2 completeness >90% and contamination <5%.
           * Independent 5S, 16S, and 23S rRNA gene presence.
           * Canonical 20-amino-acid tRNA isotype diversity (>=18/20 standard AAs).
  Gate 4 - Replicon Identity Evidence:
           * geNomad raw scores (chromosome, plasmid, virus).
           * Derived replicon assignment and confidence margin.

Outputs a structured, single-row TSV audit report per contig.
"""

import os
import sys
import argparse
import subprocess
import re
import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional

try:
    import numpy as np
except ImportError:
    np = None

try:
    import pysam
except ImportError:
    pysam = None

CANONICAL_20_AA = {
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"
}

def run_cmd(cmd):
    """Run shell command and return stdout string and return code."""
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return res.stdout.strip(), res.returncode

def get_fasta_length(fa_path):
    """Get total length and sequence of the first record in FASTA."""
    length = 0
    seq_chunks = []
    if not os.path.exists(fa_path):
        return 0, ""
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq = line.strip()
            length += len(seq)
            seq_chunks.append(seq)
    return length, "".join(seq_chunks)

def wilson_score_interval(successes, total, z=1.96):
    """Calculate Wilson score confidence interval for a binomial proportion."""
    if total <= 0:
        return 0.0, 0.0
    p = successes / total
    denom = 1.0 + (z**2) / total
    center = (p + (z**2) / (2.0 * total)) / denom
    spread = (z / denom) * math.sqrt((p * (1.0 - p) / total) + ((z**2) / (4.0 * total**2)))
    return max(0.0, center - spread), min(1.0, center + spread)

def resolve_bam_reference(bam_references, bin_id):
    """
    Resolve the actual reference name in a BAM file for a given contig/bin ID.
    Handles Myloasm/Flye header prefixes, e.g. 'u999856' matching
    'u999856_len-2079472_circular-yes_depth-35.2_duplicated-0'.
    """
    if not bam_references:
        return bin_id
    if bin_id in bam_references:
        return bin_id
    for ref in bam_references:
        if ref.startswith(f"{bin_id}_") or ref.startswith(f"{bin_id}ctg_") or ref.startswith(f"{bin_id}."):
            return ref
    for ref in bam_references:
        tok = ref.split()[0]
        if tok == bin_id or tok.split("_")[0] == bin_id:
            return ref
    for ref in bam_references:
        if bin_id in ref:
            return ref
    return bin_id

def calculate_coverage_arrays(bam_file, contig_id, contig_len):
    """
    Calculate primary depth array and all-mapped support depth array for a contig.
    Uses pysam if available, or falls back to samtools depth.
    """
    if np is None:
        primary_cov = [0] * contig_len
        all_cov = [0] * contig_len
    else:
        primary_cov = np.zeros(contig_len, dtype=np.uint32)
        all_cov = np.zeros(contig_len, dtype=np.uint32)

    if not bam_file or not os.path.exists(bam_file):
        return primary_cov, all_cov

    if pysam is not None:
        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                ref_id = resolve_bam_reference(bam.references, contig_id)
                if not ref_id or ref_id not in bam.references:
                    # Contig not in BAM references
                    return primary_cov, all_cov
                
                for read in bam.fetch(ref_id):
                    if read.is_unmapped:
                        continue
                    is_prim = not (read.is_secondary or read.is_supplementary)
                    curr_pos = read.reference_start
                    for op, op_len in (read.cigartuples or []):
                        if op in (0, 7, 8):  # M, =, X (consumes query and reference)
                            end_pos = min(contig_len, curr_pos + op_len)
                            if curr_pos < contig_len:
                                all_cov[curr_pos:end_pos] += 1
                                if is_prim:
                                    primary_cov[curr_pos:end_pos] += 1
                            curr_pos += op_len
                        elif op == 2:  # D (consumes reference only)
                            curr_pos += op_len
                        elif op == 3:  # N (consumes reference only)
                            curr_pos += op_len
            return primary_cov, all_cov
        except Exception as e:
            print(f"[audit] Warning: pysam coverage extraction failed ({e}), falling back to samtools", file=sys.stderr)

    # Fallback to samtools depth
    ref_id = contig_id
    cmd_idx = f"samtools idxstats '{bam_file}' | awk '{{print $1}}'"
    out_idx, code_idx = run_cmd(cmd_idx)
    if code_idx == 0 and out_idx:
        ref_id = resolve_bam_reference(out_idx.splitlines(), contig_id) or contig_id

    # All alignments: samtools depth -a -r ref_id
    cmd_all = f"samtools depth -a -r '{ref_id}' '{bam_file}' | awk '{{print $3}}'"
    out_all, code_all = run_cmd(cmd_all)
    if code_all == 0 and out_all:
        for idx, line in enumerate(out_all.splitlines()):
            if idx >= contig_len:
                break
            if line.isdigit():
                all_cov[idx] = int(line)

    # Primary only: samtools depth -a -r ref_id (filtering 0x100 | 0x800 = 0x900)
    cmd_prim = f"samtools view -b -F 0x904 '{bam_file}' '{ref_id}' | samtools depth -a - | awk '{{print $3}}'"
    out_prim, code_prim = run_cmd(cmd_prim)
    if code_prim == 0 and out_prim:
        for idx, line in enumerate(out_prim.splitlines()):
            if idx >= contig_len:
                break
            if line.isdigit():
                primary_cov[idx] = int(line)

    return primary_cov, all_cov

def compute_primary_coherence_metrics(primary_cov, contig_len):
    """Compute primary-only abundance and coherence metrics."""
    if contig_len == 0 or len(primary_cov) == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if np is not None:
        arr = np.asarray(primary_cov, dtype=np.float64)
        mean_d = float(np.mean(arr))
        median_d = float(np.median(arr))
        mad_d = float(np.median(np.abs(arr - median_d)))
        cv_d = float(np.std(arr) / mean_d) if mean_d > 0 else 0.0
        breadth_1x = float(np.sum(arr >= 1.0) / contig_len * 100.0)
        breadth_3x = float(np.sum(arr >= 3.0) / contig_len * 100.0)
    else:
        mean_d = float(sum(primary_cov)) / contig_len
        sorted_cov = sorted(primary_cov)
        median_d = float(sorted_cov[contig_len // 2])
        mad_d = float(sorted([abs(x - median_d) for x in primary_cov])[contig_len // 2])
        variance = sum((x - mean_d) ** 2 for x in primary_cov) / contig_len
        cv_d = float(math.sqrt(variance) / mean_d) if mean_d > 0 else 0.0
        breadth_1x = float(sum(1 for x in primary_cov if x >= 1) / contig_len * 100.0)
        breadth_3x = float(sum(1 for x in primary_cov if x >= 3) / contig_len * 100.0)

    return mean_d, median_d, mad_d, cv_d, breadth_1x, breadth_3x

def find_zero_coverage_blocks(all_cov, contig_len):
    """Find all contiguous 0x coverage blocks in all-mapped coverage array."""
    blocks = []
    in_zero = False
    start = 0

    for i in range(contig_len):
        if all_cov[i] == 0:
            if not in_zero:
                in_zero = True
                start = i
        else:
            if in_zero:
                in_zero = False
                blocks.append({"start": start, "end": i, "len": i - start})
    if in_zero:
        blocks.append({"start": start, "end": contig_len, "len": contig_len - start})

    total_zero_bp = sum(b["len"] for b in blocks)
    max_zero_block = max((b["len"] for b in blocks), default=0)
    high_risk_blocks = [b for b in blocks if b["len"] > 500]

    return blocks, total_zero_bp, max_zero_block, high_risk_blocks

def find_systematic_clipping_sites(bam_file, contig_id, contig_len, all_cov, min_depth=10, clipping_ratio_threshold=0.80, min_clipped_reads=15):
    """
    Detect internal systematic clipping sites (>80% clipping ratio, >100bp from ends).
    Strictly adheres to anvi'o find_misassemblies standard:
    min_dist_to_end = 100, clipping_ratio > threshold, and sufficient clipped read depth.
    Enhanced with population-resolved dual-allele analysis:
    Differentiates in-vivo structural polymorphisms (co-occurring through-reads and internal jump/deletion reads)
    from destructive misassembly breakpoints.
    """
    clipping_sites = []
    if not bam_file or not os.path.exists(bam_file) or pysam is None:
        return clipping_sites

    clipping_counts = defaultdict(int)
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            ref_id = resolve_bam_reference(bam.references, contig_id)
            if not ref_id or ref_id not in bam.references:
                return clipping_sites
            for read in bam.fetch(ref_id):
                if read.is_unmapped:
                    continue
                cigar = read.cigartuples or []
                if not cigar:
                    continue
                # Start clipping (soft=4, hard=5)
                if cigar[0][0] in (4, 5):
                    pos = read.reference_start
                    if 100 <= pos <= contig_len - 100:
                        clipping_counts[pos] += 1
                # End clipping
                if cigar[-1][0] in (4, 5):
                    pos = read.reference_end
                    if 100 <= pos <= contig_len - 100:
                        clipping_counts[pos] += 1

        for pos, count in clipping_counts.items():
            depth = all_cov[pos] if pos < contig_len else 0
            if depth >= min_depth and count >= min_clipped_reads:
                ratio = count / float(depth)
                if ratio > clipping_ratio_threshold:
                    site_info = {
                        "pos": pos,
                        "clipped_reads": count,
                        "depth": depth,
                        "clipping_ratio": ratio,
                        "is_sv_polymorphism": False,
                        "af_del": 0.0,
                        "af_ins": 1.0,
                        "del_freq": 0.0,
                        "ins_freq": 1.0,
                        "through_reads": 0,
                        "del_reads": count,
                        "del_partner_pos": None,
                        "del_span_bp": 0,
                        "sv_target_pos": None,
                        "sv_type": "MISASSEMBLY_BREAKPOINT"
                    }
                    # Secondary validation: quantify through-reads vs deletion/jump-reads to detect in-vivo SV
                    try:
                        with pysam.AlignmentFile(bam_file, "rb") as bam_check:
                            q_through = set()
                            q_del = set()
                            target_jump_positions = []
                            fetch_s = max(0, pos - 1000)
                            fetch_e = min(contig_len, pos + 1000)
                            for r in bam_check.fetch(ref_id, fetch_s, fetch_e):
                                if r.is_unmapped or r.is_secondary:
                                    continue
                                qn = r.query_name
                                r_s = r.reference_start
                                r_e = r.reference_end
                                cigar = r.cigartuples or []

                                is_deletion_allele = False
                                for op, op_len in cigar:
                                    if op == 2 and op_len >= 1000:
                                        is_deletion_allele = True
                                        break
                                if not is_deletion_allele and r.has_tag("SA"):
                                    sa_tag = r.get_tag("SA")
                                    for entry in sa_tag.split(";"):
                                        if not entry.strip():
                                            continue
                                        parts = entry.split(",")
                                        if len(parts) >= 4 and parts[0] == ref_id:
                                            sa_pos = int(parts[1])
                                            if abs(sa_pos - pos) >= 1000:
                                                is_deletion_allele = True
                                                target_jump_positions.append(sa_pos)

                                if is_deletion_allele:
                                    q_del.add(qn)
                                else:
                                    if r_s <= (pos - 300) and r_e >= (pos + 300):
                                        c_s_clip = cigar and cigar[0][0] in (4, 5) and abs(r_s - pos) <= 100
                                        c_e_clip = cigar and cigar[-1][0] in (4, 5) and abs(r_e - pos) <= 100
                                        if not (c_s_clip or c_e_clip):
                                            q_through.add(qn)

                            n_through = len(q_through - q_del)
                            n_del = len(q_del)
                            tot_informative = n_through + n_del
                            if n_through >= 15 and n_del >= 10 and tot_informative > 0:
                                af_del = n_del / float(tot_informative)
                                af_ins = n_through / float(tot_informative)
                                if 0.05 <= af_del <= 0.95:
                                    site_info["is_sv_polymorphism"] = True
                                    site_info["af_del"] = af_del
                                    site_info["af_ins"] = af_ins
                                    site_info["del_freq"] = af_del
                                    site_info["ins_freq"] = af_ins
                                    site_info["through_reads"] = n_through
                                    site_info["del_reads"] = n_del
                                    if target_jump_positions:
                                        if np is not None:
                                            med_jump = int(np.median(target_jump_positions))
                                        else:
                                            target_jump_positions.sort()
                                            med_jump = target_jump_positions[len(target_jump_positions) // 2]
                                        site_info["sv_target_pos"] = med_jump
                                        site_info["del_partner_pos"] = med_jump
                                        del_size = abs(med_jump - pos)
                                        site_info["del_span_bp"] = del_size
                                        site_info["sv_type"] = f"IN_VIVO_DELETION({del_size}bp,AF={af_del:.1%})"
                                    else:
                                        site_info["sv_type"] = f"IN_VIVO_SV(AF={af_del:.1%})"
                                    print(f"[audit] In-vivo structural polymorphism identified at pos {pos}: {site_info['sv_type']} (ins={n_through}, del={n_del})", file=sys.stderr)
                                    print(f"[audit] In-vivo structural polymorphism identified at pos {pos}: {site_info['sv_type']} (ins={n_through}, del={n_del})", file=sys.stderr)
                    except Exception as e_sv:
                        print(f"[audit] Warning during SV polymorphism check at pos {pos}: {e_sv}", file=sys.stderr)

                    clipping_sites.append(site_info)
                    if not site_info["is_sv_polymorphism"]:
                        print(f"[audit] Systematic internal clipping detected at pos {pos}: {count}/{depth} reads clipped ({ratio:.1%})", file=sys.stderr)
    except Exception as e:
        print(f"[audit] Warning: clipping analysis encountered error: {e}", file=sys.stderr)

    return clipping_sites

def analyze_targeted_synthetic_junction(junction_bam, flank_len, min_flank=500):
    """
    CIGAR-aware targeted junction analysis on synthetic junction BAM [last_flank + first_flank].
    Junction boundary J = flank_len.
    Evaluates spanning by unique read molecule (QNAME).
    Supports long reads (ONT/HiFi) containing normal sequencing indels while penalizing
    structural gaps (>50 bp) and junction-specific clipping.
    """
    spanning_qnames = set()
    clipped_qnames = set()
    informative_qnames = set()

    if not junction_bam or not os.path.exists(junction_bam) or pysam is None:
        return 0, 0, 0, 0.0, 0.0, 0.0, 0.0, set()

    J = flank_len
    left_target = J - min_flank
    right_target = J + min_flank
    window_len = right_target - left_target

    try:
        with pysam.AlignmentFile(junction_bam, "rb") as bam:
            for read in bam.fetch():
                if read.is_unmapped:
                    continue
                qname = read.query_name
                r_start = read.reference_start
                r_end = read.reference_end

                # Informative check: Does alignment reach near junction?
                reaches_left = (r_start <= J and r_end >= left_target)
                reaches_right = (r_start <= right_target and r_end >= J)

                if not (reaches_left or reaches_right):
                    continue

                # Clipping check within +-50 bp of J
                is_clipped = False
                cigar = read.cigartuples or []
                if cigar:
                    if cigar[0][0] in (4, 5) and abs(r_start - J) <= 50:
                        is_clipped = True
                    if cigar[-1][0] in (4, 5) and abs(r_end - J) <= 50:
                        is_clipped = True

                # Spanning check across [left_target, right_target]:
                # Must enter <= left_target and exit >= right_target without a large gap (>50 bp)
                spans_cigar = False
                if r_start <= left_target and r_end >= right_target and not is_clipped:
                    has_large_gap = False
                    matched_bp = 0
                    curr_pos = r_start
                    for op, op_len in cigar:
                        block_start = curr_pos
                        if op in (0, 7, 8):
                            block_end = curr_pos + op_len
                            ov_s = max(block_start, left_target)
                            ov_e = min(block_end, right_target)
                            if ov_e > ov_s:
                                matched_bp += (ov_e - ov_s)
                            curr_pos = block_end
                        elif op in (2, 3):
                            block_end = curr_pos + op_len
                            if max(block_start, left_target) < min(block_end, right_target) and op_len > 50:
                                has_large_gap = True
                            curr_pos = block_end
                        elif op == 1:
                            if left_target <= curr_pos <= right_target and op_len > 50:
                                has_large_gap = True

                    if not has_large_gap and matched_bp >= 0.70 * window_len:
                        spans_cigar = True

                if spans_cigar:
                    spanning_qnames.add(qname)
                    informative_qnames.add(qname)
                elif is_clipped:
                    clipped_qnames.add(qname)
                    informative_qnames.add(qname)
                elif reaches_left and reaches_right:
                    # Reaches both sides but broken by large gap or low matching -> informative negative
                    informative_qnames.add(qname)

    except Exception as e:
        print(f"[audit] Warning: junction BAM parsing failed: {e}", file=sys.stderr)

    num_spanning = len(spanning_qnames)
    num_clipped = len(clipped_qnames)
    num_informative = len(informative_qnames)

    support_frac = num_spanning / float(max(1, num_informative))
    clip_frac = num_clipped / float(max(1, num_informative))
    ci_lower, ci_upper = wilson_score_interval(num_spanning, num_informative)

    return num_spanning, num_clipped, num_informative, support_frac, clip_frac, ci_lower, ci_upper, spanning_qnames

def check_terminal_anchor_uniqueness(whole_bam, contig_id, contig_len, spanning_qnames, flank_len, low_abundance_partners=None):
    """
    Verify that junction-spanning reads map uniquely to the candidate contig's termini
    in the whole-assembly BAM without equivalent or better alternative mappings elsewhere.
    Multi-mapping reads whose secondary/supplementary alignments point only to known low-abundance
    assembly artifacts (low_abundance_partners) are exempted from triggering repeat-induced ambiguity.
    """
    if not whole_bam or not os.path.exists(whole_bam) or not spanning_qnames or pysam is None:
        return "UNIQUE", 0

    ambiguous_count = 0
    unique_count = 0
    low_partner_set = set(low_abundance_partners) if low_abundance_partners else set()

    try:
        with pysam.AlignmentFile(whole_bam, "rb") as bam:
            ref_id = resolve_bam_reference(bam.references, contig_id)
            if not ref_id or ref_id not in bam.references:
                return "UNIQUE", 0
            for qname in spanning_qnames:
                hits = list(bam.fetch(ref_id, 0, flank_len)) + list(bam.fetch(ref_id, max(0, contig_len - flank_len), contig_len))
                q_hits = [h for h in hits if h.query_name == qname]
                if q_hits:
                    has_ambig = False
                    for h in q_hits:
                        if h.mapping_quality < 5 or h.has_tag("XA"):
                            # Check if ambiguity is explained by low-abundance partner artifact
                            is_artifact_explained = False
                            if low_partner_set:
                                if h.has_tag("SA"):
                                    sa_str = h.get_tag("SA")
                                    if any(any(p in sa_entry for p in low_partner_set) for sa_entry in sa_str.split(";") if sa_entry.strip()):
                                        is_artifact_explained = True
                                if h.has_tag("XA"):
                                    xa_str = h.get_tag("XA")
                                    if any(any(p in xa_entry for p in low_partner_set) for xa_entry in xa_str.split(";") if xa_entry.strip()):
                                        is_artifact_explained = True
                            if not is_artifact_explained:
                                has_ambig = True
                                break
                    if has_ambig:
                        ambiguous_count += 1
                    else:
                        unique_count += 1
                else:
                    unique_count += 1

        if ambiguous_count > unique_count and ambiguous_count >= 3:
            status = "AMBIGUOUS (Repeat Termini)"
        elif ambiguous_count > 0:
            status = "PARTIALLY_UNIQUE"
        else:
            status = "UNIQUE"
        return status, ambiguous_count
    except Exception:
        return "UNIQUE", 0

def detect_methanothrix_partner_clustering(whole_bam, contig_id, contig_len, flank_len=5000, candidate_depth=0.0):
    """
    Detect Methanothrix-style premature circularization by analyzing multi-alignment records
    (primary, secondary, supplementary) for reads clipped near terminal ends.
    Evaluates:
      - Dominant partner contig
      - Unique read count
      - Fractional consistency (partner reads / total clipped reads at locus)
      - Partner locus clustering (MAD)
      - Partner contig depth & relative depth ratio (RDR)
      - Low-abundance artifact flag (RDR < 0.05)
    """
    if not whole_bam or not os.path.exists(whole_bam) or pysam is None:
        return 0, "None", 0.0, 0.0, 0, 0.0, 0.0, False

    partner_records = defaultdict(list)
    total_clipped_molecules = set()

    try:
        with pysam.AlignmentFile(whole_bam, "rb") as bam:
            ref_id = resolve_bam_reference(bam.references, contig_id)
            if not ref_id or ref_id not in bam.references:
                return 0, "None", 0.0, 0.0, 0, 0.0, 0.0, False

            # Scan terminal flanks of candidate contig
            terminal_reads = list(bam.fetch(ref_id, 0, flank_len)) + list(bam.fetch(ref_id, max(0, contig_len - flank_len), contig_len))
            for r in terminal_reads:
                cigar = r.cigartuples or []
                if not cigar:
                    continue
                # If clipped at either end
                if cigar[0][0] in (4, 5) or cigar[-1][0] in (4, 5):
                    total_clipped_molecules.add(r.query_name)

                # Check SA tag for supplementary alignments
                if r.has_tag("SA"):
                    sa_str = r.get_tag("SA")
                    for entry in sa_str.split(";"):
                        if not entry.strip():
                            continue
                        sa_parts = entry.split(",")
                        if len(sa_parts) >= 4:
                            sa_rname = sa_parts[0]
                            sa_pos = int(sa_parts[1])
                            sa_strand = sa_parts[2]
                            if sa_rname != contig_id and sa_rname != ref_id:
                                partner_records[sa_rname].append({
                                    "qname": r.query_name,
                                    "pos": sa_pos,
                                    "strand": sa_strand
                                })

            if not partner_records:
                return 0, "None", 0.0, 0.0, len(total_clipped_molecules), 0.0, 0.0, False

            # Find dominant partner
            dominant_partner = None
            max_partner_reads = 0
            best_records = []

            for p_name, p_list in partner_records.items():
                uniq_qnames = len(set(x["qname"] for x in p_list))
                if uniq_qnames > max_partner_reads:
                    max_partner_reads = uniq_qnames
                    dominant_partner = p_name
                    best_records = p_list

            total_clipped_count = max(1, len(total_clipped_molecules))
            partner_frac = max_partner_reads / float(total_clipped_count)

            # Calculate coordinate clustering MAD
            positions = [x["pos"] for x in best_records]
            if np is not None and positions:
                pos_med = float(np.median(positions))
                pos_mad = float(np.median(np.abs(np.array(positions) - pos_med)))
            elif positions:
                positions.sort()
                pos_med = positions[len(positions) // 2]
                pos_mad = sorted([abs(p - pos_med) for p in positions])[len(positions) // 2]
            else:
                pos_med = 0.0
                pos_mad = 0.0

            # Sample partner contig depth to evaluate relative depth ratio
            partner_depth = 0.0
            rel_depth_ratio = 1.0
            is_low_abundance_artifact = False

            if dominant_partner and dominant_partner in bam.references:
                try:
                    p_len = bam.get_reference_length(dominant_partner)
                    if p_len > 0:
                        # Check header for assembler-reported depth first (e.g. depth-2-2-2 from Myloasm)
                        m_depth = re.search(r'depth-([0-9.]+)', dominant_partner)
                        if m_depth:
                            partner_depth = float(m_depth.group(1))
                        else:
                            # To avoid measuring the shared repeat / chimeric junction flank (which can have spillover depth from candidate),
                            # measure depth across the whole contig body or interior windows away from breakpoints
                            if p_len <= 100000:
                                aligned_bases = sum(
                                    r.query_alignment_length for r in bam.fetch(dominant_partner)
                                    if not r.is_secondary
                                )
                                partner_depth = float(aligned_bases) / float(p_len)
                            else:
                                sample_windows = [(int(p_len * f) - 5000, int(p_len * f) + 5000) for f in (0.25, 0.50, 0.75)]
                                win_depths = []
                                for ws, we in sample_windows:
                                    ws = max(0, ws)
                                    we = min(p_len, we)
                                    wlen = max(1, we - ws)
                                    ab = sum(r.query_alignment_length for r in bam.fetch(dominant_partner, ws, we) if not r.is_secondary)
                                    win_depths.append(float(ab) / float(wlen))
                                partner_depth = float(np.median(win_depths)) if np is not None and win_depths else (sorted(win_depths)[len(win_depths)//2] if win_depths else 0.0)

                        if candidate_depth > 0:
                            rel_depth_ratio = partner_depth / float(candidate_depth)
                        else:
                            rel_depth_ratio = 1.0
                        if rel_depth_ratio < 0.05:
                            is_low_abundance_artifact = True
                except Exception:
                    pass

            return max_partner_reads, dominant_partner, partner_frac, pos_mad, len(total_clipped_molecules), partner_depth, rel_depth_ratio, is_low_abundance_artifact

    except Exception as e:
        print(f"[audit] Warning: partner contig tracing encountered error: {e}", file=sys.stderr)
        return 0, "None", 0.0, 0.0, 0, 0.0, 0.0, False

def parse_blastn_repeats(blast_tsv, contig_len, min_identity=80.0, min_len=200):
    """
    Parse BLASTN self-alignment tabular output (-outfmt 6) faithful to Trigodet et al.
    Filters: identity >= 80%, length > 200 bp, excluding exact self-diagonal.
    Calculates repeat breadth percentage and longest repeat length.
    """
    if not blast_tsv or not os.path.exists(blast_tsv) or contig_len == 0:
        return 0.0, 0, 0

    intervals = []
    longest_repeat = 0
    hit_count = 0

    try:
        with open(blast_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                try:
                    pident = float(parts[2])
                    aln_len = int(parts[3])
                    qstart = int(parts[6])
                    qend = int(parts[7])
                    sstart = int(parts[8])
                    send = int(parts[9])

                    # Exclude trivial diagonal
                    if qstart == sstart and qend == send:
                        continue

                    if pident >= min_identity and aln_len > min_len:
                        start_pos = min(qstart, qend)
                        end_pos = max(qstart, qend)
                        intervals.append((start_pos, end_pos))
                        if aln_len > longest_repeat:
                            longest_repeat = aln_len
                        hit_count += 1
                except ValueError:
                    continue

        if not intervals:
            return 0.0, 0, 0

        # Merge intervals to compute non-redundant repeat breadth
        intervals.sort(key=lambda x: x[0])
        merged = []
        cur_start, cur_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_end:
                cur_end = max(cur_end, e)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = s, e
        merged.append((cur_start, cur_end))

        total_repeat_bp = sum(e - s for s, e in merged)
        repeat_breadth_pct = float(total_repeat_bp) / contig_len * 100.0
        return repeat_breadth_pct, longest_repeat, hit_count

    except Exception as e:
        print(f"[audit] Warning: failed to parse BLAST TSV {blast_tsv}: {e}", file=sys.stderr)
        return 0.0, 0, 0

def parse_bakta_biological_markers(bakta_tsv):
    """
    Parse Bakta output TSV to extract:
      - Independent full-length 5S, 16S, and 23S rRNA gene presence.
      - Canonical 20 amino acid tRNA isotype diversity:
        Normalizes fMet -> Met, Ile2 -> Ile.
        Discards Sec, Pyl, Sup, Undet, pseudo.
    Supports both official 9-column Bakta TSV format:
      (#Sequence Id, Type, Start, Stop, Strand, Locus Tag, Gene, Product, DbXrefs)
    and custom/simplified TSVs with header lines.
    """
    r5s, r16s, r23s = 0, 0, 0
    canonical_isotypes = set()

    if not bakta_tsv or not os.path.exists(bakta_tsv):
        return 0, 0, 0, 0

    try:
        header_map = None
        with open(bakta_tsv) as f:
            for line in f:
                line_str = line.strip()
                if not line_str:
                    continue
                parts = line_str.split("\t")
                if line_str.startswith("#"):
                    clean_cols = [p.strip().lower().lstrip("#").replace(" ", "_") for p in parts]
                    if "type" in clean_cols:
                        header_map = {col: i for i, col in enumerate(clean_cols)}
                    continue

                if len(parts) < 2:
                    continue

                if header_map is not None:
                    type_idx = header_map.get("type", 1)
                    product_idx = header_map.get("product", -1)
                    gene_idx = header_map.get("gene", -1)
                else:
                    type_idx = 1
                    product_idx = 7 if len(parts) > 7 else (2 if len(parts) <= 4 else -1)
                    gene_idx = 6 if len(parts) > 6 else -1

                feat_type = parts[type_idx].lower() if len(parts) > type_idx else ""
                product = parts[product_idx] if (0 <= product_idx < len(parts)) else ""
                gene = parts[gene_idx] if (0 <= gene_idx < len(parts)) else ""
                search_text = f"{product} {gene}".lower()

                if feat_type == "rrna":
                    if "5s" in search_text:
                        r5s += 1
                    elif "16s" in search_text:
                        r16s += 1
                    elif "23s" in search_text:
                        r23s += 1

                elif feat_type == "trna":
                    # e.g., tRNA-Ala(tgc) or tRNA-fMet(cat)
                    match = re.search(r'trna-([A-Za-z0-9]+)', search_text, re.IGNORECASE)
                    if match:
                        raw_aa = match.group(1).capitalize()
                        # Normalization
                        if raw_aa in ("Fmet", "F-met"):
                            norm_aa = "Met"
                        elif raw_aa.startswith("Ile"):
                            norm_aa = "Ile"
                        else:
                            norm_aa = raw_aa

                        if norm_aa in CANONICAL_20_AA:
                            canonical_isotypes.add(norm_aa)

    except Exception as e:
        print(f"[audit] Warning: failed to parse Bakta TSV {bakta_tsv}: {e}", file=sys.stderr)

    return r5s, r16s, r23s, len(canonical_isotypes)

def parse_checkm2(checkm2_tsv, bin_id):
    """Parse CheckM2 summary TSV for Completeness and Contamination."""
    if not checkm2_tsv or not os.path.exists(checkm2_tsv):
        return None, None
    try:
        with open(checkm2_tsv) as f:
            header = None
            for line in f:
                parts = line.strip().split("\t")
                if not parts:
                    continue
                if header is None:
                    header = [p.lower() for p in parts]
                    continue
                row = dict(zip(header, parts))
                name = row.get("name", "")
                if (name == bin_id or name == f"{bin_id}.fa" or
                    name.startswith(f"{bin_id}_") or name.startswith(f"{bin_id}.")):
                    comp = float(row.get("completeness", 0.0))
                    cont = float(row.get("contamination", 0.0))
                    return comp, cont
    except Exception as e:
        print(f"[audit] Warning: failed to parse CheckM2 TSV {checkm2_tsv}: {e}", file=sys.stderr)
    return None, None

def parse_gtdb(gtdb_tsv, bin_id):
    """Parse GTDB-Tk classification summary."""
    if not gtdb_tsv or not os.path.exists(gtdb_tsv):
        return "Unclassified"
    try:
        with open(gtdb_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    user_genome = parts[0]
                    if (user_genome == bin_id or user_genome == f"{bin_id}.fa" or
                        user_genome.startswith(f"{bin_id}_") or user_genome.startswith(f"{bin_id}.")):
                        return parts[1]
    except Exception:
        pass
    return "Unclassified"

def parse_genomad(genomad_tsv, bin_id):
    """
    Parse geNomad summary TSV.
    Extracts raw scores and computes derived replicon assignment and score margin.
    """
    raw_scores = {"chromosome": 0.0, "plasmid": 0.0, "virus": 0.0}
    if not genomad_tsv or not os.path.exists(genomad_tsv):
        return raw_scores, "AMBIGUOUS", 0.0

    try:
        with open(genomad_tsv) as f:
            header = None
            for line in f:
                parts = line.strip().split("\t")
                if not parts:
                    continue
                if header is None:
                    header = [p.lower() for p in parts]
                    continue
                row = dict(zip(header, parts))
                seq_name = row.get("seq_name", "")
                if (seq_name == bin_id or seq_name == f"{bin_id}.fa" or
                    seq_name.startswith(f"{bin_id}_") or bin_id in seq_name):
                    c_s = float(row.get("chromosome_score", 0.0))
                    p_s = float(row.get("plasmid_score", 0.0))
                    v_s = float(row.get("virus_score", 0.0))
                    raw_scores = {"chromosome": c_s, "plasmid": p_s, "virus": v_s}
                    break
    except Exception as e:
        print(f"[audit] Warning: failed to parse geNomad TSV {genomad_tsv}: {e}", file=sys.stderr)

    scores_list = [("CHROMOSOME", raw_scores["chromosome"]),
                   ("PLASMID", raw_scores["plasmid"]),
                   ("VIRUS", raw_scores["virus"])]
    scores_list.sort(key=lambda x: x[1], reverse=True)
    top_label, top_score = scores_list[0]
    second_label, second_score = scores_list[1]
    margin = top_score - second_score

    if top_score < 0.50 or margin < 0.15:
        derived_state = "AMBIGUOUS"
    elif top_label == "VIRUS" and top_score >= 0.70:
        derived_state = "VIRUS_SUPPORTED"
    elif top_label == "PLASMID" and top_score >= 0.70:
        derived_state = "PLASMID_SUPPORTED"
    elif top_label == "CHROMOSOME" and top_score >= 0.70:
        derived_state = "CHROMOSOME_SUPPORTED"
    else:
        derived_state = "CONFLICTING"

    return raw_scores, derived_state, margin

@dataclass
class EvidenceRecord:
    metric: str
    value: Any
    threshold: str
    source: str
    passed: bool
    description: str = ""


@dataclass
class CircularityResult:
    tier: str                      # VERIFIED_CIRCULAR, SUPPORTED_CIRCULAR, CANDIDATE_CIRCULAR, AMBIGUOUS_CIRCULAR, REJECTED_LINEAR_OR_BROKEN
    evidence_type: str             # CONCORDANT_DUAL_OMICS, LONG_READ_DEFINITIVE, ORTHOGONAL_REINFORCEMENT, MOLECULE_LENGTH_SUPPORTED, LONG_READ_ONLY, CANDIDATE_MARGINAL, CONFLICTING_JUNCTION_SIGNALS, NO_JUNCTION_SUPPORT, PHYSICAL_BREAKPOINT_GAP
    spanning_reads: int
    support_fraction: float
    clip_fraction: float
    anchor_status: str
    illumina_status: str
    notes: str
    records: List[EvidenceRecord] = field(default_factory=list)


@dataclass
class ComplexityResult:
    bio_tier: str                  # CLONAL_SIMPLE, POPULATION_SV_POLYMORPHIC, STRUCTURALLY_HETEROGENEOUS, REPEAT_DOMINATED, MULTI_COMPLEX
    graph_tier: str                # GRAPH_SIMPLE, LOW_FREQUENCY_ALTERNATIVE_PATH, COMPETING_GRAPH_BRANCH, UNRESOLVED_BRANCH
    alt_path_contig: str           # e.g. "u1804999" or "None"
    rdr: float                     # Continuous float
    rdr_signal: str                # RARE_GRAPH_SIGNAL (<0.01), LOW_FREQUENCY_GRAPH_SIGNAL (0.01-0.10), POPULATION_COMPETITOR_SIGNAL (>0.10), NO_ALTERNATIVE_PATH
    notes: str                     # e.g. "SV(1.48Mb:20.3% AF_del)+Alt_Path(u1804999,RDR=0.002,rare)"
    sv_sites: List[Dict[str, Any]] = field(default_factory=list)
    clipping_sites: List[Dict[str, Any]] = field(default_factory=list)
    records: List[EvidenceRecord] = field(default_factory=list)


@dataclass
class AuditResult:
    assembly_grade: str            # HIGH_CONFIDENCE_CIRCULAR, MODERATE_CONFIDENCE_CIRCULAR, PUTATIVE_CIRCULAR, REJECTED_NON_CIRCULAR
    circularity: CircularityResult
    complexity: ComplexityResult
    gate1_status: str
    gate1_notes: str
    gate2_status: str
    gate3_status: str
    gate3_notes: str
    gate4_status: str
    gate5_status: str
    final_grade: str               # Human-readable & machine-friendly composite string

    def __iter__(self):
        """Enable tuple unpacking for 100% backward compatibility with legacy 7-tuple returns."""
        yield self.gate1_status
        yield self.gate1_notes
        yield self.gate2_status
        yield self.gate3_status
        yield self.gate3_notes
        yield self.gate4_status
        yield self.final_grade

    def __getitem__(self, index):
        return list(self)[index]


def evaluate_four_gates(
    contig_len,
    # Gate 1 physical support
    spanning_unique_reads,
    support_fraction,
    clip_fraction,
    unresolved_zero_bp,
    has_high_risk_zero_block,
    clipping_sites_count,
    partner_reads,
    partner_fraction,
    partner_name,
    anchor_status,
    # Gate 2 coherence
    primary_mean_depth,
    primary_breadth_1x,
    repeat_breadth_pct,
    longest_repeat_bp,
    # Gate 3 biological markers
    comp,
    cont,
    r5s,
    r16s,
    r23s,
    trna_aa_count,
    # Gate 4 replicon state
    replicon_state,
    # Gate 5 orthogonal short-read state
    gate5_status="NOT_APPLICABLE",
    # Tuning params
    min_spanning_reads=5,
    min_support_fraction=0.50,
    max_clip_fraction=0.15,
    min_coverage_depth=10,
    clipping_sites_detail=None,
    zero_gaps_detail=None,
    is_low_abundance_partner=False,
    partner_rel_depth_ratio=1.0,
    effective_flank=None,
    gate5_metrics=None,
):
    """
    Arbitrates the Two-Dimensional Decoupled Framework:
    - Axis 1: Circularity Confidence (Topology)
    - Axis 2: Population & Graph Complexity (Clonality & Structural Heterogeneity)
    - Assembly Usability: Downstream application clearance
    """
    # 1. Classify Internal Structural Sites (SV vs Misassembly)
    sv_sites = []
    misassembly_sites = []
    if clipping_sites_detail and isinstance(clipping_sites_detail, list):
        for cs in clipping_sites_detail:
            if cs.get("is_sv_polymorphism", False):
                sv_sites.append(cs)
            else:
                misassembly_sites.append(cs)
    elif clipping_sites_count > 0:
        misassembly_sites = [{"pos": "unknown", "clipped_reads": clipping_sites_count, "depth": 0, "clipping_ratio": 1.0}]

    # 2. Axis 2 - Level 2: Assembly Graph State
    alt_contig = partner_name if partner_name and partner_name != "None" else "None"
    rdr = float(partner_rel_depth_ratio)

    if alt_contig != "None" and partner_reads >= 3:
        if rdr < 0.01:
            rdr_signal = "RARE_GRAPH_SIGNAL"
            graph_tier = "LOW_FREQUENCY_ALTERNATIVE_PATH"
        elif 0.01 <= rdr <= 0.10:
            rdr_signal = "LOW_FREQUENCY_GRAPH_SIGNAL"
            graph_tier = "LOW_FREQUENCY_ALTERNATIVE_PATH"
        else:
            rdr_signal = "POPULATION_COMPETITOR_SIGNAL"
            graph_tier = "COMPETING_GRAPH_BRANCH"
    elif anchor_status.startswith("AMBIGUOUS"):
        rdr_signal = "NO_ALTERNATIVE_PATH"
        graph_tier = "UNRESOLVED_BRANCH"
    else:
        rdr_signal = "NO_ALTERNATIVE_PATH"
        graph_tier = "GRAPH_SIMPLE"

    # 3. Axis 2 - Level 1: Biological Complexity State
    has_sv = len(sv_sites) > 0
    has_clipping = len(misassembly_sites) > 0
    has_repeat = repeat_breadth_pct >= 70.0 or (contig_len < 50000 and longest_repeat_bp >= 5000)
    has_alt_graph = (graph_tier != "GRAPH_SIMPLE") and not is_low_abundance_partner

    comp_factors = sum([1 for f in [has_sv, has_clipping, has_repeat, (graph_tier != "GRAPH_SIMPLE")] if f])

    if comp_factors >= 2:
        bio_tier = "MULTI_COMPLEX"
    elif has_sv:
        bio_tier = "POPULATION_SV_POLYMORPHIC"
    elif has_clipping:
        bio_tier = "STRUCTURALLY_HETEROGENEOUS"
    elif has_repeat:
        bio_tier = "REPEAT_DOMINATED"
    elif has_alt_graph:
        bio_tier = "STRUCTURALLY_HETEROGENEOUS"
    else:
        bio_tier = "CLONAL_SIMPLE"

    # Complexity Summary Notes
    notes_parts = []
    if sv_sites:
        sv_strs = [f"{s['pos']}(AF={s.get('del_freq', 0.0):.1%})" for s in sv_sites[:3]]
        notes_parts.append(f"SV({len(sv_sites)}_loci:{','.join(sv_strs)})")
    if misassembly_sites:
        notes_parts.append(f"Internal_Clipping({len(misassembly_sites)}_loci)")
    if alt_contig != "None" and partner_reads >= 3:
        sig_label = "rare" if rdr < 0.01 else ("low_freq" if rdr <= 0.10 else "competitor")
        notes_parts.append(f"Alt_Path({alt_contig},RDR={rdr:.3f},{sig_label})")
    if has_repeat:
        notes_parts.append(f"Repeat({repeat_breadth_pct:.1f}%)")
    comp_notes = "+".join(notes_parts) if notes_parts else "None(Clonal_Pure)"

    # 4. Axis 1: Circularity Confidence
    gate5_is_pass = gate5_status in ("PASS", "PASS_HIGH_CONFIDENCE")
    gate5_is_conflict = (gate5_status == "JUNCTION_CONFLICT")
    gate5_is_mod = gate5_status in ("PASS_MODERATE", "RECRUITMENT_CAUTION")
    gate5_has_reads = False
    if gate5_metrics and isinstance(gate5_metrics, dict):
        g5_s = 0
        g5_b = 0
        try:
            val_s = gate5_metrics.get("Illumina_Rotated_Spanning_Reads_Unique", 0)
            if val_s not in ("NA", "None", None, ""):
                g5_s = int(float(val_s))
            val_b = gate5_metrics.get("Illumina_Rotated_Bridging_Pairs_Unique", 0)
            if val_b not in ("NA", "None", None, ""):
                g5_b = int(float(val_b))
        except Exception:
            pass
        gate5_has_reads = (g5_s > 0 or g5_b > 0)

    is_zero_gap_reject = has_high_risk_zero_block and unresolved_zero_bp >= 500
    is_zero_spanning_reject = (spanning_unique_reads == 0 and not gate5_has_reads and not gate5_is_pass)
    is_partner_conflict = (partner_reads >= 3 and partner_fraction >= 0.50 and rdr > 0.10)
    is_orthogonal_conflict = gate5_is_conflict

    if is_zero_gap_reject or is_zero_spanning_reject:
        circ_tier = "REJECTED_LINEAR_OR_BROKEN"
        circ_evidence = "PHYSICAL_BREAKPOINT_GAP" if is_zero_gap_reject else "NO_JUNCTION_SUPPORT"
    elif is_partner_conflict or is_orthogonal_conflict:
        circ_tier = "AMBIGUOUS_CIRCULAR"
        circ_evidence = "CONFLICTING_JUNCTION_SIGNALS"
    else:
        # Evaluate long-read physical junction spanning strength
        is_long_read_verified = (
            (spanning_unique_reads >= min_spanning_reads and support_fraction >= min_support_fraction and clip_fraction <= 0.25 and unresolved_zero_bp == 0)
            or (spanning_unique_reads >= 20 and unresolved_zero_bp == 0)
        )
        is_long_read_supported = (
            (spanning_unique_reads >= 3 and support_fraction >= 0.40 and unresolved_zero_bp == 0)
            or (spanning_unique_reads >= 3 and ((effective_flank and effective_flank >= 5000) or contig_len <= 100000) and unresolved_zero_bp == 0)
        )

        if is_long_read_verified:
            circ_tier = "VERIFIED_CIRCULAR"
            if gate5_is_pass:
                circ_evidence = "CONCORDANT_DUAL_OMICS"
            elif gate5_status in ("NOT_APPLICABLE", "NOT_EVALUATED", "INSUFFICIENT_DEPTH"):
                circ_evidence = "LONG_READ_DEFINITIVE"
            else:
                circ_evidence = "LONG_READ_DEFINITIVE"
        elif (gate5_is_pass or gate5_is_mod or gate5_has_reads) and (spanning_unique_reads >= 1 or primary_mean_depth >= 5.0):
            # Orthogonal Reinforcement (independent evidence confirmation, e.g. u333097)
            circ_tier = "SUPPORTED_CIRCULAR"
            circ_evidence = "ORTHOGONAL_REINFORCEMENT"
        elif is_long_read_supported:
            circ_tier = "SUPPORTED_CIRCULAR"
            if spanning_unique_reads >= 3 and ((effective_flank and effective_flank >= 5000) or contig_len <= 100000):
                circ_evidence = "MOLECULE_LENGTH_SUPPORTED"
            else:
                circ_evidence = "LONG_READ_ONLY"
        elif spanning_unique_reads >= 1:
            circ_tier = "CANDIDATE_CIRCULAR"
            circ_evidence = "CANDIDATE_MARGINAL"
        else:
            circ_tier = "REJECTED_LINEAR_OR_BROKEN"
            circ_evidence = "NO_JUNCTION_SUPPORT"

    # 5. Legacy Gate 1-4 arbitration for backward compatibility
    gate1_failures = []
    gate1_notes_list = []
    if is_partner_conflict:
        gate1_failures.append(f"Premature_Circularization_Risk({partner_reads}_reads->{partner_name})")
    elif alt_contig != "None" and partner_reads >= 3:
        if is_low_abundance_partner:
            gate1_notes_list.append(f"Low_Abundance_Partner_Artifact({alt_contig},RDR={rdr:.3f})")
        else:
            gate1_notes_list.append(f"Alternative_Assembly_Path({alt_contig},RDR={rdr:.3f})")

    if spanning_unique_reads < min_spanning_reads:
        gate1_failures.append(f"Low_Junction_Spanning({spanning_unique_reads}<{min_spanning_reads})")
    if support_fraction < min_support_fraction:
        gate1_failures.append(f"Low_Junction_Support({support_fraction:.2f}<{min_support_fraction})")
    if clip_fraction > max_clip_fraction:
        gate1_failures.append(f"High_Junction_Clipping({clip_fraction:.2f}>{max_clip_fraction})")
    if anchor_status.startswith("AMBIGUOUS"):
        gate1_failures.append("Ambiguous_Terminal_Anchors")
    if unresolved_zero_bp > 0:
        gap_str = ""
        if zero_gaps_detail:
            gap_locs = [f"{b['start']}-{b['end']}" for b in zero_gaps_detail[:3] if isinstance(b, dict) and 'start' in b and 'end' in b]
            if gap_locs:
                gap_str = f"[{','.join(gap_locs)}]"
        if has_high_risk_zero_block:
            gate1_failures.append(f"High_Risk_Zero_Coverage_Gap(>{500}bp,total={unresolved_zero_bp}bp{gap_str})")
        else:
            gate1_failures.append(f"Unresolved_Zero_Coverage({unresolved_zero_bp}bp{gap_str})")
    if misassembly_sites:
        clip_locs = [f"pos={s['pos']}" for s in misassembly_sites if isinstance(s, dict) and 'pos' in s]
        if clip_locs:
            gate1_failures.append(f"Systematic_Internal_Clipping_Sites({len(misassembly_sites)}:{','.join(clip_locs[:3])})")
        else:
            gate1_failures.append(f"Systematic_Internal_Clipping_Sites({len(misassembly_sites)})")
    if sv_sites:
        sv_strs = [f"pos={s['pos']}(AF_del={s.get('del_freq', 0.0):.1%})" for s in sv_sites[:3]]
        gate1_notes_list.append(f"Population_SV_Polymorphism({len(sv_sites)}_loci:{','.join(sv_strs)})")

    has_sv_heterogeneity = len(sv_sites) > 0 and len(misassembly_sites) == 0
    if gate1_failures:
        gate1_status = "FAIL"
    elif has_sv_heterogeneity:
        gate1_status = "PASS_WITH_SV_HETEROGENEITY"
    else:
        gate1_status = "PASS"

    gate1_summary_parts = gate1_failures + gate1_notes_list
    gate1_summary = "; ".join(gate1_summary_parts) if gate1_summary_parts else "Supported"

    # Gate 2 Coherence
    gate2_notes = []
    if primary_mean_depth < min_coverage_depth:
        gate2_notes.append(f"Low_Depth({primary_mean_depth:.1f}<{min_coverage_depth}x)")
    if primary_breadth_1x < 99.0:
        gate2_notes.append(f"Imperfect_Breadth({primary_breadth_1x:.2f}%)")
    if contig_len < 50000 and repeat_breadth_pct >= 70.0:
        gate2_notes.append(f"Repeat_Dominated_Risk({repeat_breadth_pct:.1f}%)")
    gate2_status = "COHERENT" if not gate2_notes else "; ".join(gate2_notes)

    # Gate 3 Biological Completeness
    gate3_failures = []
    is_chromosome_size = contig_len >= 500000 and replicon_state != "PLASMID_SUPPORTED" and replicon_state != "VIRUS_SUPPORTED"
    if is_chromosome_size:
        if comp is not None and comp <= 90.0:
            gate3_failures.append(f"Completeness<=90%({comp}%)")
        if cont is not None and cont >= 5.0:
            gate3_failures.append(f"Contamination>=5%({cont}%)")
        if not (r5s >= 1 and r16s >= 1 and r23s >= 1):
            gate3_failures.append(f"Incomplete_rRNA(5S:{r5s},16S:{r16s},23S:{r23s})")
        if trna_aa_count < 18:
            gate3_failures.append(f"Low_tRNA_Isotypes({trna_aa_count}<18/20)")
        gate3_status = "PASS" if not gate3_failures else "INCOMPLETE"
    else:
        gate3_status = "NOT_APPLICABLE (Non-Chromosomal/Small Replicon)"
    gate3_summary = "; ".join(gate3_failures) if gate3_failures else "Complete"

    # Gate 4 Replicon State
    gate4_status = replicon_state

    # 6. Assembly Usability Grade
    if circ_tier == "REJECTED_LINEAR_OR_BROKEN":
        assembly_grade = "REJECTED_NON_CIRCULAR"
    elif circ_tier == "AMBIGUOUS_CIRCULAR":
        assembly_grade = "PUTATIVE_CIRCULAR"
    elif circ_tier == "CANDIDATE_CIRCULAR":
        assembly_grade = "PUTATIVE_CIRCULAR"
    else:
        # Verified or Supported
        if bio_tier == "CLONAL_SIMPLE" and gate3_status in ("PASS", "NOT_APPLICABLE (Non-Chromosomal/Small Replicon)"):
            assembly_grade = "HIGH_CONFIDENCE_CIRCULAR"
        elif bio_tier in ("POPULATION_SV_POLYMORPHIC", "MULTI_COMPLEX", "STRUCTURALLY_HETEROGENEOUS"):
            assembly_grade = "HIGH_CONFIDENCE_CIRCULAR"
        elif gate3_status == "INCOMPLETE" and is_chromosome_size:
            assembly_grade = "MODERATE_CONFIDENCE_CIRCULAR"
        else:
            assembly_grade = "HIGH_CONFIDENCE_CIRCULAR"

    # 7. Final Composite Grade (Preserving backward compatibility)
    if circ_tier == "REJECTED_LINEAR_OR_BROKEN":
        if is_zero_gap_reject or misassembly_sites:
            final_grade = "POTENTIAL_MISASSEMBLY (Physical Breakpoint / Unresolved Gap)"
        elif spanning_unique_reads == 0:
            final_grade = "INSUFFICIENT_PHYSICAL_EVIDENCE (Low Depth / Non-Informative)"
        else:
            final_grade = "REJECTED_NON_CIRCULAR (Physical Breakpoint / Zero Coverage Gap)"
    elif is_partner_conflict:
        final_grade = "PREMATURE_CIRCULARIZATION_SUSPECT (Consistent Partner Alignment)"
    elif circ_tier == "AMBIGUOUS_CIRCULAR":
        if anchor_status.startswith("AMBIGUOUS"):
            final_grade = "SUSPECT_CIRCULAR_JUNCTION (Ambiguous Terminal Anchors)"
        elif is_orthogonal_conflict:
            final_grade = "SUSPECT_CIRCULAR_JUNCTION (Illumina Orthogonal Conflict)"
        else:
            final_grade = "AMBIGUOUS_CIRCULAR_REPLICON (Conflicting Graph / Orthogonal Signals)"
    elif circ_tier == "CANDIDATE_CIRCULAR":
        final_grade = "PUTATIVE_CIRCULAR_REPLICON (Manual Review Needed)"
    else:
        if replicon_state == "VIRUS_SUPPORTED":
            replicon_label = "VIRAL_REPLICON"
        elif replicon_state == "PLASMID_SUPPORTED" or (contig_len < 300000 and replicon_state != "CHROMOSOME_SUPPORTED"):
            replicon_label = "PLASMID"
        else:
            replicon_label = "CHROMOSOME"
        prefix = "ORTHOGONALLY_VERIFIED_" if circ_evidence in ("CONCORDANT_DUAL_OMICS", "ORTHOGONAL_REINFORCEMENT") else "READ_SUPPORTED_"

        if replicon_label == "CHROMOSOME":
            base_desc = f"{prefix}HIGH_QUALITY_CIRCULAR_CHROMOSOME" if gate3_status == "PASS" else f"{prefix}CIRCULAR_CHROMOSOME"
        elif replicon_label == "PLASMID":
            base_desc = f"{prefix}CIRCULAR_PLASMID" if prefix == "ORTHOGONALLY_VERIFIED_" else "HIGH_CONFIDENCE_CIRCULAR_PLASMID"
        else:
            base_desc = f"{prefix}CIRCULAR_VIRAL_REPLICON"

        if bio_tier == "MULTI_COMPLEX":
            final_grade = f"{base_desc}_WITH_POPULATION_COMPLEXITY"
        elif bio_tier == "POPULATION_SV_POLYMORPHIC" or has_sv_heterogeneity:
            final_grade = f"{base_desc}_WITH_POPULATION_HETEROGENEITY"
        elif bio_tier == "STRUCTURALLY_HETEROGENEOUS":
            final_grade = f"{base_desc}_WITH_STRUCTURAL_HETEROGENEITY"
        elif bio_tier == "REPEAT_DOMINATED":
            final_grade = f"{base_desc}_REPEAT_DOMINATED"
        else:
            final_grade = base_desc

        if gate5_status == "RECRUITMENT_CAUTION":
            final_grade = f"{final_grade} (Short-Read Recruitment Caution)"

    circ_res = CircularityResult(
        tier=circ_tier,
        evidence_type=circ_evidence,
        spanning_reads=spanning_unique_reads,
        support_fraction=support_fraction,
        clip_fraction=clip_fraction,
        anchor_status=anchor_status,
        illumina_status=gate5_status,
        notes=circ_evidence,
        records=[]
    )
    comp_res = ComplexityResult(
        bio_tier=bio_tier,
        graph_tier=graph_tier,
        alt_path_contig=alt_contig,
        rdr=rdr,
        rdr_signal=rdr_signal,
        notes=comp_notes,
        sv_sites=sv_sites,
        clipping_sites=misassembly_sites,
        records=[]
    )
    return AuditResult(
        assembly_grade=assembly_grade,
        circularity=circ_res,
        complexity=comp_res,
        gate1_status=gate1_status,
        gate1_notes=gate1_summary,
        gate2_status=gate2_status,
        gate3_status=gate3_status,
        gate3_notes=gate3_summary,
        gate4_status=gate4_status,
        gate5_status=gate5_status,
        final_grade=final_grade
    )


evaluate_five_gates = evaluate_four_gates


def parse_illumina_qc_row(illumina_qc_tsv, bin_id):
    """
    Parses circular_illumina_qc_summary.tsv to find metrics for bin_id.
    Matches bin_id with or without .fa/.fasta/.fna extension.
    Returns dict of metrics or default NAs.
    """
    default_res = {
        "Illumina_Relative_Abundance_Pct": "NA",
        "Illumina_Catalogue_Mapped_Pct": "NA",
        "Illumina_Catalogue_Unmapped_Pct": "NA",
        "Illumina_Genome_Trimmed_Depth": "NA",
        "Illumina_Genome_Covered_Fraction": "NA",
        "Illumina_Aligned_Read_Count": "NA",
        "Illumina_ANIr_Pct": "NA",
        "Illumina_Competitive_Trimmed_Depth": "NA",
        "Illumina_Competitive_Depth_Ratio_C": "NA",
        "Illumina_Strict_Depth_Ratio_S": "NA",
        "Illumina_Rotated_Junction_Depth_100bp": "NA",
        "Illumina_Rotated_Junction_Local_Ratio": "NA",
        "Illumina_Rotated_Zero_Coverage_Run_bp": "NA",
        "Illumina_Rotated_Spanning_Reads_Unique": "NA",
        "Illumina_Rotated_Bridging_Pairs_Unique": "NA",
        "Gate5_Illumina_Orthogonal_Status": "NOT_APPLICABLE",
        "Gate5_Notes": "No_Short_Reads"
    }

    if not illumina_qc_tsv or not os.path.exists(illumina_qc_tsv):
        return default_res

    clean_bin = bin_id
    for ext in (".fa", ".fasta", ".fna"):
        if clean_bin.endswith(ext):
            clean_bin = clean_bin[:-len(ext)]
            break

    target_row = None
    try:
        with open(illumina_qc_tsv, "r") as f:
            header_line = f.readline()
            if not header_line:
                return default_res
            headers = [h.strip() for h in header_line.split("\t")]
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = [p.strip() for p in line.split("\t")]
                row = dict(zip(headers, parts))
                cg = row.get("circular_genome", "")
                for ext in (".fa", ".fasta", ".fna"):
                    if cg.endswith(ext):
                        cg = cg[:-len(ext)]
                        break
                if cg == clean_bin or row.get("circular_genome", "") == bin_id:
                    target_row = row
                    break
    except Exception as e:
        print(f"[audit] Warning reading illumina-qc-tsv: {e}", file=sys.stderr)
        return default_res

    if not target_row:
        default_res["Gate5_Notes"] = "Contig_Not_In_QC_Summary"
        return default_res

    ov_status = target_row.get("overall_status", "UNKNOWN")
    flag_gap = target_row.get("flag_junction_gap", "FALSE") == "TRUE"
    flag_recruit = target_row.get("flag_possible_homologous_recruitment", "FALSE") == "TRUE"

    if flag_gap or ov_status == "JUNCTION_SUSPECT":
        gate5_status = "JUNCTION_CONFLICT"
    elif flag_recruit or ov_status == "RECRUITMENT_CAUTION":
        gate5_status = "RECRUITMENT_CAUTION"
    elif ov_status == "PASS_HIGH_CONFIDENCE":
        gate5_status = "PASS"
    elif ov_status == "PASS_MODERATE":
        gate5_status = "PASS_MODERATE"
    elif ov_status in ("LOW_DEPTH_INSUFFICIENT", "LOW_INFORMATION"):
        gate5_status = "INSUFFICIENT_DEPTH"
    elif ov_status == "UNVERIFIED_NO_ROTATED_DATA":
        gate5_status = "NOT_EVALUATED"
    else:
        gate5_status = ov_status

    flags = []
    if ov_status == "UNVERIFIED_NO_ROTATED_DATA":
        flags.append("No_Rotated_Junction_Data")
    if target_row.get("flag_low_relative_abundance") == "TRUE":
        flags.append(f"Low_RA({target_row.get('relative_abundance_percent')}%)")
    if target_row.get("flag_low_depth") == "TRUE":
        flags.append(f"Low_Depth({target_row.get('genome_trimmed_mean_depth')}x)")
    if target_row.get("flag_low_breadth") == "TRUE":
        flags.append(f"Low_Breadth({target_row.get('genome_covered_fraction')})")
    if target_row.get("flag_possible_homologous_recruitment") == "TRUE":
        flags.append(f"Recruitment_Drop(C={target_row.get('competitive_depth_ratio')})")
    if target_row.get("flag_identity_sensitive") == "TRUE":
        flags.append(f"Identity_Sensitive(S={target_row.get('strict_depth_ratio')})")
    if target_row.get("flag_junction_depth_drop") == "TRUE":
        flags.append(f"Junction_Depth_Drop(R_local={target_row.get('junction_local_ratio')})")
    if target_row.get("flag_junction_gap") == "TRUE":
        flags.append(f"Junction_Zero_Gap({target_row.get('junction_zero_coverage_run_bp')}bp)")
    if target_row.get("flag_no_junction_span") == "TRUE":
        flags.append(f"No_Junction_Span(Span={target_row.get('junction_spanning_reads_unique')},Bridge={target_row.get('junction_bridging_pairs_unique')})")

    if not flags:
        gate5_notes = f"Orthogonally_Supported(Span={target_row.get('junction_spanning_reads_unique')},Bridge={target_row.get('junction_bridging_pairs_unique')},Depth={target_row.get('genome_trimmed_mean_depth')}x)"
    else:
        gate5_notes = "; ".join(flags)

    return {
        "Illumina_Relative_Abundance_Pct": target_row.get("relative_abundance_percent", "NA"),
        "Illumina_Catalogue_Mapped_Pct": target_row.get("catalogue_mapped_percent", "NA"),
        "Illumina_Catalogue_Unmapped_Pct": target_row.get("catalogue_unmapped_percent", "NA"),
        "Illumina_Genome_Trimmed_Depth": target_row.get("genome_trimmed_mean_depth", "NA"),
        "Illumina_Genome_Covered_Fraction": target_row.get("genome_covered_fraction", "NA"),
        "Illumina_Aligned_Read_Count": target_row.get("genome_aligned_read_count", "NA"),
        "Illumina_ANIr_Pct": target_row.get("genome_anir_percent", "NA"),
        "Illumina_Competitive_Trimmed_Depth": target_row.get("contig_trimmed_mean_depth", "NA"),
        "Illumina_Competitive_Depth_Ratio_C": target_row.get("competitive_depth_ratio", "NA"),
        "Illumina_Strict_Depth_Ratio_S": target_row.get("strict_depth_ratio", "NA"),
        "Illumina_Rotated_Junction_Depth_100bp": target_row.get("junction_median_depth_100bp", "NA"),
        "Illumina_Rotated_Junction_Local_Ratio": target_row.get("junction_local_ratio", "NA"),
        "Illumina_Rotated_Zero_Coverage_Run_bp": target_row.get("junction_zero_coverage_run_bp", "NA"),
        "Illumina_Rotated_Spanning_Reads_Unique": target_row.get("junction_spanning_reads_unique", "NA"),
        "Illumina_Rotated_Bridging_Pairs_Unique": target_row.get("junction_bridging_pairs_unique", "NA"),
        "Gate5_Illumina_Orthogonal_Status": gate5_status,
        "Gate5_Notes": gate5_notes
    }

def main():
    parser = argparse.ArgumentParser(description="Trigodet-informed assembly verification for circular contigs.")
    parser.add_argument("--bin", required=True, help="Contig / Bin ID")
    parser.add_argument("--contig-fa", required=True, help="Path to original contig FASTA")
    parser.add_argument("--whole-bam", required=True, help="BAM of ONT reads mapped to whole assembly")
    parser.add_argument("--junction-bam", help="BAM of terminal reads mapped to synthetic junction")
    parser.add_argument("--self-blast", help="BLASTN tabular output (-outfmt 6) of contig against itself")
    parser.add_argument("--raw-reads", help="Path to raw ONT FASTQ for zero-coverage rescue")
    parser.add_argument("--bakta-tsv", help="Bakta output TSV")
    parser.add_argument("--checkm2-tsv", help="CheckM2 summary TSV")
    parser.add_argument("--gtdb-tsv", help="GTDB-Tk summary TSV")
    parser.add_argument("--genomad-tsv", help="geNomad summary TSV")
    parser.add_argument("--illumina-qc-tsv", help="Summary TSV from compile_circular_illumina_qc.py")

    # Tuning Parameters (ont_provisional_v1 profile)
    parser.add_argument("--terminal-flank", type=int, default=5000, help="Synthetic junction flank length (bp)")
    parser.add_argument("--min-junction-flank", type=int, default=500, help="Minimum read extension on both sides of J")
    parser.add_argument("--min-junction-spanning-reads", type=int, default=5, help="Absolute spanning read threshold")
    parser.add_argument("--min-junction-support-fraction", type=float, default=0.50, help="Junction support fraction threshold")
    parser.add_argument("--max-junction-clip-fraction", type=float, default=0.15, help="Junction clipping fraction threshold")
    parser.add_argument("--clipping-ratio-threshold", type=float, default=0.80, help="anvi'o clipping threshold (> 0.80)")
    parser.add_argument("--min-coverage-depth", type=int, default=10, help="Minimum depth for high-confidence clipping")
    parser.add_argument("--zero-cov-high-risk-threshold", type=int, default=500, help="bp threshold for high-risk zero coverage")
    parser.add_argument("--out-tsv", required=True, help="Output TSV summary path")

    args = parser.parse_args()

    # 1. Contig basic length & effective flank
    contig_len, _ = get_fasta_length(args.contig_fa)
    effective_flank = min(args.terminal_flank, contig_len // 2) if contig_len > 0 else 0

    # 2. Coverage arrays (Primary vs All-alignment support)
    primary_cov, all_cov = calculate_coverage_arrays(args.whole_bam, args.bin, contig_len)

    # 3. Primary coherence metrics
    (p_mean, p_median, p_mad, p_cv,
     p_breadth_1x, p_breadth_3x) = compute_primary_coherence_metrics(primary_cov, contig_len)

    # 4. All-alignment zero coverage
    zero_blocks, total_zero_bp, max_zero_block, high_risk_blocks = find_zero_coverage_blocks(all_cov, contig_len)
    all_support_breadth = float(contig_len - total_zero_bp) / contig_len * 100.0 if contig_len > 0 else 0.0
    has_high_risk = len(high_risk_blocks) > 0
    unresolved_zero_bp = total_zero_bp  # Default unresolved; could be rescued by raw-read search if available

    # 5. Whole-assembly systematic clipping
    clipping_sites = find_systematic_clipping_sites(
        args.whole_bam, args.bin, contig_len, all_cov,
        min_depth=args.min_coverage_depth,
        clipping_ratio_threshold=args.clipping_ratio_threshold
    )
    clipping_count = len(clipping_sites)

    # 6. Targeted synthetic junction alignment
    (span_reads, clip_reads, inf_reads,
     supp_frac, j_clip_frac,
     ci_lower, ci_upper, spanning_qnames) = analyze_targeted_synthetic_junction(
        args.junction_bam, effective_flank, min_flank=args.min_junction_flank
    )

    # 7. Methanothrix partner clustering & low-abundance chimera detection
    (partner_reads, partner_name,
     partner_frac, partner_mad,
     total_clipped_mol, partner_depth,
     rel_depth_ratio, is_low_abundance_artifact) = detect_methanothrix_partner_clustering(
        args.whole_bam, args.bin, contig_len, flank_len=effective_flank, candidate_depth=p_mean
    )

    # 8. Terminal anchor uniqueness in whole assembly (exempting low-abundance artifacts)
    low_partners = [partner_name] if is_low_abundance_artifact and partner_name != "None" else []
    anchor_status, anchor_ambig_count = check_terminal_anchor_uniqueness(
        args.whole_bam, args.bin, contig_len, spanning_qnames, effective_flank,
        low_abundance_partners=low_partners
    )

    # 9. BLASTN repeat breadth
    rep_breadth_pct, longest_rep_bp, rep_hit_count = parse_blastn_repeats(
        args.self_blast, contig_len, min_identity=80.0, min_len=200
    )

    # 10. Biological markers
    r5s, r16s, r23s, trna_aa_count = parse_bakta_biological_markers(args.bakta_tsv)
    comp, cont = parse_checkm2(args.checkm2_tsv, args.bin)
    lineage = parse_gtdb(args.gtdb_tsv, args.bin)
    genomad_raw, replicon_state, replicon_margin = parse_genomad(args.genomad_tsv, args.bin)

    # 11. Parse optional Gate 5 Illumina short-read metrics
    ill_metrics = parse_illumina_qc_row(args.illumina_qc_tsv, args.bin)

    # 12. Arbitrate Gates
    sv_sites_list = [cs for cs in clipping_sites if cs.get("is_sv_polymorphism", False)]
    sv_detail_str = ";".join([f"{cs['pos']}:AF_del={cs.get('del_freq', 0.0):.1%}(del_bp={cs.get('del_span_bp', 0)})" for cs in sv_sites_list]) if sv_sites_list else "None"

    clipping_detail_str = ";".join([f"{cs['pos']}:{cs['clipped_reads']}/{cs['depth']}:{cs['clipping_ratio']:.1%}" for cs in clipping_sites]) if clipping_sites else "None"
    gaps_to_report = high_risk_blocks if high_risk_blocks else [b for b in zero_blocks if b["len"] >= 100]
    zero_gaps_detail_str = ";".join([f"{b['start']}-{b['end']}:{b['len']}bp" for b in gaps_to_report]) if gaps_to_report else "None"

    audit_res = evaluate_four_gates(
        contig_len,
        span_reads, supp_frac, j_clip_frac,
        unresolved_zero_bp, has_high_risk, clipping_count,
        partner_reads, partner_frac, partner_name, anchor_status,
        p_mean, p_breadth_1x, rep_breadth_pct, longest_rep_bp,
        comp, cont, r5s, r16s, r23s, trna_aa_count,
        replicon_state,
        gate5_status=ill_metrics["Gate5_Illumina_Orthogonal_Status"],
        min_spanning_reads=args.min_junction_spanning_reads,
        min_support_fraction=args.min_junction_support_fraction,
        max_clip_fraction=args.max_junction_clip_fraction,
        min_coverage_depth=args.min_coverage_depth,
        clipping_sites_detail=clipping_sites,
        zero_gaps_detail=gaps_to_report,
        is_low_abundance_partner=is_low_abundance_artifact,
        partner_rel_depth_ratio=rel_depth_ratio
    )
    gate1_status = audit_res.gate1_status
    gate1_notes = audit_res.gate1_notes
    gate2_status = audit_res.gate2_status
    gate3_status = audit_res.gate3_status
    gate3_notes = audit_res.gate3_notes
    gate4_status = audit_res.gate4_status
    final_grade = audit_res.final_grade

    # 13. Write single-row TSV report
    headers = [
        "Contig_ID", "Length_bp", "Effective_Junction_Flank_bp",
        "Primary_Mean_Depth", "Primary_Median_Depth", "Primary_Depth_MAD", "Primary_Depth_CV",
        "Primary_Coverage_Breadth_1x_Pct", "Primary_Coverage_Breadth_3x_Pct",
        "All_Alignment_Support_Breadth_Pct", "Total_Zero_Coverage_bp", "Max_Zero_Coverage_Block_bp",
        "Unresolved_Zero_Coverage_bp", "Has_High_Risk_Zero_Block_gt500bp", "Zero_Coverage_Gaps_Detail",
        "Junction_Spanning_Unique_Reads", "Junction_Informative_Unique_Reads",
        "Junction_Support_Fraction", "Junction_Support_CI_95_Lower", "Junction_Support_CI_95_Upper",
        "Junction_Clip_Unique_Reads", "Junction_Clip_Fraction", "Terminal_Anchor_Uniqueness",
        "Systematic_Internal_Clipping_Sites", "Internal_Clipping_Sites_Detail",
        "Population_SV_Polymorphisms_Count", "Population_SV_Polymorphisms_Detail",
        "Supplementary_Partner_Contig", "Supplementary_Partner_Unique_Reads",
        "Supplementary_Partner_Fraction", "Supplementary_Partner_Locus_MAD",
        "Supplementary_Partner_Depth", "Supplementary_Partner_Relative_Depth_Ratio", "Supplementary_Partner_Is_Artifact",
        "Repeat_Breadth_Pct", "Longest_Repeat_bp", "Repeat_Hit_Count",
        "CheckM2_Completeness", "CheckM2_Contamination",
        "rRNA_5S_Count", "rRNA_16S_Count", "rRNA_23S_Count", "tRNA_Canonical_20_AA_Count",
        "GTDB_Lineage",
        "geNomad_Chromosome_Score", "geNomad_Plasmid_Score", "geNomad_Virus_Score",
        "geNomad_Assignment_Derived", "geNomad_Confidence_Margin_Derived",
        "Gate1_Physical_Structural_Support", "Gate1_Notes",
        "Gate2_Coverage_Coherence",
        "Gate3_Biological_Completeness", "Gate3_Notes",
        "Gate4_Replicon_Evidence_State",
        "Illumina_Relative_Abundance_Pct", "Illumina_Catalogue_Mapped_Pct", "Illumina_Catalogue_Unmapped_Pct",
        "Illumina_Genome_Trimmed_Depth", "Illumina_Genome_Covered_Fraction", "Illumina_Aligned_Read_Count", "Illumina_ANIr_Pct",
        "Illumina_Competitive_Trimmed_Depth", "Illumina_Competitive_Depth_Ratio_C", "Illumina_Strict_Depth_Ratio_S",
        "Illumina_Rotated_Junction_Depth_100bp", "Illumina_Rotated_Junction_Local_Ratio",
        "Illumina_Rotated_Zero_Coverage_Run_bp",
        "Illumina_Rotated_Spanning_Reads_Unique", "Illumina_Rotated_Bridging_Pairs_Unique",
        "Gate5_Illumina_Orthogonal_Status", "Gate5_Notes",
        "Assembly_Usability_Grade",
        "Circularity_Confidence_Tier", "Circularity_Evidence_Type",
        "Biological_Complexity_Tier", "Graph_Complexity_Tier",
        "Alternative_Path_Contig", "Alternative_Path_Relative_Depth_Ratio", "Alternative_Path_Graph_Signal",
        "Complexity_Notes",
        "Final_Confidence_Grade"
    ]

    values = [
        args.bin, str(contig_len), str(effective_flank),
        f"{p_mean:.1f}", f"{p_median:.1f}", f"{p_mad:.1f}", f"{p_cv:.2f}",
        f"{p_breadth_1x:.2f}", f"{p_breadth_3x:.2f}",
        f"{all_support_breadth:.2f}", str(total_zero_bp), str(max_zero_block),
        str(unresolved_zero_bp), "YES" if has_high_risk else "NO", zero_gaps_detail_str,
        str(span_reads), str(inf_reads),
        f"{supp_frac:.3f}", f"{ci_lower:.3f}", f"{ci_upper:.3f}",
        str(clip_reads), f"{j_clip_frac:.3f}", anchor_status,
        str(clipping_count), clipping_detail_str,
        str(len(sv_sites_list)), sv_detail_str,
        partner_name, str(partner_reads),
        f"{partner_frac:.2f}", f"{partner_mad:.1f}",
        f"{partner_depth:.1f}", f"{rel_depth_ratio:.4f}", "YES" if is_low_abundance_artifact else "NO",
        f"{rep_breadth_pct:.2f}", str(longest_rep_bp), str(rep_hit_count),
        f"{comp:.2f}" if comp is not None else "NA",
        f"{cont:.2f}" if cont is not None else "NA",
        str(r5s), str(r16s), str(r23s), str(trna_aa_count),
        lineage,
        f"{genomad_raw['chromosome']:.3f}", f"{genomad_raw['plasmid']:.3f}", f"{genomad_raw['virus']:.3f}",
        replicon_state, f"{replicon_margin:.3f}",
        gate1_status, gate1_notes,
        gate2_status,
        gate3_status, gate3_notes,
        gate4_status,
        str(ill_metrics["Illumina_Relative_Abundance_Pct"]),
        str(ill_metrics["Illumina_Catalogue_Mapped_Pct"]),
        str(ill_metrics["Illumina_Catalogue_Unmapped_Pct"]),
        str(ill_metrics["Illumina_Genome_Trimmed_Depth"]),
        str(ill_metrics["Illumina_Genome_Covered_Fraction"]),
        str(ill_metrics["Illumina_Aligned_Read_Count"]),
        str(ill_metrics["Illumina_ANIr_Pct"]),
        str(ill_metrics["Illumina_Competitive_Trimmed_Depth"]),
        str(ill_metrics["Illumina_Competitive_Depth_Ratio_C"]),
        str(ill_metrics["Illumina_Strict_Depth_Ratio_S"]),
        str(ill_metrics["Illumina_Rotated_Junction_Depth_100bp"]),
        str(ill_metrics["Illumina_Rotated_Junction_Local_Ratio"]),
        str(ill_metrics["Illumina_Rotated_Zero_Coverage_Run_bp"]),
        str(ill_metrics["Illumina_Rotated_Spanning_Reads_Unique"]),
        str(ill_metrics["Illumina_Rotated_Bridging_Pairs_Unique"]),
        ill_metrics["Gate5_Illumina_Orthogonal_Status"],
        ill_metrics["Gate5_Notes"],
        audit_res.assembly_grade,
        audit_res.circularity.tier,
        audit_res.circularity.evidence_type,
        audit_res.complexity.bio_tier,
        audit_res.complexity.graph_tier,
        audit_res.complexity.alt_path_contig,
        f"{audit_res.complexity.rdr:.4f}",
        audit_res.complexity.rdr_signal,
        audit_res.complexity.notes,
        final_grade
    ]

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    with open(args.out_tsv, "w") as f:
        f.write("\t".join(headers) + "\n")
        f.write("\t".join(values) + "\n")

    print(f"[audit] Generated audit report for {args.bin}: Grade = {final_grade}")

if __name__ == "__main__":
    main()
