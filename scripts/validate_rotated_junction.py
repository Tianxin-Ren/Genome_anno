#!/usr/bin/env python3
"""
validate_rotated_junction.py
============================
Performs multi-scale Illumina short-read physical validation across circular contig junctions:
1. Mode 'rotate': Constructs a joint rotated assembly where all circular contigs are rotated at L/2,
   shifting the assembly junction L -> 1 into the midpoint, while non-circular contigs remain untouched.
2. Mode 'evaluate': Analyzes the rotated BAM to measure:
   - Ultra-narrow, medium, and wide junction depths (DJ100, DJ300, DJ1kb)
   - Local flank baseline (Dlocal) across [J-5kb:J-1kb] + [J+1kb:J+5kb]
   - Junction local depth ratio (Rlocal = DJ100 / Dlocal) and genome ratio (Rgenome)
   - Contiguous zero-coverage run (junction_zero_coverage_run_bp)
   - Exact single-read spanning count (>= 30 bp anchors on both sides)
   - Exact paired-fragment bridging count (proper pairs spanning J)
   - Both total and MAPQ >= 20 unique counts
"""

import os
import sys
import argparse
import statistics
from collections import defaultdict

try:
    import pysam
except ImportError:
    pysam = None


def load_fasta_dict(fasta_path):
    """Load FASTA into dict {seq_id: seq} preserving order."""
    sequences = {}
    current_id = None
    current_chunks = []
    if not os.path.exists(fasta_path):
        return sequences

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_chunks)
                header = line[1:].split()[0]
                current_id = header
                current_chunks = []
            else:
                current_chunks.append(line)
        if current_id:
            sequences[current_id] = "".join(current_chunks)
    return sequences


def rotate_sequence(seq, fraction=0.5):
    """
    Rotate circular sequence at fraction * len:
    M = floor(L * fraction).
    Rotated = seq[M:] + seq[:M].
    The original junction (end -> start) is located at 1-based coordinate:
    J_rot = len(seq) - M.
    """
    seq_len = len(seq)
    if seq_len <= 1:
        return seq, 0, 0
    m = int(seq_len * fraction)
    m = max(1, min(seq_len - 1, m))
    rotated_seq = seq[m:] + seq[:m]
    j_rot = seq_len - m  # 1-based coordinate of original end base
    return rotated_seq, m, j_rot


import re


def resolve_circular_id(seq_id, circular_ids):
    """
    Matches a sequence header (e.g. 'u377048ctg_len-...' or 'u377048_...')
    against a set of circular IDs (e.g. {'u377048'}).
    Returns the matched circular ID if found, else None.
    """
    if not circular_ids:
        return None
    if seq_id in circular_ids:
        return seq_id
    for cid in circular_ids:
        if seq_id.startswith(f"{cid}_") or seq_id.startswith(f"{cid}ctg_") or seq_id.startswith(f"{cid}.") or seq_id == cid:
            return cid
    for cid in circular_ids:
        tok = seq_id.split()[0].split("_")[0]
        if tok == cid:
            return cid
    for cid in circular_ids:
        # Word boundary regex check to prevent substring collisions (e.g. 'u499' matching 'u499975')
        if re.search(rf"(^|[^0-9a-zA-Z]){re.escape(cid)}([^0-9a-zA-Z]|$)", seq_id):
            return cid
    return None


def run_rotate_assembly(args):
    """Build rotated assembly and write coordinate mapping."""
    assembly = load_fasta_dict(args.assembly)
    circular_ids = set()

    if args.circular_ids:
        if os.path.isfile(args.circular_ids):
            with open(args.circular_ids) as f:
                for line in f:
                    c = line.strip().split()[0] if line.strip() else ""
                    if c:
                        if c.endswith(".fa") or c.endswith(".fasta"):
                            c = os.path.splitext(c)[0]
                        circular_ids.add(c)
        else:
            for c in args.circular_ids.split(","):
                c = c.strip()
                if c:
                    circular_ids.add(c)

    if args.circular_dir and os.path.isdir(args.circular_dir):
        for fname in os.listdir(args.circular_dir):
            if fname.endswith((".fa", ".fasta", ".fna")):
                cid = os.path.splitext(fname)[0]
                circular_ids.add(cid)

    coords_data = []
    os.makedirs(os.path.dirname(os.path.abspath(args.out_assembly)), exist_ok=True)
    if args.out_coords:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_coords)), exist_ok=True)

    with open(args.out_assembly, "w") as out_fa:
        for seq_id, seq in assembly.items():
            matched_cid = resolve_circular_id(seq_id, circular_ids)
            if matched_cid:
                rot_seq, offset, j_rot = rotate_sequence(seq, fraction=args.rotation_fraction)
                out_fa.write(f">{seq_id}\n")
                for i in range(0, len(rot_seq), 80):
                    out_fa.write(rot_seq[i:i + 80] + "\n")
                coords_data.append({
                    "contig": matched_cid,
                    "seq_id": seq_id,
                    "length": len(seq),
                    "rotation_offset": offset,
                    "rotated_junction_coord": j_rot,
                    "is_circular": "YES"
                })
            else:
                out_fa.write(f">{seq_id}\n")
                for i in range(0, len(seq), 80):
                    out_fa.write(seq[i:i + 80] + "\n")

    if args.out_coords:
        with open(args.out_coords, "w") as out_tsv:
            out_tsv.write("contig\tseq_id\tlength\trotation_offset\trotated_junction_coord\tis_circular\n")
            for row in coords_data:
                out_tsv.write(f"{row['contig']}\t{row['seq_id']}\t{row['length']}\t{row['rotation_offset']}\t{row['rotated_junction_coord']}\t{row['is_circular']}\n")

    print(f"[rotate] Wrote rotated assembly to {args.out_assembly} ({len(coords_data)} circular contigs rotated)", file=sys.stderr)


def get_per_base_coverage(bam, ref_id, start_pos, end_pos, contig_len):
    """
    Get per-base coverage in 0-based range [start_pos, end_pos) using pysam pileup or read span.
    """
    start_pos = max(0, start_pos)
    end_pos = min(contig_len, end_pos)
    region_len = max(0, end_pos - start_pos)
    if region_len == 0:
        return []

    cov = [0] * region_len
    try:
        for col in bam.pileup(ref_id, start_pos, end_pos, truncate=True, min_base_quality=0):
            idx = col.reference_pos - start_pos
            if 0 <= idx < region_len:
                cov[idx] = col.nsegments
    except Exception:
        # Fallback: compute coverage by reading fetch alignments
        for read in bam.fetch(ref_id, start_pos, end_pos):
            if read.is_unmapped:
                continue
            r_s = max(start_pos, read.reference_start)
            r_e = min(end_pos, read.reference_end or (read.reference_start + read.query_length))
            for p in range(r_s, r_e):
                cov[p - start_pos] += 1
    return cov


def evaluate_single_junction(
    bam,
    ref_id,
    contig_len,
    j_1based,
    narrow_w=100,
    med_w=300,
    wide_w=1000,
    flank_inner=1000,
    flank_outer=5000,
    min_anchor=30,
    min_mapq=20
):
    """
    Evaluates junction continuity and physical molecular spanning.
    j_1based: 1-based coordinate of the assembly junction.
    """
    j_0 = j_1based - 1  # 0-based index of junction boundary

    # 1. Narrow, Medium, Wide window depths
    narrow_cov = get_per_base_coverage(bam, ref_id, j_0 - narrow_w, j_0 + narrow_w, contig_len)
    med_cov = get_per_base_coverage(bam, ref_id, j_0 - med_w, j_0 + med_w, contig_len)
    wide_cov = get_per_base_coverage(bam, ref_id, j_0 - wide_w, j_0 + wide_w, contig_len)

    d_j100 = statistics.median(narrow_cov) if narrow_cov else 0.0
    d_j300 = statistics.median(med_cov) if med_cov else 0.0
    d_j1k = statistics.median(wide_cov) if wide_cov else 0.0

    # 2. Local Flank Baseline (Dlocal)
    left_flank_cov = get_per_base_coverage(bam, ref_id, j_0 - flank_outer, j_0 - flank_inner, contig_len)
    right_flank_cov = get_per_base_coverage(bam, ref_id, j_0 + flank_inner, j_0 + flank_outer, contig_len)
    flank_pool = left_flank_cov + right_flank_cov
    d_local = statistics.median(flank_pool) if flank_pool else d_j1k

    # 3. Whole Contig Sampled Median
    sample_step = max(1, contig_len // 1000)
    sample_points = list(range(75, contig_len - 75, sample_step))
    sampled_covs = []
    for sp in sample_points[:500]:
        c = get_per_base_coverage(bam, ref_id, sp, sp + 10, contig_len)
        if c:
            sampled_covs.extend(c)
    d_genome = statistics.median(sampled_covs) if sampled_covs else d_local

    r_local = d_j100 / max(0.1, d_local)
    r_genome = d_j100 / max(0.1, d_genome)

    # 4. Contiguous zero coverage run around J
    zero_run_bp = 0
    if wide_cov:
        j_rel = min(wide_w, len(wide_cov) // 2)
        if wide_cov[j_rel] == 0 or (j_rel > 0 and wide_cov[j_rel - 1] == 0) or (j_rel + 1 < len(wide_cov) and wide_cov[j_rel + 1] == 0):
            l_idx = j_rel
            while l_idx >= 0 and wide_cov[l_idx] == 0:
                l_idx -= 1
            r_idx = j_rel
            while r_idx < len(wide_cov) and wide_cov[r_idx] == 0:
                r_idx += 1
            zero_run_bp = (r_idx - 1) - (l_idx + 1) + 1
            zero_run_bp = max(0, zero_run_bp)

    # 5. Physical Spanning Evidence
    span_reads_all = 0
    span_reads_unique = 0

    bridging_qnames_all = set()
    bridging_qnames_unique = set()

    fetch_start = max(0, j_0 - 1500)
    fetch_end = min(contig_len, j_0 + 1500)

    for read in bam.fetch(ref_id, fetch_start, fetch_end):
        if read.is_unmapped:
            continue
        r_s = read.reference_start
        r_e = read.reference_end or (r_s + read.query_length)
        qname = read.query_name
        is_uniq = (read.mapping_quality >= min_mapq) and not (read.is_secondary or read.is_supplementary)

        # Evidence 1: Single read spans J with >= min_anchor on both sides
        if r_s <= (j_0 - min_anchor) and r_e >= (j_0 + min_anchor):
            span_reads_all += 1
            if is_uniq:
                span_reads_unique += 1

        # Evidence 2: Paired fragment bridges J
        if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped:
            m_s = read.next_reference_start
            tlen = read.template_length
            if (r_s < j_0 and m_s >= j_0) or (m_s < j_0 and r_s >= j_0):
                if 50 <= abs(tlen) <= 2000 or (abs(tlen) == 0 and abs(r_s - m_s) <= 1500):
                    bridging_qnames_all.add(qname)
                    if is_uniq:
                        bridging_qnames_unique.add(qname)

    return {
        "junction_median_depth_100bp": round(d_j100, 2),
        "junction_median_depth_300bp": round(d_j300, 2),
        "junction_median_depth_1kb": round(d_j1k, 2),
        "local_flank_median_depth": round(d_local, 2),
        "whole_contig_median_depth": round(d_genome, 2),
        "junction_local_ratio": round(r_local, 3),
        "junction_genome_ratio": round(r_genome, 3),
        "junction_zero_coverage_run_bp": zero_run_bp,
        "junction_spanning_reads": span_reads_all,
        "junction_spanning_reads_unique": span_reads_unique,
        "junction_bridging_pairs": len(bridging_qnames_all),
        "junction_bridging_pairs_unique": len(bridging_qnames_unique),
    }


def run_evaluate_junctions(args):
    """Evaluate all circular contig junctions in rotated BAM."""
    if not os.path.exists(args.rotated_bam):
        print(f"[eval_junction] Error: Rotated BAM not found: {args.rotated_bam}", file=sys.stderr)
        sys.exit(1)

    coords_map = {}
    if not os.path.exists(args.coords_tsv):
        print(f"[eval_junction] Error: Coords TSV not found: {args.coords_tsv}", file=sys.stderr)
        sys.exit(1)

    with open(args.coords_tsv) as f:
        header = [h.strip() for h in f.readline().split("\t")]
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split("\t")]
            row = dict(zip(header, parts))
            cid = row["contig"]
            coords_map[cid] = {
                "seq_id": row.get("seq_id", ""),
                "length": int(row["length"]),
                "rotation_offset": int(row["rotation_offset"]),
                "rotated_junction_coord": int(row["rotated_junction_coord"]),
            }

    results = []
    if pysam is not None:
        with pysam.AlignmentFile(args.rotated_bam, "rb") as bam:
            for cid, info in coords_map.items():
                ref_id = info.get("seq_id") or cid
                if ref_id not in bam.references:
                    matched = [r for r in bam.references if r == cid or r.startswith(f"{cid}_") or r.startswith(f"{cid}ctg_") or r.startswith(f"{cid}.") or re.search(rf"(^|[^0-9a-zA-Z]){re.escape(cid)}([^0-9a-zA-Z]|$)", r)]
                    if matched:
                        ref_id = matched[0]
                    else:
                        print(f"[eval_junction] Warning: {cid} not in BAM references", file=sys.stderr)
                        continue

                metrics = evaluate_single_junction(
                    bam,
                    ref_id,
                    info["length"],
                    info["rotated_junction_coord"],
                    narrow_w=args.narrow_window_bp,
                    med_w=args.medium_window_bp,
                    wide_w=args.wide_window_bp,
                    flank_inner=args.flank_inner_bp,
                    flank_outer=args.flank_outer_bp,
                    min_anchor=args.min_spanning_anchor_bp,
                    min_mapq=args.min_mapq
                )
                metrics["contig"] = cid
                metrics["length"] = info["length"]
                metrics["rotation_offset"] = info["rotation_offset"]
                metrics["rotated_junction_coord"] = info["rotated_junction_coord"]
                results.append(metrics)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)), exist_ok=True)
    fields = [
        "contig", "length", "rotation_offset", "rotated_junction_coord",
        "junction_median_depth_100bp", "junction_median_depth_300bp", "junction_median_depth_1kb",
        "local_flank_median_depth", "whole_contig_median_depth",
        "junction_local_ratio", "junction_genome_ratio",
        "junction_zero_coverage_run_bp",
        "junction_spanning_reads", "junction_spanning_reads_unique",
        "junction_bridging_pairs", "junction_bridging_pairs_unique"
    ]
    with open(args.out_tsv, "w") as out_f:
        out_f.write("\t".join(fields) + "\n")
        for res in results:
            vals = [str(res.get(col, "")) for col in fields]
            out_f.write("\t".join(vals) + "\n")

    print(f"[eval_junction] Successfully evaluated {len(results)} junctions -> {args.out_tsv}", file=sys.stderr)


def build_rotated_assembly(assembly_fa, circular_dir, out_assembly, out_coords, fraction=0.5):
    """Programmatic API to build rotated assembly."""
    class MockArgs:
        def __init__(self):
            self.assembly = assembly_fa
            self.circular_ids = None
            self.circular_dir = circular_dir
            self.rotation_fraction = fraction
            self.out_assembly = out_assembly
            self.out_coords = out_coords
    args = MockArgs()
    run_rotate_assembly(args)
    if os.path.exists(out_coords):
        with open(out_coords) as f:
            return max(0, sum(1 for line in f if line.strip()) - 1)
    return 0


def main():
    parser = argparse.ArgumentParser(description="Rotate circular contigs in assembly and evaluate Illumina junctions.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: rotate
    p_rot = subparsers.add_parser("rotate", help="Build rotated full assembly")
    p_rot.add_argument("--assembly", required=True, help="Path to original full ONT assembly")
    p_rot.add_argument("--circular-ids", help="File with circular contig IDs or comma-separated list")
    p_rot.add_argument("--circular-dir", "--circular-bins-dir", dest="circular_dir", help="Directory containing circular contig FASTA files")
    p_rot.add_argument("--rotation-fraction", "--fraction", dest="rotation_fraction", type=float, default=0.5, help="Fraction of length to rotate (default 0.5)")
    p_rot.add_argument("--out-assembly", required=True, help="Path to output rotated assembly FASTA")
    p_rot.add_argument("--out-coords", required=True, help="Path to output coordinates TSV")

    # Subcommand: evaluate
    p_eval = subparsers.add_parser("evaluate", help="Evaluate rotated BAM junction metrics")
    p_eval.add_argument("--rotated-bam", "--bam", dest="rotated_bam", required=True, help="Path to rotated assembly BAM")
    p_eval.add_argument("--coords-tsv", "--coords", dest="coords_tsv", required=True, help="Coordinates TSV generated by rotate")
    p_eval.add_argument("--narrow-window-bp", "--narrow-window", dest="narrow_window_bp", type=int, default=100)
    p_eval.add_argument("--medium-window-bp", "--medium-window", dest="medium_window_bp", type=int, default=300)
    p_eval.add_argument("--wide-window-bp", "--wide-window", dest="wide_window_bp", type=int, default=1000)
    p_eval.add_argument("--flank-inner-bp", "--flank-inner", dest="flank_inner_bp", type=int, default=1000)
    p_eval.add_argument("--flank-outer-bp", "--flank-outer", dest="flank_outer_bp", type=int, default=5000)
    p_eval.add_argument("--min-spanning-anchor-bp", "--min-anchor", dest="min_spanning_anchor_bp", type=int, default=30)
    p_eval.add_argument("--min-mapq", dest="min_mapq", type=int, default=20)
    p_eval.add_argument("--contig-coverage-tsv", dest="contig_coverage_tsv", default=None, help="Optional Track B coverage TSV")
    p_eval.add_argument("--out-tsv", required=True, help="Path to output junction metrics TSV")

    args = parser.parse_args()
    if args.command == "rotate":
        run_rotate_assembly(args)
    elif args.command == "evaluate":
        run_evaluate_junctions(args)


if __name__ == "__main__":
    main()
