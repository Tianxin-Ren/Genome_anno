#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import gzip
import re
import statistics
import shlex
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Stage 3: Physical locus alignment, structural event classification, quantitative read junction verification, OFI physical metrics and locus_state determination.")
    parser.add_argument("--target_name", required=True, help="Target genome sample ID")
    parser.add_argument("--target_locus_fa", required=True, help="Target locus DNA FASTA")
    parser.add_argument("--ref_locus_fa", default="", help="Nearest confirmed positive reference locus DNA FASTA")
    parser.add_argument("--refs_locus_dir", default="", help="Stage 1 reference locus directory")
    parser.add_argument("--topology_tsv", required=True, help="Stage 1 anchor topology TSV")
    parser.add_argument("--cluster_tsv", required=True, help="Stage 2 cluster rescue summary TSV")
    parser.add_argument("--reads_bam", default="", help="Optional BAM file of reads aligned to target")
    parser.add_argument("--minimap2_env", default="", help="Conda environment for minimap2 and samtools")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPU threads")
    parser.add_argument("--module_id", default="M00159", help="ATP synthase module ID")
    parser.add_argument("--module_name", default="", help="ATP synthase module human-readable name")
    parser.add_argument("--out_paf", required=True, help="Output PAF alignment file from minimap2")
    parser.add_argument("--out_locus_tsv", required=True, help="Output TSV with physical locus metrics and locus_state")
    return parser.parse_args()

def run_cmd(cmd):
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return res

def safe_int(val, default=0):
    if pd.isna(val) or val is None or val == "":
        return default
    try:
        return int(float(val))
    except Exception:
        return default

def safe_str(val, default=""):
    if pd.isna(val) or val is None:
        return default
    return str(val)

def read_fasta_seqs(fa_path):
    """Reads fasta sequence without requiring biopython."""
    if not os.path.exists(fa_path) or os.path.getsize(fa_path) == 0:
        return ""
    seq_chunks = []
    opener = gzip.open if fa_path.endswith(".gz") else open
    with opener(fa_path, "rt", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq_chunks.append(line)
    return "".join(seq_chunks)

def parse_paf_alignments(paf_file):
    """Parses PAF file to check continuous anchor alignment and deletion gaps."""
    records = []
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return records
    with open(paf_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 12:
                try:
                    records.append({
                        "qname": parts[0],
                        "qlen": int(parts[1]),
                        "qstart": int(parts[2]),
                        "qend": int(parts[3]),
                        "strand": parts[4],
                        "tname": parts[5],
                        "tlen": int(parts[6]),
                        "tstart": int(parts[7]),
                        "tend": int(parts[8]),
                        "matches": int(parts[9]),
                        "aln_len": int(parts[10]),
                        "mapq": int(parts[11]),
                        "cigar": next((tag[5:] for tag in parts[12:] if tag.startswith("cg:Z:")), "")
                    })
                except ValueError:
                    pass
    return records

def find_deletion_junction(records, minimum_deletion=1000, locus_start=None, locus_end=None):
    """Find a deletion in the query genome relative to the reference locus."""
    candidates = []
    def inside_locus(breakpoint):
        return locus_start is None or locus_end is None or locus_start <= breakpoint <= locus_end
    records = sorted((r for r in records if r["aln_len"] >= 500 and r["mapq"] >= 20),
                     key=lambda r: (r["qname"], r["qstart"]))
    for record in records:
        cigar = record.get("cigar", "")
        if not cigar:
            continue
        qpos = record["qstart"] if record["strand"] == "+" else record["qend"]
        for count_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            count = int(count_text)
            if op in "M=XI":
                qpos += count if record["strand"] == "+" else -count
            elif op == "D" and count >= minimum_deletion:
                left_aligned = qpos - record["qstart"] if record["strand"] == "+" else record["qend"] - qpos
                right_aligned = record["qend"] - qpos if record["strand"] == "+" else qpos - record["qstart"]
                if left_aligned >= 500 and right_aligned >= 500 and inside_locus(qpos):
                    candidates.append((count, 0, qpos))
    for left, right in zip(records, records[1:]):
        if left["qname"] != right["qname"] or left["tname"] != right["tname"] or left["strand"] != right["strand"]:
            continue
        query_gap = right["qstart"] - left["qend"]
        if left["strand"] == "+":
            ref_gap = right["tstart"] - left["tend"]
        else:
            ref_gap = left["tstart"] - right["tend"]
        breakpoint = (left["qend"] + right["qstart"]) // 2
        if ref_gap >= 0 and query_gap >= -50 and ref_gap - query_gap >= minimum_deletion and inside_locus(breakpoint):
            candidates.append((ref_gap - query_gap, query_gap, breakpoint))
    return max(candidates) if candidates else None

def read_aligns_across_breakpoint(start_one_based, cigar, breakpoint_zero_based, flank_bases=20):
    """Require aligned read bases on both sides in one uninterrupted reference block."""
    ref_pos = start_one_based - 1
    block_start = None
    block_end = None
    for count_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        count = int(count_text)
        if op in "M=X":
            if block_start is None:
                block_start = ref_pos
            ref_pos += count
            block_end = ref_pos
        elif op in "DN":
            if block_start is not None and block_start <= breakpoint_zero_based - flank_bases and block_end >= breakpoint_zero_based + flank_bases:
                return True
            block_start = None
            block_end = None
            ref_pos += count
        elif op in "SHP":
            block_start = None
            block_end = None
        # Insertions consume read bases but do not interrupt reference coverage.
    return block_start is not None and block_start <= breakpoint_zero_based - flank_bases and block_end >= breakpoint_zero_based + flank_bases

def measure_junction_reads(bam_path, contig, breakpoint, slice_start, env_prefix):
    """Count observed alignments crossing the breakpoint in original coordinates."""
    genomic_pos = slice_start + breakpoint + 1
    region = f"{contig}:{max(1, genomic_pos - 250)}-{genomic_pos + 250}"
    cmd = f"{env_prefix}samtools view -q 30 {shlex.quote(bam_path)} {shlex.quote(region)}"
    result = run_cmd(cmd)
    if result.returncode != 0:
        return "AMBIGUOUS_MAPPING", 0, 0.0, "NA"
    spanning = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6 or fields[2] != contig:
            continue
        flag = int(fields[1])
        if flag & (0x4 | 0x100 | 0x400 | 0x800):
            continue
        if read_aligns_across_breakpoint(int(fields[3]), fields[5], genomic_pos - 1):
            spanning.append(int(fields[4]))
    depth_result = run_cmd(f"{env_prefix}samtools depth -a -r {shlex.quote(region)} {shlex.quote(bam_path)}")
    depth = {}
    if depth_result.returncode == 0:
        for line in depth_result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) >= 3 and fields[1].isdigit() and fields[2].isdigit():
                depth[int(fields[1])] = int(fields[2])
    near = [value for pos, value in depth.items() if abs(pos - genomic_pos) <= 50]
    flank = [value for pos, value in depth.items() if 150 <= abs(pos - genomic_pos) <= 250]
    coverage_ratio = round(statistics.mean(near) / statistics.mean(flank), 3) if near and flank and statistics.mean(flank) else "NA"
    status = "JUNCTION_SUPPORTED" if len(spanning) >= 3 else "JUNCTION_NOT_SUPPORTED"
    return status, len(spanning), statistics.mean(spanning) if spanning else 0.0, coverage_ratio

def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out_paf), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_locus_tsv), exist_ok=True)

    cmd_prefix = f"conda run -p {args.minimap2_env} " if args.minimap2_env and os.path.isabs(args.minimap2_env) else (f"conda run -n {args.minimap2_env} " if args.minimap2_env else "")

    # 1. Load upstream metrics
    topo = {}
    if os.path.exists(args.topology_tsv) and os.path.getsize(args.topology_tsv) > 0:
        topo_df = pd.read_csv(args.topology_tsv, sep="\t")
        topo = topo_df.iloc[0].to_dict() if not topo_df.empty else {}

    cluster = {}
    if os.path.exists(args.cluster_tsv) and os.path.getsize(args.cluster_tsv) > 0:
        cluster_df = pd.read_csv(args.cluster_tsv, sep="\t")
        cluster = cluster_df.iloc[0].to_dict() if not cluster_df.empty else {}

    # Select the reference named by Stage 1; never align the target to itself.
    ref_fa = args.ref_locus_fa
    if not ref_fa and args.refs_locus_dir:
        ref_acc = safe_str(topo.get("nearest_positive_ref", ""))
        candidate = os.path.join(args.refs_locus_dir, f"{ref_acc}_locus.fa")
        if ref_acc not in ("", "NONE") and os.path.isfile(candidate):
            ref_fa = candidate
    with open(args.out_paf, "w"):
        pass
    if ref_fa and os.path.isfile(ref_fa) and os.path.getsize(ref_fa) and os.path.getsize(args.target_locus_fa):
        cmd_mm = f"{cmd_prefix}minimap2 -cx asm5 {shlex.quote(ref_fa)} {shlex.quote(args.target_locus_fa)} -o {shlex.quote(args.out_paf)}"
        result = run_cmd(cmd_mm)
        if result.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {result.stderr}")

    paf_records = parse_paf_alignments(args.out_paf)
    has_reference_locus = bool(ref_fa and os.path.isfile(ref_fa) and os.path.getsize(ref_fa))

    # 3. Analyze structural event type
    nucleotide_deletion_support = "NO"
    structural_event_type = "UNRESOLVED"

    observed_gap = safe_int(topo.get("observed_anchor_gap", -1), -1)
    expected_len = safe_int(topo.get("expected_locus_length", 15000), 15000)
    same_contig = safe_str(topo.get("same_contig", "NO"))
    syntenic_geometry = safe_str(topo.get("syntenic_geometry", "BILATERAL_SAME_CONTIG"))
    gate1_status = safe_str(topo.get("gate1_status", "ASSEMBLY_UNRESOLVED"))

    # OFI Quantitative Metrics:
    # B_bridge: 1.0 if bilateral collinear anchors reside on same unbroken contig
    if same_contig == "YES" and (syntenic_geometry == "BILATERAL_SAME_CONTIG" or gate1_status == "PASS"):
        B_bridge = 1.0
    else:
        B_bridge = 0.0

    # Q_continuity: 1.0 - gap_fraction (assembly Ns in target locus sequence)
    target_dna = read_fasta_seqs(args.target_locus_fa)
    total_bp = len(target_dna)
    if total_bp > 0:
        n_count = target_dna.upper().count("N")
        gap_frac = n_count / total_bp
        Q_continuity = round(max(0.0, min(1.0, 1.0 - gap_frac)), 4)
    else:
        Q_continuity = 0.0

    # D_closure: Physical locus compression against dereplicated homologous reference span
    obs_span_raw = topo.get("observed_interanchor_span", None)
    if obs_span_raw is None or pd.isna(obs_span_raw) or safe_int(obs_span_raw, -1) < 0:
        obs_span = observed_gap if observed_gap >= 0 else 0
    else:
        obs_span = safe_int(obs_span_raw, 0)

    exp_span_raw = topo.get("expected_interanchor_span", None)
    if exp_span_raw is None or pd.isna(exp_span_raw) or safe_int(exp_span_raw, 0) <= 0:
        exp_span = expected_len if expected_len > 0 else 15000
    else:
        exp_span = safe_int(exp_span_raw, 15000)

    if exp_span > 0:
        D_closure = round(max(0.0, min(1.0, 1.0 - (obs_span / exp_span))), 4)
    else:
        D_closure = 0.0

    c_state = cluster.get("cluster_state", "FULLY_ABSENT")

    slice_start = safe_int(topo.get("target_slice_start", 0))
    interval_start = safe_int(topo.get("ancestral_interval_start", -1), -1)
    interval_end = safe_int(topo.get("ancestral_interval_end", -1), -1)
    junction = find_deletion_junction(
        paf_records,
        locus_start=interval_start - slice_start if interval_start >= 0 and interval_end > interval_start else None,
        locus_end=interval_end - slice_start if interval_start >= 0 and interval_end > interval_start else None)
    if junction:
        _, inserted_bp, _ = junction
        structural_event_type = "CLEAN_DELETION" if inserted_bp <= 200 else "DELETION_WITH_INSERTION"
        nucleotide_deletion_support = "YES"
    elif paf_records and c_state in ["COMPLETE", "COMPLETE_AFTER_RESCUE"]:
        structural_event_type = "LOCUS_INTACT"
    elif paf_records and c_state == "FULLY_ABSENT":
        structural_event_type = "POSSIBLE_LOCUS_REPLACEMENT"

    # 4. Quantitative Read-Level Junction Support
    junction_support_status = "NO_READS"
    junction_spanning_reads = 0
    junction_mean_mapq = 0.0
    junction_coverage_ratio = "NA"

    if args.reads_bam and os.path.exists(args.reads_bam) and os.path.getsize(args.reads_bam) > 0 and junction:
        junction_support_status, junction_spanning_reads, junction_mean_mapq, junction_coverage_ratio = measure_junction_reads(
            args.reads_bam, safe_str(topo.get("target_contig", "")), junction[2],
            safe_int(topo.get("target_slice_start", 0)), cmd_prefix)

    # 5. Pure Physical locus_state Determination
    if gate1_status == "ASSEMBLY_UNRESOLVED" or same_contig == "NO":
        locus_state = "ASSEMBLY_UNRESOLVED"
    elif c_state in ["COMPLETE", "COMPLETE_AFTER_RESCUE"]:
        locus_state = "PRESENT_COMPLETE"
    elif c_state in ["PARTIAL_INTACT", "PARTIAL_PSEUDOGENIZED", "PARTIAL_ABSENT", "MIXED"]:
        locus_state = "PRESENT_PARTIAL"
    elif c_state in ("ANCESTRAL_ABSENT", "UNRESOLVED_NO_QUERY"):
        locus_state = "UNRESOLVED_NO_QUERY"
        structural_event_type = "UNRESOLVED"
        nucleotide_deletion_support = "NO"
    elif c_state == "FULLY_ABSENT":
        if gate1_status == "PASS" and nucleotide_deletion_support == "YES" and (obs_span < exp_span * 0.5) and (junction_support_status in ("NO_READS", "JUNCTION_SUPPORTED")):
            locus_state = "SYNTENIC_ABSENT_HIGH_CONF"
        elif not has_reference_locus:
            locus_state = "UNRESOLVED_NO_REFERENCE"
        else:
            locus_state = "SYNTENIC_ABSENT_MEDIUM_CONF"
    else:
        locus_state = "SYNTENIC_ABSENT_MEDIUM_CONF"

    # Write output TSV
    out_row = {
        "target_id": args.target_name,
        "structural_event_type": structural_event_type,
        "nucleotide_deletion_support": nucleotide_deletion_support,
        "junction_support_status": junction_support_status,
        "junction_spanning_reads": junction_spanning_reads,
        "junction_mean_mapq": junction_mean_mapq,
        "junction_coverage_ratio": junction_coverage_ratio,
        "locus_state": locus_state,
        "B_bridge": B_bridge,
        "Q_continuity": Q_continuity,
        "D_closure": D_closure,
        "observed_interanchor_span": obs_span,
        "expected_interanchor_span": exp_span
    }

    pd.DataFrame([out_row]).to_csv(args.out_locus_tsv, sep="\t", index=False)
    print(f"Stage 3 Physical Locus Validation Completed: locus_state={locus_state}, B_bridge={B_bridge}, Q_continuity={Q_continuity}, D_closure={D_closure}, structural_event={structural_event_type}", flush=True)

if __name__ == "__main__":
    main()
