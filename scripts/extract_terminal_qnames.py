#!/usr/bin/env python3
"""
extract_terminal_qnames.py
--------------------------------------------------------------------------------
Extract terminal flanking sequences from a contig to construct a targeted
synthetic junction reference [last_flank + first_flank], and extract unique
read query names (QNAMEs) mapping to the terminal flanks from whole_assembly.bam.
"""

import os
import sys
import argparse
import subprocess

try:
    import pysam
except ImportError:
    pysam = None


def get_fasta_sequence(fa_path):
    seq_chunks = []
    if not os.path.exists(fa_path):
        return ""
    with open(fa_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_chunks.append(line.strip())
    return "".join(seq_chunks)


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


def extract_qnames_pysam(bam_path, contig_id, eff_flank, total_len):
    qnames = set()
    if not os.path.exists(bam_path) or pysam is None:
        return None
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            ref_id = resolve_bam_reference(bam.references, contig_id)
            if not ref_id or ref_id not in bam.references:
                return qnames
            # Left flank: 0 to eff_flank
            for read in bam.fetch(ref_id, 0, eff_flank):
                if not read.is_unmapped and read.query_name:
                    qnames.add(read.query_name)
            # Right flank: max(0, total_len - eff_flank) to total_len
            for read in bam.fetch(ref_id, max(0, total_len - eff_flank), total_len):
                if not read.is_unmapped and read.query_name:
                    qnames.add(read.query_name)
    except Exception as e:
        print(f"[extract_terminal_qnames] Warning: pysam fetch failed ({e}), falling back to samtools", file=sys.stderr)
        return None
    return qnames


def extract_qnames_samtools(bam_path, contig_id, eff_flank, total_len):
    qnames = set()
    if not os.path.exists(bam_path) or eff_flank <= 0:
        return qnames

    ref_id = contig_id
    cmd_idx = f"samtools idxstats '{bam_path}' | awk '{{print $1}}'"
    try:
        out_idx = subprocess.check_output(cmd_idx, shell=True, text=True, stderr=subprocess.DEVNULL)
        ref_id = resolve_bam_reference(out_idx.splitlines(), contig_id) or contig_id
    except Exception:
        pass

    # Left flank: 1 to eff_flank
    cmd1 = f"samtools view '{bam_path}' '{ref_id}:1-{eff_flank}' | cut -f1"
    try:
        out1 = subprocess.check_output(cmd1, shell=True, text=True, stderr=subprocess.DEVNULL)
        for q in out1.splitlines():
            q_str = q.strip()
            if q_str:
                qnames.add(q_str)
    except Exception:
        pass

    # Right flank: start_last to total_len
    start_last = max(1, total_len - eff_flank + 1)
    cmd2 = f"samtools view '{bam_path}' '{ref_id}:{start_last}-{total_len}' | cut -f1"
    try:
        out2 = subprocess.check_output(cmd2, shell=True, text=True, stderr=subprocess.DEVNULL)
        for q in out2.splitlines():
            q_str = q.strip()
            if q_str:
                qnames.add(q_str)
    except Exception:
        pass

    return qnames


def main():
    parser = argparse.ArgumentParser(description="Extract synthetic junction and terminal QNAMEs.")
    parser.add_argument("--contig", required=True, help="Path to contig FASTA")
    parser.add_argument("--bam", required=True, help="Path to whole assembly BAM")
    parser.add_argument("--bin", required=True, help="Contig / Bin ID")
    parser.add_argument("--flank", type=int, default=5000, help="Requested terminal flank size in bp")
    parser.add_argument("--out-junction", required=True, help="Output junction reference FASTA path")
    parser.add_argument("--out-qnames", required=True, help="Output terminal QNAMEs text file path")

    args = parser.parse_args()

    seq = get_fasta_sequence(args.contig)
    L = len(seq)
    eff_flank = min(args.flank, L // 2) if L > 0 else 0

    last_flank = seq[-eff_flank:] if eff_flank > 0 else ""
    first_flank = seq[:eff_flank] if eff_flank > 0 else ""
    junction_seq = last_flank + first_flank

    os.makedirs(os.path.dirname(os.path.abspath(args.out_junction)), exist_ok=True)
    with open(args.out_junction, "w") as f:
        header = f">{os.path.basename(args.contig)}_junction_len_{len(junction_seq)}_eff_{eff_flank}\n"
        f.write(header)
        f.write(junction_seq + "\n")

    qnames = None
    if eff_flank > 0:
        qnames = extract_qnames_pysam(args.bam, args.bin, eff_flank, L)
        if qnames is None:
            qnames = extract_qnames_samtools(args.bam, args.bin, eff_flank, L)

    if qnames is None:
        qnames = set()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_qnames)), exist_ok=True)
    with open(args.out_qnames, "w") as f:
        for q in sorted(qnames):
            f.write(q + "\n")

    print(f"[extract_terminal_qnames] {args.bin}: Length={L}bp, EffFlank={eff_flank}bp, TerminalQNAMEs={len(qnames)}")


if __name__ == "__main__":
    main()
