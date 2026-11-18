#!/usr/bin/env bash
# ==============================================================================
# generate_sample_list.sh
# 
# Automatically scans clean read directories (metagenomic mt_clean and/or
# metatranscriptomic mtt_clean) and generates paired sample tables (TSV / YAML)
# for abundance quantification.
#
# Usage:
#   bash scripts/generate_sample_list.sh pair <mt_clean_dir> <mtt_clean_dir> [output.tsv] [output.yaml]
#   bash scripts/generate_sample_list.sh single <reads_dir> [output.tsv] [output.yaml]
#
# Examples:
#   bash scripts/generate_sample_list.sh pair \
#       /dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/Mpf_2022/mt_clean/ \
#       /dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/Mpf_2022/mtt_clean/ \
#       /dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/Mpf_2022/samples_four_site.tsv
# ==============================================================================

set -euo pipefail

MODE="${1:-pair}"

if [ "$MODE" != "pair" ] && [ "$MODE" != "single" ]; then
    echo "Usage: $0 [pair|single] <dna_dir> [rna_dir] [output_tsv] [output_yaml]"
    exit 1
fi

DNA_DIR="${2:-}"
RNA_DIR=""
OUT_TSV=""
OUT_YAML=""

if [ "$MODE" == "pair" ]; then
    RNA_DIR="${3:-}"
    OUT_TSV="${4:-}"
    OUT_YAML="${5:-}"
    if [ -z "$DNA_DIR" ] || [ -z "$RNA_DIR" ]; then
        echo "Error: Both DNA directory and RNA directory are required for 'pair' mode."
        exit 1
    fi
else
    OUT_TSV="${3:-}"
    OUT_YAML="${4:-}"
    if [ -z "$DNA_DIR" ]; then
        echo "Error: Reads directory is required."
        exit 1
    fi
fi

# Python helper for robust filename parsing and pairing
python3 - << EOF
import os
import sys
import glob
import re

mode = "$MODE"
dna_dir = "$DNA_DIR"
rna_dir = "$RNA_DIR"
out_tsv = "$OUT_TSV"
out_yaml = "$OUT_YAML"

def scan_paired_files(directory):
    """Discovers R1 and R2 files in a directory and groups them by sample name."""
    samples = {}
    if not directory or not os.path.isdir(directory):
        return samples

    patterns = ["*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"]
    all_files = []
    for pat in patterns:
        all_files.extend(glob.glob(os.path.join(directory, pat)))

    for f in sorted(all_files):
        fname = os.path.basename(f)
        # Match R1 / fwd
        r1_match = re.search(r'([._-](?:R1|1|fwd|forward))(?:\.clean)?\.(?:fq|fastq)(?:\.gz)?$', fname, re.IGNORECASE)
        # Match R2 / rev
        r2_match = re.search(r'([._-](?:R2|2|rev|reverse))(?:\.clean)?\.(?:fq|fastq)(?:\.gz)?$', fname, re.IGNORECASE)

        if r1_match:
            stem = fname[:r1_match.start()]
            stem = re.sub(r'^(?:mRNA_)', '', stem)
            stem = re.sub(r'[._-]trim$', '', stem)
            if stem not in samples:
                samples[stem] = [None, None]
            samples[stem][0] = f
        elif r2_match:
            stem = fname[:r2_match.start()]
            stem = re.sub(r'^(?:mRNA_)', '', stem)
            stem = re.sub(r'[._-]trim$', '', stem)
            if stem not in samples:
                samples[stem] = [None, None]
            samples[stem][1] = f

    return samples

def normalize_sample_id(s):
    """Normalizes sample identifier across depth formats (e.g. C-1A-0 -> C1_0-0.5)."""
    clean = s.replace("-", "_").replace(".", "_")
    # Match standard patterns like C1, C2, W07, N30B
    clean = re.sub(r'^C_1A', 'C1', clean)
    clean = re.sub(r'^C_1', 'C1', clean)
    clean = re.sub(r'^C_2', 'C2', clean)
    return clean

dna_samples = scan_paired_files(dna_dir)
rna_samples = scan_paired_files(rna_dir) if mode == "pair" else {}

tsv_lines = ["sample_id\tdna_r1\tdna_r2\trna_r1\trna_r2"]

# Try to pair DNA and RNA samples
all_stems = sorted(list(set(dna_samples.keys()) | set(rna_samples.keys())))

matched_records = []
for stem in all_stems:
    dna_pair = dna_samples.get(stem, [None, None])
    rna_pair = rna_samples.get(stem, [None, None])

    # If stem didn't match directly, try fuzzy/normalized match
    if mode == "pair" and (None in dna_pair or None in rna_pair):
        norm_stem = normalize_sample_id(stem)
        for r_stem, r_pair in rna_samples.items():
            if normalize_sample_id(r_stem) == norm_stem:
                rna_pair = r_pair
                break

    d_r1 = dna_pair[0] or ""
    d_r2 = dna_pair[1] or ""
    r_r1 = rna_pair[0] or ""
    r_r2 = rna_pair[1] or ""

    if (d_r1 and d_r2) or (r_r1 and r_r2):
        matched_records.append((stem, d_r1, d_r2, r_r1, r_r2))
        tsv_lines.append(f"{stem}\t{d_r1}\t{d_r2}\t{r_r1}\t{r_r2}")

if out_tsv:
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\n".join(tsv_lines) + "\n")
    print(f"[generate_sample_list] Saved {len(matched_records)} paired samples to TSV: {out_tsv}")
else:
    print("\n".join(tsv_lines))

if out_yaml:
    os.makedirs(os.path.dirname(os.path.abspath(out_yaml)), exist_ok=True)
    with open(out_yaml, "w", encoding="utf-8") as f:
        f.write("# Auto-generated samples configuration\n")
        f.write("samples:\n")
        for s, d1, d2, r1, r2 in matched_records:
            f.write(f'  "{s}":\n')
            if d1 and d2:
                f.write('    short_read:\n')
                f.write(f'      external_clean_r1: "{d1}"\n')
                f.write(f'      external_clean_r2: "{d2}"\n')
            if r1 and r2:
                f.write('    transcriptome:\n')
                f.write(f'      mRNA_1: "{r1}"\n')
                f.write(f'      mRNA_2: "{r2}"\n')
    print(f"[generate_sample_list] Saved samples YAML: {out_yaml}")

EOF
