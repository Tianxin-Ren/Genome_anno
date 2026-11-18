#!/usr/bin/env python3
"""
summarize_abundance_by_gene.py

Integrates MAG annotation (PyKofamSearch / COG / master annotations) with
genomic & transcriptomic abundance quantification across multiple samples.

Primary Output:
  - Master Table (`gene_abundance_and_annotation.tsv`):
    Each row is one gene.
    Columns include:
      1. Identifiers & Coordinates: gene_id, MAG, Contig, Start, End, Strand, Length
      2. Annotations: KO, KO_Definition, COG, COG_Name, COG_Category, Annotation_Status
      3. Sample Abundances: {sample}_DNA_Counts, {sample}_DNA_TPM, {sample}_RNA_Counts, {sample}_RNA_TPM, {sample}_Ratio
      4. Summary Metrics: Mean_DNA_TPM, Mean_RNA_TPM, Mean_Ratio, Max_Ratio

Secondary Derivatives (Backward-Compatible):
  - by_mag/{sample}/{MAG}.tsv: Target gene metrics by MAG
  - by_type/{sample}/{type}.tsv: Aggregated metrics by MAG taxonomic type
  - master_by_mag.tsv & master_by_type.tsv

Supports both:
  1. Named CLI options: --merged-ratio, --annotation-tsv, --output, etc.
  2. Positional legacy invocation: python summarize_abundance_by_gene.py <base_dir> <output_dir>
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path
import pandas as pd
import numpy as np


def open_file(path, mode="rt"):
    """Transparently opens plain text or gzip files."""
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None, errors="replace" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None, errors="replace" if "t" in mode else None)


def parse_args():
    # Detect positional legacy mode vs named CLI flags
    if len(sys.argv) == 3 and not sys.argv[1].startswith("-") and not sys.argv[2].startswith("-"):
        parser = argparse.ArgumentParser(description="Legacy mode")
        parser.add_argument("base_dir", help="Directory containing abundance/ and annotation/ or TSVs")
        parser.add_argument("output_dir", help="Output directory")
        args = parser.parse_args()
        return "legacy", args

    parser = argparse.ArgumentParser(
        description="Integrate multi-sample MAG abundance metrics with mag_annotation."
    )
    parser.add_argument("--inputs", nargs="+", default=[], help="List of sample .ratio.tsv or .final.tsv files")
    parser.add_argument("--abundance-dir", default="", help="Directory containing sample .ratio.tsv files")
    parser.add_argument("--merged-ratio", default="", help="Path to merged_ratio.tsv")
    parser.add_argument("--merged-dna", default="", help="Path to merged_dna_tpm.tsv")
    parser.add_argument("--merged-rna", default="", help="Path to merged_rna_tpm.tsv")
    parser.add_argument("--annotation-tsv", default="", help="Path to mag_annotation output table")
    parser.add_argument("--annotation-dir", default="", help="Directory containing per-MAG annotation TSVs")
    parser.add_argument("--target-genes", default="", help="Optional path to target genes TSV (knum, Symbol)")
    parser.add_argument("--mag-types", default="", help="Optional path to MAG type mapping TSV (MAG, type)")
    parser.add_argument("--output", default="", help="Destination path for master table (TSV)")
    parser.add_argument("--output-dir", default="", help="Output directory for legacy/derivative summaries")
    args = parser.parse_args()
    return "named", args


def load_mag_annotations(anno_path, anno_dir=None):
    """
    Parses mag_annotation assets into a dictionary keyed by gene_id / protein_id / LocalGene:
      anno_map[gene_id] = {
          'KO': 'K03531',
          'KO_Definition': 'cell division protein FtsZ',
          'COG': 'COG0558',
          'COG_Name': '...',
          'COG_Category': '...',
          'Annotation_Status': '...'
      }
    """
    anno_map = {}
    if not anno_path and not anno_dir:
        return anno_map

    candidate_files = []
    if anno_path and os.path.exists(anno_path):
        candidate_files.append(anno_path)
    if anno_dir and os.path.isdir(anno_dir):
        for f in glob.glob(os.path.join(anno_dir, "*.tsv*")):
            candidate_files.append(f)

    for c_file in candidate_files:
        print(f"[summarize_abundance] Parsing annotation file: {c_file}")
        with open_file(c_file, "rt") as f:
            header_line = f.readline()
            if not header_line:
                continue
            headers = [h.strip() for h in header_line.split("\t")]
            col_map = {h.lower(): idx for idx, h in enumerate(headers)}

            # Case A: protein_annotation_master.tsv (Dual KEGG + COG)
            if "gene_id" in col_map and ("ko_id" in col_map or "cog_id" in col_map):
                for line in f:
                    parts = [p.strip() for p in line.split("\t")]
                    if len(parts) <= max(col_map.values()):
                        continue
                    gid = parts[col_map["gene_id"]]
                    pid = parts[col_map.get("protein_id", col_map["gene_id"])]
                    koid = parts[col_map["ko_id"]] if "ko_id" in col_map else "-"
                    kdef = parts[col_map["ko_definition"]] if "ko_definition" in col_map else "-"
                    cid = parts[col_map["cog_id"]] if "cog_id" in col_map else "-"
                    cname = parts[col_map["cog_name"]] if "cog_name" in col_map else "-"
                    ccat = parts[col_map["cog_category"]] if "cog_category" in col_map else "-"
                    status = parts[col_map["annotation_status"]] if "annotation_status" in col_map else "-"

                    info = {
                        "KO": koid if koid and koid != "-" else "Unannotated",
                        "KO_Definition": kdef if kdef and kdef != "-" else "-",
                        "COG": cid if cid and cid != "-" else "-",
                        "COG_Name": cname if cname and cname != "-" else "-",
                        "COG_Category": ccat if ccat and ccat != "-" else "-",
                        "Annotation_Status": status,
                    }
                    anno_map[gid] = info
                    anno_map[pid] = info

            # Case B: protein_KO_best.tsv or protein_KO_all_strict.tsv (PyKofamSearch)
            elif "gene_id" in col_map and "ko" in col_map:
                for line in f:
                    parts = [p.strip() for p in line.split("\t")]
                    if len(parts) <= max(col_map.values()):
                        continue
                    gid = parts[col_map["gene_id"]]
                    pid = parts[col_map.get("protein_id", col_map["gene_id"])]
                    ko = parts[col_map["ko"]]
                    kdef = parts[col_map["ko_definition"]] if "ko_definition" in col_map else "-"

                    # If already present from higher-priority best hit, keep first
                    if gid not in anno_map:
                        info = {
                            "KO": ko if ko and ko != "-" else "Unannotated",
                            "KO_Definition": kdef if kdef and kdef != "-" else "-",
                            "COG": "-",
                            "COG_Name": "-",
                            "COG_Category": "-",
                            "Annotation_Status": "ko_annotated",
                        }
                        anno_map[gid] = info
                        anno_map[pid] = info

            # Case C: Legacy format (e.g. knum, ?, LocalGene)
            elif len(headers) >= 3 and headers[0].lower() in ["knum", "ko", "koid"]:
                for line in f:
                    parts = [p.strip() for p in line.split("\t")]
                    if len(parts) >= 3:
                        ko = parts[0]
                        lg = parts[2]
                        if lg not in anno_map:
                            info = {
                                "KO": ko,
                                "KO_Definition": parts[3] if len(parts) > 3 else "-",
                                "COG": "-",
                                "COG_Name": "-",
                                "COG_Category": "-",
                                "Annotation_Status": "ko_annotated",
                            }
                            anno_map[lg] = info

    print(f"[summarize_abundance] Loaded annotations for {len(anno_map)} gene/protein entries.")
    return anno_map


def load_target_genes(target_path):
    """Loads target gene symbols if specified."""
    if not target_path or not os.path.exists(target_path):
        return {}, {}
    try:
        df = pd.read_csv(target_path, sep="\t")
        k_col = next((c for c in df.columns if c.lower() in ["knum", "ko", "koid"]), df.columns[0])
        s_col = next((c for c in df.columns if c.lower() in ["symbol", "gene", "name"]), df.columns[1] if len(df.columns) > 1 else k_col)
        knum_to_symbol = dict(zip(df[k_col].astype(str), df[s_col].astype(str)))
        target_knums = set(knum_to_symbol.keys())
        return target_knums, knum_to_symbol
    except Exception as e:
        print(f"Warning: Failed to load target genes from {target_path}: {e}", file=sys.stderr)
        return {}, {}


def load_mag_types(type_path):
    """Loads MAG -> type mapping if specified."""
    if not type_path or not os.path.exists(type_path):
        return {}
    try:
        df = pd.read_csv(type_path, sep="\t")
        m_col = next((c for c in df.columns if c.lower() in ["mag", "genome", "genome_id"]), df.columns[0])
        t_col = next((c for c in df.columns if c.lower() in ["type", "group", "clade", "lineage"]), df.columns[1] if len(df.columns) > 1 else m_col)
        return dict(zip(df[m_col].astype(str), df[t_col].astype(str)))
    except Exception as e:
        print(f"Warning: Failed to load MAG types from {type_path}: {e}", file=sys.stderr)
        return {}


def extract_sample_name(filepath):
    """Derives sample name from a file path."""
    base = os.path.basename(filepath)
    for ext in [".ratio.tsv", "_ratio.tsv", ".final.tsv", ".coverm.tsv", ".tsv"]:
        if base.endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]


def detect_sample_cols_in_header(fieldnames):
    """
    Examines fieldnames from a TSV file. If it represents a merged table
    with per-sample prefixes (e.g. C1_0-0.5_DNA_TPM_Lib), returns a mapping:
      sample_name -> {
         'dna_tpm': col_name,
         'rna_tpm': col_name,
         'dna_cnt': col_name,
         'rna_cnt': col_name,
         'ratio': col_name
      }
    Returns empty dict if no per-sample column prefixes are detected.
    """
    sample_cols = {}
    suffixes = [
        ("_DNA_TPM_Lib", "dna_tpm"), ("_DNA_TPM_Mapped", "dna_tpm"), ("_DNA_TPM", "dna_tpm"),
        ("_RNA_TPM_Lib", "rna_tpm"), ("_RNA_TPM_Mapped", "rna_tpm"), ("_RNA_TPM", "rna_tpm"),
        ("_DNA_Counts", "dna_cnt"), ("_RNA_Counts", "rna_cnt"),
        ("_Ratio", "ratio"),
        ("_TPM_Lib", "dna_tpm"), ("_TPM_Mapped", "dna_tpm"), ("_TPM", "dna_tpm"),
        ("_Counts", "dna_cnt"),
    ]
    suffixes.sort(key=lambda x: -len(x[0]))  # Longest suffix first
    for fn in fieldnames:
        for suf, m_type in suffixes:
            if fn.endswith(suf):
                prefix = fn[:-len(suf)]
                if prefix and prefix.lower() not in ["mean", "max", "merged", "local", "mag", "dna", "rna"]:
                    if prefix not in sample_cols:
                        sample_cols[prefix] = {}
                    if m_type not in sample_cols[prefix]:
                        sample_cols[prefix][m_type] = fn
                break
    return sample_cols


def to_num(val):
    if val is None or val == "" or val == "NA":
        return 0.0
    try:
        v = float(val)
        return v if not np.isnan(v) else 0.0
    except (ValueError, TypeError):
        return 0.0


def assemble_master_table(abundance_files, anno_map):
    """
    Merges abundance metrics across all sample ratio/final files and attaches
    gene annotations to build the comprehensive master table:
      Each gene is ONE ROW.
      Columns:
        1. Identification & Coordinates: gene_id, MAG, Contig, Start, End, Strand, Length
        2. Annotations: KO, KO_Definition, COG, COG_Name, COG_Category, Annotation_Status
        3. Per-Sample Metrics: {sample}_DNA_TPM, {sample}_RNA_TPM, {sample}_Ratio, {sample}_DNA_Counts, {sample}_RNA_Counts
        4. Summaries: Mean_DNA_TPM, Mean_RNA_TPM, Mean_Ratio, Max_Ratio
    """
    print(f"[summarize_abundance] Assembling master abundance table from {len(abundance_files)} file(s)...")
    genes_data = {}
    sample_names = []

    for f_path in abundance_files:
        if not os.path.exists(f_path):
            print(f"[summarize_abundance] Warning: File not found: {f_path}", file=sys.stderr)
            continue
        print(f"[summarize_abundance] Reading abundance file: {f_path}")
        with open_file(f_path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if not reader.fieldnames:
                continue

            # Identify key metadata columns
            fn_map = {fn.lower(): fn for fn in reader.fieldnames}
            gene_col = fn_map.get("localgene") or fn_map.get("gene_id") or fn_map.get("geneid") or "LocalGene"
            mag_col = fn_map.get("mag") or fn_map.get("genome_id") or "MAG"
            chr_col = fn_map.get("chr") or fn_map.get("contig") or "Chr"
            start_col = fn_map.get("start", "Start")
            end_col = fn_map.get("end", "End")
            strand_col = fn_map.get("strand", "Strand")
            len_col = fn_map.get("length", "Length")

            # Check if this file is a merged multi-sample file
            sample_col_map = detect_sample_cols_in_header(reader.fieldnames)

            if sample_col_map:
                # Case 1: Multi-sample merged file
                for s in sample_col_map.keys():
                    if s not in sample_names:
                        sample_names.append(s)

                for row in reader:
                    gid = row.get(gene_col)
                    if not gid:
                        continue
                    if gid not in genes_data:
                        raw_mag = row.get(mag_col, "")
                        raw_chr = row.get(chr_col, "")
                        if not raw_mag or raw_mag == "MAG":
                            if "___" in gid:
                                raw_mag = gid.split("___")[0]
                            elif "|" in gid:
                                raw_mag = gid.split("|")[0]
                        if not raw_chr:
                            if "___" in gid and len(gid.split("___")) > 1:
                                raw_chr = gid.split("___")[1]

                        genes_data[gid] = {
                            "gene_id": gid,
                            "MAG": raw_mag or "unknown",
                            "Contig": raw_chr or "unknown",
                            "Start": row.get(start_col, ""),
                            "End": row.get(end_col, ""),
                            "Strand": row.get(strand_col, ""),
                            "Length": row.get(len_col, ""),
                            "samples": {},
                        }

                    for s, cols in sample_col_map.items():
                        d_cnt_c = cols.get("dna_cnt")
                        d_tpm_c = cols.get("dna_tpm")
                        r_cnt_c = cols.get("rna_cnt")
                        r_tpm_c = cols.get("rna_tpm")
                        rat_c = cols.get("ratio")

                        genes_data[gid]["samples"][s] = {
                            "DNA_Counts": to_num(row.get(d_cnt_c)) if d_cnt_c else 0.0,
                            "DNA_TPM": to_num(row.get(d_tpm_c)) if d_tpm_c else 0.0,
                            "RNA_Counts": to_num(row.get(r_cnt_c)) if r_cnt_c else 0.0,
                            "RNA_TPM": to_num(row.get(r_tpm_c)) if r_tpm_c else 0.0,
                            "Ratio": row.get(rat_c, "NA") if rat_c and row.get(rat_c) not in ["", "NA"] else "NA",
                        }

            else:
                # Case 2: Single-sample file
                sample = extract_sample_name(f_path)
                if sample not in sample_names:
                    sample_names.append(sample)

                dna_cnt_col = fn_map.get("dna_counts") or fn_map.get("counts")
                dna_tpm_col = fn_map.get("dna_tpm_lib") or fn_map.get("dna_tpm_mapped") or fn_map.get("tpm_lib") or fn_map.get("tpm_mapped")
                rna_cnt_col = fn_map.get("rna_counts")
                rna_tpm_col = fn_map.get("rna_tpm_lib") or fn_map.get("rna_tpm_mapped")
                ratio_col = fn_map.get("ratio")

                for row in reader:
                    gid = row.get(gene_col)
                    if not gid:
                        continue
                    if gid not in genes_data:
                        raw_mag = row.get(mag_col, "")
                        raw_chr = row.get(chr_col, "")
                        if not raw_mag or raw_mag == "MAG":
                            if "___" in gid:
                                raw_mag = gid.split("___")[0]
                            elif "|" in gid:
                                raw_mag = gid.split("|")[0]
                        if not raw_chr:
                            if "___" in gid and len(gid.split("___")) > 1:
                                raw_chr = gid.split("___")[1]

                        genes_data[gid] = {
                            "gene_id": gid,
                            "MAG": raw_mag or "unknown",
                            "Contig": raw_chr or "unknown",
                            "Start": row.get(start_col, ""),
                            "End": row.get(end_col, ""),
                            "Strand": row.get(strand_col, ""),
                            "Length": row.get(len_col, ""),
                            "samples": {},
                        }

                    genes_data[gid]["samples"][sample] = {
                        "DNA_Counts": to_num(row.get(dna_cnt_col)) if dna_cnt_col else 0.0,
                        "DNA_TPM": to_num(row.get(dna_tpm_col)) if dna_tpm_col else 0.0,
                        "RNA_Counts": to_num(row.get(rna_cnt_col)) if rna_cnt_col else 0.0,
                        "RNA_TPM": to_num(row.get(rna_tpm_col)) if rna_tpm_col else 0.0,
                        "Ratio": row.get(ratio_col, "NA") if ratio_col and row.get(ratio_col) not in ["", "NA"] else "NA",
                    }

    # Build DataFrame
    sample_names = sorted(list(set(sample_names)))
    has_any_rna = any(
        s_data.get("RNA_Counts", 0.0) > 0.0 or s_data.get("RNA_TPM", 0.0) > 0.0
        for g in genes_data.values() for s_data in g.get("samples", {}).values()
    )
    has_any_ratio = any(
        s_data.get("Ratio") not in ["NA", None, ""]
        for g in genes_data.values() for s_data in g.get("samples", {}).values()
    )

    rows = []
    for gid, g_info in genes_data.items():
        # Match annotation
        anno = anno_map.get(gid, {})
        if not anno and "___" in gid:
            core_id = gid.split("___")[-1]
            anno = anno_map.get(core_id, {})

        row_dict = {
            "gene_id": g_info["gene_id"],
            "MAG": g_info["MAG"],
            "Contig": g_info["Contig"],
            "Start": g_info["Start"],
            "End": g_info["End"],
            "Strand": g_info["Strand"],
            "Length": g_info["Length"],
            "KO": anno.get("KO", "Unannotated"),
            "KO_Definition": anno.get("KO_Definition", "-"),
            "COG": anno.get("COG", "-"),
            "COG_Name": anno.get("COG_Name", "-"),
            "COG_Category": anno.get("COG_Category", "-"),
            "Annotation_Status": anno.get("Annotation_Status", "unannotated"),
        }

        dna_tpms = []
        rna_tpms = []
        ratios = []

        # 1. All Sample DNA TPMs (grouped together so they form a continuous block)
        for s in sample_names:
            s_data = g_info["samples"].get(s, {})
            d_tpm = s_data.get("DNA_TPM", 0.0)
            row_dict[f"{s}_DNA_TPM"] = round(d_tpm, 4)
            dna_tpms.append(d_tpm)

        # 2. All Sample RNA TPMs (if RNA present)
        if has_any_rna:
            for s in sample_names:
                s_data = g_info["samples"].get(s, {})
                r_tpm = s_data.get("RNA_TPM", 0.0)
                row_dict[f"{s}_RNA_TPM"] = round(r_tpm, 4)
                rna_tpms.append(r_tpm)

        # 3. All Sample Ratios (if Ratio present)
        if has_any_ratio:
            for s in sample_names:
                s_data = g_info["samples"].get(s, {})
                rat = s_data.get("Ratio", "NA")
                row_dict[f"{s}_Ratio"] = rat
                try:
                    rf = float(rat)
                    if not np.isnan(rf):
                        ratios.append(rf)
                except (ValueError, TypeError):
                    pass

        # 4. All Sample DNA Counts
        for s in sample_names:
            s_data = g_info["samples"].get(s, {})
            d_cnt = s_data.get("DNA_Counts", 0.0)
            row_dict[f"{s}_DNA_Counts"] = int(d_cnt) if d_cnt.is_integer() else d_cnt

        # 5. All Sample RNA Counts (if RNA present)
        if has_any_rna:
            for s in sample_names:
                s_data = g_info["samples"].get(s, {})
                r_cnt = s_data.get("RNA_Counts", 0.0)
                row_dict[f"{s}_RNA_Counts"] = int(r_cnt) if r_cnt.is_integer() else r_cnt

        # Summary statistics
        row_dict["Mean_DNA_TPM"] = round(float(np.mean(dna_tpms)), 4) if dna_tpms else 0.0
        row_dict["Mean_RNA_TPM"] = round(float(np.mean(rna_tpms)), 4) if rna_tpms else 0.0
        row_dict["Mean_Ratio"] = round(float(np.mean(ratios)), 4) if ratios else "NA"
        row_dict["Max_Ratio"] = round(float(np.max(ratios)), 4) if ratios else "NA"

        rows.append(row_dict)

    df_master = pd.DataFrame(rows)
    print(f"[summarize_abundance] Master table constructed: {len(df_master)} genes x {len(df_master.columns)} columns across {len(sample_names)} sample(s).")
    return df_master, sample_names


def main():
    mode, args = parse_args()

    if mode == "legacy":
        base_dir = Path(args.base_dir)
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        abundance_files = []
        for pat in ["abundance/*.ratio.tsv", "abundance/*_ratio.tsv", "*.ratio.tsv", "*_ratio.tsv", "merged_ratio.tsv"]:
            abundance_files.extend(list(base_dir.glob(pat)))
        abundance_files = sorted(list(set(abundance_files)))

        anno_file = None
        for cand in ["protein_annotation_master.tsv", "protein_KO_best.tsv", "protein_KO_all_strict.tsv"]:
            if (base_dir / cand).exists():
                anno_file = str(base_dir / cand)
                break
            if (base_dir / "annotation" / cand).exists():
                anno_file = str(base_dir / "annotation" / cand)
                break

        anno_dir = str(base_dir / "annotation") if (base_dir / "annotation").is_dir() else None
        target_genes = str(base_dir / "target.tsv") if (base_dir / "target.tsv").exists() else (
            str(base_dir / "all_lab.tsv") if (base_dir / "all_lab.tsv").exists() else ""
        )
        mag_types = str(base_dir / "type.tsv") if (base_dir / "type.tsv").exists() else ""
        master_out = str(output_dir / "gene_abundance_and_annotation.tsv")

    else:
        abundance_files = []
        if args.inputs:
            for inp in args.inputs:
                if os.path.exists(inp):
                    abundance_files.append(inp)
        if args.merged_ratio and os.path.exists(args.merged_ratio):
            abundance_files.append(args.merged_ratio)
        if args.merged_dna and os.path.exists(args.merged_dna):
            abundance_files.append(args.merged_dna)
        if args.merged_rna and os.path.exists(args.merged_rna):
            abundance_files.append(args.merged_rna)
        if args.abundance_dir and os.path.isdir(args.abundance_dir):
            for ext in ["*.ratio.tsv", "*_ratio.tsv", "*.final.tsv"]:
                for root, _, files in os.walk(args.abundance_dir):
                    for f in files:
                        if f.endswith((".ratio.tsv", "_ratio.tsv")) and not f.startswith("merged_"):
                            abundance_files.append(os.path.join(root, f))
        abundance_files = sorted(list(set(abundance_files)))

        anno_file = args.annotation_tsv
        anno_dir = args.annotation_dir
        target_genes = args.target_genes
        mag_types = args.mag_types
        master_out = args.output
        output_dir = Path(args.output_dir) if args.output_dir else (
            Path(os.path.dirname(master_out)) if master_out else Path(".")
        )
        if not master_out:
            master_out = str(output_dir / "gene_abundance_and_annotation.tsv")

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load annotations
    anno_map = load_mag_annotations(anno_file, anno_dir)

    # 2. Build master table (each gene one row)
    if not abundance_files:
        print("[summarize_abundance] Warning: No abundance files discovered to summarize.", file=sys.stderr)
        return

    df_master, sample_names = assemble_master_table(abundance_files, anno_map)

    # Save primary master table
    df_master.to_csv(master_out, sep="\t", index=False)
    print(f"[summarize_abundance] Master gene table successfully saved to: {master_out}")

    # Save dedicated matrix tables for convenience
    meta_sub = [c for c in ["gene_id", "MAG", "Contig", "KO", "KO_Definition", "Annotation_Status"] if c in df_master.columns]

    dna_tpm_cols = [f"{s}_DNA_TPM" for s in sample_names if f"{s}_DNA_TPM" in df_master.columns]
    if dna_tpm_cols:
        stat_cols = [c for c in ["Mean_DNA_TPM"] if c in df_master.columns]
        df_dna_tpm = df_master[meta_sub + dna_tpm_cols + stat_cols].copy()
        dna_tpm_out = output_dir / "gene_dna_tpm_matrix.tsv"
        df_dna_tpm.to_csv(dna_tpm_out, sep="\t", index=False)
        print(f"[summarize_abundance] Saved DNA TPM matrix to: {dna_tpm_out}")

    rna_tpm_cols = [f"{s}_RNA_TPM" for s in sample_names if f"{s}_RNA_TPM" in df_master.columns]
    if rna_tpm_cols:
        stat_cols = [c for c in ["Mean_RNA_TPM"] if c in df_master.columns]
        df_rna_tpm = df_master[meta_sub + rna_tpm_cols + stat_cols].copy()
        rna_tpm_out = output_dir / "gene_rna_tpm_matrix.tsv"
        df_rna_tpm.to_csv(rna_tpm_out, sep="\t", index=False)
        print(f"[summarize_abundance] Saved RNA TPM matrix to: {rna_tpm_out}")

    ratio_cols = [f"{s}_Ratio" for s in sample_names if f"{s}_Ratio" in df_master.columns]
    if ratio_cols:
        stat_cols = [c for c in ["Mean_Ratio", "Max_Ratio"] if c in df_master.columns]
        df_ratio = df_master[meta_sub + ratio_cols + stat_cols].copy()
        ratio_out = output_dir / "gene_ratio_matrix.tsv"
        df_ratio.to_csv(ratio_out, sep="\t", index=False)
        print(f"[summarize_abundance] Saved Ratio matrix to: {ratio_out}")

    # 3. Optional: Target genes filtering & legacy summaries
    target_knums, knum_to_symbol = load_target_genes(target_genes)
    type_map = load_mag_types(mag_types)

    if target_knums:
        # Filter master table for target genes
        df_target = df_master[df_master["KO"].isin(target_knums)].copy()
        if not df_target.empty:
            df_target["Target_Symbol"] = df_target["KO"].map(knum_to_symbol)
            target_out = output_dir / "master_target_genes_abundance.tsv"
            df_target.to_csv(target_out, sep="\t", index=False)
            print(f"[summarize_abundance] Saved target gene table to: {target_out} ({len(df_target)} records)")

    print("[summarize_abundance] Complete!")


if __name__ == "__main__":
    main()

