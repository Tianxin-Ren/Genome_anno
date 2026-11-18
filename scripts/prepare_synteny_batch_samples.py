#!/usr/bin/env python3
import os
import sys
import argparse
import yaml
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch generate synteny_samples YAML configuration for target genomes from GTDB r232 database."
    )
    parser.add_argument(
        "--target_list",
        required=True,
        help="Path to file containing target genome IDs (e.g. genome_atp_synthase_modules.long.no_M00159.free_living.genome_id.tsv)"
    )
    parser.add_argument(
        "--gtdb_metadata",
        default="/dssg/home/acct-clswfp/clswfp-icdli1/Database/GTDB/archive/metadata/232.0/ar53_metadata_r232.tsv",
        help="GTDB archaeal metadata TSV (ar53_metadata_r232.tsv)"
    )
    parser.add_argument(
        "--genomes_dir",
        default="/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/results/data_prep/archaea_all/genomes",
        help="Directory containing genome nucleotide FASTAs (*.fna, *.fna.gz, *.fasta)"
    )
    parser.add_argument(
        "--proteins_dir",
        default="/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/results/mag_annotation/gtdb_r232_archaea_all/proteins/all",
        help="Directory containing protein FASTAs (*.faa, *.faa.gz)"
    )
    parser.add_argument(
        "--out_yaml",
        default="projects/config_synteny/samples.yaml",
        help="Output YAML file path (default: projects/config_synteny/samples.yaml)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional limit on number of genomes to process (0 = all)"
    )
    return parser.parse_args()

def extract_accession_keys(name):
    base = os.path.basename(name)
    for ext in [".fna.gz", ".fna", ".fa.gz", ".fa", ".fasta.gz", ".fasta", ".faa.gz", ".faa"]:
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    keys = {base}
    clean = base.replace("GB_", "").replace("RS_", "")
    keys.add(clean)
    if "_" in clean:
        parts = clean.split("_")
        if len(parts) >= 2 and parts[0] in ["GCA", "GCF"]:
            acc_part = f"{parts[0]}_{parts[1]}"
            keys.add(acc_part)
            keys.add(f"GB_{acc_part}")
            keys.add(f"RS_{acc_part}")
            keys.add(parts[1])
    return keys

def index_dir(directory, valid_exts):
    idx = {}
    if not directory or not os.path.isdir(directory):
        print(f"Warning: Directory not found: {directory}", file=sys.stderr)
        return idx
    print(f"Scanning directory: {directory} ...", flush=True)
    for entry in os.scandir(directory):
        if entry.is_file():
            fname = entry.name
            if any(fname.endswith(ext) for ext in valid_exts):
                for k in extract_accession_keys(fname):
                    if k not in idx:
                        idx[k] = entry.path
    print(f"  -> Indexed {len(idx)} unique accession keys in {directory}", flush=True)
    return idx

def find_file(acc, idx):
    for k in extract_accession_keys(acc):
        if k in idx:
            return idx[k]
    return ""

def main():
    args = parse_args()

    # 1. Read Target Genome IDs
    if not os.path.exists(args.target_list):
        print(f"Error: Target list file not found: {args.target_list}", file=sys.stderr)
        sys.exit(1)

    raw_ids = []
    with open(args.target_list, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            first_token = line.split()[0].split("\t")[0].strip()
            if first_token and first_token.lower() not in ["genome_id", "accession", "id"]:
                raw_ids.append(first_token)

    if args.limit > 0:
        raw_ids = raw_ids[:args.limit]
    print(f"Found {len(raw_ids)} target genome IDs in {args.target_list}", flush=True)

    # 2. Index Genomes and Proteomes
    genome_idx = index_dir(args.genomes_dir, (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz"))
    protein_idx = index_dir(args.proteins_dir, (".faa", ".faa.gz"))

    # 3. Read Metadata
    meta_dict = {}
    if os.path.exists(args.gtdb_metadata):
        print(f"Reading GTDB metadata: {args.gtdb_metadata} ...", flush=True)
        try:
            with open(args.gtdb_metadata, "r", encoding="utf-8", errors="ignore") as _f:
                header_cols = [c.strip() for c in _f.readline().split("\t")]
            candidate_cols = {
                "accession", "gtdb_taxonomy", "taxonomy", 
                "checkm_completeness", "completeness", 
                "checkm_contamination", "contamination"
            }
            use_cols = [c for c in header_cols if c in candidate_cols]
            if not use_cols or len(use_cols) < 2:
                use_cols = None
            m_df = pd.read_csv(args.gtdb_metadata, sep="\t", usecols=use_cols, low_memory=False)
            acc_col = "accession" if "accession" in m_df.columns else m_df.columns[0]
            tax_col = "gtdb_taxonomy" if "gtdb_taxonomy" in m_df.columns else ("taxonomy" if "taxonomy" in m_df.columns else m_df.columns[1])
            comp_col = "checkm_completeness" if "checkm_completeness" in m_df.columns else ("completeness" if "completeness" in m_df.columns else "")
            contam_col = "checkm_contamination" if "checkm_contamination" in m_df.columns else ("contamination" if "contamination" in m_df.columns else "")

            for _, r in m_df.iterrows():
                acc_val = str(r[acc_col])
                rec = {
                    "taxonomy": str(r[tax_col]),
                    "completeness": float(r[comp_col]) if comp_col and pd.notnull(r[comp_col]) else 100.0,
                    "contamination": float(r[contam_col]) if contam_col and pd.notnull(r[contam_col]) else 0.0
                }
                for k in extract_accession_keys(acc_val):
                    if k not in meta_dict:
                        meta_dict[k] = rec
            print(f"  -> Loaded metadata for {len(meta_dict)} accession variants", flush=True)
        except Exception as e:
            print(f"Warning reading metadata: {e}", file=sys.stderr)
    else:
        print(f"Notice: GTDB metadata file {args.gtdb_metadata} not found, will use fallback taxonomy.", file=sys.stderr)

    # 4. Construct samples dictionary
    samples_dict = {}
    matched_genomes = 0
    matched_proteins = 0
    matched_meta = 0

    for gid in raw_ids:
        clean_name = gid.replace("GB_", "").replace("RS_", "")
        g_fa = find_file(gid, genome_idx)
        p_faa = find_file(gid, protein_idx)

        m_info = None
        for k in extract_accession_keys(gid):
            if k in meta_dict:
                m_info = meta_dict[k]
                break

        if g_fa:
            matched_genomes += 1
        if p_faa:
            matched_proteins += 1
        if m_info:
            matched_meta += 1

        sample_entry = {
            "fasta": g_fa if g_fa else os.path.join(args.genomes_dir, f"{clean_name}.fna.gz"),
            "taxonomy": m_info["taxonomy"] if m_info else "d__Archaea;p__Unknown;c__Unknown;o__Unknown;f__Unknown;g__Unknown;s__",
        }
        if p_faa:
            sample_entry["faa"] = p_faa
        if m_info:
            sample_entry["completeness"] = m_info["completeness"]
            sample_entry["contamination"] = m_info["contamination"]

        samples_dict[clean_name] = sample_entry

    print("\n--- Summary Audit ---")
    print(f"Total Target Genomes: {len(raw_ids)}")
    print(f"  Matched Genome FASTA on disk: {matched_genomes}/{len(raw_ids)}")
    print(f"  Matched Proteome FAA on disk: {matched_proteins}/{len(raw_ids)}")
    print(f"  Matched GTDB Taxonomy in metadata: {matched_meta}/{len(raw_ids)}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_yaml)), exist_ok=True)
    out_obj = {"synteny_samples": samples_dict}
    with open(args.out_yaml, "w", encoding="utf-8") as yf:
        yaml.dump(out_obj, yf, default_flow_style=False, sort_keys=False, allow_unicode=True)

    print(f"\nSuccessfully written synteny batch configuration to: {args.out_yaml}")

if __name__ == "__main__":
    main()
