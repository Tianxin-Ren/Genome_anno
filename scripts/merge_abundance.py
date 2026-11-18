#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Merge abundance/ratio metrics across multiple samples")
    parser.add_argument("--type", choices=["dna", "rna", "ratio", "coverm"], required=True, help="Type of merge to perform")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input files to merge")
    parser.add_argument("--output", required=True, help="Output merged TSV path")
    return parser.parse_args()

def extract_sample_name(filepath):
    base = os.path.basename(filepath)
    for ext in [".final.tsv", ".ratio.tsv", ".coverm.tsv", ".tsv"]:
        if base.endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]

def main():
    args = parse_args()
    
    if not args.inputs:
        print("Error: No input files provided.")
        sys.exit(1)
        
    print(f"Starting merge for type: {args.type} with {len(args.inputs)} files")
    
    dfs = []
    if args.type in ["dna", "rna"]:
        metadata_cols = ["MAG", "LocalGene", "Chr", "Contig", "Start", "End", "Strand", "Length"]
        value_cols = ["Counts", "TPM_Lib", "RPKM_Lib", "TPM_Mapped", "RPKM_Mapped"]
        
        for f in args.inputs:
            sample = extract_sample_name(f)
            if not os.path.exists(f):
                print(f"Warning: File {f} not found, skipping.")
                continue
            df = pd.read_csv(f, sep="\t")
            cols_to_keep = [c for c in metadata_cols if c in df.columns] + [c for c in value_cols if c in df.columns]
            df_sub = df[cols_to_keep].copy()
            rename_dict = {col: f"{sample}_{col}" for col in value_cols if col in df_sub.columns}
            df_sub = df_sub.rename(columns=rename_dict)
            
            if not dfs:
                dfs.append(df_sub)
            else:
                drop_cols = [c for c in metadata_cols if c in df_sub.columns and c != "LocalGene"]
                dfs.append(df_sub.drop(columns=drop_cols))
                
        if dfs:
            merged_df = dfs[0]
            for df_next in dfs[1:]:
                merged_df = pd.merge(merged_df, df_next, on="LocalGene", how="outer")
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
            merged_df.to_csv(args.output, sep="\t", index=False)
            print(f"Successfully merged to {args.output}")
        else:
            print("Error: No dataframes to merge.")
            sys.exit(1)
            
    elif args.type == "ratio":
        metadata_cols = ["MAG", "LocalGene", "Chr", "Contig", "Start", "End", "Strand", "Length"]
        value_cols = ["Ratio", "DNA_Counts", "DNA_TPM_Lib", "RNA_Counts", "RNA_TPM_Lib"]
        
        for f in args.inputs:
            sample = extract_sample_name(f)
            if not os.path.exists(f):
                print(f"Warning: File {f} not found, skipping.")
                continue
            df = pd.read_csv(f, sep="\t")
            present_val_cols = [c for c in value_cols if c in df.columns]
            cols_to_keep = [c for c in metadata_cols if c in df.columns] + present_val_cols
            df_sub = df[cols_to_keep].copy()
            rename_dict = {col: f"{sample}_{col}" for col in present_val_cols}
            df_sub = df_sub.rename(columns=rename_dict)
            
            if not dfs:
                dfs.append(df_sub)
            else:
                drop_cols = [c for c in metadata_cols if c in df_sub.columns and c != "LocalGene"]
                dfs.append(df_sub.drop(columns=drop_cols))
                
        if dfs:
            merged_df = dfs[0]
            for df_next in dfs[1:]:
                merged_df = pd.merge(merged_df, df_next, on="LocalGene", how="outer")
            out_dir = os.path.dirname(args.output)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            merged_df.to_csv(args.output, sep="\t", index=False)
            print(f"Successfully merged to {args.output}")
        else:
            print("Error: No dataframes to merge.")
            sys.exit(1)

    elif args.type == "coverm":
        for f in args.inputs:
            sample = os.path.basename(f).replace(".coverm.tsv", "")
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                print(f"Warning: File {f} not found or empty, skipping.")
                continue
            df = pd.read_csv(f, sep="\t")
            if df.empty or len(df.columns) < 2:
                continue
            genome_col = df.columns[0]
            val_col = df.columns[1]
            df_sub = df[[genome_col, val_col]].copy()
            df_sub = df_sub.rename(columns={genome_col: "Genome", val_col: f"{sample}_{val_col}"})
            dfs.append(df_sub)

        if dfs:
            merged_df = dfs[0]
            for df_next in dfs[1:]:
                merged_df = pd.merge(merged_df, df_next, on="Genome", how="outer")
            out_dir = os.path.dirname(args.output)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            merged_df.to_csv(args.output, sep="\t", index=False)
            print(f"Successfully merged to {args.output}")
        else:
            print("Error: No dataframes to merge.")
            sys.exit(1)

if __name__ == "__main__":
    main()
