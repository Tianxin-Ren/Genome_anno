#!/usr/bin/env python3
"""
Calculate RNA/DNA TPM-lib ratio and merge quantification results.

Usage:
    python calculate_tpm_ratio.py --dna <dna_file> --rna <rna_file> --output <output_file>

Input:
    Standard .final.tsv files from the pipeline.
    Expected columns: MAG, LocalGene, Chr, Start, End, Strand, Length, Counts, RPK, TPM_Mapped, RPKM_Mapped, RPKM_Lib, TPM_Lib

Output:
    Merged TSV with DNA_ and RNA_ prefixes for quantification columns, plus a 'Ratio' column.
    Ratio = RNA_TPM_Lib / DNA_TPM_Lib
"""

import argparse
import csv
import sys
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate RNA/DNA TPM ratio.')
    parser.add_argument('--dna', required=True, help='Path to DNA final.tsv file')
    parser.add_argument('--rna', required=True, help='Path to RNA final.tsv file')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    return parser.parse_args()

def read_tsv(filepath):
    """Read TSV file into a dictionary keyed by LocalGene."""
    data = {}
    header = []
    
    if not os.path.exists(filepath):
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        return None, None

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames
        for row in reader:
            if 'LocalGene' not in row:
                print(f"Error: 'LocalGene' column missing in {filepath}", file=sys.stderr)
                return None, None
            data[row['LocalGene']] = row
            
    return header, data

def main():
    args = parse_arguments()
    
    # 1. Read Input Files
    print(f"Reading DNA file: {args.dna}")
    dna_header, dna_data = read_tsv(args.dna)
    if not dna_data:
        sys.exit(1)
        
    print(f"Reading RNA file: {args.rna}")
    rna_header, rna_data = read_tsv(args.rna)
    if not rna_data:
        sys.exit(1)

    # 2. Define Columns
    # Metadata columns (shared)
    meta_cols = ['MAG', 'LocalGene', 'Chr', 'Start', 'End', 'Strand', 'Length']
    # Value columns to prefix
    val_cols = ['Counts', 'RPK', 'TPM_Mapped', 'RPKM_Mapped', 'RPKM_Lib', 'TPM_Lib']
    
    # Validate headers
    for col in meta_cols + val_cols:
        if col not in dna_header:
            print(f"Warning: Column '{col}' missing in DNA file. Output might be incomplete.")
        if col not in rna_header:
            print(f"Warning: Column '{col}' missing in RNA file. Output might be incomplete.")

    # 3. Prepare Output Header
    out_header = list(meta_cols)
    for col in val_cols:
        out_header.append(f"DNA_{col}")
    for col in val_cols:
        out_header.append(f"RNA_{col}")
    out_header.append("Ratio")

    # 4. Process and Merge
    output_rows = []
    
    # Iterate through all unique genes (union of keys)
    all_genes = set(dna_data.keys()) | set(rna_data.keys())
    
    print(f"Processing {len(all_genes)} genes...")
    
    count_dna_only = 0
    count_rna_only = 0
    count_shared = 0

    for gene in sorted(all_genes):
        row = {}
        
        # Get source rows
        dna_row = dna_data.get(gene)
        rna_row = rna_data.get(gene)
        
        # Fill Metadata (prefer DNA, fallback to RNA)
        source_row = dna_row if dna_row else rna_row
        for col in meta_cols:
            row[col] = source_row.get(col, "NA")

        # Fill DNA Values
        dna_tpm = 0.0
        if dna_row:
            for col in val_cols:
                val = dna_row.get(col, "0")
                row[f"DNA_{col}"] = val
                if col == "TPM_Lib":
                    try:
                        dna_tpm = float(val)
                    except ValueError:
                        dna_tpm = 0.0
        else:
            for col in val_cols:
                row[f"DNA_{col}"] = "0"
            count_rna_only += 1

        # Fill RNA Values
        rna_tpm = 0.0
        if rna_row:
            for col in val_cols:
                val = rna_row.get(col, "0")
                row[f"RNA_{col}"] = val
                if col == "TPM_Lib":
                    try:
                        rna_tpm = float(val)
                    except ValueError:
                        rna_tpm = 0.0
        else:
            for col in val_cols:
                row[f"RNA_{col}"] = "0"
            count_dna_only += 1
            
        if dna_row and rna_row:
            count_shared += 1

        # Calculate Ratio
        if dna_tpm > 0:
            ratio = rna_tpm / dna_tpm
            row["Ratio"] = f"{ratio:.6g}"
        else:
            row["Ratio"] = "NA" # Using NA for zero denominator
        
        output_rows.append(row)

    # 5. Write Output
    print(f"Writing output to: {args.output}")
    try:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=out_header, delimiter='\t')
            writer.writeheader()
            writer.writerows(output_rows)
            
        print("Done.")
        print(f"Stats: Shared={count_shared}, DNA-only={count_dna_only}, RNA-only={count_rna_only}")
            
    except Exception as e:
        print(f"Error writing output: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
