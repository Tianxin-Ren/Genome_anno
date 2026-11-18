import pandas as pd
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Aggregate multiple genome/plasmid annotation tables.")
    parser.add_argument("--bin-tables", nargs='*', help="List of master tables for bins")
    parser.add_argument("--plasmid-tables", nargs='*', help="List of master tables for plasmids")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    all_dfs = []

    # Helper to check if a file is a valid non-empty TSV
    def is_valid_table(f):
        return f and os.path.exists(f) and os.path.getsize(f) > 0

    # Process Bin Tables
    if args.bin_tables:
        for f in args.bin_tables:
            if is_valid_table(f):
                try:
                    df = pd.read_csv(f, sep='\t')
                    if df.empty:
                        continue
                    # Extract genome name from filename or path
                    genome_name = os.path.basename(f).replace("_Master_Annotation_Table.tsv", "")
                    # Insert Source information at the beginning
                    df.insert(0, 'Source_ID', genome_name)
                    df.insert(1, 'Genome_Type', 'Bin/Chromosome')
                    all_dfs.append(df)
                except Exception as e:
                    print(f"Warning: Failed to read bin table {f}: {e}", file=sys.stderr)

    # Process Plasmid Tables
    if args.plasmid_tables:
        for f in args.plasmid_tables:
            if is_valid_table(f):
                try:
                    df = pd.read_csv(f, sep='\t')
                    if df.empty:
                        continue
                    plasmid_name = os.path.basename(f).replace("_Plasmid_Annotation.tsv", "")
                    df.insert(0, 'Source_ID', plasmid_name)
                    df.insert(1, 'Genome_Type', 'Plasmid')
                    all_dfs.append(df)
                except Exception as e:
                    print(f"Warning: Failed to read plasmid table {f}: {e}", file=sys.stderr)

    if all_dfs:
        final_df = pd.concat(all_dfs, ignore_index=True)
        # Handle potential column name mismatches if schemas differ
        final_df.to_csv(args.output, sep='\t', index=False)
        print(f"Successfully aggregated {len(final_df)} genes from {len(all_dfs)} tables into {args.output}")
    else:
        print("No valid data found to aggregate. Creating empty template.", file=sys.stderr)
        # Use a more comprehensive fallback header matching merge_annotations.py
        cols = ['Source_ID', 'Genome_Type', 'Gene_ID', 'Contig', 'Start', 'End', 'Strand', 'Type', 'Gene', 'Product', 'DbXrefs']
        pd.DataFrame(columns=cols).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
