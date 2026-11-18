#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
import os
import re

def main():
    parser = argparse.ArgumentParser(description="Filter GTDB metadata and extract accessions in chunks")
    parser.add_argument("-i", "--input", required=True, help="Path to GTDB metadata TSV file")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory to save accession chunks")
    parser.add_argument("-q", "--query", required=False, help="Pandas query string", default="")
    parser.add_argument("--accessions", required=False, help="Optional file containing accessions to retain")
    parser.add_argument("-t", "--output-tsv", required=False, help="Path to save the filtered metadata TSV file", default="")
    parser.add_argument("--chunk-size", type=int, default=100, help="Number of accessions per chunk file (default: 100)")
    args = parser.parse_args()

    try:
        # Load the metadata
        print(f"Loading GTDB metadata from: {args.input}", file=sys.stderr)
        # GTDB metadata might be mixed types, set low_memory=False
        df = pd.read_csv(args.input, sep='\t', low_memory=False)
        print(f"Loaded {len(df)} records.", file=sys.stderr)
        
        # Check if accession column exists
        if 'accession' not in df.columns:
            print("Error: 'accession' column not found in the input file.", file=sys.stderr)
            sys.exit(1)

        # Apply metadata filter if provided
        if args.query and args.query.strip():
            print(f"Applying filter query: {args.query}", file=sys.stderr)
            try:
                filtered_df = df.query(args.query, engine='python')
            except Exception as e:
                print(f"Error evaluating query '{args.query}': {e}", file=sys.stderr)
                sys.exit(1)
            print(f"Filtered down to {len(filtered_df)} records.", file=sys.stderr)
        else:
            print("No filter query provided. Extracting all records.", file=sys.stderr)
            filtered_df = df

        # Optionally restrict the metadata to an external accession list.
        if args.accessions:
            requested = set()
            with open(args.accessions) as handle:
                for line in handle:
                    accession = line.strip().split("\t", 1)[0]
                    accession = re.sub(r"^(RS_|GB_)", "", accession)
                    if re.fullmatch(r"GC[AF]_\d+(\.\d+)?", accession):
                        requested.add(accession)

            metadata_accessions = (
                filtered_df["accession"].astype(str)
                .str.replace(r"^(RS_|GB_)", "", regex=True)
            )
            filtered_df = filtered_df[metadata_accessions.isin(requested)]
            print(
                f"Matched {len(filtered_df)} metadata rows from {len(requested)} requested accessions.",
                file=sys.stderr,
            )

        # Save the filtered metadata TSV if requested
        if args.output_tsv:
            print(f"Saving filtered metadata to: {args.output_tsv}", file=sys.stderr)
            # Ensure the directory for output tsv exists
            tsv_dir = os.path.dirname(os.path.abspath(args.output_tsv))
            os.makedirs(tsv_dir, exist_ok=True)
            filtered_df.to_csv(args.output_tsv, sep='\t', index=False)
            print(f"Saved {len(filtered_df)} rows to {args.output_tsv}", file=sys.stderr)

        # Extract and clean accessions
        accessions = filtered_df['accession'].dropna().astype(str).tolist()
        
        # GTDB accessions are often prefixed with 'RS_' or 'GB_'
        cleaned_accessions = []
        for acc in accessions:
            if acc.startswith('RS_'):
                cleaned_accessions.append(acc[3:])
            elif acc.startswith('GB_'):
                cleaned_accessions.append(acc[3:])
            else:
                cleaned_accessions.append(acc)
        cleaned_accessions = list(dict.fromkeys(cleaned_accessions))

        # Ensure output directory exists
        os.makedirs(args.output_dir, exist_ok=True)

        if not cleaned_accessions:
            print("No accessions matched the filter criteria.", file=sys.stderr)
            # Create an empty chunk file so Snakemake sees an output
            with open(os.path.join(args.output_dir, "chunk_000.txt"), 'w') as f:
                pass
            sys.exit(0)

        # Write to chunk files
        chunk_size = args.chunk_size if args.chunk_size > 0 else len(cleaned_accessions)
        num_chunks = (len(cleaned_accessions) + chunk_size - 1) // chunk_size

        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min(start_idx + chunk_size, len(cleaned_accessions))
            chunk_accs = cleaned_accessions[start_idx:end_idx]
            
            chunk_file = os.path.join(args.output_dir, f"chunk_{i+1:03d}.txt")
            with open(chunk_file, 'w') as f:
                for acc in chunk_accs:
                    f.write(f"{acc}\n")
            
            print(f"Wrote {len(chunk_accs)} accessions to {chunk_file}", file=sys.stderr)

        print(f"Successfully split {len(cleaned_accessions)} accessions into {num_chunks} chunks inside {args.output_dir}", file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
