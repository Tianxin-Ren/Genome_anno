#!/usr/bin/env python3
import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Partition assembly proteome into binned and unbinned sequences")
    parser.add_argument("--assembly-faa", required=True, help="Global assembly FAA file")
    parser.add_argument("--bin-faas", nargs="*", default=[], help="List of individual bin FAA files")
    parser.add_argument("--out-binned", required=True, help="Output file for all binned protein sequences")
    parser.add_argument("--out-unbinned", required=True, help="Output file for all unbinned protein sequences")
    return parser.parse_args()

def main():
    args = parse_args()
    
    binned_headers = set()
    
    # 1. Merge all binned sequences and collect their headers
    print(f"Merging {len(args.bin_faas)} bin FAA files into {args.out_binned}...", flush=True)
    os.makedirs(os.path.dirname(args.out_binned), exist_ok=True)
    
    with open(args.out_binned, "w") as f_out:
        for faa_path in args.bin_faas:
            if not os.path.exists(faa_path):
                print(f"Warning: Bin FAA file {faa_path} does not exist. Skipping.", flush=True)
                continue
            with open(faa_path) as f_in:
                for line in f_in:
                    if line.startswith(">"):
                        # Extract the sequence ID (first word after '>')
                        seq_id = line[1:].strip().split()[0]
                        binned_headers.add(seq_id)
                    f_out.write(line)
                    
    print(f"Total unique binned protein sequences: {len(binned_headers)}", flush=True)
    
    # 2. Parse global assembly and output unbinned sequences
    print(f"Filtering global assembly {args.assembly_faa} -> {args.out_unbinned}...", flush=True)
    os.makedirs(os.path.dirname(args.out_unbinned), exist_ok=True)
    
    unbinned_count = 0
    with open(args.assembly_faa) as f_in, open(args.out_unbinned, "w") as f_out:
        keep = False
        for line in f_in:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id not in binned_headers:
                    keep = True
                    unbinned_count += 1
                    f_out.write(line)
                else:
                    keep = False
            else:
                if keep:
                    f_out.write(line)
                    
    print(f"Total unique unbinned protein sequences: {unbinned_count}", flush=True)
    print("Partitioning completed successfully!", flush=True)

if __name__ == "__main__":
    main()
