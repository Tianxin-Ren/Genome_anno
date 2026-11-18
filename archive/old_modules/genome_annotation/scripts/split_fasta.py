import os
import argparse
from Bio import SeqIO
import math

def split_fasta(input_fasta, out_dir, num_chunks):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    # First pass: count total records
    total_records = 0
    with open(input_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                total_records += 1
                
    records_per_chunk = math.ceil(total_records / num_chunks)
    print(f"Splitting {total_records} records into {num_chunks} chunks (~{records_per_chunk} per chunk)")
    
    # Second pass: split
    chunk_id = 1
    current_records = 0
    out_f = open(os.path.join(out_dir, f"chunk_{chunk_id:03d}.fasta"), "w")
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        SeqIO.write(record, out_f, "fasta")
        current_records += 1
        
        if current_records >= records_per_chunk and chunk_id < num_chunks:
            out_f.close()
            chunk_id += 1
            current_records = 0
            out_f = open(os.path.join(out_dir, f"chunk_{chunk_id:03d}.fasta"), "w")
            
    out_f.close()
    
    # Touch empty files for remaining chunks if any (to keep Snakemake happy)
    for i in range(chunk_id + 1, num_chunks + 1):
        with open(os.path.join(out_dir, f"chunk_{i:03d}.fasta"), "w") as f:
            pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("-n", "--num", type=int, default=10)
    args = parser.parse_args()
    
    split_fasta(args.input, args.outdir, args.num)

if __name__ == "__main__":
    main()
