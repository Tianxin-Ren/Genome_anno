#!/usr/bin/env python3
import sys
import os
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Extract MAG/bin annotation from global assembly annotation")
    parser.add_argument("--bin_fa", required=True, help="Input bin FASTA file")
    parser.add_argument("--assembly_gff", required=True, help="Global GFF3 file")
    parser.add_argument("--assembly_faa", required=True, help="Global FAA file")
    parser.add_argument("--assembly_ffn", required=True, help="Global FFN file")
    parser.add_argument("--assembly_tsv", required=True, help="Global TSV file")
    parser.add_argument("--assembly_gbk", required=True, help="Global GBK file")
    parser.add_argument("--out_gff", required=True, help="Output GFF3 file")
    parser.add_argument("--out_faa", required=True, help="Output FAA file")
    parser.add_argument("--out_tsv", required=True, help="Output TSV file")
    parser.add_argument("--out_gbk", required=True, help="Output GBK file")
    return parser.parse_args()

def extract_contig_ids(bin_fa_path):
    contig_ids = set()
    with open(bin_fa_path) as f:
        for line in f:
            if line.startswith('>'):
                # Extract the first word of the header as the contig ID
                contig_id = line[1:].strip().split()[0]
                contig_ids.add(contig_id)
    return contig_ids

def filter_gff(assembly_gff_path, out_gff_path, contig_ids):
    gene_ids = set()
    in_fasta_section = False
    current_fasta_header = None
    current_fasta_lines = []

    with open(assembly_gff_path) as f_in, open(out_gff_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('##FASTA'):
                in_fasta_section = True
                f_out.write(line)
                continue

            if in_fasta_section:
                if line.startswith('>'):
                    if current_fasta_header and current_fasta_header in contig_ids:
                        f_out.write(">" + current_fasta_header + "\n" + "".join(current_fasta_lines))
                    current_fasta_header = line[1:].strip().split()[0]
                    current_fasta_lines = []
                else:
                    current_fasta_lines.append(line)
                continue

            # Standard GFF header/comment lines
            if line.startswith('#'):
                if line.startswith('##sequence-region'):
                    parts = line.strip().split()
                    if len(parts) >= 2 and parts[1] in contig_ids:
                        f_out.write(line)
                elif not line.startswith('##sequence-region'):
                    f_out.write(line)
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            seqid = parts[0]
            if seqid in contig_ids:
                f_out.write(line)
                # Parse gene/protein IDs from attributes column
                attributes = parts[8]
                id_match = re.search(r'ID=([^; \t\n]+)', attributes)
                if id_match:
                    gene_ids.add(id_match.group(1))
                lt_match = re.search(r'locus_tag=([^; \t\n]+)', attributes)
                if lt_match:
                    gene_ids.add(lt_match.group(1))

        # Handle the very last FASTA record
        if in_fasta_section and current_fasta_header and current_fasta_header in contig_ids:
            f_out.write(">" + current_fasta_header + "\n" + "".join(current_fasta_lines))

    return gene_ids

def filter_fasta(in_path, out_path, allowed_ids):
    with open(in_path) as f_in, open(out_path, 'w') as f_out:
        keep = False
        for line in f_in:
            if line.startswith('>'):
                header_id = line[1:].strip().split()[0]
                if header_id in allowed_ids:
                    keep = True
                    f_out.write(line)
                else:
                    keep = False
            else:
                if keep:
                    f_out.write(line)

def filter_tsv(in_path, out_path, contig_ids, gene_ids):
    with open(in_path) as f_in, open(out_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                continue
            parts = line.strip().split('\t')
            if not parts:
                continue
            
            # Check if first column (Sequence Id) is in contig_ids
            # Or if there is a locus tag in the columns that matches gene_ids
            is_match = False
            if parts[0] in contig_ids:
                is_match = True
            else:
                for part in parts:
                    if part in gene_ids:
                        is_match = True
                        break
            if is_match:
                f_out.write(line)

def filter_gbk(in_path, out_path, contig_ids):
    if not os.path.exists(in_path):
        return
    
    with open(in_path) as f_in, open(out_path, 'w') as f_out:
        current_record = []
        for line in f_in:
            current_record.append(line)
            if line.startswith('//'):
                record_str = "".join(current_record).strip()
                current_record = []  # Free memory immediately
                if not record_str:
                    continue
                
                lines = record_str.split('\n')
                locus_line = lines[0]
                if not locus_line.startswith('LOCUS'):
                    continue
                
                parts = locus_line.split()
                if len(parts) < 2:
                    continue
                locus_id = parts[1]
                
                matched = False
                if locus_id in contig_ids:
                    matched = True
                else:
                    # Check for 16-character truncation or general prefix matching
                    for cid in contig_ids:
                        if cid.startswith(locus_id) or locus_id.startswith(cid[:16]):
                            matched = True
                            break
                
                # Double check ACCESSION / VERSION lines in the record
                if not matched:
                    for r_line in lines:
                        if r_line.startswith('ACCESSION') or r_line.startswith('VERSION'):
                            acc_parts = r_line.split()
                            if len(acc_parts) >= 2 and acc_parts[1] in contig_ids:
                                matched = True
                                break
                
                if matched:
                    f_out.write(record_str + '\n//\n')

def main():
    args = parse_args()
    print(f"Loading contig IDs from {args.bin_fa}...", flush=True)
    contig_ids = extract_contig_ids(args.bin_fa)
    print(f"Found {len(contig_ids)} contig IDs in bin.", flush=True)

    print(f"Filtering GFF3: {args.assembly_gff} -> {args.out_gff}...", flush=True)
    os.makedirs(os.path.dirname(args.out_gff), exist_ok=True)
    gene_ids = filter_gff(args.assembly_gff, args.out_gff, contig_ids)
    print(f"Found {len(gene_ids)} annotated gene IDs in bin contigs.", flush=True)

    print(f"Filtering FAA: {args.assembly_faa} -> {args.out_faa}...", flush=True)
    filter_fasta(args.assembly_faa, args.out_faa, gene_ids)

    print(f"Filtering TSV: {args.assembly_tsv} -> {args.out_tsv}...", flush=True)
    filter_tsv(args.assembly_tsv, args.out_tsv, contig_ids, gene_ids)

    print(f"Filtering GBK: {args.assembly_gbk} -> {args.out_gbk}...", flush=True)
    filter_gbk(args.assembly_gbk, args.out_gbk, contig_ids)

    print("Subsetting completed successfully!", flush=True)

if __name__ == "__main__":
    main()
