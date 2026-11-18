#!/usr/bin/env python3
import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Generate global multi-genome summaries of target gene alignments")
    parser.add_argument("--query_faa", required=True, help="Target interest query protein sequences (FASTA)")
    parser.add_argument("--reports", nargs="+", required=True, help="List of individual genome target TSV files")
    parser.add_argument("--genomes", required=True, help="Comma-separated list of genome names corresponding to the reports")
    parser.add_argument("--out_tidy", required=True, help="Output global tidy TSV summary file")
    parser.add_argument("--out_matrix", required=True, help="Output global matrix TSV summary file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Parse active genome names
    genome_names = [g.strip() for g in args.genomes.split(",") if g.strip()]
    if len(genome_names) != len(args.reports):
        print(f"Error: Number of genomes ({len(genome_names)}) does not match number of report files ({len(args.reports)})!", file=sys.stderr)
        sys.exit(1)
        
    print(f"Compiling global summaries across {len(genome_names)} genomes...", flush=True)
    
    # 2. Extract target genes from the query FASTA file (preserving order)
    target_genes = []
    if os.path.exists(args.query_faa):
        print(f"Reading target genes from query FASTA {args.query_faa}...", flush=True)
        with open(args.query_faa) as f:
            for line in f:
                if line.startswith(">"):
                    target_id = line[1:].strip().split()[0]
                    if target_id not in target_genes:
                        target_genes.append(target_id)
    else:
        print(f"Warning: Query FASTA file {args.query_faa} does not exist.", flush=True)
        
    print(f"Total target genes of interest registered: {len(target_genes)}", flush=True)
    
    # 3. Read and store report data for each genome
    # Structure: database[genome][target_gene] = list of matched records (hits)
    database = {g: {t: [] for t in target_genes} for g in genome_names}
    
    headers_in = [
        "Query_Target_Gene", "Subject_Homolog_ID", "Identity_Pct", 
        "Alignment_Length", "Query_Coverage_Pct", "E_value", 
        "Bit_Score", "Source_Bin_Or_Unbinned", "Annotated_Gene", "Annotated_Product"
    ]
    
    for genome, report_path in zip(genome_names, args.reports):
        if not os.path.exists(report_path):
            print(f"Warning: Report for genome {genome} at {report_path} does not exist. Skipping.", flush=True)
            continue
            
        with open(report_path) as f:
            lines = f.readlines()
            if len(lines) <= 1:
                continue # Only header
            for line in lines[1:]:
                if not line.strip():
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 10:
                    continue
                
                qseq = cols[0]
                sseq = cols[1]
                pident = cols[2]
                length = cols[3]
                cov = cols[4]
                evalue = cols[5]
                bitscore = cols[6]
                source_bin = cols[7]
                gene = cols[8]
                prod = cols[9]
                
                # If target was not in query, dynamically add it
                if qseq not in target_genes:
                    target_genes.append(qseq)
                    for g in genome_names:
                        if qseq not in database[g]:
                            database[g][qseq] = []
                            
                database[genome][qseq].append({
                    "Subject_Homolog_ID": sseq,
                    "Identity_Pct": pident,
                    "Alignment_Length": length,
                    "Query_Coverage_Pct": cov,
                    "E_value": evalue,
                    "Bit_Score": bitscore,
                    "Source_Bin_Or_Unbinned": source_bin,
                    "Annotated_Gene": gene,
                    "Annotated_Product": prod
                })
                
    # 4. Generate the Tidy Long Summary Report
    print(f"Generating Tidy Summary -> {args.out_tidy}...", flush=True)
    os.makedirs(os.path.dirname(args.out_tidy), exist_ok=True)
    
    headers_out = ["Genome", "Query_Target_Gene"] + headers_in[1:]
    
    with open(args.out_tidy, "w") as f_tidy:
        f_tidy.write("\t".join(headers_out) + "\n")
        
        # Iterate systematically by Genome and Target Gene
        for genome in genome_names:
            for target in target_genes:
                hits = database[genome][target]
                if hits:
                    for h in hits:
                        row = [genome, target] + [h[col] for col in headers_in[1:]]
                        f_tidy.write("\t".join(str(x) for x in row) + "\n")
                else:
                    # Write blank fields for missing gene
                    row = [genome, target] + [""] * (len(headers_out) - 2)
                    f_tidy.write("\t".join(str(x) for x in row) + "\n")
                    
    # 5. Generate the Presence/Absence Matrix Summary Report
    print(f"Generating Presence/Absence Matrix -> {args.out_matrix}...", flush=True)
    os.makedirs(os.path.dirname(args.out_matrix), exist_ok=True)
    
    with open(args.out_matrix, "w") as f_matrix:
        # Header: Target_Gene_ID followed by genome names as columns
        f_matrix.write("\t".join(["Target_Gene_ID"] + genome_names) + "\n")
        
        for target in target_genes:
            row = [target]
            for genome in genome_names:
                hits = database[genome][target]
                if hits:
                    # Find the best hit based on Bit Score (highest is best)
                    try:
                        best_hit = max(hits, key=lambda x: float(x["Bit_Score"]))
                    except ValueError:
                        best_hit = hits[0]
                    
                    cell_val = f"{best_hit['Subject_Homolog_ID']} (Ident:{best_hit['Identity_Pct']}%, Cov:{best_hit['Query_Coverage_Pct']}%)"
                    row.append(cell_val)
                else:
                    row.append("") # Blank for missing gene
            f_matrix.write("\t".join(row) + "\n")
            
    print("Multi-genome global summary compilation completed successfully!", flush=True)

if __name__ == "__main__":
    main()
