#!/usr/bin/env python3
import os
import sys
import re
import argparse

def get_ko(val):
    if not isinstance(val, str):
        return None
    m = re.search(r'(K\d{5})', val)
    return m.group(1) if m else None

def load_gene_kos(file_path):
    genes_order = []
    ko_to_genes = {}
    
    if not os.path.exists(file_path):
        print(f"Error: Gene TSV file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
        
    with open(file_path, 'r') as f:
        # Filter empty lines
        lines = [line.strip().split('\t') for line in f if line.strip()]
        
    if not lines:
        print(f"Error: Gene TSV file is empty", file=sys.stderr)
        sys.exit(1)
        
    # Check if first line is a header
    # If the first line contains a valid KO (like K00436), it is NOT a header.
    first_row = lines[0]
    has_ko_in_first_row = any(re.search(r'K\d{5}', val) for val in first_row)
    
    header = None
    data_rows = lines
    if not has_ko_in_first_row:
        # First row is likely a header
        header = [h.lower().strip() for h in first_row]
        data_rows = lines[1:]
        
    # Determine indices for Gene Name and KO (default to 0 and 1)
    gene_idx = 0
    ko_idx = 1
    
    if header:
        # Try to find indices by keyword matching
        found_gene = False
        found_ko = False
        for i, h in enumerate(header):
            if any(k in h for k in ['gene', 'name', 'id', 'target', 'alias']):
                gene_idx = i
                found_gene = True
            elif any(k in h for k in ['ko', 'kegg', 'kofam', 'knum']):
                ko_idx = i
                found_ko = True
        
        # Fallback if headers exist but keywords not matched
        if not found_gene and found_ko:
            gene_idx = 1 - ko_idx if len(header) == 2 else (0 if ko_idx != 0 else 1)
        elif not found_ko and found_gene:
            ko_idx = 1 - gene_idx if len(header) == 2 else (1 if gene_idx != 1 else 0)
    else:
        # No header, auto-detect by scanning the first row
        for i, val in enumerate(first_row):
            if re.search(r'K\d{5}', val):
                ko_idx = i
        # The other index is the gene name
        gene_idx = 1 - ko_idx if len(first_row) == 2 else (0 if ko_idx != 0 else 1)
        
    print(f"Gene TSV auto-detected: Gene Name column index = {gene_idx}, KO column index = {ko_idx}")
    
    seen_genes = set()
    for row in data_rows:
        if len(row) <= max(gene_idx, ko_idx):
            continue
        gene_name = row[gene_idx].strip()
        ko_raw = row[ko_idx].strip()
        
        # Extract clean KO
        ko = get_ko(ko_raw)
        if not ko:
            continue
        
        if gene_name not in seen_genes:
            seen_genes.add(gene_name)
            genes_order.append(gene_name)
            
        ko_to_genes.setdefault(ko, []).append(gene_name)
        
    print(f"Loaded {len(genes_order)} unique genes and {len(ko_to_genes)} KO mapping relationships from {file_path}.")
    return genes_order, ko_to_genes

def process_summary(summary_path, ko_to_genes):
    genome_counts = {}
    
    if not os.path.exists(summary_path):
        print(f"Error: Summary file not found at {summary_path}", file=sys.stderr)
        sys.exit(1)
        
    with open(summary_path, 'r') as f:
        header_line = f.readline()
        header = [h.strip() for h in header_line.split('\t')]
        
        # Determine column indices (Default: Name=1st column, Evalue_Annotation=10th column)
        name_idx = 0
        ko_idx = 9  # default 10th column (Index 9)
        
        if 'Name' in header:
            name_idx = header.index('Name')
        if 'Evalue_Annotation' in header:
            ko_idx = header.index('Evalue_Annotation')
        elif 'Score_Annotation' in header:
            ko_idx = header.index('Score_Annotation')
            
        print(f"Summary file auto-detected: Genome Name column index = {name_idx}, KO column index = {ko_idx}")
        
        row_count = 0
        match_count = 0
        for line in f:
            parts = line.split('\t')
            if len(parts) <= max(name_idx, ko_idx):
                continue
            
            row_count += 1
            genome_name = parts[name_idx].strip()
            ko_val = parts[ko_idx].strip()

            # Keep genomes with zero target hits in the final matrix. Previously
            # they silently disappeared because a row was only created after a
            # matching KO had been found.
            genome_counts.setdefault(genome_name, {})
            
            # Extract clean KO
            ko = get_ko(ko_val)
            if not ko:
                continue
            
            # If this KO is mapped to one or more target genes
            if ko in ko_to_genes:
                match_count += 1
                for gene in ko_to_genes[ko]:
                    genome_counts[genome_name][gene] = genome_counts[genome_name].get(gene, 0) + 1
                    
    print(f"Processed {row_count} annotation rows. Found {match_count} hits matching the KO list.")
    return genome_counts

def write_output(output_path, genes_order, genome_counts):
    # Sort genome names alphabetically for clean output
    sorted_genomes = sorted(genome_counts.keys())
    
    output_dir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as out:
        # Write header
        header = ["Genome_Name"] + genes_order
        out.write('\t'.join(header) + '\n')
        
        for genome in sorted_genomes:
            row = [genome]
            gene_hits = genome_counts.get(genome, {})
            for gene in genes_order:
                count = gene_hits.get(gene, 0)
                row.append(str(count))
            out.write('\t'.join(row) + '\n')
            
    print(f"Success! Output written to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Map functional summary annotations into gene count matrix per genome.")
    parser.add_argument("-s", "--summary", required=True, help="Path to functional summary TSV file (e.g. bathyarchaeia_gtdb_all_functional_summary.tsv)")
    parser.add_argument("-g", "--gene", required=True, help="Path to gene TSV file (contains gene names and KO mappings)")
    parser.add_argument("-o", "--output", required=True, help="Path to write the output matrix TSV file")
    
    args = parser.parse_args()
    
    # 1. Load gene-to-KO mapping
    genes_order, ko_to_genes = load_gene_kos(args.gene)
    
    # 2. Process functional summary annotation file
    genome_counts = process_summary(args.summary, ko_to_genes)
    
    # 3. Write final matrix
    write_output(args.output, genes_order, genome_counts)

if __name__ == '__main__':
    main()
