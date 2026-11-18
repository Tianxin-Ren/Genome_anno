#!/usr/bin/env python3
import os
import sys
import argparse
import re

# Define the two operons and their properties
OPERONS = {
    "BRAO": {
        "genes": [f"BRAO_{i}" for i in range(1, 12)], # BRAO_1 to BRAO_11 (11 genes)
        "threshold": 6 # Ceiling(11/2) = 6
    },
    "FECOMJ": {
        "genes": [f"FECOMJ_{i:04d}" for i in range(94, 104)], # FECOMJ_0094 to FECOMJ_0103 (10 genes)
        "threshold": 5 # Ceiling(10/2) = 5
    }
}

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze double Operon structural integrity and physical synteny across genomes.")
    parser.add_argument("--reports", nargs="+", required=True, help="List of annotated target gene TSV reports from filter_blastp_results.py")
    parser.add_argument("--max_span", type=int, default=50000, help="Maximum physical distance span (bp) between matched genes to be considered closely grouped")
    parser.add_argument("--out_summary", required=True, help="Path to write the global operon prediction summary TSV")
    return parser.parse_args()

def parse_tsv_report(filepath):
    """Parses filter_blastp_results.py TSV output to extract hits with genomic positions."""
    hits = []
    if not os.path.exists(filepath):
        print(f"Warning: File {filepath} does not exist.", file=sys.stderr)
        return hits
        
    with open(filepath) as f:
        header = f.readline().strip().split("\t")
        
        # Ensure mandatory columns exist
        required = ["Query_Target_Gene", "Contig", "Start_Pos"]
        indices = {}
        for r in required:
            if r not in header:
                print(f"Warning: Missing required column {r} in {filepath}. Skipping.", file=sys.stderr)
                return hits
            indices[r] = header.index(r)
            
        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) <= max(indices.values()):
                continue
                
            q_gene = cols[indices["Query_Target_Gene"]].strip()
            contig = cols[indices["Contig"]].strip()
            
            try:
                start = int(cols[indices["Start_Pos"]])
            except ValueError:
                start = 0
                
            if q_gene and contig and start > 0:
                hits.append({
                    "query_gene": q_gene,
                    "contig": contig,
                    "start": start
                })
    return hits

def evaluate_operon_synteny(hits, operon_name, operon_cfg, max_span):
    """Evaluates whether the given operon is present in the genome based on coverage & physical span."""
    target_genes = set(operon_cfg["genes"])
    threshold = operon_cfg["threshold"]
    
    # 1. Filter hits that belong to this operon
    operon_hits = [h for h in hits if h["query_gene"] in target_genes]
    if not operon_hits:
        return {
            "is_present": "No",
            "matched_unique": 0,
            "contig": "",
            "physical_span": "",
            "matched_list": ""
        }
        
    # 2. Group by Contig
    contig_groups = {}
    for h in operon_hits:
        c = h["contig"]
        if c not in contig_groups:
            contig_groups[c] = []
        contig_groups[c].append(h)
        
    # 3. Analyze each contig
    best_eval = {
        "is_present": "No",
        "matched_unique": 0,
        "contig": "",
        "physical_span": "",
        "matched_list": ""
    }
    
    for contig, c_hits in contig_groups.items():
        # Get unique matched genes on this contig
        unique_genes = set(h["query_gene"] for h in c_hits)
        k = len(unique_genes)
        
        # Sort matched genes by start position
        sorted_hits = sorted(c_hits, key=lambda x: x["start"])
        
        # Check if the number of unique matched genes meets the threshold
        if k >= threshold:
            # We need to find if there is a subset of at least 'threshold' UNIQUE genes
            # that fit within 'max_span' (physical distance span).
            # To be highly robust and tolerant to duplicated paralogs:
            # We check sliding windows of size 'threshold' on the sorted coordinates.
            # However, a sliding window of coordinates might contain duplicates of the same query gene.
            # So, we check all contiguous blocks of hits that contain at least 'threshold' UNIQUE genes.
            
            n = len(sorted_hits)
            for i in range(n):
                for j in range(i + threshold - 1, n):
                    sub_hits = sorted_hits[i:j+1]
                    sub_unique = set(h["query_gene"] for h in sub_hits)
                    
                    if len(sub_unique) >= threshold:
                        span = sub_hits[-1]["start"] - sub_hits[0]["start"]
                        if span <= max_span:
                            # Found a valid physical cluster!
                            matched_list_sorted = sorted(list(sub_unique), key=lambda x: operon_cfg["genes"].index(x))
                            
                            # We keep the one with the maximum unique genes, or minimal span as the best hit
                            if k > best_eval["matched_unique"] or (k == best_eval["matched_unique"] and best_eval["is_present"] == "No"):
                                best_eval = {
                                    "is_present": "Yes",
                                    "matched_unique": k,
                                    "contig": contig,
                                    "physical_span": f"{sub_hits[0]['start']}-{sub_hits[-1]['start']}",
                                    "matched_list": ",".join(matched_list_sorted)
                                }
                            break # Found for this start position i, break inner loop to check next i
                            
        # If not present but we want to keep track of the highest partial hit
        if best_eval["is_present"] == "No" and k > best_eval["matched_unique"]:
            matched_list_sorted = sorted(list(unique_genes), key=lambda x: operon_cfg["genes"].index(x))
            best_eval = {
                "is_present": "No",
                "matched_unique": k,
                "contig": contig,
                "physical_span": f"{sorted_hits[0]['start']}-{sorted_hits[-1]['start']}",
                "matched_list": ",".join(matched_list_sorted)
            }
            
    return best_eval

def main():
    args = parse_args()
    
    # Structure to hold global records
    # records list of dict: [ {Genome, Operon, Total_Genes, Matched_Unique, Is_Present, Contig, Span, Matched_List} ]
    global_results = []
    
    print(f"Analyzing {len(args.reports)} annotated reports for operon synteny...", flush=True)
    
    for report_path in args.reports:
        # Extract genome name from directory structure
        # E.g. results/external_annotation/batch_name/genome_name/blastp/task_name/genome_target_genes_annotated.tsv
        # Or results/external_annotation/genome_name/blastp/task_name/genome_target_genes_annotated.tsv
        dir_parts = report_path.split("/")
        
        # Search backward to find genome name (the folder that has blastp as its immediate child)
        genome_name = "unknown"
        for idx in range(len(dir_parts) - 1):
            if dir_parts[idx] == "blastp":
                genome_name = dir_parts[idx - 1]
                break
        
        # Fallback if structure is different
        if genome_name == "unknown" or not genome_name:
            basename = os.path.basename(report_path)
            genome_name = basename.split("_target_genes_annotated")[0]
            
        # Parse TSV to get annotated query gene hits with coordinates
        hits = parse_tsv_report(report_path)
        
        # Evaluate each operon
        for op_name, op_cfg in OPERONS.items():
            eval_res = evaluate_operon_synteny(hits, op_name, op_cfg, args.max_span)
            
            global_results.append({
                "Genome": genome_name,
                "Operon": op_name,
                "Total_Genes": len(op_cfg["genes"]),
                "Matched_Unique_Genes": eval_res["matched_unique"],
                "Is_Present": eval_res["is_present"],
                "Contig": eval_res["contig"],
                "Physical_Span": eval_res["physical_span"],
                "Matched_Gene_List": eval_res["matched_list"]
            })
            
    # Sort results for readability
    global_results.sort(key=lambda x: (x["Genome"], x["Operon"]))
    
    # Write TSV Summary Report
    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)
    headers = [
        "Genome", "Operon", "Total_Genes", "Matched_Unique_Genes",
        "Is_Present", "Contig", "Physical_Span", "Matched_Gene_List"
    ]
    
    with open(args.out_summary, "w") as f_out:
        f_out.write("\t".join(headers) + "\n")
        for rec in global_results:
            f_out.write("\t".join(str(rec[h]) for h in headers) + "\n")
            
    print(f"Operon synteny analysis completed! Summary written to {args.out_summary}", flush=True)

if __name__ == "__main__":
    main()
