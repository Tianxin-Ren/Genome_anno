#!/usr/bin/env python3
import os
import sys
import json
import argparse
import urllib.request
import urllib.parse
import time

def parse_args():
    parser = argparse.ArgumentParser(description="Collect high-confidence gene sequences from UniProt")
    
    # Input sources
    parser.add_argument("--genes-json", help="JSON string representing the list of genes to collect")
    parser.add_argument("--kegg-tsv", help="Path to a KEGG annotated TSV file containing KO and gene_name columns")
    
    parser.add_argument("--min-sequences", type=int, default=5, help="Minimum number of sequences to collect per gene")
    parser.add_argument("--max-sequences", type=int, default=7, help="Maximum number of sequences to collect per gene")
    parser.add_argument("--output-dir", required=True, help="Output directory to save FASTA files")
    parser.add_argument("--done-flag", required=True, help="Path to sentinel file to create upon completion")
    return parser.parse_args()

def query_uniprot_api(query, size=50):
    """Query the UniProt Search API and return the list of matching entries."""
    params = {
        "query": query,
        "format": "json",
        "size": size
    }
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
    
    max_retries = 3
    retry_delay = 2
    
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "GenomeAnno-DataPrepAgent/1.0"})
            with urllib.request.urlopen(req, timeout=30) as response:
                if response.status == 200:
                    data = json.loads(response.read().decode("utf-8"))
                    return data.get("results", [])
                else:
                    print(f"Warning: UniProt API returned status {response.status} (attempt {attempt+1}/{max_retries})", file=sys.stderr)
        except Exception as e:
            print(f"Warning: UniProt query error: {e} (attempt {attempt+1}/{max_retries})", file=sys.stderr)
        
        if attempt < max_retries - 1:
            time.sleep(retry_delay)
            retry_delay *= 2
            
    print(f"Error: Failed to query UniProt after {max_retries} attempts.", file=sys.stderr)
    return []

def extract_organism_name(entry):
    org = entry.get("organism", {})
    return org.get("scientificName", "")

def extract_annotation_score(entry):
    return entry.get("annotationScore", 0)

def extract_sequence(entry):
    return entry.get("sequence", {}).get("value", "")

def is_reviewed(entry):
    entry_type = entry.get("entryType", "")
    return "reviewed" in entry_type.lower() and "unreviewed" not in entry_type.lower()

# =====================================================================
# Validation Helpers
# =====================================================================

def verify_gene_name(entry, alternatives):
    """Verify that the entry's actual gene name list contains at least one of our alternative query names (case-insensitive exact match)."""
    entry_genes = entry.get("genes", [])
    entry_names = []
    for g in entry_genes:
        if "geneName" in g:
            entry_names.append(g["geneName"].get("value", "").lower())
        for syn in g.get("synonyms", []):
            entry_names.append(syn.get("value", "").lower())
        for orf in g.get("orfNames", []):
            entry_names.append(orf.get("value", "").lower())
            
    for alt in alternatives:
        if alt.lower() in entry_names:
            return True
            
    return False

def extract_ec_numbers(entry):
    """Extract all EC numbers from entry proteinDescription."""
    ecs = []
    desc = entry.get("proteinDescription", {})
    
    def walk_name(name_dict):
        if not name_dict:
            return
        for ec in name_dict.get("ecNumbers", []):
            val = ec.get("value")
            if val:
                ecs.append(val.strip())
                
    walk_name(desc.get("recommendedName"))
    for alt in desc.get("alternativeNames", []):
        walk_name(alt)
    for sub in desc.get("submissionNames", []):
        walk_name(sub)
    return ecs

def verify_ec_number(entry, expected_ecs):
    """Verify that the entry's EC number list matches one of our expected EC numbers."""
    if not expected_ecs:
        return True # Skip EC validation if not provided
    entry_ecs = extract_ec_numbers(entry)
    for exp in expected_ecs:
        if exp in entry_ecs:
            return True
    return False

def collect_sequences_for_gene(gene_info, min_seqs, max_seqs):
    gene_name = gene_info.get("name")
    alternatives = gene_info.get("alternatives", [gene_name])
    organism = gene_info.get("organism")
    custom_query = gene_info.get("query")
    ko_id = gene_info.get("ko", "")
    expected_ecs = gene_info.get("ec_numbers", [])
    
    # Validation strategy
    # Step 1: Strict Validation (Exact Gene Name Match & EC match)
    # Step 2: Medium Validation (Exact Gene Name Match only)
    # Step 3: Loose Validation (No verification - fallback)
    
    def collect_with_validation_level(strict_gene, strict_ec):
        reviewed_entries = []
        unreviewed_entries = []
        seen_accessions = set()
        
        # A. Query Reviewed (Swiss-Prot)
        if custom_query:
            reviewed_query = f"({custom_query}) AND reviewed:true"
            results = query_uniprot_api(reviewed_query, size=50)
            for entry in results:
                acc = entry.get("primaryAccession")
                if acc and acc not in seen_accessions:
                    # Validate
                    if strict_gene and not verify_gene_name(entry, alternatives):
                        continue
                    if strict_ec and not verify_ec_number(entry, expected_ecs):
                        continue
                    reviewed_entries.append(entry)
                    seen_accessions.add(acc)
        else:
            for name in alternatives:
                if name.isalnum():
                    base_query = f"gene:{name}"
                else:
                    base_query = f'gene:"{name}"'
                    
                if organism:
                    if organism.isalnum():
                        base_query += f" AND taxonomy_name:{organism}"
                    else:
                        base_query += f' AND taxonomy_name:"{organism}"'
                        
                reviewed_query = f"({base_query}) AND reviewed:true"
                results = query_uniprot_api(reviewed_query, size=50)
                
                for entry in results:
                    acc = entry.get("primaryAccession")
                    if acc and acc not in seen_accessions:
                        # Validate
                        if strict_gene and not verify_gene_name(entry, alternatives):
                            continue
                        if strict_ec and not verify_ec_number(entry, expected_ecs):
                            continue
                        reviewed_entries.append(entry)
                        seen_accessions.add(acc)
                
                if len(reviewed_entries) >= max_seqs:
                    break
                    
        # If we have enough reviewed, return them capped
        if len(reviewed_entries) >= max_seqs:
            return reviewed_entries[:max_seqs]
            
        # B. Query Unreviewed (TrEMBL) to reach min_seqs
        if len(reviewed_entries) < min_seqs:
            if custom_query:
                results = query_uniprot_api(custom_query, size=100)
                for entry in results:
                    acc = entry.get("primaryAccession")
                    if acc and acc not in seen_accessions:
                        if strict_gene and not verify_gene_name(entry, alternatives):
                            continue
                        if strict_ec and not verify_ec_number(entry, expected_ecs):
                            continue
                        unreviewed_entries.append(entry)
                        seen_accessions.add(acc)
            else:
                for name in alternatives:
                    if name.isalnum():
                        base_query = f"gene:{name}"
                    else:
                        base_query = f'gene:"{name}"'
                        
                    if organism:
                        if organism.isalnum():
                            base_query += f" AND taxonomy_name:{organism}"
                        else:
                            base_query += f' AND taxonomy_name:"{organism}"'
                            
                    results = query_uniprot_api(base_query, size=100)
                    for entry in results:
                        acc = entry.get("primaryAccession")
                        if acc and acc not in seen_accessions:
                            if strict_gene and not verify_gene_name(entry, alternatives):
                                continue
                            if strict_ec and not verify_ec_number(entry, expected_ecs):
                                continue
                            unreviewed_entries.append(entry)
                            seen_accessions.add(acc)
                            
                    if len(reviewed_entries) + len(unreviewed_entries) >= min_seqs:
                        break
                        
        # Sort unreviewed results by annotation score (descending), then sequence length (descending)
        unreviewed_entries.sort(key=lambda x: (
            -float(extract_annotation_score(x)),
            -len(extract_sequence(x))
        ))
        
        selected_entries = reviewed_entries + unreviewed_entries
        
        # If we have between min_seqs and max_seqs, keep them all. Otherwise cap at max_seqs or min_seqs
        if len(reviewed_entries) >= min_seqs:
            return reviewed_entries[:max_seqs]
        else:
            return selected_entries[:min_seqs]

    # Try Step 1: Strict Validation (Gene name + EC number)
    print(f"Gene '{gene_name}' ({ko_id}): running STRICT validation (Gene name & EC: {expected_ecs})...", file=sys.stderr)
    selected = collect_with_validation_level(strict_gene=True, strict_ec=True)
    if len(selected) >= min_seqs:
        return selected
        
    # Try Step 2: Medium Validation (Gene name only, if EC was specified but made us fail)
    if expected_ecs:
        print(f"Warning: Strict validation returned only {len(selected)} sequences for gene '{gene_name}'. Falling back to MEDIUM validation (Gene name matching only).", file=sys.stderr)
        selected = collect_with_validation_level(strict_gene=True, strict_ec=False)
        if len(selected) >= min_seqs:
            return selected
            
    # Try Step 3: Loose Validation (No matching checks - fallback query results)
    print(f"Warning: Validation returned only {len(selected)} sequences for gene '{gene_name}'. Falling back to LOOSE validation (no matching checks).", file=sys.stderr)
    return collect_with_validation_level(strict_gene=False, strict_ec=False)

def parse_kegg_tsv(filepath):
    """Parse TSV file containing KEGG annotations and extract genes list."""
    genes = []
    if not os.path.exists(filepath):
        print(f"Error: KEGG TSV file not found: {filepath}", file=sys.stderr)
        sys.exit(1)
        
    with open(filepath, "r", encoding="utf-8") as f:
        # Read and parse header
        header_line = f.readline()
        if not header_line:
            print("Error: Empty TSV file.", file=sys.stderr)
            sys.exit(1)
            
        header = [col.strip() for col in header_line.strip().split('\t')]
        try:
            ko_idx = header.index("KO")
            gene_name_idx = header.index("gene_name")
            # ec_number is optional
            ec_idx = header.index("ec_number") if "ec_number" in header else -1
        except ValueError:
            print("Error: TSV file must contain 'KO' and 'gene_name' columns in the header.", file=sys.stderr)
            sys.exit(1)
            
        for line_num, line in enumerate(f, start=2):
            parts = [p.strip() for p in line.strip().split('\t')]
            if len(parts) <= max(ko_idx, gene_name_idx):
                continue
            
            ko = parts[ko_idx]
            raw_gene_name = parts[gene_name_idx]
            
            if not raw_gene_name or raw_gene_name == "-" or raw_gene_name.lower() == "nan":
                continue
                
            # frdA, fdrA; succinate dehydrogenase flavoprotein subunit -> names_part = "frdA, fdrA"
            names_part = raw_gene_name.split(';')[0]
            names = [n.strip() for n in names_part.split(',') if n.strip()]
            
            # Extract expected EC numbers if present
            expected_ecs = []
            if ec_idx != -1 and len(parts) > ec_idx:
                raw_ec = parts[ec_idx]
                expected_ecs = [e.strip() for e in raw_ec.split() if e.strip() and e.strip() != "-"]
            
            if names:
                genes.append({
                    "name": names[0],
                    "alternatives": names,
                    "ko": ko,
                    "ec_numbers": expected_ecs
                })
                
    return genes

def main():
    args = parse_args()
    
    # Validate sequence range
    if args.min_sequences > args.max_sequences:
        print(f"Error: --min-sequences ({args.min_sequences}) cannot be greater than --max-sequences ({args.max_sequences})", file=sys.stderr)
        sys.exit(1)
        
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Resolve genes list from input source
    if not args.genes_json and not args.kegg_tsv:
        print("Error: Either --genes-json or --kegg-tsv must be provided.", file=sys.stderr)
        sys.exit(1)
        
    if args.genes_json:
        try:
            genes = json.loads(args.genes_json)
        except Exception as e:
            print(f"Error parsing --genes-json: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Parsing KEGG TSV file from: {args.kegg_tsv}", file=sys.stderr)
        genes = parse_kegg_tsv(args.kegg_tsv)
        
    print(f"Starting collection of sequences for {len(genes)} genes in output directory: {args.output_dir}", file=sys.stderr)
    print(f"Sequence count constraint per gene: [{args.min_sequences}, {args.max_sequences}]", file=sys.stderr)
    
    for gene_info in genes:
        gene_name = gene_info.get("name")
        if not gene_name:
            print("Error: Gene name not specified in entry:", gene_info, file=sys.stderr)
            continue
            
        entries = collect_sequences_for_gene(gene_info, args.min_sequences, args.max_sequences)
        
        # Save to FASTA (name file by KO if available, otherwise by gene name)
        file_id = gene_info.get("ko") if gene_info.get("ko") else gene_name
        fasta_path = os.path.join(args.output_dir, f"{file_id}.fasta")
        with open(fasta_path, "w") as f:
            for entry in entries:
                acc = entry.get("primaryAccession")
                seq = extract_sequence(entry)
                org_name = extract_organism_name(entry)
                ann_score = extract_annotation_score(entry)
                status = "reviewed" if is_reviewed(entry) else "unreviewed"
                
                if not acc or not seq:
                    continue
                    
                # Format header (prefix sequence ID with KO if available)
                if gene_info.get("ko"):
                    header = f">{gene_info['ko']}_{acc} {gene_name} {org_name} [annotation_score={ann_score}] [status={status}]"
                else:
                    header = f">{acc} {gene_name} {org_name} [annotation_score={ann_score}] [status={status}]"
                    
                entry_ecs = extract_ec_numbers(entry)
                if entry_ecs:
                    header += f" [ec={','.join(entry_ecs)}]"
                    
                # Write to FASTA (split sequence lines at 60 chars)
                f.write(f"{header}\n")
                for i in range(0, len(seq), 60):
                    f.write(f"{seq[i:i+60]}\n")
                    
        print(f"Saved {len(entries)} sequences for gene '{gene_name}' ({file_id}) to {fasta_path}", file=sys.stderr)
        # Sleep briefly between genes to avoid API rate limits
        time.sleep(0.5)
        
    # Write done flag
    with open(args.done_flag, "w") as f:
        f.write("done\n")
        
    print(f"Finished collecting genes. Sentinel file written to {args.done_flag}", file=sys.stderr)

if __name__ == "__main__":
    main()
