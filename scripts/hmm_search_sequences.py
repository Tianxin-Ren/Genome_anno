#!/usr/bin/env python3
import sys
import os
import argparse
import time
import requests

# Global headers to prevent requests from being blocked/dropped by remote servers
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36"
}

def write_combined_fasta(query_sequences, fetched_fasta_str, output_path):
    written_ids = set()
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
        
    with open(output_path, "w") as out_f:
        # Write query sequences first
        for q_id, q_seq in query_sequences:
            clean_id = q_id.split()[0]
            out_f.write(f">{clean_id}\n{q_seq}\n")
            written_ids.add(clean_id)
            
        # Parse and write fetched sequences
        current_id = None
        current_seq = []
        
        lines = fetched_fasta_str.split("\n")
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id and current_id not in written_ids:
                    out_f.write(f">{current_id}\n" + "".join(current_seq) + "\n")
                    written_ids.add(current_id)
                header = line[1:]
                current_id = header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id and current_id not in written_ids:
            out_f.write(f">{current_id}\n" + "".join(current_seq) + "\n")
            written_ids.add(current_id)
            
    print(f"Successfully wrote {len(written_ids)} unique sequences (including query) to '{output_path}'.")

def parse_fasta(filepath):
    sequences = []
    current_id = None
    current_seq = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences.append((current_id, "".join(current_seq)))
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences.append((current_id, "".join(current_seq)))
    return sequences

def extract_accession(sseqid):
    if "|" in sseqid:
        parts = sseqid.split("|")
        if len(parts) >= 2:
            return parts[1]
    if ":" in sseqid:
        parts = sseqid.split(":")
        return parts[-1]
    return sseqid

def search_uniprot(query_id, query_seq, database, min_identity, min_coverage, max_seqs):
    print(f"Submitting query '{query_id}' to EBI NCBI BLAST using database '{database}'...")
    url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
    payload = {
        "email": "blast_pipeline@example.com",
        "program": "blastp",
        "stype": "protein",
        "database": database,
        "sequence": query_seq
    }
    
    try:
        r = requests.post(url, data=payload, headers=HEADERS)
        if r.status_code != 200:
            print(f"Error submitting BLAST for {query_id}: {r.status_code} {r.text}", file=sys.stderr)
            return []
        job_id = r.text.strip()
    except Exception as e:
        print(f"Exception submitting BLAST: {e}", file=sys.stderr)
        return []
        
    if not job_id:
        print(f"Error: No jobId returned for {query_id}", file=sys.stderr)
        return []
    
    print(f"BLAST job submitted successfully. Job ID: {job_id}. Polling status...")
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
    
    attempts = 0
    max_attempts = 90  # 90 * 10s = 900s (15 minutes)
    finished = False
    while attempts < max_attempts:
        time.sleep(10)
        try:
            status_r = requests.get(status_url, headers=HEADERS)
            if status_r.status_code != 200:
                print(f"Warning: Failed to get job status (attempt {attempts+1}): {status_r.status_code}", file=sys.stderr)
                attempts += 1
                continue
            status = status_r.text.strip()
            print(f"Job status: {status}")
            if status == "FINISHED":
                finished = True
                break
            elif status in ["FAILED", "ERROR", "NOT_FOUND"]:
                print(f"Error: Job {job_id} failed with status: {status}", file=sys.stderr)
                return []
        except Exception as e:
            print(f"Warning: Exception while polling status: {e}", file=sys.stderr)
        attempts += 1
        
    if not finished:
        print(f"Error: Job {job_id} timed out.", file=sys.stderr)
        return []
        
    print("Retrieving BLAST results (JSON format)...")
    results_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/json"
    
    try:
        results_r = requests.get(results_url, headers=HEADERS)
        if results_r.status_code != 200:
            print(f"Error retrieving results: {results_r.status_code} {results_r.text}", file=sys.stderr)
            return []
        data = results_r.json()
    except Exception as e:
        print(f"Error parsing JSON results: {e}", file=sys.stderr)
        return []
        
    hits = data.get("hits", [])
    print(f"Found {len(hits)} raw hits in BLAST result.")
    
    filtered_accessions = []
    query_len = len(query_seq)
    
    for hit in hits:
        accession = hit.get("hit_acc")
        if not accession:
            continue
        hsps = hit.get("hit_hsps", [])
        if not hsps:
            continue
        best_hsp = hsps[0]
        
        identity = best_hsp.get("hsp_identity", 0.0)
        # Normalize identity to 0.0 - 1.0 range
        identity_frac = identity / 100.0 if identity > 1.0 else identity
        
        q_from = best_hsp.get("hsp_query_from", 1)
        q_to = best_hsp.get("hsp_query_to", 1)
        coverage = (q_to - q_from + 1) / query_len
        
        if identity_frac >= min_identity and coverage >= min_coverage:
            filtered_accessions.append(accession)
            if len(filtered_accessions) >= max_seqs:
                break
                
    print(f"Filtered {len(filtered_accessions)} accessions matching identity >= {min_identity} and coverage >= {min_coverage}.")
    return filtered_accessions

def fetch_uniprot_sequences(accessions):
    if not accessions:
        return ""
    print(f"Fetching sequences for {len(accessions)} accessions from UniProt...")
    batch_size = 50
    all_fasta = []
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        query_str = " OR ".join([f"accession:{acc}" for acc in batch])
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": query_str,
            "format": "fasta",
            "size": batch_size
        }
        try:
            r = requests.get(url, params=params, headers=HEADERS)
            if r.status_code == 200:
                all_fasta.append(r.text)
            else:
                print(f"Warning: Failed to fetch batch {i//batch_size + 1}: {r.status_code} {r.text}", file=sys.stderr)
            time.sleep(1)  # Be polite to API
        except Exception as e:
            print(f"Error fetching batch {i//batch_size + 1}: {e}", file=sys.stderr)
            
    return "\n".join(all_fasta)

def search_ncbi(query_id, query_seq, min_identity, min_coverage, max_seqs):
    try:
        from Bio.Blast import NCBIWWW
        from Bio.Blast import NCBIXML
    except ImportError:
        print("Error: Biopython is required for NCBI BLAST search. Please install it or use UniProt databases.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Submitting query '{query_id}' to NCBI BLAST (blastp against 'nr')... This may take several minutes.")
    try:
        result_handle = NCBIWWW.qblast("blastp", "nr", query_seq, hitlist_size=max_seqs)
        blast_record = NCBIXML.read(result_handle)
    except Exception as e:
        print(f"Error running NCBI BLAST: {e}", file=sys.stderr)
        return []
        
    filtered_accessions = []
    query_len = len(query_seq)
    
    for alignment in blast_record.alignments:
        if not alignment.hsps:
            continue
        best_hsp = alignment.hsps[0]
        
        identity = best_hsp.identities / best_hsp.align_length if best_hsp.align_length > 0 else 0
        coverage = (best_hsp.query_end - best_hsp.query_start + 1) / query_len
        
        if identity >= min_identity and coverage >= min_coverage:
            filtered_accessions.append(alignment.accession)
            if len(filtered_accessions) >= max_seqs:
                break
                
    print(f"Filtered {len(filtered_accessions)} accessions matching identity >= {min_identity} and coverage >= {min_coverage} from NCBI.")
    return filtered_accessions

def fetch_ncbi_sequences(accessions):
    if not accessions:
        return ""
    try:
        from Bio import Entrez
    except ImportError:
        print("Error: Biopython is required for fetching NCBI sequences.", file=sys.stderr)
        sys.exit(1)
        
    Entrez.email = "blast_pipeline@example.com"
    print(f"Fetching sequences for {len(accessions)} accessions from NCBI Entrez...")
    batch_size = 50
    all_fasta = []
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        id_list = ",".join(batch)
        try:
            handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
            all_fasta.append(handle.read())
            handle.close()
            time.sleep(1)  # Be polite to Entrez API
        except Exception as e:
            print(f"Error fetching NCBI batch {i//batch_size + 1}: {e}", file=sys.stderr)
            
    return "\n".join(all_fasta)

def main():
    parser = argparse.ArgumentParser(description="Search homologs for HMM building using NCBI BLAST or UniProt.")
    parser.add_argument("--query", required=True, help="Input query FASTA file (one or more reference sequences)")
    parser.add_argument("--database", default="uniprotkb_swissprot", 
                        choices=["uniprotkb", "uniprotkb_swissprot", "uniprotkb_reference_proteomes", "ncbi"], 
                        help="Database to search against")
    parser.add_argument("--identity", type=float, default=0.3, help="Minimum sequence identity threshold (0.0 to 1.0)")
    parser.add_argument("--coverage", type=float, default=0.5, help="Minimum query coverage threshold (0.0 to 1.0)")
    parser.add_argument("--max_seqs", type=int, default=200, help="Maximum number of sequences to collect")
    parser.add_argument("--output", required=True, help="Output combined FASTA file path")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.query):
        print(f"Error: Query file '{args.query}' does not exist.", file=sys.stderr)
        sys.exit(1)
        
    query_sequences = parse_fasta(args.query)
    if not query_sequences:
        print(f"Error: No sequences found in query file '{args.query}'.", file=sys.stderr)
        sys.exit(1)
        
    all_hit_accessions = []
    
    # Search for each query sequence
    for q_id, q_seq in query_sequences:
        if args.database == "ncbi":
            hits = search_ncbi(q_id, q_seq, args.identity, args.coverage, args.max_seqs)
        else:
            hits = search_uniprot(q_id, q_seq, args.database, args.identity, args.coverage, args.max_seqs)
        all_hit_accessions.extend(hits)
        
    # De-duplicate accessions
    unique_accessions = list(dict.fromkeys(all_hit_accessions))[:args.max_seqs]
    
    # Fetch sequences
    if args.database == "ncbi":
        fetched_fasta = fetch_ncbi_sequences(unique_accessions)
    else:
        fetched_fasta = fetch_uniprot_sequences(unique_accessions)
        
    # Write combined output
    write_combined_fasta(query_sequences, fetched_fasta, args.output)

if __name__ == "__main__":
    main()
