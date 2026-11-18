#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
import argparse
import re

def check_executables(execs):
    """Check if required executables are on the PATH."""
    missing = []
    for name in execs:
        if not shutil.which(name):
            missing.append(name)
    if missing:
        print(f"Error: Missing required executables on PATH: {', '.join(missing)}")
        print("Please activate your Conda environment or install them first.")
        sys.exit(1)

def run_command(cmd, log_file=None):
    """Run a shell command, printing output on error or saving to log."""
    print(f"Running: {cmd}")
    if log_file:
        with open(log_file, "w") as log:
            res = subprocess.run(cmd, shell=True, stdout=log, stderr=subprocess.STDOUT)
    else:
        res = subprocess.run(cmd, shell=True)
    
    if res.returncode != 0:
        print(f"Error: Command failed with exit code {res.returncode}: {cmd}")
        sys.exit(1)

def sanitize_id(seq_id):
    """Sanitize sequence ID to make it safe for Newick tree and EPA-ng."""
    # Replace any character that is not alphanumeric, underscore, hyphen, or dot with underscore
    return re.sub(r'[^a-zA-Z0-9_\-\.]', '_', seq_id)

def sanitize_fasta(input_file, output_file):
    """Sanitize sequence IDs inside a FASTA file and replace non-standard U/u with X/x."""
    print(f"Sanitizing FASTA headers: {input_file} -> {output_file}")
    seen_ids = set()
    with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                header = line.strip()[1:]
                # Extract first word as sequence ID
                seq_id = header.split()[0]
                clean_id = sanitize_id(seq_id)
                
                # Prevent duplicate names by appending a suffix if needed
                final_id = clean_id
                counter = 1
                while final_id in seen_ids:
                    final_id = f"{clean_id}_{counter}"
                    counter += 1
                seen_ids.add(final_id)
                
                f_out.write(f">{final_id}\n")
            else:
                # Replace Selenocysteine 'U'/'u' with 'X'/'x' as it is not supported in phylogenetic engines
                clean_line = line.replace('U', 'X').replace('u', 'x')
                f_out.write(clean_line)

def parse_metadata(metadata_file, fasta_file):
    """
    Parse subtype information. 
    If metadata_file is provided, parse it (expected columns: SeqID \t Subtype).
    Otherwise, attempt to extract subtype from reference FASTA header.
    """
    tax_map = {}
    
    # Option A: Parse from metadata file
    if metadata_file and os.path.exists(metadata_file):
        print(f"Parsing metadata from: {metadata_file}")
        with open(metadata_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    # We sanitize the metadata ID as well to match the sanitized tree
                    seq_id = sanitize_id(parts[0])
                    subtype = parts[1]
                    tax_map[seq_id] = subtype
        return tax_map
    
    # Option B: Fallback to parsing from FASTA headers
    print(f"No metadata file provided. Attempting to parse subtypes from FASTA headers of: {fasta_file}")
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()[1:] # Remove '>'
                seq_id = header.split()[0]
                clean_id = sanitize_id(seq_id)
                
                # Pattern 1: Look for bracketed subtype e.g., >SeqID [Subtype]
                bracket_match = re.search(r'\[([^\]]+)\]', header)
                if bracket_match:
                    subtype = bracket_match.group(1).strip()
                    tax_map[clean_id] = subtype
                    continue
                
                # Pattern 2: Split by whitespace. If second field exists, treat it as subtype
                parts = header.split()
                if len(parts) >= 2:
                    subtype = parts[1]
                    tax_map[clean_id] = subtype
                    continue
                
                # Pattern 3: Look for last underscore e.g., >SeqID_Subtype
                if "_" in header:
                    parts = header.rsplit("_", 1)
                    subtype = parts[1]
                    tax_map[clean_id] = subtype
                    continue
                
                # Fallback: assign to a default group if nothing is found
                tax_map[clean_id] = "UnknownSubtype"
                
    return tax_map

def write_taxonomy_file(tax_map, output_path):
    """Write taxonomy map file in Gappa format (SeqID \t Lineage)."""
    with open(output_path, "w") as f:
        for seq_id, subtype in tax_map.items():
            # Clean subtype string for taxonomic lineage (semicolons are taxonomy level separators)
            lineage = f"Root;Subtypes;{subtype}"
            f.write(f"{seq_id}\t{lineage}\n")
    print(f"Taxonomy map file written to: {output_path}")

def get_fasta_ids(fasta_file):
    """Extract all sequence IDs from a FASTA file."""
    ids = set()
    if os.path.exists(fasta_file):
        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Get ID before first space
                    ids.add(line.strip()[1:].split()[0])
    return ids

def convert_stockholm_to_clean_fasta(sto_file, fasta_file, ref_ids=None):
    """
    Parse Stockholm file and output cleaned aligned FASTA:
    1. Parse the alignment columns using the #=GC RF line.
    2. Keep only columns that correspond to match states (GC RF is non-dot).
    3. Filter out reference sequences.
    """
    if ref_ids is None:
        ref_ids = set()
        
    seq_data = {} # seq_name -> list of sequence chunks
    rf_line = []  # list of RF chunks
    
    with open(sto_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("//"):
                continue
            if line.startswith("# STOCKHOLM"):
                continue
            
            if line.startswith("#=GC RF"):
                # RF annotation line: #=GC RF   xxx.x.x.x.x
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    rf_line.append(parts[2].strip().replace(" ", ""))
                continue
            
            if line.startswith("#"):
                # Other annotation lines, skip
                continue
                
            # Sequence line: name \t sequence
            parts = line.split(None, 1)
            if len(parts) == 2:
                name, seq = parts[0], parts[1].strip().replace(" ", "")
                if name not in seq_data:
                    seq_data[name] = []
                seq_data[name].append(seq)
                
    # Combine chunks
    full_rf = "".join(rf_line)
    
    # In HMMER Stockholm, match states in RF are marked with non-dot characters
    match_indices = []
    if full_rf:
        for i, char in enumerate(full_rf):
            if char != '.' and char != ' ':
                match_indices.append(i)
    
    print(f"Stockholm parsing: RF length = {len(full_rf)}, match states = {len(match_indices)}")
    
    with open(fasta_file, "w") as f_out:
        for name, chunks in seq_data.items():
            if name in ref_ids:
                continue
                
            full_seq = "".join(chunks)
            
            # Keep only match states
            if match_indices:
                if len(full_seq) == len(full_rf):
                    cleaned_chars = [full_seq[i] for i in match_indices]
                else:
                    # Fallback if length mismatch
                    cleaned_chars = [c for c in full_seq if not c.islower() and c != '.']
            else:
                # Fallback: remove lowercase letters (insertions) and dots
                cleaned_chars = [c for c in full_seq if not c.islower() and c != '.']
                
            # Normalize gaps to '-' and make uppercase, and replace 'U' with 'X'
            cleaned_seq = "".join(cleaned_chars).replace('.', '-').upper().replace('U', 'X')
            
            f_out.write(f">{name}\n{cleaned_seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Phylogenetic Placement Pipeline for Protein Subtyping")
    parser.add_argument("--ref-fasta", required=True, help="Path to reference protein FASTA file (approx 1,000 seqs)")
    parser.add_argument("--query-fasta", required=True, help="Path to query protein FASTA file (approx 80,000 seqs)")
    parser.add_argument("--out-dir", required=True, help="Output directory for results")
    parser.add_argument("--metadata", help="Optional TSV mapping reference SeqID to Subtype (format: SeqID \\t Subtype)")
    parser.add_argument("--tree-method", default="iqtree", choices=["iqtree", "fasttree"], help="Tree building tool")
    parser.add_argument("--threads", type=int, default=8, help="Number of CPU threads to use")
    parser.add_argument("--force", action="store_true", help="Force re-run of all steps, ignoring existing intermediate files")
    
    args = parser.parse_args()
    
    # 0. Check dependencies
    required_execs = ["mafft", "hmmbuild", "hmmalign", "epa-ng", "gappa"]
    if args.tree_method == "iqtree":
        required_execs.append("iqtree")
    else:
        required_execs.append("FastTree")
    check_executables(required_execs)
    
    # Create working directories
    os.makedirs(args.out_dir, exist_ok=True)
    tmp_dir = os.path.join(args.out_dir, "temp_run")
    os.makedirs(tmp_dir, exist_ok=True)
    
    # File paths for sanitized inputs
    ref_sanitized = os.path.join(tmp_dir, "ref.sanitized.faa")
    query_sanitized = os.path.join(tmp_dir, "query.sanitized.faa")
    
    ref_aln = os.path.join(tmp_dir, "ref.aln.faa")
    ref_tree = os.path.join(tmp_dir, "ref_tree")
    ref_tree_file = ref_tree + (".treefile" if args.tree_method == "iqtree" else ".tree")
    ref_model_file = ref_tree + ".log"
    ref_hmm = os.path.join(tmp_dir, "ref.hmm")
    query_sto = os.path.join(tmp_dir, "query.sto")
    query_aln = os.path.join(tmp_dir, "query.aln.faa")
    taxonomy_file = os.path.join(tmp_dir, "taxonomy_map.tsv")
    epa_out_dir = os.path.join(tmp_dir, "epa_out")
    jplace_file = os.path.join(epa_out_dir, "epa_result.jplace")
    final_output = os.path.join(args.out_dir, "subtype_classification_per_query.tsv")
    
    print("\n=== STEP 1: Aligning Reference Sequences (MAFFT) ===")
    if not args.force and os.path.exists(ref_aln) and os.path.getsize(ref_aln) > 0:
        print(f"Output file {ref_aln} already exists. Skipping MAFFT alignment.")
    else:
        # First, sanitize the raw reference FASTA to avoid Newick syntax errors
        sanitize_fasta(args.ref_fasta, ref_sanitized)
        run_command(f"mafft --auto --anysymbol --thread {args.threads} {ref_sanitized} > {ref_aln}", 
                    log_file=os.path.join(tmp_dir, "mafft_ref.log"))
    
    print("\n=== STEP 2: Building Reference Phylogenetic Tree ===")
    skip_step2 = not args.force and os.path.exists(ref_tree_file) and os.path.getsize(ref_tree_file) > 0
    if skip_step2:
        print(f"Output tree {ref_tree_file} already exists. Skipping tree building.")
        
    if args.tree_method == "iqtree":
        if not skip_step2:
            # IQ-TREE automatically performs ModelFinder to pick best substitution model
            run_command(f"iqtree -s {ref_aln} -m MFP -nt {args.threads} -pre {ref_tree} -quiet",
                        log_file=os.path.join(tmp_dir, "iqtree_ref.log"))
        model_string = "LG+G4" # Fallback model if parsing log fails
        # Extract best-fit model from iqtree log
        if os.path.exists(ref_model_file):
            with open(ref_model_file, "r") as f:
                for line in f:
                    if "Best-fit model:" in line:
                        model_string = line.split("Best-fit model:")[1].split()[0]
                        print(f"Detected best-fit model: {model_string}")
                        break
    else:
        if not skip_step2:
            # FastTree (Approximate Maximum Likelihood) - runs very fast
            run_command(f"FastTree -lg {ref_aln} > {ref_tree_file}",
                        log_file=os.path.join(tmp_dir, "fasttree_ref.log"))
        model_string = "LG"
        
    print("\n=== STEP 3: Preparing Subtype Taxonomy Information ===")
    if not args.force and os.path.exists(taxonomy_file) and os.path.getsize(taxonomy_file) > 0:
        print(f"Taxonomy mapping file {taxonomy_file} already exists. Skipping metadata preparation.")
    else:
        tax_map = parse_metadata(args.metadata, args.ref_fasta)
        write_taxonomy_file(tax_map, taxonomy_file)
    
    print("\n=== STEP 4: Building Reference Profile HMM ===")
    if not args.force and os.path.exists(ref_hmm) and os.path.getsize(ref_hmm) > 0:
        print(f"HMM profile {ref_hmm} already exists. Skipping hmmbuild.")
    else:
        # Use --symfrac 0.0 to ensure every column in the reference alignment is a match state,
        # making the HMM match states match the reference alignment length exactly.
        run_command(f"hmmbuild --symfrac 0.0 {ref_hmm} {ref_aln}",
                    log_file=os.path.join(tmp_dir, "hmmbuild.log"))
    
    print("\n=== STEP 5: Aligning Query Sequences to Reference Alignment ===")
    if not args.force and os.path.exists(query_sto) and os.path.getsize(query_sto) > 0:
        print(f"Stockholm alignment {query_sto} already exists. Skipping hmmalign.")
    else:
        # Sanitize query FASTA first
        sanitize_fasta(args.query_fasta, query_sanitized)
        # hmmalign keeps the reference alignment intact and maps queries onto the exact same columns.
        # We use --trim to trim off any unaligned prefix/suffix regions to ensure the queries align
        # exactly within the HMM boundaries.
        run_command(f"hmmalign --trim --mapali {ref_aln} -o {query_sto} {ref_hmm} {query_sanitized}",
                    log_file=os.path.join(tmp_dir, "hmmalign.log"))
    
    print("\n=== STEP 6: Converting Stockholm alignment to FASTA format ===")
    if not args.force and os.path.exists(query_aln) and os.path.getsize(query_aln) > 0:
        print(f"Cleaned FASTA alignment {query_aln} already exists. Skipping Stockholm conversion.")
    else:
        print("Parsing Stockholm and cleaning alignment: keeping match states and filtering out reference sequences...")
        # Parse reference sequence IDs from ref.aln.faa to filter them out of query.aln.faa
        ref_ids = get_fasta_ids(ref_aln)
        convert_stockholm_to_clean_fasta(query_sto, query_aln, ref_ids=ref_ids)
    
    print("\n=== STEP 7: Phylogenetic Placement of Queries on Reference Tree (EPA-ng) ===")
    if not args.force and os.path.exists(jplace_file) and os.path.getsize(jplace_file) > 0:
        print(f"Placement file {jplace_file} already exists. Skipping EPA-ng placement.")
    else:
        os.makedirs(epa_out_dir, exist_ok=True)
        # EPA-ng places aligned queries into the reference tree topology
        run_command(f"epa-ng --ref-msa {ref_aln} --tree {ref_tree_file} --model {model_string} --query {query_aln} --outdir {epa_out_dir} --threads {args.threads}",
                    log_file=os.path.join(tmp_dir, "epa_ng.log"))
    
    print("\n=== STEP 8: Taxonomy Assignment and Classification (Gappa) ===")
    if not args.force and os.path.exists(final_output) and os.path.getsize(final_output) > 0:
        print(f"Final output classification {final_output} already exists. Skipping Gappa classification.")
    else:
        # Gappa examine assign classifies query placement using taxon file mapping
        run_command(f"gappa examine assign --jplace-path {jplace_file} --taxon-file {taxonomy_file} --out-dir {args.out_dir} --file-prefix subtype_classification_ --per-query-results",
                    log_file=os.path.join(tmp_dir, "gappa.log"))
    
    if os.path.exists(final_output):
        print(f"\nSUCCESS! Subtype classification finished.")
        print(f"Final classifications saved to: {final_output}")
    else:
        print("\nError: Gappa run completed but output file not found. Check log files under temp_run/ directory.")

if __name__ == "__main__":
    main()
