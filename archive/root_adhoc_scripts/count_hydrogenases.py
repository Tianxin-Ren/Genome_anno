#!/usr/bin/env python3
import os
import sys
import re
from multiprocessing import Pool, cpu_count

def get_ko(val):
    if not isinstance(val, str):
        return None
    # Extract KO ID (e.g. K00436)
    m = re.search(r'(K\d{5})', val)
    return m.group(1) if m else None

def extract_fasta_records(faa_path, matched_proteins, mag_name):
    if not matched_proteins or not os.path.exists(faa_path):
        return []
    
    records = []
    targets = set(matched_proteins.keys())
    
    current_id = None
    current_lines = []
    
    with open(faa_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id in targets:
                    original_header = current_lines[0].strip()
                    group_annotation = matched_proteins[current_id]
                    modified_header = f"{original_header} [Genome:{mag_name}] [Group:{group_annotation}]\n"
                    current_lines[0] = modified_header
                    records.append("".join(current_lines))
                
                header = line.strip()
                current_id = header[1:].split()[0]
                if current_id in targets:
                    current_lines = [line]
                else:
                    current_lines = []
            else:
                if current_id in targets:
                    current_lines.append(line)
                    
        if current_id in targets and current_lines:
            original_header = current_lines[0].strip()
            group_annotation = matched_proteins[current_id]
            modified_header = f"{original_header} [Genome:{mag_name}] [Group:{group_annotation}]\n"
            current_lines[0] = modified_header
            records.append("".join(current_lines))
            
    return records

def process_single_file(args_tuple):
    file_path, mag_name, filter_set = args_tuple
    counts = {}
    matched_fasta_records = []
    if not os.path.exists(file_path):
        return mag_name, counts, matched_fasta_records

    try:
        matched_proteins = {}  # protein_id -> target_val
        with open(file_path, 'r') as f:
            header = f.readline()
            cols = [c.strip() for c in header.split('\t')]
            
            # Default column indices
            protein_idx = 0
            score_idx = 1
            evalue_idx = 2
            target_idx = 4
            
            if 'Protein_ID' in cols:
                protein_idx = cols.index('Protein_ID')
            if 'Score_Annotation' in cols:
                score_idx = cols.index('Score_Annotation')
            if 'Evalue_Annotation' in cols:
                evalue_idx = cols.index('Evalue_Annotation')
            if 'Target_Annotation' in cols:
                target_idx = cols.index('Target_Annotation')
                
            for line in f:
                parts = line.split('\t')
                if len(parts) <= max(protein_idx, score_idx, evalue_idx, target_idx):
                    continue
                
                # Extract KO from Score_Annotation or Evalue_Annotation
                ko_candidates = []
                ko_score = get_ko(parts[score_idx])
                if ko_score: 
                    ko_candidates.append(ko_score)
                ko_eval = get_ko(parts[evalue_idx])
                if ko_eval: 
                    ko_candidates.append(ko_eval)
                
                # Check if this row's KO is in the filter set
                is_match = False
                for ko in ko_candidates:
                    if ko in filter_set:
                        is_match = True
                        break
                        
                if is_match:
                    target_val = parts[target_idx].strip()
                    if target_val and target_val != '-':
                        counts[target_val] = counts.get(target_val, 0) + 1
                        protein_id = parts[protein_idx].strip()
                        matched_proteins[protein_id] = target_val
                        
        if matched_proteins:
            mag_dir = os.path.dirname(file_path)
            faa_path = os.path.join(mag_dir, f"{mag_name}.faa")
            if os.path.exists(faa_path):
                matched_fasta_records = extract_fasta_records(faa_path, matched_proteins, mag_name)
            else:
                print(f"WARNING: Protein file not found for {mag_name} at {faa_path}", file=sys.stderr)
    except Exception as e:
        print(f"Error processing {mag_name}: {e}", file=sys.stderr)
        
    return mag_name, counts, matched_fasta_records

def main():
    base_dir = "/dssg/home/acct-schhui/gses0/rentianxin/results/external_annotation/gses_all_rawMAGs_fa_V0"
    filter_path = "/dssg/home/acct-schhui/gses0/rentianxin/results/external_annotation/ICDLI_1953_F13_1_3/filter_Hydrogenase.tsv"
    output_path = "/dssg/home/acct-schhui/gses0/rentianxin/results/external_annotation/gses_all_rawMAGs_fa_V0_hydrogenase_counts.tsv"

    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    if len(sys.argv) > 2:
        filter_path = sys.argv[2]
    if len(sys.argv) > 3:
        output_path = sys.argv[3]

    output_faa_path = output_path.replace(".tsv", ".faa")
    if len(sys.argv) > 4:
        output_faa_path = sys.argv[4]

    print(f"Loading filter list from {filter_path}...", flush=True)
    filter_set = set()
    if os.path.exists(filter_path):
        with open(filter_path) as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    ko = get_ko(parts[0])
                    if ko:
                        filter_set.add(ko)
    else:
        print(f"ERROR: Filter list file not found at {filter_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(filter_set)} KO filters.", flush=True)

    print(f"Scanning directories in {base_dir}...", flush=True)
    tasks = []
    if os.path.exists(base_dir):
        with os.scandir(base_dir) as entries:
            for entry in entries:
                if entry.is_dir():
                    mag_name = entry.name
                    file_path = os.path.join(base_dir, mag_name, f"{mag_name}_advanced_hmm_resolved.tsv")
                    tasks.append((file_path, mag_name, filter_set))
    else:
        print(f"ERROR: Base directory not found at {base_dir}", file=sys.stderr)
        sys.exit(1)

    total_mags = len(tasks)
    print(f"Found {total_mags} MAG directories to process.", flush=True)
    if total_mags == 0:
        print("No subdirectories found. Exiting.", flush=True)
        sys.exit(0)

    # Use multiprocessing to maximize cores
    workers = min(cpu_count(), 64)
    print(f"Starting processing with {workers} parallel workers...", flush=True)

    results = []
    all_groups = set(['NiFe_G1', 'NiFe_G2', 'NiFe_G3', 'NiFe_G4', 'FeFe_GA', 'FeFe_GB', 'FeFe_GC', 'Fe_G1'])
    with Pool(workers) as pool:
        for idx, (mag_name, counts, fasta_records) in enumerate(pool.imap_unordered(process_single_file, tasks, chunksize=100)):
            results.append((mag_name, counts, fasta_records))
            for group in counts.keys():
                all_groups.add(group)
            if (idx + 1) % 5000 == 0 or (idx + 1) == total_mags:
                print(f"Processed {idx + 1}/{total_mags} genomes...", flush=True)

    print(f"Writing results to {output_path}...", flush=True)
    # Sort results by MAG name for clean output
    results.sort(key=lambda x: x[0])
    
    # Sort groups logically: NiFe first, then FeFe, then Fe, then others
    def group_sort_key(group_name):
        gl = group_name.lower()
        if 'nife' in gl:
            category_order = 0
        elif 'fefe' in gl:
            category_order = 1
        elif 'fe' in gl:
            category_order = 2
        else:
            category_order = 3
        return (category_order, gl)

    sorted_groups = sorted(all_groups, key=group_sort_key)
    
    with open(output_path, 'w') as out:
        header_cols = ["Genome_Name"] + sorted_groups
        out.write("\t".join(header_cols) + "\n")
        for mag_name, counts, _ in results:
            row_vals = [mag_name]
            for group in sorted_groups:
                row_vals.append(str(counts.get(group, 0)))
            out.write("\t".join(row_vals) + "\n")

    print(f"Writing annotated proteins to {output_faa_path}...", flush=True)
    with open(output_faa_path, 'w') as out_faa:
        for _, _, fasta_records in results:
            for record in fasta_records:
                out_faa.write(record)

    print(f"Success! TSV written to {output_path}, FAA written to {output_faa_path}", flush=True)

if __name__ == '__main__':
    main()
