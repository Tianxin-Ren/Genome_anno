#!/usr/bin/env python3
import sys
import os
import argparse
import glob
import urllib.parse

def parse_args():
    parser = argparse.ArgumentParser(description="Filter raw BLASTp results and extract matched sequences")
    parser.add_argument("--raw_tsv", required=True, help="Raw BLAST tabular fmt 6 file")
    parser.add_argument("--gff3", required=True, help="Merged global GFF3 file")
    parser.add_argument("--faa", required=True, help="Global protein FAA file")
    parser.add_argument("--ffn", required=True, help="Global nucleotide FFN file")
    parser.add_argument("--binned", required=True, help="Path to all_bins.faa")
    parser.add_argument("--evalue", type=float, default=1e-5, help="Maximum E-value threshold")
    parser.add_argument("--identity", type=float, default=30.0, help="Minimum identity percentage")
    parser.add_argument("--coverage", type=float, default=50.0, help="Minimum query coverage percentage")
    parser.add_argument("--out_tsv", required=True, help="Output filtered TSV report")
    parser.add_argument("--out_faa", required=True, help="Output extracted matched FAA sequences")
    parser.add_argument("--out_ffn", required=True, help="Output extracted matched FFN sequences")
    return parser.parse_args()

def load_protein_to_bin_mapping(genome_dir):
    """Scan circular_refined directory to map protein IDs to MAG bins."""
    protein_to_bin = {}
    circular_refined_dir = os.path.join(genome_dir, "circular_refined")
    if os.path.exists(circular_refined_dir):
        print(f"Scanning bin FAA files in {circular_refined_dir} to map protein IDs...", flush=True)
        for bin_dir in glob.glob(os.path.join(circular_refined_dir, "*")):
            if os.path.isdir(bin_dir):
                bin_name = os.path.basename(bin_dir)
                for faa_file in glob.glob(os.path.join(bin_dir, "*.faa")):
                    with open(faa_file) as f:
                        for line in f:
                            if line.startswith(">"):
                                seq_id = line[1:].strip().split()[0]
                                protein_to_bin[seq_id] = bin_name
    return protein_to_bin

def load_gff3_annotations(gff3_path):
    """Parse GFF3 file to construct a dictionary mapping protein IDs to annotations."""
    protein_to_annotation = {}
    if not os.path.exists(gff3_path):
        print(f"Warning: GFF3 file {gff3_path} does not exist.", flush=True)
        return protein_to_annotation

    print(f"Parsing annotations from GFF3 {gff3_path}...", flush=True)
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            feature_type = cols[2]
            if feature_type == "CDS":
                attributes = cols[8]
                attr_dict = {}
                for item in attributes.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attr_dict[k] = v
                
                seq_id = attr_dict.get("ID")
                if seq_id:
                    gene = attr_dict.get("gene", "")
                    product = attr_dict.get("product", "")
                    protein_to_annotation[seq_id] = {
                        "gene": urllib.parse.unquote(gene),
                        "product": urllib.parse.unquote(product),
                        "contig": cols[0],
                        "start": int(cols[3])
                    }
    return protein_to_annotation


def stream_and_enrich_fasta(in_fasta, out_fasta, matched_subjects, genome_name):
    """Streams a FASTA file and writes matched records with enriched headers."""
    if not os.path.exists(in_fasta):
        print(f"Warning: FASTA file {in_fasta} does not exist. Skipping extraction.", flush=True)
        return
        
    print(f"Extracting matched sequences from {in_fasta} -> {out_fasta}...", flush=True)
    os.makedirs(os.path.dirname(out_fasta), exist_ok=True)
    
    with open(in_fasta) as f_in, open(out_fasta, "w") as f_out:
        keep = False
        seq_id = None
        for line in f_in:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id in matched_subjects:
                    keep = True
                    # Construct enriched header
                    metadata = matched_subjects[seq_id]
                    header = f">{seq_id} [Target={metadata['target']}, Genome={genome_name}, Identity={metadata['identity']}%, Coverage={metadata['coverage']}%, Product={metadata['product']}]\n"
                    f_out.write(header)
                else:
                    keep = False
            else:
                if keep:
                    f_out.write(line)

def main():
    args = parse_args()
    
    # Extract genome root dir
    genome_dir = os.path.dirname(os.path.dirname(args.raw_tsv))
    genome_name = os.path.basename(genome_dir)
    
    # 1. Load mappings
    protein_to_bin = load_protein_to_bin_mapping(genome_dir)
    protein_to_annotation = load_gff3_annotations(args.gff3)
    
    # 2. Filter raw BLAST alignment table
    print(f"Filtering raw alignments {args.raw_tsv}...", flush=True)
    matched_subjects = {}  # subject_id -> metadata_dict
    
    # Format of blastfmt 6:
    # 0: qseqid, 1: sseqid, 2: pident, 3: length, 4: mismatch, 5: gapopen,
    # 6: qstart, 7: qend, 8: sstart, 9: send, 10: evalue, 11: bitscore, 12: qlen, 13: slen
    
    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    
    filtered_records = []
    
    if os.path.exists(args.raw_tsv) and os.path.getsize(args.raw_tsv) > 0:
        with open(args.raw_tsv) as f_in:
            for line in f_in:
                if not line.strip():
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 14:
                    continue
                
                qseqid = cols[0]
                sseqid = cols[1]
                pident = float(cols[2])
                alignment_length = int(cols[3])
                evalue = float(cols[10])
                bitscore = float(cols[11])
                qlen = int(cols[12])
                
                # Calculate dynamic query coverage percentage
                coverage = (alignment_length / qlen) * 100.0
                
                # Check user thresholds
                if evalue <= args.evalue and pident >= args.identity and coverage >= args.coverage:
                    # Retrieve metadata
                    bin_name = protein_to_bin.get(sseqid, "unbinned")
                    anno = protein_to_annotation.get(sseqid, {"gene": "", "product": "", "contig": "", "start": 0})
                    
                    record = {
                        "Query_Target_Gene": qseqid,
                        "Subject_Homolog_ID": sseqid,
                        "Identity_Pct": pident,
                        "Alignment_Length": alignment_length,
                        "Query_Coverage_Pct": round(coverage, 2),
                        "E_value": evalue,
                        "Bit_Score": bitscore,
                        "Source_Bin_Or_Unbinned": bin_name,
                        "Annotated_Gene": anno.get("gene", ""),
                        "Annotated_Product": anno.get("product", ""),
                        "Contig": anno.get("contig", ""),
                        "Start_Pos": anno.get("start", 0)
                    }
                    filtered_records.append(record)
                    
                    # Store for FASTA sequence extraction (keep best hit if duplicate hits exist)
                    if sseqid not in matched_subjects or pident > matched_subjects[sseqid]["identity"]:
                        matched_subjects[sseqid] = {
                            "target": qseqid,
                            "identity": pident,
                            "coverage": round(coverage, 2),
                            "product": anno.get("product", "")
                        }
                        
    # Write TSV report
    headers = [
        "Query_Target_Gene", "Subject_Homolog_ID", "Identity_Pct", 
        "Alignment_Length", "Query_Coverage_Pct", "E_value", 
        "Bit_Score", "Source_Bin_Or_Unbinned", "Annotated_Gene", "Annotated_Product",
        "Contig", "Start_Pos"
    ]

    
    with open(args.out_tsv, "w") as f_out:
        f_out.write("\t".join(headers) + "\n")
        for rec in filtered_records:
            f_out.write("\t".join(str(rec[h]) for h in headers) + "\n")
            
    print(f"Found and annotated {len(filtered_records)} high-confidence homologous genes.", flush=True)
    
    # 3. Stream-extract sequences
    stream_and_enrich_fasta(args.faa, args.out_faa, matched_subjects, genome_name)
    stream_and_enrich_fasta(args.ffn, args.out_ffn, matched_subjects, genome_name)
    
    # Create empty files if no matched subjects to avoid Snakemake errors
    if not matched_subjects:
        print("No matches found. Creating empty FAA and FFN files.", flush=True)
        open(args.out_faa, "w").close()
        open(args.out_ffn, "w").close()
        
    print("Alignment filtering and extraction completed successfully!", flush=True)

if __name__ == "__main__":
    main()
