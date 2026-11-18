import os
import argparse
import pandas as pd
import shutil
import re
from Bio import SeqIO

def get_attributes(attr_str):
    """Parse GFF3 attributes string into a dictionary."""
    attrs = {}
    for item in attr_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs

def subset_results(bin_fasta, assembly_gff, assembly_faa, assembly_tsv, assembly_eggnog,
                   assembly_kofam, assembly_dbcan,
                   out_gff, out_faa, out_tsv, out_eggnog, out_kofam, out_dbcan):
    
    # 1. Get simplified contig IDs from the bin
    bin_contig_ids = set()
    for record in SeqIO.parse(bin_fasta, "fasta"):
        # Unified simplified ID: u124927ctg... -> u124927
        simple_id = record.id.split('ctg')[0]
        bin_contig_ids.add(simple_id)
    
    print(f"Found {len(bin_contig_ids)} simplified contigs in bin.")

    # 2. Parse GFF3 to get the mapping of Protein IDs to our target Contigs
    target_protein_ids = set()
    if os.path.exists(assembly_gff):
        with open(assembly_gff, 'r') as f_in, open(out_gff, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    # Also simplify ##sequence-region in output
                    if line.startswith('##sequence-region'):
                        parts = line.split()
                        if len(parts) > 1:
                            parts[1] = parts[1].split('ctg')[0]
                            line = " ".join(parts) + "\n"
                    f_out.write(line)
                    continue
                parts = line.split('\t')
                if len(parts) >= 9:
                    contig_id_full = parts[0]
                    contig_id_simple = contig_id_full.split('ctg')[0]
                    
                    if contig_id_simple in bin_contig_ids:
                        # Use simplified ID in the subsetted GFF
                        parts[0] = contig_id_simple
                        f_out.write("\t".join(parts))
                        # Extract Protein ID (usually ID=...)
                        attrs = get_attributes(parts[8])
                        p_id = attrs.get('ID')
                        if p_id:
                            target_protein_ids.add(p_id)
    
    print(f"Found {len(target_protein_ids)} proteins belonging to target contigs.")

    # 3. Subset FAA based on found Protein IDs
    if os.path.exists(assembly_faa):
        faa_records = []
        for record in SeqIO.parse(assembly_faa, "fasta"):
            # Bakta/Prodigal ID might be the first part of the record ID
            p_id = record.id.split()[0]
            if p_id in target_protein_ids:
                # Ensure contig in description is simplified
                if "[contig=" in record.description:
                    record.description = re.sub(r'\[contig=[^\]]+ctg[^\]]*\]', 
                                              lambda m: m.group(0).split('ctg')[0] + ']', 
                                              record.description)
                faa_records.append(record)
        with open(out_faa, "w") as f_out:
            SeqIO.write(faa_records, f_out, "fasta")

    # 4. Subset TSV (Bakta summary or similar)
    if os.path.exists(assembly_tsv):
        try:
            # Bakta TSV usually has Locus Tag in a column
            df = pd.read_csv(assembly_tsv, sep='\t', comment='#')
            # Look for Locus Tag or ID column
            id_cols = [c for c in df.columns if 'locus' in c.lower() or 'id' in c.lower()]
            if id_cols:
                id_col = id_cols[0]
                subset_df = df[df[id_col].isin(target_protein_ids)]
                subset_df.to_csv(out_tsv, sep='\t', index=False)
            else:
                # Fallback: if we can't find ID, maybe it's subsetted by Contig?
                # But usually TSV is gene-based.
                df.to_csv(out_tsv, sep='\t', index=False)
        except Exception as e:
            print(f"Warning subsetting TSV: {e}")
            with open(out_tsv, 'w') as f: f.write("# Error or no results\n")

    # 5. Subset eggNOG
    if os.path.exists(assembly_eggnog):
        with open(assembly_eggnog, 'r') as f_in, open(out_eggnog, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                    continue
                parts = line.split('\t')
                query = parts[0]
                if query in target_protein_ids:
                    f_out.write(line)

    # 6. Subset Kofam
    if os.path.exists(assembly_kofam):
        with open(assembly_kofam, 'r') as f_in, open(out_kofam, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                    continue
                parts = line.split('\t')
                if len(parts) > 1:
                    query = parts[1] # Kofam format: * gene_id KO ...
                    if query in target_protein_ids:
                        f_out.write(line)
                elif len(parts) > 0:
                    query = parts[0]
                    if query in target_protein_ids:
                        f_out.write(line)

    # 7. Subset dbCAN (directory)
    if os.path.exists(assembly_dbcan):
        if os.path.exists(out_dbcan): shutil.rmtree(out_dbcan)
        os.makedirs(out_dbcan, exist_ok=True)
        for f_name in ["overview.txt", "cgc.out", "dbsub.out", "hmmer.out", "diamond.out"]:
            src = os.path.join(assembly_dbcan, f_name)
            dst = os.path.join(out_dbcan, f_name)
            if os.path.exists(src):
                try:
                    # Check first line for separator
                    with open(src, 'r') as f:
                        header = f.readline()
                    sep = '\t' if '\t' in header else ','
                    
                    df = pd.read_csv(src, sep=sep)
                    id_col = df.columns[0]
                    subset_df = df[df[id_col].isin(target_protein_ids)]
                    subset_df.to_csv(dst, sep=sep, index=False)
                except Exception as e:
                    print(f"Warning subsetting dbCAN {f_name}: {e}")
                    with open(dst, 'w') as f: f.write("# No results\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin_fasta", required=True)
    parser.add_argument("--assembly_gff", required=True)
    parser.add_argument("--assembly_faa", required=True)
    parser.add_argument("--assembly_tsv", required=True)
    parser.add_argument("--assembly_eggnog", required=True)
    parser.add_argument("--assembly_kofam", required=True)
    parser.add_argument("--assembly_dbcan", required=True)
    parser.add_argument("--out_gff", required=True)
    parser.add_argument("--out_faa", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--out_eggnog", required=True)
    parser.add_argument("--out_kofam", required=True)
    parser.add_argument("--out_dbcan", required=True)
    args = parser.parse_args()
    
    subset_results(args.bin_fasta, args.assembly_gff, args.assembly_faa, args.assembly_tsv, 
                   args.assembly_eggnog, args.assembly_kofam, args.assembly_dbcan,
                   args.out_gff, args.out_faa, args.out_tsv, args.out_eggnog, 
                   args.out_kofam, args.out_dbcan)
