import os
import argparse
import pandas as pd
from Bio import SeqIO

def get_attributes(attr_str):
    attrs = {}
    for item in attr_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs

def subset_functional(bin_fasta, assembly_gff, assembly_faa, assembly_tsv, 
                      assembly_eggnog, assembly_kofam, 
                      out_gff, out_faa, out_tsv, out_eggnog, out_kofam):
    
    # 1. Get exact contig IDs from the bin
    bin_headers = set()
    for record in SeqIO.parse(bin_fasta, "fasta"):
        bin_headers.add(record.id)
    
    # Also support simplified IDs (splitting at 'ctg') for flexibility
    simple_headers = {h.split('ctg')[0] for h in bin_headers}
    
    print(f"Targeting {len(bin_headers)} contigs based on exact ID matches.")

    # 2. Collect Protein IDs and subset GFF
    target_protein_ids = set()
    with open(assembly_gff, 'r') as f_in, open(out_gff, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                continue
            
            parts = line.split('\t')
            if len(parts) >= 9:
                contig_name = parts[0]
                contig_simple = contig_name.split('ctg')[0]
                
                # Direct match or simple ID match
                if contig_name in bin_headers or contig_simple in simple_headers:
                    f_out.write(line)
                    attrs = get_attributes(parts[8])
                    p_id = attrs.get('ID')
                    if p_id: target_protein_ids.add(p_id)
                    lt = attrs.get('locus_tag')
                    if lt: target_protein_ids.add(lt)
    
    print(f"Found {len(target_protein_ids)} protein IDs linked to these contigs.")

    # 3. Subset FAA
    if os.path.exists(assembly_faa):
        faa_records = []
        for record in SeqIO.parse(assembly_faa, "fasta"):
            p_id = record.id.split()[0]
            if p_id in target_protein_ids:
                faa_records.append(record)
        with open(out_faa, "w") as f_out:
            SeqIO.write(faa_records, f_out, "fasta")

    # 4. Subset TSV
    if os.path.exists(assembly_tsv):
        try:
            df = pd.read_csv(assembly_tsv, sep='\t', comment='#')
            lt_col = [c for c in df.columns if 'locus' in c.lower() or 'id' in c.lower()]
            if lt_col:
                subset_df = df[df[lt_col[0]].isin(target_protein_ids)]
                subset_df.to_csv(out_tsv, sep='\t', index=False)
        except: pass

    # 5. Subset eggNOG
    if os.path.exists(assembly_eggnog):
        with open(assembly_eggnog, 'r') as f_in, open(out_eggnog, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                    continue
                parts = line.split('\t')
                if parts[0] in target_protein_ids:
                    f_out.write(line)

    # 6. Subset Kofam
    if os.path.exists(assembly_kofam):
        with open(assembly_kofam, 'r') as f_in, open(out_kofam, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                    continue
                parts = line.split('\t')
                query = parts[1] if len(parts) > 1 else parts[0]
                if query in target_protein_ids:
                    f_out.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin_fasta", required=True)
    parser.add_argument("--assembly_gff", required=True)
    parser.add_argument("--assembly_faa", required=True)
    parser.add_argument("--assembly_tsv", required=True)
    parser.add_argument("--assembly_eggnog", required=True)
    parser.add_argument("--assembly_kofam", required=True)
    parser.add_argument("--out_gff", required=True)
    parser.add_argument("--out_faa", required=True)
    parser.add_argument("--out_tsv", required=True)
    parser.add_argument("--out_eggnog", required=True)
    parser.add_argument("--out_kofam", required=True)
    args = parser.parse_args()
    
    subset_functional(args.bin_fasta, args.assembly_gff, args.assembly_faa, args.assembly_tsv, 
                      args.assembly_eggnog, args.assembly_kofam, 
                      args.out_gff, args.out_faa, args.out_tsv, args.out_eggnog, args.out_kofam)
