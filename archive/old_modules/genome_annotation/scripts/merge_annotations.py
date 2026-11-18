import pandas as pd
import argparse
import os
import re
import numpy as np

# Fixed schema for Master Table
GENE_COLS = ['Gene_ID', 'Contig', 'Start', 'End', 'Strand', 'Type', 'Gene', 'Product', 'DbXrefs']
EGGNOG_COLS = ['seed_ortholog', 'evalue', 'COG_category', 'Description', 'Preferred_name', 'GOs', 'KEGG_ko']
KOFAM_COLS = ['Kofam_KO', 'Kofam_Definition']
DBCAN_COLS = ['CAZy_HMMER', 'CAZy_dbCAN', 'CAZy_DIAMOND']

def parse_gff(gff_file):
    data = []
    cols = ['Gene_ID', 'Contig', 'Start', 'End', 'Strand', 'Type']
    if not os.path.exists(gff_file): 
        return pd.DataFrame(columns=cols)
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('##FASTA'): break
            if line.startswith('#') or not line.strip(): continue
            parts = line.split('\t')
            if len(parts) < 9: continue
            if parts[2] == 'region': continue
            attr = parts[8]
            gene_id = None
            id_match = re.search(r'ID=([^;]+)', attr)
            if id_match: gene_id = id_match.group(1)
            if gene_id:
                data.append({'Gene_ID': gene_id, 'Contig': parts[0], 'Start': int(parts[3]), 'End': int(parts[4]), 'Strand': parts[6], 'Type': parts[2]})
    
    df = pd.DataFrame(data)
    if df.empty:
        return pd.DataFrame(columns=cols)
    return df

def parse_bakta_tsv(tsv_file):
    cols = ['Gene_ID', 'Gene', 'Product', 'DbXrefs']
    if not os.path.exists(tsv_file) or os.path.getsize(tsv_file) == 0: 
        return pd.DataFrame(columns=cols)
    try:
        df = pd.read_csv(tsv_file, sep='\t', comment='#')
        if 'Locus Tag' in df.columns: 
            df = df.rename(columns={'Locus Tag': 'Gene_ID'})
        elif 'Sequence Id' in df.columns:
            df = df.rename(columns={'Sequence Id': 'Gene_ID'})
            
        # Ensure all required columns exist
        for c in cols:
            if c not in df.columns: df[c] = "-"
            
        return df[cols]
    except:
        return pd.DataFrame(columns=cols)

def parse_eggnog(eggnog_file):
    cols = ['Gene_ID'] + EGGNOG_COLS
    if not eggnog_file or not os.path.exists(eggnog_file): 
        return pd.DataFrame(columns=cols)
    try:
        header = None
        with open(eggnog_file, 'r') as f:
            for line in f:
                if line.startswith('# query'):
                    header = line.lstrip('# ').strip().split('\t')
                    break
        if not header: return pd.DataFrame(columns=cols)
        df = pd.read_csv(eggnog_file, sep='\t', comment='#', names=header)
        df = df.rename(columns={'query': 'Gene_ID'})
        valid_cols = [c for c in cols if c in df.columns]
        res = df[valid_cols]
        for c in cols:
            if c not in res.columns: res[c] = "-"
        return res
    except:
        return pd.DataFrame(columns=cols)

def parse_kofam(kofam_file):
    cols = ['Gene_ID', 'Kofam_KO', 'Kofam_Definition']
    if not kofam_file or not os.path.exists(kofam_file): return pd.DataFrame(columns=cols)
    data = []
    try:
        with open(kofam_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    # Kofam detail-tsv format:
                    # [tag] [gene_id] [ko] [threshold] [score] [e-value] [definition]
                    # 0:tag, 1:gene_id, 2:ko, 3:threshold, 4:score, 5:e-value, 6:definition
                    gene_id = parts[1]
                    ko = parts[2]
                    # Definition is usually at index 6 for 'detail' format
                    if len(parts) > 6:
                        definition = parts[6]
                    elif len(parts) > 3:
                        # Fallback for 'mapper' format: gene_id, ko, definition
                        # In mapper format, parts[2] is definition if len is small
                        definition = parts[2] if len(parts) == 3 else parts[3]
                    else:
                        definition = "-"
                    
                    data.append({'Gene_ID': gene_id, 'Kofam_KO': ko, 'Kofam_Definition': definition.strip()})
        return pd.DataFrame(data)
    except:
        return pd.DataFrame(columns=cols)

def parse_dbcan(dbcan_dir):
    cols = ['Gene_ID'] + DBCAN_COLS
    if not dbcan_dir or not os.path.exists(dbcan_dir): return pd.DataFrame(columns=cols)
    try:
        overview = os.path.join(dbcan_dir, 'overview.txt')
        if not os.path.exists(overview): return pd.DataFrame(columns=cols)
        df = pd.read_csv(overview, sep='\t')
        df = df.rename(columns={df.columns[0]: 'Gene_ID', 'HMMER': 'CAZy_HMMER', 'dbCAN': 'CAZy_dbCAN', 'DIAMOND': 'CAZy_DIAMOND'})
        return df[cols]
    except:
        return pd.DataFrame(columns=cols)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--bakta", required=True)
    parser.add_argument("--eggnog")
    parser.add_argument("--kofam")
    parser.add_argument("--dbcan")
    parser.add_argument("--genomespot")
    parser.add_argument("--checkm2")
    parser.add_argument("--gtdbtk")
    parser.add_argument("--quast")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # 1. Base annotation
    df_final = parse_gff(args.gff)
    if df_final.empty:
        pd.DataFrame(columns=GENE_COLS).to_csv(args.output, sep='\t', index=False)
        return

    # 2. Bakta
    df_bakta = parse_bakta_tsv(args.bakta)
    df_final = df_final.merge(df_bakta, on='Gene_ID', how='left')

    # 3. Functional Modules
    if args.eggnog:
        df_final = df_final.merge(parse_eggnog(args.eggnog), on='Gene_ID', how='left')
    if args.kofam:
        df_final = df_final.merge(parse_kofam(args.kofam), on='Gene_ID', how='left')
    if args.dbcan:
        df_final = df_final.merge(parse_dbcan(args.dbcan), on='Gene_ID', how='left')

    # 4. Add Genome-level info as metadata columns if needed (Optional)
    # GenomeSpot, CheckM2, GTDB-Tk info can be merged as constant columns or handled separately
    
    df_final.fillna("-").to_csv(args.output, sep='\t', index=False)
    print(f"Master Table saved to {args.output}")

if __name__ == "__main__":
    main()
