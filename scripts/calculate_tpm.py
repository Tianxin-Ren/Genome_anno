#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate TPM/RPKM from featureCounts output (Pandas Version)")
    parser.add_argument("--counts", required=True, help="Input featureCounts TSV")
    parser.add_argument("--summary", required=True, help="Input featureCounts summary file")
    parser.add_argument("--output", required=True, help="Output TSV")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Parse Total Library Size from summary
    summary_df = pd.read_csv(args.summary, sep='\t', index_col=0)
    total_library_reads = summary_df.iloc[:, 0].sum()
    
    if total_library_reads == 0:
        total_library_reads = 1
        print("Warning: Total library reads is 0, setting to 1.")

    # 2. Read featureCounts output
    df = pd.read_csv(args.counts, sep='\t', comment='#')
    
    # 3. Identifiers
    df['LocalGene'] = df['Geneid']
    def extract_mag(gid, chr_str):
        gid_s = str(gid)
        chr_s = str(chr_str)
        if '___' in gid_s:
            return gid_s.split('___')[0]
        if '|' in gid_s:
            return gid_s.split('|')[0]
        if '___' in chr_s:
            return chr_s.split('___')[0]
        if '|' in chr_s:
            return chr_s.split('|')[0]
        return "MAG"
    df['MAG'] = [extract_mag(g, c) for g, c in zip(df['Geneid'], df['Chr'])]
    df['Contig'] = df['Chr'].apply(lambda x: str(x).split(';')[-1])
    # Keep Chr column as string
    df['Chr'] = df['Contig']

    # 4. Calculate RPK, TPM, RPKM
    counts_col = df.columns[6] 
    length_col = 'Length'
    
    # RPK = Counts / (Length / 1000.0)
    df['RPK'] = df[counts_col] / (df[length_col] / 1000.0)
    
    sum_rpk = df['RPK'].sum()
    total_mapped = df[counts_col].sum()
    
    # TPM_Mapped
    if sum_rpk > 0:
        df['TPM_Mapped'] = (df['RPK'] / sum_rpk) * 1000000.0
    else:
        df['TPM_Mapped'] = 0.0
        
    # RPKM_Mapped
    if total_mapped > 0:
        df['RPKM_Mapped'] = df['RPK'] / (total_mapped / 1000000.0)
    else:
        df['RPKM_Mapped'] = 0.0
        
    # RPKM_Lib
    df['RPKM_Lib'] = df['RPK'] / (total_library_reads / 1000000.0)
    
    # TPM_Lib
    mapping_rate = total_mapped / total_library_reads
    df['TPM_Lib'] = df['TPM_Mapped'] * mapping_rate

    # 5. Save final table
    final_cols = ['MAG', 'LocalGene', 'Chr', 'Contig', 'Start', 'End', 'Strand', 'Length', counts_col, 
                  'RPK', 'TPM_Mapped', 'RPKM_Mapped', 'RPKM_Lib', 'TPM_Lib']
    df_final = df[final_cols].rename(columns={counts_col: 'Counts'})
    
    df_final.to_csv(args.output, sep='\t', index=False)
    print(f"Abundance metrics (Pandas) saved to {args.output}")

if __name__ == "__main__":
    main()
