#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate TPM/RPKM from featureCounts output")
    parser.add_argument("--counts", required=True, help="Input featureCounts TSV")
    parser.add_argument("--summary", required=True, help="Input featureCounts summary file")
    parser.add_argument("--output", required=True, help="Output TSV")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # 1. Parse Total Library Size from summary
    # summary file format: Status\tSampleName
    # Row 1: Assigned\t123
    # Row 2: Unassigned_Ambiguity\t45
    # ...
    summary_df = pd.read_csv(args.summary, sep='\t', index_col=0)
    total_library_reads = summary_df.iloc[:, 0].sum()
    
    if total_library_reads == 0:
        total_library_reads = 1
        print("Warning: Total library reads is 0, setting to 1 to avoid division by zero.")

    # 2. Read featureCounts output
    # Skip the first comment line
    df = pd.read_csv(args.counts, sep='\t', comment='#')
    
    # 3. Parse MAG and LocalGene information
    # Column 1 (Geneid): e.g. "5_2"
    # Column 2 (Chr): e.g. "Bins.sr.10|Bins.sr.10_contig_5"
    
    def parse_row(row):
        gene_id = str(row['Geneid'])
        chr_field = str(row['Chr'])
        
        # Split chr_field by "|"
        chr_parts = chr_field.split('|')
        mag = chr_parts[0]
        contig_part = chr_parts[1] if len(chr_parts) > 1 else chr_field
        
        # Extract gene number (last part after _)
        gene_num = gene_id.split('_')[-1]
        
        local_gene = f"{contig_part}_{gene_num}"
        return mag, local_gene

    info_df = df.apply(parse_row, axis=1, result_type='expand')
    df.insert(0, 'MAG', info_df[0])
    df.insert(1, 'LocalGene', info_df[1])

    # 4. Calculate RPK, TPM, RPKM
    # Column mapping (assuming single sample BAM input to featureCounts)
    # featureCounts cols: MAG, LocalGene, Geneid, Chr, Start, End, Strand, Length, Counts
    # index of Counts is 8
    counts_col = df.columns[8]
    length_col = 'Length'
    
    # RPK = Counts / (Length / 1000)
    df['RPK'] = df[counts_col] / (df[length_col] / 1000.0)
    
    sum_rpk = df['RPK'].sum()
    total_mapped = df[counts_col].sum()
    
    # TPM_Mapped = (RPK / sum_RPK) * 1,000,000
    if sum_rpk > 0:
        df['TPM_Mapped'] = (df['RPK'] / sum_rpk) * 1000000.0
    else:
        df['TPM_Mapped'] = 0.0
        
    # RPKM_Mapped = (Counts / (Length/1000)) / (total_mapped / 1,000,000)
    if total_mapped > 0:
        df['RPKM_Mapped'] = df['RPK'] / (total_mapped / 1000000.0)
    else:
        df['RPKM_Mapped'] = 0.0
        
    # RPKM_Lib = (Counts / (Length/1000)) / (total_library_reads / 1,000,000)
    df['RPKM_Lib'] = df['RPK'] / (total_library_reads / 1000000.0)
    
    # TPM_Lib = TPM_Mapped * (total_mapped / total_library_reads)
    mapping_rate = total_mapped / total_library_reads
    df['TPM_Lib'] = df['TPM_Mapped'] * mapping_rate

    # 5. Save final table
    # Standardize columns to match user's final.tsv expectation
    final_cols = ['MAG', 'LocalGene', 'Chr', 'Start', 'End', 'Strand', 'Length', counts_col, 
                  'RPK', 'TPM_Mapped', 'RPKM_Mapped', 'RPKM_Lib', 'TPM_Lib']
    df_final = df[final_cols].rename(columns={counts_col: 'Counts'})
    
    df_final.to_csv(args.output, sep='\t', index=False)
    print(f"Abundance metrics saved to {args.output}")

if __name__ == "__main__":
    main()
