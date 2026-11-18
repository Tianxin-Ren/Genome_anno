import pandas as pd
import argparse
import os

MISSING_VAL = "Not_Run"

def main():
    parser = argparse.ArgumentParser(description="Generate Global Summary with Not_Run logic")
    parser.add_argument("--checkm2", required=False)
    parser.add_argument("--gtdbtk", required=False)
    parser.add_argument("--quast-reports", nargs='*', required=False)
    parser.add_argument("--genomespot-reports", nargs='*', required=False)
    parser.add_argument("--genomes", nargs='+', required=True)
    parser.add_argument("--output", required=True)
    
    args = parser.parse_args()
    
    # Initialize with all genomes
    summary = pd.DataFrame({'Genome': args.genomes})
    
    # 1. Load CheckM2
    if args.checkm2 and os.path.exists(args.checkm2):
        df_qc = pd.read_csv(args.checkm2, sep='\t')
        df_qc['Genome'] = df_qc['Name'].apply(lambda x: x.split('.')[0])
        df_qc = df_qc[['Genome', 'Completeness', 'Contamination']].rename(columns={'Completeness': 'QC_Completeness', 'Contamination': 'QC_Contamination'})
        summary = summary.merge(df_qc, on='Genome', how='left')
    else:
        for c in ['QC_Completeness', 'QC_Contamination']: summary[c] = MISSING_VAL
    
    # 2. Load GTDB-Tk
    if args.gtdbtk and os.path.exists(args.gtdbtk) and os.path.getsize(args.gtdbtk) > 0:
        df_tax = pd.read_csv(args.gtdbtk, sep='\t')
        df_tax['Genome'] = df_tax['user_genome'].apply(lambda x: str(x).split('.')[0])
        
        # Robustly handle potentially missing columns (like fastani_ani in skip-ANI mode)
        keep_cols = ['Genome', 'classification']
        rename_dict = {'classification': 'Tax_Classification'}
        if 'fastani_ani' in df_tax.columns:
            keep_cols.append('fastani_ani')
            rename_dict['fastani_ani'] = 'Tax_ANI'
        else:
            df_tax['Tax_ANI'] = MISSING_VAL
            
        df_tax = df_tax[keep_cols].rename(columns=rename_dict)
        summary = summary.merge(df_tax, on='Genome', how='left')
    else:
        for c in ['Tax_Classification', 'Tax_ANI']: summary[c] = MISSING_VAL
    
    # 3. Integrate QUAST
    quast_data = []
    if args.quast_reports:
        for q_dir in args.quast_reports:
            if not os.path.exists(q_dir): continue
            # Extract genome name from path results/genome_annotation/genome/01_qc/quast
            genome = os.path.basename(os.path.dirname(os.path.dirname(q_dir)))
            report = os.path.join(q_dir, "report.tsv")
            if os.path.exists(report):
                d = {'Genome': genome}
                with open(report, 'r') as f:
                    for line in f:
                        p = line.strip().split('\t')
                        if len(p) < 2: continue
                        if p[0] == 'Total length': d['Assm_Total_Length'] = p[1]
                        elif p[0] == 'N50': d['Assm_N50'] = p[1]
                        elif p[0] == 'GC (%)': d['Assm_GC_Percent'] = p[1]
                quast_data.append(d)
    
    if quast_data:
        df_quast = pd.DataFrame(quast_data)
        summary = summary.merge(df_quast, on='Genome', how='left')
    else:
        for c in ['Assm_Total_Length', 'Assm_N50', 'Assm_GC_Percent']: summary[c] = MISSING_VAL

    # 4. Integrate GenomeSpot OGT
    gs_data = []
    if args.genomespot_reports:
        mapping = {
            'optimum_temperature': 'Pheno_OGT',
            'ogt': 'Pheno_OGT',
            'optimum_salinity': 'Pheno_Salinity',
            'salinity': 'Pheno_Salinity'
        }
        for gs_file in args.genomespot_reports:
            if not os.path.exists(gs_file): continue
            genome = os.path.basename(gs_file).split('.')[0]
            try:
                df_gs = pd.read_csv(gs_file, sep='\t')
                if len(df_gs) >= 1:
                    d = {'Genome': genome}
                    for old_col, new_col in mapping.items():
                        matching = [c for c in df_gs.columns if old_col in str(c).lower()]
                        if matching:
                            d[new_col] = df_gs[matching[0]].iloc[0]
                    gs_data.append(d)
            except:
                pass
            
    if gs_data:
        df_gs_all = pd.DataFrame(gs_data)
        summary = summary.merge(df_gs_all, on='Genome', how='left')
    else:
        for c in ['Pheno_OGT', 'Pheno_Salinity']: summary[c] = MISSING_VAL

    # Final cleanup
    summary = summary.fillna(MISSING_VAL)
    summary.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
