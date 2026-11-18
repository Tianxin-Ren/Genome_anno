import pandas as pd
import sys
import os
import re

def aggregate_functional(faa_path, egg_path, kofam_path, advanced_hmm_path, gff_path, sample_name, alias, output_path):
    print(f"[agg] Processing sample: {sample_name} (Alias: {alias})", file=sys.stderr)
    
    # 1. 从 GFF 提取基础坐标信息 (Contig, Start, End, Strand)
    # ... (keeping existing logic)
    gff_data = []
    if os.path.exists(gff_path):
        with open(gff_path) as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.split('\t')
                if len(parts) == 9 and parts[2] == 'CDS':
                    attr = parts[8]
                    m_id = re.search(r'ID=([^;]+)', attr)
                    if m_id:
                        clean_id = m_id.group(1)
                        gff_data.append({
                            'Protein_ID': clean_id,
                            'Contig': parts[0],
                            'Start': parts[3],
                            'End': parts[4],
                            'Strand': parts[6]
                        })
    df_base = pd.DataFrame(gff_data)
    if df_base.empty:
        print(f"WARNING: No CDS entries found in GFF: {gff_path}. Writing empty summary table.", file=sys.stderr)
        df_empty = pd.DataFrame(columns=[
            'Name', 'Alias', 'Protein_ID', 'Contig', 'Start', 'End', 'Strand', 'Gene_Index_on_Contig',
            'Score_Annotation', 'Evalue_Annotation', 'Best_Evalue', 'Target_Annotation', 'Pfam_Domains',
            'All_Similar_Hits', 'Decision_Logic', 'EggNOG_Description', 'EggNOG_Name', 'EggNOG_KO',
            'EggNOG_Pfam', 'KOfam_KO', 'KOfam_Hit'
        ])
        df_empty.to_csv(output_path, sep='\t', index=False)
        return

    # 新增：计算在每个 Contig 上的基因序号 (从 1 开始)
    df_base['Gene_Index_on_Contig'] = df_base.groupby('Contig').cumcount() + 1

    # 2. 读取高级 HMM 仲裁结果 (优先放置)
    if advanced_hmm_path != "none" and os.path.exists(advanced_hmm_path) and os.path.getsize(advanced_hmm_path) > 0:
        try:
            df_hmm = pd.read_csv(advanced_hmm_path, sep='\t')
            if 'Protein' in df_hmm.columns:
                df_hmm = df_hmm.rename(columns={'Protein': 'Protein_ID'})
            if 'Protein_ID' in df_hmm.columns:
                df_base = pd.merge(df_base, df_hmm, on='Protein_ID', how='left')
        except Exception as e:
            print(f"Warning: Failed to load advanced HMM results: {e}", file=sys.stderr)

    # 3. 读取 EggNOG (仅保留 KO, Pfam, Description, Preferred_name)
    if egg_path != "none" and os.path.exists(egg_path) and os.path.getsize(egg_path) > 0:
        try:
            df_egg = pd.read_csv(egg_path, sep='\t', comment='#', header=None)
            if not df_egg.empty:
                # EggNOG-mapper v2 columns (0-indexed): 
                # 0:query, 7:Description, 8:Preferred_name, 11:KEGG_ko, 20:PFAMs
                cols_to_idx = {
                    'Protein_ID': 0,
                    'EggNOG_Description': 7,
                    'EggNOG_Name': 8,
                    'EggNOG_KO': 11,
                    'EggNOG_Pfam': 20
                }
                
                max_col = df_egg.shape[1]
                valid_indices = [idx for name, idx in cols_to_idx.items() if idx < max_col]
                valid_names = [name for name, idx in cols_to_idx.items() if idx < max_col]
                
                df_egg = df_egg.iloc[:, valid_indices]
                df_egg.columns = valid_names
                if 'Protein_ID' in df_egg.columns:
                    df_base = pd.merge(df_base, df_egg, on='Protein_ID', how='left')
        except Exception as e:
            print(f"Warning: Failed to load EggNOG annotations: {e}", file=sys.stderr)

    # 4. 读取 KofamScan (如果有)
    if kofam_path != "none" and os.path.exists(kofam_path):
        kofam_data = []
        with open(kofam_path) as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split()
                if not parts: continue
                
                # KofamScan detail format starts with '*' for hits above threshold, or space/nothing for others
                if parts[0] == '*':
                    if len(parts) >= 3:
                        gene_id = parts[1]
                        ko_id = parts[2]
                        kofam_data.append({'Protein_ID': gene_id, 'KOfam_KO': ko_id, 'KOfam_Hit': True})
                else:
                    if len(parts) >= 2:
                        gene_id = parts[0]
                        ko_id = parts[1]
                        kofam_data.append({'Protein_ID': gene_id, 'KOfam_KO': ko_id, 'KOfam_Hit': False})
                        
        if kofam_data:
            df_ko = pd.DataFrame(kofam_data).drop_duplicates('Protein_ID')
            df_base = pd.merge(df_base, df_ko, on='Protein_ID', how='left')
        else:
            # 物理上 kofam 结果没有命中任何基因，所以不含有有效行。我们直接手动补齐这一空列以保持表头结构一致！
            df_base['KOfam_KO'] = ""
            df_base['KOfam_Hit'] = False

    # 5. 插入 Name 和 Alias 列
    df_base.insert(0, 'Alias', alias)
    df_base.insert(0, 'Name', sample_name)
    
    # 保存结果
    df_base.to_csv(output_path, sep='\t', index=False)
    print(f"[agg] Final table written to: {output_path} (Total proteins: {len(df_base)})", file=sys.stderr)

if __name__ == "__main__":
    # 参数顺序由 functional.smk 决定
    # python3 scripts/aggregate_functional.py {faa} {egg} {kofam} {advanced_hmm} {gff} {sample} {alias} {output}
    if len(sys.argv) < 9:
        print("Usage: python3 aggregate_functional.py <faa> <egg> <kofam> <adv_hmm> <gff> <sample> <alias> <output>")
        sys.exit(1)
        
    aggregate_functional(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
