#!/usr/bin/env python3
import sys
import os
import pandas as pd

def evaluate_oxphos(input_path, output_path):
    print(f"正在读取输入文件: {input_path} ...")
    # 支持 TSV 格式读取
    df_input = pd.read_csv(input_path, sep="\t")
    
    results = []
    columns = df_input.columns[1:]
    
    # 建立基因名清理与匹配映射
    def get_presence_dict(row):
        presence = {}
        for col in columns:
            # 提取标准的基因简称，例如 "nuoA; NADH-..." -> "nuoA"
            gene_id = col.split(';')[0].split(',')[0].strip().split(' ')[0]
            # Multiple KOs can map to the same canonical gene name. Preserve a
            # positive hit instead of allowing a later zero-valued column to
            # overwrite it.
            is_present = 1 if row[col] > 0 else 0
            presence[gene_id] = max(presence.get(gene_id, 0), is_present)
        return presence

    print("开始针对 Step 1 - Step 5 进行定量评估...")
    for index, row in df_input.iterrows():
        genome = row.iloc[0]
        p = get_presence_dict(row)
        
        # ----------------------------------------------------
        # STEP 1 EVALUATION (NADH 脱氢酶 / 复合物 I)
        # ----------------------------------------------------
        nuo_subunits = ['nuoA', 'nuoB', 'nuoC', 'nuoD', 'nuoE', 'nuoF', 'nuoG', 'nuoH', 'nuoI', 'nuoJ', 'nuoK', 'nuoL', 'nuoM', 'nuoN', 'nuoCD', 'nuoBCD', 'nuoLM']
        nuo_count = sum([p.get(g, 0) for g in nuo_subunits])
        if p.get('nuoBCD', 0): nuo_count += 2
        if p.get('nuoCD', 0): nuo_count += 1
        if p.get('nuoLM', 0): nuo_count += 1
        nuo_score = min(nuo_count / 14.0, 1.0)
        
        ndh_score = 1.0 if (p.get('ndh', 0) or any(p.get(f'ndh{c}', 0) for c in 'ABCDEFGHIJKLMN')) else 0.0
        hox_score = sum([p.get('hoxE', 0), p.get('hoxF', 0), p.get('hoxU', 0)]) / 3.0
        
        step1_score = max(nuo_score, ndh_score, hox_score)
        if step1_score >= 0.8:   step1_status = 'Full'
        elif step1_score >= 0.3: step1_status = 'Partial'
        else:                    step1_status = 'Absent'
            
        # ----------------------------------------------------
        # STEP 2 EVALUATION (琥珀酸脱氢酶 / 复合物 II)
        # ----------------------------------------------------
        pos1 = max(p.get('sdhA',0), p.get('frdA',0), p.get('SDHA',0))
        pos2 = max(p.get('sdhB',0), p.get('frdB',0), p.get('SDHB',0))
        pos3 = max(p.get('sdhC',0), p.get('frdC',0), p.get('SDHC',0))
        pos4 = max(p.get('sdhD',0), p.get('frdD',0), p.get('SDHD',0), p.get('SDH4',0))
        step2_score = (pos1 + pos2 + pos3 + pos4) / 4.0
        if step2_score >= 0.75:  step2_status = 'Full'
        elif step2_score >= 0.25: step2_status = 'Partial'
        else:                     step2_status = 'Absent'

        # ----------------------------------------------------
        # STEP 3 EVALUATION (细胞色素 bc1 / 复合物 III)
        # ----------------------------------------------------
        bc1_core = [max(p.get('qcrA',0), p.get('UQCRFS1',0)), 
                    max(p.get('qcrB',0), p.get('CYTB',0), p.get('fbcH',0)), 
                    max(p.get('qcrC',0), p.get('CYC1',0))]
        step3_score = sum(bc1_core) / 3.0
        if step3_score >= 0.66:  step3_status = 'Full'
        elif step3_score >= 0.33: step3_status = 'Partial'
        else:                     step3_status = 'Absent'

        # ----------------------------------------------------
        # STEP 4 EVALUATION (末端氧化酶 / 复合物 IV)
        # ----------------------------------------------------
        cox_type = (max(p.get('COX1',0), p.get('coxA',0)) + max(p.get('COX2',0), p.get('coxB',0))) / 2.0
        cco_type = (p.get('ccoN', 0) + p.get('ccoO', 0) + p.get('ccoNO', 0)) / 2.0
        cyd_type = (p.get('cydA', 0) + p.get('cydB', 0)) / 2.0
        cyo_type = (p.get('cyoB', 0) + p.get('cyoA', 0)) / 2.0
        sox_type = (p.get('soxB', 0) + p.get('soxA', 0)) / 2.0
        
        step4_score = max(cox_type, cco_type, cyd_type, cyo_type, sox_type)
        step4_score = min(step4_score, 1.0)
        if step4_score >= 1.0:   step4_status = 'Full'
        elif step4_score >= 0.5: step4_status = 'Partial'
        else:                    step4_status = 'Absent'

        # ----------------------------------------------------
        # STEP 5 EVALUATION (ATP 合成酶 / 复合物 V)
        # ----------------------------------------------------
        f_type_subunits = ['ATPF0A', 'ATPF0B', 'ATPF0C', 'ATPF1A', 'ATPF1B', 'ATPF1D', 'ATPF1E', 'ATPF1G']
        f_type_score = sum([p.get(g, 0) for g in f_type_subunits]) / 8.0
        
        va_type_subunits = ['ATPVA', 'ATPVB', 'ATPVC', 'ATPVD', 'ATPVE', 'ATPVF', 'ATPVG', 'ATPVI', 'ATPVK']
        va_type_score = sum([p.get(g, 0) for g in va_type_subunits]) / 9.0
        
        step5_score = max(f_type_score, va_type_score)
        if step5_score >= 0.75:  step5_status = 'Full'
        elif step5_score >= 0.35: step5_status = 'Partial'
        else:                     step5_status = 'Absent'
        
        results.append({
            'Genome_Name': genome,
            'Step1_Score': round(step1_score, 2), 'Step1_Status': step1_status,
            'Step2_Score': round(step2_score, 2), 'Step2_Status': step2_status,
            'Step3_Score': round(step3_score, 2), 'Step3_Status': step3_status,
            'Step4_Score': round(step4_score, 2), 'Step4_Status': step4_status,
            'Step5_Score': round(step5_score, 2), 'Step5_Status': step5_status,
        })
        
    df_out = pd.DataFrame(results)
    output_dir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_dir, exist_ok=True)
    df_out.to_csv(output_path, sep="\t", index=False)
    print(f"评估完成！结果已成功写入 TSV 文件: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("使用错误！正确运行命令格式为:")
        print("python evaluate_oxphos.py <输入注释结果.tsv> <输出评估结果.tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"错误：找不到输入文件 {input_file}")
        sys.exit(1)
        
    evaluate_oxphos(input_file, output_file)
