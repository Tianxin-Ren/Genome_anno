#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import sys
from multiprocessing import Pool

# 导入配置
try:
    import config
except ImportError:
    print("[-] 错误: 找不到 config.py 配置文件。")
    sys.exit(1)

import logging
logger = logging.getLogger("Pipeline")

def build_rpob_mapping(fasta_path):
    mapping = {}
    if not os.path.exists(fasta_path):
        logger.warning(f"未找到 FASTA 文件 {fasta_path}，无法进行 ID 还原映射。")
        return mapping
        
    logger.info(f"正在从 {fasta_path} 构建 ID 映射字典...")
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split(maxsplit=1)
                seqid = parts[0][1:]
                desc = parts[1].lower() if len(parts) > 1 else ""
                
                if 'nrfa' in desc or 'nitrite reductase' in desc:
                    mapping[seqid] = 'nrfA'
                else:
                    mapping[seqid] = 'rpoB'
    return mapping

def process_single_sample(args):
    """
    单样本处理函数，供多进程调用
    """
    sample, rpob_map = args
    nif_file = os.path.join(config.OUT_DIR, f"{sample}__combined_db_nif_r207__result_bsr.tsv")
    rpob_file = os.path.join(config.OUT_DIR, f"{sample}__rpoB_final__result_bsr.tsv")

    sample_dfs = []
    files_found = 0
    
    for fpath in [nif_file, rpob_file]:
        if not os.path.exists(fpath):
            continue
        files_found += 1
        try:
            # 仅读取必要的列
            df = pd.read_csv(fpath, sep='\t', usecols=['sseqid_db', 'score_db', 'bsr'])
            # 应用全局配置中的阈值
            df_filtered = df[(df['score_db'] >= config.MIN_SCORE) & (df['bsr'] >= config.MIN_BSR)].copy()
            
            if not df_filtered.empty:
                if "nif_r207" in fpath:
                    # NIF 库按下划线提取基因名
                    df_filtered['GeneFamily'] = df_filtered['sseqid_db'].apply(lambda x: str(x).split('_')[0])
                else:
                    # RPOB 库使用字典映射
                    df_filtered['GeneFamily'] = df_filtered['sseqid_db'].map(rpob_map).fillna('rpoB')
                
                # 仅保留目标列表中的基因
                target_df = df_filtered[df_filtered['GeneFamily'].isin(config.TARGET_GENES)]
                sample_dfs.append(target_df)
        except Exception:
            continue

    if sample_dfs:
        combined = pd.concat(sample_dfs)
        counts = combined['GeneFamily'].value_counts().reset_index()
        counts.columns = ['GeneFamily', 'Count']
        counts['SampleID'] = sample
        return counts
    elif files_found == 0:
        return f"MISSING:{sample}"
    else:
        return f"EMPTY:{sample}"

def main():
    logger.info("-" * 40)
    logger.info("模块：后期整合统计处理")
    logger.info("-" * 40)
    
    # 1. 扫描样本
    if config.TARGET_SAMPLES:
        samples = config.TARGET_SAMPLES
    else:
        if not os.path.exists(config.METAGENOMICS_DIR):
            logger.error(f"输入目录 {config.METAGENOMICS_DIR} 不存在。")
            return
        # 扫描比对结果目录
        files = [f for f in os.listdir(config.OUT_DIR) if f.endswith("__result_bsr.tsv")]
        samples = sorted(list(set([f.split("__")[0] for f in files])))
    
    if not samples:
        logger.error("未发现待处理样本。请检查 config.py 或输出目录。")
        return

    logger.info(f"待处理样本总数: {len(samples)}")
    
    # 2. 并行提取矩阵
    rpob_map = build_rpob_mapping(config.RPOB_FASTA)
    
    logger.info(f"正在并行处理各样本 (Multiprocessing)...")
    with Pool() as pool:
        results = pool.map(process_single_sample, [(s, rpob_map) for s in samples])

    all_counts = []
    missing = []
    empty = []
    
    for res in results:
        if isinstance(res, pd.DataFrame):
            all_counts.append(res)
        elif isinstance(res, str):
            if res.startswith("MISSING"): missing.append(res.split(":")[1])
            else: empty.append(res.split(":")[1])

    if missing: logger.warning(f"有 {len(missing)} 个样本未找到结果文件。")
    if empty: logger.warning(f"有 {len(empty)} 个样本在过滤后无有效基因。")

    if not all_counts:
        logger.error("所有样本均无有效结果数据。")
        return

    # 3. 生成原始矩阵
    logger.info("正在聚合数据生成原始矩阵...")
    raw_df = pd.concat(all_counts, ignore_index=True)
    matrix = raw_df.pivot_table(index='SampleID', columns='GeneFamily', values='Count', aggfunc='sum').fillna(0).astype(int)
    
    # 强制对齐目标基因列
    for gene in config.TARGET_GENES:
        if gene not in matrix.columns: matrix[gene] = 0
    matrix = matrix[config.TARGET_GENES]
    
    raw_out = os.path.join(config.OUT_DIR, 'targeted_gene_matrix_raw.tsv')
    matrix.to_csv(raw_out, sep='\t')
    logger.info(f"原始矩阵已保存: {raw_out}")

    # 4. 长度标准化 (Step 1)
    df_step1 = matrix.copy().astype(float)
    for gene in df_step1.columns:
        if gene in config.GENE_LENGTHS:
            df_step1[gene] = df_step1[gene] / (config.GENE_LENGTHS[gene] * 3)
    
    step1_out = os.path.join(config.OUT_DIR, 'matrix_step1_length_normalized.tsv')
    df_step1.round(6).to_csv(step1_out, sep='\t')
    logger.info(f"长度标准化矩阵已保存: {step1_out}")

    # 5. rpoB 标准化 (Step 2)
    if 'rpoB' in df_step1.columns:
        rpob_vals = df_step1['rpoB']
        df_step2 = df_step1.div(rpob_vals, axis=0) * 100
        df_step2 = df_step2.replace([np.inf, -np.inf], np.nan).fillna(0)
        
        step2_out = os.path.join(config.OUT_DIR, 'matrix_step2_rpob_normalized.tsv')
        df_step2.round(4).to_csv(step2_out, sep='\t')
        logger.info(f"rpoB 标准化矩阵已保存: {step2_out}")
    else:
        logger.warning("未找到 rpoB 列，跳过相对丰度标准化。")

    logger.info("后期整合处理顺利完成。")

if __name__ == "__main__":
    main()
