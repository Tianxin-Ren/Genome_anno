#!/usr/bin/env python3
# ==============================================================================
# GTDB 细菌参考基因组按需分层下载与操纵子检出工具 (fetch_bacterial_references.py)
# 作用：
#   1. 根据待测 MAG 分类学，按亲缘阶元（同属 -> 同科 -> 同目 -> 同纲）分层构建候选代表菌池；
#   2. 对多个目标基因组的亲缘关系进行全局去冗余，统一集中下载到同一目录，杜绝重复下载；
#   3. 逐个下载并现场调用 Diamond 比对检查目标 operon（如 sox 探针）；
#   4. 检出完整 operon 的确立为阳性参考，无目标 operon 的确立为近缘阴性参考；
#   5. 达到设置的参考数量要求（--max_positive_refs, --max_negative_refs）立即停止下载，绝不多下。
# ==============================================================================
import os
import sys
import argparse
import subprocess
import shutil
import zipfile
import re
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Tiered, deduplicated on-demand fetcher & operon validator for bacterial reference genomes."
    )
    parser.add_argument(
        "--gtdb_metadata",
        default="/dssg/home/acct-clswfp/clswfp-icdli1/Database/GTDB/archive/metadata/232.0/bac120_metadata_r232.tsv",
        help="Path to GTDB bac120 metadata TSV"
    )
    parser.add_argument(
        "--target_classification",
        default="",
        help="Path to mag_classification file (MAG_id \\t taxonomy)"
    )
    parser.add_argument(
        "--target_taxonomy",
        default="",
        help="Fallback target taxonomy string if no classification file provided"
    )
    parser.add_argument(
        "--query_operon_faa",
        default="",
        help="Path to query operon FAA (e.g. Reinekea_forsetii_sox_operon.faa)"
    )
    parser.add_argument(
        "--seed_positive_accs",
        nargs="*",
        default=["GCF_000223525.1"],
        help="Known ground-truth positive reference accessions (e.g. Reinekea forsetii GCF_000223525.1)"
    )
    parser.add_argument(
        "--out_genomes_dir",
        default="results/synteny_validation/database/bac120_references/genomes",
        help="Unified destination directory for reference .fna files"
    )
    parser.add_argument(
        "--out_proteins_dir",
        default="results/synteny_validation/database/bac120_references/proteins",
        help="Unified destination directory for reference .faa files"
    )
    parser.add_argument(
        "--out_manifest",
        default="results/synteny_validation/database/bac120_references/downloaded_references_manifest.tsv",
        help="Path to write reference summary manifest TSV"
    )
    parser.add_argument(
        "--max_positive_refs",
        type=int,
        default=3,
        help="Maximum positive reference genomes to acquire (default: 3)"
    )
    parser.add_argument(
        "--max_negative_refs",
        type=int,
        default=2,
        help="Maximum negative sister reference genomes to acquire (default: 2)"
    )
    parser.add_argument(
        "--positive_gene_threshold",
        type=int,
        default=4,
        help="Minimum number of operon genes matched to qualify as positive reference (default: 4)"
    )
    parser.add_argument(
        "--conda_env",
        default="bio",
        help="Conda environment name for datasets/diamond/prodigal"
    )
    return parser.parse_args()

def extract_clean_acc(acc_str):
    s = str(acc_str).strip()
    s = s.replace("GB_", "").replace("RS_", "")
    m = re.search(r"(GC[AF]_\d+\.\d+)", s)
    return m.group(1) if m else s

def download_genome_via_datasets(acc, out_fna, out_faa, env_name="bio", tmp_dir="/tmp"):
    clean = extract_clean_acc(acc)
    zip_path = os.path.join(tmp_dir, f"{clean}_dataset.zip")
    extract_dir = os.path.join(tmp_dir, f"{clean}_extracted")
    
    cmd = [
        "datasets", "download", "genome", "accession", clean,
        "--include", "genome,protein",
        "--filename", zip_path
    ]
    if env_name:
        cmd = ["conda", "run", "-n", env_name] + cmd
        
    print(f"  [NCBI Datasets] 正在下载 {clean} ...", flush=True)
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=300)
        if res.returncode == 0 and os.path.exists(zip_path) and os.path.getsize(zip_path) > 100:
            if os.path.exists(extract_dir):
                shutil.rmtree(extract_dir)
            with zipfile.ZipFile(zip_path, 'r') as zf:
                zf.extractall(extract_dir)
            
            found_fna, found_faa = "", ""
            for root, _, files in os.walk(extract_dir):
                for f in files:
                    if f.endswith(".fna") and not f.startswith("cds_") and not f.startswith("rna_"):
                        found_fna = os.path.join(root, f)
                    elif f.endswith("protein.faa") or f.endswith(".faa"):
                        found_faa = os.path.join(root, f)
            
            if found_fna:
                shutil.copyfile(found_fna, out_fna)
            if found_faa:
                shutil.copyfile(found_faa, out_faa)
            
            if os.path.exists(zip_path): os.remove(zip_path)
            if os.path.exists(extract_dir): shutil.rmtree(extract_dir)
            
            if os.path.exists(out_fna) and os.path.getsize(out_fna) > 0:
                print(f"  [+] 下载成功: {clean} -> {os.path.basename(out_fna)}", flush=True)
                return True
    except Exception as e:
        print(f"  [-] NCBI datasets 下载提示: {clean} 失败 ({e})", file=sys.stderr)
        
    return False

def predict_proteins_if_needed(fna_path, faa_path, env_name="bio"):
    if os.path.exists(faa_path) and os.path.getsize(faa_path) > 100:
        return True
    print(f"  [Prodigal] 为 {os.path.basename(fna_path)} 现场预测蛋白质序列...", flush=True)
    cmd = ["prodigal", "-i", fna_path, "-a", faa_path, "-p", "single", "-q"]
    if env_name:
        cmd = ["conda", "run", "-n", env_name] + cmd
    try:
        subprocess.run(cmd, check=True, timeout=180)
        return os.path.exists(faa_path) and os.path.getsize(faa_path) > 100
    except Exception as e:
        print(f"  [-] Prodigal 预测失败: {e}", file=sys.stderr)
        return False

def verify_with_diamond(query_faa, target_faa, env_name="bio", tmp_dir="/tmp"):
    if not query_faa or not os.path.exists(query_faa) or not target_faa or not os.path.exists(target_faa):
        return 0, []
    
    clean_target = os.path.splitext(os.path.basename(target_faa))[0]
    db_prefix = os.path.join(tmp_dir, f"dmnd_{clean_target}")
    out_tsv = os.path.join(tmp_dir, f"dmnd_out_{clean_target}.tsv")
    
    cmd_mkdb = ["diamond", "makedb", "--in", target_faa, "-d", db_prefix, "--quiet"]
    cmd_blastp = [
        "diamond", "blastp",
        "-q", query_faa,
        "-d", db_prefix,
        "-o", out_tsv,
        "--id", "25",
        "--query-cover", "30",
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue",
        "--quiet"
    ]
    if env_name:
        cmd_mkdb = ["conda", "run", "-n", env_name] + cmd_mkdb
        cmd_blastp = ["conda", "run", "-n", env_name] + cmd_blastp
        
    matched_queries = set()
    try:
        subprocess.run(cmd_mkdb, check=True, timeout=120)
        subprocess.run(cmd_blastp, check=True, timeout=120)
        if os.path.exists(out_tsv):
            with open(out_tsv, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if parts:
                        matched_queries.add(parts[0])
            os.remove(out_tsv)
        if os.path.exists(f"{db_prefix}.dmnd"):
            os.remove(f"{db_prefix}.dmnd")
    except Exception as e:
        print(f"  [-] Diamond 探针比对提示: {e}", file=sys.stderr)
        
    return len(matched_queries), sorted(list(matched_queries))

def main():
    args = parse_args()
    os.makedirs(args.out_genomes_dir, exist_ok=True)
    os.makedirs(args.out_proteins_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_manifest)), exist_ok=True)
    
    print("========================================================================")
    print(">>> GTDB 细菌参考基因组按需分层下载与验证引擎")
    print(f">>> 统一参考核酸存放目录: {args.out_genomes_dir}")
    print(f">>> 统一参考蛋白存放目录: {args.out_proteins_dir}")
    print(f">>> 设定数量目标: 阳性参考 {args.max_positive_refs} 株 | 阴性姐妹参考 {args.max_negative_refs} 株")
    print("========================================================================")
    
    # 1. 解析所有待测 MAG 的分类学，全局去冗余
    target_taxonomies = []
    target_accessions = set()
    
    if args.target_classification and os.path.exists(args.target_classification):
        with open(args.target_classification, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0].lower() not in ["mag_id", "genome_id", "accession"]:
                    target_taxonomies.append(parts[1].strip())
                    clean_t = extract_clean_acc(parts[0])
                    target_accessions.add(clean_t)
                    
    if not target_taxonomies and args.target_taxonomy:
        target_taxonomies.append(args.target_taxonomy.strip())
        
    if not target_taxonomies:
        target_taxonomies.append("d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Methylococcales;f__Methylomonadaceae;g__Methyloprofundus;s__")

    unique_genera, unique_families, unique_orders, unique_classes = set(), set(), set(), set()
    for t_str in target_taxonomies:
        for p in t_str.split(";"):
            p = p.strip()
            if p.startswith("g__") and len(p) > 3: unique_genera.add(p)
            elif p.startswith("f__") and len(p) > 3: unique_families.add(p)
            elif p.startswith("o__") and len(p) > 3: unique_orders.add(p)
            elif p.startswith("c__") and len(p) > 3: unique_classes.add(p)

    print(f"[+] 待测样本群亲缘阶元覆盖: {len(target_taxonomies)} 个目标样本")
    print(f"    - 涉及属 (Genera): {', '.join(sorted(unique_genera)) or 'None'}")
    print(f"    - 涉及科 (Families): {', '.join(sorted(unique_families)) or 'None'}")
    print(f"    - 涉及目 (Orders): {', '.join(sorted(unique_orders)) or 'None'}")
    print(f"    - 涉及纲 (Classes): {', '.join(sorted(unique_classes)) or 'None'}")

    # 2. 检查本地已存在的参考基因组（避免任何重复下载）
    confirmed_positive = []
    confirmed_negative = []
    already_evaluated_accs = set()
    
    existing_fnas = [f for f in os.listdir(args.out_genomes_dir) if any(f.endswith(ext) for ext in [".fna", ".fa", ".fasta"])]
    if existing_fnas:
        print(f"\n[+] 扫描本地统一参考目录，发现 {len(existing_fnas)} 个已有基因组，优先就地筛查（0 下载开销）...")
        for ef in sorted(existing_fnas):
            acc = extract_clean_acc(ef)
            if acc in already_evaluated_accs: continue
            out_fna = os.path.join(args.out_genomes_dir, ef)
            out_faa = os.path.join(args.out_proteins_dir, f"{acc}.faa")
            
            predict_proteins_if_needed(out_fna, out_faa, env_name=args.conda_env)
            hits, genes = verify_with_diamond(args.query_operon_faa, out_faa, env_name=args.conda_env)
            already_evaluated_accs.add(acc)
            
            rec = {
                "accession": acc,
                "organism_name": f"Local cache ({acc})",
                "taxonomy": "Cached reference",
                "sox_genes_detected": hits,
                "matched_genes": ",".join(genes),
                "genome_fna": out_fna,
                "protein_faa": out_faa
            }
            if hits >= args.positive_gene_threshold and len(confirmed_positive) < args.max_positive_refs:
                rec["verdict"] = "CONFIRMED_POSITIVE"
                confirmed_positive.append(rec)
                print(f"  [+] [本地缓存阳性] {acc}: 检出 {hits} 个 Sox 基因 ({', '.join(genes)}) [阳性: {len(confirmed_positive)}/{args.max_positive_refs}]")
            elif hits == 0 and len(confirmed_negative) < args.max_negative_refs:
                rec["verdict"] = "CONFIRMED_NEGATIVE"
                confirmed_negative.append(rec)
                print(f"  [-] [本地缓存阴性] {acc}: 检出 0 个 Sox 基因 [阴性: {len(confirmed_negative)}/{args.max_negative_refs}]")
                
        if len(confirmed_positive) >= args.max_positive_refs and len(confirmed_negative) >= args.max_negative_refs:
            print(f"\n[+] 本地已有参考基因组已完全满足设定要求（阳性 {len(confirmed_positive)}/{args.max_positive_refs}, 阴性 {len(confirmed_negative)}/{args.max_negative_refs}），无需下载！")
            pd.DataFrame(confirmed_positive + confirmed_negative).to_csv(args.out_manifest, sep="\t", index=False)
            return

    # 3. 读取 GTDB bac120 元数据表构建候选池
    meta_df = pd.DataFrame()
    if os.path.exists(args.gtdb_metadata):
        print(f"\n[+] 正在载入 GTDB 细菌元数据表: {args.gtdb_metadata} ...", flush=True)
        try:
            with open(args.gtdb_metadata, "r", encoding="utf-8", errors="ignore") as _f:
                header_cols = [c.strip() for c in _f.readline().split("\t")]
            candidate_cols = {
                "accession", "gtdb_taxonomy", "taxonomy", 
                "checkm_completeness", "completeness", 
                "checkm_contamination", "contamination",
                "gtdb_representative", "representative",
                "ncbi_organism_name", "ncbi_assembly_level"
            }
            use_cols = [c for c in header_cols if c in candidate_cols]
            if not use_cols or len(use_cols) < 2:
                use_cols = None
            meta_df = pd.read_csv(args.gtdb_metadata, sep="\t", usecols=use_cols, low_memory=False)
            print(f"[+] 成功读取 {len(meta_df)} 条 GTDB 细菌记录。", flush=True)
        except Exception as e:
            print(f"[-] 读取元数据表失败: {e}", file=sys.stderr)
            
    candidate_queue = []
    
    # 优先压入已知的 Seed 阳性菌（如 Reinekea forsetii，探针直属来源）
    for s_acc in args.seed_positive_accs:
        clean_s = extract_clean_acc(s_acc)
        if clean_s not in already_evaluated_accs and clean_s not in target_accessions:
            candidate_queue.append({
                "accession": clean_s,
                "tier": "Tier0_Seed_Positive",
                "organism_name": "Reinekea forsetii Hel1_31_1 (Seed)",
                "taxonomy": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Oceanospirillales;f__Saccharospirillaceae;g__Reinekea;s__Reinekea forsetii"
            })

    if not meta_df.empty:
        acc_col = "accession" if "accession" in meta_df.columns else meta_df.columns[0]
        tax_col = "gtdb_taxonomy" if "gtdb_taxonomy" in meta_df.columns else "taxonomy"
        comp_col = "checkm_completeness" if "checkm_completeness" in meta_df.columns else ""
        contam_col = "checkm_contamination" if "checkm_contamination" in meta_df.columns else ""
        rep_col = "gtdb_representative" if "gtdb_representative" in meta_df.columns else ""
        org_col = "ncbi_organism_name" if "ncbi_organism_name" in meta_df.columns else acc_col
        lvl_col = "ncbi_assembly_level" if "ncbi_assembly_level" in meta_df.columns else ""

        # 仅筛选 GTDB 代表基因组 + 高质控标准（代表基因组天然杜绝同一物种/菌株的重复冗余）
        sub = meta_df.copy()
        if rep_col and rep_col in sub.columns:
            sub = sub[sub[rep_col].astype(str).str.lower().isin(["t", "true", "1"])]
        if comp_col and contam_col:
            sub = sub[(sub[comp_col] >= 90.0) & (sub[contam_col] <= 5.0)]
            sub["q_score"] = sub[comp_col].astype(float) - 5.0 * sub[contam_col].astype(float)
        else:
            sub["q_score"] = 100.0

        if lvl_col and lvl_col in sub.columns:
            level_weights = {"complete genome": 100, "chromosome": 80, "scaffold": 50, "contig": 10}
            sub["lvl_score"] = sub[lvl_col].astype(str).str.lower().map(lambda x: level_weights.get(x, 20))
            sub["rank_score"] = sub["q_score"] + sub["lvl_score"]
        else:
            sub["rank_score"] = sub["q_score"]

        sub = sub.sort_values(by="rank_score", ascending=False)

        # 分层挑选候选菌：
        # Tier 1: 同属其他代表物种
        for g in unique_genera:
            g_hits = sub[sub[tax_col].str.contains(g, na=False, regex=False)]
            for _, r in g_hits.head(3).iterrows():
                c_acc = extract_clean_acc(r[acc_col])
                if c_acc not in already_evaluated_accs and c_acc not in target_accessions and c_acc not in [x["accession"] for x in candidate_queue]:
                    candidate_queue.append({"accession": c_acc, "tier": f"Tier1_{g}", "organism_name": str(r.get(org_col, c_acc)), "taxonomy": str(r[tax_col])})

        # Tier 2: 同科其他代表属
        for f in unique_families:
            f_hits = sub[sub[tax_col].str.contains(f, na=False, regex=False)]
            for _, r in f_hits.head(6).iterrows():
                c_acc = extract_clean_acc(r[acc_col])
                if c_acc not in already_evaluated_accs and c_acc not in target_accessions and c_acc not in [x["accession"] for x in candidate_queue]:
                    candidate_queue.append({"accession": c_acc, "tier": f"Tier2_{f}", "organism_name": str(r.get(org_col, c_acc)), "taxonomy": str(r[tax_col])})

        # Tier 3: 同目其他代表科
        for o in unique_orders:
            o_hits = sub[sub[tax_col].str.contains(o, na=False, regex=False)]
            for _, r in o_hits.head(6).iterrows():
                c_acc = extract_clean_acc(r[acc_col])
                if c_acc not in already_evaluated_accs and c_acc not in target_accessions and c_acc not in [x["accession"] for x in candidate_queue]:
                    candidate_queue.append({"accession": c_acc, "tier": f"Tier3_{o}", "organism_name": str(r.get(org_col, c_acc)), "taxonomy": str(r[tax_col])})

        # Tier 4: 同纲（如 Gammaproteobacteria）其他代表目（用于在近缘科全缺失时，扩展捕获阳性参考）
        for c in unique_classes:
            c_hits = sub[sub[tax_col].str.contains(c, na=False, regex=False)]
            # 优先包含已知硫代硫酸盐氧化近缘目代表
            sox_clades = ["Halothiobacillaceae", "Hydrogenovibrionaceae", "Thioalkalivibrionaceae", "Chromatiaceae", "Oceanospirillales"]
            for s_clade in sox_clades:
                sc_hits = c_hits[c_hits[tax_col].str.contains(s_clade, na=False, regex=False)]
                for _, r in sc_hits.head(2).iterrows():
                    c_acc = extract_clean_acc(r[acc_col])
                    if c_acc not in already_evaluated_accs and c_acc not in target_accessions and c_acc not in [x["accession"] for x in candidate_queue]:
                        candidate_queue.append({"accession": c_acc, "tier": f"Tier4_{s_clade}", "organism_name": str(r.get(org_col, c_acc)), "taxonomy": str(r[tax_col])})

    print(f"\n[+] 构建分层下载待选队列: 共 {len(candidate_queue)} 株高质量代表菌")
    print(f"    - 严格遵循：按亲缘由近至远逐个尝试，一旦达到配额（阳性={args.max_positive_refs}，阴性={args.max_negative_refs}）立即停止！")

    # 4. 逐株下载并现场比对判定
    manifest_records = confirmed_positive + confirmed_negative
    
    for cand in candidate_queue:
        if len(confirmed_positive) >= args.max_positive_refs and len(confirmed_negative) >= args.max_negative_refs:
            print(f"\n========================================================================")
            print(f"[+] 完美达到设置要求（阳性参考: {len(confirmed_positive)}/{args.max_positive_refs}，阴性参考: {len(confirmed_negative)}/{args.max_negative_refs}）！")
            print(f"[+] 触发早停机制，停止后续候选菌下载。")
            print(f"========================================================================")
            break
            
        acc = cand["accession"]
        tier = cand["tier"]
        org_name = cand["organism_name"]
        out_fna = os.path.join(args.out_genomes_dir, f"{acc}.fna")
        out_faa = os.path.join(args.out_proteins_dir, f"{acc}.faa")
        
        print(f"\n>>> 考察候选参考 [{tier}]: {acc} ({org_name})")
        
        # 1. 尝试下载
        dl_ok = os.path.exists(out_fna) and os.path.getsize(out_fna) > 0
        if not dl_ok:
            dl_ok = download_genome_via_datasets(acc, out_fna, out_faa, env_name=args.conda_env)
            
        if not dl_ok:
            print(f"  [-] 下载失败或暂不可用，跳过该菌株。")
            continue
            
        # 2. 预测蛋白质
        predict_proteins_if_needed(out_fna, out_faa, env_name=args.conda_env)
        
        # 3. 现场 Diamond 比对判断是否具有目标 operon
        hits, genes = verify_with_diamond(args.query_operon_faa, out_faa, env_name=args.conda_env)
        
        rec = {
            "accession": acc,
            "organism_name": org_name,
            "taxonomy": cand["taxonomy"],
            "tier": tier,
            "sox_genes_detected": hits,
            "matched_genes": ",".join(genes),
            "genome_fna": out_fna,
            "protein_faa": out_faa
        }
        
        if hits >= args.positive_gene_threshold:
            if len(confirmed_positive) < args.max_positive_refs:
                rec["verdict"] = "CONFIRMED_POSITIVE"
                confirmed_positive.append(rec)
                manifest_records.append(rec)
                print(f"  [+] [锁定阳性参考] 检出 {hits} 个 Sox 基因 ({', '.join(genes)}) -> 阳性参考进度: {len(confirmed_positive)}/{args.max_positive_refs}")
            else:
                print(f"  [~] 该菌具有 Sox 基因 ({hits}个)，但阳性参考配额已满 ({len(confirmed_positive)}/{args.max_positive_refs})，保留文件并继续。")
        elif hits == 0:
            if len(confirmed_negative) < args.max_negative_refs:
                rec["verdict"] = "CONFIRMED_NEGATIVE"
                confirmed_negative.append(rec)
                manifest_records.append(rec)
                print(f"  [-] [锁定阴性参考] 检出 0 个 Sox 基因 -> 阴性参考进度: {len(confirmed_negative)}/{args.max_negative_refs}")
            else:
                print(f"  [~] 该菌无 Sox 基因，但阴性参考配额已满 ({len(confirmed_negative)}/{args.max_negative_refs})，保留文件并继续。")
        else:
            rec["verdict"] = "PARTIAL_MATCH"
            manifest_records.append(rec)
            print(f"  [!] 提示: 该菌仅有部分亚基残余 ({hits}个)，不作为严格参考。")

    # 5. 输出汇总清单
    manifest_df = pd.DataFrame(manifest_records)
    manifest_df.to_csv(args.out_manifest, sep="\t", index=False)
    
    print("\n========================================================================")
    print(f"[+] 参考基因组按需下载与验证汇总清单 -> {args.out_manifest}")
    print("========================================================================")
    if not manifest_df.empty:
        cols_to_show = [c for c in ["accession", "organism_name", "verdict", "sox_genes_detected", "tier"] if c in manifest_df.columns]
        print(manifest_df[cols_to_show].to_string(index=False))

if __name__ == "__main__":
    main()
