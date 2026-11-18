#!/usr/bin/env python3
import argparse
import sys
import subprocess
from pathlib import Path
import pandas as pd

# 尝试导入本地配置
try:
    import config
except ImportError:
    # 如果没有 config.py，则使用内部默认值
    class config:
        MIN_SCORE = 100
        MIN_BSR = 0.9
        LOOKUP_PL = "./blast_based_read_lookup_new.pl"

def main():
    parser = argparse.ArgumentParser(description="工业级宏基因组功能基因定量比对脚本 (BSR 计算)")

    # 核心输入输出参数
    parser.add_argument("-q", "--query", required=True, type=Path, help="输入 Fastq 文件 (支持 .gz)")
    parser.add_argument("-f", "--fasta_db", required=True, type=Path, help="目标功能基因库 FASTA")
    parser.add_argument("-d", "--outgroup_db", required=True, type=Path, help="外群库 DIAMOND DB (.dmnd)")
    parser.add_argument("-o", "--output", required=True, type=Path, help="输出 BSR 结果路径")

    # 高级参数
    parser.add_argument("-t", "--threads", default=64, type=int, help="线程数 (默认 64)")
    parser.add_argument("-k", "--max_target_seqs", default=500, type=int, help="DIAMOND 最大命中数")
    parser.add_argument("--keep_temp", action="store_true", help="保留中间文件")
    parser.add_argument("--use_perl", action="store_true", help="使用外部 Perl 脚本提取序列 (默认使用 Python 内置)")
    parser.add_argument("--dry_run", action="store_true", help="空运行模式")

    args = parser.parse_args()

    # 检查必需文件
    for file_path, name in [(args.query, "Query"), (args.fasta_db, "Fasta DB"), (args.outgroup_db, "Outgroup DB")]:
        if not file_path.is_file():
            sys.exit(f"[-] 错误: {name} 文件不存在: {file_path}")

    outdir = args.output.parent
    if not args.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.output.stem

    # 中间文件定义
    gene_dmnd_db = outdir / f"{prefix}_gene_temp"
    actual_db_file = gene_dmnd_db.with_suffix('.dmnd')
    query_fasta = outdir / f"{prefix}_query_temp.fasta"  
    unfiltered_1st = outdir / f"{prefix}_1st_unfiltered.tsv"
    filtered_1st = outdir / f"{prefix}_1st_filtered.tsv"
    retrieved_reads = outdir / f"{prefix}_1st_hits.faa"
    unfiltered_2nd = outdir / f"{prefix}_2nd_unfiltered.tsv"
    filtered_2nd = outdir / f"{prefix}_2nd_filtered.tsv"

    # DIAMOND 参数配置
    blast_col = ["qseqid", "sseqid", "pident", "qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score"]
    scoring = ["--matrix", "blosum45", "--gapopen", "14", "--gapextend", "2", "--masking", "0", "--min-score", "10"]
    performance = ["-b", "12", "-c", "1"]

    def run_cmd(cmd, shell=False):
        if args.dry_run:
            print(f"[DRY-RUN] {' '.join(map(str, cmd)) if not shell else cmd}")
        else:
            if shell:
                subprocess.run(cmd, shell=True, check=True)
            else:
                subprocess.run([str(x) for x in cmd], check=True)

    import time
    start_time = time.time()

    try:
        if args.output.exists():
            print(f"[*] 结果文件 {args.output.name} 已存在，跳过该任务。")
            return

        print(f"\n{'='*60}\n[*] 开始处理样本: {args.output.name}\n{'='*60}")

        # 1. 构建参考库
        t0 = time.time()
        if not actual_db_file.exists():
            print(f"[*] [Step 1/7] 正在构建 DIAMOND 库: {args.fasta_db.name}...")
            run_cmd(["diamond", "makedb", "--in", args.fasta_db, "-d", gene_dmnd_db, "-p", args.threads])
        print(f"[+] Step 1 完成, 耗时: {time.time()-t0:.2f}s")

        # 2. FASTQ 转 FASTA
        t0 = time.time()
        if not query_fasta.exists():
            print("[*] [Step 2/7] 正在转换序列格式 (FASTQ -> FASTA)...")
            cat_cmd = "zcat" if str(args.query).endswith(".gz") else "cat"
            awk_cmd = f"{cat_cmd} {args.query} | awk 'NR%4==1{{print \">\" substr($0,2)}} NR%4==2{{print}}' > {query_fasta}"
            run_cmd(awk_cmd, shell=True)
        print(f"[+] Step 2 完成, 耗时: {time.time()-t0:.2f}s")

        # 3. 第一轮靶向比对
        t0 = time.time()
        if not unfiltered_1st.exists():
            print("[*] [Step 3/7] 启动第一轮靶向比对 (k={})...".format(args.max_target_seqs))
            run_cmd(["diamond", "blastx", "-q", query_fasta, "-d", gene_dmnd_db, "-o", unfiltered_1st, 
                     "-p", args.threads, "-k", args.max_target_seqs, "--outfmt", "6", *blast_col, *scoring, *performance])
        print(f"[+] Step 3 完成, 耗时: {time.time()-t0:.2f}s")

        # 4. 去重过滤
        t0 = time.time()
        if not filtered_1st.exists():
            print("[*] [Step 4/7] 正在对命中结果进行去重处理 (Top 1 Hit)...")
            if not args.dry_run:
                chunk_iter = pd.read_csv(unfiltered_1st, sep="\t", names=blast_col, chunksize=5000000)
                chunks = [chunk.sort_values('score', ascending=False).drop_duplicates(['qseqid']) for chunk in chunk_iter]
                if chunks:
                    pd.concat(chunks).sort_values('score', ascending=False).drop_duplicates(['qseqid']).to_csv(filtered_1st, sep="\t", index=False)
                else:
                    print("[-] 警告: 第一轮比对无命中结果。")
                    return
        print(f"[+] Step 4 完成, 耗时: {time.time()-t0:.2f}s")

        # 5. 提取 Reads
        t0 = time.time()
        if not retrieved_reads.exists():
            if args.use_perl:
                print(f"[*] [Step 5/7] 正在使用 Perl 脚本提取序列...")
                run_cmd(["perl", config.LOOKUP_PL, filtered_1st, query_fasta, retrieved_reads])
            else:
                print("[*] [Step 5/7] 正在使用 Python 内置逻辑提取序列...")
                if not args.dry_run:
                    target_ids = set(pd.read_csv(filtered_1st, sep="\t", usecols=["qseqid"])["qseqid"])
                    with open(query_fasta, 'r') as f_in, open(retrieved_reads, 'w') as f_out:
                        write_flag = False
                        for line in f_in:
                            if line.startswith('>'):
                                seq_id = line[1:].split()[0]
                                write_flag = (seq_id in target_ids)
                            if write_flag: f_out.write(line)
        print(f"[+] Step 5 完成, 耗时: {time.time()-t0:.2f}s")

        # 6. 第二轮外群比对 (重大优化：强制使用 k=1 减少临时文件)
        t0 = time.time()
        if not unfiltered_2nd.exists():
            print("[*] [Step 6/7] 启动第二轮外群比对 (优化版: k=1)...")
            # 强制将 k 设为 1，因为我们只需要 Best Outgroup Hit
            run_cmd(["diamond", "blastx", "-q", retrieved_reads, "-d", args.outgroup_db, "-o", unfiltered_2nd, 
                     "-p", args.threads, "-k", 1, "--outfmt", "6", *blast_col, *scoring, *performance])
        print(f"[+] Step 6 完成, 耗时: {time.time()-t0:.2f}s")

        # 7. 计算 BSR
        t0 = time.time()
        if not args.output.exists():
            print("[*] [Step 7/7] 正在计算最终 BSR 矩阵...")
            if not args.dry_run:
                # 即使 k=1，也使用分块读取以防万一
                chunk_iter_og = pd.read_csv(unfiltered_2nd, sep="\t", names=blast_col, usecols=["qseqid", "sseqid", "pident", "score"], chunksize=5000000)
                chunks_og = [chunk.sort_values('score', ascending=False).drop_duplicates(['qseqid']) for chunk in chunk_iter_og]
                derep_og = pd.concat(chunks_og).sort_values('score', ascending=False).drop_duplicates(['qseqid'])
                
                df1 = pd.read_csv(filtered_1st, sep="\t", usecols=["qseqid", "sseqid", "pident", "score"])
                res = df1.rename(columns={'sseqid': 'sseqid_db', 'score': 'score_db', 'pident': 'pident_db'}).merge(
                      derep_og.rename(columns={'sseqid': 'sseqid_og', 'score': 'score_og', 'pident': 'pident_og'}), on="qseqid")
                res["bsr"] = res["score_db"] / res["score_og"]
                res.to_csv(args.output, sep="\t", index=False)
        print(f"[+] Step 7 完成, 耗时: {time.time()-t0:.2f}s")

        print(f"\n{'-'*60}\n[SUCCESS] 样本 {args.output.name} 处理完成! 总耗时: {time.time()-start_time:.2f}s\n{'-'*60}")

    finally:
        if not args.keep_temp and not args.dry_run and args.output.exists():
            print("[*] 正在清理临时文件...")
            for tmp in [actual_db_file, query_fasta, unfiltered_1st, filtered_1st, retrieved_reads, unfiltered_2nd]:
                if tmp.exists(): tmp.unlink()

if __name__ == "__main__":
    main()
