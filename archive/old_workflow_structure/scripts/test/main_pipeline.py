#!/usr/bin/env python3
import os
import sys
import subprocess
import logging
from multiprocessing import Pool
from pathlib import Path
from datetime import datetime

# 导入配置
try:
    import config
    import process_quant
except ImportError:
    print("[-] 错误: 找不到 config.py 或 process_quant.py。")
    sys.exit(1)

def setup_logging():
    """
    设置双通道日志系统：控制台 + 文件
    """
    logger = logging.getLogger("Pipeline")
    logger.setLevel(logging.INFO)
    
    # 格式化
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    
    # 控制台 Handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    # 文件 Handler
    fh = logging.FileHandler(config.LOG_FILE, mode='a', encoding='utf-8')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    return logger

def run_alignment_job(args):
    """
    并行子任务：比对
    注意：在 multiprocessing 中，子进程不会自动继承父进程的 logger 实例，
    但由于我们是调用外部命令，其 stdout 仍会流向父进程。
    """
    sample_id, fastq_path, db_key = args
    db_info = config.DB_CONFIG[db_key]
    output_path = os.path.join(config.OUT_DIR, f"{sample_id}{db_info['suffix']}")
    
    cmd = [
        "python3", "read_trawler_V1.py",
        "-q", fastq_path,
        "-f", db_info['faa'],
        "-d", db_info['dmnd'],
        "-o", output_path,
        "-t", str(config.THREADS_PER_JOB)
    ]
    
    # 手动构建日志消息
    msg = f"[*] [Launcher] 启动比对: {sample_id} -> {db_key}"
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} [INFO] {msg}")
    
    try:
        # 将子进程的输出合并到主进程的 stdout/stderr
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError:
        return False

def main():
    logger = setup_logging()
    
    logger.info("="*60)
    logger.info("🚀 宏基因组功能基因定量：全自动工业级流水线启动")
    logger.info("="*60)

    # 1. 确定样本列表
    if config.TARGET_SAMPLES:
        sample_ids = config.TARGET_SAMPLES
    else:
        fastq_files = list(Path(config.METAGENOMICS_DIR).glob("**/*_R1.clean.fastq.gz"))
        sample_ids = sorted(list(set([f.name.replace("_R1.clean.fastq.gz", "") for f in fastq_files])))

    if not sample_ids:
        logger.error("未发现待处理样本，请检查配置或目录结构。")
        return

    logger.info(f"发现待处理样本: {len(sample_ids)} 个")
    os.makedirs(config.OUT_DIR, exist_ok=True)

    # 2. 构建任务队列
    tasks = []
    for sid in sample_ids:
        fastq = os.path.join(config.METAGENOMICS_DIR, sid, "01_qc", f"{sid}_R1.clean.fastq.gz")
        if not os.path.exists(fastq):
            logger.warning(f"找不到样本 {sid} 的序列文件 (路径: {fastq})，已跳过。")
            continue
        
        for db_key in config.DB_CONFIG.keys():
            tasks.append((sid, fastq, db_key))

    # 3. 并行比对
    logger.info(f"正在并行对比 (并发任务数: {config.CONCURRENT_JOBS})...")
    with Pool(config.CONCURRENT_JOBS) as pool:
        results = pool.map(run_alignment_job, tasks)

    success_count = sum(results)
    logger.info(f"比对阶段结束。成功: {success_count}/{len(tasks)}")

    if success_count == 0:
        logger.error("所有比对任务均失败，停止后续统计。")
        return

    # 4. 后期整合
    logger.info("正在启动后期整合统计 (process_quant.py)...")
    # 为了让 process_quant 的输出也能进日志，我们稍微包装一下
    process_quant.main()

    logger.info("="*60)
    logger.info(f"✨ 流程执行完毕！最终矩阵见 {config.OUT_DIR}")
    logger.info("="*60)

if __name__ == "__main__":
    main()
