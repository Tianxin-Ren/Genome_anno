# Metagenomic Read Profiling (config_read)

本配置工程专用于宏基因组测序短读长（Short-read reads）的分类丰度（Taxonomic Profiling）分析，使用超快速宏基因组丰度估算工具 **Sylph**。

## 目录结构与配置文件
- [pipeline.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_read/pipeline.yaml): 流水线主控配置文件，负责平台挂载（`gses0`）与模块开关（纯净架构，不包含任何 `project_root` 数据路径）。
- [samples.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_read/samples.yaml): 样本定义核心文件，包含全局 `project_root` 以及显式列出的各个 `{sample}` 名称与洁净 reads 物理路径。
- [read_profiling.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_read/read_profiling.yaml): Sylph 模块配置，定义相对输出路径 `results_dir`、`logs_dir`、数据库路径及单样本/批次线程数。
- [run_cluster.sh](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_read/scripts/run_cluster.sh): 预置 SLURM 作业调度提交脚本（保留原始 `#SBATCH` 资源参数）。

---

## 样本定义与 TSV 自动导入

当前 `samples.yaml` 的 `data.samples` 指向包含 `RUN_id` 和 `Store_address` 的 TSV。每行按 `Store_address/00.QC/{RUN_id}.clean.{1,2}.fq.gz` 查找 R1/R2；`read_profiling.yaml` 中的 `skip_missing_inputs: true` 使缺失任一文件的样本不进入调度和汇总。

每次运行（包括 `-n` dry-run）都会更新 `{project_root}/logs/read_profiling/input_status.tsv`，记录每个样本的原始目录、解析后的目录、R1/R2 路径、文件是否存在、`00.QC` 符号链接目标及 `status`。`status=ready` 表示进入调度；`missing_r1`、`missing_r2`、`missing_both` 表示被跳过。可用下面的命令查看跳过的样本：

```bash
awk -F '\t' 'NR == 1 || $NF != "ready"' /dssg/home/acct-schhui/gses0/rentianxin/logs/read_profiling/input_status.tsv
```

### fastp 质控损坏文件容错与报告机制
对于需要通过 fastp 进行质控的原始测序文件（`status=ready_needs_qc`）：
- 若某个样本的原始 `.fq.gz` 文件损坏（如意外截断、CRC32 校验失败、解压异常或 0 字节空文件），**流水线不会中断报错退出**。
- 流程将自动捕获异常并执行以下操作：
  1. 记录损坏样本信息到全局报告表：`{project_root}/logs/read_profiling/corrupted_samples.tsv`
  2. 在该样本的质控输出目录下生成标记文件：`{store}/00.QC/{sample}.corrupted` 和 `CORRUPTED.txt`
  3. 创建洁净占位数据，下游 Sylph 物种丰度分析与合并步骤会自动识别并安全跳过该样本
  4. 保证其余健康样本不受任何影响，顺利完成分析与最终多样本矩阵聚合（`merged_sylph_profiling.tsv` 与 `merged_*_abundance.tsv`）

查看所有损坏样本清单：
```bash
cat /dssg/home/acct-schhui/gses0/rentianxin/logs/read_profiling/corrupted_samples.tsv
```

也可以在 [samples.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_read/samples.yaml) 中显式提供单个样本的 R1/R2 路径：

```yaml
project_root: "/dssg/home/acct-schhui/gses0/liliuyang/osma/test"

samples:
  "DP8480005241BL_L01_102":
    short_read:
      r1: "/union/home/acct-schhui/gses0/output_all/DP8480005241BL_L01_102/00.QC/DP8480005241BL_L01_102.clean.1.fq.gz"
      r2: "/union/home/acct-schhui/gses0/output_all/DP8480005241BL_L01_102/00.QC/DP8480005241BL_L01_102.clean.2.fq.gz"
  "DP8480005241BL_L01_60":
    short_read:
      r1: "/union/home/acct-schhui/gses0/output_all/DP8480005241BL_L01_60/00.QC/DP8480005241BL_L01_60.clean.1.fq.gz"
      r2: "/union/home/acct-schhui/gses0/output_all/DP8480005241BL_L01_60/00.QC/DP8480005241BL_L01_60.clean.2.fq.gz"
```

如果您手头有包含 `RUN_id` 和 `Store_address` 的 TSV 文件（例如 `samples.tsv`），可直接使用配套工具一键更新生成：
```bash
python3 scripts/tsv_to_samples_yaml.py \
  --tsv projects/config_read/samples.tsv \
  --project-root "/dssg/home/acct-schhui/gses0/liliuyang/osma/test" \
  --output projects/config_read/samples.yaml
```

---

## 运行方式

### 归档到 `/union`（在登录节点手动执行）

分析完成、确认没有正在写入这些文件的作业后，在登录节点重新生成 `/dssg` → `/union` 归档清单：

```bash
python3 projects/config_read/scripts/archive_plan.py \
  --tsv /dssg/home/acct-schhui/gses0/rentianxin/data/osma_progress/sylph_batch2.tsv \
  --results-dir /dssg/home/acct-schhui/gses0/rentianxin/results/MetaG-reads/sylph \
  --output /dssg/home/acct-schhui/gses0/rentianxin/logs/read_profiling/archive_plan.tsv
```

清单逐文件记录原路径、符号链接解析后的真实路径、对应的 `/union` 目标、大小及候选动作。`move_to_union_and_link` 表示实体文件位于 `/dssg`；`already_on_union` 表示文件已在 `/union`；`delete_intermediate` 表示已确认可删除的 fastp 中间文件。盘点不会搬动或删除文件，也不会被默认的 `rule all` 调用。

检查清单后，在**可写 `/union` 的登录节点**执行：

```bash
python3 projects/config_read/scripts/archive_plan.py \
  --apply-plan /dssg/home/acct-schhui/gses0/rentianxin/logs/read_profiling/archive_plan.tsv
```

执行时先核对清单里的文件大小与修改时间，逐个复制到 `/union` 并校验 SHA-256，最后把 `/dssg` 原实体文件替换为指向 `/union` 的软链接。所有保留文件迁移成功后，才删除清单中的 fastp merged/unmerged FASTQ 与 JSON/HTML；若某样本缺少原始成对 reads、两份 clean reads 或 Sylph 结果，该样本的中间文件继续保留。日志、审计表、`done` 和损坏标记保留在 `/dssg`。结果写入 `logs/read_profiling/archive_apply_report.tsv`。发现目标同名不同内容、源文件变化或路径异常时会停止，不会覆盖目标文件；未识别的文件不会自动处理。

### 方式 1: 直接通过 SLURM 提交集群调度运行（推荐）
```bash
sbatch projects/config_read/scripts/run_cluster.sh
```

### 方式 2: 使用 Snakemake 命令行运行
```bash
snakemake \
  --workflow-profile profiles/slurm \
  --config conf=projects/config_read \
  --rerun-incomplete \
  --printshellcmds
```

### 方式 3: 本地/交互节点测试
```bash
snakemake \
  --config conf=projects/config_read \
  -c 16 \
  --printshellcmds
```

---

## 结果与日志层级规范（参照既有设计）

产物严格按照仓库标准隔离在 `results_dir` 和 `logs_dir` 子目录下：

- **单样本丰度结果**：
  `{project_root}/result/MetaG-reads/sylph/{sample}/{sample}.sylph.tsv`
- **单样本运行日志**：
  `{project_root}/logs/read_profiling/sylph/{sample}.log`
- **全局汇总丰度表**：
  `{project_root}/result/MetaG-reads/sylph/merged_sylph_profiling.tsv`
- **汇总步骤日志**：
  `{project_root}/logs/read_profiling/merge_sylph.log`
