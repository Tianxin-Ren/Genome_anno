# Unified bacterial/archaeal annotation task

`config_En` 现在只描述任务和科学参数；Conda 环境及数据库路径由 platform 提供。

## 文件职责

- `pipeline.yaml`：样本文件、模块开关、默认 platform 和输出根目录；
- `samples.yaml`：基因组 FASTA 与样本属性；
- `unified_annotation.yaml`：已经实现的步骤开关及科学参数；
- `check_software_working.sh`：只读软件检查；
- `check_existing_databases.sh`：只读数据库检查。

默认部署是 `config/platforms/slurm.yaml`，因此原命令保持有效：

```bash
snakemake --workflow-profile profiles/slurm --config conf=config_En
```

切换到 PI 集群：

```bash
snakemake \
  --workflow-profile profiles/slurm_pi \
  --config conf=config_En platform=config/platforms/slurm_pi.yaml
```

如需同时覆盖输出根目录，可继续使用原有的扁平命令行参数：

```bash
snakemake \
  --workflow-profile profiles/slurm_pi \
  --config \
  conf=config_En \
  platform=config/platforms/slurm_pi.yaml \
  results_dir=/lustre/path/to/project/results
```

配置解析检查：

```bash
python scripts/validate_config.py --project config_En
```

## Step 语义

以下步骤开关均已接入 DAG：

- `genome_qc`、`taxonomy`、`cds_prediction`、`rna_prediction`；
- `protein_similarity`、`orthology`、`ko_profiles`、`protein_domains`；
- `localization_topology`、`signal_peptide_prediction`、`biob_hmm`；
- `candidate_integration`。

关闭某个蛋白证据通道时，整合流程读取明确的空证据文件，不会重新触发该通道。所有蛋白证据与整合依赖 `cds_prediction: true`；违反依赖关系会在构建 DAG 时直接报错。

CRISPR、deep disambiguation、Bakta/Prokka reference check 尚未实现，因此不再出现在活动配置中。如果旧配置将这些步骤设为 `true`，工作流会明确报错，而不是静默忽略。

## 输出

默认结果根目录来自 `pipeline.yaml`：

```text
/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/results/unified_annotation/
```

规则日志和临时文件分别位于：

```text
<project_root>/logs/unified_annotation/
<project_root>/work/                     # 统一预留工作目录
```
