# Run preparation

## 1. 准备样本

编辑 `config_En/samples.yaml`，每个样本至少提供 FASTA：

```yaml
samples:
  sample_A:
    fasta: "/path/to/project/genomes/sample_A.fna"
    assembly_type: MAG
    expected_domain: auto
    genetic_code: auto
```

`assembly_type` 可使用 `isolate`、`MAG`、`SAG` 或 `metagenome_assembly`。已知古菌/细菌时可设置 `expected_domain: Archaea/Bacteria`；否则使用 `auto`。

## 2. 选择集群部署

DSSG：

```bash
python scripts/validate_config.py \
  --project config_En \
  --platform config/platforms/slurm.yaml
```

PI/Lustre：

```bash
python scripts/validate_config.py \
  --project config_En \
  --platform config/platforms/slurm_pi.yaml
```

登录对应集群后可增加 `--check-paths`，检查配置中的环境根目录和数据库文件是否存在。

## 3. 软件和数据库检查

DSSG 默认路径：

```bash
ENVS=/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs \
  bash config_En/check_software_working.sh

DBROOT=/dssg/home/acct-clswfp/clswfp-icdli1/Database \
  bash config_En/check_existing_databases.sh
```

PI 集群运行检查时，把 `ENVS` 替换为 `config/platforms/slurm_pi.yaml` 中的 `paths.conda_envs`。

## 4. 当前默认分析

当前 `unified_annotation.yaml` 开启：

- Pyrodigal CDS prediction；
- tRNAscan-SE、Infernal/Rfam、ARAGORN RNA prediction；
- MMseqs2/Swiss-Prot；
- EggNOG-mapper；
- KOfamScan；
- InterProScan/Pfam；
- 候选证据整合和跨基因组汇总。

QC、taxonomy、SignalP、DeepTMHMM 和 BioB HMM 默认关闭，可直接通过对应 step 开关启用。CRISPR、deep disambiguation、Bakta/Prokka review 尚未实现，配置为开启时会明确报错。

## 5. 运行

DSSG dry-run：

```bash
snakemake -n --workflow-profile profiles/slurm --config conf=config_En
```

DSSG 正式运行：

```bash
snakemake --workflow-profile profiles/slurm --config conf=config_En
```

PI/Lustre：

```bash
snakemake \
  --workflow-profile profiles/slurm_pi \
  --config conf=config_En platform=config/platforms/slurm_pi.yaml
```

## 6. 主要输出

```text
results/unified_annotation/<sample>/structural/
results/unified_annotation/<sample>/function/
results/unified_annotation/<sample>/tables/
results/unified_annotation/<sample>/report/
results/unified_annotation/summary/
```

实际 `results` 根目录由 `config_En/pipeline.yaml` 的 `paths.results` 决定。

