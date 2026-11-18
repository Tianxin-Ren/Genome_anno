# Genome annotation Snakemake workflow

这是一个面向本地工作站和多 SLURM 集群的模块化生物信息工作流。核心原则是把“分析任务”“机器部署”和“调度后端”分开：

- `projects/<task>/`：样本、启用模块和分析参数；
- `config/platforms/<site>.yaml`：该机器已有的 Conda 环境和数据库路径；
- `profiles/<backend>/`：Snakemake 本地或集群调度参数；
- `modules/<module>/`：可组合、可独立扩展的规则实现；
- `results/<module>/` 与 `logs/<module>/`：项目级产物和日志。

详细设计与迁移说明见 [docs/architecture.md](docs/architecture.md)。

## 快速开始

1. 创建任务配置：

   ```bash
   cp -R projects/template projects/my_task
   ```

2. 编辑 `projects/my_task/pipeline.yaml` 和 `samples.yaml`。项目数据及输出根目录由 `paths.project_root` 指定。

3. 本地部署使用环境变量，或直接编辑 `config/platforms/local.yaml`：

   ```bash
   export GENOME_ANNO_PROJECT_ROOT=/absolute/path/to/my_task
   export GENOME_ANNO_CONDA_ENVS=/absolute/path/to/conda/envs
   export GENOME_ANNO_DATABASES=/absolute/path/to/databases
   ```

4. 校验配置：

   ```bash
   python scripts/validate_config.py --project projects/my_task
   ```

5. 本地 dry-run 和执行：

   ```bash
   scripts/run_pipeline.sh --project projects/my_task --profile profiles/local -n
   scripts/run_pipeline.sh --project projects/my_task --profile profiles/local
   ```

6. 切换到 SLURM 集群只替换 platform 与 profile：

   ```bash
   scripts/run_pipeline.sh \
     --project projects/my_task \
     --platform config/platforms/my-cluster.yaml \
     --profile profiles/slurm \
     -n
   ```

   仓库内现有两套集群配对为：`profiles/slurm` + `config/platforms/slurm.yaml`，以及 `profiles/slurm_pi` + `config/platforms/slurm_pi.yaml`。

已有 `projects/config_*` 继续兼容；也仍可直接运行：

```bash
snakemake -s Snakefile --profile profiles/slurm --config conf=projects/config_an
```

## 新增模块

一个模块应包含 `modules/<name>/<name>.smk`、可选 `scripts/`，并只读取 `config["<name>"]`。在项目 `pipeline.yaml` 的 `modules` 中注册并启用它，模块把最终产物加入全局 `TARGETS`。机器路径不要写进模块或项目参数，应放到 platform 的 `environments`/`databases` 中。
