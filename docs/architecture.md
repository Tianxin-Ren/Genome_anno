# 架构设计与重构说明

## 1. 现有代码库的核心设计思想

原仓库已经形成了四个有价值的基础设计：

1. 顶层 `Snakefile` 是工作流注册器，通过 `modules.*.enable` 动态 include 功能模块。
2. 每个模块拥有独立配置和命名空间，例如 `annotation`、`binning`、`abundance`。
3. `projects/config_*` 代表不同任务实例，同一套规则可服务不同样本和分析组合。
4. `profiles/slurm`、`profiles/slurm_pi` 把部分集群资源参数从规则中抽离。

此前的主要耦合点是：各模块重复实现路径/环境解析；项目配置含有具体集群的绝对 Conda 和数据库路径；部分模块仍写死 `results/...`；日志散落在结果目录；根目录同时存在正式脚本、调试脚本、历史副本和归档代码。

## 2. 目标架构：project × platform × profile

三个配置层的职责必须保持正交：

| 层 | 描述 | 可复用范围 |
|---|---|---|
| project | 样本、模块开关、科学参数、项目输出根目录 | 一个分析任务 |
| platform | 已部署 Conda 环境和本地数据库的逻辑别名 | 一台工作站或一个集群 |
| profile | executor、分区、CPU、内存、并发数 | 一类调度后端/队列 |

同一项目从本地切到集群时，不修改样本和分析参数：

```text
projects/my_task + config/platforms/local.yaml + profiles/local
projects/my_task + config/platforms/cluster-a.yaml + profiles/slurm
projects/my_task + config/platforms/cluster-b.yaml + profiles/slurm_pi
```

配置优先级为：platform 基础值 < project 分析配置 < 命令行显式参数。环境与数据库属于基础设施，platform 中的逻辑别名会在最后映射到模块旧键，从而覆盖旧项目中写死的机器路径。

## 3. Conda 与数据库部署

工作流目前采用已经部署好的命名 Conda 环境，通过 `conda run -p <env-path>` 调用。platform 中只需声明一个根目录和环境目录名：

```yaml
paths:
  conda_envs: /opt/conda/envs

environments:
  bakta: bakta_1.12.0
  bio: bio
  gtdbtk: gtdbtk-2.4.0
```

环境值也可以直接是绝对路径。`common/scripts/runtime.py` 会把全局别名注入每个模块的 `conda_envs`，模块仍可保留自己的默认环境名。

数据库同样使用逻辑别名：

```yaml
databases:
  bakta: /db/bakta/db-light
  gtdbtk: /db/gtdbtk/release232
  checkm2: /db/checkm2/uniref100.KO.1.dmnd
```

目前兼容映射覆盖 annotation 的 Bakta、GTDB-Tk、CheckM2、EggNOG、KOfam、dbCAN、Pfam/HMM，binning 的 CheckM，abundance 的 SortMeRNA，以及 plasmid_analysis 的 geNomad。

本地模板支持环境变量：

```bash
export GENOME_ANNO_CONDA_ENVS=/opt/conda/envs
export GENOME_ANNO_DATABASES=/data/databases
```

在真正运行前可检查变量和物理路径：

```bash
python scripts/validate_config.py \
  --project projects/my_task \
  --platform config/platforms/local.yaml \
  --check-paths
```

## 4. 多集群运行

一个集群通常需要两个文件：

- `config/platforms/<cluster>.yaml`：共享存储上的环境和数据库；
- `profiles/<cluster>/config.yaml`：executor、partition、内存、线程和并发数。

仓库已有两组可直接配对的集群配置：

| 调度 profile | 部署 platform | 主要存储 |
|---|---|---|
| `profiles/slurm` | `config/platforms/slurm.yaml` | DSSG 环境与数据库 |
| `profiles/slurm_pi` | `config/platforms/slurm_pi.yaml` | Lustre 项目环境 + DSSG 共享数据库 |

新集群可复制其中之一，并从 `config/platforms/slurm-example.yaml` 建立对应 platform。科学参数不应放进 profile，分区和调度参数不应放进 project。

## 5. 项目与输出组织

推荐每个任务使用独立的 `paths.project_root`：

```text
<project_root>/
├── data/                    # 原始或外部输入
├── results/
│   ├── preprocess/
│   ├── assembly/
│   ├── binning/
│   ├── annotation/
│   └── ...
├── logs/
│   ├── preprocess/
│   ├── assembly/
│   ├── annotation/
│   └── ...
└── work/                    # 临时/中间工作区的统一预留位置
```

`paths.results`、`paths.logs`、`paths.work` 均可自定义。兼容层会把旧的 `results/<module>` 自动重定位到新的 results 根目录；模块间引用应从其他模块的 `results_dir` 读取，不再自行拼写 `results/...`。

## 6. 多任务复用与扩展约束

- 新任务：复制 `projects/template`，只改项目配置，不复制规则和通用脚本。
- 新机器：新增 platform，不复制项目配置。
- 新队列/集群：新增 profile，不修改规则。
- 新功能：新增模块目录，并在顶层注册；模块只暴露最终 targets。
- 模块路径：使用本模块的 `results_dir`/`logs_dir`，跨模块输入读取对方配置。
- 外部资源：使用 platform 逻辑别名，不提交用户或集群绝对路径。
- 活动配置只声明已经接入规则的参数；未实现功能必须在解析时明确报错，不能静默忽略。
- 可选证据通道通过条件输入和空证据占位文件组合，关闭通道不会被下游整合重新触发。
- 调试与一次性分析：保留在 `scratch/`；稳定后迁移到对应模块的 `scripts/`。
- 历史实现：只放 `archive/`，不被顶层工作流 include。

## 7. 后续精简建议

本次改动保留旧项目兼容，没有删除任何用户文件。下一阶段可安全分批处理：

1. 删除或归档 `*冲突副本*`、根目录调试脚本和已确认不用的旧配置（已完成第一批清理，移入 `archive/`）。
2. 将所有规则从手写 `conda run -p` 逐步迁移到统一 helper，或标准 Snakemake `conda:` 声明。
3. 为每个模块增加最小 dry-run fixture 与配置 schema。
4. 将规则中的脚本路径改为基于 workflow root 的绝对路径，随后可让每个任务拥有独立 Snakemake working directory，支持完全并发运行。

## 8. 冗余模块收敛约定（避免重复造轮子）

仓库内部分功能模块因历史演进存在能力重叠，新增或修改功能时遵循以下收敛约定：

- **功能注释主干**：新任务优先使用 `modules/unified_annotation/`；`modules/annotation/` 仅保留兼容旧项目及 `circular_evidence` 等特有子流程；`modules/mag_annotation/`（高通量 PyKofamSearch）与 `modules/cog_annotation/` 仅作为专项目标模块维护，不再向旧版 `annotation` 重复堆叠新工具。
- **Hydrogenase 注释**：`tools/hydrogenase_annotation/` 作为底层分析内核与参考库保留，调度与路径解析统一由 `modules/hydrogenase_annotation/` 经 `common/scripts/runtime.py`（`MODULE_NAMES`）接管，不再扩展外部独立调度入口。
- **环境与配置入口**：Conda 环境定义统一收口至 `envs/`（重建参考）与 `config/platforms/*.yaml`（运行时挂载），`modules/backup_envs/` 与重复的 `config_mag_annotation/` 已归档至 `archive/`。
- **单模块脚本归属**：新增仅服务于单一模块的脚本应直接放入 `modules/<module>/scripts/`，避免继续向根目录 `scripts/` 堆积一次性脚本。
