# 高通量古菌与 MAG 功能注释流水线 (config_mag_annotation)

本配置工程 `projects/config_mag_annotation/` 专为 **GTDB r232 全量古菌（22,343 基因组，35–55M 蛋白）** 及超大规模 MAG 集合量身定制，通过 **PyKofamSearch + 全量序列化 KOfam 库 (`database.pkl`)**，在 HPC 集群 **8 节点 × 64 核心 × 512 GB RAM** 条件下实现生产级、方法学发表标准的高通量注释。

---

## 🔬 核心数据链与设计原则

$$\text{Genome} \longrightarrow \text{Gene} \longrightarrow \text{Contig / Strand / Coordinate Order} \longrightarrow \text{Protein} \longrightarrow \text{KO}$$

1. **分片级原生 SLURM 调度 (Shard-Level Independent SLURM Scheduling)**：
   - 适配共享资源池模式（按实际申请 CPU/Mem 计费，多任务共享计算节点）。
   - 每个 Shard 作为一个独立的 SLURM 作业（默认 12 核心，~96 GB 内存）由 Snakemake 原生提交与调度。
   - **最高弹性与容错隔离**：最大并发度达 ~42 个作业并发（42 × 12 = 504 核心），任意单个 Shard 的偶发重试或失败不会影响其他 Shard；节点有几核空闲就跑几个作业，资源利用率最大化。
   - 支持动态自动暂存：若节点具备本地 `$SLURM_TMPDIR`，单 Shard 作业自动复用本地缓存的 `database.pkl`，规避网络存储并发 I/O 抖动。
2. **残基复杂度感知均衡切片 (Complexity-Aware Residue-Balanced Sharding)**：
   - 以 **~100–150 M amino acids (aa)** 为主要划分中枢（测试集为 ~15M aa），**500,000 proteins** 为上限保护，输出包含平均长度、中位数长度与最大长度统计的 `shards_manifest.tsv`，消除由于个别超长蛋白引发的拖尾效应。
3. **防碰撞精确去重 (Collision-Checked Exact Dedup)**：
   - 采用 SHA256 初筛 + 真实序列比对二级防御，杜绝极小概率哈希碰撞；在 `run_metadata.json` 中完整记录原始规模与搜索宇宙 $Z$。
4. **一级核心资产定义**：
   - `protein_KO_all_strict.tsv` 作为一级核心资产，保留所有满足 KOfam curated adaptive threshold（或 80% 相对阈值放宽）的候选命中（含 bitscore, score_type, 阈值差额 `score_margin`、命中位次 `ko_rank` 与严格判定 `is_strict`）。
   - 结合 `protein_metadata.tsv.gz` 记录的 Contig、起止坐标与链信息，为后续 ATP 合成酶、氢化酶等操纵子共线性与基因邻接分析奠定完整数据基座。


---

## 📂 配置文件索引

| 文件 | 路径与说明 |
| :--- | :--- |
| **`pipeline.yaml`** | [`pipeline.yaml`](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_mag_annotation/pipeline.yaml) - 主控配置，调度 `mag_annotation` 模块与 `config/platforms/slurm.yaml`。 |
| **`samples.yaml`** | [`samples.yaml`](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_mag_annotation/samples.yaml) - 输入声明，提供 `genome_dir`（DNA FASTA）或 `protein_dir`（预预测蛋白）。 |
| **`mag_annotation.yaml`** | [`mag_annotation.yaml`](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/projects/config_mag_annotation/mag_annotation.yaml) - 科学参数、残基切片阈值、Node 捆绑核心数与数据库配置。 |

---

## 💻 集群运行指令

### 1. 投递前预处理：数据库序列化
若尚未生成单文件序列化数据库 `database.pkl`，执行：
```bash
bash projects/config_mag_annotation/scripts/serialize_db.sh
```
该脚本将基于 `profiles/` 与 `ko_list` 生成 `database.pkl`，并伴随产出记录发布日期与所有 SHA256 校验码的 `kofam_manifest.json` 溯源表。

### 2. 投递前基准测试与验证门禁（推荐在单个计算节点上先跑）
- **6 点三维基准测试**（测定 20×3、10×6、8×8、5×12、4×16、2×32 在真实 CPU 上的峰值吞吐）：
  ```bash
  bash projects/config_mag_annotation/scripts/run_benchmark.sh
  ```
- **三级科学验证门禁**（检验分片大小不变性与古菌标志物召回）：
  ```bash
  bash projects/config_mag_annotation/scripts/run_validation.sh
  ```

### 3. 正式提交 SLURM 集群生产任务
```bash
sbatch projects/config_mag_annotation/scripts/run_cluster.sh
```

### 4. 超大规模输入的两阶段 Snakemake 运行方式

第一阶段按每批 100 个 genome 打包成 Snakemake job，只构建 gene-calling 任务：

```bash
snakemake --workflow-profile profiles/slurm \
  --config conf=projects/config_mag_annotation stage=gene_calling \
  platform=config/platforms/icdli1.yaml \
  --rerun-incomplete --keep-going
```

每批完成后会写入 `results/mag_annotation/bacteria_pure_culture_filter/gene_calling_batches/`；
如需调整每个 job 的规模，修改 `mag_annotation.yaml` 中的 `gene_calling.batch_size`。

确认第一阶段的 `proteins/all/` 已全部生成后，第二阶段只读取预先生成的蛋白目录：

```bash
snakemake --workflow-profile profiles/slurm \
  --config conf=projects/config_mag_annotation stage=annotation \
  platform=config/platforms/icdli1.yaml \
  --rerun-incomplete --keep-going
```

### 4. 完整闭合基因组的可选基因预测策略

当前配置使用 `auto`，适合混合的完整基因组和 draft genome。对于 FASTA header 明确带有
`circular-yes`、`closed-yes`、`complete-yes` 或 `topology=circular` 的完整闭合基因组，
可以启用自动策略：

```yaml
mag_annotation:
  gene_calling:
    mode: "auto"
    trans_table: 11
    auto:
      complete_mode: "single"
      fallback_mode: "meta"
```

`auto` 只有在所有保留 contig 都带有上述显式闭合标记时才选择 `single`；否则回退到
`meta`。闭合模式使用 genome-specific training、开放序列边界和 N 区域 masking，避免
完整闭合基因组在起始密码子选择上出现不必要的短 CDS。每个基因组的实际选择会写入
`pyrodigal` 日志，metadata 中记录最终解析后的 `single` 或 `meta` 模式。

如果一个输入目录中的所有基因组都已确认完整闭合，也可以直接使用：

```yaml
mag_annotation:
  gene_calling:
    mode: "single"
```

---

## 📊 核心产物文件结构

所有结果统一输出至 `results/mag_annotation/{sample}/`（当前示例为 `hadar_Ytt`）目录下：

```text
results/mag_annotation/hadar_Ytt/
├── proteins/
│   ├── all/                                      # Pyrodigal 统一预测的原始全量蛋白 (.faa.gz)
│   └── unique/
│       ├── unique_proteins.faa.gz                # 100% 严格去重后的唯一蛋白序列集
│       ├── protein_to_unique.tsv.gz              # 原始蛋白 -> 唯一蛋白映射表
│       ├── run_metadata.json                     # 统计背景空间 Z 与双规模溯源记录
│       └── dedup_stats.tsv                       # 去重比例与残基统计表
├── genes/                                        # 规范 GFF3 注释文件 (.gff.gz)
├── metadata/                                     # 包含 10 核心拓扑字段的单蛋白结构元数据表
├── shards/
│   ├── shard_0001.faa.gz                         # ~120M aa 残基均衡切片
│   ├── shard_list.txt                            # 切片编号列表
│   └── shards_manifest.tsv                       # 包含 mean, median, max 长度统计的分片清单
├── results/                                      # 原生并发计算结果（1 个 Shard 对应 1 个独立结果）
│   ├── shard_0001.tsv                            # PyKofamSearch 独立搜索输出
│   └── shard_0002.tsv
├── validation/
│   └── validation_report.json                    # 三级科学验证门禁评估报告
└── merged/
    ├── protein_KO_all_strict.tsv                 # 【核心一级资产】全部严格达标 KO 命中与位次表
    ├── protein_KO_best.tsv                       # 过滤单最佳 KO 衍生结果
    └── genome_KO_counts.long.tsv                 # 基因组级稀疏长表计数矩阵
```

同时，四层功能证据体系重构产物输出至 `results/pathway_reconstruction/{sample}/`：

```text
results/pathway_reconstruction/hadar_Ytt/
├── modules/                                      # [Tier 1] 模块重构 (complete / near_complete / partial)
│   ├── genome_module_completeness.tsv            # 全量 574 官方 + 自定义模块得分表 (含 definition_satisfied)
│   ├── genome_module_presence.tsv                # 模块定性有无判定
│   ├── module_evidence.tsv                       # 坐标级全链路证据追踪表
│   └── genome_module_matrix.tsv                  # 基因组 x 模块完整度矩阵大表
├── pathways/                                     # [Tier 2] 通路图谱网络投影
│   ├── genome_pathway_footprint.tsv              # map_ko_fraction, reaction_node_fraction 及 NA 标记长表
│   └── genome_pathway_matrix.tsv                 # 原核相关通路网络投影宽表
├── brite/                                        # [Tier 3] 分子机器与系统画像
│   ├── genome_brite_counts.tsv                   # 转运体/双组分/防御系统归一化拷贝数表 (copies_per_1000_proteins)
│   └── genome_brite_matrix.tsv                   # 核心系统丰度宽表
├── audit/                                        # [Global Audit] 覆盖度审计
│   └── genome_ko_coverage_audit.tsv              # 蛋白注释率、层级追踪率与真实孤儿 KO 审计表
└── capabilities/                                 # [Tier 4 & Unified Profiles]
    ├── genome_metabolism_capabilities.json       # [Schema v1.0] 核心能量代谢画像 (100% 格式路径零破坏)
    └── genome_functional_profile.json            # [Schema v2.0] 四层体系全域画像 (含 ED 家族与古菌特殊系统)
```

