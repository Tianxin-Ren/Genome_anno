# MAG High-Throughput Functional Annotation Module (`mag_annotation`)

## 1. 架构与数据链设计

本模块专为 **GTDB r232 全量古菌（~22,343 基因组，35–55M 蛋白）**及超大规模环境 MAG 集合打造，基于 **PyKofamSearch + 序列化全量 KOfam 数据库 (`database.pkl`)**，在 HPC 集群 **8 节点 × 64 核心 × 512 GB RAM** 条件下实现生产级功能注释。

### 核心生物信息数据链路
$$\text{Genome} \longrightarrow \text{Gene} \longrightarrow \text{Contig / Strand / Coordinate Order} \longrightarrow \text{Protein} \longrightarrow \text{KO}$$

### 资产层级定义
- **一级核心资产 (Primary Artifact)**：`protein_KO_all_strict.tsv`
  保留所有通过 KOfam 自适应得分阈值（或 80% 相对放宽阈值）的候选 KO 命中，记录得分与自适应阈值的相对差额 (`score_margin`)、候选命中位次 (`ko_rank`)、严格判定标志 (`is_strict`) 以及阈值评分类型 (`score_type`)。与 `protein_metadata.tsv.gz` 结合，为后续 ATP 合成酶多亚基化学计量、操纵子共线性及基因邻接 (Gene Neighborhood) 分析提供确定性微观拓扑数据。
- **二级衍生矩阵 (Secondary Derivative)**：`protein_KO_best.tsv.gz` 及 `genome_KO_counts.long.tsv.gz`
  提供过滤为唯一最佳命中的基因组级稀疏计数长表，避免 5.75 亿单元格的超大稠密矩阵造成 I/O 负担。

---

## 2. 生产级优化机理

1. **节点整作业打包 (Node-Level Bundling)**：
   - 1 个 SLURM Job = 1 个完整的 `64c512g` 节点 (`threads: 64`, `mem_mb: ~500000`)；
   - 节点主进程将 `database.pkl` 从共享存储单次阶段化至 `$SLURM_TMPDIR/kofam.pkl` 并校验 SHA256；
   - 利用 Linux `multiprocessing.set_start_method('fork')`，**通过 Copy-on-Write (CoW) 共享只读数据库物理内存页面**，并发启动 5 个 12-thread worker，彻底保证 8 节点 × 5 进程 = 40 并发 = 480 核心主力计算。
2. **残基复杂度感知均衡切片 (Complexity-Aware Residue-Balanced Sharding)**：
   - 以 **~100–150 M amino acids (aa)** 为主要划分基准，**500,000 proteins** 为上限保护；
   - 输出包含 `mean_length`, `median_length`, `max_length` 的 `shards_manifest.tsv`，消除因超长蛋白造成的单 Shard 拖尾长尾效应。
3. **防碰撞精确去重 (Collision-Checked Exact Dedup)**：
   - 采用 SHA256 初筛 + 真实氨基酸字符二次比对，100% 保真且杜绝碰撞；
   - 在 `run_metadata.json` 中记录双统计背景：原始全量蛋白规模与 HMM 检索宇宙 $Z$，彻底解决审稿人对 E-value 背景空间的质疑。
4. **松弛预过滤参数与科学真实性**：
   - 设置 `prefilter_evalue: 10.0`，最大限度降低 HMMER 前序启发式过滤（MSV/Bias/Viterbi/Forward）对边界 bitscore 的截断风险。

---

## 3. 论文方法学标准文本 (Methods Template)

> **Methods: High-Throughput Functional Annotation of Archaeal Genomes**
> 
> "Predicted proteomes from GTDB release 232 archaeal genomes were deduplicated at 100% sequence identity using collision-checked exact hashing to establish a non-redundant protein search universe while preserving complete protein-to-genome topological mappings. The resulting unique sequences were partitioned into residue-balanced shards (~120M amino acids per shard). Functional annotations were assigned using PyKofamSearch against the full KOfam database (profiles and adaptive score thresholds derived from KEGG release 2026-08-25). Candidate hits were retained using relaxed HMMER heuristic pre-filtering (E-value ≤ 10.0) followed by KO-specific curated score threshold evaluation to minimize potential loss introduced by statistical pre-filtering. All hits satisfying the curated adaptive thresholds were retained in the primary annotation asset (`protein_KO_all_strict.tsv.gz`) with explicit reporting of bitscore, domain/full score threshold, margin, and candidate rank. Node-level worker placement utilized memory page-sharing under Linux fork semantics to optimize memory efficiency across 64-core compute allocations."

---

## 4. 核心产物文件结构

```text
results/mag_annotation/
├── proteins/
│   ├── all/                                      # [可选] Pyrodigal 预测原始全量蛋白 (.faa.gz)
│   └── unique/
│       ├── unique_proteins.faa.gz                # 100% 严格去重后的唯一蛋白序列集合
│       ├── protein_to_unique.tsv.gz              # 原始蛋白 -> 唯一蛋白全量映射表
│       ├── run_metadata.json                     # 统计背景空间 Z 与去重比例溯源
│       └── dedup_stats.tsv                       # 去重统计概要
├── genes/                                        # 基因结构注释 (.gff.gz)
├── metadata/                                     # 单蛋白规范拓扑元数据 (*_protein_metadata.tsv.gz)
├── shards/
│   ├── shard_0001.faa.gz                         # 残基均衡切片
│   ├── shard_list.txt                            # 切片清单
│   └── shards_manifest.tsv                       # 包含长度分布统计的分片元数据
├── node_batches/                                 # 单节点 5x12c 打包作业目录
│   ├── node_batch_0001.tsv
│   └── node_batch_0001.done
├── results/                                      # Shard 检索中间输出
│   ├── shard_0001.tsv
│   └── ...
├── validation/                                   # 三级验证门禁报告
│   └── validation_report.json
└── merged/
    ├── protein_KO_all_strict.tsv.gz              # 【核心一级资产】全部严格达标 KO 命中与位次
    ├── protein_KO_best.tsv.gz                    # 过滤最佳单个 KO 衍生结果
    ├── genome_KO_counts.long.tsv.gz              # 基因组级长格式稀疏计数矩阵
    └── genome_KO_counts.dense.tsv.gz             # [可选] 稠密计数宽表
```
