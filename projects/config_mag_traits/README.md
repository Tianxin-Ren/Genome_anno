# 批量 MAG 生态表型与代谢特征分析流程 (config_mag_traits)

本配置工程 `projects/config_mag_traits/` 专为**超大批量 MAG（数千至上万个基因组）**量身打造，旨在以极高算力效率开展两轮去冗余、生态表型预测与特征代谢通路挖掘。

核心集成工具：
0. **[前置两轮去冗余] dRep** (`dRep dereplicate`)：
   - **第一轮 (Round 1)**：默认 **ANI 0.99**（亚种/菌株株系水平聚类），剔除近乎完全相同的克隆株。
   - **第二轮 (Round 2)**：默认 **ANI 0.95**（物种水平聚类），从 Round 1 胜出株中筛选出代表性菌株（SRGs）。
   - **CheckM 规避优化**：支持传入一次性预先计算的 `genomeInfo.csv`（`genome,completeness,contamination`），两轮均直接复用，完全跳过耗时极长的 CheckM 计算；未提供时自动回退为 `--ignoreGenomeQuality`。
1. **分批批次并行 (Batch Chunking)**：
   - 筛选出的代表性 MAGs 自动按每批 50 个 MAG 进行切割分发（`batch_size: 50`）。
   - **Prodigal (2 CPU)**：批次内串行生成 `.faa`、`.ffn`、`.gff`，避免成千上万小作业对 SLURM 调度系统的挤压。
   - **dbCAN (4 CPU)**：将批次内 50 个 MAG 的 `.faa` 合并为单个多序列文件（加头 `>{mag}___{gene_id}`），一次性调用 `run_dbcan`（Diamond + PyHMMER），避免 50 次载入 10GB+ 数据库的 IO 开销。
   - **CAMPER (4 CPU)**：多酚降解与酚酸转化酶系统（Polyphenol metabolism）特征注释，批次化 `hmmsearch` 统一检索。
   - **gRodon2 (1 CPU)**：启动单一 R 进程在内存中常驻，连续遍历处理批次内的 50 个 `.ffn` 文件，采用 `mode="partial"` 模式，彻底避免频繁启动 R 解释器的延时。
   - **GenomeSPOT (2 CPU)**：原核微生物理化生态表型机器学习预测（最适/最高/最低温度 OGT/Tmin/Tmax、最适/耐受 pH、最适/耐受盐度 Salinity、氧气耐受性及概率 Oxygen Tolerance/Prob），批次化执行与单菌容错隔离。
2. **终极多维特征矩阵大表整合 (Summary)**：
   - 自动映射还原两轮 dRep 聚类层级（原始 MAG -> R1 菌株聚类 -> R1 代表株 -> R2 种聚类 -> 最终代表株）。
   - 汇总 CAZymes 统计、多酚酶命中、gRodon2 倍增时间与生长预测、GenomeSPOT 生态温盐 pH 氧气耐受表型特征。

---

## 🔬 分析流程与数据流向

```text
[海量原始 MAG DNA 序列 (FASTA)]
         │
         ├─── (可选前置: steps.drep = true, two_rounds: true)
         │       │
         │       ▼
         │   0.1 Round 1 dRep (ANI 0.99 菌株水平去冗余, 复用 genomeInfo.csv)
         │       │
         │       ▼
         │   0.2 Round 2 dRep (ANI 0.95 种水平去冗余, 复用 winners genomeInfo)
         │       ├──> 获得最终代表性基因组 (SRGs, N_final)
         │       └──> 产出两轮聚类映射拓扑 (R1_Cdb.csv, R2_Cdb.csv)
         │               │
         ├───────────────┘
         ▼
  1. 批次分区 (Partition into 50-MAG batches)
         │
         ▼
  2. Batch Prodigal (2 CPU per batch)
         ├──> 批次合并 FAA (>{mag}___{gene}) ────────┬────────────────────────┬────────────────────────┐
         │                                          │                        │                        │
         ├──> 批次 FFN 目录 (50x {mag}.ffn)          ▼                        ▼                        ▼
         │         │                          3. Batch dbCAN (4 CPU)    4. Batch CAMPER (4 CPU)  5. Batch GenomeSPOT (2 CPU)
         │         │                          - 单次加载大库检索        - 多酚转化酶检索         - 温/盐/pH/氧气耐受性预测
         │         │                          - 自动拆分还原 MAG 计数   - 代谢模块汇总           - 单菌容错隔离机制
         │         ▼                                │                        │                        │
         │   6. Batch gRodon2 (1 CPU)               │                        │                        │
         │      - 单一 R 进程循环处理 50 个 FFN     │                        │                        │
         │      - mode="partial" 偏差矫正           │                        │                        │
         │      - 预测最大生长速率与倍增时间        │                        │                        │
         │         │                                │                        │                        │
         └─────────┴────────────────────────────────┴────────────────────────┴────────────────────────┘
                                           │
                                           ▼
                    7. 终极多维特征矩阵大表整合 (merge_mag_traits.py)
                       - 自动关联 R1/R2 两轮聚类归属与成员覆盖
                       results/mag_traits/summary/all_mag_traits_matrix.tsv
```

---

## 📂 配置文件与参数说明

| 文件 | 职责与配置要点 |
| :--- | :--- |
| **`pipeline.yaml`** | 主控文件。调度 `mag_traits` 模块开关及运行平台（已设为 `config/platforms/gses0.yaml`）。 |
| **`samples.yaml`** | 输入声明。提供 `project_root`、`mag_dir` 和可选的 `genome_info`。 |
| **`mag_traits.yaml`** | 工具科学参数、两轮去冗余阈值、批次切分大小（`batch_size: 50`）与各工具线程配置。 |

### 快速接入指南：
1. **配置项目根目录与 MAG 输入路径**：
   打开 `projects/config_mag_traits/samples.yaml`：
   ```yaml
   # 指定项目主目录（若提供，下方相对路径自动以此为基准拼接）
   project_root: "/dssg/home/acct-schhui/gses0/rentianxin"

   # 指定 MAG 存放目录（支持绝对路径，或基于 project_root 的相对路径）
   mag_dir: "data/gses_all_rawMAGs_fa_V0"

   # 可选：CheckM 质量评估结果（提供即可跳过耗时极长的内置 CheckM 计算）
   # CSV 格式要求头：genome,completeness,contamination
   genome_info: "data/checkm_summary.csv"
   ```

2. **配置两轮 dRep 与批次大小**：
   在 `projects/config_mag_traits/mag_traits.yaml` 中：
   ```yaml
   mag_traits:
     steps:
       drep: true            # 开启 dRep 去冗余
     drep:
       two_rounds: true      # 开启两轮去冗余
       round1_ani: 0.99      # 第一轮 0.99 菌株去冗余
       round2_ani: 0.95      # 第二轮 0.95 种水平去冗余
       threads: 64
     batching:
       batch_size: 50        # 每批 50 个 MAG 进行打包批次并行
     prodigal:
       threads: 2
     dbcan:
       threads: 4
     camper:
       threads: 4
     grodon2:
       threads: 1
       mode: "partial"       # 针对 MAG 的 partial 模式
       temperature: 20.0
     genomespot:
       threads: 2           # 每个批次 2 CPU 进行表型推断
       models_dir: ""        # 留空自动定位内置或平台 models
   ```

---

## 💻 集群提交与运行指令

### 提交到 gses0 平台 (搭配 profiles/slurm_pi)
```bash
snakemake --workflow-profile profiles/slurm_pi \
  --config conf=projects/config_mag_traits \
  --rerun-incomplete --keep-going
```

---

## 📊 核心产物文件结构

所有结果统一输出至 `results/mag_traits/` 目录下：

```text
results/mag_traits/
├── drep/                                         # 两轮去冗余产物
│   ├── round1/
│   │   ├── dereplicated_genomes/                 # Round 1 胜出 MAG (0.99 ANI)
│   │   ├── genomeInfo.csv                        # 标准化质量信息表
│   │   └── data_tables/Cdb.csv                   # Round 1 聚类拓扑
│   └── round2/
│       ├── dereplicated_genomes/                 # 最终代表性基因组 (SRGs, 0.95 ANI)
│       ├── genomeInfo.csv                        # 自动筛选过滤出的胜出株质量表
│       └── data_tables/Cdb.csv                   # Round 2 聚类拓扑
├── batches/                                      # 批次并行工作区 (每批 50 个 MAG)
│   ├── batch_list.txt                            # 批次编号列表 (batch_0001, ...)
│   └── {batch_id}/
│       ├── mag_list.tsv                          # 该批次包含的 50 个 MAG 映射
│       ├── faa/ & ffn/                           # 预测蛋白与编码核酸
│       ├── dbcan/overview.txt                    # 批次合并比对 CAZymes 概览
│       ├── camper/camper_hits.tsv                # 批次合并比对多酚酶命中
│       ├── grodon2/grodon2_growth.tsv            # 单一 R 进程连续估算批次生长速率
│       └── genomespot/                           # 批次 GenomeSPOT 理化表型预测
│           ├── {mag}.predictions.tsv             # 单菌预测详情 (OGT, pH, 盐度, 氧气)
│           └── genomespot_summary.tsv            # 批次 50 个 MAG 汇总表
└── summary/
    └── all_mag_traits_matrix.tsv                 # 【核心终极产物】全 MAG 多维特征宽表矩阵 (包含 CAZyme, 多酚, gRodon2, GenomeSPOT)
```
