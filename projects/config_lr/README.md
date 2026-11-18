# 三代长读段 (ONT / HiFi) 组装与后续分析工程配置指引 (config_lr)

本配置文件夹 `projects/config_lr/` 专门针对 **Oxford Nanopore (ONT)** 及 **PacBio HiFi** 第三代长读段测序数据量身定制，集成了从原始读段质控、Myloasm 高精度从头组装、环状染色体/质粒自动提取，到全基因组与闭环染色体双轨功能注释、质粒病毒鉴定的全套标准分析流程。

针对 Bathy/u999856 与 Hadar/u499975 的实测长、短读长环化证据导出，见 [环化证据数据指南](../../docs/circular_evidence.md)。该指南区分论文方法与项目自定义判断阈值。

---

## 🔬 流程全景与核心机制 (Workflow Architecture)

本流水线整合了学术界三代测序通用分析范式与实验室自研的高精度闭环染色体解析技术，主要包含以下六大分析阶段：

```
[原始 ONT / HiFi Reads]
       │
       ▼
1. 读段质控与长度评估 (NanoPlot)
       │  產物: results/preprocess/{sample}/long_read/nanoplot/NanoPlot-report.html
       ▼
2. 长读段高分辨率从头组装 (Myloasm)
       │  產物: results/assembly/{sample}/long_read/myloasm/{sample}/{sample}.myloasm.fa
       │
       ├───────────────────────────────────────────────────────┐
       ▼                                                       ▼
3A. 全组装通用结构与功能注释 (Full Assembly Track)         3B. 闭环重叠群自动解析 (Circular Extraction Track)
    - 组装体分块并行 Bakta 基因预测                            - 正则解析 Myloasm 序列头成环状态与测序深度
    - EggNOG 进化谱系 (COG/KEGG/GO) 注释                       - 独立分离环状 Contigs (染色体/质粒)
    - Unified HMM 实验室特征能量代谢通路仲裁                     │  產物: results/binning/{sample}/circular_info.tsv
    - 综合功能大表 (*_functional_summary.tsv)                  │        results/binning/{sample}/circular_contigs/*.fa
       │                                                       ▼
       │                                              4. 闭环染色体分类与质量评估 (Taxonomy & QC)
       │                                                  - GTDB-Tk 分类学定种定属 (release232)
       │                                                  - CheckM2 预测完整度与污染度
       │                                                       │
       │                                                       ▼
       │                                              5. 闭环染色体精细化注释 (Refined Complete Annotation)
       │                                                  - Bakta (--complete 模式独立基因预测)
       │                                                  - 单菌级 EggNOG + HMM 映射与功能汇总
       │                                                  - 样本终极大表 (circular_summary.tsv)
       ▼                                                       ▼
6. 质粒与可移动元件端到端预测 (geNomad) ──────────────────────────┘
    產物: results/plasmid_analysis/{sample}/{sample}_genomad_summary.tsv
```

---

## 📂 配置文件结构说明

| 配置文件 | 功能职责与关键参数 |
| :--- | :--- |
| **`pipeline.yaml`** | **主控注册中心**：调度各模块的开启 (`enable: True/False`) 及全局 Conda/数据根目录定义。 |
| **`samples.yaml`** | **样本数据表**：声明 ONT 样本的长读段文件路径 (`long_read: lr: ...`)，可指定人工标记的 `selected_contigs`。 |
| **`preprocess.yaml`** | **三代质控配置**：启用 `long_qc: True`，自动执行 NanoPlot 质量评估。 |
| **`assembly.yaml`** | **Myloasm 组装配置**：启用 `myloasm: True`，支持传入 `extra` 测序化学版本参数（如 `--nano-r10` 等）。 |
| **`binning.yaml`** | **成环提取配置**：启用 `extract_circular: True`，根据 Myloasm 头信息自动抽取成环染色体与质粒。 |
| **`annotation.yaml`** | **双轨注释配置**：控制 Bakta、GTDB-Tk、CheckM2、EggNOG、Unified HMM 开关与数据库路径。 |
| **`plasmid_analysis.yaml`**| **质粒与病毒预测**：配置 geNomad 端到端模型和输出阈值。 |
| **`abundance.yaml`** | **定量配置 (可选)**：针对有配套二代 short reads 或 RNA-seq 的样本，计算基因 TPM 与覆盖度。 |

---

## 🚀 两种推荐分析模式与配置建议

### 模式 A：宏基因组 / 未完全成环全景分析 (Full Assembly Focus)
- **适用场景**：复杂群落宏基因组长读段、低覆盖度测序或包含较多线性片段（Linear Contigs）的样本。
- **配置建议**：
  - 在 `annotation.yaml` 中保持 `bakta: True`, `eggnog: True`, `advanced_hmm: True`。
  - 将 `circular: False`，或者保留 `circular: True` 但重点关注全组装输出的 `results/annotation/{sample}/bakta/` 与 `{sample}_functional_summary.tsv`。

### 模式 B：纯培养单菌 / 富集物高深度闭环精细化分析 (Circular Refinement Focus)
- **适用场景**：单菌分离株、高深度富集物（如 Bathyarchaeia 富集培养物），期望获得完整无缝闭环的染色体。
- **配置建议**：
  - 在 `annotation.yaml` 中保持 `circular: True`, `gtdbtk: True`, `checkm2: True`。
  - 可在 `samples.yaml` 中通过 `selected_contigs` 显式标记重点 contig（如 `u1287049: { complete: true }`），强制跳过长度筛选并赋予官方 `--complete` 标志。

---

## 💻 运行指令与集群调度

### 1. 提交至 SLURM 集群后台运行 (生产环境推荐)
通过 `--config conf=projects/config_lr` 指定本配置工程，并挂载对应的平台资源 profile：

```bash
snakemake --workflow-profile profiles/slurm \
  --config conf=projects/config_lr platform=config/platforms/slurm.yaml \
  --rerun-incomplete --keep-going
```

*（若在 `slurm_pi` 分区运行，可指定 `platform=config/platforms/slurm_pi.yaml`）*

### 2. 语法依赖编译校验 (本地 MacBook / 登录节点)
> [!IMPORTANT]
> 本地开发环境仅包含代码，不包含真实原始测序数据，严禁运行 `-n` 或 `--dry-run`。
> 本地语法审查请使用：
```bash
conda run -n snakemake_9.21.0 snakemake --lint
```

---

## 📊 核心产物文件路径一览

所有输出产物均保存在全局 `project_root`（默认为 `/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/results/`）下：

| 分析阶段 | 产物物理路径 | 内容说明 |
| :--- | :--- | :--- |
| **质控报告** | `results/preprocess/{sample}/long_read/nanoplot/NanoPlot-report.html` | 交互式长读段长度分布、N50 及 Q 值统计 |
| **三代组装体** | `results/assembly/{sample}/long_read/myloasm/{sample}/{sample}.myloasm.fa` | Myloasm 组装得到的全套 Primary Contigs (免截断) |
| **环状列表** | `results/binning/{sample}/circular_info.tsv` | 提取出的所有成环 Contig 元数据 (ID、长度、深度、多重度) |
| **单条环状 Contigs**| `results/binning/{sample}/circular_contigs/{contig_id}.fa` | 独立拆分出的单条闭环染色体/大质粒 FASTA |
| **物种分类鉴定** | `results/annotation/{sample}/circular_taxonomy/gtdbtk.summary.tsv` | GTDB-Tk 分类学谱系推断与参考基因组距离比对 |
| **完整度/污染度** | `results/annotation/{sample}/circular_quality/checkm2_summary.tsv` | CheckM2 评估的 Completeness (%) 与 Contamination (%) |
| **闭环基因组注释** | `results/annotation/{sample}/circular_refined/{contig_id}/` | 针对单条环状染色体的 Bakta GFF3/FAA/GBK/TSV 全套文件 |
| **KEGG 功能注释** | `results/annotation/{sample}/circular_refined/{contig_id}/functional/{contig_id}.kofam.tsv` | KofamScan 基于 Bakta 蛋白产出的官方 KO 命中表 |
| **物理比对质控** | `results/annotation/{sample}/circular_audit/{contig_id}/{contig_id}_audit.tsv` | Trigodet et al. 2026 四门控物理支撑与结构校验单 contig 报告 |
| **全样本闭环总表** | `results/annotation/{sample}/circular_summary.tsv` | 整合了 Gate 1~4 物理质检、分类、质量、功能注释率的终极大表 |
| **全组装功能表** | `results/annotation/{sample}/{sample}_functional_summary.tsv` | 全组装水平的 Unified HMM 与 EggNOG 功能整合表 |
| **质粒/病毒预测**| `results/plasmid_analysis/{sample}/{sample}_genomad_summary.tsv` | geNomad 预测出的质粒和噬菌体序列评分与分类 |

---

### 🛡️ Trigodet et al. 2026 (Nat Biotechnol) 四门控环状组装物理审计体系 (Version 2.1)

在长读长宏基因组组装中，组装软件（如 Flye、metaMDBG、HiCanu）常因转座子、长重复序列或端粒产生**过早假性成环（Premature Circularization）**、内部剪切断裂（Clipping）或重复嵌合。根据 *Trigodet et al. 2026 (Nature Biotechnology)* 的前沿评测与指导原则，本流水线彻底摒弃单纯依赖组装软件成环标签的做法，建立基于**双比对架构（Dual-Mapping Architecture）**与**原始三代读段（Raw ONT Reads）分子级物理支撑**的四门控（Four-Gate）科学验证体系：

### 1. 核心架构设计 (Dual-Mapping Architecture)
- **全局组装未旋转比对 (`whole_assembly.bam`)**：
  将所有 ONT 读段一次性比对回原始组装体（保留所有 linear 与 circular contigs 的真实坐标）。
  用于全长系统性截断排查（Systematic Clipping）、0× 物理空洞识别、多重映射追踪（`--secondary-seq`）以及伴侣序列捕获（`detect_methanothrix_partner_clustering`，识别如 *Methanothrix* 般与线性 contig 嵌合断裂的假环）。
- **靶向人造成环接头比对 (`junction.bam`)**：
  从候选 contig 两端提取侧翼序列（动态自适应：$\min(5000, L // 2)$，小质粒完全免疫假重叠），拼接为靶向接头 `[last_flank] + [first_flank]`。从原始测序数据中提取唯一分子标识读段（Deduplicated QNAMEs），利用 CIGAR 算子精准评估跨越点 $J$ 的物理跨越读段数（$\ge 5$ 条）、支持率（$\ge 0.50$，Wilson score 95% CI）及接头截断率（$\le 0.15$）。
- **覆盖度统计严谨解耦（Primary-only vs. All-alignment Support）**：
  - 测序深度、MAD、CV、1×/3× 覆盖广度采用 **Primary Alignment Only**，杜绝多重比对导致的深度虚高；
  - 0× 覆盖空洞检测采用 **All Valid Alignments**，只要存在真实比对即认定具有物理核苷酸支撑（忠实还原 Trigodet / anvi'o 定义），并支持原始读段 gap rescue；
- **BLASTN 全长自比对重复分析**：
  严格执行 Trigodet et al. 过滤标准：$\ge 80\%$ 相似度、长度 $>200$ bp，过滤自对角线后合并区间，精准量化非对角线重复覆盖比例（Repeat Breadth Pct）。
- **MIMAG 生物学标志物严谨判定**：
  - 必须同时独立检出全长 5S、16S 和 23S rRNA 基因；
  - 标准 20 种氨基酸同工受体归一化统计（`fMet` $\to$ `Met`, `Ile2` $\to$ `Ile`，剔除 `Sec`/`Pyl`/`Sup`/`Undet`），要求覆盖 $\ge 18/20$ 种标准氨基酸。

```
                              [待检环状 Contig]
                                      │
              ┌───────────────────────┴───────────────────────┐
              ▼                                               ▼
   [全局组装比对 whole_assembly.bam]              [靶向接头比对 junction.bam]
   - 全局深度/MAD/CV (Primary only)              - CIGAR 算子跨接头读段 (>= 5 reads)
   - 0x 覆盖空洞排查 (All alignments)            - 接头支持比例 (>= 50%, Wilson CI)
   - 系统性截断位点 (>80% clipping)              - 接头剪切率 (< 15%)
   - 伴侣 Contig 补充比对 (Methanothrix 排查)     - 末端比对全组装唯一性审查
              │                                               │
              └───────────────────────┬───────────────────────┘
                                      ▼
             [Gate 1: 物理结构支撑 (Physical Structural Support)]
              * 跨接头唯一分子读段 >= 5 且 支持率 >= 50%
              * 接头剪切率 <= 15% 且 末端比对具有全组装唯一性
              * 无 >500 bp 高危 0x 空洞且未解决 0x bp == 0
              * 系统性截断位点 == 0
              * 无统计学显著的多重伴侣 Contig 嵌合断裂 (Partner reads < 3 或 < 50%)
                                      │
                                      ├── [FAIL] ──► POTENTIAL_MISASSEMBLY /
                                      │              PREMATURE_CIRCULARIZATION_SUSPECT /
                                      │              INSUFFICIENT_PHYSICAL_EVIDENCE
                                      ▼ [PASS]
             [Gate 2: 覆盖度一致性与重复特征 (Structural & Coverage Coherence)]
              * 测序深度均值 (Mean Depth >= 10x) 与 1x 覆盖广度 (>= 99%)
              * BLASTN 自比对重复广度 (< 70% for <50kb)
                                      │
                                      ▼
             [Gate 3: 生物学完整性 (Biological Completeness, 染色体适用)]
              * CheckM2 完整度 > 90%, 污染度 < 5%
              * rRNA 独立基因完整性: 5S >= 1, 16S >= 1, 23S >= 1
              * 标准 20 种氨基酸 tRNA 同工受体多样性: >= 18 / 20
              *(质粒与病毒复制子自主豁免 Gate 3)*
                                      │
                                      ▼
             [Gate 4: 复制子类型证据 (Replicon Identity Evidence)]
              * geNomad 端到端染色体、质粒与病毒原始评分
              * 复制子置信判定与 Margin 置信差距
                                      │
                                      ▼
                       【最终置信度评级 (Confidence Grade)】
              1. READ_SUPPORTED_HIGH_QUALITY_CIRCULAR_CHROMOSOME (顶尖完整染色体)
              2. PHYSICALLY_VERIFIED_CIRCULAR_CHROMOSOME (物理验证染色体，标记物待补)
              3. HIGH_CONFIDENCE_CIRCULAR_PLASMID (高置信度物理成环质粒)
              4. READ_SUPPORTED_CIRCULAR_VIRAL_REPLICON (物理成环病毒/噬菌体复制子)
              5. PHYSICALLY_VERIFIED_CIRCULAR_REPLICON (物理验证中等环状复制子)
              6. PREMATURE_CIRCULARIZATION_SUSPECT (过早闭环高危嫌疑，与线性重叠群嵌合)
              7. POTENTIAL_MISASSEMBLY (存在系统性物理断裂或高危未解决空洞)
              8. INSUFFICIENT_PHYSICAL_EVIDENCE (深度过低或接头缺乏跨越读段)
              9. PUTATIVE_CIRCULAR_REPLICON (人工审查级备选)
```

---

## 📦 环境配置与数据库部署指引 (DSSG HPC)

### 1. Conda 虚拟环境清单与创建命令

在集群管理节点或交互节点中，请确保已创建以下 Conda 环境（环境目录统一位于 `conda_env_dir`，如 `/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/`）：

```bash
# 1. 物理审计比对、BLAST与解析环境 (必需)
conda create -y -n longread_audit -c bioconda -c conda-forge \
    python=3.10 minimap2 samtools blast seqkit pysam pandas numpy matplotlib
# (可选，用于高级双圈 Circos 渲染): pip install pycirclize

# 2. anvi'o 9 错误扫描独立验证环境 (推荐，Trigodet et al. 2026 同款架构)
# 支持 anvi-script-find-misassemblies 自动化结构错误侦测
conda create -y -n anvio_9 -c bioconda -c conda-forge \
    python=3.10 anvio=9.0

# 3. 基础与已有环境 (DSSG 集群已预装确认)
# - bakta_1.12.0: Bakta 全量基因预测
# - checkm2_1.0.2: CheckM2 机器学习质量预测
# - gtdbtk-2.4.0: GTDB-Tk 物种分类
# - kofamscan: KEGG 官方 HMM 扫描
# - biopython_1.85: 序列操作与格式转换
```

> [!TIP]
> **双引擎容错与独立验证设计**：
> - 流程核心裁决引擎为 `scripts/audit_circular_assembly.py`，依托 `longread_audit` 环境运行，具备自给自足的 CIGAR 跨接头分析、伴侣 contig 聚类与 0x 空洞扫描能力；
> - `anvio_9` 环境并行运行 `anvi-script-find-misassemblies --just-do-it` 生成独立校验产物 `results/annotation/{sample}/circular_audit/{contig_id}/anvio_misassemblies.txt`，为科研复核提供双重证据支持。

### 2. 数据库配置检查

物理审计 Gate 1（Trigodet 框架）以待测环状 Contig 为参考序列、原始 ONT 读段为比对源，**无需额外下载外部数据库**！
下游 Gate 3 与注释使用的数据库已挂载在共享存储：
- **Bakta DB**: `/dssg/home/acct-clswfp/clswfp-icdli1/Database/Bakta/db-light`
- **CheckM2 DB**: `/dssg/home/acct-clswfp/clswfp-icdli1/Database/CheckM2_database/uniref100.KO.1.dmnd`
- **GTDB-Tk DB**: `/dssg/home/acct-clswfp/clswfp-icdli1/Database/GTDB/archive/gtdb-tk/release232`
- **KofamScan DB & ko_list**: `/dssg/home/acct-clswfp/clswfp-icdli1/Database/kofam`
