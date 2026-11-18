# 基因簇保留与候选 Operon 结构分析 (`config_synteny`)

## 1. 流程背景与科学定位
在比较基因组学中，传统共线性展示（如 Clinker/Cblaster）主要用于“发现共线性”，但在回答核心进化问题：
> **“目标基因簇没注释出来，究竟是注释遗漏 (Annotation Miss)、组装断裂 (Assembly Break)、基因退化 (Pseudogenization)，还是发生了真正的演化丢失 (Evolutionary Gene Loss)？”**

`config_synteny` 结合参考基因组、侧翼同源基因、蛋白与 DNA 搜索及组装比对，评估目标基因簇在同源位点的保留状态。`operon_prediction=CANDIDATE_OPERON` 表示编码基因的排列支持同一转录单元假说；单凭基因组序列不能确认共转录，最终确认需转录组或转录边界证据。

---

## 2. 核心架构与 4 层门控 (Gates)

1. **Stage 0: 近缘参考筛选与验证**
   - 阳性参考初筛后，至少 40%（且至少 3 个）目标成员须有真实坐标并落在同一 50 kb 区域；KO 计数只用于筛选候选，不单独构成阳性证据。
   - 从经验证的阳性参考提取侧翼锚点。阴性参考还须通过双侧锚点闭合及全基因组 DNA 翻译同源搜索，才标记 `CONFIRMED_NEGATIVE`。
   - 参考位点共定位用于建立锚点；operon 的链向、顺序、间距及中间基因状态在 Stage 2 单独判断。缺少坐标或搜索失败时记为 provisional；没有确认阳性参考时，位点定位阶段输出 `NO_VALIDATED_POSITIVE_REF`，并保留候选参考的失败说明。
2. **Stage 1: 目标基因组组装可解析性检验 (Gate 1)**
   - 在目标基因组定位共识锚点，检测是否位于同一 Contig、计算观测锚点间隙 $\text{observed\_anchor\_gap}$ 与预期长度 $\text{expected\_locus\_length}$；
   - 依据基因簇尺寸动态计算边缘安全余量：$\text{required\_edge\_margin} = \max(\text{min\_edge\_bp}, \text{expected\_locus\_length} \times 1.5)$。
3. **Stage 2: 基因挽救与候选 Operon 结构**
   - 以 `blastp` 搜索蛋白，以位点 DNA / 全基因组 DNA 为 `blastx` 查询、目标蛋白为数据库进行翻译搜索。
   - 蛋白命中的局部性由真实 CDS 坐标决定；搜索失败或坐标未知返回未解析状态。
   - DIAMOND 2.0.7 之前不支持 `qseq_translated`；旧版会自动重试并将 DNA 命中的移码/终止密码子判定记为未评估，不据此宣称完整 ORF 或伪基因。
   - `cluster_state` 描述成员保留，`operon_structure_state` 和 `operon_prediction` 描述同向性、顺序、基因间距与中间基因。间距阈值由 `max_operon_intergenic_bp` 控制，建议用研究类群的已知 operon 校准。
   - 未核验密码子重编码或测序错误时，不自动抹去内部终止或移码证据。
4. **Stage 3: 位点物理证据**
   - 缺失断点须由位点内的 PAF 共线分段或单条比对中的大段 `D` CIGAR 操作支持；锚点距离缩短本身不算核酸断点证据。
   - 提供 BAM 时按实际位置统计跨候选断点且 MAPQ≥30 的 reads；没有断点时 read 支撑保持未取得。
5. **Stage 4: 演化解释**
   - 可为样本提供独立推断、已定根的 `species_tree_nwk`，或同时提供 `species_tree_outgroup`。未提供真实树时输出 `UNRESOLVED_PHYLOGENY`，不会构造固定树来宣称高置信度演化丢失。
   - `I_novel` 当前仅表示位点中未解释基因的占位比例。没有跨参考正交比较时 `novel_insertion_evidence=NOT_ASSESSED`，不能据此判为新插入。
   - `CLEAN_LOSS` 需要真实比对断点；只有距离缩短时使用 `LOSS_WITH_UNRESOLVED_STRUCTURE`，避免把结构变化当成已证实的整段缺失。

---

## 3. 配置文件目录结构 (`projects/config_synteny/`)

```text
projects/config_synteny/
├── pipeline.yaml            # 主调度开关（启用 synteny_validation）
├── synteny_validation.yaml  # 算法参数、数据库路径与环境配置
├── samples.yaml             # 目标样本元数据表
└── README.md                # 流程说明文档
```

---

## 4. 运行调度指令

在集群或工作站执行 Snakemake 调度：

```bash
# 激活 snakemake 环境
conda activate snakemake_9.21.0

# 启动 config_synteny 流程
# 集群统一使用 Slurm workflow profile 和平台配置，不再使用本地 -j 参数。
snakemake --workflow-profile profiles/slurm --config conf=projects/config_synteny platform=config/platforms/slurm.yaml

# 强制重跑 Clinker 可视化
snakemake --workflow-profile profiles/slurm --config conf=projects/config_synteny platform=config/platforms/slurm.yaml -R synteny_clinker_visualization
```

---

## 5. 核心输出文件说明

分析结果将统一保存在 `results/synteny_validation/` 下：

1. **最终全景评估总表** (`results/synteny_validation/locus_synteny_validation_summary.tsv`)
   包含成员基因保留、候选 operon 结构、位点断点证据和系统树状态；请结合 `mutation_evidence_status`、`species_tree_status`、`replacement_search_status` 与 `novel_insertion_evidence` 解读演化标签。
2. **各样本阶段明细表**：
   - `stage0/selected_validated_references.tsv`：通过验证的近缘阳性与阴性参考。
   - `stage0/consensus_anchors.tsv`：正交同源锚点。
   - `stage1/anchor_topology.tsv`：锚点空间拓扑与动态边缘余量。
   - `stage2/gene_level_rescue_summary.tsv`：逐基因三级挽救明细。
   - `stage3/locus_absence_validation.tsv`：物理位点状态判定。
   - `stage4/{sample}_local_tree.nwk`：局部系统发育树。
3. **交互式局部共线性矢量图** (`visualization/{sample}_clinker.html`)
   调用 Clinker 绘制目标位点与阳性参考的蛋白质同源比对连线图：先确定 target genes / 模块核心基因的完整范围，再从其左右边界各向外保留 `clinker_flank_genes` 个 CDS（默认 25 个）；可用 `clinker_context_bp` 添加额外 bp 上下文（默认 0）。统一侧翼基因数有助于让 ref/sample 轨道长度接近；实际长度仍受基因大小、核心区域跨度和 contig 边界影响。完整目标 DNA 切片仍供 Stage 2/3 使用，聚焦只影响可视化；核心区域之间的整段序列会保留。无连线的基因表示当前 Clinker 阈值下未检出可绘制的同源关系，不单独构成基因缺失证据。
