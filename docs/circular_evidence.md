# Bathy / Hadar 环化证据数据导出

目标在 [config_lr/annotation.yaml](../projects/config_lr/annotation.yaml) 的 `annotation.circular_evidence.targets` 中指定：

| 样本 | contig |
| --- | --- |
| Bathy_enrichment_2025 | u999856 |
| Hadar_Ytt | u499975 |

本目标导出实测数据供后续绘图和人工判断，不自动认定基因组环化可信。

用户已确认 Bathy 长/短 reads 来自**同一富集体系、不同时间取样**；已记录在 samples.yaml 并传入 manifest。短读长属于跨时间补充证据，缺少支持可能反映群体变化，不能单独否定长读长组装。Hadar 的取样关系暂未确认。

## 方法依据与适用范围

参考 Trigodet 等的 [Troubleshooting common errors in assemblies of long-read metagenomes](https://doi.org/10.1038/s41587-025-02971-8) 及[作者复现流程](https://merenlab.org/data/benchmarking-long-read-assemblers/)。保留原始长 reads 对完整组装的比对，检查零覆盖、系统性 clipping 与重复序列；原始组装中其他 contig 作为竞争参考，不能先只取目标基因组再作全长比对。

论文使用 HiFi 数据。本项目使用配置中的 ONT `map-ont`，不照搬 HiFi 错误率。论文复现流程保留同等得分的 secondary 比对（`-p 1 --secondary-seq`）：这里同时导出 **all-alignments 支持**与 **MAPQ ≥20 的 primary 深度**。前者用来寻找真正没有比对支持的位置，不能当作独立分子丰度；仅在过滤后出现的空洞不能直接判为组装错误。

clipping 表保留所有有末端 soft/hard clip 的位置及比例，不先按阈值删除。论文高置信事件重点为 100% clipping 且至少 10× 局部覆盖；[anvi'o 工具](https://anvio.org/help/main/programs/anvi-script-find-misassemblies/) 默认筛查比例 >0.8，并排除距 contig 端部 100 bp 内的位置。本导出器按左右端分别计数，坐标锚定在首/末实际匹配碱基，**不是 anvi'o 的逐字段复刻**。不要将相邻或左右 clipping 比例直接相加。

论文对 >500 bp 无 apparent coverage 的区间另作原始 reads BLAST 回溯；>1 kb 是文中另一统计尺度。这里只输出候选区间；发现长读长 all-alignments 内部空洞后，仍需按作者流程用 `blastn -dust no` 回查原始 reads，排除 mapper 失败。重复比对采用 BLASTN self-alignment；绘图时排除自对角线，筛选 identity ≥80%、alignment length >200 bp。重复本身不等于错误。

合成接头长读长回贴、旋转参考上的 Illumina 跨接/成对桥接是本项目补充。长读长默认两侧 anchor ≥500 bp，短读长 ≥30 bp；这些值以及旧 audit 的“至少5条 reads/支持率50%”均为项目阈值，不是论文给出的通用环化标准。

## 在服务器运行

在 Genome_anno 代码目录执行。配置的结果根目录为 `/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin`，如实际部署不同，应修改 `samples.yaml` 或通过 `project_root` 覆盖。

```bash
# 先检查任务清单
snakemake annotation_circular_evidence_targets \
  --profile profiles/slurm --dry-run --notemp \
  --config conf=projects/config_lr platform=config/platforms/icdli1.yaml

# 两个目标一起执行
snakemake annotation_circular_evidence_targets \
  --profile profiles/slurm --notemp \
  --config conf=projects/config_lr platform=config/platforms/icdli1.yaml
```

`--notemp` 保留短读长 BAM/BAI 和旋转参考，便于 IGV 复核。已有且被 Snakemake 判为最新的中间产物会复用；缺失的比对需要重跑。该目标只依赖提取、比对、self-BLAST 和数据导出，不依赖 Bakta、GTDB-Tk、CheckM2 或旧绘图脚本。若上游组装本身尚未生成，Snakemake 仍会追溯到组装依赖。

运行环境复用 `annotation.conda_envs.longread_audit`（需 Python、pysam、minimap2、samtools、BLAST）。self-BLAST 缺软件时直接失败，不用其他算法或空文件冒充完成。

## 结果及绘图用途

两个结果目录分别为项目根下：

```text
results/annotation/Bathy_enrichment_2025/circular_evidence/u999856/
results/annotation/Hadar_Ytt/circular_evidence/u499975/
```

| 文件 | 内容与用途 |
| --- | --- |
| `genome_coverage.tsv.gz` | 每1 kb长/短 reads 深度、零覆盖碱基数、覆盖比例；同时含 primary/MAPQ过滤与 all-alignments 两种口径，画全长覆盖轨道 |
| `junction_depth.tsv.gz` | 原始 L→1 接头附近默认 ±5 kb 逐碱基深度；`track` 区分 long/short，画接头覆盖曲线 |
| `junction_reads.tsv.gz` | 窗口内全部 alignment 的 read名、长度、坐标、方向、MAPQ、CIGAR、flags、mate/TLEN、过滤与跨接/桥接状态，画 reads 排列图；不是抽样数据 |
| `evidence_summary.tsv` | 全长平均深度/覆盖比例、接头窗口深度/覆盖比例、跨接独立 QNAME 数、短读 fragment 桥接数 |
| `zero_coverage_regions.tsv` | 全部精确零覆盖区间，按长短 reads 和过滤口径分列，画缺口轨道或提取待回查序列 |
| `long_clipping.tsv` | 长 reads 全组装比对的左右 clipping 位点、alignment数、当地支持深度、比例、离末端距离 |
| `repeat_self_blast.tsv` | BLAST outfmt 6 的原始12列，无表头；可筛选非对角线重复并作 dotplot |
| `manifest.json` | 参数、过滤口径、输入路径/大小/mtime、FASTA序列SHA256、坐标定义和限制；成功导出后最后写入 |

建议下载整个目录。BAM 仍留服务器供复核；需要 IGV 时还应使用相应参考 FASTA 和 BAM/BAI。

## 坐标和计数口径

- 导出的区间统一为 **0-based、左闭右开**。`relative_pos=0` 是接头右侧第一个碱基，`-1` 是原始 contig 最后一个碱基；`original_pos_0based` 可映射回原始 FASTA。
- 例外是原样保留的 `repeat_self_blast.tsv`：BLAST 的 qstart/qend/sstart/send 为 **1-based、两端包含**，反向比对的起止坐标可能递减；转换为区间时使用 `[min(start,end)-1, max(start,end))`。
- 长读接头是 `[last_flank]+[first_flank]`，边界取实际 synthetic reference 长度的一半。短读参考为 `seq[offset:]+seq[:offset]`，接头边界为 `L-offset`；包括奇数长度，不再减1。
- depth 只累计 CIGAR `M/=/X`，不将 `D/N` 当作碱基覆盖。长读跨接允许通常的测序 indel：每侧anchor至少70%实际匹配，窗口内不含N或>50 bp的indel；这是项目补充口径。另导出 `spans_junction_strict` / `strict_spanning_unique_qnames`，要求两个anchor处于同一参考连续匹配块（允许I，不允许D/N）；短读长使用这一严格口径。
- primary 深度排除 unmapped/secondary/supplementary/QC-fail/duplicate，并要求 MAPQ ≥20。paired reads 重叠区按两个 reads 计深度；跨接按 QNAME 去重。
- 短读桥接要求实际观察到两端、两端通过过滤、proper pair、同参考、FR方向、互相一致的mate位置和TLEN，fragment长50–2000 bp。非标准文库需重新设定口径；同一fragment的spanning与bridging计数不能相加。
- `query_length` 是 BAM 中 read 序列长度，不含硬剪切部分，不一定等于 FASTQ 原始长度。
- 全长轨道使用原始线性坐标，其端点可能有截断效应，应结合接头轨道解释。长读合成接头 BAM 的 MAPQ 只反映该局部参考，不能证明其在全组装中的唯一性。

## 判断与后续复核

优先检查接头是否有多个独立长 reads 连续跨越、两侧 anchor 是否落在重复区，再看全长内部是否出现一致 clipping、all-alignments 零覆盖或深度突变。短读长应检查接头深度相对邻近侧翼是否骤降，以及跨接 read / 成对 fragment 是否支持。短读长低深度或不同时间取样的群体变化，应记为证据不足或样本差异，不能直接推翻长读长结果。

目前旧 `plot_circular_audit.py` 在缺失 BAM 且 TSV 有均值时存在模拟覆盖 fallback。本导出目标不调用它；后续绘图应使用这里的实测 TSV。缺失输入会报错，真实空比对会导出0及空reads表，不能画成均匀覆盖或自动写“verified”。
