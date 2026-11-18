# Read Profiling Module (Sylph Taxonomic Profiling)

本模块专为高通量宏基因组测序短读长（Short-read Reads）群落物种丰度分析而设计，集成了超快速、低内存占用的物种分类丰度分析工具 **Sylph**。

## 功能特性
1. **多样本批次一键分析 (Batch/Dataset Mode)**：
   直接支持通配符模式匹配输入目录下所有样本质控后的洁净 reads（例如 `*/00.QC/*.clean.1.fq.gz` 与 `*/00.QC/*.clean.2.fq.gz`），复用加载至内存的 GTDB `.syldb` 数据库，极大减少 I/O 与重复建构开销。
2. **单样本独立分析与自动合并 (Per-sample Mode)**：
   支持按样本配置 `samples.yaml`，分别生成各个样本的丰度文件并自动合并为全局丰度矩阵。
3. **完全兼容 SLURM 调度规范**：
   动态资源申请，自动配置 `7d` 超时保护，支持统一 Conda 环境执行（`sylph`）。

## 配置文件说明
- [config.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/modules/read_profiling/config.yaml): 模块默认参数与数据库路径。
- 项目级配置位于 `projects/config_read/`。
