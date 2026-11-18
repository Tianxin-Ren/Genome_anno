# KEGG 信息抓取与 KO 注释工具使用指南

本项目已整合了 `KEGG` 信息检索与注释工具，您可以在当前仓库目录下直接使用。该工具支持根据 KO 号或 Pathway Map 快速抓取对应的基因名、功能描述和 EC 编号，并支持智能本地缓存（7 天自动过期）以减少重复请求。

---

## 1. 结构与位置

- **主程序**: [kegg_ko_annotate.py](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/scripts/kegg_ko_annotate.py)
- **快捷包装脚本**: [kegg](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/scripts/kegg)
- **缓存目录**: `scripts/cache/` (运行后自动创建)
- **默认输出目录**: `scripts/data/output/` (运行后自动创建)

---

## 2. 三种运行模式

您可以在项目根目录下使用 `./scripts/kegg` 快捷启动（已设置可执行权限）。

### A. 通路下载模式 (`--map`)
获取指定的 KEGG Pathway Map 关联的所有 KO，并自动输出包含完整注释（定义、基因名、EC 编号）的表格。如果不指定 `--output`，默认会根据通路名称自动命名并保存至 `scripts/data/output/`。

```bash
# 下载并注释氧化磷酸化通路 (map00190) 的所有 KO
./scripts/kegg --map 00190

# 指定输出 file 路径
./scripts/kegg --map map00190 --output results/data_prep/phospho_kos.tsv
```

### B. 命令行直接检索模式 (`--ko`)
直接在命令行输入逗号分隔的 KO 号进行快速检索。

```bash
./scripts/kegg --ko K00399,K00398,K00397 --output results/data_prep/mcr_kos.tsv
```

### C. 传统文件注释模式 (`--input` 与 `--ko_column`)
读取已有的包含 KO 列的 TSV/CSV/Excel 文件，并在右侧追加 `definition`、`gene_name` 和 `ec_number` 注释列。

```bash
./scripts/kegg --input results/annotation/my_sample_kos.tsv --ko_column KO_ID
```

---

## 3. 智能缓存与更新

为保障请求速度并遵循 API 礼仪，该工具包含本地 JSON 缓存：
- 缓存有效期默认为 **7 天**（可以在运行命令时使用 `--expire-days` 自定义失效时间，设为 `-1` 永不过期）。
- 若想强制重新从 KEGG 服务器上同步最新的数据并更新本地缓存，可以在命令行加上 `--update-cache` 参数：

```bash
./scripts/kegg --map 00190 --update-cache
```

---

## 4. 与数据搜集模块 (data_prep) 联动

您可以先使用此工具检索出感兴趣的 KEGG 生物通路的所有基因名称，然后将这些基因配置到 [data_prep.yaml](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/config/data_prep.yaml) 中进行高可信度序列的自动搜集。

1. **查找通路基因**：
   ```bash
   ./scripts/kegg --map 00190 --output results/data_prep/photosynthesis_kos.tsv
   ```
2. **提取基因名填入配置**：
   打开生成的 TSV 文件，挑选其中需要的基因名称（`gene_name`）并配置到 `config/data_prep.yaml`：
   ```yaml
       "collect_photosynthesis_genes":
         genes:
           - name: "psbA"
             organism: "Cyanobacteria"
           - name: "psbD"
             organism: "Cyanobacteria"
         min_sequences: 5
   ```
3. **运行搜集流程**：
   ```bash
   conda run -n snakemake_9.21.0 snakemake -s Snakefile --config conf=config results/data_prep/collect_photosynthesis_genes/genes_collected.done
   ```
