#!/usr/bin/env bash
# ==============================================================================
# diagnose_abundance_counts.sh
#
# 一键排查 DNA / RNA 计数全为 0 的全流程诊断脚本
# 用法:
#   bash scripts/diagnose_abundance_counts.sh [样本名] [结果目录]
# 示例:
#   bash scripts/diagnose_abundance_counts.sh
#   bash scripts/diagnose_abundance_counts.sh C1_0-0.5
# ==============================================================================

set -u

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

echo -e "${BOLD}${BLUE}================================================================${NC}"
echo -e "${BOLD}${BLUE}       MAG Abundance 基因丰度与计数全为 0 一键排查诊断工具       ${NC}"
echo -e "${BOLD}${BLUE}================================================================${NC}"

# 1. 定位结果目录
RES_DIR="${2:-}"
if [ -z "$RES_DIR" ]; then
    if [ -d "results/abundance/NS_2022" ]; then
        RES_DIR="results/abundance/NS_2022"
    elif [ -d "results/abundance" ]; then
        RES_DIR=$(find results/abundance -mindepth 1 -maxdepth 1 -type d | head -n 1)
    else
        RES_DIR="."
    fi
fi

echo -e "📂 分析结果主目录: ${BOLD}${RES_DIR}${NC}"

# 2. 定位样本目录
SAMPLE="${1:-}"
if [ -z "$SAMPLE" ]; then
    # 自动发现第一个样本文件夹
    SAMPLE_DIR=$(find "$RES_DIR" -mindepth 1 -maxdepth 1 -type d ! -name "references" ! -name "genome_abundance" ! -name "logs" | head -n 1)
    if [ -n "$SAMPLE_DIR" ] && [ -d "$SAMPLE_DIR" ]; then
        SAMPLE=$(basename "$SAMPLE_DIR")
    fi
fi

if [ -z "$SAMPLE" ]; then
    echo -e "${RED}[ERROR] 未在 ${RES_DIR} 下找到任何样本目录！请检查是否已运行流程。${NC}"
    exit 1
fi

SAMPLE_DIR="${RES_DIR}/${SAMPLE}"
echo -e "🧪 目标排查样本: ${BOLD}${SAMPLE}${NC} (路径: ${SAMPLE_DIR})"
echo ""

# ------------------------------------------------------------------------------
# 检查项 1: 上游单样本文件与中间文件是否存在真实的非零计数
# ------------------------------------------------------------------------------
echo -e "${BOLD}▶ [排查 1/5] 上游单样本结果 (final.tsv / featureCounts) 是否存在非零计数？${NC}"
DNA_FINAL="${SAMPLE_DIR}/dna_quant/${SAMPLE}.final.tsv"
RNA_FINAL="${SAMPLE_DIR}/rna_quant/${SAMPLE}.final.tsv"
DNA_COUNTS_FILE="${SAMPLE_DIR}/dna_quant/${SAMPLE}.featureCounts.tsv"

NONZERO_DNA=0
NONZERO_RNA=0

if [ -f "$DNA_FINAL" ]; then
    NONZERO_DNA=$(awk -F'\t' 'NR>1 && $8>0 {cnt++} END {print cnt+0}' "$DNA_FINAL" 2>/dev/null || echo 0)
    echo -e "  - DNA final.tsv 统计: 共有 ${BOLD}${NONZERO_DNA}${NC} 个基因 Counts > 0"
    if [ "$NONZERO_DNA" -gt 0 ]; then
        echo -e "    ${GREEN}✔ DNA 计数正常！非零条目示例:${NC}"
        awk -F'\t' 'NR>1 && $8>0 {print "      " $2 "\tCounts=" $8 "\tTPM=" $13}' "$DNA_FINAL" | head -n 3
    fi
else
    echo -e "  - ${YELLOW}未找到 DNA final.tsv: ${DNA_FINAL}${NC}"
fi

if [ -f "$RNA_FINAL" ]; then
    NONZERO_RNA=$(awk -F'\t' 'NR>1 && $8>0 {cnt++} END {print cnt+0}' "$RNA_FINAL" 2>/dev/null || echo 0)
    echo -e "  - RNA final.tsv 统计: 共有 ${BOLD}${NONZERO_RNA}${NC} 个基因 Counts > 0"
    if [ "$NONZERO_RNA" -gt 0 ]; then
        echo -e "    ${GREEN}✔ RNA 计数正常！非零条目示例:${NC}"
        awk -F'\t' 'NR>1 && $8>0 {print "      " $2 "\tCounts=" $8 "\tTPM=" $13}' "$RNA_FINAL" | head -n 3
    fi
fi

if [ "$NONZERO_DNA" -gt 0 ] || [ "$NONZERO_RNA" -gt 0 ]; then
    echo -e "\n  ${GREEN}【关键发现】${BOLD}上游单样本 featureCounts / final.tsv 已经成功比对并计数！${NC}"
    echo -e "  之所以汇总大表（gene_abundance_and_annotation.tsv）显示为 0，是因为旧版汇总脚本读取合并表列名时的 bug（现已修复）。"
    echo -e "  👉 ${BOLD}解决方案：无需重新比对，仅需重新运行汇总规则即可恢复数据：${NC}"
    echo -e "     snakemake --snakefile modules/abundance/abundance.smk --configfile projects/config_abun/pipeline.yaml -R abundance_summarize_by_gene -j 1\n"
else
    echo -e "  ${RED}✖ 警告: 上游单样本 final.tsv 中计数确实全为 0，继续向下排查...${NC}\n"
fi

# ------------------------------------------------------------------------------
# 检查项 2: featureCounts 运行摘要统计 (Assigned vs Unassigned)
# ------------------------------------------------------------------------------
echo -e "${BOLD}▶ [排查 2/5] featureCounts 读段分配状态 (.featureCounts.tsv.summary)${NC}"

check_summary() {
    local TYPE="$1"
    local SUM_FILE="${SAMPLE_DIR}/${TYPE}_quant/${SAMPLE}.featureCounts.tsv.summary"
    if [ -f "$SUM_FILE" ]; then
        echo -e "  ${BOLD}[${TYPE} 模式分配状态]:${NC}"
        cat "$SUM_FILE" | sed 's/^/    /'
        local ASSIGNED=$(awk '$1=="Assigned" {print $2}' "$SUM_FILE")
        local NO_FEAT=$(awk '$1=="Unassigned_NoFeatures" {print $2}' "$SUM_FILE")
        local UNMAPPED=$(awk '$1=="Unassigned_Unmapped" {print $2}' "$SUM_FILE")
        
        if [ "${ASSIGNED:-0}" -eq 0 ]; then
            if [ "${NO_FEAT:-0}" -gt 0 ]; then
                echo -e "    ${RED}✖ 诊断提示: Assigned 为 0，而 Unassigned_NoFeatures 极高！${NC}"
                echo -e "      原因很可能是 BAM 参考序列名称与 GFF 中的 Contig 名称不一致，或者 GFF 缺少 featureCounts 所需的特征/标签！"
            elif [ "${UNMAPPED:-0}" -gt 0 ]; then
                echo -e "    ${RED}✖ 诊断提示: Assigned 为 0，全部落在 Unassigned_Unmapped！${NC}"
                echo -e "      原因通常是比对率极低，或者 Reads 与 MAG 基因组完全没有同源性！"
            fi
        else
            echo -e "    ${GREEN}✔ 成功分配了 ${ASSIGNED} 条 Reads！${NC}"
        fi
    else
        echo -e "  ${YELLOW}未找到 ${TYPE} summary 文件: ${SUM_FILE}${NC}"
    fi
}

check_summary "dna"
check_summary "rna"
echo ""

# ------------------------------------------------------------------------------
# 检查项 3: BAM 与 GFF 染色体 / Contig 名称是否一致
# ------------------------------------------------------------------------------
echo -e "${BOLD}▶ [排查 3/5] BAM 头部参考序列名称 与 GFF 第一列是否精准吻合？${NC}"
BAM_FILE="${SAMPLE_DIR}/dna_quant/${SAMPLE}.sorted.bam"
[ ! -f "$BAM_FILE" ] && BAM_FILE="${SAMPLE_DIR}/rna_quant/${SAMPLE}.sorted.bam"

GFF_FILE=$(find "$RES_DIR" -maxdepth 1 -name "*.gff" | head -n 1)

if [ -f "$BAM_FILE" ] && command -v samtools &>/dev/null; then
    echo -e "  - 从 BAM 头部提取的前 5 个参考序列名 (@SQ):"
    samtools view -H "$BAM_FILE" 2>/dev/null | grep '^@SQ' | head -n 5 | awk '{print "      " $2}' | sed 's/SN://'
    
    if [ -n "$GFF_FILE" ] && [ -f "$GFF_FILE" ]; then
        echo -e "  - 从 GFF 文件 (${GFF_FILE}) 提取的前 5 个参考序列名 (第 1 列):"
        grep -v '^#' "$GFF_FILE" | head -n 5 | awk '{print "      " $1}'
        
        # 交叉验证
        SAMPLE_BAM_CHR=$(samtools view -H "$BAM_FILE" 2>/dev/null | grep '^@SQ' | head -n 1 | awk '{print $2}' | sed 's/SN://')
        if grep -q -w "$SAMPLE_BAM_CHR" "$GFF_FILE" 2>/dev/null; then
            echo -e "  ${GREEN}✔ 命名一致性验证通过！BAM 参考序列 ($SAMPLE_BAM_CHR) 在 GFF 中存在。${NC}"
        else
            echo -e "  ${RED}✖ 命名不匹配！BAM 参考序列 ($SAMPLE_BAM_CHR) 未在 GFF 第 1 列中找到！${NC}"
            echo -e "    这会导致 featureCounts 无法将 Reads 与基因关联，从而导致 Reads 计数全为 0。"
        fi
    fi
else
    echo -e "  ${YELLOW}无法解析 BAM (请确认已安装 samtools 且 BAM 文件存在)${NC}"
fi
echo ""

# ------------------------------------------------------------------------------
# 检查项 4: GFF 第三列特征类型与第九列属性标签
# ------------------------------------------------------------------------------
echo -e "${BOLD}▶ [排查 4/5] GFF 特征类型 (第 3 列) 与属性标签 (第 9 列)${NC}"
if [ -n "$GFF_FILE" ] && [ -f "$GFF_FILE" ]; then
    echo -e "  - GFF 中出现的特征类型统计 (第 3 列):"
    grep -v '^#' "$GFF_FILE" | awk '{print $3}' | sort | uniq -c | head -n 5 | sed 's/^/    /'
    
    echo -e "  - GFF 前 3 行属性列预览 (第 9 列):"
    grep -v '^#' "$GFF_FILE" | head -n 3 | cut -f 9 | sed 's/^/      /'
    
    if grep -v '^#' "$GFF_FILE" | awk '{print $3}' | grep -q -w "CDS"; then
        echo -e "  ${GREEN}✔ 存在 'CDS' 特征，与 featureCounts 参数 (-t CDS) 匹配。${NC}"
    else
        echo -e "  ${RED}✖ 未发现 'CDS' 特征！featureCounts 默认查找 -t CDS 会匹配不到任何基因！${NC}"
    fi
else
    echo -e "  ${YELLOW}未找到 GFF 文件。${NC}"
fi
echo ""

# ------------------------------------------------------------------------------
# 检查项 5: 比对日志与比对率 (Alignment Rate)
# ------------------------------------------------------------------------------
echo -e "${BOLD}▶ [排查 5/5] Mapping 日志中的整体比对率 (Alignment Rate)${NC}"
DNA_MAP_LOG="logs/abundance/NS_2022/dna_mapping/${SAMPLE}.log"
[ ! -f "$DNA_MAP_LOG" ] && DNA_MAP_LOG=$(find logs/abundance -name "*${SAMPLE}*.log" | grep -i "dna_mapping" | head -n 1)

if [ -f "$DNA_MAP_LOG" ]; then
    echo -e "  - DNA Bowtie2 比对日志结尾摘要 (${DNA_MAP_LOG}):"
    tail -n 10 "$DNA_MAP_LOG" | sed 's/^/      /'
else
    echo -e "  ${YELLOW}未找到 DNA mapping 日志: ${DNA_MAP_LOG}${NC}"
fi

RNA_MAP_LOG="logs/abundance/NS_2022/rna_mapping/${SAMPLE}.log"
[ ! -f "$RNA_MAP_LOG" ] && RNA_MAP_LOG=$(find logs/abundance -name "*${SAMPLE}*.log" | grep -i "rna_mapping" | head -n 1)

if [ -f "$RNA_MAP_LOG" ]; then
    echo -e "  - RNA HISAT2 比对日志结尾摘要 (${RNA_MAP_LOG}):"
    tail -n 10 "$RNA_MAP_LOG" | sed 's/^/      /'
fi

echo -e "\n${BOLD}${BLUE}================================================================${NC}"
echo -e "${BOLD}${BLUE}                         排查总结与建议                         ${NC}"
echo -e "${BOLD}${BLUE}================================================================${NC}"
echo -e "1. 若 [排查 1] 中已显示有非零条目，说明比对与定量均已成功，直接重新执行：\n"
echo -e "   ${BOLD}snakemake --snakefile modules/abundance/abundance.smk \\"
echo -e "             --configfile projects/config_abun/pipeline.yaml \\"
echo -e "             -R abundance_summarize_by_gene -j 1${NC}\n"
echo -e "2. 若 [排查 2] 中 Assigned 为 0 且 Unassigned_NoFeatures 极高，请根据 [排查 3/4] 检查 Contig 命名一致性。"
echo -e "3. 若 [排查 5] 中整体比对率接近 0.00%，请检查干净 reads 输入文件路径及数据质量。"
echo -e "${BOLD}${BLUE}================================================================${NC}"
