#!/usr/bin/env bash
# ==============================================================================
# 单样本独立共线性验证与基因丢失推断脚本 (run_single_synteny.sh)
# 作用：脱离 Snakemake，直接调用各 Stage Python 脚本完成单样本四大 ATP 合成酶模块端到端验证与 Clinker 绘图
# 用法: bash scripts/run_single_synteny.sh <TARGET_ACC> [MODULE(default:all)|M00159|M00157|M00158|M00160] [CUSTOM_TAXONOMY] [--clinker]
# 选项:
#   --force, -f      强制重跑所有步骤 (默认: 仅运行缺失/未完成的步骤)
#   --clinker, -c    开启 Clinker 可视化 HTML 对比图渲染 (默认: 关闭)
#   --no-clinker     显式关闭 Clinker 渲染
#   --threshold=NUM  设置 Contig 共有性收敛阈值 (默认: 0.75)
#   环境变量: FORCE=true/false, RUN_CLINKER=true/false, CONSENSUS_THRESHOLD=NUM
# ==============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FORCE="${FORCE:-false}"
RUN_CLINKER="${RUN_CLINKER:-false}"
CONSENSUS_THRESHOLD="${CONSENSUS_THRESHOLD:-0.75}"
SIGNIFICANCE_MODE="${SIGNIFICANCE_MODE:-empirical_max}"
SIGNIFICANCE_ALPHA="${SIGNIFICANCE_ALPHA:-0.01}"
N_PERMUTATIONS="${N_PERMUTATIONS:-1000}"
ALLOW_INSUFFICIENT="${ALLOW_INSUFFICIENT:-false}"
POSITIONAL_ARGS=()
DOMAIN_ARG=""
PRESET_ARG=""

PYTHON_BIN="python3"
if ! command -v python3 &>/dev/null || ! python3 -c "import yaml" &>/dev/null; then
    if [ -x "/opt/homebrew/Caskroom/miniconda/base/bin/python" ]; then
        PYTHON_BIN="/opt/homebrew/Caskroom/miniconda/base/bin/python"
    elif [ -x "$HOME/miniconda3/bin/python" ]; then
        PYTHON_BIN="$HOME/miniconda3/bin/python"
    fi
fi

for arg in "$@"; do
    case "$arg" in
        --force|-f)
            FORCE="true"
            ;;
        --clinker|-c)
            RUN_CLINKER="true"
            ;;
        --no-clinker)
            RUN_CLINKER="false"
            ;;
        --threshold=*|--consensus_contig_threshold=*)
            CONSENSUS_THRESHOLD="${arg#*=}"
            ;;
        --significance_mode=*|--mode=*)
            SIGNIFICANCE_MODE="${arg#*=}"
            ;;
        --alpha=*|--significance_alpha=*)
            SIGNIFICANCE_ALPHA="${arg#*=}"
            ;;
        --preset=*|--active_preset=*)
            PRESET_ARG="${arg#*=}"
            ;;
        --domain=*|--organism_domain=*)
            DOMAIN_ARG="${arg#*=}"
            ;;
        --perms=*|--n_permutations=*)
            N_PERMUTATIONS="${arg#*=}"
            ;;
        --allow_insufficient_evidence_locus)
            ALLOW_INSUFFICIENT="true"
            ;;
        *)
            POSITIONAL_ARGS+=("$arg")
            ;;
    esac
done

TARGET="${POSITIONAL_ARGS[0]:-GCA_016839045.1}"
CLEAN_TARGET=$(echo "$TARGET" | sed 's/^GB_//' | sed 's/^RS_//')

# 读取配置文件与样本清单 (优先从 samples.yaml 读取目标样本的独立 domain、presets 与文件路径)
CONFIG_YAML="${SCRIPT_DIR}/../projects/config_synteny/synteny_validation.yaml"
SAMPLES_YAML="${SCRIPT_DIR}/../projects/config_synteny/samples.yaml"

SAMPLE_PARSED=$($PYTHON_BIN -c "
import yaml, os

target = '$CLEAN_TARGET'
raw_target = '$TARGET'
samples_yaml = '$SAMPLES_YAML'
config_yaml = '$CONFIG_YAML'

s_dom, s_presets = '', ''
s_fa, s_faa, s_gff, s_tax, s_anno = '', '', '', '', ''
s_comp, s_contam = '', ''
def_dom, def_presets = '', ''

if os.path.exists(samples_yaml):
    try:
        with open(samples_yaml, 'r', encoding='utf-8') as f:
            y = yaml.safe_load(f) or {}
        def_dom = str(y.get('default_domain', '')).strip()
        dp = y.get('default_presets', [])
        if isinstance(dp, list):
            def_presets = ','.join(str(x).strip() for x in dp if str(x).strip())
        elif dp:
            def_presets = str(dp).strip()
        samples = y.get('synteny_samples', {})
        cand = samples.get(target) or samples.get(raw_target) or samples.get(raw_target.replace('GB_', '').replace('RS_', ''))
        if not cand:
            for k, v in samples.items():
                if k == target or k == raw_target or k.replace('GB_', '').replace('RS_', '') == target:
                    cand = v
                    break
        if isinstance(cand, dict):
            s_dom = str(cand.get('domain', '')).strip()
            pr = cand.get('presets') or cand.get('preset') or cand.get('active_presets') or cand.get('active_preset')
            if isinstance(pr, list):
                s_presets = ','.join(str(x).strip() for x in pr if str(x).strip())
            elif pr:
                s_presets = str(pr).strip()
            s_fa = str(cand.get('fasta') or cand.get('path') or '').strip()
            s_faa = str(cand.get('faa') or '').strip()
            s_gff = str(cand.get('gff') or cand.get('gff3') or '').strip()
            s_tax = str(cand.get('taxonomy') or '').strip()
            s_anno = str(cand.get('annotation') or cand.get('kofam_tsv') or '').strip()
            if 'completeness' in cand:
                s_comp = str(cand['completeness']).strip()
            if 'contamination' in cand:
                s_contam = str(cand['contamination']).strip()
    except Exception:
        pass

cfg_dom, cfg_presets = '', ''
if os.path.exists(config_yaml):
    try:
        with open(config_yaml, 'r', encoding='utf-8') as f:
            cy = yaml.safe_load(f).get('synteny_validation', {})
        cfg_dom = str(cy.get('domain', '')).strip()
        cp = cy.get('active_presets') or cy.get('active_preset')
        if isinstance(cp, list):
            cfg_presets = ','.join(str(x).strip() for x in cp if str(x).strip())
        elif cp:
            cfg_presets = str(cp).strip()
    except Exception:
        pass

print(f'{s_dom}\t{s_presets}\t{s_fa}\t{s_faa}\t{s_gff}\t{s_tax}\t{s_anno}\t{s_comp}\t{s_contam}\t{def_dom}\t{def_presets}\t{cfg_dom}\t{cfg_presets}')
" 2>/dev/null || true)

SAMPLE_DOMAIN=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $1}')
SAMPLE_PRESETS=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $2}')
SAMPLE_FA=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $3}')
SAMPLE_FAA=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $4}')
SAMPLE_GFF=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $5}')
SAMPLE_TAX=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $6}')
SAMPLE_ANNO=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $7}')
SAMPLE_COMP=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $8}')
SAMPLE_CONTAM=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $9}')
DEFAULT_DOMAIN=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $10}')
DEFAULT_PRESETS=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $11}')
CFG_DOMAIN=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $12}')
CFG_PRESETS=$(echo "$SAMPLE_PARSED" | awk -F'\t' '{print $13}')

ARG2="${POSITIONAL_ARGS[1]:-all}"
CUSTOM_TAXONOMY=""
MODULE_REQ="preset_default"

# 智能识别第2个参数是模块、预设快捷别名还是分类学字符串
if [[ "$ARG2" == d__* ]] || [[ "$ARG2" == *";"* ]]; then
    CUSTOM_TAXONOMY="$ARG2"
    MODULE_REQ="preset_default"
elif [ "$ARG2" = "sox" ] || [ "$ARG2" = "sox_system" ]; then
    MODULE_REQ="M00595"
    PRESET_ARG="sox_system"
elif [ "$ARG2" = "atp" ] || [ "$ARG2" = "atp_synthase" ]; then
    MODULE_REQ="atp"
    PRESET_ARG="atp_synthase"
elif [[ "$ARG2" =~ ^M[0-9]{5}$ ]]; then
    MODULE_REQ="$ARG2"
    PRESET_ARG="$ARG2"
elif [ "$ARG2" = "all" ] || [ "$ARG2" = "ALL" ]; then
    MODULE_REQ="preset_default"
else
    MODULE_REQ="$ARG2"
    CUSTOM_TAXONOMY="${POSITIONAL_ARGS[2]:-}"
fi

# 1. 判定运行生物域 (优先级: 命令行参数 > 样本配置 domain > 分类学推导 > 全局默认 domain > 细菌)
if [ -n "$DOMAIN_ARG" ]; then
    DOMAIN="$DOMAIN_ARG"
elif [ -n "$SAMPLE_DOMAIN" ]; then
    DOMAIN="$SAMPLE_DOMAIN"
elif [[ "$CUSTOM_TAXONOMY" == *"d__Archaea"* ]] || [[ "$ARG2" == *"d__Archaea"* ]] || [[ "$SAMPLE_TAX" == *"d__Archaea"* ]]; then
    DOMAIN="archaea"
elif [[ "$CUSTOM_TAXONOMY" == *"d__Bacteria"* ]] || [[ "$ARG2" == *"d__Bacteria"* ]] || [[ "$SAMPLE_TAX" == *"d__Bacteria"* ]]; then
    DOMAIN="bacteria"
elif [ -n "$DEFAULT_DOMAIN" ]; then
    DOMAIN="$DEFAULT_DOMAIN"
elif [ -n "$CFG_DOMAIN" ]; then
    DOMAIN="$CFG_DOMAIN"
else
    DOMAIN="bacteria"
fi

# 2. 判定运行预设 (优先级: 命令行参数 > 样本配置 presets > 全局默认 presets > 配置文件 active_presets > domain默认)
if [ -n "$PRESET_ARG" ]; then
    PRESET="$PRESET_ARG"
elif [ -n "$SAMPLE_PRESETS" ]; then
    PRESET="$SAMPLE_PRESETS"
elif [ -n "$DEFAULT_PRESETS" ]; then
    PRESET="$DEFAULT_PRESETS"
elif [ -n "$CFG_PRESETS" ]; then
    PRESET="$CFG_PRESETS"
elif [ "$DOMAIN" = "archaea" ]; then
    PRESET="M00159"
else
    PRESET="sox_system"
fi

echo "========================================================================"
echo ">>> 开始单样本共线性验证: $CLEAN_TARGET"
echo ">>> 运行生物域 (Domain): $DOMAIN"
echo ">>> 激活分析预设 (Preset): $PRESET"
echo ">>> 请求分析模块: $MODULE_REQ"
echo ">>> Clinker 可视化开关: $RUN_CLINKER"
echo ">>> 强制重跑所有步骤 (--force): $FORCE"
echo ">>> Contig 共有性收敛阈值: $CONSENSUS_THRESHOLD"
echo "========================================================================"

# 1. 基础路径定义与配置文件同步解析 (支持 domain_presets 独立分发)
BASE_DIR="/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin"

CFG_PATHS=$($PYTHON_BIN -c "
import yaml, os
cfg = '$CONFIG_YAML'
domain = '$DOMAIN'
base = '${BASE_DIR}'

if domain == 'archaea':
    g_meta = f'{base}/Database/GTDB/archive/metadata/232.0/ar53_metadata_r232.tsv'
    g_gen = f'{base}/results/data_prep/archaea_all/genomes'
    g_prot = f'{base}/results/mag_annotation/gtdb_r232_archaea_all/proteins/all'
    g_ko = f'{base}/results/mag_annotation/gtdb_r232_archaea_all/merged/protein_KO_all_strict.tsv'
    t_gen = ''
    t_prot = ''
    t_gene = ''
    t_class = ''
    q_faa = ''
    only_rep = 'false'
    max_pos = 5
    max_neg = 2
else:
    g_meta = f'{base}/Database/GTDB/archive/metadata/232.0/bac120_metadata_r232.tsv'
    g_gen = f'{base}/results/data_prep/bacteria_all/genomes'
    g_prot = f'{base}/results/data_prep/bacteria_all/proteins'
    g_ko = ''
    t_gen = f'{base}/Mpf_2022/MOB_abundance/data/MAGs'
    t_prot = f'{base}/results/mag_annotation/NS_MOB/proteins/all'
    t_gene = f'{base}/results/mag_annotation/NS_MOB/genes'
    t_class = f'{base}/results/mag_annotation/NS_MOB/mag_classification'
    q_faa = f'{base}/results/mag_annotation/NS_MOB/Reinekea_forsetii_sox_operon.faa'
    only_rep = 'true'
    max_pos = 10
    max_neg = 1

if os.path.exists(cfg):
    try:
        with open(cfg) as f:
            y = yaml.safe_load(f).get('synteny_validation', {})
            d_presets = y.get('domain_presets', {})
            d_cfg = d_presets.get(domain, {})
            g_meta = d_cfg.get('gtdb_metadata', y.get('gtdb_metadata', g_meta))
            g_gen = d_cfg.get('gtdb_genomes_dir', y.get('gtdb_genomes_dir', g_gen))
            g_prot = d_cfg.get('gtdb_proteins_dir', y.get('gtdb_proteins_dir', g_prot))
            g_ko = d_cfg.get('gtdb_ko_table', y.get('gtdb_ko_table', g_ko))
            t_gen = d_cfg.get('target_genomes_dir', y.get('target_genomes_dir', t_gen))
            t_prot = d_cfg.get('target_proteins_dir', y.get('target_proteins_dir', t_prot))
            t_gene = d_cfg.get('target_genes_dir', y.get('target_genes_dir', t_gene))
            t_class = d_cfg.get('target_classification', y.get('target_classification', t_class))
            q_faa = d_cfg.get('query_operon_faa', y.get('query_operon_faa', q_faa))
            only_rep = str(d_cfg.get('only_representative', y.get('only_representative', only_rep))).lower()
            max_pos = d_cfg.get('max_positive_refs', y.get('max_positive_refs', max_pos))
            max_neg = d_cfg.get('max_negative_refs', y.get('max_negative_refs', max_neg))
    except Exception:
        pass
print(f'{g_meta}\t{g_gen}\t{g_prot}\t{g_ko}\t{t_gen}\t{t_prot}\t{t_gene}\t{t_class}\t{q_faa}\t{only_rep}\t{max_pos}\t{max_neg}')
" 2>/dev/null || true)

GTDB_META=$(echo "$CFG_PATHS" | awk -F'\t' '{print $1}')
GENOMES_DIR=$(echo "$CFG_PATHS" | awk -F'\t' '{print $2}')
PROTEINS_DIR=$(echo "$CFG_PATHS" | awk -F'\t' '{print $3}')
GTDB_KO_TABLE=$(echo "$CFG_PATHS" | awk -F'\t' '{print $4}')
TARGET_GENOMES_DIR=$(echo "$CFG_PATHS" | awk -F'\t' '{print $5}')
TARGET_PROTEINS_DIR=$(echo "$CFG_PATHS" | awk -F'\t' '{print $6}')
TARGET_GENES_DIR=$(echo "$CFG_PATHS" | awk -F'\t' '{print $7}')
TARGET_CLASS_FILE=$(echo "$CFG_PATHS" | awk -F'\t' '{print $8}')
QUERY_OPERON_FAA=$(echo "$CFG_PATHS" | awk -F'\t' '{print $9}')
ONLY_REPRESENTATIVE=$(echo "$CFG_PATHS" | awk -F'\t' '{print $10}')
MAX_POS_REFS=$(echo "$CFG_PATHS" | awk -F'\t' '{print $11}')
MAX_NEG_REFS=$(echo "$CFG_PATHS" | awk -F'\t' '{print $12}')

# 2. 解析目标文件 (优先使用 samples.yaml 声明的路径，次选在待测 MAG 数据集中检索，最后在 GTDB 库中检索)
TARGET_FA=""
if [ -n "$SAMPLE_FA" ] && [ -f "$SAMPLE_FA" ]; then
    TARGET_FA="$SAMPLE_FA"
fi
if [ -z "$TARGET_FA" ] && [ -n "$TARGET_GENOMES_DIR" ] && [ -d "$TARGET_GENOMES_DIR" ]; then
    TARGET_FA=$(ls -1 "${TARGET_GENOMES_DIR}"/*"${CLEAN_TARGET}"*.fna* 2>/dev/null | head -n 1 || ls -1 "${TARGET_GENOMES_DIR}"/*"${CLEAN_TARGET}"*.fa* 2>/dev/null | head -n 1 || true)
fi
if [ -z "$TARGET_FA" ] && [ -n "$GENOMES_DIR" ] && [ -d "$GENOMES_DIR" ]; then
    TARGET_FA=$(ls -1 "${GENOMES_DIR}"/*"${CLEAN_TARGET}"*.fna* 2>/dev/null | head -n 1 || ls -1 "${GENOMES_DIR}"/*"${CLEAN_TARGET}"*.fa* 2>/dev/null | head -n 1 || true)
fi

TARGET_FAA=""
if [ -n "$SAMPLE_FAA" ] && [ -f "$SAMPLE_FAA" ]; then
    TARGET_FAA="$SAMPLE_FAA"
fi
if [ -z "$TARGET_FAA" ] && [ -n "$TARGET_PROTEINS_DIR" ] && [ -d "$TARGET_PROTEINS_DIR" ]; then
    TARGET_FAA=$(ls -1 "${TARGET_PROTEINS_DIR}"/*"${CLEAN_TARGET}"*.faa* 2>/dev/null | head -n 1 || true)
fi
if [ -z "$TARGET_FAA" ] && [ -n "$PROTEINS_DIR" ] && [ -d "$PROTEINS_DIR" ]; then
    TARGET_FAA=$(ls -1 "${PROTEINS_DIR}"/*"${CLEAN_TARGET}"*.faa* 2>/dev/null | head -n 1 || true)
fi

TARGET_GFF=""
if [ -n "$SAMPLE_GFF" ] && [ -f "$SAMPLE_GFF" ]; then
    TARGET_GFF="$SAMPLE_GFF"
fi
if [ -z "$TARGET_GFF" ] && [ -n "$TARGET_GENES_DIR" ] && [ -d "$TARGET_GENES_DIR" ]; then
    TARGET_GFF=$(ls -1 "${TARGET_GENES_DIR}"/*"${CLEAN_TARGET}"*.gff* 2>/dev/null | head -n 1 || true)
fi

if [ -z "$TARGET_FA" ] || [ ! -f "$TARGET_FA" ]; then
    echo "[-] 错误: 未找到目标核酸基因组文件: ${CLEAN_TARGET} (搜索目录: ${TARGET_GENOMES_DIR}, ${GENOMES_DIR})"
    exit 1
fi
echo "[+] 目标核酸: $TARGET_FA"
echo "[+] 目标蛋白: ${TARGET_FAA:-'未提供 (将由 Prodigal 现场直译)'}"
echo "[+] 目标 GFF: ${TARGET_GFF:-'未提供'}"

# 3. 提取 Taxonomy 与 CheckM 指标
# 优先使用 samples.yaml 声明的 taxonomy
if [ -z "$CUSTOM_TAXONOMY" ] && [ -n "$SAMPLE_TAX" ]; then
    CUSTOM_TAXONOMY="$SAMPLE_TAX"
    echo "[+] 从 samples.yaml 提取分类成功: $CUSTOM_TAXONOMY"
fi
# 其次从 target_classification (mag_classification) 提取分类
if [ -z "$CUSTOM_TAXONOMY" ] && [ -n "$TARGET_CLASS_FILE" ] && [ -f "$TARGET_CLASS_FILE" ]; then
    FOUND_TAX=$(awk -F'\t' -v t="$CLEAN_TARGET" '($1 == t || $1 ~ t) {print $2; exit}' "$TARGET_CLASS_FILE" || true)
    if [ -n "$FOUND_TAX" ]; then
        CUSTOM_TAXONOMY="$FOUND_TAX"
        echo "[+] 从待测分类表提取成功: $CUSTOM_TAXONOMY"
    fi
fi
META_INFO=$($PYTHON_BIN -c "
import sys, gzip, os

meta_file = '$GTDB_META'
target = '$CLEAN_TARGET'
tax = '$CUSTOM_TAXONOMY'
comp, contam = '100.0', '0.0'

if os.path.exists(meta_file):
    opener = gzip.open if meta_file.endswith('.gz') else open
    try:
        with opener(meta_file, 'rt', encoding='utf-8', errors='ignore') as f:
            header_line = f.readline()
            if header_line:
                header = [c.strip().lower() for c in header_line.strip().split('\t')]
                acc_idx = header.index('accession') if 'accession' in header else 0
                tax_idx = header.index('gtdb_taxonomy') if 'gtdb_taxonomy' in header else (header.index('taxonomy') if 'taxonomy' in header else -1)
                comp_idx = -1
                for col in ['checkm_completeness', 'checkm2_completeness', 'completeness']:
                    if col in header:
                        comp_idx = header.index(col)
                        break
                contam_idx = -1
                for col in ['checkm_contamination', 'checkm2_contamination', 'contamination']:
                    if col in header:
                        contam_idx = header.index(col)
                        break

                for line in f:
                    if target in line:
                        parts = line.strip().split('\t')
                        if len(parts) > acc_idx:
                            acc_val = parts[acc_idx].replace('GB_', '').replace('RS_', '').strip()
                            if acc_val == target or parts[acc_idx].strip() == target:
                                if not tax and tax_idx >= 0 and tax_idx < len(parts):
                                    tax = parts[tax_idx].strip()
                                if comp_idx >= 0 and comp_idx < len(parts):
                                    comp = parts[comp_idx].strip()
                                if contam_idx >= 0 and contam_idx < len(parts):
                                    contam = parts[contam_idx].strip()
                                break
    except Exception as e:
        print(f'Warning reading metadata: {e}', file=sys.stderr)

if not tax:
    tax = 'd__Archaea;'
print(f'{tax}\t{comp}\t{contam}')
")

TAXONOMY=$(echo "$META_INFO" | awk -F'\t' '{print $1}')
COMPLETENESS=$(echo "$META_INFO" | awk -F'\t' '{print $2}')
CONTAMINATION=$(echo "$META_INFO" | awk -F'\t' '{print $3}')

COMPLETENESS=${SAMPLE_COMP:-${COMPLETENESS:-100.0}}
CONTAMINATION=${SAMPLE_CONTAM:-${CONTAMINATION:-0.0}}

echo "[+] 目标分类: $TAXONOMY"
echo "[+] 完整度: ${COMPLETENESS}% | 污染度: ${CONTAMINATION}%"

# 4. 确定待运行模块列表
MODULES_TO_RUN=()
if [ "$MODULE_REQ" = "preset_default" ] || [ "$MODULE_REQ" = "all" ] || [ "$MODULE_REQ" = "ALL" ]; then
    if [ "$PRESET" = "atp_synthase" ] || [ "$PRESET" = "atp" ]; then
        MODULES_TO_RUN=("M00159" "M00157" "M00158" "M00160")
    elif [ "$PRESET" = "sox_system" ] || [ "$PRESET" = "sox" ]; then
        MODULES_TO_RUN=("M00595")
    else
        # 尝试从 YAML 的 operon_presets 中动态提取所有激活预设的模块列表 (支持多预设如 M00159,M00157)
        YAML_MODS=$($PYTHON_BIN -c "
import yaml, sys
try:
    with open('$CONFIG_YAML', 'r', encoding='utf-8') as f:
        y = yaml.safe_load(f).get('synteny_validation', {})
    presets = y.get('operon_presets', {})
    raw_p = [p.strip() for p in '$PRESET'.split(',') if p.strip()]
    all_mods = []
    for p in raw_p:
        if p in ['atp_synthase', 'atp']:
            all_mods.extend(['M00159', 'M00157', 'M00158', 'M00160'])
        elif p in ['sox_system', 'sox']:
            all_mods.append('M00595')
        elif p in presets:
            p_data = presets[p]
            mods = p_data.get('enabled_modules', [])
            if not mods and 'modules' in p_data:
                mods = list(p_data['modules'].keys())
            all_mods.extend(mods)
        elif p.startswith('M') and len(p) == 6:
            all_mods.append(p)
    seen = set()
    res = [m for m in all_mods if not (m in seen or seen.add(m))]
    if res:
        print(' '.join(res))
except Exception:
    pass
" 2>/dev/null || true)
        if [ -n "$YAML_MODS" ]; then
            read -r -a MODULES_TO_RUN <<< "$YAML_MODS"
        else
            MODULES_TO_RUN=("$PRESET")
        fi
    fi
elif [ "$MODULE_REQ" = "atp" ] || [ "$MODULE_REQ" = "atp_synthase" ]; then
    MODULES_TO_RUN=("M00159" "M00157" "M00158" "M00160")
elif [ "$MODULE_REQ" = "sox" ] || [ "$MODULE_REQ" = "sox_system" ]; then
    MODULES_TO_RUN=("M00595")
else
    MODULES_TO_RUN=("$MODULE_REQ")
fi

GFF_ARG=""
if [ -n "$TARGET_GFF" ] && [ -f "$TARGET_GFF" ]; then
    GFF_ARG="--target_gff $TARGET_GFF"
fi

# 5. 循环执行各模块独立多轨分析
ANY_STAGE_RERUN=false

for MOD in "${MODULES_TO_RUN[@]}"; do
    case "$MOD" in
        "M00595")
            MOD_NAME="Sulfur oxidation, SOX system (soxBAXZYCD), thiosulfate => sulfate"
            TARGET_KOS="K17224,K17222,K17223,K17227,K17226,K17225,K22622"
            MOD_EXP_LEN=10000
            ;;
        "M00159")
            MOD_NAME="V/A-type ATPase, prokaryotes"
            TARGET_KOS="K02117,K02118,K02119,K02120,K02121,K02122,K02107,K02123,K02124"
            MOD_EXP_LEN=15000
            ;;
        "M00157")
            MOD_NAME="F-type ATPase, prokaryotes and chloroplasts"
            TARGET_KOS="K02111,K02112,K02113,K02114,K02115,K02108,K02109,K02110"
            MOD_EXP_LEN=12000
            ;;
        "M00158")
            MOD_NAME="F-type ATPase, eukaryotes"
            TARGET_KOS="K02126,K02127,K02128,K02129,K02130,K02131,K02132,K02133,K02134,K02135,K02136,K02137,K02138"
            MOD_EXP_LEN=20000
            ;;
        "M00160")
            MOD_NAME="V-type ATPase, eukaryotes"
            TARGET_KOS="K02145,K02146,K02147,K02148,K02149,K02150,K02151,K02152,K02153,K02154,K02155"
            MOD_EXP_LEN=18000
            ;;
        *)
            echo "[-] 提示: 使用自定义模块: $MOD，尝试从预设或知识库动态解析参数..."
            DYN_INFO=$($PYTHON_BIN -c "
import yaml, os, sys
m = '$MOD'
kos = ''
name = m
exp_len = 15000
if os.path.exists('$CONFIG_YAML'):
    try:
        with open('$CONFIG_YAML') as f:
            y = yaml.safe_load(f).get('synteny_validation', {})
            for p in y.get('operon_presets', {}).values():
                mods = p.get('modules', {})
                if m in mods:
                    minf = mods[m]
                    name = minf.get('name', name)
                    exp_len = minf.get('expected_locus_length', exp_len)
                    tk = minf.get('target_kos', [])
                    kos = ','.join(tk) if isinstance(tk, list) else str(tk)
                    break
    except Exception:
        pass
if not kos:
    k_path = os.path.join('$SCRIPT_DIR', '..', 'knowledge', 'kegg', 'modules', f'{m}.yaml')
    if os.path.exists(k_path):
        try:
            with open(k_path) as f:
                ky = yaml.safe_load(f)
                name = ky.get('name', name)
                tk = ky.get('key_kos', [])
                kos = ','.join(tk) if isinstance(tk, list) else str(tk)
        except Exception:
            pass
print(f'{name}\t{kos}\t{exp_len}')
" 2>/dev/null || true)
            if [ -n "$DYN_INFO" ]; then
                MOD_NAME=$(echo "$DYN_INFO" | awk -F'\t' '{print $1}')
                DYN_KOS=$(echo "$DYN_INFO" | awk -F'\t' '{print $2}')
                DYN_LEN=$(echo "$DYN_INFO" | awk -F'\t' '{print $3}')
                if [ -n "$DYN_KOS" ]; then
                    TARGET_KOS="$DYN_KOS"
                else
                    TARGET_KOS="K02117,K02118,K02119,K02120,K02121,K02122,K02107,K02123,K02124"
                fi
                MOD_EXP_LEN=${DYN_LEN:-15000}
            else
                MOD_NAME="$MOD"
                TARGET_KOS="K02117,K02118,K02119,K02120,K02121,K02122,K02107,K02123,K02124"
                MOD_EXP_LEN=15000
            fi
            ;;
    esac

    MOD_OUT_DIR="${BASE_DIR}/results/synteny_validation/${CLEAN_TARGET}/${MOD}"
    mkdir -p "${MOD_OUT_DIR}/stage0" "${MOD_OUT_DIR}/stage1" "${MOD_OUT_DIR}/stage2" "${MOD_OUT_DIR}/stage3" "${MOD_OUT_DIR}/stage4" "${MOD_OUT_DIR}/visualization"
    KO_CACHE="${BASE_DIR}/results/synteny_validation/cache/${MOD}_gtdb_index.tsv"
    mkdir -p "$(dirname "$KO_CACHE")"

    echo -e "\n========================================================================"
    echo ">>> [模块轨道: ${MOD}] 检查全流程分析: ${MOD_NAME}"
    echo ">>> 输出目录: ${MOD_OUT_DIR}"
    echo "========================================================================"

    STAGE_RERUN=false

    # Stage 0: 参考基因组筛选
    STAGE0_OUT="${MOD_OUT_DIR}/stage0/selected_validated_references.tsv"
    STAGE0_ANCHORS="${MOD_OUT_DIR}/stage0/consensus_anchors.tsv"
    STAGE0_SEED="${MOD_OUT_DIR}/stage0/${CLEAN_TARGET}_seed_query_operon.faa"

    if [ "$FORCE" = "true" ] || [ ! -s "$STAGE0_OUT" ] || [ ! -s "$STAGE0_ANCHORS" ] || [ ! -s "$STAGE0_SEED" ]; then
        echo -e "\n------------------------------------------------------------------------"
        echo ">>> [${MOD} | Stage 0] 筛选阳性与阴性参考基因组..."
        echo "------------------------------------------------------------------------"

        # 检查参考基因组库是否需要按需下载 (仅在细菌模式下触发)
        mkdir -p "$GENOMES_DIR" "$PROTEINS_DIR"
        REF_COUNT=$(ls -1 "$GENOMES_DIR"/*.fna* 2>/dev/null | wc -l || true)
        if [ "$REF_COUNT" -eq 0 ] && [ "$DOMAIN" = "bacteria" ]; then
            echo -e "\n[!] 提示: 本地近缘细菌参考基因组库为空 (${GENOMES_DIR})。"
            echo ">>> 正在启动按需下载程序 (fetch_bacterial_references.py) 获取近缘阳性与阴性参考..."
            python3 "${SCRIPT_DIR}/fetch_bacterial_references.py" \
                --gtdb_metadata "$GTDB_META" \
                --target_taxonomy "${CUSTOM_TAXONOMY:-}" \
                --target_classification "${TARGET_CLASS_FILE:-}" \
                --query_operon_faa "${QUERY_OPERON_FAA:-}" \
                --out_genomes_dir "$GENOMES_DIR" \
                --out_proteins_dir "$PROTEINS_DIR" \
                --max_positive_refs "$MAX_POS_REFS" \
                --max_negative_refs "$MAX_NEG_REFS" || true
        fi

        QUERY_ARG=""
        if [ -n "$QUERY_OPERON_FAA" ] && [ -f "$QUERY_OPERON_FAA" ]; then
            QUERY_ARG="--query_operon_faa $QUERY_OPERON_FAA"
        fi

        KO_ARG=""
        if [ -n "$GTDB_KO_TABLE" ] && [ -f "$GTDB_KO_TABLE" ]; then
            KO_ARG="--gtdb_ko_table $GTDB_KO_TABLE"
        fi

        AUTO_FETCH_FLAG=""
        if [ "$DOMAIN" = "bacteria" ]; then
            AUTO_FETCH_FLAG="--auto_fetch"
        else
            AUTO_FETCH_FLAG="--no_auto_fetch"
        fi

        ONLY_REP_FLAG=""
        if [ "$ONLY_REPRESENTATIVE" = "true" ]; then
            ONLY_REP_FLAG="--only_representative"
        fi

        python3 "${SCRIPT_DIR}/select_and_validate_refs.py" \
            --domain "$DOMAIN" \
            $AUTO_FETCH_FLAG \
            $ONLY_REP_FLAG \
            --target_name "$CLEAN_TARGET" \
            --target_fa "$TARGET_FA" \
            --target_faa "${TARGET_FAA:-}" \
            --target_taxonomy "$TAXONOMY" \
            --target_completeness "$COMPLETENESS" \
            --target_contamination "$CONTAMINATION" \
            --module_id "$MOD" \
            --module_name "$MOD_NAME" \
            --gtdb_metadata "$GTDB_META" \
            --gtdb_genomes_dir "$GENOMES_DIR" \
            --gtdb_proteins_dir "$PROTEINS_DIR" \
            --target_kos "$TARGET_KOS" \
            $KO_ARG \
            $QUERY_ARG \
            --ko_index_cache "$KO_CACHE" \
            --min_af 0.50 \
            --adaptive_af \
            --min_ani 80.0 \
            --max_positive_refs "$MAX_POS_REFS" \
            --max_negative_refs "$MAX_NEG_REFS" \
            --threads 4 \
            --out_refs_tsv "$STAGE0_OUT" \
            --out_anchors_tsv "$STAGE0_ANCHORS" \
            --out_seed_operon_faa "$STAGE0_SEED" \
            --out_dir "${MOD_OUT_DIR}/stage0"
        STAGE_RERUN=true
        ANY_STAGE_RERUN=true
    else
        echo -e "\n[-] [${MOD} | Stage 0] 产物已存在，跳过参考筛选 (使用 --force 可重新运行)"
    fi

    # Stage 1: 保守侧翼锚点定位与局部基因座切片
    STAGE1_TOPO="${MOD_OUT_DIR}/stage1/anchor_topology.tsv"
    STAGE1_FA="${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.fa"
    STAGE1_GBK="${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.gbk"
    STAGE1_PROTS="${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_100kb_locus_proteins.tsv"
    STAGE1_PROTS_LINES=0
    if [ -f "$STAGE1_PROTS" ]; then
        STAGE1_PROTS_LINES=$(wc -l < "$STAGE1_PROTS" | tr -d ' ')
    fi

    if [ "$FORCE" = "true" ] || [ "$STAGE_RERUN" = "true" ] || \
       [ ! -s "$STAGE1_TOPO" ] || [ ! -s "$STAGE1_FA" ] || [ ! -s "$STAGE1_GBK" ] || \
       [ ! -s "$STAGE1_PROTS" ] || [ "$STAGE1_PROTS_LINES" -le 1 ]; then
        echo -e "\n------------------------------------------------------------------------"
        echo ">>> [${MOD} | Stage 1] 保守侧翼锚点定位与局部切片提取..."
        echo "------------------------------------------------------------------------"
        ANNO_ARG=""
        if [ -n "$SAMPLE_ANNO" ] && [ -f "$SAMPLE_ANNO" ]; then
            ANNO_ARG="--target_annotation_tsv $SAMPLE_ANNO"
        fi

        python3 "${SCRIPT_DIR}/extract_anchor_loci.py" \
            --target_name "$CLEAN_TARGET" \
            --target_fa "$TARGET_FA" \
            --target_faa "${TARGET_FAA:-}" \
            --module_id "$MOD" \
            --module_name "$MOD_NAME" \
            --target_kos "$TARGET_KOS" \
            $GFF_ARG \
            $ANNO_ARG \
            --gtdb_ko_table "$GTDB_KO_TABLE" \
            --ko_index_cache "$KO_CACHE" \
            --threads 4 \
            --anchors_tsv "$STAGE0_ANCHORS" \
            --refs_tsv "$STAGE0_OUT" \
            --seed_operon_faa "$STAGE0_SEED" \
            --expected_locus_length "$MOD_EXP_LEN" \
            --min_edge_bp 10000 \
            --flank_bp 20000 \
            --n_flank_genes 25 \
            --max_visualized_refs 10 \
            --consensus_contig_threshold "$CONSENSUS_THRESHOLD" \
            --significance_mode "$SIGNIFICANCE_MODE" \
            --significance_alpha "$SIGNIFICANCE_ALPHA" \
            --n_permutations "$N_PERMUTATIONS" \
            $([ "$ALLOW_INSUFFICIENT" = "true" ] && echo "--allow_insufficient_evidence_locus") \
            --extended_window_bp 100000 \
            --out_topology_tsv "$STAGE1_TOPO" \
            --out_target_locus_fa "$STAGE1_FA" \
            --out_target_locus_gbk "$STAGE1_GBK" \
            --out_refs_locus_dir "${MOD_OUT_DIR}/stage1/refs_loci" \
            --out_extended_proteins_tsv "$STAGE1_PROTS"
        STAGE_RERUN=true
        ANY_STAGE_RERUN=true
    else
        echo -e "\n[-] [${MOD} | Stage 1] 产物已存在且完整，跳过切片提取 (使用 --force 可重新运行)"
    fi

    # Stage 2: 序列三级挽救与假基因化判定
    STAGE2_GENE="${MOD_OUT_DIR}/stage2/gene_level_rescue_summary.tsv"
    STAGE2_CLUST="${MOD_OUT_DIR}/stage2/cluster_state.tsv"

    if [ "$FORCE" = "true" ] || [ "$STAGE_RERUN" = "true" ] || [ ! -s "$STAGE2_GENE" ] || [ ! -s "$STAGE2_CLUST" ]; then
        echo -e "\n------------------------------------------------------------------------"
        echo ">>> [${MOD} | Stage 2] 序列三级挽救与假基因化判定..."
        echo "------------------------------------------------------------------------"
        STAGE1_PROTS_ARG=""
        if [ -f "$STAGE1_PROTS" ] && [ -s "$STAGE1_PROTS" ]; then
            STAGE1_PROTS_ARG="--target_annotation_tsv $STAGE1_PROTS"
        fi

        python3 "${SCRIPT_DIR}/rescue_and_classify_gene_state.py" \
            --target_name "$CLEAN_TARGET" \
            --target_locus_fa "${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.fa" \
            --topology_tsv "${MOD_OUT_DIR}/stage1/anchor_topology.tsv" \
            --target_faa "${TARGET_FAA:-}" \
            --query_operon_faa "${MOD_OUT_DIR}/stage0/${CLEAN_TARGET}_seed_query_operon.faa" \
            --module_id "$MOD" \
            --module_name "$MOD_NAME" \
            $GFF_ARG \
            $STAGE1_PROTS_ARG \
            --blastp_identity 25.0 \
            --blastp_coverage 30.0 \
            --blastx_evalue 1e-3 \
            --out_gene_tsv "$STAGE2_GENE" \
            --out_cluster_tsv "$STAGE2_CLUST"
        STAGE_RERUN=true
        ANY_STAGE_RERUN=true
    else
        echo -e "\n[-] [${MOD} | Stage 2] 产物已存在，跳过三级挽救 (使用 --force 可重新运行)"
    fi

    # Stage 3: 物理对齐与断裂点验证
    STAGE3_OUT="${MOD_OUT_DIR}/stage3/locus_absence_validation.tsv"
    STAGE3_PAF="${MOD_OUT_DIR}/stage3/${CLEAN_TARGET}_vs_ref.paf"

    if [ "$FORCE" = "true" ] || [ "$STAGE_RERUN" = "true" ] || [ ! -s "$STAGE3_OUT" ]; then
        echo -e "\n------------------------------------------------------------------------"
        echo ">>> [${MOD} | Stage 3] 物理对齐与断裂点验证..."
        echo "------------------------------------------------------------------------"
        REF_FA=$(ls -1 "${MOD_OUT_DIR}/stage1/refs_loci"/*.fa 2>/dev/null | head -n 1 || echo "${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.fa")

        python3 "${SCRIPT_DIR}/validate_syntenic_absence.py" \
            --target_name "$CLEAN_TARGET" \
            --target_locus_fa "${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.fa" \
            --ref_locus_fa "$REF_FA" \
            --topology_tsv "${MOD_OUT_DIR}/stage1/anchor_topology.tsv" \
            --cluster_tsv "${MOD_OUT_DIR}/stage2/cluster_state.tsv" \
            --module_id "$MOD" \
            --module_name "$MOD_NAME" \
            --out_paf "$STAGE3_PAF" \
            --out_locus_tsv "$STAGE3_OUT"
        STAGE_RERUN=true
        ANY_STAGE_RERUN=true
    else
        echo -e "\n[-] [${MOD} | Stage 3] 产物已存在，跳过断裂点对齐 (使用 --force 可重新运行)"
    fi

    # Stage 4: Fitch Parsimony 演化推断与全量汇总表
    STAGE4_OUT="${MOD_OUT_DIR}/stage4/${CLEAN_TARGET}_locus_validation_summary.tsv"
    STAGE4_TREE="${MOD_OUT_DIR}/stage4/${CLEAN_TARGET}_local_tree.nwk"

    if [ "$FORCE" = "true" ] || [ "$STAGE_RERUN" = "true" ] || [ ! -s "$STAGE4_OUT" ]; then
        echo -e "\n------------------------------------------------------------------------"
        echo ">>> [${MOD} | Stage 4] 演化推断与报告生成..."
        echo "------------------------------------------------------------------------"
        python3 "${SCRIPT_DIR}/infer_evolutionary_loss.py" \
            --target_name "$CLEAN_TARGET" \
            --target_taxonomy "$TAXONOMY" \
            --target_completeness "$COMPLETENESS" \
            --target_contamination "$CONTAMINATION" \
            --module_id "$MOD" \
            --module_name "$MOD_NAME" \
            --refs_tsv "${MOD_OUT_DIR}/stage0/selected_validated_references.tsv" \
            --topology_tsv "${MOD_OUT_DIR}/stage1/anchor_topology.tsv" \
            --cluster_tsv "${MOD_OUT_DIR}/stage2/cluster_state.tsv" \
            --locus_tsv "${MOD_OUT_DIR}/stage3/locus_absence_validation.tsv" \
            --out_summary_tsv "$STAGE4_OUT" \
            --out_tree_nwk "$STAGE4_TREE"
        STAGE_RERUN=true
        ANY_STAGE_RERUN=true
    else
        echo -e "\n[-] [${MOD} | Stage 4] 产物已存在，跳过演化推断 (使用 --force 可重新运行)"
    fi

    # Clinker 可视化绘图与 Label 优化
    if [ "$RUN_CLINKER" = "true" ] || [ "$RUN_CLINKER" = "TRUE" ] || [ "$RUN_CLINKER" = "1" ]; then
        MULTI_HTML="${MOD_OUT_DIR}/visualization/${CLEAN_TARGET}_clinker.html"
        if [ "$FORCE" = "true" ] || [ "$STAGE_RERUN" = "true" ] || [ ! -s "$MULTI_HTML" ]; then
            echo -e "\n------------------------------------------------------------------------"
            echo ">>> [${MOD} | Clinker] 可视化多轨图与两两对比图渲染 (RUN_CLINKER=true)..."
            echo "------------------------------------------------------------------------"
            ORDER_FILE="${MOD_OUT_DIR}/stage1/refs_loci/ordered_clinker_files.txt"

            # 多轨图
            if [ -f "$ORDER_FILE" ] && [ -s "$ORDER_FILE" ]; then
                REF_ARGS=$(cat "$ORDER_FILE" | tr '\n' ' ')
                conda run -n clinker_0.0.32 clinker $REF_ARGS "${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.gbk" -p "$MULTI_HTML" -i 0.25 || true
            else
                conda run -n clinker_0.0.32 clinker "${MOD_OUT_DIR}/stage1/refs_loci"/*.gbk "${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.gbk" -p "$MULTI_HTML" -i 0.25 || true
            fi
            if [ -f "$MULTI_HTML" ]; then
                python3 "${SCRIPT_DIR}/refine_clinker_labels.py" "$MULTI_HTML" || true
            fi

            # 两两对比图 (按 max 设置遍历所有筛选到的 Ref)
            REF_LIST=""
            if [ -f "$ORDER_FILE" ] && [ -s "$ORDER_FILE" ]; then
                REF_LIST=$(cat "$ORDER_FILE")
            else
                REF_LIST=$(ls -1 "${MOD_OUT_DIR}/stage1/refs_loci"/*_locus.gbk 2>/dev/null || true)
            fi

            for REF_GBK in $REF_LIST; do
                if [ -f "$REF_GBK" ]; then
                    REF_BASE=$(basename "$REF_GBK" _locus.gbk)
                    PAIR_HTML="${MOD_OUT_DIR}/visualization/${CLEAN_TARGET}_vs_${REF_BASE}_clinker.html"
                    REF_TARGET_GBK="${MOD_OUT_DIR}/stage1/refs_loci/${CLEAN_TARGET}_vs_${REF_BASE}_target_locus.gbk"
                    if [ ! -f "$REF_TARGET_GBK" ]; then
                        REF_TARGET_GBK="${MOD_OUT_DIR}/stage1/${CLEAN_TARGET}_locus.gbk"
                    fi
                    echo "[+] 绘制两两对比: $PAIR_HTML"
                    conda run -n clinker_0.0.32 clinker "$REF_GBK" "$REF_TARGET_GBK" -p "$PAIR_HTML" -i 0.25 || true
                    if [ -f "$PAIR_HTML" ]; then
                        python3 "${SCRIPT_DIR}/refine_clinker_labels.py" "$PAIR_HTML" || true
                    fi
                fi
            done
        else
            echo -e "\n[-] [${MOD} | Clinker] HTML 对比图已存在，跳过渲染 (使用 --force 可重新渲染)"
        fi
    else
        echo -e "\n------------------------------------------------------------------------"
        echo ">>> [${MOD} | Clinker] 跳过 Clinker 绘图 (RUN_CLINKER=false, 默认关闭)"
        echo ">>> 提示: 如需生成 Clinker 对比图，可在命令行中追加 --clinker 或设置 RUN_CLINKER=true"
        echo "------------------------------------------------------------------------"
    fi

    echo "[+] 模块 ${MOD} 完成: ${MOD_OUT_DIR}/stage4/${CLEAN_TARGET}_locus_validation_summary.tsv"
done

# 6. 汇总结算全模块汇总表与四大 ATP 综合透视矩阵
MATRIX_TSV="${BASE_DIR}/results/synteny_validation/${CLEAN_TARGET}/atp_synthase_synteny_matrix.tsv"
STACKED_TSV="${BASE_DIR}/results/synteny_validation/${CLEAN_TARGET}/locus_synteny_all_modules_summary.tsv"

if [ "$FORCE" = "true" ] || [ "$ANY_STAGE_RERUN" = "true" ] || [ ! -s "$MATRIX_TSV" ] || [ ! -s "$STACKED_TSV" ]; then
    echo -e "\n========================================================================"
    echo ">>> [Summary] 编译综合全模块汇总表与四模块透视矩阵..."
    echo "========================================================================"

    STAGE4_FILES=$(ls -1 "${BASE_DIR}/results/synteny_validation/${CLEAN_TARGET}"/*/stage4/${CLEAN_TARGET}_locus_validation_summary.tsv 2>/dev/null || true)

    if [ -n "$STAGE4_FILES" ]; then
        python3 "${SCRIPT_DIR}/compile_atp_synteny_matrix.py" \
            --input_tsvs $STAGE4_FILES \
            --out_stacked_tsv "$STACKED_TSV" \
            --out_matrix_tsv "$MATRIX_TSV"
        
        echo "[+] 全模块堆叠总表: $STACKED_TSV"
        echo "[+] 四模块透视矩阵: $MATRIX_TSV"
        echo -e "\n>>> 矩阵快速预览 (atp_synthase_synteny_matrix.tsv):"
        head -n 5 "$MATRIX_TSV"
    fi
else
    echo -e "\n[-] [Summary] 综合全模块矩阵已存在且无上游模块更新，跳过重新编译 (使用 --force 可重新编译)"
fi

echo -e "\n========================================================================"
echo ">>> 目标样本 [${CLEAN_TARGET}] 分析流程全部完成！"
echo "========================================================================"
