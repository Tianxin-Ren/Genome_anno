# ==========================================
# 宏基因组功能基因定量流程 - 全局配置
# ==========================================

# 1. 核心过滤阈值
MIN_SCORE = 100
MIN_BSR = 0.9

# 2. 目标样本列表
# 如果留空，统计脚本将自动扫描输出目录中的所有样本
TARGET_SAMPLES = [
    'SK1712_6961','SK1712_6964','SK1712_6967','SK1712_6970','SK1715_6498','SK1715_6506','SK1715_6512','SK1715_6520','SK1712_6962','SK1712_6965','SK1712_6968','SK1712_6971','SK1715_6500','SK1715_6508','SK1715_6517','SK1715_6522','SK1712_6963','SK1712_6966','SK1712_6969','SK1712_6972','SK1715_6502','SK1715_6510','SK1715_6518','SK1715_6959'
]

# 3. 目标功能基因列表
TARGET_GENES = [
    'amoA', 'bNOR', 'cNOR', 'dur3', 'eNOR', 'gNOR', 'hdh', 'hzsA', 
    'nNOR', 'napAB', 'napA', 'narG', 'nirK', 'nirS', 'nosZ', 'qNOR', 
    'sNOR', 'ureC', 'nrfA', 'rpoB'
]

# 4. 基因长度字典 (氨基酸长度，计算时会自动乘以 3)
GENE_LENGTHS = {
    'amoA': 250, 'bNOR': 560, 'cNOR': 460, 'dur3': 660, 'eNOR': 564, 
    'gNOR': 544, 'hdh': 558, 'hzsA': 809, 'nNOR': 504, 'napAB': 834, 
    'napA': 776, 'narG': 1223, 'nirK': 413, 'nirS': 583, 'nosZ': 670, 
    'qNOR': 762, 'sNOR': 489, 'ureC': 572, 'nrfA': 474, 'rpoB': 1172
}

# 5. 路径配置
METAGENOMICS_DIR = '/home/data-1/rtx/Mariana/results/metagenomics'
OUT_DIR = './read_function'
LOOKUP_PL = './blast_based_read_lookup_new.pl'
LOG_FILE = './pipeline.log'

# 6. 数据库配置
DB_CONFIG = {
    'nif': {
        'faa': '../../nitrogen/combined/combined_db_r207.faa',
        'dmnd': '../../nitrogen/outgroup/combined_db_r207_outgroup.dmnd',
        'suffix': '__combined_db_nif_r207__result_bsr.tsv'
    },
    'rpob': {
        'faa': '../../nitrogen/combined/combined_nrfA-rpoB.faa',
        'dmnd': '../../nitrogen/outgroup/combined_nrfA-rpoB_outgroup.dmnd',
        'suffix': '__rpoB_final__result_bsr.tsv'
    }
}
RPOB_FASTA = DB_CONFIG['rpob']['faa']

# 7. 并行执行配置
CONCURRENT_JOBS = 4  # 同时并行的样本数
THREADS_PER_JOB = 64 # 每个比对任务使用的线程数
