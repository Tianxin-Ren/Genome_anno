# Genome Annotation Module
# ------------------------------------------------------------------------------

import os
import glob

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
AN_ROOT = config.get("project_root", ".")
AN_ENV_DIR = config.get("conda_env_dir", "")

# Dynamically retrieve annotation config to ensure registry-loaded values are present
if "annotation" not in config:
    config["annotation"] = {}
AN_C = config["annotation"]

# Fallback to advanced_analysis settings or root config if annotation config is absent/minimal
# We merge keys that are not present in AN_C to avoid overwriting custom values (like unified_hmm_db)
if not AN_C.get("bakta_db"):
    fallback_src = config.get("advanced_analysis", config)
    if isinstance(fallback_src, dict):
        for k, v in fallback_src.items():
            if k not in AN_C and k != "annotation":
                AN_C[k] = v

if "steps" not in AN_C:
    AN_C["steps"] = {"bakta": "True", "cdhit": "False"}




def an_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(AN_ROOT, rel_path))

DEFAULT_AN_ENVS = {
    "gtdbtk": "gtdbtk_2.7.1",
    "checkm2": "checkm2_1.0.2",
    "checkm": "checkm_1.2.4",
    "bakta": "bakta_1.12.0",
    "prokka": "prokka_1.14.6",
    "biopython": "biopython_1.85",
    "cdhit": "cd-hit_4.8.1",
    "quast": "quast_5.3.0",
    "anvio": "anvio_9",
    "longread_audit": "longread_audit",
    "coverm": "coverm_0.7.0",
    "mapping": "longread_audit",
    "minimap2": "longread_audit",
    "samtools": "longread_audit",
}

def an_get_env_dir():
    dir_candidate = config.get("conda_env_dir", "") or config.get("paths", {}).get("conda_envs", "")
    if dir_candidate and os.path.isabs(dir_candidate):
        return dir_candidate
    for fallback in [
        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs",
        os.path.expanduser("~/.conda/envs"),
    ]:
        if os.path.isdir(fallback):
            return fallback
    return dir_candidate or os.path.expanduser("~/.conda/envs")

def an_get_env(name):
    default_name = DEFAULT_AN_ENVS.get(name, name)
    env_name = AN_C.get("conda_envs", {}).get(name) or config.get("environments", {}).get(name) or default_name
    if os.path.isabs(env_name):
        return env_name
    env_dir = an_get_env_dir()
    if env_dir:
        return os.path.join(env_dir, env_name)
    return env_name

def an_get_locus_tag(genome, bin_or_chunk):
    # Bakta limit is 24 chars. 
    # 优先检查 bin_or_chunk 本身是否为符合规范的合法 Locus Tag (例如 u999856)
    # Locus Tag 规范：英文字母开头，只含英文字母、数字和下划线，至少 3 个字符。
    import re
    clean_bin = re.sub(r'[^A-Za-z0-9_]', '', bin_or_chunk)
    
    # 如果 clean_bin 已经是字母开头，且长度 >= 3 字符，则直接使用它作为 locus_tag！
    # 这将完美解决类似 "u999856" 带有 "Bathy_enrich_u999856" 前缀的问题，使其直降为极简的 "u999856"
    if len(clean_bin) >= 3 and clean_bin[0].isalpha():
        return clean_bin[:24]
        
    # 如果 clean_bin 不是字母开头（例如是纯数字或数字开头的简写 "1" 等），则加上极简前缀 "c" (形成 c1, c10 这种 cxx 风格)，绝不包含 Bathyxxx
    if clean_bin:
        return f"c{clean_bin}"[:24]

        
    # 终极兜底方案：
    clean_genome = genome.replace("/", "_")
    tag = f"l_{clean_genome[:8]}_{bin_or_chunk}"
    return tag[:24]


AN_RES = an_full_path(AN_C.get("results_dir", "results/annotation"))

AN_LOGS = an_full_path(AN_C.get("logs_dir", os.path.join("logs", "annotation")))
AN_IN_DIR = an_full_path(AN_C.get("input_dir", "results/assembly"))
AN_SAMPLES = config.get("samples", {})
if AN_SAMPLES is None: AN_SAMPLES = {}
AN_BN_RES = an_full_path(config.get("binning", {}).get("results_dir", "results/binning"))
AN_AS_RES = an_full_path(config.get("assembly", {}).get("results_dir", "results/assembly"))

EXT_GENOMES = config.get("external_genomes", {})
if EXT_GENOMES is None:
    EXT_GENOMES = {}
EXT_RES = an_full_path(config.get("annotation", {}).get("external_results_dir", "results/external_annotation"))

# 2. Scatter-Gather Configuration
# ------------------------------------------------------------------------------
AN_SPLIT_COUNT = AN_C.get("split_chunks", 10)
if AN_SPLIT_COUNT < 1: AN_SPLIT_COUNT = 1

scattergather:
    an_chunks = AN_SPLIT_COUNT,
    hmm_chunks = AN_C.get("hmm_split_chunks", 20)

# 3. Target Registration
# ------------------------------------------------------------------------------
AN_GENOMES = [
    g for g in AN_SAMPLES.keys()
    if "external_assembly" in AN_SAMPLES[g] 
    or AN_SAMPLES[g].get("type") == "protein"
    or config.get("modules", {}).get("assembly", {}).get("enable", False)
    or config.get("assembly", {}).get("enable", False)
]
if not AN_GENOMES:
    # Fallback to standard flat directory scanning if no active samples configured.
    # We only scan direct child files ending in .fa to prevent matching subdirectories.
    if os.path.exists(AN_IN_DIR):
        AN_GENOMES = [
            os.path.splitext(f)[0]
            for f in os.listdir(AN_IN_DIR)
            if f.endswith(".fa") and os.path.isfile(os.path.join(AN_IN_DIR, f))
        ]
    else:
        AN_GENOMES = []

# --- A. Normal Assembly Annotation ---
for g in AN_GENOMES:
    steps = config.get("annotation", {}).get("steps", {})
    do_bakta = str(steps.get("bakta", "False")).lower() == "true"
    do_eggnog = str(steps.get("eggnog", "False")).lower() == "true"
    do_kofam = str(steps.get("kofamscan", "False")).lower() == "true"
    do_adv_hmm = str(steps.get("advanced_hmm", "False")).lower() == "true"

    # Support protein type (bypassing Bakta/Cd-hit/Circular assembly rules)
    sample_info = AN_SAMPLES.get(g, {})
    is_protein = isinstance(sample_info, dict) and sample_info.get("type") == "protein"
    if is_protein:
        if do_eggnog:
            TARGETS.append(os.path.join(AN_RES, g, f"eggnog/{g}.emapper.annotations"))
        if do_kofam:
            TARGETS.append(os.path.join(AN_RES, g, f"kofamscan/{g}_kofam.tsv"))
        if do_adv_hmm:
            TARGETS.append(os.path.join(AN_RES, g, "advanced_hmm/unified_all_raw.tbl"))
            do_pfam = str(steps.get("pfam", "False")).lower() == "true"
            if do_pfam:
                TARGETS.append(os.path.join(AN_RES, g, "advanced_hmm/pfam_all_raw.domtbl"))
            do_resolve = str(steps.get("resolve_hmm", "True")).lower() == "true"
            if do_resolve:
                TARGETS.append(os.path.join(AN_RES, g, f"{g}_advanced_hmm_resolved.tsv"))
        continue

    # 1. Bakta
    if do_bakta:
        TARGETS.append(os.path.join(AN_RES, g, "bakta", f"{g}.gff3"))
        do_cdhit = str(steps.get("cdhit", "False")).lower() == "true"
        if do_cdhit:
            c_val = AN_C.get("cdhit_c", 0.95)
            TARGETS.append(os.path.join(AN_RES, g, "cdhit", f"{g}_proteomics_database_c{c_val}.faa"))
    
    # 2. EggNOG & KofamScan (Assembly-level: only if assembly-level Bakta is enabled)
    if do_bakta:
        if do_eggnog:
            TARGETS.append(os.path.join(AN_RES, g, f"eggnog/{g}.emapper.annotations"))
        if do_kofam:
            TARGETS.append(os.path.join(AN_RES, g, f"kofamscan/{g}_kofam.tsv"))
    
    # 3. Refined Annotation (Auto-disabled for standalone external assemblies)
    is_external_asm = "external_assembly" in AN_SAMPLES.get(g, {})
    do_circular = str(steps.get("circular", "False")).lower() == "true"
    if do_circular and not is_external_asm:
        TARGETS.append(os.path.join(AN_RES, g, "circular_summary.tsv"))
    
    # 4. Advanced HMM (Intermediate, integrated into Summary)
    do_hmm = str(steps.get("advanced_hmm", "False")).lower() == "true"
    do_resolve = str(steps.get("resolve_hmm", "True")).lower() == "true"
    if do_hmm and not do_resolve:
        TARGETS.append(os.path.join(AN_RES, g, "advanced_hmm/unified_all_raw.tbl"))
        do_pfam = str(steps.get("pfam", "False")).lower() == "true"
        if do_pfam:
            TARGETS.append(os.path.join(AN_RES, g, "advanced_hmm/pfam_all_raw.domtbl"))

    # 5. Functional Summary
    if do_resolve:
        TARGETS.append(os.path.join(AN_RES, g, f"{g}_functional_summary.tsv"))


# --- B. External Genomes (Logic Reuse) ---
if EXT_GENOMES:
    _batch_groups = sorted(set(
        info.get("batch_group") for info in EXT_GENOMES.values()
        if isinstance(info, dict) and "batch_group" in info
    ))
    _do_resolve_hmm = str(config.get("annotation", {}).get("steps", {}).get("resolve_hmm", "True")).lower() == "true"
    if _do_resolve_hmm:
        # Instead of a global summary, merge each batch group/sample separately
        for bg in _batch_groups:
            TARGETS.append(os.path.join(EXT_RES, f"{bg}_functional_summary.tsv"))

    # Register batch sentinels to TARGETS
    for bg in _batch_groups:
        TARGETS.append(os.path.join(EXT_RES, f"{bg}_batch_hmm_completed.txt"))
        dedup_cfg = config.get("annotation", {}).get("protein_dedup", {}) or {}
        if str(dedup_cfg.get("enabled", False)).lower() == "true":
            TARGETS.append(os.path.join(EXT_RES, f"batch_{bg}_protein_dedup/unique_proteins.faa"))
            TARGETS.append(os.path.join(EXT_RES, f"batch_{bg}_protein_dedup/protein_map.tsv"))
            TARGETS.append(os.path.join(EXT_RES, f"batch_{bg}_protein_dedup/catalog.sqlite"))
            if str(dedup_cfg.get("cluster_99", True)).lower() == "true":
                TARGETS.append(os.path.join(EXT_RES, f"batch_{bg}_protein_dedup/cluster_99/representatives.faa"))

for name, info in EXT_GENOMES.items():
    if isinstance(info, dict) and info.get("type") == "dna_dir":
        continue
    # If this genome belongs to a batch group, its targets are represented by the batch sentinel
    if isinstance(info, dict) and "batch_group" in info:
        continue
    steps = config.get("annotation", {}).get("steps", {})
    do_bakta = str(steps.get("bakta", "False")).lower() == "true"
    do_eggnog = str(steps.get("eggnog", "False")).lower() == "true"
    do_kofam = str(steps.get("kofamscan", "False")).lower() == "true"

    # Determine path prefix: batch genomes nest under {batch_group}/
    batch_group = info.get("batch_group") if isinstance(info, dict) else None
    if batch_group:
        res_prefix = os.path.join(EXT_RES, batch_group, name)
    else:
        res_prefix = os.path.join(EXT_RES, name)

    # Support protein type (bypassing Bakta/Cd-hit/Aggregation rules)
    is_protein = isinstance(info, dict) and info.get("type") == "protein"
    if is_protein:
        if do_eggnog:
            TARGETS.append(os.path.join(res_prefix, f"eggnog/{name}.emapper.annotations"))
        if do_kofam:
            TARGETS.append(os.path.join(res_prefix, f"kofamscan/{name}_kofam.tsv"))
        do_adv_hmm = str(steps.get("advanced_hmm", "False")).lower() == "true"
        if do_adv_hmm:
            TARGETS.append(os.path.join(res_prefix, "advanced_hmm/unified_all_raw.tbl"))
            do_pfam = str(steps.get("pfam", "False")).lower() == "true"
            if do_pfam:
                TARGETS.append(os.path.join(res_prefix, "advanced_hmm/pfam_all_raw.domtbl"))
            do_resolve = str(steps.get("resolve_hmm", "True")).lower() == "true"
            if do_resolve:
                TARGETS.append(os.path.join(res_prefix, f"{name}_advanced_hmm_resolved.tsv"))
        continue

    # 1. Bakta (Reference)
    if do_bakta:
        TARGETS.append(os.path.join(res_prefix, f"{name}.gff3"))
        do_cdhit = str(steps.get("cdhit", "False")).lower() == "true"
        if do_cdhit:
            c_val = AN_C.get("cdhit_c", 0.95)
            TARGETS.append(os.path.join(res_prefix, f"cdhit/{name}_proteomics_database_c{c_val}.faa"))
    
    # 2. Functional (Reference)
    if do_eggnog:
        TARGETS.append(os.path.join(res_prefix, f"eggnog/{name}.emapper.annotations"))
    if do_kofam:
        TARGETS.append(os.path.join(res_prefix, f"kofamscan/{name}_kofam.tsv"))
    
    # 3. Advanced HMM (Intermediate, integrated into Summary)
    do_hmm = str(steps.get("advanced_hmm", "False")).lower() == "true"
    do_resolve = str(steps.get("resolve_hmm", "True")).lower() == "true"
    if do_hmm and not do_resolve:
        TARGETS.append(os.path.join(res_prefix, "advanced_hmm/unified_all_raw.tbl"))
        do_pfam = str(steps.get("pfam", "False")).lower() == "true"
        if do_pfam:
            TARGETS.append(os.path.join(res_prefix, "advanced_hmm/pfam_all_raw.domtbl"))
    
    # 4. Aggregate (Reference)
    if do_resolve:
        TARGETS.append(os.path.join(res_prefix, f"{name}_functional_summary.tsv"))






# 4. Helpers
# ------------------------------------------------------------------------------
def get_sample_name_from_genome(genome):
    samples = config.get("samples", {})
    if not samples:
        return genome.split("/")[0]
    
    # Sort sample names by length descending
    for s in sorted(samples.keys(), key=len, reverse=True):
        if genome.startswith(s + "/") or genome == s:
            return s
        if genome.startswith(s + "_"):
            return s
            
    return genome.split("/")[0]

def an_get_asm_path(genome):
    if "get_filtered_assembly" in globals():
        return get_filtered_assembly(genome)
    else:
        return os.path.join(AN_IN_DIR, f"{genome}.fa")

def get_anno_faa_input(wildcards):
    # Determine sample name from wildcards (genome or bin)
    name = wildcards.genome if hasattr(wildcards, "genome") else (wildcards.bin if hasattr(wildcards, "bin") else None)
    if not name:
        return ""
    
    # 1. Check in AN_SAMPLES
    sample_info = AN_SAMPLES.get(name, {})
    if isinstance(sample_info, dict) and sample_info.get("type") == "protein":
        return an_full_path(sample_info.get("path"))
        
    # 2. Check in EXT_GENOMES
    ext_info = EXT_GENOMES.get(name, {})
    if isinstance(ext_info, dict) and ext_info.get("type") == "protein":
        return an_full_path(ext_info.get("path"))
        
    # 3. Standard Bakta FAA output
    is_external = name in EXT_GENOMES
    if is_external:
        return os.path.join(EXT_RES, f"{name}/{name}.faa")
    else:
        return os.path.join(AN_RES, f"{name}/bakta/{name}.faa")

def get_selected_contigs_info(sample):
    """Returns a dict of {contig_id: {complete: bool}} from config."""
    sample_info = AN_SAMPLES.get(sample, {})
    # Support both old 'selected_circulars' (as true) and new 'selected_contigs'
    contigs = sample_info.get("selected_contigs", {})
    if not contigs:
        old_circulars = sample_info.get("selected_circulars", [])
        contigs = {c: {"complete": True} for c in old_circulars}
    return contigs

def get_complete_contigs(sample):
    """Returns a list of contig IDs explicitly marked as complete in config."""
    info = get_selected_contigs_info(sample)
    return [c for c, cfg in info.items() if cfg.get("complete") is True]

def get_non_complete_contigs(sample):
    """Returns a list of contig IDs explicitly marked as NOT complete in config."""
    info = get_selected_contigs_info(sample)
    return [c for c, cfg in info.items() if cfg.get("complete") is False]

def get_cdhit_n_helper(c_str):
    """Helper to dynamically calculate word size (n) for cdhit based on similarity threshold (c)."""
    try:
        c = float(c_str)
    except ValueError:
        return 5
    if c >= 0.90:
        return 5
    elif c >= 0.85:
        return 4
    elif c >= 0.80:
        return 3
    else:
        return 2

# 5. Rule Modules Inclusion
# ------------------------------------------------------------------------------
include: "rules/bakta.smk"
include: "rules/functional.smk"
include: "rules/circular.smk"
include: "rules/circular_illumina.smk"
include: "rules/circular_evidence.smk"
include: "rules/external.smk"
