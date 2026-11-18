import copy
import glob
import os
import re
import sys

# ------------------------------------------------------------------------------
# 0. Global Targets (Must be defined first for rule all)
# ------------------------------------------------------------------------------
TARGETS = []


rule all:
    input:
        lambda wildcards: TARGETS

# ------------------------------------------------------------------------------
# 1. Master Configuration & Registry Loading
# ------------------------------------------------------------------------------
import yaml
from snakemake.utils import update_config

# Preserve command-line values so they keep the highest precedence after all
# project files have been loaded.
CLI_CONFIG = copy.deepcopy(dict(config))
WORKFLOW_ROOT = os.path.abspath(workflow.current_basedir or os.getcwd())
sys.path.insert(0, WORKFLOW_ROOT)

from common.scripts.resources import calc_threads, calc_mem_mb
from common.scripts.runtime import (
    apply_runtime_config,
    deep_merge,
    expand_path,
    unresolved_variables,
)

# Can be overridden via: --config conf=projects/<name> or cnf=projects/<name>
CONF_DIR = os.path.normpath(str(expand_path(CLI_CONFIG.get("conf") or CLI_CONFIG.get("cnf") or "config")))
if not os.path.isabs(CONF_DIR):
    CONF_DIR = os.path.join(WORKFLOW_ROOT, CONF_DIR)

def get_threads(requested):
    return calc_threads(requested, config, workflow)

def get_mem_mb(threads):
    return calc_mem_mb(threads, workflow)

def get_slurm_account():
    # 1. Check loaded config hierarchy (runtime -> cluster -> default_resources -> root)
    for src in (config.get("runtime", {}), config.get("cluster", {}), config.get("default_resources", {}), config):
        if isinstance(src, dict):
            for k in ("slurm_account", "account"):
                val = src.get(k)
                if val:
                    return str(val).strip()
    # 2. Check environment variable
    env_acc = os.environ.get("SLURM_ACCOUNT") or os.environ.get("SNAKEMAKE_SLURM_ACCOUNT")
    if env_acc:
        return env_acc.strip()
    # 3. Environment heuristic based on current cluster user or directory
    user = os.environ.get("USER", "").lower()
    cwd = os.getcwd().lower()
    if any(k in user or k in cwd for k in ("clswfp", "icdli")):
        return "acct-clswfp"
    if any(k in user or k in cwd for k in ("schhui", "gses0")):
        return "acct-schhui"
    return ""

def get_slurm_partition(default_val="64c512g"):
    for src in (config.get("runtime", {}), config.get("cluster", {}), config.get("default_resources", {}), config):
        if isinstance(src, dict):
            for k in ("slurm_partition", "partition"):
                val = src.get(k)
                if val:
                    return str(val).strip()
    env_part = os.environ.get("SLURM_PARTITION")
    if env_part:
        return env_part.strip()
    return default_val


# 1.1 Configuration helpers
def load_yaml(path):
    with open(path) as handle:
        return yaml.safe_load(handle) or {}


def remap_config_paths(data, conf_dir):
    """Resolve project-local config references without touching data paths."""
    if isinstance(data, dict):
        for key, value in data.items():
            data[key] = remap_config_paths(value, conf_dir)
    elif isinstance(data, list):
        for index, value in enumerate(data):
            data[index] = remap_config_paths(value, conf_dir)
    elif isinstance(data, str):
        if data.startswith("./"):
            return os.path.normpath(os.path.join(conf_dir, data))
        if data.startswith("config/") and os.path.basename(conf_dir) != "config":
            project_candidate = os.path.normpath(os.path.join(conf_dir, data[len("config/"):]))
            # Legacy module references such as config/annotation.yaml should
            # follow a custom project directory only when that file exists.
            # Shared references such as config/platforms/slurm.yaml remain
            # relative to the workflow root.
            if os.path.exists(project_candidate):
                return project_candidate
        if data == "config" and os.path.basename(conf_dir) != "config":
            return conf_dir
    return data


def resolve_platform_path(value):
    value = str(expand_path(value))
    if os.path.isabs(value) and os.path.exists(value):
        return value

    candidates = [
        value,
        os.path.join(CONF_DIR, value),
        os.path.join(WORKFLOW_ROOT, value),
        os.path.join(WORKFLOW_ROOT, "config", "platforms", value),
    ]
    if not value.endswith((".yaml", ".yml")):
        candidates.extend([
            f"{value}.yaml",
            os.path.join(CONF_DIR, f"{value}.yaml"),
            os.path.join(WORKFLOW_ROOT, f"{value}.yaml"),
            os.path.join(WORKFLOW_ROOT, "config", "platforms", f"{value}.yaml"),
        ])
    for c in candidates:
        if os.path.exists(c):
            return c

    # Fuzzy matching for platform aliases (e.g., icdli <-> icdli1)
    base_name = os.path.basename(value)
    stem = re.sub(r'\.ya?ml$', '', base_name)
    alt_stems = [stem[:-1]] if stem.endswith("1") else [stem + "1"]
    for alt in alt_stems:
        for folder in (os.path.join(WORKFLOW_ROOT, "config", "platforms"), CONF_DIR, WORKFLOW_ROOT):
            for ext in (".yaml", ".yml"):
                alt_path = os.path.join(folder, alt + ext)
                if os.path.exists(alt_path):
                    return alt_path

    return os.path.join(WORKFLOW_ROOT, value)


# 1.2 Load project + optional platform configuration.
PIPELINE_YAML = os.path.join(CONF_DIR, "pipeline.yaml")
if not os.path.exists(PIPELINE_YAML):
    raise FileNotFoundError(
        f"ERROR: Configuration directory '{CONF_DIR}' or its master configuration '{PIPELINE_YAML}' "
        "could not be found. Pass --config conf=projects/<name>."
    )
PROJECT_CONFIG = remap_config_paths(load_yaml(PIPELINE_YAML), CONF_DIR)
PLATFORM_REF = CLI_CONFIG.get("platform")
if not PLATFORM_REF:
    # 1. Environment heuristic: auto-detect cluster platform from user or working directory
    user = os.environ.get("USER", "").lower()
    cwd = os.getcwd().lower()
    if any(k in user or k in cwd for k in ("schhui", "gses0")):
        gses0_cand = os.path.join(WORKFLOW_ROOT, "config", "platforms", "gses0.yaml")
        if os.path.exists(gses0_cand):
            PLATFORM_REF = gses0_cand
    elif any(k in user or k in cwd for k in ("clswfp", "icdli")):
        icdli_cand = os.path.join(WORKFLOW_ROOT, "config", "platforms", "icdli1.yaml")
        if os.path.exists(icdli_cand):
            PLATFORM_REF = icdli_cand

if not PLATFORM_REF:
    PLATFORM_REF = PROJECT_CONFIG.get("runtime", {}).get("platform")

if not PLATFORM_REF:
    # 2. Auto-detect platform from --workflow-profile or --profile CLI arguments
    for idx, arg in enumerate(sys.argv):
        if arg in ("--workflow-profile", "--profile") and idx + 1 < len(sys.argv):
            prof_name = sys.argv[idx + 1].rstrip("/").split("/")[-1]
            candidate = os.path.join(WORKFLOW_ROOT, "config", "platforms", f"{prof_name}.yaml")
            if os.path.exists(candidate):
                PLATFORM_REF = candidate
                break
    # 3. Universal fallback: default to slurm on HPC or local platform
    if not PLATFORM_REF:
        for fallback_candidate in ("config/platforms/slurm.yaml", "config/platforms/local.yaml"):
            candidate = os.path.join(WORKFLOW_ROOT, fallback_candidate)
            if os.path.exists(candidate):
                PLATFORM_REF = candidate
                break

PLATFORM_CONFIG = {}
if PLATFORM_REF:
    PLATFORM_YAML = resolve_platform_path(PLATFORM_REF)
    if not os.path.exists(PLATFORM_YAML):
        raise FileNotFoundError(f"ERROR: Platform configuration not found: {PLATFORM_YAML}")
    PLATFORM_CONFIG = load_yaml(PLATFORM_YAML)

config.clear()
update_config(config, deep_merge(PLATFORM_CONFIG, PROJECT_CONFIG, CLI_CONFIG))

# 1.3 Load sample configuration.
SAMPLES_PATH = config.get("data", {}).get("samples", os.path.join(CONF_DIR, "samples.yaml"))
if os.path.exists(SAMPLES_PATH):
    sample_cfg = load_yaml(SAMPLES_PATH)
    update_config(config, sample_cfg)
    if sample_cfg.get("project_root"):
        config.setdefault("paths", {})["project_root"] = sample_cfg["project_root"]

# 1.4 Load enabled module configuration.
MODULES = config.get("modules", {})
for mod_name, mod_info in MODULES.items():
    if mod_info.get("enable", False) and mod_info.get("path"):
        if os.path.exists(mod_info["path"]):
            update_config(config, load_yaml(mod_info["path"]))
        else:
            print(f"WARNING: Module {mod_name} enabled but config file missing: {mod_info['path']}", file=sys.stderr)

# Explicit CLI values win. The runtime adapter then maps the portable schema to
# legacy module keys (project_root, conda_env_dir, per-module database paths).
update_config(config, CLI_CONFIG)
batch_name = config.get("batch_name")
if batch_name is not None and (not isinstance(batch_name, str) or not batch_name or batch_name in (".", "..") or "/" in batch_name or "\\" in batch_name):
    raise ValueError("batch_name must be a non-empty directory name without path separators")
if CLI_CONFIG.get("project_root"):
    config.setdefault("paths", {})["project_root"] = CLI_CONFIG["project_root"]
if CLI_CONFIG.get("conda_env_dir"):
    config.setdefault("paths", {})["conda_envs"] = CLI_CONFIG["conda_env_dir"]
if CLI_CONFIG.get("results_dir"):
    config.setdefault("paths", {})["results"] = CLI_CONFIG["results_dir"]
if CLI_CONFIG.get("logs_dir"):
    config.setdefault("paths", {})["logs"] = CLI_CONFIG["logs_dir"]
if CLI_CONFIG.get("work_dir"):
    config.setdefault("paths", {})["work"] = CLI_CONFIG["work_dir"]
apply_runtime_config(config)

# Propagate platform cluster account/partition to workflow default_resources
_platform_acc = get_slurm_account()
_platform_part = get_slurm_partition("64c512g")

def _update_default_res_obj(target_dr, acc, part):
    if target_dr is None:
        return
    # 1. Update .args list (evaluated by Snakemake when constructing job resources)
    if hasattr(target_dr, "args") and isinstance(target_dr.args, list):
        if acc:
            target_dr.args = [a for a in target_dr.args if not a.startswith("slurm_account=")]
            target_dr.args.append(f'slurm_account="{acc}"')
        if part:
            target_dr.args = [a for a in target_dr.args if not a.startswith("slurm_partition=")]
            target_dr.args.append(f'slurm_partition="{part}"')
    # 2. Update .parsed dict if present
    if hasattr(target_dr, "parsed") and isinstance(target_dr.parsed, dict):
        if acc:
            target_dr.parsed["slurm_account"] = acc
        if part:
            target_dr.parsed["slurm_partition"] = part
    # 3. Direct mapping if supported
    if hasattr(target_dr, "__setitem__"):
        try:
            if acc:
                target_dr["slurm_account"] = acc
            if part:
                target_dr["slurm_partition"] = part
        except Exception:
            pass

if hasattr(workflow, "default_resources") and workflow.default_resources is not None:
    try:
        _update_default_res_obj(workflow.default_resources, _platform_acc, _platform_part)
    except Exception:
        pass

if hasattr(workflow, "resource_settings") and hasattr(workflow.resource_settings, "default_resources"):
    try:
        _update_default_res_obj(workflow.resource_settings.default_resources, _platform_acc, _platform_part)
    except Exception:
        pass

UNRESOLVED = unresolved_variables(config)
if UNRESOLVED:
    print(
        "WARNING: Unresolved environment variables in configuration: " + ", ".join(UNRESOLVED),
        file=sys.stderr,
    )


# 1.5 Dynamic Genome Auto-Discovery (HPC Reference/Isolate Genomes Folder Scanning)
genomes_dir = config.get("genomes_dir")
ext_genomes = config.setdefault("external_genomes", {})
if ext_genomes is None:
    ext_genomes = {}

# If genomes_dir is specified, we create a default batch entry for it in ext_genomes
if genomes_dir:
    if isinstance(genomes_dir, str):
        genomes_dirs = [genomes_dir]
    elif isinstance(genomes_dir, list):
        genomes_dirs = genomes_dir
    else:
        genomes_dirs = []
    
    for i, g_dir in enumerate(genomes_dirs):
        discovered_batch_name = f"default_batch_{i+1}" if len(genomes_dirs) > 1 else "default_batch"
        if discovered_batch_name not in ext_genomes:
            ext_genomes[discovered_batch_name] = {
                "path": g_dir,
                "type": "dna_dir",
                # Inherit global settings
                "hmm_db": config.get("annotation", {}).get("unified_hmm_db", ""),
                "ko_list": config.get("annotation", {}).get("ko_list", ""),
                "batch_size": config.get("annotation", {}).get("batch_size", 500)
            }

# Now we scan all "dna_dir" entries in ext_genomes
discovered_count = 0
discovered_genomes = {}
for discovered_batch_name, batch_info in list(ext_genomes.items()):
    if isinstance(batch_info, dict) and batch_info.get("type") == "dna_dir":
        batch_dir = batch_info.get("path")
        # Resolve path relative to project_root if relative and proj_root exists
        proj_root = config.get("project_root", "")
        if not os.path.isabs(batch_dir) and not batch_dir.startswith("./") and proj_root:
            abs_batch_dir = os.path.normpath(os.path.join(proj_root, batch_dir))
        else:
            abs_batch_dir = os.path.abspath(batch_dir)
            
        if os.path.exists(abs_batch_dir):
            for ext in ["*.fna", "*.fasta", "*.fa", "*.fna.gz", "*.fasta.gz", "*.fa.gz"]:
                for f_path in glob.glob(os.path.join(abs_batch_dir, ext)):
                    basename = os.path.basename(f_path)
                    # Skip protein files
                    if basename.endswith((".faa", ".faa.gz")):
                        continue
                    name = re.sub(r'\.(fna|fasta|fa)(\.gz)?$', '', basename)
                    if name not in ext_genomes:
                        discovered_genomes[name] = {
                            "path": f_path,
                            "complete": False,
                            "type": "dna",
                            "batch_group": discovered_batch_name
                        }
                        discovered_count += 1

ext_genomes.update(discovered_genomes)
config["external_genomes"] = ext_genomes

if discovered_count > 0:
    print(f"INFO: Dynamically discovered {discovered_count} external genomes via batch folders.", file=sys.stderr)






# ------------------------------------------------------------------------------
# 2. Global Setup & Helpers
# ------------------------------------------------------------------------------
MODULES_DIR = os.path.join(WORKFLOW_ROOT, "modules")


# Helper to resolve conda environments
def resolve_env(name_or_path):
    if not name_or_path: return ""
    if os.path.isabs(name_or_path): return name_or_path
    local_path = os.path.join(config.get("conda_env_dir", ""), name_or_path)
    if os.path.exists(local_path): return local_path
    return name_or_path

# ------------------------------------------------------------------------------
# 3. Module Inclusion (Driven by Master Config)
# ------------------------------------------------------------------------------

# ENGINE 1: Preprocess
if MODULES.get("preprocess", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "preprocess/preprocess.smk")

# ENGINE 2: Assembly
if MODULES.get("assembly", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "assembly/assembly.smk")

# ENGINE 3: Binning
if MODULES.get("binning", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "binning/binning.smk")

# ENGINE 4: Unified Annotation
if MODULES.get("annotation", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "annotation/annotation.smk")
    TARGETS.extend(AN_CIRCULAR_EVIDENCE_TARGETS)

# ENGINE 4a: Reusable annotation post-processing
if MODULES.get("postprocess", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "postprocess/postprocess.smk")

# ENGINE 4b: Unified Bacteria/Archaea BioB Annotation
if MODULES.get("unified_annotation", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "unified_annotation/unified_annotation.smk")

# ENGINE 5: Advanced Analysis
if MODULES.get("advanced_analysis", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "advanced_analysis/advanced_analysis.smk")
    TARGETS.extend(ADV_TARGETS)

# ENGINE 6: HMM Build
if MODULES.get("hmm_build", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "advanced_analysis/hmm_build.smk")
    TARGETS.extend(HMM_TARGETS)

# ENGINE 7: Plasmid Analysis
if MODULES.get("plasmid_analysis", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "plasmid_analysis/plasmid_analysis.smk")

# ENGINE 8: Abundance Analysis
if MODULES.get("abundance", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "abundance/abundance.smk")

# ENGINE 9: Data Preparation
if MODULES.get("data_prep", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "data_prep/data_prep.smk")

# ENGINE 10: Phylogeny Analysis
if MODULES.get("phylogeny", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "phylogeny/phylogeny.smk")

# ENGINE 11: Visualization Module
if MODULES.get("visualize", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "visualize/visualize.smk")

# ENGINE 12: Protein Phylogeny Module
if MODULES.get("protein_phylogeny", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "protein_phylogeny/protein_phylogeny.smk")

# ENGINE 13: MAG Ecological & Metabolic Traits Module (dbCAN, CAMPER, gRodon2)
if MODULES.get("mag_traits", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "mag_traits/mag_traits.smk")
    TARGETS.extend(MAG_TRAITS_TARGETS)

# ENGINE 14: MAG High-Throughput Functional Annotation (PyKofamSearch + Full KOfam)
if MODULES.get("mag_annotation", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "mag_annotation/mag_annotation.smk")
    TARGETS.extend(MAG_ANNO_TARGETS)

# ENGINE 15: KEGG Module Reconstruction & Bioenergetic Capability Engine
if MODULES.get("pathway_reconstruction", {}).get("enable") and str(config.get("stage", "all")).lower() != "gene_calling":
    include: os.path.join(MODULES_DIR, "pathway_reconstruction/pathway_reconstruction.smk")
    TARGETS.extend(PATHWAY_TARGETS)

# ENGINE 16: Synteny & Operon Gene-Loss Validation Engine
if MODULES.get("synteny_validation", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "synteny_validation/synteny_validation.smk")
    TARGETS.extend(SYNTENY_TARGETS)

# ENGINE 17: COG High-Throughput Footprint-Aware Annotation Engine (COG2024 / Cog_LE)
if MODULES.get("cog_annotation", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "cog_annotation/cog_annotation.smk")
    TARGETS.extend(COG_ANNO_TARGETS)

# ENGINE 18: Metagenomic Read Taxonomic Profiling Engine (Sylph)
if MODULES.get("read_profiling", {}).get("enable") or MODULES.get("sylph", {}).get("enable"):
    include: os.path.join(MODULES_DIR, "read_profiling/read_profiling.smk")
    TARGETS.extend(READ_PROFILING_TARGETS)

# ENGINE 19: Hydrogenase Functional Annotation Engine (34-Class Mechanistic Engine)
if MODULES.get("hydrogenase_annotation", {}).get("enable") and str(config.get("stage", "all")).lower() != "gene_calling":
    include: os.path.join(MODULES_DIR, "hydrogenase_annotation/hydrogenase_annotation.smk")
    TARGETS.extend(HYD_ANNO_TARGETS)





# ------------------------------------------------------------------------------
# 4. Global Resource Attachment (Ensures SLURM account/partition propagation)
# ------------------------------------------------------------------------------
_final_acc = get_slurm_account()
_final_part = get_slurm_partition("64c512g")
if hasattr(workflow, "rules"):
    for r in workflow.rules:
        if hasattr(r, "resources") and isinstance(r.resources, dict):
            if _final_acc and ("slurm_account" not in r.resources or not r.resources["slurm_account"]):
                r.resources["slurm_account"] = _final_acc
            if _final_part and ("slurm_partition" not in r.resources or not r.resources["slurm_partition"]):
                r.resources["slurm_partition"] = _final_part
