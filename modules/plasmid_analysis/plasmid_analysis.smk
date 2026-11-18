# Plasmid Analysis Module: geNomad and MOB-suite integration
import os

# 1. Global Path Logic
# ------------------------------------------------------------------------------
ROOT = config.get("project_root", ".")
ENV_DIR = config.get("conda_env_dir", "")
c = config.get("plasmid_analysis", {})

def full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(ROOT, rel_path))

def get_env(name):
    env_name = config.get("environments", {}).get(name) or c.get("conda_envs", {}).get(name, name)
    if not ENV_DIR: return env_name
    return os.path.join(ENV_DIR, env_name)

PL_RES = full_path(c.get("results_dir", "results/plasmid_analysis"))
PL_LOGS = full_path(c.get("logs_dir", os.path.join("logs", "plasmid_analysis")))

# Input Handling (Typically from Long-read Assembly)
IN_DIR = full_path(c.get("input_dir", "results/assembly"))
SAMPLES_DICT = config.get("samples", {})

# 2. Target Registration (Smart Filtering)
# ------------------------------------------------------------------------------
if c.get("enable", False):
    for sample in SAMPLES_DICT:
        # Only analyze if sample has long-read data (as geNomad is usually for 3rd gen genomes here)
        if "long_read" in SAMPLES_DICT[sample]:
            TARGETS.append(os.path.join(PL_RES, sample, f"{sample}_genomad_summary.tsv"))

# 3. Rules
# ------------------------------------------------------------------------------

rule plasmid_genomad:
    input:
        fa = lambda wildcards: full_path(os.path.join(IN_DIR, wildcards.sample, "long_read/myloasm", wildcards.sample, f"{wildcards.sample}.myloasm.fa")),
        db = c.get("genomad", {}).get("db", "")
    output:
        summary = os.path.join(PL_RES, "{sample}", "{sample}_genomad_summary.tsv"),
        outdir = directory(os.path.join(PL_RES, "{sample}", "genomad_out"))
    log:
        os.path.join(PL_LOGS, "genomad/{sample}.log")
    threads: get_threads(c.get("genomad", {}).get("threads", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = get_env("genomad"),
        min_score = c.get("genomad", {}).get("min_score", 0.7)
    shell:
        """
        set -euo pipefail
        rm -rf {output.outdir}
        conda run -p {params.env} genomad end-to-end \
            --threads {threads} \
            --min-score {params.min_score} \
            --cleanup \
            {input.fa} {output.outdir} {input.db} 2>&1 | tee {log}
        # Extract summary to fixed location
        cp {output.outdir}/*_summary/*_summary.tsv {output.summary} 2>/dev/null || touch {output.summary}
        
        """
