# Preprocess Module: Quality Control (SR and LR)
import os

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
PP_ROOT = config.get("project_root", ".")
PP_ENV_DIR = config.get("conda_env_dir", "")
PP_C = config.get("preprocess", {})

def pp_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(PP_ROOT, rel_path))

def pp_get_env(name):
    env_name = config.get("environments", {}).get(name) or PP_C.get("conda_envs", {}).get(name, name)
    if not PP_ENV_DIR or os.path.isabs(env_name): return env_name
    return os.path.join(PP_ENV_DIR, env_name)

PP_RES = pp_full_path(PP_C.get("results_dir", "results/preprocess"))
PP_LOGS = pp_full_path(PP_C.get("logs_dir", os.path.join("logs", "preprocess")))

# Sample discovery
PP_SAMPLES = config.get("samples", {})

# 2. Target Registration
# ------------------------------------------------------------------------------
for sample in PP_SAMPLES:
    info = PP_SAMPLES[sample]
    if PP_C.get("steps", {}).get("short_qc", True) and "short_read" in info:
        TARGETS.append(os.path.join(PP_RES, f"{sample}/short_read/{sample}_R1.clean.fastq.gz"))
    
    if PP_C.get("steps", {}).get("long_qc", True) and "long_read" in info:
        TARGETS.append(os.path.join(PP_RES, f"{sample}/long_read/nanoplot/NanoPlot-report.html"))

# 3. Rules
# ------------------------------------------------------------------------------

ruleorder: preprocess_import_external > preprocess_fastp

def get_external_r1(wildcards):
    path = PP_SAMPLES[wildcards.sample].get("short_read", {}).get("external_clean_r1")
    return pp_full_path(path) if path else "EXTERNAL_R1_NOT_PROVIDED"

def get_external_r2(wildcards):
    path = PP_SAMPLES[wildcards.sample].get("short_read", {}).get("external_clean_r2")
    return pp_full_path(path) if path else "EXTERNAL_R2_NOT_PROVIDED"

rule preprocess_import_external:
    input:
        r1 = get_external_r1,
        r2 = get_external_r2
    output:
        r1 = os.path.join(PP_RES, "{sample}/short_read/{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(PP_RES, "{sample}/short_read/{sample}_R2.clean.fastq.gz"),
        json = os.path.join(PP_RES, "{sample}/short_read/{sample}.fastp.json"),
        html = os.path.join(PP_RES, "{sample}/short_read/{sample}.fastp.html")
    log:
        os.path.join(PP_LOGS, "import_external/{sample}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname {output.r1})
        ln -sf {input.r1} {output.r1}
        ln -sf {input.r2} {output.r2}
        touch {output.json} {output.html}
        }} 2>&1 | tee {log}
        """
def get_raw_r1(wildcards):
    path = PP_SAMPLES[wildcards.sample].get("short_read", {}).get("r1")
    return pp_full_path(path) if path else "RAW_R1_NOT_PROVIDED"

def get_raw_r2(wildcards):
    path = PP_SAMPLES[wildcards.sample].get("short_read", {}).get("r2")
    return pp_full_path(path) if path else "RAW_R2_NOT_PROVIDED"

rule preprocess_fastp:
    input:
        r1 = get_raw_r1,
        r2 = get_raw_r2
    output:
        r1 = os.path.join(PP_RES, "{sample}/short_read/{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(PP_RES, "{sample}/short_read/{sample}_R2.clean.fastq.gz"),
        json = os.path.join(PP_RES, "{sample}/short_read/{sample}.fastp.json"),
        html = os.path.join(PP_RES, "{sample}/short_read/{sample}.fastp.html")
    log:
        os.path.join(PP_LOGS, "fastp/{sample}.log")
    threads: get_threads(PP_C.get("threads", {}).get("fastp", 4))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = pp_get_env("bio")
    shell:
        """
        set -euo pipefail
        if [ -z "{input.r1}" ]; then echo "Error: No short_read input for {wildcards.sample}"; exit 1; fi
        conda run -p {params.env} fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} \
            --thread {threads} 2>&1 | tee {log}
        
        """
rule preprocess_nanoplot:
    input:
        lr = lambda wildcards: pp_full_path(PP_SAMPLES[wildcards.sample].get("long_read", {}).get("lr") or 
                                            PP_SAMPLES[wildcards.sample].get("long_read", {}).get("fastq")) or []
    output:
        report = os.path.join(PP_RES, "{sample}/long_read/nanoplot/NanoPlot-report.html")
    log:
        os.path.join(PP_LOGS, "nanoplot/{sample}.log")
    threads: get_threads(PP_C.get("threads", {}).get("nanoplot", 4))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = pp_get_env("nanoplot"),
        outdir = lambda wildcards: os.path.join(PP_RES, wildcards.sample, "long_read/nanoplot")
    shell:
        """
        set -euo pipefail
        if [ -z "{input.lr}" ]; then echo "Error: No long_read input for {wildcards.sample}"; exit 1; fi
        FORMAT_ARG="--fastq {input.lr}"
        if [[ "{input.lr}" == *.bam ]]; then
            FORMAT_ARG="--bam {input.lr}"
        fi
        conda run -p {params.env} NanoPlot \
            --threads {threads} \
            $FORMAT_ARG \
            -o {params.outdir} 2>&1 | tee {log}
        
        """
