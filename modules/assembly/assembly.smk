# Assembly Module: MEGAHIT (SR) and Myloasm (LR)
import os

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
# Sample discovery
AS_SAMPLES = config.get("samples", {})
if AS_SAMPLES is None: AS_SAMPLES = {}

# Robust constraint: Only allow actual sample names
import re
AS_SAMPLE_LIST = [re.escape(s) for s in AS_SAMPLES.keys()]
wildcard_constraints:
    sample = "|".join(AS_SAMPLE_LIST) if AS_SAMPLE_LIST else "[^/]+"

AS_ROOT = config.get("project_root", ".")
AS_ENV_DIR = config.get("conda_env_dir", "")
AS_C = config.get("assembly", {})

def as_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(AS_ROOT, rel_path))

def as_get_env(name):
    env_name = config.get("environments", {}).get(name) or AS_C.get("conda_envs", {}).get(name, name)
    if not AS_ENV_DIR or os.path.isabs(env_name): return env_name
    return os.path.join(AS_ENV_DIR, env_name)

AS_RES = as_full_path(AS_C.get("results_dir", "results/assembly"))
AS_LOGS = as_full_path(AS_C.get("logs_dir", os.path.join("logs", "assembly")))

# Dependencies from Preprocess
AS_PP_RES = as_full_path(config.get("preprocess", {}).get("results_dir", "results/preprocess"))

def as_get_reads(wildcards, side):
    return os.path.join(AS_PP_RES, f"{wildcards.sample}/short_read/{wildcards.sample}_R{side}.clean.fastq.gz")


def get_raw_assembly(wildcards):
    info = AS_SAMPLES.get(wildcards.sample, {})
    if "external_assembly" in info:
        return os.path.join(AS_RES, f"{wildcards.sample}/external/{wildcards.sample}.fa")
    if "long_read" in info:
        return os.path.join(AS_RES, f"{wildcards.sample}/long_read/myloasm/{wildcards.sample}/{wildcards.sample}.myloasm.fa")
    return os.path.join(AS_RES, f"{wildcards.sample}/megahit/{wildcards.sample}.contigs.fa")

def get_filtered_assembly(sample):
    min_len = AS_C.get("min_len", 1000)
    class DummyWildcards: pass
    w = DummyWildcards()
    w.sample = sample
    raw_fa = get_raw_assembly(w)
    
    # [NEW] Skip filtering for long read assemblies (myloasm)
    if "myloasm" in raw_fa:
        return raw_fa
        
    base, ext = os.path.splitext(raw_fa)
    if not ext: ext = ".fa"
    return f"{base}_filtered_{min_len}{ext}"

# 2. Target Registration
# ------------------------------------------------------------------------------
for sample in AS_SAMPLES:
    TARGETS.append(get_filtered_assembly(sample))

# 3. Rules
# ------------------------------------------------------------------------------

rule assembly_import_external:
    input:
        fa = lambda wildcards: as_full_path(AS_SAMPLES[wildcards.sample]["external_assembly"])
    output:
        fa = os.path.join(AS_RES, "{sample}/external/{sample}.fa")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        "mkdir -p $(dirname {output.fa}) && ln -sf {input.fa} {output.fa}"
rule assembly_filter:
    input:
        fa = "{prefix}{ext}"
    output:
        fa = "{prefix}_filtered_{min_len}{ext}"
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = as_get_env("bio")
    wildcard_constraints:
        ext=r"\.(fa|fasta|fna)",
        min_len=r"\d+"
    shell:
        "conda run -p {params.env} seqkit seq -m {wildcards.min_len} {input.fa} -o {output.fa}"
rule assembly_megahit:
    input:
        r1 = lambda wildcards: as_get_reads(wildcards, "1"),
        r2 = lambda wildcards: as_get_reads(wildcards, "2")
    output:
        fa = os.path.join(AS_RES, "{sample}/megahit/{sample}.contigs.fa")
    log:
        os.path.join(AS_LOGS, "megahit/{sample}.log")
    threads: get_threads(AS_C.get("threads", {}).get("megahit", AS_C.get("threads", 64)))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = as_get_env("megahit"),
        out_dir = lambda wildcards: os.path.join(AS_RES, wildcards.sample, "megahit")
    shell:
        """
        set -euo pipefail
        rm -rf "{params.out_dir}"
        conda run -p {params.env} megahit \
            -1 {input.r1} -2 {input.r2} \
            -o "{params.out_dir}" \
            -t {threads} --out-prefix {wildcards.sample} 2>&1 | tee "{log}"
        
        """
rule assembly_myloasm:
    input:
        lr = lambda wildcards: as_full_path(AS_SAMPLES[wildcards.sample].get("long_read", {}).get("lr") or 
                                            AS_SAMPLES[wildcards.sample].get("long_read", {}).get("fastq", "")) or []
    output:
        fa = os.path.join(AS_RES, "{sample}/long_read/myloasm/{sample}/{sample}.myloasm.fa")
    log:
        os.path.join(AS_LOGS, "myloasm/{sample}.log")
    threads: get_threads(AS_C.get("threads", {}).get("myloasm", AS_C.get("threads", 64)))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = as_get_env("myloasm"),
        extra = AS_C.get("extra", ""),
        out_base = lambda wildcards: os.path.join(AS_RES, wildcards.sample, "long_read/myloasm", wildcards.sample)
    shell:
        """
        set -euo pipefail
        if [ -z "{input.lr}" ]; then echo "Error: No long_read input for {wildcards.sample}"; exit 1; fi
        mkdir -p "{params.out_base}"
        conda run -p {params.env} myloasm \
            -t {threads} \
            -o "{params.out_base}" \
            {params.extra} \
            "{input.lr}" 2>&1 | tee "{log}"
        if [ -f "{params.out_base}/assembly_primary.fa" ]; then
            mv "{params.out_base}/assembly_primary.fa" "{output.fa}"
        elif [ -f "{params.out_base}/final.assembly.fa" ]; then
            mv "{params.out_base}/final.assembly.fa" "{output.fa}"
        else
            echo "Error: myloasm failed to produce assembly_primary.fa or final.assembly.fa"
            exit 1
        fi
        
        """
