# ==============================================================================
# Advanced Analysis Module: Standalone HMM Construction Pipeline
# ==============================================================================
import os
import glob

# Load independent configuration (Managed by Snakefile registry)
# configfile: "config/hmm_build.yaml"
HMM_C = config.get("hmm_build", {})

def hmm_get_env(name):
    env_name = HMM_C.get("conda_envs", {}).get(name, name)
    env_dir = config.get("conda_env_dir", "")
    if not env_dir or os.path.isabs(env_name):
        return env_name
    return os.path.join(env_dir, env_name)

def hmm_full_path(path):
    if not path or os.path.isabs(path):
        return path
    return os.path.join(config.get("project_root", "."), path)

# Directories
HMM_IN_DIR = hmm_full_path(HMM_C.get("input_dir", "data/hmm_raw_inputs"))
HMM_OUT_DIR = hmm_full_path(HMM_C.get("output_dir") or HMM_C.get("results_dir", "results/hmm_build"))
HMM_LOG_DIR = hmm_full_path(HMM_C.get("logs_dir", "logs/hmm_build"))
HMM_TMP_DIR = os.path.join(HMM_OUT_DIR, "temp_processing")

# Target collection
HMM_TARGETS = []
HMM_SUFFIX = HMM_C.get("input_suffix", "_combined.faa")

# Detect input files and create targets
if os.path.exists(HMM_IN_DIR):
    raw_files = glob.glob(os.path.join(HMM_IN_DIR, f"*{HMM_SUFFIX}"))
    basenames = [os.path.basename(f).replace(HMM_SUFFIX, "") for f in raw_files]
    for b in basenames:
        HMM_TARGETS.append(os.path.join(HMM_OUT_DIR, f"{b}.hmm"))

# Constraints
wildcard_constraints:
    gene="[^/]+"

# --- Rules ---

# Step 0: Homolog Sequence Search (UniProt or NCBI)
rule hmm_search_sequences:
    input:
        ref_faa = os.path.join(HMM_IN_DIR, f"{{gene}}{HMM_SUFFIX}")
    output:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_searched.faa")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_search.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = hmm_get_env("biopython"),
        enable = HMM_C.get("params", {}).get("search", {}).get("enable", True),
        database = HMM_C.get("params", {}).get("search", {}).get("database", "uniprotkb_swissprot"),
        identity = HMM_C.get("params", {}).get("search", {}).get("identity", 0.3),
        coverage = HMM_C.get("params", {}).get("search", {}).get("coverage", 0.5),
        max_seqs = HMM_C.get("params", {}).get("search", {}).get("max_seqs", 200)
    shell:
        """
        set -euo pipefail
        if [ "{params.enable}" = "True" ]; then
            echo "Searching homologs for {wildcards.gene} using database {params.database}..." | tee {log}
            conda run -p {params.env} python3 scripts/hmm_search_sequences.py \
                --query {input.ref_faa} \
                --database {params.database} \
                --identity {params.identity} \
                --coverage {params.coverage} \
                --max_seqs {params.max_seqs} \
                --output {output.faa} 2>&1 | tee -a {log}
        else
            echo "Sequence search disabled. Copying reference input to output..." | tee {log}
            cp {input.ref_faa} {output.faa}
        fi
        
        """
rule hmm_clean_header:
    input:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_searched.faa")
    output:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_clean.faa")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_clean_header.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = hmm_get_env("seqkit")
    shell:
        """
        set -euo pipefail
        # Remove everything after the first space in the header
        conda run -p {params.env} seqkit replace -p " .*" -r "" {input.faa} -o {output.faa} 2>&1 | tee {log}
        
        """
rule hmm_filter_length:
    input:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_clean.faa")
    output:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_len.faa")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_filter_length.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = hmm_get_env("seqkit"),
        enable = HMM_C.get("params", {}).get("length_filter", {}).get("enable", False),
        adaptive = HMM_C.get("params", {}).get("length_filter", {}).get("adaptive", True),
        margin = HMM_C.get("params", {}).get("length_filter", {}).get("margin", 0.25),
        min_len = HMM_C.get("params", {}).get("length_filter", {}).get("min_len", 0),
        max_len = HMM_C.get("params", {}).get("length_filter", {}).get("max_len", 100000)
    shell:
        """
        set -euo pipefail
        if [ "{params.enable}" = "True" ]; then
            if [ "{params.adaptive}" = "True" ]; then
                # Extract the length of the first sequence (considered the reference sequence)
                REF_LEN=$(conda run -p {params.env} seqkit fx2tab -l -n -i {input.faa} | head -n 1 | cut -f 2)
                # Calculate adaptive thresholds (1-margin to 1+margin)
                MIN_L=$(python3 -c "print(int($REF_LEN * (1 - {params.margin})))")
                MAX_L=$(python3 -c "print(int($REF_LEN * (1 + {params.margin})))")
                echo "Adaptive Filter for {wildcards.gene}: Ref Length=$REF_LEN, Min=$MIN_L, Max=$MAX_L" | tee {log}
                conda run -p {params.env} seqkit seq -m $MIN_L -M $MAX_L {input.faa} -o {output.faa} 2>&1 | tee -a {log}
            else
                echo "Static Filter for {wildcards.gene}: Min={params.min_len}, Max={params.max_len}" | tee {log}
                conda run -p {params.env} seqkit seq -m {params.min_len} -M {params.max_len} {input.faa} -o {output.faa} 2>&1 | tee -a {log}
            fi
        else
            cp {input.faa} {output.faa}
            echo "Length filtering disabled. Copied input to output." > {log}
        fi
        
        """
rule hmm_cdhit_cluster:
    input:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_len.faa")
    output:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_nr.faa")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_cdhit.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = hmm_get_env("cd_hit"),
        enable = HMM_C.get("params", {}).get("cd_hit", {}).get("enable", False),
        ident = HMM_C.get("params", {}).get("cd_hit", {}).get("identity", 0.95),
        word = HMM_C.get("params", {}).get("cd_hit", {}).get("word_size", 5),
        mem = HMM_C.get("params", {}).get("cd_hit", {}).get("memory_limit", 8000)
    shell:
        """
        set -euo pipefail
        if [ "{params.enable}" = "True" ]; then
            conda run -p {params.env} cd-hit -i {input.faa} -o {output.faa} -c {params.ident} -n {params.word} -d 0 -M {params.mem} -T {threads} 2>&1 | tee {log}
        else
            cp {input.faa} {output.faa}
            echo "CD-HIT clustering disabled. Copied input to output." > {log}
        fi
        
        """
rule hmm_mafft_align:
    input:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_nr.faa")
    output:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_aligned.faa")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_mafft.log")
    threads: get_threads(8)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = hmm_get_env("mafft")
    shell:
        """
        set -euo pipefail
        # Align with the user specified localpair alignment options
        conda run -p {params.env} mafft --anysymbol --localpair --maxiterate 1000 --thread {threads} {input.faa} > {output.faa} 2> {log}
        
        """
rule hmm_clean_u_to_c:
    input:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_aligned.faa")
    output:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_aligned_clean.faa")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_clean_u_to_c.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        """
        set -euo pipefail
        # Replace Selenocysteine (U) with Cysteine (C) as it is invalid in HMMER
        sed 's/[uU]/C/g' {input.faa} > {output.faa} 2> {log}
        
        """
rule hmm_trimal_trim:
    input:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_aligned_clean.faa")
    output:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_trimmed.faa"),
        html = os.path.join(HMM_TMP_DIR, "{gene}_trimmed.html")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_trimal.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = hmm_get_env("trimal"),
        enable = HMM_C.get("params", {}).get("trimal", {}).get("enable", True),
        gt = HMM_C.get("params", {}).get("trimal", {}).get("gt", 0.5),
        st = HMM_C.get("params", {}).get("trimal", {}).get("st", 0.001),
        cons = HMM_C.get("params", {}).get("trimal", {}).get("cons", 60)
    shell:
        """
        set -euo pipefail
        if [ "{params.enable}" = "True" ]; then
            conda run -p {params.env} trimal -in {input.faa} -out {output.faa} -gt {params.gt} -st {params.st} -cons {params.cons} -htmlout {output.html} 2>&1 | tee {log}
        else
            echo "trimAl trimming disabled. Copying alignment to output..." | tee {log}
            cp {input.faa} {output.faa}
            touch {output.html}
        fi
        
        """
rule hmm_build_model:
    input:
        faa = os.path.join(HMM_TMP_DIR, "{gene}_trimmed.faa")
    output:
        hmm = os.path.join(HMM_OUT_DIR, "{gene}.hmm")
    log:
        os.path.join(HMM_LOG_DIR, "{gene}_hmmbuild.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = hmm_get_env("hmmer")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} hmmbuild -n {wildcards.gene} --cpu {threads} {output.hmm} {input.faa} 2>&1 | tee {log}
        
        """
