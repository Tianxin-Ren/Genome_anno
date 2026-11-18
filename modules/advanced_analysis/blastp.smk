# BLASTp Targeted Gene Discovery in Advanced Analysis
# ------------------------------------------------------------------------------

# Retrieve settings from ADV_C namespaced under 'blastp'
BLASTP_C = ADV_C.get("blastp", {})

# Setup task dictionary (multiple search tasks)
BLASTP_TASKS = BLASTP_C.get("tasks", {})
if not BLASTP_TASKS:
    default_query = config.get("target_protein") or BLASTP_C.get("blastp_query") or "config/target_genes.faa"
    BLASTP_TASKS = {
        "default": {
            "query": default_query,
            "blastp_tool": BLASTP_C.get("blastp_tool", "blastp"),
            "blastp_evalue": BLASTP_C.get("blastp_evalue", 1e-5),
            "blastp_identity": BLASTP_C.get("blastp_identity", 30.0),
            "blastp_coverage": BLASTP_C.get("blastp_coverage", 50.0)
        }
    }

# Setup paths and environment resolution
ext_batch = config.get("external_genomes_batch", "")
EXT_BASE = config.get("annotation", {}).get("external_results_dir", "results/external_annotation")
if ext_batch:
    EXT_RES = os.path.join(config.get("project_root", "."), EXT_BASE, ext_batch)
else:
    EXT_RES = os.path.join(config.get("project_root", "."), EXT_BASE)
EXT_GENOMES = config.get("external_genomes", {})
if EXT_GENOMES is None: EXT_GENOMES = {}

def adv_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(config.get("project_root", "."), rel_path))

def adv_get_env(name):
    # Try advanced conda environments first, fallback to annotation conda environments
    env_name = ADV_C.get("conda_envs", {}).get(name, config.get("annotation", {}).get("conda_envs", {}).get(name, name))
    env_dir = config.get("conda_env_dir", "")
    if not env_dir: return env_name
    return os.path.join(env_dir, env_name)

def get_batch_sentinel_if_applicable(genome):
    if "/" in genome:
        bg = genome.split("/")[0]
        return [os.path.join(EXT_RES, f"{bg}_batch_hmm_completed.txt")]
    return []


def get_blastp_faa_input(wildcards):
    genome = wildcards.genome
    if "/" in genome:
        bg = genome.split("/")[0]
        checkpoints.ext_annotation_batch_prodigal_combine.get(batch_group=bg)
    return os.path.join(EXT_RES, f"{genome}/{genome}.faa")


def get_blastp_gff_input(wildcards):
    genome = wildcards.genome
    if "/" in genome:
        bg = genome.split("/")[0]
        checkpoints.ext_annotation_batch_prodigal_combine.get(batch_group=bg)
    return os.path.join(EXT_RES, f"{genome}/{genome}.gff3")


def get_db_done(wildcards):
    tool = str(BLASTP_TASKS[wildcards.task_name].get("blastp_tool", "blastp")).lower()
    return os.path.join(EXT_RES, wildcards.genome, "blastp/db", f"{wildcards.genome}_{tool}_db.done")

rule advanced_blastp_makeblastdb:
    input:
        faa = get_blastp_faa_input,
        sentinel = lambda wildcards: get_batch_sentinel_if_applicable(wildcards.genome)
    output:
        done = os.path.join(EXT_RES, "{genome}/blastp/db/{genome}_{tool}_db.done")
    log:
        os.path.join(ADV_LOGS, "blastp_db/{genome}_{tool}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        tool = lambda wildcards: wildcards.tool,
        db_prefix = lambda wildcards: os.path.join(EXT_RES, wildcards.genome, "blastp/db", f"{wildcards.genome}_{wildcards.tool}_db"),
        env = lambda wildcards: adv_get_env("blast") if wildcards.tool == "blastp" else adv_get_env("diamond")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname "{output.done}")
        
        if [ "{params.tool}" = "diamond" ]; then
            echo "Compiling Diamond database..."
            conda run -p {params.env} diamond makedb --in {input.faa} --db {params.db_prefix} 2>&1 | tee {log}
        else
            echo "Compiling NCBI BLAST database..."
            conda run -p {params.env} makeblastdb -in {input.faa} -dbtype prot -out {params.db_prefix} 2>&1 | tee {log}
        fi
        
        touch {output.done}
        }} 2>&1 | tee {log}
        """

rule advanced_blastp_search:
    input:
        db_done = get_db_done,
        query = lambda wildcards: adv_full_path(BLASTP_TASKS[wildcards.task_name]["query"])
    output:
        raw_tsv = os.path.join(EXT_RES, "{genome}/blastp/{task_name}/{genome}_vs_targets_raw.tsv")
    log:
        os.path.join(ADV_LOGS, "blastp_search/{task_name}/{genome}.log")
    threads: get_threads(16)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        tool = lambda wildcards: str(BLASTP_TASKS[wildcards.task_name].get("blastp_tool", "blastp")).lower(),
        db_prefix = lambda wildcards: os.path.join(EXT_RES, wildcards.genome, "blastp/db", f"{wildcards.genome}_{str(BLASTP_TASKS[wildcards.task_name].get('blastp_tool', 'blastp')).lower()}_db"),
        evalue = lambda wildcards: BLASTP_TASKS[wildcards.task_name].get("blastp_evalue", 1e-5),
        env = lambda wildcards: adv_get_env("blast") if str(BLASTP_TASKS[wildcards.task_name].get("blastp_tool", "blastp")).lower() == "blastp" else adv_get_env("diamond")
    shell:
        r"""
        {{
        set -euo pipefail
        
        # Check if query file exists and has size
        if [ ! -s "{input.query}" ]; then
            echo "Query file {input.query} is missing or empty. Creating empty output."
            touch {output.raw_tsv}
            exit 0
        fi

        if [ "{params.tool}" = "diamond" ]; then
            echo "Running Diamond blastp..."
            conda run -p {params.env} diamond blastp \
                -q {input.query} \
                -d {params.db_prefix} \
                -o {output.raw_tsv} \
                -e {params.evalue} \
                -p {threads} \
                --sensitive \
                --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen 2>&1 | tee {log}
        else
            echo "Running NCBI blastp..."
            conda run -p {params.env} blastp \
                -query {input.query} \
                -db {params.db_prefix} \
                -out {output.raw_tsv} \
                -evalue {params.evalue} \
                -num_threads {threads} \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" 2>&1 | tee {log}
        fi
        }} 2>&1 | tee {log}
        """

rule advanced_blastp_filter:
    input:
        raw_tsv = os.path.join(EXT_RES, "{genome}/blastp/{task_name}/{genome}_vs_targets_raw.tsv"),
        gff3 = get_blastp_gff_input,
        faa = get_blastp_faa_input,
        ffn = os.path.join(EXT_RES, "{genome}/{genome}.ffn"),
        sentinel = lambda wildcards: get_batch_sentinel_if_applicable(wildcards.genome)
    output:
        annotated_tsv = os.path.join(EXT_RES, "{genome}/blastp/{task_name}/{genome}_target_genes_annotated.tsv"),
        matched_faa = os.path.join(EXT_RES, "{genome}/blastp/{task_name}/{genome}_matched_target_genes.faa"),
        matched_ffn = os.path.join(EXT_RES, "{genome}/blastp/{task_name}/{genome}_matched_target_genes.ffn")
    log:
        os.path.join(ADV_LOGS, "blastp_filter/{task_name}/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_get_env("biopython"),
        identity = lambda wildcards: BLASTP_TASKS[wildcards.task_name].get("blastp_identity", 30.0),
        coverage = lambda wildcards: BLASTP_TASKS[wildcards.task_name].get("blastp_coverage", 50.0),
        evalue = lambda wildcards: BLASTP_TASKS[wildcards.task_name].get("blastp_evalue", 1e-5)
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} python3 scripts/filter_blastp_results.py \
            --raw_tsv {input.raw_tsv} \
            --gff3 {input.gff3} \
            --faa {input.faa} \
            --ffn {input.ffn} \
            --binned none \
            --evalue {params.evalue} \
            --identity {params.identity} \
            --coverage {params.coverage} \
            --out_tsv {output.annotated_tsv} \
            --out_faa {output.matched_faa} \
            --out_ffn {output.matched_ffn}
        }} 2>&1 | tee {log}
        """

rule advanced_global_blastp_summary:
    input:
        query_faa = lambda wildcards: adv_full_path(BLASTP_TASKS[wildcards.task_name]["query"]),
        reports = lambda wildcards: expand(os.path.join(EXT_RES, "{genome}/blastp/{task_name}/{genome}_target_genes_annotated.tsv"), genome=EXT_GENOMES.keys(), task_name=wildcards.task_name)
    output:
        tidy = os.path.join(EXT_RES, "blastp/{task_name}/global_target_genes_tidy_summary.tsv"),
        matrix = os.path.join(EXT_RES, "blastp/{task_name}/global_target_genes_matrix_summary.tsv")
    log:
        os.path.join(ADV_LOGS, "global_blastp_summary/{task_name}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_get_env("biopython"),
        genomes = lambda wildcards: ",".join(EXT_GENOMES.keys())
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} python3 scripts/generate_global_blastp_summary.py \
            --query_faa {input.query_faa} \
            --reports {input.reports} \
            --genomes "{params.genomes}" \
            --out_tidy {output.tidy} \
            --out_matrix {output.matrix}
        }} 2>&1 | tee {log}
        """

rule advanced_operon_summary:
    input:
        reports = lambda wildcards: expand(os.path.join(EXT_RES, "{genome}/blastp/{task_name}/{genome}_target_genes_annotated.tsv"), genome=EXT_GENOMES.keys(), task_name=wildcards.task_name)
    output:
        summary = os.path.join(EXT_RES, "blastp/{task_name}/global_operon_prediction_summary.tsv")
    log:
        os.path.join(ADV_LOGS, "operon_summary/{task_name}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_get_env("biopython"),
        max_span = 50000
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} python3 scripts/analyze_operons.py \
            --reports {input.reports} \
            --max_span {params.max_span} \
            --out_summary {output.summary}
        }} 2>&1 | tee {log}
        """
