# Data Preparation Workflow
# Parses GTDB metadata and downloads corresponding genomes via NCBI Datasets CLI in chunks

import os
import json

# --- Configuration Setup ---
DATA_PREP_CONF = config.get("data_prep", {})
PROJ = os.path.abspath(config.get("project_root", "."))
DATA_PREP_TASKS = DATA_PREP_CONF.get("tasks", {})

# Environments
DATA_PREP_ENV_DIR = config.get("conda_env_dir", "")

def data_prep_get_env(name):
    env_name = DATA_PREP_CONF.get("conda_envs", {}).get(name, name)
    if not DATA_PREP_ENV_DIR: return env_name
    return os.path.normpath(os.path.join(DATA_PREP_ENV_DIR, env_name))

ENV_DATASETS = data_prep_get_env("datasets_cli")

# Paths
DATA_PREP_RES = os.path.join(PROJ, DATA_PREP_CONF.get("results_dir", "results/data_prep"))
DATA_PREP_LOGS = os.path.join(
    PROJ,
    DATA_PREP_CONF.get("logs_dir", "logs/data_prep"),
    "{task_name}",
)

# Register global ncbi_conn resource limit (default: 5)
NCBI_CONN_LIMIT = int(DATA_PREP_CONF.get("ncbi_conn", config.get("ncbi_conn", 5)))
if hasattr(workflow, "global_resources"):
    if "ncbi_conn" not in workflow.global_resources:
        if hasattr(workflow, "register_resource"):
            workflow.register_resource("ncbi_conn", NCBI_CONN_LIMIT)
        else:
            workflow.global_resources["ncbi_conn"] = NCBI_CONN_LIMIT

# Register targets for all tasks
for task_name, task_info in DATA_PREP_TASKS.items():
    if "genes" in task_info:
        TARGETS.append(os.path.join(DATA_PREP_RES, f"{task_name}/genes_collected.done"))
    elif not task_info.get("download", task_info.get("download_genomes", True)):
        TARGETS.append(os.path.join(DATA_PREP_RES, f"{task_name}/01_filter/filtered_metadata.tsv"))
    else:
        TARGETS.append(os.path.join(DATA_PREP_RES, f"{task_name}/genomes_downloaded.done"))

# ------------------------------------------------------------------------------
# Checkpoint: Filter GTDB Metadata and split into chunks
# ------------------------------------------------------------------------------
checkpoint filter_gtdb_metadata:
    output:
        chunks_dir = directory(os.path.join(DATA_PREP_RES, "{task_name}/01_filter/chunks")),
        filtered_tsv = os.path.join(DATA_PREP_RES, "{task_name}/01_filter/filtered_metadata.tsv")
    log:
        os.path.join(DATA_PREP_LOGS, "filter_gtdb.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ENV_DATASETS,
        metadata = lambda wildcards: DATA_PREP_TASKS.get(wildcards.task_name, {}).get("gtdb_metadata", ""),
        query = lambda wildcards: DATA_PREP_TASKS.get(wildcards.task_name, {}).get("filter_query", ""),
        accessions = lambda wildcards: DATA_PREP_TASKS.get(wildcards.task_name, {}).get("accession_file", ""),
        chunk_size = lambda wildcards: DATA_PREP_TASKS.get(wildcards.task_name, {}).get("chunk_size", 100),
        script = os.path.join(MODULES_DIR, "data_prep", "scripts", "filter_gtdb.py")
    shell:
        """
        set -euo pipefail
        ACCESSION_ARG=""
        if [ -n "{params.accessions}" ]; then
            ACCESSION_ARG="--accessions {params.accessions}"
        fi
        conda run -p {params.env} python {params.script} \
            -i {params.metadata} \
            -o {output.chunks_dir} \
            -t {output.filtered_tsv} \
            -q "{params.query}" \
            $ACCESSION_ARG \
            --chunk-size {params.chunk_size} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# 2. Download Genomes with NCBI Datasets (Runs in parallel for each chunk)
# ------------------------------------------------------------------------------
rule download_ncbi_genomes:
    input:
        chunk = os.path.join(DATA_PREP_RES, "{task_name}/01_filter/chunks/chunk_{chunk_id}.txt")
    output:
        dataset_zip = os.path.join(DATA_PREP_RES, "{task_name}/02_download/chunk_{chunk_id}.zip")
    log:
        os.path.join(DATA_PREP_LOGS, "download_chunk_{chunk_id}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads),
        ncbi_conn=1
    params:
        env = ENV_DATASETS,
        api_key = DATA_PREP_CONF.get("ncbi_api_key", os.environ.get("NCBI_API_KEY", ""))
    shell:
        """
        set -euo pipefail

        # Check if the chunk file is empty (e.g. chunk_000.txt for empty search results)
        if [ ! -s {input.chunk} ]; then
            echo "Chunk file is empty. Skipping download." 2>&1 | tee {log}
            touch {output.dataset_zip}
            exit 0
        fi

        API_ARG=""
        if [ -n "{params.api_key}" ]; then
            API_ARG="--api-key {params.api_key}"
        fi

        # Jitter: sleep random seconds (1-5s) to avoid API burst when multiple parallel jobs start
        sleep $(( (RANDOM % 5) + 1 ))

        # Download with exponential backoff retry loop (handles NCBI gateway timeouts / transient 502/504)
        MAX_RETRIES=4
        DELAY=10
        SUCCESS=0

        for attempt in $(seq 1 $MAX_RETRIES); do
            echo "Download attempt $attempt/$MAX_RETRIES for {input.chunk}..." 2>&1 | tee -a {log}
            if conda run -p {params.env} datasets download genome accession \
                --inputfile {input.chunk} \
                --include genome \
                $API_ARG \
                --filename {output.dataset_zip} 2>&1 | tee -a {log}; then
                
                # Check that output zip was created and is non-empty
                if [ -s {output.dataset_zip} ]; then
                    echo "Download successful on attempt $attempt." 2>&1 | tee -a {log}
                    SUCCESS=1
                    break
                fi
            fi

            echo "Attempt $attempt failed. Cleaning up partial download and waiting $DELAY seconds before retry..." 2>&1 | tee -a {log}
            rm -f {output.dataset_zip}
            sleep $DELAY
            DELAY=$((DELAY * 2))
        done

        if [ $SUCCESS -ne 1 ]; then
            echo "Error: Failed to download {input.chunk} after $MAX_RETRIES attempts." 2>&1 | tee -a {log}
            rm -f {output.dataset_zip}
            exit 1
        fi
        """
def get_chunk_zips(wildcards):
    # Obtain the dynamic chunks directory once checkpoint filter_gtdb_metadata is done
    chunks_dir = checkpoints.filter_gtdb_metadata.get(**wildcards).output.chunks_dir
    
    # Extract chunk ids (e.g. chunk_001.txt -> 001)
    chunk_ids = glob_wildcards(os.path.join(chunks_dir, "chunk_{chunk_id}.txt")).chunk_id
    
    # Return the expected zip file paths for each chunk
    return [os.path.join(DATA_PREP_RES, f"{wildcards.task_name}/02_download/chunk_{chunk_id}.zip") for chunk_id in chunk_ids]

# ------------------------------------------------------------------------------
# 3. Extract and Aggregate all chunk zip packages
# ------------------------------------------------------------------------------
rule extract_genomes:
    input:
        dataset_zips = get_chunk_zips
    output:
        done_flag = os.path.join(DATA_PREP_RES, "{task_name}/genomes_downloaded.done")
    log:
        os.path.join(DATA_PREP_LOGS, "extract_genomes.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        out_dir = os.path.join(DATA_PREP_RES, "{task_name}/genomes")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.out_dir}
        # Clear out log
        echo "Starting extraction..." > {log}
        # Iterate and extract each non-empty zip file
        for zip_file in {input.dataset_zips}; do
            if [ ! -s "$zip_file" ]; then
                echo "Skipping empty or missing zip: $zip_file" 2>&1 | tee -a {log}
                continue
            fi
            echo "Extracting chunk: $zip_file" 2>&1 | tee -a {log}
            TMP_DIR=$(mktemp -d -t ncbi_extract_XXXXXX)
            unzip -q "$zip_file" -d $TMP_DIR 2>&1 | tee -a {log}
            # Move all .fna files from the nested ncbi_dataset/data/GC*_*/ folders to the final output dir
            find $TMP_DIR -name "*.fna" -exec cp {{}} {params.out_dir}/ \\; 2>&1 | tee -a {log}
            # Clean up temporary space
            rm -rf $TMP_DIR
        done
        # Touch sentinel file
        touch {output.done_flag}
        
        """

# ------------------------------------------------------------------------------
# Rule to collect high-confidence gene sequences from UniProt
# ------------------------------------------------------------------------------
rule collect_genes:
    output:
        done_flag = os.path.join(DATA_PREP_RES, "{task_name}/genes_collected.done")
    log:
        os.path.join(DATA_PREP_LOGS, "collect_genes.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ENV_DATASETS,
        genes_json = lambda wildcards: json.dumps(DATA_PREP_TASKS.get(wildcards.task_name, {}).get("genes", [])),
        kegg_tsv = lambda wildcards: DATA_PREP_TASKS.get(wildcards.task_name, {}).get("kegg_tsv", ""),
        min_sequences = lambda wildcards: DATA_PREP_TASKS.get(wildcards.task_name, {}).get("min_sequences", 5),
        max_sequences = lambda wildcards: DATA_PREP_TASKS.get(wildcards.task_name, {}).get("max_sequences", 7),
        out_dir = lambda wildcards: os.path.join(DATA_PREP_RES, wildcards.task_name, "genes"),
        script = os.path.join(MODULES_DIR, "data_prep", "scripts", "collect_genes.py")
    shell:
        """
        set -euo pipefail
        
        INPUT_ARGS=""
        if [ -n "{params.kegg_tsv}" ]; then
            INPUT_ARGS="--kegg-tsv {params.kegg_tsv}"
        else
            INPUT_ARGS="--genes-json '{params.genes_json}'"
        fi
        
        conda run -p {params.env} python {params.script} \
            $INPUT_ARGS \
            --min-sequences {params.min_sequences} \
            --max-sequences {params.max_sequences} \
            --output-dir {params.out_dir} \
            --done-flag {output.done_flag} 2>&1 | tee {log}
        """
