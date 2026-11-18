# Batch mode rules for external genome annotation (Sentinel-file based redirection)
# Output structure: {EXT_RES}/{batch_group}/{bin}/{bin}.faa  (nested under sample name)
import os
import re
import shlex

DEDUP_C = AN_C.get("protein_dedup", {}) or {}
DEDUP_ENABLED = str(DEDUP_C.get("enabled", False)).lower() == "true"
CLUSTER_99_ENABLED = DEDUP_ENABLED and str(DEDUP_C.get("cluster_99", True)).lower() == "true"


def batch_dedup_faa(wildcards):
    return os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_protein_dedup/unique_proteins.faa")


def batch_dedup_mapping(wildcards):
    return os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_protein_dedup/protein_map.tsv")


def batch_hmm_source(wildcards):
    if CLUSTER_99_ENABLED:
        return os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_protein_dedup/cluster_99/representatives.faa")
    if DEDUP_ENABLED:
        return batch_dedup_faa(wildcards)
    return os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_combined.faa")


def batch_cluster_mapping(wildcards):
    if CLUSTER_99_ENABLED:
        return os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_protein_dedup/cluster_99/clusters.tsv")
    return []


def batch_hmm_resolved_tbl(wildcards):
    if DEDUP_ENABLED:
        return os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_combined_hmm_expanded.tbl")
    return os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_combined_hmm.tbl")

# Constrain wildcard 'bin' to only match genomes belonging to a batch group
batch_bins = [
    name for name, info in EXT_GENOMES.items()
    if isinstance(info, dict) and "batch_group" in info
]

# Collect unique batch group names
_batch_groups = sorted(set(
    info.get("batch_group") for info in EXT_GENOMES.values()
    if isinstance(info, dict) and "batch_group" in info
))

wildcard_constraints:
    batch_group = "|".join(re.escape(bg) for bg in _batch_groups) if _batch_groups else "none_matched_placeholder"


def get_hmm_gather_inputs(wildcards):
    checkpoint_output = checkpoints.ext_annotation_batch_hmm_split.get(**wildcards).output[0]
    import glob
    chunk_files = glob.glob(os.path.join(checkpoint_output, "*.faa"))
    if len(chunk_files) == 0:
        # If the split folder exists but is empty (e.g., from a previously aborted run or cluster kill),
        # we delete the empty directory to force Snakemake to reschedule and rerun the checkpoint.
        try:
            if os.path.exists(checkpoint_output):
                import shutil
                shutil.rmtree(checkpoint_output)
        except Exception:
            pass
        # Re-trigger evaluation to let Snakemake raise HelperNotReady and queue the split job.
        checkpoint_output = checkpoints.ext_annotation_batch_hmm_split.get(**wildcards).output[0]
    tbls = []
    for f in chunk_files:
        chunk_name = os.path.splitext(os.path.basename(f))[0]
        tbls.append(os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_combined_chunks/{chunk_name}.tbl"))
    return sorted(tbls)


def get_batch_hmm_db_input(wildcards):
    batch_name = wildcards.batch_group
    batch_info = EXT_GENOMES.get(batch_name, {})
    db_param = batch_info.get("hmm_db", "")
    resolved = resolve_hmm_inputs(db_param)
    if not resolved:
        raise ValueError(f"No HMM files resolved for batch '{batch_name}' from HMM parameter '{db_param}'")
    if len(resolved) > 1:
        return os.path.join(EXT_RES, f"batch_{batch_name}_combined_unified.hmm")
    else:
        return resolved[0]


# Rule to combine global HMM database if there are multiple HMM profiles for a batch
rule ext_annotation_batch_cat_hmm:
    input:
        hmms = lambda wildcards: resolve_hmm_inputs(EXT_GENOMES.get(wildcards.batch_group, {}).get("hmm_db", ""))
    output:
        db = os.path.join(EXT_RES, "batch_{batch_group}_combined_unified.hmm")
    log:
        os.path.join(AN_LOGS, "batch_cat_unified_hmm_{batch_group}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    run:
        if len(input.hmms) > 1:
            shell("cat {input.hmms} > {output.db}")
        else:
            shell("cp {input.hmms[0]} {output.db}")


# 1. Prodigal + Combine FAA step (Checkpoint)
checkpoint ext_annotation_batch_prodigal_combine:
    input:
        # No direct file dependency to dynamically discover
    output:
        faa = os.path.join(EXT_RES, "batch_{batch_group}_combined.faa")
    log:
        os.path.join(AN_LOGS, "batch_prodigal_combine_{batch_group}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        genomes_dirs = lambda wildcards: f"--genomes-dir {an_full_path(EXT_GENOMES.get(wildcards.batch_group, {}).get('path'))}",
        output_dir = EXT_RES,
        batch_name = lambda wildcards: wildcards.batch_group,
        conda_env_dir = AN_ENV_DIR,
        bakta_env = AN_C.get("conda_envs", {}).get("bakta", "bakta_1.12.0"),
        biopython_env = AN_C.get("conda_envs", {}).get("biopython", "biopython_1.85"),
        prodigal_env_flag = lambda wildcards: f"--prodigal-env {AN_C.get('conda_envs', {}).get('prodigal')}" if AN_C.get("conda_envs", {}).get("prodigal") else "",
        output_ffn_flag = lambda wildcards: "--output-ffn" if str(EXT_GENOMES.get(wildcards.batch_group, {}).get("output_ffn", "False")).lower() == "true" else "",
        batch_size = 100000
    shell:
        r"""
        {{
            set -euo pipefail
            python3 scripts/batch_prodigal_hmm.py \
                {params.genomes_dirs} \
                --output-dir {params.output_dir} \
                --batch-name {params.batch_name} \
                --hmm-db none \
                --ko-list none \
                --conda-env-dir {params.conda_env_dir} \
                --bakta-env {params.bakta_env} \
                --biopython-env {params.biopython_env} \
                {params.prodigal_env_flag} \
                {params.output_ffn_flag} \
                --batch-size {params.batch_size} \
                --threads {threads} \
                --only-predict-combine \
                --combined-faa {output.faa}
        }} 2>&1 | tee {log}
        """


# 1b. Lossless protein deduplication. The target HMM database is searched on
# all unique proteins, while protein_map.tsv allows results to be expanded back
# to every original genome/protein ID before arbitration.
if DEDUP_ENABLED:
    rule ext_annotation_batch_protein_deduplicate:
        input:
            faa = os.path.join(EXT_RES, "batch_{batch_group}_combined.faa")
        output:
            unique_faa = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/unique_proteins.faa"),
            mapping = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/protein_map.tsv"),
            stats = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/stats.tsv"),
            database = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/catalog.sqlite")
        log:
            os.path.join(AN_LOGS, "batch_protein_deduplicate_{batch_group}.log")
        threads: get_threads(8)
        resources:
            runtime="7d",
            mem_mb=lambda wildcards, threads: get_mem_mb(threads)
        params:
            script = os.path.join(WORKFLOW_ROOT, "scripts/deduplicate_proteins.py"),
            python_env = an_get_env("python"),
            gff_root = lambda wildcards: os.path.join(EXT_RES, wildcards.batch_group)
        shell:
            r"""
            {{
                set -euo pipefail
                PYTHON_BIN="{params.python_env}/bin/python3"
                if [ ! -x "$PYTHON_BIN" ]; then PYTHON_BIN="python3"; fi
                "$PYTHON_BIN" {params.script} \
                    --input-faa {input.faa} \
                    --output-faa {output.unique_faa} \
                    --mapping {output.mapping} \
                    --stats {output.stats} \
                    --database {output.database} \
                    --gff-root {params.gff_root}
            }} 2>&1 | tee {log}
            """


if CLUSTER_99_ENABLED:
    rule ext_annotation_batch_cluster_99:
        input:
            unique_faa = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/unique_proteins.faa")
        output:
            representatives = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/cluster_99/representatives.faa"),
            clusters = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/cluster_99/clusters.tsv"),
            parameters = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/cluster_99/parameters.tsv")
        log:
            os.path.join(AN_LOGS, "batch_protein_cluster_99_{batch_group}.log")
        threads: get_threads(8)
        resources:
            runtime="7d",
            mem_mb=lambda wildcards, threads: get_mem_mb(threads)
        params:
            script = os.path.join(WORKFLOW_ROOT, "scripts/cluster_proteins.py"),
            env = an_get_env("mmseqs2"),
            python_env = an_get_env("python"),
            identity = DEDUP_C.get("identity", 0.99),
            coverage = DEDUP_C.get("coverage", 0.95),
            coverage_mode = DEDUP_C.get("coverage_mode", 0)
        shell:
            r"""
            {{
                set -euo pipefail
                mkdir -p "$(dirname {output.representatives})"
                MMSEQS_BIN="mmseqs"
                if [ -x "{params.env}/bin/mmseqs" ]; then MMSEQS_BIN="{params.env}/bin/mmseqs"; fi
                PYTHON_BIN="{params.python_env}/bin/python3"
                if [ ! -x "$PYTHON_BIN" ]; then PYTHON_BIN="python3"; fi
                "$PYTHON_BIN" {params.script} \
                    --input-faa {input.unique_faa} \
                    --output-dir "$(dirname {output.representatives})" \
                    --mmseqs "$MMSEQS_BIN" \
                    --identity {params.identity} \
                    --coverage {params.coverage} \
                    --coverage-mode {params.coverage_mode} \
                    --threads {threads}
            }} 2>&1 | tee {log}
            """


# 2. Split HMM step using seqkit (Dynamic checkpoint)
checkpoint ext_annotation_batch_hmm_split:
    input:
        faa = batch_hmm_source
    output:
        outdir = directory(os.path.join(EXT_RES, "batch_{batch_group}_combined_chunks"))
    log:
        os.path.join(AN_LOGS, "batch_hmm_split_{batch_group}.log")
    threads: get_threads(32)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        outdir = lambda wildcards: os.path.join(EXT_RES, f"batch_{wildcards.batch_group}_combined_chunks"),
        split_count = lambda wildcards: EXT_GENOMES.get(wildcards.batch_group, {}).get("hmm_split_chunks") or AN_C.get("hmm_split_chunks", 20),
        env = an_get_env("seqkit")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p "{params.outdir}/tmp"
        conda run -p {params.env} seqkit split "{input.faa}" -p {params.split_count} -j {threads} --two-pass -O "{params.outdir}/tmp" --force
        i=1
        for file in $(ls "{params.outdir}/tmp"/* | sort -V); do
            mv "$file" "{params.outdir}/chunk_${{i}}-of-{params.split_count}.faa"
            ((i++))
        done
        rm -rf "{params.outdir}/tmp"
        }} 2>&1 | tee {log}
        """


# 3. Running hmmsearch on each chunk
rule ext_annotation_batch_hmmsearch_chunk:
    input:
        faa = os.path.join(EXT_RES, "batch_{batch_group}_combined_chunks/{chunk}.faa"),
        db = get_batch_hmm_db_input
    output:
        tbl = os.path.join(EXT_RES, "batch_{batch_group}_combined_chunks/{chunk}.tbl")
    log:
        os.path.join(AN_LOGS, "batch_hmm_unified/{batch_group}_{chunk}.log")
    threads: get_threads(8)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("hmmer"),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5)
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} hmmsearch --cpu {threads} --tblout {output.tbl} -E {params.evalue} {input.db} {input.faa} > {log} 2>&1
        }} 2>&1 | tee {log}
        """


# 4. Gather chunk tables into a single tbl
rule ext_annotation_batch_hmm_gather:
    input:
        unified = get_hmm_gather_inputs
    output:
        unified_raw = os.path.join(EXT_RES, "batch_{batch_group}_combined_hmm.tbl")
    log:
        os.path.join(AN_LOGS, "batch_hmm_gather_{batch_group}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        chunk_dir = os.path.join(EXT_RES, "batch_{batch_group}_combined_chunks")
    shell:
        r"""
        {{
        set -euo pipefail
        echo "# hmmsearch" > {output.unified_raw}
        if [ -n "{input.unified}" ]; then
            cat {input.unified} | grep -v '^#' >> {output.unified_raw} || true
        fi
        # Clean up chunk directory to save space
        rm -rf "{params.chunk_dir}"
        }} 2>&1 | tee {log}
        """


if DEDUP_ENABLED:
    rule ext_annotation_batch_hmm_expand:
        input:
            raw = os.path.join(EXT_RES, "batch_{batch_group}_combined_hmm.tbl"),
            mapping = batch_dedup_mapping,
            database = os.path.join(EXT_RES, "batch_{batch_group}_protein_dedup/catalog.sqlite"),
            cluster_mapping = batch_cluster_mapping
        output:
            expanded = os.path.join(EXT_RES, "batch_{batch_group}_combined_hmm_expanded.tbl")
        log:
            os.path.join(AN_LOGS, "batch_hmm_expand_{batch_group}.log")
        threads: 1
        resources:
            runtime="7d",
            mem_mb=lambda wildcards, threads: get_mem_mb(threads)
        params:
            script = os.path.join(WORKFLOW_ROOT, "scripts/expand_dedup_hmm.py"),
            python_env = an_get_env("python"),
            cluster_mapping_flag = lambda wildcards: f"--cluster-mapping {shlex.quote(batch_cluster_mapping(wildcards))}" if CLUSTER_99_ENABLED else ""
        shell:
            r"""
            {{
                set -euo pipefail
                PYTHON_BIN="{params.python_env}/bin/python3"
                if [ ! -x "$PYTHON_BIN" ]; then PYTHON_BIN="python3"; fi
                "$PYTHON_BIN" {params.script} \
                    --input {input.raw} \
                    --database {input.database} \
                    {params.cluster_mapping_flag} \
                    --output {output.expanded}
            }} 2>&1 | tee {log}
            """


# 5. Split tbl back to MAGs, run resolve/aggregation in parallel (Checkpoint)
checkpoint ext_annotation_batch_hmm_resolve_agg:
    input:
        tbl = batch_hmm_resolved_tbl,
        ko_list = lambda wildcards: EXT_GENOMES.get(wildcards.batch_group, {}).get("ko_list", "")
    output:
        sentinel = os.path.join(EXT_RES, "{batch_group}_batch_hmm_completed.txt")
    log:
        os.path.join(AN_LOGS, "batch_hmm_resolve_agg_{batch_group}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        genomes_dirs = lambda wildcards: f"--genomes-dir {an_full_path(EXT_GENOMES.get(wildcards.batch_group, {}).get('path'))}",
        output_dir = EXT_RES,
        batch_name = lambda wildcards: wildcards.batch_group,
        conda_env_dir = AN_ENV_DIR,
        bakta_env = AN_C.get("conda_envs", {}).get("bakta", "bakta_1.12.0"),
        biopython_env = AN_C.get("conda_envs", {}).get("biopython", "biopython_1.85"),
        run_agg_flag = lambda wildcards: "--run-aggregation" if str(config.get("annotation", {}).get("steps", {}).get("resolve_hmm", "True")).lower() == "true" else "",
        target_profiles_flag = lambda wildcards: f"--target-profiles {shlex.quote(get_target_profiles_for_genome(wildcards.batch_group))}" if get_target_profiles_for_genome(wildcards.batch_group) else "",
        negative_map_flag = lambda wildcards: f"--negative-map {shlex.quote(AN_C.get('negative_map', ''))}" if AN_C.get('negative_map', '') else "",
        competition_fraction_flag = lambda wildcards: f"--competition-fraction {AN_C.get('competition_fraction', 0.10)}",
        overwrite_flag = lambda wildcards: "--force-overwrite" if str(EXT_GENOMES.get(wildcards.batch_group, {}).get("force_overwrite", config.get("force_overwrite_hmm", AN_C.get("force_overwrite_hmm", "False")))).lower() == "true" else ""
    shell:
        r"""
        {{
            set -euo pipefail
            BIOPYTHON_ENV="{params.biopython_env}"
            if [ -n "{params.conda_env_dir}" ]; then
                BIOPYTHON_ENV="{params.conda_env_dir}/{params.biopython_env}"
            fi
            
            # Direct python binary execution to bypass conda run lock issues on HPC
            if [ -f "$BIOPYTHON_ENV/bin/python3" ]; then
                echo "Running directly using python binary: $BIOPYTHON_ENV/bin/python3"
                "$BIOPYTHON_ENV/bin/python3" scripts/batch_prodigal_hmm.py \
                    {params.genomes_dirs} \
                    --output-dir {params.output_dir} \
                    --batch-name {params.batch_name} \
                    --hmm-db none \
                    --ko-list {input.ko_list} \
                    --conda-env-dir {params.conda_env_dir} \
                    --bakta-env {params.bakta_env} \
                    --biopython-env {params.biopython_env} \
                    --threads {threads} \
                    {params.run_agg_flag} \
                    {params.target_profiles_flag} \
                    {params.negative_map_flag} \
                    {params.competition_fraction_flag} \
                    {params.overwrite_flag} \
                    --only-resolve-aggregate \
                    --combined-tbl {input.tbl}
            else
                echo "Falling back to conda run..."
                conda run -p "$BIOPYTHON_ENV" python3 scripts/batch_prodigal_hmm.py \
                    {params.genomes_dirs} \
                    --output-dir {params.output_dir} \
                    --batch-name {params.batch_name} \
                    --hmm-db none \
                    --ko-list {input.ko_list} \
                    --conda-env-dir {params.conda_env_dir} \
                    --bakta-env {params.bakta_env} \
                    --biopython-env {params.biopython_env} \
                    --threads {threads} \
                    {params.run_agg_flag} \
                    {params.target_profiles_flag} \
                    {params.negative_map_flag} \
                    {params.competition_fraction_flag} \
                    {params.overwrite_flag} \
                    --only-resolve-aggregate \
                    --combined-tbl {input.tbl}
            fi
            touch {output.sentinel}
        }} 2>&1 | tee {log}
        """
