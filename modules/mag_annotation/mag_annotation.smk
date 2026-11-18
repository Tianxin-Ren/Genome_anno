# ==============================================================================
# MAG High-Throughput Functional Annotation Module (PyKofamSearch + Full KOfam)
# ==============================================================================
import os
import glob
import json
import re
import sys

# 1. Global Path Logic (Namespaced)
MA_C = config.get("mag_annotation", {})
RUN_STAGE = str(config.get("stage", "all")).strip().lower()
if RUN_STAGE not in {"all", "gene_calling", "annotation"}:
    raise ValueError(f"Unsupported mag_annotation stage: {RUN_STAGE}")
MA_ROOT = config.get("project_root", "") or config.get("paths", {}).get("project_root", "")
if MA_ROOT:
    MA_ROOT = os.path.expanduser(os.path.expandvars(str(MA_ROOT)))
    if not os.path.isabs(MA_ROOT):
        MA_ROOT = os.path.abspath(MA_ROOT)
else:
    MA_ROOT = os.getcwd()


def ma_full_path(path):
    if not path:
        return ""
    path_str = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path_str):
        return os.path.normpath(path_str)
    if MA_ROOT:
        return os.path.normpath(os.path.join(MA_ROOT, path_str))
    return os.path.normpath(path_str)


def get_active_samples():
    cli_val = config.get("sample") or config.get("sample_name") or MA_C.get("sample") or MA_C.get("sample_name")
    if isinstance(cli_val, str) and cli_val.strip() and cli_val.strip().lower() not in ["all", "*"]:
        return [cli_val.strip()]
    for s_candidate in [config.get("samples"), config.get("sample"), MA_C.get("samples"), MA_C.get("sample")]:
        if isinstance(s_candidate, dict) and s_candidate:
            return list(s_candidate.keys())
        elif isinstance(s_candidate, list) and s_candidate:
            return [str(x).strip() for x in s_candidate if str(x).strip()]
        elif isinstance(s_candidate, str) and s_candidate.strip() and s_candidate.strip().lower() not in ["all", "*"]:
            return [s_candidate.strip()]
    return ["default"]

ACTIVE_SAMPLES = get_active_samples()
MA_BASE_RES = ma_full_path(MA_C.get("results_dir", "results/mag_annotation"))
MA_BASE_LOGS = ma_full_path(MA_C.get("logs_dir", "logs/mag_annotation"))

SCRIPTS_DIR = os.path.join(WORKFLOW_ROOT, "modules", "mag_annotation", "scripts")


def ma_env(name):
    env_name = config.get("environments", {}).get(name) or MA_C.get("conda_envs", {}).get(name)
    if not env_name:
        env_name = name
    env_dir = config.get("conda_env_dir", "") or config.get("paths", {}).get("conda_envs", "")
    if not env_dir or os.path.isabs(env_name):
        return env_name
    candidate = os.path.join(env_dir, env_name)
    if os.path.exists(candidate):
        return candidate
    # Automatic fallback for common versioned suffixes (e.g. pykofamsearch -> pykofamsearch_26.7.21)
    if os.path.exists(env_dir):
        try:
            for entry in os.listdir(env_dir):
                if entry == env_name or entry.startswith(f"{env_name}_"):
                    return os.path.join(env_dir, entry)
        except Exception:
            pass
    return candidate



# 2. Input Registration & Target Definitions
MAG_ANNO_TARGETS = []
steps = MA_C.get("steps", {})
do_merge_summary = str(steps.get("merge_summary", True)).lower() == "true"
do_validation = str(steps.get("validation", False)).lower() == "true"
do_dense = str(steps.get("export_dense", False)).lower() == "true"
do_postprocess = str(steps.get("postprocess_atp", True)).lower() == "true"
do_correlation = str(steps.get("analyze_atp_correlation", True)).lower() == "true"
_cog_enabled_in_pipeline = config.get("modules", {}).get("cog_annotation", {}).get("enable", False)
do_combine = str(steps.get("combine_kegg_cog", _cog_enabled_in_pipeline)).lower() == "true" and (
    _cog_enabled_in_pipeline or str(steps.get("combine_kegg_cog", "false")).lower() == "true"
)
COG_RES_FOR_MA = ma_full_path(config.get("cog_annotation", {}).get("results_dir", MA_BASE_RES))


def get_ma_strict_hits(wildcards):
    s = wildcards.sample
    base = os.path.join(MA_BASE_RES, s, "merged", "protein_KO_all_strict.tsv")
    if os.path.exists(base):
        return base
    if os.path.exists(base + ".gz"):
        return base + ".gz"
    return base


def get_ma_long_matrix(wildcards):
    s = wildcards.sample
    base = os.path.join(MA_BASE_RES, s, "merged", "genome_KO_counts.long.tsv")
    if os.path.exists(base):
        return base
    if os.path.exists(base + ".gz"):
        return base + ".gz"
    return base


for s in ACTIVE_SAMPLES:
    if RUN_STAGE != "gene_calling" and do_merge_summary:
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "merged", "protein_KO_all_strict.tsv"))
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "merged", "protein_KO_best.tsv"))
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "merged", "genome_KO_counts.long.tsv"))
    if RUN_STAGE != "gene_calling" and do_validation:
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "validation", "validation_report.json"))
    if RUN_STAGE != "gene_calling" and do_dense:
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "merged", "genome_KO_counts.dense.tsv"))
    if RUN_STAGE != "gene_calling" and do_postprocess:
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "postprocess", "genome_atp_synthase_master_table.tsv"))
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "postprocess", "genome_atp_synthase_modules.long.tsv"))
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "postprocess", "genome_atp_synthase_subunits.tsv"))
    if RUN_STAGE != "gene_calling" and do_correlation:
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "postprocess", "atp_gene_correlations_all_kos.tsv"))
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "postprocess", "atp_gene_correlations_candidates.tsv"))
    if RUN_STAGE != "gene_calling" and do_combine:
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "combined", "protein_annotation_master.tsv"))
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "combined", "genome_annotation_summary.tsv"))
    if RUN_STAGE != "gene_calling" and str(steps.get("split_genome_annotations", False)).lower() == "true":
        MAG_ANNO_TARGETS.append(os.path.join(MA_BASE_RES, s, "per_genome", "split_complete.json"))



# 3. Discovery of Genomes or Proteins
SAMPLE_GENOMES = {}
SAMPLE_PROTEIN_DIR = {}
GENE_BATCHES = {}

def get_sample_info(sample_tag):
    if not sample_tag:
        return {}
    for s_dict in [config.get("samples"), config.get("sample"), MA_C.get("samples"), MA_C.get("sample")]:
        if isinstance(s_dict, dict) and sample_tag in s_dict:
            val = s_dict[sample_tag]
            if isinstance(val, dict):
                return val
            elif isinstance(val, str):
                return {"genome_dir": val}
    return {}

for s in ACTIVE_SAMPLES:
    s_info = get_sample_info(s)
    g_dir = s_info.get("genome_dir") or s_info.get("path") or config.get("genome_dir") or MA_C.get("genome_dir", "")
    p_dir = s_info.get("protein_dir") or s_info.get("faa") or config.get("protein_dir") or MA_C.get("protein_dir", "")

    SAMPLE_PROTEIN_DIR[s] = ma_full_path(p_dir) if p_dir else ""
    SAMPLE_GENOMES[s] = {}
    GENE_BATCHES[s] = {}

    # In stage 2, pre-called proteins are authoritative. Avoid scanning the
    # original genome directory, which can contain tens of thousands of files.
    if g_dir and (not p_dir or RUN_STAGE == "gene_calling"):
        abs_gdir = ma_full_path(g_dir)
        if os.path.exists(abs_gdir):
            valid_exts = {".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz", ".fa.gz"}
            for root, _, files in os.walk(abs_gdir):
                for f in files:
                    for ext in valid_exts:
                        if f.endswith(ext) and not f.endswith((".faa", ".faa.gz", ".ffn", ".ffn.gz")):
                            stem = f[:-len(ext)]
                            candidate = os.path.join(root, f)
                            previous = SAMPLE_GENOMES[s].get(stem)
                            if previous and os.path.normpath(previous) != os.path.normpath(candidate):
                                raise ValueError(
                                    f"[mag_annotation ERROR] Duplicate genome ID '{stem}' for sample '{s}': "
                                    f"{previous} and {candidate}"
                                )
                            SAMPLE_GENOMES[s][stem] = candidate
                            break

    # Also check external_genomes in config
    ext_genomes = config.get("external_genomes", {})
    if ext_genomes and isinstance(ext_genomes, dict):
        for name, info in ext_genomes.items():
            if isinstance(info, dict) and "path" in info:
                SAMPLE_GENOMES[s][name] = ma_full_path(info["path"])


# Stage 1 target registration happens after genome discovery so the input
# directory is not expanded when stage 2 uses pre-called proteins.
if RUN_STAGE == "gene_calling":
    for s in ACTIVE_SAMPLES:
        genomes = SAMPLE_GENOMES.get(s, {})
        if not genomes:
            raise ValueError(f"[mag_annotation ERROR] No genomes found for gene-calling stage: {s}")
        batch_size = int(MA_C.get("gene_calling", {}).get("batch_size", 100))
        if batch_size < 1:
            raise ValueError("mag_annotation.gene_calling.batch_size must be >= 1")
        items = sorted(genomes.items())
        for offset in range(0, len(items), batch_size):
            batch_id = f"{offset // batch_size:05d}"
            GENE_BATCHES[s][batch_id] = items[offset:offset + batch_size]
            MAG_ANNO_TARGETS.append(
                os.path.join(MA_BASE_RES, s, "gene_calling_batches", f"{batch_id}.complete.json")
            )



# 4. Rules Chain
# ------------------------------------------------------------------------------

# 4.1 Batched gene calling for very large genome collections. The marker is
# written only after every genome in the batch has produced its FAA/GFF/meta.
rule mag_annotation_call_genes_batch:
    input:
        fna=lambda wildcards: [path for _, path in GENE_BATCHES[wildcards.sample][wildcards.batch]]
    output:
        marker=os.path.join(MA_BASE_RES, "{sample}", "gene_calling_batches", "{batch}.complete.json")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "pyrodigal_batches", "{batch}.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("pyrodigal"),
        script = os.path.join(SCRIPTS_DIR, "run_pyrodigal_batch.py"),
        out_faa = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "proteins", "all"),
        out_gff = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "genes"),
        out_meta = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "metadata"),
        mode = MA_C.get("gene_calling", {}).get("mode", "auto"),
        trans_table = MA_C.get("gene_calling", {}).get("trans_table", 11),
        auto_complete_mode = MA_C.get("gene_calling", {}).get("auto", {}).get("complete_mode", "single"),
        auto_fallback_mode = MA_C.get("gene_calling", {}).get("auto", {}).get("fallback_mode", "meta")
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_faa}" "{params.out_gff}" "{params.out_meta}" "$(dirname "{output.marker}")" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi
        $PY_BIN "{params.script}" \
          --genomes {input.fna} \
          --out-faa "{params.out_faa}" --out-gff "{params.out_gff}" --out-meta "{params.out_meta}" \
          --marker "{output.marker}" --mode "{params.mode}" --trans-table {params.trans_table} \
          --auto-complete-mode "{params.auto_complete_mode}" \
          --auto-fallback-mode "{params.auto_fallback_mode}" 2>&1 | tee "{log}"
        """

# 4.2 Gene Calling (Only if starting from DNA genomes)
rule mag_annotation_call_genes:
    input:
        fna = lambda wildcards: SAMPLE_GENOMES[wildcards.sample][wildcards.genome]
    output:
        faa = os.path.join(MA_BASE_RES, "{sample}", "proteins", "all", "{genome}.faa.gz"),
        gff = os.path.join(MA_BASE_RES, "{sample}", "genes", "{genome}.gff.gz"),
        meta = os.path.join(MA_BASE_RES, "{sample}", "metadata", "{genome}_protein_metadata.tsv.gz")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "pyrodigal", "{genome}.log")
    threads: get_threads(1)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("pyrodigal"),
        script = os.path.join(SCRIPTS_DIR, "run_pyrodigal.py"),
        mode = MA_C.get("gene_calling", {}).get("mode", "meta"),
        trans_table = MA_C.get("gene_calling", {}).get("trans_table", 11),
        auto_complete_mode = MA_C.get("gene_calling", {}).get("auto", {}).get("complete_mode", "single"),
        auto_fallback_mode = MA_C.get("gene_calling", {}).get("auto", {}).get("fallback_mode", "meta")
    shell:
        """
        set -euo pipefail
        ulimit -c 0
        mkdir -p "$(dirname "{output.faa}")" "$(dirname "{output.gff}")" "$(dirname "{output.meta}")" "$(dirname "{log}")"

        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        $PY_BIN "{params.script}" \
            -i "{input.fna}" \
            -g "{wildcards.genome}" \
            --out_faa "{output.faa}" \
            --out_gff "{output.gff}" \
            --out_meta "{output.meta}" \
            -m "{params.mode}" \
            -t {params.trans_table} \
            --auto-complete-mode "{params.auto_complete_mode}" \
            --auto-fallback-mode "{params.auto_fallback_mode}" 2>&1 | tee "{log}"
        """


def get_all_protein_inputs(wildcards):
    s = wildcards.sample
    prot_dir = SAMPLE_PROTEIN_DIR.get(s, "")
    if prot_dir:
        if os.path.isdir(prot_dir):
            expected_dir = os.path.normpath(os.path.join(MA_BASE_RES, s, "proteins", "all"))
            if RUN_STAGE == "annotation" and os.path.normpath(prot_dir) == expected_dir:
                batch_dir = os.path.join(MA_BASE_RES, s, "gene_calling_batches")
                markers = sorted(glob.glob(os.path.join(batch_dir, "*.complete.json")))
                if not markers:
                    raise ValueError(
                        f"[mag_annotation ERROR] No gene-calling completion markers found for sample '{s}': {batch_dir}. "
                        "Finish stage 1 gene calling before starting stage 2."
                    )
                missing = []
                for marker in markers:
                    try:
                        with open(marker, "r", encoding="utf-8") as handle:
                            payload = json.load(handle)
                    except Exception as exc:
                        raise ValueError(f"[mag_annotation ERROR] Invalid gene-calling marker: {marker}: {exc}") from exc
                    for genome in payload.get("genomes", []):
                        faa = os.path.join(prot_dir, f"{genome}.faa.gz")
                        if not os.path.isfile(faa) or os.path.getsize(faa) == 0:
                            missing.append(faa)
                if missing:
                    preview = ", ".join(missing[:3])
                    suffix = f" (+{len(missing) - 3} more)" if len(missing) > 3 else ""
                    raise ValueError(
                        f"[mag_annotation ERROR] Gene-calling markers exist but {len(missing)} FAA files are missing or empty: "
                        f"{preview}{suffix}. Complete stage 1 before starting stage 2."
                    )
            return [prot_dir]
        raise ValueError(
            f"[mag_annotation ERROR] Configured protein_dir does not exist for sample '{s}': {prot_dir}. "
            "Complete stage 1 gene calling before starting stage 2."
        )
    genomes = sorted(SAMPLE_GENOMES.get(s, {}).keys())
    if not genomes:
        s_info = get_sample_info(s)
        g_dir = s_info.get("genome_dir") or s_info.get("path") or ""
        raise ValueError(
            f"\n[mag_annotation ERROR] No input genomes or proteins found for sample '{s}'!\n"
            f"  - Configured genome_dir: {ma_full_path(g_dir) if g_dir else '(not set)'}\n"
            f"  - Configured protein_dir: {prot_dir or '(not set)'}\n"
            f"Please verify that 'genome_dir' exists and contains FASTA files (.fna/.fasta/.fa/.gz) on the cluster."
        )
    return [os.path.join(MA_BASE_RES, s, "proteins", "all", f"{g}.faa.gz") for g in genomes]


# 4.2 Exact Sequence Deduplication & Dual Statistical Provenance
rule mag_annotation_dedup_proteins:
    input:
        prots = get_all_protein_inputs
    output:
        unique_faa = os.path.join(MA_BASE_RES, "{sample}", "proteins", "unique", "unique_proteins.faa.gz"),
        mapping = os.path.join(MA_BASE_RES, "{sample}", "proteins", "unique", "protein_to_unique.tsv.gz"),
        metadata = os.path.join(MA_BASE_RES, "{sample}", "proteins", "unique", "run_metadata.json"),
        stats = os.path.join(MA_BASE_RES, "{sample}", "proteins", "unique", "dedup_stats.tsv")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "dedup", "dedup_proteins.log")
    threads: get_threads(8)
    resources:
        runtime = "7d",
        # The exact-dedup index retains every unique sequence in RAM. Keep
        # this job on the 512 GB partition instead of the profile default.
        mem_mb = 500000,
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("python"),
        script = os.path.join(SCRIPTS_DIR, "dedup_proteins.py"),
        out_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "proteins", "unique"),
        prots_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "proteins", "all"),
        context = lambda wildcards: MA_C.get("proteome_context") or wildcards.sample or "GTDB_r232_archaea",
        first_prot = lambda wildcards, input: str(input.prots[0]) if input.prots else ""
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        LIST_FILE="{params.out_dir}/input_protein_files.txt"
        if [ -n "{params.first_prot}" ] && [ -d "{params.first_prot}" ]; then
            echo "{params.first_prot}" > "$LIST_FILE"
        else
            echo "{params.prots_dir}" > "$LIST_FILE"
        fi

        $PY_BIN "{params.script}" \
            -i "$LIST_FILE" \
            -o "{params.out_dir}" \
            -c "{params.context}" 2>&1 | tee "{log}"
        """


# 4.3 Checkpoint: Complexity-Aware Residue-Balanced Sharding
checkpoint mag_annotation_partition_shards:
    input:
        unique_faa = os.path.join(MA_BASE_RES, "{sample}", "proteins", "unique", "unique_proteins.faa.gz")
    output:
        shard_list = os.path.join(MA_BASE_RES, "{sample}", "shards", "shard_list.txt"),
        manifest = os.path.join(MA_BASE_RES, "{sample}", "shards", "shards_manifest.tsv")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "shards", "partition_shards.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("python"),
        script = os.path.join(SCRIPTS_DIR, "shard_fasta.py"),
        shards_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "shards"),
        target_aa = MA_C.get("sharding", {}).get("target_aa", 120000000),
        target_proteins = MA_C.get("sharding", {}).get("target_proteins", 500000)
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.shards_dir}" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        $PY_BIN "{params.script}" \
            -i "{input.unique_faa}" \
            -o "{params.shards_dir}" \
            --target_aa {params.target_aa} \
            --target_proteins {params.target_proteins} 2>&1 | tee "{log}"
        """


# 4.4 Shard-Level Independent PyKofamSearch Job (1 Shard = 1 Slurm Job, 12 Cores)
rule mag_annotation_search_shard:
    input:
        shard = os.path.join(MA_BASE_RES, "{sample}", "shards", "{shard}.faa.gz"),
        db = ma_full_path(MA_C.get("database_pkl") or config.get("databases", {}).get("database_pkl") or config.get("databases", {}).get("kofam_database_pkl") or config.get("databases", {}).get("kofam_db", "kofam_db/database.pkl"))
    output:
        tsv = os.path.join(MA_BASE_RES, "{sample}", "results", "{shard}.tsv")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "shards", "{shard}.log")
    threads: get_threads(MA_C.get("pykofam", {}).get("threads_per_shard", 12))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("pykofamsearch"),
        prefilter_evalue = MA_C.get("pykofam", {}).get("prefilter_evalue", 10.0),
        all_hits_flag = "--all_hits" if float(MA_C.get("pykofam", {}).get("threshold_factor", 0.80)) < 1.0 else ""
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.tsv}")" "$(dirname "{log}")"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ] && [ ! -d "$ENV_DIR" ] && [ -d "$ENV_DIR"_26.7.21 ]; then
            ENV_DIR="$ENV_DIR"_26.7.21
        fi

        PYKOFAM_CMD="pykofamsearch"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/pykofamsearch" ]; then
                PYKOFAM_CMD="$ENV_DIR/bin/pykofamsearch"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PYKOFAM_CMD="conda run --no-capture-output -p $ENV_DIR pykofamsearch"
            else
                PYKOFAM_CMD="conda run --no-capture-output -n $ENV_DIR pykofamsearch"
            fi
        fi

        DB_PATH="{input.db}"
        if [ -n "${{SLURM_TMPDIR:-}}" ] && [ -d "${{SLURM_TMPDIR:-}}" ]; then
            LOCAL_DB="${{SLURM_TMPDIR}}/$(basename "{input.db}")"
            if [ ! -f "$LOCAL_DB" ]; then
                TMP_LOCAL_DB="${{LOCAL_DB}}.tmp.$$"
                cp "{input.db}" "$TMP_LOCAL_DB"
                mv -n "$TMP_LOCAL_DB" "$LOCAL_DB" 2>/dev/null || rm -f "$TMP_LOCAL_DB"
            fi
            if [ -f "$LOCAL_DB" ]; then
                DB_PATH="$LOCAL_DB"
            fi
        fi

        TMP_OUT="{output.tsv}.tmp.$$"
        trap 'rm -f "$TMP_OUT"' EXIT ERR

        $PYKOFAM_CMD \
            -i "{input.shard}" \
            -o "$TMP_OUT" \
            -b "$DB_PATH" \
            -p {threads} \
            --evalue {params.prefilter_evalue} \
            {params.all_hits_flag} 2>&1 | tee "{log}"

        mv "$TMP_OUT" "{output.tsv}"
        trap - EXIT ERR
        """


# 4.5 Dynamic Aggregation of Shards
def get_all_shard_tsvs(wildcards):
    shard_list_file = checkpoints.mag_annotation_partition_shards.get(sample=wildcards.sample).output.shard_list
    shards = []
    if os.path.exists(shard_list_file):
        with open(shard_list_file, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s:
                    shards.append(s)
    return [os.path.join(MA_BASE_RES, wildcards.sample, "results", f"{s}.tsv") for s in shards]


# 4.6 Merge Primary Artifacts & Secondary Matrix Aggregation
rule mag_annotation_merge_summary:
    input:
        shard_tsvs = get_all_shard_tsvs,
        mapping = os.path.join(MA_BASE_RES, "{sample}", "proteins", "unique", "protein_to_unique.tsv.gz")
    output:
        strict_hits = os.path.join(MA_BASE_RES, "{sample}", "merged", "protein_KO_all_strict.tsv"),
        best_hits = os.path.join(MA_BASE_RES, "{sample}", "merged", "protein_KO_best.tsv"),
        long_matrix = os.path.join(MA_BASE_RES, "{sample}", "merged", "genome_KO_counts.long.tsv")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "summary", "merge_summary.log")
    threads: get_threads(MA_C.get("merge_summary", {}).get("threads", 8))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("python"),
        script = os.path.join(SCRIPTS_DIR, "merge_kofam_results.py"),
        shards_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "results"),
        out_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "merged"),
        ko_count_strategy = MA_C.get("ko_count_strategy", "best_hit"),
        threshold_factor = MA_C.get("pykofam", {}).get("threshold_factor", 0.80),
        evalue_fallback_cutoff = MA_C.get("pykofam", {}).get("evalue_fallback_cutoff", 1e-3),
        ko_list = ma_full_path(MA_C.get("ko_list") or config.get("databases", {}).get("ko_list") or MA_C.get("database_pkl") or config.get("databases", {}).get("database_pkl") or config.get("databases", {}).get("kofam_database_pkl") or config.get("databases", {}).get("kofam_db", "")),
        export_dense_flag = "--export_dense" if do_dense else ""
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        KO_LIST_FLAG=""
        if [ -n "{params.ko_list}" ] && [ -f "{params.ko_list}" ]; then
            KO_LIST_FLAG="--ko_list {params.ko_list}"
        fi

        $PY_BIN "{params.script}" \
            -s "{params.shards_dir}" \
            -m "{input.mapping}" \
            -o "{params.out_dir}" \
            -p {threads} \
            --threshold_factor {params.threshold_factor} \
            --evalue_fallback_cutoff {params.evalue_fallback_cutoff} \
            --ko_count_strategy "{params.ko_count_strategy}" \
            $KO_LIST_FLAG \
            {params.export_dense_flag} 2>&1 | tee "{log}"
        """



# 4.7 Optional Scientific Validation Gate (Level 1, 2, 3)
rule mag_annotation_scientific_validation:
    input:
        strict_hits = get_ma_strict_hits,
        long_matrix = get_ma_long_matrix
    output:
        report = os.path.join(MA_BASE_RES, "{sample}", "validation", "validation_report.json")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "validation", "validate_invariance.log")
    threads: get_threads(1)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("python"),
        script = os.path.join(SCRIPTS_DIR, "validate_invariance.py")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.report}")" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        $PY_BIN "{params.script}" \
            -o "{output.report}" \
            --annotation_tsv "{input.long_matrix}" 2>&1 | tee "{log}"
        """


# 4.8 Post-Annotation ATP Synthase Structural Gating & Metadata Integration
rule mag_annotation_postprocess_atp_synthase:
    input:
        strict_hits = get_ma_strict_hits
    output:
        master = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "genome_atp_synthase_master_table.tsv"),
        modules_long = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "genome_atp_synthase_modules.long.tsv"),
        subunits = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "genome_atp_synthase_subunits.tsv"),
        pruned_meta = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "pruned_metadata.tsv")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "postprocess", "postprocess_atp_synthase.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("python"),
        script = os.path.join(SCRIPTS_DIR, "postprocess_atp_synthase.py"),
        metadata = ma_full_path(MA_C.get("postprocess", {}).get("metadata_tsv", "")),
        proteins_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "proteins", "all"),
        out_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "postprocess"),
        operon_gap = MA_C.get("postprocess", {}).get("operon_max_gene_gap", 3),
        w_borderline = MA_C.get("postprocess", {}).get("weight_borderline", 0.5),
        w_supported = MA_C.get("postprocess", {}).get("weight_supported_borderline", 0.75),
        hq_comp = MA_C.get("postprocess", {}).get("high_quality_genome_completeness", 90.0),
        hq_contam = MA_C.get("postprocess", {}).get("high_quality_contamination_limit", 5.0),
        lq_comp = MA_C.get("postprocess", {}).get("low_quality_genome_completeness", 70.0)
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        META_FLAG=""
        if [ -n "{params.metadata}" ] && [ -f "{params.metadata}" ]; then
            META_FLAG="-m {params.metadata}"
        fi

        PROT_FLAG=""
        if [ -d "{params.proteins_dir}" ]; then
            PROT_FLAG="--proteins_dir {params.proteins_dir}"
        fi

        $PY_BIN "{params.script}" \
            -k "{input.strict_hits}" \
            $META_FLAG \
            $PROT_FLAG \
            -o "{params.out_dir}" \
            --operon_max_gene_gap {params.operon_gap} \
            --weight_borderline {params.w_borderline} \
            --weight_supported_borderline {params.w_supported} \
            --high_quality_genome_completeness {params.hq_comp} \
            --high_quality_contamination_limit {params.hq_contam} \
            --low_quality_genome_completeness {params.lq_comp} 2>&1 | tee "{log}"
        """


# 4.9 Whole-Genome KO & ATP Synthase Co-occurrence / Correlation Analysis
rule mag_annotation_analyze_atp_correlation:
    input:
        master = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "genome_atp_synthase_master_table.tsv"),
        ko_counts = get_ma_long_matrix
    output:
        all_kos = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "atp_gene_correlations_all_kos.tsv"),
        pos_ctrls = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "atp_gene_correlations_positive_controls.tsv"),
        candidates = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "atp_gene_correlations_candidates.tsv"),
        top_associated = os.path.join(MA_BASE_RES, "{sample}", "postprocess", "atp_gene_correlations_top_associated.tsv")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "postprocess", "analyze_atp_correlation.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("python"),
        script = os.path.join(SCRIPTS_DIR, "analyze_atp_correlation.py"),
        out_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "postprocess"),
        min_prevalence = MA_C.get("correlation_analysis", {}).get("min_ko_prevalence_count", 5),
        fdr_alpha = MA_C.get("correlation_analysis", {}).get("fdr_alpha", 0.05),
        delta_prev = MA_C.get("correlation_analysis", {}).get("delta_prevalence_cutoff", 0.20),
        min_stratum = MA_C.get("correlation_analysis", {}).get("min_stratum_group_size", 3),
        sampling_mode = MA_C.get("correlation_analysis", {}).get("sampling_mode_primary", "all")
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        $PY_BIN "{params.script}" \
            -m "{input.master}" \
            -k "{input.ko_counts}" \
            -o "{params.out_dir}" \
            --min_ko_prevalence_count {params.min_prevalence} \
            --fdr_alpha {params.fdr_alpha} \
            --delta_prevalence_cutoff {params.delta_prev} \
            --min_stratum_group_size {params.min_stratum} \
            --sampling_mode_primary "{params.sampling_mode}" 2>&1 | tee "{log}"
        """


# 4.10 Unified KEGG + COG Dual Functional Annotation Engine
rule mag_annotation_combine_kegg_cog:
    input:
        protein_mapping = os.path.join(MA_BASE_RES, "{sample}", "proteins", "unique", "protein_to_unique.tsv.gz"),
        ko_best = os.path.join(MA_BASE_RES, "{sample}", "merged", "protein_KO_best.tsv"),
        cog_summary = os.path.join(COG_RES_FOR_MA, "{sample}", "cog", "master", "protein_COG_summary.tsv"),
        cat_def = os.path.join(COG_RES_FOR_MA, "{sample}", "cog", "audit", "cog_category_definition.tsv"),
        pssm_mapping = os.path.join(COG_RES_FOR_MA, "{sample}", "cog", "audit", "cdd_pssm_mapping.tsv")
    output:
        master = os.path.join(MA_BASE_RES, "{sample}", "combined", "protein_annotation_master.tsv"),
        summary = os.path.join(MA_BASE_RES, "{sample}", "combined", "genome_annotation_summary.tsv")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "combined", "combine_kegg_cog.log")
    threads: get_threads(4)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = ma_env("python"),
        script = os.path.join(SCRIPTS_DIR, "combine_kegg_cog_annotations.py"),
        out_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "combined")
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        PY_BIN="python3"
        ENV_DIR="{params.env}"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi

        $PY_BIN "{params.script}" \
            -p "{input.protein_mapping}" \
            -k "{input.ko_best}" \
            -c "{input.cog_summary}" \
            --cat_def "{input.cat_def}" \
            -m "{input.pssm_mapping}" \
            -o "{params.out_dir}" 2>&1 | tee "{log}"
        """


# 4.11 Split the all-hit KEGG table into one file per genome.
rule mag_annotation_split_per_genome:
    input:
        kegg = os.path.join(MA_BASE_RES, "{sample}", "merged", "protein_KO_all_strict.tsv")
    output:
        marker = os.path.join(MA_BASE_RES, "{sample}", "per_genome", "split_complete.json")
    log:
        os.path.join(MA_BASE_LOGS, "{sample}", "per_genome", "split_annotations.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        script = os.path.join(SCRIPTS_DIR, "split_annotations_by_genome.py"),
        out_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "per_genome"),
        proteins_dir = lambda wildcards: os.path.join(MA_BASE_RES, wildcards.sample, "proteins", "all")
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        python3 "{params.script}" \
            --kegg "{input.kegg}" \
            --proteins_dir "{params.proteins_dir}" \
            --out_dir "{params.out_dir}" --marker "{output.marker}" 2>&1 | tee "{log}"
        """


