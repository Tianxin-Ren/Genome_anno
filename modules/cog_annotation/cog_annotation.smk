# ==============================================================================
# COG2024-Compatible Footprint-Aware Annotation Module
# ==============================================================================
import os
import glob
import re
import sys

# 1. Global Path Logic (Namespaced)
COG_C = config.get("cog_annotation", {})
COG_ROOT = config.get("project_root", "") or config.get("paths", {}).get("project_root", "")
if COG_ROOT:
    COG_ROOT = os.path.expanduser(os.path.expandvars(str(COG_ROOT)))
    if not os.path.isabs(COG_ROOT):
        COG_ROOT = os.path.abspath(COG_ROOT)
else:
    COG_ROOT = os.getcwd()


def cog_full_path(path):
    if not path:
        return ""
    path_str = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path_str):
        return os.path.normpath(path_str)
    if COG_ROOT:
        return os.path.normpath(os.path.join(COG_ROOT, path_str))
    return os.path.normpath(path_str)


def get_active_samples_cog():
    cli_val = config.get("sample") or config.get("sample_name") or COG_C.get("sample") or COG_C.get("sample_name")
    if isinstance(cli_val, str) and cli_val.strip() and cli_val.strip().lower() not in ["all", "*"]:
        return [cli_val.strip()]
    for s_candidate in [config.get("samples"), config.get("sample"), COG_C.get("samples"), COG_C.get("sample")]:
        if isinstance(s_candidate, dict) and s_candidate:
            return list(s_candidate.keys())
        elif isinstance(s_candidate, list) and s_candidate:
            return [str(x).strip() for x in s_candidate if str(x).strip()]
        elif isinstance(s_candidate, str) and s_candidate.strip() and s_candidate.strip().lower() not in ["all", "*"]:
            return [s_candidate.strip()]
    return ["default"]


ACTIVE_SAMPLES_COG = get_active_samples_cog()
COG_BASE_RES = cog_full_path(COG_C.get("results_dir", "results/mag_annotation"))
COG_BASE_LOGS = cog_full_path(COG_C.get("logs_dir", "logs/mag_annotation"))
SCRIPTS_COG_DIR = os.path.join(WORKFLOW_ROOT, "modules", "cog_annotation", "scripts")


def cog_env(name):
    env_name = config.get("environments", {}).get(name) or COG_C.get("conda_envs", {}).get(name)
    if not env_name:
        env_name = name
    env_dir = config.get("conda_env_dir", "") or config.get("paths", {}).get("conda_envs", "")
    if not env_dir or os.path.isabs(env_name):
        return env_name
    candidate = os.path.join(env_dir, env_name)
    if os.path.exists(candidate):
        return candidate
    if os.path.exists(env_dir):
        try:
            for entry in os.listdir(env_dir):
                if entry == env_name or entry.startswith(f"{env_name}_"):
                    return os.path.join(env_dir, entry)
        except Exception:
            pass
    return candidate


def _cluster_cog_db_fallback(item):
    user = os.environ.get("USER", "").lower()
    cwd = os.getcwd().lower()
    is_gses0 = any(k in user or k in cwd for k in ("schhui", "gses0"))
    if is_gses0:
        mapping = {
            "cogle_db": "/dssg/home/acct-schhui/gses0/Database/COG_2024/CDD/Cog",
            "cddid": "/dssg/home/acct-schhui/gses0/Database/COG_2024/CDD/cddid.tbl",
            "cog_def": "/dssg/home/acct-schhui/gses0/Database/COG_2024/cog-24.def.tab",
            "cog_fun": "/dssg/home/acct-schhui/gses0/Database/COG_2024/cog-24.fun.tab",
            "bitscore_specific": "/dssg/home/acct-schhui/gses0/Database/COG_2024/CDD/bitscore_specific.txt",
        }
    else:
        mapping = {
            "cogle_db": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/CDD/Cog",
            "cddid": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/CDD/cddid.tbl",
            "cog_def": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/cog-24.def.tab",
            "cog_fun": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/cog-24.fun.tab",
            "bitscore_specific": "/dssg/home/acct-clswfp/clswfp-icdli1/Database/COG_2024/CDD/bitscore_specific.txt",
        }
    return mapping.get(item, "")


def get_cogle_db():
    val = (
        COG_C.get("databases", {}).get("cogle_db") or
        config.get("databases", {}).get("cogle_db") or
        config.get("databases", {}).get("cdd_cogle") or
        config.get("databases", {}).get("cdd") or
        _cluster_cog_db_fallback("cogle_db")
    )
    return cog_full_path(val)


def get_cddid_tbl():
    val = (
        COG_C.get("databases", {}).get("cddid") or
        config.get("databases", {}).get("cddid") or
        config.get("databases", {}).get("cddid_tbl") or
        _cluster_cog_db_fallback("cddid")
    )
    return cog_full_path(val)


def get_cog_def_tbl():
    val = (
        COG_C.get("databases", {}).get("cog_def") or
        config.get("databases", {}).get("cog_def") or
        config.get("databases", {}).get("cog2024_def") or
        _cluster_cog_db_fallback("cog_def")
    )
    return cog_full_path(val)


def get_cog_fun_tbl():
    val = (
        COG_C.get("databases", {}).get("cog_fun") or
        config.get("databases", {}).get("cog_fun") or
        config.get("databases", {}).get("cog2024_fun") or
        _cluster_cog_db_fallback("cog_fun")
    )
    return cog_full_path(val)


def get_bitscore_specific_txt():
    val = (
        COG_C.get("databases", {}).get("bitscore_specific") or
        config.get("databases", {}).get("bitscore_specific") or
        config.get("databases", {}).get("cdd_bitscore_specific") or
        _cluster_cog_db_fallback("bitscore_specific")
    )
    return cog_full_path(val)


# 2. Target Definitions
COG_ANNO_TARGETS = []
steps = COG_C.get("steps", {})
do_audit = str(steps.get("audit_provenance", True)).lower() == "true"
do_master = str(steps.get("merge_summary", True)).lower() == "true"
do_matrices = str(steps.get("matrices", True)).lower() == "true"
do_qc = str(steps.get("qc", True)).lower() == "true"
do_dense = str(steps.get("export_dense", False)).lower() == "true"

for s in ACTIVE_SAMPLES_COG:
    if do_audit:
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "audit", "cog_database_provenance_report.json"))
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "audit", "cdd_pssm_mapping.tsv"))
    if do_master:
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "master", "protein_COG_master.tsv"))
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "master", "protein_COG_summary.tsv"))
    if do_matrices:
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "matrices", "genome_COG_protein_counts.long.tsv"))
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "matrices", "genome_COG_footprint_counts.long.tsv"))
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "matrices", "genome_COG_category_protein_normalized.tsv"))
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "matrices", "genome_COG_category_footprint_normalized.tsv"))
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "matrices", "genome_COG_category_multilabel.tsv"))
        if do_dense:
            COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "matrices", "genome_COG_protein_counts.dense.tsv"))
    if do_qc:
        COG_ANNO_TARGETS.append(os.path.join(COG_BASE_RES, s, "cog", "validation", "cog_qc_report.json"))


# 3. Input Synergy Helper: Reuse mag_annotation deduplication if present
def get_unique_proteins_input(wildcards):
    s = wildcards.sample
    # 1. Check mag_annotation output first
    mag_anno_res = cog_full_path(config.get("mag_annotation", {}).get("results_dir", "results/mag_annotation"))
    cand_mag = os.path.join(mag_anno_res, s, "proteins", "unique", "unique_proteins.faa.gz")
    if os.path.exists(cand_mag):
        return cand_mag
    # 2. Check cog_annotation local unique
    cand_local = os.path.join(COG_BASE_RES, s, "proteins", "unique", "unique_proteins.faa.gz")
    if os.path.exists(cand_local):
        return cand_local
    return cand_mag


def get_protein_mapping_input(wildcards):
    s = wildcards.sample
    mag_anno_res = cog_full_path(config.get("mag_annotation", {}).get("results_dir", "results/mag_annotation"))
    cand_mag = os.path.join(mag_anno_res, s, "proteins", "unique", "protein_to_unique.tsv.gz")
    if os.path.exists(cand_mag):
        return cand_mag
    cand_local = os.path.join(COG_BASE_RES, s, "proteins", "unique", "protein_to_unique.tsv.gz")
    if os.path.exists(cand_local):
        return cand_local
    return cand_mag


# 4. Rules

# 4.1 Gate 0: Database Provenance & Compatibility Audit
rule cog_annotation_audit_provenance:
    output:
        report = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "cog_database_provenance_report.json"),
        mapping = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "cdd_pssm_mapping.tsv"),
        only_cog = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "only_cog2024.tsv"),
        only_cogle = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "only_cogle.tsv"),
        cat_def = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "cog_category_definition.tsv")
    log:
        os.path.join(COG_BASE_LOGS, "{sample}", "cog", "audit", "audit_provenance.log")
    threads: get_threads(1)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = cog_env("python"),
        script = os.path.join(SCRIPTS_COG_DIR, "audit_database_provenance.py"),
        cddid = lambda wildcards: get_cddid_tbl(),
        cog_def = lambda wildcards: get_cog_def_tbl(),
        cog_fun = lambda wildcards: get_cog_fun_tbl(),
        bitscore_specific = lambda wildcards: get_bitscore_specific_txt(),
        cogle_db = lambda wildcards: get_cogle_db(),
        out_dir = lambda wildcards: os.path.join(COG_BASE_RES, wildcards.sample, "cog", "audit")
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

        if [ -z "{params.cddid}" ] || [ ! -f "{params.cddid}" ]; then
            echo "[cog_annotation ERROR] cddid.tbl not found at: '{params.cddid}'. Please verify databases.cddid in platform/project config." | tee -a "{log}" >&2
            exit 1
        fi
        if [ -z "{params.cog_def}" ] || [ ! -f "{params.cog_def}" ]; then
            echo "[cog_annotation ERROR] cog-24.def.tab not found at: '{params.cog_def}'. Please verify databases.cog_def in platform/project config." | tee -a "{log}" >&2
            exit 1
        fi

        COG_FUN_FLAG=""
        if [ -n "{params.cog_fun}" ] && [ -f "{params.cog_fun}" ]; then
            COG_FUN_FLAG="--cog_fun {params.cog_fun}"
        fi

        BITSCORE_FLAG=""
        if [ -n "{params.bitscore_specific}" ] && [ -f "{params.bitscore_specific}" ]; then
            BITSCORE_FLAG="--bitscore_specific {params.bitscore_specific}"
        fi

        $PY_BIN "{params.script}" \
            --cddid "{params.cddid}" \
            --cog_def "{params.cog_def}" \
            $COG_FUN_FLAG \
            $BITSCORE_FLAG \
            --cogle_db "{params.cogle_db}" \
            -o "{params.out_dir}" 2>&1 | tee "{log}"
        """


# 4.2 Dynamic Residue-Balanced Sharding (Checkpoint)
checkpoint cog_annotation_partition_shards:
    input:
        unique_faa = get_unique_proteins_input
    output:
        shard_list = os.path.join(COG_BASE_RES, "{sample}", "cog", "shards", "shard_list.txt"),
        manifest = os.path.join(COG_BASE_RES, "{sample}", "cog", "shards", "shards_manifest.tsv")
    log:
        os.path.join(COG_BASE_LOGS, "{sample}", "cog", "shards", "partition_shards.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = cog_env("python"),
        script = os.path.join(WORKFLOW_ROOT, "modules", "mag_annotation", "scripts", "shard_fasta.py"),
        shards_dir = lambda wildcards: os.path.join(COG_BASE_RES, wildcards.sample, "cog", "shards"),
        target_aa = COG_C.get("sharding", {}).get("target_aa", 100000000),
        target_proteins = COG_C.get("sharding", {}).get("target_proteins", 400000)
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


# 4.3 Search and Assign Shard (Per Shard RPS-BLAST + Footprint Assignment)
rule cog_annotation_search_and_assign_shard:
    input:
        shard = os.path.join(COG_BASE_RES, "{sample}", "cog", "shards", "{shard}.faa.gz"),
        pssm_mapping = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "cdd_pssm_mapping.tsv")
    output:
        raw_hits = os.path.join(COG_BASE_RES, "{sample}", "cog", "raw_hits", "{shard}.hits.tsv.gz"),
        assigned = os.path.join(COG_BASE_RES, "{sample}", "cog", "assigned", "{shard}.assigned.tsv.gz")
    log:
        os.path.join(COG_BASE_LOGS, "{sample}", "cog", "shards", "{shard}.log")
    threads: get_threads(COG_C.get("rpsblast", {}).get("threads_per_shard", 12))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        python_env = cog_env("python"),
        blast_env = cog_env("blast"),
        assign_script = os.path.join(SCRIPTS_COG_DIR, "assign_cog_footprints.py"),
        cogle_db = lambda wildcards: get_cogle_db(),
        evalue = COG_C.get("rpsblast", {}).get("evalue", 1e-4),
        min_aln_len = COG_C.get("rpsblast", {}).get("min_aln_len", 40),
        min_prof_cov = COG_C.get("rpsblast", {}).get("min_prof_cov", 0.30),
        ambiguity_ratio = COG_C.get("rpsblast", {}).get("ambiguity_ratio_cutoff", 0.15),
        ambiguity_delta = COG_C.get("rpsblast", {}).get("ambiguity_delta_cutoff", 20.0)
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.raw_hits}")" "$(dirname "{output.assigned}")" "$(dirname "{log}")"
        : > "{log}"

        # 1. Resolve RPS-BLAST binary
        RPS_BIN="rpsblast"
        BLAST_ENV="{params.blast_env}"
        if [ -n "$BLAST_ENV" ]; then
            if [ -x "$BLAST_ENV/bin/rpsblast" ]; then
                RPS_BIN="$BLAST_ENV/bin/rpsblast"
            elif [[ "$BLAST_ENV" == *"/"* ]]; then
                RPS_BIN="conda run --no-capture-output -p $BLAST_ENV rpsblast"
            else
                RPS_BIN="conda run --no-capture-output -n $BLAST_ENV rpsblast"
            fi
        fi

        # 2. Auto-detect valid RPS database prefix (supports Cog, Cog_LE, Cog_LE/Cog, etc.)
        RAW_DB_PATH="{params.cogle_db}"
        DB_DIR="$(dirname "$RAW_DB_PATH")"
        RESOLVED_DB=""

        for cand in \
            "$RAW_DB_PATH" \
            "$DB_DIR/Cog" \
            "$DB_DIR/Cog_LE" \
            "$DB_DIR/Cog_LE/Cog" \
            "$RAW_DB_PATH/Cog" \
            "$DB_DIR/Cdd" \
            "$DB_DIR/cdd"; do
            if [ -f "$cand.rps" ] || [ -f "$cand.pal" ] || [ -f "$cand.00.rps" ] || [ -f "$cand.pn" ]; then
                RESOLVED_DB="$cand"
                break
            fi
        done

        if [ -n "$RESOLVED_DB" ]; then
            DB_PATH="$RESOLVED_DB"
            DB_DIR="$(dirname "$DB_PATH")"
            DB_PREFIX="$(basename "$DB_PATH")"
        else
            echo "[cog_annotation ERROR] Cannot find valid RPS-BLAST database (.rps/.pal/.pn) for '{params.cogle_db}'!" | tee -a "{log}" >&2
            echo "[cog_annotation ERROR] Searched candidate directory: $DB_DIR" | tee -a "{log}" >&2
            echo "[cog_annotation ERROR] Directory contents:" | tee -a "{log}" >&2
            ls -la "$DB_DIR" 2>&1 | head -n 30 | tee -a "{log}" >&2
            exit 1
        fi

        echo "[cog_annotation] Resolved RPS database prefix: $DB_PATH" | tee -a "{log}"
        ls -lh "$DB_PATH".* 2>&1 | tee -a "{log}"

        # A .pal prefix is an alias; copying only files named after the alias
        # can omit the RPS volumes listed inside it. Keep aliases on shared storage.
        IS_RPS_ALIAS=0
        if [ -f "$DB_PATH.pal" ]; then
            IS_RPS_ALIAS=1
            echo "[cog_annotation] Keeping alias database on shared storage: $DB_PATH.pal" | tee -a "{log}"
        fi

        # Stage direct RPS database files to local SSD scratch if available.
        if [ "$IS_RPS_ALIAS" -eq 0 ] && [ -n "${{SLURM_TMPDIR:-}}" ] && [ -d "${{SLURM_TMPDIR:-}}" ]; then
            LOCAL_DIR="${{SLURM_TMPDIR}}/cog_db"
            mkdir -p "$LOCAL_DIR"
            (
                flock -x 200
                if [ ! -f "$LOCAL_DIR/.synced" ]; then
                    if ! cp "$DB_DIR"/$DB_PREFIX.* "$LOCAL_DIR/"; then
                        echo "[cog_annotation ERROR] Failed to stage RPS database files from $DB_DIR for prefix $DB_PREFIX" >&2
                        exit 1
                    fi
                    touch "$LOCAL_DIR/.synced"
                fi
            ) 200>"$LOCAL_DIR/.lock"
            if [ -f "$LOCAL_DIR/$DB_PREFIX.pal" ] || [ -f "$LOCAL_DIR/$DB_PREFIX.rps" ] || [ -f "$LOCAL_DIR/$DB_PREFIX.00.rps" ] || [ -f "$LOCAL_DIR/$DB_PREFIX.pn" ]; then
                DB_PATH="$LOCAL_DIR/$DB_PREFIX"
                DB_DIR="$LOCAL_DIR"
                echo "[cog_annotation] Using staged RPS database prefix: $DB_PATH" | tee -a "{log}"
            else
                echo "[cog_annotation ERROR] Staged RPS database prefix is incomplete: $LOCAL_DIR/$DB_PREFIX" | tee -a "{log}" >&2
                ls -lh "$LOCAL_DIR" 2>&1 | tee -a "{log}" >&2
                exit 1
            fi
        fi

        # Alias DBLIST entries may be relative to their source directory.
        export BLASTDB="$DB_DIR${{BLASTDB:+:$BLASTDB}}"

        # 3. Decompress query shard and execute RPS-BLAST with exact slen output
        TMP_DIR="${{SLURM_TMPDIR:-/tmp}}"
        TMP_RAW="$TMP_DIR/{wildcards.shard}.raw.$$.tsv"
        TMP_QUERY="$TMP_DIR/{wildcards.shard}.$$.faa"
        trap 'rm -f "$TMP_RAW" "$TMP_QUERY"' EXIT ERR

        if [[ "{input.shard}" == *.gz ]]; then
            gzip -dc "{input.shard}" > "$TMP_QUERY"
        else
            cp "{input.shard}" "$TMP_QUERY"
        fi

        $RPS_BIN \
            -db "$DB_PATH" \
            -query "$TMP_QUERY" \
            -out "$TMP_RAW" \
            -evalue {params.evalue} \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
            -num_threads {threads} \
            -seg no \
            -comp_based_stats 1 2>&1 | tee -a "{log}"

        # Compress raw hits
        gzip -c "$TMP_RAW" > "{output.raw_hits}"

        # 4. Execute Footprint Assignment & Arbitration
        PY_BIN="python3"
        P_ENV="{params.python_env}"
        if [ -n "$P_ENV" ]; then
            if [ -x "$P_ENV/bin/python3" ]; then
                PY_BIN="$P_ENV/bin/python3"
            elif [ -x "$P_ENV/bin/python" ]; then
                PY_BIN="$P_ENV/bin/python"
            elif [[ "$P_ENV" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $P_ENV python3"
            else
                PY_BIN="conda run --no-capture-output -n $P_ENV python3"
            fi
        fi

        $PY_BIN "{params.assign_script}" \
            -i "$TMP_RAW" \
            -m "{input.pssm_mapping}" \
            -o "{output.assigned}" \
            --evalue {params.evalue} \
            --min_aln_len {params.min_aln_len} \
            --min_prof_cov {params.min_prof_cov} \
            --ambiguity_ratio_cutoff {params.ambiguity_ratio} \
            --ambiguity_delta_cutoff {params.ambiguity_delta} 2>&1 | tee -a "{log}"

        rm -f "$TMP_RAW"
        trap - EXIT ERR
        """


# 4.4 Dynamic Aggregation of Assigned Shards
def get_all_assigned_shards(wildcards):
    shard_list_file = checkpoints.cog_annotation_partition_shards.get(sample=wildcards.sample).output.shard_list
    shards = []
    if os.path.exists(shard_list_file):
        with open(shard_list_file, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s:
                    shards.append(s)
    return [os.path.join(COG_BASE_RES, wildcards.sample, "cog", "assigned", f"{s}.assigned.tsv.gz") for s in shards]


# 4.5 Expand Unique Assignments to Full Genome Master Tables
rule cog_annotation_expand_master:
    input:
        assigned_shards = get_all_assigned_shards,
        mapping = get_protein_mapping_input,
        provenance = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "cog_database_provenance_report.json")
    output:
        master_tsv = os.path.join(COG_BASE_RES, "{sample}", "cog", "master", "protein_COG_master.tsv"),
        summary_tsv = os.path.join(COG_BASE_RES, "{sample}", "cog", "master", "protein_COG_summary.tsv")
    log:
        os.path.join(COG_BASE_LOGS, "{sample}", "cog", "master", "expand_master.log")
    threads: get_threads(8)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = cog_env("python"),
        script = os.path.join(SCRIPTS_COG_DIR, "expand_dedup_cog.py"),
        assigned_dir = lambda wildcards: os.path.join(COG_BASE_RES, wildcards.sample, "cog", "assigned"),
        out_dir = lambda wildcards: os.path.join(COG_BASE_RES, wildcards.sample, "cog", "master")
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
            -s "{params.assigned_dir}" \
            -m "{input.mapping}" \
            -p "{input.provenance}" \
            -o "{params.out_dir}" 2>&1 | tee "{log}"
        """


# 4.6 Comparative Matrices Generation
rule cog_annotation_generate_matrices:
    input:
        master_tsv = os.path.join(COG_BASE_RES, "{sample}", "cog", "master", "protein_COG_master.tsv"),
        cat_def = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "cog_category_definition.tsv")
    output:
        protein_counts = os.path.join(COG_BASE_RES, "{sample}", "cog", "matrices", "genome_COG_protein_counts.long.tsv"),
        footprint_counts = os.path.join(COG_BASE_RES, "{sample}", "cog", "matrices", "genome_COG_footprint_counts.long.tsv"),
        protein_norm = os.path.join(COG_BASE_RES, "{sample}", "cog", "matrices", "genome_COG_category_protein_normalized.tsv"),
        footprint_norm = os.path.join(COG_BASE_RES, "{sample}", "cog", "matrices", "genome_COG_category_footprint_normalized.tsv"),
        multilabel = os.path.join(COG_BASE_RES, "{sample}", "cog", "matrices", "genome_COG_category_multilabel.tsv"),
        cat_def = os.path.join(COG_BASE_RES, "{sample}", "cog", "matrices", "cog_category_definition.tsv")
    log:
        os.path.join(COG_BASE_LOGS, "{sample}", "cog", "matrices", "generate_matrices.log")
    threads: get_threads(4)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = cog_env("python"),
        script = os.path.join(SCRIPTS_COG_DIR, "generate_cog_matrices.py"),
        out_dir = lambda wildcards: os.path.join(COG_BASE_RES, wildcards.sample, "cog", "matrices"),
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

        $PY_BIN "{params.script}" \
            -m "{input.master_tsv}" \
            --category_def "{input.cat_def}" \
            -o "{params.out_dir}" \
            {params.export_dense_flag} 2>&1 | tee "{log}"
        """


# 4.7 6-Gate Quality Control Validator
rule cog_annotation_validate_qc:
    input:
        master_tsv = os.path.join(COG_BASE_RES, "{sample}", "cog", "master", "protein_COG_master.tsv"),
        summary_tsv = os.path.join(COG_BASE_RES, "{sample}", "cog", "master", "protein_COG_summary.tsv"),
        mapping = get_protein_mapping_input,
        provenance = os.path.join(COG_BASE_RES, "{sample}", "cog", "audit", "cog_database_provenance_report.json"),
        protein_norm = os.path.join(COG_BASE_RES, "{sample}", "cog", "matrices", "genome_COG_category_protein_normalized.tsv")
    output:
        report = os.path.join(COG_BASE_RES, "{sample}", "cog", "validation", "cog_qc_report.json"),
        summary = os.path.join(COG_BASE_RES, "{sample}", "cog", "validation", "cog_qc_summary.tsv")
    log:
        os.path.join(COG_BASE_LOGS, "{sample}", "cog", "validation", "validate_qc.log")
    threads: get_threads(1)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = cog_env("python"),
        script = os.path.join(SCRIPTS_COG_DIR, "validate_cog_qc.py"),
        matrices_dir = lambda wildcards: os.path.join(COG_BASE_RES, wildcards.sample, "cog", "matrices"),
        out_dir = lambda wildcards: os.path.join(COG_BASE_RES, wildcards.sample, "cog", "validation")
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
            -m "{input.master_tsv}" \
            -s "{input.summary_tsv}" \
            -p "{input.mapping}" \
            -r "{input.provenance}" \
            -d "{params.matrices_dir}" \
            -o "{params.out_dir}" 2>&1 | tee "{log}"
        """
