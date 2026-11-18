# ==============================================================================
# MAG Ecological & Metabolic Traits Module (dRep, dbCAN, CAMPER, gRodon2)
# ==============================================================================
import os
import glob
import re
import sys

# 1. Global Path Logic (Namespaced)
MT_C = config.get("mag_traits", {})
MT_ROOT = config.get("project_root", "") or config.get("paths", {}).get("project_root", "")
if MT_ROOT:
    MT_ROOT = os.path.expanduser(os.path.expandvars(str(MT_ROOT)))
    if not os.path.isabs(MT_ROOT):
        MT_ROOT = os.path.abspath(MT_ROOT)
else:
    MT_ROOT = os.getcwd()


def mt_full_path(path):
    if not path:
        return ""
    path_str = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path_str):
        return os.path.normpath(path_str)
    if MT_ROOT:
        return os.path.normpath(os.path.join(MT_ROOT, path_str))
    return os.path.normpath(path_str)


def get_active_samples():
    cli_val = config.get("sample") or config.get("sample_name") or MT_C.get("sample") or MT_C.get("sample_name")
    if isinstance(cli_val, str) and cli_val.strip() and cli_val.strip().lower() not in ["all", "*"]:
        return [cli_val.strip()]
    for s_candidate in [config.get("samples"), config.get("sample"), MT_C.get("samples"), MT_C.get("sample")]:
        if isinstance(s_candidate, dict) and s_candidate:
            return list(s_candidate.keys())
        elif isinstance(s_candidate, list) and s_candidate:
            return [str(x).strip() for x in s_candidate if str(x).strip()]
        elif isinstance(s_candidate, str) and s_candidate.strip() and s_candidate.strip().lower() not in ["all", "*"]:
            return [s_candidate.strip()]
    # Fallback to directory basename of mag_dir if available, else "default"
    m_dir = config.get("mag_dir") or MT_C.get("mag_dir", "")
    if m_dir:
        base = os.path.basename(os.path.normpath(str(m_dir).strip()))
        if base:
            return [base]
    return ["default"]


ACTIVE_SAMPLES = get_active_samples()
MT_BASE_RES = mt_full_path(MT_C.get("results_dir", "results/mag_traits"))
MT_BASE_LOGS = mt_full_path(MT_C.get("logs_dir", "logs/mag_traits"))

# Retain MT_RES and MT_LOGS for backward compatibility if referenced elsewhere
MT_RES = MT_BASE_RES
MT_LOGS = MT_BASE_LOGS


def get_sample_info(sample_tag):
    if not sample_tag:
        return {}
    for s_dict in [config.get("samples"), config.get("sample"), MT_C.get("samples"), MT_C.get("sample")]:
        if isinstance(s_dict, dict) and sample_tag in s_dict:
            val = s_dict[sample_tag]
            if isinstance(val, dict):
                return val
            elif isinstance(val, str):
                return {"mag_dir": val}
    return {}


def mt_env(name):
    env_map = MT_C.get("conda_envs", {})
    env_name = env_map.get(name) or config.get("environments", {}).get(name)
    if not env_name and name == "grodon2":
        env_name = (
            env_map.get("grodon")
            or env_map.get("phydon")
            or config.get("environments", {}).get("grodon")
            or config.get("environments", {}).get("phydon")
        )
    if not env_name and name == "genomespot":
        env_name = (
            env_map.get("genome_spot")
            or config.get("environments", {}).get("genome_spot")
            or "GenomeSPOT_1.0.1"
        )
    if not env_name:
        env_name = name
    env_dir = config.get("conda_env_dir", "")
    if not env_dir or os.path.isabs(env_name):
        return env_name
    return os.path.join(env_dir, env_name)


# 2. Config & Parameter Toggles
# ------------------------------------------------------------------------------
drep_cfg = MT_C.get("drep", {})
two_rounds_drep = str(drep_cfg.get("two_rounds", True)).lower() == "true"
r1_ani = drep_cfg.get("round1_ani", 0.99)
r2_ani = drep_cfg.get("round2_ani", 0.95)
single_ani = drep_cfg.get("ani_threshold", 0.95)
drep_comp = drep_cfg.get("min_completeness", drep_cfg.get("completeness", drep_cfg.get("comp", 50.0)))
drep_cont = drep_cfg.get("max_contamination", drep_cfg.get("contamination", drep_cfg.get("cont", 10.0)))
drep_len = drep_cfg.get("min_length", drep_cfg.get("length", drep_cfg.get("l", 50000)))
batch_size = MT_C.get("batching", {}).get("batch_size", 50)
drep_block_size = int(drep_cfg.get("block_size", 10000))

# 3. MAG Discovery and Input Registration per Sample
# ------------------------------------------------------------------------------
SAMPLE_MAG_FASTAS = {}
SAMPLE_ALL_MAGS = {}
SAMPLE_GENOME_INFO = {}
SAMPLE_TAXONOMY = {}

refresh_mags = str(config.get("refresh_mags", False)).lower() in ("true", "1")

for s in ACTIVE_SAMPLES:
    s_info = get_sample_info(s)
    SAMPLE_MAG_FASTAS[s] = {}
    s_mag_dir = s_info.get("mag_dir") or s_info.get("genome_dir") or s_info.get("path") or config.get("mag_dir") or MT_C.get("mag_dir", "")
    abs_s_mag_dir = mt_full_path(s_mag_dir) if s_mag_dir else ""

    s_ginfo = s_info.get("genome_info") or config.get("genome_info") or drep_cfg.get("genome_info", "")
    SAMPLE_GENOME_INFO[s] = mt_full_path(s_ginfo) if s_ginfo else ""

    s_tax = s_info.get("taxonomy") or s_info.get("taxonomy_tsv") or config.get("taxonomy") or drep_cfg.get("taxonomy", "")
    if not s_tax:
        s_tax = s_ginfo
    SAMPLE_TAXONOMY[s] = mt_full_path(s_tax) if s_tax else ""

    manifest_cache = os.path.join(MT_BASE_RES, s, "drep", "raw_mags_manifest.tsv")

    # Fast-path: read directly from manifest_cache if already generated (<0.05s vs 60s Lustre directory scan)
    if not refresh_mags and os.path.isfile(manifest_cache) and os.path.getsize(manifest_cache) > 0:
        try:
            with open(manifest_cache, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2 and parts[0] and parts[1]:
                        SAMPLE_MAG_FASTAS[s][parts[0]] = parts[1]
        except Exception as e:
            print(f"WARNING: [mag_traits] Failed reading cached manifest {manifest_cache}: {e}", file=sys.stderr)

    # Fallback: scan mag_dir if SAMPLE_MAG_FASTAS[s] is still empty
    if not SAMPLE_MAG_FASTAS[s] and abs_s_mag_dir and os.path.exists(abs_s_mag_dir):
        valid_exts = {".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz", ".fa.gz"}
        try:
            with os.scandir(abs_s_mag_dir) as entries:
                for entry in entries:
                    name = entry.name
                    if name.startswith(".") or name.endswith((".faa", ".faa.gz", ".ffn", ".ffn.gz")):
                        continue
                    for ext in valid_exts:
                        if name.endswith(ext):
                            stem = name[:-len(ext)]
                            if stem not in SAMPLE_MAG_FASTAS[s]:
                                SAMPLE_MAG_FASTAS[s][stem] = entry.path
                            break
        except Exception as e:
            print(f"WARNING: [mag_traits] Failed scanning mag_dir {abs_s_mag_dir}: {e}", file=sys.stderr)

    # Check sample-specific external_genomes
    s_ext = s_info.get("external_genomes", {})
    if s_ext and isinstance(s_ext, dict):
        for name, info in s_ext.items():
            if isinstance(info, dict) and "path" in info:
                SAMPLE_MAG_FASTAS[s][name] = mt_full_path(info["path"])

    # Check global external_genomes declaration
    ext_genomes = config.get("external_genomes", {})
    if ext_genomes and isinstance(ext_genomes, dict):
        for name, info in ext_genomes.items():
            if isinstance(info, dict) and "path" in info and name not in SAMPLE_MAG_FASTAS[s]:
                SAMPLE_MAG_FASTAS[s][name] = mt_full_path(info["path"])

    SAMPLE_ALL_MAGS[s] = sorted(SAMPLE_MAG_FASTAS[s].keys())

# Backward compatibility accessors (first active sample)
ALL_MAGS = SAMPLE_ALL_MAGS.get(ACTIVE_SAMPLES[0], []) if ACTIVE_SAMPLES else []
MAG_FASTAS = SAMPLE_MAG_FASTAS.get(ACTIVE_SAMPLES[0], {}) if ACTIVE_SAMPLES else {}
GENOME_INFO = SAMPLE_GENOME_INFO.get(ACTIVE_SAMPLES[0], "") if ACTIVE_SAMPLES else ""

# 4. Target Registration & Step Toggles
# ------------------------------------------------------------------------------
MAG_TRAITS_TARGETS = []
steps = MT_C.get("steps", {})
do_drep = str(steps.get("drep", False)).lower() == "true"
do_prodigal = str(steps.get("prodigal", True)).lower() == "true"
do_pool_genes = str(steps.get("pool_genes", True)).lower() == "true"
do_dbcan = str(steps.get("dbcan", True)).lower() == "true"
do_camper = str(steps.get("camper", True)).lower() == "true"
do_grodon2 = str(steps.get("grodon2", steps.get("phydon", True))).lower() == "true"
do_genomespot = str(steps.get("genomespot", True)).lower() == "true"
do_summary = str(steps.get("summary", True)).lower() == "true"

for s in ACTIVE_SAMPLES:
    if do_summary:
        MAG_TRAITS_TARGETS.append(os.path.join(MT_BASE_RES, s, "summary", "all_mag_traits_matrix.tsv"))
    else:
        if do_drep:
            if two_rounds_drep:
                MAG_TRAITS_TARGETS.append(os.path.join(MT_BASE_RES, s, "drep", "round2", "data_tables", "Widb.csv"))
            else:
                MAG_TRAITS_TARGETS.append(os.path.join(MT_BASE_RES, s, "drep", "single", "data_tables", "Widb.csv"))
        else:
            MAG_TRAITS_TARGETS.append(os.path.join(MT_BASE_RES, s, "batches", "batch_list.txt"))

    if do_prodigal and do_pool_genes:
        MAG_TRAITS_TARGETS.append(os.path.join(MT_BASE_RES, s, "proteins", ".pool.done"))

TARGETS.extend(MAG_TRAITS_TARGETS)


# 5. Analysis Rules
# ------------------------------------------------------------------------------

# --- 5.0 Helper: Sample MAG Directory Lookup ---
def get_sample_mag_dir(wildcards):
    s = wildcards.sample
    s_info = get_sample_info(s)
    s_mag_dir = s_info.get("mag_dir") or s_info.get("genome_dir") or s_info.get("path") or config.get("mag_dir") or MT_C.get("mag_dir", "")
    abs_d = mt_full_path(s_mag_dir) if s_mag_dir else ""
    return [abs_d] if (abs_d and os.path.isdir(abs_d)) else []


def get_block_threads(wildcards):
    s = wildcards.sample
    b_id = wildcards.block_id
    manifest_file = os.path.join(MT_BASE_RES, s, "drep", "round1_blocks", "partitions", b_id + ".tsv")
    n_mags = 0
    if os.path.isfile(manifest_file):
        try:
            with open(manifest_file, "r", encoding="utf-8") as f:
                n_mags = sum(1 for line in f if line.strip())
        except Exception:
            pass
    max_threads = int(drep_cfg.get("threads", 64))
    if n_mags <= 1:
        return 1
    elif n_mags <= 50:
        return get_threads(min(4, max_threads))
    elif n_mags <= 500:
        return get_threads(min(8, max_threads))
    elif n_mags <= 2000:
        return get_threads(min(16, max_threads))
    elif n_mags <= 5000:
        return get_threads(min(32, max_threads))
    else:
        return get_threads(max_threads)


# --- 5.1 Optional Pre-step: Two-Round or Single-Round dRep Dereplication ---
if do_drep:
    rule mag_traits_drep_manifest:
        input:
            mag_dir = get_sample_mag_dir
        output:
            manifest = os.path.join(MT_BASE_RES, "{sample}", "drep", "raw_mags_manifest.tsv")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: get_mem_mb(threads)
        wildcard_constraints:
            sample = r"[^/]+"
        run:
            import os
            s = wildcards.sample
            os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
            mags_dict = dict(SAMPLE_MAG_FASTAS.get(s, {}))
            if not mags_dict and input.mag_dir:
                d = input.mag_dir[0] if isinstance(input.mag_dir, list) else str(input.mag_dir)
                if os.path.isdir(d):
                    valid_exts = {".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz", ".fa.gz"}
                    for entry in os.scandir(d):
                        name = entry.name
                        if name.startswith(".") or name.endswith((".faa", ".faa.gz", ".ffn", ".ffn.gz")):
                            continue
                        for ext in valid_exts:
                            if name.endswith(ext):
                                stem = name[:-len(ext)]
                                if stem not in mags_dict:
                                    mags_dict[stem] = entry.path
                                break
            with open(output.manifest, "w", encoding="utf-8") as f:
                for m in sorted(mags_dict.keys()):
                    p = mags_dict[m]
                    if p:
                        f.write(f"{m}\t{os.path.abspath(p)}\n")

    if two_rounds_drep or drep_block_size > 0:
        if drep_block_size > 0:
            checkpoint mag_traits_drep_round1_partition:
                input:
                    manifest = os.path.join(MT_BASE_RES, "{sample}", "drep", "raw_mags_manifest.tsv")
                output:
                    blocks_dir = directory(os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "partitions")),
                    block_list = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "block_list.txt")
                log:
                    os.path.join(MT_BASE_LOGS, "{sample}", "drep", "round1_partition.log")
                threads: 1
                resources:
                    runtime = "7d",
                    mem_mb = lambda wildcards, threads: get_mem_mb(threads)
                params:
                    block_size = drep_block_size,
                    min_phylum_size = drep_cfg.get("min_phylum_size", 500),
                    taxonomy = lambda wildcards: SAMPLE_TAXONOMY.get(wildcards.sample, ""),
                    partition_script = os.path.join(WORKFLOW_ROOT, "scripts", "partition_manifest.py")
                wildcard_constraints:
                    sample = r"[^/]+"
                shell:
                    """
                    set -euo pipefail
                    mkdir -p "{output.blocks_dir}" "$(dirname "{output.block_list}")" "$(dirname "{log}")"
                    TAX_ARG=""
                    if [ -n "{params.taxonomy}" ] && [ -f "{params.taxonomy}" ]; then
                        TAX_ARG="--taxonomy {params.taxonomy}"
                    fi
                    python3 "{params.partition_script}" \
                        --manifest "{input.manifest}" \
                        $TAX_ARG \
                        --block-size {params.block_size} \
                        --min-phylum-size {params.min_phylum_size} \
                        --out-dir "{output.blocks_dir}" \
                        --out-list "{output.block_list}" 2>&1 | tee "{log}"
                    """

            rule mag_traits_drep_round1_block:
                input:
                    manifest = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "partitions", "{block_id}.tsv")
                output:
                    block_derep = directory(os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "{block_id}", "dereplicated_genomes")),
                    block_cdb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "{block_id}", "data_tables", "Cdb.csv"),
                    block_widb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "{block_id}", "data_tables", "Widb.csv"),
                    block_ginfo = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "{block_id}", "genomeInfo.csv")
                log:
                    os.path.join(MT_BASE_LOGS, "{sample}", "drep", "round1_blocks", "{block_id}.log")
                threads: get_block_threads
                resources:
                    runtime = "7d",
                    mem_mb = lambda wildcards, threads: get_mem_mb(threads)
                params:
                    env = mt_env("drep"),
                    ani = r1_ani if two_rounds_drep else single_ani,
                    comp = drep_comp,
                    cont = drep_cont,
                    length = drep_len,
                    genome_info = lambda wildcards: SAMPLE_GENOME_INFO.get(wildcards.sample, ""),
                    extra_flags = drep_cfg.get("extra_flags", ""),
                    prepare_script = os.path.join(WORKFLOW_ROOT, "scripts", "prepare_drep_inputs.py"),
                    out_dir = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "{block_id}")
                wildcard_constraints:
                    sample = r"[^/]+",
                    block_id = r"[^/]+"
                shell:
                    """
                    set -euo pipefail
                    mkdir -p "{params.out_dir}" "$(dirname "{log}")"
                    TMP_DIR="{params.out_dir}/tmp_input_genomes"
                    rm -rf "$TMP_DIR"
                    mkdir -p "$TMP_DIR"
                    trap 'rm -rf "$TMP_DIR"' EXIT

                    python3 "{params.prepare_script}" \
                        --manifest "{input.manifest}" \
                        --stage-dir "$TMP_DIR" \
                        --out-list "$TMP_DIR/genomes_list.txt" \
                        --genome-info-src "{params.genome_info}" \
                        --genome-info-dst "{output.block_ginfo}"

                    NUM_GENOMES=$(grep -cve '^\\s*$' "$TMP_DIR/genomes_list.txt" || true)

                    if [ "$NUM_GENOMES" -le 1 ]; then
                        echo "[round1_block] Block '{wildcards.block_id}' has $NUM_GENOMES input genome(s). Direct bypass without dRep..." | tee -a "{log}"
                        rm -rf "{output.block_derep}"
                        mkdir -p "{output.block_derep}" "$(dirname "{output.block_cdb}")"

                        if [ "$NUM_GENOMES" -eq 1 ]; then
                            SINGLE_PATH=$(head -n 1 "$TMP_DIR/genomes_list.txt" | tr -d '[:space:]')
                            SINGLE_NAME=$(basename "$SINGLE_PATH")
                            cp -f "$SINGLE_PATH" "{output.block_derep}/$SINGLE_NAME"
                            echo "genome,secondary_cluster,primary_cluster" > "{output.block_cdb}"
                            echo "$SINGLE_NAME,1_0,1" >> "{output.block_cdb}"
                            echo "genome,score,cluster" > "{output.block_widb}"
                            echo "$SINGLE_NAME,100.0,1_0" >> "{output.block_widb}"
                        else
                            echo "genome,secondary_cluster,primary_cluster" > "{output.block_cdb}"
                            echo "genome,score,cluster" > "{output.block_widb}"
                        fi

                        if [ ! -f "{output.block_ginfo}" ] || [ ! -s "{output.block_ginfo}" ]; then
                            touch "{output.block_ginfo}"
                        fi
                        echo "[round1_block] Bypass completed successfully." | tee -a "{log}"
                    else
                        CMD_PREFIX=""
                        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                            CMD_PREFIX="conda run -p {params.env}"
                        fi

                        DREP_BIN="dRep"
                        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                            if [ -x "{params.env}/bin/drep" ] && [ ! -x "{params.env}/bin/dRep" ]; then
                                DREP_BIN="drep"
                            fi
                        else
                            if command -v drep >/dev/null 2>&1 && ! command -v dRep >/dev/null 2>&1; then
                                DREP_BIN="drep"
                            fi
                        fi

                        GINFO_ARG=""
                        if [ -s "{output.block_ginfo}" ]; then
                            GINFO_ARG="--genomeInfo {output.block_ginfo}"
                        else
                            if [[ "{params.extra_flags}" != *"--ignoreGenomeQuality"* ]]; then
                                GINFO_ARG="--ignoreGenomeQuality"
                            fi
                        fi

                        set +e
                        $CMD_PREFIX $DREP_BIN dereplicate "{params.out_dir}" \
                            -g "$TMP_DIR/genomes_list.txt" \
                            -p {threads} \
                            -sa {params.ani} \
                            -comp {params.comp} \
                            -con {params.cont} \
                            -l {params.length} \
                            $GINFO_ARG \
                            {params.extra_flags} > "{log}" 2>&1
                        DREP_EXIT=$?
                        set -e
                        cat "{log}"

                        if [ "$DREP_EXIT" -ne 0 ]; then
                            PASS_COUNT=0
                            if [ -f "{params.out_dir}/data_tables/Chdb.csv" ]; then
                                PASS_COUNT=$(tail -n +2 "{params.out_dir}/data_tables/Chdb.csv" | grep -cve '^\\s*$' || true)
                            elif [ -f "{params.out_dir}/data_tables/genomeInfo.csv" ]; then
                                PASS_COUNT=$(tail -n +2 "{params.out_dir}/data_tables/genomeInfo.csv" | grep -cve '^\\s*$' || true)
                            fi

                            if [ "$PASS_COUNT" -eq 1 ]; then
                                echo "[round1_block] Caught dRep <=1 genome filtering edge case in block '{wildcards.block_id}'. Recovering single survivor..." | tee -a "{log}"
                                SURVIVOR=""
                                if [ -f "{params.out_dir}/data_tables/Chdb.csv" ]; then
                                    SURVIVOR=$(tail -n +2 "{params.out_dir}/data_tables/Chdb.csv" | head -n 1 | cut -d',' -f1 | tr -d '[:space:]')
                                fi
                                if [ -z "$SURVIVOR" ] && [ -f "{params.out_dir}/data_tables/genomeInfo.csv" ]; then
                                    SURVIVOR=$(tail -n +2 "{params.out_dir}/data_tables/genomeInfo.csv" | head -n 1 | cut -d',' -f1 | tr -d '[:space:]')
                                fi
                                SURV_NAME=$(basename "$SURVIVOR")
                                rm -rf "{output.block_derep}"
                                mkdir -p "{output.block_derep}" "$(dirname "{output.block_cdb}")"
                                if [ -f "$TMP_DIR/$SURV_NAME" ]; then
                                    cp -f "$TMP_DIR/$SURV_NAME" "{output.block_derep}/$SURV_NAME"
                                elif [ -f "$SURVIVOR" ]; then
                                    cp -f "$SURVIVOR" "{output.block_derep}/$SURV_NAME"
                                fi
                                echo "genome,secondary_cluster,primary_cluster" > "{output.block_cdb}"
                                echo "$SURV_NAME,1_0,1" >> "{output.block_cdb}"
                                echo "genome,score,cluster" > "{output.block_widb}"
                                echo "$SURV_NAME,100.0,1_0" >> "{output.block_widb}"
                                echo "[round1_block] Recovery completed successfully." | tee -a "{log}"
                            else
                                echo "ERROR: dRep dereplicate failed with exit code $DREP_EXIT" >&2
                                exit "$DREP_EXIT"
                            fi
                        fi

                        if [ ! -f "{output.block_widb}" ]; then
                            if [ -f "{params.out_dir}/data_tables/Wdb.csv" ]; then
                                cp -f "{params.out_dir}/data_tables/Wdb.csv" "{output.block_widb}"
                            else
                                touch "{output.block_widb}"
                            fi
                        fi
                    fi
                    """

            def get_round1_blocks(wildcards):
                s = wildcards.sample
                block_list_file = checkpoints.mag_traits_drep_round1_partition.get(sample=s).output.block_list
                blocks = []
                if os.path.exists(block_list_file):
                    with open(block_list_file, "r", encoding="utf-8") as f:
                        for line in f:
                            b = line.strip()
                            if b:
                                blocks.append(b)
                return blocks

            def get_round1_block_dereps(wildcards):
                s = wildcards.sample
                blocks = get_round1_blocks(wildcards)
                return [os.path.join(MT_BASE_RES, s, "drep", "round1_blocks", b, "dereplicated_genomes") for b in blocks]

            def get_round1_block_cdbs(wildcards):
                s = wildcards.sample
                blocks = get_round1_blocks(wildcards)
                return [os.path.join(MT_BASE_RES, s, "drep", "round1_blocks", b, "data_tables", "Cdb.csv") for b in blocks]

            def get_round1_block_widbs(wildcards):
                s = wildcards.sample
                blocks = get_round1_blocks(wildcards)
                return [os.path.join(MT_BASE_RES, s, "drep", "round1_blocks", b, "data_tables", "Widb.csv") for b in blocks]

            def get_round1_block_ginfos(wildcards):
                s = wildcards.sample
                blocks = get_round1_blocks(wildcards)
                return [os.path.join(MT_BASE_RES, s, "drep", "round1_blocks", b, "genomeInfo.csv") for b in blocks]

            checkpoint mag_traits_drep_round1:
                input:
                    manifest = os.path.join(MT_BASE_RES, "{sample}", "drep", "raw_mags_manifest.tsv"),
                    block_list = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1_blocks", "block_list.txt"),
                    block_dereps = get_round1_block_dereps,
                    block_cdbs = get_round1_block_cdbs,
                    block_widbs = get_round1_block_widbs,
                    block_ginfos = get_round1_block_ginfos
                output:
                    r1_derep = directory(os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "dereplicated_genomes")),
                    r1_cdb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "data_tables", "Cdb.csv"),
                    r1_widb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "data_tables", "Widb.csv"),
                    r1_ginfo = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "genomeInfo.csv")
                log:
                    os.path.join(MT_BASE_LOGS, "{sample}", "drep", "round1_consolidation.log")
                threads: get_threads(drep_cfg.get("threads", 64))
                resources:
                    runtime = "7d",
                    mem_mb = lambda wildcards, threads: get_mem_mb(threads)
                params:
                    env = mt_env("drep"),
                    ani = r1_ani if two_rounds_drep else single_ani,
                    comp = drep_comp,
                    cont = drep_cont,
                    length = drep_len,
                    genome_info = lambda wildcards: SAMPLE_GENOME_INFO.get(wildcards.sample, ""),
                    extra_flags = drep_cfg.get("extra_flags", ""),
                    prepare_script = os.path.join(WORKFLOW_ROOT, "scripts", "prepare_drep_inputs.py"),
                    chain_script = os.path.join(WORKFLOW_ROOT, "scripts", "chain_drep_clusters.py"),
                    drep_utils = os.path.join(WORKFLOW_ROOT, "scripts", "mag_traits_drep_utils.py"),
                    out_dir = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1")
                wildcard_constraints:
                    sample = r"[^/]+"
                shell:
                    """
                    set -euo pipefail
                    mkdir -p "{params.out_dir}" "$(dirname "{log}")"
                    NUM_BLOCKS=$(grep -cve '^\\s*$' "{input.block_list}" || true)

                    if [ "$NUM_BLOCKS" -le 1 ]; then
                        echo "[round1_consolidation] Single block detected ($NUM_BLOCKS block). Direct pass-through..."
                        FIRST_BLOCK=$(head -n 1 "{input.block_list}" | tr -d '[:space:]')
                        BLOCK_DIR="$(dirname "{input.block_list}")/$FIRST_BLOCK"

                        rm -rf "{output.r1_derep}"
                        mkdir -p "{output.r1_derep}"
                        if [ -d "$BLOCK_DIR/dereplicated_genomes" ]; then
                            cp -a "$BLOCK_DIR/dereplicated_genomes"/* "{output.r1_derep}/" 2>/dev/null || true
                        fi

                        mkdir -p "$(dirname "{output.r1_cdb}")"
                        cp -f "$BLOCK_DIR/data_tables/Cdb.csv" "{output.r1_cdb}"
                        if [ -f "$BLOCK_DIR/data_tables/Widb.csv" ]; then
                            cp -f "$BLOCK_DIR/data_tables/Widb.csv" "{output.r1_widb}"
                        elif [ -f "$BLOCK_DIR/data_tables/Wdb.csv" ]; then
                            cp -f "$BLOCK_DIR/data_tables/Wdb.csv" "{output.r1_widb}"
                        else
                            touch "{output.r1_widb}"
                        fi
                        cp -f "$BLOCK_DIR/genomeInfo.csv" "{output.r1_ginfo}"
                        echo "[round1_consolidation] Pass-through completed successfully."
                    elif grep -q '^tax_' "{input.block_list}"; then
                        echo "[round1_consolidation] Taxonomy blocks detected ($NUM_BLOCKS blocks). Performing direct union..."
                        python3 "{params.drep_utils}" union-blocks \
                            --block-list "{input.block_list}" \
                            --blocks-dir "$(dirname "{input.block_list}")" \
                            --out-derep "{output.r1_derep}" \
                            --out-cdb "{output.r1_cdb}" \
                            --out-widb "{output.r1_widb}" \
                            --out-ginfo "{output.r1_ginfo}"
                        echo "[round1_consolidation] Taxonomy direct union completed successfully."
                    else
                        echo "[round1_consolidation] Multi-block detected ($NUM_BLOCKS blocks). Staging pooled winners..."
                        TMP_DIR="{params.out_dir}/tmp_pool"
                        rm -rf "$TMP_DIR"
                        mkdir -p "$TMP_DIR/staged" "$TMP_DIR/pool_drep"
                        trap 'rm -rf "$TMP_DIR"' EXIT

                        POOLED_MANIFEST="$TMP_DIR/pooled_manifest.tsv"
                        python3 "{params.drep_utils}" pool-manifest \
                            --out-manifest "$POOLED_MANIFEST" \
                            --block-dirs {input.block_dereps}

                        N_POOLED=$(wc -l < "$POOLED_MANIFEST" | tr -d '[:space:]')
                        echo "[round1_consolidation] Total pooled candidates from all blocks: $N_POOLED"

                        if [ "$N_POOLED" -eq 0 ]; then
                            echo "ERROR: No dereplicated genomes found in block outputs!" >&2
                            exit 1
                        fi

                        python3 "{params.prepare_script}" \
                            --manifest "$POOLED_MANIFEST" \
                            --stage-dir "$TMP_DIR/staged" \
                            --out-list "$TMP_DIR/genomes_list.txt" \
                            --genome-info-src "{params.genome_info}" \
                            --genome-info-dst "$TMP_DIR/pool_genomeInfo.csv"

                        CMD_PREFIX=""
                        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                            CMD_PREFIX="conda run -p {params.env}"
                        fi

                        DREP_BIN="dRep"
                        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                            if [ -x "{params.env}/bin/drep" ] && [ ! -x "{params.env}/bin/dRep" ]; then
                                DREP_BIN="drep"
                            fi
                        else
                            if command -v drep >/dev/null 2>&1 && ! command -v dRep >/dev/null 2>&1; then
                                DREP_BIN="drep"
                            fi
                        fi

                        GINFO_ARG=""
                        if [ -s "$TMP_DIR/pool_genomeInfo.csv" ]; then
                            GINFO_ARG="--genomeInfo $TMP_DIR/pool_genomeInfo.csv"
                        else
                            if [[ "{params.extra_flags}" != *"--ignoreGenomeQuality"* ]]; then
                                GINFO_ARG="--ignoreGenomeQuality"
                            fi
                        fi

                        $CMD_PREFIX $DREP_BIN dereplicate "$TMP_DIR/pool_drep" \
                            -g "$TMP_DIR/genomes_list.txt" \
                            -p {threads} \
                            -sa {params.ani} \
                            -comp {params.comp} \
                            -con {params.cont} \
                            -l {params.length} \
                            $GINFO_ARG \
                            {params.extra_flags} 2>&1 | tee "{log}"

                        POOL_WIDB="$TMP_DIR/pool_drep/data_tables/Widb.csv"
                        if [ ! -f "$POOL_WIDB" ]; then
                            if [ -f "$TMP_DIR/pool_drep/data_tables/Wdb.csv" ]; then
                                cp -f "$TMP_DIR/pool_drep/data_tables/Wdb.csv" "$POOL_WIDB"
                            else
                                touch "$POOL_WIDB"
                            fi
                        fi

                        PROV_TSV="$(dirname "{output.r1_cdb}")/provenance.tsv"
                        python3 "{params.chain_script}" \
                            --block-cdbs {input.block_cdbs} \
                            --block-widbs {input.block_widbs} \
                            --pool-cdb "$TMP_DIR/pool_drep/data_tables/Cdb.csv" \
                            --pool-widb "$POOL_WIDB" \
                            --out-cdb "{output.r1_cdb}" \
                            --out-widb "{output.r1_widb}" \
                            --out-provenance "$PROV_TSV"

                        rm -rf "{output.r1_derep}"
                        mkdir -p "{output.r1_derep}"
                        if [ -d "$TMP_DIR/pool_drep/dereplicated_genomes" ]; then
                            cp -a "$TMP_DIR/pool_drep/dereplicated_genomes"/* "{output.r1_derep}/" 2>/dev/null || true
                        fi

                        python3 "{params.drep_utils}" merge-genome-info \
                            --out "{output.r1_ginfo}" \
                            --inputs {input.block_ginfos}
                        echo "[round1_consolidation] Multi-block consolidation completed successfully."
                    fi
                    """
        else:
            checkpoint mag_traits_drep_round1:
                input:
                    manifest = os.path.join(MT_BASE_RES, "{sample}", "drep", "raw_mags_manifest.tsv")
                output:
                    r1_derep = directory(os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "dereplicated_genomes")),
                    r1_cdb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "data_tables", "Cdb.csv"),
                    r1_widb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "data_tables", "Widb.csv"),
                    r1_ginfo = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "genomeInfo.csv")
                log:
                    os.path.join(MT_BASE_LOGS, "{sample}", "drep", "round1.log")
                threads: get_threads(drep_cfg.get("threads", 64))
                resources:
                    runtime = "7d",
                    mem_mb = lambda wildcards, threads: get_mem_mb(threads)
                params:
                    env = mt_env("drep"),
                    ani = r1_ani,
                    comp = drep_comp,
                    cont = drep_cont,
                    length = drep_len,
                    genome_info = lambda wildcards: SAMPLE_GENOME_INFO.get(wildcards.sample, ""),
                    extra_flags = drep_cfg.get("extra_flags", ""),
                    prepare_script = os.path.join(WORKFLOW_ROOT, "scripts", "prepare_drep_inputs.py"),
                    out_dir = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1")
                wildcard_constraints:
                    sample = r"[^/]+"
                shell:
                    """
                    set -euo pipefail
                    mkdir -p "{params.out_dir}" "$(dirname "{log}")"
                    TMP_DIR="{params.out_dir}/tmp_input_genomes"
                    rm -rf "$TMP_DIR"
                    mkdir -p "$TMP_DIR"
                    trap 'rm -rf "$TMP_DIR"' EXIT

                    python3 "{params.prepare_script}" \
                        --manifest "{input.manifest}" \
                        --stage-dir "$TMP_DIR" \
                        --out-list "$TMP_DIR/genomes_list.txt" \
                        --genome-info-src "{params.genome_info}" \
                        --genome-info-dst "{output.r1_ginfo}"

                    CMD_PREFIX=""
                    if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                        CMD_PREFIX="conda run -p {params.env}"
                    fi

                    DREP_BIN="dRep"
                    if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                        if [ -x "{params.env}/bin/drep" ] && [ ! -x "{params.env}/bin/dRep" ]; then
                            DREP_BIN="drep"
                        fi
                    else
                        if command -v drep >/dev/null 2>&1 && ! command -v dRep >/dev/null 2>&1; then
                            DREP_BIN="drep"
                        fi
                    fi

                    GINFO_ARG=""
                    if [ -s "{output.r1_ginfo}" ]; then
                        GINFO_ARG="--genomeInfo {output.r1_ginfo}"
                    else
                        if [[ "{params.extra_flags}" != *"--ignoreGenomeQuality"* ]]; then
                            GINFO_ARG="--ignoreGenomeQuality"
                        fi
                    fi

                    $CMD_PREFIX $DREP_BIN dereplicate "{params.out_dir}" \
                        -g "$TMP_DIR/genomes_list.txt" \
                        -p {threads} \
                        -sa {params.ani} \
                        -comp {params.comp} \
                        -con {params.cont} \
                        -l {params.length} \
                        $GINFO_ARG \
                        {params.extra_flags} 2>&1 | tee "{log}"

                    if [ ! -f "{output.r1_widb}" ]; then
                        if [ -f "{params.out_dir}/data_tables/Wdb.csv" ]; then
                            cp -f "{params.out_dir}/data_tables/Wdb.csv" "{output.r1_widb}"
                        else
                            touch "{output.r1_widb}"
                        fi
                    fi
                    """

        if two_rounds_drep:
            rule mag_traits_prep_round2_genome_info:
                input:
                    r1_derep = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "dereplicated_genomes"),
                    r1_ginfo = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "genomeInfo.csv")
                output:
                    r2_ginfo = os.path.join(MT_BASE_RES, "{sample}", "drep", "round2", "genomeInfo.csv")
                log:
                    os.path.join(MT_BASE_LOGS, "{sample}", "drep", "round2_prep_ginfo.log")
                threads: 1
                resources:
                    runtime = "7d",
                    mem_mb = lambda wildcards, threads: get_mem_mb(threads)
                params:
                    drep_utils = os.path.join(WORKFLOW_ROOT, "scripts", "mag_traits_drep_utils.py")
                wildcard_constraints:
                    sample = r"[^/]+"
                shell:
                    """
                    (
                    set -euo pipefail
                    mkdir -p "$(dirname "{output.r2_ginfo}")" "$(dirname "{log}")"
                    python3 "{params.drep_utils}" filter-genome-info \
                        --derep-dir "{input.r1_derep}" \
                        --src "{input.r1_ginfo}" \
                        --dst "{output.r2_ginfo}"
                    ) > "{log}" 2>&1
                    """

            checkpoint mag_traits_drep_round2:
                input:
                    r1_derep = os.path.join(MT_BASE_RES, "{sample}", "drep", "round1", "dereplicated_genomes"),
                    r2_ginfo = os.path.join(MT_BASE_RES, "{sample}", "drep", "round2", "genomeInfo.csv")
                output:
                    r2_derep = directory(os.path.join(MT_BASE_RES, "{sample}", "drep", "round2", "dereplicated_genomes")),
                    r2_cdb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round2", "data_tables", "Cdb.csv"),
                    r2_widb = os.path.join(MT_BASE_RES, "{sample}", "drep", "round2", "data_tables", "Widb.csv")
                log:
                    os.path.join(MT_BASE_LOGS, "{sample}", "drep", "round2.log")
                threads: get_threads(drep_cfg.get("threads", 64))
                resources:
                    runtime = "7d",
                    mem_mb = lambda wildcards, threads: get_mem_mb(threads)
                params:
                    env = mt_env("drep"),
                    ani = r2_ani,
                    comp = drep_comp,
                    cont = drep_cont,
                    length = drep_len,
                    extra_flags = drep_cfg.get("extra_flags", ""),
                    out_dir = os.path.join(MT_BASE_RES, "{sample}", "drep", "round2")
                wildcard_constraints:
                    sample = r"[^/]+"
                shell:
                    """
                    set -euo pipefail
                    mkdir -p "{params.out_dir}" "$(dirname "{log}")"

                    CMD_PREFIX=""
                    if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                        CMD_PREFIX="conda run -p {params.env}"
                    fi

                    DREP_BIN="dRep"
                    if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                        if [ -x "{params.env}/bin/drep" ] && [ ! -x "{params.env}/bin/dRep" ]; then
                            DREP_BIN="drep"
                        fi
                    else
                        if command -v drep >/dev/null 2>&1 && ! command -v dRep >/dev/null 2>&1; then
                            DREP_BIN="drep"
                        fi
                    fi

                    R1_LIST="{params.out_dir}/r1_derep_list.txt"
                    find "$(cd "{input.r1_derep}" && pwd)" -maxdepth 1 -name "*.fa" > "$R1_LIST"
                    if [ ! -s "$R1_LIST" ]; then
                        find "$(cd "{input.r1_derep}" && pwd)" -maxdepth 1 -type f > "$R1_LIST"
                    fi
                    if [ ! -s "$R1_LIST" ]; then
                        echo "ERROR: No dereplicated genomes found in {input.r1_derep}" >&2
                        exit 1
                    fi
                    trap 'rm -f "$R1_LIST"' EXIT

                    GINFO_ARG=""
                    if [ -s "{input.r2_ginfo}" ]; then
                        GINFO_ARG="--genomeInfo {input.r2_ginfo}"
                    else
                        if [[ "{params.extra_flags}" != *"--ignoreGenomeQuality"* ]]; then
                            GINFO_ARG="--ignoreGenomeQuality"
                        fi
                    fi

                    $CMD_PREFIX $DREP_BIN dereplicate "{params.out_dir}" \
                        -g "$R1_LIST" \
                        -p {threads} \
                        -sa {params.ani} \
                        -comp {params.comp} \
                        -con {params.cont} \
                        -l {params.length} \
                        $GINFO_ARG \
                        {params.extra_flags} 2>&1 | tee "{log}"

                    if [ ! -f "{output.r2_widb}" ]; then
                        if [ -f "{params.out_dir}/data_tables/Wdb.csv" ]; then
                            cp -f "{params.out_dir}/data_tables/Wdb.csv" "{output.r2_widb}"
                        else
                            touch "{output.r2_widb}"
                        fi
                    fi
                    """
    else:
        checkpoint mag_traits_drep_single:
            input:
                manifest = os.path.join(MT_BASE_RES, "{sample}", "drep", "raw_mags_manifest.tsv")
            output:
                derep_dir = directory(os.path.join(MT_BASE_RES, "{sample}", "drep", "single", "dereplicated_genomes")),
                cdb = os.path.join(MT_BASE_RES, "{sample}", "drep", "single", "data_tables", "Cdb.csv"),
                widb = os.path.join(MT_BASE_RES, "{sample}", "drep", "single", "data_tables", "Widb.csv")
            log:
                os.path.join(MT_BASE_LOGS, "{sample}", "drep", "single.log")
            threads: get_threads(drep_cfg.get("threads", 16))
            resources:
                runtime = "7d",
                mem_mb = lambda wildcards, threads: get_mem_mb(threads)
            params:
                env = mt_env("drep"),
                ani = single_ani,
                comp = drep_comp,
                cont = drep_cont,
                length = drep_len,
                genome_info = lambda wildcards: SAMPLE_GENOME_INFO.get(wildcards.sample, ""),
                extra_flags = drep_cfg.get("extra_flags", "--ignoreGenomeQuality"),
                prepare_script = os.path.join(WORKFLOW_ROOT, "scripts", "prepare_drep_inputs.py"),
                out_dir = os.path.join(MT_BASE_RES, "{sample}", "drep", "single")
            wildcard_constraints:
                sample = r"[^/]+"
            shell:
                """
                set -euo pipefail
                mkdir -p "{params.out_dir}" "$(dirname "{log}")"
                TMP_DIR="{params.out_dir}/tmp_input_genomes"
                rm -rf "$TMP_DIR"
                mkdir -p "$TMP_DIR"
                trap 'rm -rf "$TMP_DIR"' EXIT

                python3 "{params.prepare_script}" \
                    --manifest "{input.manifest}" \
                    --stage-dir "$TMP_DIR" \
                    --out-list "$TMP_DIR/genomes_list.txt" \
                    --genome-info-src "{params.genome_info}" \
                    --genome-info-dst "$TMP_DIR/genomeInfo.csv"

                CMD_PREFIX=""
                if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                    CMD_PREFIX="conda run -p {params.env}"
                fi

                DREP_BIN="dRep"
                if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
                    if [ -x "{params.env}/bin/drep" ] && [ ! -x "{params.env}/bin/dRep" ]; then
                        DREP_BIN="drep"
                    fi
                else
                    if command -v drep >/dev/null 2>&1 && ! command -v dRep >/dev/null 2>&1; then
                        DREP_BIN="drep"
                    fi
                fi

                GINFO_ARG=""
                if [ -s "$TMP_DIR/genomeInfo.csv" ]; then
                    GINFO_ARG="--genomeInfo $TMP_DIR/genomeInfo.csv"
                else
                    if [[ "{params.extra_flags}" != *"--ignoreGenomeQuality"* ]]; then
                        GINFO_ARG="--ignoreGenomeQuality"
                    fi
                fi

                $CMD_PREFIX $DREP_BIN dereplicate "{params.out_dir}" \
                    -g "$TMP_DIR/genomes_list.txt" \
                    -p {threads} \
                    -sa {params.ani} \
                    -comp {params.comp} \
                    -con {params.cont} \
                    -l {params.length} \
                    $GINFO_ARG \
                    {params.extra_flags} 2>&1 | tee "{log}"

                if [ ! -f "{output.widb}" ]; then
                    if [ -f "{params.out_dir}/data_tables/Wdb.csv" ]; then
                        cp -f "{params.out_dir}/data_tables/Wdb.csv" "{output.widb}"
                    else
                        touch "{output.widb}"
                    fi
                fi
                """


# --- 5.2 Batch Partitioning (Groups MAGs into Batches per Sample) ---
def get_partition_input(wildcards):
    s = wildcards.sample
    if do_drep:
        if two_rounds_drep:
            return checkpoints.mag_traits_drep_round2.get(sample=s).output.r2_derep
        elif drep_block_size > 0:
            return checkpoints.mag_traits_drep_round1.get(sample=s).output.r1_derep
        else:
            return checkpoints.mag_traits_drep_single.get(sample=s).output.derep_dir
    else:
        return os.path.join(MT_BASE_RES, s, "drep", "raw_mags_manifest.tsv")


checkpoint mag_traits_partition_batches:
    input:
        mags_source = get_partition_input
    output:
        batch_list = os.path.join(MT_BASE_RES, "{sample}", "batches", "batch_list.txt")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "batches", "partition.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        batch_size = batch_size,
        batches_dir = os.path.join(MT_BASE_RES, "{sample}", "batches")
    wildcard_constraints:
        sample = r"[^/]+"
    run:
        import os, sys
        s = wildcards.sample
        os.makedirs(params.batches_dir, exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        mags = []
        src = input.mags_source
        if isinstance(src, str) and os.path.isdir(src):
            for f in sorted(os.listdir(src)):
                for ext in ['.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz']:
                    if f.endswith(ext):
                        stem = f[:-len(ext)]
                        mags.append((stem, os.path.abspath(os.path.join(src, f))))
                        break
        elif isinstance(src, str) and os.path.isfile(src):
            with open(src, 'r', encoding='utf-8') as mf:
                for line in mf:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2 and parts[0] and parts[1]:
                        mags.append((parts[0], os.path.abspath(parts[1])))
        else:
            for stem in SAMPLE_ALL_MAGS.get(s, []):
                fpath = SAMPLE_MAG_FASTAS.get(s, {}).get(stem, "")
                if fpath:
                    mags.append((stem, os.path.abspath(fpath)))

        seen = set()
        dedup_mags = [m for m in mags if not (m[0] in seen or seen.add(m[0]))]

        total = len(dedup_mags)
        batch_ids = []
        bsize = int(params.batch_size)
        if total > 0:
            num_batches = (total + bsize - 1) // bsize
            for i in range(num_batches):
                b_id = 'batch_%04d' % (i + 1)
                batch_ids.append(b_id)
                chunk = dedup_mags[i * bsize : (i + 1) * bsize]
                b_dir = os.path.join(params.batches_dir, b_id)
                os.makedirs(b_dir, exist_ok=True)
                with open(os.path.join(b_dir, 'mag_list.tsv'), 'w', encoding='utf-8') as out:
                    for stem, fpath in chunk:
                        out.write(stem + '\t' + fpath + '\n')

        with open(output.batch_list, 'w', encoding='utf-8') as out:
            for b_id in batch_ids:
                out.write(b_id + '\n')

        msg = 'INFO: Successfully partitioned %d MAGs into %d batches (size=%d) at %s\n' % (total, len(batch_ids), bsize, params.batches_dir)
        print(msg)
        with open(log[0], 'w', encoding='utf-8') as lf:
            lf.write(msg)


# --- 5.3 Batch Prodigal Gene Prediction ---
rule mag_traits_batch_prodigal:
    input:
        mag_list = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "mag_list.tsv")
    output:
        done = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "prodigal.done")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "batches", "{batch_id}", "prodigal.log")
    threads: get_threads(MT_C.get("prodigal", {}).get("threads", 2))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = mt_env("prodigal"),
        mode = MT_C.get("prodigal", {}).get("mode", "meta"),
        batch_dir = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}")
    wildcard_constraints:
        sample = r"[^/]+",
        batch_id = r"batch_\d+"
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.batch_dir}/faa" "{params.batch_dir}/ffn" "{params.batch_dir}/gff" "$(dirname "{log}")"
        PROD_BIN="prodigal"
        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
            if [ -x "{params.env}/bin/pyrodigal" ] && [ ! -x "{params.env}/bin/prodigal" ]; then
                PROD_BIN="pyrodigal"
            fi
            EXEC_CMD="conda run -p {params.env} $PROD_BIN"
        else
            if command -v pyrodigal >/dev/null 2>&1 && ! command -v prodigal >/dev/null 2>&1; then
                PROD_BIN="pyrodigal"
            fi
            EXEC_CMD="$PROD_BIN"
        fi

        while IFS=$'\t' read -r mag fna_path; do
            [ -z "$mag" ] && continue
            FAA_OUT="{params.batch_dir}/faa/$mag.faa"
            FFN_OUT="{params.batch_dir}/ffn/$mag.ffn"
            GFF_OUT="{params.batch_dir}/gff/$mag.gff"

            INPUT_FNA="$fna_path"
            TMP_FNA=""
            if [[ "$fna_path" == *.gz ]]; then
                TMP_FNA="$(mktemp --suffix=.fna)"
                gzip -dc "$fna_path" > "$TMP_FNA"
                INPUT_FNA="$TMP_FNA"
            fi

            $EXEC_CMD -i "$INPUT_FNA" -a "$FAA_OUT" -d "$FFN_OUT" -o "$GFF_OUT" -f gff -p "{params.mode}" >> "{log}" 2>&1 || true

            if [ -n "$TMP_FNA" ] && [ -f "$TMP_FNA" ]; then
                rm -f "$TMP_FNA"
            fi
            touch "$FAA_OUT" "$FFN_OUT" "$GFF_OUT"
        done < "{input.mag_list}"

        touch "{output.done}"
        """


# --- 5.4 Pool Predicted Genes and Proteins Across Batches ---
def get_all_prodigal_done(wildcards):
    s = wildcards.sample
    batch_list_file = checkpoints.mag_traits_partition_batches.get(sample=s).output.batch_list
    batches = []
    if os.path.exists(batch_list_file):
        with open(batch_list_file, "r", encoding="utf-8") as f:
            for line in f:
                b = line.strip()
                if b:
                    batches.append(b)
    return [os.path.join(MT_BASE_RES, s, "batches", b, "prodigal.done") for b in batches]


rule mag_traits_pool_genes_and_proteins:
    input:
        prodigal_dones = get_all_prodigal_done,
        batch_list = os.path.join(MT_BASE_RES, "{sample}", "batches", "batch_list.txt")
    output:
        proteins_dir = directory(os.path.join(MT_BASE_RES, "{sample}", "proteins", "all")),
        genes_dir = directory(os.path.join(MT_BASE_RES, "{sample}", "genes")),
        done = os.path.join(MT_BASE_RES, "{sample}", "proteins", ".pool.done")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "pool_genes_and_proteins.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        batches_dir = os.path.join(MT_BASE_RES, "{sample}", "batches"),
        script = os.path.join(WORKFLOW_ROOT, "scripts", "pool_mag_genes_and_proteins.py")
    wildcard_constraints:
        sample = r"[^/]+"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        python3 "{params.script}" \
            --batches_dir "{params.batches_dir}" \
            --batch_list "{input.batch_list}" \
            --proteins_dir "{output.proteins_dir}" \
            --genes_dir "{output.genes_dir}" \
            --done_file "{output.done}" 2>&1 | tee "{log}"
        """


# --- 5.5 Batch dbCAN (Concatenated FAA, 4 CPUs) ---
rule mag_traits_batch_dbcan:
    input:
        prodigal_done = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "prodigal.done")
    output:
        overview = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "dbcan", "overview.txt")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "batches", "{batch_id}", "dbcan.log")
    threads: get_threads(MT_C.get("dbcan", {}).get("threads", 4))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = mt_env("dbcan"),
        db = mt_full_path(MT_C.get("dbcan", {}).get("db") or config.get("databases", {}).get("dbcan", "")),
        tools_arg = (
            f"--tools {' '.join(MT_C.get('dbcan', {}).get('tools'))}"
            if (isinstance(MT_C.get("dbcan", {}).get("tools"), list) and MT_C.get("dbcan", {}).get("tools"))
            else (f"--tools {MT_C.get('dbcan', {}).get('tools')}" if MT_C.get("dbcan", {}).get("tools") else "")
        ),
        extra = MT_C.get("dbcan", {}).get("extra_flags", ""),
        batch_dir = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}")
    wildcard_constraints:
        sample = r"[^/]+",
        batch_id = r"batch_\d+"
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.batch_dir}/dbcan" "$(dirname "{log}")"
        CONCAT_FAA="{params.batch_dir}/batch_concat_dbcan.faa"
        > "$CONCAT_FAA"

        for faa_file in "{params.batch_dir}"/faa/*.faa; do
            [ ! -f "$faa_file" ] && continue
            bname="$(basename "$faa_file")"
            mag="${{bname%.faa}}"
            awk -v mag="$mag" '/^>/ {{print ">" mag "___" substr($0, 2); next}} {{print}}' "$faa_file" >> "$CONCAT_FAA"
        done

        CMD_PREFIX=""
        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
            CMD_PREFIX="conda run -p {params.env}"
        fi

        DB_ARG=""
        if [ -n "{params.db}" ] && [ -d "{params.db}" ]; then
            DB_ARG="--db_dir {params.db}"
        fi

        if [ -n "$DB_ARG" ] && [ -s "$CONCAT_FAA" ]; then
            if $CMD_PREFIX run_dbcan --help 2>&1 | grep -q "CAZyme_annotation"; then
                $CMD_PREFIX run_dbcan CAZyme_annotation \
                    --input_raw_data "$CONCAT_FAA" \
                    --mode protein \
                    --output_dir "{params.batch_dir}/dbcan" \
                    $DB_ARG \
                    --threads {threads} \
                    {params.extra} 2>&1 | tee "{log}" || true
            else
                DBCAN_EXTRA="{params.extra}"
                if [[ "$DBCAN_EXTRA" != *"--dbCANFile"* ]] && [ -n "{params.db}" ] && [ -d "{params.db}" ]; then
                    if [ -f "{params.db}/dbCAN-HMMdb-V12.txt" ]; then
                        DBCAN_EXTRA="$DBCAN_EXTRA --dbCANFile dbCAN-HMMdb-V12.txt"
                    elif [ -f "{params.db}/dbCAN-HMMdb-V11.txt" ]; then
                        DBCAN_EXTRA="$DBCAN_EXTRA --dbCANFile dbCAN-HMMdb-V11.txt"
                    elif [ -f "{params.db}/dbCAN.hmm" ]; then
                        DBCAN_EXTRA="$DBCAN_EXTRA --dbCANFile dbCAN.hmm"
                    fi
                fi

                $CMD_PREFIX run_dbcan "$CONCAT_FAA" protein \
                    --out_dir "{params.batch_dir}/dbcan" \
                    $DB_ARG \
                    {params.tools_arg} \
                    $DBCAN_EXTRA \
                    --dbcan_thread {threads} \
                    --hmm_cpu {threads} \
                    --dia_cpu {threads} 2>&1 | tee "{log}" || true
            fi
        else
            echo "INFO: dbCAN database not configured or concatenated FAA is empty." >> "{log}"
        fi

        if [ -f "{params.batch_dir}/dbcan/overview.tsv" ] && [ -s "{params.batch_dir}/dbcan/overview.tsv" ]; then
            cp -f "{params.batch_dir}/dbcan/overview.tsv" "{output.overview}"
        elif [ -f "{params.batch_dir}/dbcan/overview.txt" ] && [ -s "{params.batch_dir}/dbcan/overview.txt" ]; then
            cp -f "{params.batch_dir}/dbcan/overview.txt" "{output.overview}"
        fi

        if [ ! -f "{output.overview}" ] || [ ! -s "{output.overview}" ]; then
            echo -e "Gene ID\\tHMMER\\tDiamond\\tdbCAN_sub\\t#ofTools" > "{output.overview}"
        fi

        rm -f "$CONCAT_FAA"
        """


# --- 5.6 Batch CAMPER (4 CPUs) ---
rule mag_traits_batch_camper:
    input:
        prodigal_done = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "prodigal.done")
    output:
        tsv = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "camper", "camper_hits.tsv")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "batches", "{batch_id}", "camper.log")
    threads: get_threads(MT_C.get("camper", {}).get("threads", 4))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = mt_env("camper"),
        hmmer_env = mt_env("hmmer") or mt_env("bakta"),
        hmm_db = mt_full_path(MT_C.get("camper", {}).get("hmm_db") or config.get("databases", {}).get("camper", "")),
        evalue = MT_C.get("camper", {}).get("evalue", 1e-5),
        batch_dir = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}")
    wildcard_constraints:
        sample = r"[^/]+",
        batch_id = r"batch_\d+"
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.batch_dir}/camper" "$(dirname "{log}")"
        CONCAT_FAA="{params.batch_dir}/camper/batch_concat.faa"
        > "$CONCAT_FAA"

        for faa_file in "{params.batch_dir}"/faa/*.faa; do
            [ ! -f "$faa_file" ] && continue
            bname="$(basename "$faa_file")"
            mag="${{bname%.faa}}"
            awk -v mag="$mag" '/^>/ {{print ">" mag "___" substr($0, 2); next}} {{print}}' "$faa_file" >> "$CONCAT_FAA"
        done

        CMD_PREFIX=""
        if [ -n "{params.env}" ] && [ -d "{params.env}" ] && [ -x "{params.env}/bin/hmmsearch" ]; then
            CMD_PREFIX="conda run -p {params.env}"
        elif [ -n "{params.hmmer_env}" ] && [ -d "{params.hmmer_env}" ]; then
            CMD_PREFIX="conda run -p {params.hmmer_env}"
        elif [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
            CMD_PREFIX="conda run -p {params.env}"
        fi

        echo -e "MAG\\tGene_ID\\tTarget_HMM\\tE_value\\tScore" > "{output.tsv}"

        if [ -n "{params.hmm_db}" ] && [ -e "{params.hmm_db}" ] && [ -s "$CONCAT_FAA" ]; then
            TBL_OUT="{params.batch_dir}/camper/camper.tbl"
            if [ -d "{params.hmm_db}" ]; then
                TMP_HMM="$(mktemp --suffix=.hmm)"
                find "{params.hmm_db}" -maxdepth 2 -name "*.hmm" -type f -exec cat {{}} + > "$TMP_HMM" 2>/dev/null || touch "$TMP_HMM"
                if [ ! -s "$TMP_HMM" ]; then
                    cat "{params.hmm_db}"/*.hmm > "$TMP_HMM" 2>/dev/null || touch "$TMP_HMM"
                fi
                $CMD_PREFIX hmmsearch --cpu {threads} -E {params.evalue} --tblout "$TBL_OUT" "$TMP_HMM" "$CONCAT_FAA" 2>&1 | tee "{log}" || true
                rm -f "$TMP_HMM"
            else
                $CMD_PREFIX hmmsearch --cpu {threads} -E {params.evalue} --tblout "$TBL_OUT" "{params.hmm_db}" "$CONCAT_FAA" 2>&1 | tee "{log}" || true
            fi

            if [ -f "$TBL_OUT" ]; then
                awk '!/^#/ && NF>=5 {{
                    n = split($1, a, "___");
                    mag = (n >= 2 ? a[1] : "Unknown");
                    gene = (n >= 2 ? a[2] : $1);
                    print mag "\\t" gene "\\t" $3 "\\t" $5 "\\t" $6
                }}' "$TBL_OUT" >> "{output.tsv}" || true
            fi
        else
            echo "INFO: CAMPER database missing or empty proteins." >> "{log}"
        fi
        """


# --- 5.7 Batch Ribosomal Protein Identification & gRodon2 Growth Estimation ---
rule mag_traits_batch_ribosomal:
    input:
        prodigal_done = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "prodigal.done")
    output:
        tsv = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "grodon2", "ribosomal_genes.tsv")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "batches", "{batch_id}", "ribosomal.log")
    threads: get_threads(MT_C.get("grodon2", {}).get("ribo_threads", 4))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = mt_env("camper") or mt_env("prodigal") or mt_env("drep") or mt_env("dbcan"),
        hmmer_env = mt_env("hmmer") or mt_env("bakta"),
        ribo_db = mt_full_path(MT_C.get("grodon2", {}).get("ribo_db") or config.get("databases", {}).get("checkm", "")),
        evalue = MT_C.get("grodon2", {}).get("ribo_evalue", 1e-5),
        batch_dir = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}")
    wildcard_constraints:
        sample = r"[^/]+",
        batch_id = r"batch_\d+"
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.batch_dir}/grodon2" "$(dirname "{log}")"
        echo -e "MAG\\tGene_ID\\tTarget_HMM\\tE_value\\tScore" > "{output.tsv}"

        CONCAT_FAA="{params.batch_dir}/grodon2/batch_concat.faa"
        > "$CONCAT_FAA"

        for faa_file in "{params.batch_dir}"/faa/*.faa; do
            [ ! -f "$faa_file" ] && continue
            bname="$(basename "$faa_file")"
            mag="${{bname%.faa}}"
            awk -v mag="$mag" '/^>/ {{print ">" mag "___" substr($0, 2); next}} {{print}}' "$faa_file" >> "$CONCAT_FAA"
        done

        CMD_PREFIX=""
        if [ -n "{params.env}" ] && [ -d "{params.env}" ] && [ -x "{params.env}/bin/hmmsearch" ]; then
            CMD_PREFIX="conda run -p {params.env}"
        elif [ -n "{params.hmmer_env}" ] && [ -d "{params.hmmer_env}" ] && [ -x "{params.hmmer_env}/bin/hmmsearch" ]; then
            CMD_PREFIX="conda run -p {params.hmmer_env}"
        elif command -v hmmsearch >/dev/null 2>&1; then
            CMD_PREFIX=""
        elif [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
            CMD_PREFIX="conda run -p {params.env}"
        fi

        HMM_FILE=""
        if [ -n "{params.ribo_db}" ] && [ -e "{params.ribo_db}" ]; then
            if [ -f "{params.ribo_db}" ]; then
                HMM_FILE="{params.ribo_db}"
            elif [ -d "{params.ribo_db}" ]; then
                if [ -f "{params.ribo_db}/hmms/phylo.hmm" ]; then
                    HMM_FILE="{params.ribo_db}/hmms/phylo.hmm"
                elif [ -f "{params.ribo_db}/phylo.hmm" ]; then
                    HMM_FILE="{params.ribo_db}/phylo.hmm"
                elif [ -d "{params.ribo_db}/hmms" ]; then
                    FIRST_HMM="$(find "{params.ribo_db}/hmms" -maxdepth 1 -name "*.hmm" | head -n 1 || true)"
                    if [ -n "$FIRST_HMM" ]; then
                        COUNT="$(find "{params.ribo_db}/hmms" -maxdepth 1 -name "*.hmm" | wc -l)"
                        if [ "$COUNT" -eq 1 ]; then
                            HMM_FILE="$FIRST_HMM"
                        else
                            TMP_HMM="$(mktemp --suffix=.hmm)"
                            cat "{params.ribo_db}/hmms"/*.hmm > "$TMP_HMM" 2>/dev/null || touch "$TMP_HMM"
                            if [ -s "$TMP_HMM" ]; then
                                HMM_FILE="$TMP_HMM"
                            fi
                        fi
                    fi
                elif [ -d "{params.ribo_db}/pfam" ]; then
                    PFAM_HMM="$(find "{params.ribo_db}/pfam" -maxdepth 1 -name "*.hmm" | head -n 1 || true)"
                    if [ -n "$PFAM_HMM" ]; then
                        HMM_FILE="$PFAM_HMM"
                    fi
                fi
            fi
        fi

        if [ -n "$HMM_FILE" ] && [ -s "$HMM_FILE" ] && [ -s "$CONCAT_FAA" ]; then
            TBL_OUT="{params.batch_dir}/grodon2/ribosomal.tbl"
            $CMD_PREFIX hmmsearch --cpu {threads} -E {params.evalue} --tblout "$TBL_OUT" "$HMM_FILE" "$CONCAT_FAA" 2>&1 | tee "{log}" || true

            if [ -f "$TBL_OUT" ]; then
                awk '!/^#/ && NF>=5 {{
                    n = split($1, a, "___");
                    mag = (n >= 2 ? a[1] : "Unknown");
                    gene = (n >= 2 ? a[2] : $1);
                    print mag "\\t" gene "\\t" $3 "\\t" $5 "\\t" $6
                }}' "$TBL_OUT" >> "{output.tsv}" || true
            fi
            if [[ "$HMM_FILE" == *"tmp."* ]]; then
                rm -f "$HMM_FILE"
            fi
        else
            echo "INFO: Ribosomal HMM database missing or empty protein input." >> "{log}"
        fi

        rm -f "$CONCAT_FAA"
        """


rule mag_traits_batch_grodon2:
    input:
        prodigal_done = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "prodigal.done"),
        ribo_tsv = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "grodon2", "ribosomal_genes.tsv")
    output:
        tsv = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "grodon2", "grodon2_growth.tsv")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "batches", "{batch_id}", "grodon2.log")
    threads: get_threads(MT_C.get("grodon2", {}).get("threads", 1))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = mt_env("grodon2"),
        temp = MT_C.get("grodon2", {}).get("temperature", MT_C.get("phydon", {}).get("temperature", 20.0)),
        mode = MT_C.get("grodon2", {}).get("mode", "partial"),
        script = os.path.join(WORKFLOW_ROOT, "scripts", "run_grodon2.R"),
        batch_dir = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}")
    wildcard_constraints:
        sample = r"[^/]+",
        batch_id = r"batch_\d+"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.tsv}")" "$(dirname "{log}")"

        RIBO_ARG=""
        if [ -s "{input.ribo_tsv}" ]; then
            RIBO_ARG="--ribo_file {input.ribo_tsv}"
        fi

        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
            conda run -p "{params.env}" Rscript "{params.script}" \
                --batch_dir "{params.batch_dir}/ffn" \
                --out "{output.tsv}" \
                --temp "{params.temp}" \
                --mode "{params.mode}" \
                $RIBO_ARG 2>&1 | tee "{log}" || true
        else
            Rscript "{params.script}" \
                --batch_dir "{params.batch_dir}/ffn" \
                --out "{output.tsv}" \
                --temp "{params.temp}" \
                --mode "{params.mode}" \
                $RIBO_ARG 2>&1 | tee "{log}" || true
        fi

        if [ ! -f "{output.tsv}" ]; then
            echo -e "MAG\\tDoublingTime_hours\\tMaxGrowthRate_per_hour\\tCUB_HEG\\tGrowthCategory\\tMethod\\tMode\\tTemperature\\tStatus" > "{output.tsv}"
        fi
        """


# --- 5.8 Batch GenomeSPOT Phenotypic Trait Prediction (2 CPUs) ---
rule mag_traits_batch_genomespot:
    input:
        prodigal_done = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "prodigal.done"),
        mag_list = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "mag_list.tsv")
    output:
        summary = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}", "genomespot", "genomespot_summary.tsv")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "batches", "{batch_id}", "genomespot.log")
    threads: get_threads(MT_C.get("genomespot", {}).get("threads", 2))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = mt_env("genomespot"),
        models = mt_full_path(MT_C.get("genomespot", {}).get("models_dir") or config.get("databases", {}).get("genomespot", "")),
        save_feat = "--save_features" if MT_C.get("genomespot", {}).get("save_features", False) else "",
        script = os.path.join(WORKFLOW_ROOT, "scripts", "run_genomespot.py"),
        batch_dir = os.path.join(MT_BASE_RES, "{sample}", "batches", "{batch_id}")
    wildcard_constraints:
        sample = r"[^/]+",
        batch_id = r"batch_\d+"
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.batch_dir}/genomespot" "$(dirname "{log}")"

        CMD_PREFIX=""
        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
            CMD_PREFIX="conda run -p {params.env}"
        fi

        MODELS_ARG=""
        if [ -n "{params.models}" ]; then
            MODELS_ARG="--models {params.models}"
        fi

        $CMD_PREFIX python3 "{params.script}" \
            --mag_list "{input.mag_list}" \
            --faa_dir "{params.batch_dir}/faa" \
            --out_dir "{params.batch_dir}/genomespot" \
            --threads {threads} \
            {params.save_feat} \
            $MODELS_ARG 2>&1 | tee "{log}" || true

        if [ ! -f "{output.summary}" ]; then
            echo -e "MAG\\tOGT\\tTmin\\tTmax\\tpH_opt\\tpH_min\\tpH_max\\tSalinity_opt\\tSalinity_min\\tSalinity_max\\tOxygen_tolerance\\tOxygen_prob\\tStatus" > "{output.summary}"
        fi
        """


# --- 5.9 Summary Aggregation Across Batches & dRep per Sample ---
def get_all_batch_outputs(wildcards):
    s = wildcards.sample
    batch_list_file = checkpoints.mag_traits_partition_batches.get(sample=s).output.batch_list
    batches = []
    if os.path.exists(batch_list_file):
        with open(batch_list_file, "r", encoding="utf-8") as f:
            for line in f:
                b = line.strip()
                if b:
                    batches.append(b)

    outputs = []
    for b in batches:
        if do_dbcan:
            outputs.append(os.path.join(MT_BASE_RES, s, "batches", b, "dbcan", "overview.txt"))
        if do_camper:
            outputs.append(os.path.join(MT_BASE_RES, s, "batches", b, "camper", "camper_hits.tsv"))
        if do_grodon2:
            outputs.append(os.path.join(MT_BASE_RES, s, "batches", b, "grodon2", "grodon2_growth.tsv"))
        if do_genomespot:
            outputs.append(os.path.join(MT_BASE_RES, s, "batches", b, "genomespot", "genomespot_summary.tsv"))
    return outputs


def get_summary_cdb_r1(wildcards):
    if do_drep and two_rounds_drep:
        return os.path.join(MT_BASE_RES, wildcards.sample, "drep", "round1", "data_tables", "Cdb.csv")
    return []


def get_summary_cdb_r2(wildcards):
    if do_drep and two_rounds_drep:
        return os.path.join(MT_BASE_RES, wildcards.sample, "drep", "round2", "data_tables", "Cdb.csv")
    return []


def get_summary_cdb_single(wildcards):
    if do_drep and not two_rounds_drep:
        if drep_block_size > 0:
            return os.path.join(MT_BASE_RES, wildcards.sample, "drep", "round1", "data_tables", "Cdb.csv")
        else:
            return os.path.join(MT_BASE_RES, wildcards.sample, "drep", "single", "data_tables", "Cdb.csv")
    return []


def get_summary_pool_done(wildcards):
    if do_prodigal and do_pool_genes:
        return os.path.join(MT_BASE_RES, wildcards.sample, "proteins", ".pool.done")
    return []


rule mag_traits_summary:
    input:
        batch_outputs = get_all_batch_outputs,
        pool_done = get_summary_pool_done,
        cdb_r1 = get_summary_cdb_r1,
        cdb_r2 = get_summary_cdb_r2,
        cdb = get_summary_cdb_single
    output:
        matrix = os.path.join(MT_BASE_RES, "{sample}", "summary", "all_mag_traits_matrix.tsv")
    log:
        os.path.join(MT_BASE_LOGS, "{sample}", "summary", "merge_traits.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = mt_env("python"),
        res_dir = os.path.join(MT_BASE_RES, "{sample}"),
        script = os.path.join(WORKFLOW_ROOT, "scripts", "merge_mag_traits.py"),
        cdb_arg = lambda wildcards, input: (
            f'--cdb-r1 "{input.cdb_r1}" --cdb-r2 "{input.cdb_r2}"'
            if (do_drep and two_rounds_drep and input.cdb_r1 and input.cdb_r2)
            else (f'--cdb "{input.cdb}"' if (do_drep and input.cdb) else "")
        )
    wildcard_constraints:
        sample = r"[^/]+"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.matrix}")" "$(dirname "{log}")"

        CMD_PREFIX=""
        if [ -n "{params.env}" ] && [ -d "{params.env}" ]; then
            CMD_PREFIX="conda run -p {params.env}"
        fi

        $CMD_PREFIX python3 "{params.script}" \
            --results-dir "{params.res_dir}" \
            --output "{output.matrix}" \
            {params.cdb_arg} 2>&1 | tee "{log}"
        """
