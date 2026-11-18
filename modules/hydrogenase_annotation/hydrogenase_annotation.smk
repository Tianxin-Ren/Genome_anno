# ==============================================================================
# Hydrogenase Functional Annotation Module (34-Class Mechanistic Engine)
# ==============================================================================
import csv
import os


HYD_C = config.get("hydrogenase_annotation", {})
HYD_ROOT = config.get("project_root", "") or config.get("paths", {}).get("project_root", "")
if HYD_ROOT:
    HYD_ROOT = os.path.expanduser(os.path.expandvars(str(HYD_ROOT)))
    if not os.path.isabs(HYD_ROOT):
        HYD_ROOT = os.path.abspath(HYD_ROOT)
else:
    HYD_ROOT = os.getcwd()


def hyd_full_path(path):
    if not path:
        return ""
    path_str = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path_str):
        return os.path.normpath(path_str)
    if HYD_ROOT:
        return os.path.normpath(os.path.join(HYD_ROOT, path_str))
    return os.path.normpath(path_str)


def hyd_active_samples():
    cli_val = config.get("sample") or config.get("sample_name") or HYD_C.get("sample") or HYD_C.get("sample_name")
    if isinstance(cli_val, str) and cli_val.strip() and cli_val.strip().lower() not in ["all", "*"]:
        return [cli_val.strip()]
    for s_cand in [config.get("samples"), config.get("sample"), HYD_C.get("samples"), HYD_C.get("sample")]:
        if isinstance(s_cand, dict) and s_cand:
            return list(s_cand.keys())
        if isinstance(s_cand, list) and s_cand:
            return [str(x).strip() for x in s_cand if str(x).strip()]
        if isinstance(s_cand, str) and s_cand.strip() and s_cand.strip().lower() not in ["all", "*"]:
            return [s_cand.strip()]
    return ["default"]


def hyd_sample_info(sample_tag):
    for s_dict in [config.get("samples"), config.get("sample"), HYD_C.get("samples"), HYD_C.get("sample")]:
        if isinstance(s_dict, dict) and sample_tag in s_dict:
            val = s_dict[sample_tag]
            if isinstance(val, dict):
                return val
            if isinstance(val, str):
                return {"genome_dir": val}
    return {}


def hyd_resolve_pipeline_dir():
    p = HYD_C.get("pipeline_dir", "tools/hydrogenase_annotation")
    if os.path.isabs(p):
        return p
    return os.path.normpath(os.path.join(WORKFLOW_ROOT, p))


def hyd_env_target():
    env_name = config.get("environments", {}).get("hydrogenase") or HYD_C.get("conda_envs", {}).get("hydrogenase", "")
    if env_name:
        env_dir = config.get("conda_env_dir", "") or config.get("paths", {}).get("conda_envs", "")
        return env_name if os.path.isabs(env_name) else (os.path.join(env_dir, env_name) if env_dir else env_name)
    local_env = os.path.join(hyd_resolve_pipeline_dir(), ".env")
    if os.path.isdir(local_env):
        return local_env
    return ""


HYD_ACTIVE_SAMPLES = hyd_active_samples()
HYD_BASE_RES = hyd_full_path(HYD_C.get("results_dir", "results/hydrogenase_annotation"))
HYD_BASE_LOGS = hyd_full_path(HYD_C.get("logs_dir", "logs/hydrogenase_annotation"))
HYD_PIPELINE_DIR = hyd_resolve_pipeline_dir()
HYD_MERGE_SCRIPT = os.path.join(WORKFLOW_ROOT, "modules", "hydrogenase_annotation", "scripts", "merge_hydrogenase_batches.py")

HYD_BATCH_SIZE = max(1, int(HYD_C.get("batch_size", 200)))
HYD_THREADS_PER_BATCH = max(1, int(HYD_C.get("threads_per_batch", 40)))
HYD_THREADS_PER_GENOME = max(1, int(HYD_C.get("threads_per_genome", 1)))
HYD_MODE = str(HYD_C.get("mode", "genome"))
HYD_TOPOLOGY = str(HYD_C.get("topology", "auto"))
HYD_TRANS_TABLE = int(HYD_C.get("translation_table", 11))
HYD_RESUME = str(HYD_C.get("resume", True)).lower() == "true"
HYD_INTEGRATE_MA = str(HYD_C.get("integrate_mag_annotation", True)).lower() == "true"
HYD_REUSE_MA_GENES = str(HYD_C.get("reuse_mag_gene_calling", True)).lower() == "true"

_MA_ENABLED = bool(config.get("modules", {}).get("mag_annotation", {}).get("enable", False))
_MA_MERGE = str(config.get("mag_annotation", {}).get("steps", {}).get("merge_summary", True)).lower() == "true"
_MA_BASE_RES = hyd_full_path(config.get("mag_annotation", {}).get("results_dir", "results/mag_annotation"))
_MA_GENE_BATCHES = globals().get("GENE_BATCHES", {})

HYD_BATCHES = {}
HYD_ANNO_TARGETS = []

for s in HYD_ACTIVE_SAMPLES:
    s_info = hyd_sample_info(s)
    g_dir = s_info.get("genome_dir") or s_info.get("path") or config.get("genome_dir") or HYD_C.get("genome_dir", "")
    HYD_BATCHES[s] = []
    if not g_dir:
        continue
    abs_gdir = hyd_full_path(g_dir)
    if not os.path.isdir(abs_gdir):
        continue

    genomes = {}
    valid_exts = (".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz", ".fa.gz")
    for root, _, files in os.walk(abs_gdir):
        for fname in sorted(files):
            if fname.endswith((".faa", ".faa.gz", ".ffn", ".ffn.gz")):
                continue
            for ext in valid_exts:
                if fname.endswith(ext):
                    stem = fname[:-len(ext)]
                    genomes.setdefault(stem, os.path.join(root, fname))
                    break

    sorted_items = sorted(genomes.items())
    if not sorted_items:
        continue

    for idx in range(0, len(sorted_items), HYD_BATCH_SIZE):
        bid = f"{(idx // HYD_BATCH_SIZE) + 1:04d}"
        chunk = sorted_items[idx:idx + HYD_BATCH_SIZE]
        HYD_BATCHES[s].append(bid)
        manifest_path = os.path.join(HYD_BASE_RES, s, "batches", f"batch_{bid}", "manifest.tsv")
        os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
        if HYD_REUSE_MA_GENES and _MA_ENABLED:
            p_base = hyd_full_path(s_info.get("protein_dir") or "") or os.path.join(_MA_BASE_RES, s, "proteins", "all")
            g_base = os.path.join(_MA_BASE_RES, s, "genes")
            lines = ["sample\tproteins\tgff\tpath\n"] + [
                f"{stem}\t{os.path.join(p_base, stem + '.faa.gz')}\t{os.path.join(g_base, stem + '.gff.gz')}\t{fpath}\n"
                for stem, fpath in chunk
            ]
        else:
            lines = ["sample\tpath\n"] + [f"{stem}\t{fpath}\n" for stem, fpath in chunk]
        new_text = "".join(lines)
        if not os.path.isfile(manifest_path) or open(manifest_path, encoding="utf-8").read() != new_text:
            with open(manifest_path, "w", encoding="utf-8") as fh:
                fh.write(new_text)

    HYD_ANNO_TARGETS.append(os.path.join(HYD_BASE_RES, s, "cohort_summary.tsv"))
    HYD_ANNO_TARGETS.append(os.path.join(HYD_BASE_RES, s, "cohort_functional_annotations.tsv"))
    HYD_ANNO_TARGETS.append(os.path.join(HYD_BASE_RES, s, "cohort_class_matrix.tsv"))


def hyd_batch_inputs(wildcards):
    s = wildcards.sample
    inp = {
        "manifest": os.path.join(HYD_BASE_RES, s, "batches", f"batch_{wildcards.bid}", "manifest.tsv"),
    }
    if HYD_REUSE_MA_GENES and _MA_ENABLED and _MA_GENE_BATCHES.get(s):
        inp["gene_calls"] = [
            os.path.join(_MA_BASE_RES, s, "gene_calling_batches", f"{gb_id}.complete.json")
            for gb_id in sorted(_MA_GENE_BATCHES[s])
        ]
    return inp


rule hydrogenase_annotation_batch:
    input:
        unpack(hyd_batch_inputs)
    output:
        summary = os.path.join(HYD_BASE_RES, "{sample}", "batches", "batch_{bid}", "cohort_summary.tsv"),
        done = os.path.join(HYD_BASE_RES, "{sample}", "batches", "batch_{bid}", "batch.done")
    log:
        os.path.join(HYD_BASE_LOGS, "{sample}", "batch_{bid}.log")
    threads: get_threads(HYD_THREADS_PER_BATCH)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = hyd_env_target(),
        runner = os.path.join(HYD_PIPELINE_DIR, "run_parallel.py"),
        out_dir = lambda wildcards: os.path.join(HYD_BASE_RES, wildcards.sample, "batches", f"batch_{wildcards.bid}"),
        tpw = HYD_THREADS_PER_GENOME,
        mode = HYD_MODE,
        topology = HYD_TOPOLOGY,
        table = HYD_TRANS_TABLE,
        resume_flag = "--resume" if HYD_RESUME else ""
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        ENV_DIR="{params.env}"
        PY_BIN="python3"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                export PATH="$ENV_DIR/bin:$PATH"
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                export PATH="$ENV_DIR/bin:$PATH"
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi
        if ! $PY_BIN -c "import Bio, openpyxl" 2>/dev/null; then
            echo "[hydrogenase_annotation ERROR] Python environment '$ENV_DIR' is missing or lacks Bio/openpyxl. Create it with: mamba env create -p '$ENV_DIR' -f envs/hydrogenase_annotation.yaml" >&2
            exit 1
        fi
        $PY_BIN -u "{params.runner}" \
            --manifest "{input.manifest}" \
            --out "{params.out_dir}" \
            --workers {threads} \
            --threads-per-worker {params.tpw} \
            --mode "{params.mode}" \
            --topology "{params.topology}" \
            --translation-table {params.table} \
            {params.resume_flag} 2>&1 | tee "{log}"
        touch "{output.done}"
        """


def hyd_merge_inputs(wildcards):
    s = wildcards.sample
    inp = {
        "done": [os.path.join(HYD_BASE_RES, s, "batches", f"batch_{bid}", "batch.done") for bid in HYD_BATCHES.get(s, [])],
    }
    if HYD_INTEGRATE_MA and _MA_ENABLED and _MA_MERGE:
        inp["ko_best"] = os.path.join(_MA_BASE_RES, s, "merged", "protein_KO_best.tsv")
    return inp


rule hydrogenase_annotation_merge:
    input:
        unpack(hyd_merge_inputs)
    output:
        summary = os.path.join(HYD_BASE_RES, "{sample}", "cohort_summary.tsv"),
        annotations = os.path.join(HYD_BASE_RES, "{sample}", "cohort_functional_annotations.tsv"),
        matrix = os.path.join(HYD_BASE_RES, "{sample}", "cohort_class_matrix.tsv"),
        manifest = os.path.join(HYD_BASE_RES, "{sample}", "cohort_manifest.json")
    log:
        os.path.join(HYD_BASE_LOGS, "{sample}", "merge.log")
    threads: get_threads(4)
    resources:
        runtime = "1d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = hyd_env_target(),
        script = HYD_MERGE_SCRIPT,
        out_dir = lambda wildcards: os.path.join(HYD_BASE_RES, wildcards.sample),
        batch_dirs = lambda wildcards: [os.path.join(HYD_BASE_RES, wildcards.sample, "batches", f"batch_{bid}") for bid in HYD_BATCHES.get(wildcards.sample, [])],
        mag_dir_arg = lambda wildcards: f'--mag-sample-dir "{os.path.join(_MA_BASE_RES, wildcards.sample)}"' if HYD_INTEGRATE_MA else ""
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname "{log}")"
        ENV_DIR="{params.env}"
        PY_BIN="python3"
        if [ -n "$ENV_DIR" ]; then
            if [ -x "$ENV_DIR/bin/python3" ]; then
                export PATH="$ENV_DIR/bin:$PATH"
                PY_BIN="$ENV_DIR/bin/python3"
            elif [ -x "$ENV_DIR/bin/python" ]; then
                export PATH="$ENV_DIR/bin:$PATH"
                PY_BIN="$ENV_DIR/bin/python"
            elif [[ "$ENV_DIR" == *"/"* ]]; then
                PY_BIN="conda run --no-capture-output -p $ENV_DIR python3"
            else
                PY_BIN="conda run --no-capture-output -n $ENV_DIR python3"
            fi
        fi
        $PY_BIN -u "{params.script}" \
            --batches {params.batch_dirs} \
            --out "{params.out_dir}" \
            {params.mag_dir_arg} 2>&1 | tee "{log}"
        """
