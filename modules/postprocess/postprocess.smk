# Reusable downstream integration of batch annotation tables.
#
# Current pipeline:
#   batch functional summary -> per-genome gene counts -> OXPHOS evaluation

import os
import re


PP_C = config.get("postprocess", {}) or {}
PP_TASKS = PP_C.get("tasks", {}) or {}
PP_ROOT = config.get("project_root", ".")
PP_RES = os.path.join(PP_ROOT, PP_C.get("results_dir", "results/postprocess"))
PP_LOGS = os.path.join(PP_ROOT, PP_C.get("logs_dir", "logs/postprocess"))
PP_ENV_ROOT = config.get("conda_env_dir", "")


def pp_full_path(value):
    if not value:
        return ""
    if os.path.isabs(value):
        return value
    return os.path.normpath(os.path.join(PP_ROOT, value))


def pp_task(wildcards):
    task = PP_TASKS.get(wildcards.task)
    if not isinstance(task, dict):
        raise ValueError(f"Unknown postprocess task: {wildcards.task}")
    return task


def pp_batch_summary(wildcards):
    task = pp_task(wildcards)
    explicit = task.get("functional_summary")
    if explicit:
        return pp_full_path(explicit)
    batch = task.get("annotation_batch", wildcards.task)
    annotation_results = config.get("annotation", {}).get(
        "external_results_dir", "results/external_annotation"
    )
    return os.path.join(pp_full_path(annotation_results), f"{batch}_functional_summary.tsv")


def pp_gene_map(wildcards):
    value = pp_task(wildcards).get("gene_map")
    if not value:
        raise ValueError(f"postprocess.tasks.{wildcards.task}.gene_map is required")
    return pp_full_path(value)


def pp_python_env():
    env_name = PP_C.get("conda_envs", {}).get("python", "biopython_1.85")
    if not PP_ENV_ROOT or os.path.isabs(env_name):
        return env_name
    return os.path.join(PP_ENV_ROOT, env_name)


_pp_task_names = sorted(PP_TASKS)
wildcard_constraints:
    task = "|".join(re.escape(name) for name in _pp_task_names) if _pp_task_names else "none_matched_placeholder"


for _task_name, _task_info in PP_TASKS.items():
    if not isinstance(_task_info, dict):
        continue
    TARGETS.append(os.path.join(PP_RES, _task_name, "gene_counts.tsv"))
    if str(_task_info.get("evaluate_oxphos", False)).lower() == "true":
        TARGETS.append(os.path.join(PP_RES, _task_name, "oxphos_evaluation.tsv"))


rule postprocess_gene_counts:
    input:
        summary = pp_batch_summary,
        gene_map = pp_gene_map
    output:
        matrix = os.path.join(PP_RES, "{task}/gene_counts.tsv")
    log:
        os.path.join(PP_LOGS, "{task}/gene_counts.log")
    threads: get_threads(2)
    resources:
        runtime = "1d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = pp_python_env(),
        script = os.path.join(WORKFLOW_ROOT, "map_gene_counts.py")
    shell:
        r"""
        {{
            set -euo pipefail
            mkdir -p "$(dirname "{output.matrix}")" "$(dirname "{log}")"
            PYTHON_BIN="python3"
            if [ -x "{params.env}/bin/python3" ]; then
                PYTHON_BIN="{params.env}/bin/python3"
            elif [ -x "{params.env}/bin/python" ]; then
                PYTHON_BIN="{params.env}/bin/python"
            fi
            "$PYTHON_BIN" "{params.script}" \
                --summary "{input.summary}" \
                --gene "{input.gene_map}" \
                --output "{output.matrix}"
        }} 2>&1 | tee "{log}"
        """


rule postprocess_oxphos_evaluation:
    input:
        matrix = os.path.join(PP_RES, "{task}/gene_counts.tsv")
    output:
        evaluation = os.path.join(PP_RES, "{task}/oxphos_evaluation.tsv")
    log:
        os.path.join(PP_LOGS, "{task}/oxphos_evaluation.log")
    threads: get_threads(2)
    resources:
        runtime = "1d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = pp_python_env(),
        script = os.path.join(WORKFLOW_ROOT, "OP_evaluate.py")
    shell:
        r"""
        {{
            set -euo pipefail
            mkdir -p "$(dirname "{output.evaluation}")" "$(dirname "{log}")"
            PYTHON_BIN="python3"
            if [ -x "{params.env}/bin/python3" ]; then
                PYTHON_BIN="{params.env}/bin/python3"
            elif [ -x "{params.env}/bin/python" ]; then
                PYTHON_BIN="{params.env}/bin/python"
            fi
            "$PYTHON_BIN" "{params.script}" \
                "{input.matrix}" "{output.evaluation}"
        }} 2>&1 | tee "{log}"
        """
