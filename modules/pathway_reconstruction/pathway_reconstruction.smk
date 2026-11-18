# ==============================================================================
# KEGG Module Reconstruction & Bioenergetic Capability Engine Module
# ==============================================================================
import os
import sys

# 1. Global Path Logic (Namespaced)
PR_C = config.get("pathway_reconstruction", {})
PR_ROOT = config.get("project_root", "") or config.get("paths", {}).get("project_root", "")
if PR_ROOT:
    PR_ROOT = os.path.expanduser(os.path.expandvars(str(PR_ROOT)))
    if not os.path.isabs(PR_ROOT):
        PR_ROOT = os.path.abspath(PR_ROOT)
else:
    PR_ROOT = os.getcwd()


def pr_full_path(path):
    if not path:
        return ""
    path_str = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path_str):
        return os.path.normpath(path_str)
    if PR_ROOT:
        return os.path.normpath(os.path.join(PR_ROOT, path_str))
    return os.path.normpath(path_str)


def resolve_module_catalog(path):
    if not path:
        return os.path.join(WORKFLOW_ROOT, "knowledge", "kegg", "all_kegg_modules.tsv")
    path_str = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path_str):
        if os.path.exists(path_str):
            return path_str
        wf_cand = os.path.join(WORKFLOW_ROOT, "knowledge", "kegg", os.path.basename(path_str))
        if os.path.exists(wf_cand):
            return wf_cand
        return path_str
    wf_path = os.path.normpath(os.path.join(WORKFLOW_ROOT, path_str))
    if os.path.exists(wf_path):
        return wf_path
    if PR_ROOT:
        pr_path = os.path.normpath(os.path.join(PR_ROOT, path_str))
        if os.path.exists(pr_path):
            return pr_path
    return wf_path


def get_pr_active_samples():
    cli_val = config.get("sample") or config.get("sample_name") or PR_C.get("sample") or PR_C.get("sample_name")
    if isinstance(cli_val, str) and cli_val.strip() and cli_val.strip().lower() not in ["all", "*"]:
        return [cli_val.strip()]

    for s_candidate in [
        config.get("samples"),
        config.get("sample"),
        PR_C.get("samples"),
        PR_C.get("sample"),
        config.get("mag_annotation", {}).get("samples"),
        config.get("mag_annotation", {}).get("sample"),
    ]:
        if isinstance(s_candidate, dict) and s_candidate:
            return list(s_candidate.keys())
        elif isinstance(s_candidate, list) and s_candidate:
            return [str(x).strip() for x in s_candidate if str(x).strip()]
        elif isinstance(s_candidate, str) and s_candidate.strip() and s_candidate.strip().lower() not in ["all", "*"]:
            return [s_candidate.strip()]

    return ["default"]


ACTIVE_PR_SAMPLES = get_pr_active_samples()
PR_BASE_RES = pr_full_path(PR_C.get("results_dir", "results/pathway_reconstruction"))
PR_BASE_LOGS = pr_full_path(PR_C.get("logs_dir", "logs/pathway_reconstruction"))

PR_SCRIPTS = os.path.join(WORKFLOW_ROOT, "modules", "pathway_reconstruction", "scripts")


def pr_env(name):
    env_name = config.get("environments", {}).get(name) or PR_C.get("conda_envs", {}).get(name)
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


# 2. Target Definitions
ENABLE_PW_FOOTPRINT = PR_C.get("enable_pathway_footprint", True)
ENABLE_BRITE = PR_C.get("enable_brite_profiling", True)
ENABLE_AUDIT = PR_C.get("enable_coverage_audit", True)

PATHWAY_TARGETS = []
for s in ACTIVE_PR_SAMPLES:
    # Tier 1 (v1 targets)
    PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "modules", "genome_module_completeness.tsv"))
    PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "modules", "genome_module_presence.tsv"))
    PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "modules", "module_evidence.tsv"))
    PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "modules", "genome_module_matrix.tsv"))
    PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "capabilities", "genome_metabolism_capabilities.json"))

    # Tier 2: Pathway Footprint
    if ENABLE_PW_FOOTPRINT:
        PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "pathways", "genome_pathway_footprint.tsv"))
        PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "pathways", "genome_pathway_matrix.tsv"))

    # Tier 3: BRITE System Profiling
    if ENABLE_BRITE:
        PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "brite", "genome_brite_counts.tsv"))
        PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "brite", "genome_brite_matrix.tsv"))

    # Global Coverage Audit
    if ENABLE_AUDIT:
        PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "audit", "genome_ko_coverage_audit.tsv"))

    # Tier 4 & Unified v2.0 Profile
    if ENABLE_PW_FOOTPRINT and ENABLE_BRITE and ENABLE_AUDIT:
        PATHWAY_TARGETS.append(os.path.join(PR_BASE_RES, s, "capabilities", "genome_functional_profile.json"))


# 3. Input Dependency Resolution from Annotation Module
def get_input_ko_strict(wildcards):
    custom_ko = PR_C.get("ko_strict_tsv")
    if custom_ko:
        return pr_full_path(custom_ko)
    # Default to mag_annotation output in the matching sample subfolder
    ma_res_base = pr_full_path(config.get("mag_annotation", {}).get("results_dir", "results/mag_annotation"))
    p = os.path.join(ma_res_base, wildcards.sample, "merged", "protein_KO_all_strict.tsv")
    if os.path.exists(p):
        return p
    if os.path.exists(p + ".gz"):
        return p + ".gz"
    return p



# 4. Rules
# ------------------------------------------------------------------------------

# 4.1 Evaluate KEGG Modules against Genome KO Evidence
rule pathway_evaluate_kegg_modules:
    input:
        ko_strict = get_input_ko_strict
    output:
        completeness = os.path.join(PR_BASE_RES, "{sample}", "modules", "genome_module_completeness.tsv"),
        presence = os.path.join(PR_BASE_RES, "{sample}", "modules", "genome_module_presence.tsv"),
        evidence = os.path.join(PR_BASE_RES, "{sample}", "modules", "module_evidence.tsv")
    log:
        os.path.join(PR_BASE_LOGS, "{sample}", "evaluate_kegg_modules.log")
    threads: get_threads(PR_C.get("threads", 4))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = pr_env("python"),
        script = os.path.join(PR_SCRIPTS, "evaluate_module.py"),
        modules_dir = resolve_module_catalog(PR_C.get("modules_catalog")),
        out_dir = lambda wildcards: os.path.join(PR_BASE_RES, wildcards.sample, "modules"),
        threshold = PR_C.get("completion_threshold", 0.75),
        max_missing = PR_C.get("max_missing_steps", 1),
        evidence_strategy = PR_C.get("ko_evidence_strategy", "all_strict")
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
            -k "{input.ko_strict}" \
            -d "{params.modules_dir}" \
            -o "{params.out_dir}" \
            --ko_evidence_strategy "{params.evidence_strategy}" \
            --completion_threshold {params.threshold} \
            --max_missing_steps {params.max_missing} 2>&1 | tee "{log}"
        """


# 4.2 Aggregate Module Matrices & Level 4 Bioenergetics Profile
rule pathway_generate_metabolism_capabilities:
    input:
        completeness = os.path.join(PR_BASE_RES, "{sample}", "modules", "genome_module_completeness.tsv")
    output:
        matrix = os.path.join(PR_BASE_RES, "{sample}", "modules", "genome_module_matrix.tsv"),
        capabilities = os.path.join(PR_BASE_RES, "{sample}", "capabilities", "genome_metabolism_capabilities.json")
    log:
        os.path.join(PR_BASE_LOGS, "{sample}", "generate_capabilities.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = pr_env("python"),
        script = os.path.join(PR_SCRIPTS, "generate_module_matrix.py"),
        modules_out = lambda wildcards: os.path.join(PR_BASE_RES, wildcards.sample, "modules"),
        cap_out = lambda wildcards: os.path.join(PR_BASE_RES, wildcards.sample, "capabilities")
    shell:
        """
        set -euo pipefail
        mkdir -p "{params.modules_out}" "{params.cap_out}" "$(dirname "{log}")"
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
            -c "{input.completeness}" \
            -o "{params.modules_out}" 2>&1 | tee "{log}"

        if [ -f "{params.modules_out}/genome_metabolism_capabilities.json" ]; then
            mv "{params.modules_out}/genome_metabolism_capabilities.json" "{output.capabilities}"
        fi
        """


# 4.3 Tier 2: Evaluate Pathway Map Footprints & Reaction Node Coverages
rule pathway_evaluate_pathway_footprint:
    input:
        ko_strict = get_input_ko_strict,
        completeness = os.path.join(PR_BASE_RES, "{sample}", "modules", "genome_module_completeness.tsv")
    output:
        footprint = os.path.join(PR_BASE_RES, "{sample}", "pathways", "genome_pathway_footprint.tsv"),
        matrix = os.path.join(PR_BASE_RES, "{sample}", "pathways", "genome_pathway_matrix.tsv")
    log:
        os.path.join(PR_BASE_LOGS, "{sample}", "evaluate_pathway_footprint.log")
    threads: get_threads(PR_C.get("threads", 4))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = pr_env("python"),
        script = os.path.join(PR_SCRIPTS, "evaluate_pathway_footprint.py"),
        pw_tsv = os.path.join(WORKFLOW_ROOT, PR_C.get("knowledge_dir", "knowledge/kegg"), "all_kegg_pathways.tsv"),
        triples_tsv = os.path.join(WORKFLOW_ROOT, PR_C.get("knowledge_dir", "knowledge/kegg"), "pathway_ko_reaction.tsv"),
        out_dir = lambda wildcards: os.path.join(PR_BASE_RES, wildcards.sample, "pathways"),
        scope_filter = PR_C.get("pathway_scope_filter", "prokaryote_relevant")
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
            -k "{input.ko_strict}" \
            -p "{params.pw_tsv}" \
            -t "{params.triples_tsv}" \
            -m "{input.completeness}" \
            -o "{params.out_dir}" \
            --scope_filter "{params.scope_filter}" 2>&1 | tee "{log}"
        """


# 4.4 Tier 3: Evaluate BRITE Functional Hierarchy & Molecular Systems
rule pathway_evaluate_brite_hierarchy:
    input:
        ko_strict = get_input_ko_strict
    output:
        counts = os.path.join(PR_BASE_RES, "{sample}", "brite", "genome_brite_counts.tsv"),
        matrix = os.path.join(PR_BASE_RES, "{sample}", "brite", "genome_brite_matrix.tsv")
    log:
        os.path.join(PR_BASE_LOGS, "{sample}", "evaluate_brite_hierarchy.log")
    threads: get_threads(PR_C.get("threads", 4))
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = pr_env("python"),
        script = os.path.join(PR_SCRIPTS, "evaluate_brite_hierarchy.py"),
        brite_tsv = os.path.join(WORKFLOW_ROOT, PR_C.get("knowledge_dir", "knowledge/kegg"), "ko_hierarchy_membership.tsv"),
        out_dir = lambda wildcards: os.path.join(PR_BASE_RES, wildcards.sample, "brite")
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
            -k "{input.ko_strict}" \
            -b "{params.brite_tsv}" \
            -o "{params.out_dir}" 2>&1 | tee "{log}"
        """


# 4.5 Global Coverage Audit (Separating Protein->KO from KO->Hierarchy)
rule pathway_audit_genome_coverage:
    input:
        ko_strict = get_input_ko_strict
    output:
        audit = os.path.join(PR_BASE_RES, "{sample}", "audit", "genome_ko_coverage_audit.tsv")
    log:
        os.path.join(PR_BASE_LOGS, "{sample}", "audit_genome_coverage.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = pr_env("python"),
        script = os.path.join(PR_SCRIPTS, "audit_genome_coverage.py"),
        modules_tsv = os.path.join(WORKFLOW_ROOT, PR_C.get("knowledge_dir", "knowledge/kegg"), "all_kegg_modules.tsv"),
        pw_tsv = os.path.join(WORKFLOW_ROOT, PR_C.get("knowledge_dir", "knowledge/kegg"), "all_kegg_pathways.tsv"),
        brite_tsv = os.path.join(WORKFLOW_ROOT, PR_C.get("knowledge_dir", "knowledge/kegg"), "ko_hierarchy_membership.tsv")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.audit}")" "$(dirname "{log}")"
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
            -k "{input.ko_strict}" \
            -m "{params.modules_tsv}" \
            -p "{params.pw_tsv}" \
            -b "{params.brite_tsv}" \
            -o "{output.audit}" 2>&1 | tee "{log}"
        """


# 4.6 Tier 4: Unified Schema v2.0 Functional Profile
rule pathway_generate_functional_profile:
    input:
        module_completeness = os.path.join(PR_BASE_RES, "{sample}", "modules", "genome_module_completeness.tsv"),
        footprint = os.path.join(PR_BASE_RES, "{sample}", "pathways", "genome_pathway_footprint.tsv"),
        brite_counts = os.path.join(PR_BASE_RES, "{sample}", "brite", "genome_brite_counts.tsv"),
        coverage_audit = os.path.join(PR_BASE_RES, "{sample}", "audit", "genome_ko_coverage_audit.tsv")
    output:
        profile = os.path.join(PR_BASE_RES, "{sample}", "capabilities", "genome_functional_profile.json")
    log:
        os.path.join(PR_BASE_LOGS, "{sample}", "generate_functional_profile.log")
    threads: get_threads(2)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads),
        slurm_account = lambda wildcards: get_slurm_account(),
        slurm_partition = lambda wildcards: get_slurm_partition("64c512g")
    params:
        env = pr_env("python"),
        script = os.path.join(PR_SCRIPTS, "generate_functional_profile.py"),
        manifest = os.path.join(WORKFLOW_ROOT, PR_C.get("knowledge_dir", "knowledge/kegg"), "kegg_snapshot_manifest.json")
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.profile}")" "$(dirname "{log}")"
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
            -c "{input.module_completeness}" \
            -f "{input.footprint}" \
            -b "{input.brite_counts}" \
            -a "{input.coverage_audit}" \
            -m "{params.manifest}" \
            -o "{output.profile}" 2>&1 | tee "{log}"
        """

