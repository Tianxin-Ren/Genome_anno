# ==============================================================================
# Visualize Module: Pipeline Visualization & Interactive HTML Dashboard Reporting
# ==============================================================================
import os

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
VI_ROOT = config.get("project_root", "")
VI_ENV_DIR = config.get("conda_env_dir", "")
VI_C = config.get("visualize", {})

def vi_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    # Prefix relative paths with the global project_root
    return os.path.normpath(os.path.join(VI_ROOT, rel_path))

def vi_get_env(name):
    env_name = VI_C.get("conda_envs", {}).get(name, name)
    if not VI_ENV_DIR: return env_name
    return os.path.join(VI_ENV_DIR, env_name)

# Results directories
VI_RES = vi_full_path(VI_C.get("results_dir", "results/visualize"))
VI_LOGS = vi_full_path(VI_C.get("logs_dir", os.path.join("logs", "visualize")))

def serialize_selected_genes(w):
    genes_config = (config.get("samples", {}).get(w.sample, {}) or {}).get("selected_genes", {})
    if genes_config is None:
        return ""
    serialized = []
    if isinstance(genes_config, dict):
        for k, v in genes_config.items():
            serialized.append(f"{k}:{v or ''}")
    elif isinstance(genes_config, list):
        for item in genes_config:
            if isinstance(item, dict):
                for k, v in item.items():
                    serialized.append(f"{k}:{v or ''}")
            else:
                serialized.append(f"{item or ''}")
    return ",".join(serialized)

# 2. Target Registration
# ------------------------------------------------------------------------------
VI_TARGETS = []

VI_SAMPLES = config.get("samples", {})
if VI_SAMPLES is None: VI_SAMPLES = {}

if VI_C.get("enable", False):
    # 1. Circular Genome Targets (Dynamic based on samples.yaml GFF configurations)
    for sample, s_info in VI_SAMPLES.items():
        if "gff" in s_info:
            VI_TARGETS.append(os.path.join(VI_RES, f"{sample}/{sample}_circular_map.html"))
            VI_TARGETS.append(os.path.join(VI_RES, f"{sample}/{sample}_circular_map.pdf"))
            VI_TARGETS.append(os.path.join(VI_RES, f"{sample}/{sample}_selected_genes.pdf"))

    # 2. Abundance Heatmap Target
    if VI_C.get("plots", {}).get("abundance_heatmap", {}).get("enable", False):
        VI_TARGETS.append(os.path.join(VI_RES, "abundance/abundance_heatmap.pdf"))

    # 3. Quality Assessment Target
    if VI_C.get("plots", {}).get("quality_assessment", {}).get("enable", False):
        VI_TARGETS.append(os.path.join(VI_RES, "quality/quality_scatter_plot.png"))

    # 4. Taxonomic Profile Target
    if VI_C.get("plots", {}).get("taxonomic_profile", {}).get("enable", False):
        VI_TARGETS.append(os.path.join(VI_RES, "taxonomy/taxonomic_distribution.html"))

    # 5. Global Summary Report (always built if enabled)
    VI_TARGETS.append(os.path.join(VI_RES, "summary_report.html"))

TARGETS.extend(VI_TARGETS)


# 3. Rules
# ------------------------------------------------------------------------------

rule visualize_abundance_heatmap:
    output:
        pdf = os.path.join(VI_RES, "abundance/abundance_heatmap.pdf")
    log:
        os.path.join(VI_LOGS, "abundance_heatmap.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = vi_get_env("python"),
        palette = VI_C.get("plots", {}).get("abundance_heatmap", {}).get("palette", "viridis"),
        input_type = VI_C.get("plots", {}).get("abundance_heatmap", {}).get("input_type", "coverm"),
        abundance_dir = vi_full_path(config.get("abundance", {}).get("results_dir", "results/abundance"))
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.pdf}")"
        conda run -p {params.env} python3 scripts/generate_vi_report.py \
            --mode heatmap \
            --palette {params.palette} \
            --input-type {params.input_type} \
            --abundance-dir {params.abundance_dir} \
            --output {output.pdf} 2>&1 | tee {log}
        """


rule visualize_quality_assessment:
    output:
        png = os.path.join(VI_RES, "quality/quality_scatter_plot.png")
    log:
        os.path.join(VI_LOGS, "quality_assessment.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = vi_get_env("python"),
        input_type = VI_C.get("plots", {}).get("quality_assessment", {}).get("input_type", "checkm"),
        binning_dir = vi_full_path(config.get("binning", {}).get("results_dir", "results/binning"))
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.png}")"
        conda run -p {params.env} python3 scripts/generate_vi_report.py \
            --mode quality \
            --input-type {params.input_type} \
            --binning-dir {params.binning_dir} \
            --output {output.png} 2>&1 | tee {log}
        """


rule visualize_taxonomic_profile:
    output:
        html = os.path.join(VI_RES, "taxonomy/taxonomic_distribution.html")
    log:
        os.path.join(VI_LOGS, "taxonomic_profile.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = vi_get_env("python"),
        input_type = VI_C.get("plots", {}).get("taxonomic_profile", {}).get("input_type", "gtdbtk"),
        annotation_dir = vi_full_path(config.get("annotation", {}).get("results_dir", "results/annotation"))
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.html}")"
        conda run -p {params.env} python3 scripts/generate_vi_report.py \
            --mode taxonomy \
            --input-type {params.input_type} \
            --annotation-dir {params.annotation_dir} \
            --output {output.html} 2>&1 | tee {log}
        """


def get_circular_inputs(wildcards):
    s_info = config.get("samples", {}).get(wildcards.sample, {})
    inputs = {}
    gff_path = s_info.get("gff", "")
    if gff_path:
        inputs["gff"] = vi_full_path(gff_path)
    tpm_path = s_info.get("tpm", "")
    if tpm_path:
        inputs["tpm"] = vi_full_path(tpm_path)
    return inputs

rule visualize_circular_genome:
    input:
        unpack(get_circular_inputs)
    output:
        html = os.path.join(VI_RES, "{sample}/{sample}_circular_map.html"),
        pdf = os.path.join(VI_RES, "{sample}/{sample}_circular_map.pdf"),
        bar_pdf = os.path.join(VI_RES, "{sample}/{sample}_selected_genes.pdf")
    log:
        os.path.join(VI_LOGS, "circular_map_{sample}.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = vi_get_env("python"),
        gff_path = lambda w, input: getattr(input, "gff", ""),
        tpm_path = lambda w, input: getattr(input, "tpm", ""),
        abundance_metric = VI_C.get("circular_genome", {}).get("abundance_metric", "auto"),
        theme = VI_C.get("circular_genome", {}).get("theme", VI_C.get("report", {}).get("theme", "dark")),
        selected_genes = serialize_selected_genes
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.html}")"
        
        conda run -p {params.env} python3 scripts/generate_vi_report.py \
            --mode circular_genome \
            --gff "{params.gff_path}" \
            --tpm "{params.tpm_path}" \
            --theme "{params.theme}" \
            --sample "{wildcards.sample}" \
            --selected-genes "{params.selected_genes}" \
            --title "{wildcards.sample} Circular Map" \
            --abundance-metric "{params.abundance_metric}" \
            --output {output.html} \
            --output-pdf {output.pdf} \
            --output-bar-pdf {output.bar_pdf} \
            --interactive "True" 2>&1 | tee {log}
        """


def get_report_inputs(wildcards):
    inputs = []
    # 1. Abundance Heatmap
    if VI_C.get("plots", {}).get("abundance_heatmap", {}).get("enable", False):
        inputs.append(os.path.join(VI_RES, "abundance/abundance_heatmap.pdf"))
    # 2. Quality Assessment
    if VI_C.get("plots", {}).get("quality_assessment", {}).get("enable", False):
        inputs.append(os.path.join(VI_RES, "quality/quality_scatter_plot.png"))
    # 3. Taxonomic Profile
    if VI_C.get("plots", {}).get("taxonomic_profile", {}).get("enable", False):
        inputs.append(os.path.join(VI_RES, "taxonomy/taxonomic_distribution.html"))
    # 4. Circular Genome Maps
    for s, info in config.get("samples", {}).items():
        if "gff" in info:
            inputs.append(os.path.join(VI_RES, f"{s}/{s}_circular_map.html"))
    return inputs

rule visualize_generate_report:
    input:
        plots = get_report_inputs
    output:
        report = os.path.join(VI_RES, "summary_report.html")
    log:
        os.path.join(VI_LOGS, "summary_report.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = vi_get_env("python"),
        title = VI_C.get("report", {}).get("title", "Genome_anno Analysis Report"),
        theme = VI_C.get("report", {}).get("theme", "dark"),
        interactive = str(VI_C.get("report", {}).get("interactive", True)),
        abundance_dir = vi_full_path(config.get("abundance", {}).get("results_dir", "results/abundance")),
        binning_dir = vi_full_path(config.get("binning", {}).get("results_dir", "results/binning")),
        annotation_dir = vi_full_path(config.get("annotation", {}).get("results_dir", "results/annotation"))
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname "{output.report}")"
        conda run -p {params.env} python3 scripts/generate_vi_report.py \
            --mode report \
            --title "{params.title}" \
            --theme {params.theme} \
            --interactive {params.interactive} \
            --abundance-dir {params.abundance_dir} \
            --binning-dir {params.binning_dir} \
            --annotation-dir {params.annotation_dir} \
            --output {output.report} 2>&1 | tee {log}
        """
