# ==============================================================================
# Advanced Analysis Module: Standalone Custom Tasks
# ==============================================================================
import os

# Load independent configuration (Managed by Snakefile registry)
# configfile: "config/advanced_analysis.yaml"
ADV_C = config.get("advanced_analysis", config)

# Results directory
ADV_RES = os.path.join(
    config.get("project_root", ""),
    ADV_C.get("results_dir", "results/advanced_analysis"),
)
ADV_LOGS = os.path.join(
    config.get("project_root", ""),
    ADV_C.get("logs_dir", "logs/advanced_analysis"),
)

def adv_env(name):
    env_name = config.get("environments", {}).get(name) or ADV_C.get("conda_envs", {}).get(name, name)
    env_dir = config.get("conda_env_dir", "")
    if not env_dir or os.path.isabs(env_name):
        return env_name
    return os.path.join(env_dir, env_name)

# Helper for cblaster parameter encoding
def get_cblaster_param_str(p):
    return "g{}_i{}_h{}_ig{}_d{}".format(
        p.get("gap_size", 20000),
        p.get("min_identity", 0.3),
        p.get("min_hits", 3),
        "T" if p.get("intermediate_genes", True) else "F",
        p.get("intermediate_max_distance", 5000)
    )

# --- Target and Path Mapping ---
ADV_TARGETS = []

# Setup paths and environment resolution for external genome BLASTp
ext_batch = config.get("external_genomes_batch", "")
EXT_BASE = config.get("annotation", {}).get("external_results_dir", "results/external_annotation")
if ext_batch:
    EXT_RES = os.path.join(config.get("project_root", "."), EXT_BASE, ext_batch)
else:
    EXT_RES = os.path.join(config.get("project_root", "."), EXT_BASE)
EXT_GENOMES = config.get("external_genomes", {})
if EXT_GENOMES is None: EXT_GENOMES = {}

BLASTP_C = ADV_C.get("blastp", {})
do_blastp = str(BLASTP_C.get("enable", "False")).lower() == "true"
if do_blastp:
    tasks_dict = BLASTP_C.get("tasks") or {}
    if not tasks_dict:
        default_query = config.get("target_protein") or BLASTP_C.get("blastp_query") or "config/target_genes.faa"
        tasks_dict = {
            "default": {
                "query": default_query,
                "blastp_tool": BLASTP_C.get("blastp_tool", "blastp"),
                "blastp_evalue": BLASTP_C.get("blastp_evalue", 1e-5),
                "blastp_identity": BLASTP_C.get("blastp_identity", 30.0),
                "blastp_coverage": BLASTP_C.get("blastp_coverage", 50.0)
            }
        }
    for task_name in tasks_dict.keys():
        for name in EXT_GENOMES.keys():
            ADV_TARGETS.append(os.path.join(EXT_RES, name, f"blastp/{task_name}/{name}_target_genes_annotated.tsv"))
            ADV_TARGETS.append(os.path.join(EXT_RES, name, f"blastp/{task_name}/{name}_matched_target_genes.faa"))
            ADV_TARGETS.append(os.path.join(EXT_RES, name, f"blastp/{task_name}/{name}_matched_target_genes.ffn"))
        
        if EXT_GENOMES:
            ADV_TARGETS.append(os.path.join(EXT_RES, f"blastp/{task_name}/global_target_genes_tidy_summary.tsv"))
            ADV_TARGETS.append(os.path.join(EXT_RES, f"blastp/{task_name}/global_target_genes_matrix_summary.tsv"))
            ADV_TARGETS.append(os.path.join(EXT_RES, f"blastp/{task_name}/global_operon_prediction_summary.tsv"))


SYNTENY_MAP = {} # name -> (path1, path2)

# Resolve ambiguity: local search and hmm search have same output patterns
ruleorder: advanced_cblaster_search_hmm > advanced_cblaster_search_local

if ADV_C.get("synteny_comparison", {}).get("enable"):
    pairs_dict = ADV_C.get("synteny_comparison", {}).get("pairs") or {}
    for name, paths in pairs_dict.items():
        path1 = paths.get("table1")
        path2 = paths.get("table2")
        if path1 and path2:
            out_file = os.path.join(ADV_RES, f"synteny/{name}_alignment.tsv")
            ADV_TARGETS.append(out_file)
            SYNTENY_MAP[name] = (path1, path2)

CLINKER_MAP = {} # seq_id -> absolute_gbk_path

if ADV_C.get("clinker_synteny", {}).get("enable"):
    clinker_pairs = ADV_C.get("clinker_synteny", {}).get("pairs") or {}
    for name, paths in clinker_pairs.items():
        gbk1 = paths.get("gbk1")
        gbk2 = paths.get("gbk2")
        if gbk1 and gbk2:
            ADV_TARGETS.append(os.path.join(ADV_RES, f"clinker/{name}.html"))
            CLINKER_MAP[f"{name}_gbk1"] = gbk1
            CLINKER_MAP[f"{name}_gbk2"] = gbk2

CBLASTER_TASKS = {}
CBLASTER_GBK_MAP = {}   # target_id -> gbk_path
CBLASTER_QUERY_MAP = {} # task_name -> final_query_path (possibly cropped)

if ADV_C.get("cblaster_search", {}).get("enable"):
    tasks = ADV_C.get("cblaster_search", {}).get("tasks") or {}
    for t_name, t_info in tasks.items():
        if t_info is None: t_info = {}
        # Skip example task
        if t_name == "BGC_Discovery_Task1" and t_info.get("query") == "/path/to/your/100kb_query.gbk":
            continue
            
        q_raw = t_info.get("query_gbk") or t_info.get("query")
        q_hmm = t_info.get("query_hmm")
        region = t_info.get("region")
        targets = t_info.get("targets") or []
        
        # 1. Resolve Query Path
        final_q = q_raw
        if q_raw and region:
            # If region is specified, we need to extract it first
            final_q = os.path.join(ADV_RES, f"cblaster/queries/{t_name}_region.gbk")
            # Register this query for the prep rule
            CBLASTER_QUERY_MAP[t_name] = {
                "raw": q_raw,
                "region": region,
                "final": final_q
            }
        
        # 2. Resolve Targets
        target_ids = []
        if isinstance(targets, dict):
            for t_id, info in targets.items():
                target_ids.append(t_id)
                CBLASTER_GBK_MAP[t_id] = info.get("path") if isinstance(info, dict) else info
        elif isinstance(targets, list):
            for gbk_path in targets:
                if gbk_path:
                    t_id = os.path.basename(gbk_path).rsplit(".", 1)[0]
                    target_ids.append(t_id)
                    CBLASTER_GBK_MAP[t_id] = gbk_path
        else:
            # fallback/single string target
            gbk_path = str(targets)
            if gbk_path and gbk_path != "[]":
                t_id = os.path.basename(gbk_path).rsplit(".", 1)[0]
                target_ids.append(t_id)
                CBLASTER_GBK_MAP[t_id] = gbk_path
            
        if (q_raw or q_hmm) and targets:
            mode = "hmm" if q_hmm else "search"
            
            hmm_names = ""
            if mode == "hmm" and q_hmm and os.path.isdir(q_hmm):
                h_files = [os.path.join(q_hmm, f) for f in os.listdir(q_hmm) if f.endswith(".hmm")]
                hmm_names = " ".join(h_files)

            CBLASTER_TASKS[t_name] = {
                "mode": mode,
                "query": final_q,
                "query_hmm": q_hmm,
                "hmm_names": hmm_names,
                "target_ids": target_ids,
                "params": t_info.get("params", {}),
                "param_str": get_cblaster_param_str(t_info.get("params", {}))
            }
            
            # Register targets
            for tid in target_ids:
                p_str = CBLASTER_TASKS[t_name]["param_str"]
                ADV_TARGETS.append(os.path.join(ADV_RES, f"cblaster/{t_name}/{tid}/{p_str}/{tid}_clinker.html"))

# --- Rules ---

rule advanced_analysis_all:
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
rule advanced_cblaster_prep_query:
    input:
        gbk = lambda w: CBLASTER_QUERY_MAP[w.task_name]["raw"]
    output:
        gbk = os.path.join(ADV_RES, "cblaster/queries/{task_name}_region.gbk")
    log:
        os.path.join(ADV_LOGS, "cblaster_prep_query/{task_name}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        region = lambda w: CBLASTER_QUERY_MAP[w.task_name]["region"]
    shell:
        r"""
        {{
        set -euo pipefail
        python3 scripts/extract_gbk_region.py {input.gbk} {output.gbk} {params.region}
        }} 2>&1 | tee {log}
        """
rule advanced_synteny_standalone:
    output:
        alignment = os.path.join(ADV_RES, "synteny/{name}_alignment.tsv")
    log:
        os.path.join(ADV_LOGS, "synteny/{name}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        paths = lambda wildcards: SYNTENY_MAP.get(wildcards.name, ("", "")),
        anchor = ADV_C.get("synteny_comparison", {}).get("anchor_ko", "K03040")
    shell:
        r"""
        {{
        set -euo pipefail
        python3 scripts/compare_synteny.py \
            {params.paths[0]} \
            {params.paths[1]} \
            {output.alignment} \
            {params.anchor}
        }} 2>&1 | tee {log}
        """
rule advanced_clinker_run:
    input:
        gbk1 = lambda wildcards: CLINKER_MAP.get(f"{wildcards.name}_gbk1", ""),
        gbk2 = lambda wildcards: CLINKER_MAP.get(f"{wildcards.name}_gbk2", "")
    output:
        html = os.path.join(ADV_RES, f"clinker/{{name}}.html")
    log:
        os.path.join(ADV_LOGS, "clinker_run/{name}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_env("clinker")
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} clinker {input.gbk1} {input.gbk2} -p {output.html}
        }} 2>&1 | tee {log}
        """
rule advanced_cblaster_makedb:
    input:
        gbk = lambda wildcards: CBLASTER_GBK_MAP.get(wildcards.target_id, "")
    output:
        sqlite = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.sqlite3"),
        dmnd = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.dmnd"),
        fasta = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.fasta")
    log:
        os.path.join(ADV_LOGS, "cblaster_makedb/{target_id}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_env("cblaster"),
        prefix = lambda w: os.path.join(ADV_RES, f"cblaster/db_prep/{w.target_id}.sqlite")
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} cblaster makedb {input.gbk} -n {params.prefix} --force
        }} 2>&1 | tee {log}
        """
rule advanced_cblaster_search_local:
    input:
        query = lambda w: [CBLASTER_TASKS[w.task_name]["query"]] if CBLASTER_TASKS[w.task_name]["mode"] == "search" else [],
        sqlite = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.sqlite3"),
        dmnd = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.dmnd")
    output:
        session = os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/session.json"),
        report = os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/search_report.txt")
    log:
        os.path.join(ADV_LOGS, "cblaster_search_local/{task_name}/{target_id}_{param_str}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_env("cblaster"),
        gap = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("gap_size", 20000),
        ident = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("min_identity", 0.3),
        hits = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("min_hits", 3),
        ig = lambda w: "-ig" if CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("intermediate_genes", True) else "",
        ig_dist = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("intermediate_max_distance", 5000)
    wildcard_constraints:
        task_name = "|".join([name for name, info in CBLASTER_TASKS.items() if info["mode"] == "search"]) or "MATCH_NOTHING"
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} cblaster search -m local -qf {input.query} -db {input.dmnd} \
            --gap {params.gap} --min_identity {params.ident} --min_hits {params.hits} \
            {params.ig} --max_distance {params.ig_dist} \
            -o {output.report} -s {output.session}
        }} 2>&1 | tee {log}
        """
rule advanced_cblaster_search_hmm:
    input:
        sqlite = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.sqlite3"),
        dmnd = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.dmnd"),
        fasta = os.path.join(ADV_RES, "cblaster/db_prep/{target_id}.sqlite.fasta")
    output:
        session = os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/session.json"),
        report = os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/search_report.txt")
    log:
        os.path.join(ADV_LOGS, "cblaster_search_hmm/{task_name}/{target_id}_{param_str}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_env("cblaster"),
        hmm_profiles = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("hmm_names", ""),
        gap = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("gap_size", 20000),
        ident = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("min_identity", 0.3),
        hits = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("min_hits", 3),
        ig = lambda w: "-ig" if CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("intermediate_genes", True) else "",
        ig_dist = lambda w: CBLASTER_TASKS.get(w.task_name, {}).get("params", {}).get("intermediate_max_distance", 5000)
    wildcard_constraints:
        task_name = "|".join([name for name, info in CBLASTER_TASKS.items() if info["mode"] == "hmm"]) or "MATCH_NOTHING"
    shell:
        r"""
        {{
        set -euo pipefail
        # -qp takes full paths to local .hmm files (cblaster checks Path.exists() to decide local vs Pfam)
        # -db takes the fasta produced by makedb (cblaster auto-finds the matching .sqlite3 by stripping extension)
        conda run -p {params.env} cblaster search -m hmm -qp {params.hmm_profiles} -db {input.fasta} \
            --gap {params.gap} --min_identity {params.ident} --min_hits {params.hits} \
            {params.ig} --max_distance {params.ig_dist} \
            -o {output.report} -s {output.session}
        }} 2>&1 | tee {log}
        """
rule advanced_cblaster_extract:
    input:
        session = os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/session.json")
    output:
        cluster_dir = directory(os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/extracted_clusters"))
    log:
        os.path.join(ADV_LOGS, "cblaster_extract/{task_name}/{target_id}_{param_str}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_env("cblaster")
    shell:
        r"""
        {{
        set -euo pipefail
        rm -rf {output.cluster_dir}
        mkdir -p {output.cluster_dir}
        conda run -p {params.env} cblaster extract_clusters {input.session} -o {output.cluster_dir} || true
        }} 2>&1 | tee {log}
        """
rule advanced_cblaster_clinker:
    input:
        query = lambda w: [CBLASTER_TASKS[w.task_name]["query"]] if CBLASTER_TASKS[w.task_name]["mode"] == "search" else [],
        cluster_dir = os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/extracted_clusters")
    output:
        html = os.path.join(ADV_RES, "cblaster/{task_name}/{target_id}/{param_str}/{target_id}_clinker.html")
    log:
        os.path.join(ADV_LOGS, "cblaster_clinker/{task_name}/{target_id}_{param_str}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = adv_env("clinker")
    shell:
        r"""
        {{
        set -euo pipefail
        count=$(ls -1 {input.cluster_dir}/*.gbk 2>/dev/null | wc -l || echo 0)
        if [ "$count" -eq 0 ]; then
            echo "<html><body><h1>No homologous clusters were found by cblaster in this genome under current parameters.</h1></body></html>" > {output.html}
            exit 0
        fi
        conda run -p {params.env} clinker {input.query} {input.cluster_dir}/*.gbk -p {output.html}
        }} 2>&1 | tee {log}
        """

include: "blastp.smk"
