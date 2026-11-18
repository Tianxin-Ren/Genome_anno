# Functional Annotation Rules
# ------------------------------------------------------------------------------

# ==============================================================================
# Helper for Dynamic Unified HMM Database (File, List, or Directory Concatenation)
# ==============================================================================
def resolve_hmm_inputs(db_param):
    if not db_param:
        return []
    
    # Normalize to list
    if isinstance(db_param, str):
        if "," in db_param:
            raw_items = [s.strip() for s in db_param.split(",")]
        else:
            raw_items = [db_param]
    elif isinstance(db_param, list):
        raw_items = db_param
    else:
        raw_items = []
        
    resolved_files = []
    for item in raw_items:
        path = an_full_path(item)
        if not path:
            continue
        if os.path.isdir(path):
            resolved_files.extend(sorted(glob.glob(os.path.join(path, "*.hmm"))))
        elif os.path.isfile(path):
            resolved_files.append(path)
        else:
            if path.endswith(".hmm"):
                resolved_files.append(path)
    return resolved_files

def get_hmm_files_for_genome(genome):
    sample_info = AN_SAMPLES.get(genome, {})
    if not isinstance(sample_info, dict):
        sample_info = {}
    ext_info = EXT_GENOMES.get(genome, {})
    if not isinstance(ext_info, dict):
        ext_info = {}
        
    hmm_param = sample_info.get("hmm_db") or ext_info.get("hmm_db") or AN_C.get("unified_hmm_db", "")
    return resolve_hmm_inputs(hmm_param)

def get_hmm_tool_for_genome(genome):
    sample_info = AN_SAMPLES.get(genome, {})
    if not isinstance(sample_info, dict):
        sample_info = {}
    ext_info = EXT_GENOMES.get(genome, {})
    if not isinstance(ext_info, dict):
        ext_info = {}
    return sample_info.get("hmm_tool") or ext_info.get("hmm_tool") or AN_C.get("hmm_tool", "hmmsearch")

def get_target_profiles_for_genome(genome):
    sample_info = AN_SAMPLES.get(genome, {})
    if not isinstance(sample_info, dict):
        sample_info = {}
    ext_info = EXT_GENOMES.get(genome, {})
    if not isinstance(ext_info, dict):
        ext_info = {}
    val = sample_info.get("target_profiles") or ext_info.get("target_profiles") or AN_C.get("target_profiles", "")
    if isinstance(val, list):
        return ",".join(val)
    return str(val) if val is not None else ""

def is_hmm_split(genome):
    info = AN_SAMPLES.get(genome)
    if isinstance(info, dict) and "split_hmm" in info:
        return str(info.get("split_hmm")).lower() == "true"
    return str(AN_C.get("split_hmm", "True")).lower() == "true"

def get_hmm_gather_inputs(wildcards):
    genome = wildcards.genome
    if is_hmm_split(genome):
        return {
            "unified": expand(os.path.join(AN_RES, f"{genome}/advanced_hmm/unified/chunk_{{sid}}.tbl"), 
                           sid=[f"{i+1}-of-{AN_C.get('hmm_split_chunks', 20)}" for i in range(AN_C.get('hmm_split_chunks', 20))]),
            "pfam": get_pfam_chunks(wildcards),
            "faa": [],
            "db": [],
            "pfam_db": []
        }
    else:
        steps = config.get("annotation", {}).get("steps", {})
        do_pfam = str(steps.get("pfam", "True")).lower() == "true"
        inputs = {
            "unified": [],
            "pfam": [],
            "faa": get_anno_faa_input(wildcards),
            "db": get_unified_hmm_db(wildcards)
        }
        if do_pfam:
            inputs["pfam_db"] = get_pfam_db_input(wildcards)
        else:
            inputs["pfam_db"] = []
        return inputs

def get_unified_hmm_db(wildcards):
    genome = wildcards.genome if hasattr(wildcards, "genome") else wildcards.bin
    resolved = get_hmm_files_for_genome(genome)
    if not resolved:
        sample_info = AN_SAMPLES.get(genome, {})
        ext_info = EXT_GENOMES.get(genome, {})
        hmm_param = (sample_info.get("hmm_db") if isinstance(sample_info, dict) else None) or \
                    (ext_info.get("hmm_db") if isinstance(ext_info, dict) else None) or \
                    AN_C.get("unified_hmm_db", "")
        raise ValueError(
            f"ERROR: No HMM database files could be resolved for sample '{genome}' "
            f"from HMM parameter '{hmm_param}'. Please check that the file/folder path is correct "
            "and accessible on this cluster account."
        )
    tool = get_hmm_tool_for_genome(genome)
    if len(resolved) > 1 or tool == "hmmscan":
        res_dir = EXT_RES if genome in EXT_GENOMES else AN_RES
        return os.path.join(res_dir, genome, "advanced_hmm/combined_unified.hmm")
    else:
        return resolved[0]

def get_cat_hmm_inputs(wildcards):
    genome = wildcards.genome if hasattr(wildcards, "genome") else wildcards.bin
    return get_hmm_files_for_genome(genome)

def get_pfam_db_input(wildcards):
    db_path = AN_C.get("pfam_hmm_db", "")
    return [db_path] if db_path else []

rule annotation_cat_unified_hmm:
    input:
        hmms = get_cat_hmm_inputs
    output:
        db = os.path.join(AN_RES, "{genome}/advanced_hmm/combined_unified.hmm")
    log:
        os.path.join(AN_LOGS, "cat_unified_hmm/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        tool = lambda wildcards: get_hmm_tool_for_genome(wildcards.genome)
    shell:
        """
        set -euo pipefail
        {{
        cat {input.hmms} > {output.db}
        if [ "{params.tool}" = "hmmscan" ]; then
            conda run -p {params.env} hmmpress -f {output.db}
        fi
        }} 2>&1 | tee {log}
        """
rule annotation_eggnog:
    input:
        faa = get_anno_faa_input
    output:
        anno = os.path.join(AN_RES, "{genome}/eggnog/{genome}.emapper.annotations")
    log:
        os.path.join(AN_LOGS, "eggnog/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("eggnog"),
        db = AN_C.get("eggnog_db", ""),
        prefix = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "eggnog", wildcards.genome)
    shell:
        r"""
        set -euo pipefail
        conda run -p {params.env} emapper.py -m mmseqs --dbmem -i {input.faa} --output {params.prefix} --cpu {threads} --data_dir {params.db} --override 2>&1 | tee {log}
        
        """
rule annotation_kofamscan:
    input:
        faa = get_anno_faa_input
    output:
        tsv = os.path.join(AN_RES, "{genome}/kofamscan/{genome}_kofam.tsv")
    log:
        os.path.join(AN_LOGS, "kofamscan/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("kofamscan"),
        db = AN_C.get("kofam_db", ""),
        tmpdir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "kofamscan/tmp")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.tmpdir}
        conda run -p {params.env} exec_annotation \
            -o {output.tsv} \
            --cpu {threads} \
            --tmp {params.tmpdir} \
            -p {params.db}/profiles \
            -k {params.db}/ko_list \
            {input.faa} 2>&1 | tee {log}
        rm -rf {params.tmpdir}
        
        """
def get_aggregate_inputs(wildcards):
    """Dynamically determine inputs for aggregation based on config."""
    import sys
    anno_cfg = config.get("annotation", {})
    steps = anno_cfg.get("steps", {})
    
    # 打印到屏幕，这样你在 snakemake -n 时能看到
    print(f"--- [DEBUG] Genome: {wildcards.genome} ---", file=sys.stderr)
    print(f"--- [DEBUG] Steps dict: {steps} ---", file=sys.stderr)
    
    do_eggnog = str(steps.get("eggnog", "False")).lower() == "true"
    do_kofam = str(steps.get("kofamscan", "False")).lower() == "true"
    do_adv_hmm = str(steps.get("advanced_hmm", "False")).lower() == "true"

    print(f"--- [DEBUG] do_kofam resolved to: {do_kofam} ---", file=sys.stderr)

    inputs = {
        "faa": os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.faa"),
        "gff": os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.gff3"),
        "script": "scripts/aggregate_functional.py"
    }
    
    if do_eggnog:
        inputs["eggnog"] = os.path.join(AN_RES, f"{wildcards.genome}/eggnog/{wildcards.genome}.emapper.annotations")
        
    if do_kofam:
        inputs["kofam"] = os.path.join(AN_RES, f"{wildcards.genome}/kofamscan/{wildcards.genome}_kofam.tsv")
        
    if do_adv_hmm:
        inputs["advanced_hmm"] = os.path.join(AN_RES, f"{wildcards.genome}/{wildcards.genome}_advanced_hmm_resolved.tsv")
    
    return inputs



# ==========================================
# Advanced HMM Arbitration (Scatter/Gather)
# ==========================================

# 必须对 merge 后且重命名过的 FAA 进行切分，以包含环状 Contig 和正确 ID
rule annotation_hmm_split:
    input:
        faa = get_anno_faa_input
    output:
        # temp(expand(os.path.join(AN_RES, "{{genome}}/advanced_hmm/chunks/chunk_{sid}.faa"), 
        #      sid=[f"{i+1}-of-{AN_C.get('hmm_split_chunks', 20)}" for i in range(AN_C.get('hmm_split_chunks', 20))]))
        expand(os.path.join(AN_RES, "{{genome}}/advanced_hmm/chunks/chunk_{sid}.faa"), 
             sid=[f"{i+1}-of-{AN_C.get('hmm_split_chunks', 20)}" for i in range(AN_C.get('hmm_split_chunks', 20))])
    log:
        os.path.join(AN_LOGS, "hmm_split/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        outdir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "advanced_hmm/chunks"),
        split_count = AN_C.get("hmm_split_chunks", 20),
        env = an_get_env("bakta")
    shell:
        """
            set -euo pipefail
        {{
        mkdir -p {params.outdir}/tmp
        # Clean gap characters '-' and stop codons '*' from sequence lines in FAA
        CLEAN_FAA="{params.outdir}/cleaned_input.faa"
        sed '/^[^>]/s/[-*]//g' "{input.faa}" > "$CLEAN_FAA"
        conda run -p {params.env} seqkit split "$CLEAN_FAA" -p {params.split_count} -O {params.outdir}/tmp
        rm -f "$CLEAN_FAA"
        # 重命名以匹配 sid (i-of-N)
        i=1
        for file in $(ls {params.outdir}/tmp/*.faa | sort -V); do
            mv "$file" "{params.outdir}/chunk_${{i}}-of-{params.split_count}.faa"
            ((i++))
        done
        rm -rf {params.outdir}/tmp
        }} 2>&1 | tee {log}
        
        """
rule annotation_hmmsearch_unified:
    input:
        faa = os.path.join(AN_RES, "{genome}/advanced_hmm/chunks/chunk_{scatteritem}.faa"),
        db = get_unified_hmm_db
    output:
        tbl = os.path.join(AN_RES, "{genome}/advanced_hmm/unified/chunk_{scatteritem}.tbl")
    log:
        os.path.join(AN_LOGS, "hmm_unified/{genome}_{scatteritem}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        tool = lambda wildcards: get_hmm_tool_for_genome(wildcards.genome),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5)
    shell:
        "conda run -p {params.env} {params.tool} --cpu {threads} --tblout {output.tbl} -E {params.evalue} {input.db} {input.faa} > {log} 2>&1"
rule annotation_hmmsearch_pfam:
    input:
        faa = os.path.join(AN_RES, "{genome}/advanced_hmm/chunks/chunk_{scatteritem}.faa"),
        db = get_pfam_db_input
    output:
        tbl = os.path.join(AN_RES, "{genome}/advanced_hmm/pfam/chunk_{scatteritem}.domtbl")
    log:
        os.path.join(AN_LOGS, "hmm_pfam/{genome}_{scatteritem}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5)
    shell:
        "conda run -p {params.env} hmmsearch --cpu {threads} --domtblout {output.tbl} -E {params.evalue} {input.db} {input.faa} > {log} 2>&1"
def get_pfam_chunks(wildcards):
    """Dynamically return Pfam chunk files if steps.pfam is True, else empty list."""
    steps = config.get("annotation", {}).get("steps", {})
    do_pfam = str(steps.get("pfam", "True")).lower() == "true"
    genome = wildcards.genome if hasattr(wildcards, "genome") else wildcards.bin
    if do_pfam:
        return expand(os.path.join(AN_RES, f"{genome}/advanced_hmm/pfam/chunk_{{sid}}.domtbl"), 
                      sid=[f"{i+1}-of-{AN_C.get('hmm_split_chunks', 20)}" for i in range(AN_C.get('hmm_split_chunks', 20))])
    else:
        return []

rule annotation_hmm_gather:
    input:
        unpack(get_hmm_gather_inputs)
    output:
        unified_raw = os.path.join(AN_RES, "{genome}/advanced_hmm/unified_all_raw.tbl"),
        pfam_raw = os.path.join(AN_RES, "{genome}/advanced_hmm/pfam_all_raw.domtbl")
    log:
        os.path.join(AN_LOGS, "hmm_gather/{genome}.log")
    threads: lambda wildcards: 1 if is_hmm_split(wildcards.genome) else get_threads(16)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        is_split = lambda wildcards: "true" if is_hmm_split(wildcards.genome) else "false",
        tool = lambda wildcards: get_hmm_tool_for_genome(wildcards.genome if hasattr(wildcards, "genome") else wildcards.bin),
        env = an_get_env("bakta"),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5),
        run_pfam = lambda wildcards: str(config.get("annotation", {}).get("steps", {}).get("pfam", "True")).lower() == "true"
    shell:
        """
        set -euo pipefail
        {{
        if [ "{params.is_split}" = "true" ]; then
            echo "# {params.tool}" > {output.unified_raw}
            cat {input.unified} | grep -v '^#' >> {output.unified_raw} || true
            if [ -n "{input.pfam}" ]; then
                cat {input.pfam} | grep -v '^#' > {output.pfam_raw} || true
            else
                echo "# Skipping Pfam because steps.pfam is set to False" > {output.pfam_raw}
            fi
        else
            conda run -p {params.env} hmmsearch --cpu {threads} --tblout {output.unified_raw} -E {params.evalue} {input.db} {input.faa} > {log} 2>&1
            if [ "{params.run_pfam}" = "true" ]; then
                conda run -p {params.env} hmmsearch --cpu {threads} --domtblout {output.pfam_raw} -E {params.evalue} {input.pfam_db} {input.faa} >> {log} 2>&1
            else
                echo "# Skipping Pfam because steps.pfam is set to False" > {output.pfam_raw}
            fi
        fi
        }} 2>&1 | tee {log}
        """
rule annotation_hmm_resolve:
    input:
        unified = os.path.join(AN_RES, "{genome}/advanced_hmm/unified_all_raw.tbl"),
        pfam = os.path.join(AN_RES, "{genome}/advanced_hmm/pfam_all_raw.domtbl"),
        ko_list = AN_C.get("ko_list", ""),
        faa = get_anno_faa_input
    output:
        final_table = os.path.join(AN_RES, "{genome}/{genome}_advanced_hmm_resolved.tsv")
    log:
        os.path.join(AN_LOGS, "hmm_resolve/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("biopython"),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5),
        log_factor = AN_C.get("hmm_log_factor", 0.2),
        negative_map = AN_C.get("negative_map", ""),
        competition_fraction = AN_C.get("competition_fraction", 0.10),
        target_profiles = lambda wildcards: get_target_profiles_for_genome(wildcards.genome if hasattr(wildcards, "genome") else wildcards.bin)
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 scripts/resolve_annotations.py \
            --unified {input.unified} \
            --pfam {input.pfam} \
            --ko_list {input.ko_list} \
            --faa {input.faa} \
            --output {output.final_table} \
            --evalue_cutoff {params.evalue} \
            --log_factor {params.log_factor} \
            --negative_map {params.negative_map} \
            --competition_fraction {params.competition_fraction} \
            --target_profiles "{params.target_profiles}" 2>&1 | tee {log}
        
        """
rule annotation_aggregate:
    input:
        unpack(get_aggregate_inputs)
    output:
        summary = os.path.join(AN_RES, "{genome}/{genome}_functional_summary.tsv")
    log:
        os.path.join(AN_LOGS, "aggregate/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("biopython"),
        egg_in = lambda w, input: input.eggnog if hasattr(input, "eggnog") else "none",
        ko_in = lambda w, input: input.kofam if hasattr(input, "kofam") else "none",
        adv_hmm = lambda w, input: input.advanced_hmm if hasattr(input, "advanced_hmm") else "none"
    shell:
        r"""
        set -euo pipefail
        conda run -p {params.env} python3 {input.script} \
            {input.faa} {params.egg_in} {params.ko_in} {params.adv_hmm} \
            {input.gff} {wildcards.genome} {wildcards.genome} {output.summary} 2>&1 | tee {log}
        
        """
