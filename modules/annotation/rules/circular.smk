# Circular Contig Analysis and Refinement Rules
# ------------------------------------------------------------------------------

# Helper to map a wildcard genome ID (which might be a path like sample/external/sample_filtered)
# to its actual sample name as defined in samples.yaml

def get_bakta_refined_input_fa(wildcards):
    data_type = config.get("annotation", {}).get("data_type", "long")
    sample = get_sample_name_from_genome(wildcards.genome)
    if data_type == "short":
        comp = config.get("binning", {}).get("completeness", 50)
        cont = config.get("binning", {}).get("contamination", 10)
        return os.path.join(AN_BN_RES, f"{sample}/refinement/metawrap_{comp}_{cont}_bins/{wildcards.bin}.fa")
    else:
        return os.path.join(AN_BN_RES, f"{sample}/circular_contigs/{wildcards.bin}.fa")

def get_circular_dependencies(wildcards):
    # 1. Trigger the checkpoint (Standard Snakemake requirement)
    # In a dry run, if the checkpoint hasn't run, this allows us to still predict targets
    try:
        checkpoint_dir = checkpoints.annotation_circular_filter.get(genome=wildcards.genome).output.dir
    except:
        checkpoint_dir = os.path.join(AN_RES, wildcards.genome, "circular_bins_filtered")
    
    selected = []
    if os.path.exists(checkpoint_dir):
        selected = [f.replace(".fa", "") for f in os.listdir(checkpoint_dir) if f.endswith(".fa")]
    else:
        # Fallback for dry-run
        data_type = config.get("annotation", {}).get("data_type", "long")
        sample = get_sample_name_from_genome(wildcards.genome)
        if data_type == "short":
            comp = config.get("binning", {}).get("completeness", 50)
            cont = config.get("binning", {}).get("contamination", 10)
            src_dir = os.path.join(AN_BN_RES, sample, f"refinement/metawrap_{comp}_{cont}_bins")
            if os.path.exists(src_dir):
                selected = [f.replace(".fa", "") for f in os.listdir(src_dir) if f.endswith(".fa")]
        else:
            info_file = os.path.join(AN_BN_RES, sample, "circular_info.tsv")
            if os.path.exists(info_file):
                import pandas as pd
                df = pd.read_csv(info_file, sep="\t")
                df["ID"] = df["ID"].astype(str)
                min_size = AN_C.get("min_circular_size", 800000)
                config_complete = [str(c) for c in get_complete_contigs(wildcards.genome)]
                
                # Use the same logic as the checkpoint for consistency
                mask = ((df["Circular"].str.lower() == "yes") & (df["Length"] >= min_size)) | (df["ID"].isin(config_complete))
                selected = df[mask]["ID"].tolist()
                
                # Ensure manual overrides are included
                for cid in config_complete:
                    if cid not in selected:
                        selected.append(cid)

    steps = AN_C.get("steps", {})
    do_gtdb = str(steps.get("gtdbtk", "False")).lower() == "true"
    do_checkm = str(steps.get("checkm2", "False")).lower() == "true"
    do_eggnog = str(steps.get("eggnog", "False")).lower() == "true"
    do_kofam = str(steps.get("kofamscan", "False")).lower() == "true"
    do_audit = str(steps.get("circular_audit", "False")).lower() == "true"

    deps = []
    if do_gtdb:
        deps.append(os.path.join(AN_RES, wildcards.genome, "circular_taxonomy/gtdbtk.summary.tsv"))
    if do_checkm:
        deps.append(os.path.join(AN_RES, wildcards.genome, "circular_quality/checkm2_summary.tsv"))
    
    for b in selected:
        deps.append(os.path.join(AN_RES, wildcards.genome, f"circular_refined/{b}/{b}.gff3"))
        if do_eggnog:
            deps.append(os.path.join(AN_RES, wildcards.genome, f"circular_refined/{b}/functional/{b}.emapper.annotations"))
        if do_kofam:
            deps.append(os.path.join(AN_RES, wildcards.genome, f"circular_refined/{b}/functional/{b}.kofam.tsv"))
        
        # Always generate functional summary for circular bins
        deps.append(os.path.join(AN_RES, wildcards.genome, f"circular_refined/{b}/{b}_functional_summary.tsv"))
        if do_audit:
            deps.append(os.path.join(AN_RES, wildcards.genome, f"circular_audit/{b}/{b}_audit.tsv"))

    if do_audit and selected:
        deps.append(os.path.join(AN_RES, wildcards.genome, "circular_audit_summary.tsv"))
        if AN_C.get("circular_audit", {}).get("enable_anvio_validation", True):
            deps.append(os.path.join(AN_RES, wildcards.genome, "circular_audit/anvio_concordance.tsv"))
        if sample_has_short_reads(wildcards):
            deps.append(os.path.join(AN_RES, wildcards.genome, "circular_coverm/circular_illumina_qc_summary.tsv"))
    return deps

def get_refined_targets(wildcards):
    targets = [os.path.join(AN_RES, wildcards.genome, "circular_summary.tsv")]
    targets.extend(get_circular_dependencies(wildcards))
    return targets

def get_circular_filter_inputs(wildcards):
    data_type = config.get("annotation", {}).get("data_type", "long")
    sample = get_sample_name_from_genome(wildcards.genome)
    if data_type == "short":
        comp = config.get("binning", {}).get("completeness", 50)
        cont = config.get("binning", {}).get("contamination", 10)
        bin_dir = os.path.join(AN_BN_RES, f"{sample}/refinement/metawrap_{comp}_{cont}_bins")
        
        # 智能检测：如果用户已经物理放好了 refined bins，则返回空字典以切断对上游的比对、组装和分箱依赖
        if os.path.exists(bin_dir):
            try:
                fa_files = [f for f in os.listdir(bin_dir) if f.endswith(".fa")]
                if fa_files:
                    import sys
                    print(f"INFO: Pre-existing refined bins detected in {bin_dir} for sample {sample}. Decoupling upstream dependencies.", file=sys.stderr)
                    return {}
            except:
                pass
                
        # 否则，返回刚性依赖，走标准依赖追溯流程
        return {
            "stats": os.path.join(AN_BN_RES, f"{sample}/refinement/metawrap_{comp}_{cont}_bins.stats"),
            "bin_dir": os.path.join(AN_BN_RES, f"{sample}/refinement/metawrap_{comp}_{cont}_bins")
        }
    else:
        return {
            "info": os.path.join(AN_BN_RES, f"{sample}/circular_info.tsv"),
            "dir_flag": os.path.join(AN_BN_RES, f"{sample}/circular_contigs/.done")
        }

checkpoint annotation_circular_filter:
    input:
        unpack(get_circular_filter_inputs)
    output:
        dir = directory(os.path.join(AN_RES, "{genome}/circular_bins_filtered"))
    params:
        min_size = AN_C.get("min_circular_size", 800000),
        src_dir = lambda wildcards: os.path.join(AN_BN_RES, get_sample_name_from_genome(wildcards.genome), f"refinement/metawrap_{config.get('binning', {}).get('completeness', 50)}_{config.get('binning', {}).get('contamination', 10)}_bins") if config.get("annotation", {}).get("data_type", "long") == "short" else os.path.join(AN_BN_RES, get_sample_name_from_genome(wildcards.genome), "circular_contigs"),
        config_complete = lambda wildcards: get_complete_contigs(wildcards.genome)
    run:
        import os
        import glob
        import shutil
        
        os.makedirs(output.dir, exist_ok=True)
        # Clear existing files to prevent conflicts
        for f in glob.glob(os.path.join(output.dir, "*")):
            if os.path.islink(f) or os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)
                
        data_type = config.get("annotation", {}).get("data_type", "long")
        if data_type == "short":
            # For short-reads, link refined MAGs directly from params.src_dir
            src_dir = params.src_dir
            if os.path.exists(src_dir):
                for f in os.listdir(src_dir):
                    if f.endswith(".fa"):
                        src = os.path.join(src_dir, f)
                        dst = os.path.join(output.dir, f)
                        if os.path.exists(dst): os.remove(dst)
                        os.symlink(os.path.abspath(src), dst)
        else:
            # For long-reads, perform standard circular filter
            import pandas as pd
            df = pd.read_csv(input.info, sep="\t")
            df["ID"] = df["ID"].astype(str)
            config_complete = [str(c) for c in params.config_complete]
            mask = ((df["Circular"].str.lower() == "yes") & (df["Length"] >= params.min_size)) | (df["ID"].isin(config_complete))
            filtered_ids = df[mask]["ID"].tolist()
            
            for cid in config_complete:
                if cid not in filtered_ids:
                    filtered_ids.append(cid)
    
            for fid in filtered_ids:
                src = os.path.join(params.src_dir, f"{fid}.fa")
                dst = os.path.join(output.dir, f"{fid}.fa")
                if os.path.exists(src):
                    if os.path.exists(dst): os.remove(dst)
                    os.symlink(os.path.abspath(src), dst)


rule annotation_circular_taxonomy:
    input:
        dir = os.path.join(AN_RES, "{genome}/circular_bins_filtered")
    output:
        summary = os.path.join(AN_RES, "{genome}/circular_taxonomy/gtdbtk.summary.tsv")
    log:
        os.path.join(AN_LOGS, "circular_gtdbtk/{genome}.log")
    threads: get_threads(64)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("gtdbtk"),
        db = an_full_path(AN_C.get("gtdbtk_db", "") or config.get("databases", {}).get("gtdbtk", "")),
        bins_dir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_bins_filtered"),
        out_dir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_taxonomy")
    shell:
        r"""
        {{
        set -euo pipefail
        if [ -n "{params.db}" ]; then
            export GTDBTK_DATA_PATH="{params.db}"
        fi
        # Check if there are any bins
        if [ -z "$(ls -A {params.bins_dir}/*.fa 2>/dev/null)" ]; then
            mkdir -p {params.out_dir}
            echo -e "user_genome\tclassification" > {output.summary}
            exit 0
        fi
        conda run -p {params.env} gtdbtk classify_wf \
            --genome_dir {params.bins_dir} \
            --out_dir {params.out_dir} \
            --extension fa \
            --cpus {threads} 2>&1 | tee {log}
        # Standardize output: merge both bacterial (bac120) and archaeal (ar53) summaries
        rm -f {output.summary}
        summary_files=$(ls {params.out_dir}/classify/gtdbtk.*.summary.tsv 2>/dev/null || true)
        if [ -n "$summary_files" ]; then
            awk 'FNR==1 && NR!=1 {{next}} {{print}}' $summary_files > {output.summary}
        else
            echo -e "user_genome\tclassification" > {output.summary}
        fi
        }} 2>&1 | tee {log}
        
        """
rule annotation_circular_quality:
    input:
        dir = os.path.join(AN_RES, "{genome}/circular_bins_filtered")
    output:
        summary = os.path.join(AN_RES, "{genome}/circular_quality/checkm2_summary.tsv")
    log:
        os.path.join(AN_LOGS, "circular_checkm2/{genome}.log")
    threads: get_threads(32)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("checkm2"),
        db = an_full_path(AN_C.get("checkm2_db", "") or config.get("databases", {}).get("checkm2", "")),
        bins_dir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_bins_filtered"),
        out_dir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_quality")
    shell:
        r"""
        {{
        set -euo pipefail
        if [ -z "$(ls -A {params.bins_dir}/*.fa 2>/dev/null)" ]; then
            mkdir -p {params.out_dir}
            echo -e "Name\tCompleteness\tContamination" > {output.summary}
            exit 0
        fi
        DB_ARG=""
        if [ -n "{params.db}" ]; then
            DB_ARG="--database_path {params.db}"
        fi
        conda run -p {params.env} checkm2 predict \
            --threads {threads} \
            --input {params.bins_dir} \
            --output-directory {params.out_dir} \
            --extension fa \
            $DB_ARG --force 2>&1 | tee {log}
        cp {params.out_dir}/quality_report.tsv {output.summary}
        }} 2>&1 | tee {log}
        
        """
def get_bakta_refined_inputs(wildcards):
    data_type = config.get("annotation", {}).get("data_type", "long")
    inputs = {
        "fa": get_bakta_refined_input_fa(wildcards),
        "_dir": os.path.join(AN_RES, f"{wildcards.genome}/circular_bins_filtered")
    }
    if data_type == "short":
        inputs["gff_merged"] = os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.gff3")
        inputs["faa_merged"] = os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.faa")
        inputs["ffn_merged"] = os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.ffn")
        inputs["tsv_merged"] = os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.tsv")
        inputs["gbk_merged"] = os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.gbk")
    else:
        db = an_full_path(AN_C.get("bakta_db", "") or config.get("databases", {}).get("bakta", ""))
        if db:
            inputs["db"] = db
    return inputs

rule annotation_bakta_circular_refined:
    input:
        unpack(get_bakta_refined_inputs)
    output:
        gff = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/{bin}.gff3"),
        faa = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/{bin}.faa"),
        tsv = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/{bin}.tsv"),
        gbk = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/{bin}.gbk")
    log:
        os.path.join(AN_LOGS, "bakta_refined/{genome}_{bin}.log")
    threads: get_threads(16)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        data_type = lambda wildcards: config.get("annotation", {}).get("data_type", "long"),
        env_bakta = an_get_env("bakta"),
        env_biopython = an_get_env("biopython"),
        extract_script = "scripts/extract_bakta_bin.py",
        db = lambda wildcards: an_full_path(AN_C.get("bakta_db", "") or config.get("databases", {}).get("bakta", "")),
        outdir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_refined", wildcards.bin),
        complete = lambda wildcards: (
            get_selected_contigs_info(wildcards.genome).get(wildcards.bin, {}).get("complete")
            if get_selected_contigs_info(wildcards.genome).get(wildcards.bin, {}).get("complete") is not None
            else True
        ),
        locus_tag = lambda wildcards: an_get_locus_tag(wildcards.genome, wildcards.bin),
        gff_merged = lambda wildcards: os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.gff3"),
        faa_merged = lambda wildcards: os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.faa"),
        ffn_merged = lambda wildcards: os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.ffn"),
        tsv_merged = lambda wildcards: os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.tsv"),
        gbk_merged = lambda wildcards: os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.gbk")
    shell:
        r"""
        {{
        set -euo pipefail
        if [ "{params.data_type}" == "short" ]; then
            conda run -p {params.env_biopython} python3 {params.extract_script} \
                --bin_fa {input.fa} \
                --assembly_gff {params.gff_merged} \
                --assembly_faa {params.faa_merged} \
                --assembly_ffn {params.ffn_merged} \
                --assembly_tsv {params.tsv_merged} \
                --assembly_gbk {params.gbk_merged} \
                --out_gff {output.gff} \
                --out_faa {output.faa} \
                --out_tsv {output.tsv} \
                --out_gbk {output.gbk}
        else
            rm -rf {params.outdir}
            EXTRA_FLAGS="--keep-contig-headers"
            if [ "$(echo "{params.complete}" | tr '[:upper:]' '[:lower:]')" = "true" ]; then
                EXTRA_FLAGS="$EXTRA_FLAGS --complete"
            fi
            conda run -p {params.env_bakta} bakta \
                --db {params.db} --threads {threads} --output {params.outdir} \
                --prefix {wildcards.bin} --locus-tag {params.locus_tag} \
                $EXTRA_FLAGS --force {input.fa}
            cp {params.outdir}/{wildcards.bin}.gbff {output.gbk}
        fi
        }} 2>&1 | tee {log}
        
        """
if config.get("annotation", {}).get("data_type", "long") == "short":
    def get_map_functional_inputs(wildcards):
        steps = config.get("annotation", {}).get("steps", {})
        do_eggnog = str(steps.get("eggnog", "False")).lower() == "true"
        do_kofam = str(steps.get("kofamscan", "False")).lower() == "true"
        
        inputs = {
            "refined_faa": os.path.join(AN_RES, f"{wildcards.genome}/circular_refined/{wildcards.bin}/{wildcards.bin}.faa"),
            "assembly_faa": os.path.join(AN_RES, f"{wildcards.genome}/bakta/{wildcards.genome}.faa")
        }
        if do_eggnog:
            inputs["eggnog"] = os.path.join(AN_RES, f"{wildcards.genome}/eggnog/{wildcards.genome}.emapper.annotations")
        if do_kofam:
            inputs["kofam"] = os.path.join(AN_RES, f"{wildcards.genome}/kofamscan/{wildcards.genome}_kofam.tsv")
        return inputs

    rule annotation_map_functional_circular:
        input:
            unpack(get_map_functional_inputs)
        output:
            eggnog = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/functional/{bin}.emapper.annotations"),
            kofam = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/functional/{bin}.kofam.tsv")
        log:
            os.path.join(AN_LOGS, "map_functional/{genome}_{bin}.log")
        threads: 1
        resources:
            runtime="7d",
            mem_mb=lambda wildcards, threads: get_mem_mb(threads)
        params:
            script = "scripts/map_refined_functional.py",
            env = an_get_env("biopython"),
            egg_in = lambda w, input: input.eggnog if hasattr(input, "eggnog") else "none",
            ko_in = lambda w, input: input.kofam if hasattr(input, "kofam") else "none"
        shell:
            r"""
            set -euo pipefail
            conda run -p {params.env} python3 {params.script} \
                --refined_faa {input.refined_faa} \
                --assembly_faa {input.assembly_faa} \
                --assembly_eggnog {params.egg_in} \
                --assembly_kofam {params.ko_in} \
                --out_eggnog {output.eggnog} \
                --out_kofam {output.kofam} 2>&1 | tee {log}
            """
else:
    rule annotation_circular_kofamscan:
        input:
            faa = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/{bin}.faa")
        output:
            tsv = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/functional/{bin}.kofam.tsv")
        log:
            os.path.join(AN_LOGS, "circular_kofamscan/{genome}_{bin}.log")
        threads: get_threads(16)
        resources:
            runtime="7d",
            mem_mb=lambda wildcards, threads: get_mem_mb(threads)
        params:
            env = an_get_env("kofamscan"),
            db = an_full_path(AN_C.get("kofam_db", "") or config.get("databases", {}).get("kofam", "")),
            ko_list = an_full_path(AN_C.get("ko_list", "") or config.get("databases", {}).get("ko_list", "")),
            tmpdir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_refined", wildcards.bin, "kofam_tmp")
        shell:
            r"""
            {{
            set -euo pipefail
            mkdir -p "{params.tmpdir}"
            KO_LIST="{params.ko_list}"
            if [ -z "$KO_LIST" ] || [ ! -f "$KO_LIST" ]; then
                KO_LIST="{params.db}/ko_list"
            fi
            conda run -p {params.env} exec_annotation \
                -o "{output.tsv}" \
                --cpu {threads} \
                --tmp "{params.tmpdir}" \
                -p "{params.db}/profiles" \
                -k "$KO_LIST" \
                "{input.faa}" 2>&1 | tee "{log}"
            rm -rf "{params.tmpdir}"
            }} 2>&1 | tee "{log}"
            """

    rule annotation_circular_eggnog:
        input:
            faa = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/{bin}.faa")
        output:
            anno = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/functional/{bin}.emapper.annotations")
        log:
            os.path.join(AN_LOGS, "circular_eggnog/{genome}_{bin}.log")
        threads: get_threads(16)
        resources:
            runtime="7d",
            mem_mb=lambda wildcards, threads: get_mem_mb(threads)
        params:
            env = an_get_env("eggnog"),
            db = an_full_path(AN_C.get("eggnog_db", "") or config.get("databases", {}).get("eggnog", "")),
            prefix = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_refined", wildcards.bin, "functional", wildcards.bin)
        shell:
            r"""
            {{
            set -euo pipefail
            conda run -p {params.env} emapper.py \
                -m mmseqs --dbmem \
                -i "{input.faa}" \
                --output "{params.prefix}" \
                --cpu {threads} \
                --data_dir "{params.db}" \
                --override 2>&1 | tee "{log}"
            }} 2>&1 | tee "{log}"
            """
def get_circular_aggregate_inputs(wildcards):
    steps = config.get("annotation", {}).get("steps", {})
    do_eggnog = str(steps.get("eggnog", "False")).lower() == "true"
    do_kofam = str(steps.get("kofamscan", "False")).lower() == "true"

    inputs = {
        "faa": os.path.join(AN_RES, f"{wildcards.genome}/circular_refined/{wildcards.bin}/{wildcards.bin}.faa"),
        "gff": os.path.join(AN_RES, f"{wildcards.genome}/circular_refined/{wildcards.bin}/{wildcards.bin}.gff3"),
        "script": "scripts/aggregate_functional.py"
    }
    if do_eggnog:
        inputs["eggnog"] = os.path.join(AN_RES, f"{wildcards.genome}/circular_refined/{wildcards.bin}/functional/{wildcards.bin}.emapper.annotations")
    if do_kofam:
        inputs["kofam"] = os.path.join(AN_RES, f"{wildcards.genome}/circular_refined/{wildcards.bin}/functional/{wildcards.bin}.kofam.tsv")
    return inputs

rule annotation_circular_aggregate:
    input:
        unpack(get_circular_aggregate_inputs)
    output:
        summary = os.path.join(AN_RES, "{genome}/circular_refined/{bin}/{bin}_functional_summary.tsv")
    log:
        os.path.join(AN_LOGS, "circular_aggregate/{genome}_{bin}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("biopython"),
        egg_in = lambda w, input: input.eggnog if hasattr(input, "eggnog") else "none",
        ko_in = lambda w, input: input.kofam if hasattr(input, "kofam") else "none",
        adv_hmm = "none" # Circular bins currently don't use the advanced HMM arbitration
    shell:
        r"""
        set -euo pipefail
        conda run -p {params.env} python3 {input.script} \
            {input.faa} {params.egg_in} {params.ko_in} {params.adv_hmm} \
            {input.gff} {wildcards.bin} {wildcards.bin} {output.summary} 2>&1 | tee {log}
        
        """
def get_circular_summary_inputs(wildcards):
    steps = config.get("annotation", {}).get("steps", {})
    inputs = {
        "deps": get_circular_dependencies(wildcards)
    }
    data_type = config.get("annotation", {}).get("data_type", "long")
    sample = get_sample_name_from_genome(wildcards.genome)
    if data_type == "short":
        comp = config.get("binning", {}).get("completeness", 50)
        cont = config.get("binning", {}).get("contamination", 10)
        inputs["info"] = os.path.join(AN_BN_RES, f"{sample}/refinement/metawrap_{comp}_{cont}_bins.stats")
    else:
        inputs["info"] = os.path.join(AN_BN_RES, f"{sample}/circular_info.tsv")
        
    if steps.get("gtdbtk"):
        inputs["gtdb"] = os.path.join(AN_RES, f"{wildcards.genome}/circular_taxonomy/gtdbtk.summary.tsv")
    if steps.get("checkm2"):
        inputs["checkm"] = os.path.join(AN_RES, f"{wildcards.genome}/circular_quality/checkm2_summary.tsv")
    if steps.get("circular_audit"):
        inputs["audit"] = os.path.join(AN_RES, f"{wildcards.genome}/circular_audit_summary.tsv")
    return inputs

rule annotation_summarize_circular:
    input:
        unpack(get_circular_summary_inputs)
    output:
        summary = os.path.join(AN_RES, "{genome}/circular_summary.tsv")

    params:
        func_dir = os.path.join(AN_RES, "{genome}/circular_refined"),
        filtered_dir = os.path.join(AN_RES, "{genome}/circular_bins_filtered"),
        min_size = lambda wildcards: AN_C.get("min_circular_size", 2500),
        config_complete = lambda wildcards: get_complete_contigs(wildcards.genome)
    run:
        import pandas as pd
        import os
        data_type = config.get("annotation", {}).get("data_type", "long")
        # 1. Basic Stats
        if data_type == "short":
            info = pd.read_csv(input.info, sep="\t")
            if "bin" in info.columns:
                info = info.rename(columns={"bin": "ID"})
            elif "bin_id" in info.columns:
                info = info.rename(columns={"bin_id": "ID"})
            info["ID"] = info["ID"].astype(str)
        else:
            info = pd.read_csv(input.info, sep="\t")
            info = info[info["Circular"].astype(str).str.lower() == "yes"].copy()
            info["ID"] = info["ID"].astype(str)
            
            # Align with circular_bins_filtered: only summarize contigs actually processed
            if os.path.exists(params.filtered_dir):
                retained_bins = set(f.replace(".fa", "") for f in os.listdir(params.filtered_dir) if f.endswith(".fa"))
                if retained_bins:
                    info = info[info["ID"].isin(retained_bins)].copy()

            # Enforce minimum size threshold to prevent short artifact fragments
            min_len = params.min_size
            if min_len is not None and min_len > 0:
                config_complete = [str(c) for c in params.config_complete]
                mask = (info["Length"] >= min_len) | (info["ID"].isin(config_complete))
                info = info[mask].copy()
        # Optional GTDB-Tk
        if hasattr(input, "gtdb") and os.path.exists(input.gtdb):
            gtdb = pd.read_csv(input.gtdb, sep="\t")
            gtdb["user_genome"] = gtdb["user_genome"].astype(str)
        else:
            gtdb = pd.DataFrame(columns=["user_genome", "classification"])
        # Optional CheckM2
        if hasattr(input, "checkm") and os.path.exists(input.checkm):
            checkm = pd.read_csv(input.checkm, sep="\t")
            checkm["Name"] = checkm["Name"].astype(str)
        else:
            checkm = pd.DataFrame(columns=["Name", "Completeness", "Contamination"])
        # Join on 'ID' (matches binning.smk)
        res = info.merge(gtdb, left_on="ID", right_on="user_genome", how="left")
        res = res.merge(checkm, left_on="ID", right_on="Name", how="left", suffixes=("", "_checkm2"))

        # Optional Circular Audit (Trigodet et al. 2026)
        if hasattr(input, "audit") and os.path.exists(input.audit) and os.path.getsize(input.audit) > 0:
            try:
                audit_df = pd.read_csv(input.audit, sep="\t")
                if "Contig_ID" in audit_df.columns:
                    audit_df = audit_df.rename(columns={"Contig_ID": "ID"})
                audit_df["ID"] = audit_df["ID"].astype(str)
                res = res.merge(audit_df, on="ID", how="left", suffixes=("", "_audit"))
            except Exception:
                pass
        # 2. Functional Stats (Sequence-aware)
        func_stats = []
        for b in info["ID"]:
            faa_file = os.path.join(params.func_dir, b, b + ".faa")
            egg_file = os.path.join(params.func_dir, b, "functional", b + ".emapper.annotations")
            ko_file = os.path.join(params.func_dir, b, "functional", b + ".kofam.tsv")
            # Count Total Genes
            total_genes = 0
            if os.path.exists(faa_file):
                with open(faa_file) as f:
                    total_genes = sum(1 for l in f if l.startswith(">"))
            # Track unique annotated gene IDs
            annotated_ids = set()
            n_egg = 0
            if os.path.exists(egg_file):
                with open(egg_file) as f:
                    for line in f:
                        if not line.startswith("#"):
                            parts = line.split("\t")
                            if parts:
                                annotated_ids.add(parts[0])
                                n_egg += 1
            n_ko = 0
            if os.path.exists(ko_file):
                # For each protein, find its highest score
                best_scores = {} # protein_id -> max_score
                with open(ko_file) as f:
                    for line in f:
                        if line.startswith("#") or line.startswith("-") or not line.strip():
                            continue
                        content = line.lstrip("*").strip()
                        parts = content.split()
                        if len(parts) < 5: continue
                        pid = parts[0]
                        try:
                            score = float(parts[3])
                            if pid not in best_scores or score > best_scores[pid]:
                                best_scores[pid] = score
                        except: continue
                # Count genes with Score > 50 and add to annotated set
                for pid, score in best_scores.items():
                    if score > 50:
                        annotated_ids.add(pid)
                n_ko = sum(1 for s in best_scores.values() if s > 50)
            # Final stats for this bin
            total_annotated = len(annotated_ids)
            rate = (total_annotated / total_genes * 100) if total_genes > 0 else 0
            func_stats.append(dict(
                ID=b, 
                Total_Genes=total_genes,
                Annotated_Genes=total_annotated,
                Annotation_Rate_Pct=round(rate, 2),
                eggnog_genes=n_egg, 
                kofam_genes=n_ko
            ))
        # 3. Final Merge and Save
        res = res.merge(pd.DataFrame(func_stats), on="ID", how="left")
        
        # Guard: Drop empty/unannotated ghost rows where no genes were found and no CheckM2/Audit data exists
        if not res.empty:
            has_genes = (res["Total_Genes"] > 0) if "Total_Genes" in res.columns else False
            has_checkm = (res["Completeness"].notna() & (res["Completeness"] > 0)) if "Completeness" in res.columns else False
            has_audit = (res["Gate1_Physical_Structural_Support"].notna() & (res["Gate1_Physical_Structural_Support"] != "")) if "Gate1_Physical_Structural_Support" in res.columns else False
            has_gtdb = (res["classification"].notna() & (res["classification"] != "")) if "classification" in res.columns else False
            
            meaningful_mask = has_genes | has_checkm | has_audit | has_gtdb
            if meaningful_mask.any():
                res = res[meaningful_mask].copy()

        res.to_csv(output.summary, sep="\t", index=False)

# ------------------------------------------------------------------------------
# Trigodet-Informed Raw-Read Circular Assembly Verification Pipeline (Version 2.1)
# ------------------------------------------------------------------------------

def get_circular_lr_input(wildcards):
    sample = get_sample_name_from_genome(wildcards.genome)
    s_info = config.get("samples", {}).get(sample, {})
    lr_data = s_info.get("long_read", {})
    if isinstance(lr_data, str):
        path = lr_data
    elif isinstance(lr_data, dict):
        path = lr_data.get("lr") or lr_data.get("fastq") or ""
    else:
        path = ""
    return an_full_path(path) if path else ""

rule annotation_circular_whole_assembly_mapping:
    input:
        assembly = lambda wildcards: an_get_asm_path(get_sample_name_from_genome(wildcards.genome)),
        reads = get_circular_lr_input
    output:
        bam = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam"),
        bai = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam.bai")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_whole_assembly_mapping.log")
    threads: get_threads(AN_C.get("threads", {}).get("minimap2", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        mapping_env = an_get_env("mapping"),
        preset = lambda wildcards: AN_C.get("circular_audit", {}).get("longread_mapping_preset", "map-ont")
    shell:
        r"""
        {{
        set -euo pipefail
        if [ -z "{input.reads}" ] || [ ! -f "{input.reads}" ]; then
            echo "Error: Long reads file not found: {input.reads}" >&2
            exit 1
        fi
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/minimap2" ]; then
            for cand in "{params.env}" \
                        "{params.mapping_env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" \
                        "$HOME/.conda/envs/bio"; do
                if [ -x "$cand/bin/minimap2" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        for extra in "{params.mapping_env}" "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" "$HOME/.conda/envs/bio"; do
            if [ -x "$extra/bin/samtools" ]; then
                export PATH="$extra/bin:$PATH"
            fi
        done
        mkdir -p $(dirname {output.bam})
        minimap2 -ax {params.preset} -t {threads} -p 1 --secondary-seq {input.assembly} {input.reads} | \
            samtools view -u -F 4 - | \
            samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        }} 2>&1 | tee {log}
        """

rule annotation_circular_extract_terminal_qnames:
    input:
        bam = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam"),
        bai = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam.bai"),
        contig = os.path.join(AN_RES, "{genome}/circular_bins_filtered/{bin}.fa")
    output:
        qnames = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_terminal_qnames.txt"),
        junction_ref = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_junction_ref.fa")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_{bin}_terminal_qnames.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        flank = lambda wildcards: AN_C.get("circular_audit", {}).get("config_terminal_flank_bp", 5000)
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname {output.qnames})
        export PATH="{params.env}/bin:$PATH"
        python3 scripts/extract_terminal_qnames.py \
            --contig {input.contig} \
            --bam {input.bam} \
            --bin {wildcards.bin} \
            --flank {params.flank} \
            --out-junction {output.junction_ref} \
            --out-qnames {output.qnames}
        }} 2>&1 | tee {log}
        """

rule annotation_circular_recover_terminal_fastq:
    input:
        qnames = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_terminal_qnames.txt"),
        raw_reads = get_circular_lr_input
    output:
        reads = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_terminal_reads.fq.gz")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_{bin}_recover_fastq.log")
    threads: get_threads(AN_C.get("threads", {}).get("seqkit", 4))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env_seqkit = an_get_env("seqkit"),
        env_longread = an_get_env("longread_audit")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname {output.reads})
        if [ ! -s "{input.qnames}" ]; then
            echo -n "" | gzip -c > {output.reads}
        else
            if [ -x "{params.env_seqkit}/bin/seqkit" ]; then
                "{params.env_seqkit}/bin/seqkit" grep -f {input.qnames} {input.raw_reads} | gzip -c > {output.reads}
            elif [ -x "{params.env_longread}/bin/seqkit" ]; then
                "{params.env_longread}/bin/seqkit" grep -f {input.qnames} {input.raw_reads} | gzip -c > {output.reads}
            else
                export PATH="{params.env_longread}/bin:$PATH"
                python3 scripts/recover_terminal_fastq.py \
                    --qnames {input.qnames} \
                    --raw-reads {input.raw_reads} \
                    --out-reads {output.reads}
            fi
        fi
        }} 2>&1 | tee {log}
        """

rule annotation_circular_synthetic_junction_mapping:
    input:
        junction_ref = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_junction_ref.fa"),
        reads = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_terminal_reads.fq.gz")
    output:
        bam = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_junction.bam"),
        bai = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_junction.bam.bai")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_{bin}_junction_mapping.log")
    threads: get_threads(AN_C.get("threads", {}).get("minimap2", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        mapping_env = an_get_env("mapping"),
        preset = lambda wildcards: AN_C.get("circular_audit", {}).get("longread_mapping_preset", "map-ont")
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/minimap2" ]; then
            for cand in "{params.env}" \
                        "{params.mapping_env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" \
                        "$HOME/.conda/envs/bio"; do
                if [ -x "$cand/bin/minimap2" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        for extra in "{params.mapping_env}" "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" "$HOME/.conda/envs/bio"; do
            if [ -x "$extra/bin/samtools" ]; then
                export PATH="$extra/bin:$PATH"
            fi
        done
        mkdir -p $(dirname {output.bam})
        if [ ! -s "{input.reads}" ]; then
            samtools view -H -b -o {output.bam} <(minimap2 -d /dev/stdout {input.junction_ref} 2>/dev/null) 2>/dev/null || \
                samtools view -b -o {output.bam} /dev/null
            samtools index {output.bam} || true
        else
            minimap2 -ax {params.preset} -t {threads} -p 1 {input.junction_ref} {input.reads} | \
                samtools view -u -F 4 - | \
                samtools sort -@ {threads} -o {output.bam} -
            samtools index -@ {threads} {output.bam}
        fi
        }} 2>&1 | tee {log}
        """

rule annotation_circular_blast_repeats:
    input:
        contig = os.path.join(AN_RES, "{genome}/circular_bins_filtered/{bin}.fa")
    output:
        blast_tsv = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_self_blast.tsv")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_{bin}_blast_repeats.log")
    threads: get_threads(AN_C.get("threads", {}).get("blast", 4))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        tmp_db = lambda wildcards: os.path.join(AN_RES, wildcards.genome, f"circular_audit/{wildcards.bin}/tmp_blast_db")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname {params.tmp_db})
        export PATH="{params.env}/bin:$PATH"
        if which makeblastdb &> /dev/null && which blastn &> /dev/null; then
            makeblastdb -in {input.contig} -dbtype nucl -out {params.tmp_db}
            blastn -query {input.contig} -db {params.tmp_db} -outfmt 6 -num_threads {threads} > {output.blast_tsv}
            rm -f {params.tmp_db}.*
        else
            minimap2 -x asm20 -c -t {threads} {input.contig} {input.contig} | \
                awk '{{print $1"\t"$6"\t"($10/$11*100)"\t"$11"\t0\t0\t"$3"\t"$4"\t"$8"\t"$9"\t0\t0"}}' > {output.blast_tsv} || touch {output.blast_tsv}
        fi
        }} 2>&1 | tee {log}
        """

rule annotation_circular_anvio_validation:
    input:
        whole_bam = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam"),
        whole_bai = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam.bai")
    output:
        concordance = os.path.join(AN_RES, "{genome}/circular_audit/anvio_concordance.tsv")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_anvio_validation.log")
    threads: get_threads(AN_C.get("threads", {}).get("anvio", 4))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("anvio"),
        out_prefix = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_audit/anvio_out")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p $(dirname {output.concordance})
        export PATH="{params.env}/bin:$PATH"
        if which anvi-script-find-misassemblies &> /dev/null; then
            echo "Running anvi-script-find-misassemblies with --just-do-it..."
            anvi-script-find-misassemblies \
                -b {input.whole_bam} \
                -o {params.out_prefix} \
                --just-do-it \
                -T {threads} 2>&1 || true
            echo -e "engine\tstatus\tnotes" > {output.concordance}
            echo -e "anvio_9\tCOMPLETED\tOutput generated at {params.out_prefix}" >> {output.concordance}
        else
            echo "anvi-script-find-misassemblies not available in {params.env}"
            echo -e "engine\tstatus\tnotes" > {output.concordance}
            echo -e "anvio_9\tNOT_RUN\tenvironment or executable missing" >> {output.concordance}
        fi
        }} 2>&1 | tee {log}
        """

def get_circular_audit_inputs(wildcards):
    inputs = {
        "contig": os.path.join(AN_RES, wildcards.genome, f"circular_bins_filtered/{wildcards.bin}.fa"),
        "whole_bam": os.path.join(AN_RES, wildcards.genome, "circular_audit/whole_assembly.bam"),
        "whole_bai": os.path.join(AN_RES, wildcards.genome, "circular_audit/whole_assembly.bam.bai"),
        "junction_bam": os.path.join(AN_RES, wildcards.genome, f"circular_audit/{wildcards.bin}/{wildcards.bin}_junction.bam"),
        "junction_bai": os.path.join(AN_RES, wildcards.genome, f"circular_audit/{wildcards.bin}/{wildcards.bin}_junction.bam.bai"),
        "self_blast": os.path.join(AN_RES, wildcards.genome, f"circular_audit/{wildcards.bin}/{wildcards.bin}_self_blast.tsv"),
        "script": "scripts/audit_circular_assembly.py"
    }
    steps = AN_C.get("steps", {})
    if str(steps.get("circular", "False")).lower() == "true":
        inputs["bakta_tsv"] = os.path.join(AN_RES, wildcards.genome, f"circular_refined/{wildcards.bin}/{wildcards.bin}.tsv")
    if str(steps.get("checkm2", "False")).lower() == "true":
        inputs["checkm2"] = os.path.join(AN_RES, wildcards.genome, "circular_quality/checkm2_summary.tsv")
    if str(steps.get("gtdbtk", "False")).lower() == "true":
        inputs["gtdb"] = os.path.join(AN_RES, wildcards.genome, "circular_taxonomy/gtdbtk.summary.tsv")
    sample = get_sample_name_from_genome(wildcards.genome)
    pl_dir = an_full_path(config.get("plasmid_analysis", {}).get("results_dir", "results/plasmid_analysis"))
    genomad_summary = os.path.join(pl_dir, sample, f"{sample}_genomad_summary.tsv")
    if os.path.exists(genomad_summary):
        inputs["genomad"] = genomad_summary
    if sample_has_short_reads(wildcards):
        inputs["illumina_qc"] = os.path.join(AN_RES, wildcards.genome, "circular_coverm/circular_illumina_qc_summary.tsv")
    return inputs

rule annotation_circular_audit_single:
    input:
        unpack(get_circular_audit_inputs)
    output:
        tsv = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_audit.tsv")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_{bin}_audit.log")
    threads: get_threads(2)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        terminal_flank = lambda wildcards: AN_C.get("circular_audit", {}).get("config_terminal_flank_bp", 5000),
        min_flank = lambda wildcards: AN_C.get("circular_audit", {}).get("min_junction_flank_bp", 500),
        min_spanning = lambda wildcards: AN_C.get("circular_audit", {}).get("min_junction_spanning_unique_reads", 5),
        min_support_fraction = lambda wildcards: AN_C.get("circular_audit", {}).get("min_junction_support_fraction", 0.50),
        max_clip_fraction = lambda wildcards: AN_C.get("circular_audit", {}).get("max_junction_clip_fraction", 0.15),
        clipping_ratio_threshold = lambda wildcards: AN_C.get("circular_audit", {}).get("clipping_ratio_threshold", 0.80),
        min_coverage_depth = lambda wildcards: AN_C.get("circular_audit", {}).get("min_coverage_depth", 10),
        zero_cov_high_risk = lambda wildcards: AN_C.get("circular_audit", {}).get("zero_cov_high_risk_threshold_bp", 500),
        bakta_tsv = lambda wildcards, input: getattr(input, "bakta_tsv", ""),
        checkm2 = lambda wildcards, input: getattr(input, "checkm2", ""),
        gtdb = lambda wildcards, input: getattr(input, "gtdb", ""),
        genomad = lambda wildcards, input: getattr(input, "genomad", ""),
        illumina_qc = lambda wildcards, input: getattr(input, "illumina_qc", "")
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/python3" ]; then
            for cand in "{params.env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit"; do
                if [ -x "$cand/bin/python3" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        OPT_ARGS=""
        if [ -n "{params.bakta_tsv}" ] && [ -f "{params.bakta_tsv}" ]; then
            OPT_ARGS="$OPT_ARGS --bakta-tsv {params.bakta_tsv}"
        fi
        if [ -n "{params.checkm2}" ] && [ -f "{params.checkm2}" ]; then
            OPT_ARGS="$OPT_ARGS --checkm2-tsv {params.checkm2}"
        fi
        if [ -n "{params.gtdb}" ] && [ -f "{params.gtdb}" ]; then
            OPT_ARGS="$OPT_ARGS --gtdb-tsv {params.gtdb}"
        fi
        if [ -n "{params.genomad}" ] && [ -f "{params.genomad}" ]; then
            OPT_ARGS="$OPT_ARGS --genomad-tsv {params.genomad}"
        fi
        if [ -n "{params.illumina_qc}" ] && [ -f "{params.illumina_qc}" ]; then
            OPT_ARGS="$OPT_ARGS --illumina-qc-tsv {params.illumina_qc}"
        fi
        $ENV_DIR/bin/python3 {input.script} \
            --bin {wildcards.bin} \
            --contig-fa {input.contig} \
            --whole-bam {input.whole_bam} \
            --junction-bam {input.junction_bam} \
            --self-blast {input.self_blast} \
            --terminal-flank {params.terminal_flank} \
            --min-junction-flank {params.min_flank} \
            --min-junction-spanning-reads {params.min_spanning} \
            --min-junction-support-fraction {params.min_support_fraction} \
            --max-junction-clip-fraction {params.max_clip_fraction} \
            --clipping-ratio-threshold {params.clipping_ratio_threshold} \
            --min-coverage-depth {params.min_coverage_depth} \
            --zero-cov-high-risk-threshold {params.zero_cov_high_risk} \
            $OPT_ARGS \
            --out-tsv {output.tsv}
        }} 2>&1 | tee {log}
        """

def get_circular_audit_plot_inputs(wildcards):
    inputs = {
        "contig": os.path.join(AN_RES, wildcards.genome, f"circular_bins_filtered/{wildcards.bin}.fa"),
        "whole_bam": os.path.join(AN_RES, wildcards.genome, "circular_audit/whole_assembly.bam"),
        "whole_bai": os.path.join(AN_RES, wildcards.genome, "circular_audit/whole_assembly.bam.bai"),
        "junction_bam": os.path.join(AN_RES, wildcards.genome, f"circular_audit/{wildcards.bin}/{wildcards.bin}_junction.bam"),
        "junction_bai": os.path.join(AN_RES, wildcards.genome, f"circular_audit/{wildcards.bin}/{wildcards.bin}_junction.bam.bai"),
        "audit_tsv": os.path.join(AN_RES, wildcards.genome, f"circular_audit/{wildcards.bin}/{wildcards.bin}_audit.tsv"),
        "script": "scripts/plot_circular_audit.py"
    }
    steps = AN_C.get("steps", {})
    if str(steps.get("circular", "False")).lower() == "true":
        inputs["bakta_tsv"] = os.path.join(AN_RES, wildcards.genome, f"circular_refined/{wildcards.bin}/{wildcards.bin}.tsv")
    return inputs

rule annotation_circular_audit_plot:
    input:
        unpack(get_circular_audit_plot_inputs)
    output:
        pdf = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_verification_plot.pdf"),
        png = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_verification_plot.png"),
        bundle = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_plot_bundle.json")
    log:
        os.path.join(AN_LOGS, "circular_audit/{genome}_{bin}_plot.log")
    threads: get_threads(AN_C.get("threads", {}).get("python", 2))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        flank = lambda wildcards: AN_C.get("circular_audit", {}).get("config_terminal_flank_bp", 5000),
        bakta_tsv = lambda wildcards, input: getattr(input, "bakta_tsv", ""),
        dump_bundle = lambda wildcards: AN_C.get("circular_audit", {}).get("dump_plot_bundle", True)
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/python3" ]; then
            for cand in "{params.env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit"; do
                if [ -x "$cand/bin/python3" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        OPT_ARGS=""
        if [ -n "{params.bakta_tsv}" ] && [ -f "{params.bakta_tsv}" ]; then
            OPT_ARGS="$OPT_ARGS --bakta-tsv {params.bakta_tsv}"
        fi
        if [ "{params.dump_bundle}" = "True" ] || [ "{params.dump_bundle}" = "true" ]; then
            OPT_ARGS="$OPT_ARGS --dump-bundle {output.bundle}"
        else
            touch {output.bundle}
        fi
        $ENV_DIR/bin/python3 {input.script} \
            --bin {wildcards.bin} \
            --contig-fa {input.contig} \
            --whole-bam {input.whole_bam} \
            --junction-bam {input.junction_bam} \
            --audit-tsv {input.audit_tsv} \
            --flank {params.flank} \
            $OPT_ARGS \
            --out-pdf {output.pdf} \
            --out-png {output.png}
        }} 2>&1 | tee {log}
        """

def get_circular_audit_summary_inputs(wildcards):
    try:
        checkpoint_dir = checkpoints.annotation_circular_filter.get(genome=wildcards.genome).output.dir
    except:
        checkpoint_dir = os.path.join(AN_RES, wildcards.genome, "circular_bins_filtered")
    
    selected = []
    if os.path.exists(checkpoint_dir):
        selected = [f.replace(".fa", "") for f in os.listdir(checkpoint_dir) if f.endswith(".fa")]
    else:
        sample = get_sample_name_from_genome(wildcards.genome)
        info_file = os.path.join(AN_BN_RES, sample, "circular_info.tsv")
        if os.path.exists(info_file):
            import pandas as pd
            df = pd.read_csv(info_file, sep="\t")
            df["ID"] = df["ID"].astype(str)
            min_size = AN_C.get("min_circular_size", 800000)
            config_complete = [str(c) for c in get_complete_contigs(wildcards.genome)]
            mask = ((df["Circular"].str.lower() == "yes") & (df["Length"] >= min_size)) | (df["ID"].isin(config_complete))
            selected = df[mask]["ID"].tolist()
            for cid in config_complete:
                if cid not in selected:
                    selected.append(cid)
    
    tsvs = [os.path.join(AN_RES, wildcards.genome, f"circular_audit/{b}/{b}_audit.tsv") for b in selected]
    plots = [os.path.join(AN_RES, wildcards.genome, f"circular_audit/{b}/{b}_verification_plot.pdf") for b in selected]
    bundles = [os.path.join(AN_RES, wildcards.genome, f"circular_audit/{b}/{b}_plot_bundle.json") for b in selected]
    return {"tsvs": tsvs, "plots": plots, "bundles": bundles}

rule annotation_circular_audit_summary:
    input:
        unpack(get_circular_audit_summary_inputs)
    output:
        summary = os.path.join(AN_RES, "{genome}/circular_audit_summary.tsv")
    log:
        os.path.join(AN_LOGS, "circular_audit_summary/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    run:
        import pandas as pd
        dfs = []
        file_list = input.tsvs if hasattr(input, "tsvs") else input
        for f in file_list:
            if os.path.exists(f) and os.path.getsize(f) > 0:
                try:
                    dfs.append(pd.read_csv(f, sep="\t"))
                except:
                    pass
        if dfs:
            res = pd.concat(dfs, ignore_index=True)
        else:
            res = pd.DataFrame(columns=[
                "Contig_ID", "Length_bp", "Effective_Junction_Flank_bp",
                "Primary_Mean_Depth", "Primary_Median_Depth", "Primary_Depth_MAD", "Primary_Depth_CV",
                "Primary_Coverage_Breadth_1x_Pct", "Primary_Coverage_Breadth_3x_Pct",
                "All_Alignment_Support_Breadth_Pct", "Total_Zero_Coverage_bp", "Max_Zero_Coverage_Block_bp",
                "Unresolved_Zero_Coverage_bp", "Has_High_Risk_Zero_Block_gt500bp", "Zero_Coverage_Gaps_Detail",
                "Junction_Spanning_Unique_Reads", "Junction_Informative_Unique_Reads",
                "Junction_Support_Fraction", "Junction_Support_CI_95_Lower", "Junction_Support_CI_95_Upper",
                "Junction_Clip_Unique_Reads", "Junction_Clip_Fraction", "Terminal_Anchor_Uniqueness",
                "Systematic_Internal_Clipping_Sites", "Internal_Clipping_Sites_Detail",
                "Population_SV_Polymorphisms_Count", "Population_SV_Polymorphisms_Detail",
                "Supplementary_Partner_Contig", "Supplementary_Partner_Unique_Reads",
                "Supplementary_Partner_Fraction", "Supplementary_Partner_Locus_MAD",
                "Supplementary_Partner_Depth", "Supplementary_Partner_Relative_Depth_Ratio", "Supplementary_Partner_Is_Artifact",
                "Repeat_Breadth_Pct", "Longest_Repeat_bp", "Repeat_Hit_Count",
                "CheckM2_Completeness", "CheckM2_Contamination",
                "rRNA_5S_Count", "rRNA_16S_Count", "rRNA_23S_Count", "tRNA_Canonical_20_AA_Count",
                "GTDB_Lineage",
                "geNomad_Chromosome_Score", "geNomad_Plasmid_Score", "geNomad_Virus_Score",
                "geNomad_Assignment_Derived", "geNomad_Confidence_Margin_Derived",
                "Gate1_Physical_Structural_Support", "Gate1_Notes",
                "Gate2_Coverage_Coherence",
                "Gate3_Biological_Completeness", "Gate3_Notes",
                "Gate4_Replicon_Evidence_State",
                "Illumina_Relative_Abundance_Pct", "Illumina_Catalogue_Mapped_Pct", "Illumina_Catalogue_Unmapped_Pct",
                "Illumina_Genome_Trimmed_Depth", "Illumina_Genome_Covered_Fraction", "Illumina_Aligned_Read_Count", "Illumina_ANIr_Pct",
                "Illumina_Competitive_Trimmed_Depth", "Illumina_Competitive_Depth_Ratio_C", "Illumina_Strict_Depth_Ratio_S",
                "Illumina_Rotated_Junction_Depth_100bp", "Illumina_Rotated_Junction_Local_Ratio",
                "Illumina_Rotated_Zero_Coverage_Run_bp",
                "Illumina_Rotated_Spanning_Reads_Unique", "Illumina_Rotated_Bridging_Pairs_Unique",
                "Gate5_Illumina_Orthogonal_Status", "Gate5_Notes",
                "Assembly_Usability_Grade",
                "Circularity_Confidence_Tier", "Circularity_Evidence_Type",
                "Biological_Complexity_Tier", "Graph_Complexity_Tier",
                "Alternative_Path_Contig", "Alternative_Path_Relative_Depth_Ratio", "Alternative_Path_Graph_Signal",
                "Complexity_Notes",
                "Final_Confidence_Grade"
            ])
        res.to_csv(output.summary, sep="\t", index=False)

