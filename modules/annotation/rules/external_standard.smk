# Standard rules for external genome annotation (non-batch mode)
import os
import re

# Constrain wildcard 'bin' to only match individual genomes that do not belong to any batch group
# when batch HMM mode is active. This avoids Snakemake creating standard jobs for the batch MAGs.
if AN_C.get("batch_hmm", False):
    individual_bins = [
        name for name, info in EXT_GENOMES.items() 
        if not (isinstance(info, dict) and "batch_group" in info) and not (isinstance(info, dict) and info.get("type") == "dna_dir")
    ]
    wildcard_constraints:
        bin = "|".join(re.escape(b) for b in individual_bins) if individual_bins else "none_matched_placeholder"

def is_ext_split_bakta(bin_key):
    """Checks if the external genome is configured to run bakta in split mode."""
    info = EXT_GENOMES.get(bin_key)
    if isinstance(info, dict):
        return str(info.get("split_bakta", "False")).lower() == "true"
    return False

def get_ext_bakta_inputs(wildcards):
    """Determines whether to use single or split/merged bakta outputs based on split_bakta."""
    bin_key = wildcards.bin
    if is_ext_split_bakta(bin_key):
        return {
            "gff": os.path.join(EXT_RES, f"{bin_key}/merged/{bin_key}.gff3"),
            "faa": os.path.join(EXT_RES, f"{bin_key}/merged/{bin_key}.faa"),
            "ffn": os.path.join(EXT_RES, f"{bin_key}/merged/{bin_key}.ffn"),
            "tsv": os.path.join(EXT_RES, f"{bin_key}/merged/{bin_key}.tsv"),
            "gbk": os.path.join(EXT_RES, f"{bin_key}/merged/{bin_key}.gbk")
        }
    else:
        return {
            "gff": os.path.join(EXT_RES, f"{bin_key}/single/{bin_key}.gff3"),
            "faa": os.path.join(EXT_RES, f"{bin_key}/single/{bin_key}.faa"),
            "ffn": os.path.join(EXT_RES, f"{bin_key}/single/{bin_key}.ffn"),
            "tsv": os.path.join(EXT_RES, f"{bin_key}/single/{bin_key}.tsv"),
            "gbk": os.path.join(EXT_RES, f"{bin_key}/single/{bin_key}.gbk")
        }

rule ext_annotation_bakta:
    input:
        unpack(get_ext_bakta_inputs)
    output:
        gff = os.path.join(EXT_RES, "{bin}/{bin}.gff3"),
        faa = os.path.join(EXT_RES, "{bin}/{bin}.faa"),
        ffn = os.path.join(EXT_RES, "{bin}/{bin}.ffn"),
        tsv = os.path.join(EXT_RES, "{bin}/{bin}.tsv"),
        gbk = os.path.join(EXT_RES, "{bin}/{bin}.gbk")
    log: os.path.join(AN_LOGS, "external_bakta_final/{bin}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        {{
        set -euo pipefail
        cp "{input.gff}" "{output.gff}"
        cp "{input.faa}" "{output.faa}"
        cp "{input.ffn}" "{output.ffn}"
        cp "{input.tsv}" "{output.tsv}"
        cp "{input.gbk}" "{output.gbk}"
        }} 2>&1 | tee {log}
        """

rule ext_annotation_bakta_single:
    input:
        fa = lambda wildcards: EXT_GENOMES[wildcards.bin] if isinstance(EXT_GENOMES[wildcards.bin], str) else EXT_GENOMES[wildcards.bin].get("path")
    output:
        gff = os.path.join(EXT_RES, "{bin}/single/{bin}.gff3"),
        faa = os.path.join(EXT_RES, "{bin}/single/{bin}.faa"),
        ffn = os.path.join(EXT_RES, "{bin}/single/{bin}.ffn"),
        tsv = os.path.join(EXT_RES, "{bin}/single/{bin}.tsv"),
        gbk = os.path.join(EXT_RES, "{bin}/single/{bin}.gbk")
    log: os.path.join(AN_LOGS, "external_bakta_single/{bin}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        db = an_full_path(AN_C.get("bakta_db", "")),
        outdir = os.path.join(EXT_RES, "{bin}/single"),
        complete = lambda wildcards: (EXT_GENOMES[wildcards.bin] if isinstance(EXT_GENOMES[wildcards.bin], dict) else {}).get("complete", False),
        locus_tag = lambda wildcards: get_ext_inherent_name(wildcards.bin)[:24],
        extra_flags = lambda wildcards: AN_C.get("ext_bakta_extra_flags", "")
    shell:
        r"""
        {{
        set -euo pipefail
        
        ulimit -s unlimited || true
        export OMP_STACKSIZE=128M || true
        
        rm -rf "{params.outdir}"
        EXTRA_FLAGS="--keep-contig-headers"
        if [ "{params.complete}" == "True" ]; then
            EXTRA_FLAGS="$EXTRA_FLAGS --complete"
        fi
        conda run -p {params.env} bakta \
            --db {params.db} --threads {threads} --output "{params.outdir}" \
            --prefix {wildcards.bin} --locus-tag {params.locus_tag} \
            {params.extra_flags} \
            $EXTRA_FLAGS --force {input.fa} 2>&1 | tee {log}
        cp "{params.outdir}/{wildcards.bin}.gbff" "{output.gbk}"
        }} 2>&1 | tee {log}
        """

rule ext_annotation_bakta_split:
    input:
        fa = lambda wildcards: EXT_GENOMES[wildcards.bin] if isinstance(EXT_GENOMES[wildcards.bin], str) else EXT_GENOMES[wildcards.bin].get("path")
    output:
        chunks = scatter.an_chunks(os.path.join(EXT_RES, "{{bin}}/chunks/chunk_{scatteritem}.fasta")),
        meta = os.path.join(EXT_RES, "{bin}/chunks/chunk_metadata.tsv")
    log:
        os.path.join(AN_LOGS, "external_bakta_split/{bin}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        outdir = lambda wildcards: os.path.join(EXT_RES, wildcards.bin, "chunks"),
        num = AN_SPLIT_COUNT,
        env = an_get_env("seqkit")
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p "{params.outdir}"
        echo -e "chunk\tis_complete" > "{output.meta}"
        conda run -p {params.env} seqkit split2 --by-part {params.num} {input.fa} -O "{params.outdir}/tmp_remain" --force
        for i in $(seq 1 {params.num}); do
            TARGET_FILE="{params.outdir}/chunk_$i-of-{params.num}.fasta"
            SUFFIX=$(printf "%03d" $i)
            SRC_FILE=$(ls {params.outdir}/tmp_remain/*.part_$SUFFIX.* 2>/dev/null || true)
            if [ -n "$SRC_FILE" ]; then
                mv "$SRC_FILE" "$TARGET_FILE"
            else
                touch "$TARGET_FILE"
            fi
            echo -e "chunk_$i-of-{params.num}\tfalse" >> "{output.meta}"
        done
        rm -rf "{params.outdir}/tmp_remain"
        }} 2>&1 | tee {log}
        """

rule ext_annotation_bakta_chunk:
    input:
        fa = os.path.join(EXT_RES, "{bin}/chunks/chunk_{scatteritem}.fasta"),
        db = AN_C.get("bakta_db", ""),
        meta = os.path.join(EXT_RES, "{bin}/chunks/chunk_metadata.tsv")
    output:
        gff = os.path.join(EXT_RES, "{bin}/chunks/bakta_{scatteritem}/{bin}.gff3"),
        faa = os.path.join(EXT_RES, "{bin}/chunks/bakta_{scatteritem}/{bin}.faa"),
        ffn = os.path.join(EXT_RES, "{bin}/chunks/bakta_{scatteritem}/{bin}.ffn"),
        tsv = os.path.join(EXT_RES, "{bin}/chunks/bakta_{scatteritem}/{bin}.tsv"),
        gbk = os.path.join(EXT_RES, "{bin}/chunks/bakta_{scatteritem}/{bin}.gbk")
    log:
        os.path.join(AN_LOGS, "external_bakta_chunk/{bin}_{scatteritem}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        db = an_full_path(AN_C.get("bakta_db", "")),
        outdir = lambda wildcards, output: os.path.dirname(output.gff),
        locus_tag = lambda wildcards: f"c{wildcards.scatteritem.split('-of-')[0]}",
        extra_flags = lambda wildcards: AN_C.get("ext_bakta_extra_flags", "")
    shell:
        r"""
        {{
        set -euo pipefail
        
        ulimit -s unlimited || true
        export OMP_STACKSIZE=128M || true
        
        IS_COMPLETE=$(grep "chunk_{wildcards.scatteritem}" {input.meta} | cut -f2 || echo "false")
        EXTRA_FLAGS=""
        if [ "$IS_COMPLETE" == "true" ]; then
            EXTRA_FLAGS="--complete"
        fi
        mkdir -p "{params.outdir}"
        conda run -p {params.env} bakta \
            --db {params.db} --threads {threads} --output "{params.outdir}" \
            --prefix "{wildcards.bin}" --locus-tag {params.locus_tag} \
            {params.extra_flags} \
            --keep-contig-headers $EXTRA_FLAGS --force {input.fa} 2>&1 | tee {log}
        cp "{params.outdir}/{wildcards.bin}.gbff" "{output.gbk}"
        }} 2>&1 | tee {log}
        """

def get_ext_bakta_merge_inputs(wildcards):
    """Gathers all chunks to merge for the external genome."""
    bin_key = wildcards.bin
    all_gffs, all_faas, all_ffns, all_gbks = [], [], [], []
    for i in range(1, AN_SPLIT_COUNT + 1):
        sid = f"{i}-of-{AN_SPLIT_COUNT}"
        base = os.path.join(EXT_RES, bin_key, "chunks", f"bakta_{sid}")
        all_gffs.append(os.path.join(base, f"{bin_key}.gff3"))
        all_faas.append(os.path.join(base, f"{bin_key}.faa"))
        all_ffns.append(os.path.join(base, f"{bin_key}.ffn"))
        all_gbks.append(os.path.join(base, f"{bin_key}.gbk"))
    
    return {
        "gffs": all_gffs,
        "faas": all_faas,
        "ffns": all_ffns,
        "gbks": all_gbks,
        "script": "scripts/merge_bakta.py"
    }

rule ext_annotation_bakta_merge:
    input:
        unpack(get_ext_bakta_merge_inputs)
    output:
        gff = os.path.join(EXT_RES, "{bin}/merged/{bin}.gff3"),
        faa = os.path.join(EXT_RES, "{bin}/merged/{bin}.faa"),
        ffn = os.path.join(EXT_RES, "{bin}/merged/{bin}.ffn"),
        tsv = os.path.join(EXT_RES, "{bin}/merged/{bin}.tsv"),
        gbk = os.path.join(EXT_RES, "{bin}/merged/{bin}.gbk")
    log:
        os.path.join(AN_LOGS, "external_bakta_merge/{bin}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("biopython")
    shell:
        r"""
        {{
        set -euo pipefail
        conda run -p {params.env} python3 {input.script} \
            --gffs {input.gffs} \
            --faas {input.faas} \
            --ffns {input.ffns} \
            --gbks {input.gbks} \
            --out_gff {output.gff} \
            --out_faa {output.faa} \
            --out_ffn {output.ffn} \
            --out_tsv {output.tsv} \
            --out_gbk {output.gbk}
        }} 2>&1 | tee {log}
        """

# Advanced HMM for External Genomes (Chunked Job Splitting)
rule ext_annotation_hmm_split:
    input:
        faa = get_anno_faa_input
    output:
        expand(os.path.join(EXT_RES, "{{bin}}/advanced_hmm/chunks/chunk_{sid}.faa"),
             sid=[f"{i+1}-of-{AN_C.get('hmm_split_chunks', 20)}" for i in range(AN_C.get('hmm_split_chunks', 20))])
    log:
        os.path.join(AN_LOGS, "external_hmm_split/{bin}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        outdir = lambda wildcards: os.path.join(EXT_RES, wildcards.bin, "advanced_hmm/chunks"),
        split_count = AN_C.get("hmm_split_chunks", 20),
        env = an_get_env("bakta")
    shell:
        """
        set -euo pipefail
        {{
        mkdir -p {params.outdir}/tmp
        CLEAN_FAA="{params.outdir}/cleaned_input.faa"
        sed '/^[^>]/s/[-*]//g' "{input.faa}" > "$CLEAN_FAA"
        conda run -p {params.env} seqkit split "$CLEAN_FAA" -p {params.split_count} -O {params.outdir}/tmp
        rm -f "$CLEAN_FAA"
        i=1
        for file in $(ls {params.outdir}/tmp/*.faa | sort -V); do
            mv "$file" "{params.outdir}/chunk_${{i}}-of-{params.split_count}.faa"
            ((i++))
        done
        rm -rf {params.outdir}/tmp
        }} 2>&1 | tee {log}
        """

rule ext_annotation_hmmsearch_unified:
    input:
        faa = os.path.join(EXT_RES, "{bin}/advanced_hmm/chunks/chunk_{scatteritem}.faa"),
        db = get_unified_hmm_db
    output:
        tbl = os.path.join(EXT_RES, "{bin}/advanced_hmm/unified/chunk_{scatteritem}.tbl")
    log: os.path.join(AN_LOGS, "external_hmm_unified/{bin}_{scatteritem}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5)
    shell:
        "conda run -p {params.env} hmmsearch --cpu {threads} --tblout {output.tbl} -E {params.evalue} {input.db} {input.faa} > {log} 2>&1"

rule ext_annotation_hmmsearch_pfam:
    input:
        faa = os.path.join(EXT_RES, "{bin}/advanced_hmm/chunks/chunk_{scatteritem}.faa"),
        db = get_pfam_db_input
    output:
        tbl = os.path.join(EXT_RES, "{bin}/advanced_hmm/pfam/chunk_{scatteritem}.domtbl")
    log: os.path.join(AN_LOGS, "external_hmm_pfam/{bin}_{scatteritem}.log")
    threads: get_threads(4)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5),
        run_pfam = lambda wildcards: str(config.get("annotation", {}).get("steps", {}).get("pfam", "True")).lower() == "true"
    shell:
        """
        set -euo pipefail
        if [ "{params.run_pfam}" = "true" ]; then
            conda run -p {params.env} hmmsearch --cpu {threads} --domtblout {output.tbl} -E {params.evalue} {input.db} {input.faa} > {log} 2>&1
        else
            echo "# Skipping Pfam because steps.pfam is set to False" > {output.tbl}
        fi
        """

def is_ext_hmm_split(bin_key):
    """Checks if the external genome is configured to run HMM in split mode."""
    info = EXT_GENOMES.get(bin_key)
    if isinstance(info, dict) and "split_hmm" in info:
        return str(info.get("split_hmm")).lower() == "true"
    return str(AN_C.get("ext_split_hmm", "False")).lower() == "true"

def get_ext_pfam_chunks(wildcards):
    """Dynamically return Pfam chunk files for external genome if steps.pfam is True, else empty list."""
    steps = config.get("annotation", {}).get("steps", {})
    do_pfam = str(steps.get("pfam", "True")).lower() == "true"
    bin_key = wildcards.bin
    if do_pfam:
        return expand(os.path.join(EXT_RES, f"{bin_key}/advanced_hmm/pfam/chunk_{{sid}}.domtbl"), 
                      sid=[f"{i+1}-of-{AN_C.get('hmm_split_chunks', 20)}" for i in range(AN_C.get('hmm_split_chunks', 20))])
    else:
        return []

def get_ext_hmm_gather_inputs(wildcards):
    """Determines whether to gather HMM search chunks or run unified HMM search on full FAA."""
    bin_key = wildcards.bin
    if is_ext_hmm_split(bin_key):
        return {
            "unified": expand(os.path.join(EXT_RES, f"{bin_key}/advanced_hmm/unified/chunk_{{sid}}.tbl"), 
                           sid=[f"{i+1}-of-{AN_C.get('hmm_split_chunks', 20)}" for i in range(AN_C.get('hmm_split_chunks', 20))]),
            "pfam": get_ext_pfam_chunks(wildcards),
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

rule ext_annotation_hmm_gather:
    input:
        unpack(get_ext_hmm_gather_inputs)
    output:
        unified_raw = os.path.join(EXT_RES, "{bin}/advanced_hmm/unified_all_raw.tbl"),
        pfam_raw = os.path.join(EXT_RES, "{bin}/advanced_hmm/pfam_all_raw.domtbl")
    log:
        os.path.join(AN_LOGS, "external_hmm_gather/{bin}.log")
    threads: lambda wildcards: 1 if is_ext_hmm_split(wildcards.bin) else get_threads(16)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        is_split = lambda wildcards: "true" if is_ext_hmm_split(wildcards.bin) else "false",
        env = an_get_env("bakta"),
        evalue = AN_C.get("hmm_evalue_threshold", 1e-5),
        run_pfam = lambda wildcards: str(config.get("annotation", {}).get("steps", {}).get("pfam", "True")).lower() == "true"
    shell:
        """
        set -euo pipefail
        {{
        if [ "{params.is_split}" = "true" ]; then
            echo "# hmmsearch" > {output.unified_raw}
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

use rule annotation_hmm_resolve as ext_annotation_hmm_resolve with:
    input:
        unified = os.path.join(EXT_RES, "{bin}/advanced_hmm/unified_all_raw.tbl"),
        pfam = os.path.join(EXT_RES, "{bin}/advanced_hmm/pfam_all_raw.domtbl"),
        ko_list = AN_C.get("ko_list", ""),
        faa = get_anno_faa_input
    output:
        final_table = os.path.join(EXT_RES, "{bin}/{bin}_advanced_hmm_resolved.tsv")
    log: os.path.join(AN_LOGS, "external_hmm_resolve/{bin}.log")

def get_ext_aggregate_inputs(wildcards):
    steps = config.get("annotation", {}).get("steps", {})
    inputs = {
        "faa": os.path.join(EXT_RES, f"{wildcards.bin}/{wildcards.bin}.faa"),
        "gff": os.path.join(EXT_RES, f"{wildcards.bin}/{wildcards.bin}.gff3"),
        "script": "scripts/aggregate_functional.py"
    }
    if str(steps.get("eggnog", "False")).lower() == "true":
        inputs["eggnog"] = os.path.join(EXT_RES, f"{wildcards.bin}/eggnog/{wildcards.bin}.emapper.annotations")
    if str(steps.get("kofamscan", "False")).lower() == "true":
        inputs["kofam"] = os.path.join(EXT_RES, f"{wildcards.bin}/kofamscan/{wildcards.bin}_kofam.tsv")
    if str(steps.get("advanced_hmm", "False")).lower() == "true":
        inputs["advanced_hmm"] = os.path.join(EXT_RES, f"{wildcards.bin}/{wildcards.bin}_advanced_hmm_resolved.tsv")
    return inputs

rule ext_annotation_aggregate:
    input:
        unpack(get_ext_aggregate_inputs)
    output:
        summary = os.path.join(EXT_RES, "{bin}/{bin}_functional_summary.tsv")
    log:
        os.path.join(AN_LOGS, "external_aggregate/{bin}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("biopython"),
        egg_in = lambda w, input: input.eggnog if hasattr(input, "eggnog") else "none",
        ko_in = lambda w, input: input.kofam if hasattr(input, "kofam") else "none",
        adv_hmm = lambda w, input: input.advanced_hmm if hasattr(input, "advanced_hmm") else "none",
        inherent_name = lambda wildcards: get_ext_inherent_name(wildcards.bin)
    shell:
        r"""
        set -euo pipefail
        conda run -p {params.env} python3 {input.script} \
            {input.faa} {params.egg_in} {params.ko_in} {params.adv_hmm} \
            {input.gff} {params.inherent_name} {wildcards.bin} {output.summary} 2>&1 | tee {log}
        """

rule ext_annotation_cat_unified_hmm:
    input:
        hmms = get_cat_hmm_inputs
    output:
        db = os.path.join(EXT_RES, "{bin}/advanced_hmm/combined_unified.hmm")
    log:
        os.path.join(AN_LOGS, "ext_cat_unified_hmm/{bin}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        tool = lambda wildcards: get_hmm_tool_for_genome(wildcards.bin)
    shell:
        r"""
        set -euo pipefail
        {{
        cat {input.hmms} > {output.db}
        if [ "{params.tool}" = "hmmscan" ]; then
            conda run -p {params.env} hmmpress -f {output.db}
        fi
        }} 2>&1 | tee {log}
        """
