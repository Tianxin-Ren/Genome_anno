# Bakta Annotation Rules
# ------------------------------------------------------------------------------



rule annotation_bakta_split:
    input:
        fa = lambda wildcards: an_get_asm_path(wildcards.genome)
    output:
        chunks = scatter.an_chunks(os.path.join(AN_RES, "{{genome}}/bakta/chunks/chunk_{scatteritem}.fasta")),
        meta = os.path.join(AN_RES, "{genome}/bakta/chunks/chunk_metadata.tsv")
    log:
        os.path.join(AN_LOGS, "bakta_split/{genome}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        outdir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "bakta/chunks"),
        num = AN_SPLIT_COUNT,
        env = an_get_env("seqkit"),
        min_size = AN_C.get("min_circular_size", 800000),
        config_complete = lambda wildcards: ",".join(get_complete_contigs(wildcards.genome)),
        config_excluded = lambda wildcards: ",".join(get_non_complete_contigs(wildcards.genome))
    shell:
        r"""
        {{
        set -euo pipefail
        mkdir -p {params.outdir}
        echo -e "chunk\tis_complete" > {output.meta}
        # 1. Detect Auto-Complete Contigs (circular=yes and length > 800kb)
        AUTO_COMPLETE=""
        if [ -f "{input.fa}" ]; then
            if grep -E -q "circular=yes|circular-yes" "{input.fa}"; then
                AUTO_COMPLETE=$(grep ">" "{input.fa}" | grep -E "circular=yes|circular-yes" | while read -r line; do
                    ID=$(echo "$line" | cut -d' ' -f1 | sed 's/>//')
                    LEN=$(echo "$line" | sed -n 's/.*len[-=]\([0-9]*\).*/\1/p')
                    if [ -n "$LEN" ] && [ "$LEN" -gt {params.min_size} ]; then
                        echo "$ID"
                    fi
                done | tr '\n' ',' | sed 's/,$//' || true)
            fi
        fi
        # 2. Merge with config and handle overrides
        FINAL_COMPLETE_FILE="{params.outdir}/final_complete_list.txt"
        touch "$FINAL_COMPLETE_FILE"
        # Add Auto
        IFS=',' read -ra ADDR_AUTO <<< "$AUTO_COMPLETE"
        for c in "${{ADDR_AUTO[@]}}"; do [ -n "$c" ] && echo "$c" >> "$FINAL_COMPLETE_FILE" || true; done
        # Add Config_True
        IFS=',' read -ra ADDR_TRUE <<< "{params.config_complete}"
        for c in "${{ADDR_TRUE[@]}}"; do [ -n "$c" ] && echo "$c" >> "$FINAL_COMPLETE_FILE" || true; done
        # Remove Config_False (Excluded)
        TEMP_LIST="{params.outdir}/temp_list.txt"
        IFS=',' read -ra ADDR_FALSE <<< "{params.config_excluded}"
        cp "$FINAL_COMPLETE_FILE" "$TEMP_LIST"
        > "$FINAL_COMPLETE_FILE"
        while read -r c; do
            EXCLUDE=0
            for ex in "${{ADDR_FALSE[@]}}"; do
                if [ "$c" == "$ex" ] || [[ "$c" == "$ex"* ]]; then EXCLUDE=1; break; fi
            done
            if [ $EXCLUDE -eq 0 ]; then echo "$c" >> "$FINAL_COMPLETE_FILE"; fi
        done < "$TEMP_LIST"
        rm -f "$TEMP_LIST"
        # Deduplicate
        sort -u "$FINAL_COMPLETE_FILE" -o "$FINAL_COMPLETE_FILE"
        # 3. Extract remaining and split to fill exactly {params.num} chunks
        REMAINING="{params.outdir}/remaining_assembly.fa"
        if [ -s "$FINAL_COMPLETE_FILE" ]; then
            conda run -p {params.env} seqkit grep -v -f "$FINAL_COMPLETE_FILE" {input.fa} > "$REMAINING"
        else
            cp {input.fa} "$REMAINING"
        fi
        # Split the REMAINDER into exactly {params.num} parts
        conda run -p {params.env} seqkit split2 --by-part {params.num} "$REMAINING" -O {params.outdir}/tmp_remain --force
        for i in $(seq 1 {params.num}); do
            TARGET_FILE="{params.outdir}/chunk_$i-of-{params.num}.fasta"
            SUFFIX=$(printf "%03d" $i)
            SRC_FILE=$(ls {params.outdir}/tmp_remain/*.part_$SUFFIX.* 2>/dev/null || true)
            if [ -n "$SRC_FILE" ]; then
                mv "$SRC_FILE" "$TARGET_FILE"
            else
                touch "$TARGET_FILE"
            fi
            echo -e "chunk_$i-of-{params.num}\tfalse" >> {output.meta}
        done
        rm -rf {params.outdir}/tmp_remain
        rm -f "$REMAINING" "$FINAL_COMPLETE_FILE"
        }} 2>&1 | tee {log}
        
        """
rule annotation_bakta_chunk:
    input:
        fa = os.path.join(AN_RES, "{genome}/bakta/chunks/chunk_{scatteritem}.fasta"),
        db = AN_C.get("bakta_db", ""),
        meta = os.path.join(AN_RES, "{genome}/bakta/chunks/chunk_metadata.tsv")
    output:
        gff = os.path.join(AN_RES, "{genome}/bakta/chunks/bakta_{scatteritem}/{genome}.gff3"),
        faa = os.path.join(AN_RES, "{genome}/bakta/chunks/bakta_{scatteritem}/{genome}.faa"),
        ffn = os.path.join(AN_RES, "{genome}/bakta/chunks/bakta_{scatteritem}/{genome}.ffn"),
        tsv = os.path.join(AN_RES, "{genome}/bakta/chunks/bakta_{scatteritem}/{genome}.tsv"),
        gbk = os.path.join(AN_RES, "{genome}/bakta/chunks/bakta_{scatteritem}/{genome}.gbk")
    log:
        os.path.join(AN_LOGS, "bakta_chunk/{genome}_{scatteritem}.log")
    threads: get_threads(16)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("bakta"),
        db = an_full_path(AN_C.get("bakta_db", "")),
        outdir = lambda wildcards, output: os.path.dirname(output.gff),
        locus_tag = lambda wildcards: f"c{wildcards.scatteritem.split('-of-')[0]}"
    shell:
        r"""
        {{
        set -euo pipefail
        # Check if this chunk is marked as complete in metadata
        IS_COMPLETE=$(grep "chunk_{wildcards.scatteritem}" {input.meta} | cut -f2 || echo "false")
        EXTRA_FLAGS=""
        if [ "$IS_COMPLETE" == "true" ]; then
            EXTRA_FLAGS="--complete"
        fi
        mkdir -p "{params.outdir}"
        conda run -p {params.env} bakta \
            --db {params.db} --threads {threads} --output "{params.outdir}" \
            --prefix "{wildcards.genome}" --locus-tag {params.locus_tag} \
            --keep-contig-headers $EXTRA_FLAGS --force {input.fa} 2>&1 | tee {log}
        cp "{params.outdir}/{wildcards.genome}.gbff" "{output.gbk}"
        }} 2>&1 | tee {log}
        
        """
def get_bakta_merge_inputs(wildcards):
    """Gather ALL Bakta outputs: chunks + circular refined."""
    genome = wildcards.genome
    all_gffs, all_faas, all_ffns, all_gbks = [], [], [], []

    # 1. Chunk files (replicate gather.an_chunks)
    for i in range(1, AN_SPLIT_COUNT + 1):
        sid = f"{i}-of-{AN_SPLIT_COUNT}"
        base = os.path.join(AN_RES, genome, "bakta", "chunks", f"bakta_{sid}")
        all_gffs.append(os.path.join(base, f"{genome}.gff3"))
        all_faas.append(os.path.join(base, f"{genome}.faa"))
        all_ffns.append(os.path.join(base, f"{genome}.ffn"))
        all_gbks.append(os.path.join(base, f"{genome}.gbk"))

    # 2. Circular refined files (checkpoint-dependent)
    # Auto-disabled for standalone external assemblies (which skip binning)
    steps = AN_C.get("steps", {})
    do_circular = str(steps.get("circular", "False")).lower() == "true"
    is_external_asm = "external_assembly" in AN_SAMPLES.get(genome, {})
    data_type = config.get("annotation", {}).get("data_type", "long")
    if do_circular and not is_external_asm and data_type != "short":
        try:
            checkpoints.annotation_circular_filter.get(genome=genome)
        except:
            pass
        cp_dir = os.path.join(AN_RES, genome, "circular_bins_filtered")
        if os.path.exists(cp_dir):
            for fa in sorted(os.listdir(cp_dir)):
                if fa.endswith(".fa"):
                    b = fa.replace(".fa", "")
                    all_gffs.append(os.path.join(AN_RES, genome, "circular_refined", b, f"{b}.gff3"))
                    all_faas.append(os.path.join(AN_RES, genome, "circular_refined", b, f"{b}.faa"))
                    all_ffns.append(os.path.join(AN_RES, genome, "circular_refined", b, f"{b}.ffn"))
                    all_gbks.append(os.path.join(AN_RES, genome, "circular_refined", b, f"{b}.gbk"))

    return {
        "gffs": all_gffs,
        "faas": all_faas,
        "ffns": all_ffns,
        "gbks": all_gbks,
        "script": "scripts/merge_bakta.py"
    }

rule annotation_bakta_merge:
    input:
        unpack(get_bakta_merge_inputs)
    output:
        gff = os.path.join(AN_RES, "{genome}/bakta/{genome}.gff3"),
        faa = os.path.join(AN_RES, "{genome}/bakta/{genome}.faa"),
        ffn = os.path.join(AN_RES, "{genome}/bakta/{genome}.ffn"),
        tsv = os.path.join(AN_RES, "{genome}/bakta/{genome}.tsv"),
        gbk = os.path.join(AN_RES, "{genome}/bakta/{genome}.gbk")
    log:
        os.path.join(AN_LOGS, "bakta_merge/{genome}.log")
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
def get_cdhit_partition_inputs(wildcards):
    genome = wildcards.genome
    inputs = {
        "assembly_faa": os.path.join(AN_RES, f"{genome}/bakta/{genome}.faa")
    }
    
    # Trigger the checkpoint (Standard Snakemake requirement)
    try:
        checkpoint_dir = checkpoints.annotation_circular_filter.get(genome=genome).output.dir
    except Exception:
        checkpoint_dir = os.path.join(AN_RES, genome, "circular_bins_filtered")
        
    selected = []
    if os.path.exists(checkpoint_dir):
        selected = [f.replace(".fa", "") for f in os.listdir(checkpoint_dir) if f.endswith(".fa")]
    else:
        # Fallback for dry-run
        data_type = config.get("annotation", {}).get("data_type", "long")
        sample = get_sample_name_from_genome(genome)
        if data_type == "short":
            comp = config.get("binning", {}).get("completeness", 50)
            cont = config.get("binning", {}).get("contamination", 10)
            src_dir = os.path.join(AN_BN_RES, sample, f"refinement/metawrap_{comp}_{cont}_bins")
            if os.path.exists(src_dir):
                selected = [f.replace(".fa", "") for f in os.listdir(src_dir) if f.endswith(".fa")]
                
    inputs["bin_faas"] = [os.path.join(AN_RES, f"{genome}/circular_refined/{b}/{b}.faa") for b in selected]
    return inputs

rule annotation_cdhit_partition:
    input:
        unpack(get_cdhit_partition_inputs)
    output:
        binned = os.path.join(AN_RES, "{genome}/cdhit/all_bins.faa"),
        unbinned = os.path.join(AN_RES, "{genome}/cdhit/all_unbinned.faa")
    log:
        os.path.join(AN_LOGS, "cdhit_partition/{genome}.log")
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
        conda run -p {params.env} python3 scripts/extract_binned_unbinned.py \
            --assembly-faa {input.assembly_faa} \
            --bin-faas {input.bin_faas} \
            --out-binned {output.binned} \
            --out-unbinned {output.unbinned}
        }} 2>&1 | tee {log}
        
        """
rule annotation_cdhit_binned:
    input:
        binned = os.path.join(AN_RES, "{genome}/cdhit/all_bins.faa")
    output:
        faa = os.path.join(AN_RES, "{genome}/cdhit/all_bins_c{c}.faa"),
        clstr = os.path.join(AN_RES, "{genome}/cdhit/all_bins_c{c}.faa.clstr")
    log:
        os.path.join(AN_LOGS, "cdhit_binned/{genome}_c{c}.log")
    threads: get_threads(12)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("cdhit"),
        c = lambda wildcards: wildcards.c,
        d = 400,
        n = lambda wildcards: get_cdhit_n_helper(wildcards.c)
    shell:
        r"""
        {{
        set -euo pipefail
        if [ ! -s "{input.binned}" ]; then
            echo "No binned sequences found. Creating empty outputs."
            touch {output.faa} {output.clstr}
        else
            conda run -p {params.env} cd-hit \
                -i {input.binned} \
                -o {output.faa} \
                -c {params.c} \
                -d {params.d} \
                -T {threads} \
                -M 20000 \
                -n {params.n}
        fi
        }} 2>&1 | tee {log}
        
        """
rule annotation_cdhit_unbinned:
    input:
        unbinned = os.path.join(AN_RES, "{genome}/cdhit/all_unbinned.faa")
    output:
        faa = os.path.join(AN_RES, "{genome}/cdhit/all_unbinned_c{c}.faa"),
        clstr = os.path.join(AN_RES, "{genome}/cdhit/all_unbinned_c{c}.faa.clstr")
    log:
        os.path.join(AN_LOGS, "cdhit_unbinned/{genome}_c{c}.log")
    threads: get_threads(12)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("cdhit"),
        c = lambda wildcards: wildcards.c,
        d = 400,
        n = lambda wildcards: get_cdhit_n_helper(wildcards.c)
    shell:
        r"""
        {{
        set -euo pipefail
        if [ ! -s "{input.unbinned}" ]; then
            echo "No unbinned sequences found. Creating empty outputs."
            touch {output.faa} {output.clstr}
        else
            conda run -p {params.env} cd-hit \
                -i {input.unbinned} \
                -o {output.faa} \
                -c {params.c} \
                -d {params.d} \
                -T {threads} \
                -M 20000 \
                -n {params.n}
        fi
        }} 2>&1 | tee {log}
        
        """
rule annotation_cdhit_2d:
    input:
        binned = os.path.join(AN_RES, "{genome}/cdhit/all_bins_c{c}.faa"),
        unbinned = os.path.join(AN_RES, "{genome}/cdhit/all_unbinned_c{c}.faa")
    output:
        faa = os.path.join(AN_RES, "{genome}/cdhit/unbinned_minus_binned_c{c}.faa")
    log:
        os.path.join(AN_LOGS, "cdhit_2d/{genome}_c{c}.log")
    threads: get_threads(12)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("cdhit"),
        c = lambda wildcards: wildcards.c,
        n = lambda wildcards: get_cdhit_n_helper(wildcards.c)
    shell:
        r"""
        {{
        set -euo pipefail
        if [ ! -s "{input.binned}" ]; then
            echo "Binned clustered file is empty. Copying unbinned clustered directly."
            cp {input.unbinned} {output.faa}
        elif [ ! -s "{input.unbinned}" ]; then
            echo "Unbinned clustered file is empty. Creating empty output."
            touch {output.faa}
        else
            conda run -p {params.env} cd-hit-2d \
                -i {input.binned} \
                -i2 {input.unbinned} \
                -o {output.faa} \
                -c {params.c} \
                -n {params.n} \
                -S2 30 \
                -M 0 \
                -T {threads}
        fi
        }} 2>&1 | tee {log}
        
        """
rule annotation_cdhit_final:
    input:
        binned = os.path.join(AN_RES, "{genome}/cdhit/all_bins_c{c}.faa"),
        unbinned_minus_binned = os.path.join(AN_RES, "{genome}/cdhit/unbinned_minus_binned_c{c}.faa")
    output:
        db = os.path.join(AN_RES, "{genome}/cdhit/{genome}_proteomics_database_c{c}.faa")
    log:
        os.path.join(AN_LOGS, "cdhit_final/{genome}_c{c}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        {{
        set -euo pipefail
        cat {input.binned} {input.unbinned_minus_binned} > {output.db}
        }} 2>&1 | tee {log}
        
        """
