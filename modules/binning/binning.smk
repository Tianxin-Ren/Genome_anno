# Binning Module: MAG Generation & Circular Contig Extraction
import os
import re

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
BN_ROOT = config.get("project_root", ".")
BN_ENV_DIR = config.get("conda_env_dir", "")
BN_C = config.get("binning", {})

def bn_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(BN_ROOT, rel_path))

def bn_get_env(name):
    env_name = config.get("environments", {}).get(name) or BN_C.get("conda_envs", {}).get(name, name)
    if not BN_ENV_DIR or os.path.isabs(env_name): return env_name
    return os.path.join(BN_ENV_DIR, env_name)

BN_RES = bn_full_path(BN_C.get("results_dir", "results/binning"))
BN_LOGS = bn_full_path(BN_C.get("logs_dir", os.path.join("logs", "binning")))
BN_BATCH_RES = os.path.join(BN_RES, "batches", config["batch_name"]) if config.get("batch_name") else BN_RES
BN_BATCH_LOGS = os.path.join(BN_LOGS, "batches", config["batch_name"]) if config.get("batch_name") else BN_LOGS

# Dependencies
BN_AS_RES = bn_full_path(config.get("assembly", {}).get("results_dir", "results/assembly"))
BN_PP_RES = bn_full_path(config.get("preprocess", {}).get("results_dir", "results/preprocess"))

# Sample discovery
BN_SAMPLES = config.get("samples", {})

# 2. Target Registration
# ------------------------------------------------------------------------------
steps = BN_C.get("steps", {})
for sample in BN_SAMPLES:
    info = BN_SAMPLES[sample]
    # Traditional Binning
    if "short_read" in info:
        if steps.get("metawrap"):
            t = os.path.join(BN_RES, f"{sample}/metawrap")
            TARGETS.append(t)
        if steps.get("semibin2"):
            TARGETS.append(os.path.join(BN_RES, f"{sample}/semibin2"))
        if steps.get("refinement"):
            TARGETS.append(os.path.join(BN_RES, f"{sample}/refinement/metawrap_{BN_C.get('completeness',50)}_{BN_C.get('contamination',10)}_bins.stats"))
    
    # [NEW] Circular Extraction
    if steps.get("extract_circular") and "long_read" in info:
        TARGETS.append(os.path.join(BN_RES, f"{sample}/circular_info.tsv"))
        
if steps.get("drep") and steps.get("refinement"):
    TARGETS.append(os.path.join(BN_BATCH_RES, "drep/dereplicated_genomes"))
if steps.get("gtdbtk"):
    TARGETS.append(os.path.join(BN_BATCH_RES, "gtdbtk"))

# 3. Rules
# ------------------------------------------------------------------------------

checkpoint binning_circular_filter:
    input:
        fasta = get_raw_assembly
    output:
        tsv = os.path.join(BN_RES, "{sample}/circular_info.tsv"),
        dir_flag = touch(os.path.join(BN_RES, "{sample}/circular_contigs/.done"))
    log:
        os.path.join(BN_LOGS, "circular_filter/{sample}.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        out_dir = lambda wildcards: os.path.join(BN_RES, wildcards.sample, "circular_contigs")
    run:
        os.makedirs(params.out_dir, exist_ok=True)
        header_pattern = re.compile(
             r"^>?(?P<id>[A-Za-z0-9_]+?)(?:ctg)?_len-(?P<len>[0-9]+)_circular-(?P<status>[^_]+)_depth-(?P<depth>[^_]+)_duplicated-(?P<dup>[^_ ]+)(?: mult=(?P<mult>[^ ]+))?"
        )
        sample_info = config.get("samples", {}).get(wildcards.sample, {})
        selected = set(str(k) for k in sample_info.get("selected_contigs", {}).keys())
        selected.update(str(k) for k in sample_info.get("selected_circulars", []))
        
        results = []
        current_seq = []
        current_header = None
        
        def save_contig(header, seq_lines):
            m = header_pattern.match(header)
            if m:
                info = m.groupdict()
                results.append(info)
                if info["status"].lower() in ["yes", "possibly"] or info["id"] in selected:
                    with open(os.path.join(params.out_dir, f"{info['id']}.fa"), "w") as f:
                        f.write(f">{header}\n")
                        f.write("".join(seq_lines) + "\n")
            else:
                first_tok = header.split()[0].lstrip(">")
                if first_tok in selected:
                    info = {"id": first_tok, "len": str(sum(len(l) for l in seq_lines)), "status": "yes", "depth": "NA", "dup": "NA", "mult": "1.00"}
                    results.append(info)
                    with open(os.path.join(params.out_dir, f"{first_tok}.fa"), "w") as f:
                        f.write(f">{header}\n")
                        f.write("".join(seq_lines) + "\n")
        
        with open(input.fasta, "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if current_header:
                        save_contig(current_header, current_seq)
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header:
                save_contig(current_header, current_seq)
        
        with open(output.tsv, "w") as f:
            f.write("ID\tLength\tCircular\tDepth\tDuplicated\tMultiplicity\n")
            for res in results:
                f.write(f"{res['id']}\t{res['len']}\t{res['status']}\t{res['depth']}\t{res['dup']}\t{res.get('mult', '1.00')}\n")

rule binning_mapping:
    input:
        r1 = os.path.join(BN_PP_RES, "{sample}/short_read/{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(BN_PP_RES, "{sample}/short_read/{sample}_R2.clean.fastq.gz"),
        asm = lambda w: get_filtered_assembly(w.sample)
    output:
        bam = os.path.join(BN_RES, "{sample}/mapping/{sample}.sorted.bam"),
        depth = os.path.join(BN_RES, "{sample}/mapping/{sample}.depth.txt")
    log:
        os.path.join(BN_LOGS, "mapping/{sample}.log")
    threads: get_threads(32)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = bn_get_env("mapping")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} bowtie2-build {input.asm} {input.asm}.idx > /dev/null 2>&1
        conda run -p {params.env} bowtie2 -p {threads} -x {input.asm}.idx -1 {input.r1} -2 {input.r2} | \
        conda run -p {params.env} samtools sort -@ {threads} -o {output.bam} -
        conda run -p {params.env} samtools index {output.bam}
        conda run -p {params.env} jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam} 2>&1 | tee {log}
        
        """
rule binning_metawrap:
    input:
        r1 = os.path.join(BN_PP_RES, "{sample}/short_read/{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(BN_PP_RES, "{sample}/short_read/{sample}_R2.clean.fastq.gz"),
        asm = lambda w: get_filtered_assembly(w.sample)
    output:
        directory(os.path.join(BN_RES, "{sample}/metawrap"))
    log:
        os.path.join(BN_LOGS, "metawrap/{sample}.log")
    threads: get_threads(BN_C.get("threads", {}).get("metawrap", 64))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = bn_get_env("metawrap"),
        outdir = lambda wildcards: os.path.join(BN_RES, wildcards.sample, "metawrap")
    shell:
        """
        set -euo pipefail
        TMP_READS="$(dirname "{params.outdir}")/tmp_reads_{wildcards.sample}"
        mkdir -p "$TMP_READS"
        ln -sf {input.r1} "$TMP_READS/{wildcards.sample}_1.fastq"
        ln -sf {input.r2} "$TMP_READS/{wildcards.sample}_2.fastq"
        conda run -p {params.env} metawrap binning \
            -o "{params.outdir}" \
            -t {threads} \
            -a {input.asm} \
            --metabat2 --maxbin2 --concoct \
            "$TMP_READS/{wildcards.sample}_1.fastq" "$TMP_READS/{wildcards.sample}_2.fastq" 2>&1 | tee "{log}"
        rm -rf "$TMP_READS"
        
        """
rule binning_semibin2:
    input:
        asm = lambda w: get_filtered_assembly(w.sample),
        bam = os.path.join(BN_RES, "{sample}/mapping/{sample}.sorted.bam")
    output:
        bins = directory(os.path.join(BN_RES, "{sample}/semibin2"))
    log:
        os.path.join(BN_LOGS, "semibin2/{sample}.log")
    threads: get_threads(BN_C.get("threads", {}).get("semibin2", 32))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = bn_get_env("semibin2")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} SemiBin2 single_easy_bin \
            -i {input.asm} -b {input.bam} \
            -o "{output.bins}" \
            -t {threads} --environment global 2>&1 | tee "{log}"
        
        """
rule binning_refinement:
    input:
        metawrap = os.path.join(BN_RES, "{sample}/metawrap")
    output:
        stats = os.path.join(BN_RES, "{sample}/refinement/metawrap_{comp}_{cont}_bins.stats"),
        bin_dir = directory(os.path.join(BN_RES, "{sample}/refinement/metawrap_{comp}_{cont}_bins"))
    log:
        os.path.join(BN_LOGS, "refinement/{sample}_{comp}_{cont}.log")
    threads: get_threads(BN_C.get("threads", {}).get("refinement", 64))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = bn_get_env("metawrap"),
        checkm_db = bn_full_path(BN_C.get("checkm_db", "") or config.get("databases", {}).get("checkm", ""))
    shell:
        """
        set -euo pipefail
        if [ -n "{params.checkm_db}" ] && [ -d "{params.checkm_db}" ]; then
            echo "{params.checkm_db}" | conda run -p {params.env} checkm data setRoot "{params.checkm_db}"
        fi
        conda run -p {params.env} metawrap bin_refinement \
            -o "$(dirname "{output.stats}")" \
            -t {threads} \
            -A "{input.metawrap}/metabat2_bins" \
            -B "{input.metawrap}/maxbin2_bins" \
            -C "{input.metawrap}/concoct_bins" \
            -c {wildcards.comp} -x {wildcards.cont} 2>&1 | tee "{log}"
        
        """
rule binning_gather_all:
    input:
        bin_dirs = expand(os.path.join(BN_RES, "{sample}/refinement/metawrap_{comp}_{cont}_bins"), 
                          sample=config.get("samples", {}).keys(), comp=BN_C.get("completeness", 50), cont=BN_C.get("contamination", 10))
    output:
        all_bins_dir = directory(os.path.join(BN_BATCH_RES, "all_bins_combined"))
    log:
        os.path.join(BN_BATCH_LOGS, "binning_gather_all.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        """
        set -euo pipefail
        mkdir -p "{output.all_bins_dir}"
        for d in {input.bin_dirs}; do
            sample="$(basename "$(dirname "$(dirname "$d")")")"
            for f in "$d"/*.fa; do
                [ ! -f "$f" ] && continue
                fname="$(basename "$f")"
                [[ "$fname" == "${{sample}}_"* ]] || fname="${{sample}}_$fname"
                [ ! -e "{output.all_bins_dir}/$fname" ] || {{ echo "Duplicate bin: $fname" >&2; exit 1; }}
                cp "$f" "{output.all_bins_dir}/$fname"
            done
        done 2>&1 | tee "{log}"
        """
rule binning_drep:
    input:
        bins = os.path.join(BN_BATCH_RES, "all_bins_combined")
    output:
        derep_dir = directory(os.path.join(BN_BATCH_RES, "drep/dereplicated_genomes")),
        summary = os.path.join(BN_BATCH_RES, "drep/data_tables/Widb.csv")
    log:
        os.path.join(BN_BATCH_LOGS, "drep.log")
    threads: get_threads(BN_C.get("threads", {}).get("drep", 64))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = bn_get_env("drep"),
        ani = BN_C.get("ani_threshold", 0.95),
        comp = BN_C.get("completeness", 50),
        cont = BN_C.get("contamination", 10)
    shell:
        """

        set -euo pipefail
        rm -rf "{output.derep_dir}"
        conda run -p {params.env} dRep dereplicate "$(dirname "{output.derep_dir}")" \
            -g "{input.bins}"/*.fa \
            -p {threads} \
            -sa {params.ani} \
            -comp {params.comp} \
            -con {params.cont} 2>&1 | tee "{log}"
        
        """


rule binning_gtdbtk:
    input:
        genomes = os.path.join(BN_BATCH_RES, "drep/dereplicated_genomes")
    output:
        directory(os.path.join(BN_BATCH_RES, "gtdbtk"))
    log:
        os.path.join(BN_BATCH_LOGS, "gtdbtk.log")
    threads: get_threads(BN_C.get("threads", {}).get("gtdbtk", 64))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = bn_get_env("gtdbtk"),
        db = bn_full_path(BN_C.get("gtdbtk_db", "") or config.get("databases", {}).get("gtdbtk", ""))
    shell:
        """
        set -euo pipefail
        export GTDBTK_DATA_PATH={params.db:q}
        conda run -p {params.env:q} gtdbtk classify_wf \
            --genome_dir {input.genomes:q} --extension fa \
            --out_dir {output:q} --cpus {threads} 2>&1 | tee {log:q}
        """
