# Abundance Module: Genomics & Transcriptomics Quantification
import os

# 1. Global Path Logic (Namespaced)
# ------------------------------------------------------------------------------
AB_ROOT = config.get("project_root", ".")
AB_ENV_DIR = config.get("conda_env_dir", "")
AB_C = config.get("abundance", {})

def ab_full_path(rel_path):
    if not rel_path: return ""
    if os.path.isabs(rel_path): return rel_path
    return os.path.normpath(os.path.join(AB_ROOT, rel_path))

def ab_get_env(name):
    env_name = config.get("environments", {}).get(name) or AB_C.get("conda_envs", {}).get(name, name)
    if not AB_ENV_DIR: return env_name
    return os.path.join(AB_ENV_DIR, env_name)

AB_RES = ab_full_path(AB_C.get("results_dir", "results/abundance"))
AB_LOGS = ab_full_path(AB_C.get("logs_dir", os.path.join("logs", "abundance")))
AB_COVERM_RES = os.path.join(AB_RES, "batches", config["batch_name"]) if config.get("batch_name") else AB_RES
AB_COVERM_LOGS = os.path.join(AB_LOGS, "batches", config["batch_name"]) if config.get("batch_name") else AB_LOGS
AB_COVERM_MERGED = os.path.join(AB_COVERM_RES, f"{config['batch_name']}_coverm_abundance.tsv" if config.get("batch_name") else "merged_coverm_abundance.tsv")
AB_BIN_RES = ab_full_path(config.get("binning", {}).get("results_dir", "results/binning"))
AB_BIN_BATCH_RES = os.path.join(AB_BIN_RES, "batches", config["batch_name"]) if config.get("batch_name") else AB_BIN_RES
AB_COVERM_REF = os.path.join(AB_BIN_BATCH_RES, "drep/dereplicated_genomes") if AB_C.get("coverm_use_drep") else ab_full_path(AB_C.get("coverm_ref", ""))

# Sample discovery
AB_SAMPLES = config.get("samples", {})
AB_DNA_SAMPLES = [s for s, info in AB_SAMPLES.items() if "short_read" in info]
AB_RNA_SAMPLES = [s for s, info in AB_SAMPLES.items() if "transcriptome" in info]
AB_RATIO_SAMPLES = [s for s, info in AB_SAMPLES.items() if "short_read" in info and "transcriptome" in info]

# Reference Genome Resolution Helpers
def get_ref_name_helper(sample):
    info = AB_SAMPLES.get(sample, {})
    ref_path = info.get("abundance_ref") or info.get("external_assembly")
    if not ref_path:
        return "default_ref"
    basename = os.path.basename(os.path.normpath(ref_path))
    for ext in [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"]:
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
            break
    # Add a prefix "ref_" to distinguish from regular sample folders in results/annotation
    return f"ref_{basename}"

def get_ref_dir_by_name(ref_name):
    # Strip "ref_" prefix if present
    pure_name = ref_name[4:] if ref_name.startswith("ref_") else ref_name
    for sample, info in AB_SAMPLES.items():
        ref_path = info.get("abundance_ref") or info.get("external_assembly")
        if ref_path:
            resolved_path = ab_full_path(ref_path)
            name = os.path.basename(os.path.normpath(resolved_path))
            for ext in [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"]:
                if name.endswith(ext):
                    name = name[:-len(ext)]
                    break
            if name == pure_name:
                return resolved_path
    raise ValueError(f"Could not resolve reference source path for name: {ref_name}")

def get_bowtie2_index_prefix(sample):
    info = AB_SAMPLES.get(sample, {})
    if "abundance_ref" in info or "external_assembly" in info:
        ref_name = get_ref_name_helper(sample)
        return os.path.join(AB_RES, f"references/{ref_name}/bowtie2_index/genome")
    else:
        return os.path.join(AB_RES, f"{sample}/reference/bowtie2_index/genome")

def get_bowtie2_index_file(sample):
    return get_bowtie2_index_prefix(sample) + ".1.bt2"

def get_hisat2_index_prefix(sample):
    info = AB_SAMPLES.get(sample, {})
    if "abundance_ref" in info or "external_assembly" in info:
        ref_name = get_ref_name_helper(sample)
        return os.path.join(AB_RES, f"references/{ref_name}/hisat2_index/genome")
    else:
        return os.path.join(AB_RES, f"{sample}/reference/hisat2_index/genome")

def get_hisat2_index_file(sample):
    return get_hisat2_index_prefix(sample) + ".1.ht2"

def get_gff3_path(sample):
    info = AB_SAMPLES.get(sample, {})
    if "abundance_ref" in info or "external_assembly" in info:
        ref_name = get_ref_name_helper(sample)
        return os.path.join(AB_RES, f"{ref_name}.gff")
    else:
        annotation_dir = config.get("annotation", {}).get("results_dir", "results/annotation")
        return ab_full_path(os.path.join(annotation_dir, sample, "bakta", f"{sample}.gff3"))

def get_mag_annotation_tsv():
    cfg_tsv = AB_C.get("mag_annotation", {}).get("annotation_tsv", "")
    if cfg_tsv:
        resolved = ab_full_path(cfg_tsv)
        if os.path.exists(resolved):
            return resolved
    # Heuristic fallback discovery for primary/secondary mag_annotation assets
    candidates = [
        ab_full_path("protein_KO_all_strict.tsv"),
        ab_full_path("protein_KO_best.tsv"),
        ab_full_path("protein_annotation_master.tsv"),
        os.path.join(WORKFLOW_ROOT, "protein_KO_all_strict.tsv")
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return ""

# Mapping and quantification read resolvers (Bypasses fastp/sortmerna if clean/mRNA inputs are provided)
def get_dna_mapping_r1(sample):
    sr = AB_SAMPLES.get(sample, {}).get("short_read", {})
    if "external_clean_r1" in sr:
        return ab_full_path(sr["external_clean_r1"])
    if AB_C.get("use_preprocessed_reads", False):
        return ab_full_path(os.path.join(config.get("preprocess", {}).get("results_dir", "results/preprocess"), f"{sample}/short_read/{sample}_R1.clean.fastq.gz"))
    return os.path.join(AB_RES, f"{sample}/preprocess/genomics/{sample}_R1.clean.fastq.gz")

def get_dna_mapping_r2(sample):
    sr = AB_SAMPLES.get(sample, {}).get("short_read", {})
    if "external_clean_r2" in sr:
        return ab_full_path(sr["external_clean_r2"])
    if AB_C.get("use_preprocessed_reads", False):
        return ab_full_path(os.path.join(config.get("preprocess", {}).get("results_dir", "results/preprocess"), f"{sample}/short_read/{sample}_R2.clean.fastq.gz"))
    return os.path.join(AB_RES, f"{sample}/preprocess/genomics/{sample}_R2.clean.fastq.gz")

def get_rna_mapping_r1(sample):
    info = AB_SAMPLES.get(sample, {})
    tr = info.get("transcriptome", {})
    if "mRNA_1" in tr or "external_mrna_r1" in tr or "mrna_r1" in tr:
        path = tr.get("mRNA_1") or tr.get("external_mrna_r1") or tr.get("mrna_r1")
        return ab_full_path(path)
    return os.path.join(AB_RES, f"{sample}/preprocess/transcriptomics/{sample}.mRNA_1.fastq.gz")

def get_rna_mapping_r2(sample):
    info = AB_SAMPLES.get(sample, {})
    tr = info.get("transcriptome", {})
    if "mRNA_2" in tr or "external_mrna_r2" in tr or "mrna_r2" in tr:
        path = tr.get("mRNA_2") or tr.get("external_mrna_r2") or tr.get("mrna_r2")
        return ab_full_path(path)
    return os.path.join(AB_RES, f"{sample}/preprocess/transcriptomics/{sample}.mRNA_2.fastq.gz")

# 2. Target Registration
# ------------------------------------------------------------------------------
AB_TARGETS = []

RUN_COVERM = AB_C.get("run_coverm", False)
RUN_GENE_QUANT = AB_C.get("run_gene_quant", True)
RUN_SUMMARIZE_BY_GENE = AB_C.get("run_summarize_by_gene", True)

for sample in AB_SAMPLES:
    info = AB_SAMPLES[sample]
    has_dna = "short_read" in info
    has_rna = "transcriptome" in info
    has_ref = "abundance_ref" in info
    
    if has_ref or AB_COVERM_REF:
        if has_dna:
            if RUN_COVERM:
                AB_TARGETS.append(os.path.join(AB_COVERM_RES, f"{sample}/genome_abundance/{sample}.coverm.tsv"))
            if RUN_GENE_QUANT and has_ref:
                AB_TARGETS.append(os.path.join(AB_RES, f"{sample}/dna_quant/{sample}.final.tsv"))
        if has_rna and RUN_GENE_QUANT and has_ref:
            AB_TARGETS.append(os.path.join(AB_RES, f"{sample}/rna_quant/{sample}.final.tsv"))
        if has_dna and has_rna and RUN_GENE_QUANT and has_ref:
            AB_TARGETS.append(os.path.join(AB_RES, f"{sample}/ratio/{sample}.ratio.tsv"))

# Add merged target files if we have active samples
if RUN_GENE_QUANT:
    if AB_DNA_SAMPLES:
        AB_TARGETS.append(os.path.join(AB_RES, "merged_dna_tpm.tsv"))
    if AB_RNA_SAMPLES:
        AB_TARGETS.append(os.path.join(AB_RES, "merged_rna_tpm.tsv"))
    if AB_RATIO_SAMPLES:
        AB_TARGETS.append(os.path.join(AB_RES, "merged_ratio.tsv"))
    if RUN_SUMMARIZE_BY_GENE and (AB_RATIO_SAMPLES or AB_DNA_SAMPLES or AB_RNA_SAMPLES):
        AB_TARGETS.append(os.path.join(AB_RES, "gene_abundance_and_annotation.tsv"))

if RUN_COVERM and AB_DNA_SAMPLES:
    AB_TARGETS.append(AB_COVERM_MERGED)

TARGETS.extend(AB_TARGETS)


# 3. Rules
# ------------------------------------------------------------------------------

# --- 3.1 Preprocessing ---

def get_genomics_r1(wildcards):
    sr = AB_SAMPLES[wildcards.sample].get("short_read", {})
    path = sr.get("r1") or sr.get("external_clean_r1")
    return ab_full_path(path) if path else "SHORT_READ_R1_NOT_PROVIDED"

def get_genomics_r2(wildcards):
    sr = AB_SAMPLES[wildcards.sample].get("short_read", {})
    path = sr.get("r2") or sr.get("external_clean_r2")
    return ab_full_path(path) if path else "SHORT_READ_R2_NOT_PROVIDED"

rule abundance_preprocess_genomics_fastp:
    input:
        r1 = get_genomics_r1,
        r2 = get_genomics_r2
    output:
        r1 = os.path.join(AB_RES, "{sample}/preprocess/genomics/{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(AB_RES, "{sample}/preprocess/genomics/{sample}_R2.clean.fastq.gz"),
        json = os.path.join(AB_RES, "{sample}/preprocess/genomics/{sample}.fastp.json"),
        html = os.path.join(AB_RES, "{sample}/preprocess/genomics/{sample}.fastp.html")
    log: os.path.join(AB_LOGS, "fastp_genomics/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("fastp", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("fastp"),
        extra = AB_C.get("params", {}).get("fastp", "")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} --thread {threads} {params.extra} 2>&1 | tee {log}
        
        """
def get_transcriptomics_r1(wildcards):
    path = AB_SAMPLES[wildcards.sample].get("transcriptome", {}).get("r1")
    return ab_full_path(path) if path else "TRANSCRIPTOMICS_R1_NOT_PROVIDED"

def get_transcriptomics_r2(wildcards):
    path = AB_SAMPLES[wildcards.sample].get("transcriptome", {}).get("r2")
    return ab_full_path(path) if path else "TRANSCRIPTOMICS_R2_NOT_PROVIDED"

rule abundance_preprocess_transcriptomics_fastp:
    input:
        r1 = get_transcriptomics_r1,
        r2 = get_transcriptomics_r2
    output:
        r1 = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}_R1.trim.fastq.gz"),
        r2 = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}_R2.trim.fastq.gz"),
        json = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}.fastp.json"),
        html = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}.fastp.html")
    log: os.path.join(AB_LOGS, "fastp_transcriptomics/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("fastp", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("fastp"),
        extra = AB_C.get("params", {}).get("fastp", "")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} --thread {threads} {params.extra} 2>&1 | tee {log}
        
        """
rule abundance_preprocess_transcriptomics_sortmerna:
    input:
        r1 = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}_R1.trim.fastq.gz"),
        r2 = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}_R2.trim.fastq.gz")
    output:
        mRNA_r1 = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}.mRNA_1.fastq.gz"),
        mRNA_r2 = os.path.join(AB_RES, "{sample}/preprocess/transcriptomics/{sample}.mRNA_2.fastq.gz")
    log: os.path.join(AB_LOGS, "sortmerna/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("sortmerna", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("sortmerna"),
        db = AB_C.get("sortmerna_db", ""),
        workdir = lambda w: os.path.join(AB_RES, w.sample, "preprocess/transcriptomics/sortmerna_work"),
        out_prefix = lambda w: os.path.join(AB_RES, w.sample, "preprocess/transcriptomics", w.sample)
    shell:
        """
        set -euo pipefail
        rm -rf {params.workdir} && mkdir -p {params.workdir}
        conda run -p {params.env} sortmerna --threads {threads} --ref {params.db} \
            --reads {input.r1} --reads {input.r2} --workdir {params.workdir} \
            --aligned {params.out_prefix}.rRNA --other {params.out_prefix}.mRNA \
            --paired_out --fastx --out2 2>&1 | tee {log}
        # Robustly handle output naming (SortMeRNA 4.x can be inconsistent)
        # Check for _1/_2 or _fwd/_rev or .mRNA (interleaved)
        if [ -f {params.out_prefix}.mRNA_1.fastq ]; then
            gzip -c {params.out_prefix}.mRNA_1.fastq > {output.mRNA_r1}
            gzip -c {params.out_prefix}.mRNA_2.fastq > {output.mRNA_r2}
        elif [ -f {params.out_prefix}.mRNA_fwd.fastq ]; then
            gzip -c {params.out_prefix}.mRNA_fwd.fastq > {output.mRNA_r1}
            gzip -c {params.out_prefix}.mRNA_rev.fastq > {output.mRNA_r2}
        elif [ -f {params.out_prefix}.mRNA.fastq ]; then
            # If interleaved, we try to de-interleave or just use it as r1
            gzip -c {params.out_prefix}.mRNA.fastq > {output.mRNA_r1}
            touch {output.mRNA_r2}
        elif ls {params.out_prefix}.mRNA*.gz 1> /dev/null 2>&1; then
            # If already compressed, just move them
            mv $(ls {params.out_prefix}.mRNA*_1.f*gz | head -n 1) {output.mRNA_r1} || mv $(ls {params.out_prefix}.mRNA*fwd.f*gz | head -n 1) {output.mRNA_r1}
            mv $(ls {params.out_prefix}.mRNA*_2.f*gz | head -n 1) {output.mRNA_r2} || mv $(ls {params.out_prefix}.mRNA*rev.f*gz | head -n 1) {output.mRNA_r2}
        else
            echo "Error: SortMeRNA output files not found with expected patterns in $(dirname {params.out_prefix})"
            ls -lh $(dirname {params.out_prefix})
            exit 1
        fi
        # Clean up temporary files
        rm -f {params.out_prefix}.mRNA*.fastq {params.out_prefix}.rRNA*.fastq
        rm -rf {params.workdir}
        
        """
rule abundance_coverm_genome:
    input:
        r1 = lambda w: get_dna_mapping_r1(w.sample),
        r2 = lambda w: get_dna_mapping_r2(w.sample),
        ref_dir = lambda w: ab_full_path(AB_SAMPLES[w.sample].get("abundance_ref")) or AB_COVERM_REF
    output:
        tsv = os.path.join(AB_COVERM_RES, "{sample}/genome_abundance/{sample}.coverm.tsv")
    log: os.path.join(AB_COVERM_LOGS, "coverm/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("coverm", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("coverm"),
        ext = AB_C.get("mag_extension", "fa"),
        extra = AB_C.get("params", {}).get("coverm", "--min-read-percent-identity 95")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} coverm genome \
            --coupled {input.r1} {input.r2} \
            -d {input.ref_dir} -x {params.ext} \
            -p bwa-mem2 {params.extra} \
            -t {threads} -m relative_abundance \
            --output-file {output.tsv} 2>&1 | tee {log}
        
        """
rule abundance_merge_coverm:
    input:
        tsvs = lambda wildcards: expand(os.path.join(AB_COVERM_RES, "{sample}/genome_abundance/{sample}.coverm.tsv"), sample=AB_DNA_SAMPLES)
    output:
        merged = AB_COVERM_MERGED
    log:
        os.path.join(AB_COVERM_LOGS, "merge_coverm.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("python")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/merge_abundance.py" --type coverm --inputs {input.tsvs} --output {output.merged} 2>&1 | tee {log}
        
        """
rule abundance_combine_reference:
    output:
        fa = os.path.join(AB_RES, "{ref_name}.fa"),
        gff = os.path.join(AB_RES, "{ref_name}.gff")
    log:
        os.path.join(AB_LOGS, "combine_reference/{ref_name}.log")
    threads: get_threads(AB_C.get("threads", {}).get("combine_ref", 4))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("python"),
        ref_dir = lambda w: get_ref_dir_by_name(w.ref_name),
        genes_dir = ab_full_path(AB_C.get("mag_annotation", {}).get("genes_dir", "")),
        annotation_tsv = lambda w: get_mag_annotation_tsv(),
        mag_ext = AB_C.get("mag_extension", "fa")
    shell:
        """
        set -euo pipefail
        GENES_FLAG=""
        if [ -n "{params.genes_dir}" ] && [ -d "{params.genes_dir}" ]; then
            GENES_FLAG="--genes-dir {params.genes_dir}"
        fi
        ANNO_FLAG=""
        if [ -n "{params.annotation_tsv}" ] && [ -f "{params.annotation_tsv}" ]; then
            ANNO_FLAG="--annotation-tsv {params.annotation_tsv}"
        fi
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/combine_mag_reference.py" \
            --ref-dir "{params.ref_dir}" \
            --out-fa "{output.fa}" \
            --out-gff "{output.gff}" \
            --mag-ext "{params.mag_ext}" \
            $GENES_FLAG $ANNO_FLAG 2>&1 | tee {log}
        """

rule abundance_build_bt2_index_shared:
    input:
        fa = os.path.join(AB_RES, "{ref_name}.fa")
    output:
        os.path.join(AB_RES, "references/{ref_name}/bowtie2_index/genome.1.bt2")
    log:
        os.path.join(AB_LOGS, "bowtie2_build_shared/{ref_name}.log")
    threads: get_threads(AB_C.get("threads", {}).get("mapping", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("mapping"),
        prefix = lambda w: os.path.join(AB_RES, f"references/{w.ref_name}/bowtie2_index/genome")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} bowtie2-build --threads {threads} -f {input.fa} {params.prefix} 2>&1 | tee {log}
        
        """
rule abundance_build_hisat2_index_shared:
    input:
        fa = os.path.join(AB_RES, "{ref_name}.fa")
    output:
        os.path.join(AB_RES, "references/{ref_name}/hisat2_index/genome.1.ht2")
    log:
        os.path.join(AB_LOGS, "hisat2_build_shared/{ref_name}.log")
    threads: get_threads(AB_C.get("threads", {}).get("mapping", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("mapping"),
        prefix = lambda w: os.path.join(AB_RES, f"references/{w.ref_name}/hisat2_index/genome")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} hisat2-build -p {threads} {input.fa} {params.prefix} 2>&1 | tee {log}
        
        """
def get_assembly_fasta(wildcards):
    info = AB_SAMPLES[wildcards.sample]
    assembly_dir = config.get("assembly", {}).get("results_dir", "results/assembly")
    if "external_assembly" in info:
        return ab_full_path(info["external_assembly"])
    if "long_read" in info:
        return ab_full_path(os.path.join(assembly_dir, wildcards.sample, "long_read/myloasm", wildcards.sample, f"{wildcards.sample}.myloasm.fa"))
    return ab_full_path(os.path.join(assembly_dir, wildcards.sample, f"{wildcards.sample}_megahit", f"{wildcards.sample}.contigs.fa"))

rule abundance_build_bt2_index_standard:
    input:
        fa = get_assembly_fasta
    output:
        os.path.join(AB_RES, "{sample}/reference/bowtie2_index/genome.1.bt2")
    log:
        os.path.join(AB_LOGS, "bowtie2_build/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("mapping", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("mapping"),
        prefix = lambda w: os.path.join(AB_RES, w.sample, "reference/bowtie2_index/genome")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} bowtie2-build --threads {threads} -f {input.fa} {params.prefix} 2>&1 | tee {log}
        
        """
rule abundance_build_hisat2_index_standard:
    input:
        fa = get_assembly_fasta
    output:
        os.path.join(AB_RES, "{sample}/reference/hisat2_index/genome.1.ht2")
    log:
        os.path.join(AB_LOGS, "hisat2_build/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("mapping", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("mapping"),
        prefix = lambda w: os.path.join(AB_RES, w.sample, "reference/hisat2_index/genome")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} hisat2-build -p {threads} {input.fa} {params.prefix} 2>&1 | tee {log}
        
        """
rule abundance_dna_mapping:
    input:
        r1 = lambda w: get_dna_mapping_r1(w.sample),
        r2 = lambda w: get_dna_mapping_r2(w.sample),
        idx = lambda w: get_bowtie2_index_file(w.sample)
    output: bam = os.path.join(AB_RES, "{sample}/dna_quant/{sample}.sorted.bam")
    log: os.path.join(AB_LOGS, "dna_mapping/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("mapping", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("mapping"),
        prefix = lambda w: get_bowtie2_index_prefix(w.sample)
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} bash -c "
        bowtie2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} | \
        samtools view -@ {threads} -bS - | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}" 2>&1 | tee {log}
        
        """
rule abundance_rna_mapping:
    input:
        r1 = lambda w: get_rna_mapping_r1(w.sample),
        r2 = lambda w: get_rna_mapping_r2(w.sample),
        idx = lambda w: get_hisat2_index_file(w.sample)
    output: bam = os.path.join(AB_RES, "{sample}/rna_quant/{sample}.sorted.bam")
    log: os.path.join(AB_LOGS, "rna_mapping/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("mapping", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("mapping"),
        prefix = lambda w: get_hisat2_index_prefix(w.sample)
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} bash -c "
        hisat2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} --no-spliced-alignment | \
        samtools view -@ {threads} -bS - | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}" 2>&1 | tee {log}
        
        """
rule abundance_feature_counts:
    input:
        bam = os.path.join(AB_RES, "{sample}/{type}_quant/{sample}.sorted.bam"),
        gff = lambda w: get_gff3_path(w.sample)
    output:
        counts = os.path.join(AB_RES, "{sample}/{type}_quant/{sample}.featureCounts.tsv"),
        summary = os.path.join(AB_RES, "{sample}/{type}_quant/{sample}.featureCounts.tsv.summary")
    log: os.path.join(AB_LOGS, "{type}_counts/{sample}.log")
    threads: get_threads(AB_C.get("threads", {}).get("mapping", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("subread"),
        extra = AB_C.get("params", {}).get("featureCounts", "-p -t CDS -g ID -M --primary")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} featureCounts -T {threads} {params.extra} -a {input.gff} -o {output.counts} {input.bam} 2>&1 | tee {log}
        
        """
rule abundance_calculate_tpm:
    input:
        counts = os.path.join(AB_RES, "{sample}/{type}_quant/{sample}.featureCounts.tsv"),
        summary = os.path.join(AB_RES, "{sample}/{type}_quant/{sample}.featureCounts.tsv.summary")
    output: final = os.path.join(AB_RES, "{sample}/{type}_quant/{sample}.final.tsv")
    log: os.path.join(AB_LOGS, "{type}_tpm/{sample}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params: env = ab_get_env("python")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/calculate_tpm.py" --counts {input.counts} --summary {input.summary} --output {output.final} 2>&1 | tee {log}
        
        """
rule abundance_calculate_ratio:
    input:
        dna = os.path.join(AB_RES, "{sample}/dna_quant/{sample}.final.tsv"),
        rna = os.path.join(AB_RES, "{sample}/rna_quant/{sample}.final.tsv")
    output:
        ratio = os.path.join(AB_RES, "{sample}/ratio/{sample}.ratio.tsv")
    log:
        os.path.join(AB_LOGS, "ratio/{sample}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("python")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/calculate_tpm_ratio.py" --dna {input.dna} --rna {input.rna} --output {output.ratio} 2>&1 | tee {log}
        
        """
rule abundance_merge_dna:
    input:
        files = lambda wildcards: expand(os.path.join(AB_RES, "{sample}/dna_quant/{sample}.final.tsv"), sample=AB_DNA_SAMPLES)
    output:
        merged = os.path.join(AB_RES, "merged_dna_tpm.tsv")
    log:
        os.path.join(AB_LOGS, "merge_dna.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("python")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/merge_abundance.py" --type dna --inputs {input.files} --output {output.merged} 2>&1 | tee {log}
        
        """
rule abundance_merge_rna:
    input:
        files = lambda wildcards: expand(os.path.join(AB_RES, "{sample}/rna_quant/{sample}.final.tsv"), sample=AB_RNA_SAMPLES)
    output:
        merged = os.path.join(AB_RES, "merged_rna_tpm.tsv")
    log:
        os.path.join(AB_LOGS, "merge_rna.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("python")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/merge_abundance.py" --type rna --inputs {input.files} --output {output.merged} 2>&1 | tee {log}
        
        """
rule abundance_merge_ratio:
    input:
        files = lambda wildcards: expand(os.path.join(AB_RES, "{sample}/ratio/{sample}.ratio.tsv"), sample=AB_RATIO_SAMPLES)
    output:
        merged = os.path.join(AB_RES, "merged_ratio.tsv")
    log:
        os.path.join(AB_LOGS, "merge_ratio.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("python")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/merge_abundance.py" --type ratio --inputs {input.files} --output {output.merged} 2>&1 | tee {log}
        
        """

rule abundance_summarize_by_gene:
    input:
        ratio = lambda wildcards: (
            os.path.join(AB_RES, "merged_ratio.tsv") if AB_RATIO_SAMPLES else (
                os.path.join(AB_RES, "merged_dna_tpm.tsv") if AB_DNA_SAMPLES else os.path.join(AB_RES, "merged_rna_tpm.tsv")
            )
        )
    output:
        master = os.path.join(AB_RES, "gene_abundance_and_annotation.tsv")
    log:
        os.path.join(AB_LOGS, "summarize_abundance_by_gene.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = ab_get_env("python"),
        anno = lambda wildcards: get_mag_annotation_tsv(),
        target_genes = ab_full_path(AB_C.get("mag_annotation", {}).get("target_genes", "")),
        mag_types = ab_full_path(AB_C.get("mag_annotation", {}).get("mag_types", ""))
    shell:
        """
        set -euo pipefail
        TARGET_FLAG=""
        if [ -n "{params.target_genes}" ] && [ -f "{params.target_genes}" ]; then
            TARGET_FLAG="--target-genes {params.target_genes}"
        fi
        TYPE_FLAG=""
        if [ -n "{params.mag_types}" ] && [ -f "{params.mag_types}" ]; then
            TYPE_FLAG="--mag-types {params.mag_types}"
        fi
        ANNO_FLAG=""
        if [ -n "{params.anno}" ] && [ -f "{params.anno}" ]; then
            ANNO_FLAG="--annotation-tsv {params.anno}"
        fi
        conda run -p {params.env} python3 "{WORKFLOW_ROOT}/scripts/summarize_abundance_by_gene.py" \
            --merged-ratio "{input.ratio}" \
            $ANNO_FLAG \
            --output "{output.master}" \
            $TARGET_FLAG $TYPE_FLAG 2>&1 | tee {log}
        """

