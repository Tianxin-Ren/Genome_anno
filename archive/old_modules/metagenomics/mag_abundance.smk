# MAG Abundance Pipeline Module (Sub-module of Metagenomics)
import os

# 1. Setup & Path Logic
# ------------------------------------------------------------------------------
c_abund = config.get("mag_abundance", {})
MA_PROJ = os.path.abspath(c_abund.get("project_dir", "."))

def norm_ma_path(base, rel):
    if not rel: return base
    if os.path.isabs(rel): return os.path.normpath(rel)
    return os.path.normpath(os.path.join(base, rel))

MA_RES  = norm_ma_path(MA_PROJ, c_abund.get("results_dir", "results/mag_abundance"))
MA_DATA = norm_ma_path(MA_PROJ, c_abund.get("data_dir", "data"))
MA_OUT  = MA_RES 
MA_LOGS = os.path.join(MA_RES, "logs")
COMBINED_REF_DIR = os.path.join(MA_OUT, "combined_reference")

DNA_FQ_DIR = norm_ma_path(MA_DATA, c_abund.get("dna_fq_dir", "dna_fq"))
RNA_FQ_DIR = norm_ma_path(MA_DATA, c_abund.get("rna_fq_dir", "rna_fq"))

# Resolve MAG Directory
MAG_DIR = c_abund.get("mag_dir")
if not MAG_DIR:
    # Fallback/Default
    MAG_DIR = "MAGs"
elif not os.path.isabs(MAG_DIR):
    MAG_DIR = norm_ma_path(MA_DATA, MAG_DIR)

# Sample Handling
SAMPLES_MAP = config.get("samples", {})
if SAMPLES_MAP:
    DNA_SAMPLES = [s for s, info in SAMPLES_MAP.items() if "dna" in info]
    RNA_SAMPLES = [s for s, info in SAMPLES_MAP.items() if "rna" in info]
    SAMPLES_WITH_RATIO = [s for s, info in SAMPLES_MAP.items() if "dna" in info and "rna" in info]
else:
    DNA_SAMPLES = c_abund.get("dna_samples", [])
    RNA_SAMPLES = c_abund.get("rna_samples", [])
    SAMPLES_WITH_RATIO = []

# Helper functions
def get_dna_prefix(wildcards):
    return config.get("samples", {}).get(wildcards.sample, {}).get("dna", wildcards.sample)

def get_rna_prefix(wildcards):
    return config.get("samples", {}).get(wildcards.sample, {}).get("rna", wildcards.sample)

# 2. Target Registration
# ------------------------------------------------------------------------------
if DNA_SAMPLES:
    TARGETS.append(os.path.join(MA_OUT, "abundance", "merged_dna_tpm_rpkm.final.tsv"))
if RNA_SAMPLES:
    TARGETS.append(os.path.join(MA_OUT, "abundance", "merged_rna_tpm_rpkm.final.tsv"))
if SAMPLES_WITH_RATIO:
    TARGETS.extend(expand(os.path.join(MA_OUT, "abundance", "{sample}.ratio.tsv"), sample=SAMPLES_WITH_RATIO))
    if c_abund.get("summarize"):
        TARGETS.append(os.path.join(MA_OUT, "summary", "master_by_mag.tsv"))

# 3. Rules
# ------------------------------------------------------------------------------

# 1. Prepare Combined Reference
rule prepare_combined_ref:
    input:
        mag_dir = MAG_DIR
    output:
        fa = os.path.join(COMBINED_REF_DIR, "combined_mags.fa")
    log:
        os.path.join(MA_LOGS, "prepare_ref.log")
    shell:
        """
        mkdir -p {COMBINED_REF_DIR}
        rm -f {output.fa}
        # Combine all MAGs, prefixing contig names with MAG name
        find {input.mag_dir} -maxdepth 1 -name "*.fa" -type f | while read mag_fa; do
            mag_name=$(basename "$mag_fa" .fa)
            sed "s/^>/>$mag_name|/" "$mag_fa" >> {output.fa}
        done 2>&1 | tee {log}
        """

rule combined_prodigal:
    input:
        fa = os.path.join(COMBINED_REF_DIR, "combined_mags.fa")
    output:
        gff = os.path.join(COMBINED_REF_DIR, "combined_mags.protein.gff"),
        faa = os.path.join(COMBINED_REF_DIR, "combined_mags.protein.faa")
    log:
        os.path.join(MA_LOGS, "combined_prodigal.log")
    params:
        env = c_abund.get("conda_envs", {}).get("bio", "bio")
    shell:
        "conda run -p {params.env} prodigal -i {input.fa} -a {output.faa} -f gff -o {output.gff} -p single -q 2>&1 | tee {log}"

rule combined_gtf:
    input:
        gff = os.path.join(COMBINED_REF_DIR, "combined_mags.protein.gff")
    output:
        gtf = os.path.join(COMBINED_REF_DIR, "combined_mags.protein.gtf")
    log:
        os.path.join(MA_LOGS, "combined_gtf.log")
    params:
        env = c_abund.get("conda_envs", {}).get("bio", "bio")
    shell:
        "conda run -p {params.env} gffread -E {input.gff} -T -o {output.gtf} 2>&1 | tee {log}"

rule index_bt2:
    input:
        fa = os.path.join(COMBINED_REF_DIR, "combined_mags.fa")
    output:
        os.path.join(COMBINED_REF_DIR, "bowtie2_index/combined_mags.1.bt2")
    log:
        os.path.join(MA_LOGS, "index_bt2.log")
    threads: c_abund.get("threads", 8)
    params:
        env = c_abund.get("conda_envs", {}).get("mapping", "mapping")
    shell:
        "conda run -p {params.env} bowtie2-build --threads {threads} -f {input.fa} {COMBINED_REF_DIR}/bowtie2_index/combined_mags 2>&1 | tee {log}"

rule index_hisat2:
    input:
        fa = os.path.join(COMBINED_REF_DIR, "combined_mags.fa")
    output:
        os.path.join(COMBINED_REF_DIR, "hisat2_index/combined_mags.1.ht2")
    log:
        os.path.join(MA_LOGS, "index_hisat2.log")
    threads: c_abund.get("threads", 8)
    params:
        env = c_abund.get("conda_envs", {}).get("mapping", "mapping")
    shell:
        "conda run -p {params.env} hisat2-build -p {threads} {input.fa} {COMBINED_REF_DIR}/hisat2_index/combined_mags 2>&1 | tee {log}"

# 2. Mapping Abundance (DNA - Metagenomics)
rule dna_abundance:
    input:
        fq1 = lambda wildcards: os.path.join(DNA_FQ_DIR, f"{get_dna_prefix(wildcards)}_trim.R1.fq.gz"),
        fq2 = lambda wildcards: os.path.join(DNA_FQ_DIR, f"{get_dna_prefix(wildcards)}_trim.R2.fq.gz"),
        u1 = lambda wildcards: os.path.join(DNA_FQ_DIR, f"{get_dna_prefix(wildcards)}_unpair.R1.fq.gz"),
        u2 = lambda wildcards: os.path.join(DNA_FQ_DIR, f"{get_dna_prefix(wildcards)}_unpair.R2.fq.gz"),
        idx = os.path.join(COMBINED_REF_DIR, "bowtie2_index/combined_mags.1.bt2"),
        gtf = os.path.join(COMBINED_REF_DIR, "combined_mags.protein.gtf")
    output:
        bam = os.path.join(MA_OUT, "dna_mapping/{sample}.sorted.bam"),
        counts = os.path.join(MA_OUT, "dna_mapping/{sample}.featureCounts.tsv"),
        summary = os.path.join(MA_OUT, "dna_mapping/{sample}.featureCounts.tsv.summary")
    log:
        os.path.join(MA_LOGS, "dna_mapping/{sample}.log")
    threads: c_abund.get("threads", 16)
    params:
        env = c_abund.get("conda_envs", {}).get("mapping", "mapping")
    shell:
        """
        conda run -p {params.env} bash -c '
        mkdir -p "$(dirname {output.bam})"
        bowtie2 -q --phred33 -p {threads} -x {COMBINED_REF_DIR}/bowtie2_index/combined_mags \
            -1 {input.fq1} -2 {input.fq2} -U {input.u1},{input.u2} | \
        samtools view -@ {threads} -bS - | \
        samtools sort -@ {threads} -o {output.bam} - 
        samtools index {output.bam}
        ' 2> {log}
        
        conda run -p {params.env} featureCounts -T {threads} -p -t CDS -f -g transcript_id -a {input.gtf} -o {output.counts} {output.bam} 2>&1 | tee -a {log}
        """

# 3. Mapping Abundance (RNA - Metatranscriptomics)
rule rna_abundance:
    input:
        fwd = lambda wildcards: os.path.join(RNA_FQ_DIR, f"{get_rna_prefix(wildcards)}_trim_fwd.fq"),
        rev = lambda wildcards: os.path.join(RNA_FQ_DIR, f"{get_rna_prefix(wildcards)}_trim_rev.fq"),
        u = lambda wildcards: os.path.join(RNA_FQ_DIR, f"{get_rna_prefix(wildcards)}_unpair.fq"),
        idx = os.path.join(COMBINED_REF_DIR, "hisat2_index/combined_mags.1.ht2"),
        gtf = os.path.join(COMBINED_REF_DIR, "combined_mags.protein.gtf")
    output:
        bam = os.path.join(MA_OUT, "rna_mapping/{sample}.sorted.bam"),
        counts = os.path.join(MA_OUT, "rna_mapping/{sample}.featureCounts.tsv"),
        summary = os.path.join(MA_OUT, "rna_mapping/{sample}.featureCounts.tsv.summary")
    log:
        os.path.join(MA_LOGS, "rna_mapping/{sample}.log")
    threads: c_abund.get("threads", 16)
    params:
        env = c_abund.get("conda_envs", {}).get("mapping", "mapping")
    shell:
        """
        conda run -p {params.env} bash -c '
        mkdir -p "$(dirname {output.bam})"
        hisat2 -p {threads} -x {COMBINED_REF_DIR}/hisat2_index/combined_mags \
            -1 {input.fwd} -2 {input.rev} -U {input.u} --no-spliced-alignment | \
        samtools view -@ {threads} -bS - | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        ' 2> {log}

        conda run -p {params.env} featureCounts -T {threads} -p -t CDS -f -g transcript_id -M --primary -a {input.gtf} -o {output.counts} {output.bam} 2>&1 | tee -a {log}
        """

# 4. Metrics Calculation (TPM/RPKM)
rule calculate_mag_metrics:
    input:
        counts = os.path.join(MA_OUT, "{type}_mapping/{sample}.featureCounts.tsv"),
        summary = os.path.join(MA_OUT, "{type}_mapping/{sample}.featureCounts.tsv.summary")
    output:
        final = os.path.join(MA_OUT, "{type}_mapping/{sample}.final.tsv")
    shell:
        """
        printf "MAG\\tLocalGene\\tChr\\tStart\\tEnd\\tStrand\\tLength\\tCounts\\tRPK\\tTPM_Mapped\\tRPKM_Mapped\\tRPKM_Lib\\tTPM_Lib\\n" > {output.final}
        
        lib_size=$(awk 'NR>1 {{sum += $2}} END {{printf "%d", sum}}' {input.summary})
        if [ -z "$lib_size" ] || [ "$lib_size" -eq 0 ]; then
            lib_size=1
        fi

        awk -v lib_size="$lib_size" 'BEGIN{{OFS="\\t"}} 
            NR<=2 {{next}} 
            {{ 
                gene_suffix = $1
                chr_field = $2
                split(chr_field, chr_arr, "|")
                mag = chr_arr[1]
                contig_part = chr_arr[2]
                split(gene_suffix, suffix_arr, "_")
                gene_num = suffix_arr[length(suffix_arr)]
                local_gene = contig_part "_" gene_num
                gene_len = $6
                counts = $7
                if (gene_len > 0) {{ rpk = counts / (gene_len / 1000) }} else {{ rpk = 0 }}
                mags[NR] = mag; local_genes[NR] = local_gene; chrs[NR] = $2; starts[NR] = $3; ends[NR] = $4;
                strands[NR] = $5; lens[NR] = $6; cnts[NR] = $7; rpks[NR] = rpk;
                sum_rpk += rpk; total_mapped += counts; n = NR
            }} 
            END {{ 
                mapping_rate = (lib_size > 0) ? (total_mapped / lib_size) : 0
                for(i=3; i<=n; i++) {{ 
                    tpm_mapped = (sum_rpk > 0) ? (rpks[i]/sum_rpk*1000000) : 0
                    tpm_lib = tpm_mapped * mapping_rate
                    rpkm_mapped = (total_mapped > 0) ? (rpks[i] / (total_mapped / 1000000)) : 0
                    rpkm_lib = (lib_size > 0) ? (rpks[i] / (lib_size / 1000000)) : 0
                    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%d\\t%d\\t%.6f\\t%.6f\\t%.6f\\t%.6f\\t%.6f\\n", \
                        mags[i], local_genes[i], chrs[i], starts[i], ends[i], strands[i], lens[i], cnts[i], rpks[i], \
                        tpm_mapped, rpkm_mapped, rpkm_lib, tpm_lib
                }} 
            }}' {input.counts} >> {output.final}
        """

# 5. Merge Abundance results
rule calculate_sample_ratio:
    input:
        dna = os.path.join(MA_OUT, "dna_mapping/{sample}.final.tsv"),
        rna = os.path.join(MA_OUT, "rna_mapping/{sample}.final.tsv")
    output:
        ratio = os.path.join(MA_OUT, "abundance/{sample}.ratio.tsv")
    log:
        os.path.join(MA_LOGS, "{sample}_ratio.log")
    shell:
        """
        python3 modules/metagenomics/scripts/calculate_tpm_ratio.py \
            --dna {input.dna} \
            --rna {input.rna} \
            --output {output.ratio} > {log} 2>&1
        """

rule merge_abundance:
    input:
        files = lambda wildcards: expand(os.path.join(MA_OUT, wildcards.type + "_mapping/{sample}.final.tsv"),
                                         sample=DNA_SAMPLES if wildcards.type == "dna" else RNA_SAMPLES),
        ratios = lambda wildcards: expand(os.path.join(MA_OUT, "abundance/{sample}.ratio.tsv"), 
                                         sample=SAMPLES_WITH_RATIO) if wildcards.type == "rna" else []
    output:
        merged = os.path.join(MA_OUT, "abundance/merged_{type}_tpm_rpkm.final.tsv")
    shell:
        """
        python3 - <<EOF
import pandas as pd
import os
import sys
files = "{input.files}".split()
ratios = "{input.ratios}".split()
if not files: exit(0)
dfs = []
metric_cols = ["TPM_Mapped", "RPKM_Mapped", "RPKM_Lib", "TPM_Lib"]
for f in files:
    sample = os.path.basename(f).replace(".final.tsv", "")
    df = pd.read_csv(f, sep=None, engine='python')
    cols_to_keep = ["MAG", "LocalGene"] + metric_cols
    df_sub = df[cols_to_keep]
    rename_dict = {{col: f"{{sample}}_{{col}}" for col in metric_cols}}
    df_sub = df_sub.rename(columns=rename_dict)
    if not dfs: dfs.append(df_sub)
    else: dfs.append(df_sub.drop(columns=["MAG", "LocalGene"]))

if ratios:
    ratio_dfs = []
    for f in ratios:
        sample = os.path.basename(f).replace(".ratio.tsv", "")
        df_ratio = pd.read_csv(f, sep=None, engine='python')
        if "Ratio" in df_ratio.columns:
            df_ratio = df_ratio[["MAG", "LocalGene", "Ratio"]]
            df_ratio = df_ratio.rename(columns={{"Ratio": f"{{sample}}_Ratio"}})
            ratio_dfs.append(df_ratio.drop(columns=["MAG", "LocalGene"]))
    if ratio_dfs:
        ratio_combined = pd.concat(ratio_dfs, axis=1)
        dfs.append(ratio_combined)

merged_df = pd.concat(dfs, axis=1)
os.makedirs(os.path.dirname("{output.merged}"), exist_ok=True)
merged_df.to_csv("{output.merged}", sep="\\t", index=False)
EOF
        """

# 6. Global Summarization (by Target Genes)
rule summarize_abundance_by_gene:
    input:
        ratios = expand(os.path.join(MA_OUT, "abundance/{sample}.ratio.tsv"), sample=SAMPLES_WITH_RATIO),
        target_genes = os.path.join(MA_DATA, c_abund.get("summarize", {}).get("target_genes", "target.tsv")),
        mag_types = os.path.join(MA_DATA, c_abund.get("summarize", {}).get("mag_types", "type.tsv")),
        annotation_dir = os.path.join(MA_DATA, c_abund.get("summarize", {}).get("annotation_dir", "annotation"))
    output:
        master_mag = os.path.join(MA_OUT, "summary/master_by_mag.tsv"),
        master_mag_wide = os.path.join(MA_OUT, "summary/master_by_mag_wide.tsv"),
        master_type = os.path.join(MA_OUT, "summary/master_by_type.tsv"),
        master_type_wide = os.path.join(MA_OUT, "summary/master_by_type_wide.tsv")
    params:
        abundance_dir = os.path.join(MA_OUT, "abundance"),
        output_dir = os.path.join(MA_OUT, "summary")
    shell:
        """
        python3 modules/metagenomics/scripts/summarize_abundance_by_gene.py \
            --abundance-dir {params.abundance_dir} \
            --annotation-dir {input.annotation_dir} \
            --target-genes {input.target_genes} \
            --mag-types {input.mag_types} \
            --output-dir {params.output_dir}
        """
