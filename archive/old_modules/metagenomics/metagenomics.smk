# Metagenomics Pipeline Module
# Key Steps: Read QC, Binning, Refinement, dRep, GTDB-Tk, CoverM.

import os

wildcard_constraints:
    sample="[^/]+",
    genome="[^/]+"

# 1. Setup & Path Logic
# ------------------------------------------------------------------------------
c_meta = config.get("metagenomics", {})
META_PROJ = os.path.abspath(c_meta.get("project_dir", "."))

# Helper to normalize paths
def norm_path(base, rel):
    if not rel: return base
    if os.path.isabs(rel): return os.path.normpath(rel)
    return os.path.normpath(os.path.join(base, rel))

proj_name = c_meta.get("project_name", "")
base_res = c_meta.get("results_dir", "results")
META_RES = norm_path(META_PROJ, os.path.join(base_res, proj_name)) if proj_name else norm_path(META_PROJ, base_res)

META_DATA = norm_path(META_PROJ, c_meta.get("data_dir", "data"))
RAW_FQ_DIR  = norm_path(META_DATA, c_meta.get("raw_fq_dir", "raw_fq"))
SAMPLES_OUT = norm_path(META_RES,  c_meta.get("samples_out_dir", "samples"))
META_OUT    = norm_path(META_RES,  c_meta.get("metagenomics_out_dir", "MAGs"))
META_LOGS   = os.path.join(META_RES, "logs")

# Global Suffixes for raw fastq files
R1_SUF = c_meta.get("r1_suffix", ".R1.raw.fastq.gz")
R2_SUF = c_meta.get("r2_suffix", ".R2.raw.fastq.gz")

# Resolve samples: prioritize config list, fallback to glob
if c_meta.get("samples"):
    SAMPLES = c_meta["samples"]
    if isinstance(SAMPLES, str):
        SAMPLES = [SAMPLES]
else:
    SAMPLES, = glob_wildcards(os.path.join(RAW_FQ_DIR, f"{{sample}}{R1_SUF}"))

# 2. Target Registration
# ------------------------------------------------------------------------------
if c_meta.get("run_drep", True):
    TARGETS.append(os.path.join(META_OUT, "drep_results", "dereplicated_genomes"))
    TARGETS.append(os.path.join(META_OUT, "drep_results", "data_tables", "Widb.csv"))
else:
    comp = c_meta.get("completeness", 50)
    cont = c_meta.get("contamination", 10)
    TARGETS.extend(expand(os.path.join(SAMPLES_OUT, "{sample}", "04_refinement", "metawrap", "metawrap_" + str(comp) + "_" + str(cont) + "_bins.stats"), sample=SAMPLES))
    TARGETS.append(os.path.join(META_OUT, "all_bins"))

if c_meta.get("run_gtdbtk", True):
    TARGETS.append(os.path.join(META_OUT, "gtdbtk_results", "merged_gtdbtk_summary.tsv"))

if c_meta.get("run_coverm", True):
    TARGETS.append(os.path.join(META_OUT, "abundance", "merged_coverm_abundance.tsv"))

TARGETS.append(os.path.join(META_OUT, "all_bins_quality_summary.tsv"))

# 3. Rules
# ------------------------------------------------------------------------------

# 1. Quality Control using fastp
rule run_fastp:
    input:
        r1 = os.path.join(RAW_FQ_DIR, "{sample}" + R1_SUF),
        r2 = os.path.join(RAW_FQ_DIR, "{sample}" + R2_SUF)
    output:
        r1 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R2.clean.fastq.gz"),
        html = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_fastp.html"),
        json = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_fastp.json")
    log:
        os.path.join(META_LOGS, "qc_fastp", "{sample}.log")
    threads: 16
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c_meta.get("conda_envs", {}).get("bio", "bio"),
        qos_extra = lambda wildcards: "--cut_front --cut_right --cut_window_size 4 --cut_mean_quality 20 --trim_poly_g" \
                    if c_meta.get("fastp", {}).get("quality_filter", True) else "",
        length_req = c_meta.get("fastp", {}).get("length_required", 75)
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} fastp \\
            -i {input.r1} -I {input.r2} \\
            -o {output.r1} -O {output.r2} \\
            --thread {threads} \\
            {params.qos_extra} \\
            --length_required {params.length_req} \\
            --html {output.html} --json {output.json} 2>&1 | tee {log}
        """

# 1.4 Interleave R1 and R2
rule run_reformat_interleave:
    input:
        r1 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R2.clean.fastq.gz")
    output:
        interleaved = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "interleaved", "{sample}_interleaved.fastq.gz")
    log:
        os.path.join(META_LOGS, "interleave", "{sample}.log")
    threads: 8
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c_meta.get("conda_envs", {}).get("bbmap", "bbmap")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} reformat.sh \\
            in1={input.r1} in2={input.r2} \\
            out={output.interleaved} \\
            threads={threads} 2>&1 | tee {log}
        """

# 1.5 Digital Normalization on interleaved file
rule run_bbnorm:
    input:
        interleaved = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "interleaved", "{sample}_interleaved.fastq.gz")
    output:
        norm_interleaved = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "norm", "{sample}_interleaved.norm.fastq.gz")
    log:
        os.path.join(META_LOGS, "bbnorm", "{sample}.log")
    threads: 64
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000,
        slurm_extra = "'--exclusive'"
    params:
        env = c_meta.get("conda_envs", {}).get("bbmap", "bbmap"),
        bbtools_jar = lambda wildcards: os.path.join(c_meta.get("conda_envs", {}).get("bbmap", "bbmap"), "opt", "bbmap-39.81-1", "bbtools.jar")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} java -ea --add-modules jdk.incubator.vector \\
             -Xmx480000M -Xms480000M \\
             -cp {params.bbtools_jar} \\
             jgi.KmerNormalize \\
             in={input.interleaved} \\
             out={output.norm_interleaved} \\
             target=100 threads=32 bits=16 interleaved=t \\
             prefilter=t 2>&1 | tee {log}
        """

# 1.6 De-interleave
rule run_reformat_split:
    input:
        norm_interleaved = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "norm", "{sample}_interleaved.norm.fastq.gz")
    output:
        r1 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "norm", "{sample}_R1.norm.fastq.gz"),
        r2 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "norm", "{sample}_R2.norm.fastq.gz")
    log:
        os.path.join(META_LOGS, "reformat_split", "{sample}.log")
    threads: 8
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c_meta.get("conda_envs", {}).get("bbmap", "bbmap")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} reformat.sh \\
            in={input.norm_interleaved} \\
            out1={output.r1} \\
            out2={output.r2} \\
            threads={threads} 2>&1 | tee {log}
        """

# 2. Assembly using MEGAHIT
rule run_megahit:
    input:
        r1 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "norm", "{sample}_R1.norm.fastq.gz"),
        r2 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "norm", "{sample}_R2.norm.fastq.gz")
    output:
        contigs = os.path.join(SAMPLES_OUT, "{sample}", "02_assembly", "megahit", "{sample}.contigs.fa")
    log:
        os.path.join(META_LOGS, "megahit", "{sample}.log")
    threads: 64
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000,
        slurm_extra = "'--exclusive'"
    params:
        env = c_meta.get("conda_envs", {}).get("megahit", "bio"),
        outdir = os.path.join(SAMPLES_OUT, "{sample}", "02_assembly", "megahit")
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        conda run -p {params.env} megahit \\
            -1 {input.r1} -2 {input.r2} \\
            -o {params.outdir} \\
            -t {threads} \\
            -m $(({resources.mem_mb} * 1024 * 1024)) \\
            --out-prefix {wildcards.sample} 2>&1 | tee {log}
        """

# 2.5 Filter contigs using seqkit
rule run_seqkit_filter:
    input:
        fa = os.path.join(SAMPLES_OUT, "{sample}", "02_assembly", "megahit", "{sample}.contigs.fa")
    output:
        fa = os.path.join(SAMPLES_OUT, "{sample}", "02_assembly", "megahit", "{sample}.filtered.fa")
    log:
        os.path.join(META_LOGS, "seqkit", "{sample}.log")
    threads: 4
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: threads * 8000
    params:
        env = c_meta.get("conda_envs", {}).get("bio", "bio")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} seqkit seq -m 1000 {input.fa} > {output.fa} 2>&1 | tee {log}
        """

# 3. Binning using MetaWRAP
rule run_metawrap_binning:
    input:
        r1 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R2.clean.fastq.gz"),
        asm = os.path.join(SAMPLES_OUT, "{sample}", "02_assembly", "megahit", "{sample}.filtered.fa")
    output:
        directory(os.path.join(SAMPLES_OUT, "{sample}", "03_binning", "metawrap"))
    log:
        os.path.join(META_LOGS, "metawrap_binning", "{sample}.log")
    threads: 64
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000,
        slurm_extra = "'--exclusive'"
    params:
        env = c_meta.get("conda_envs", {}).get("metawrap", "metawrap"),
        outdir = os.path.join(SAMPLES_OUT, "{sample}", "03_binning", "metawrap"),
        tmp_dir = os.path.join(SAMPLES_OUT, "{sample}", "03_binning", "tmp")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.tmp_dir}
        ln -sf {input.r1} {wildcards.sample}_1.fastq
        ln -sf {input.r2} {wildcards.sample}_2.fastq
        
        conda run -p {params.env} metawrap binning \\
            -o {params.outdir} \\
            -t {threads} \\
            -a {input.asm} \\
            --metabat2 --maxbin2 --concoct \\
            {wildcards.sample}_1.fastq {wildcards.sample}_2.fastq 2>&1 | tee {log}
        rm -rf {params.tmp_dir}
        rm {wildcards.sample}_1.fastq {wildcards.sample}_2.fastq
        """

# 4. Bin Refinement using MetaWRAP
COMP = c_meta.get("completeness", 50)
CONT = c_meta.get("contamination", 10)

rule run_metawrap_refinement:
    input:
        bins = os.path.join(SAMPLES_OUT, "{sample}", "03_binning", "metawrap")
    output:
        stats = os.path.join(SAMPLES_OUT, "{sample}", "04_refinement", "metawrap", "metawrap_" + str(COMP) + "_" + str(CONT) + "_bins.stats"),
        bin_dir = directory(os.path.join(SAMPLES_OUT, "{sample}", "04_refinement", "metawrap", "metawrap_" + str(COMP) + "_" + str(CONT) + "_bins"))
    log:
        os.path.join(META_LOGS, "metawrap_refinement", "{sample}.log")
    threads: 64
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000,
        slurm_extra = "'--exclusive'"
    params:
        env = c_meta.get("conda_envs", {}).get("metawrap", "metawrap"),
        outdir = os.path.join(SAMPLES_OUT, "{sample}", "04_refinement", "metawrap"),
        comp = COMP,
        cont = CONT,
        checkm_db = c_meta.get("checkm_db", "")
    shell:
        """
        set -euo pipefail
        if [ -d "{params.checkm_db}" ]; then
            echo "{params.checkm_db}" | conda run -p {params.env} checkm data setRoot {params.checkm_db}
        fi
        
        conda run -p {params.env} metawrap bin_refinement \\
            -o {params.outdir} \\
            -t {threads} \\
            -A {input.bins}/metabat2_bins \\
            -B {input.bins}/maxbin2_bins \\
            -C {input.bins}/concoct_bins \\
            -c {params.comp} -x {params.cont} 2>&1 | tee {log}
        """

# 5. Gathering and Renaming MAGs for dRep
rule gather_bins:
    input:
        stats = expand(os.path.join(SAMPLES_OUT, "{sample}", "04_refinement", "metawrap", "metawrap_" + str(COMP) + "_" + str(CONT) + "_bins.stats"), sample=SAMPLES)
    output:
        bin_dir = directory(os.path.join(META_OUT, "all_bins"))
    log:
        os.path.join(META_LOGS, "utils", "gather_bins.log")
    threads: 4
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env_name = c_meta.get("env_name", "Metagenomics_Project"),
        bin_suffix = "metawrap_" + str(COMP) + "_" + str(CONT) + "_bins"
    run:
        import shutil
        import os
        os.makedirs(output.bin_dir, exist_ok=True)
        for sample in SAMPLES:
            bin_dir = os.path.join(SAMPLES_OUT, sample, "04_refinement", "metawrap", params.bin_suffix)
            if os.path.exists(bin_dir):
                for f in os.listdir(bin_dir):
                    if f.endswith(".fa"):
                        new_name = f"{params.env_name}_{sample}_{f}"
                        shutil.copy2(os.path.join(bin_dir, f), os.path.join(output.bin_dir, new_name))

checkpoint run_drep:
    input:
        bins = os.path.join(META_OUT, "all_bins")
    output:
        derep_dir = directory(os.path.join(META_OUT, "drep_results", "dereplicated_genomes")),
        summary = os.path.join(META_OUT, "drep_results", "data_tables", "Widb.csv")
    log:
        os.path.join(META_LOGS, "run_drep", "drep.log")
    threads: 64
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000,
        slurm_extra = "'--exclusive'"
    params:
        env = c_meta.get("conda_envs", {}).get("drep", "bio"),
        outdir = os.path.join(META_OUT, "drep_results"),
        comp = c_meta.get("completeness", 50),
        cont = c_meta.get("contamination", 10),
        ani = c_meta.get("ani_threshold", 0.95),
        checkm_db = c_meta.get("checkm_db", "")
    shell:
        """
        set -euo pipefail
        if [ -d "{params.checkm_db}" ]; then
            echo "{params.checkm_db}" | conda run -p {params.env} checkm data setRoot {params.checkm_db}
        fi
        
        conda run -p {params.env} dRep dereplicate "{params.outdir}" \\
            -g "{input.bins}"/*.fa \\
            -p {threads} \\
            -comp {params.comp} \\
            -con {params.cont} \\
            -sa {params.ani} 2>&1 | tee {log}
        """

# 8. CoverM Abundance Calculation (Per Sample)
rule run_coverm_per_sample:
    input:
        mag_dir = os.path.join(SAMPLES_OUT, "{sample}", "04_refinement", "metawrap", "metawrap_" + str(COMP) + "_" + str(CONT) + "_bins"),
        r1 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(SAMPLES_OUT, "{sample}", "01_qc", "fastp", "{sample}_R2.clean.fastq.gz")
    output:
        abundance = os.path.join(SAMPLES_OUT, "{sample}", "06_abundance", "{sample}.coverm.tsv")
    log:
        os.path.join(META_LOGS, "coverm", "{sample}.log")
    threads: 16
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c_meta.get("conda_envs", {}).get("coverm", "coverm"),
        tmp_dir = os.path.join(SAMPLES_OUT, "{sample}", "06_abundance", "tmp")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.tmp_dir}
        if [ ! -d "{input.mag_dir}" ] || [ $(ls {input.mag_dir}/*.fa 2>/dev/null | wc -l) -eq 0 ]; then
            echo "Warning: No MAGs found for {wildcards.sample}. Creating empty abundance." | tee {log}
            touch {output.abundance}
            exit 0
        fi

        conda run -p {params.env} bash -c "
            export TMPDIR={params.tmp_dir}
            coverm genome \
                --coupled {input.r1} {input.r2} \
                --genome-fasta-directory {input.mag_dir} \
                -x fa \
                -p bwa-mem2 --min-read-percent-identity 95 --min-read-aligned-percent 75 \
                -m relative_abundance \
                --threads {threads} \
                -o {output.abundance}
        " 2>&1 | tee {log}
            
        rm -rf {params.tmp_dir}
        """

rule merge_coverm_abundance:
    input:
        tsvs = expand(os.path.join(SAMPLES_OUT, "{sample}", "06_abundance", "{sample}.coverm.tsv"), sample=SAMPLES)
    output:
        merged = os.path.join(META_OUT, "abundance", "merged_coverm_abundance.tsv")
    log:
        os.path.join(META_LOGS, "merge", "coverm_abundance.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    run:
        import pandas as pd
        import os
        dfs = []
        for f in input.tsvs:
            if not os.path.exists(f) or os.path.getsize(f) == 0: continue
            sample = os.path.basename(f).replace(".coverm.tsv", "")
            df = pd.read_csv(f, sep='\t')
            df.insert(0, 'Sample', sample)
            dfs.append(df)
        if dfs:
            merged_df = pd.concat(dfs, ignore_index=True)
            merged_df.to_csv(output.merged, sep='\t', index=False)
        else:
            with open(output.merged, 'w') as f:
                pass

# 9. Merge Refinement Quality Stats
rule merge_refinement_stats:
    input:
        stats_files = expand(os.path.join(SAMPLES_OUT, "{sample}", "04_refinement", "metawrap", "metawrap_" + str(COMP) + "_" + str(CONT) + "_bins.stats"), sample=SAMPLES)
    output:
        merged_stats = os.path.join(META_OUT, "all_bins_quality_summary.tsv")
    log:
        os.path.join(META_LOGS, "utils", "merge_stats.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    run:
        import pandas as pd
        import os
        dfs = []
        for f in input.stats_files:
            if not os.path.exists(f) or os.path.getsize(f) == 0: continue
            sample = f.split("/")[-4] # Extract sample name from path
            df = pd.read_csv(f, sep='\t')
            df.insert(0, 'Sample', sample)
            dfs.append(df)
        if dfs:
            merged_df = pd.concat(dfs, ignore_index=True)
            merged_df.to_csv(output.merged_stats, sep='\t', index=False)
        else:
            with open(output.merged_stats, 'w') as f:
                f.write("Sample\tbin\tcompleteness\tcontamination\tGC\tsize\n")
