# MAG Annotation Pipeline Module
# Independent module: reads .fa files from data/mag_fasta/ (MAG_ANNOT_IN)

# 1. Gene Prediction
rule mag_prodigal:
    input:
        os.path.join(MAG_ANNOT_IN, "{mag}.fa")
    output:
        gff = os.path.join(MAG_ANNOT_OUT, "{mag}/01_prodigal/{mag}.gff"),
        faa = os.path.join(MAG_ANNOT_OUT, "{mag}/01_prodigal/{mag}.faa")
    log:
        os.path.join(MAG_ANNOT_LOGS, "mag_prodigal/{mag}.log")
    threads: 8
    resources:
        runtime="7d"
    conda:
        "../envs/bio.yaml"
    shell:
        "prodigal -i {input} -a {output.faa} -f gff -o {output.gff} -p single -q 2>&1 | tee {log}"

# 2. RNA Prediction
rule mag_rna_prediction:
    input:
        os.path.join(MAG_ANNOT_IN, "{mag}.fa")
    output:
        trna = os.path.join(MAG_ANNOT_OUT, "{mag}/02_rna/{mag}.trna.txt"),
        rrna = os.path.join(MAG_ANNOT_OUT, "{mag}/02_rna/{mag}.rrna.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/mag_rna/{mag}.log")
    threads: 8
    resources:
        runtime="7d"
    conda:
        "../envs/bio.yaml"
    shell:
        """
        aragorn -t -gc1 -w -o {output.trna} {input} 2>> {log}
        barrnap --threads {threads} --quiet {input} > {output.rrna} 2>> {log}
        """

# 3. Quality Assessment
rule mag_checkm2:
    input:
        mag_dir = MAG_ANNOT_IN
    output:
        summary = os.path.join(MAG_ANNOT_OUT, "quality/checkm2_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/mag_checkm2/checkm2.log")
    threads: 64
    resources:
        runtime="7d",
        slurm_extra="'--exclusive'"
    conda:
        "../envs/checkm2.yaml"
    params:
        outdir = os.path.join(MAG_ANNOT_OUT, "quality/checkm2_out")
    shell:
        """
        rm -rf {params.outdir}
        checkm2 predict --threads {threads} --input {input.mag_dir} --output-directory {params.outdir} --extension fa 2>&1 | tee {log}
        cp {params.outdir}/quality_report.tsv {output.summary}
        """

# 4. Taxonomic Assignment
rule mag_gtdbtk:
    input:
        mag_dir = MAG_ANNOT_IN
    output:
        summary = os.path.join(MAG_ANNOT_OUT, "taxonomy/gtdbtk.summary.tsv")
    log:
        os.path.join(MAG_ANNOT_LOGS, "mag_gtdbtk/gtdbtk.log")
    threads: 16
    resources:
        runtime="7d"
    params:
        env = config["mag_annotation"]["conda_envs"]["gtdbtk"],
        outdir = os.path.join(MAG_ANNOT_OUT, "taxonomy/gtdbtk_out")
    shell:
        """
        rm -rf {params.outdir}
        conda run -p {params.env} gtdbtk classify_wf --cpus {threads} --genome_dir {input.mag_dir} --out_dir {params.outdir} --extension fa 2>&1 | tee {log}
        # Correctly handle both Bacteria and Archaea summaries
        if [ -f "{params.outdir}/gtdbtk.bac120.summary.tsv" ] && [ -f "{params.outdir}/gtdbtk.ar122.summary.tsv" ]; then
            head -n 1 {params.outdir}/gtdbtk.bac120.summary.tsv > {output.summary}
            tail -n +2 {params.outdir}/gtdbtk.bac120.summary.tsv >> {output.summary}
            tail -n +2 {params.outdir}/gtdbtk.ar122.summary.tsv >> {output.summary}
        elif [ -f "{params.outdir}/gtdbtk.bac120.summary.tsv" ]; then
            cp {params.outdir}/gtdbtk.bac120.summary.tsv {output.summary}
        elif [ -f "{params.outdir}/gtdbtk.ar122.summary.tsv" ]; then
            cp {params.outdir}/gtdbtk.ar122.summary.tsv {output.summary}
        else
            touch {output.summary}
        fi
        """

# 5. Functional Annotation (HMM based)
rule mag_hmm_annotation:
    input:
        faa = os.path.join(MAG_ANNOT_OUT, "{mag}/01_prodigal/{mag}.faa")
    output:
        anno = os.path.join(MAG_ANNOT_OUT, "{mag}/03_hmm/{mag}.hmm.tsv")
    log:
        os.path.join(MAG_ANNOT_LOGS, "mag_hmm/{mag}.log")
    conda:
        "../envs/bio.yaml"
    shell:
        "hmmsearch --cpu {threads} --tblout {output.anno} {params.db} {input.faa} 2>&1 | tee {log}"

# 6. Phenotype Prediction
rule mag_genomespot:
    input:
        faa = os.path.join(MAG_ANNOT_OUT, "{mag}/01_prodigal/{mag}.faa")
    output:
        pheno = os.path.join(MAG_ANNOT_OUT, "{mag}/04_phenotype/{mag}.phenotype.tsv")
    log:
        os.path.join(MAG_ANNOT_LOGS, "mag_genomespot/{mag}.log")
    threads: 8
    resources:
        runtime="7d"
    params:
        env = config["mag_annotation"]["conda_envs"].get("genomespot", "genomespot")
    shell:
        "conda run -p {params.env} genomespot {input.faa} --threads {threads} > {output.pheno} 2> {log}"

# 7. Consolidation
rule mag_summary:
    input:
        mag_dir = MAG_ANNOT_IN,
        quality = os.path.join(MAG_ANNOT_OUT, "quality/checkm2_summary.tsv"),
        taxonomy = os.path.join(MAG_ANNOT_OUT, "taxonomy/gtdbtk.summary.tsv"),
        anno = expand(os.path.join(MAG_ANNOT_OUT, "{mag}/03_hmm/{mag}.hmm.tsv"), mag=MAGS)
    output:
        summary = os.path.join(MAG_ANNOT_OUT, "merged_summary.tsv")
    log:
        os.path.join(MAG_ANNOT_LOGS, "mag_summary/merge.log")
    threads: 1
    resources:
        runtime="7d"
    params:
        running_sh = config.get("mag_annotation", {}).get("running_sh", "scripts/mag_utils.sh")
    shell:
        """
        # Example call to user script
        source {params.running_sh}
        run_mags_summarize "all_mags" 1 {input.mag_dir} {MAG_ANNOT_OUT} {MAG_ANNOT_OUT} {input.quality} 2>&1 | tee {log}
        """
