# Rules for circular DNA/mags annotation and analysis

# 1. Prodigal for gene prediction
rule prodigal_circular:
    input:
        "results/circular_mags/{mag}.fa"
    output:
        faa="results/circular_annotation/{mag}/{mag}.faa",
        gff="results/circular_annotation/{mag}/{mag}.gff"
    log:
        "logs/prodigal/{mag}.log"
    threads: 8
    resources:
        runtime="7d"
    shell:
        "prodigal -i {input} -a {output.faa} -o {output.gff} -f gff -p meta 2> {log}"

# 2. HMM search for specific markers
rule hmmsearch_circular:
    input:
        faa="results/circular_annotation/{mag}/{mag}.faa"
    output:
        table="results/circular_annotation/{mag}/{mag}_{db}.hmm.tsv"
    log:
        "logs/hmmsearch/{mag}_{db}.log"
    params:
        db_path=lambda wildcards: config["hmm_databases"][wildcards.db]
    threads: 16
    resources:
        runtime="7d"
    shell:
        "hmmsearch --cpu {threads} --tblout {output.table} {params.db_path} {input.faa} > /dev/null 2> {log}"

# 3. BLAST comparison against reference
rule blast_circular:
    input:
        query="results/circular_mags/{mag}.fa",
        db=config["blast_db"]
    output:
        fmt6="results/circular_annotation/{mag}/{mag}_blast.tsv"
    log:
        "logs/blast/{mag}.log"
    threads: 16
    resources:
        runtime="7d"
    shell:
        "blastn -query {input.query} -db {input.db} -out {output.fmt6} -outfmt 6 -num_threads {threads} 2> {log}"

# 4. tRNA and rRNA prediction
rule rna_prediction_circular:
    input:
        "results/circular_mags/{mag}.fa"
    output:
        trna="results/circular_annotation/{mag}/{mag}.trna.txt",
        rrna="results/circular_annotation/{mag}/{mag}.rrna.txt"
    log:
        "logs/rna_pred/{mag}.log"
    threads: 8
    resources:
        runtime="7d"
    shell:
        """
        aragorn -t -gc1 -w -o {output.trna} {input} 2>> {log}
        barrnap --threads {threads} --quiet {input} > {output.rrna} 2>> {log}
        """

# 5. Bakta functional annotation (Integrated)
rule bakta_circular:
    input:
        fasta = "results/circular_mags/{mag}.fa",
        db = config["annotation"]["parameters"]["bakta_db"],
        custom_hmms = config["annotation"]["parameters"]["hmms"]
    output:
        gff = "results/circular_annotation/{mag}/{mag}.bakta.gff",
        tsv = "results/circular_annotation/{mag}/{mag}.bakta.tsv",
        faa = "results/circular_annotation/{mag}/{mag}.bakta.faa",
        dir = directory("results/circular_annotation/{mag}/bakta_out")
    log:
        "logs/bakta/{mag}.log"
    threads: 16
    resources:
        runtime="7d"
    conda:
        "../envs/bakta.yaml"
    shell:
        """
        bakta --db {input.db} \
              --threads {threads} \
              --output {output.dir} \
              --prefix {wildcards.mag} \
              --hmms {input.custom_hmms} \
              --force \
              {input.fasta} 2>&1 | tee {log}
        
        cp {output.dir}/{wildcards.mag}.gff {output.gff}
        cp {output.dir}/{wildcards.mag}.tsv {output.tsv}
        cp {output.dir}/{wildcards.mag}.faa {output.faa}
        """

# 6. Summarize circular annotations
rule summarize_circular:
    input:
        faa=expand("results/circular_annotation/{mag}/{mag}.faa", mag=CIRCULAR_MAGS),
        bakta=expand("results/circular_annotation/{mag}/{mag}.bakta.tsv", mag=CIRCULAR_MAGS),
        hmm=expand("results/circular_annotation/{mag}/{mag}_{db}.hmm.tsv", mag=CIRCULAR_MAGS, db=config["hmm_databases"].keys()),
        blast=expand("results/circular_annotation/{mag}/{mag}_blast.tsv", mag=CIRCULAR_MAGS),
        rna=expand("results/circular_annotation/{mag}/{mag}.trna.txt", mag=CIRCULAR_MAGS)
    output:
        report="results/circular_annotation/circular_mags_summary.tsv"
    log:
        "logs/summarize_circular.log"
    threads: 1
    resources:
        runtime="7d"
    shell:
        "python modules/metagenomics/scripts/summarize_circular.py --input_dir results/circular_annotation --output {output.report} 2> {log}"
