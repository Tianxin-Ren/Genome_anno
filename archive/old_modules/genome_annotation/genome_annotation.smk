# Integrated Genome Annotation Pipeline (Prokaryotic)
# -----------------------------------------------
# Author: Antigravity Architecture Steward
# Function: Comprehensive structural and functional annotation of prokaryotic genomes.
# Key Steps: Bakta (Structural), eggNOG (Functional), KofamScan (KO), dbCAN (CAZy).

import os

wildcard_constraints:
    sample="[^/]+",
    genome="[^/]+",
    plasmid="[^/]+"

# ------------------------------------------------------------------------------
# 1. Setup & Wildcards
# ------------------------------------------------------------------------------
c = config.get("genome_annotation", {})
GA_PROJ = os.path.abspath(c.get("project_dir", "."))

# All paths below are resolved relative to GA_PROJ
MODE = c.get("mode", "standard") 

if MODE == "hybrid":
    GA_IN_DIR  = os.path.join(GA_PROJ, c.get("bins_dir", "data/bins"))
    GA_ASSEMBLY = os.path.join(GA_PROJ, c.get("assembly_fasta", "data/assembly.fasta"))
    
    # Robust output directory handling: prevent results/results/
    out_val = c.get("output_dir", "results/annotation")
    if out_val.startswith("results/") and GA_PROJ.rstrip("/").endswith("/results"):
        GA_OUT_DIR = os.path.join(GA_PROJ, out_val.replace("results/", "", 1))
    else:
        GA_OUT_DIR = os.path.join(GA_PROJ, out_val)
    GA_ASSEMBLY_OUT = os.path.join(GA_OUT_DIR, "assembly_annotation")
else:
    GA_IN_DIR  = os.path.join(GA_PROJ, c.get("input_dir", "data/genomes"))
    out_val = c.get("output_dir", "results/genome_annotation")
    GA_OUT_DIR = os.path.join(GA_PROJ, out_val)

# Identify genomes (bins) and their extensions
_genomes_fa, = glob_wildcards(os.path.join(GA_IN_DIR, "{genome}.fa"))
_genomes_fasta, = glob_wildcards(os.path.join(GA_IN_DIR, "{genome}.fasta"))
GENOMES = list(set(_genomes_fa) | set(_genomes_fasta))

# Helper to get the correct input path for a genome
def get_genome_input(wildcards):
    path_fa = os.path.join(GA_IN_DIR, f"{wildcards.genome}.fa")
    if os.path.exists(path_fa):
        return path_fa
    return os.path.join(GA_IN_DIR, f"{wildcards.genome}.fasta")

# Identify standalone plasmids
GA_PLASMID_DIR = os.path.join(GA_PROJ, c.get("plasmids_dir", "data/plasmids"))
_plasmids_fa, = glob_wildcards(os.path.join(GA_PLASMID_DIR, "{plasmid}.fa"))
_plasmids_fasta, = glob_wildcards(os.path.join(GA_PLASMID_DIR, "{plasmid}.fasta"))
PLASMIDS = list(set(_plasmids_fa) | set(_plasmids_fasta))

# Helper to get the correct input path for a plasmid
def get_plasmid_input(wildcards):
    path_fa = os.path.join(GA_PLASMID_DIR, f"{wildcards.plasmid}.fa")
    if os.path.exists(path_fa):
        return path_fa
    return os.path.join(GA_PLASMID_DIR, f"{wildcards.plasmid}.fasta")

# ------------------------------------------------------------------------------
# Module 0: Assembly-level Annotation (Hybrid & Scatter-Gather Mode)
# ------------------------------------------------------------------------------
if MODE == "hybrid":
    import math
    
    # Dynamic Split Calculation
    MB_PER_CHUNK = c.get("mb_per_chunk", 0)
    FIXED_CHUNKS = c.get("split_chunks", 0)
    
    if MB_PER_CHUNK > 0 and os.path.exists(GA_ASSEMBLY):
        size_mb = os.path.getsize(GA_ASSEMBLY) / (1024 * 1024)
        SPLIT_COUNT = max(1, math.ceil(size_mb / MB_PER_CHUNK))
    elif FIXED_CHUNKS > 0:
        SPLIT_COUNT = FIXED_CHUNKS
    else:
        SPLIT_COUNT = 1 # Fallback to no splitting

    CHUNKS = [f"{i:03d}" for i in range(1, SPLIT_COUNT + 1)]

    rule split_assembly:
        input:
            GA_ASSEMBLY
        output:
            expand(os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.fasta"), chunk=CHUNKS)
        log:
            os.path.join(GA_OUT_DIR, "logs", "utils", "split_assembly.log")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("utils", "utils"),
            outdir = os.path.join(GA_ASSEMBLY_OUT, "chunks"),
            num = SPLIT_COUNT
        shell:
            """
            set -euo pipefail
            conda run -p {params.env} python3 modules/genome_annotation/scripts/split_fasta.py \
                -i {input} \
                -o {params.outdir} \
                -n {params.num} 2>&1 | tee {log}
            """

    rule run_bakta_chunk:
        input:
            fasta = os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.fasta"),
            db = c.get("bakta", {}).get("db", "")
        output:
            gff = os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.gff3"),
            faa = os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.faa"),
            tsv = os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.tsv")
        log:
            os.path.join(GA_OUT_DIR, "logs", "bakta", "chunk_{chunk}.log")
        threads: 16 # Chunks are smaller, so we use fewer threads per chunk
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("bakta", "bakta"),
            outdir = os.path.join(GA_ASSEMBLY_OUT, "chunks", "bakta_{chunk}")
        shell:
            """
            set -euo pipefail
            rm -rf {params.outdir}
            conda run -p {params.env} bakta \
                --db {input.db} \
                --threads {threads} \
                --output {params.outdir} \
                --prefix chunk_{wildcards.chunk} \
                --keep-contig-headers \
                --force {input.fasta} 2>&1 | tee {log}
            
            cp {params.outdir}/chunk_{wildcards.chunk}.gff3 {output.gff}
            cp {params.outdir}/chunk_{wildcards.chunk}.faa {output.faa}
            cp {params.outdir}/chunk_{wildcards.chunk}.tsv {output.tsv}
            rm -rf {params.outdir}
            """

    rule merge_bakta_assembly:
        input:
            gffs = expand(os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.gff3"), chunk=CHUNKS),
            faas = expand(os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.faa"), chunk=CHUNKS),
            tsvs = expand(os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.tsv"), chunk=CHUNKS)
        output:
            gff = os.path.join(GA_ASSEMBLY_OUT, "assembly.gff3"),
            faa = os.path.join(GA_ASSEMBLY_OUT, "assembly.faa"),
            tsv = os.path.join(GA_ASSEMBLY_OUT, "assembly.tsv")
        log:
            os.path.join(GA_OUT_DIR, "logs", "utils", "merge_bakta.log")
        threads: 8
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("utils", "utils"),
            chunks_dir = os.path.join(GA_ASSEMBLY_OUT, "chunks"),
            fixed_dir = os.path.join(GA_ASSEMBLY_OUT, "fixed_chunks")
        shell:
            """
            set -euo pipefail
            mkdir -p {params.fixed_dir}
            
            # 1. Fix IDs for each chunk
            for chunk in {CHUNKS}; do
                conda run -p {params.env} python3 modules/genome_annotation/scripts/fix_bakta_ids.py \
                    {params.chunks_dir}/chunk_$chunk.gff3 \
                    {params.chunks_dir}/chunk_$chunk.faa \
                    {params.chunks_dir}/chunk_$chunk.tsv \
                    {params.chunks_dir}/chunk_$chunk.fasta \
                    {params.fixed_dir}/chunk_$chunk.gff3 \
                    {params.fixed_dir}/chunk_$chunk.faa \
                    {params.fixed_dir}/chunk_$chunk.tsv 2>&1 | tee -a {log}
            done

            # 2. Merge GFF (Keep header from first, skip for others)
            FIXED_GFFS=$(ls {params.fixed_dir}/*.gff3 | sort)
            head -n 5 $(echo $FIXED_GFFS | awk '{{print $1}}') > {output.gff}
            for f in $FIXED_GFFS; do grep -v '^#' $f >> {output.gff}; done
            
            # 3. Merge FAA
            FIXED_FAAS=$(ls {params.fixed_dir}/*.faa | sort)
            cat $FIXED_FAAS > {output.faa}
            
            # 4. Merge TSV
            FIXED_TSVS=$(ls {params.fixed_dir}/*.tsv | sort)
            head -n 6 $(echo $FIXED_TSVS | awk '{{print $1}}') > {output.tsv}
            for f in $FIXED_TSVS; do tail -n +7 $f >> {output.tsv}; done
            
            rm -rf {params.fixed_dir}
            """

    rule run_eggnog_chunk:
        input: faa = os.path.join(GA_ASSEMBLY_OUT, "chunks/chunk_{chunk}.faa")
        output: anno = os.path.join(GA_ASSEMBLY_OUT, "chunks/chunk_{chunk}.emapper.annotations")
        log: os.path.join(GA_OUT_DIR, "logs/eggnog/chunk_{chunk}.log")
        threads: 16
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("eggnog", "eggnog"),
            db = c.get("eggnog", {}).get("db", ""),
            prefix = os.path.join(GA_ASSEMBLY_OUT, "chunks/chunk_{chunk}")
        shell:
            "conda run -p {params.env} emapper.py -m mmseqs --dbmem -i {input.faa} --output {params.prefix} --cpu {threads} --data_dir {params.db} --override 2>&1 | tee {log}"

    rule merge_eggnog_assembly:
        input:
            annos = expand(os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.emapper.annotations"), chunk=CHUNKS)
        output:
            anno = os.path.join(GA_ASSEMBLY_OUT, "assembly.emapper.annotations")
        log:
            os.path.join(GA_OUT_DIR, "logs", "utils", "merge_eggnog.log")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        shell:
            """
            set -euo pipefail
            grep '^#' {input.annos[0]} > {output.anno}
            for f in {input.annos}; do grep -v '^#' $f >> {output.anno}; done
            """

    rule run_kofamscan_chunk:
        input: faa = os.path.join(GA_ASSEMBLY_OUT, "chunks/chunk_{chunk}.faa")
        output: txt = os.path.join(GA_ASSEMBLY_OUT, "chunks/chunk_{chunk}.kofam.txt")
        log: os.path.join(GA_OUT_DIR, "logs/kofamscan/chunk_{chunk}.log")
        threads: 16
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("kofamscan", "kofamscan"),
            db_profiles = os.path.join(c.get("kofamscan", {}).get("db", ""), "profiles"),
            db_kolist = os.path.join(c.get("kofamscan", {}).get("db", ""), "ko_list"),
            tmp_dir = os.path.join(GA_ASSEMBLY_OUT, "chunks/kofam_tmp_{chunk}")
        shell:
            """
            export LANG=C; export LC_ALL=C
            rm -rf {params.tmp_dir} && mkdir -p {params.tmp_dir}
            conda run -p {params.env} exec_annotation -f detail-tsv -o {output.txt} --cpu {threads} \
                --tmp {params.tmp_dir} -p {params.db_profiles} -k {params.db_kolist} {input.faa} 2>&1 | tee {log}
            rm -rf {params.tmp_dir}
            """

    rule merge_kofam_assembly:
        input: txts = expand(os.path.join(GA_ASSEMBLY_OUT, "chunks/chunk_{chunk}.kofam.txt"), chunk=CHUNKS)
        output: txt = os.path.join(GA_ASSEMBLY_OUT, "assembly.kofam.txt")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        shell:
            """
            set -euo pipefail
            grep '^#' {input.txts[0]} > {output.txt}
            for f in {input.txts}; do grep -v '^#' $f >> {output.txt}; done
            """

    rule run_dbcan_chunk:
        input:
            faa = os.path.join(GA_ASSEMBLY_OUT, "chunks", "chunk_{chunk}.faa")
        output:
            dir = directory(os.path.join(GA_ASSEMBLY_OUT, "chunks", "dbcan_{chunk}"))
        log:
            os.path.join(GA_OUT_DIR, "logs", "dbcan", "chunk_{chunk}.log")
        threads: 8
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("dbcan", "dbcan"),
            db = c.get("dbcan", {}).get("db", "")
        shell:
            """
            set -euo pipefail
            conda run -p {params.env} run_dbcan {input.faa} protein \\
                --out_dir {output.dir} \\
                --db_dir {params.db} \\
                --dbCANFile dbCAN.hmm \\
                --dbcan_thread {threads} \\
                --hmm_cpu {threads} \\
                --dia_cpu {threads} 2>&1 | tee {log}
            """

    rule merge_dbcan_assembly:
        input:
            dirs = expand(os.path.join(GA_ASSEMBLY_OUT, "chunks", "dbcan_{chunk}"), chunk=CHUNKS)
        output:
            dir = directory(os.path.join(GA_ASSEMBLY_OUT, "dbcan"))
        log:
            os.path.join(GA_OUT_DIR, "logs", "utils", "merge_dbcan.log")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        shell:
            """
            set -euo pipefail
            mkdir -p {output.dir}
            # For each expected output file in dbcan results, merge them
            for f_name in overview.txt cgc.out dbsub.out hmmer.out diamond.out; do
                head -n 1 {input.dirs[0]}/$f_name > {output.dir}/$f_name 2>/dev/null || true
                for d in {input.dirs}; do
                    tail -n +2 $d/$f_name >> {output.dir}/$f_name 2>/dev/null || true
                done
            done
            """

    rule subset_assembly_to_bins:
        input:
            bin_fa = get_genome_input,
            gff = os.path.join(GA_ASSEMBLY_OUT, "assembly.gff3"),
            faa = os.path.join(GA_ASSEMBLY_OUT, "assembly.faa"),
            tsv = os.path.join(GA_ASSEMBLY_OUT, "assembly.tsv"),
            eggnog = os.path.join(GA_ASSEMBLY_OUT, "assembly.emapper.annotations"),
            kofam = os.path.join(GA_ASSEMBLY_OUT, "assembly.kofam.txt"),
            dbcan = os.path.join(GA_ASSEMBLY_OUT, "dbcan")
        output:
            gff = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.gff3"),
            faa = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.faa"),
            tsv = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.tsv"),
            eggnog = os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "{genome}.emapper.annotations"),
            kofam = os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "{genome}.kofam.txt"),
            dbcan = directory(os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "dbcan"))
        log:
            os.path.join(GA_OUT_DIR, "logs", "subset", "{genome}.log")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("utils", "utils")
        shell:
            """
            set -euo pipefail
            conda run -p {params.env} python3 modules/genome_annotation/scripts/subset_assembly_results.py \\
                --bin_fasta {input.bin_fa} \\
                --assembly_gff {input.gff} \\
                --assembly_faa {input.faa} \\
                --assembly_tsv {input.tsv} \\
                --assembly_eggnog {input.eggnog} \\
                --assembly_kofam {input.kofam} \\
                --assembly_dbcan {input.dbcan} \\
                --out_gff {output.gff} \\
                --out_faa {output.faa} \\
                --out_tsv {output.tsv} \\
                --out_eggnog {output.eggnog} \\
                --out_kofam {output.kofam} \\
                --out_dbcan {output.dbcan} 2>&1 | tee {log}
            """

# ------------------------------------------------------------------------------
# Module 1: Quality Assessment (QC)
# ------------------------------------------------------------------------------
rule run_checkm2:
    input:
        mag_dir = GA_IN_DIR
    output:
        summary = os.path.join(GA_OUT_DIR, "quality", "checkm2_summary.tsv")
    log:
        os.path.join(GA_OUT_DIR, "logs", "checkm2", "checkm2.log")
    threads: c.get("checkm2", {}).get("threads", 64)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("checkm2", ""),
        outdir = os.path.join(GA_OUT_DIR, "quality", "checkm2_out"),
        db = c.get("checkm2", {}).get("db", "")
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        conda run -p {params.env} checkm2 predict \\
            --threads {threads} \\
            --input {input.mag_dir} \\
            --output-directory {params.outdir} \\
            --extension fa \\
            --database {params.db} 2>&1 | tee {log}
        cp {params.outdir}/quality_report.tsv {output.summary}
        """

rule run_quast:
    input:
        f = os.path.join(GA_IN_DIR, "{genome}.fa")
    output:
        dir = directory(os.path.join(GA_OUT_DIR, "{genome}", "01_qc", "quast"))
    log:
        os.path.join(GA_OUT_DIR, "logs", "quast", "{genome}.log")
    threads: 8
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("quast", ""),
        extra = c.get("quast", {}).get("params", "")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} quast.py {params.extra} \\
            -t {threads} \\
            -o {output.dir} \\
            {input.f} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Module 2: Taxonomic Assignment
# ------------------------------------------------------------------------------
rule run_gtdbtk:
    input:
        mag_dir = GA_IN_DIR
    output:
        summary = os.path.join(GA_OUT_DIR, "taxonomy", "gtdbtk.summary.tsv")
    log:
        os.path.join(GA_OUT_DIR, "logs", "gtdbtk", "gtdbtk.log")
    threads: c.get("gtdbtk", {}).get("threads", 64)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("gtdbtk", "gtdbtk"),
        outdir = os.path.join(GA_OUT_DIR, "taxonomy", "gtdbtk_out"),
        db = c.get("gtdbtk", {}).get("db", "")
    shell:
        """
        set -euo pipefail
        if [ -n "{params.db}" ]; then export GTDBTK_DATA_PATH={params.db}; fi
        rm -rf {params.outdir}
        conda run -p {params.env} gtdbtk classify_wf \\
            --cpus {threads} \\
            --genome_dir {input.mag_dir} \\
            --out_dir {params.outdir} \\
            --extension fa \\
            --skip_ani_screen 2>&1 | tee {log}
        
        # Merge bacterial and archaeal summaries if they exist
        SUMMARIES=$(ls {params.outdir}/gtdbtk.*.summary.tsv 2>/dev/null || true)
        if [ -n "$SUMMARIES" ]; then
            head -n 1 $(echo $SUMMARIES | awk '{{print $1}}') > {output.summary}
            for s in $SUMMARIES; do
                tail -n +2 $s >> {output.summary}
            done
        else
            touch {output.summary}
        fi
        """

# ------------------------------------------------------------------------------
# Module 3: Structural Annotation (Bakta Hub)
# ------------------------------------------------------------------------------
if MODE == "standard":
    rule run_bakta:
        doc:
            """
            Rapid and standard structural annotation using Bakta.
            Includes tRNA, rRNA, and protein-coding gene identification.
            """
        input:
            f = os.path.join(GA_IN_DIR, "{genome}.fa"),
            db = c.get("bakta", {}).get("db", "")
        output:
            gff = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.gff3"),
            faa = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.faa"),
            tsv = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.tsv")
        log:
            os.path.join(GA_OUT_DIR, "logs", "bakta", "{genome}.log")
        threads: c.get("bakta", {}).get("threads", 32)
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("bakta", "bakta"),
            outdir = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "bakta_out")
        shell:
            """
            set -euo pipefail
            rm -rf {params.outdir}
            conda run -p {params.env} bakta \\
                --db {params.db} \\
                --threads {threads} \\
                --output {params.outdir} \\
                --prefix {wildcards.genome} \\
                --keep-contig-headers \\
                --force {input.f} 2>&1 | tee {log}
            
            cp {params.outdir}/{wildcards.genome}.gff3 {output.gff}
            cp {params.outdir}/{wildcards.genome}.faa {output.faa}
            cp {params.outdir}/{wildcards.genome}.tsv {output.tsv}
            """

# ------------------------------------------------------------------------------
# Module 4: Functional Annotation (using Bakta .faa)
# ------------------------------------------------------------------------------
if MODE == "standard":
    rule run_eggnog:
        doc:
            """
            Functional annotation of predicted protein sequences using eggNOG-mapper.
            Assigns Orthologous Groups, COG categories, and Gene Ontology terms.
            """
        input:
            faa = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.faa")
        output:
            anno = os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "{genome}.emapper.annotations")
        log:
            os.path.join(GA_OUT_DIR, "logs", "eggnog", "{genome}.log")
        threads: c.get("eggnog", {}).get("threads", 16)
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("eggnog", "eggnog"),
            db = c.get("eggnog", {}).get("db", ""),
            out_prefix = os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "{genome}")
        shell:
            """
            set -euo pipefail
            conda run -p {params.env} emapper.py -m mmseqs --dbmem -i {input.faa} --output {params.out_prefix} --cpu {threads} --data_dir {params.db} 2>&1 | tee {log}
            """

    rule run_kofamscan:
        input:
            faa = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.faa")
        output:
            txt = os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "{genome}.kofam.txt")
        log:
            os.path.join(GA_OUT_DIR, "logs", "kofamscan", "{genome}.log")
        threads: c.get("kofamscan", {}).get("threads", 16)
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("kofamscan", "kofamscan"),
            db_profiles = os.path.join(c.get("kofamscan", {}).get("db", ""), "profiles"),
            db_kolist = os.path.join(c.get("kofamscan", {}).get("db", ""), "ko_list"),
            tmp_dir = os.path.join(GA_OUT_DIR, "{genome}/03_functional/kofam_tmp")
        shell:
            """
            export LANG=C
            export LC_ALL=C
            rm -rf {params.tmp_dir} && mkdir -p {params.tmp_dir}
            conda run -p {params.env} exec_annotation \
                -f detail-tsv \
                -o {output.txt} \
                --cpu {threads} \
                --tmp {params.tmp_dir} \
                -p {params.db_profiles} \
                -k {params.db_kolist} \
                {input.faa} 2>&1 | tee {log}
            rm -rf {params.tmp_dir}
            """

    rule run_dbcan:
        input:
            faa = os.path.join(GA_OUT_DIR, "{genome}/02_bakta/{genome}.faa")
        output:
            dir = directory(os.path.join(GA_OUT_DIR, "{genome}/03_functional/dbcan"))
        log:
            os.path.join(GA_OUT_DIR, "logs/dbcan/{genome}.log")
        threads: c.get("dbcan", {}).get("threads", 8)
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("dbcan", "dbcan"),
            db = c.get("dbcan", {}).get("db", "")
        shell:
            """
            set -euo pipefail
            conda run -p {params.env} run_dbcan {input.faa} protein --out_dir {output.dir} --db_dir {params.db} --dbCANFile dbCAN.hmm --dbcan_thread {threads} --hmm_cpu {threads} --dia_cpu {threads} 2>&1 | tee {log}
            """

# ------------------------------------------------------------------------------
# Module 5: Phenotype Prediction (Always Bin-level, depends on Bakta proteins)
# ------------------------------------------------------------------------------
rule run_metabolic:
    input:
        faa = os.path.join(GA_OUT_DIR, "{genome}/02_bakta/{genome}.faa")
    output:
        dir = directory(os.path.join(GA_OUT_DIR, "{genome}/04_phenotype/metabolic"))
    log:
        os.path.join(GA_OUT_DIR, "logs/metabolic/{genome}.log")
    threads: c.get("metabolic", {}).get("threads", 32)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("metabolic", "metabolic")
    shell:
        "conda run -p {params.env} METABOLIC-g.pl -in-gn {input.faa} -o {output.dir} -t {threads} 2>&1 | tee {log}"

rule run_genomespot:
    input:
        fna = os.path.join(GA_IN_DIR, "{genome}.fa"),
        faa = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.faa")
    output:
        pheno = os.path.join(GA_OUT_DIR, "{genome}", "04_phenotype", "{genome}.genomespot.tsv")
    log:
        os.path.join(GA_OUT_DIR, "logs", "genomespot", "{genome}.log")
    threads: 8
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("genomespot", "genomespot"),
        models = c.get("genomespot", {}).get("db", ""),
        out_prefix = os.path.join(GA_OUT_DIR, "{genome}", "04_phenotype", "{genome}"),
        tmp_fna = os.path.join(GA_OUT_DIR, "{genome}", "04_phenotype", "{genome}.tmp_simple.fa")
    shell:
        """
        set -euo pipefail
        # 1. Safety check
        if [ ! -s {input.faa} ] || [ $(grep -c ">" {input.faa}) -eq 0 ]; then
            echo "Warning: No proteins found for {wildcards.genome}." | tee {log}
            echo -e "genome\tprediction\tprobability" > {output.pheno}
            exit 0
        fi

        # 2. Simplify FASTA headers
        sed 's/ctg.*//' {input.fna} > {params.tmp_fna}
        
        # 3. Run GenomeSpot
        conda run -p {params.env} python -m genome_spot.genome_spot \\
            --models {params.models} \\
            --contigs {params.tmp_fna} \\
            --proteins {input.faa} \\
            --output {params.out_prefix} 2>&1 | tee {log}
            
        # 4. Cleanup
        rm {params.tmp_fna}
        """

# ------------------------------------------------------------------------------
# Module 6: Integration & Reporting
# ------------------------------------------------------------------------------

def get_annot_input(wildcards):
    """Dynamically build input for merge_annotations based on enabled steps"""
    steps = config.get("genome_annotation", {}).get("steps", {})
    inputs = {
        "gff": os.path.join(GA_OUT_DIR, "{genome}/02_bakta/{genome}.gff3").format(genome=wildcards.genome),
        "bakta_tsv": os.path.join(GA_OUT_DIR, "{genome}/02_bakta/{genome}.tsv").format(genome=wildcards.genome)
    }
    
    # Mapping of step name to (arg_name, path_pattern)
    optional_steps = {
        "eggnog": ("eggnog", os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "{genome}.emapper.annotations")),
        "kofamscan": ("kofam", os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "{genome}.kofam.txt")),
        "dbcan": ("dbcan", os.path.join(GA_OUT_DIR, "{genome}", "03_functional", "dbcan")),
        "genomespot": ("genomespot", os.path.join(GA_OUT_DIR, "{genome}", "04_phenotype", "{genome}.genomespot.tsv")),
        "checkm2": ("checkm2", os.path.join(GA_OUT_DIR, "quality", "checkm2_summary.tsv")),
        "gtdbtk": ("gtdbtk", os.path.join(GA_OUT_DIR, "taxonomy", "gtdbtk.summary.tsv")),
        "quast": ("quast", os.path.join(GA_OUT_DIR, "{genome}", "01_qc", "quast"))
    }
    
    for step, (key, pattern) in optional_steps.items():
        if steps.get(step, True):
            inputs[key] = pattern.format(genome=wildcards.genome)
        else:
            pass
            
    return inputs

rule merge_annotations:
    input:
        gff = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.gff3"),
        bakta_tsv = os.path.join(GA_OUT_DIR, "{genome}", "02_bakta", "{genome}.tsv")
    output:
        master = os.path.join(GA_OUT_DIR, "{genome}", "05_integrated", "{genome}_ordered_gene_set.tsv")
    log:
        os.path.join(GA_OUT_DIR, "logs", "merge", "{genome}.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("utils", "utils"),
        params_str = get_annot_input
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 modules/genome_annotation/scripts/merge_annotations.py \\
            --genome {wildcards.genome} \\
            --gff {input.gff} \\
            --bakta {input.bakta_tsv} \\
            {params.params_str} \\
            --output {output.master} 2>&1 | tee {log}
        """

rule generate_genome_summary:
    input:
        checkm2 = os.path.join(GA_OUT_DIR, "quality/checkm2_summary.tsv") if c.get("steps", {}).get("checkm2", True) else [],
        gtdbtk = os.path.join(GA_OUT_DIR, "taxonomy/gtdbtk.summary.tsv") if c.get("steps", {}).get("gtdbtk", True) else [],
        quast_reports = expand(os.path.join(GA_OUT_DIR, "{genome}/01_qc/quast"), genome=GENOMES) if c.get("steps", {}).get("quast", True) else [],
        gs_reports = expand(os.path.join(GA_OUT_DIR, "{genome}/04_phenotype/{genome}.predictions.tsv"), genome=GENOMES) if c.get("steps", {}).get("genomespot", True) else []
    output:
        summary = os.path.join(GA_OUT_DIR, "Batch_Genome_Summary.tsv")
    log:
        os.path.join(GA_OUT_DIR, "logs/summary/batch_summary.log")
    params:
        env = c.get("conda_envs", {}).get("utils", "utils"),
        # Re-construct flags for the script
        checkm2_arg = lambda wildcards, input: f"--checkm2 {input.checkm2}" if hasattr(input, "checkm2") and input.checkm2 else "",
        gtdbtk_arg = lambda wildcards, input: f"--gtdbtk {input.gtdbtk}" if hasattr(input, "gtdbtk") and input.gtdbtk else "",
        quast_arg = lambda wildcards, input: f"--quast-reports {' '.join(input.quast_reports)}" if hasattr(input, "quast_reports") and input.quast_reports else "",
        gs_arg = lambda wildcards, input: f"--genomespot-reports {' '.join(input.gs_reports)}" if hasattr(input, "gs_reports") and input.gs_reports else ""
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 modules/genome_annotation/scripts/summarize_genome_qc_tax.py \
            --genomes {GENOMES} \
            {params.checkm2_arg} \
            {params.gtdbtk_arg} \
            {params.quast_arg} \
            {params.gs_arg} \
            --output {output.summary} 2>&1 | tee {log}
        """

rule run_multiqc:
    input:
        # Wait for all genome master tables
        expand(os.path.join(GA_OUT_DIR, "{genome}", "05_integrated", "{genome}_ordered_gene_set.tsv"), genome=GENOMES),
        # Wait for assembly level results if in hybrid mode
        assembly_gff = os.path.join(GA_ASSEMBLY_OUT, "assembly.gff3") if MODE == "hybrid" else []
    output:
        report = os.path.join(GA_OUT_DIR, "report", "multiqc_report.html")
    log:
        os.path.join(GA_OUT_DIR, "logs", "multiqc", "run.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("bio", "bio"),
        outdir = os.path.join(GA_OUT_DIR, "report")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        conda run -p {params.env} multiqc {GA_OUT_DIR} \\
            -o {params.outdir} \\
            -n multiqc_report 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Module 6: Standalone Plasmid Annotation (via Subsetting)
# ------------------------------------------------------------------------------
if PLASMIDS:
    # We use the same subsetting logic for plasmids to avoid re-running Bakta/eggNOG
    rule subset_assembly_to_plasmids:
        input:
            plasmid_fa = get_plasmid_input,
            gff = os.path.join(GA_ASSEMBLY_OUT, "assembly.gff3"),
            faa = os.path.join(GA_ASSEMBLY_OUT, "assembly.faa"),
            tsv = os.path.join(GA_ASSEMBLY_OUT, "assembly.tsv"),
            eggnog = os.path.join(GA_ASSEMBLY_OUT, "assembly.emapper.annotations"),
            kofam = os.path.join(GA_ASSEMBLY_OUT, "assembly.kofam.txt"),
            dbcan = os.path.join(GA_ASSEMBLY_OUT, "dbcan")
        output:
            gff = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.gff3"),
            faa = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.faa"),
            tsv = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.tsv"),
            eggnog = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.emapper.annotations"),
            kofam = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.kofam.txt"),
            dbcan = directory(os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/dbcan"))
        log: os.path.join(GA_OUT_DIR, "logs/subset/plasmid_{plasmid}.log")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("utils", "utils")
        shell:
            """
            set -euo pipefail
            conda run -p {params.env} python3 modules/genome_annotation/scripts/subset_assembly_results.py \
                --bin_fasta {input.plasmid_fa} \
                --assembly_gff {input.gff} \
                --assembly_faa {input.faa} \
                --assembly_tsv {input.tsv} \
                --assembly_eggnog {input.eggnog} \
                --assembly_kofam {input.kofam} \
                --assembly_dbcan {input.dbcan} \
                --out_gff {output.gff} \
                --out_faa {output.faa} \
                --out_tsv {output.tsv} \
                --out_eggnog {output.eggnog} \
                --out_kofam {output.kofam} \
                --out_dbcan {output.dbcan} 2>&1 | tee {log}
            """

    rule merge_plasmid_annotations:
        input:
            gff = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.gff3"),
            bakta_tsv = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.tsv"),
            eggnog = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}.emapper.annotations")
        output:
            master = os.path.join(GA_OUT_DIR, "plasmids/{plasmid}/{plasmid}_Plasmid_Annotation.tsv")
        log: os.path.join(GA_OUT_DIR, "logs/merge/plasmid_{plasmid}.log")
        threads: 1
        resources:
            runtime = "7d",
            mem_mb = lambda wildcards, threads: threads * 8000
        params:
            env = c.get("conda_envs", {}).get("utils", "utils")
        shell:
            """
            set -euo pipefail
            conda run -p {params.env} python3 modules/genome_annotation/scripts/merge_annotations.py \
                --genome {wildcards.plasmid} --gff {input.gff} --bakta {input.bakta_tsv} \
                --eggnog {input.eggnog} --output {output.master} 2>&1 | tee {log}
            """

rule aggregate_all_annotations:
    input:
        bin_tables = expand(os.path.join(GA_OUT_DIR, "{genome}", "05_integrated", "{genome}_ordered_gene_set.tsv"), genome=GENOMES),
        plasmid_tables = expand(os.path.join(GA_OUT_DIR, "plasmids", "{plasmid}", "{plasmid}_Plasmid_Annotation.tsv"), plasmid=PLASMIDS) if PLASMIDS else []
    output:
        all_genes = os.path.join(GA_OUT_DIR, "All_Genes_Master_Table.tsv")
    log:
        os.path.join(GA_OUT_DIR, "logs", "utils", "aggregate_all.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c.get("conda_envs", {}).get("utils", "utils")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} python3 modules/genome_annotation/scripts/aggregate_all_results.py \\
            --bin-tables {input.bin_tables} \\
            --plasmid-tables {input.plasmid_tables} \\
            --output {output.all_genes} 2>&1 | tee {log}
        """

# Add all targets
TARGETS.extend(expand(os.path.join(GA_OUT_DIR, "{genome}", "05_integrated", "{genome}_ordered_gene_set.tsv"), genome=GENOMES))

if PLASMIDS:
    TARGETS.extend(expand(os.path.join(GA_OUT_DIR, "plasmids", "{plasmid}", "{plasmid}_Plasmid_Annotation.tsv"), plasmid=PLASMIDS))

TARGETS.append(os.path.join(GA_OUT_DIR, "All_Genes_Master_Table.tsv"))
TARGETS.append(os.path.join(GA_OUT_DIR, "Batch_Genome_Summary.tsv"))
TARGETS.append(os.path.join(GA_OUT_DIR, "report", "multiqc_report.html"))

# Add assembly level results to targets if in hybrid mode
if MODE == "hybrid":
    TARGETS.append(os.path.join(GA_ASSEMBLY_OUT, "assembly.gff3"))
    TARGETS.append(os.path.join(GA_ASSEMBLY_OUT, "assembly.emapper.annotations"))

# Optional standalone targets (if someone wants just taxonomy)
if c.get("steps", {}).get("gtdbtk", True):
    TARGETS.append(os.path.join(GA_OUT_DIR, "taxonomy", "gtdbtk.summary.tsv"))
