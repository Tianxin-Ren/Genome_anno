# Unified bacterial/archaeal annotation module.
# First runnable implementation for config_En.

import os

UA_ROOT = config.get("project_root", ".") or "."
UA_ENV_DIR = config.get("conda_env_dir", "")
UA_C = config.get("unified_annotation", {})
UA_STEPS = UA_C.get("steps", {})
# samples.yaml is intentionally allowed to be either a plain mapping or the
# conventional {samples: {...}} document.  The latter is used by config_En.
UA_SAMPLES = config.get("samples", {}) or {}
if isinstance(UA_SAMPLES, dict) and isinstance(UA_SAMPLES.get("samples"), dict):
    UA_SAMPLES = UA_SAMPLES["samples"]

def ua_full_path(path):
    if not path:
        return ""
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(UA_ROOT, path))


UA_RESULTS_ROOT = config.get("results_dir", "results")
# The top-level pipeline result root is the normal control point.  A module
# level results_dir remains supported for backwards compatibility.
UA_RES = ua_full_path(UA_C.get("results_dir", os.path.join(UA_RESULTS_ROOT, "unified_annotation")))
UA_LOGS = ua_full_path(UA_C.get("logs_dir", os.path.join(UA_RES, "logs")))
UA_TMP = ua_full_path(UA_C.get("tmp_dir", os.path.join(UA_RES, "tmp")))


def ua_env(name):
    env_name = config.get("environments", {}).get(name) or UA_C.get("conda_envs", {}).get(name, name)
    if not UA_ENV_DIR:
        return env_name
    if os.path.isabs(env_name):
        return env_name
    return os.path.join(UA_ENV_DIR, env_name)


def ua_sample_fasta(wildcards):
    info = UA_SAMPLES[wildcards.sample]
    return ua_full_path(info.get(UA_C.get("input", {}).get("fasta_key", "fasta")))


def ua_assembly_type(sample):
    return (UA_SAMPLES.get(sample, {}) or {}).get("assembly_type", "MAG")


def ua_pyrodigal_mode(wildcards):
    modes = UA_C.get("cds_prediction", {}).get("mode_by_assembly_type", {})
    return modes.get(ua_assembly_type(wildcards.sample), "meta")


def ua_translation_table(wildcards):
    sample_info = UA_SAMPLES.get(wildcards.sample, {}) or {}
    value = sample_info.get("genetic_code", UA_C.get("cds_prediction", {}).get("translation_table", "auto"))
    if str(value).lower() in {"", "auto", "none", "auto_from_sample"}:
        return "auto"
    return str(value)


def ua_domain_mode(wildcards):
    domain = (UA_SAMPLES.get(wildcards.sample, {}) or {}).get("expected_domain", "auto")
    if str(domain).lower().startswith("arch"):
        return "archaea"
    if str(domain).lower().startswith("bac"):
        return "bacteria"
    return "bacteria"


UA_GENOMES = sorted([
    name for name, info in UA_SAMPLES.items()
    if isinstance(info, dict) and info.get(UA_C.get("input", {}).get("fasta_key", "fasta"))
])

UA_TARGETS = []
for sample in UA_GENOMES:
    if UA_STEPS.get("genome_qc", True):
        UA_TARGETS.extend([
            os.path.join(UA_RES, sample, "quality/quast/report.tsv"),
            os.path.join(UA_RES, sample, "quality/checkm2/quality_report.tsv"),
        ])
    if UA_STEPS.get("taxonomy", True):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "taxonomy/gtdbtk.done"))
    if UA_STEPS.get("cds_prediction", True):
        UA_TARGETS.extend([
            os.path.join(UA_RES, sample, "structural/features.gff3"),
            os.path.join(UA_RES, sample, "structural/cds.fna"),
            os.path.join(UA_RES, sample, "structural/proteins.faa"),
            os.path.join(UA_RES, sample, "structural/gene_prediction_evidence.jsonl"),
        ])
    if UA_STEPS.get("rna_prediction", True):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "structural/rna_features.gff3"))
    if UA_STEPS.get("protein_similarity", True):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "function/mmseqs/swissprot_hits.m8"))
    if UA_STEPS.get("orthology", True):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "function/eggnog/{sample}.emapper.annotations".format(sample=sample)))
    if UA_STEPS.get("ko_profiles", True):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "function/kofam/kofam.tsv"))
    if UA_STEPS.get("protein_domains", True):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "function/interpro/interpro.tsv"))
    if UA_STEPS.get("localization_topology", False):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "function/topology/deeptmhmm.tsv"))
    if UA_STEPS.get("signal_peptide_prediction", False):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "function/localization/signalp.tsv"))
    if UA_STEPS.get("biob_hmm", False):
        UA_TARGETS.append(os.path.join(UA_RES, sample, "function/biob_hmm/biob_hmm.tbl"))
    if UA_STEPS.get("candidate_integration", True):
        UA_TARGETS.extend([
            os.path.join(UA_RES, sample, "function/protein_function_candidates.tsv"),
            os.path.join(UA_RES, sample, "function/protein_function_consensus.jsonl"),
            os.path.join(UA_RES, sample, "function/annotation_conflicts.tsv"),
            os.path.join(UA_RES, sample, "function/unresolved_proteins.faa"),
            os.path.join(UA_RES, sample, "genome/genome_function_inventory.json"),
            os.path.join(UA_RES, sample, "tables/genes.tsv"),
            os.path.join(UA_RES, sample, "tables/evidence.tsv"),
            os.path.join(UA_RES, sample, "tables/modules.tsv"),
            os.path.join(UA_RES, sample, "tables/unknown_clusters.tsv"),
            os.path.join(UA_RES, sample, "report/annotation_report.html"),
        ])

if UA_GENOMES and UA_STEPS.get("candidate_integration", True):
    UA_TARGETS.append(os.path.join(UA_RES, "summary/genomes.tsv"))
    UA_TARGETS.append(os.path.join(UA_RES, "summary/unknown_clusters.tsv"))
    UA_TARGETS.append(os.path.join(UA_RES, "summary/annotation_report.html"))

TARGETS.extend(UA_TARGETS)


rule ua_normalize_fasta:
    input:
        fasta=ua_sample_fasta
    output:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    log:
        os.path.join(UA_LOGS, "normalize/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fasta}) $(dirname {log})
        if [[ "{input.fasta}" == *.gz ]]; then
            gzip -cd {input.fasta} > {output.fasta}
        else
            cp {input.fasta} {output.fasta}
        fi
        echo "normalized {input.fasta}" > {log}
        """


rule ua_quast:
    input:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    output:
        report=os.path.join(UA_RES, "{sample}", "quality/quast/report.tsv")
    log:
        os.path.join(UA_LOGS, "quast/{sample}.log")
    threads: get_threads(4)
    resources:
        runtime="4h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("quality_and_taxonomy", {}).get("assembly_stats", {}).get("env", "quast_5.3.0")),
        outdir=lambda wildcards: os.path.join(UA_RES, wildcards.sample, "quality/quast")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        conda run -p {params.env} quast.py {input.fasta} -o {params.outdir} -t {threads} 2>&1 | tee {log}
        """


rule ua_checkm2:
    input:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    output:
        report=os.path.join(UA_RES, "{sample}", "quality/checkm2/quality_report.tsv")
    log:
        os.path.join(UA_LOGS, "checkm2/{sample}.log")
    threads: get_threads(8)
    resources:
        runtime="12h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("quality_and_taxonomy", {}).get("completeness_contamination", {}).get("env", "checkm2_1.0.2")),
        db=UA_C.get("quality_and_taxonomy", {}).get("completeness_contamination", {}).get("database", ""),
        indir=lambda wildcards: os.path.join(UA_TMP, "checkm2_input", wildcards.sample),
        outdir=lambda wildcards: os.path.join(UA_RES, wildcards.sample, "quality/checkm2")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.indir} {params.outdir} $(dirname {log})
        ln -sf $(readlink -f {input.fasta}) {params.indir}/{wildcards.sample}.fna
        conda run -p {params.env} checkm2 predict \
            --threads {threads} \
            --input {params.indir} \
            --output-directory {params.outdir} \
            --database_path {params.db} \
            -x fna 2>&1 | tee {log}
        """


rule ua_gtdbtk:
    input:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    output:
        done=os.path.join(UA_RES, "{sample}", "taxonomy/gtdbtk.done")
    log:
        os.path.join(UA_LOGS, "gtdbtk/{sample}.log")
    threads: get_threads(16)
    resources:
        runtime="24h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("quality_and_taxonomy", {}).get("taxonomy", {}).get("env", config.get("environments", {}).get("gtdbtk", "gtdbtk_2.7.1"))),
        db=UA_C.get("quality_and_taxonomy", {}).get("taxonomy", {}).get("database", ""),
        indir=lambda wildcards: os.path.join(UA_TMP, "gtdbtk_input", wildcards.sample),
        outdir=lambda wildcards: os.path.join(UA_RES, wildcards.sample, "taxonomy/gtdbtk")
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.indir} {params.outdir} $(dirname {output.done}) $(dirname {log})
        ln -sf $(readlink -f {input.fasta}) {params.indir}/{wildcards.sample}.fna
        GTDBTK_DATA_PATH={params.db} conda run -p {params.env} gtdbtk classify_wf \
            --genome_dir {params.indir} \
            --out_dir {params.outdir} \
            --extension fna \
            --cpus {threads} 2>&1 | tee {log}
        touch {output.done}
        """


rule ua_pyrodigal:
    input:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    output:
        gff=os.path.join(UA_RES, "{sample}", "structural/features.gff3"),
        cds=os.path.join(UA_RES, "{sample}", "structural/cds.fna"),
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa"),
        evidence=os.path.join(UA_RES, "{sample}", "structural/gene_prediction_evidence.jsonl")
    log:
        os.path.join(UA_LOGS, "pyrodigal/{sample}.log")
    threads: 1
    resources:
        runtime="4h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("cds_prediction", {}).get("env", "pyrodigal_3.7.1")),
        mode=ua_pyrodigal_mode,
        translation_table=ua_translation_table
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gff}) $(dirname {log})
        GENETIC_CODE_ARG=""
        if [ "{params.translation_table}" != "auto" ]; then
            GENETIC_CODE_ARG="-g {params.translation_table}"
        fi
        conda run -p {params.env} pyrodigal \
            -i {input.fasta} \
            -p {params.mode} \
            $GENETIC_CODE_ARG \
            -f gff \
            -o {output.gff} \
            -d {output.cds} \
            -a {output.faa} 2>&1 | tee {log}
        sed -i 's/\*//g' {output.faa}
        python scripts/unified_annotation_integrate.py gene-evidence \
            --sample {wildcards.sample} \
            --gff {output.gff} \
            --faa {output.faa} \
            --mode {params.mode} \
            --translation-table {params.translation_table} \
            --output {output.evidence}
        """


rule ua_trnascan:
    input:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    output:
        tbl=os.path.join(UA_RES, "{sample}", "rna/trnascan/trnascan.tbl")
    log:
        os.path.join(UA_LOGS, "trnascan/{sample}.log")
    threads: get_threads(4)
    resources:
        runtime="4h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("rna_prediction", {}).get("trna", {}).get("env", "trnascan-se_2.0.12")),
        mode=lambda wildcards: (
            UA_C.get("rna_prediction", {}).get("trna", {}).get("archaeal_mode", "-A")
            if ua_domain_mode(wildcards) == "archaea"
            else UA_C.get("rna_prediction", {}).get("trna", {}).get("bacterial_mode", "-B")
        )
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tbl}) $(dirname {log})
        conda run -p {params.env} tRNAscan-SE {params.mode} --thread {threads} -o {output.tbl} {input.fasta} 2>&1 | tee {log}
        """


rule ua_rfam:
    input:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    output:
        tbl=os.path.join(UA_RES, "{sample}", "rna/rfam/rfam.tbl")
    log:
        os.path.join(UA_LOGS, "rfam/{sample}.log")
    threads: get_threads(8)
    resources:
        runtime="12h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("rna_prediction", {}).get("rfam", {}).get("env", "infernal_1.1.5")),
        db=UA_C.get("rna_prediction", {}).get("rfam", {}).get("database", ""),
        clanin=lambda wildcards: os.path.join(os.path.dirname(UA_C.get("rna_prediction", {}).get("rfam", {}).get("database", "")), "Rfam.clanin")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tbl}) $(dirname {log})
        CLANIN_ARG=""
        if [ -s "{params.clanin}" ]; then
            CLANIN_ARG="--fmt 2 --clanin {params.clanin}"
        fi
        if conda run -p {params.env} cmscan --cpu {threads} --cut_ga --rfam --nohmmonly --tblout {output.tbl} $CLANIN_ARG {params.db} {input.fasta} > {log} 2>&1; then
            :
        else
            status=$?
            {{
                echo "# cmscan failed with exit status ${{status}}; preserving empty Rfam table so optional RNA evidence does not block the first workflow run."
                echo "# See log: {log}"
            }} > {output.tbl}
            exit 0
        fi
        if [ ! -s {output.tbl} ]; then
            echo "# cmscan completed but produced no table rows" > {output.tbl}
        fi
        """


rule ua_aragorn:
    input:
        fasta=os.path.join(UA_RES, "{sample}", "input/genome.fna")
    output:
        txt=os.path.join(UA_RES, "{sample}", "rna/aragorn/aragorn.txt")
    log:
        os.path.join(UA_LOGS, "aragorn/{sample}.log")
    threads: 1
    resources:
        runtime="2h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("rna_prediction", {}).get("tmrna", {}).get("env", "aragorn_1.2.41"))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.txt}) $(dirname {log})
        if conda run -p {params.env} aragorn -w -m -t {input.fasta} > {output.txt} 2> {log}; then
            :
        else
            status=$?
            {{
                echo "# ARAGORN failed with exit status ${{status}}"
                echo "# See log: {log}"
            }} > {output.txt}
            exit 0
        fi
        if [ ! -s {output.txt} ]; then
            echo "# ARAGORN completed but produced no rows" > {output.txt}
        fi
        """


rule ua_rna_gff:
    input:
        trna=os.path.join(UA_RES, "{sample}", "rna/trnascan/trnascan.tbl"),
        rfam=os.path.join(UA_RES, "{sample}", "rna/rfam/rfam.tbl"),
        aragorn=os.path.join(UA_RES, "{sample}", "rna/aragorn/aragorn.txt")
    output:
        gff=os.path.join(UA_RES, "{sample}", "structural/rna_features.gff3")
    log:
        os.path.join(UA_LOGS, "rna_gff/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gff}) $(dirname {log})
        python scripts/unified_annotation_integrate.py rna-gff \
            --trna {input.trna} \
            --rfam {input.rfam} \
            --aragorn {input.aragorn} \
            --output {output.gff}
        echo "Wrote parsed RNA feature GFF: {output.gff}" > {log}
        """


rule ua_mmseqs_swissprot:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa")
    output:
        hits=os.path.join(UA_RES, "{sample}", "function/mmseqs/swissprot_hits.m8")
    log:
        os.path.join(UA_LOGS, "mmseqs/{sample}.log")
    threads: get_threads(8)
    resources:
        runtime="12h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("protein_function", {}).get("similarity", {}).get("env", "mmseqs2_17.b804f")),
        db=UA_C.get("protein_function", {}).get("similarity", {}).get("databases", {}).get("swissprot", ""),
        tmp=lambda wildcards: os.path.join(UA_TMP, "mmseqs", wildcards.sample),
        min_identity=UA_C.get("protein_function", {}).get("similarity", {}).get("min_identity", 0.30),
        min_coverage=UA_C.get("protein_function", {}).get("similarity", {}).get("min_query_coverage", 0.70),
        max_hits=UA_C.get("protein_function", {}).get("similarity", {}).get("max_hits_per_protein", 25)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.hits}) {params.tmp} $(dirname {log})
        conda run -p {params.env} mmseqs easy-search \
            {input.faa} {params.db} {output.hits} {params.tmp} \
            --threads {threads} \
            --min-seq-id {params.min_identity} \
            -c {params.min_coverage} \
            --max-seqs {params.max_hits} 2>&1 | tee {log}
        """


rule ua_mmseqs_disabled:
    output:
        hits=os.path.join(UA_RES, "{sample}", "function/mmseqs/swissprot.disabled.m8")
    log:
        os.path.join(UA_LOGS, "mmseqs_disabled/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.hits}) $(dirname {log})
        echo "# MMseqs/Swiss-Prot disabled by config" > {output.hits}
        echo "MMseqs/Swiss-Prot disabled by config" > {log}
        """


rule ua_eggnog:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa")
    output:
        anno=os.path.join(UA_RES, "{sample}", "function/eggnog/{sample}.emapper.annotations")
    log:
        os.path.join(UA_LOGS, "eggnog/{sample}.log")
    threads: get_threads(8)
    resources:
        runtime="24h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("protein_function", {}).get("orthology", {}).get("env", "eggnog-mapper")),
        db=UA_C.get("protein_function", {}).get("orthology", {}).get("database", ""),
        mode=UA_C.get("protein_function", {}).get("orthology", {}).get("search_mode", "diamond"),
        outdir=lambda wildcards: os.path.join(UA_RES, wildcards.sample, "function/eggnog"),
        prefix=lambda wildcards: wildcards.sample,
        scratch=lambda wildcards: os.path.join(UA_TMP, "eggnog", wildcards.sample)
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} {params.scratch} $(dirname {log})
        if conda run -p {params.env} emapper.py \
            -m {params.mode} \
            -i {input.faa} \
            --output {params.prefix} \
            --output_dir {params.outdir} \
            --cpu {threads} \
            --scratch_dir {params.scratch} \
            --data_dir {params.db} \
            --override > {log} 2>&1; then
            :
        else
            status=$?
            {{
                echo "# eggNOG-mapper failed with exit status ${{status}}"
                echo "# See log: {log}"
            }} > {output.anno}
            exit 0
        fi
        if [ ! -s {output.anno} ]; then
            echo "# eggNOG-mapper completed but did not create annotations" > {output.anno}
        fi
        """


rule ua_eggnog_disabled:
    output:
        anno=os.path.join(UA_RES, "{sample}", "function/eggnog/{sample}.disabled.annotations")
    log:
        os.path.join(UA_LOGS, "eggnog_disabled/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.anno}) $(dirname {log})
        echo "# EggNOG disabled by config" > {output.anno}
        echo "EggNOG disabled by config" > {log}
        """


rule ua_kofam:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa")
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/kofam/kofam.tsv")
    log:
        os.path.join(UA_LOGS, "kofam/{sample}.log")
    threads: get_threads(8)
    resources:
        runtime="12h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("protein_function", {}).get("ko_profiles", {}).get("env", "kofamscan")),
        profiles=UA_C.get("protein_function", {}).get("ko_profiles", {}).get("profiles", ""),
        ko_list=UA_C.get("protein_function", {}).get("ko_profiles", {}).get("ko_list", ""),
        tmp=lambda wildcards: os.path.join(UA_TMP, "kofam", wildcards.sample)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) {params.tmp} $(dirname {log})
        conda run -p {params.env} exec_annotation \
            -o {output.tsv} \
            --cpu {threads} \
            --tmp {params.tmp} \
            -p {params.profiles} \
            -k {params.ko_list} \
            {input.faa} 2>&1 | tee {log}
        rm -rf {params.tmp}
        """


rule ua_kofam_disabled:
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/kofam/kofam.disabled.tsv")
    log:
        os.path.join(UA_LOGS, "kofam_disabled/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        echo "# KOfam disabled by config" > {output.tsv}
        echo "KOfam disabled by config" > {log}
        """


rule ua_interpro:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa")
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/interpro/interpro.tsv")
    log:
        os.path.join(UA_LOGS, "interpro/{sample}.log")
    threads: get_threads(16)
    resources:
        runtime="48h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        script=os.path.join(ua_env(UA_C.get("protein_function", {}).get("domains", {}).get("env", "interproscan_5.77-108.0")), "interproscan.sh"),
        java_env=ua_env("interproscan_java11"),
        apps=lambda wildcards: ",".join(UA_C.get("protein_function", {}).get("domains", {}).get("applications", ["Pfam"]))
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        if conda run -p {params.java_env} {params.script} \
            -i {input.faa} \
            -f TSV \
            -appl {params.apps} \
            -cpu {threads} \
            -o {output.tsv} \
            -dp > {log} 2>&1; then
            :
        else
            status=$?
            {{
                echo "# InterProScan failed with exit status ${{status}}"
                echo "# See log: {log}"
            }} > {output.tsv}
            exit 0
        fi
        if [ ! -s {output.tsv} ]; then
            echo "# InterProScan completed but produced no TSV rows" > {output.tsv}
        fi
        """


rule ua_interpro_disabled:
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/interpro/interpro.disabled.tsv")
    log:
        os.path.join(UA_LOGS, "interpro_disabled/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        echo "# InterProScan disabled by config" > {output.tsv}
        echo "InterProScan disabled by config" > {log}
        """


rule ua_signalp_disabled:
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/localization/signalp.disabled.tsv")
    log:
        os.path.join(UA_LOGS, "signalp_disabled/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        printf "feature_id\tprediction\treason\n" > {output.tsv}
        echo "SignalP disabled by config" > {log}
        """


rule ua_signalp:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa")
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/localization/signalp.tsv")
    log:
        os.path.join(UA_LOGS, "signalp/{sample}.log")
    threads: get_threads(4)
    resources:
        runtime="8h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("protein_function", {}).get("signal_peptide", {}).get("env", "signalp6_6.0h")),
        outdir=lambda wildcards: os.path.join(UA_RES, wildcards.sample, "function/localization/signalp_raw")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) {params.outdir} $(dirname {log})
        if conda run -p {params.env} signalp6 --fastafile {input.faa} --output_dir {params.outdir} --format none > {log} 2>&1; then
            FOUND="$(find {params.outdir} -type f -name "prediction_results.txt" -print -quit)"
            if [ -n "$FOUND" ]; then
                cp "$FOUND" {output.tsv}
            else
                printf "feature_id\tprediction\treason\n" > {output.tsv}
            fi
        else
            status=$?
            {{
                echo "feature_id	prediction	reason"
                echo "NA	unresolved	SignalP failed with exit status ${{status}}; see {log}"
            }} > {output.tsv}
            exit 0
        fi
        """


rule ua_deeptmhmm:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa")
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/topology/deeptmhmm.tsv")
    log:
        os.path.join(UA_LOGS, "deeptmhmm/{sample}.log")
    threads: get_threads(4)
    resources:
        runtime="12h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("protein_function", {}).get("topology", {}).get("env", "deeptmhmm_1.0")),
        outdir=lambda wildcards: os.path.join(UA_RES, wildcards.sample, "function/topology/deeptmhmm_raw")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) {params.outdir} $(dirname {log})
        if conda run -p {params.env} biolib run DTU/DeepTMHMM --fasta {input.faa} --outdir {params.outdir} > {log} 2>&1; then
            if [ -s {params.outdir}/TMRs.gff3 ]; then
                awk 'BEGIN{{FS=OFS="\t"; print "feature_id","prediction","source"}} !/^#/ {{print $1,$3,"DeepTMHMM"}}' {params.outdir}/TMRs.gff3 > {output.tsv}
            else
                printf "feature_id\tprediction\tsource\n" > {output.tsv}
            fi
        else
            status=$?
            {{
                echo "feature_id	prediction	reason"
                echo "NA	unresolved	DeepTMHMM failed with exit status ${{status}}; see {log}"
            }} > {output.tsv}
            exit 0
        fi
        """


rule ua_deeptmhmm_disabled:
    output:
        tsv=os.path.join(UA_RES, "{sample}", "function/topology/deeptmhmm.disabled.tsv")
    log:
        os.path.join(UA_LOGS, "deeptmhmm_disabled/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        printf "feature_id\tprediction\treason\n" > {output.tsv}
        echo "DeepTMHMM disabled by config" > {log}
        """


rule ua_biob_hmm_disabled:
    output:
        tbl=os.path.join(UA_RES, "{sample}", "function/biob_hmm/biob_hmm.disabled.tbl")
    log:
        os.path.join(UA_LOGS, "biob_hmm_disabled/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tbl}) $(dirname {log})
        echo "# BioB HMM disabled by config" > {output.tbl}
        echo "BioB HMM disabled by config" > {log}
        """


rule ua_biob_hmm:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa")
    output:
        tbl=os.path.join(UA_RES, "{sample}", "function/biob_hmm/biob_hmm.tbl")
    log:
        os.path.join(UA_LOGS, "biob_hmm/{sample}.log")
    threads: get_threads(8)
    resources:
        runtime="12h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env=ua_env(UA_C.get("protein_function", {}).get("biob_hmm", {}).get("env", "hmmer_3.4")),
        db=UA_C.get("protein_function", {}).get("biob_hmm", {}).get("hmm_db", "")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tbl}) $(dirname {log})
        if [ ! -s "{params.db}" ]; then
            {{
                echo "# BioB HMM database is missing or empty: {params.db}"
                echo "# This is an attempted formal BioB HMM channel, not a disabled channel."
            }} > {log}
            {{
                echo "# BioB HMM attempted but database is missing or empty: {params.db}"
            }} > {output.tbl}
            exit 0
        fi
        if conda run -p {params.env} hmmsearch --cpu {threads} --tblout {output.tbl} {params.db} {input.faa} > {log} 2>&1; then
            if [ ! -s {output.tbl} ]; then
                echo "# BioB HMM completed but produced no table rows" > {output.tbl}
            fi
        else
            status=$?
            {{
                echo "# BioB HMM attempted but hmmsearch failed with exit status ${{status}}"
                echo "# See log: {log}"
            }} > {output.tbl}
            exit 0
        fi
        """


def ua_signalp_input(wildcards):
    if UA_STEPS.get("signal_peptide_prediction", False):
        return os.path.join(UA_RES, wildcards.sample, "function/localization/signalp.tsv")
    return os.path.join(UA_RES, wildcards.sample, "function/localization/signalp.disabled.tsv")


def ua_biob_hmm_input(wildcards):
    if UA_STEPS.get("biob_hmm", False):
        return os.path.join(UA_RES, wildcards.sample, "function/biob_hmm/biob_hmm.tbl")
    return os.path.join(UA_RES, wildcards.sample, "function/biob_hmm/biob_hmm.disabled.tbl")


def ua_topology_input(wildcards):
    if UA_STEPS.get("localization_topology", False):
        return os.path.join(UA_RES, wildcards.sample, "function/topology/deeptmhmm.tsv")
    return os.path.join(UA_RES, wildcards.sample, "function/topology/deeptmhmm.disabled.tsv")


def ua_mmseqs_input(wildcards):
    if UA_STEPS.get("protein_similarity", True):
        return os.path.join(UA_RES, wildcards.sample, "function/mmseqs/swissprot_hits.m8")
    return os.path.join(UA_RES, wildcards.sample, "function/mmseqs/swissprot.disabled.m8")


def ua_eggnog_input(wildcards):
    if UA_STEPS.get("orthology", True):
        return os.path.join(UA_RES, wildcards.sample, f"function/eggnog/{wildcards.sample}.emapper.annotations")
    return os.path.join(UA_RES, wildcards.sample, f"function/eggnog/{wildcards.sample}.disabled.annotations")


def ua_kofam_input(wildcards):
    if UA_STEPS.get("ko_profiles", True):
        return os.path.join(UA_RES, wildcards.sample, "function/kofam/kofam.tsv")
    return os.path.join(UA_RES, wildcards.sample, "function/kofam/kofam.disabled.tsv")


def ua_interpro_input(wildcards):
    if UA_STEPS.get("protein_domains", True):
        return os.path.join(UA_RES, wildcards.sample, "function/interpro/interpro.tsv")
    return os.path.join(UA_RES, wildcards.sample, "function/interpro/interpro.disabled.tsv")


rule ua_consensus_integrate:
    input:
        faa=os.path.join(UA_RES, "{sample}", "structural/proteins.faa"),
        gff=os.path.join(UA_RES, "{sample}", "structural/features.gff3"),
        mmseqs=ua_mmseqs_input,
        eggnog=ua_eggnog_input,
        kofam=ua_kofam_input,
        interpro=ua_interpro_input,
        signalp=ua_signalp_input,
        topology=ua_topology_input,
        biob_hmm=ua_biob_hmm_input
    output:
        candidates=os.path.join(UA_RES, "{sample}", "function/protein_function_candidates.tsv"),
        consensus=os.path.join(UA_RES, "{sample}", "function/protein_function_consensus.jsonl"),
        conflicts=os.path.join(UA_RES, "{sample}", "function/annotation_conflicts.tsv"),
        unresolved=os.path.join(UA_RES, "{sample}", "function/unresolved_proteins.faa"),
        inventory=os.path.join(UA_RES, "{sample}", "genome/genome_function_inventory.json"),
        genes=os.path.join(UA_RES, "{sample}", "tables/genes.tsv"),
        evidence_table=os.path.join(UA_RES, "{sample}", "tables/evidence.tsv"),
        modules=os.path.join(UA_RES, "{sample}", "tables/modules.tsv"),
        unknown_clusters=os.path.join(UA_RES, "{sample}", "tables/unknown_clusters.tsv"),
        report=os.path.join(UA_RES, "{sample}", "report/annotation_report.html")
    log:
        os.path.join(UA_LOGS, "consensus/{sample}.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        fasta = lambda wildcards: os.path.join(UA_RES, wildcards.sample, "input/genome.fna")
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.candidates}) $(dirname {output.inventory}) $(dirname {log})
        python scripts/unified_annotation_integrate.py consensus \
            --sample {wildcards.sample} \
            --fasta {params.fasta} \
            --faa {input.faa} \
            --gff {input.gff} \
            --mmseqs {input.mmseqs} \
            --eggnog {input.eggnog} \
            --kofam {input.kofam} \
            --interpro {input.interpro} \
            --signalp {input.signalp} \
            --topology {input.topology} \
            --biob-hmm {input.biob_hmm} \
            --candidates {output.candidates} \
            --consensus {output.consensus} \
            --conflicts {output.conflicts} \
            --unresolved {output.unresolved} \
            --inventory {output.inventory} \
            --genes {output.genes} \
            --evidence-table {output.evidence_table} \
            --modules {output.modules} \
            --unknown-clusters {output.unknown_clusters} \
            --report {output.report} > {log} 2>&1
        echo "Consensus integration complete for {wildcards.sample}" >> {log}
        """


rule ua_cross_genome_summary:
    input:
        genes=expand(os.path.join(UA_RES, "{sample}", "tables/genes.tsv"), sample=UA_GENOMES),
        inventories=expand(os.path.join(UA_RES, "{sample}", "genome/genome_function_inventory.json"), sample=UA_GENOMES)
    output:
        genomes=os.path.join(UA_RES, "summary/genomes.tsv"),
        unknown=os.path.join(UA_RES, "summary/unknown_clusters.tsv"),
        report=os.path.join(UA_RES, "summary/annotation_report.html")
    log:
        os.path.join(UA_LOGS, "summary/cross_genome.log")
    threads: 1
    resources:
        runtime="1h",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        root = UA_RES
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.genomes}) $(dirname {log})
        python scripts/unified_annotation_integrate.py cross-genome \
            --root {params.root} --genomes {output.genomes} --unknown {output.unknown} \
            --report {output.report} > {log} 2>&1
        """
