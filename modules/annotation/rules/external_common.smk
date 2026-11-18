# Common rules and helper functions for external genome annotation
import os
import re

_batch_groups = sorted(set(
    info.get("batch_group") for info in EXT_GENOMES.values()
    if isinstance(info, dict) and "batch_group" in info
))


def get_ext_batch_summary_paths(wildcards):
    bg = wildcards.batch_group
    # Trigger resolution checkpoint
    checkpoints.ext_annotation_batch_hmm_resolve_agg.get(batch_group=bg)
    
    # Retrieve all MAG names belonging to this batch group
    mag_names = [
        name for name, info in EXT_GENOMES.items()
        if isinstance(info, dict) and info.get("batch_group") == bg
    ]
    
    paths = []
    for name in sorted(mag_names):
        paths.append(os.path.join(EXT_RES, bg, name, f"{name}_functional_summary.tsv"))
    return paths


def get_ext_inherent_name(bin_key):
    """Extracts the filename without extension to use as an inherent identifier."""
    path = EXT_GENOMES[bin_key] if isinstance(EXT_GENOMES[bin_key], str) else EXT_GENOMES[bin_key].get("path")
    base = os.path.basename(path)
    # Remove common extensions
    name = re.sub(r'\.(fna|fasta|fa|gz|genomic)$', '', base)
    return name

use rule annotation_eggnog as ext_annotation_eggnog with:
    input:
        faa = get_anno_faa_input
    output:
        anno = os.path.join(EXT_RES, "{bin}/eggnog/{bin}.emapper.annotations")
    log: os.path.join(AN_LOGS, "external_eggnog/{bin}.log")
    params:
        env = an_get_env("eggnog"),
        db = AN_C.get("eggnog_db", ""),
        prefix = lambda wildcards: os.path.join(EXT_RES, wildcards.bin, "eggnog", wildcards.bin)

use rule annotation_kofamscan as ext_annotation_kofamscan with:
    input:
        faa = get_anno_faa_input
    output:
        tsv = os.path.join(EXT_RES, "{bin}/kofamscan/{bin}_kofam.tsv")
    log: os.path.join(AN_LOGS, "external_kofamscan/{bin}.log")
    params:
        env = an_get_env("kofamscan"),
        db = AN_C.get("kofam_db", ""),
        tmpdir = os.path.join(EXT_RES, "{bin}/kofamscan/kofam_tmp")

rule ext_annotation_merge_batch:
    input:
        summaries = get_ext_batch_summary_paths,
        sentinel = os.path.join(EXT_RES, "{batch_group}_batch_hmm_completed.txt")
    output:
        master = os.path.join(EXT_RES, "{batch_group}_functional_summary.tsv")
    log: os.path.join(AN_LOGS, "external_merge_batch_{batch_group}.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    shell:
        """
        {{
        set -euo pipefail
        first_file=""
        for f in {input.summaries}; do
            first_file="$f"
            break
        done
        if [ -n "$first_file" ] && [ -f "$first_file" ]; then
            head -n 1 "$first_file" > {output.master}
            for f in {input.summaries}; do
                if [ -f "$f" ]; then
                    tail -n +2 "$f" >> {output.master}
                fi
            done
        else
            touch {output.master}
        fi
        }} 2>&1 | tee {log}
        """

rule ext_annotation_cdhit:
    input:
        faa = os.path.join(EXT_RES, "{bin}/{bin}.faa")
    output:
        faa = os.path.join(EXT_RES, "{bin}/cdhit/{bin}_proteomics_database_c{c}.faa"),
        clstr = os.path.join(EXT_RES, "{bin}/cdhit/{bin}_proteomics_database_c{c}.faa.clstr")
    log:
        os.path.join(AN_LOGS, "external_cdhit/{bin}_c{c}.log")
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
        conda run -p {params.env} cd-hit \
            -i {input.faa} \
            -o {output.faa} \
            -c {params.c} \
            -d {params.d} \
            -T {threads} \
            -M {resources.mem_mb} \
            -n {params.n}
        }} 2>&1 | tee {log}
        """
