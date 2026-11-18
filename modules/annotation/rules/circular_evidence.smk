AN_CIRCULAR_EVIDENCE_TARGETS = [
    os.path.join(AN_RES, sample, "circular_evidence", contig, "manifest.json")
    for sample, contigs in AN_C.get("circular_evidence", {}).get("targets", {}).items()
    for contig in contigs
]


rule annotation_circular_evidence_targets:
    input:
        AN_CIRCULAR_EVIDENCE_TARGETS


rule annotation_circular_evidence_self_blast:
    input:
        contig = os.path.join(AN_RES, "{genome}/circular_bins_filtered/{bin}.fa")
    output:
        tsv = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/repeat_self_blast.tsv")
    log:
        os.path.join(AN_LOGS, "circular_evidence/{genome}_{bin}_self_blast.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit")
    shell:
        """
        {{
        set -euo pipefail
        mkdir -p $(dirname {output.tsv:q})
        conda run -p {params.env:q} blastn \
            -query {input.contig:q} -subject {input.contig:q} \
            -outfmt 6 -out {output.tsv:q}
        }} 2>&1 | tee {log:q}
        """


rule annotation_circular_evidence_export:
    input:
        contig = os.path.join(AN_RES, "{genome}/circular_bins_filtered/{bin}.fa"),
        long_bam = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam"),
        long_bai = os.path.join(AN_RES, "{genome}/circular_audit/whole_assembly.bam.bai"),
        junction_bam = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_junction.bam"),
        junction_bai = os.path.join(AN_RES, "{genome}/circular_audit/{bin}/{bin}_junction.bam.bai"),
        short_bam = os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam"),
        short_bai = os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam.bai"),
        rotated_bam = os.path.join(AN_RES, "{genome}/circular_coverm/rotated.raw.bam"),
        rotated_bai = os.path.join(AN_RES, "{genome}/circular_coverm/rotated.raw.bam.bai"),
        coords = os.path.join(AN_RES, "{genome}/circular_coverm/rotated_junction_coords.tsv"),
        self_blast = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/repeat_self_blast.tsv"),
        script = os.path.join(WORKFLOW_ROOT, "scripts/export_circular_evidence.py")
    output:
        manifest = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/manifest.json"),
        coverage = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/genome_coverage.tsv.gz"),
        depth = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/junction_depth.tsv.gz"),
        reads = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/junction_reads.tsv.gz"),
        summary = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/evidence_summary.tsv"),
        gaps = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/zero_coverage_regions.tsv"),
        clipping = os.path.join(AN_RES, "{genome}/circular_evidence/{bin}/long_clipping.tsv")
    log:
        os.path.join(AN_LOGS, "circular_evidence/{genome}_{bin}_export.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        out_dir = lambda wildcards: os.path.join(AN_RES, wildcards.genome, "circular_evidence", wildcards.bin),
        long_anchor = AN_C.get("circular_audit", {}).get("min_junction_flank_bp", 500),
        short_anchor = AN_C.get("circular_illumina", {}).get("junction", {}).get("min_spanning_anchor_bp", 30),
        min_mapq = AN_C.get("circular_illumina", {}).get("junction", {}).get("min_mapq", 20),
        short_context = lambda wildcards: config.get("samples", {}).get(wildcards.genome, {}).get("short_read", {}).get("sampling_context", "unspecified")
    wildcard_constraints:
        genome = "[^/]+",
        bin = "[^/]+"
    shell:
        """
        {{
        set -euo pipefail
        conda run -p {params.env:q} python {input.script:q} \
            --contig {wildcards.bin:q} --contig-fa {input.contig:q} \
            --long-bam {input.long_bam:q} --long-junction-bam {input.junction_bam:q} \
            --short-bam {input.short_bam:q} --short-rotated-bam {input.rotated_bam:q} \
            --rotation-coords {input.coords:q} --self-blast {input.self_blast:q} \
            --out-dir {params.out_dir:q} --min-mapq {params.min_mapq} \
            --long-anchor {params.long_anchor} --short-anchor {params.short_anchor} \
            --short-read-context {params.short_context:q}
        }} 2>&1 | tee {log:q}
        """
