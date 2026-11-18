# Circular Contigs Illumina Short-Read Orthogonal Abundance & Verification Rules
# ------------------------------------------------------------------------------

def get_sample_short_reads_dict(wildcards):
    sample = get_sample_name_from_genome(wildcards.genome)
    s_cfg = config.get("samples", {}).get(sample, {})
    r1, r2 = None, None
    if "short_read" in s_cfg:
        sr = s_cfg["short_read"]
        if isinstance(sr, dict):
            r1 = sr.get("r1")
            r2 = sr.get("r2")
        elif isinstance(sr, (list, tuple)) and len(sr) >= 2:
            r1 = sr[0]
            r2 = sr[1]
    if not r1 and "short" in s_cfg and isinstance(s_cfg["short"], dict):
        r1 = s_cfg["short"].get("r1")
        r2 = s_cfg["short"].get("r2")
    if not r1 and "r1" in s_cfg and "r2" in s_cfg:
        r1 = s_cfg["r1"]
        r2 = s_cfg["r2"]

    if r1 and r2:
        return {"r1": an_full_path(r1), "r2": an_full_path(r2)}
    return {"r1": "", "r2": ""}

def get_sample_short_r1(wildcards):
    return get_sample_short_reads_dict(wildcards).get("r1", "")

def get_sample_short_r2(wildcards):
    return get_sample_short_reads_dict(wildcards).get("r2", "")

def sample_has_short_reads(wildcards):
    sr = get_sample_short_reads_dict(wildcards)
    return bool(sr.get("r1") and sr.get("r2"))


# Track A: CoverM Genome Mode Abundance Across Non-Redundant Circular Genomes
# ------------------------------------------------------------------------------
rule annotation_circular_genome_abundance:
    input:
        bins_dir = os.path.join(AN_RES, "{genome}/circular_bins_filtered"),
        r1 = get_sample_short_r1,
        r2 = get_sample_short_r2
    output:
        tsv = os.path.join(AN_RES, "{genome}/circular_coverm/circular_genome_abundance.tsv")
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_genome_abundance.log")
    threads: get_threads(AN_C.get("threads", {}).get("coverm", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("coverm"),
        mapper = lambda wildcards: AN_C.get("circular_illumina", {}).get("mapping", {}).get("mapper", "minimap2-sr"),
        min_pid = lambda wildcards: AN_C.get("circular_illumina", {}).get("standard", {}).get("min_identity", 95),
        min_aln_pct = lambda wildcards: AN_C.get("circular_illumina", {}).get("standard", {}).get("min_aligned_percent", 75),
        end_exclusion = lambda wildcards: AN_C.get("circular_illumina", {}).get("coverm", {}).get("contig_end_exclusion", 75),
        min_cov_frac = lambda wildcards: AN_C.get("circular_illumina", {}).get("abundance", {}).get("min_covered_fraction", 0)
    shell:
        r"""
        {{
        set -euo pipefail
        export PATH="{params.env}/bin:$PATH"
        mkdir -p $(dirname {output.tsv})

        shopt -s nullglob
        fa_files=({input.bins_dir}/*.fa)
        if [ ${{#fa_files[@]}} -eq 0 ]; then
            echo -e "Genome\tRelative Abundance (%)\tMean\tTrimmed Mean\tCovered Fraction\tRead Count" > {output.tsv}
            exit 0
        fi

        coverm genome \
            --coupled {input.r1} {input.r2} \
            --genome-fasta-directory {input.bins_dir} \
            -x fa \
            -p {params.mapper} \
            --min-read-percent-identity {params.min_pid} \
            --min-read-aligned-percent {params.min_aln_pct} \
            --contig-end-exclusion {params.end_exclusion} \
            --min-covered-fraction {params.min_cov_frac} \
            --exclude-supplementary \
            -m relative_abundance mean trimmed_mean covered_fraction count \
            -t {threads} > {output.tsv}
        }} 2>&1 | tee {log}
        """


# Track B: Raw Competitive Mapping on Whole Assembly (Single Mapping Job)
# ------------------------------------------------------------------------------
rule annotation_circular_original_mapping:
    input:
        assembly = lambda wildcards: an_get_asm_path(get_sample_name_from_genome(wildcards.genome)),
        r1 = get_sample_short_r1,
        r2 = get_sample_short_r2
    output:
        bam = temp(os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam")),
        bai = temp(os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam.bai"))
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_original_mapping.log")
    threads: get_threads(AN_C.get("threads", {}).get("minimap2", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        mapping_env = an_get_env("mapping")
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/minimap2" ]; then
            for cand in "{params.env}" \
                        "{params.mapping_env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" \
                        "$HOME/.conda/envs/bio"; do
                if [ -x "$cand/bin/minimap2" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        for extra in "{params.mapping_env}" "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" "$HOME/.conda/envs/bio"; do
            if [ -x "$extra/bin/samtools" ]; then
                export PATH="$extra/bin:$PATH"
            fi
        done
        mkdir -p $(dirname {output.bam})
        minimap2 -ax sr -t {threads} -p 1 {input.assembly} {input.r1} {input.r2} | \
            samtools view -u -F 4 - | \
            samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        }} 2>&1 | tee {log}
        """


# Track B: CoverM Contig Standard (95/75)
# ------------------------------------------------------------------------------
rule annotation_circular_coverm_standard:
    input:
        bam = os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam"),
        bai = os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam.bai")
    output:
        tsv = os.path.join(AN_RES, "{genome}/circular_coverm/circular_contigs.coverage.tsv")
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_coverm_standard.log")
    threads: get_threads(AN_C.get("threads", {}).get("coverm", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("coverm"),
        min_pid = lambda wildcards: AN_C.get("circular_illumina", {}).get("standard", {}).get("min_identity", 95),
        min_aln_pct = lambda wildcards: AN_C.get("circular_illumina", {}).get("standard", {}).get("min_aligned_percent", 75),
        end_exclusion = lambda wildcards: AN_C.get("circular_illumina", {}).get("coverm", {}).get("contig_end_exclusion", 75)
    shell:
        r"""
        {{
        set -euo pipefail
        export PATH="{params.env}/bin:$PATH"
        coverm contig \
            -b {input.bam} \
            --min-read-percent-identity {params.min_pid} \
            --min-read-aligned-percent {params.min_aln_pct} \
            --contig-end-exclusion {params.end_exclusion} \
            --exclude-supplementary \
            -m mean trimmed_mean covered_fraction variance count \
            -t {threads} > {output.tsv}
        }} 2>&1 | tee {log}
        """


# Track B: CoverM Contig Strict (98/90)
# ------------------------------------------------------------------------------
rule annotation_circular_coverm_strict:
    input:
        bam = os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam"),
        bai = os.path.join(AN_RES, "{genome}/circular_coverm/original.raw.bam.bai")
    output:
        tsv = os.path.join(AN_RES, "{genome}/circular_coverm/circular_contigs.coverage.strict.tsv")
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_coverm_strict.log")
    threads: get_threads(AN_C.get("threads", {}).get("coverm", 8))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("coverm"),
        strict_pid = lambda wildcards: AN_C.get("circular_illumina", {}).get("strict", {}).get("min_identity", 98),
        strict_aln_pct = lambda wildcards: AN_C.get("circular_illumina", {}).get("strict", {}).get("min_aligned_percent", 90),
        end_exclusion = lambda wildcards: AN_C.get("circular_illumina", {}).get("coverm", {}).get("contig_end_exclusion", 75)
    shell:
        r"""
        {{
        set -euo pipefail
        export PATH="{params.env}/bin:$PATH"
        coverm contig \
            -b {input.bam} \
            --min-read-percent-identity {params.strict_pid} \
            --min-read-aligned-percent {params.strict_aln_pct} \
            --contig-end-exclusion {params.end_exclusion} \
            --exclude-supplementary \
            -m trimmed_mean covered_fraction count \
            -t {threads} > {output.tsv}
        }} 2>&1 | tee {log}
        """


# Track C: Joint Rotated Assembly Generation
# ------------------------------------------------------------------------------
rule annotation_circular_build_rotated_assembly:
    input:
        assembly = lambda wildcards: an_get_asm_path(get_sample_name_from_genome(wildcards.genome)),
        bins_dir = os.path.join(AN_RES, "{genome}/circular_bins_filtered")
    output:
        rotated_fa = temp(os.path.join(AN_RES, "{genome}/circular_coverm/rotated_assembly.fasta")),
        coords_tsv = os.path.join(AN_RES, "{genome}/circular_coverm/rotated_junction_coords.tsv")
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_build_rotated_assembly.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        fraction = lambda wildcards: AN_C.get("circular_illumina", {}).get("junction", {}).get("rotation_fraction", 0.5)
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/python3" ]; then
            for cand in "{params.env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit"; do
                if [ -x "$cand/bin/python3" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        mkdir -p $(dirname {output.rotated_fa})
        $ENV_DIR/bin/python3 scripts/validate_rotated_junction.py rotate \
            --assembly {input.assembly} \
            --circular-bins-dir {input.bins_dir} \
            --fraction {params.fraction} \
            --out-assembly {output.rotated_fa} \
            --out-coords {output.coords_tsv}
        }} 2>&1 | tee {log}
        """


# Track C: Rotated Assembly Short-Read Alignment
# ------------------------------------------------------------------------------
rule annotation_circular_rotated_mapping:
    input:
        assembly = os.path.join(AN_RES, "{genome}/circular_coverm/rotated_assembly.fasta"),
        r1 = get_sample_short_r1,
        r2 = get_sample_short_r2
    output:
        bam = temp(os.path.join(AN_RES, "{genome}/circular_coverm/rotated.raw.bam")),
        bai = temp(os.path.join(AN_RES, "{genome}/circular_coverm/rotated.raw.bam.bai"))
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_rotated_mapping.log")
    threads: get_threads(AN_C.get("threads", {}).get("minimap2", 16))
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        mapping_env = an_get_env("mapping")
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/minimap2" ]; then
            for cand in "{params.env}" \
                        "{params.mapping_env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" \
                        "$HOME/.conda/envs/bio"; do
                if [ -x "$cand/bin/minimap2" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        for extra in "{params.mapping_env}" "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio" "$HOME/.conda/envs/bio"; do
            if [ -x "$extra/bin/samtools" ]; then
                export PATH="$extra/bin:$PATH"
            fi
        done
        mkdir -p $(dirname {output.bam})
        minimap2 -ax sr -t {threads} -p 1 {input.assembly} {input.r1} {input.r2} | \
            samtools view -u -F 4 - | \
            samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        }} 2>&1 | tee {log}
        """


# Track C: Junction Depth, Ratio, and Physical Spanning Evaluation
# ------------------------------------------------------------------------------
rule annotation_circular_junction_metrics:
    input:
        bam = os.path.join(AN_RES, "{genome}/circular_coverm/rotated.raw.bam"),
        bai = os.path.join(AN_RES, "{genome}/circular_coverm/rotated.raw.bam.bai"),
        coords = os.path.join(AN_RES, "{genome}/circular_coverm/rotated_junction_coords.tsv"),
        contig_coverage = os.path.join(AN_RES, "{genome}/circular_coverm/circular_contigs.coverage.tsv")
    output:
        tsv = os.path.join(AN_RES, "{genome}/circular_coverm/rotated_junction_qc.tsv")
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_junction_metrics.log")
    threads: get_threads(2)
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        narrow_window = lambda wildcards: AN_C.get("circular_illumina", {}).get("junction", {}).get("narrow_window_bp", 100),
        flank_inner = lambda wildcards: AN_C.get("circular_illumina", {}).get("junction", {}).get("flank_inner_bp", 1000),
        flank_outer = lambda wildcards: AN_C.get("circular_illumina", {}).get("junction", {}).get("flank_outer_bp", 5000),
        min_anchor = lambda wildcards: AN_C.get("circular_illumina", {}).get("junction", {}).get("min_spanning_anchor_bp", 30),
        min_mapq = lambda wildcards: AN_C.get("circular_illumina", {}).get("junction", {}).get("min_mapq", 20)
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/python3" ]; then
            for cand in "{params.env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit"; do
                if [ -x "$cand/bin/python3" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        $ENV_DIR/bin/python3 scripts/validate_rotated_junction.py evaluate \
            --bam {input.bam} \
            --coords {input.coords} \
            --contig-coverage-tsv {input.contig_coverage} \
            --narrow-window {params.narrow_window} \
            --flank-inner {params.flank_inner} \
            --flank-outer {params.flank_outer} \
            --min-anchor {params.min_anchor} \
            --min-mapq {params.min_mapq} \
            --out-tsv {output.tsv}
        }} 2>&1 | tee {log}
        """


# Synthesis: Compile Circular Illumina QC Summary
# ------------------------------------------------------------------------------
rule annotation_circular_illumina_qc_summary:
    input:
        genome_abundance = os.path.join(AN_RES, "{genome}/circular_coverm/circular_genome_abundance.tsv"),
        contig_coverage = os.path.join(AN_RES, "{genome}/circular_coverm/circular_contigs.coverage.tsv"),
        strict_coverage = os.path.join(AN_RES, "{genome}/circular_coverm/circular_contigs.coverage.strict.tsv"),
        junction_qc = os.path.join(AN_RES, "{genome}/circular_coverm/rotated_junction_qc.tsv")
    output:
        summary = os.path.join(AN_RES, "{genome}/circular_coverm/circular_illumina_qc_summary.tsv")
    log:
        os.path.join(AN_LOGS, "circular_coverm/{genome}_illumina_qc_summary.log")
    threads: 1
    resources:
        runtime="7d",
        mem_mb=lambda wildcards, threads: get_mem_mb(threads)
    params:
        env = an_get_env("longread_audit"),
        sample = lambda wildcards: get_sample_name_from_genome(wildcards.genome),
        low_ra = lambda wildcards: AN_C.get("circular_illumina", {}).get("thresholds", {}).get("low_relative_abundance", 0.1),
        low_depth = lambda wildcards: AN_C.get("circular_illumina", {}).get("thresholds", {}).get("low_depth", 5.0),
        low_breadth = lambda wildcards: AN_C.get("circular_illumina", {}).get("thresholds", {}).get("low_breadth", 0.95),
        comp_warn = lambda wildcards: AN_C.get("circular_illumina", {}).get("thresholds", {}).get("competitive_warning_ratio", 0.50),
        id_sens = lambda wildcards: AN_C.get("circular_illumina", {}).get("thresholds", {}).get("identity_sensitive_ratio", 0.50),
        junc_warn = lambda wildcards: AN_C.get("circular_illumina", {}).get("thresholds", {}).get("junction_warning_ratio", 0.70),
        junc_fail = lambda wildcards: AN_C.get("circular_illumina", {}).get("thresholds", {}).get("junction_fail_ratio", 0.30)
    shell:
        r"""
        {{
        set -euo pipefail
        ENV_DIR="{params.env}"
        if [ ! -x "$ENV_DIR/bin/python3" ]; then
            for cand in "{params.env}" \
                        "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/longread_audit" \
                        "$HOME/.conda/envs/longread_audit"; do
                if [ -x "$cand/bin/python3" ]; then
                    ENV_DIR="$cand"
                    break
                fi
            done
        fi
        export PATH="$ENV_DIR/bin:$PATH"
        $ENV_DIR/bin/python3 scripts/compile_circular_illumina_qc.py \
            --sample "{params.sample}" \
            --genome-abundance-tsv {input.genome_abundance} \
            --contig-coverage-tsv {input.contig_coverage} \
            --strict-coverage-tsv {input.strict_coverage} \
            --junction-qc-tsv {input.junction_qc} \
            --low-ra-threshold {params.low_ra} \
            --low-depth-threshold {params.low_depth} \
            --low-breadth-threshold {params.low_breadth} \
            --competitive-warning-threshold {params.comp_warn} \
            --identity-sensitive-threshold {params.id_sens} \
            --junction-warning-threshold {params.junc_warn} \
            --junction-fail-threshold {params.junc_fail} \
            --out-tsv {output.summary}
        }} 2>&1 | tee {log}
        """
