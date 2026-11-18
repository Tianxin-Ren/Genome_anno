# Long-Read Assembly Pipeline
# Key Steps: ONT/HiFi assembly using Myloasm & QC.

import os

wildcard_constraints:
    sample="[^/]+",
    genome="[^/]+"

# 1. Setup & Path Logic
# ------------------------------------------------------------------------------
c_lr = config.get("longread_assembly", {})
LR_PROJ = os.path.abspath(c_lr.get("project_dir", "."))

def norm_lr_path(base, rel):
    if not rel: return base
    if os.path.isabs(rel): return os.path.normpath(rel)
    return os.path.normpath(os.path.join(base, rel))

LR_RES = norm_lr_path(LR_PROJ, c_lr.get("results_dir", "results/longread_assembly"))
LR_DATA = norm_lr_path(LR_PROJ, c_lr.get("data_dir", "data"))
LR_IN_DIR = norm_lr_path(LR_DATA, c_lr.get("input_dir", "raw_fq"))
LR_LOGS = os.path.join(LR_RES, "logs")

# Resolve samples
if c_lr.get("samples"):
    LR_SAMPLES = c_lr["samples"]
    if isinstance(LR_SAMPLES, str):
        LR_SAMPLES = [LR_SAMPLES]
else:
    LR_SAMPLES, = glob_wildcards(os.path.join(LR_IN_DIR, "{sample}.fastq.gz"))

# 2. Target Registration
# ------------------------------------------------------------------------------
TARGETS.extend(expand(os.path.join(LR_RES, "{sample}", "01_qc", "nanoplot"), sample=LR_SAMPLES))
TARGETS.extend(expand(os.path.join(LR_RES, "{sample}", "02_assembly", "myloasm", "{sample}.myloasm.fa"), sample=LR_SAMPLES))
TARGETS.extend(expand(os.path.join(LR_RES, "{sample}", "02_assembly", "myloasm", "circular_info.tsv"), sample=LR_SAMPLES))

# 3. Rules
# ------------------------------------------------------------------------------

rule run_nanoplot:
    input:
        fq = os.path.join(LR_IN_DIR, "{sample}.fastq.gz")
    output:
        dir = directory(os.path.join(LR_RES, "{sample}", "01_qc", "nanoplot"))
    log:
        os.path.join(LR_LOGS, "nanoplot", "{sample}.log")
    threads: 16
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        env = c_lr.get("conda_envs", {}).get("nanoplot", "nanoplot")
    shell:
        """
        set -euo pipefail
        conda run -p {params.env} NanoPlot \\
            -t {threads} \\
            --fastq {input.fq} \\
            -o {output.dir} 2>&1 | tee {log}
        """

rule run_myloasm:
    input:
        fq = os.path.join(LR_IN_DIR, "{sample}.fastq.gz")
    output:
        fa = os.path.join(LR_RES, "{sample}", "02_assembly", "myloasm", "{sample}.myloasm.fa")
    log:
        os.path.join(LR_LOGS, "myloasm", "{sample}.log")
    threads: c_lr.get("threads", 64)
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000,
        slurm_extra = "'--exclusive'"
    params:
        env = c_lr.get("conda_envs", {}).get("myloasm", "myloasm"),
        out_dir = os.path.join(LR_RES, "{sample}", "02_assembly", "myloasm_out"),
        platform = lambda wildcards: "--nano-r10" if "ont" in wildcards.sample.lower() else "--hifi",
        extra = c_lr.get("params", "")
    shell:
        """
        set -euo pipefail
        rm -rf {params.out_dir}
        conda run -p {params.env} myloasm \\
            -t {threads} \\
            -i {input.fq} \\
            -o {params.out_dir} \\
            {params.platform} \\
            {params.extra} 2>&1 | tee {log}
        cp {params.out_dir}/assembly.fasta {output.fa}
        """

rule extract_circularity_advanced:
    input:
        fasta = os.path.join(LR_RES, "{sample}", "02_assembly", "myloasm_out", "assembly.fasta")
    output:
        tsv = os.path.join(LR_RES, "{sample}", "02_assembly", "myloasm", "circular_info.tsv"),
        dir_flag = touch(os.path.join(LR_RES, "{sample}", "02_assembly", "myloasm", "circular_contigs", ".done"))
    log:
        os.path.join(LR_LOGS, "extract_circularity", "{sample}.log")
    threads: 1
    resources:
        runtime = "7d",
        mem_mb = lambda wildcards, threads: threads * 8000
    params:
        out_dir = lambda wildcards: os.path.join(LR_RES, wildcards.sample, "02_assembly", "myloasm", "circular_contigs")
    run:
        import os
        import re
        
        os.makedirs(params.out_dir, exist_ok=True)
        
        header_pattern = re.compile(
             r"^>?(?P<id>u[0-9]+)ctg_len-(?P<len>[0-9]+)_circular-(?P<status>[^_]+)_depth-(?P<depth>[^_]+)_duplicated-(?P<dup>[^_ ]+)( mult=(?P<mult>[^ ]+))?"
        )
        
        results = []
        current_seq = []
        current_header = None
        
        def save_contig(header, seq_lines):
            m = header_pattern.match(header)
            if m:
                info = m.groupdict()
                results.append(info)
                if info["status"].lower() in ["yes", "possibly"]:
                    with open(os.path.join(params.out_dir, f"{info['id']}.fa"), "w") as f:
                        f.write(f">{header}\n")
                        f.write("".join(seq_lines) + "\n")
        
        with open(input.fasta, "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if current_header:
                        save_contig(current_header, current_seq)
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header:
                save_contig(current_header, current_seq)
        
        with open(output.tsv, "w") as f:
            f.write("ID\tLength\tCircular\tDepth\tDuplicated\tMultiplicity\n")
            for res in results:
                f.write(f"{res['id']}\t{res['len']}\t{res['status']}\t{res['depth']}\t{res['dup']}\t{res.get('mult', '1.00')}\n")
