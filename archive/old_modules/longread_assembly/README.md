# Long-Read Assembly Pipeline Execution Guide

This guide covers the assembly of ONT (Oxford Nanopore) or HiFi (PacBio) reads using the `longread_assembly.smk` pipeline.

## 1. Configure the Pipeline

Customize the assembly parameters and paths in the YAML configuration.

👉 **[Click here to edit: `modules/longread_assembly/config.yaml`](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/modules/longread_assembly/config.yaml)**

### Critical Parameters:
- **`project_dir`**: Absolute path to your project root.
- **`input_dir`**: Folder containing your `*.fastq.gz` files.
- **Naming convention**: Files should be named `{sample}_{platform}.fastq.gz` (e.g., `SampleA_ont.fastq.gz` or `SampleB_hifi.fastq.gz`).

## 2. Result Directory Structure (Skill 9C)

The pipeline organizes results into numbered steps for clear lineage:

- `{sample}/01_qc/nanoplot/`: NanoPlot reports and statistics.
- `{sample}/02_assembly/myloasm/`: 
  - `{sample}.myloasm.fa`: The final assembled primary contigs.
  - `circular_info.tsv`: Information on circularized contigs (plasmids/chromosomes).

## 3. Execution Commands

### Dry Run (Recommended)
```bash
snakemake -n --workflow-profile workflow/profiles/slurm --configfile modules/longread_assembly/config.yaml
```

### Full Run on SLURM
```bash
snakemake --workflow-profile workflow/profiles/slurm --configfile modules/longread_assembly/config.yaml
```

> [!IMPORTANT]
> **Resource Intensive**: Myloasm is highly multi-threaded and memory-intensive. It uses `--exclusive` nodes by default for maximum performance.
