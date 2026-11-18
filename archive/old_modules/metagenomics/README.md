# Metagenomics Pipeline Execution Guide

This guide will walk you through configuring and running the Metagenomics pipeline (`metagenomics.smk`).

## 1. Configure the Pipeline

Before running the pipeline, you need to set up your project paths and choose your steps in the configuration file.

👉 **[Click here to edit: `modules/metagenomics/config.yaml`](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/modules/metagenomics/config.yaml)**

### Key Settings to Check:
- **`project_dir`**: Must be the absolute path to your workspace.
- **`data_dir` / `raw_fq_dir`**: Where your raw reads are located (relative to project).
- **`steps`**: Toggle modules like `run_drep`, `run_gtdbtk`, `run_coverm`.

## 2. Result Directory Structure (Skill 9C)

Outputs are organized hierarchically to keep your project directory clean:

- `{sample}/01_qc/fastp/`: Cleaned reads and QC reports.
- `{sample}/02_assembly/megahit/`: Raw assembly contigs.
- `{sample}/03_binning/metawrap/`: Initial bin sets (Metabat2, Maxbin2, Concoct).
- `{sample}/04_refinement/metawrap/`: High-quality refined MAGs and CheckM stats.
- `{sample}/05_taxonomy/`: GTDB-Tk classification per sample.
- `{sample}/06_abundance/`: CoverM relative abundance tables.
- `MAGs/drep_results/`: Dereplicated global genome set.

## 3. Dry Run (Test the Configuration)

Always perform a dry run first to ensure Snakemake can resolve the DAG (Directed Acyclic Graph) without errors.

```bash
snakemake -n --workflow-profile workflow/profiles/slurm --configfile modules/metagenomics/config.yaml
```

## 4. Submit to SLURM

Once the dry run is successful, submit the full pipeline to the cluster. The `slurm` profile will automatically handle job submission, memory, and threads based on the rules.

```bash
snakemake --workflow-profile workflow/profiles/slurm --configfile modules/metagenomics/config.yaml
```

> [!TIP]
> **Running Multiple Projects?**
> You don't need to change the YAML file every time. You can override the project directory directly from the command line:
> ```bash
> snakemake --workflow-profile workflow/profiles/slurm \
>           --configfile modules/metagenomics/config.yaml \
>           --config project_dir=/absolute/path/to/new_project
> ```
