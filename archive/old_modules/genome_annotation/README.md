# Genome Annotation Pipeline Execution Guide

This guide will walk you through configuring and running the Genome Annotation pipeline (`genome_annotation.smk`).

## 1. Configure the Pipeline

This pipeline supports multiple modes (`standard`, `hybrid`). Ensure you have selected the correct mode and provided the right inputs.

👉 **[Click here to edit: `modules/genome_annotation/config.yaml`](file:///Users/rentianxin/Documents/Mac_FN/code/Genome_anno/modules/genome_annotation/config.yaml)**

### Key Settings to Check:
- **`mode`**: Choose `"standard"` (for MAGs) or `"hybrid"` (for full assemblies).
- **`steps`**: Toggle `bakta`, `eggnog`, `dbcan`, `genomespot`.

## 2. Result Directory Structure (Skill 9C)

- `{genome}/01_qc/quast/`: Assembly quality metrics.
- `{genome}/02_bakta/`: Structural annotation (GFF, GBK, FAA).
- `{genome}/03_functional/`: eggNOG, KofamScan, and dbCAN results.
- `{genome}/04_phenotype/`: GenomeSpot predictions (Topt, pH, etc.).
- `results/annotation/All_Genes_Master_Table.tsv`: The final unified annotation matrix.

## 3. Dry Run (Test the Configuration)

Verify the execution plan and ensure all input genomes are correctly detected.

```bash
snakemake -n --workflow-profile workflow/profiles/slurm --configfile modules/genome_annotation/config.yaml
```

## 4. Submit to SLURM

Submit the pipeline to the cluster. 

```bash
snakemake --workflow-profile workflow/profiles/slurm --configfile modules/genome_annotation/config.yaml
```

> [!NOTE]
> All results will be seamlessly merged into a final `All_Genes_Master_Table.tsv` in your designated `output_dir`.
