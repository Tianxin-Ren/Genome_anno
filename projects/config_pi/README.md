# Snakemake Configuration Guide

This pipeline has been modularized to allow calling specific workflows independently.

## Modular Config Files

The configurations are now split into functional modules:

- `config/genome_annotation.yaml`: [NEW] Integrated Prokaryotic Genome Annotation (Bakta, Functional, Phenotype, etc.).
- `config/metagenomics.yaml`: Quality control, assembly, and binning.
- `config/mag_abundance.yaml`: Abundance calculation (TPM/RPKM) for MAGs.
- `config/longread_assembly.yaml`: ONT/HiFi assembly and plasmid analysis.
- `config/misc.yaml`: Homolog search and standalone plasmid analysis.

> [!NOTE]
> The separate `config/annotation.yaml` and `config/mag_annotation.yaml` are now integrated into `config/genome_annotation.yaml`.

## Usage

To run a specific workflow, pass its configuration file using the `--configfile` flag.

### Run the integrated annotation workflow
```bash
snakemake --use-conda --configfile config/genome_annotation.yaml
```

### Run a single legacy workflow (if still needed)
```bash
snakemake --use-conda --configfile config/metagenomics.yaml
```

### Run multiple workflows
```bash
snakemake --use-conda --configfile config/mag_annotation.yaml config/mag_abundance.yaml
```

### Samples Mapping
Sample mappings should still be defined in `config/samples.yaml`. This file is automatically loaded if it exists.

## Benefits
- Clearer separation of concerns.
- No need to toggle `enable: True/False` flags in a giant file.
- Individual modules can be versioned or shared more easily.
