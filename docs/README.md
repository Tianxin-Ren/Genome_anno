# Bioinformatics Pipeline Framework Documentation

Welcome to the documentation for the "Four Engines" bioinformatics pipeline. This repository is designed for high-performance, modular analysis of 2nd gen (Illumina) and 3rd gen (ONT/HiFi) sequencing data.

## Architecture: The Four Engines

The pipeline is organized into four core functional engines and optional add-on modules:

1.  **[Preprocess](modules/preprocess.md)**: Quality control for both short and long reads.
2.  **[Assembly](modules/assembly.md)**: De novo assembly using MEGAHIT (SR) or Myloasm (LR).
3.  **[Binning](modules/binning.md)**: Metagenomic MAG generation, supporting both SR-only and Hybrid (LR+SR) modes.
4.  **[Annotation](modules/annotation.md)**: Unified functional and structural annotation factory with special enhancements for 3rd gen genomes.

## Key Design Principles

-   **Modularity**: Each engine is self-contained in its own `modules/` directory.
-   **Config-Driven**: All logic is controlled via `config/pipeline.yaml`.
-   **HPC Ready**: Built-in support for SLURM with dynamic resource allocation.
-   **Data-Type Aware**: Intelligent detection of 2nd vs 3rd gen data to apply the best analytical tools.

---
*For quick start instructions, please refer to the [Usage Guide](usage.md).*
