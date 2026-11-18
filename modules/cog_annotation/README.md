# COG2024-Compatible Footprint-Aware Annotation Engine

## Overview
This module provides a production-grade, highly scalable COG functional annotation pipeline tailored for large-scale bacterial and archaeal datasets (~22,343 Archaea genomes, ~40 million proteins).

Unlike naive best-hit approaches, this engine enforces official COG footprint principles:
1. **Gate 0 Provenance Audit**: Verifies database consistency between NCBI CDD Cog_LE PSSMs (`cddid.tbl`) and COG2024 definitions (`cog-24.def.tab`). Results are conservatively designated as:
   `CDD Cog_LE-based COG annotation mapped to COG2024 definitions`.
2. **HSP-to-Footprint Collinear Aggregation**: Aggregates multi-HSP alignments on the same query $\times$ PSSM before computing coverage.
3. **Three Orthogonal Classification Dimensions**:
   - **Match Quality**: Evaluates COG2024 75% coverage rule on aggregated footprints (`cog_match_code` 0–3, `profile_coverage_class`, `query_coverage_class`).
   - **Region Evidence**: Distinguishes unique, competing, and ambiguous profiles using continuous score margins ($\Delta\text{bitscore}$, $R_\Delta = \frac{S_1 - S_2}{S_1}$, `bitscore_ratio`).
   - **Protein Architecture**: Classifies `single_cog`, `multi_cog`, `repeated_cog`, and `mixed_multi_repeated`, alongside boolean flags (`is_multi_cog`, `is_repeated_cog`, `has_competing_profiles`, `has_ambiguous_region`).
4. **Decoupled Gene Truncation Vector**: Evaluates Prodigal partial flags, contig edge distance, and length ratios independently.
5. **Comparative Genomics Matrices**: Produces both Protein Incidence counts and Footprint Frequency counts, alongside Protein-normalized ($\sum = N_{\text{proteins}}$), Footprint-normalized, and Multi-label category matrices.
6. **HPC Optimization**: Node-local scratch caching (`$SLURM_TMPDIR`) with atomic `flock` and node-worker bundling.

## Directory Structure
```
modules/cog_annotation/
├── cog_annotation.smk         # Main Snakemake workflow module
├── config.yaml                # Default module configuration
├── README.md                  # This documentation
└── scripts/
    ├── audit_database_provenance.py   # Gate 0 Provenance Audit
    ├── assign_cog_footprints.py       # Core Footprint Arbitration Engine
    ├── expand_dedup_cog.py            # Streaming Master Table Expander
    ├── generate_cog_matrices.py       # Comparative Matrices Generator
    ├── validate_cog_qc.py             # 6-Gate Scientific Invariant Validator
    ├── benchmark_cog_recovery.py      # Ground-truth Footprint IoU Benchmark
    └── run_cog_node_worker.py         # SLURM 64-core Node Worker Manager
```

## Running the Pipeline
Enable the module in your project configuration (e.g., `projects/config_mag_annotation/pipeline.yaml`):
```yaml
modules:
  cog_annotation:
    path: "./cog_annotation.yaml"
    enable: true
```

Run Snakemake:
```bash
snakemake --workflow-profile profiles/slurm --config conf=projects/config_mag_annotation
```
