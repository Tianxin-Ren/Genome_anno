#!/bin/bash

# Export script to capture existing Conda environments from absolute paths
# Run this on the OLD server where the absolute paths exist.

set -euo pipefail

# 1. Configuration
ENV_OUT_DIR="workflow/envs"
mkdir -p "$ENV_OUT_DIR"

# 2. Define environments (Key:Path)
# Keys will be used as the YAML filename
declare -A ENVS=(
    ["bio"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bio"
    ["metawrap_1.3.2"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/metawrap_1.3.2"
    ["bbmap_39.81"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/bbmap_39.81"
    ["megahit_1.2.9"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/megahit_1.2.9"
    ["drep_3.5.0"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/drep_3.5.0"
    ["gtdbtk-2.4.0"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/gtdbtk-2.4.0"
    ["coverm_0.7.0"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/coverm_0.7.0"
    ["nanoplot_1.20.0"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/nanoplot_1.20.0"
    ["myloasm"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/myloasm"
    ["checkm2_1.0.2"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/checkm2_1.0.2"
    ["quast_5.3.0"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/quast_5.3.0"
    ["bakta_1.12.0"]="/dssg/home/acct-clswfp/clswfp-icdli1/software/miniconda3/envs/bakta_1.12.0"
    ["eggnog-mapper"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/eggnog-mapper"
    ["kofamscan"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/kofamscan"
    ["dbcan"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/dbcan"
    ["multiqc_1.33"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/multiqc_1.33"
    ["biopython_1.85"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/biopython_1.85"
    ["genomespot_1.0.1"]="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/genomespot_1.0.1"
)

echo ">>> Starting environment export from absolute paths..."

for name in "${!ENVS[@]}"; do
    path="${ENVS[$name]}"
    output="$ENV_OUT_DIR/${name}.yaml"
    
    echo "----------------------------------------------------------------"
    if [ -d "$path" ]; then
        echo ">>> Exporting '$name' from $path..."
        # Using --no-builds for better portability between different HPC clusters
        conda env export -p "$path" --no-builds | grep -v "^prefix: " > "$output"
        echo ">>> Saved to $output"
    else
        echo ">>> WARNING: Path not found for '$name': $path"
        echo ">>> Skipping export for this environment."
    fi
done

echo "----------------------------------------------------------------"
echo ">>> Export completed. You can now use scripts/setup_envs.sh on the new server."
