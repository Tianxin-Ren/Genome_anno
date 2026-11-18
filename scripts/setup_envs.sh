#!/bin/bash

# Bootstrap script to setup Conda environments for Genome_anno pipeline
# This script installs environments into the project-local .conda_envs/ directory.

set -euo pipefail

# 1. Configuration
ENV_DIR="workflow/envs"
TARGET_DIR=".conda_envs"
CONDA_CMD="conda"

# 2. Check for mamba (recommended for speed)
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo ">>> Found mamba, using it for faster installation."
else
    echo ">>> mamba not found, falling back to conda."
fi

mkdir -p "$TARGET_DIR"

# 3. Iterate over YAML files
echo ">>> Starting environment bootstrap..."
for yaml_file in "$ENV_DIR"/*.yaml; do
    env_name=$(basename "$yaml_file" .yaml)
    prefix="$TARGET_DIR/$env_name"
    
    echo "----------------------------------------------------------------"
    if [ -d "$prefix" ]; then
        echo ">>> Environment '$env_name' already exists at $prefix. Skipping."
        echo ">>> To force reinstall, delete the directory and run again."
    else
        echo ">>> Installing '$env_name' into $prefix..."
        $CONDA_CMD env create -f "$yaml_file" -p "$prefix" --quiet
        echo ">>> Successfully installed '$env_name'."
    fi
done

echo "----------------------------------------------------------------"
echo ">>> All environments processed."
echo ">>> IMPORTANT: Update your config.yaml files to point to the paths in $TARGET_DIR/"
echo ">>> Example: gtdbtk: \"$(pwd)/$TARGET_DIR/gtdbtk\""
