#!/bin/bash
#SBATCH --job-name=gtdb_pykofam_master
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/pykofam_master_%j.out
#SBATCH --error=logs/pykofam_master_%j.err

set -euo pipefail

# Navigate to workflow root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$WORKFLOW_ROOT"

echo "=== Launching Production PyKofamSearch Pipeline ==="
echo "Workflow Root: $WORKFLOW_ROOT"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Time: $(date)"

# Launch Snakemake master process with profiles/slurm
snakemake \
  --workflow-profile profiles/slurm \
  --config conf=projects/config_mag_annotation stage=annotation \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds

echo "=== Pipeline Completed at $(date) ==="
