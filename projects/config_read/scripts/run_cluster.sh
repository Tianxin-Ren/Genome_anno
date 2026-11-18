#!/bin/bash
#SBATCH -J sylph.OSMA_3nd_4and5_20241024.20260306
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -o sylph.OSMA_3nd_4and5_20241024.20260306.%A.out
#SBATCH -e sylph.OSMA_3nd_4and5_20241024.20260306.%A.err
#SBATCH --mail-user=liuyangli@sjtu.edu.cn
#SBATCH --time=7-00:00:00

set -euo pipefail

ulimit -s unlimited || true
ulimit -l unlimited || true

eval "$(conda shell.bash hook 2>/dev/null)" || true

# Navigate to workflow root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$WORKFLOW_ROOT"

echo "=== Launching Sylph Read Profiling Pipeline (config_read) ==="
echo "Workflow Root: $WORKFLOW_ROOT"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Time: $(date)"

# Execute Snakemake pipeline with platform profile or local profile
snakemake \
  --workflow-profile profiles/slurm \
  --config conf=projects/config_read \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds

echo "=== Pipeline Completed at $(date) ==="
