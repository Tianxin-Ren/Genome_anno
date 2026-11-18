#!/usr/bin/env bash
# Collect the most relevant logs after a unified_annotation Snakemake failure.

set -euo pipefail

ROOT="${ROOT:-results/unified_annotation}"
SNAKE_LOG_DIR="${SNAKE_LOG_DIR:-.snakemake/slurm_logs}"
LINES="${LINES:-80}"

show_tail() {
  local label="$1"
  local path="$2"
  echo
  echo "===== ${label}: ${path} ====="
  if [ -s "${path}" ]; then
    tail -n "${LINES}" "${path}"
  elif [ -e "${path}" ]; then
    echo "[EMPTY] ${path}"
  else
    echo "[MISSING] ${path}"
  fi
}

echo "Collecting logs under: ${ROOT}"

for subdir in rfam eggnog interpro pyrodigal checkm2 gtdbtk quast trnascan; do
  if [ -d "${ROOT}/logs/${subdir}" ]; then
    while IFS= read -r log; do
      show_tail "rule log" "${log}"
    done < <(find "${ROOT}/logs/${subdir}" -maxdepth 1 -type f -name "*.log" | sort | head -n 12)
  fi
done

if [ -d "${SNAKE_LOG_DIR}" ]; then
  echo
  echo "===== Recent SLURM wrapper logs ====="
  find "${SNAKE_LOG_DIR}" -type f -name "*.log" -printf "%T@ %p\n" 2>/dev/null \
    | sort -nr \
    | head -n 20 \
    | cut -d' ' -f2- \
    | while IFS= read -r log; do show_tail "slurm log" "${log}"; done
fi

