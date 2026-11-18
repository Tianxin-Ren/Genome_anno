#!/bin/bash
# ==============================================================================
# Run Three-Level Scientific Invariance & Biological Consistency Validation Gate
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$WORKFLOW_ROOT"

VALIDATION_SCRIPT="$WORKFLOW_ROOT/modules/mag_annotation/scripts/validate_invariance.py"
OUT_REPORT="results/mag_annotation/validation/validation_report.json"
ANNOTATION_TSV="${1:-results/mag_annotation/merged/genome_KO_counts.long.tsv.gz}"

mkdir -p "results/mag_annotation/validation"

echo "=== Running Scientific Validation Gate (Level 1, 2, 3) ==="
python3 "$VALIDATION_SCRIPT" \
  -o "$OUT_REPORT" \
  --annotation_tsv "$ANNOTATION_TSV"

echo "Validation report generated at: $OUT_REPORT"
