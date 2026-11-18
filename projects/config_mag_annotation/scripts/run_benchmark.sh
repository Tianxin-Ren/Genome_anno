#!/bin/bash
# ==============================================================================
# Run 6-Point 3D Benchmark Matrix & CoW Memory Gate on Single 64c512g Node
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$WORKFLOW_ROOT"

BENCHMARK_SCRIPT="$WORKFLOW_ROOT/modules/mag_annotation/scripts/benchmark_matrix.py"
COW_SCRIPT="$WORKFLOW_ROOT/modules/mag_annotation/scripts/benchmark_cow.py"

# Configure test inputs (default to representative slice if available)
TEST_FASTA="${1:-results/mag_annotation/shards/shard_0001.faa.gz}"
DATABASE="${2:-/dssg/home/acct-clswfp/clswfp-icdli1/Database/kofam/database.pkl}"
OUT_DIR="results/mag_annotation/benchmarks"

echo "=== Starting 6-Point 3D Benchmark Matrix ==="
echo "Test FASTA: $TEST_FASTA"
echo "Database:   $DATABASE"
echo "Output:     $OUT_DIR"

python3 "$BENCHMARK_SCRIPT" \
  -i "$TEST_FASTA" \
  -d "$DATABASE" \
  -o "$OUT_DIR"

echo "=== Running CoW Memory Sharing Gate ==="
python3 "$COW_SCRIPT" \
  -d "$DATABASE" \
  -o "$OUT_DIR/cow_benchmark.json"

echo "=== Benchmark Suite Finished Successfully ==="
