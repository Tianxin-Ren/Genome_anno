#!/bin/bash
# ==============================================================================
# Serialize KOfam Models and Generate Provenance Manifest
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$WORKFLOW_ROOT"

SERIALIZE_SCRIPT="$WORKFLOW_ROOT/modules/mag_annotation/scripts/serialize_kofam.py"

PROFILES_DIR="${1:-/dssg/home/acct-clswfp/clswfp-icdli1/Database/kofam/profiles}"
KO_LIST="${2:-/dssg/home/acct-clswfp/clswfp-icdli1/Database/kofam/ko_list}"
OUT_PKL="${3:-/dssg/home/acct-clswfp/clswfp-icdli1/Database/kofam/database.pkl}"
MANIFEST="${4:-/dssg/home/acct-clswfp/clswfp-icdli1/Database/kofam/kofam_manifest.json}"
PYKOFAM_ENV="${5:-pykofamsearch_26.7.21}"

echo "=== Serializing KOfam Models ==="
echo "Profiles: $PROFILES_DIR"
echo "KO List:  $KO_LIST"
echo "Output:   $OUT_PKL"
echo "Manifest: $MANIFEST"
echo "Conda Env: $PYKOFAM_ENV"

python3 "$SERIALIZE_SCRIPT" \
  -d "$PROFILES_DIR" \
  -k "$KO_LIST" \
  -o "$OUT_PKL" \
  -m "$MANIFEST" \
  --env "$PYKOFAM_ENV"

echo "Serialization completed successfully."
