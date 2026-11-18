#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

project="projects/template"
profile="profiles/local"
platform=""

usage() {
    echo "Usage: scripts/run_pipeline.sh --project DIR [--profile DIR] [--platform FILE] [snakemake arguments...]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project)
            project="$2"
            shift 2
            ;;
        --profile)
            profile="$2"
            shift 2
            ;;
        --platform)
            platform="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            break
            ;;
    esac
done

[[ "${project}" = /* ]] || project="${WORKFLOW_ROOT}/${project}"
[[ "${profile}" = /* ]] || profile="${WORKFLOW_ROOT}/${profile}"

if [[ ! -f "${project}/pipeline.yaml" ]]; then
    echo "Project configuration not found: ${project}/pipeline.yaml" >&2
    exit 2
fi
if [[ ! -f "${profile}/config.yaml" ]]; then
    echo "Snakemake profile not found: ${profile}/config.yaml" >&2
    exit 2
fi

command=(
    snakemake
    "$@"
    --snakefile "${WORKFLOW_ROOT}/Snakefile"
    --profile "${profile}"
    --config "conf=${project}"
)
if [[ -n "${platform}" ]]; then
    [[ "${platform}" = /* ]] || platform="${WORKFLOW_ROOT}/${platform}"
    command+=("platform=${platform}")
fi

cd "${WORKFLOW_ROOT}"
exec "${command[@]}"
