#!/usr/bin/env bash
# Report placeholder outputs in results/unified_annotation.

set -euo pipefail

ROOT="${ROOT:-results/unified_annotation}"

if [ ! -d "${ROOT}" ]; then
  echo "Missing results directory: ${ROOT}" >&2
  exit 1
fi

is_sample_dir() {
  local dir="$1"
  [ -d "${dir}/input" ] || [ -d "${dir}/structural" ] || [ -d "${dir}/function" ] || [ -d "${dir}/genome" ]
}

status_file() {
  local label="$1"
  local path="$2"

  case "${label}" in
    rna_features_gff3)
      if [ -e "${path}" ] && grep -q "tRNAscan-SE table:" "${path}" && grep -q "Infernal/Rfam table:" "${path}"; then
        echo "POINTER_PLACEHOLDER"
        return
      fi
      ;;
    unresolved_proteins)
      if [ -s "${path}" ]; then
        echo "UNFILTERED_PLACEHOLDER"
        return
      fi
      ;;
  esac

  if [ ! -e "${path}" ]; then
    echo "MISSING"
  elif [ ! -s "${path}" ]; then
    echo "EMPTY_PLACEHOLDER"
  elif grep -qE 'failed with exit status|completed but did not|completed but produced no|pending_|See log:' "${path}"; then
    echo "PLACEHOLDER"
  elif [[ "${label}" == "protein_function_candidates" || "${label}" == "annotation_conflicts" ]] && [ "$(wc -l < "${path}")" -le 1 ]; then
    echo "HEADER_ONLY_PLACEHOLDER"
  else
    echo "REAL_OR_NONPLACEHOLDER"
  fi
}

print_check() {
  local sample="$1"
  local label="$2"
  local path="$3"
  local status
  status="$(status_file "${label}" "${path}")"
  printf "%s\t%s\t%s\t%s\n" "${sample}" "${label}" "${status}" "${path}"
}

echo -e "sample\toutput\tstatus\tpath"

for sample_dir in "${ROOT}"/*; do
  [ -d "${sample_dir}" ] || continue
  sample="$(basename "${sample_dir}")"
  case "${sample}" in
    logs|tmp) continue ;;
  esac
  is_sample_dir "${sample_dir}" || continue

  print_check "${sample}" "rfam_tbl" "${sample_dir}/rna/rfam/rfam.tbl"
  print_check "${sample}" "rna_features_gff3" "${sample_dir}/structural/rna_features.gff3"
  print_check "${sample}" "eggnog_annotations" "${sample_dir}/function/eggnog/${sample}.emapper.annotations"
  print_check "${sample}" "interpro_tsv" "${sample_dir}/function/interpro/interpro.tsv"
  print_check "${sample}" "protein_function_candidates" "${sample_dir}/function/protein_function_candidates.tsv"
  print_check "${sample}" "protein_function_consensus" "${sample_dir}/function/protein_function_consensus.jsonl"
  print_check "${sample}" "annotation_conflicts" "${sample_dir}/function/annotation_conflicts.tsv"
  print_check "${sample}" "unresolved_proteins" "${sample_dir}/function/unresolved_proteins.faa"
  print_check "${sample}" "genome_function_inventory" "${sample_dir}/genome/genome_function_inventory.json"
  print_check "${sample}" "complex_reconstruction" "${sample_dir}/genome/complex_reconstruction.json"
  print_check "${sample}" "reaction_inventory" "${sample_dir}/genome/reaction_inventory.json"
  print_check "${sample}" "transport_inventory" "${sample_dir}/genome/transport_inventory.json"
  print_check "${sample}" "pathway_closure" "${sample_dir}/genome/pathway_closure.json"
done
