#!/usr/bin/env bash
# Read-only database checker for the unified BioB workflow.
# It never downloads, deletes, or modifies database files.

if [ -z "${BASH_VERSION:-}" ]; then
  echo "This script requires bash. Run it as:"
  echo "  bash config_En/check_existing_databases.sh"
  exit 2
fi

set -euo pipefail

DBROOT="${DBROOT:-/dssg/home/acct-clswfp/clswfp-icdli1/Database}"
MAX_SHOW="${MAX_SHOW:-5}"

ok=0
miss=0

print_hit() {
  local label="$1"
  local path="$2"
  echo "[OK]   ${label} -> ${path}"
  ok=$((ok + 1))
}

print_miss() {
  local label="$1"
  local target="$2"
  echo "[MISS] ${label} -> ${target}"
  miss=$((miss + 1))
}

check_path() {
  local label="$1"
  local path="$2"
  if [ -e "${path}" ]; then
    print_hit "${label}" "${path}"
  else
    print_miss "${label}" "${path}"
  fi
}

check_find() {
  local label="$1"
  local dir="$2"
  shift 2

  if [ ! -d "${dir}" ]; then
    print_miss "${label}" "${dir}"
    return
  fi

  mapfile -t hits < <(find "${dir}" "$@" -print 2>/dev/null | head -n "${MAX_SHOW}")
  if [ "${#hits[@]}" -gt 0 ]; then
    print_hit "${label}" "${dir}"
    printf '       %s\n' "${hits[@]}"
  else
    print_miss "${label}" "${dir} :: find $*"
  fi
}

check_find_limited() {
  local label="$1"
  shift
  local found=0
  local dir

  for dir in "$@"; do
    if [ -d "${dir}" ]; then
      mapfile -t hits < <(find "${dir}" -maxdepth 4 -type f \( -name "Rfam.cm" -o -name "Rfam.cm.i1*" \) -print 2>/dev/null | head -n "${MAX_SHOW}")
      if [ "${#hits[@]}" -gt 0 ]; then
        print_hit "${label}" "${dir}"
        printf '       %s\n' "${hits[@]}"
        found=1
      fi
    fi
  done

  if [ "${found}" -eq 0 ]; then
    print_miss "${label}" "checked only likely Rfam dirs: ${*}"
  fi
}

echo "Database root: ${DBROOT}"
echo

echo "== Core QC / taxonomy =="
check_path "CheckM2 DIAMOND DB" "${DBROOT}/CheckM2_database/uniref100.KO.1.dmnd"
check_find "GTDB bac120 metadata" "${DBROOT}/GTDB" -type f -name "bac120_metadata*.tsv.gz"
check_find "GTDB ar53 metadata" "${DBROOT}/GTDB" -type f -name "ar53_metadata*.tsv.gz"
check_find "GTDB marker/pplacer files" "${DBROOT}/GTDB" -type f \( -name "gtdbtk.*.json" -o -name "*.hmm" -o -name "*.refpkg" \)

echo
echo "== Functional annotation =="
check_path "eggNOG directory" "${DBROOT}/eggnog"
check_find "eggNOG DB files" "${DBROOT}/eggnog" -type f \( -name "*.db" -o -name "*.dmnd" -o -name "*mmseqs*" \)
check_path "KOfam profiles" "${DBROOT}/kofam/profiles"
check_path "KOfam ko_list" "${DBROOT}/kofam/ko_list"
check_path "Pfam directory" "${DBROOT}/Pfam"
check_find "Pfam HMM files" "${DBROOT}/Pfam" -type f \( -name "Pfam-A.hmm" -o -name "Pfam-A.hmm.h3*" \)
check_path "CDD directory" "${DBROOT}/CDD"
check_find "CDD files" "${DBROOT}/CDD" -type f

echo
echo "== RNA annotation =="
check_find_limited "Rfam covariance model / pressed index" \
  "${DBROOT}/Rfam" \
  "${DBROOT}/rfam" \
  "${DBROOT}/Infernal" \
  "${DBROOT}/infernal" \
  "${DBROOT}/RNA" \
  "${DBROOT}/rna"

echo
echo "== Protein similarity =="
check_path "Swiss-Prot directory" "${DBROOT}/swiss-prot"
check_find "Swiss-Prot FASTA or DB files" "${DBROOT}/swiss-prot" -maxdepth 3 -type f \( -name "*.faa" -o -name "*.fa" -o -name "*.fasta" -o -name "*.faa.gz" -o -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.dmnd" \)
check_path "NR directory" "${DBROOT}/NR"
check_find "NR BLAST/DIAMOND/FASTA files" "${DBROOT}/NR" -maxdepth 3 -type f \( -name "*.pin" -o -name "*.psq" -o -name "*.phr" -o -name "*.dmnd" -o -name "*.faa" -o -name "*.faa.gz" \)
check_path "MMseqs2 DB directory" "${DBROOT}/mmseqs2_db"
check_find "Any MMseqs2 DB-like files" "${DBROOT}/mmseqs2_db" -maxdepth 3 -type f
check_path "Swiss-Prot MMseqs2 DB optional build target" "${DBROOT}/mmseqs2_db/swiss-prot"
check_path "UniRef90 MMseqs2 DB optional" "${DBROOT}/mmseqs2_db/uniref90"
check_path "RefSeq bacteria/archaea MMseqs2 DB optional" "${DBROOT}/RefSeq/protein/refseq_bacteria_archaea.mmseqs"

echo
echo "== Optional reference / review =="
check_path "Bakta DB" "${DBROOT}/Bakta"
check_find "Bakta DB files" "${DBROOT}/Bakta" -maxdepth 4 -type f
check_path "BioB curated HMM DB" "${DBROOT}/BioB/hmm/biob_core_functions.hmm"
check_find "BioB curated HMM pressed index" "${DBROOT}/BioB/hmm" -maxdepth 1 -type f -name "biob_core_functions.hmm.h3*"
check_find "BioB source HMM files" "${DBROOT}/BioB/hmm/source_hmms" -type f -name "*.hmm"
check_path "BioB reference proteins" "${DBROOT}/BioB/reference_proteins"
check_path "Foldseek DB/resources" "${DBROOT}/foldseek"
check_find "Foldseek DB files" "${DBROOT}/foldseek" -maxdepth 3 -type f
check_path "KEGG" "${DBROOT}/KEGG"
check_path "arCOG" "${DBROOT}/arCOG"

echo
echo "== Summary =="
echo "OK:   ${ok}"
echo "MISS: ${miss}"
