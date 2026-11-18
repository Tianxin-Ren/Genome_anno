#!/usr/bin/env bash
# Read-only software checker for the unified BioB workflow.
# It tests whether required conda envs/commands can start and report help/version.

if [ -z "${BASH_VERSION:-}" ]; then
  echo "This script requires bash. Run it as:"
  echo "  bash config_En/check_software_working.sh"
  exit 2
fi

set -u

ENVS="${ENVS:-/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs}"

ok=0
fail=0
warn=0

pass() {
  echo "[OK]   $1"
  ok=$((ok + 1))
}

bad() {
  echo "[FAIL] $1"
  fail=$((fail + 1))
}

note() {
  echo "[WARN] $1"
  warn=$((warn + 1))
}

env_exists() {
  local env_name="$1"
  if [[ "${env_name}" = /* ]]; then
    [ -d "${env_name}" ]
  else
    [ -d "${ENVS}/${env_name}" ]
  fi
}

env_path() {
  local env_name="$1"
  if [[ "${env_name}" = /* ]]; then
    printf '%s\n' "${env_name}"
  else
    printf '%s\n' "${ENVS}/${env_name}"
  fi
}

check_env() {
  local env_name="$1"
  if env_exists "${env_name}"; then
    pass "env exists: ${env_name}"
    return 0
  fi
  bad "env missing: ${env_name}"
  return 1
}

check_cmd() {
  local label="$1"
  local env_name="$2"
  shift 2

  if ! env_exists "${env_name}"; then
    bad "${label}: env missing (${env_name})"
    return 1
  fi

  local tmp
  tmp="$(mktemp)"
  if conda run -p "$(env_path "${env_name}")" "$@" >"${tmp}" 2>&1; then
    pass "${label}: $*"
    head -n 3 "${tmp}" | sed 's/^/       /'
  else
    bad "${label}: $*"
    head -n 8 "${tmp}" | sed 's/^/       /'
  fi
  rm -f "${tmp}"
}

check_cmd_optional() {
  local label="$1"
  local env_name="$2"
  shift 2

  if ! env_exists "${env_name}"; then
    note "${label}: optional env missing (${env_name})"
    return 0
  fi

  local tmp
  tmp="$(mktemp)"
  if conda run -p "$(env_path "${env_name}")" "$@" >"${tmp}" 2>&1; then
    pass "${label}: $*"
    head -n 3 "${tmp}" | sed 's/^/       /'
  else
    note "${label}: optional command did not start cleanly ($*)"
    head -n 8 "${tmp}" | sed 's/^/       /'
  fi
  rm -f "${tmp}"
}

check_cmd_exists() {
  local label="$1"
  local env_name="$2"
  local cmd_name="$3"

  if ! env_exists "${env_name}"; then
    bad "${label}: env missing (${env_name})"
    return 1
  fi

  local tmp
  tmp="$(mktemp)"
  if conda run -p "$(env_path "${env_name}")" which "${cmd_name}" >"${tmp}" 2>&1; then
    pass "${label}: command found (${cmd_name})"
    head -n 2 "${tmp}" | sed 's/^/       /'
  else
    bad "${label}: command not found (${cmd_name})"
    head -n 5 "${tmp}" | sed 's/^/       /'
  fi
  rm -f "${tmp}"
}

check_cmd_any_version() {
  local label="$1"
  local env_name="$2"
  shift 2

  if ! env_exists "${env_name}"; then
    bad "${label}: env missing (${env_name})"
    return 1
  fi

  local cmd_name
  local tmp
  for cmd_name in "$@"; do
    tmp="$(mktemp)"
    if conda run -p "$(env_path "${env_name}")" "${cmd_name}" --version >"${tmp}" 2>&1; then
      pass "${label}: ${cmd_name} --version"
      head -n 3 "${tmp}" | sed 's/^/       /'
      rm -f "${tmp}"
      return 0
    fi
    rm -f "${tmp}"
  done

  bad "${label}: none of these version commands worked: $*"
}

check_file_exec() {
  local label="$1"
  local path="$2"
  if [ -x "${path}" ]; then
    pass "${label}: ${path}"
  elif [ -e "${path}" ]; then
    bad "${label}: exists but is not executable: ${path}"
  else
    bad "${label}: missing: ${path}"
  fi
}

echo "== Conda =="
if command -v conda >/dev/null 2>&1; then
  pass "conda found: $(command -v conda)"
  conda --version | sed 's/^/       /'
else
  bad "conda not found in PATH"
fi

echo
echo "== Quality and taxonomy =="
check_cmd "QUAST" "quast_5.3.0" quast.py --version
check_cmd "CheckM2" "checkm2_1.0.2" checkm2 --version
check_cmd "GTDB-Tk" "gtdbtk_2.7.1" gtdbtk --version

echo
echo "== CDS and RNA prediction =="
check_cmd "Pyrodigal" "pyrodigal_3.7.1" pyrodigal --version
check_cmd "tRNAscan-SE" "trnascan-se_2.0.12" tRNAscan-SE -h
check_cmd "Infernal cmscan" "infernal_1.1.5" cmscan -h
check_cmd "Infernal cmpress" "infernal_1.1.5" cmpress -h
check_cmd "ARAGORN" "aragorn_1.2.41" aragorn -h
note "CRISPR/MinCED checks skipped: the CRISPR workflow step is not implemented yet."

echo
echo "== Protein functional annotation =="
check_cmd "MMseqs2" "mmseqs2_17.b804f" mmseqs version
check_cmd "eggNOG-mapper" "eggnog-mapper" emapper.py --version
check_cmd "KOfamScan" "kofamscan" exec_annotation -h
check_cmd "HMMER hmmsearch" "hmmer_3.4" hmmsearch -h

echo
echo "== Localization and topology =="
check_cmd_optional "SignalP6 optional" "signalp6_6.0h" signalp6 --help
check_cmd_optional "DeepTMHMM BioLib package optional" "deeptmhmm_1.0" python -c "import biolib; print('pybiolib import OK')"
note "DeepTMHMM may still require BioLib login/network/model access for real predictions."

echo
echo "== Future deep disambiguation (optional) =="
check_cmd_optional "MAFFT optional" "mafft_7.525" mafft --version
check_cmd_optional "IQ-TREE optional" "iqtree_3.1.2" iqtree2 --version
check_cmd_optional "Foldseek optional" "foldseek_5.53465f0" foldseek version

echo
echo "== Optional review tools =="
check_cmd_optional "Bakta optional review" "bakta_1.12.0" bakta --version
check_cmd_optional "Prokka optional review" "prokka" prokka --version

echo
echo "== InterProScan =="
check_file_exec "InterProScan executable" "${ENVS}/interproscan_5.77-108.0/interproscan.sh"
if [ -x "${ENVS}/interproscan_5.77-108.0/interproscan.sh" ]; then
  tmp="$(mktemp)"
  if env_exists interproscan_java11; then
    interpro_cmd=(conda run -p "$(env_path interproscan_java11)" "${ENVS}/interproscan_5.77-108.0/interproscan.sh" --version)
  else
    interpro_cmd=("${ENVS}/interproscan_5.77-108.0/interproscan.sh" --version)
  fi
  if "${interpro_cmd[@]}" >"${tmp}" 2>&1; then
    pass "InterProScan version command"
    head -n 3 "${tmp}" | sed 's/^/       /'
  else
    bad "InterProScan version command"
    head -n 8 "${tmp}" | sed 's/^/       /'
  fi
  rm -f "${tmp}"
fi

echo
echo "== Summary =="
echo "OK:   ${ok}"
echo "WARN: ${warn}"
echo "FAIL: ${fail}"

if [ "${fail}" -gt 0 ]; then
  exit 1
fi
