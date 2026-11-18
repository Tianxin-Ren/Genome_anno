#!/usr/bin/env bash
# Environment installation commands for the unified BioB workflow.
# This script intentionally does not download or build databases.

set -euo pipefail

ENVS="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs"

env_exists() {
  conda env list | awk '{print $1}' | grep -qx "$1"
}

create_env_if_missing() {
  local env_name="$1"
  shift
  if env_exists "${env_name}"; then
    echo "[SKIP] conda env exists: ${env_name}"
  else
    conda create -y -n "${env_name}" "$@"
  fi
}

# Use strict channel priority for Bioconda installs.
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Already present in your HPC env list and therefore not recreated here:
# quast_5.3.0, checkm2_1.0.2/checkm2, gtdbtk_2.7.1, aragorn_1.2.41,
# mmseqs2_17.b804f, eggnog-mapper, kofamscan, hmmer_3.4, mafft_7.525,
# iqtree_3.1.2, foldseek_5.53465f0, bakta_1.12.0, prokka,
# crisprcasfinder_20260706.

create_env_if_missing pyrodigal_3.7.1 -c conda-forge -c bioconda pyrodigal=3.7.1
create_env_if_missing trnascan-se_2.0.12 -c conda-forge -c bioconda trnascan-se=2.0.12
create_env_if_missing infernal_1.1.5 -c conda-forge -c bioconda infernal=1.1.5
create_env_if_missing minced_0.4.2 -c conda-forge -c bioconda minced=0.4.2

# SignalP 6.0 is licence-gated. This conda package is a helper/placeholder;
# follow the SignalP licence/download instructions if the real model files are
# not bundled in your installation.
create_env_if_missing signalp6_6.0h -c predector -c conda-forge signalp6=6.0h

# Java 11 runtime for InterProScan.
create_env_if_missing interproscan_java11 -c conda-forge openjdk=11

# DeepTMHMM is usually run through BioLib. It needs network access on first use.
if env_exists deeptmhmm_1.0; then
  echo "[SKIP] conda env exists: deeptmhmm_1.0"
else
  conda create -y -n deeptmhmm_1.0 -c conda-forge python=3.11 pip
  conda run -n deeptmhmm_1.0 python -m pip install pybiolib
fi

# InterProScan is not normally installed as a small conda package. Keep it in
# one env-like directory named with the software version.
if [ -x "${ENVS}/interproscan_5.77-108.0/interproscan.sh" ]; then
  echo "[SKIP] InterProScan exists: ${ENVS}/interproscan_5.77-108.0"
else
  mkdir -p "${ENVS}/interproscan_5.77-108.0"
  cd "${ENVS}/interproscan_5.77-108.0"
  curl -L -O https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz
  curl -L -O https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz.md5
  md5sum -c interproscan-5.77-108.0-64-bit.tar.gz.md5
  tar -xzf interproscan-5.77-108.0-64-bit.tar.gz --strip-components=1
  python3 setup.py -f interproscan.properties
fi
