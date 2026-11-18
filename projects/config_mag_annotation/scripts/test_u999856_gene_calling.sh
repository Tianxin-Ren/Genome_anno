#!/usr/bin/env bash
set -euo pipefail

# One-genome gene-calling comparison for u999856.
#
# Compared modes:
#   1) current run_pyrodigal.py --mode meta
#   2) current run_pyrodigal.py --mode single
#   3) Bakta-equivalent Pyrodigal default trained mode, circular/open edge
#   4) Bakta-equivalent Pyrodigal --meta mode, circular/open edge
#
# The Bakta-equivalent tests call Pyrodigal from the Bakta environment
# directly, so CDS prediction does not require AMRFinderPlus or its database.
#
# Usage:
#   bash projects/config_mag_annotation/scripts/test_u999856_gene_calling.sh
#   bash projects/config_mag_annotation/scripts/test_u999856_gene_calling.sh /path/to/output

GENOME_ID="${GENOME_ID:-u999856}"
THREADS="${THREADS:-8}"
OUT_DIR="${1:-$PWD/work/u999856_gene_calling_test}"

FA="${FA:-/dssg/home/acct-clswfp/clswfp-icdli1/RenTianxin/results/binning/Bathy_enrichment_2025/circular_contigs/u999856.fa}"
BAKTA_DB="${BAKTA_DB:-/dssg/home/acct-clswfp/clswfp-icdli1/Database/Bakta/db-light}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
RUN_PYRO="$REPO_ROOT/modules/mag_annotation/scripts/run_pyrodigal.py"

CONTIG="u999856ctg_len-2079472_circular-yes_depth-2443-2393-2315_duplicated-no"
ENV_ROOTS=(
  "/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs"
  "/dssg/home/acct-clswfp/clswfp-icdli1/software/miniconda3/envs"
)

resolve_bin() {
  local env_name="$1"
  local bin_name="$2"
  local root
  for root in "${ENV_ROOTS[@]}"; do
    if [[ -x "$root/$env_name/bin/$bin_name" ]]; then
      printf '%s\n' "$root/$env_name/bin/$bin_name"
      return 0
    fi
  done
  if command -v "$bin_name" >/dev/null 2>&1; then
    command -v "$bin_name"
    return 0
  fi
  echo "ERROR: cannot find $bin_name in environment $env_name" >&2
  return 1
}

PYTHON="${PYTHON:-$(resolve_bin pyrodigal_3.7.1 python)}"
BAKTA_PYTHON="${BAKTA_PYTHON:-$(resolve_bin bakta_1.12.0 python)}"
BAKTA="${BAKTA:-$(resolve_bin bakta_1.12.0 bakta || true)}"

[[ -r "$FA" ]] || { echo "ERROR: genome FASTA not readable: $FA" >&2; exit 1; }
[[ -r "$RUN_PYRO" ]] || { echo "ERROR: run_pyrodigal.py not found: $RUN_PYRO" >&2; exit 1; }

mkdir -p "$OUT_DIR" "$OUT_DIR/current_meta" "$OUT_DIR/current_single" \
  "$OUT_DIR/bakta_cds_default" "$OUT_DIR/bakta_cds_meta"

printf 'original sequence id\tnew sequence id\ttype\ttopology\tname\n%s\t-\tchromosome\tcircular\t-\n' \
  "$CONTIG" > "$OUT_DIR/u999856.replicons.tsv"

echo "[CONFIG]"
echo "FA=$FA"
echo "BAKTA_DB=$BAKTA_DB"
echo "PYTHON=$PYTHON"
echo "BAKTA_PYTHON=$BAKTA_PYTHON"
echo "BAKTA=${BAKTA:-NOT_FOUND}"
echo "OUT_DIR=$OUT_DIR"
"$PYTHON" -c 'import pyrodigal; print("pyrodigal=" + str(getattr(pyrodigal, "__version__", "unknown")))'
"$BAKTA_PYTHON" -c 'import pyrodigal; print("bakta_env_pyrodigal=" + str(getattr(pyrodigal, "__version__", "unknown")))'
if [[ -n "$BAKTA" && -x "$BAKTA" ]]; then
  "$BAKTA" --version
else
  echo "bakta_executable=NOT_USED_FOR_CDS_ONLY_TEST"
fi

echo "[1/4] current run_pyrodigal: meta"
"$PYTHON" "$RUN_PYRO" \
  -i "$FA" -g "$GENOME_ID" -m meta -t 11 \
  --out_faa "$OUT_DIR/current_meta/u999856.faa.gz" \
  --out_gff "$OUT_DIR/current_meta/u999856.gff.gz" \
  --out_meta "$OUT_DIR/current_meta/u999856_protein_metadata.tsv.gz" \
  2>&1 | tee "$OUT_DIR/current_meta/run.log"

echo "[2/4] current run_pyrodigal: single"
"$PYTHON" "$RUN_PYRO" \
  -i "$FA" -g "$GENOME_ID" -m single -t 11 \
  --out_faa "$OUT_DIR/current_single/u999856.faa.gz" \
  --out_gff "$OUT_DIR/current_single/u999856.gff.gz" \
  --out_meta "$OUT_DIR/current_single/u999856_protein_metadata.tsv.gz" \
  2>&1 | tee "$OUT_DIR/current_single/run.log"

run_bakta_cds_only() {
  local mode="$1"
  local out_dir="$2"
  mkdir -p "$out_dir"
  echo "mode=$mode; closed=False; mask=True; environment=bakta_1.12.0" | tee "$out_dir/run.log"

  "$BAKTA_PYTHON" - "$FA" "$CONTIG" "$out_dir" "$mode" <<'PY' 2>&1 | tee -a "$out_dir/run.log"
import gzip
import os
import sys

import pyrodigal

fasta, target_contig, out_dir, mode = sys.argv[1:]

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

seqs = {}
name = None
buf = []
with open_text(fasta) as handle:
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
if name is not None:
    seqs[name] = "".join(buf).upper()

if target_contig not in seqs:
    raise SystemExit(f"target contig not found: {target_contig}")

# Bakta's CDS stage uses Pyrodigal with mask=True. For a complete circular
# replicon Bakta deliberately leaves the sequence open to preserve edge CDSs.
if mode == "meta":
    finder = pyrodigal.GeneFinder(
        meta=True,
        metagenomic_bins=None,
        closed=False,
        mask=True,
    )
else:
    finder = pyrodigal.GeneFinder(
        meta=False,
        closed=False,
        mask=True,
    )
    finder.train(
        *[seq.encode("ascii") for seq in seqs.values()],
        translation_table=11,
    )

genes = finder.find_genes(seqs[target_contig].encode("ascii"))
faa_path = os.path.join(out_dir, "u999856.faa")
gff_path = os.path.join(out_dir, "u999856.gff3")
meta_path = os.path.join(out_dir, "u999856_protein_metadata.tsv")

def esc(value):
    return str(value).replace("%", "%25").replace(";", "%3B").replace("=", "%3D")

with open(faa_path, "w") as faa, open(gff_path, "w") as gff, open(meta_path, "w") as meta:
    gff.write("##gff-version 3\n")
    meta.write("protein_id\tstart\tend\tstrand\taa_length\tstart_type\trbs_motif\tpartial_begin\tpartial_end\tscore\n")
    for index, gene in enumerate(genes, start=1):
        protein_id = f"u999856___{target_contig}___bakta_cds_{index:05d}"
        translation = gene.translate(
            translation_table=None if mode == "meta" else 11
        ).rstrip("*").upper()
        strand = "+" if gene.strand == 1 else "-"
        attrs = {
            "ID": protein_id,
            "protein_id": protein_id,
            "start_type": getattr(gene, "start_type", "-"),
            "rbs_motif": getattr(gene, "rbs_motif", "-"),
        }
        attr_text = ";".join(f"{key}={esc(value)}" for key, value in attrs.items())
        gff.write(
            f"{target_contig}\tpyrodigal_bakta_mode\tCDS\t{gene.begin}\t{gene.end}\t{gene.score:.3f}\t"
            f"{strand}\t0\t{attr_text}\n"
        )
        faa.write(f">{protein_id}\n{translation}\n")
        meta.write(
            f"{protein_id}\t{gene.begin}\t{gene.end}\t{strand}\t{len(translation)}\t"
            f"{getattr(gene, 'start_type', '-') }\t{getattr(gene, 'rbs_motif', '-') }\t"
            f"{getattr(gene, 'partial_begin', False)}\t{getattr(gene, 'partial_end', False)}\t{gene.score:.3f}\n"
        )

print(f"mode={mode}")
print(f"predicted={len(genes)}")
print(f"faa={faa_path}")
print(f"gff={gff_path}")
print(f"metadata={meta_path}")
PY
}

echo "[3/4] Bakta CDS engine equivalent: default trained mode, circular/open edge"
run_bakta_cds_only default "$OUT_DIR/bakta_cds_default"

echo "[4/4] Bakta CDS engine equivalent: meta mode, circular/open edge"
run_bakta_cds_only meta "$OUT_DIR/bakta_cds_meta"

echo "[SUMMARY] CDSs overlapping the prot01279/prot01280 test region"
python3 - "$OUT_DIR" "$CONTIG" <<'PY'
import gzip
import os
import sys

out, contig = sys.argv[1:]
region_start, region_end = 1284830, 1286696

sources = {
    "current_meta": os.path.join(out, "current_meta", "u999856.gff.gz"),
    "current_single": os.path.join(out, "current_single", "u999856.gff.gz"),
    "bakta_cds_default": os.path.join(out, "bakta_cds_default", "u999856.gff3"),
    "bakta_cds_meta": os.path.join(out, "bakta_cds_meta", "u999856.gff3"),
}

def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def attrs(text):
    result = {}
    for item in text.rstrip(";\n").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result

print("source\tstart\tend\tstrand\taa_no_stop\tfeature_id\tstart_type\trbs_motif")
for source, path in sources.items():
    if not os.path.exists(path):
        print(f"{source}\tMISSING")
        continue
    with opener(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2].lower() != "cds":
                continue
            if fields[0] != contig:
                continue
            start, end = int(fields[3]), int(fields[4])
            if end < region_start or start > region_end:
                continue
            a = attrs(fields[8])
            feature_id = a.get("protein_id") or a.get("locus_tag") or a.get("ID", "-")
            start_type = a.get("start_type", "-")
            rbs_motif = a.get("rbs_motif", "-")
            nt_len = end - start + 1
            aa_len = nt_len // 3 - 1 if nt_len % 3 == 0 else nt_len // 3
            print(f"{source}\t{start}\t{end}\t{fields[6]}\t{aa_len}\t{feature_id}\t{start_type}\t{rbs_motif}")
PY

echo "[DONE] Results are in: $OUT_DIR"
