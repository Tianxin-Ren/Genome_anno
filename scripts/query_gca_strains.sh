#!/usr/bin/env bash
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  query_gca_strains.sh -i INPUT [options]

Options:
  -i, --input FILE       TSV containing GCA/GCF accessions; required
  -o, --output-dir DIR   Output directory; default: current directory
  -e, --env DIR          ncbi-datasets environment; default: icdli1 path
  -c, --chunk-size N     Accessions per NCBI request; default: 200
      --reuse-summary FILE  Reprocess an existing JSONL without querying NCBI
  -h, --help             Show this help
EOF
}

INPUT=""
OUT_DIR="$(pwd)"
ENV_DIR="/dssg/home/acct-clswfp/clswfp-icdli1/.conda/envs/ncbi-datasets-cli_18.26.0"
CHUNK_SIZE=200
REUSE_SUMMARY=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input|-o|--output-dir|-e|--env|-c|--chunk-size|--reuse-summary)
            if [[ $# -lt 2 || -z "$2" || "$2" == -* ]]; then
                echo "Missing value for $1" >&2; exit 2
            fi ;;
    esac
    case "$1" in
        -i|--input) INPUT="${2:-}"; shift 2 ;;
        -o|--output-dir) OUT_DIR="${2:-}"; shift 2 ;;
        -e|--env) ENV_DIR="${2:-}"; shift 2 ;;
        -c|--chunk-size) CHUNK_SIZE="${2:-}"; shift 2 ;;
        --reuse-summary) REUSE_SUMMARY="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$INPUT" || ! -f "$INPUT" ]]; then
    echo "Input file is required and must exist: $INPUT" >&2
    usage >&2
    exit 2
fi
if ! [[ "$CHUNK_SIZE" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid chunk size: $CHUNK_SIZE" >&2
    exit 2
fi

mkdir -p "$OUT_DIR"
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
OUT_DIR="$(cd "$OUT_DIR" && pwd)"
if [[ -n "$REUSE_SUMMARY" ]]; then
    [[ -s "$REUSE_SUMMARY" ]] || { echo "JSONL missing or empty: $REUSE_SUMMARY" >&2; exit 2; }
    REUSE_SUMMARY="$(cd "$(dirname "$REUSE_SUMMARY")" && pwd)/$(basename "$REUSE_SUMMARY")"
fi

DATASETS="$ENV_DIR/bin/datasets"
DATAFORMAT="$ENV_DIR/bin/dataformat"
if [[ ! -x "$DATAFORMAT" || ( -z "$REUSE_SUMMARY" && ! -x "$DATASETS" ) ]]; then
    echo "datasets/dataformat not found in: $ENV_DIR" >&2
    exit 1
fi
command -v python >/dev/null || { echo "python is required" >&2; exit 1; }

cd "$OUT_DIR"
RUN_DIR="gca_query_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN_DIR/chunks"

# Header rows fail the accession pattern; headerless lists keep their first entry.
LC_ALL=C awk -F'\t' '
{
    for (i = 1; i <= NF; i++) {
        x = $i
        gsub(/\357\273\277/, "", x)
        sub(/^[[:space:]]+/, "", x)
        sub(/[[:space:]]+$/, "", x)
        gsub(/^(RS_|GB_)/, "", x)
        if (x ~ /^GC[AF]_[0-9]+(\.[0-9]+)?$/) print x
    }
}' "$INPUT" | sort -u > "$RUN_DIR/accessions.txt"
cp "$RUN_DIR/accessions.txt" gca_input.txt

if [[ ! -s gca_input.txt ]]; then
    echo "No GCA/GCF accession found in: $INPUT" >&2
    exit 1
fi

if [[ -n "$REUSE_SUMMARY" ]]; then
    if [[ ! gca_summary.jsonl -ef "$REUSE_SUMMARY" ]]; then
        cp "$REUSE_SUMMARY" gca_summary.jsonl
    fi
else
    split -l "$CHUNK_SIZE" gca_input.txt "$RUN_DIR/chunks/chunk_"
    : > gca_summary.jsonl
    : > gca_query.err
    for chunk in "$RUN_DIR"/chunks/chunk_*; do
        "$DATASETS" summary genome accession \
            --inputfile "$chunk" \
            --as-json-lines \
            >> gca_summary.jsonl \
            2>> gca_query.err
    done
fi
[[ -s gca_summary.jsonl ]] || { echo "No JSONL returned; see gca_query.err" >&2; exit 1; }

"$DATAFORMAT" tsv genome \
    --fields accession,organism-name,organism-infraspecific-strain,organism-infraspecific-isolate,assminfo-biosample-accession,assminfo-biosample-strain,assminfo-biosample-isolate,assminfo-biosample-isolate-name-alias,assminfo-biosample-sample-name,assminfo-biosample-ids-label,assminfo-biosample-ids-value,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value,assminfo-biosample-description-comment,assminfo-biosample-isolation-source \
    < gca_summary.jsonl > gca_strain_metadata.tsv

python - <<'PY' | tee gca_statistics.log
import csv
import json
import re
from collections import defaultdict

BAD = {"", "missing", "not applicable", "not available", "unknown", "none", "na", "n/a", ".", "-"}
LABEL_RE = re.compile(r"^(strain|isolate|isolate_name_alias|culture_collection|culture collection)$", re.I)
COLLECTION_RE = re.compile(
    r"(?<![A-Z0-9])(?P<code>ATCC|DSMZ?|JCM|CGMCC|NCTC|CIP|LMG|NBRC|KCTC|"
    r"BCRC|CCUG|NRRL|MTCC|VKPM|NCIMB|NCIM|NCIB|IFO|CCTCC|CICC|CBS|IAM|"
    r"CCAP|CCMP|UTEX|BCCM)[\s:_-]*(?P<id>(?:[A-Z]+[: _-]*)*\d[A-Z0-9_.:/-]*)",
    re.I,
)

def clean(value):
    value = re.sub(r"\s+", " ", str(value or "")).strip()
    return "" if value.lower() in BAD else value.upper()

def first(row, *names):
    for name in names:
        if name in row:
            return row[name]
    return ""

records = defaultdict(lambda: {"organism": "", "ids": set(), "collections": set()})
with open("gca_input.txt", encoding="utf-8") as handle:
    input_accs = {line.strip() for line in handle if line.strip()}

with open("gca_strain_metadata.tsv", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        acc = clean(row.get("Assembly Accession"))
        if acc not in input_accs:
            continue
        records[acc]["organism"] = clean(row.get("Organism Name"))

        for field in (
            "Organism Infraspecific Names Strain",
            "Organism Infraspecific Names Isolate",
            "Assembly BioSample Strain",
            "Assembly BioSample Isolate",
            "Assembly BioSample Isolate name alias",
        ):
            value = clean(row.get(field))
            if value:
                records[acc]["ids"].add(value)

        label = clean(first(row, "Assembly BioSample Sample Identifiers Label", "Assembly BioSample IDs Label"))
        value = clean(first(row, "Assembly BioSample Sample Identifiers Value", "Assembly BioSample IDs Value"))
        if value and LABEL_RE.fullmatch(label):
            records[acc]["ids"].add(value)

        attr = clean(row.get("Assembly BioSample Attribute Name"))
        value = clean(row.get("Assembly BioSample Attribute Value"))
        if value and LABEL_RE.fullmatch(attr):
            records[acc]["ids"].add(value)


# Scan original fields independently, retaining evidence lost by TSV projection.
def text_fields(obj, path="", label=""):
    if isinstance(obj, dict):
        attribute = obj.get("name") or obj.get("harmonizedName") or obj.get("attributeName")
        for key, value in obj.items():
            yield from text_fields(value, f"{path}.{key}".lstrip("."),
                                   attribute if key == "value" and attribute else key)
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            yield from text_fields(value, f"{path}[{index}]", label)
    elif isinstance(obj, str):
        yield path, re.sub(r"[^a-z0-9]", "", label.lower()), obj

evidence = set()
with open("gca_summary.jsonl", encoding="utf-8") as handle:
    for line in handle:
        if not line.strip():
            continue
        report = json.loads(line)
        acc = report.get("accession", "")
        if not acc:
            raise ValueError("JSONL record lacks accession")
        if acc not in input_accs:
            continue
        rec = records[acc]
        for path, label, text in text_fields(report):
            value = clean(text)
            if not value:
                continue
            if label in {"strain", "isolate", "isolatenamealias"}:
                rec["ids"].add(value)
            if label == "culturecollection":
                rec["collections"].add(value)
                evidence.add((acc, "declared_culture_collection", value, path, text))
            else:
                # shortcut: free-text matches are candidates; verify them in collection catalogues.
                for match in COLLECTION_RE.finditer(text):
                    code = match["code"].upper()
                    if code == "DSMZ":
                        code = "DSM"
                    identifier = re.sub(r"\s+", "", match["id"].upper()).rstrip(".")
                    candidate = f"{code}:{identifier}"
                    rec["collections"].add(candidate)
                    evidence.add((acc, "text_candidate", candidate, path, text))

with open("gca_collection_evidence.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession", "evidence_type", "collection_id_or_value", "source_field", "source_text"])
    writer.writerows(sorted(evidence))

with open("gca_strain_presence.tsv", "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["accession", "organism", "strain_or_isolate_ids", "culture_collection_ids"])
    for acc in sorted(records):
        writer.writerow([
            acc,
            records[acc]["organism"],
            ";".join(sorted(records[acc]["ids"])) or "NA",
            ";".join(sorted(records[acc]["collections"])) or "NA",
        ])

owners = defaultdict(set)
for acc, record in records.items():
    for strain_id in record["ids"]:
        owners[strain_id].add(acc)

with open("duplicated_strain_ids.tsv", "w", encoding="utf-8") as handle:
    handle.write("strain_id\tgca_count\tgca_accessions\n")
    for strain_id, gcas in sorted(owners.items()):
        if len(gcas) > 1:
            handle.write(f"{strain_id}\t{len(gcas)}\t{';'.join(sorted(gcas))}\n")

input_count = len(input_accs)

with_id = sum(bool(record["ids"]) for record in records.values())
unique_ids = len(owners)
single_gca_ids = sum(len(gcas) == 1 for gcas in owners.values())
with_collection = sum(bool(record["collections"]) for record in records.values())

print(f"Input GCA/GCF: {input_count}")
print(f"Returned metadata GCA: {len(records)}")
print(f"GCA without returned metadata: {len(input_accs - records.keys())}")
print(f"GCA with >=1 strain/isolate ID: {with_id}")
print(f"Returned GCA without strain/isolate ID: {len(records) - with_id}")
print(f"Distinct strain/isolate IDs: {unique_ids}")
print(f"IDs belonging to exactly one GCA: {single_gca_ids}")
print(f"GCA with a culture-collection-like ID: {with_collection}")
print(f"GCA with declared culture_collection: {len({acc for acc, kind, *_ in evidence if kind == 'declared_culture_collection' and acc in input_accs})}")
print("Distinct IDs are text identifiers, not a count of biological strains.")
print("Wrote: gca_strain_metadata.tsv, gca_strain_presence.tsv, duplicated_strain_ids.tsv, gca_collection_evidence.tsv")
PY

echo "Query run: $OUT_DIR/$RUN_DIR"
