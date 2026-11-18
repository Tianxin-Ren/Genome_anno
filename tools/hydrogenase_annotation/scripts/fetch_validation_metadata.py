#!/usr/bin/env python3
"""Cache public accession/taxonomy metadata; product names are not activity truth."""
import csv
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio import SeqIO
from annotate import ROOT


RAW = ROOT / "references/raw/validation_corpus"


def accession(value):
    match = re.match(r"^([A-Z]+_?\d+\.\d+)(?:$|[ _|])", value.strip())
    if not match:
        raise ValueError(f"Invalid source protein accession: {value}")
    return match[1]


def main():
    output = RAW / "ncbi_metadata"
    output.mkdir(exist_ok=True)
    manifest_file = output / "download_manifest.json"
    manifest = json.loads(manifest_file.read_text()) if manifest_file.exists() else {}

    def fetch(name, db, operation, ids, fmt):
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + operation + ".fcgi?" + urllib.parse.urlencode(
            dict(db=db, id=",".join(map(str, ids)), retmode=fmt, tool="hydrogenase_validation_audit"))
        path = output / name
        if path.exists():
            entry = manifest[name]
            if entry["url"] != url or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
                raise ValueError(f"Changed cache: {path}")
            return path.read_bytes()
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=90) as response:
                    data = response.read()
                if fmt == "json":
                    result = json.loads(data)
                    if "error" in result:
                        raise ValueError(result["error"])
                else:
                    result = ET.fromstring(data)
                    if result.tag != "TaxaSet" or result.find("ERROR") is not None:
                        raise ValueError("NCBI taxonomy did not return a TaxaSet")
                path.write_bytes(data)
                manifest[name] = dict(url=url, sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
                manifest_file.write_text(json.dumps(manifest, indent=2) + "\n")
                print(name, len(data), flush=True)
                time.sleep(0.4)
                return data
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)

    ids = {accession(r["ncbi_accession"]) for r in csv.DictReader((RAW / "HydDB_original_dataset.csv").open(), delimiter=";")}
    ids.update(accession(r["accession"]) for r in csv.DictReader((RAW / "HydDB_original_non_hydrogenase_references.csv").open(), delimiter=";"))
    ids.update(accession(r.id) for r in SeqIO.parse(RAW / "HydDB_all_hydrogenases.faa", "fasta"))
    refs = json.loads((ROOT / "references/processed/references.json").read_text())
    ids.update(accession(r["protein_accession"]) for r in refs)
    ids = sorted(ids)
    protein_rows = []
    for index in range(0, len(ids), 100):
        result = json.loads(fetch(f"protein_{index // 100:03}.json", "protein", "esummary", ids[index:index+100], "json"))["result"]
        protein_rows.extend(result[uid] for uid in result["uids"])
    taxids = {int(r["taxid"]) for r in protein_rows if r.get("taxid")}
    for path in sorted(RAW.glob("reviewed_*.json")):
        if not path.name.endswith(".headers.json"):
            taxids.update(r["organism"]["taxonId"] for r in json.loads(path.read_text())["results"])
    for sample in {r["sample"] for r in refs}:
        rec = SeqIO.read(ROOT / "references/raw" / (sample + ".gb"), "genbank")
        for feature in rec.features:
            if feature.type == "source":
                taxids.update(int(v.split(":")[1]) for v in feature.qualifiers.get("db_xref", []) if v.startswith("taxon:"))
    # Previously examined development structures/cases cannot become final tests through a new name.
    for path in sorted((ROOT / "references/raw/structures").glob("*.cif")):
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        cif = MMCIF2Dict(str(path))
        for field in ["_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id", "_entity_src_nat.pdbx_ncbi_taxonomy_id"]:
            taxids.update(int(v) for v in cif.get(field, []) if v.isdigit())
    extra = ROOT / "references/raw/independent/tkivui_DSM2030.gb"
    if extra.exists():
        rec = SeqIO.read(extra, "genbank")
        for feature in rec.features:
            if feature.type == "source":
                taxids.update(int(v.split(":")[1]) for v in feature.qualifiers.get("db_xref", []) if v.startswith("taxon:"))
    taxa = sorted(taxids)
    for index in range(0, len(taxa), 100):
        fetch(f"taxonomy_{index // 100:03}.xml", "taxonomy", "efetch", taxa[index:index+100], "xml")
    retrieved = {r.get("accessionversion") for r in protein_rows}
    report = dict(requested_proteins=len(ids), retrieved_proteins=len(protein_rows), taxonomy_ids=len(taxa),
                  missing_accessions=sorted(set(ids)-retrieved), limitation="NCBI taxonomy and products are metadata, not catalytic/function experimental evidence.")
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    assert accession("WP_005026529.1_Bilophila_wadsworthia|[FeFe]_Group_A1") == "WP_005026529.1"
    assert accession("WP_005026529.1 ") == "WP_005026529.1"
    main()
