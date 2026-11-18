#!/usr/bin/env python3
"""Download versioned public genomes; preserve bytes and record SHA256 provenance."""
import hashlib
import gzip
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENOMES = {
    "ecoli_MG1655": "NC_000913.3",
    "msmegmatis_mc2155": "NC_008596.1",
    "swolfei_Goettingen": "NC_008346.1",
    "tmaritima_MSB8": "NC_000853.1",
    "mmarburgensis_Marburg": "NC_014408.1",
    "mbarkeri_Fusaro": "NC_007355.1",
    "pfuriosus_DSM3638": "NC_003413.1",
    "ecoli_W3110_holdout": "NC_007779.1",
    "dvulgaris_Hildenborough": "NC_002937.3",
    "wsuccinogenes_DSM1740": "NC_005090.1",
    "npunctiforme_ATCC29133": "NC_010628.1",
    "cnecator_H16_pHG1": "NC_005241.1",
    "rcapsulatus_SB1003": "NC_014034.1",
    "tonnurineus_NA1": "NC_011529.1",
    "chydrogenoformans_Z2901": "NC_007503.1",
    "cpasteurianum_DSM525": "CP009268.1",
    "awoodii_DSM1030": "NC_016894.1",
    "cautoethanogenum_DSM10061": "NC_022592.1",
    "creinhardtii_HydA1_cDNA": "AF289201.1",
}


def fetch(url, path, required):
    if path.exists() and required in path.read_bytes():
        return
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "hydrogenase-reference-workflow/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response:
                data = response.read()
            if required not in data:
                raise ValueError(f"Unexpected response from {url}")
            temporary = path.with_suffix(path.suffix + ".part")
            temporary.write_bytes(data)
            temporary.replace(path)
            return
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))


def main():
    raw = ROOT / "references/raw"
    raw.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, accession in GENOMES.items():
        params = urllib.parse.urlencode(dict(db="nuccore", id=accession, rettype="gbwithparts", retmode="text"))
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
        path = raw / f"{name}.gb"
        fetch(url, path, b"LOCUS ")
        version=next(line.split()[1] for line in path.read_text().splitlines() if line.startswith("VERSION"))
        if version!=accession:
            raise ValueError(f"Requested {accession}, received {version}")
        manifest.append(dict(sample=name, accession=accession, url=url,
                             bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                             downloaded_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(path.stat().st_mtime))))
        print(name, accession, path.stat().st_size, flush=True)
        time.sleep(0.4)
    (raw / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    experimental=[]
    for accession in ["P29166","Q9FYU1"]:
        url=f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        path=raw/f"uniprot_{accession}.json"
        fetch(url,path,b'"primaryAccession"')
        if json.loads(path.read_text())["primaryAccession"]!=accession:
            raise ValueError(f"Unexpected UniProt accession for {accession}")
        experimental.append(dict(accession=accession,url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (raw/"experimental_sequence_manifest.json").write_text(json.dumps(experimental,indent=2)+"\n")
    policy=json.loads((ROOT/"config/maturation_domain_policy.json").read_text())
    domain_dir=raw/"domain_review"
    domain_dir.mkdir(exist_ok=True)
    for model in policy["models"]:
        path=domain_dir/model["file"]
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==model["sha256"]:
            continue
        with urllib.request.urlopen(model["url"],timeout=45) as response:
            data=response.read()
        if data[:2]==b'\x1f\x8b':
            data=gzip.decompress(data)
        if hashlib.sha256(data).hexdigest()!=model["sha256"] or not data.startswith(b"HMMER"):
            raise ValueError(f"Pinned domain model changed: {model['accession']}; review its version before updating")
        temporary=path.with_suffix(".hmm.part")
        temporary.write_bytes(data)
        temporary.replace(path)
    structure_policy=json.loads((ROOT/"config/structure_sources.json").read_text())
    structures=raw/"structures"
    structures.mkdir(exist_ok=True)
    for source in structure_policy["sources"]:
        path=structures/(source["pdb_id"]+".cif")
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==source["sha256"]:
            continue
        with urllib.request.urlopen(source["download_url"],timeout=60) as response:
            data=response.read()
        if hashlib.sha256(data).hexdigest()!=source["sha256"]:
            raise ValueError(f"Pinned structural source changed: {source['pdb_id']}; review the revision before updating")
        temporary=path.with_suffix(".cif.part")
        temporary.write_bytes(data)
        temporary.replace(path)
    # Preserve an original download record; source-policy hashes also permit fresh setup.
    if not (structures/"download_manifest.json").exists():
        (structures/"download_manifest.json").write_text(json.dumps(structure_policy["sources"],indent=2)+"\n")
    split_policy=json.loads((ROOT/'config/split_centers.json').read_text())
    for name,source in split_policy['source_files'].items():
        path=ROOT/name
        if not path.resolve().is_relative_to(raw.resolve()):raise ValueError('Split-center download path is outside raw sources')
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256']:continue
        with urllib.request.urlopen(source['url'],timeout=60) as response:data=response.read()
        if hashlib.sha256(data).hexdigest()!=source['sha256']:raise ValueError('Pinned split-center source changed: '+name)
        path.parent.mkdir(parents=True,exist_ok=True);temporary=path.with_suffix(path.suffix+'.part');temporary.write_bytes(data);temporary.replace(path)


if __name__ == "__main__":
    main()
