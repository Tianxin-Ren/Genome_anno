#!/usr/bin/env python3
"""Freeze source-labelled catalytic identity partitions before examining predictions."""
import argparse
import collections
import csv
import hashlib
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, write_tsv
from fetch_validation_metadata import RAW, accession


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()


def genus_ownership(genus, policy):
    value = int(digest(policy["genus_split_salt"] + ":" + genus), 16) % 100
    return next(key for key, (low, high) in policy["genus_hash_percentiles"].items() if low <= value <= high)


def component_split(genera, exposed, valid, policy):
    if exposed:
        return "reference_overlap"
    if not valid or not genera:
        return "source_review"
    owners = {genus_ownership(g, policy) for g in genera}
    return owners.pop() if len(owners) == 1 else "split_bridge"


def curated_role(entry):
    name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    domains = {x["id"] for x in entry.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam"}
    lower = name.lower()
    if lower.startswith("nadh-quinone oxidoreductase subunit"):
        return "negative", "Nuo", "reviewed respiratory complex I subunit D-containing product"
    if any(v in lower for v in ["nuclear prelamin a recognition", "cytosolic iron-sulfur assembly component 3", "cytosolic fe-s cluster assembly", "protein nar1"]):
        return "negative", "NARF", "reviewed Fe-S assembly / NARF homolog, not catalytic hydrogenase"
    if "hydrogenase expression/formation protein" in lower:
        return "negative", "HupK_HoxV", "reviewed maturation/chaperone hydrogenase homolog"
    if lower == "f(420)h(2) dehydrogenase subunit d":
        return "negative", "FpoD", "reviewed F420H2 dehydrogenase, a non-hydrogenase catalytic-core homolog"
    if "methylenetetrahydromethanopterin dehydrogenase-related" in lower:
        return "excluded", "HmdII", "Hmd paralog name does not establish absence of all H2 catalysis; exact experimental review required"
    if "hydrogenase" in lower and not any(v in lower for v in ["uncharacterized", "probable", "formation", "expression", "related"]):
        if "PF00374" in domains:
            return "positive", "NiFe", "reviewed hydrogenase catalytic-chain product plus NiFe core domain"
        if "PF03201" in domains and lower == "5,10-methenyltetrahydromethanopterin hydrogenase":
            return "positive", "Fe", "reviewed Hmd catalytic product plus Hmd domain"
        if {"PF02906", "PF02256"} <= domains:
            return "positive", "FeFe", "reviewed FeFe catalytic-chain product plus H-domain regions"
    return "excluded", "", "keyword/domain search result is insufficient catalytic identity truth"


def taxonomy():
    result = {}
    for path in sorted((RAW / "ncbi_metadata").glob("taxonomy_*.xml")):
        for node in ET.fromstring(path.read_bytes()).findall("Taxon"):
            lineage = [*node.findall("LineageEx/Taxon"), node]
            genus = next((n for n in lineage if n.findtext("Rank") == "genus"), None)
            row = dict(taxid=node.findtext("TaxId"), name=node.findtext("ScientificName"), rank=node.findtext("Rank"),
                       genus_taxid=genus.findtext("TaxId") if genus is not None else "",
                       genus_name=genus.findtext("ScientificName") if genus is not None else "")
            for key in [row["taxid"], *[n.text for n in node.findall("AkaTaxIds/TaxId")]]:
                result[key] = row
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if (out / "freeze.json").exists():
        parser.error("Frozen partitions must not be overwritten; use a new directory")
    out.mkdir(parents=True, exist_ok=True)
    policy_file = ROOT / "config/identity_corpus_policy.json"
    policy = json.loads(policy_file.read_text())
    taxa = taxonomy()
    if not taxa or not (RAW / "ncbi_metadata/summary.json").exists():
        raise ValueError("Complete the pinned NCBI metadata download first")
    ncbi = {}
    for path in sorted((RAW / "ncbi_metadata").glob("protein_*.json")):
        data = json.loads(path.read_text())["result"]
        ncbi.update({data[uid]["accessionversion"]: data[uid] for uid in data["uids"] if "accessionversion" in data[uid]})
    units, audit, variants = {}, [], []

    def add(sequence, source, acc, label, category, tier, organism, taxid, raw_file, review_note):
        normalized = sequence.strip().upper()
        if not normalized or re.search(r"[^ACDEFGHIKLMNPQRSTVWYBXZJUO]", normalized):
            raise ValueError(f"Non-protein letters: {source}: {acc}")
        key = digest(normalized)
        record = units.setdefault(key, dict(unit_id="seq_"+key[:20], sequence=normalized, sequence_sha256=key, sources=[]))
        row = dict(source=source, accession=acc, label=label, category=category, evidence_tier=tier, organism=organism,
                   taxid=str(taxid or ""), source_file=raw_file, source_sha256=sha(RAW/raw_file),
                   raw_sequence_sha256=digest(sequence), normalized_case=sequence != normalized, review_note=review_note)
        record["sources"].append(row)
        audit.append(dict(unit_id=record["unit_id"], **row))
        return record["unit_id"]

    original = {}
    for row in csv.DictReader((RAW / "HydDB_original_dataset.csv").open(), delimiter=";"):
        acc = accession(row["ncbi_accession"])
        original[acc] = row
        add(row["protein_sequence"], "HydDB_original", acc, "positive", row["hydrogenase"], "author_curated_catalytic_identity",
            row["organism"], ncbi.get(acc, {}).get("taxid"), "HydDB_original_dataset.csv",
            "Source classification, not per-protein experimental activity or a confirmed molecular reaction")
    for row in csv.DictReader((RAW / "HydDB_original_non_hydrogenase_references.csv").open(), delimiter=";"):
        acc = accession(row["accession"])
        add(row["protein_sequence"], "HydDB_negative", acc, "negative", row["domain"], "author_curated_non_hydrogenase_homolog",
            row["organism"], ncbi.get(acc, {}).get("taxid"), "HydDB_original_non_hydrogenase_references.csv",
            "Negative catalytic identity; subtype label is author curation, not new biochemical measurement")
    for row in SeqIO.parse(RAW / "HydDB_all_hydrogenases.faa", "fasta"):
        acc, category, seq = accession(row.id), row.description.split("|")[-1].replace("_", " "), str(row.seq).upper()
        old = original.get(acc)
        oldseq = old["protein_sequence"].strip().upper() if old else ""
        status = "identical" if seq == oldseq else "exact_subsequence" if oldseq and seq in oldseq else "different_sequence" if old else "new_accession"
        variants.append(dict(accession=acc, raw_header=row.description, category=category, sequence_sha256=digest(seq), length=len(seq),
            relationship=status, original_unit_id="seq_"+digest(oldseq)[:20] if old else "", original_length=len(oldseq),
            original_start=oldseq.find(seq)+1 if seq in oldseq else "", class_conflict=bool(old and old["hydrogenase"] != category)))
        if not old:
            meta = ncbi.get(acc, {})
            add(seq, "HydDB_current_new", acc, "positive", category, "author_curated_catalytic_identity", meta.get("organism", ""),
                meta.get("taxid"), "HydDB_all_hydrogenases.faa", "New author-curated accession; source sequence scope retained")
    reviewed = {}
    for path in sorted(RAW.glob("reviewed_*.json")):
        if not path.name.endswith(".headers.json"):
            for entry in json.loads(path.read_text())["results"]:
                reviewed.setdefault(entry["primaryAccession"], (path.name, entry))
    for acc, (filename, entry) in sorted(reviewed.items()):
        label, category, note = curated_role(entry)
        name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        if label == "excluded":
            audit.append(dict(unit_id="", source="UniProt", accession=acc, label=label, category=category,
                evidence_tier="not_catalytic_truth", organism=entry["organism"]["scientificName"], taxid=entry["organism"]["taxonId"],
                source_file=filename, source_sha256=sha(RAW/filename), raw_sequence_sha256=digest(entry["sequence"]["value"]),
                normalized_case=False, review_note=note+": "+name))
            continue
        experimental = sorted({e.get("id", "") for c in entry.get("comments", []) if c["commentType"] == "CATALYTIC ACTIVITY"
            for e in c.get("reaction", {}).get("evidences", []) if e["evidenceCode"] == "ECO:0000269"})
        tier = "reviewed_catalytic_identity_with_direct_reaction_evidence" if label == "positive" and experimental else "reviewed_catalytic_identity" if label == "positive" else "reviewed_non_hydrogenase_identity"
        add(entry["sequence"]["value"], "UniProt", acc, label, category, tier, entry["organism"]["scientificName"],
            entry["organism"]["taxonId"], filename, note+": "+name+"; reaction_ECO0000269_PubMed="+",".join(experimental))
    rows = sorted(units.values(), key=lambda r: r["unit_id"])
    if len({r["unit_id"] for r in rows}) != len(rows):
        raise ValueError("Truncated sequence digest collision")
    for row in rows:
        sources = row["sources"]
        labels = {s["label"] for s in sources}
        genus_rows = [taxa.get(s["taxid"], {}) for s in sources]
        row.update(label=next(iter(labels)) if len(labels) == 1 else "conflict",
            categories=sorted({s["category"] for s in sources}), evidence_tiers=sorted({s["evidence_tier"] for s in sources}),
            genus_taxids=sorted({g["genus_taxid"] for g in genus_rows if g.get("genus_taxid")}),
            genus_names=sorted({g["genus_name"] for g in genus_rows if g.get("genus_name")}))
        problems = []
        if row["label"] == "conflict":
            problems.append("positive_negative_source_conflict")
        if any(not g.get("genus_taxid") for g in genus_rows):
            problems.append("missing_genus_taxonomy")
        for s in sources:
            if s["source"].startswith("HydDB"):
                meta = ncbi.get(s["accession"], {})
                if not meta:
                    problems.append("accession_summary_missing")
                elif int(meta["slen"]) != len(row["sequence"]):
                    problems.append("source_vs_current_accession_length_discrepancy")
        row["source_review_flags"] = sorted(set(problems))
    SeqIO.write([SeqRecord(Seq(r["sequence"]), id=r["unit_id"], description="") for r in rows], out/"identity_units.faa", "fasta")
    SeqIO.write([SeqRecord(Seq(r["sequence"].replace("U", "C")), id=r["unit_id"], description="") for r in rows], out/"leakage_units.faa", "fasta")
    # Include prior-query exposure as well as training seeds, so a past diagnostic cannot reappear as a blind test.
    exposure, exposed_taxids, inputs = {}, set(), []
    def expose(seq, source, role):
        key = digest(seq)
        item = exposure.setdefault(key, dict(sequence=seq, sources=[], roles=set()))
        item["sources"].append(source)
        item["roles"].add(role)
    for filename in ["references.faa", "reference_source_proteins.faa"]:
        path = ROOT/"references/processed"/filename
        inputs.append(path)
        for rec in SeqIO.parse(path, "fasta"):
            expose(str(rec.seq).upper(), filename+":"+rec.id, "training_reference")
    refs = json.loads((ROOT/"references/processed/references.json").read_text())
    genome_paths = [ROOT/"references/raw"/(s+".gb") for s in sorted({r["sample"] for r in refs})]
    genome_paths += sorted((ROOT/"references/raw/independent").glob("*.gb"))
    for path in genome_paths:
        inputs.append(path)
        for rec in SeqIO.parse(path, "genbank"):
            for feature in rec.features:
                if feature.type == "source":
                    exposed_taxids.update(v.split(":")[1] for v in feature.qualifiers.get("db_xref", []) if v.startswith("taxon:"))
                elif path.parent.name == "independent" and feature.type == "CDS":
                    for seq in feature.qualifiers.get("translation", []):
                        expose(seq.upper(), path.name+":"+feature.qualifiers.get("locus_tag", [""])[0], "previous_development")
    for path in sorted((ROOT/"references/raw/structures").glob("*.cif")):
        inputs.append(path)
        cif = MMCIF2Dict(str(path))
        for seq in cif.get("_entity_poly.pdbx_seq_one_letter_code_can", []):
            sequence = re.sub(r"\s+", "", seq).upper()
            if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+", sequence):
                expose(sequence, path.name, "previous_structure_or_template")
            else:
                raise ValueError(f"Unparsed canonical polymer in {path}")
        for field in ["_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id", "_entity_src_nat.pdbx_ncbi_taxonomy_id"]:
            exposed_taxids.update(v for v in cif.get(field, []) if v.isdigit())
    for path in sorted((ROOT/"results/hadar_pilot_round08/current").glob("*/proteins.faa")):
        inputs.append(path)
        for rec in SeqIO.parse(path, "fasta"):
            expose(str(rec.seq).upper(), path.parent.name+":"+rec.id, "previous_pilot")
    exposed_genera = {taxa[t]["genus_taxid"] for t in exposed_taxids if t in taxa and taxa[t]["genus_taxid"]}
    if not exposed_genera or len(list((ROOT/"results/hadar_pilot_round08/current").glob("*/proteins.faa"))) != 14:
        raise ValueError("Expected prior taxon and all 14-genome exposure controls")
    exposures = []
    for key, item in sorted(exposure.items()):
        item.update(exposure_id="prior_"+key[:20], sequence_sha256=key, roles=sorted(item["roles"]))
        exposures.append(item)
    SeqIO.write([SeqRecord(Seq(r["sequence"].replace("U", "C")), id=r["exposure_id"], description="") for r in exposures], out/"prior_exposure.faa", "fasta")
    (out/"prior_exposure.json").write_text(json.dumps(exposures, indent=2)+"\n")
    commands = []
    def run(cmd):
        commands.append(cmd)
        (out/"commands.json").write_text(json.dumps(commands, indent=2)+"\n")
        with (out/"blast.log").open("a") as handle:
            subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
    for name in ["leakage_units", "prior_exposure"]:
        run(["makeblastdb", "-in", str(out/(name+".faa")), "-dbtype", "prot", "-out", str(out/name)])
    fmt = "6 qseqid sseqid pident qlen slen qstart qend sstart send length bitscore evalue"
    for label, db, count in [("all_vs_all", "leakage_units", len(rows)), ("vs_prior", "prior_exposure", len(exposures))]:
        print("leakage_alignment", label, len(rows), count, flush=True)
        run(["blastp", "-query", str(out/"leakage_units.faa"), "-db", str(out/db), "-evalue", str(policy["blast_evalue"]),
             "-max_target_seqs", str(count), "-max_hsps", str(policy["blast_max_hsps"]), "-num_threads", "4", "-outfmt", fmt,
             "-out", str(out/(label+".tsv"))])
    parent = {r["unit_id"]: r["unit_id"] for r in rows}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[max(a, b)] = min(a, b)
    def alignment_rows(path):
        for line in path.read_text().splitlines():
            q, s, identity, qlen, slen, qs, qe, ss, se, length, bits, ev = line.split("\t")
            qcov, scov = (int(qe)-int(qs)+1)/int(qlen), (int(se)-int(ss)+1)/int(slen)
            yield q, s, float(identity), qcov, scov, qcov if int(qlen) <= int(slen) else scov
    edges = 0
    for q, s, identity, qcov, scov, shorter_cov in alignment_rows(out/"all_vs_all.tsv"):
        if q != s and identity >= policy["homolog_identity_percent"] and shorter_cov >= policy["homolog_shorter_sequence_coverage"]:
            union(q, s)
            edges += 1
    # Exact accession identity is also a biological link even if distinct source sequence versions disagree.
    alias_owner = {}
    for row in rows:
        for s in row["sources"]:
            key = s["source"].split("_")[0]+":"+s["accession"]
            if key in alias_owner:
                union(row["unit_id"], alias_owner[key])
            alias_owner[key] = row["unit_id"]
    prior = {r["exposure_id"]: r for r in exposures}
    overlap, comparable_max = set(), {}
    for q, s, identity, qcov, scov, shorter_cov in alignment_rows(out/"vs_prior.tsv"):
        if identity >= policy["homolog_identity_percent"] and shorter_cov >= policy["homolog_shorter_sequence_coverage"]:
            overlap.add(q)
        if "training_reference" in prior[s]["roles"] and min(qcov, scov) >= policy["remote_group4_both_sequence_coverage"]:
            comparable_max[q] = max(identity, comparable_max.get(q, 0))
    groups = collections.defaultdict(list)
    for row in rows:
        groups[find(row["unit_id"])].append(row)
    clusters = []
    for cid, members in sorted(groups.items()):
        genera = {g for row in members for g in row["genus_taxids"]}
        exposed = bool(genera & exposed_genera or {r["unit_id"] for r in members} & overlap)
        split = component_split(genera, exposed, all(not r["source_review_flags"] for r in members), policy)
        for row in members:
            row.update(homolog_component=cid, split=split, prior_homolog_overlap=row["unit_id"] in overlap,
                max_comparable_training_identity=comparable_max.get(row["unit_id"]),
                remote_group4=bool(any(v.startswith("[NiFe] Group 4") for v in row["categories"]) and
                    row["unit_id"] in comparable_max and comparable_max[row["unit_id"]] <= policy["remote_group4_max_training_identity_percent"]))
        clusters.append(dict(component=cid, split=split, units=len(members), genus_taxids=sorted(genera),
            exposed_genus_taxids=sorted(genera & exposed_genera), prior_homolog_units=sum(r["unit_id"] in overlap for r in members)))
    byid = {r["unit_id"]: r for r in rows}
    for variant in variants:
        unit = byid.get(variant["original_unit_id"])
        variant["inherited_split"] = unit["split"] if unit else next(r["split"] for r in rows if any(s["source"].startswith("HydDB") and s["accession"] == variant["accession"] for s in r["sources"]))
    # Invariant: accepted partitions share neither genus ownership nor homolog components.
    accepted = {name: [r for r in rows if r["split"] == name] for name in policy["genus_hash_percentiles"]}
    for a in accepted:
        for b in accepted:
            if a != b:
                assert not {g for r in accepted[a] for g in r["genus_taxids"]} & {g for r in accepted[b] for g in r["genus_taxids"]}
                assert not {r["homolog_component"] for r in accepted[a]} & {r["homolog_component"] for r in accepted[b]}
    for split, selected in accepted.items():
        SeqIO.write([SeqRecord(Seq(r["sequence"]), id=r["unit_id"], description="") for r in selected], out/(split+".faa"), "fasta")
    (out/"units.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2)+"\n")
    (out/"components.json").write_text(json.dumps(clusters, indent=2)+"\n")
    write_tsv(out/"source_selection.tsv", audit)
    write_tsv(out/"sequence_variant_audit.tsv", variants)
    write_tsv(out/"partition_summary.tsv", [dict(split=s, label=label, units=sum(r["label"] == label for r in selected),
        components=len({r["homolog_component"] for r in selected if r["label"] == label}),
        genera=len({g for r in selected if r["label"] == label for g in r["genus_taxids"]}))
        for s in sorted({r["split"] for r in rows}) for selected in [[r for r in rows if r["split"] == s]] for label in ["positive", "negative", "conflict"]])
    source_files = sorted(p for p in RAW.rglob("*") if p.is_file())
    frozen = dict(status="frozen_before_first_prediction", policy_sha256=sha(policy_file), builder_sha256=sha(Path(__file__)),
        sources={str(p.relative_to(ROOT)): sha(p) for p in [*source_files, *inputs]},
        prior_exposed_genus_taxids=sorted(exposed_genera), sequences=len(rows), prior_exposure_sequences=len(exposures),
        directed_homolog_edges=edges, partitions=dict(collections.Counter(r["split"] for r in rows)),
        prediction_models={str(p.relative_to(ROOT)): sha(p) for p in [ROOT/"references/processed/catalytic_families.hmm", ROOT/"references/processed/profile_manifest.json", ROOT/"scripts/annotate.py"]},
        artifacts={p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.suffix in [".json", ".faa", ".tsv"]},
        limitations=["Author-curated identity labels are not new biochemical measurements.", "No complete-complex or reaction truth is created by this corpus.",
            "Operational local-HSP clustering does not exclude every distant domain relationship.", "Reserved final validation has not been predicted."])
    (out/"freeze.json").write_text(json.dumps(frozen, indent=2)+"\n")
    print(json.dumps({k: frozen[k] for k in ["status", "sequences", "prior_exposure_sequences", "partitions"]}), flush=True)


if __name__ == "__main__":
    main()
