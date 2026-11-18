#!/usr/bin/env python3
"""Build the reviewed seed library, family HMMs, and auditable Excel import."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import openpyxl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from reference_config import SYSTEMS, DECOYS, UNRESOLVED, SOURCE_REVIEW_NOTES
from annotate import ROOT, extract_genbank, write_tsv
from profile_models import build_profile, profile_training
from audit_maturation_domains import audit_review
from extract_structure_sites import main as extract_structural_sources
from fes_sites import build_model as build_fes_model


def pick(genes, selector):
    matches = [g for g in genes if selector in [g["id"], g["gene"], *g["aliases"].split(";")]]
    if len(matches) != 1:
        raise ValueError(f"Reference selector {selector!r}: expected one translated CDS, found {len(matches)}")
    return matches[0]


def ligand_positions(sequence, family):
    if family == "NiFe":
        pairs = list(re.finditer("[CU]..[CU]", sequence))
        assert len(pairs) >= 2 and pairs[-1].start() - pairs[0].start() > 150
        return [pairs[0].start()+1, pairs[0].start()+4, pairs[-1].start()+1, pairs[-1].start()+4]
    if family == "FeFe":
        match = re.search(r"C...C[IV][GNS]GGG", sequence)
        assert match is not None, "FeFe H-domain signature missing in curated reference"
        return [match.start()+1, match.start()+5]
    assert sequence[171] == "C", "M. marburgensis Hmd FeGP ligand: C172; numbering is species-specific"
    return [172]


def reviewed_maturation_entries(review, genomes, occupied):
    """An explicit, sequence-pinned review replaces gene-name based import."""
    seen = set()
    for row in review["records"]:
        key = (row["sample"], row["locus_tag"])
        if key in seen:
            raise ValueError(f"Duplicate maturation review entry: {key}")
        seen.add(key)
        gene = pick(genomes[row["sample"]], row["locus_tag"])
        if (gene["protein_accession"] != row["protein_accession"] or
                hashlib.sha256(gene["sequence"].encode()).hexdigest() != row["sequence_sha256"]):
            raise ValueError(f"Maturation review sequence/accession changed: {key}")
        if row["status"] == "exclude_non_maturase_homonym":
            continue
        if row["status"] != "reviewed_database_maturation_family":
            raise ValueError(f"Unreviewed maturation entry: {key}")
        if key in occupied or row["sequence_sha256"] in {r["sequence_sha256"] for r in occupied.values()}:
            raise ValueError(f"Maturation entry overlaps a functional subunit: {key}")
        if not row.get("source_url") or not row.get("review_note") or not row.get("source_product"):
            raise ValueError(f"Missing maturation review evidence: {key}")
        yield row


def reference_region(gene, review=None):
    sequence=gene["sequence"]
    metadata=dict(source_protein_length=len(sequence),source_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        protein_region_start=1,protein_region_end=len(sequence),coverage_scope="full_reference_protein")
    if review is not None:
        if gene["protein_accession"]!=review["source_protein_accession"] or metadata["source_sequence_sha256"]!=review["source_sequence_sha256"]:
            raise ValueError("Reviewed reference domain source changed")
        start,end=review["protein_region_start"],review["protein_region_end"]
        if not 1<=start<=end<=len(sequence):
            raise ValueError("Reviewed reference domain outside source protein")
        sequence=sequence[start-1:end]
        if hashlib.sha256(sequence.encode()).hexdigest()!=review["region_sequence_sha256"]:
            raise ValueError("Reviewed reference domain sequence changed")
        metadata.update(protein_region_start=start,protein_region_end=end,coverage_scope="reviewed_reference_domain",
            domain_evidence_url=review["evidence_url"],domain_experimental_accession=review["experimental_accession"])
    return sequence,metadata


def reviewed_function_records(records):
    """Correct the reviewed H01 scope while retaining the original workbook rows."""
    functions = {r["功能ID"]: dict(r) for r in records}
    hmd = functions["H01"]
    field = "大亚基特征或必须保留的分型证据"
    original = "辨认有活性的Hmd-I样序列与FeGP结合环境，排除HmdII类失活同源物"
    corrected = "辨认Hmd-I样分支与FeGP结合环境；经典H4MPT反应须与HmdII/III等近缘解释区分，不能由家族名称推定普遍失活"
    if hmd["结构系统ID"] != "Fe_Hmd" or hmd[field] not in {original, corrected}:
        raise ValueError("H01 source metadata differs from its reviewed scope correction")
    hmd[field] = corrected
    return functions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    args = parser.parse_args()
    output = ROOT / "references/processed"
    output.mkdir(parents=True, exist_ok=True)
    (output/"build_manifest.json").write_text(json.dumps(dict(status="running"))+"\n")
    workbook = openpyxl.load_workbook(args.workbook, data_only=True, read_only=True)
    imported = {}
    for sheet in workbook:
        rows = list(sheet.values)
        if not rows:
            continue
        headers = [str(x or f"column_{i+1}") for i,x in enumerate(rows[0])]
        records = [dict(zip(headers, r)) for r in rows[1:] if any(x is not None for x in r)]
        imported[sheet.title] = records
        write_tsv(output / f"excel_{sheet.title}.tsv", records, headers)
    (output / "excel_import.json").write_text(json.dumps(imported, ensure_ascii=False, indent=2)+"\n")
    functions = reviewed_function_records(imported["实证功能表"])
    samples = {s["sample"] for s in SYSTEMS} | set(DECOYS) | set(UNRESOLVED)
    genomes = {name:extract_genbank(ROOT / "references/raw" / f"{name}.gb")[0] for name in sorted(samples)}
    source_records = {name:SeqIO.read(ROOT/"references/raw"/f"{name}.gb", "genbank") for name in samples}
    source_features = {name:{f.qualifiers.get("locus_tag", f.qualifiers.get("protein_id", [""]))[0]:f
        for f in record.features if f.type == "CDS"} for name,record in source_records.items()}
    maturation_review = json.loads((ROOT/"config/maturation_review.json").read_text())
    domain_audit=audit_review(maturation_review,genomes,ROOT,output/"maturation_domain_audit")
    if domain_audit["status"]!="passed":
        raise ValueError("Maturation domain audit failed; inspect maturation_domain_audit/summary.json")
    audited={(r["sample"],r["locus_tag"]):r for r in domain_audit["rows"]}
    references, fasta = [], []
    domain_review=json.loads((ROOT/"config/reference_domains.json").read_text())["domains"]
    domains={r["reference_id"]:r for r in domain_review}
    if len(domains)!=len(domain_review):
        raise ValueError("Duplicate reviewed reference domain")
    full_sources=[]
    operons=ROOT/"references/operons"
    operons.mkdir(exist_ok=True)
    operon_manifest=[]

    def add(sample, selector, sid, role, kind, family=""):
        gene = pick(genomes[sample], selector)
        rid = f"{sid}__{role}"
        review=domains.get(rid)
        if review and (review["sample"]!=sample or review["locus_tag"]!=gene["id"]):
            raise ValueError("Reviewed reference domain CDS mapping changed")
        sequence,region=reference_region(gene,review)
        full_sources.append(SeqRecord(Seq(gene["sequence"]),id=rid,description="complete_source_CDS_translation"))
        entry = dict(reference_id=rid, structure_id=sid, role=role, kind=kind, center=family,
                     sample=sample, locus_tag=gene["id"], aliases=gene["aliases"],
                     protein_accession=gene["protein_accession"], contig=gene["contig"],
                     start=gene["start"], end=gene["end"], strand=gene["strand"],
                     sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                     cysteine_positions=[i+1 for i,a in enumerate(sequence) if a=="C"],
                     selenocysteine_positions=[i+1 for i,a in enumerate(sequence) if a=="U"],
                     diagnostic_positions=ligand_positions(sequence,family) if kind=="catalytic" else [])
        entry.update(region)
        source = source_records[sample]
        feature = source_features[sample][gene["id"]]
        system = next((s for s in SYSTEMS if s["structure_id"] == sid), None)
        entry.update(source_accession=source.id,
                     source_url="https://www.ncbi.nlm.nih.gov/nuccore/"+source.id,
                     protein_source_url="https://www.ncbi.nlm.nih.gov/protein/"+gene["protein_accession"],
                     source_product=feature.qualifiers.get("product", [""])[0],
                     source_inference=feature.qualifiers.get("inference", []),
                     evidence_url=system["evidence_url"] if system else "https://www.ncbi.nlm.nih.gov/nuccore/"+source.id,
                     evidence_scope="experimental_reference_system_CDS_mapping" if system and kind != "context" else "database_annotation_not_direct_biochemical_truth",
                     domain_validation="pending_independent_profile_check")
        entry["diagnostic_residues"]={str(i):sequence[i-1] for i in entry["diagnostic_positions"]}
        references.append(entry)
        fasta.append(SeqRecord(Seq(sequence), id=rid, description=""))
        return rid

    for system in SYSTEMS:
        source=source_records[system["sample"]]
        local=[pick(genomes[system["sample"]],selector) for selector in system["local"].values()]
        lo,hi=min(g["start"] for g in local),max(g["end"] for g in local)
        region_start,region_end=max(1,lo-500),min(len(source),hi+500)
        region=source[region_start-1:region_end]
        region.id=source.id
        transcript=source.annotations.get("molecule_type")=="mRNA"
        region_kind="experimental cDNA CDS, not an operon" if transcript else "local reference gene cluster"
        region.description=f"{system['structure_id']}; {source.id}:{region_start}..{region_end}; {region_kind}"
        SeqIO.write(region,operons/f"{system['structure_id']}.local.gb","genbank")
        SeqIO.write(region,operons/f"{system['structure_id']}.local.fna","fasta")
        SeqIO.write([SeqRecord(Seq(g["sequence"]),id=g["id"],description="") for g in local],
                    operons/f"{system['structure_id']}.local.faa","fasta")
        operon_manifest.append(dict(structure_id=system["structure_id"],sample=system["sample"],accession=source.id,
            organism=source.annotations.get("organism",""),core_start_1based=lo,core_end_1based=hi,
            extracted_region_start_1based=region_start,extracted_region_end_1based=region_end,
            strand=local[0]["strand"],local_roles=";".join(system["local"]),
            remote_roles=";".join(system["remote"]),context_roles=";".join(system.get("context",{})),
            reference_record_type="cDNA_single_CDS_not_operon" if transcript else "genomic_gene_cluster",
            evidence_url=system["evidence_url"],
            evidence_scope="experimental enzyme system; CDS mapping; not a new transcription experiment"))
        for role,selector in {**system["local"], **system["remote"], **system.get("context",{})}.items():
            kind = "catalytic" if role==system["catalytic_role"] else "context" if role in system.get("context",{}) else "module"
            add(system["sample"], selector, system["structure_id"], role, kind, system["center"])
        for fid in system["function_ids"]:
            assert functions[fid]["结构系统ID"] == system["structure_id"], fid
        system["reference_function_records"] = [functions[fid] for fid in system["function_ids"]]
    for sample,selectors in DECOYS.items():
        for selector in selectors:
            add(sample,selector,"decoy",sample+"_"+selector,"decoy")
    for sample,selectors in UNRESOLVED.items():
        for selector,family in selectors.items():
            add(sample,selector,"unresolved",sample+"_"+selector,"unresolved",family)
    occupied = {(r["sample"], r["locus_tag"]):r for r in references if r["kind"] in ["catalytic", "module"]}
    for row in reviewed_maturation_entries(maturation_review, genomes, occupied):
        add(row["sample"], row["locus_tag"], "maturation", row["sample"]+"_"+row["role"], "maturation")
        evidence=audited[(row["sample"],row["locus_tag"])]
        references[-1].update(evidence_scope="database_CDS_review_plus_independent_profile_architecture", review_note=row["review_note"],
            evidence_url=row["source_url"],domain_validation=evidence["domain_validation"],
            domain_specificity_scope=evidence["specificity_scope"],domain_evidence=evidence["domain_evidence"])
    assert len({r["reference_id"] for r in references}) == len(references)
    if domains.keys()-{r["reference_id"] for r in references}:
        raise ValueError("Reviewed reference domain not used in reference build")
    write_tsv(operons/"operon_manifest.tsv",operon_manifest)
    SeqIO.write(fasta, output / "references.faa", "fasta")
    SeqIO.write(full_sources,output/"reference_source_proteins.faa","fasta")
    (output / "references.json").write_text(json.dumps(references, indent=2)+"\n")
    (ROOT / "config/systems.json").write_text(json.dumps(SYSTEMS, ensure_ascii=False, indent=2)+"\n")
    from audit_reference_domains import audit as audit_reference_roles, write_reference_tables
    references,reference_audit=audit_reference_roles(references,ROOT,output/'reference_domain_audit')
    if reference_audit['status']!='review_complete':raise ValueError('Reference role audit incomplete')
    (output/'references.json').write_text(json.dumps(references,indent=2)+'\n')
    reviewed = {fid for s in SYSTEMS for fid in s["function_ids"]}
    write_tsv(output / "function_rule_status.tsv", [dict(function_id=fid, structure_id=f["结构系统ID"],
        status="reviewed_seed_rule" if fid in reviewed else "imported_not_enabled",
        excel_source=f["关键机制或来源URL"]) for fid,f in functions.items()])
    write_tsv(output/"source_review.tsv", [dict(function_id=fid,structure_id=s["structure_id"],
        imported_primary_url=functions[fid]["关键机制或来源URL"],active_primary_url=s["evidence_url"],
        review_note=SOURCE_REVIEW_NOTES.get(fid,"Reference enzyme system reviewed; genome CDS mapping does not establish function in a new query"))
        for s in SYSTEMS for fid in s["function_ids"]])
    write_tsv(output/"maturation_review.tsv", maturation_review["records"])
    write_reference_tables(references,output)
    subprocess.run(["makeblastdb", "-in",str(output/"references.faa"),"-dbtype","prot","-out",str(output/"seed")], check=True)
    profiles = []
    profile_manifest=[]
    for family in ["NiFe","FeFe","Fe"]:
        ids = {r["reference_id"] for r in references if r["center"]==family and r["kind"] in ["catalytic","unresolved"]}
        unaligned, aligned, profile = [output/f"{family}.{suffix}" for suffix in ["faa","aligned.faa","hmm"]]
        SeqIO.write([SeqRecord(Seq(str(s.seq).replace("U","X")),id=s.id,description="") for s in fasta if s.id in ids], unaligned, "fasta")
        if len(ids)>1:
            with aligned.open("w") as handle:
                subprocess.run(["mafft","--quiet","--auto",str(unaligned)],stdout=handle,check=True)
        else:
            shutil.copyfile(unaligned,aligned)
        with (ROOT / "logs" / f"hmmbuild_{family}.log").open("w") as handle:
            subprocess.run(["hmmbuild","--amino","-n",family,str(profile),str(aligned)],stdout=handle,check=True)
        profiles.append(profile.read_text())
        profile_manifest.append(dict(profile=family,family=family,sequences=sorted(ids),
            model_sha256=hashlib.sha256(profile.read_bytes()).hexdigest(),scope="legacy_seed_family_model"))
    branch_config=json.loads((ROOT/"config/discovery_profiles.json").read_text())
    records={s.id:s for s in fasta}
    for config in branch_config["profiles"]:
        training,provenance=profile_training(ROOT,config,records)
        profile,detail=build_profile(training,output,config["name"],config["alignment_occupancy"],config.get('reference_core'))
        profiles.append(profile.read_text())
        profile_manifest.append(dict(**detail,**provenance,family=config["family"],scope=config["scope"],evidence=config["evidence"]))
    (output/"profile_manifest.json").write_text(json.dumps(profile_manifest,indent=2)+"\n")
    (output/"catalytic_families.hmm").write_text("".join(profiles))
    extract_structural_sources()
    fes_manifest=build_fes_model()
    from active_centers import build_model as build_active_center_model
    active_manifest=build_active_center_model()
    from split_centers import build_model as build_split_center_model
    build_split_center_model()
    from catalytic_confirmation import build_model as build_confirmation_model
    build_confirmation_model()
    provenance = dict(status="complete",workbook_name=args.workbook.name, workbook_sha256=hashlib.sha256(args.workbook.read_bytes()).hexdigest(),
                      function_records=len(functions), reviewed_structures=len(SYSTEMS),
                      enabled_function_records=sum(len(s["function_ids"]) for s in SYSTEMS),
                      references=len(references), profile_scope="small curated seed library; not HydDB/Pfam classifier",
                      reference_domain_audit_sha256=hashlib.sha256((output/'reference_domain_audit/summary.json').read_bytes()).hexdigest(),
                      maturation_review_sha256=hashlib.sha256((ROOT/"config/maturation_review.json").read_bytes()).hexdigest(),
                      reference_domains_sha256=hashlib.sha256((ROOT/"config/reference_domains.json").read_bytes()).hexdigest(),
                      reference_source_proteins_sha256=hashlib.sha256((output/"reference_source_proteins.faa").read_bytes()).hexdigest(),
                      maturation_domain_policy_sha256=domain_audit["policy_sha256"],
                      maturation_domain_audit_sha256=hashlib.sha256((output/"maturation_domain_audit/summary.json").read_bytes()).hexdigest(),
                      excluded_maturation_homonyms=[r["locus_tag"] for r in maturation_review["records"] if r["status"].startswith("exclude")],
                      branch_profile_config_sha256=hashlib.sha256((ROOT/"config/discovery_profiles.json").read_bytes()).hexdigest(),
                      profile_manifest_sha256=hashlib.sha256((output/"profile_manifest.json").read_bytes()).hexdigest(),
                      fes_model_manifest_sha256=hashlib.sha256((output/"fes/manifest.json").read_bytes()).hexdigest(),
                      fes_site_policy_sha256=fes_manifest['policy_sha256'],
                      active_center_model_manifest_sha256=hashlib.sha256((output/'active_centers/manifest.json').read_bytes()).hexdigest(),
                      active_center_policy_sha256=active_manifest['policy_sha256'],
                      confirmation_model_manifest_sha256=hashlib.sha256((output/'catalytic_confirmation/manifest.json').read_bytes()).hexdigest(),
                      split_center_model_manifest_sha256=hashlib.sha256((output/'split_centers/manifest.json').read_bytes()).hexdigest(),
                      combined_hmm_sha256=hashlib.sha256((output/"catalytic_families.hmm").read_bytes()).hexdigest(),
                      alignment_mask="Sec U -> X only in MSA/HMM/trees; sequence database and residue evidence preserve U")
    (output/"build_manifest.json").write_text(json.dumps(provenance,indent=2)+"\n")
    print(json.dumps(provenance,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
