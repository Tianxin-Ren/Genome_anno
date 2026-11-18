#!/usr/bin/env python3
"""Check accession-pinned maturation candidates against independent curated profiles."""
import argparse
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, extract_genbank, write_tsv


def architecture_supported(accessions,rule):
    return set(rule.get("all_of",[]))<=set(accessions) and (not rule.get("any_of") or bool(set(rule["any_of"])&set(accessions)))


def scan_profile_domains(profile,query,models,min_coverage,output,threads=2):
    output.mkdir(parents=True,exist_ok=True)
    dom=output/"domains.domtbl"
    target=output/"domains.tblout"
    with (output/"hmmsearch.log").open("w") as log:
        subprocess.run(["hmmsearch","--cpu",str(threads),"--cut_ga","--noali","--domtblout",str(dom),"--tblout",str(target),str(profile),str(query)],stdout=log,check=True)
    domains=defaultdict(list)
    for line in dom.read_text().splitlines():
        if line.startswith("#"):
            continue
        c=line.split();model=models[c[3]]
        coverage=(int(c[16])-int(c[15])+1)/int(c[5])
        passed=float(c[7])>=model["sequence_ga"] and float(c[13])>=model["domain_ga"] and coverage>=min_coverage
        domains[c[0]].append(dict(accession=model["accession"],profile=c[3],domain_bitscore=float(c[13]),
            sequence_bitscore=float(c[7]),independent_evalue=float(c[12]),model_coverage=coverage,
            model_start=int(c[15]),model_end=int(c[16]),model_length=int(c[5]),query_length=int(c[2]),
            query_start=int(c[17]),query_end=int(c[18]),passes=passed,
            sequence_ga=model["sequence_ga"],domain_ga=model["domain_ga"],source=model["url"],model_sha256=model["sha256"]))
    return domains


def audit_review(review,genomes,root,output):
    output.mkdir(parents=True,exist_ok=True)
    policy_path=root/"config/maturation_domain_policy.json"
    policy=json.loads(policy_path.read_text())
    models={}; texts=[]
    for model in policy["models"]:
        path=root/"references/raw/domain_review"/model["file"]
        data=path.read_bytes()
        if hashlib.sha256(data).hexdigest()!=model["sha256"]:
            raise ValueError(f"Domain profile hash changed: {path}")
        text=data.decode();ga=next((line.split()[1:] for line in text.splitlines() if line.startswith("GA ")),None)
        if ga is None:
            raise ValueError(f"Domain profile has no curated gathering threshold: {path}")
        models[model["profile"]]={**model,"sequence_ga":float(ga[0].strip(';')),"domain_ga":float(ga[1].strip(';'))}
        texts.append(text)
    combined=output/"maturation_domains.hmm";combined.write_text("".join(texts))
    records=[];by_key={}
    for row in review["records"]:
        gene=next(g for g in genomes[row["sample"]] if g["id"]==row["locus_tag"])
        if hashlib.sha256(gene["sequence"].encode()).hexdigest()!=row["sequence_sha256"]:
            raise ValueError(f"Maturation audit input changed: {row['locus_tag']}")
        key=row["sample"]+"__"+row["locus_tag"]
        by_key[key]=row
        records.append(SeqRecord(Seq(gene["sequence"].replace('U','X')),id=key,description=""))
    query=output/"reviewed_candidates.faa";SeqIO.write(records,query,"fasta")
    domains=scan_profile_domains(combined,query,models,policy["min_model_coverage"],output)
    audit=[]
    for key,row in by_key.items():
        accessions={d["accession"].split('.')[0] for d in domains[key] if d["passes"]}
        expected=architecture_supported(accessions,policy["rules"][row["role"]])
        excluded=row["status"]=="exclude_non_maturase_homonym"
        alternative=bool({"PF00970","PF00175"}&accessions)
        passed=(not expected and alternative) if excluded else expected
        audit.append(dict(sample=row["sample"],locus_tag=row["locus_tag"],role=row["role"],
            protein_accession=row["protein_accession"],sequence_sha256=row["sequence_sha256"],
            review_status=row["status"],passed=passed,
            domain_validation="excluded_homonym_alternative_domain_confirmed" if excluded and passed else
                "profile_architecture_supported" if passed else "requires_manual_review",
            accepted_accessions=";".join(sorted(accessions)),specificity_scope=policy["rules"][row["role"]]["scope"],
            cofactor_assembly="unknown",domain_evidence=domains[key]))
    write_tsv(output/"audit.tsv",[{**r,"domain_evidence":json.dumps(r["domain_evidence"])} for r in audit])
    summary=dict(status="passed" if all(r["passed"] for r in audit) else "failed",total=len(audit),passed=sum(r["passed"] for r in audit),
        failed=[{k:v for k,v in r.items() if k!="domain_evidence"} for r in audit if not r["passed"]],
        rows=audit,profile_models=models,min_model_coverage=policy["min_model_coverage"],
        combined_model_sha256=hashlib.sha256(combined.read_bytes()).hexdigest(),policy_sha256=hashlib.sha256(policy_path.read_bytes()).hexdigest(),
        review_sha256=hashlib.sha256(json.dumps(review,sort_keys=True).encode()).hexdigest(),
        caveat=policy["caveat"])
    (output/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=ROOT/"results/maturation_domain_audit")
    args=parser.parse_args()
    review=json.loads((ROOT/"config/maturation_review.json").read_text())
    genomes={s:extract_genbank(ROOT/"references/raw"/f"{s}.gb")[0] for s in {r["sample"] for r in review["records"]}}
    result=audit_review(review,genomes,ROOT,args.out)
    print(json.dumps({k:v for k,v in result.items() if k!="rows"},ensure_ascii=False,indent=2))
    if result["status"]!="passed":
        raise SystemExit(1)


if __name__=="__main__":
    main()
