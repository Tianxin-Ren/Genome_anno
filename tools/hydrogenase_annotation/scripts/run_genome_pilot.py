#!/usr/bin/env python3
"""Run real genome FASTAs separately, preserve inputs, and compare a pinned baseline."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, write_tsv
from run_acceptance_checks import state_hashes


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle,delimiter="\t"))


def inventory(directory):
    rows=[]
    for path in sorted(directory.glob("*.fna")):
        records=list(SeqIO.parse(path,"fasta"))
        if not records or len({r.id for r in records})!=len(records):
            raise ValueError(f"Empty FASTA or repeated contig IDs: {path}")
        dna="".join(str(r.seq).upper() for r in records)
        if set(dna)-set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Non-IUPAC DNA characters: {path}")
        rows.append(dict(sample=path.stem,path=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            contig_count=len(records),bp=len(dna),ambiguous_bp=sum(c not in "ACGT" for c in dna),
            gc_fraction=(dna.count("G")+dna.count("C"))/max(1,sum(dna.count(c) for c in "ACGT")),
            headers=[r.description for r in records],
            source_qc_flags=["source_header_duplicated_possibly"] if any("duplicated-possibly" in r.description for r in records) else [],
            taxonomy="not_supplied",completeness="not_independently_verified"))
    if not rows:
        raise ValueError(f"No .fna inputs found: {directory}")
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--baseline-script",type=Path)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--threads",type=int,default=1)
    args=parser.parse_args()
    if args.workers<1 or args.threads<1:
        parser.error("workers and threads must be positive")
    inputs=inventory(args.inputs)
    output=args.out.resolve()
    if (output/"execution.json").exists():
        parser.error("Use a new pilot directory; existing execution evidence is preserved")
    output.mkdir(parents=True,exist_ok=True)
    (output/"input_manifest.json").write_text(json.dumps(inputs,ensure_ascii=False,indent=2)+"\n")
    before=state_hashes()
    scripts={"current":ROOT/"scripts/annotate.py"}
    if args.baseline_script:
        scripts["baseline"]=args.baseline_script.resolve()
    script_hashes={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in scripts.items()}
    env=dict(os.environ);env["PATH"]=str(ROOT/".env/bin")+os.pathsep+env.get("PATH","")
    jobs=[]
    report=dict(status="running",input_genomes=len(inputs),total_bp=sum(x["bp"] for x in inputs),
        before_hashes=before,script_hashes=script_hashes,workers=args.workers,threads_per_run=args.threads,executions=jobs,
        scope="User-supplied real genomes, not an experimental truth set or a representative 1000-MAG sample")
    def save():
        (output/"execution.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    save()

    def run(sample):
        executions=[]
        for version,script in scripts.items():
            directory=output/version/sample["sample"];directory.mkdir(parents=True,exist_ok=True)
            cmd=[sys.executable,str(script),"--genome",sample["path"],"--out",str(directory),"--threads",str(args.threads),"--translation-table","11"]
            if version=="current":
                cmd.extend(["--mode","genome"])
            start=time.monotonic()
            with (directory/"command.log").open("w") as log:
                process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
                _,status,usage=os.wait4(process.pid,0)
                process.returncode=os.waitstatus_to_exitcode(status)
            row=dict(sample=sample["sample"],version=version,command=cmd,exit_code=process.returncode,
                elapsed_seconds=round(time.monotonic()-start,3),wait4_cpu_user_seconds=usage.ru_utime,
                wait4_cpu_system_seconds=usage.ru_stime,
                wait4_peak_rss_bytes=usage.ru_maxrss if sys.platform=="darwin" else usage.ru_maxrss*1024,
                resource_scope="POSIX wait4 child rusage; not simultaneous process-tree peak memory",
                input_sha256=sample["sha256"],script_sha256=script_hashes[version])
            (directory/"pilot_execution.json").write_text(json.dumps(row,indent=2)+"\n")
            executions.append(row)
        return executions

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run,s):s for s in inputs}
        for future in as_completed(futures):
            jobs.extend(future.result());save()
            print(futures[future]["sample"],[(r["version"],r["exit_code"]) for r in jobs if r["sample"]==futures[future]["sample"]],flush=True)
    summary=[];differences=[]
    for sample in inputs:
        versions={}
        for version in scripts:
            d=output/version/sample["sample"]
            execution=next(j for j in jobs if j["sample"]==sample["sample"] and j["version"]==version)
            if execution["exit_code"]:
                continue
            rows=read_tsv(d/"complexes.tsv");versions[version]={r["protein_id"]:r for r in rows}
            manifest=json.loads((d/"run_manifest.json").read_text())
            summary.append(dict(sample=sample["sample"],version=version,bp=sample["bp"],proteins=manifest["protein_count"],
                candidates=len(rows),high_structure=sum(r["structure_confidence"]=="high" for r in rows),
                potential_complexes=sum(bool(r.get("potential_function_ids")) for r in rows),
                reaction_rows=len(read_tsv(d/"potential_functions.tsv")),unresolved=sum(r["structure_id"]=="unresolved" for r in rows),
                group4_hmm_anchors=sum("NiFe_Group4" in r.get("hmm_support","") for r in rows),
                elapsed_seconds=execution["elapsed_seconds"],wait4_peak_rss_bytes=execution["wait4_peak_rss_bytes"]))
        if "baseline" in versions and "current" in versions:
            for gid in sorted(set(versions["baseline"])|set(versions["current"])):
                old,new=versions["baseline"].get(gid,{}),versions["current"].get(gid,{})
                keys=["structure_id","structure_confidence","potential_function_ids","hmm_support"]
                if any(old.get(k)!=new.get(k) for k in keys):
                    differences.append(dict(sample=sample["sample"],protein_id=gid,
                        **{"baseline_"+k:old.get(k,"") for k in keys},**{"current_"+k:new.get(k,"") for k in keys},
                        current_flags=new.get("flags",""),current_missing=new.get("missing_roles",""),
                        review_status="requires_biological_review_not_accuracy_truth"))
    write_tsv(output/"genome_summary.tsv",summary)
    write_tsv(output/"baseline_differences.tsv",differences,["sample","protein_id",*[prefix+k for prefix in ["baseline_","current_"]
        for k in ["structure_id","structure_confidence","potential_function_ids","hmm_support"]],"current_flags","current_missing","review_status"])
    report["after_hashes"]=state_hashes()
    report["inputs_unchanged"]=all(hashlib.sha256(Path(s["path"]).read_bytes()).hexdigest()==s["sha256"] for s in inputs)
    report["status"]="complete" if all(j["exit_code"]==0 for j in jobs) and report["inputs_unchanged"] and before==report["after_hashes"] else "failed_or_state_changed"
    report["failed"]=[j for j in jobs if j["exit_code"]]
    report["summary"]=summary;report["difference_rows"]=len(differences)
    report["limitations"]=["Genetic code 11 explicitly used, source taxonomy not supplied.",
        "Circularity follows source headers, not independent closure validation; origin-spanning CDS/operons still require review.",
        "Header duplication warning retained for u1320898; no sample silently excluded.",
        "Candidate/function counts are predictions, not true enzyme counts or precision estimates."]
    save()
    print(json.dumps(dict(status=report["status"],genomes=len(inputs),runs=len(jobs),differences=len(differences))))
    if report["status"]!="complete":
        sys.exit(1)


if __name__=="__main__":
    main()
