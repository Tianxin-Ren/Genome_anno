#!/usr/bin/env python3
"""Verify branch-model coverage and negative competition, separately from blind acceptance."""
import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, annotate, extract_genbank, write_tsv
from evaluate_profile_design import search
from discovery_evidence import search_discovery


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=ROOT/"results/round02_discovery_checks")
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    refs=json.loads((ROOT/"references/processed/references.json").read_text())
    systems=json.loads((ROOT/"config/systems.json").read_text())
    profiles=json.loads((ROOT/"references/processed/profile_manifest.json").read_text())
    fasta={s.id:s for s in SeqIO.parse(ROOT/"references/processed/references.faa","fasta")}
    branch=next(p for p in profiles if p["profile"]=="NiFe_Group4")
    target_ids={r["reference_id"] for r in refs if r["kind"]=="catalytic" and r["structure_id"] in
        {s["structure_id"] for s in systems if s["reference_clade"].startswith("NiFe 4")}}
    query=args.out/"reference_and_negative_controls.faa"
    ids=target_ids | {r["reference_id"] for r in refs if r["kind"]=="decoy"}
    SeqIO.write([fasta[rid] for rid in sorted(ids)],query,"fasta")
    domains=search(ROOT/"references/processed/catalytic_families.hmm",query,args.out,"production_reference_scan")
    supported={d["protein_id"] for d in domains if d["profile"]=="NiFe_Group4" and d["passes"]}
    checks=[]
    def check(name,passed,observed):
        checks.append(dict(test=name,passed=bool(passed),observed=observed))
    for rid in sorted(target_ids):
        check("HMM_only_seed:"+rid,rid in supported,"No BLAST evidence used")
    genome=ROOT/"references/raw/independent/tkivui_DSM2030.gb"
    genes=extract_genbank(genome)[0]
    catalytic=next(g for g in genes if "TKV_c19740" in [g["id"],*g["aliases"].split(";")])
    digest=hashlib.sha256(catalytic["sequence"].encode()).hexdigest()
    check("development_case_excluded_from_training",digest not in branch["training_sequence_sha256"].values(),digest)
    found=annotate(copy.deepcopy(genes),args.out/"development_genome",input_files=[genome],mode="genome")
    candidates=[]
    for line in (args.out/"development_genome/family.domtbl").read_text().splitlines():
        if line.startswith("#"):
            continue
        c=line.split()
        if c[0]==catalytic["id"] and c[3]=="NiFe_Group4":
            candidates.append(dict(evalue=float(c[12]),coverage=(int(c[16])-int(c[15])+1)/int(c[5])))
    check("HMM_only_external_development_case",any(d["evalue"]<=1e-5 and d["coverage"]>=0.55 for d in candidates),candidates)
    anchor=next(a for a in found if a["protein_id"]==catalytic["id"])
    check("profile_name_not_misreported_as_family",anchor["catalytic_family"]=="NiFe",anchor["catalytic_family"])
    check("no_function_from_discovery_only",not anchor["potential_function_ids"],anchor["function_status"])
    negatives=[]
    for ref in [r for r in refs if r["kind"]=="decoy"]:
        source=extract_genbank(ROOT/"references/raw"/f"{ref['sample']}.gb")[0]
        g=copy.deepcopy(next(g for g in source if g["id"]==ref["locus_tag"]))
        g["id"]=ref["reference_id"]
        negatives.append(g)
    found=annotate(negatives,args.out/"negative_competition",mode="assembly")
    check("hard_negative_no_high_function",not any(a.get("potential_function_ids") for a in found),found)
    from identity_truth import adjudications
    reviewed=adjudications(ROOT)
    hmd_ids={r['reference_id'] for r in refs if r['kind']=='decoy' and reviewed.get(r['sequence_sha256'],{}).get('family')=='Fe'}
    check("known_decoys_rejected_despite_branch_HMM",{a['protein_id'] for a in found}==hmd_ids and
        all(a['catalytic_family']=='Fe' and a['catalytic_status']=='candidate' and a['structure_id']=='unresolved' and
            'Hmd_paralog_competitor_blocks_function_transfer_not_discovery' in a['flags'] for a in found),
        'Nuo/Mbs/other decoys rejected; only two reviewed Hmd paralogs retained as unresolved Fe candidates, never assigned Hmd function')
    corpus=json.loads((ROOT/'references/validation/round09_identity_corpus/units.json').read_text())
    unit=next(r for r in corpus if any(s['accession']=='WP_028309992.1' for s in r['sources']))
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    probe=SeqRecord(Seq(unit['sequence']),id=unit['unit_id'],description='')
    thresholds=json.loads((ROOT/'config/thresholds.json').read_text());batches=[]
    for label,records in [('single',[probe]),('context',[probe,*[fasta[r['reference_id']] for r in refs if r['kind']=='decoy']])]:
        directory=args.out/('fixed_search_space_'+label);directory.mkdir()
        query=directory/'query.faa';SeqIO.write(records,query,'fasta')
        hits,_,_=search_discovery(query,directory,ROOT,2,thresholds)
        batches.append([h for h in hits if h['protein_id']==probe.id])
    check('fixed_search_space_batch_invariance',batches[0]==batches[1],batches)
    check('collinear_sensory_development_candidate_recovered',any(h['profile']=='NiFe' and h['evidence_kind']=='collinear_domain_chain' for h in batches[0]),batches[0])
    write_tsv(args.out/"checks.tsv",checks)
    report=dict(status="passed" if all(c["passed"] for c in checks) else "failed",passed=sum(c["passed"] for c in checks),
        total=len(checks),failed=[c for c in checks if not c["passed"]],
        independent_Group4_acceptance="not_met: one development case; requires frozen independent sample count and clade coverage",
        source_withheld_diagnostics=str(ROOT/"results/profile_design_round02/summary.json"),
        model_sha256=hashlib.sha256((ROOT/"references/processed/catalytic_families.hmm").read_bytes()).hexdigest())
    (args.out/"summary.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if report["failed"]:
        sys.exit(1)


if __name__=="__main__":
    main()
