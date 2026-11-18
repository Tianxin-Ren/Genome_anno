#!/usr/bin/env python3
"""Compare seed model designs without promoting a development case to blind-test truth."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, extract_genbank, write_tsv
from profile_models import build_profile


def search(profile,query,directory,label,reference_sequences=None):
    dom=directory/(label+".domtbl")
    with (directory/(label+".log")).open("w") as handle:
        subprocess.run(["hmmsearch","--cpu","2",*(['-Z',str(reference_sequences)] if reference_sequences is not None else []),"--domtblout",str(dom),str(profile),str(query)],stdout=handle,check=True)
    domains=[]
    for line in dom.read_text().splitlines():
        if line.startswith("#"):
            continue
        c=line.split();coverage=(int(c[16])-int(c[15])+1)/int(c[5])
        domains.append(dict(protein_id=c[0],profile=c[3],independent_evalue=float(c[12]),domain_score=float(c[13]),
            model_coverage=coverage,query_start=int(c[17]),query_end=int(c[18]),
            passes=float(c[12])<=1e-5 and coverage>=0.55))
    return domains


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=ROOT/"results/profile_design_round02")
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    refs=json.loads((ROOT/"references/processed/references.json").read_text())
    systems=json.loads((ROOT/"config/systems.json").read_text())
    ids={s["structure_id"] for s in systems if s["reference_clade"].startswith("NiFe 4")}
    selected=[r for r in refs if r["kind"]=="catalytic" and r["structure_id"] in ids or
        r["kind"]=="unresolved" and r["sample"]=="ecoli_MG1655"]
    fasta={s.id:s for s in SeqIO.parse(ROOT/"references/processed/references.faa","fasta")}
    training=[fasta[r["reference_id"]] for r in selected]
    # This case is a development challenge, never a final test used to justify performance.
    genome=ROOT/"references/raw/independent/tkivui_DSM2030.gb"
    genes=extract_genbank(genome)[0]
    cat=next(g for g in genes if "TKV_c19740" in [g["id"],*g["aliases"].split(";")])
    challenge=args.out/"development_challenge.faa"
    SeqIO.write(SeqRecord(Seq(cat["sequence"]),id=cat["id"],description=""),challenge,"fasta")
    frozen=dict(split="development_challenge_not_independent_final_test",sample="T. kivui DSM2030 Ech2E",
        genome_version="CP009170.1",protein_accession=cat["protein_accession"],selector="TKV_c19740",protein_id=cat["id"],
        sequence_sha256=hashlib.sha256(cat["sequence"].encode()).hexdigest(),evidence="https://pmc.ncbi.nlm.nih.gov/articles/PMC9356269/",
        rule="Not added to model training or positive sequence reference. Changes chosen using this case invalidate it as final blind-test evidence.")
    (args.out/"challenge_manifest.json").write_text(json.dumps(frozen,indent=2)+"\n")
    models=[];rows=[]
    original=search(ROOT/"references/processed/NiFe.hmm",challenge,args.out,"original_challenge")
    rows.extend(dict(experiment="old_global_model_development_challenge",**d) for d in original)
    hmm,manifest=build_profile(training,args.out,"NiFe_Group4")
    models.append(manifest)
    rows.extend(dict(experiment="branch_core_model_development_challenge",**d) for d in search(hmm,challenge,args.out,"branch_challenge"))
    queries=args.out/"group4_and_decoys.faa"
    SeqIO.write([*training,*[fasta[r["reference_id"]] for r in refs if r["kind"]=="decoy"]],queries,"fasta")
    rows.extend(dict(experiment="branch_core_in_library_and_decoys",**d) for d in search(hmm,queries,args.out,"branch_seeds"))
    for ref in selected:
        if ref["kind"]!="catalytic":
            continue
        kept=[fasta[r["reference_id"]] for r in selected if r["sample"]!=ref["sample"]]
        label="leave_source_out_"+ref["structure_id"]
        hmm,manifest=build_profile(kept,args.out/label,"Group4_LOSO")
        models.append(dict(excluded_sample=ref["sample"],**manifest))
        query=args.out/label/"query.faa"
        SeqIO.write(fasta[ref["reference_id"]],query,"fasta")
        rows.extend(dict(experiment=label,**d) for d in search(hmm,query,args.out/label,"heldout"))
    blast=args.out/"challenge_vs_all_references.tsv"
    subprocess.run(["blastp","-query",str(challenge),"-subject",str(ROOT/"references/processed/references.faa"),
        "-max_hsps","1","-outfmt","6 qseqid sseqid pident qcovhsp qlen slen length bitscore evalue","-out",str(blast)],check=True)
    write_tsv(args.out/"domains.tsv",rows)
    report=dict(status="diagnostic_complete",challenge=frozen,profiles=models,domains=rows,
        limitations="Five source-withheld seed controls and one external development challenge are not a 20-sample independent Group4 acceptance set. HMM hit alone cannot establish catalytic activity or function.")
    (args.out/"summary.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(dict(models=len(models),domains=len(rows),challenge=cat["id"])))


if __name__=="__main__":
    main()
