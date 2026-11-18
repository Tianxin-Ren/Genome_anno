#!/usr/bin/env python3
"""Catalytic reference/query trees for review, without automatic function upgrades."""
import argparse
import csv
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from annotate import ROOT,write_tsv


def trim_to_reference_core(records,reference_ids):
    lengths={len(r.seq) for r in records}
    if len(lengths)!=1:raise ValueError("Unequal sequence lengths in alignment")
    refs=[str(r.seq) for r in records if r.id in reference_ids]
    if not refs:raise ValueError("No reference sequences in tree alignment")
    columns=[i for i in range(next(iter(lengths))) if sum(s[i]!="-" for s in refs)/len(refs)>=0.6]
    return [SeqRecord(Seq("".join(str(r.seq)[i] for i in columns)),id=r.id,description="") for r in records],len(columns)


def build_trees(result,threads=2):
    result=Path(result).resolve()
    with (result/"complexes.tsv").open() as handle:
        candidates=list(csv.DictReader(handle,delimiter="\t"))
    proteins=SeqIO.to_dict(SeqIO.parse(result/"proteins.faa","fasta"))
    out=result/"phylogeny"
    out.mkdir(exist_ok=True)
    statuses=[]
    for family in ["NiFe","FeFe","Fe"]:
        genes=[r["protein_id"] for r in candidates if r["catalytic_family"]==family]
        if not genes:continue
        ref=ROOT/"references/processed"/f"{family}.aligned.faa"
        ref_records=list(SeqIO.parse(ref,"fasta"))
        if len({str(r.seq) for r in ref_records})<4:
            statuses.append(dict(family=family,status="not_run_fewer_than_4_distinct_reference_sequences",tree=""))
            continue
        query=out/f"{family}.queries.faa"
        labels=[dict(tree_label=f"query_{i+1}",protein_id=gid) for i,gid in enumerate(genes)]
        SeqIO.write([SeqRecord(Seq(str(proteins[r["protein_id"]].seq).replace("U","X")),id=r["tree_label"],description="") for r in labels],query,"fasta")
        write_tsv(out/f"{family}.query_labels.tsv",labels)
        alignment=out/f"{family}.reference_plus_query.faa"
        with alignment.open("w") as handle, (out/f"{family}.mafft.log").open("w") as log:
            subprocess.run(["mafft","--quiet","--add",str(query),"--keeplength",str(ref)],stdout=handle,stderr=log,check=True)
        trimmed,sites=trim_to_reference_core(list(SeqIO.parse(alignment,"fasta")),{r.id for r in ref_records})
        if sites<100:raise ValueError(f"Only {sites} reference-core sites for {family}; review alignment")
        core=out/f"{family}.reference_core.faa"
        SeqIO.write(trimmed,core,"fasta")
        prefix=out/family
        with (out/f"{family}.iqtree.stdout.log").open("w") as log:
            subprocess.run(["iqtree3","-s",str(core),"-m","LG+G4","-B","1000","-T",str(threads),
                            "--seed","1729","--prefix",str(prefix),"--redo"],stdout=log,stderr=log,check=True)
        statuses.append(dict(family=family,status=f"review_tree_only; LG+G4; 1000_UFBoot; unrooted; sparse_reference_sampling; reference_occupancy>=60%; sites={sites}",tree=str(prefix)+".treefile"))
    write_tsv(out/"status.tsv",statuses,["family","status","tree"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=2)
    args=parser.parse_args()
    build_trees(args.result,args.threads)
