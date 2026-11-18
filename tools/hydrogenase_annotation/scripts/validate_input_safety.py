#!/usr/bin/env python3
"""Execute real searches on reference-derived boundary controls; not an accuracy benchmark."""
import argparse
import copy
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, extract_genbank, write_tsv


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=ROOT/"results/input_safety")
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    refs=json.loads((ROOT/"references/processed/references.json").read_text())
    genomes={}
    checks=[]
    label_tables=["complexes.tsv","catalytic_anchors.tsv","catalytic_sites.tsv","complex_members.tsv",
        "reaction_capabilities.tsv","maturation.tsv","domain_hits.tsv","annotation_evidence.tsv",
        "fes_sites.tsv","fes_ligands.tsv","fes_template_evidence.tsv","active_center_summary.tsv",
        "active_center_template_evidence.tsv","seed_diagnostic_sites.tsv","catalytic_confirmation.tsv",
        "family_candidates.tsv","sequence_evidence.tsv","rejected.tsv","genes.tsv","module_graph_nodes.tsv",
        "module_graph_edges.tsv","cross_anchor_associations.tsv","catalytic_center_candidates.tsv","catalytic_center_components.tsv"]

    def subset(sid):
        selected=[r for r in refs if r["structure_id"]==sid]
        sample=selected[0]["sample"]
        if sample not in genomes:
            genomes[sample]=extract_genbank(ROOT/"references/raw"/f"{sample}.gb")[0]
        ids={r["locus_tag"] for r in selected}
        return copy.deepcopy([g for g in genomes[sample] if g["id"] in ids])

    def run(name,genes,mode,ownership=None,transcript=False):
        directory=args.out/name
        directory.mkdir(parents=True,exist_ok=True)
        faa,gff=directory/"input.faa",directory/"input.gff"
        SeqIO.write([SeqRecord(Seq(g["sequence"]),id=g["id"],description="") for g in genes],faa,"fasta")
        cmd=[sys.executable,str(ROOT/"scripts/annotate.py"),"--proteins",str(faa),"--mode",mode,"--out",str(directory)]
        if not transcript:
            gff.write_text("##gff-version 3\n"+"".join(
                f"{g['contig']}\tboundary_control\tCDS\t{g['start']}\t{g['end']}\t.\t{'+' if g['strand']==1 else '-'}\t0\tID={g['id']};partial={'11' if g['partial'] else '00'}\n"
                for g in genes))
            cmd.extend(["--gff",str(gff)])
        if ownership is not None:
            bins=directory/"bins.tsv"
            write_tsv(bins,[dict(contig=c,bin_id=b) for c,b in ownership.items()],["contig","bin_id"])
            cmd.extend(["--bin-map",str(bins)])
        with (directory/"command.log").open("w") as log:
            result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Boundary input failed: {name}; see {directory/'command.log'}")
        return directory,read_tsv(directory/"complexes.tsv")

    def record(name,condition,observed):
        checks.append(dict(test=name,passed=bool(condition),observed=observed,
            evidence_type="reference_derived_control_not_independent_performance"))

    for sid,remote_roles in [("NiFe_FHL1_Hyc",["FdhF"]),("NiFe_MvhHdr",["HdrA","HdrB","HdrC"])]:
        genes=subset(sid)
        remotes={r["locus_tag"] for r in refs if r["structure_id"]==sid and r["role"] in remote_roles}
        for g in genes:
            g["contig"]="partner_contig" if g["id"] in remotes else "core_contig"
        cases=[("genome",None,True),("assembly",None,False),
               ("same_bin",{"core_contig":"binA","partner_contig":"binA"},True),
               ("cross_bin",{"core_contig":"binA","partner_contig":"binB"},False),
               ("unbinned",{"core_contig":".","partner_contig":"."},False)]
        for label,bins,expected in cases:
            directory,annotations=run(sid+"_"+label,genes,"bin" if bins is not None else label,bins)
            a=next(a for a in annotations if a["structure_id"]==sid)
            assigned=bool(a["potential_function_ids"])
            record(sid+":"+label,assigned==expected and a["catalytic_status"]=="supported",
                dict(function_assigned=assigned,structure=a["structure_status"],missing=a["missing_roles"]))
            if label=="cross_bin":
                members=read_tsv(directory/"complex_members.tsv")
                record(sid+":no_cross_bin_members",not any(m["role"] in remote_roles for m in members),
                    [m["role"] for m in members])
            repeated,_=run(sid+"_"+label+"_repeat",genes,"bin" if bins is not None else label,bins)
            record(sid+":"+label+":repeatability",all((directory/t).read_bytes()==(repeated/t).read_bytes() for t in label_tables),label_tables)

    hya=subset("NiFe_Ecoli_Hyd1")
    ref=next(r for r in refs if r["structure_id"]=="NiFe_Ecoli_Hyd1" and r["kind"]=="catalytic")
    catalytic=next(g for g in hya if g["id"]==ref["locus_tag"])
    cutoff=ref["diagnostic_positions"][1]+5
    catalytic["sequence"]=catalytic["sequence"][cutoff:]
    catalytic["partial"]=True
    directory,annotations=run("N_terminal_truncation",hya,"genome")
    sites=read_tsv(directory/"catalytic_sites.tsv")
    record("truncation_not_substitution",sum(s["status"]=="not_covered" for s in sites)>=2 and
        not any(s["status"]=="substituted" for s in sites) and not any(a["potential_function_ids"] for a in annotations),
        [s["status"] for s in sites])

    hys=subset("NiFe_periplasmic_Hys")
    ref=next(r for r in refs if r["structure_id"]=="NiFe_periplasmic_Hys" and r["kind"]=="catalytic")
    catalytic=next(g for g in hys if g["id"]==ref["locus_tag"])
    original,original_annotations=run('real_Sec_preserved',copy.deepcopy(hys),'genome')
    original_sites=read_tsv(original/'catalytic_sites.tsv')
    sec_positions={str(p) for p in ref['selenocysteine_positions']}
    record('real_Sec_coordinates_preserved_in_complete_model',any(s['residue_role']=='direct_metal_ligand' and
        s['query_position'] in sec_positions and s['observed_residue']=='U' and s['status']=='present' for s in original_sites),
        [dict(query_position=s['query_position'],observed=s['observed_residue'],status=s['status']) for s in original_sites if s['residue_role']=='direct_metal_ligand'])
    catalytic["sequence"]=catalytic["sequence"].replace("U","C")
    directory,annotations=run("Sec_to_Cys",hys,"genome")
    seed_sites=read_tsv(directory/'seed_diagnostic_sites.tsv');sites=read_tsv(directory/"catalytic_sites.tsv")
    record("Sec_Cys_distinction",any(s["expected_residue"]=="U" and s["status"]=="substituted" for s in seed_sites)
        and any(s['residue_role']=='direct_metal_ligand' and s['query_position'] in sec_positions and s['observed_residue']=='C' for s in sites)
        and not any(a["potential_function_ids"] for a in annotations),dict(reference_specific=seed_sites,
            complete_model=[s for s in sites if s['residue_role']=='direct_metal_ligand']))

    hya=subset("NiFe_Ecoli_Hyd1")
    ref=next(r for r in refs if r["structure_id"]=="NiFe_Ecoli_Hyd1" and r["kind"]=="catalytic")
    catalytic=next(g for g in hya if g["id"]==ref["locus_tag"])
    seq=list(catalytic["sequence"]);seq[ref["diagnostic_positions"][0]-1]="X";catalytic["sequence"]="".join(seq)
    directory,annotations=run("unknown_ligand",hya,"genome")
    sites=read_tsv(directory/"catalytic_sites.tsv")
    record("unknown_not_substitution",any(s["status"]=="ambiguous" for s in sites) and
        not any(s["status"]=="substituted" for s in sites) and not any(a["potential_function_ids"] for a in annotations),sites)

    for sample,selector in [("pfuriosus_DSM3638","PF_RS04490"),("tonnurineus_NA1","TON_RS02775")]:
        genes=[g for g in extract_genbank(ROOT/"references/raw"/f"{sample}.gb")[0] if g["id"]==selector]
        directory,annotations=run("homonym_"+selector,genes,"genome")
        maturation=read_tsv(directory/"maturation.tsv")
        record("homonym_"+selector,not any(m["maturation_gene"]=="hydG" for m in maturation) and not annotations,maturation)

    directory,annotations=run("algal_transcript",subset("FeFe_CrHydA1"),"transcript",transcript=True)
    a=next(a for a in annotations if a["structure_id"]=="FeFe_CrHydA1")
    record("transcript_no_genomic_coordinates",a["catalytic_status"]=="supported" and a["start"]=="" and
        a["end"]=="" and a["operon_status"]=="not_applicable_translated_transcript",a)
    directory2,annotations2=run("algal_transcript_repeat",subset("FeFe_CrHydA1"),"transcript",transcript=True)
    record("fixed_input_repeatability",all((directory/t).read_bytes()==(directory2/t).read_bytes() for t in label_tables),label_tables)

    write_tsv(args.out/"checks.tsv",checks)
    summary=dict(status="passed" if all(c["passed"] for c in checks) else "failed",total=len(checks),
        passed=sum(c["passed"] for c in checks),failed=[c for c in checks if not c["passed"]],
        annotation_script_sha256=hashlib.sha256((ROOT/"scripts/annotate.py").read_bytes()).hexdigest(),
        limitations="Reference-derived boundary controls, repeated across four supported modes; not independent accuracy or a large-MAG benchmark.")
    (args.out/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    if summary["failed"]:
        sys.exit(1)


if __name__=="__main__":
    main()
