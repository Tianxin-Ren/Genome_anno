#!/usr/bin/env python3
"""Real reference retrieval, near-strain holdout, and explicit perturbation controls."""
import copy
import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from Bio import SeqIO

from annotate import ROOT, annotate, extract_genbank, write_tsv


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,default=ROOT/"results/validation")
    args=parser.parse_args()
    systems=json.loads((ROOT/"config/systems.json").read_text())
    refs=json.loads((ROOT/"references/processed/references.json").read_text())
    output=args.out
    output.mkdir(parents=True,exist_ok=True)
    (output/"summary.json").write_text(json.dumps(dict(status="running"))+"\n")
    genomes={name:extract_genbank(ROOT/"references/raw"/f"{name}.gb")[0] for name in sorted({s["sample"] for s in systems})}
    genomes["ecoli_W3110_holdout"]=extract_genbank(ROOT/"references/raw/ecoli_W3110_holdout.gb")[0]

    def run(item):
        name,genes=item
        # Strip reference names and descriptions: classification must use sequence and coordinates.
        genes=copy.deepcopy(genes)
        for g in genes:g.update(gene="",aliases="",protein_accession="")
        result=annotate(genes,output/name,input_files=[ROOT/"references/raw"/f"{name}.gb"],mode="genome")
        print(name,len(result),flush=True)
        return name,result
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=dict(pool.map(run,genomes.items()))
    checks=[]
    for label,sample,gid,accession in [("CpI","cpasteurianum_DSM525","CLPA_c00280","P29166"),
                                       ("HydA1_cDNA","creinhardtii_HydA1_cDNA","AAG00591.1","Q9FYU1")]:
        protein=next(g for g in genomes[sample] if g["id"]==gid)
        experimental=json.loads((ROOT/"references/raw"/f"uniprot_{accession}.json").read_text())
        checks.append(dict(test="experimental_sequence_match:"+label,type="reference_identity_mapping",passed=protein["sequence"]==experimental["sequence"]["value"],
            observed=accession,limitation="Confirms reference sequence identity; not an independent prediction-accuracy test"))
    for system in systems:
        catalytic=next(r for r in refs if r["structure_id"]==system["structure_id"] and r["kind"]=="catalytic")
        found=next((r for r in results[system["sample"]] if r["protein_id"]==catalytic["locus_tag"]),{})
        ok=(found.get("structure_id")==system["structure_id"] and found.get("structure_confidence")=="high" and
            set(found.get("potential_function_ids","").split(";"))==set(system["function_ids"]))
        checks.append(dict(test="reference_retrieval:"+system["structure_id"],type="in_library_smoke_test",passed=ok,
            observed=found.get("structure_confidence","missing")+":"+found.get("potential_function_ids","") ,
            limitation="Reference sequences are in training library; not independent sensitivity/specificity"))
    for sid in ["NiFe_Ecoli_Hyd1","NiFe_Ecoli_Hyd2","NiFe_FHL1_Hyc"]:
        matched=[r for r in results["ecoli_W3110_holdout"] if r["structure_id"]==sid and r.get("structure_confidence")=="high"]
        checks.append(dict(test="near_strain_holdout:"+sid,type="held_out_genome",passed=len(matched)==1,
            observed=len(matched),limitation="W3110 genome excluded from seed library; K-12 proteins may be identical; not remote generalization"))

    def genes_for(sid):
        s=next(s for s in systems if s["structure_id"]==sid)
        ids={r["locus_tag"] for r in refs if r["structure_id"]==sid}
        return copy.deepcopy([g for g in genomes[s["sample"]] if g["id"] in ids])

    controls=[]
    huc=genes_for("NiFe_HucLMS")
    controls.append(("missing_HucM",[g for g in huc if g["id"]!="MSMEG_RS10955"]))
    hyb=genes_for("NiFe_Ecoli_Hyd2")
    controls.append(("missing_HybA",[g for g in hyb if g["gene"]!="hybA"]))
    broken=copy.deepcopy(huc)
    cat=next(g for g in broken if g["id"]=="MSMEG_RS10965")
    positions=next(r["diagnostic_positions"] for r in refs if r["structure_id"]=="NiFe_HucLMS" and r["kind"]=="catalytic")
    sequence=list(cat["sequence"])
    for p in positions:sequence[p-1]="A"
    cat["sequence"]="".join(sequence)
    controls.append(("mutated_NiFe_ligands",broken))
    small=copy.deepcopy(huc)
    subunit=next(g for g in small if g["id"]=="MSMEG_RS10960")
    subunit["sequence"]=subunit["sequence"].replace("C","A")
    controls.append(("mutated_small_subunit_cysteines",small))
    split=copy.deepcopy(huc)
    split[0]["contig"]="separate_contig"
    controls.append(("split_local_module",split))
    anti=copy.deepcopy(huc)
    anti[0]["strand"]*=-1
    controls.append(("opposite_strand_module",anti))
    truncated=copy.deepcopy(huc)
    next(g for g in truncated if g["id"]=="MSMEG_RS10965")["partial"]=True
    controls.append(("partial_catalytic_CDS",truncated))
    mvh=genes_for("NiFe_MvhHdr")
    hdrids={r["locus_tag"] for r in refs if r["structure_id"]=="NiFe_MvhHdr" and r["role"].startswith("Hdr")}
    controls.append(("missing_remote_Hdr",[g for g in mvh if g["id"] not in hdrids]))
    fhl=genes_for("NiFe_FHL1_Hyc")
    controls.append(("missing_remote_FdhF",[g for g in fhl if g["gene"]!="fdhF"]))
    tm=genes_for("FeFe_TmHydABC")
    sw=genes_for("FeFe_SwHyd1ABC")
    fefe_cys=copy.deepcopy(tm)
    cat=next(g for g in fefe_cys if g["gene"]=="hydA")
    cat["sequence"]=cat["sequence"].replace("C","A",1)
    controls.append(("mutated_FeFe_reference_cysteine",fefe_cys))
    hybrid=[g for g in tm if g["gene"]!="hydB"]
    oldb=next(g for g in tm if g["gene"]=="hydB")
    swb=next(g for g in sw if g["id"]=="SWOL_RS05170")
    hybrid.append({**oldb,"sequence":swb["sequence"]})
    controls.append(("incompatible_FeFe_beta_module",hybrid))
    for sid,role in [("NiFe_HynABC","HydC"),("NiFe_VhoABC","VhtC"),("NiFe_H2_sensor_2b","HupT"),
                     ("NiFe_HyhNADPH","SHI_gamma"),("NiFe_HoxFUYH","HoxF"),
                     ("NiFe_Tonn_Fdh2Mfh2Mnh2","Fdh2_A"),("NiFe_Tonn_Fdh2Mfh2Mnh2","Mnh2_TON1580"),
                     ("NiFe_Coo_CODH","CooS"),("FeFe_HDCR_Aw","FdhF2"),("FeFe_HytABCDE","FdhA")]:
        gid=next(r["locus_tag"] for r in refs if r["structure_id"]==sid and r["role"]==role)
        controls.append(("missing_"+sid+"_"+role,[g for g in genes_for(sid) if g["id"]!=gid]))
    sec=genes_for("NiFe_periplasmic_Hys")
    for gene in sec:gene["sequence"]=gene["sequence"].replace("U","C")
    controls.append(("Hys_Sec_changed_to_Cys_requires_review",sec))
    for kind in ["decoy","maturation"]:
        selected=[]
        for sample,genes in genomes.items():
            if sample.endswith("holdout"):continue
            ids={r["locus_tag"] for r in refs if r["sample"]==sample and r["kind"]==kind}
            selected.extend(copy.deepcopy([g for g in genes if g["id"] in ids]))
        controls.append((kind+"_only",selected))
    for name,genes in controls:
        found=annotate(genes,output/"controls"/name,mode="genome")
        high=[r for r in found if r.get("potential_function_ids")]
        ok=not high
        if name=="maturation_only":ok=ok and not found
        if name=="decoy_only":ok=ok and all(r['catalytic_family']=='Fe' and r['catalytic_status']=='candidate' and
            r['structure_id']=='unresolved' and r['structure_confidence']=='low' and
            'Hmd_paralog_competitor_blocks_function_transfer_not_discovery' in r['flags'] for r in found)
        checks.append(dict(test=name,type="perturbation_control",passed=ok,observed=len(high),
            limitation="Synthetic deletion/mutation/mixture of real reference CDSs; not a biological knockout assay"))
    hup=genes_for("NiFe_HupSL")
    context_ids={r["locus_tag"] for r in refs if r["structure_id"]=="NiFe_HupSL" and r["kind"]=="context"}
    result=annotate([g for g in hup if g["id"] not in context_ids],output/"controls/Hup_without_nitrogen_fixation_context",mode="genome")
    ok=any(a["structure_id"]=="NiFe_HupSL" and a["structure_confidence"]=="high" and not a.get("potential_function_ids") for a in result)
    checks.append(dict(test="Hup_without_nif_context",type="physiological_context_control",passed=ok,
        observed="Hup structure retained; nitrogenase-recycling function withheld",limitation="Nif presence supports context, not measured H2 recycling flux"))
    for sample,gid,label in [("mbarkeri_Fusaro","MBAR_RS09780","Vhx_not_Vht"),
                             ("cnecator_H16_pHG1","H16_RS31660","H16_sensor_not_functionally_validated")]:
        found=next((r for r in results[sample] if r["protein_id"]==gid),{})
        checks.append(dict(test=label,type="real_paralog_control",passed=found.get("structure_id")=="unresolved" and not found.get("potential_function_ids"),
            observed=found.get("structure_id","missing"),limitation="Native reference paralogs remain unresolved; no claim that the catalytic protein itself is inactive"))
    fna=output/"ecoli_MG1655.fna"
    SeqIO.write(SeqIO.read(ROOT/"references/raw/ecoli_MG1655.gb","genbank"),fna,"fasta")
    predicted=output/"from_genome"
    subprocess.run([sys.executable,str(ROOT/"scripts/annotate.py"),"--genome",str(fna),"--mode","genome","--out",str(predicted)],check=True)
    from_gff=output/"from_faa_gff"
    subprocess.run([sys.executable,str(ROOT/"scripts/annotate.py"),"--proteins",str(predicted/"prodigal.faa"),
                    "--gff",str(predicted/"prodigal.gff"),"--mode","genome","--out",str(from_gff)],check=True)
    for name,directory in [("genome_Prodigal_entry",predicted),("faa_gff_entry",from_gff)]:
        with (directory/"complexes.tsv").open() as handle:found=list(csv.DictReader(handle,delimiter="\t"))
        detected={r["structure_id"] for r in found}
        fhl=next((r for r in found if r["structure_id"]=="NiFe_FHL1_Hyc"),{})
        ok={"NiFe_Ecoli_Hyd1","NiFe_Ecoli_Hyd2","NiFe_FHL1_Hyc"}<=detected and not fhl.get("potential_function_ids")
        checks.append(dict(test=name,type="end_to_end_input",passed=ok,observed="3 structures; FHL function withheld",
            limitation="Prodigal truncates Sec-containing FdhF; supplied curated CDS proteins are preferred for recoded genes"))
    from phylogeny import build_trees
    for sample,family in [("ecoli_W3110_holdout","NiFe"),("tmaritima_MSB8","FeFe")]:
        build_trees(output/sample)
        tree=output/sample/"phylogeny"/f"{family}.treefile"
        checks.append(dict(test="phylogeny:"+family,type="tree_execution",passed=tree.exists() and tree.stat().st_size>0,
            observed="LG+G4; 1000 ultrafast bootstraps",limitation="Sparse seed sampling; review tree only; does not independently establish HydDB subgroup or function"))
    write_tsv(output/"checks.tsv",checks)
    summary=dict(status="passed" if all(r["passed"] for r in checks) else "failed",total_checks=len(checks),passed=sum(r["passed"] for r in checks),
        failed=[r for r in checks if not r["passed"]],reference_structures=len(systems),heldout_genomes=1,
        independent_generalization_accuracy="not_estimated")
    (output/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    if summary["failed"]:sys.exit(1)


if __name__=="__main__":main()
