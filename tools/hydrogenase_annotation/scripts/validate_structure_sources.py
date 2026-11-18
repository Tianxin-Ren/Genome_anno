#!/usr/bin/env python3
"""Validate structural extraction against declared chemical/site examples, not query accuracy."""
import json
import sys
from annotate import ROOT, write_tsv


def main():
    directory=ROOT/"references/structural_sites"
    templates=json.loads((directory/"templates.json").read_text())
    reserved=json.loads((directory/"reserved_validation.json").read_text())
    by_id={r["pdb_id"]:r for r in templates}
    checks=[]
    def check(name,passed,observed):
        checks.append(dict(test=name,passed=bool(passed),observed=observed,evidence_scope="structural_source_preparation_not_sequence_prediction_accuracy"))
    check("reserved_structures_never_in_templates",not {r["pdb_id"] for r in templates}&{r["pdb_id"] for r in reserved},
        dict(templates=len(templates),reserved=len(reserved)))
    for pdb,expected in {"3UQY":3,"6EHQ":3,"7UTD":3,"7P92":10,"1FEH":5,"4CI0":4,"9E1J":3,"1CC1":3}.items():
        check("canonical_site_count:"+pdb,len(by_id[pdb]["sites"])==expected,len(by_id[pdb]["sites"]))
    hya=by_id["3UQY"]["sites"]
    check("Hya_special_proximal_six_Cys",any(s["cluster_type"]=="[4Fe-3S]" and
        sorted(r["residue"] for r in s["ligand_residues"])==["C"]*6 for s in hya),"F4S with six sequence-checked cysteines")
    check("Hya_distal_His_ligand",any(s["cluster_type"]=="[4Fe-4S]" and
        sorted(r["residue"] for r in s["ligand_residues"])==["C","C","C","H"] for s in hya),"Histidine retained, not a Cys-count-only rule")
    check("Huc_three_threeFe_clusters",len(by_id["7UTD"]["sites"])==3 and
        all(s["cluster_type"]=="[3Fe-4S]" and len(s["ligand_residues"])==3 for s in by_id["7UTD"]["sites"]),
        "Three F3S components, each with three primary Cys ligands")
    check('NiFeSe_small_three_fourFe_clusters',len(by_id['1CC1']['sites'])==3 and
        all(s['cluster_type']=='[4Fe-4S]' and len(s['ligand_residues'])==4 for s in by_id['1CC1']['sites']),
        'Reviewed1CC1smallsubunit;three4Fe4Sclusters;reference replay only, not independent performance')
    frhb=next(s for s in by_id['4CI0']['sites'] if s['site_id']=='4CI0__K_1282_SF4')
    check('FrhB_four_ligands_include_coordinate_supplement',[(r['position'],r['residue']) for r in frhb['ligand_residues']]==[(104,'C'),(134,'C'),(192,'C'),(195,'C')],
        'Cys104 omitted from deposited connections is corroborated by exact-sequence 4OMF; 4CI0 short-contact warning remains in site review')
    check("long_Hyb_Asp_contact_kept_separate",not any(r["position"]==82 and r["residue"]=="D"
        for s in by_id["6EHQ"]["sites"] for r in s["ligand_residues"]) and
        any(c["position"]==82 and c["residue"]=="D" for s in by_id["6EHQ"]["sites"] for c in s["contacts"]),
        "Deposited 2.770-A Asp contact is preserved but not forced into first shell")
    check("all_query_occupancies_unknown",all(s["query_cofactor_occupancy"]=="unknown" for r in templates+reserved for s in r["sites"]),
        "Experimental reference occupancy never becomes query occupancy")
    write_tsv(directory/"checks.tsv",checks)
    summary=dict(status="passed" if all(c["passed"] for c in checks) else "failed",total=len(checks),passed=sum(c["passed"] for c in checks),
        failed=[c for c in checks if not c["passed"]],
        scope="Source extraction only. Templates still require coordination review and sequence predictor implementation; no >=95% topology precision claim.")
    (directory/"validation_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    if summary["failed"]:
        sys.exit(1)


if __name__=="__main__":
    main()
