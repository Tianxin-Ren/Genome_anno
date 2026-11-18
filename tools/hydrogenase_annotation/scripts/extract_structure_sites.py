#!/usr/bin/env python3
"""Extract experimental Fe-S components and explicit mmCIF metal connections.

This is structural source preparation, not a sequence predictor. Reserved
validation structures are written separately and never enter templates.json.
"""
import hashlib
import json
from math import dist
from collections import defaultdict
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from annotate import ROOT, write_tsv


CLUSTER_TYPES={"SF4":"[4Fe-4S]","F3S":"[3Fe-4S]","F4S":"[4Fe-3S]","FES":"[2Fe-2S]"}


def coordinate_connections(sequences,atoms):
    """Deposited metalc tables can omit bonds or entire modeled Fe-S clusters."""
    metals=defaultdict(list)
    for atom in atoms:
        if atom['component'] in CLUSTER_TYPES and atom['element']=='FE':
            metals[atom['chain'],atom['author_position'],atom['component']].append(atom)
    donors=[a for a in atoms if a['entity'] in sequences and a['position'].isdigit() and a['element'] in ['S','SE','N','O']]
    result=defaultdict(list)
    for key,cluster in metals.items():
        result[key]
        for protein in donors:
            for metal in cluster:
                distance=dist(protein['xyz'],metal['xyz'])
                if distance>(2.8 if protein['element'] in ['S','SE'] else 2.6):continue
                position=int(protein['position'])
                result[key].append(dict(chain=protein['chain'],entity_id=protein['entity'],position=position,
                    author_position=protein['author_position'],residue=sequences[protein['entity']][position-1],atom=protein['atom'],
                    metal_atom=metal['atom'],distance_angstrom=distance,symmetry='1_555',metal_symmetry='1_555',
                    evidence_origin='first_model_positive_occupancy_ASU_coordinates'))
    return result


def primary_contact(contact):
    """Keep longer deposited contacts as evidence, not as proven first-shell ligands."""
    limit=2.8 if contact["atom"].startswith("S") else 2.6
    return (contact["distance_angstrom"] is not None and contact["distance_angstrom"]<=limit and
            contact["symmetry"]==contact["metal_symmetry"]=="1_555")


def extract(path):
    d=MMCIF2Dict(str(path))
    pdb=d["_entry.id"][0]
    sequences={entity:sequence.replace("\n","").replace(" ","") for entity,sequence in
        zip(d["_entity_poly.entity_id"],d["_entity_poly.pdbx_seq_one_letter_code_can"])}
    chains=dict(zip(d["_struct_asym.id"],d["_struct_asym.entity_id"]))
    chem={cid:dict(name=name,formula=formula) for cid,name,formula in
        zip(d["_chem_comp.id"],d["_chem_comp.name"],d["_chem_comp.formula"])}
    components=defaultdict(list)
    for i,kind in enumerate(d.get("_struct_conn.conn_type_id",[])):
        if kind!="metalc":
            continue
        sides=[]
        for n in [1,2]:
            prefix=f"_struct_conn.ptnr{n}_"
            sides.append(dict(chain=d[prefix+"label_asym_id"][i],component=d[prefix+"label_comp_id"][i],
                seq_id=d[prefix+"label_seq_id"][i],atom=d[prefix+"label_atom_id"][i],
                author_seq_id=d[prefix+"auth_seq_id"][i],symmetry=d[prefix+"symmetry"][i]))
        for metal,protein in [sides,list(reversed(sides))]:
            if metal["component"] not in CLUSTER_TYPES or not protein["seq_id"].isdigit():
                continue
            entity=chains[protein["chain"]]
            if entity not in sequences:
                continue
            position=int(protein["seq_id"])
            expected=seq1(protein["component"],custom_map={"SEC":"U"})
            sequence=sequences[entity]
            if position>len(sequence) or sequence[position-1]!=expected:
                raise ValueError(f"mmCIF sequence/contact mismatch: {pdb} {protein} {expected}")
            key=(metal["chain"],metal["author_seq_id"],metal["component"])
            distance=d.get("_struct_conn.pdbx_dist_value",["?"]*len(d["_struct_conn.conn_type_id"]))[i]
            components[key].append(dict(chain=protein["chain"],entity_id=entity,position=position,
                author_position=protein["author_seq_id"],residue=expected,atom=protein["atom"],metal_atom=metal["atom"],
                distance_angstrom=float(distance) if distance not in ["?","."] else None,
                symmetry=protein["symmetry"],metal_symmetry=metal["symmetry"]))
    supplements=0;atoms=[]
    if '_atom_site.id' in d:
        from active_centers import structure_atoms
        atom_sequences,atoms=structure_atoms(path)
        # entity_poly's canonical sequence resolves modified residues such as DHI -> H; chirality is not transferred.
        if {e:len(s) for e,s in atom_sequences.items()}!={e:len(s) for e,s in sequences.items()}:
            raise ValueError('Canonical FeS polymer coordinates disagree')
        for key,contacts in coordinate_connections(sequences,atoms).items():
            existing={(c['chain'],c['position'],c['atom'],c['metal_atom']) for c in components[key] if primary_contact(c)}
            for contact in contacts:
                identity=(contact['chain'],contact['position'],contact['atom'],contact['metal_atom'])
                if identity not in existing:
                    components[key].append(contact);existing.add(identity);supplements+=1
    sites=[];raw_sites=[];seen=set()
    canonical={entity:min(c for c,e in chains.items() if e==entity) for entity in sequences}
    modeled={entity:sorted({int(a['position']) for a in atoms if a['entity']==entity and a['chain']==canonical[entity] and a['position'].isdigit()}) for entity in sequences}
    ligand_atoms={entity:{} for entity in sequences}
    for a in atoms:
        if a['entity'] in sequences and a['chain']==canonical[a['entity']] and a['position'].isdigit() and a['element'] in ['S','SE','N','O']:
            names=ligand_atoms[a['entity']].setdefault(a['position'],set());names.add(a['atom'])
    ligand_atoms={e:{p:sorted(names) for p,names in positions.items()} for e,positions in ligand_atoms.items()}
    for (chain,author_position,cid),contacts in sorted(components.items()):
        ligands=sorted({(c["chain"],c["entity_id"],c["position"],c["residue"]) for c in contacts if primary_contact(c)})
        site=dict(site_id=f"{pdb}__{chain}_{author_position}_{cid}",pdb_id=pdb,
            component_id=cid,cluster_type=CLUSTER_TYPES[cid],chemical_name=chem[cid]["name"],chemical_formula=chem[cid]["formula"],
            ligand_residues=[dict(chain=c,entity_id=e,position=p,residue=r) for c,e,p,r in ligands],
            contacts=contacts,polypeptide_chains=sorted({c[0] for c in ligands}),
            source_evidence="deposited_connections_supplemented_by_experimental_ASU_first_shell_coordinates",
            first_shell_screen="S/Se <=2.8 A; N/O <=2.6 A; identity symmetry only; all deposited contacts retained",
            coordination_review="requires_template_review_before_prediction",
            query_cofactor_occupancy="unknown")
        raw_sites.append(site)
        # Choose a declared chain per entity, not whichever copy gives the desired ligands.
        if not ligands or any(c!=canonical[e] for c,e,p,r in ligands):
            continue
        signature=(cid,tuple(sorted({(e,p,r) for c,e,p,r in ligands})))
        if signature in seen:
            continue
        seen.add(signature)
        sites.append(site)
    return dict(pdb_id=pdb,method=d.get("_exptl.method",[]),sequences=sequences,
        chain_entities={c:e for c,e in chains.items() if e in sequences},canonical_chains=canonical,sites=sites,raw_site_instances=raw_sites,
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),coordinate_supplemented_contacts=supplements,
        modeled_protein_positions=modeled if atoms else None,modeled_ligand_atoms=ligand_atoms if atoms else None)


def main():
    config=json.loads((ROOT/"config/structure_sources.json").read_text())
    raw=ROOT/"references/raw/structures"
    downloads={r["pdb_id"]:r for r in json.loads((raw/"download_manifest.json").read_text())}
    output=ROOT/"references/structural_sites";output.mkdir(exist_ok=True)
    templates=[];validation=[];summary=[]
    for source in config["sources"]:
        path=raw/(source["pdb_id"]+".cif")
        if hashlib.sha256(path.read_bytes()).hexdigest()!=downloads[source["pdb_id"]]["sha256"]:
            raise ValueError(f"Structure source changed: {path}")
        record={**extract(path),**source}
        (templates if source["split"]=="template" else validation).append(record)
        summary.append(dict(pdb_id=source["pdb_id"],split=source["split"],
            unique_sites=len(record["sites"]),ligand_positions=sum(len(s["ligand_residues"]) for s in record["sites"]),
            cluster_types=";".join(sorted({s["cluster_type"] for s in record["sites"]}))))
    assert not {r["pdb_id"] for r in templates}&{r["pdb_id"] for r in validation}
    for name,records in [("templates",templates),("reserved_validation",validation)]:
        (output/(name+".json")).write_text(json.dumps(records,indent=2)+"\n")
    write_tsv(output/"extraction_summary.tsv",summary)
    print(json.dumps(dict(template_structures=len(templates),reserved_structures=len(validation),sources=summary,
        status="source_extraction_only_not_prediction_acceptance"),indent=2))


if __name__=="__main__":
    main()
