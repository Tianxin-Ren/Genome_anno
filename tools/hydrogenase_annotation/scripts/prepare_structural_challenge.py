#!/usr/bin/env python3
"""Freeze source-selected development coordinates and template-exposure audit before prediction."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT,write_tsv
from catalytic_confirmation import sha
from extract_structure_sites import extract


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    policy_path=ROOT/'config/round19_structural_challenge.json';policy=json.loads(policy_path.read_text())
    raw=ROOT/'references/raw/round19_structural_challenge';downloads=json.loads((raw/'download_manifest.json').read_text())
    if downloads['protocol_sha256']!=sha(policy_path):raise ValueError('Challenge source protocol changed')
    files={r['pdb_id']:r for r in downloads['files']};sources=[];proteins=[];source_rows=[]
    for spec in policy['records']:
        if spec['split']!='development_challenge':continue
        pdb=spec['pdb_id'];path=raw/(pdb+'.cif')
        if sha(path)!=files[pdb]['sha256']:raise ValueError('Challenge raw structure changed')
        source=dict(extract(path),source='https://www.rcsb.org/structure/'+pdb,split='development_challenge',sha256=sha(path));sources.append(source)
        metadata={p['rcsb_polymer_entity_container_identifiers']['entity_id']:p for p in spec['metadata']['polymer_entities']}
        entities={l['entity_id'] for s in source['sites'] for l in s['ligand_residues']}
        for entity in sorted(entities):
            sequence=source['sequences'][entity];meta=metadata[entity]
            if sequence!=meta['entity_poly']['pdbx_seq_one_letter_code_can'].replace('\n','').replace(' ',''):raise ValueError('RCSB metadata/canonical polymer mismatch')
            pid=pdb+'__entity_'+entity;proteins.append(SeqRecord(Seq(sequence),id=pid,description='frozen_experimental_structure_query'))
            source_rows.append(dict(protein_id=pid,pdb_id=pdb,entity=entity,sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                taxonomy=meta['rcsb_entity_source_organism'],source=source['source'],source_sha256=sha(path)))
    (out/'sources.json').write_text(json.dumps(sources,indent=2)+'\n');SeqIO.write(proteins,out/'proteins.faa','fasta')
    with (out/'template_exposure.log').open('w') as log:
        subprocess.run(['blastp','-query',str(out/'proteins.faa'),'-db',str(ROOT/'references/processed/fes/seed'),'-evalue','1e-5','-max_target_seqs','1000','-max_hsps','1',
            '-outfmt','6 qseqid sseqid pident qlen slen qstart qend sstart send bitscore evalue','-out',str(out/'template_exposure.tsv')],stderr=log,check=True)
    exposure={p.id:dict(max_comparable_template_identity=0,homolog_overlap=False) for p in proteins}
    for line in (out/'template_exposure.tsv').read_text().splitlines():
        q,s,identity,qlen,slen,qs,qe,ss,se,score,evalue=line.split('\t');identity=float(identity);qc=(int(qe)-int(qs)+1)/int(qlen);sc=(int(se)-int(ss)+1)/int(slen)
        if min(qc,sc)>=.7:exposure[q]['max_comparable_template_identity']=max(exposure[q]['max_comparable_template_identity'],identity)
        if identity>=70 and (qc if int(qlen)<=int(slen) else sc)>=.7:exposure[q]['homolog_overlap']=True
    for row in source_rows:row.update(exposure[row['protein_id']])
    (out/'source_audit.json').write_text(json.dumps(source_rows,indent=2)+'\n')
    freeze=dict(status='development_sources_and_truth_frozen_before_prediction',protocol_sha256=sha(policy_path),
        predictor_files_before_first_prediction={n:sha(ROOT/n) for n in policy['frozen_predictor_files']},
        gold_extractor_sha256=sha(ROOT/'scripts/extract_structure_sites.py'),
        sources=len(sources),query_proteins=len(proteins),canonical_sites=sum(len(s['sites']) for s in sources),
        artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        limitations='Development sources are now exposed. Homology audit is against FeS templates only; no all-pipeline or public-model independence guarantee. Asymmetric-unit copies and same-source structures are not biological replicates.')
    (out/'freeze.json').write_text(json.dumps(freeze,indent=2)+'\n');print(json.dumps({k:freeze[k] for k in ['status','sources','query_proteins','canonical_sites']},indent=2))


if __name__=='__main__':main()
