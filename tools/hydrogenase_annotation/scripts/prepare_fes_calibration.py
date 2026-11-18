#!/usr/bin/env python3
"""Freeze two declared calibration structures and audit future sources before prediction."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT
from catalytic_confirmation import sha
from extract_structure_sites import extract


def homolog_overlap(identity, qlen, slen, qstart, qend, sstart, send, policy):
    covered=(qend-qstart+1)/qlen if qlen<=slen else (send-sstart+1)/slen
    return identity>=policy['identity_percent'] and covered>=policy['shorter_coverage']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    protocol=ROOT/'config/round20_fes_validation.json';policy=json.loads(protocol.read_text())
    raw=ROOT/'references/raw/round20_fes_calibration';download=json.loads((raw/'download_manifest.json').read_text())
    if download['protocol_sha256']!=sha(protocol):raise ValueError('Calibration source selection changed')
    downloaded={r['pdb_id']:r for r in download['files']};calibration=[]
    for spec in policy['records']:
        if spec['split']!='calibration':continue
        pdb=spec['pdb_id'];path=raw/(pdb+'.cif')
        if sha(path)!=downloaded[pdb]['sha256']:raise ValueError('Calibration coordinates changed')
        calibration.append(dict(extract(path),source='https://www.rcsb.org/structure/'+pdb,split='calibration',sha256=sha(path)))
    sources=(json.loads((ROOT/'references/structural_sites/templates.json').read_text())+
        json.loads((ROOT/'references/structural_sites/reserved_validation.json').read_text())+
        json.loads((ROOT/'references/validation/round20_structural_challenge_development_v3/sources.json').read_text()))
    metadata={}
    metadata_files=sorted((ROOT/'references/raw/round19_structural_search').glob('entries_*.json'))
    for path in metadata_files:
        for entry in json.loads(path.read_text())['data']['entries']:metadata[entry['rcsb_id']]=entry
    metadata.update({r['pdb_id']:r['metadata'] for r in policy['records']})
    def genera(pdb,entity):
        match=next(e for e in metadata[pdb]['polymer_entities'] if e['rcsb_polymer_entity_container_identifiers']['entity_id']==entity)
        return sorted({r['ncbi_scientific_name'].split()[0] for r in match['rcsb_entity_source_organism']})
    subjects=[];subject_info={};n=0
    for source in sources+calibration:
        pdb=source['pdb_id'];entities={r['entity_id'] for s in source['sites'] for r in s['ligand_residues']}
        for entity in sorted(entities):
            sequence=source['sequences'][entity];windows={(1,len(sequence))}
            # Ligand-spanning domains also guard against fusion proteins hiding a close core.
            for site in source['sites']:
                positions=[r['position'] for r in site['ligand_residues'] if r['entity_id']==entity]
                if positions:windows.add((max(1,min(positions)-15),min(len(sequence),max(positions)+15)))
            for lo,hi in sorted(windows):
                n+=1;sid='exposed_'+str(n)
                subjects.append(SeqRecord(Seq(sequence[lo-1:hi]),id=sid,description=''))
                subject_info[sid]=dict(pdb_id=pdb,entity_id=entity,start=lo,end=hi,genera=genera(pdb,entity),split='calibration' if source in calibration else 'development_or_template')
    queries=[];rows={};calibration_proteins=[]
    for spec in policy['records']:
        pdb=spec['pdb_id']
        for entity in spec['metadata']['polymer_entities']:
            eid=entity['rcsb_polymer_entity_container_identifiers']['entity_id'];pid=pdb+'__entity_'+eid
            sequence=''.join(entity['entity_poly']['pdbx_seq_one_letter_code_can'].split())
            queries.append(SeqRecord(Seq(sequence),id=pid,description='metadata_before_prediction'))
            row=dict(protein_id=pid,pdb_id=pdb,entity_id=eid,declared_split=spec['split'],genera=genera(pdb,eid),homology_overlaps=[],genus_overlaps=[],reasons=[])
            if spec['split']=='technical_same_source_control':row['reasons'].append('declared_same_source_technical_control')
            if pdb=='1HFE':row['reasons'].append('historical_publication_vs_current_metadata_taxonomy_requires_review')
            rows[pid]=row
            if spec['split']=='calibration':
                source=next(s for s in calibration if s['pdb_id']==pdb)
                if sequence!=source['sequences'][eid]:raise ValueError('Calibration metadata/coordinate sequence mismatch')
                if any(r['entity_id']==eid for s in source['sites'] for r in s['ligand_residues']):calibration_proteins.append(queries[-1])
    SeqIO.write(subjects,out/'exposed_full_and_ligand_domains.faa','fasta');SeqIO.write(queries,out/'metadata_queries.faa','fasta')
    with (out/'exposure.log').open('w') as log:
        subprocess.run(['blastp','-query',str(out/'metadata_queries.faa'),'-subject',str(out/'exposed_full_and_ligand_domains.faa'),'-evalue','1e-5',
            '-max_target_seqs','10000','-max_hsps','1','-outfmt','6 qseqid sseqid pident qlen slen qstart qend sstart send bitscore evalue','-out',str(out/'exposure.tsv')],stderr=log,check=True)
    for line in (out/'exposure.tsv').read_text().splitlines():
        q,s,identity,qlen,slen,qs,qe,ss,se,score,evalue=line.split('\t');row=rows[q];subject=subject_info[s]
        if subject['pdb_id']==row['pdb_id']:continue
        if row['declared_split']=='calibration' and subject['split']=='calibration':continue
        if homolog_overlap(float(identity),*[int(x) for x in [qlen,slen,qs,qe,ss,se]],policy['independence']):
            row['homology_overlaps'].append(dict(subject_id=s,identity=float(identity),subject=subject))
    for row in rows.values():
        comparisons=[s for s in subject_info.values() if s['pdb_id']!=row['pdb_id'] and not(row['declared_split']=='calibration' and s['split']=='calibration')]
        row['genus_overlaps']=sorted({g for s in comparisons for g in s['genera'] if g in row['genera']})
        if row['genus_overlaps']:row['reasons'].append('exposed_genus_overlap')
        if row['homology_overlaps']:row['reasons'].append('exposed_full_or_ligand_domain_homology_overlap')
        row['independent_eligible']=not row['reasons']
    SeqIO.write(calibration_proteins,out/'proteins.faa','fasta')
    (out/'sources.json').write_text(json.dumps(calibration,indent=2)+'\n')
    (out/'source_audit.json').write_text(json.dumps(list(rows.values()),indent=2)+'\n')
    (out/'exposed_subjects.json').write_text(json.dumps(subject_info,indent=2)+'\n')
    freeze=dict(status='calibration_sources_and_truth_frozen_before_prediction',protocol_sha256=sha(protocol),
        gold_extractor_sha256=sha(ROOT/'scripts/extract_structure_sites.py'),
        metadata_files={str(p.relative_to(ROOT)):sha(p) for p in metadata_files},
        source_files={str((raw/(r['pdb_id']+'.cif')).relative_to(ROOT)):r['sha256'] for r in download['files']},
        sources=len(calibration),query_proteins=len(calibration_proteins),canonical_sites=sum(len(s['sites']) for s in calibration),
        future_coordinates_downloaded=False,artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        limitations=['This audit establishes operational FeS-stage separation only, not all-pipeline independence.',
            'Genus names follow cached RCSB metadata, not an independently reconstructed taxonomy tree; 1HFE is excluded pending source review.',
            'Calibration outcomes may inform development; final structure coordinates and predictions remain unopened.',
            'Domain homology checks use every exposed ligand-spanning window, including calibration sources; not evidence of absence of evolutionary relatedness.'])
    (out/'freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
    print(json.dumps(dict(status=freeze['status'],sources=len(calibration),proteins=len(calibration_proteins),
        eligible_metadata_proteins=[r['protein_id'] for r in rows.values() if r['independent_eligible']]),indent=2))


if __name__=='__main__':main()
