"""Report a literature ligand hypothesis without modifying any annotation decision."""
import json
import subprocess
from math import dist
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, BLAST_FIELDS, load_blast, write_tsv
from active_centers import structure_atoms, map_hits
from catalytic_confirmation import sha
from extract_structure_sites import primary_contact

FIELDS=('protein_id candidate_model_id source_status pattern_status context_alignment_status template_span_status '
        'query_partial query_ligand_positions ligand_residues confidence query_activity cofactor_occupancy '
        'primary_source source_geometry limitations').split()


def load_review(root=ROOT):
    root=Path(root);path=root/'config/noncanonical_center_evidence.json';policy=json.loads(path.read_text())
    ref_fasta=root/'references/processed/noncanonical_center_literature_model.faa'
    reference=SeqIO.read(ref_fasta,'fasta')
    geometry=json.loads((root/'references/processed/noncanonical_center_geometry.json').read_text())
    model=dict(reference_id=policy['model_id'],sequence=str(reference.seq),family='FeFe',
               window_start=policy['window_start'],window_end=policy['window_end'],
               source_sha256=sha(path),evidence_url=policy['primary_source'],author_positions={},
               sites=[dict(reference_position=int(p),expected_residue=r,residue_role='direct_metal_ligand',
                           requirement='required_for_template_coordination',reference_residue=r) for p,r in policy['expected_ligands'].items()])
    parameters=json.loads((root/policy['parameters_source']).read_text())['parameters']
    return policy,model,parameters,geometry


def candidate_rows(hits,policy,model,parameters,geometry,partial_ids=(),query_sequences=None):
    selected,_=map_hits(hits,[model],parameters,query_sequences=query_sequences);rows=[]
    for row in selected:
        supported=(row['active_center_status']=='template_ligand_set_supported' and
                   row['template_span_alignment_status']=='supported' and row['protein_id'] not in partial_ids)
        rows.append(dict(protein_id=row['protein_id'],candidate_model_id=policy['model_id'],
            source_status='literature_reported_serine_ligand_geometry_not_qualified',
            pattern_status='literature_ligand_pattern_candidate' if supported else 'pattern_not_supported',
            context_alignment_status=row['context_alignment_status'],template_span_status=row['template_span_alignment_status'],
            query_partial=row['protein_id'] in partial_ids,
            query_ligand_positions=[s['query_position'] for s in row['sites']],
            ligand_residues=[{k:v for k,v in s.items() if k!='model_status'} for s in row['sites']],
            confidence='uncalibrated_literature_hypothesis_only',query_activity='not_measured',cofactor_occupancy='unknown',
            primary_source=policy['primary_source'],source_geometry=geometry,
            limitations='Alternative reference residue mapping only. Source geometry does not qualify under the fixed screen. '
                'No query Ser coordination, H-cluster assembly, catalytic identity, complete complex, reaction or physiology is established. '
                'The canonical active-center and all annotation decisions remain separately recorded.'))
    return rows


def annotate_evidence(fasta,out,eligible_ids,partial_ids=(),threads=2,root=ROOT):
    """The only annotation input is a set of IDs; there is no mutable decision object."""
    root=Path(root);out=Path(out)
    if out.name!='noncanonical_centers':raise ValueError('Expected dedicated auxiliary evidence directory')
    out.mkdir(parents=True,exist_ok=False);policy,model,parameters,geometry=load_review(root)
    if any(g['complete_primary_ligand_set'] for g in geometry):raise ValueError('Geometry now needs a fresh source review; do not auto-promote')
    records=[r for r in SeqIO.parse(fasta,'fasta') if r.id in eligible_ids]
    query=out/'queries.faa';SeqIO.write(records,query,'fasta');raw=out/'blast.tsv';raw.write_text('')
    if records:
        with (out/'blast.log').open('w') as log:
            subprocess.run(['blastp','-query',str(query),'-subject',str(root/'references/processed/noncanonical_center_literature_model.faa'),
                '-evalue',str(parameters['max_evalue']),'-max_hsps','20','-num_threads',str(threads),
                '-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(raw)],stderr=log,check=True)
    hits=load_blast(raw,dict(blast_evalue=parameters['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
    rows=candidate_rows(hits,policy,model,parameters,geometry,partial_ids,query_sequences={r.id:str(r.seq) for r in records})
    write_tsv(out/'alternatives.tsv',[{k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows],FIELDS)
    manifest=dict(status='auxiliary_literature_evidence_complete',policy_sha256=sha(root/'config/noncanonical_center_evidence.json'),
        script_sha256=sha(Path(__file__)),query_sha256=sha(query),raw_alignment_sha256=sha(raw),output_sha256=sha(out/'alternatives.tsv'),
        query_proteins=len(records),rows=len(rows),pattern_candidates=sum(r['pattern_status']=='literature_ligand_pattern_candidate' for r in rows),
        source_geometry=geometry,changes_annotation_decisions=False,independent_validation=False,scope=policy['scope'])
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest
