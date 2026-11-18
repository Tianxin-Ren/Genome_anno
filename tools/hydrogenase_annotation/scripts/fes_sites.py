#!/usr/bin/env python3
"""Map reviewed experimental Fe-S ligand sets to query coordinates."""
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, BLAST_FIELDS, load_blast, residue_states, write_tsv


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build_model(root=ROOT):
    root=Path(root)
    source_path=root/'references/structural_sites/templates.json'
    policy_path=root/'config/fes_site_review.json'
    sources=json.loads(source_path.read_text());policy=json.loads(policy_path.read_text())
    for name,digest in policy.get('source_files',{}).items():
        if sha(root/name)!=digest:raise ValueError('FeS review evidence source changed: '+name)
    directory=root/'references/processed/fes';directory.mkdir(exist_ok=True)
    (directory/'manifest.json').write_text('{"status":"running"}\n')
    review={r['site_id']:r for r in policy['sites']}
    if len(review)!=len(policy['sites']):
        raise ValueError('Duplicate FeS site review')
    proteins={};sites=[]
    for source in sources:
        if source['split']!='template':
            raise ValueError('Reserved structure cannot enter FeS templates')
        if sha(root/'references/raw/structures'/f"{source['pdb_id']}.cif")!=source['sha256']:
            raise ValueError('FeS structural source changed')
        for site in source['sites']:
            row=review.pop(site['site_id'])
            for corroboration in row.get('coordination_corroboration',[]):
                from extract_structure_sites import extract
                path=root/corroboration['file']
                if sha(path)!=corroboration['sha256']:raise ValueError('FeS corroboration source changed')
                record=extract(path);corroborating=next(s for s in record['sites'] if s['site_id']==corroboration['site_id'])
                if (record['sequences'][corroboration['entity_id']]!=source['sequences'][row['entity_id']] or
                    {str(l['position']):l['residue'] for l in corroborating['ligand_residues']}!=row['expected_ligands']):
                    raise ValueError('FeS exact-sequence coordination corroboration changed')
            ligands={str(r['position']):r['residue'] for r in site['ligand_residues']}
            entity=row['entity_id'];sequence=source['sequences'][entity]
            if (ligands!=row['expected_ligands'] or site['cluster_type']!=row['cluster_type'] or
                {r['entity_id'] for r in site['ligand_residues']}!={entity} or
                hashlib.sha256(sequence.encode()).hexdigest()!=row['sequence_sha256']):
                raise ValueError('FeS reviewed sequence/ligand set changed: '+site['site_id'])
            if not 1<=row['domain_start']<=min(map(int,ligands))<=max(map(int,ligands))<=row['domain_end']<=len(sequence):
                raise ValueError('FeS site window does not contain all reviewed ligands')
            if any(sequence[int(pos)-1]!=residue for pos,residue in row.get('alternative_coordination_positions',{}).items()):
                raise ValueError('Alternative coordination reference position changed')
            pid=source['pdb_id']+'__entity_'+entity
            proteins[pid]=SeqRecord(Seq(sequence),id=pid,description='experimental_structure_entity')
            sites.append(dict(**row,reference_id=pid,pdb_id=source['pdb_id'],source_sha256=source['sha256'],
                reference_sequence=sequence))
    if review:
        raise ValueError('FeS review contains unused sites')
    fasta=directory/'templates.faa';SeqIO.write(list(proteins.values()),fasta,'fasta')
    (directory/'sites.json').write_text(json.dumps(sites,indent=2)+'\n')
    windows={f'site_{i:04}':s for i,s in enumerate(sites)}
    window_fasta=directory/'site_windows.faa'
    SeqIO.write([SeqRecord(Seq(s['reference_sequence'][s['domain_start']-1:s['domain_end']]),id=k,description='reviewed_site_window') for k,s in windows.items()],window_fasta,'fasta')
    (directory/'site_window_index.json').write_text(json.dumps({k:s['site_id'] for k,s in windows.items()},indent=2)+'\n')
    with (directory/'makeblastdb.log').open('w') as log:
        subprocess.run(['makeblastdb','-in',str(fasta),'-dbtype','prot','-out',str(directory/'seed')],stdout=log,stderr=log,check=True)
        subprocess.run(['makeblastdb','-in',str(window_fasta),'-dbtype','prot','-out',str(directory/'site_seed')],stdout=log,stderr=log,check=True)
    manifest=dict(status='complete',sources_sha256=sha(source_path),policy_sha256=sha(policy_path),
        review_source_files=policy.get('source_files',{}),
        prediction_script_sha256=sha(Path(__file__)),
        files={p.name:sha(p) for p in [fasta,directory/'sites.json',window_fasta,directory/'site_window_index.json',*sorted(directory.glob('seed.*')),*sorted(directory.glob('site_seed.*'))]},
        template_entities=len(proteins),sites=len(sites),reviewed_sites=sum(s['status']=='reviewed_primary_ligand_set' for s in sites),
        parameters=policy['parameters'],scope=policy['scope'])
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def site_evidence(hit, site, parameters):
    lo,hi=site['domain_start'],site['domain_end']
    sequence=site['reference_sequence']
    window,_=residue_states(hit,list(range(lo,hi+1)),{str(i):sequence[i-1] for i in range(lo,hi+1)})
    covered=[r for r in window if r['query_position'] is not None]
    coverage=len(covered)/len(window)
    identity=100*sum(r['status']=='present' for r in covered)/max(1,len(covered))
    states,_=residue_states(hit,[int(p) for p in site['expected_ligands']],site['expected_ligands'])
    alternatives,_=residue_states(hit,[int(p) for p in site.get('alternative_coordination_positions',{})],site.get('alternative_coordination_positions',{}))
    extra=[];rpos=hit['sstart']-1;qpos=hit['qstart']-1
    for q,r in zip(hit['qseq'],hit['sseq']):
        if r!='-':rpos+=1
        if q!='-':qpos+=1
        if lo<=rpos<=hi and q=='C' and r not in ['C','U']:
            extra.append(qpos)
    reasons=[]
    if site['status']!='reviewed_primary_ligand_set':reasons.append('incomplete_source_ligands')
    if coverage<parameters['min_window_coverage']:reasons.append('site_window_not_covered')
    if identity<parameters['min_window_identity']:reasons.append('site_window_low_identity')
    for status in ['not_covered','substituted','ambiguous']:
        if any(r['status']==status for r in states):reasons.append('ligand_'+status)
    # Extra non-ligand Cys need not change an already complete four-ligand cluster.
    if parameters['extra_cysteine_guard'] and extra and site['cluster_type']=='[3Fe-4S]':reasons.append('additional_Cys_in_site_window')
    if any(a['status']=='substituted' and a['observed_residue'] in 'CHD' for a in alternatives):reasons.append('possible_alternative_ligand')
    if any(a['status']=='not_covered' for a in alternatives):reasons.append('alternative_coordination_position_not_covered')
    if any(a['status']=='ambiguous' for a in alternatives):reasons.append('alternative_coordination_position_ambiguous')
    status=('ambiguous' if any(r in reasons for r in ['additional_Cys_in_site_window','ligand_ambiguous','possible_alternative_ligand','alternative_coordination_position_ambiguous']) else
            'unknown' if reasons else 'template_supported')
    positions=[r['query_position'] for r in states if r['query_position'] is not None]
    return dict(protein_id=hit['qseqid'],site_position=site['site_position'],template_site_id=site['site_id'],
        structure_template=site['pdb_id'],template_entity_id=site['entity_id'],reference_id=site['reference_id'],
        template_cluster_type=site['cluster_type'],cluster_type_predicted=site['cluster_type'] if status=='template_supported' else status,
        status=status,ligand_residues=states,alternative_coordination_states=alternatives,query_ligand_positions=positions,
        motif_match='complete' if all(r['status']=='present' for r in states) else 'incomplete_or_changed',
        domain_alignment='supported' if coverage>=parameters['min_window_coverage'] and identity>=parameters['min_window_identity'] else 'insufficient',
        reference_window_start=lo,reference_window_end=hi,window_coverage=coverage,window_identity=identity,
        query_alignment_start=hit['qstart'],query_alignment_end=hit['qend'],
        bitscore=hit['bitscore'],evalue=hit['evalue'],additional_Cys_positions=extra,
        reasons=reasons,source=site['source'],source_sha256=site['source_sha256'],
        coordination_corroboration=site.get('coordination_corroboration',[]),source_geometry_warnings=site.get('source_geometry_warnings',[]),
        topology_evidence='structural_homology_with_mapped_ligands',cofactor_occupancy='unknown')


def overlap(a,b):
    return max(0,min(a['qend'],b['qend'])-max(a['qstart'],b['qstart'])+1)/min(a['qend']-a['qstart']+1,b['qend']-b['qstart']+1)


def full_reference_coordinates(hit,site):
    lo,hi=site['domain_start'],site['domain_end']
    if hit['slen']!=hi-lo+1 or not 1<=hit['sstart']<=hit['send']<=hit['slen']:
        raise ValueError('FeS window alignment/source coordinates changed')
    return dict(hit,sseqid=site['reference_id'],sstart=hit['sstart']+lo-1,send=hit['send']+lo-1,slen=len(site['reference_sequence']))


def rescue_with_site_windows(rows,evidence):
    by_query=defaultdict(list)
    for row in evidence:by_query[row['protein_id']].append(row)
    for i,row in enumerate(rows):
        row['alignment_source']='full_reference_search'
        row.setdefault('full_alignment_leader_bitscore',row['bitscore'])
        allowed={'site_window_not_covered','site_window_low_identity','ligand_substituted'}
        if row['status']!='unknown' or not row['reasons'] or not set(row['reasons'])<=allowed:continue
        coordinates={(r['query_position'],r['observed_residue']) for r in row['ligand_residues']}
        matched=[r for r in by_query[row['protein_id']] if r['status']=='template_supported' and
            {(s['query_position'],s['observed_residue']) for s in r['ligand_residues']}==coordinates]
        types={r['template_cluster_type'] for r in matched}
        if len(types)>1:
            row.update(status='ambiguous',cluster_type_predicted='ambiguous',reasons=[*row['reasons'],'conflicting_site_window_templates'])
        elif types=={row['template_cluster_type']}:
            best=dict(max(matched,key=lambda r:(r['window_coverage'],r['window_identity'],r['bitscore'],r['template_site_id'])))
            best.update(fes_site_id=row['fes_site_id'],full_alignment_leader_template_site=row['full_alignment_leader_template_site'],
                full_alignment_leader_bitscore=row['full_alignment_leader_bitscore'],site_selection_reason='reviewed_window_supports_identical_ligands_and_type',
                competing_template_sites=sorted(set(row['competing_template_sites']+[r['template_site_id'] for r in matched])),
                alignment_source='reviewed_site_window_search',confidence=row['confidence'])
            rows[i]=best
    return rows


def map_hits(hits, sites, parameters):
    by_reference=defaultdict(list);by_query=defaultdict(list)
    for site in sites:by_reference[site['reference_id']].append(site)
    for hit in hits:
        if hit['evalue']<=parameters['max_evalue']:by_query[hit['qseqid']].append(hit)
    results=[];all_evidence=[]
    for qid,query_hits in sorted(by_query.items()):
        selected=[]
        for hit in sorted(query_hits,key=lambda h:(-h['bitscore'],h['sseqid'],h['qstart'])):
            competing=[h for h in selected if overlap(hit,h)>=.5]
            if competing and hit['bitscore']<max(h['bitscore'] for h in competing)*parameters['competitive_bitscore_ratio']:
                continue
            selected.append(hit)
        candidates=[site_evidence(h,s,parameters) for h in selected for s in by_reference[h['sseqid']]]
        all_evidence.extend(candidates)
        groups=[]
        for row in sorted(candidates,key=lambda r:(-r['bitscore'],r['template_site_id'],r['query_alignment_start'])):
            positions=set(row['query_ligand_positions'])
            same=[]
            for group in groups:
                for member in group:
                    other=set(member['query_ligand_positions'])
                    # Compare ligand coordinates, not overlap of interleaved FeS domains.
                    shared=bool(positions and other) and len(positions&other)>=.5*min(len(positions),len(other))
                    absent=(not positions and not other and row['site_position']==member['site_position'] and
                        min(row['query_alignment_end'],member['query_alignment_end'])>=max(row['query_alignment_start'],member['query_alignment_start']))
                    if shared or absent:
                        same.append(group);break
            if same:
                first=same[0];first.append(row)
                for other in same[1:]:first.extend(other);groups.remove(other)
            else:groups.append([row])
        for i,group in enumerate(sorted(groups,key=lambda g:(min(g[0]['query_ligand_positions'] or [g[0]['query_alignment_start']]),g[0]['site_position'])),1):
            if not any(g['query_ligand_positions'] for g in group):
                continue
            best=dict(group[0]);leader=best['template_site_id'];selection='highest_full_alignment_score'
            # A stronger whole-protein alignment must not hide a valid local map of the identical ligand set.
            if best['status']=='unknown' and best['reasons'] and set(best['reasons'])<={'site_window_not_covered','site_window_low_identity'}:
                coordinates={(r['query_position'],r['observed_residue']) for r in best['ligand_residues']}
                qualified=[r for r in group if r['status']=='template_supported' and r['template_cluster_type']==best['template_cluster_type'] and
                    {(s['query_position'],s['observed_residue']) for s in r['ligand_residues']}==coordinates]
                if qualified:
                    best=dict(max(qualified,key=lambda r:(r['window_coverage'],r['window_identity'],r['bitscore'],r['template_site_id'])))
                    selection='qualified_site_window_with_identical_ligand_coordinates'
            types={g['template_cluster_type'] for g in group if g['status']=='template_supported'}
            if len(types)>1:
                best.update(status='ambiguous',cluster_type_predicted='ambiguous',reasons=[*best['reasons'],'conflicting_competitive_templates'])
            best.update(fes_site_id=f'{qid}::FeS{i:03}',competing_template_sites=[g['template_site_id'] for g in group],
                full_alignment_leader_template_site=leader,full_alignment_leader_bitscore=group[0]['bitscore'],site_selection_reason=selection,
                confidence='provisional_structural_homology_not_calibrated')
            results.append(best)
    return results,all_evidence


SITE_FIELDS=['fes_site_id','protein_id','complex_anchors','association_scope','site_position','cluster_type_predicted','status',
    'ligand_residues','alternative_coordination_states','motif_match','domain_alignment','structure_template','template_entity_id','template_site_id',
    'reference_window_start','reference_window_end','query_alignment_start','query_alignment_end','window_coverage','window_identity',
    'bitscore','evalue','additional_Cys_positions','reasons','competing_template_sites','topology_evidence','confidence',
    'cofactor_occupancy','source','source_sha256','coordination_corroboration','source_geometry_warnings',
    'full_alignment_leader_template_site','site_selection_reason','alignment_source','full_alignment_leader_bitscore']


def annotate_fes(fasta, output, modules=(), threads=1, root=ROOT, anchors=()):
    root=Path(root);output=Path(output);directory=root/'references/processed/fes'
    manifest=json.loads((directory/'manifest.json').read_text())
    if manifest.get('status')!='complete':
        raise ValueError('FeS model is incomplete')
    raw=output/'fes_blast.tsv'
    with (output/'fes_blast.log').open('w') as log:
        subprocess.run(['blastp','-query',str(fasta),'-db',str(directory/'seed'),'-evalue',str(manifest['parameters']['max_evalue']),
            '-max_target_seqs','1000','-max_hsps','10','-num_threads',str(threads),'-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(raw)],stderr=log,check=True)
    thresholds=dict(blast_evalue=manifest['parameters']['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0)
    sites=json.loads((directory/'sites.json').read_text());by_site={s['site_id']:s for s in sites}
    rows,evidence=map_hits(load_blast(raw,thresholds),sites,manifest['parameters'])
    for row in evidence:row['alignment_source']='full_reference_search'
    window_index=json.loads((directory/'site_window_index.json').read_text());window_raw=output/'fes_site_window_blast.tsv'
    with (output/'fes_site_window_blast.log').open('w') as log:
        subprocess.run(['blastp','-query',str(fasta),'-db',str(directory/'site_seed'),'-evalue',str(manifest['parameters']['max_evalue']),
            '-max_target_seqs','1000','-max_hsps','10','-num_threads',str(threads),'-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(window_raw)],stderr=log,check=True)
    local=[]
    for hit in load_blast(window_raw,thresholds):
        site=by_site[window_index[hit['sseqid']]]
        local.append(dict(site_evidence(full_reference_coordinates(hit,site),site,manifest['parameters']),alignment_source='reviewed_site_window_search'))
    rows=rescue_with_site_windows(rows,local);evidence.extend(local)
    links=defaultdict(set)
    for member in modules:links[member['protein_id']].add(member['anchor'])
    for qid in anchors:links[qid].add(qid)
    for row in rows:
        row['complex_anchors']=sorted(links[row['protein_id']])
        row['association_scope']=('candidate_catalytic_anchor' if row['protein_id'] in anchors else
            'candidate_complex_member' if row['complex_anchors'] else 'unlinked_FeS_homolog_not_hydrogenase_identity')
    def serialized(row):return {k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in row.items()}
    write_tsv(output/'fes_sites.tsv',[serialized(r) for r in rows],SITE_FIELDS)
    write_tsv(output/'fes_template_evidence.tsv',[serialized(r) for r in evidence],list(evidence[0]) if evidence else ['protein_id','template_site_id','status'])
    ligands=[dict(fes_site_id=r['fes_site_id'],protein_id=r['protein_id'],template_site_id=r['template_site_id'],**ligand) for r in rows for ligand in r['ligand_residues']]
    write_tsv(output/'fes_ligands.tsv',ligands,['fes_site_id','protein_id','template_site_id','reference_position','expected_residue','query_position','observed_residue','status','coverage_reason'])
    return rows,manifest
