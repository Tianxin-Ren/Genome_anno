#!/usr/bin/env python3
"""Map complete reviewed protein ligand sets; sequence evidence is not enzyme activity."""
import hashlib
import json
import subprocess
from collections import defaultdict
from math import dist
from pathlib import Path

from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1
from annotate import ROOT, BLAST_FIELDS, load_blast, residue_states, repeat_mapping_ambiguity, write_tsv


SITE_FIELDS=('protein_id reference_id family reference_position expected_residue query_position observed_residue status '
    'coverage_reason evidence_url site_model residue_role requirement reference_residue source_sha256 '
    'template_author_position model_status window_coverage window_identity alignment_components cofactor_occupancy '
    'site_context_coverage site_context_identity site_context_start site_context_end catalytic_anchor source_center_id center_component_id').split()
SUMMARY_FIELDS=('protein_id family reference_id active_center_status direct_ligand_status reaction_signature_status '
    'window_coverage window_identity competing_templates confidence cofactor_occupancy '
    'context_alignment_status min_required_context_coverage min_required_context_identity site_contexts '
    'catalytic_domain_integrity template_span_alignment_status center_organization center_component_proteins').split()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def structure_atoms(path):
    """Read canonical polymer positions and positive-occupancy first-model atoms."""
    data=MMCIF2Dict(str(path))
    sequences=defaultdict(dict)
    for e,p,r in zip(data['_entity_poly_seq.entity_id'],data['_entity_poly_seq.num'],data['_entity_poly_seq.mon_id']):
        letter=seq1(r,custom_map={'SEC':'U','MSE':'M'})
        if int(p) in sequences[e] and sequences[e][int(p)]!=letter:
            raise ValueError('Microheterogeneous polymer needs explicit review')
        sequences[e][int(p)]=letter
    for positions in sequences.values():
        if sorted(positions)!=list(range(1,len(positions)+1)):raise ValueError('Non-contiguous polymer sequence')
    sequences={e:''.join(p[n] for n in sorted(p)) for e,p in sequences.items()}
    atoms=[]
    for i in range(len(data['_atom_site.id'])):
        if data['_atom_site.pdbx_PDB_model_num'][i]!='1':continue
        if float(data['_atom_site.occupancy'][i])<=0:continue
        atoms.append(dict(entity=data['_atom_site.label_entity_id'][i],chain=data['_atom_site.label_asym_id'][i],
            component=data['_atom_site.label_comp_id'][i],position=data['_atom_site.label_seq_id'][i],
            author_position=data['_atom_site.auth_seq_id'][i],element=data['_atom_site.type_symbol'][i].upper(),
            atom=data['_atom_site.label_atom_id'][i],xyz=tuple(float(data['_atom_site.Cartn_'+axis][i]) for axis in ['x','y','z'])))
    return sequences,atoms


def nife_metals(atoms):
    metals=[a for a in atoms if a['component'] in ['NI','3NI'] and a['element'] in ['FE','NI']]
    nickel=[a for a in metals if a['element']=='NI']
    # FCO is the carbonyl/dicyano iron; only admit it when paired with a nearby nickel.
    return metals+[a for a in atoms if a['component']=='FCO' and a['element']=='FE' and any(dist(a['xyz'],n['xyz'])<=3.5 for n in nickel)]


def structural_record(path, entity, family):
    """Read canonical polymer positions, retaining SEC as U and author offsets separately."""
    sequences,atoms=structure_atoms(path);entity=str(entity);sequence=sequences[entity]
    chain=min(a['chain'] for a in atoms if a['entity']==entity)
    protein=[a for a in atoms if a['entity']==entity and a['chain']==chain and a['position'].isdigit()]
    donors=[a for a in protein if a['component'] in ['CYS','SEC'] and a['element'] in ['S','SE']]
    allowed={'NiFe':{'NI','3NI'},'FeFe':{'HC1','FE2'},'Fe':{'FE2','FE9'}}[family]
    metals=nife_metals(atoms) if family=='NiFe' else [a for a in atoms if a['component'] in allowed and a['element'] in ['FE','NI']]
    contacts=[dict(position=int(p['position']),residue=sequence[int(p['position'])-1],
        author_position=p['author_position'],metal_component=m['component'],metal_chain=m['chain'],
        metal_author_position=m['author_position'],distance_angstrom=dist(p['xyz'],m['xyz']))
        for p in donors for m in metals if dist(p['xyz'],m['xyz'])<=2.8]
    ligands={c['position']:c['residue'] for c in contacts}
    if family=='FeFe' and len(ligands)==1:
        # H-cluster's bridging Cys identifies its cubane, excluding accessory cubanes.
        bridges=[p for p in donors if int(p['position']) in ligands]
        cubanes={(m['chain'],m['author_position']) for m in atoms if m['component']=='SF4' and m['element']=='FE'
            and any(dist(p['xyz'],m['xyz'])<=2.8 for p in bridges)}
        for p in donors:
            for m in atoms:
                if m['component']!='SF4' or m['element']!='FE' or (m['chain'],m['author_position']) not in cubanes:continue
                distance=dist(p['xyz'],m['xyz'])
                if distance<=2.8:
                    pos=int(p['position']);ligands[pos]=sequence[pos-1]
                    contacts.append(dict(position=pos,residue=sequence[pos-1],author_position=p['author_position'],
                        metal_component='SF4',metal_chain=m['chain'],metal_author_position=m['author_position'],distance_angstrom=distance))
    return dict(sequence=sequence,entity_id=entity,canonical_chain=chain,family=family,
        ligands={str(p):r for p,r in sorted(ligands.items())},contacts=contacts,
        author_positions={str(int(a['position'])):a['author_position'] for a in protein},
        source_sha256=sha(path),coordinate_evidence='experimental_ASU_first_model_S_Se_distance_le_2.8A_not_query_occupancy')


def split_structural_record(path, entities, nickel_chain, nickel_author_position):
    sequences,atoms=structure_atoms(path)
    nickel=[a for a in nife_metals(atoms) if a['element']=='NI' and a['chain']==nickel_chain and a['author_position']==str(nickel_author_position)]
    if len(nickel)!=1:raise ValueError('Reviewed nickel center is missing or nonunique')
    iron=[a for a in nife_metals(atoms) if a['element']=='FE' and dist(a['xyz'],nickel[0]['xyz'])<=3.5]
    if len(iron)!=1:raise ValueError('Reviewed NiFe iron pairing is missing or ambiguous')
    metals=nickel+iron;components=[]
    for entity in map(str,entities):
        donors=[a for a in atoms if a['entity']==entity and a['component'] in ['CYS','SEC'] and a['element'] in ['S','SE'] and
            any(dist(a['xyz'],m['xyz'])<=2.8 for m in metals)]
        chains={a['chain'] for a in donors}
        if len(chains)!=1:raise ValueError('Split component chain is missing or ambiguous at the selected center')
        chain=next(iter(chains));protein=[a for a in atoms if a['entity']==entity and a['chain']==chain and a['position'].isdigit()]
        contacts=[dict(position=int(p['position']),residue=sequences[entity][int(p['position'])-1],author_position=p['author_position'],
            metal_component=m['component'],metal_chain=m['chain'],metal_author_position=m['author_position'],distance_angstrom=dist(p['xyz'],m['xyz']))
            for p in donors for m in metals if dist(p['xyz'],m['xyz'])<=2.8]
        components.append(dict(entity_id=entity,canonical_chain=chain,sequence=sequences[entity],ligands={str(c['position']):c['residue'] for c in contacts},
            contacts=contacts,author_positions={str(int(a['position'])):a['author_position'] for a in protein}))
    return dict(components=components,source_sha256=sha(path),nickel=nickel[0],iron=iron[0],
        coordinate_evidence='same_reviewed_NiFe_pair; S/Se_distance_le_2.8A; separate_polypeptide_coordinates; not_query_occupancy')


def build_model(root=ROOT):
    root=Path(root);policy_path=root/'config/active_center_review.json';policy=json.loads(policy_path.read_text())
    parameters=policy['parameters'];flank=parameters['site_context_flank']
    if (not isinstance(flank,int) or isinstance(flank,bool) or flank<1 or
        not 0<parameters['min_window_coverage']<=1 or not 0<parameters['min_window_identity']<=100 or
        not 0<parameters['competitive_bitscore_ratio']<=1 or parameters['max_evalue']<=0):
        raise ValueError('Invalid active-center alignment policy')
    plan_path=root/'config/active_center_sources.json';plan=json.loads(plan_path.read_text())
    raw=root/'references/raw/catalytic_structures';downloads=json.loads((raw/'download_manifest.json').read_text())
    if downloads['source_plan_sha256']!=sha(plan_path):raise ValueError('Catalytic source plan changed')
    source_index={s['pdb_id']:s for s in plan['sources']};pins={s['pdb_id']:s for s in downloads['sources']}
    directory=root/'references/processed/active_centers';directory.mkdir(parents=True,exist_ok=True)
    (directory/'manifest.json').write_text('{"status":"running"}\n')
    models=[];records=[];extracted=[]
    for review in policy['proteins']:
        pdb=review['pdb_id'];path=raw/(pdb+'.cif')
        if sha(path)!=pins[pdb]['sha256']:raise ValueError('Catalytic structural source changed: '+pdb)
        row=structural_record(path,review['entity_id'],review['family'])
        if row['ligands']!=review['expected_structural_ligands']:raise ValueError('Structural ligand review mismatch: '+pdb)
        row.update(pdb_id=pdb,split=source_index[pdb]['split'],evidence_url=review['evidence_url'])
        extracted.append(row)
        if row['split']!='template':continue
        expected_count=1 if row['family']=='Fe' else 4
        if len(row['ligands'])!=expected_count:raise ValueError('Incomplete catalytic template')
        model=dict(row,reference_id=pdb+'__entity_'+row['entity_id'],
            window_start=review['window_start'],window_end=review['window_end'])
        sites={}
        for p,r in row['ligands'].items():
            sites[int(p)]=dict(expected_residue='CU' if row['family']=='NiFe' else r,
                residue_role='direct_metal_ligand',requirement='required_for_template_coordination',reference_residue=r)
        for neighbor in review.get('neighbors',[]):
            p=neighbor['position'];r=row['sequence'][p-1]
            if r!=neighbor['reference_residue'] or p in sites:raise ValueError('Neighbor sequence review mismatch')
            sites[p]=dict(expected_residue=neighbor.get('allowed_residues',r),residue_role=neighbor['role'],
                requirement=neighbor['requirement'],reference_residue=r,evidence_url=neighbor['evidence_url'])
        if not 1<=model['window_start']<=min(sites)<=max(sites)<=model['window_end']<=len(row['sequence']):
            raise ValueError('Catalytic model window excludes reviewed sites')
        model['sites']=[dict(reference_position=p,**s) for p,s in sorted(sites.items())]
        models.append(model);records.append(SeqRecord(Seq(row['sequence']),id=model['reference_id'],description='reviewed_structural_catalytic_template'))
    if len({r.id for r in records})!=len(records):raise ValueError('Duplicate active-center template')
    fasta=directory/'templates.faa';SeqIO.write(records,fasta,'fasta')
    for name,rows in [('models',models),('structural_truth',extracted)]:
        (directory/(name+'.json')).write_text(json.dumps(rows,indent=2)+'\n')
    with (directory/'makeblastdb.log').open('w') as log:
        subprocess.run(['makeblastdb','-in',str(fasta),'-dbtype','prot','-out',str(directory/'seed')],stdout=log,stderr=log,check=True)
    variants=policy.get('exact_active_reference_variants',[])
    for variant in variants:
        source=root/variant['source_fasta']
        protein=next(r for r in SeqIO.parse(source,'fasta') if r.id==variant['reference_id'])
        if (hashlib.sha256(str(protein.seq).encode()).hexdigest()!=variant['sequence_sha256'] or
            str(protein.seq)[variant['query_position']-1]!=variant['query_residue']):
            raise ValueError('Reviewed active noncanonical reference changed')
    manifest=dict(status='complete',policy_sha256=sha(policy_path),source_plan_sha256=sha(plan_path),
        download_manifest_sha256=sha(raw/'download_manifest.json'),parameters=policy['parameters'],
        template_proteins=len(models),positions=sum(len(m['sites']) for m in models),exact_active_reference_variants=variants,
        files={p.name:sha(p) for p in [fasta,directory/'models.json',directory/'structural_truth.json',*sorted(directory.glob('seed.*'))]},
        scope=policy['scope'])
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def collinear_hits(hits):
    """One non-overlapping chain per template and query; no repeated-domain double count."""
    ordered=sorted(hits,key=lambda h:(h['qstart'],h['sstart'],h['qend'],h['send']));chains=[]
    for i,h in enumerate(ordered):
        candidates=[[h]]+[chains[j]+[h] for j,a in enumerate(ordered[:i])
            if a['qend']<h['qstart'] and a['send']<h['sstart'] and
            (h['qstart']-a['qend'])-(h['sstart']-a['send'])<=h['slen']]
        chains.append(max(candidates,key=lambda c:sum(x['bitscore'] for x in c)))
    return max(chains,key=lambda c:sum(x['bitscore'] for x in c)) if chains else []


def map_model(components, model, parameters, query_sequence=None):
    sequence=model['sequence'];lo,hi=model['window_start'],model['window_end'];mapped={};ambiguous_positions=set()
    for hit in components:
        states,piece=residue_states(hit,list(range(lo,hi+1)))
        ambiguous_positions.update(s['reference_position'] for s in states if s['coverage_reason']=='repeat_indel_alignment_ambiguity')
        if set(piece)&set(mapped):raise ValueError('Overlapping catalytic reference segments')
        mapped.update(piece)
    if query_sequence is not None:ambiguous_positions.update(repeat_mapping_ambiguity(mapped,sequence,query_sequence))
    covered={p:(q,r) for p,(q,r) in mapped.items() if q is not None}
    coverage=len(covered)/(hi-lo+1)
    identity=100*sum(sequence[p-1]==r for p,(q,r) in covered.items())/max(1,len(covered))
    # Site mapping needs its local context; whole-template completeness is a separate question.
    radius=parameters['site_context_flank']
    contexts=[]
    for site in model['sites']:
        if site['requirement'] not in ['required_for_template_coordination','canonical_reaction_signature']:continue
        p=site['reference_position'];start=max(lo,p-radius);end=min(hi,p+radius)
        local={n:covered[n] for n in range(start,end+1) if n in covered}
        contexts.append(dict(reference_position=p,start=start,end=end,coverage=len(local)/(end-start+1),
            identity=100*sum(sequence[n-1]==r for n,(q,r) in local.items())/max(1,len(local))))
    if not contexts:raise ValueError('Active-center model has no required sites')
    adequate=all(c['coverage']>=parameters['min_window_coverage'] and c['identity']>=parameters['min_window_identity'] for c in contexts)
    context_by_position={c['reference_position']:c for c in contexts}
    sites=[]
    coordinates=[dict(query_start=h['qstart'],query_end=h['qend'],reference_start=h['sstart'],reference_end=h['send'],
        bitscore=h['bitscore'],evalue=h['evalue']) for h in components]
    for site in model['sites']:
        p=site['reference_position'];q,r=mapped.get(p,(None,''));expected=site['expected_residue']
        state=('ambiguous' if p in ambiguous_positions else 'not_covered' if r in ['','-'] else 'ambiguous' if r in 'XBZJO' else 'present' if r in expected else 'substituted')
        if p in ambiguous_positions:q=None
        context=context_by_position.get(p)
        sites.append(dict(protein_id=components[0]['qseqid'],reference_id=model['reference_id'],family=model['family'],**{k:v for k,v in site.items() if k!='evidence_url'},
            query_position=q,observed_residue=r,status=state,coverage_reason='repeat_indel_alignment_ambiguity' if p in ambiguous_positions else 'alignment_deletion' if r=='-' else 'outside_alignment' if r=='' else 'aligned',
            evidence_url=site.get('evidence_url',model['evidence_url']),site_model='complete_template_protein_ligands_and_reviewed_neighbors',
            source_sha256=model['source_sha256'],template_author_position=model['author_positions'].get(str(p),''),
            window_coverage=coverage,window_identity=identity,alignment_components=coordinates,cofactor_occupancy='unknown',
            site_context_coverage=context['coverage'] if context else None,site_context_identity=context['identity'] if context else None,
            site_context_start=context['start'] if context else None,site_context_end=context['end'] if context else None))
    def combined(rows):
        for state in ['ambiguous','not_covered','substituted']:
            if any(r['status']==state for r in rows):return state
        return 'present' if rows else 'not_applicable'
    ligands=combined([s for s in sites if s['residue_role']=='direct_metal_ligand'])
    signature=combined([s for s in sites if s['requirement']=='canonical_reaction_signature'])
    status=('insufficient_alignment' if not adequate else
        'ligands_'+ligands if ligands!='present' else 'canonical_reaction_signature_changed' if signature not in ['present','not_applicable'] else
        'template_ligand_set_supported')
    for s in sites:s['model_status']=status
    summary=dict(protein_id=components[0]['qseqid'],family=model['family'],reference_id=model['reference_id'],
        active_center_status=status,direct_ligand_status=ligands,reaction_signature_status=signature,
        window_coverage=coverage,window_identity=identity,bitscore=sum(h['bitscore'] for h in components),sites=sites,
        context_alignment_status='supported' if adequate else 'insufficient',
        min_required_context_coverage=min(c['coverage'] for c in contexts),min_required_context_identity=min(c['identity'] for c in contexts),
        site_contexts=contexts,catalytic_domain_integrity='not_established_by_active_site_mapping',
        template_span_alignment_status='supported' if coverage>=parameters['min_window_coverage'] and identity>=parameters['min_window_identity'] else 'insufficient',
        confidence='provisional_structural_mapping_not_catalytic_activity',cofactor_occupancy='unknown')
    return summary


def map_hits(hits, models, parameters, query_sequences=None):
    by_pair=defaultdict(list);index={m['reference_id']:m for m in models}
    for h in hits:
        if h['evalue']<=parameters['max_evalue']:by_pair[(h['qseqid'],h['sseqid'])].append(h)
    by_query=defaultdict(list)
    for (qid,rid),rows in sorted(by_pair.items()):
        by_query[qid].append(map_model(collinear_hits(rows),index[rid],parameters,query_sequence=query_sequences.get(qid) if query_sequences is not None else None))
    selected=[];evidence=[]
    for qid,rows in sorted(by_query.items()):
        rows.sort(key=lambda r:(-r['bitscore'],r['reference_id']));evidence.extend(rows)
        competing=[r for r in rows if r['bitscore']>=rows[0]['bitscore']*parameters['competitive_bitscore_ratio']]
        adequate=[r for r in competing if r['context_alignment_status']=='supported']
        span_supported=[r for r in adequate if r['template_span_alignment_status']=='supported']
        best=dict((span_supported or adequate or competing)[0])
        signatures={(r['family'],tuple(s['query_position'] for s in r['sites'] if s['residue_role']=='direct_metal_ligand')) for r in (adequate or competing)}
        if len(signatures)>1:
            best['active_center_status']='ambiguous_competing_ligand_maps'
        best['competing_templates']=[r['reference_id'] for r in competing]
        selected.append(best)
    return selected,evidence


def apply_site_guard(annotations, selected):
    by_query={r['protein_id']:r for r in selected}
    for anchor in annotations:
        row=by_query.get(anchor['protein_id'])
        status=row['active_center_status'] if row else 'no_reviewed_template_mapping'
        anchor['active_center_status']=status
        anchor['active_center_template']=row['reference_id'] if row else ''
        anchor['active_center_span_status']=row['template_span_alignment_status'] if row else 'not_evaluable'
        anchor['catalytic_domain_integrity']='not_established_by_active_site_mapping'
        # A noncanonical ligand is a review boundary, not proof of an inactive enzyme.
        if status not in ['template_ligand_set_supported','experimentally_active_exact_reference_noncanonical_site','multi_polypeptide_ligand_set_supported','fused_split_template_ligand_set_supported']:
            anchor.update(catalytic_status='candidate',function_status='unassigned',function_confidence='unassigned',potential_function_ids='')
            if anchor['structure_confidence']=='high':anchor['structure_confidence']='medium'
            anchor['flags']=';'.join(filter(None,[anchor.get('flags',''),'complete_active_center_requires_review']))
        elif row['template_span_alignment_status']!='supported':
            # Local ligand mapping must not rescue a missing catalytic scaffold for function transfer.
            anchor.update(function_status='unassigned',function_confidence='unassigned',potential_function_ids='')
            if anchor['structure_confidence']=='high':anchor['structure_confidence']='medium'
            anchor['flags']=';'.join(filter(None,[anchor.get('flags',''),'catalytic_template_span_incomplete_or_divergent']))
        if row and status in ['multi_polypeptide_ligand_set_supported','fused_split_template_ligand_set_supported']:
            anchor.update(function_status='unassigned',function_confidence='unassigned',potential_function_ids='')
            if anchor['structure_confidence']=='high':anchor['structure_confidence']='medium'
            anchor['flags']=';'.join(filter(None,[anchor.get('flags',''),'split_template_center_composition_and_confirmation_not_calibrated']))


def write_active_center_tables(output, selected, evidence=None):
    def serialized(row):return {k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in row.items()}
    for row in selected:
        row.setdefault('center_organization','single_polypeptide');row.setdefault('center_component_proteins',[row['protein_id']])
        for site in row['sites']:
            site.setdefault('catalytic_anchor',row['protein_id']);site.setdefault('source_center_id',row['reference_id']);site.setdefault('center_component_id','single_polypeptide')
    write_tsv(output/'catalytic_sites.tsv',[serialized(s) for r in selected for s in r['sites']],SITE_FIELDS)
    write_tsv(output/'active_center_summary.tsv',[serialized(r) for r in selected],SUMMARY_FIELDS)
    if evidence is not None:write_tsv(output/'active_center_template_evidence.tsv',[serialized(r) for r in evidence],SUMMARY_FIELDS+['bitscore','sites'])


def annotate_active_centers(fasta, output, threads=1, root=ROOT, anchors=None):
    root=Path(root);output=Path(output);directory=root/'references/processed/active_centers'
    manifest=json.loads((directory/'manifest.json').read_text())
    if manifest.get('status')!='complete':
        raise ValueError('Active-center model is incomplete')
    raw=output/'active_center_blast.tsv'
    with (output/'active_center_blast.log').open('w') as log:
        subprocess.run(['blastp','-query',str(fasta),'-db',str(directory/'seed'),'-evalue',str(manifest['parameters']['max_evalue']),
            '-max_target_seqs','1000','-max_hsps','20','-num_threads',str(threads),'-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(raw)],stderr=log,check=True)
    thresholds=dict(blast_evalue=manifest['parameters']['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0)
    hits=load_blast(raw,thresholds)
    if anchors is not None:hits=[h for h in hits if h['qseqid'] in anchors]
    sequences={r.id:str(r.seq) for r in SeqIO.parse(fasta,'fasta')}
    selected,evidence=map_hits(hits,json.loads((directory/'models.json').read_text()),manifest['parameters'],query_sequences=sequences)
    for row in selected:
        for variant in manifest.get('exact_active_reference_variants',[]):
            if hashlib.sha256(sequences[row['protein_id']].encode()).hexdigest()!=variant['sequence_sha256']:continue
            changed=[s for s in row['sites'] if s['residue_role']=='direct_metal_ligand' and s['status']!='present']
            if (row['active_center_status']=='ligands_substituted' and row['reference_id']==variant['template_id'] and len(changed)==1 and
                changed[0]['reference_position']==variant['changed_reference_position'] and changed[0]['query_position']==variant['query_position'] and
                changed[0]['observed_residue']==variant['query_residue']):
                row.update(active_center_status='experimentally_active_exact_reference_noncanonical_site',exact_active_variant_review=variant)
    write_active_center_tables(output,selected,evidence)
    return selected,manifest


if __name__=='__main__':
    print(json.dumps(build_model(),indent=2))
