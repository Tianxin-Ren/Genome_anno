"""Source-pinned split NiFe centers; residue evidence does not establish complex activity."""
import copy
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT,BLAST_FIELDS,extract_genbank,load_blast,association_allowed,write_tsv,_circular_contig_length
from active_centers import sha,split_structural_record,map_model,collinear_hits,write_active_center_tables
from circular_coordinates import valid_cds_interval

CANDIDATE_FIELDS='catalytic_anchor center_id center_status pairing_status organization component_proteins component_status partner_window_score candidate_bitscore score_scope source source_sha256 sites'.split()
COMPONENT_FIELDS='catalytic_anchor center_id component_id protein_id component_status reference_id reference_accession window_coverage window_identity template_span_status alignment_components source source_sha256'.split()


def source_genes(root,case,resolve_reviewed=True):
    genes=extract_genbank(root/case['genome_path'])[0]
    if resolve_reviewed and case['resolved_positions']:
        gene=next(g for g in genes if g['id']==case['partner_id'])
        if hashlib.sha256(gene['sequence'].encode()).hexdigest()!=case['raw_partner_sequence_sha256']:raise ValueError('Historical split-CDS source changed')
        record=SeqIO.read(root/case['genome_path'],'genbank');sequence=list(gene['sequence'])
        for pos,review in case['resolved_positions'].items():
            i=int(pos)-1;dna=str(record.seq[gene['start']-1+3*i:gene['start']-1+3*i+3])
            if gene['strand']!=1 or sequence[i]!=review['raw'] or dna!=review['genomic_codon']:raise ValueError('Reviewed Sec correction no longer matches exact CDS/codon')
            sequence[i]=review['reviewed']
        gene['sequence']=''.join(sequence)
        if hashlib.sha256(gene['sequence'].encode()).hexdigest()!=case['reviewed_partner_sequence_sha256']:raise ValueError('Reviewed protein sequence mismatch')
    return genes


def build_model(root=ROOT):
    root=Path(root);config_path=root/'config/split_centers.json';policy=json.loads(config_path.read_text())
    for path,source in policy['source_files'].items():
        if sha(root/path)!=source['sha256']:raise ValueError('Split-center source changed: '+path)
    template=policy['template'];structure=split_structural_record(root/template['path'],[c['entity_id'] for c in template['components']],template['nickel_chain'],template['nickel_author_position'])
    if structure['iron']['chain']!=template['iron_chain'] or structure['iron']['author_position']!=template['iron_author_position']:raise ValueError('Reviewed paired iron changed')
    components=[];records=[]
    for spec,physical in zip(template['components'],structure['components']):
        gene=next(g for g in extract_genbank(root/spec['genome_path'])[0] if g['id']==spec['locus_tag'])
        if (gene['protein_accession']!=spec['protein_accession'] or hashlib.sha256(gene['sequence'].encode()).hexdigest()!=spec['sequence_sha256'] or
            physical['canonical_chain']!=spec['chain'] or physical['ligands']!=spec['expected_ligands'] or
            not gene['sequence'].startswith(physical['sequence'])):raise ValueError('Split component CDS/structure correspondence changed')
        model=dict(physical,sequence=gene['sequence'],reference_id=template['center_id']+'__'+spec['component_id'],family='NiFe',
            component_id=spec['component_id'],reference_accession=spec['protein_accession'],window_start=spec['window_start'],window_end=spec['window_end'],
            source_sha256=structure['source_sha256'],evidence_url=template['evidence_url'],
            sites=[dict(reference_position=int(p),reference_residue=residue,expected_residue='CU',residue_role='direct_metal_ligand',
                requirement='required_for_template_coordination') for p,residue in sorted(physical['ligands'].items(),key=lambda x:int(x[0]))])
        if not 1<=model['window_start']<=min(int(p) for p in physical['ligands'])<=max(int(p) for p in physical['ligands'])<=model['window_end']<=len(physical['sequence']):raise ValueError('Split template window excludes reviewed sites')
        components.append(model);records.append(SeqRecord(Seq(gene['sequence']),id=model['reference_id'],description='source_CDS_full_precursor; structural_window_separate'))
    if [len(c['sites']) for c in components]!=[2,2]:raise ValueError('Split center is not the reviewed two-plus-two ligand set')
    directory=root/'references/processed/split_centers';directory.mkdir(parents=True,exist_ok=True)
    (directory/'models.json').write_text(json.dumps(components,indent=2)+'\n');(directory/'structural_truth.json').write_text(json.dumps(structure,indent=2)+'\n')
    SeqIO.write(records,directory/'components.faa','fasta')
    with (directory/'makeblastdb.log').open('w') as log:subprocess.run(['makeblastdb','-in',str(directory/'components.faa'),'-dbtype','prot','-out',str(directory/'seed')],stdout=log,stderr=log,check=True)
    manifest=dict(status='complete',policy_sha256=sha(config_path),center_id=template['center_id'],parameters=policy['parameters'],
        files={p.name:sha(p) for p in sorted(directory.iterdir()) if p.is_file() and p.name not in ['manifest.json','makeblastdb.log']},
        source_files=policy['source_files'],scope=policy['scope'])
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');return manifest


def component_link(a,u,genes,mode,parameters):
    if a['protein_id']==u['protein_id']:
        return max(h['query_end'] for h in a['sites'][0]['alignment_components'])<min(h['query_start'] for h in u['sites'][0]['alignment_components'])
    if mode=='transcript':return False
    ga,gu=genes[a['protein_id']],genes[u['protein_id']]
    if not association_allowed(ga,gu,mode) or ga['contig']!=gu['contig'] or ga['strand']!=gu['strand']:return False
    left,right=(ga,gu) if ga['strand']==1 else (gu,ga)
    circular_length=_circular_contig_length([ga,gu])
    if circular_length is not None:
        if not circular_length or not all(valid_cds_interval(g, allow_unknown_length=False) for g in (left, right)):
            return False
        offset = (right['start'] - left['start']) % circular_length
        gap = offset - (left['end'] - left['start'] + 1)
        span = offset + right['end'] - right['start'] + 1
        return 0 <= gap <= parameters['max_intergenic_gap_bp'] and span <= min(circular_length, parameters['max_span_bp'])
    return (0<=right['start']-left['end']-1<=parameters['max_intergenic_gap_bp'] and
        right['end']-left['start']+1<=parameters['max_span_bp'])


def annotate_split_centers(fasta,output,genes,annotations,selected,mode,threads=1,root=ROOT):
    root=Path(root);output=Path(output);directory=root/'references/processed/split_centers';policy=json.loads((root/'config/split_centers.json').read_text());manifest=json.loads((directory/'manifest.json').read_text())
    if manifest.get('status')!='complete':raise ValueError('Split-center model is incomplete')
    models=json.loads((directory/'models.json').read_text());parameters=manifest['parameters'];raw=output/'split_center_blast.tsv'
    with (output/'split_center_blast.log').open('w') as log:
        subprocess.run(['blastp','-task','blastp-short','-seg','no','-query',str(fasta),'-db',str(directory/'seed'),'-evalue',str(parameters['max_evalue']),
            '-max_target_seqs','1000','-max_hsps','20','-num_threads',str(threads),'-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(raw)],stdout=log,stderr=log,check=True)
    hits=load_blast(raw,dict(blast_evalue=parameters['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0));groups=defaultdict(list)
    for h in hits:
        if h['evalue']<=parameters['max_evalue']:groups[h['qseqid'],h['sseqid']].append(h)
    sequences={g['id']:g['sequence'] for g in genes}
    mapped={key:map_model(collinear_hits(values),next(m for m in models if m['reference_id']==key[1]),parameters,query_sequence=sequences[key[0]]) for key,values in groups.items()}
    old={r['protein_id']:r for r in selected};replacements={};candidates=[];component_rows=[];by_id={g['id']:g for g in genes}
    for anchor in annotations:
        qid=anchor['protein_id'];anchor.update(split_center_status='not_detected_or_not_evaluable',split_center_component_proteins='')
        a=mapped.get((qid,models[0]['reference_id']))
        if anchor['catalytic_family']!='NiFe' or not a or a['template_span_alignment_status']!='supported':continue
        if old.get(qid,{}).get('direct_ligand_status')=='present':
            anchor['split_center_status']='single_chain_ligands_already_mapped';continue
        partners=[u for (q,rid),u in mapped.items() if rid==models[1]['reference_id'] and component_link(a,u,by_id,mode,parameters)]
        # Mature structural windows, not precursor tails or the shared large chain, rank partners.
        partner_score=lambda u:u['window_coverage']*u['window_identity']
        partners.sort(key=lambda u:(-partner_score(u),u['protein_id']))
        competing=[u for u in partners if partner_score(u)>=partner_score(partners[0])*parameters['competitive_bitscore_ratio']] if partners else []
        anchor['split_center_status']='ambiguous_component_pairing' if len(competing)>1 else 'required_component_not_detected' if not partners else 'component_ligand_evidence_incomplete'
        anchor['split_center_component_proteins']=';'.join(dict.fromkeys([qid,*[u['protein_id'] for u in competing]]))
        if not partners:
            cid=qid+':'+manifest['center_id']+':missing_component';model=models[1]
            missing=[dict(s,protein_id='',catalytic_anchor=qid,reference_id=model['reference_id'],family='NiFe',source_center_id=manifest['center_id'],
                center_component_id=model['component_id'],query_position=None,observed_residue='',status='not_covered',coverage_reason='required_component_not_detected_or_unpredicted',
                evidence_url=model['evidence_url'],source_sha256=model['source_sha256'],cofactor_occupancy='unknown') for s in model['sites']]
            present=[dict(s,catalytic_anchor=qid,source_center_id=manifest['center_id'],center_component_id=models[0]['component_id']) for s in a['sites']]
            candidates.append(dict(catalytic_anchor=qid,center_id=cid,center_status='required_component_not_detected',pairing_status='unresolved',organization='incomplete',partner_window_score=0,
                component_proteins=json.dumps([qid]),component_status=json.dumps([a['active_center_status'],'not_detected_or_unpredicted']),candidate_bitscore=a['bitscore'],
                score_scope='A_component_only; no_partner_evidence; not_combined_significance',source=model['evidence_url'],source_sha256=model['source_sha256'],sites=json.dumps([*present,*missing],sort_keys=True)))
            component_rows.append(dict(catalytic_anchor=qid,center_id=cid,component_id=models[0]['component_id'],protein_id=qid,component_status=a['active_center_status'],
                reference_id=models[0]['reference_id'],reference_accession=models[0]['reference_accession'],window_coverage=a['window_coverage'],window_identity=a['window_identity'],template_span_status=a['template_span_alignment_status'],
                alignment_components=json.dumps(a['sites'][0]['alignment_components']),source=models[0]['evidence_url'],source_sha256=models[0]['source_sha256']))
            component_rows.append(dict(catalytic_anchor=qid,center_id=cid,component_id=model['component_id'],protein_id='',component_status='not_detected_or_unpredicted',
                reference_id=model['reference_id'],reference_accession=model['reference_accession'],window_coverage=0,window_identity=0,template_span_status='not_evaluable',alignment_components='[]',source=model['evidence_url'],source_sha256=model['source_sha256']))
        for i,u in enumerate(partners):
            cid=qid+':'+manifest['center_id']+':'+str(i);parts=[a,u];sites=[];organization='single_polypeptide' if qid==u['protein_id'] else 'multiple_polypeptides'
            for model,row in zip(models,parts):
                for site in row['sites']:
                    sites.append(dict(site,catalytic_anchor=qid,source_center_id=manifest['center_id'],center_component_id=model['component_id']))
                component_rows.append(dict(catalytic_anchor=qid,center_id=cid,component_id=model['component_id'],protein_id=row['protein_id'],component_status=row['active_center_status'],
                    reference_id=model['reference_id'],reference_accession=model['reference_accession'],window_coverage=row['window_coverage'],window_identity=row['window_identity'],
                    template_span_status=row['template_span_alignment_status'],alignment_components=json.dumps(row['sites'][0]['alignment_components']),source=model['evidence_url'],source_sha256=model['source_sha256']))
            valid=all(p['active_center_status']=='template_ligand_set_supported' and p['template_span_alignment_status']=='supported' for p in parts)
            status='split_template_ligand_set_supported' if valid else 'component_ligand_evidence_incomplete'
            if len(competing)>1 and u in competing:status='ambiguous_component_pairing'
            candidates.append(dict(catalytic_anchor=qid,center_id=cid,center_status=status,pairing_status='competitive' if u in competing else 'below_competitive_score',organization=organization,
                component_proteins=json.dumps([p['protein_id'] for p in parts]),component_status=json.dumps([p['active_center_status'] for p in parts]),partner_window_score=partner_score(u),candidate_bitscore=sum(p['bitscore'] for p in parts),
                score_scope='component_bitscore_sum_not_combined_significance; pairing_compares_U_structural_window_only',source=models[0]['evidence_url'],source_sha256=models[0]['source_sha256'],sites=json.dumps(sites,sort_keys=True)))
            if not valid or len(competing)!=1 or u is not competing[0]:continue
            previous_positions={s['query_position'] for s in old.get(qid,{}).get('sites',[]) if s['residue_role']=='direct_metal_ligand' and s['status']=='present'}
            if not previous_positions<={s['query_position'] for s in a['sites']}:
                anchor['split_center_status']='conflicts_with_single_chain_ligand_mapping';continue
            row=copy.deepcopy(a);row.update(reference_id=manifest['center_id'],sites=sites,
                active_center_status='multi_polypeptide_ligand_set_supported' if organization=='multiple_polypeptides' else 'fused_split_template_ligand_set_supported',
                center_organization=organization,center_component_proteins=list(dict.fromkeys(p['protein_id'] for p in parts)),competing_templates=[manifest['center_id']],
                window_coverage=min(p['window_coverage'] for p in parts),window_identity=min(p['window_identity'] for p in parts),
                min_required_context_coverage=min(p['min_required_context_coverage'] for p in parts),min_required_context_identity=min(p['min_required_context_identity'] for p in parts),
                confidence='provisional_multi_component_mapping_not_measured_assembly_or_activity',
                template_span_alignment_status='supported',bitscore=sum(p['bitscore'] for p in parts),site_contexts=[dict(c,component_id=m['component_id'],protein_id=p['protein_id']) for m,p in zip(models,parts) for c in p['site_contexts']])
            for site in row['sites']:site['model_status']=row['active_center_status']
            replacements[qid]=row;anchor.update(split_center_status=row['active_center_status'],split_center_component_proteins=';'.join(row['center_component_proteins']))
            anchor['flags']=';'.join(filter(None,[anchor.get('flags',''),'legacy_seed_diagnostics_are_single_polypeptide; split_center_evidence_separate']))
    merged=[replacements.get(r['protein_id'],r) for r in selected]
    merged.extend(r for q,r in replacements.items() if q not in old)
    write_active_center_tables(output,merged)
    write_tsv(output/'catalytic_center_candidates.tsv',candidates,CANDIDATE_FIELDS);write_tsv(output/'catalytic_center_components.tsv',component_rows,COMPONENT_FIELDS)
    return merged,manifest


if __name__=='__main__':print(json.dumps(build_model(),indent=2))
