"""Reference-like canonical FeGP core and H4MPT potential are separate claims."""
import json,math,hashlib
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT,load_blast
from localization import sha
def _validate_hit(hit, queries, references):
    if hit['qseqid'] not in queries or hit['sseqid'] not in references: raise ValueError('Unknown BLAST sequence ownership')
    if len(hit['qseq'])!=len(hit['sseq']): raise ValueError('BLAST alignment lengths differ')
    for field, sequence, lo, hi in [('qseq',queries[hit['qseqid']],hit['qstart'],hit['qend']),('sseq',references[hit['sseqid']],hit['sstart'],hit['send'])]:
        observed=hit[field].replace('-','').replace('U','X').replace('O','X')
        if not 1<=lo<=hi<=len(sequence) or observed!=sequence[lo-1:hi].replace('U','X').replace('O','X'):
            raise ValueError('BLAST alignment differs from bound '+field+' sequence')
def _partial(gene):
    return str((gene or {}).get('partial', '')).strip().lower() in {'true', '1', 'yes'}
from coding_completeness import supports as coding_supports, supports_translated_region
from catalytic_confirmation import core_similarity
from active_centers import map_model,collinear_hits

RULE='canonical_Hmd_FeGP_H4MPT_potential_v4'


def _canonical_group(model, root, refs, prototype_records):
    prototypes={r['prototype_id']:r['sequence'] for r in prototype_records}
    ids=model.get('canonical_prototypes')
    review_path=model.get('canonical_prototype_review')
    if not isinstance(ids,list) or not ids or len(ids)!=len(set(ids)) or not review_path:
        raise ValueError('Canonical Hmd prototype group is missing or nonunique')
    review=json.loads((root/review_path).read_text()) if (root/review_path).is_file() else {}
    rows=review.get('references',review) if isinstance(review,dict) else review
    if not isinstance(rows,list):
        rows=[dict(r,sequence_sha256=sha_text(r['sequence'])) for r in prototype_records if r.get('prototype_id') in set(ids)]
    reviewed={r.get('prototype_id'):r for r in rows}
    if len(reviewed)!=len(rows) or len(prototypes)!=len(prototype_records): raise ValueError('Duplicate canonical Hmd review or prototype records')
    if set(reviewed)!=set(ids): raise ValueError('Canonical Hmd prototype review/group mismatch')
    old=model.get('canonical_reference','Fe_Hmd__Hmd')
    if old not in refs: raise ValueError('Original canonical Hmd reference is missing')
    old_matches=[rid for rid,seq in prototypes.items() if seq==refs[old]]
    if len(old_matches)!=1 or old_matches[0] not in ids: raise ValueError('Original canonical Hmd prototype was removed or nonunique')
    for rid in ids:
        if rid not in prototypes or rid not in reviewed: raise ValueError('Canonical Hmd prototype missing from bound banks')
        raw=next(r for r in prototype_records if r['prototype_id']==rid)
        if raw.get('label')!='positive' or raw.get('family')!='Fe': raise ValueError('Canonical Hmd prototype bank record is not positive Fe')
        if reviewed[rid].get('label')!='positive' or reviewed[rid].get('family')!='Fe': raise ValueError('Canonical Hmd prototype review is not positive Fe')
        if reviewed[rid].get('sequence_sha256') != sha_text(prototypes[rid]): raise ValueError('Canonical Hmd prototype sequence hash mismatch')
    return ids,reviewed


def sha_text(value):
    import hashlib
    return hashlib.sha256(value.encode()).hexdigest()

def branch_comparison(feature,hits,queries,records,canonical_ids,model):
    """Require concordant complete canonical references to beat every saved rival."""
    bank={r['prototype_id']:r['sequence'] for r in records};metadata={r['prototype_id']:r for r in records}
    positive={};rivals=[];q=feature['protein_id']
    for h in hits:
        if h['qseqid']!=q or h['evalue']>model['max_evalue']:continue
        _validate_hit(h,queries,bank)
        if h['qlen']!=len(queries[q]) or h['slen']!=len(bank[h['sseqid']]) or not all(math.isfinite(h[k]) for k in ('pident','evalue','bitscore')):raise ValueError('Canonical branch alignment length/score mismatch')
        score=core_similarity(h,feature['core_positions'])
        pairs=sum(a not in '-XBZJO' and b not in '-XBZJO' for a,b in zip(h['qseq'],h['sseq']))
        coverage=dict(query=pairs/h['qlen'],reference=pairs/h['slen'])
        row=dict(reference_id=h['sseqid'],reference_identity_label=metadata[h['sseqid']]['label'],reference_family=metadata[h['sseqid']]['family'],score=score,paired_coverage=coverage,raw_alignment=h)
        if h['sseqid'] in canonical_ids:
            if min(coverage.values())>=model['min_coverage'] and score['core_coverage']>=model['min_coverage']:
                old=positive.get(h['sseqid'])
                if old is None or row['score']['similarity']>old['score']['similarity']:positive[h['sseqid']]=row
        else:rivals.append(row)
    floor=min((r['score']['similarity'] for r in positive.values()),default=0)
    counter=max((r['score']['similarity'] for r in rivals),default=0)
    blocks=[]
    if set(positive)!=set(canonical_ids):blocks.append('complete_concordant_canonical_reference_group_not_supported')
    if floor<=counter:blocks.append('canonical_branch_does_not_beat_all_alternative_prototypes')
    return dict(blockers=blocks,canonical_reference_scores=positive,alternative_reference_scores=rivals,weakest_canonical_score=floor,strongest_alternative_score=counter,margin=floor-counter,score_scope='Matching residues per same query core; alternate/indeterminate prototypes are competitors, not assumed inactive.',reference_transfer_identity_threshold_used=False)

def prepare(annotation,features,genes,root=ROOT,coding_completeness=None):
    annotation,root=Path(annotation),Path(root)
    model_path=root/'research_workflow/canonical_hmd_model.json';model=json.loads(model_path.read_text())
    if model['version']!=RULE:raise ValueError('Canonical Hmd rule version mismatch')
    rows=features if isinstance(features,list) else list(features.values());rows=[f for f in rows if f['family']=='Fe']
    if not rows:return {},{}
    with (annotation/'proteins.faa').open() as handle:queries={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    with (root/'references/processed/references.faa').open() as handle:refs={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    prototype_records=json.loads((root/'references/processed/catalytic_confirmation/prototypes.json').read_text())
    prototypes={r['prototype_id']:r['sequence'] for r in prototype_records}
    canonical_ids,prototype_review=_canonical_group(model,root,refs,prototype_records)
    branch_review=json.loads((root/'references/processed/canonical_hmd_review.json').read_text())
    if model.get('branch_rule')!=branch_review['version'] or branch_review['version']!='canonical_Hmd_concordant_branch_v1' or branch_review['positive_source_groups']!=len(canonical_ids):raise ValueError('Canonical branch rule/review mismatch')
    site_models={r['reference_id']:r for r in json.loads((root/'references/processed/active_centers/models.json').read_text())}
    templates={rid:r['sequence'] for rid,r in site_models.items()}
    site_parameters=json.loads((root/'config/active_center_review.json').read_text())['parameters']
    genes={g['id']:g for g in genes};selected={f['protein_id']:[] for f in rows}
    paths=[annotation/'proteins.faa',model_path,Path(__file__)]
    search_verified=False;feature_manifest=annotation/'confirmation/feature_manifest.json'
    if feature_manifest.is_file():
        paths.append(feature_manifest);fm=json.loads(feature_manifest.read_text());blast=annotation/'confirmation/prototype_blast.tsv'
        search_verified=(fm.get('status')=='features_extracted' and fm.get('original_query_sha256')==sha(annotation/'proteins.faa') and fm.get('model_manifest_sha256')==sha(root/'references/processed/catalytic_confirmation/manifest.json') and blast.is_file() and fm.get('artifacts',{}).get('prototype_blast.tsv')==sha(blast))
    raw_by_file={}
    for filename,rid,bank in [('confirmation/prototype_blast.tsv',None,prototypes),('active_center_blast.tsv',model['template'],templates)]:
        path=annotation/filename
        if not path.is_file():continue
        paths.append(path)
        raw_by_file[filename]=load_blast(path,dict(blast_evalue=site_parameters['max_evalue'] if filename=='active_center_blast.tsv' else model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
        for hit in raw_by_file[filename]:
            target_id=hit['sseqid']
            if hit['qseqid'] not in selected or (rid is not None and target_id!=rid) or (rid is None and target_id not in canonical_ids):continue
            _validate_hit(hit,queries,bank)
            if hit['qlen']!=len(queries[hit['qseqid']]) or hit['slen']!=len(bank[target_id]):raise ValueError('Canonical Hmd sequence length mismatch')
            pairs=[(q,s) for q,s in zip(hit['qseq'],hit['sseq']) if q not in '-XBZJO' and s not in '-XBZJO']
            coverage=dict(query=len(pairs)/len(queries[hit['qseqid']]),reference=len(pairs)/len(bank[target_id]))
            identity=100*sum(q==s for q,s in pairs)/len(pairs) if pairs else 0
            selected[hit['qseqid']].append(dict(reference_id=target_id,identity=identity,paired_coverage=coverage,
                supported=hit['evalue']<=model['max_evalue'] and identity>=model['min_identity'] and min(coverage.values())>=model['min_coverage'],raw_alignment=hit,source_path=str(path)))
    results={}
    for f in rows:
        q=f['protein_id'];center=f.get('source_active_center',{});blocks=[]
        if f.get('positive_prototype'):
            winners=[h for h in raw_by_file.get('confirmation/prototype_blast.tsv',[]) if h['qseqid']==q and h['sseqid']==f['positive_prototype']]
            if len(winners)!=1:raise ValueError('Canonical Hmd positive core alignment missing or nonunique')
            hit=winners[0];_validate_hit(hit,queries,prototypes);calculated=core_similarity(hit,f['core_positions'])
            if any((f.get('positive_evidence') or {}).get(k)!=v for k,v in calculated.items()) or f['positive_score']!=calculated['similarity']:
                raise ValueError('Canonical Hmd saved core evidence differs from actual alignment')
        proof_used=False; coding_region_used=False
        proof_ok=bool(coding_completeness and q in coding_completeness and q in genes and coding_supports(coding_completeness[q],q,expected_protein_sha256=hashlib.sha256(queries[q].encode()).hexdigest(),feature=f,gene=genes[q]))
        region_ok=bool(coding_completeness and q in coding_completeness and q in genes and not proof_ok and
                       supports_translated_region(coding_completeness[q],q,queries[q],feature=f,gene=genes[q]))
        if region_ok:
            complete_ref=any(h.get('supported') and h.get('paired_coverage',{}).get('reference') == 1.0
                             for h in selected[q])
            if complete_ref: coding_region_used=True
            else: blocks.append('translated_region_requires_complete_canonical_reference_alignment')
        if q not in genes or (_partial(genes[q]) and not proof_ok and not coding_region_used):blocks.append('partial_or_unknown_CDS')
        elif _partial(genes[q]) and proof_ok: proof_used=True
        if not f.get('positive_prototype') or f.get('positive_score',0)<=0 or f.get('competitor_margin',0)<=0:blocks.append('positive_core_competition_not_supported')
        if f.get('positive_prototype') not in canonical_ids:blocks.append('canonical_positive_prototype_not_selected')
        allowed={'family_outside_calibration_scope','reviewed_single_chain_center_or_span_not_supported','required_domain_architecture_not_colocated_with_center'}|({'partial_CDS'} if (proof_used or coding_region_used) else set())
        if set(f.get('eligibility_reasons',[]))-allowed:blocks.append('core_evidence_conflict')
        if center.get('center_organization')!='single_polypeptide':blocks.append('single_chain_not_supported')
        active_hits=[h for h in raw_by_file.get('active_center_blast.tsv',[]) if h['qseqid']==q and h['sseqid']==model['template']]
        mapped=map_model(collinear_hits(active_hits),site_models[model['template']],site_parameters,query_sequence=queries[q]) if active_hits else {}
        if not mapped:blocks.append('raw_canonical_site_alignment_missing')
        for key in ['direct_ligand_status','reaction_signature_status']:
            if mapped and center.get('reference_id')==model['template'] and mapped[key]!=center.get(key):raise ValueError('Canonical Hmd saved ligand evidence differs from actual alignment')
        center=dict(center,**{k:v for k,v in mapped.items() if k in ['direct_ligand_status','reaction_signature_status','site_contexts']})
        contexts=center.get('site_contexts',[]);contexts=json.loads(contexts) if isinstance(contexts,str) else contexts
        contexts={c['reference_position']:c for c in contexts}
        def local(position):
            c=contexts.get(position,{})
            return c.get('coverage',0)>=model['context_min_coverage'] and c.get('identity',0)>=model['context_min_identity']
        if center.get('reference_id')!=model['template'] or center.get('direct_ligand_status')!='present' or not local(model['ligand_position']):blocks.append('FeGP_ligand_context_not_supported')
        if not any(h['supported'] for h in selected[q]):blocks.append('full_canonical_reference_match_not_supported')
        reaction_blocks=list(blocks)
        if center.get('reaction_signature_status')!='present' or not local(model['signature_position']):reaction_blocks.append('canonical_H14_reaction_context_not_supported')
        branch=branch_comparison(f,raw_by_file.get('confirmation/prototype_blast.tsv',[]),queries,prototype_records,canonical_ids,model)
        branch_blocks=[b for b in blocks if b!='full_canonical_reference_match_not_supported']+branch['blockers']
        if not search_verified:branch_blocks.append('complete_current_prototype_search_not_verified')
        if f.get('core_domain_coverage',0)<model['min_coverage']:branch_blocks.append('complete_Hmd_domain_not_supported')
        branch_reaction_blocks=list(branch_blocks)
        if 'canonical_H14_reaction_context_not_supported' in reaction_blocks:branch_reaction_blocks.append('canonical_H14_reaction_context_not_supported')
        results[q]=dict(candidate_id=q,rule_version=RULE,core_status='reference_supported_uncalibrated' if not blocks else 'unresolved',
            reaction_status='compatible_reaction_candidate_uncalibrated' if not reaction_blocks else 'unresolved',core_blockers=blocks,reaction_blockers=reaction_blocks,
            branch_core_status='sequence_supported_candidate' if not branch_blocks else 'unresolved',branch_reaction_status='compatible_reaction_candidate_uncalibrated' if not branch_reaction_blocks else 'unresolved',branch_core_blockers=branch_blocks,branch_reaction_blockers=branch_reaction_blocks,branch_comparison=branch,
            canonical_reference_evidence=selected[q],canonical_reference_group=dict(reference_ids=canonical_ids,review=prototype_review),selected_canonical_positive_prototype=f.get('positive_prototype'),source_active_center=center,remapped_sites=mapped.get('sites',[]),positive_core_reference=f.get('positive_evidence'),competitor_margin=f.get('competitor_margin'),
            coding_completeness_used=proof_used,coding_region_used=coding_region_used,coding_completeness_scope='uncalibrated coding-completeness support, not gene-caller correction or calibrated identity' if proof_used else 'translated DNA region support only; source partial/CDS status retained' if coding_region_used else 'not_used',
            limits='Reference transfer and concordant branch hypotheses are distinct, uncalibrated predictions; query activity, FeGP occupancy, substrate exclusivity and in-vivo direction untested. Full-CDS completeness is separate from a source-bound complete catalytic region; the translated-region route preserves partial and missing-stop status. His14 is claim-specific support, not a universal activity requirement. Branch evidence is limited to the searched prototype bank and two reviewed canonical source groups.')
    return results,{str(p):sha(p) for p in paths}

def attributes(result):
    return [dict(attribute=name,status=status,value=result,evidence_dependencies=dependencies,evidence_source='canonical_hmd_evidence.json',limitations=result['limits'])
        for name,status,dependencies in [
            ('canonical_Hmd_FeGP_core','sequence_supported_candidate' if result['core_status']=='reference_supported_uncalibrated' else 'unresolved',['canonical_full_sequence_reference','positive_negative_core_competition','Cys176_context','complete_CDS_or_bound_complete_catalytic_region']),
            ('FeGP_H4MPT_hydrogenation_potential',result['reaction_status'],['canonical_Hmd_FeGP_core','H14_reaction_context','experimentally_supported_H4MPT_reference_reaction']),
            ('canonical_Hmd_branch_core',result['branch_core_status'],['complete_concordant_canonical_references','all_prototype_competition','FeGP_ligand_context','complete_CDS_or_bound_complete_catalytic_region']),
            ('canonical_H4MPT_branch_reaction_potential',result['branch_reaction_status'],['canonical_Hmd_branch_core','H14_reaction_context','experimentally_supported_H4MPT_reference_reaction'])]]

def identity_evidence(annotation,features,genes,root=ROOT,coding_completeness=None):
    """Use the same Fe evidence in the core and additional-candidate entry points."""
    rows=features if isinstance(features,list) else list(features.values())
    if not any(f['family']=='Fe' for f in rows):return {},{}
    proof_rows,proof_bound=coding_completeness,{}
    if isinstance(coding_completeness,tuple) and len(coding_completeness)==2:
        proof_rows,proof_bound=coding_completeness
    elif isinstance(coding_completeness,dict) and isinstance(coding_completeness.get('rows'),dict):
        proof_rows=coding_completeness['rows'];proof_bound=coding_completeness.get('input_hashes',{})
    canonical,bound=prepare(annotation,features,genes,root,coding_completeness=proof_rows)
    hmdii={};bound.update(proof_bound)
    return {q:dict(canonical=canonical.get(q),hmdii=hmdii.get(q)) for q in set(canonical)|set(hmdii)},bound
