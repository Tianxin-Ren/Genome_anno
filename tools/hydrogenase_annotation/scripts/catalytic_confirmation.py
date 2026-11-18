"""Separate source-pinned catalytic identity evidence from reaction classification."""
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, BLAST_FIELDS, load_blast, residue_states, write_tsv
from active_centers import collinear_hits
from audit_maturation_domains import scan_profile_domains
from audit_reference_domains import model_evidence

FEATURE_FIELDS=('protein_id family eligible eligibility_reasons core_positions core_domain_coverage positive_score '
    'competitor_score competitor_margin positive_prototype competitor_prototype positive_evidence competitor_evidence '
    'domain_evidence source_active_center').split()
CONFIRMATION_FIELDS=FEATURE_FIELDS+['confirmation_status','identity_confidence','rule_parameters','rule_scope']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stage_fingerprint(root=ROOT):
    from discovery_evidence import stage_fingerprint as discovery_fingerprint
    root=Path(root)
    files=['scripts/catalytic_confirmation.py','scripts/identity_truth.py','config/catalytic_identity_adjudications.json','scripts/annotate.py','scripts/active_centers.py','scripts/split_centers.py',
        'scripts/module_graph.py','scripts/complex_selection.py','config/complex_selection.json','scripts/module_architecture.py','config/module_architecture.json','scripts/audit_reference_domains.py','scripts/audit_maturation_domains.py',
        'scripts/noncanonical_center_evidence.py','config/noncanonical_center_evidence.json',
        'scripts/complex_extensions.py','scripts/build_complex_extensions.py','config/complex_extension_sources.json','config/reaction_contracts.json','config/exact_reference_functions.json','references/processed/complex_extensions/manifest.json',
        'config/catalytic_confirmation.json','config/reference_catalytic_identity_review.json','config/reference_domain_policy.json','config/thresholds.json','config/systems.json',
        'references/processed/references.json','references/processed/references.faa',
        'references/processed/catalytic_confirmation/manifest.json','references/processed/active_centers/manifest.json','references/processed/split_centers/manifest.json']
    return dict(files={name:sha(root/name) for name in files},discovery=discovery_fingerprint(root))


def table(path):
    import csv
    with Path(path).open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def source_family(row):
    return next((f for f in ['NiFe','FeFe','Fe'] if any(c==f or c.startswith('['+f+']') for c in row['categories'])), '')


def build_model(root=ROOT):
    root=Path(root);policy_path=root/'config/catalytic_confirmation.json';policy=json.loads(policy_path.read_text())
    corpus=root/policy['corpus'];frozen=json.loads((corpus/'freeze.json').read_text())
    for name,digest in frozen['artifacts'].items():
        if sha(corpus/name)!=digest:raise ValueError('Frozen confirmation source changed: '+name)
    audit_path=root/policy['training_audit']/'summary.json';audit=json.loads(audit_path.read_text());independent=set(audit['independent_unit_ids'])
    if audit['corpus_freeze_sha256']!=sha(corpus/'freeze.json'):raise ValueError('Confirmation training audit changed')
    by_sequence={}
    from identity_truth import adjudications
    identity_review=adjudications(root)
    def add(sequence,label,family,source,component='',genera=()):
        digest=hashlib.sha256(sequence.encode()).hexdigest()
        correction=identity_review.get(digest)
        if correction:
            source=dict(source,frozen_identity_label=label,identity_adjudication=correction)
            label=correction['identity_label'];family=correction['family']
        row=by_sequence.setdefault(digest,dict(prototype_id='prototype_'+digest[:20],sequence=sequence,sequence_sha256=digest,
            label=label,family=family,sources=[],homolog_components=[],genus_taxids=[]))
        if row['label']!=label or row['family']!=family:raise ValueError('Conflicting confirmation prototype labels')
        row['sources'].append(source)
        if component:row['homolog_components']=sorted(set(row['homolog_components'])|{component})
        row['genus_taxids']=sorted(set(row['genus_taxids'])|set(genera))
    units=json.loads((corpus/'units.json').read_text());byid={r['unit_id']:r for r in units}
    for row in units:
        if row['split']==policy['training_split'] and row['unit_id'] in independent:
            add(row['sequence'],row['label'],source_family(row),dict(unit_id=row['unit_id'],tier=row['evidence_tiers'],accessions=[s['accession'] for s in row['sources']]),row['homolog_component'],row['genus_taxids'])
    sequences={r.id:str(r.seq) for r in SeqIO.parse(root/'references/processed/references.faa','fasta')}
    refs=json.loads((root/'references/processed/references.json').read_text())
    review_path=root/policy['reference_identity_review'];review=json.loads(review_path.read_text());identity_reviews={r['reference_id']:r for r in review['records']}
    if review['corpus_freeze_sha256']!=sha(corpus/'freeze.json') or set(identity_reviews)!={r['reference_id'] for r in refs if r['kind']=='unresolved'}:raise ValueError('Reference identity review changed')
    for row in refs:
        if row['kind'] not in ['catalytic','decoy','unresolved']:continue
        sequence=sequences[row['reference_id']]
        if hashlib.sha256(sequence.encode()).hexdigest()!=row['sequence_sha256']:raise ValueError('Confirmation reference changed')
        label={'catalytic':'positive','decoy':'negative','unresolved':'indeterminate'}[row['kind']];source=dict(reference_id=row['reference_id'],accession=row['protein_accession'],scope=row['evidence_scope']);genera=[]
        if row['kind']=='unresolved':
            reviewed=identity_reviews[row['reference_id']]
            if reviewed['sequence_sha256']!=row['sequence_sha256'] or reviewed['protein_accession']!=row['protein_accession'] or reviewed['original_kind']!='unresolved':raise ValueError('Reference identity correspondence changed')
            label=reviewed['identity_label']
            if label=='positive':
                unit=byid[reviewed['source_unit_id']]
                if unit['split']!='reference_overlap' or unit['label']!='positive' or unit['sequence']!=sequence or source_family(unit)!=row['center']:raise ValueError('Reference identity source is not an exact exposed positive')
                source['identity_review']=reviewed;genera=unit['genus_taxids']
            elif label!='indeterminate':raise ValueError('Unsupported reference identity review label')
        add(sequence,label,row['center'],source,genera=genera)
    directory=root/'references/processed/catalytic_confirmation';directory.mkdir(parents=True,exist_ok=True)
    prototypes=sorted(by_sequence.values(),key=lambda r:r['prototype_id'])
    (directory/'prototypes.json').write_text(json.dumps(prototypes,indent=2)+'\n')
    SeqIO.write([SeqRecord(Seq(r['sequence']),id=r['prototype_id'],description=r['label']) for r in prototypes],directory/'prototypes.faa','fasta')
    source_policy_path=root/'config/reference_domain_policy.json';source_policy=json.loads(source_policy_path.read_text());models={};texts=[]
    for model in source_policy['models']:
        if model['accession'] not in policy['domain_profiles']:continue
        path=root/'references/raw/reference_domain_review'/(model['accession']+'.hmm')
        if sha(path)!=model['sha256']:raise ValueError('Confirmation profile changed')
        fields={line.split()[0]:' '.join(line.split()[1:]) for line in path.read_text().split('HMM          ')[0].splitlines()[1:] if line.strip()}
        if fields['NAME']!=model['profile'] or fields['ACC']!=model['versioned_accession']:raise ValueError('Unexpected confirmation profile')
        ga=fields['GA'].strip(';').split();models[model['profile']]=dict(model,accession=model['versioned_accession'],sequence_ga=float(ga[0]),domain_ga=float(ga[1]));texts.append(path.read_text())
    if {m['accession'].split('.')[0] for m in models.values()}!=set(policy['domain_profiles']):raise ValueError('Confirmation profile missing')
    (directory/'domains.hmm').write_text(''.join(texts));(directory/'domain_models.json').write_text(json.dumps(models,indent=2)+'\n')
    with (directory/'makeblastdb.log').open('w') as log:
        subprocess.run(['makeblastdb','-in',str(directory/'prototypes.faa'),'-dbtype','prot','-out',str(directory/'seed')],stdout=log,stderr=log,check=True)
    manifest=dict(status='built_before_calibration',policy_sha256=sha(policy_path),corpus_freeze_sha256=sha(corpus/'freeze.json'),training_audit_sha256=sha(audit_path),
        reference_domain_policy_sha256=sha(source_policy_path),reference_identity_review_sha256=sha(review_path),identity_adjudications_sha256=sha(root/policy['identity_adjudications']),references_sha256=sha(root/'references/processed/references.json'),
        prototype_count=len(prototypes),source_scope=policy['training_selection'],
        files={p.name:sha(p) for p in sorted(directory.iterdir()) if p.is_file() and p.name not in ['manifest.json','makeblastdb.log']})
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');return manifest


def core_similarity(hit,positions):
    """Compare both classes on the same observed query-core residues."""
    positions=set(positions);qpos=hit['qstart']-1;covered=set();matched=set()
    for q,s in zip(hit['qseq'],hit['sseq']):
        if q!='-':qpos+=1
        if q!='-' and s!='-' and qpos in positions:
            covered.add(qpos)
            if q==s and q not in 'XBZJO':matched.add(qpos)
    return dict(core_coverage=len(covered)/len(positions) if positions else 0,
        core_identity=100*len(matched)/len(covered) if covered else 0,
        similarity=100*len(matched)/len(positions) if positions else 0,
        reference_start=hit['sstart'],reference_end=hit['send'],query_start=hit['qstart'],query_end=hit['qend'],
        bitscore=hit['bitscore'],evalue=hit['evalue'],score_scope='exact_matches_per_query_core_not_combined_significance')


def excluded_prototype(prototype,query_metadata):
    return bool(set(prototype['homolog_components'])&{query_metadata.get('homolog_component','')} or
        set(prototype['genus_taxids'])&set(query_metadata.get('genus_taxids',[])))


def qualifies(feature,parameters):
    return bool(feature['eligible'] and feature['positive_prototype'] and
        feature['core_domain_coverage']>=parameters['min_domain_coverage'] and
        feature['positive_score']>=parameters['min_core_similarity'] and
        feature['competitor_margin']>0 and feature['competitor_margin']>=parameters['min_competitor_margin'])


def extract_features(annotation,output,root=ROOT,metadata=None,leave_out_development=False,threads=2,anchors=None):
    root=Path(root);annotation=Path(annotation);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    policy=json.loads((root/'config/catalytic_confirmation.json').read_text());directory=root/'references/processed/catalytic_confirmation'
    manifest=json.loads((directory/'manifest.json').read_text())
    if (manifest['policy_sha256']!=sha(root/'config/catalytic_confirmation.json') or
        manifest['identity_adjudications_sha256']!=sha(root/policy['identity_adjudications']) or
        manifest['reference_identity_review_sha256']!=sha(root/policy['reference_identity_review']) or
        manifest['reference_domain_policy_sha256']!=sha(root/'config/reference_domain_policy.json') or
        any(sha(directory/n)!=digest for n,digest in manifest['files'].items())):raise ValueError('Confirmation model changed')
    query=annotation/'proteins.faa';models=json.loads((directory/'domain_models.json').read_text());profile_policy=json.loads((root/'config/reference_domain_policy.json').read_text())
    normalized=output/'query_U_to_X.faa'
    SeqIO.write([SeqRecord(Seq(str(r.seq).replace('U','X')),id=r.id,description='HMM_only_normalization') for r in SeqIO.parse(query,'fasta')],normalized,'fasta')
    domains=scan_profile_domains(directory/'domains.hmm',normalized,models,0,output,threads=threads)
    raw=output/'prototype_blast.tsv'
    with (output/'prototype_blast.log').open('w') as log:
        subprocess.run(['blastp','-query',str(query),'-db',str(directory/'seed'),'-evalue',str(policy['comparison_max_evalue']),'-max_target_seqs','10000','-max_hsps','1',
            '-num_threads',str(threads),'-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(raw)],stderr=log,check=True)
    prototypes={r['prototype_id']:r for r in json.loads((directory/'prototypes.json').read_text())}
    hits=load_blast(raw,dict(blast_evalue=policy['comparison_max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0));by_query=defaultdict(list)
    for h in hits:by_query[h['qseqid']].append(h)
    site_models={m['reference_id']:m for m in json.loads((root/'references/processed/active_centers/models.json').read_text())}
    active_hits=load_blast(annotation/'active_center_blast.tsv',dict(blast_evalue=1,min_identity=0,min_query_coverage=0,min_reference_coverage=0));by_pair=defaultdict(list)
    for h in active_hits:by_pair[h['qseqid'],h['sseqid']].append(h)
    active={r['protein_id']:r for r in table(annotation/'active_center_summary.tsv')}
    anchors={r['protein_id']:r for r in (anchors if anchors is not None else table(annotation/'catalytic_anchors.tsv'))}
    genes={r['id']:r for r in table(annotation/'genes.tsv')};features=[]
    for record in SeqIO.parse(query,'fasta'):
        qid=record.id;summary=active.get(qid,{});anchor=anchors.get(qid,{});family=anchor.get('catalytic_family','');rid=summary.get('reference_id','');positions=set();reasons=[]
        model=site_models.get(rid)
        if model:
            for h in collinear_hits(by_pair[qid,rid]):
                positions.update(q for q,r in residue_states(h,range(model['window_start'],model['window_end']+1))[1].values() if q is not None and r!='-')
        if summary.get('active_center_status')!='template_ligand_set_supported' or summary.get('template_span_alignment_status')!='supported':reasons.append('reviewed_single_chain_center_or_span_not_supported')
        if genes[qid]['partial'] not in [False,'False']:reasons.append('partial_CDS')
        if family not in policy['eligible_calibration_families']:reasons.append('family_outside_calibration_scope')
        if not positions:reasons.append('no_reviewed_query_core')
        evidence=model_evidence(domains[qid],profile_policy)
        for e in evidence:
            for p in e['components']:p['protein_id']=qid
            segments=[(p['query_start'],p['query_end']) for p in e['components']]
            length=sum(b-a+1 for a,b in segments)
            e['query_core_overlap']=sum(any(a<=p<=b for a,b in segments) for p in positions)/length if length else 0
        coverage_by_accession={e['accession']:e['model_coverage'] for e in evidence if e['query_core_overlap']>=policy['domain_core_overlap_min']}
        domain_coverage=max((min(coverage_by_accession.get(a,0) for a in option) for option in policy['architecture'].get(family,[])),default=0)
        if not domain_coverage:reasons.append('required_domain_architecture_not_colocated_with_center')
        best_positive=None;best_competitor=None
        for h in by_query[qid] if positions else []:
            prototype=prototypes[h['sseqid']]
            if leave_out_development and excluded_prototype(prototype,(metadata or {}).get(qid,{})):continue
            score=core_similarity(h,positions)
            if score['core_coverage']<policy['core_comparison_min_coverage']:continue
            candidate=dict(score,prototype_id=prototype['prototype_id'],label=prototype['label'],family=prototype['family'])
            positive=prototype['label']=='positive' and prototype['family']==family
            previous=best_positive if positive else best_competitor
            if previous is None or (candidate['similarity'],candidate['prototype_id'])>(previous['similarity'],previous['prototype_id']):
                if positive:best_positive=candidate
                else:best_competitor=candidate
        ps=best_positive['similarity'] if best_positive else 0;cs=best_competitor['similarity'] if best_competitor else 0
        features.append(dict(protein_id=qid,family=family,eligible=not reasons,eligibility_reasons=reasons,core_positions=sorted(positions),core_domain_coverage=domain_coverage,
            positive_score=ps,competitor_score=cs,competitor_margin=ps-cs,positive_prototype=best_positive['prototype_id'] if best_positive else '',
            competitor_prototype=best_competitor['prototype_id'] if best_competitor else '',positive_evidence=best_positive,competitor_evidence=best_competitor,
            domain_evidence=evidence,source_active_center=summary))
    (output/'features.json').write_text(json.dumps(features,indent=2)+'\n')
    write_tsv(output/'features.tsv',[{k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) or v is None else v for k,v in r.items()} for r in features],FEATURE_FIELDS)
    provenance=dict(status='features_extracted',policy_sha256=sha(root/'config/catalytic_confirmation.json'),model_manifest_sha256=sha(directory/'manifest.json'),
        original_query_sha256=sha(query),active_summary_sha256=sha(annotation/'active_center_summary.tsv'),active_alignment_sha256=sha(annotation/'active_center_blast.tsv'),
        query_family_inputs={q:r.get('catalytic_family','') for q,r in anchors.items()},leave_out_development=leave_out_development,
        artifacts={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file() and p.name!='feature_manifest.json'})
    (output/'feature_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n');return features


def annotate_confirmation(output,annotations,threads=2,root=ROOT):
    root=Path(root);output=Path(output);directory=root/'references/processed/catalytic_confirmation'
    features=extract_features(output,output/'confirmation',root,threads=threads,anchors=annotations)
    rule_path=root/'config/catalytic_confirmation_thresholds.json';rules=json.loads(rule_path.read_text()) if rule_path.exists() else None
    byid={r['protein_id']:r for r in annotations};rows=[]
    for feature in features:
        parameters=rules['families'].get(feature['family']) if rules else None
        supported=bool(parameters and qualifies(feature,parameters))
        status='rule_supported' if supported else 'unresolved' if rules else 'calibration_pending'
        confidence='calibrated_high' if supported else 'not_high_confidence'
        row=dict(feature,confirmation_status=status,identity_confidence=confidence,rule_parameters=parameters,
            rule_scope='Catalytic identity under locked rule; no reaction, folding, cofactor occupancy or measured activity inference')
        rows.append(row);anchor=byid.get(feature['protein_id'])
        if anchor is None:continue
        anchor.update(identity_confirmation_status=status,catalytic_identity_confidence=confidence,identity_confirmation_family=feature['family'])
        if supported:
            anchor['catalytic_status']='supported'
            anchor['flags']=';'.join(filter(None,[anchor.get('flags',''),'identity_confirmation_separate_from_nearest_function_reference']))
        elif rules and feature['eligible'] and feature['competitor_prototype'] and feature['competitor_score']>=feature['positive_score']:
            anchor.update(catalytic_status='candidate',potential_function_ids='',function_status='unassigned',function_confidence='unassigned')
            if anchor['structure_confidence']=='high':anchor['structure_confidence']='medium'
            anchor['flags']=';'.join(filter(None,[anchor.get('flags',''),'confirmation_non_catalytic_indeterminate_or_other_family_competitor']))
    write_tsv(output/'catalytic_confirmation.tsv',[{k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) or v is None else v for k,v in r.items()} for r in rows],CONFIRMATION_FIELDS)
    return dict(model_manifest_sha256=sha(directory/'manifest.json'),policy_sha256=sha(root/'config/catalytic_confirmation.json'),
        script_sha256=sha(Path(__file__)),thresholds_sha256=sha(rule_path) if rules else None,
        feature_manifest_sha256=sha(output/'confirmation/feature_manifest.json'),high_identity_outputs=sum(r['identity_confidence']=='calibrated_high' for r in rows))


if __name__=='__main__':print(json.dumps(build_model(),indent=2))
