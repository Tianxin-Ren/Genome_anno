"""Uncalibrated family-sequence support; domain completeness and ligation stay separate."""
import hashlib
import json
from pathlib import Path
PROFILES={'NiFe':{'PF00374','PF00346'},'FeFe':{'PF02906'}}
REPLACED_REASONS={'reviewed_single_chain_center_or_span_not_supported','required_domain_architecture_not_colocated_with_center'}
FEFE_PROFILE_CORE_RULE='FeFe_complete_Hdomain_uncalibrated_v1'
NIFE_BACKBONE_RULE='NiFe_reference_backbone_and_ligand_contexts_uncalibrated_v1'
FEFE_PROFILE_LIGAND_RULE='FeFe-PF02906-profile-ligand-pattern-v1'
FEFE_EXPERIMENTAL_BACKBONE_RULE='FeFe-active-backbone-homology-pattern-v2'
NIFE_CORE_REGION_RULE='NiFe_catalytic_core_region_uncalibrated_v1'
NIFE_CORE_REGION_STATUS='sequence_supported_uncalibrated_core_region'
NIFE_INDETERMINATE_CORE_RULE='NiFe_intact_core_indeterminate_competitor_uncalibrated_v1'
FEFE_INDETERMINATE_CORE_RULE='FeFe_intact_core_indeterminate_competitor_uncalibrated_v1'
FEFE_INDETERMINATE_PROFILE_RULE='FeFe_complete_profile_ligands_indeterminate_competitor_uncalibrated_v1'


def _coding_support(coding_completeness, feature, gene=None):
    if not coding_completeness or not feature: return False
    pid=feature.get('protein_id')
    if isinstance(coding_completeness,dict) and coding_completeness.get('protein_id')==pid:
        row=coding_completeness
    else:
        row=coding_completeness.get(pid) if isinstance(coding_completeness,dict) else None
    if not isinstance(row,dict): return False
    from coding_completeness import supports
    expected=(gene or {}).get('protein_sha256') or feature.get('protein_sha256') or feature.get('query_protein_sha256')
    if not expected: return False
    return supports(row,pid,expected_protein_sha256=expected,feature=feature,gene=gene)


def _eligibility_reasons(feature,coding_completeness=None,gene=None):
    reasons=set(feature.get('eligibility_reasons',[]))
    if 'partial_CDS' in reasons and _coding_support(coding_completeness,feature,gene): reasons.remove('partial_CDS')
    return reasons


def reference_backbones(annotation,features,root):
    """Bind a winning NiFe alignment to actual full sequences, excluding alignment gaps."""
    from Bio import SeqIO
    from annotate import load_blast
    annotation,root=Path(annotation),Path(root)
    targets={q:f for q,f in features.items() if f['family']=='NiFe' and f.get('positive_prototype') and not local_domain_evidence(f)}
    if not targets:return {},{}
    paths=[annotation/'proteins.faa',annotation/'confirmation/prototype_blast.tsv',root/'references/processed/catalytic_confirmation/prototypes.faa']
    with paths[0].open() as handle:queries={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    with paths[2].open() as handle:references={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    result={}
    for h in load_blast(paths[1],dict(blast_evalue=1e-10,min_identity=0,min_query_coverage=0,min_reference_coverage=0)):
        f=targets.get(h['qseqid'])
        if not f or h['sseqid']!=f['positive_prototype']:continue
        match=f['positive_evidence']
        if (h['qstart'],h['qend'],h['sstart'],h['send'])!=(match['query_start'],match['query_end'],match['reference_start'],match['reference_end']):continue
        q,s=queries[h['qseqid']],references[h['sseqid']]
        if (h['qlen']!=len(q) or h['slen']!=len(s) or len(h['qseq'])!=len(h['sseq']) or
            h['qseq'].replace('-','')!=q[h['qstart']-1:h['qend']] or h['sseq'].replace('-','')!=s[h['sstart']-1:h['send']]):
            raise ValueError('NiFe backbone alignment does not match source sequences')
        if h['qseqid'] in result:raise ValueError('Ambiguous NiFe winning backbone alignment')
        paired=sum(a not in '-XBZJO' and b not in '-XBZJO' for a,b in zip(h['qseq'],h['sseq']))
        result[h['qseqid']]=dict(reference_id=h['sseqid'],paired_reference_coverage=paired/len(s),paired_residues=paired,
            reference_length=len(s),query_start=h['qstart'],query_end=h['qend'],evalue=h['evalue'],
            query_sha256=hashlib.sha256(q.encode()).hexdigest(),reference_sha256=hashlib.sha256(s.encode()).hexdigest())
    return result,{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def nife_backbone_support(feature,alternative,parameters,source_review=None,coding_completeness=None,gene=None):
    """Full reference backbone and local ligands can substitute for an unavailable profile."""
    if not feature or not parameters or feature['family']!='NiFe' or local_domain_evidence(feature):return False
    if _eligibility_reasons(feature,coding_completeness,gene)-REPLACED_REASONS or lineage_uncertain(feature,source_review):return False
    if alternative.get('query_partial') not in [False,'False'] and not _coding_support(coding_completeness,feature,gene):return False
    if not feature['positive_prototype'] or feature['positive_score']<parameters['min_core_similarity']:return False
    if feature['competitor_margin']<=0 or feature['competitor_margin']<parameters['min_competitor_margin']:return False
    backbone=alternative.get('NiFe_reference_backbone',{})
    if backbone.get('reference_id')!=feature['positive_prototype'] or backbone.get('paired_reference_coverage',0)<parameters['min_domain_coverage']:return False
    center=feature.get('source_active_center',{})
    contexts=center.get('site_contexts',[])
    if isinstance(contexts,str):contexts=json.loads(contexts)
    return bool(center.get('center_organization')=='single_polypeptide' and center.get('direct_ligand_status')=='present' and
        center.get('context_alignment_status')=='supported' and center.get('active_center_status')=='template_ligand_set_supported' and
        len(contexts)==4 and len({c['reference_position'] for c in contexts})==4)


def load_lineage_review(root):
    root=Path(root)
    return json.loads((root/'research_workflow/identity_scope_source_review.json').read_text())


def lineage_uncertain(feature,review):
    return bool(review and feature and feature.get('positive_prototype') in {r['prototype_id'] for r in review['uncertain_positive_prototypes']})


def local_domain_evidence(feature):
    accepted=[]
    for domain in feature.get('domain_evidence',[]):
        if domain['accession'] not in PROFILES.get(feature['family'],set()) or domain.get('query_core_overlap',0)<.7:continue
        parts=domain['components']
        if not parts or not all(p['sequence_bitscore']>=p['sequence_ga'] and p['domain_bitscore']>=p['domain_ga'] for p in parts):continue
        accepted.append(dict(accession=domain['accession'],profile=domain['profile'],model_coverage=domain['model_coverage'],
            full_domain_model_supported=domain['supported'],query_core_overlap=domain['query_core_overlap'],
            evidence='GA-qualified local profile evidence; not complete-domain or biochemical-activity proof'))
    return accepted


def nife_core_regions(annotation, features, root):
    """Recompute a hash-bound NiFe active-center region proof for partial CDS."""
    from Bio import SeqIO
    from annotate import load_blast
    from active_centers import collinear_hits, map_model
    def _validate_hit(hit, queries, references):
        if hit['qseqid'] not in queries or hit['sseqid'] not in references: raise ValueError('Unknown BLAST sequence ownership')
        if len(hit['qseq'])!=len(hit['sseq']): raise ValueError('BLAST alignment lengths differ')
        for field, sequence, lo, hi in [('qseq',queries[hit['qseqid']],hit['qstart'],hit['qend']),('sseq',references[hit['sseqid']],hit['sstart'],hit['send'])]:
            observed=hit[field].replace('-','').replace('U','X').replace('O','X')
            if not 1<=lo<=hi<=len(sequence) or observed!=sequence[lo-1:hi].replace('U','X').replace('O','X'):
                raise ValueError('BLAST alignment differs from bound '+field+' sequence')
    annotation, root = Path(annotation), Path(root)
    blast = annotation/'active_center_blast.tsv'
    models_path = root/'references/processed/active_centers/models.json'
    templates_path = root/'references/processed/active_centers/templates.faa'
    targets={pid:f for pid,f in features.items() if f.get('family')=='NiFe' and 'partial_CDS' in f.get('eligibility_reasons',[])}
    if not targets or not blast.is_file() or not models_path.is_file() or not templates_path.is_file(): return {}, {}
    models = {m['reference_id']: m for m in json.loads(models_path.read_text()) if m.get('family') == 'NiFe'}
    policy = json.loads((root/'references/processed/active_centers/manifest.json').read_text())['parameters']
    source_review = load_lineage_review(root)
    hits = load_blast(blast, {'blast_evalue': policy['max_evalue'], 'min_identity': 0, 'min_query_coverage': 0, 'min_reference_coverage': 0})
    with (annotation/'proteins.faa').open() as handle: proteins = {r.id: str(r.seq) for r in SeqIO.parse(handle, 'fasta') if r.id in targets}
    with templates_path.open() as handle: templates={r.id:str(r.seq) for r in SeqIO.parse(handle, 'fasta')}
    if any(models[r]['sequence'] != templates.get(r) for r in models): raise ValueError('Active-center template FASTA does not match models')
    by_pair={}
    for h in hits:
        if h['sseqid'] in models: by_pair.setdefault((h['qseqid'], h['sseqid']), []).append(h)
    proofs={}
    for pid, feature in targets.items():
        if feature.get('family') != 'NiFe' or 'partial_CDS' not in feature.get('eligibility_reasons',[]): continue
        if pid not in proteins or lineage_uncertain(feature, source_review) or set(feature.get('eligibility_reasons',[])) - {'partial_CDS'} - REPLACED_REASONS: continue
        if not local_domain_evidence(feature) or not feature.get('positive_prototype') or feature.get('positive_score', 0) < 0 or feature.get('competitor_margin', 0) <= 0: continue
        center=feature.get('source_active_center', {}); rid=center.get('reference_id') or feature.get('active_center_template'); model=models.get(rid); rows=by_pair.get((pid, rid), [])
        if not model or not rows: continue
        try:
            for h in rows: _validate_hit(h, {pid: proteins[pid]}, {rid: model['sequence']})
        except ValueError:
            continue
        components=collinear_hits(rows); mapped=map_model(components, model, policy, query_sequence=proteins[pid]); contexts=mapped.get('site_contexts', []); sites=mapped.get('sites', [])
        direct=[s for s in sites if s.get('residue_role') == 'direct_metal_ligand']
        complete=(mapped.get('window_coverage') == 1.0 and mapped.get('template_span_alignment_status') == 'supported' and mapped.get('context_alignment_status') == 'supported' and mapped.get('active_center_status') == 'template_ligand_set_supported' and mapped.get('direct_ligand_status') == 'present' and len(direct) == 4 and len({s.get('reference_position') for s in direct}) == 4 and all(s.get('status') == 'present' for s in direct) and len(contexts) == 4 and len({c.get('reference_position') for c in contexts}) == 4)
        if not complete: continue
        feature_digest=hashlib.sha256(json.dumps(feature,sort_keys=True,separators=(',',':')).encode()).hexdigest(); protein_digest=hashlib.sha256(proteins[pid].encode()).hexdigest()
        proofs[pid]=dict(status=NIFE_CORE_REGION_STATUS,rule_version=NIFE_CORE_REGION_RULE,protein_id=pid,protein_sha256=protein_digest,feature_sha256=feature_digest,reference_id=rid,reference_sha256=hashlib.sha256(model['sequence'].encode()).hexdigest(),coordinate_source_sha256=model['source_sha256'],window_start=model['window_start'],window_end=model['window_end'],hsp_coordinates=[{k:h[k] for k in ('qstart','qend','sstart','send','evalue','bitscore')} for h in components],recomputed_metrics={k:mapped[k] for k in ('window_coverage','window_identity','template_span_alignment_status','context_alignment_status','active_center_status','direct_ligand_status') }|{'center_organization':'single_polypeptide'},sites=sites,contexts=contexts,scope='NiFe catalytic active-center region only; incomplete CDS remains incomplete and no full-function claim')
    review_paths=[root/'research_workflow/identity_scope_source_review.json']
    review=json.loads(review_paths[0].read_text())
    review_paths += [root/review[k] for k in ('source_path','prototype_path')]
    deps=[blast,models_path,templates_path,annotation/'proteins.faa',root/'references/processed/active_centers/manifest.json',root/'research_workflow/identity_scope.py',root/'scripts/active_centers.py',root/'scripts/annotate.py',root/'research_workflow/mbh_function.py',*review_paths]
    return proofs, {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in deps if p.is_file()}


def nife_core_region_support(feature, alternative, parameters, source_review=None, gene=None):
    """Validate a source-bound partial-CDS catalytic-region proof."""
    if not feature or not parameters or feature.get('family') != 'NiFe' or 'partial_CDS' not in feature.get('eligibility_reasons',[]): return False
    if set(feature.get('eligibility_reasons',[])) - {'partial_CDS'} - REPLACED_REASONS or lineage_uncertain(feature, source_review): return False
    if not local_domain_evidence(feature) or not feature.get('positive_prototype') or feature.get('positive_score', 0) < parameters.get('min_core_similarity', 0) or feature.get('competitor_margin', 0) <= 0 or feature.get('competitor_margin', 0) < parameters.get('min_competitor_margin', 0): return False
    proof=(alternative or {}).get('NiFe_core_region', {})
    if not isinstance(proof,dict) or proof.get('status') != NIFE_CORE_REGION_STATUS or proof.get('rule_version') != NIFE_CORE_REGION_RULE or proof.get('protein_id') != feature.get('protein_id'): return False
    expected=(gene or {}).get('protein_sha256') or feature.get('protein_sha256') or feature.get('query_protein_sha256')
    if not expected or proof.get('protein_sha256') != expected: return False
    digest=hashlib.sha256(json.dumps(feature,sort_keys=True,separators=(',',':')).encode()).hexdigest(); center=feature.get('source_active_center',{})
    if proof.get('feature_sha256') != digest or proof.get('reference_id') != (center.get('reference_id') or feature.get('active_center_template')): return False
    if center.get('center_organization') != 'single_polypeptide' or center.get('active_center_status') != 'template_ligand_set_supported' or center.get('template_span_alignment_status') != 'supported' or center.get('context_alignment_status') != 'supported' or center.get('direct_ligand_status') != 'present' or float(center.get('window_coverage', 0)) != 1.0: return False
    metrics=proof.get('recomputed_metrics',{}); sites=proof.get('sites',[]); contexts=proof.get('contexts',[]); direct=[s for s in sites if s.get('residue_role')=='direct_metal_ligand']
    try: feature_contexts=center.get('site_contexts',[]); feature_contexts=json.loads(feature_contexts) if isinstance(feature_contexts,str) else feature_contexts
    except (TypeError,ValueError): return False
    same_contexts=len(feature_contexts)==len(contexts)==4 and all(a.get('reference_position')==b.get('reference_position') and a.get('coverage')==b.get('coverage') and a.get('identity')==b.get('identity') for a,b in zip(feature_contexts,contexts))
    return bool(metrics.get('center_organization') == center.get('center_organization') == 'single_polypeptide' and metrics.get('window_coverage') == 1.0 and float(center.get('window_coverage',0)) == metrics.get('window_coverage') and float(metrics.get('window_identity',-1)) == float(center.get('window_identity',-2)) and metrics.get('template_span_alignment_status') == center.get('template_span_alignment_status') and metrics.get('context_alignment_status') == center.get('context_alignment_status') and metrics.get('active_center_status') == center.get('active_center_status') and metrics.get('direct_ligand_status') == center.get('direct_ligand_status') and same_contexts and proof.get('hsp_coordinates') and len(direct)==4 and len({s.get('reference_position') for s in direct})==4 and all(s.get('status')=='present' and s.get('alignment_components') for s in direct))


def complete_fefe_core_support(feature,alternative,parameters,source_review=None,coding_completeness=None,gene=None):
    """A complete H-domain and ligand pattern support family potential, not atomic ligation."""
    if not feature or not parameters or feature['family']!='FeFe':return False
    if _eligibility_reasons(feature,coding_completeness,gene)-REPLACED_REASONS or lineage_uncertain(feature,source_review):return False
    if alternative.get('query_partial') not in [False,'False'] and not _coding_support(coding_completeness,feature,gene):return False
    if not feature['positive_prototype'] or feature['positive_score']<=0 or feature['competitor_margin']<=0:return False
    match=feature.get('positive_evidence',{})
    if match.get('core_coverage',0)<.7 or match.get('evalue',1)>1e-10:return False
    center=feature['source_active_center']
    if center.get('center_organization','single_polypeptide')!='single_polypeptide' or center.get('direct_ligand_status')!='present':return False
    return any(d['accession']=='PF02906' and d.get('supported') is True and d.get('evidence_kind')=='single_domain' and
        d['model_coverage']>=parameters['min_domain_coverage'] and d.get('query_core_overlap',0)>=.7 and len(d['components'])==1 and
        all(p['sequence_bitscore']>=p['sequence_ga'] and p['domain_bitscore']>=p['domain_ga'] for p in d['components']) for d in feature.get('domain_evidence',[]))


def fefe_profile_ligand_support(feature,alternative,parameters,source_review=None,coding_completeness=None,gene=None):
    """A complete profile can resolve an uncovered sequence pattern, never a conflict."""
    if not feature or not parameters or feature['family']!='FeFe':return False
    evidence=alternative.get('FeFe_profile_ligands',{})
    if evidence.get('profile_ligand_status')!='profile_mapped_ligand_pattern_supported':return False
    if _eligibility_reasons(feature,coding_completeness,gene)-REPLACED_REASONS or lineage_uncertain(feature,source_review):return False
    if alternative.get('query_partial') not in [False,'False'] and not _coding_support(coding_completeness,feature,gene):return False
    if not feature.get('positive_prototype') or feature['positive_score']<parameters['min_core_similarity']:return False
    if feature['competitor_margin']<=0 or feature['competitor_margin']<parameters['min_competitor_margin']:return False
    match=feature.get('positive_evidence',{})
    if match.get('core_coverage',0)<.7 or match.get('evalue',1)>1e-10:return False
    center=feature.get('source_active_center',{})
    if center.get('center_organization','single_polypeptide')!='single_polypeptide' or center.get('direct_ligand_status')!='not_covered':return False
    return _fefe_profile_ligand_pattern_supported(feature,alternative,parameters)


def _fefe_profile_ligand_pattern_supported(feature,alternative,parameters):
    """Validate the current mapped four-Cys proof independently of atomic span."""
    evidence=alternative.get('FeFe_profile_ligands',{})
    if evidence.get('profile_ligand_status')!='profile_mapped_ligand_pattern_supported':return False
    from fefe_profile_ligands import _digest,_load_model,_query_domain
    if evidence.get('rule_version')!=FEFE_PROFILE_LIGAND_RULE or evidence.get('feature_sha256')!=_digest(feature):return False
    model,_,_,_=_load_model(Path(__file__).resolve().parents[1])
    if evidence.get('hmm_sha256')!=model['hmm_sha256'] or evidence.get('reference_model_sha256')!=model['reference_model_sha256']:return False
    domains=_query_domain(feature,model,parameters['min_domain_coverage'],feature.get('protein_id'),evidence.get('query_length'))
    if len(domains)!=1:return False
    lo,hi,domain=domains[0]
    if evidence.get('domain_ownership')!=dict(accession='PF02906',query_start=lo,query_end=hi,domain_sha256=_digest(domain)):return False
    sites=evidence.get('sites',[])
    if len(sites)!=4 or [(s.get('reference_position'),s.get('match_column')) for s in sites]!=list(zip(model['reference_positions'],model['expected_match_columns'])):return False
    positions=[s.get('query_position') for s in sites]
    return bool(evidence.get('protein_id')==feature.get('protein_id') and not evidence.get('blockers') and
        all(s.get('residue_status')=='present' and s.get('residue')=='C' and s.get('pp')=='*' and s.get('within_domain') is True and
            isinstance(s.get('query_position'),int) and lo<=s['query_position']<=hi for s in sites) and
        len(set(positions))==4 and positions==sorted(positions))


def fefe_indeterminate_profile_core_support(feature,alternative,parameters,source_review=None,coding_completeness=None,gene=None):
    """Complete H-domain and bound ligands support potential despite unknown reference truth."""
    if not feature or not parameters or feature.get('family')!='FeFe' or feature.get('competitor_margin',1)>0:return False
    if _eligibility_reasons(feature,coding_completeness,gene)-REPLACED_REASONS or lineage_uncertain(feature,source_review):return False
    if alternative.get('query_partial') not in [False,'False'] and not _coding_support(coding_completeness,feature,gene):return False
    positive=feature.get('positive_evidence') or {};competitor=feature.get('competitor_evidence') or {}
    if positive.get('label')!='positive' or positive.get('family')!='FeFe' or competitor.get('label')!='indeterminate' or competitor.get('family')!='FeFe':return False
    if not feature.get('positive_prototype') or positive.get('prototype_id')!=feature['positive_prototype'] or not feature.get('competitor_prototype') or competitor.get('prototype_id')!=feature['competitor_prototype']:return False
    if feature.get('positive_score',0)<parameters['min_core_similarity'] or positive.get('core_coverage',0)<parameters['min_domain_coverage'] or positive.get('evalue',1)>1e-10:return False
    center=feature.get('source_active_center',{})
    if center.get('center_organization','single_polypeptide')!='single_polypeptide' or center.get('direct_ligand_status') not in {'present','not_covered'}:return False
    evidence=alternative.get('FeFe_profile_ligands',{})
    expected=(gene or {}).get('protein_sha256') or feature.get('protein_sha256') or feature.get('query_protein_sha256')
    if not expected or evidence.get('query_sha256')!=expected:return False
    return _fefe_profile_ligand_pattern_supported(feature,alternative,parameters)


def fefe_experimental_backbone_support(feature,alternative,parameters,source_review=None,coding_completeness=None,gene=None):
    """A complete measured-active backbone can bridge a partial profile, never a ligand conflict."""
    if not feature or not parameters or feature.get('family')!='FeFe':return False
    if _eligibility_reasons(feature,coding_completeness,gene)-REPLACED_REASONS or lineage_uncertain(feature,source_review):return False
    if alternative.get('query_partial') not in [False,'False'] and not _coding_support(coding_completeness,feature,gene):return False
    if not feature.get('positive_prototype') or feature['positive_score']<parameters['min_core_similarity'] or feature['competitor_margin']<=0 or feature['competitor_margin']<parameters['min_competitor_margin']:return False
    center=feature.get('source_active_center',{})
    if center.get('center_organization','single_polypeptide')!='single_polypeptide' or center.get('direct_ligand_status') not in ['not_covered','present']:return False
    from fefe_experimental_backbone import record_supported
    return record_supported(alternative.get('FeFe_experimental_backbone',{}),feature,parameters)


def intact_core_indeterminate_competitor_support(feature,parameters):
    """Unknown same-family reference truth does not negate an intact catalytic core."""
    if not feature or not parameters or feature.get('family') not in {'NiFe','FeFe'} or feature.get('competitor_margin',1)>0:return False
    family=feature['family'];profile={'NiFe':'PF00374','FeFe':'PF02906'}[family]
    positive=feature.get('positive_evidence') or {};competitor=feature.get('competitor_evidence') or {}
    if positive.get('label')!='positive' or positive.get('family')!=family or competitor.get('label')!='indeterminate' or competitor.get('family')!=family:return False
    if not feature.get('positive_prototype') or positive.get('prototype_id')!=feature['positive_prototype'] or not feature.get('competitor_prototype') or competitor.get('prototype_id')!=feature['competitor_prototype']:return False
    coverage=parameters['min_domain_coverage']
    if feature.get('positive_score',0)<parameters['min_core_similarity'] or feature.get('core_domain_coverage',0)<coverage or positive.get('core_coverage',0)<coverage or positive.get('evalue',1)>1e-10:return False
    if not any(d['accession']==profile and d['full_domain_model_supported'] and d['model_coverage']>=coverage for d in local_domain_evidence(feature)):return False
    center=feature.get('source_active_center',{})
    contexts=center.get('site_contexts',[]);contexts=json.loads(contexts) if isinstance(contexts,str) else contexts
    return bool(center.get('center_organization')=='single_polypeptide' and center.get('active_center_status')=='template_ligand_set_supported' and
        center.get('direct_ligand_status')=='present' and center.get('template_span_alignment_status')=='supported' and center.get('context_alignment_status')=='supported' and
        len(contexts)==4 and len({c.get('reference_position') for c in contexts})==4)


def sequence_core_support(feature,alternative,parameters,keep_center_context=False,source_review=None,coding_completeness=None,gene=None):
    if not feature or not parameters or feature['family'] not in PROFILES:return False
    if _eligibility_reasons(feature,coding_completeness,gene)-REPLACED_REASONS:return False
    if lineage_uncertain(feature,source_review):return False
    if feature['family']=='FeFe' and not keep_center_context:
        from fefe_group_b_evidence import record_supported
        from fefe_profile_ligands import _digest
        proof=alternative.get('FeFe_profile_ligands',{})
        expected=(gene or {}).get('protein_sha256') or feature.get('protein_sha256') or feature.get('query_protein_sha256')
        if (expected and proof.get('query_sha256')==expected and proof.get('feature_sha256')==_digest(feature)
                and proof.get('protein_id')==feature.get('protein_id')
                and (alternative.get('query_partial') in [False,'False'] or _coding_support(coding_completeness,feature,gene))
                and record_supported(feature,proof,parameters)):
            return True
    scored=feature['positive_prototype'] and feature['positive_score']>=parameters['min_core_similarity'] and feature['competitor_margin']>0 and feature['competitor_margin']>=parameters['min_competitor_margin']
    profile_indeterminate=fefe_indeterminate_profile_core_support(feature,alternative,parameters,source_review,coding_completeness,gene)
    if not scored and not intact_core_indeterminate_competitor_support(feature,parameters) and not complete_fefe_core_support(feature,alternative,parameters,source_review,coding_completeness,gene) and not profile_indeterminate:return False
    backbone=fefe_experimental_backbone_support(feature,alternative,parameters,source_review,coding_completeness,gene)
    if not local_domain_evidence(feature) and not nife_backbone_support(feature,alternative,parameters,source_review,coding_completeness,gene) and not backbone:return False
    center=feature['source_active_center']
    if center.get('center_organization','single_polypeptide')!='single_polypeptide':return False
    if keep_center_context:return center.get('active_center_status')=='template_ligand_set_supported' and center.get('template_span_alignment_status')=='supported'
    direct=center.get('direct_ligand_status')=='present'
    pattern=(alternative.get('candidate_model_id')=='Aw_H_cubane_3C1Ser_literature' and
        alternative.get('pattern_status')=='literature_ligand_pattern_candidate' and alternative.get('query_partial') in [False,'False'])
    return bool(direct or pattern or fefe_profile_ligand_support(feature,alternative,parameters,source_review,coding_completeness,gene) or backbone or profile_indeterminate)
