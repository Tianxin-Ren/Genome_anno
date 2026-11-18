"""Export claim-specific evidence and loss reasons without changing frozen predictions."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from identity_scope import sequence_core_support,local_domain_evidence,load_lineage_review,lineage_uncertain,complete_fefe_core_support,FEFE_PROFILE_CORE_RULE,reference_backbones,nife_backbone_support,NIFE_BACKBONE_RULE,fefe_profile_ligand_support,FEFE_PROFILE_LIGAND_RULE,fefe_experimental_backbone_support,FEFE_EXPERIMENTAL_BACKBONE_RULE
from identity_scope import nife_core_regions,nife_core_region_support,NIFE_CORE_REGION_STATUS,NIFE_CORE_REGION_RULE
from identity_scope import intact_core_indeterminate_competitor_support,NIFE_INDETERMINATE_CORE_RULE,FEFE_INDETERMINATE_CORE_RULE
from identity_scope import fefe_indeterminate_profile_core_support,FEFE_INDETERMINATE_PROFILE_RULE
TYPES={'observed_conflict','sequence_incomplete','partner_not_detected','partner_ambiguous',
       'outside_calibration_scope','reference_rule_missing','truth_unavailable','not_tested'}
FIELDS=('candidate_id family discovery_pass identity_status original_identity_status identity_sequence_evidence center_model_status query_coordination_status module_graph_status '
        'complex_status reaction_claim_status first_blocking_reason all_blocking_reasons missing_evidence_type calibration_scope '
        'evaluation_truth_status source_annotation reference_reaction_status').split()
ATTR_FIELDS='candidate_id attribute status value evidence_dependencies evidence_source limitations'.split()
SUPPORTED_ATTRIBUTE_STATES={
    'reference_like_modules_detected','domain_profile_supported','profile_in_candidate_architecture',
    'locally_unambiguous_sequence_candidate','reference_like_components_cooccur','nonoverlapping_reference_domains_supported',
    'sequence_supported_candidate','local_coupling_potential_uncalibrated','reference_like_H2_Hdr_coupling_hypothesis',
    'compatible_reaction_candidate_uncalibrated','reference_like_module_uncalibrated','existing_reaction_potential',
    'reference_like_FeS_relay_candidate','reference_like_H2_Fd_coupling_potential','reference_like_nonredox_scaffold_candidate',
    'local_membrane_H2_redox_potential_uncalibrated',
    'reference_like_quinone_site','reference_like_membrane_delivery_module','local_quinone_coupling_potential_uncalibrated',
    'candidate_H2_oxidation_potential',
    'assay_supported_candidate','sequence_supported_fusion',
    'candidate_CO_linked_H2_pathway','candidate_H2_formate_coupling','candidate_genome_H2_formate_pathway',
    'candidate_H2_formate_pathway_conditional_on_recoding',
    'sequence_supported_core_region',
    'supported_role_module_homology',
    'exact_reference_sequence_potential','exact_reference_coding_locus_potential','inherited_reference_potential',
    'complete_supported','core_complete_accessory_missing','composition_supported_reorganized',
    'coding_locus_supported_conditional','candidate_H2_formate_pathway_conditional_on_recoding'}

from preferred_prediction import PREDICTION_FIELDS, VERSION as PREFERRED_POLICY_VERSION, predict as preferred_predict
from reference_reactions import TUPLE_POTENTIAL_STATUSES, coding_locus_candidate_ids
from split_core_identity import SPLIT_IDENTITY, SPLIT_LABEL, prepare as prepare_split_identity
CODING_LOCUS_STATUS='exact_reference_coding_locus_potential'


def load(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def table(path):
    if not path.exists():return []
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))
def yes(value):return value is True or value=='True'
def write(path,rows,fields):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader()
        writer.writerows({k:json.dumps(r.get(k),ensure_ascii=False) if isinstance(r.get(k),(list,dict)) else r.get(k,'') for k in fields} for r in rows)



def preferred_role_options(nodes, anchor_id, structure_id):
    """Keep partner alternatives without copying residue-level search matrices."""
    options=[]
    for node in nodes:
        if node.get('anchor')!=anchor_id or node.get('structure_id')!=structure_id:continue
        raw=node.get('candidate_realizations',[])
        if isinstance(raw,str):raw=json.loads(raw or '[]')
        candidates=[{k:c[k] for k in ('protein_ids','coverage','identity','fragmented','jointly_feasible','competitive','selected') if k in c} for c in raw]
        options.append(dict(role=node['role'],selected_protein_ids=node.get('protein_ids','').split(';'),candidates=candidates))
    return options


def identity_reasons(feature,parameters):
    reasons=[]
    def add(code,kind,detail):reasons.append(dict(stage='identity',claim='calibrated_high_identity',code=code,type=kind,detail=detail,source='confirmation/features.json'))
    if not feature:add('identity_features_missing','not_tested','No identity features available');return reasons
    center=feature.get('source_active_center',{});status=center.get('active_center_status','')
    for reason in feature['eligibility_reasons']:
        kind='sequence_incomplete' if reason=='partial_CDS' else 'reference_rule_missing' if reason=='no_reviewed_query_core' else 'outside_calibration_scope'
        if reason=='reviewed_single_chain_center_or_span_not_supported' and status=='ligands_substituted':kind='observed_conflict'
        add(reason,kind,dict(center_status=status,scope='Conflict with canonical template is not proof of inactivity') if reason.startswith('reviewed_single') else reason)
    if not parameters:add('family_not_calibrated','outside_calibration_scope',feature['family']);return reasons
    for key,threshold,code in [('core_domain_coverage','min_domain_coverage','domain_coverage_below_rule'),('positive_score','min_core_similarity','core_similarity_below_rule')]:
        if feature[key]<parameters[threshold]:add(code,'outside_calibration_scope',dict(observed=feature[key],required=parameters[threshold]))
    if not feature['positive_prototype']:add('positive_prototype_missing','reference_rule_missing','No comparable positive prototype')
    margin=feature['competitor_margin']
    if margin<=0:
        competitor_evidence=feature.get('competitor_evidence') or {}
        label=competitor_evidence.get('label')
        competitor=feature.get('competitor_prototype') or ''
        if not competitor and not competitor_evidence:
            kind='reference_rule_missing'
            scope='No comparable competitor reference; default zero margin is not an observed conflict'
        else:
            kind='truth_unavailable' if label=='indeterminate' else 'observed_conflict'
            scope='Unresolved reference truth; retained calibration margin is not biochemical evidence of inactivity' if label=='indeterminate' else 'Competing reference evidence'
        add('nonpositive_competitor_margin',kind,
            dict(margin=margin,competitor=feature['competitor_prototype'],competitor_label=label,
                 scope=scope))
    elif margin<parameters['min_competitor_margin']:add('competitor_margin_below_rule','outside_calibration_scope',dict(margin=margin,required=parameters['min_competitor_margin']))
    return reasons


def sequence_status(anchor,feature,alternative,parameters,source_review=None,coding_completeness=None,gene=None):
    original=anchor.get('catalytic_identity_confidence','not_high_confidence')
    if original=='calibrated_high':return original
    if feature and feature['family']=='Fe':
        evidence=alternative.get('Fe_reference_identity',{})
        return 'sequence_supported_uncalibrated_core' if ((evidence.get('canonical') or {}).get('core_status')=='reference_supported_uncalibrated' or (evidence.get('canonical') or {}).get('branch_core_status')=='sequence_supported_candidate' or (evidence.get('hmdii') or {}).get('hmdii_status')=='reference_supported_uncalibrated') else 'evidence_insufficient'
    if not feature or not parameters:return 'evidence_insufficient'
    if feature['family']=='FeFe' and alternative.get('Assay_sequence_identity'):
        from assay_sequence_identity import record_supported
        if record_supported(alternative['Assay_sequence_identity'],feature):return 'sequence_supported_uncalibrated_core'
    coding_completeness = coding_completeness or alternative.get('Coding_completeness')
    if sequence_core_support(feature,alternative,parameters,coding_completeness=coding_completeness,gene=gene):
        if lineage_uncertain(feature,source_review):return 'hydrogenase_like_lineage_unconfirmed'
        if alternative.get('pattern_status')=='literature_ligand_pattern_candidate':return 'sequence_supported_uncalibrated_noncanonical'
        return 'sequence_supported_uncalibrated_core'
    if nife_core_region_support(feature,alternative,parameters,source_review,gene=gene):
        return NIFE_CORE_REGION_STATUS
    from coding_completeness import supports as coding_supports
    expected=(gene or {}).get('protein_sha256') or feature.get('protein_sha256') or feature.get('query_protein_sha256')
    coding_complete=bool(expected and coding_supports(coding_completeness,feature.get('protein_id'),expected_protein_sha256=expected,feature=feature,gene=gene))
    scoped_reasons=set(feature['eligibility_reasons'])-({'partial_CDS'} if coding_complete else set())
    only_center=scoped_reasons<={'reviewed_single_chain_center_or_span_not_supported'}
    score_ok=(feature['positive_prototype'] and feature['core_domain_coverage']>=parameters['min_domain_coverage'] and
        feature['positive_score']>=parameters['min_core_similarity'] and feature['competitor_margin']>0 and feature['competitor_margin']>=parameters['min_competitor_margin'])
    if only_center and score_ok:
        if (alternative.get('candidate_model_id')=='Aw_H_cubane_3C1Ser_literature' and
            alternative.get('pattern_status')=='literature_ligand_pattern_candidate' and (alternative.get('query_partial') in [False,'False'] or coding_complete)):
            return 'sequence_supported_uncalibrated_noncanonical'
        return 'sequence_support_center_model_unresolved'
    return 'evidence_insufficient'


def _candidate_context(annotation, screen_dir, confirmation_dir):
    """Load a hash-bound, in-memory extension; never edits core tables."""
    screen_dir, confirmation_dir = Path(screen_dir).resolve(), Path(confirmation_dir).resolve()
    smp = screen_dir/'candidate_screen_manifest.json'; cmp = confirmation_dir/'candidate_confirmation_manifest.json'
    if not smp.is_file() or not cmp.is_file(): raise ValueError('Candidate screen/confirmation manifest missing')
    sm = load(smp); cm = load(cmp)
    if sm.get('status') != 'complete' or cm.get('status') != 'complete': raise ValueError('Candidate screen/confirmation incomplete')
    if str(annotation) != str(Path(sm.get('annotation','')).resolve()): raise ValueError('Candidate screen annotation mismatch')
    for name,digest in sm.get('output_sha256',{}).items():
        if not (screen_dir/name).is_file() or sha(screen_dir/name) != digest: raise ValueError('Candidate screen output missing or changed: '+name)
    if sha(screen_dir/'candidate_screen.tsv') != sm.get('output_sha256',{}).get('candidate_screen.tsv'):
        raise ValueError('Candidate screen table changed')
    report_path = confirmation_dir/'candidate_confirmation.json'
    if sha(report_path) != cm.get('report_sha256'): raise ValueError('Candidate confirmation report changed')
    report = load(report_path)
    if report.get('source_screen_manifest_sha256') != sha(smp) or str(Path(report.get('annotation','')).resolve()) != str(annotation):
        raise ValueError('Candidate confirmation source mismatch')
    if report.get('source_screen_inputs') != sm.get('input_sha256'): raise ValueError('Confirmation does not bind the complete screen input set')
    for path, digest in sm.get('input_sha256',{}).items():
        p=Path(path)
        if not p.exists() or sha(p) != digest: raise ValueError('Candidate source input missing or changed: '+path)
    for group in ['code_sha256','model_sha256']:
        if not report.get(group): raise ValueError('Candidate confirmation lacks '+group)
        for path,digest in report.get(group,{}).items():
            p=Path(path)
            if not p.is_file() or sha(p) != digest: raise ValueError('Candidate confirmation '+group+' missing or changed: '+path)
    for path,digest in report.get('backbone_input_sha256',{}).items():
        if not Path(path).is_file() or sha(Path(path))!=digest:raise ValueError('Candidate backbone source changed')
    for path,digest in report.get('Fe_input_sha256',{}).items():
        if not Path(path).is_file() or sha(Path(path))!=digest:raise ValueError('Candidate Fe identity source changed')
    for path,digest in report.get('FeFe_profile_input_sha256',{}).items():
        if not Path(path).is_file() or sha(Path(path))!=digest:raise ValueError('Candidate FeFe profile source changed')
    for path,digest in report.get('FeFe_experimental_backbone_input_sha256',{}).items():
        if not Path(path).is_file() or sha(Path(path))!=digest:raise ValueError('Candidate FeFe experimental backbone source changed')
    coding_rows={}; coding_inputs=report.get('coding_completeness_input_sha256',{})
    if coding_inputs:
        from coding_completeness import load_evidence as load_coding
        coding_dir=confirmation_dir/'coding_completeness'
        coding_manifest,coding_rows=load_coding(coding_dir,annotation=annotation,root=ROOT)
        expected_coding_inputs=dict(coding_manifest['input_hashes'])
        expected_coding_inputs.update({str(p):sha(p) for p in [coding_dir/'coding_completeness.json',coding_dir/'coding_completeness_manifest.json']})
        if coding_inputs!=expected_coding_inputs:raise ValueError('Candidate coding completeness binding is incomplete or changed')
    screened={r['protein_id'] for r in table(screen_dir/'candidate_screen.tsv')}
    core={r['protein_id'] for r in table(annotation/'catalytic_anchors.tsv')}
    decisions=report.get('decisions',[]); ids=[d.get('candidate_id') for d in decisions]
    if len(ids) != len(set(ids)) or set(ids)-screened or set(ids)&core: raise ValueError('Candidate IDs are not unique or not extra screen candidates')
    features={}; alternatives={}; accepted=[]; bound=[smp,screen_dir/'candidate_screen.tsv',cmp,report_path,annotation/'genes.tsv',annotation/'proteins.faa']
    bound.extend(Path(p) for p in sm['input_sha256'])
    bound.extend(screen_dir/p for p in sm['output_sha256'])
    bound.extend(Path(p) for group in ['code_sha256','model_sha256'] for p in report[group])
    bound.extend(Path(p) for p in report.get('backbone_input_sha256',{}))
    bound.extend(Path(p) for p in report.get('Fe_input_sha256',{}))
    bound.extend(Path(p) for p in report.get('FeFe_profile_input_sha256',{}))
    bound.extend(Path(p) for p in report.get('FeFe_experimental_backbone_input_sha256',{}))
    bound.extend(Path(p) for p in coding_inputs)
    profile_evidence={}
    if report.get('FeFe_profile_input_sha256'):
        from fefe_profile_ligands import verify_saved
        subset=confirmation_dir/'candidate_annotation'
        profile_features={f['protein_id']:f for f in load(subset/'confirmation/features.json')}
        profile_parameters=load(ROOT/'config/catalytic_confirmation_thresholds.json')['families']['FeFe']
        profile_evidence,_=verify_saved(confirmation_dir/'fefe_profile_ligands',annotation=subset,features=profile_features,parameters=profile_parameters,root=ROOT)
    experimental_backbone_evidence={}
    if report.get('FeFe_experimental_backbone_input_sha256'):
        from fefe_experimental_backbone import verify_saved
        subset=confirmation_dir/'candidate_annotation'
        backbone_features={f['protein_id']:f for f in load(subset/'confirmation/features.json')}
        backbone_parameters=load(ROOT/'config/catalytic_confirmation_thresholds.json')['families']['FeFe']
        experimental_backbone_evidence,_=verify_saved(confirmation_dir/'fefe_experimental_backbone',annotation=subset,features=backbone_features,parameters=backbone_parameters,root=ROOT)
    assay_evidence={};assay_inputs=report.get('Assay_sequence_identity_input_sha256',{})
    if assay_inputs:
        from assay_sequence_identity import verify_saved
        subset=confirmation_dir/'candidate_annotation';assay_dir=confirmation_dir/'assay_sequence_identity'
        assay_features={f['protein_id']:f for f in load(subset/'confirmation/features.json')}
        assay_evidence=verify_saved(assay_dir,subset,assay_features,ROOT)
        expected=dict(load(assay_dir/'assay_sequence_identity_manifest.json')['input_hashes'])
        expected.update({str(p.resolve()):sha(p) for p in [assay_dir/'assay_sequence_identity.json',assay_dir/'assay_sequence_identity_manifest.json']})
        if assay_inputs!=expected:raise ValueError('Candidate assayed sequence binding is incomplete or changed')
        bound.extend(Path(p) for p in assay_inputs)
    allowed={'sequence_supported_uncalibrated_core','sequence_supported_uncalibrated_noncanonical'}
    for d in decisions:
        status=d.get('sequence_status','unconfirmed')
        fam=(d.get('family') or '').strip(); normalized=d.get('normalized_families',[])
        if status in allowed and not d.get('blocking_reasons') and fam in {'NiFe','FeFe','Fe'} and normalized == [fam]: accepted.append(d)
    source_gene_rows=table(annotation/'genes.tsv'); source_genes={r['id']:r for r in source_gene_rows}
    if len(source_genes)!=len(source_gene_rows): raise ValueError('Duplicate source gene IDs')
    with (annotation/'proteins.faa').open() as handle: source_records=list(SeqIO.parse(handle,'fasta'))
    if len({r.id for r in source_records})!=len(source_records): raise ValueError('Duplicate source protein IDs')
    source_proteins={r.id:str(r.seq) for r in source_records}
    for d in decisions:
        if d not in accepted: continue
        pid=d['candidate_id']
        ev=report.get('evidence',{}).get(d.get('candidate_id'),{})
        fp=Path(ev.get('features',''))
        if not fp.is_file(): raise ValueError('Confirmed candidate feature file missing')
        feature_rows=load(fp); feature_rows=feature_rows if isinstance(feature_rows,list) else feature_rows.get('features',[])
        if len({x.get('protein_id') for x in feature_rows}) != len(feature_rows): raise ValueError('Duplicate candidate feature IDs')
        feature=next((x for x in feature_rows if x.get('protein_id')==d['candidate_id']),None)
        if feature is None or d.get('features') != feature: raise ValueError('Confirmed candidate feature evidence changed')
        if feature.get('family')!=d['family']: raise ValueError('Confirmed candidate feature family differs from its decision')
        features[d['candidate_id']]=feature; bound.append(fp)
        if d.get('FeFe_profile_ligands')!=profile_evidence.get(pid):raise ValueError('Candidate FeFe profile evidence changed')
        if pid in profile_evidence:alternatives.setdefault(pid,{})['FeFe_profile_ligands']=profile_evidence[pid]
        if d.get('FeFe_experimental_backbone')!=experimental_backbone_evidence.get(pid):raise ValueError('Candidate FeFe experimental backbone evidence changed')
        if pid in experimental_backbone_evidence:alternatives.setdefault(pid,{})['FeFe_experimental_backbone']=experimental_backbone_evidence[pid]
        if d.get('Assay_sequence_identity')!=assay_evidence.get(pid):raise ValueError('Candidate assayed sequence evidence changed')
        if pid in assay_evidence:alternatives.setdefault(pid,{})['Assay_sequence_identity']=assay_evidence[pid]
        if d.get('Coding_completeness')!=coding_rows.get(pid):raise ValueError('Candidate coding completeness evidence changed')
        if pid in coding_rows:alternatives.setdefault(pid,{})['Coding_completeness']=coding_rows[pid]
        ap=Path(ev.get('noncanonical',''))
        if ev.get('noncanonical'):
            if not ap.is_file(): raise ValueError('Declared candidate alternative file missing')
            alternatives.setdefault(d['candidate_id'],{}).update(next((x for x in table(ap) if x.get('protein_id')==d['candidate_id']),{})) ; bound.append(ap)
        subset=Path(ev.get('annotation',''))
        if not subset.is_dir(): raise ValueError('Confirmed candidate annotation subset missing')
        subset_manifest=subset/'candidate_subset_manifest.json'
        if not subset_manifest.is_file(): raise ValueError('Confirmed candidate subset manifest missing')
        subset_lock=load(subset_manifest)
        if Path(subset_lock.get('source_annotation','')).resolve() != annotation: raise ValueError('Candidate subset source annotation mismatch')
        if pid not in subset_lock.get('candidate_proteins',[]): raise ValueError('Candidate subset ownership mismatch')
        if subset_lock.get('genes_sha256') != sha(subset/'genes.tsv') or subset_lock.get('proteins_sha256') != sha(subset/'proteins.faa'):
            raise ValueError('Confirmed candidate subset artifacts changed')
        subset_gene_rows=table(subset/'genes.tsv'); subset_genes={r['id']:r for r in subset_gene_rows}
        if len(subset_genes)!=len(subset_gene_rows): raise ValueError('Duplicate candidate subset gene IDs')
        with (subset/'proteins.faa').open() as handle: subset_records=list(SeqIO.parse(handle,'fasta'))
        if len({r.id for r in subset_records})!=len(subset_records): raise ValueError('Duplicate candidate subset protein IDs')
        subset_proteins={r.id:str(r.seq) for r in subset_records}
        declared=set(subset_lock.get('candidate_proteins',[]))
        if declared != set(subset_genes) or declared != set(subset_proteins): raise ValueError('Candidate subset protein/gene ownership list mismatch')
        pid=d['candidate_id']
        if pid not in subset_genes or pid not in source_genes or subset_genes[pid] != source_genes[pid] or subset_proteins.get(pid) != source_proteins.get(pid): raise ValueError('Candidate subset differs from source annotation')
        backbone,backbone_inputs=reference_backbones(subset,{pid:feature},ROOT)
        if backbone.get(pid)!=d.get('NiFe_reference_backbone'):raise ValueError('Candidate backbone evidence changed')
        if pid in backbone:alternatives.setdefault(pid,{})['NiFe_reference_backbone']=backbone[pid]
        bound.extend(Path(p) for p in backbone_inputs)
        if feature['family']=='Fe':
            from canonical_hmd import identity_evidence
            fe,fe_inputs=identity_evidence(subset,{pid:feature},subset_gene_rows,ROOT,coding_completeness=(coding_rows,coding_inputs))
            if fe.get(pid)!=d.get('Fe_reference_identity'):raise ValueError('Candidate Fe reference evidence changed')
            alternatives.setdefault(pid,{})['Fe_reference_identity']=fe[pid]
            bound.extend(Path(p) for p in fe_inputs)
        bound.extend([subset_manifest,subset/'genes.tsv',subset/'proteins.faa'])
    return accepted, features, alternatives, bound


def explain(annotation,out,candidate_ids=None,rules=None,domain_evidence=None,candidate_screen=None,candidate_confirmation=None):
    annotation=Path(annotation).resolve();out=Path(out).resolve()
    if not (annotation/'catalytic_anchors.tsv').exists():raise ValueError('Missing candidate table')
    inputs=[annotation/n for n in ['catalytic_anchors.tsv','family_candidates.tsv','confirmation/features.json','active_center_summary.tsv',
        'noncanonical_centers/alternatives.tsv','complex_members.tsv','module_graph_nodes.tsv','domain_hits.tsv','complex_hypotheses.tsv','reaction_capabilities.tsv','genes.tsv','proteins.faa','catalytic_sites.tsv','run_manifest.json'] if (annotation/n).exists()]
    before={str(p):sha(p) for p in inputs}
    before[str(ROOT/'research_workflow/preferred_prediction.py')]=sha(ROOT/'research_workflow/preferred_prediction.py')
    before[str(ROOT/'research_workflow/mechanistic_functions.py')]=sha(ROOT/'research_workflow/mechanistic_functions.py')
    before[str(ROOT/'config/mechanistic_rules.json')]=sha(ROOT/'config/mechanistic_rules.json')
    with (annotation/'proteins.faa').open() as handle:
        protein_records=list(SeqIO.parse(handle,'fasta'))
    if len({r.id for r in protein_records})!=len(protein_records):raise ValueError('Duplicate source protein IDs')
    protein_sha256={r.id:hashlib.sha256(str(r.seq).upper().encode()).hexdigest() for r in protein_records}
    query_sequences={r.id:str(r.seq) for r in protein_records}
    with (ROOT/'references/processed/references.faa').open() as handle:
        reference_sequences={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    source_review=load_lineage_review(ROOT)
    for p in [ROOT/'research_workflow/identity_scope.py',ROOT/'research_workflow/identity_scope_source_review.json']:before[str(p)]=sha(p)
    anchors=table(annotation/'catalytic_anchors.tsv');features={r['protein_id']:r for r in load(annotation/'confirmation/features.json')}
    backbones,backbone_inputs=reference_backbones(annotation,features,ROOT);before.update(backbone_inputs)
    core_regions,region_inputs=nife_core_regions(annotation,features,ROOT);before.update(region_inputs)
    if rules is None: rules=load(ROOT/'config/catalytic_confirmation_thresholds.json')['families']
    extra_features={};extra_alternatives={};extra_context=[];candidate_bound=[]
    if (candidate_screen is None) != (candidate_confirmation is None): raise ValueError('Candidate screen and confirmation must be supplied together')
    if candidate_screen is not None:
        accepted,extra_features,extra_alternatives,candidate_bound=_candidate_context(annotation,candidate_screen,candidate_confirmation)
        source_genes={r['id']:dict(r,protein_sha256=protein_sha256[r['id']]) for r in table(annotation/'genes.tsv')}
        for d in accepted:
            pid=d['candidate_id']; fam=d.get('family',''); feature=extra_features[pid]
            alt=dict(extra_alternatives.get(pid,{})); alt['query_partial']=source_genes[pid].get('partial','False') in [True,'True']
            checked=sequence_status(dict(protein_id=pid,catalytic_identity_confidence='not_high_confidence'),feature,alt,rules.get(fam),source_review,
                                    coding_completeness=alt.get('Coding_completeness'),gene=source_genes[pid])
            if checked != d.get('sequence_status'): raise ValueError('Candidate confirmation sequence status is stale or bypassed')
            extra_alternatives[pid]=alt
            input_mode=load(annotation/'run_manifest.json')['input_mode']
            if input_mode not in {'assembly','genome','bin','transcript'}: raise ValueError('Unknown candidate input linkage mode')
            extra_context.append(dict(protein_id=pid,catalytic_family=fam,catalytic_identity_confidence=checked,structure_id='unresolved',structure_confidence='unresolved',structure_status='candidate_scoped_unresolved',function_status='unassigned',input_mode=input_mode,active_center_status=feature.get('source_active_center',{}).get('active_center_status','not_tested'),active_center_span_status=feature.get('source_active_center',{}).get('template_span_alignment_status','not_tested'),flags='partial_CDS' if alt['query_partial'] else '',missing_roles='',module_graph_search_status='not_tested',candidate_scoped_research_prediction='True'))
        features.update(extra_features)
        anchors.extend(extra_context)
        before.update({str(p):sha(p) for p in candidate_bound})
    alternatives={r['protein_id']:r for r in table(annotation/'noncanonical_centers/alternatives.tsv')}
    alternatives.update(extra_alternatives)
    source_genes={r['id']:dict(r,protein_sha256=protein_sha256[r['id']]) for r in table(annotation/'genes.tsv')}
    for a in anchors:
        qid=a['protein_id'];alternatives.setdefault(qid,{})['query_partial']=source_genes.get(qid,{}).get('partial')
        if qid in backbones:alternatives[qid]['NiFe_reference_backbone']=backbones[qid]
        if qid in core_regions:alternatives[qid]['NiFe_core_region']=core_regions[qid]
    discovered={r['protein_id'] for r in table(annotation/'family_candidates.tsv')}
    members=table(annotation/'complex_members.tsv');reactions=table(annotation/'reaction_capabilities.tsv');nodes=table(annotation/'module_graph_nodes.tsv')
    hypotheses=table(annotation/'complex_hypotheses.tsv')
    domains=table(annotation/'domain_hits.tsv')
    if domain_evidence is not None:
        from query_domain_organization import verify_completed
        domain_evidence=Path(domain_evidence).resolve();domain_manifest=verify_completed(domain_evidence)
        if str(annotation) not in domain_manifest['source_runs']:raise ValueError('Domain evidence belongs to another annotation source')
        domains=[d for d in table(domain_evidence/'query_domains.tsv') if d['source_run']==str(annotation)]
        for name in ['query_domain_manifest.json','query_domains.tsv']:before[str(domain_evidence/name)]=sha(domain_evidence/name)
    thresholds=load(annotation/'run_manifest.json')['thresholds']
    if candidate_ids is not None:
        if set(candidate_ids)-{a['protein_id'] for a in anchors}:raise ValueError('Requested candidate absent from saved prediction')
        anchors=[a for a in anchors if a['protein_id'] in candidate_ids]
    out.mkdir(parents=True,exist_ok=False);progression=[];attributes=[];frh_results={};module_role_inventory=None
    coding_rows={};coding_inputs={}
    core_ids={a['protein_id'] for a in anchors if a.get('candidate_scoped_research_prediction')!='True'}
    coding_ids=core_ids | {pid for pid,g in source_genes.items() if g.get('partial') in [True,'True']}
    coding_ids.update(coding_locus_candidate_ids(annotation, ROOT))
    if coding_ids:
        from coding_completeness import prepare as prepare_coding,load_evidence as load_coding
        coding_dir=out/'coding_completeness'
        _,coding_inputs=prepare_coding(annotation,coding_dir,candidate_ids=coding_ids,root=ROOT)
        _,coding_rows=load_coding(coding_dir,annotation=annotation,root=ROOT)
        coding_inputs.update({str(p):sha(p) for p in [coding_dir/'coding_completeness.json',coding_dir/'coding_completeness_manifest.json']})
        before.update(coding_inputs)
        for q,e in coding_rows.items():alternatives.setdefault(q,{})['Coding_completeness']=e
    function_coding={q:alt['Coding_completeness'] for q,alt in alternatives.items() if alt.get('Coding_completeness')}
    if domain_evidence is not None:
        from module_role_evidence import prepare as prepare_module_roles
        role_candidates=table(Path(candidate_screen)/'candidate_screen.tsv') if candidate_screen is not None else anchors
        if candidate_ids is not None: role_candidates=[c for c in role_candidates if c['protein_id'] in candidate_ids]
        module_role_inventory,role_inputs=prepare_module_roles(annotation,domain_evidence,role_candidates,out/'independent_roles',ROOT,coding_completeness=function_coding)
        before.update(role_inputs)
    from fefe_profile_ligands import prepare as prepare_profile_ligands
    core_features={a['protein_id']:features[a['protein_id']] for a in anchors if a.get('candidate_scoped_research_prediction')!='True' and a['protein_id'] in features}
    profile_evidence,profile_inputs=prepare_profile_ligands(annotation,core_features,out/'fefe_profile_ligands',rules.get('FeFe',{}),ROOT)
    before.update(profile_inputs)
    for q,e in profile_evidence.items():alternatives.setdefault(q,{})['FeFe_profile_ligands']=e
    from fefe_experimental_backbone import prepare as prepare_experimental_backbone
    backbone_features=core_features
    experimental_backbone_evidence,experimental_backbone_inputs=prepare_experimental_backbone(annotation,backbone_features,out/'fefe_experimental_backbone',rules.get('FeFe',{}),ROOT)
    before.update(experimental_backbone_inputs)
    for q,e in experimental_backbone_evidence.items():alternatives.setdefault(q,{})['FeFe_experimental_backbone']=e
    assay_evidence={};assay_inputs={}
    if any(f.get('family')=='FeFe' for f in core_features.values()):
        from assay_sequence_identity import prepare as prepare_assay
        assay_evidence,assay_inputs=prepare_assay(annotation,core_features,out/'assay_sequence_identity',ROOT)
        before.update(assay_inputs)
        for q,e in assay_evidence.items():alternatives.setdefault(q,{})['Assay_sequence_identity']=e
    from input_translation_evidence import audit as audit_translation
    translation_audit,translation_inputs=audit_translation(annotation,table(annotation/'genes.tsv'),table(annotation/'catalytic_sites.tsv'))
    before.update(translation_inputs);before[str(ROOT/'research_workflow/input_translation_evidence.py')]=sha(ROOT/'research_workflow/input_translation_evidence.py')
    (out/'input_translation_evidence.json').write_text(json.dumps(translation_audit,ensure_ascii=False,indent=2)+'\n')
    frh_model=None;frh_competition={}
    frh_results={}
    hmd_results={};canonical_hmd_results={}
    if any(a['catalytic_family']=='Fe' for a in anchors):
        from canonical_hmd import identity_evidence
        core_features={a['protein_id']:features[a['protein_id']] for a in anchors if a.get('candidate_scoped_research_prediction')!='True'}
        fe,fe_inputs=identity_evidence(annotation,core_features,table(annotation/'genes.tsv'),ROOT,coding_completeness=(coding_rows,coding_inputs));before.update(fe_inputs)
        for q,e in fe.items():alternatives.setdefault(q,{})['Fe_reference_identity']=e
        for a in anchors:
            q=a['protein_id'];e=alternatives.get(q,{}).get('Fe_reference_identity',{})
            if e.get('hmdii'):hmd_results[q]=e['hmdii']
            if e.get('canonical'):canonical_hmd_results[q]=e['canonical']
    identities={a['protein_id']:sequence_status(a,features.get(a['protein_id']),alternatives.get(a['protein_id'],{}),rules.get(a.get('catalytic_family','')),source_review,
        coding_completeness=alternatives.get(a['protein_id'],{}).get('Coding_completeness'),gene=source_genes.get(a['protein_id'])) for a in anchors}
    split_centers, split_inputs = prepare_split_identity(annotation, anchors, source_genes, query_sequences, ROOT)
    before.update(split_inputs)
    split_centers = {q: proof for q, proof in split_centers.items()
                     if identities[q] in {'evidence_insufficient', 'sequence_support_center_model_unresolved'}}
    identities.update(dict.fromkeys(split_centers, SPLIT_IDENTITY))
    from complementary_core_context import prepare as prepare_complementary_context, attributes as complementary_context_attributes
    complementary_anchors = [a for a in anchors if a.get('catalytic_family') == 'FeFe'
                             and identities[a['protein_id']] in {'calibrated_high', 'sequence_supported_uncalibrated_core'}]
    complementary_context, complementary_inputs = prepare_complementary_context(
        annotation, out/'complementary_core_context', complementary_anchors, source_genes, query_sequences, ROOT)
    before.update(complementary_inputs)
    fusion_results={};mvh_results={};pyridine_results={};pyridine_model=None;prosite_rows=[]
    context_results={};huc_results={};mbh_results={};ech_results={};hdcr_results={};ag_results={}
    fhl_results={};coo_results={};quinone_results={};tonn_results={};composition_results={}
    fdh_summary=dict(status='not_assessed_mechanistic_mode',candidate_count=0,supported_conditional_count=0)
    from reference_reactions import prepare as prepare_reference_reactions, VERSION as REFERENCE_REACTION_VERSION
    reference_results,reference_inputs=prepare_reference_reactions(annotation,anchors,identities,pyridine_results,ROOT,
        coding_completeness=function_coding,module_role_inventory=module_role_inventory)
    before.update(reference_inputs)
    if reference_results:(out/'reference_reaction_evidence.json').write_text(json.dumps(reference_results,ensure_ascii=False,indent=2)+'\n')
    for anchor in anchors:
        qid=anchor['protein_id'];family=anchor.get('catalytic_family','');feature=features.get(qid);parameters=rules.get(family);alt=alternatives.get(qid,{})
        coding=alt.get('Coding_completeness')
        coding_kwargs=dict(coding_completeness=coding,gene=source_genes.get(qid))
        fefe_profile_core=complete_fefe_core_support(feature,alt,parameters,source_review,**coding_kwargs)
        fefe_profile_mapped=fefe_profile_ligand_support(feature,alt,parameters,source_review,**coding_kwargs)
        fefe_experimental_backbone=fefe_experimental_backbone_support(feature,alt,parameters,source_review,**coding_kwargs)
        from assay_sequence_identity import record_supported as assay_supported
        assayed_core=assay_supported(alt.get('Assay_sequence_identity'),feature) if feature else False
        reasons=identity_reasons(feature,parameters) if anchor.get('catalytic_identity_confidence')!='calibrated_high' else []
        def add(stage,claim,code,kind,detail,source):reasons.append(dict(stage=stage,claim=claim,code=code,type=kind,detail=detail,source=source))
        selected=[m for m in members if m['anchor']==qid];flags=anchor.get('flags','').split(';');missing=[s for s in anchor.get('missing_roles','').split(';') if s]
        graph=[n for n in nodes if n['anchor']==qid and n['structure_id']==anchor['structure_id']]
        if 'no_reviewed_structure_rule' in missing:add('structure','candidate_architecture','reference_rule_missing','reference_rule_missing','No reviewed architecture rule; this is not evidence of absent query partners','catalytic_anchors.tsv:missing_roles')
        elif missing:
            structural_roles={n['role'] for n in graph if yes(n['required_for_structure'])}
            for stage,roles in [('structure',[r for r in missing if r in structural_roles]),('reaction',[r for r in missing if r not in structural_roles])]:
                if roles:add(stage,'reference_composition' if stage=='structure' else 'reference_reaction_participants','partner_not_detected','partner_not_detected',dict(roles=roles,input_mode=anchor.get('input_mode'),scope='Not detected under the saved reference rule and permitted linkage; not proven genomic absence'),'catalytic_anchors.tsv:missing_roles;module_graph_nodes.tsv')
        for stage,required,identity,coverage in [('structure','required_for_structure',thresholds['high_structure_identity'],thresholds['high_structure_coverage']),('reaction','required_for_function',thresholds['function_transfer_identity'],.8)]:
            below=[dict(role=n['role'],identity=float(n['aligned_pair_identity']),coverage=float(n['reference_coverage']),identity_required=identity,coverage_required=coverage)
                for n in graph if n[required] and (required!='required_for_structure' or yes(n[required])) and n['protein_ids'] and n['aligned_pair_identity'] and n['reference_coverage']
                and (float(n['aligned_pair_identity'])<identity or float(n['reference_coverage'])<coverage)]
            if below:add(stage,'legacy_high_structure' if stage=='structure' else 'precise_reference_reaction','module_similarity_below_rule','outside_calibration_scope',below,'module_graph_nodes.tsv;run_manifest.json:thresholds')
        if 'partial_CDS' in flags:add('structure','complete_reference_composition','partial_member_CDS','sequence_incomplete','At least one required member is explicitly partial','catalytic_anchors.tsv:flags;genes.tsv')
        if 'local_gene_cluster_incomplete_or_disrupted' in flags:add('structure','local_module_linkage','local_linkage_rule_not_met','partner_not_detected','Required local linkage/completeness constraint not satisfied; not proof of absent genes','catalytic_anchors.tsv:flags')
        if any(f.startswith(('ambiguous_modules','cross_anchor_provisional_member_overlap')) for f in flags):add('structure','partner_assignment','partner_ambiguous','partner_ambiguous',anchor.get('cross_anchor_pairing_status'),'complex_members.tsv;cross_anchor_associations.tsv')
        if 'catalytic_reference_ambiguous' in flags:add('structure','reference_architecture_choice','competing_catalytic_references','partner_ambiguous','Alternative catalytic reference assignments remain','catalytic_anchors.tsv:flags')
        if 'module_cysteine_topology_not_conserved' in flags:add('structure','legacy_high_structure_or_reaction','legacy_all_reference_Cys_gate','observed_conflict',dict(roles=sorted({m['role'] for m in selected if not yes(m.get('reference_cysteines_conserved'))}),scope='Legacy reference-Cys mismatch; not proof of metal-cluster loss or absent catalytic activity'),'complex_members.tsv:reference_cysteines_conserved')
        for flag in flags:
            if flag.startswith('module_below_function_transfer_threshold:'):add('reaction','precise_reference_reaction','legacy_module_transfer_threshold','outside_calibration_scope',dict(roles=flag.split(':',1)[1].split(','),scope='Legacy similarity/coverage rule; not an independently established biochemical boundary'),'catalytic_anchors.tsv:flags')
        if anchor.get('structure_confidence')!='high' and not any(r['stage']=='structure' for r in reasons):add('structure','high_structure','high_structure_rule_not_met','outside_calibration_scope',dict(structure_status=anchor.get('structure_status'),scope='Inspect per-module similarity, completeness, organization and center constraints'),'complex_members.tsv;complex_hypotheses.tsv')
        assigned=[r for r in reactions if r.get('complex_anchor',r.get('anchor'))==qid]
        if not assigned and reference_results.get(qid,{}).get('status') not in TUPLE_POTENTIAL_STATUSES and not any(r['stage']=='reaction' for r in reasons):add('reaction','precise_reaction','mechanism_specific_evidence_missing','not_tested','No validated discriminator or assigned reaction contract at this resolution','reaction_capabilities.tsv;reaction_contract_candidates.tsv')
        identity=identities[qid]
        if anchor.get('candidate_scoped_research_prediction') == 'True' and family != 'NiFe' and qid not in pyridine_results and qid not in hmd_results and qid not in canonical_hmd_results:
            add('reaction','precise_reaction','reference_rule_missing','reference_rule_missing',
                'Confirmed candidate family has no deployed function rule at this resolution; identity evidence is retained without a reaction claim.',
                'candidate_confirmation/candidate_confirmation.json')
        center='noncanonical_literature_pattern_supported' if alt.get('pattern_status')=='literature_ligand_pattern_candidate' else anchor.get('active_center_status','not_tested')
        translation_sites=[s for s in translation_audit['sensitive_sites'] if s['protein_id']==qid]
        if translation_sites:
            add('identity','genomic_ligand_state','translation_model_dependent_ligand','not_tested',translation_sites,'prodigal.gff;raw genomic codon;input_translation_evidence.json')
            center='translation_model_dependent_ligand_state'
        for r in reasons:assert r['type'] in TYPES
        reasons.sort(key=lambda r:['identity','structure','reaction'].index(r['stage']))
        progression.append(dict(candidate_id=qid,family=family,discovery_pass=qid in discovered,identity_status=identity,
            original_identity_status=anchor.get('catalytic_identity_confidence',''),center_model_status=center,query_coordination_status='unknown',module_graph_status=anchor.get('module_graph_search_status','not_tested'),
            identity_sequence_evidence=dict(local_profile_evidence=local_domain_evidence(feature) if feature else [],
                NiFe_indeterminate_competitor_core_support=family=='NiFe' and identity=='sequence_supported_uncalibrated_core' and intact_core_indeterminate_competitor_support(feature,parameters),
                NiFe_indeterminate_competitor_core_rule=NIFE_INDETERMINATE_CORE_RULE,
                FeFe_indeterminate_competitor_core_support=family=='FeFe' and identity=='sequence_supported_uncalibrated_core' and intact_core_indeterminate_competitor_support(feature,parameters),
                FeFe_indeterminate_competitor_core_rule=FEFE_INDETERMINATE_CORE_RULE,
                FeFe_indeterminate_profile_core_support=family=='FeFe' and identity=='sequence_supported_uncalibrated_core' and fefe_indeterminate_profile_core_support(feature,alt,parameters,source_review,function_coding,source_genes.get(qid)),
                FeFe_indeterminate_profile_core_rule=FEFE_INDETERMINATE_PROFILE_RULE,
                NiFe_reference_backbone_support=nife_backbone_support(feature,alt,parameters,source_review,**coding_kwargs),NiFe_reference_backbone_rule=NIFE_BACKBONE_RULE,
                NiFe_reference_backbone=alt.get('NiFe_reference_backbone'),
                HmdII_reference_support=hmd_results.get(qid),
                canonical_Hmd_reference_support=canonical_hmd_results.get(qid),
                complete_FeFe_Hdomain_support=fefe_profile_core,complete_FeFe_Hdomain_rule=FEFE_PROFILE_CORE_RULE,
                FeFe_experimental_backbone_support=fefe_experimental_backbone,FeFe_experimental_backbone_rule=FEFE_EXPERIMENTAL_BACKBONE_RULE,
                direct_ligand_status=feature.get('source_active_center',{}).get('direct_ligand_status','unknown') if feature else 'unknown',
                coding_completeness=coding,
                assayed_core_sequence=alt.get('Assay_sequence_identity'),
                nearest_positive_reference_lineage_unconfirmed=lineage_uncertain(feature,source_review),
                scope='Family sequence evidence only; full Pfam coverage, precise ligation and new high-confidence calibration remain separate.'),
            complex_status=anchor.get('structure_status','not_tested'),reaction_claim_status=anchor.get('function_status','unassigned'),
            first_blocking_reason=reasons[0]['code'] if reasons else '',all_blocking_reasons=reasons,missing_evidence_type=sorted({r['type'] for r in reasons}),
            calibration_scope=dict(original=parameters or 'family_not_calibrated',new_sequence_status='Uncalibrated; does not inherit canonical precision' if identity!='calibrated_high' else 'Original canonical identity scope only'),
            evaluation_truth_status='not_provided_to_this_export; never a prediction blocker',source_annotation=str(annotation),reference_reaction_status=reference_results.get(qid,{}).get('status','not_assessed')))
        if translation_sites:progression[-1]['identity_sequence_evidence']['translation_context']=dict(sites=translation_sites,called_protein_center_status=anchor.get('active_center_status'))
        if alt.get('FeFe_profile_ligands'):
            progression[-1]['identity_sequence_evidence'].update(FeFe_profile_ligand_support=fefe_profile_mapped,FeFe_profile_ligand_rule=FEFE_PROFILE_LIGAND_RULE,FeFe_profile_ligands=alt['FeFe_profile_ligands'])
        if alt.get('FeFe_experimental_backbone'):
            progression[-1]['identity_sequence_evidence'].update(FeFe_experimental_backbone_support=fefe_experimental_backbone,FeFe_experimental_backbone_rule=FEFE_EXPERIMENTAL_BACKBONE_RULE,FeFe_experimental_backbone=alt['FeFe_experimental_backbone'])
        if identity==NIFE_CORE_REGION_STATUS:
            progression[-1]['identity_sequence_evidence'].update(NiFe_core_region=alt['NiFe_core_region'],NiFe_core_region_rule=NIFE_CORE_REGION_RULE)
        def attr(name,status,value,dependencies,source,limit):attributes.append(dict(candidate_id=qid,attribute=name,status=status,value=value,evidence_dependencies=dependencies,evidence_source=source,limitations=limit))
        if qid in complementary_context:
            attributes.extend(dict(candidate_id=qid, **row) for row in complementary_context_attributes(complementary_context[qid]))
        if qid in split_centers:
            proof = split_centers[qid]
            progression[-1]['identity_sequence_evidence']['NiFe_split_core'] = proof
            attr('NiFe_split_catalytic_core_potential', 'sequence_supported_candidate', proof,
                 ['unique bound split ligand mapping', 'complete associated component CDSs'],
                 'catalytic_center_candidates.tsv;catalytic_center_components.tsv;active_center_summary.tsv',
                 'Split catalytic center sequence scope only. Assembly, occupancy, electron partner and reaction remain unresolved.')
        module_context=(module_role_inventory or {}).get('candidate_contexts',{}).get(qid,{})
        input_modules=[m for m in module_context.get('input_module_contexts',[]) if m['status']=='supported_role_module_homology']
        if input_modules:
            attr('Context_NuoEF_like_input_module_potential','supported_role_module_homology',input_modules,
                ['complete independently supported NuoE/NuoF roles','disjoint fusion or unique local same-strand role pairing','input ownership and separate anchor association tiers'],
                'independent_roles/inventory.json',
                'Electron-input role module potential in the permitted genome/contig/bin context. This does not assign it to the hydrogenase, establish physical assembly, turnover, nucleotide specificity or H2 coupling.')
        if identity==NIFE_CORE_REGION_STATUS:
            attr('NiFe_catalytic_core_region_potential','sequence_supported_core_region',alt['NiFe_core_region'],
                ['current protein and selected template alignment','complete template window and four ligand contexts','local GA profile and positive reference competition'],
                'candidate_progression.tsv;annotation/active_center_blast.tsv',
                'Uncalibrated catalytic-region sequence potential only. CDS remains incomplete; complete protein, maturation, small-subunit association, coupling, query activity and in-vivo direction remain unresolved.')
        if coding:
            attr('Protein_coding_completeness',coding['status'],coding,
                ['bound nucleotide and effective translation model','canonical start and terminal stop','exact current protein correspondence'],
                'candidate_confirmation/coding_completeness/' if anchor.get('candidate_scoped_research_prediction')=='True' else 'coding_completeness/',
                'Coding completeness under the bound translation model only. Source partial annotations, calibrated eligibility, catalytic ligands and function-specific requirements remain separate.')
        if assayed_core and not (fefe_profile_core or fefe_profile_mapped or fefe_experimental_backbone) and identity in ['calibrated_high','sequence_supported_uncalibrated_core']:
            attr('H_cluster_redox_core_potential','sequence_supported_candidate',alt['Assay_sequence_identity'],
                ['exact complete assayed reference core sequence','owned complete GA H-domain','source-bound direct activity evidence'],
                'assay_sequence_identity/assay_sequence_identity.json',
                'Uncalibrated core potential from reference sequence equivalence; tags and mature boundary are unbound. No query activity, atomic coordination, genomic CDS completeness, partner, direction or bifurcation assignment.')
        if (fefe_profile_core or fefe_profile_mapped or fefe_experimental_backbone) and identity in ['calibrated_high','sequence_supported_uncalibrated_core']:
            profile=[d for d in feature['domain_evidence'] if d['accession']=='PF02906' and d.get('supported') is True and d.get('evidence_kind')=='single_domain']
            attr('H_cluster_redox_core_potential','sequence_supported_candidate',
                dict(rule_version=FEFE_PROFILE_CORE_RULE if fefe_profile_core else FEFE_PROFILE_LIGAND_RULE if fefe_profile_mapped else FEFE_EXPERIMENTAL_BACKBONE_RULE,profile_evidence=profile,positive_reference_evidence=feature['positive_evidence'],competitor_margin=feature['competitor_margin'],
                    direct_ligand_pattern=feature['source_active_center']['direct_ligand_status'] if fefe_profile_core else alt['FeFe_profile_ligands'] if fefe_profile_mapped else alt['FeFe_experimental_backbone'],center_qualification=feature['source_active_center']['active_center_status']),
                ['experimental_backbone_complete_coverage' if fefe_experimental_backbone and not (fefe_profile_core or fefe_profile_mapped) else 'complete_GA_H_domain','trusted_positive_reference_competition','single_nonpartial_chain','direct_ligand_pattern'],
                'confirmation/features.json;active_center_template_evidence.tsv;noncanonical_centers/alternatives.tsv'+('' if fefe_profile_core else ';fefe_profile_ligands/PF02906.query_mapping.json' if fefe_profile_mapped else ';fefe_experimental_backbone/backbone_evidence.json'),
                'Uncalibrated H-cluster redox-core potential; no exact coordination, occupancy, native maturation, electron partner, direction, bifurcation or measured activity.')
            folds=[d for d in domains if d['protein_id']==qid and yes(d['passes']) and d['accession'].split('.')[0] in ['PF01058','PF02906']]
            fusion=any(a['accession'].split('.')[0]=='PF01058' and b['accession'].split('.')[0]=='PF02906' and int(a['query_end'])<int(b['query_start']) for a in folds for b in folds)
            if fusion:attr('FeFe_fused_SSU_like_relay_organization','profile_in_candidate_architecture',folds,
                ['same_polypeptide','ordered_GA_PF01058_then_PF02906'],str(domain_evidence/'query_domains.tsv') if domain_evidence else 'domain_hits.tsv',
                'Shared small-subunit-like fold fused to a FeFe H-domain; not NiFe catalytic identity, an external partner assignment or a complete electron-transfer mechanism.')
        et=[m for m in selected if m.get('module_category','').startswith('electron_transfer')]
        attr('electron_transfer_module_presence','reference_like_modules_detected' if et else 'not_assessed_or_not_detected',
            [dict(role=m['role'],protein_id=m['protein_id'],reference_protein=m['reference_protein'],query_start=m['query_domain_start'],query_end=m['query_domain_end']) for m in et],
            ['module_alignment','permitted_member_linkage'],'complex_members.tsv','Module homology only; exact FeS ligation, cofactor specificity, interaction and mechanism remain separate.')
        attr('candidate_complex_architecture',anchor.get('structure_status','unresolved'),anchor.get('structure_id','unresolved'),
            ['required_modules','resource_ownership','candidate_pairing'],'complex_members.tsv;complex_hypotheses.tsv','Candidate composition is not measured physical assembly.')
        for name,accessions in [('nucleotide_binding_domain_in_associated_module',{'PF00175','PF01512','PF13450'}),('flavin_binding_domain_in_associated_module',{'PF00970','PF00941','PF12831'})]:
            matched=[]
            for domain in domains:
                if domain['accession'].split('.')[0] not in accessions or not yes(domain['passes']):continue
                positions=set(range(int(domain['query_start']),int(domain['query_end'])+1));resources=set()
                for member in selected:
                    if member['protein_id']==domain['protein_id']:resources.update(range(int(member['query_domain_start']),int(member['query_domain_end'])+1))
                if positions and positions<=resources:matched.append({k:domain[k] for k in ['protein_id','accession','profile','query_start','query_end','model_coverage','independent_evalue']})
            attr(name,'domain_profile_supported' if matched else 'not_assessed_or_not_detected',matched,
                ['passing_domain_profile','domain_coordinates_within_assigned_module_resources'],str(domain_evidence/'query_domains.tsv')+';complex_members.tsv' if domain_evidence else 'domain_hits.tsv;complex_members.tsv',
                'A fold/profile assignment supports a partial module annotation, not substrate turnover, NAD/NADP selectivity or obligatory electron coupling.')
        partner_hypotheses={}
        for hypothesis in hypotheses:
            if hypothesis['anchor']!=qid:continue
            for protein in json.loads(hypothesis['required_member_ids']):partner_hypotheses.setdefault(protein,[]).append(hypothesis['structure_id'])
        candidate_domains=[dict(**{k:d[k] for k in ['protein_id','accession','profile','query_start','query_end','model_coverage','independent_evalue']},
            candidate_architectures=sorted(set(partner_hypotheses[d['protein_id']]))) for d in domains
            if d['accession'].split('.')[0] in ['PF00175','PF01512','PF13450','PF00970','PF00941','PF12831'] and yes(d['passes']) and d['protein_id'] in partner_hypotheses]
        attr('candidate_partner_domain_profiles','profile_in_candidate_architecture' if candidate_domains else 'not_assessed_or_not_detected',candidate_domains,
            ['passing_domain_profile','member_of_saved_permitted_architecture_hypothesis'],'domain_hits.tsv;complex_hypotheses.tsv',
            'A profile on a proposed partner is retained even when architecture selection or the center is unresolved. No domain-to-role assignment, physical interaction, specificity or precise reaction is established.')
        attr('nucleotide_selectivity','unresolved',None,['mechanism_specific_cofactor_discrimination'],'reaction_contract_candidates.tsv','A nucleotide/flavin module or positive reference reaction does not establish NAD versus NADP exclusivity.')
        for name in ['H2_reaction_role','ferredoxin_requirement','electron_coupling']:
            attr(name,'unresolved',None,['query-specific reaction or mechanism evidence'],'reference_reaction_evidence.json;reaction_capabilities.tsv',
                'Reference-assay capabilities are reported separately; they do not measure the query reaction, obligatory partners or mechanism in its biological context.')
        attr('query_activity','not_tested',None,['query_measurement'],'no query experiment supplied','No sequence output establishes occupancy, maturity, activity or in-vivo direction.')
        reference=reference_results.get(qid,{})
        if reference.get('capabilities'):
            attr('Reference_assay_reaction_capabilities',reference['status'],reference['capabilities'],
                ['bound typed reaction registry','existing assigned reference, complete exact full-CDS protein-set sequence equivalence, or bound complete CDS coding-locus equivalence','independent catalytic identity and claim-scoped ownership'],
                'reference_reaction_evidence.json;reaction_capabilities.tsv',
                'Reference-assay molecular potential only. Inherited homology, exact full-CDS protein-set sequence equivalence and complete-CDS coding-locus evidence are separate; coding-locus potential does not resolve expression, called start, mature construct or query activity, selectivity exclusivity, cofactor occupancy, physical assembly or in-vivo net direction.')

    write(out/'candidate_progression.tsv',progression,FIELDS);write(out/'function_attributes.tsv',attributes,ATTR_FIELDS)
    from mechanistic_functions import evaluate_all as evaluate_mechanistic, RULE_VERSION as MECHANISTIC_RULE_VERSION
    site_gates_by_id={}
    for r in progression:
        qid=r['candidate_id']; gates={}
        if qid in fhl_results:
            gates['fdh_site_status']='supported' if fhl_results[qid].get('coupling_status') in ['candidate_H2_formate_coupling','candidate_genome_H2_formate_pathway'] else 'unverified'
        if qid in canonical_hmd_results:
            gates['hmd_gate_status']='canonical_hmd' if canonical_hmd_results[qid].get('reaction_status')=='compatible_reaction_candidate_uncalibrated' or canonical_hmd_results[qid].get('branch_reaction_status')=='compatible_reaction_candidate_uncalibrated' else None
        elif qid in hmd_results and hmd_results[qid].get('hmdii_status')=='reference_supported_uncalibrated':
            gates['hmd_gate_status']='hmdii_paralog'
        site_gates_by_id[qid]=gates
    mechanistic_results=evaluate_mechanistic(anchors,progression,source_genes,domains,root=ROOT,site_gates_by_id=site_gates_by_id)
    (out/'mechanistic_function_evidence.json').write_text(json.dumps(mechanistic_results,ensure_ascii=False,indent=2)+'\n')
    functional=[]
    for row in progression:
        qid=row['candidate_id']; attrs=[a for a in attributes if a['candidate_id']==qid]
        supported=[a['attribute'] for a in attrs if a['status'] in SUPPORTED_ATTRIBUTE_STATES]
        mech=mechanistic_results[qid]
        label=mech['function_prediction']
        if row['identity_status']==SPLIT_IDENTITY and label=='hydrogenase_candidate_identity_unresolved':label=SPLIT_LABEL
        coarse_label=label;reference=reference_results.get(qid,{})
        for rname in mech['supported_roles']:
            if rname not in supported: supported.append(rname)
        functional.append(dict(candidate_id=qid,family=row['family'],identity_status=row['identity_status'],function_prediction=label,
            coarse_function_prediction=coarse_label,mechanistic_subtype=mech.get('mechanistic_subtype',''),reference_reaction_status=reference.get('status','not_assessed'),
            reference_contract_ids=[c['contract_id'] for c in reference.get('capabilities',[])],
            reference_reaction_keys=[c['reaction_key'] for c in reference.get('capabilities',[])],
            reference_assay_reactions=[c['molecular_reaction'] for c in reference.get('capabilities',[])],
            reference_assay_capabilities=reference.get('capabilities',[]),
            supported_attributes=supported,unresolved_attributes=[a['attribute'] for a in attrs if a['status']=='unresolved'],
            confidence='uncalibrated_core_region' if row['identity_status']==NIFE_CORE_REGION_STATUS else 'uncalibrated_split_center' if label==SPLIT_LABEL else reference.get('status') if label in {'reference_sequence_matched_reaction_potential','reference_coding_locus_matched_reaction_potential'} else mech['confidence'],
            NAD_NADP_selectivity=mech['NAD_NADP_selectivity'],in_vivo_direction='unknown',query_activity='not_tested',FeFe_rubrerythrin_assay_capabilities='not_assessed',FeFe_rubrerythrin_caveat='not_assessed',
            evidence=(str(Path(candidate_confirmation).resolve()/'candidate_confirmation.json')+';' if any(a['protein_id']==qid and a.get('candidate_scoped_research_prediction')=='True' for a in anchors) else '')+'mechanistic_function_evidence.json;function_attributes.tsv'))
    claim_results=dict(mechanistic=mechanistic_results,canonical_hmd=canonical_hmd_results,hmdii=hmd_results,reference=reference_results)
    preferred={}
    screened={r['protein_id']:r for r in table(Path(candidate_screen)/'candidate_screen.tsv')} if candidate_screen else {}
    progress_by_id={r['candidate_id']:r for r in progression}
    for record in functional:
        qid=record['candidate_id']
        prediction=preferred_predict(record,evidence=dict(
            attributes=[a for a in attributes if a['candidate_id']==qid],
            screen=screened.get(qid,dict(domain_evidence=[d for d in domains if d.get('protein_id')==qid])),
            progression=progress_by_id[qid],
            identity_feature=features.get(qid,{}),
            claims={name:results[qid] for name,results in claim_results.items() if qid in results},
            gene=source_genes.get(qid,{}),genes=source_genes,
            module_graph_nodes=preferred_role_options(nodes,qid,next((a['structure_id'] for a in anchors if a['protein_id']==qid),'')),
            role_context=(module_role_inventory or {}).get('candidate_contexts',{}).get(qid,{})))
        record.update(prediction);preferred[qid]=prediction
    (out/'preferred_predictions.json').write_text(json.dumps(preferred,ensure_ascii=False,indent=2)+'\n')
    from function_rules import FIELDS as FUNCTION_RULE_FIELDS, annotate as annotate_function_rules
    class_predictions,class_manifest=annotate_function_rules(annotation,out/'function_rules',functional)
    for record in functional:record.update(class_predictions[record['candidate_id']])
    write(out/'functional_annotations.tsv',functional,['candidate_id',*PREDICTION_FIELDS,*FUNCTION_RULE_FIELDS,'family','identity_status','function_prediction','supported_attributes','unresolved_attributes','confidence','NAD_NADP_selectivity','in_vivo_direction','query_activity','FeFe_rubrerythrin_assay_capabilities','FeFe_rubrerythrin_caveat','evidence','coarse_function_prediction','mechanistic_subtype','reference_reaction_status','reference_contract_ids','reference_reaction_keys','reference_assay_reactions','reference_assay_capabilities'])
    assert all(sha(Path(p))==digest for p,digest in before.items())
    manifest=dict(status='claim_specific_explanation_complete',candidate_count=len(progression),attribute_rows=len(attributes),input_sha256=before,
        script_sha256=sha(Path(__file__)),identity_rules=rules,changes_original_predictions=False,
        domain_evidence=str(domain_evidence) if domain_evidence else 'legacy_domain_hits_only_not_a_complete_domain_inventory',
        reason_semantics='Reasons apply to named claims, not a serial universal gate. First reason is display order, not a causal assertion that all later claims depend on it.',
        truth_semantics='Evaluation truth availability never blocks a prediction. No biological accuracy or new calibrated_high scope is asserted.',
        frh_function_scope='Intact local three-protein candidate; partial module and F420-compatible potential separately uncalibrated; no exact-FeS gate.' if frh_results else 'not_applicable',
        frh_output_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and ('frh' in str(p.relative_to(out)))},
        mvh_hdr_scope='Claim-specific local relay, nonoverlapping fusion and uncalibrated coupling hypothesis; disulfide specificity/Fd/bifurcation remain separate.' if mvh_results else 'not_applicable',
        input_translation_evidence=dict(path='input_translation_evidence.json',sha256=sha(out/'input_translation_evidence.json'),status=translation_audit['status']),
        coding_completeness_output_sha256={str(p.relative_to(out)):sha(p) for p in (out/'coding_completeness').glob('*') if p.is_file()},
        coding_completeness_scope='Bound coding integrity may support an uncalibrated claim despite a source partial flag; source flags and calibrated eligibility remain unchanged.',
        FeFe_profile_output_sha256={str(p.relative_to(out)):sha(p) for p in (out/'fefe_profile_ligands').glob('*') if p.is_file()},
        FeFe_experimental_backbone_output_sha256={str(p.relative_to(out)):sha(p) for p in (out/'fefe_experimental_backbone').glob('*') if p.is_file()},
        candidate_scoped_research_predictions=sorted(a['protein_id'] for a in anchors if a.get('candidate_scoped_research_prediction')=='True'),
        candidate_context_metadata={'screen':str(Path(candidate_screen).resolve()) if candidate_screen else None,'confirmation':str(Path(candidate_confirmation).resolve()) if candidate_confirmation else None,'confirmation_manifest_sha256':sha(Path(candidate_confirmation)/'candidate_confirmation_manifest.json') if candidate_confirmation else None,'scope':'Extra candidates are research-scoped and never written to core anchor/module tables.'},
        mechanistic_rule_version=MECHANISTIC_RULE_VERSION,
        mechanistic_output_sha256=sha(out/'mechanistic_function_evidence.json'),
        canonical_hmd_output_sha256=sha(out/'canonical_hmd_evidence.json') if canonical_hmd_results and (out/'canonical_hmd_evidence.json').is_file() else None,
        reference_reaction_rule_version=REFERENCE_REACTION_VERSION if reference_inputs else 'not_applicable',
        reference_reaction_output_sha256=sha(out/'reference_reaction_evidence.json') if reference_results else None,
        reference_reaction_scope='Source assay capabilities only; exact full-CDS protein-set potential, coding-locus potential and inherited homology reference assignments are distinct from raw query-protein identity (including shorter or longer calls and raw initiator representation), query activity, specificity exclusivity and in-vivo direction.',
        function_rules_version=class_manifest['version'],
        function_rules_manifest_sha256=sha(out/'function_rules/manifest.json'),
        function_rules_class_counts=class_manifest['predicted_class_counts'],
        functional_annotations_sha256=sha(out/'functional_annotations.tsv'),
        function_attributes_sha256=sha(out/'function_attributes.tsv'),
        complementary_core_context_output_sha256=sha(out/'complementary_core_context_evidence.json') if (out/'complementary_core_context_evidence.json').is_file() else None,
        preferred_prediction_output_sha256=sha(out/'preferred_predictions.json'),
        preferred_prediction_policy_version=PREFERRED_POLICY_VERSION,
        preferred_prediction_policy_sha256=sha(ROOT/'research_workflow/preferred_prediction.py'),
        preferred_prediction_semantics='Primary sequence-function best explanation; qualitative high/medium/low is not a calibrated probability. Experimental unknowns and independent validation do not block sequence prediction.',
        module_role_inventory_sha256=sha(out/'independent_roles/inventory.json') if module_role_inventory is not None else None,
        module_role_table_sha256=sha(out/'independent_roles/gene_roles.tsv') if module_role_inventory is not None else None,
        pyridine_role_composition_output_sha256=sha(out/'pyridine_role_composition_evidence.json') if composition_results else None,
        assayed_core_input_sha256=assay_inputs,
        pyridine_output_sha256=sha(out/'pyridine_function_evidence.json') if pyridine_results else None,
        mvh_output_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and ('mvh' in str(p.relative_to(out)))})
    (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    return progression,attributes


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--annotation',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--domain-evidence',type=Path,help='Verified full query-domain output containing this exact source annotation')
    parser.add_argument('--candidate-screen',type=Path)
    parser.add_argument('--candidate-confirmation',type=Path)
    args=parser.parse_args();rows,attributes=explain(args.annotation,args.out,domain_evidence=args.domain_evidence,candidate_screen=args.candidate_screen,candidate_confirmation=args.candidate_confirmation);print(f'{len(rows)} candidates; {len(attributes)} attribute rows')
