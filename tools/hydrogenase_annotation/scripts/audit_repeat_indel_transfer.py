"""Audit a conservative repeat guard using frozen traces; never predict final proteins."""
import ast
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from annotate import load_blast,residue_states,repeat_indel_ambiguity,repeat_mapping_ambiguity
from active_centers import collinear_hits
from catalytic_confirmation import stage_fingerprint,qualifies,source_family
from evaluate_identity_corpus import identity_metrics
from identity_truth import load_units
from audit_input_translation_transfer import normalized_annotate_source


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def require(condition,message):
    if not condition:raise ValueError(message)


def retained_bounds(truth,high_ids,possibly_changed):
    old=identity_metrics(truth,high_ids);retained=set(high_ids)-set(possibly_changed)
    lower=identity_metrics(truth,retained);positives=old['TP']+old['FN']
    return dict(TP_lower=lower['TP'],TP_upper=old['TP'],FP_upper=old['FP'],FN_lower=old['FN'],FN_upper=positives-lower['TP'],
        precision_lower=lower['TP']/(lower['TP']+old['FP']) if lower['TP'] else None,
        recall_lower=lower['TP']/positives if positives else None,truth_positives=positives,
        changed_high_upper=len((set(high_ids)&{r['unit_id'] for r in truth})&set(possibly_changed)),
        wilson95_retained=lower['precision_wilson95'])


def normalized_source(text,kind):
    if kind=='annotate':
        text=normalized_annotate_source(text)
    tree=ast.parse(text)
    if kind=='annotate':
        tree.body=[n for n in tree.body if not isinstance(n,ast.FunctionDef) or n.name not in ['residue_states','repeat_indel_ambiguity','repeat_mapping_ambiguity']]
    elif kind=='active_centers':
        tree.body=[n for n in tree.body if not isinstance(n,ast.FunctionDef) or n.name not in ['map_model','map_hits','annotate_active_centers']]
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom) and node.module=='annotate':node.names=[n for n in node.names if n.name!='repeat_mapping_ambiguity']
    elif kind=='split_centers':
        for node in ast.walk(tree):
            if isinstance(node,ast.FunctionDef) and node.name=='annotate_split_centers':
                node.body=[n for n in node.body if not isinstance(n,ast.Assign) or not any(isinstance(t,ast.Name) and t.id=='sequences' for t in n.targets)]
            if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='map_model':node.keywords=[k for k in node.keywords if k.arg!='query_sequence']
    elif kind=='noncanonical_center_evidence':
        for node in ast.walk(tree):
            if isinstance(node,ast.FunctionDef) and node.name=='candidate_rows' and node.args.args[-1].arg=='query_sequences':node.args.args.pop();node.args.defaults.pop()
            if isinstance(node,ast.Call):node.keywords=[k for k in node.keywords if k.arg!='query_sequences']
    return ast.dump(tree,include_attributes=False)


def affected_frozen_highs(root,annotation,high_ids):
    """Inspect only already deposited alignment coordinates, with the locked guard."""
    models={m['reference_id']:m for m in read(root/'references/processed/active_centers/models.json')}
    policy=read(root/'config/active_center_review.json')['parameters']
    raw=load_blast(annotation/'active_center_blast.tsv',dict(blast_evalue=policy['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
    by_pair=defaultdict(list)
    for hit in raw:
        if hit['qseqid'] in high_ids:by_pair[hit['qseqid'],hit['sseqid']].append(hit)
    with (annotation/'proteins.faa').open() as handle:sequences={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta') if r.id in high_ids}
    require(set(sequences)==set(high_ids),'Frozen high-ID sequence ownership changed')
    rows=[]
    for (qid,rid),hits in sorted(by_pair.items()):
        model=models[rid];mapped={};ambiguous=set()
        for hit in collinear_hits(hits):
            _,piece=residue_states(hit,range(model['window_start'],model['window_end']+1));mapped.update(piece)
            ambiguous.update(repeat_indel_ambiguity(hit))
        ambiguous.update(repeat_mapping_ambiguity(mapped,model['sequence'],sequences[qid]))
        required={s['reference_position'] for s in model['sites'] if s['requirement'] in ['required_for_template_coordination','canonical_reaction_signature']}
        if ambiguous&required:rows.append(dict(protein_id=qid,reference_id=rid,ambiguous_required_positions=sorted(ambiguous&required)))
    return rows


def verify_transition(root,original_verifier):
    root=Path(root);lock_path=root/'references/validation_locks/round47_repeat_indel_transition.json';lock=read(lock_path)
    require(lock['status']=='guard_frozen_before_final_trace_audit','Guard audit was not preregistered')
    require(sha(root/'experiments/round47/plan.json')==lock['plan_sha256'],'Guard transition plan changed')
    require(sha(Path(__file__))==lock['auditor_sha256'],'Guard trace auditor changed after registration')
    require(stage_fingerprint(root)==lock['feature_stage_fingerprint'],'Guard predictor changed after registration')
    baseline=Path(lock['baseline_view']);plan=read(root/'experiments/round47/plan.json')
    for name,value in plan['production_before'].items():require(sha(baseline/name)==value,'Baseline source changed: '+name)
    previous=original_verifier(baseline)
    for name in ['annotate','active_centers','split_centers','noncanonical_center_evidence']:
        require(normalized_source((root/'scripts'/(name+'.py')).read_text(),name)==normalized_source((baseline/'scripts'/(name+'.py')).read_text(),name),'Unreviewed change outside repeat guard/caller plumbing: '+name)
    old_rule=read(baseline/'config/catalytic_confirmation_thresholds.json');rule=read(root/'config/catalytic_confirmation_thresholds.json')
    require(rule['families']==old_rule['families']==lock['families'],'Guard changed calibrated decision parameters')
    cal_path=Path(rule['calibration_summary']);cal=read(cal_path);old_cal=read(Path(old_rule['calibration_summary']))
    require(sha(cal_path)==rule['calibration_summary_sha256']==lock['calibration_summary_sha256'],'New calibration integrity changed')
    require(cal['feature_stage_fingerprint']==lock['feature_stage_fingerprint'],'New calibration belongs to another predictor')
    for name,value in cal['artifacts'].items():require(sha(cal_path.parent/name)==value,'New calibration artifact changed')
    before={f['protein_id']:f for f in read(Path(old_rule['calibration_summary']).parent/'annotation/confirmation/features.json')}
    after={f['protein_id']:f for f in read(cal_path.parent/'annotation/confirmation/features.json')}
    require(set(before)==set(after),'Calibration query set changed')
    keys=['family','core_positions','core_domain_coverage','positive_score','competitor_score','competitor_margin','positive_prototype','competitor_prototype','positive_evidence','competitor_evidence','domain_evidence']
    for qid,new in after.items():
        old=before[qid]
        if new['family'] not in ['NiFe','FeFe']:continue
        require(all(new[k]==old[k] for k in keys),'Guard changed identity feature mathematics: '+qid)
        require(not new['eligible'] or old['eligible'],'Guard promoted calibration eligibility')
    final=root/'results/round17_final_confirmation';summary=read(final/'summary.json');features=read(final/'annotation/confirmation/features.json')
    for name,value in summary['artifacts'].items():require(sha(final/name)==value,'Frozen final artifact changed')
    chosen={f['protein_id'] for f in features if old_rule['families'].get(f['family']) and qualifies(f,old_rule['families'][f['family']])}
    affected=affected_frozen_highs(root,final/'annotation',chosen);possibly_changed={r['protein_id'] for r in affected}
    audit=read(Path(summary['source_audit']));require(sha(Path(summary['source_audit']))==summary['source_audit_sha256'],'Final source audit changed')
    corpus=root/read(root/'config/catalytic_confirmation.json')['corpus'];truth=[r for r in load_units(corpus) if r['split']=='final_validation' and r['unit_id'] in set(audit['independent_unit_ids'])]
    bounds={'all':retained_bounds(truth,chosen,possibly_changed)};family_checks={};selection=read(root/'config/catalytic_confirmation.json')['calibration_selection']
    byid={f['protein_id']:f for f in features}
    for family in ['NiFe','FeFe']:
        family_truth=[r for r in truth if r['label']=='negative' or source_family(r)==family]
        family_chosen={q for q in chosen if byid[q]['family']==family};bounds[family]=retained_bounds(family_truth,family_chosen,possibly_changed)
        b=bounds[family];family_checks[family]=bool(b['TP_lower']>=selection['min_true_positive_outputs_per_family'] and b['precision_lower']>=selection['min_precision'] and b['recall_lower']>=selection['min_recall'])
    checks=previous['checks']+['registered_guard_and_baseline_sources','unchanged_identity_math_outside_guard',
        'same_calibration_features_and_queries','no_calibration_eligibility_promotion','frozen_alignment_guard_impact_bounds']
    return dict(status='conservative_guard_transition_verified',validation_scope='repeat_guard_conservative_subset',parent_final=str(final/'summary.json'),
        parent_final_summary_sha256=sha(final/'summary.json'),parent_thresholds=previous['parent_thresholds'],parent_thresholds_sha256=previous['parent_thresholds_sha256'],
        current_feature_stage_fingerprint=stage_fingerprint(root),current_thresholds_sha256=sha(root/'config/catalytic_confirmation_thresholds.json'),
        baseline_transition_status=previous['status'],checks=checks,passed=len(checks),calibration_queries=len(before),
        inherited_final_metrics=summary['metrics']['source_audited_independent'],current_metric_bounds=bounds,family_checks=family_checks,
        precision_gate_supported=bounds['all']['truth_positives']>=100 and all(family_checks.values()) and bounds['all']['precision_lower']>=.99,
        possibly_changed_high_ids=sorted(possibly_changed),guard_trace_evidence=affected,new_final_prediction_runs=0,
        interpretation='Fixed prototype/core scores and thresholds; the reviewed guard only removes eligible high outputs. Frozen alignments identify a conservative superset of potentially removed high IDs. Counts are bounds, not newly measured final predictions. No final data fitted or searched.',
        limits=['Conditional finite identity evidence only; no new final test or population guarantee.','No transfer for Hmd, complex, function or FeS accuracy.','Any guard/dependency/threshold change invalidates this locked audit.'])
