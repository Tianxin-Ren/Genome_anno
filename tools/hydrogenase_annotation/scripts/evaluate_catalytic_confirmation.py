#!/usr/bin/env python3
"""Calibrate catalytic identity separately, then evaluate a locked source-separated final set."""
import argparse
import itertools
import json
import subprocess
import sys
from pathlib import Path
from annotate import ROOT,write_tsv
from catalytic_confirmation import sha,table,qualifies,stage_fingerprint,source_family
from evaluate_identity_corpus import identity_metrics
from run_acceptance_checks import state_hashes


def prediction_metrics(truth,features,rules):
    selected={r['protein_id'] for r in features if rules.get(r['family']) and qualifies(r,rules[r['family']])}
    metrics=identity_metrics(truth,selected)
    byid={r['protein_id']:r for r in features}
    wrong=[r['unit_id'] for r in truth if r['unit_id'] in selected and r['label']=='positive' and source_family(r)!=byid[r['unit_id']]['family']]
    metrics['wrong_family_high_outputs']=wrong
    return metrics,selected


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--split',choices=['calibration','final_validation'],required=True)
    parser.add_argument('--source-audit',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--lock',type=Path,help='Required preprediction lock for final validation')
    args=parser.parse_args();out=args.out.resolve()
    if (out/'summary.json').exists():parser.error('Use a fresh prediction directory')
    policy_path=ROOT/'config/catalytic_confirmation.json';policy=json.loads(policy_path.read_text());corpus=ROOT/policy['corpus'];freeze=json.loads((corpus/'freeze.json').read_text())
    for name,digest in freeze['artifacts'].items():
        if sha(corpus/name)!=digest:raise ValueError('Confirmation evaluation source changed')
    audit_path=args.source_audit.resolve()/'summary.json';audit=json.loads(audit_path.read_text())
    if (audit['status']!='source_audit_frozen_before_confirmation_prediction' or audit['corpus_freeze_sha256']!=sha(corpus/'freeze.json') or
        audit['confirmation_model_manifest_sha256']!=sha(ROOT/'references/processed/catalytic_confirmation/manifest.json')):raise ValueError('Confirmation source audit is stale')
    for name,digest in audit['artifacts'].items():
        if sha(audit_path.parent/name)!=digest:raise ValueError('Confirmation source audit artifact changed')
    from identity_truth import load_units
    all_rows=load_units(corpus);rows=[r for r in all_rows if r['split']==args.split]
    independent=[r for r in rows if r['unit_id'] in set(audit['independent_unit_ids'])]
    rule_path=ROOT/'config/catalytic_confirmation_thresholds.json';fingerprint=stage_fingerprint(ROOT)
    if args.split=='calibration':
        if rule_path.exists():parser.error('Calibration rule already exists; preserve it and its reserved final test')
        locked=None
    else:
        if not args.lock:parser.error('Final confirmation validation requires a preprediction --lock')
        locked=json.loads(args.lock.read_text());calibration=Path(locked['calibration_summary']);cal=json.loads(calibration.read_text())
        if (locked['stage']!='catalytic_identity_confirmation_only' or locked['feature_stage_fingerprint']!=fingerprint or
            locked['thresholds_sha256']!=sha(rule_path) or locked['calibration_summary_sha256']!=sha(calibration) or
            locked['source_audit_sha256']!=sha(audit_path) or locked['acceptance_targets_sha256']!=sha(ROOT/'config/acceptance_targets.json') or
            cal['status']!='confirmation_calibration_complete' or cal['feature_stage_fingerprint']!=fingerprint or cal['source_audit_sha256']!=sha(audit_path)):
            raise ValueError('Final confirmation lock/calibration is missing or changed')
    out.mkdir(parents=True,exist_ok=True);before=state_hashes()
    command=[sys.executable,str(ROOT/'scripts/annotate.py'),'--proteins',str(corpus/(args.split+'.faa')),'--mode','transcript','--threads','4','--out',str(out/'annotation')]
    with (out/'annotation.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    features=json.loads((out/'annotation/confirmation/features.json').read_text());byid={r['protein_id']:r for r in features};grid=[]
    if {r['unit_id'] for r in rows}!=set(byid):raise ValueError('Confirmation query ownership changed')
    if args.split=='calibration':
        rules={};selection=policy['calibration_selection']
        for family in policy['eligible_calibration_families']:
            truth=[r for r in independent if r['label']=='negative' or source_family(r)==family];options=[]
            for values in itertools.product(*policy['parameter_grid'].values()):
                parameters=dict(zip(policy['parameter_grid'],values));metrics,_=prediction_metrics(truth,features,{family:parameters})
                row=dict(family=family,parameters=parameters,**metrics);grid.append(row)
                if (metrics['TP']>=selection['min_true_positive_outputs_per_family'] and metrics['precision']>=selection['min_precision'] and
                    metrics['recall']>=selection['min_recall'] and not metrics['wrong_family_high_outputs']):options.append(row)
            best=max(options,key=lambda r:(r['TP'],r['parameters']['min_competitor_margin'],r['parameters']['min_domain_coverage'],r['parameters']['min_core_similarity']),default=None)
            rules[family]=best['parameters'] if best else None
        write_tsv(out/'parameter_grid.tsv',[dict(r,parameters=json.dumps(r['parameters']),wrong_family_high_outputs=json.dumps(r['wrong_family_high_outputs'])) for r in grid])
    else:
        rules=json.loads(rule_path.read_text())['families']
    metrics={};selected_by_stratum={}
    reps={}
    for r in sorted(independent,key=lambda r:r['unit_id']):reps.setdefault((r['homolog_component'],r['label'],source_family(r)),r)
    strata=dict(all_frozen_units=rows,source_audited_independent=independent,independent_component_representatives=list(reps.values()))
    for family in policy['eligible_calibration_families']:
        strata['family_'+family]=[r for r in independent if r['label']=='negative' or source_family(r)==family]
    for name,truth in strata.items():
        active_rules={name[7:]:rules[name[7:]]} if name.startswith('family_') else rules
        metrics[name],selected_by_stratum[name]=prediction_metrics(truth,features,active_rules)
    chosen=selected_by_stratum['all_frozen_units'];independent_ids={r['unit_id'] for r in independent};predictions=[]
    for row in rows:
        f=byid[row['unit_id']]
        predictions.append(dict(unit_id=row['unit_id'],truth=row['label'],truth_family=source_family(row),source_audited_independent=row['unit_id'] in independent_ids,
            homolog_component=row['homolog_component'],genera=';'.join(row['genus_names']),evidence_tiers=';'.join(row['evidence_tiers']),
            accessions=';'.join(sorted({s['accession'] for s in row['sources']})),predicted_family=f['family'],high_identity=row['unit_id'] in chosen,
            eligible=f['eligible'],eligibility_reasons=json.dumps(f['eligibility_reasons']),core_domain_coverage=f['core_domain_coverage'],positive_score=f['positive_score'],
            competitor_score=f['competitor_score'],competitor_margin=f['competitor_margin'],positive_prototype=f['positive_prototype'],competitor_prototype=f['competitor_prototype']))
    write_tsv(out/'predictions.tsv',predictions)
    after=state_hashes()
    if before!=after:raise ValueError('Confirmation state changed during prediction')
    if args.split=='final_validation':
        actual=table(out/'annotation/catalytic_confirmation.tsv');actual_high={r['protein_id'] for r in actual if r['identity_confidence']=='calibrated_high'}
        if actual_high!=chosen:raise ValueError('Pipeline confirmation labels differ from locked rule')
    report=dict(status='confirmation_calibration_complete' if args.split=='calibration' else 'final_confirmation_validation_complete',split=args.split,
        policy_sha256=sha(policy_path),corpus_freeze_sha256=sha(corpus/'freeze.json'),model_manifest_sha256=sha(ROOT/'references/processed/catalytic_confirmation/manifest.json'),
        feature_stage_fingerprint=fingerprint,source_audit=str(audit_path),source_audit_sha256=sha(audit_path),families=rules,metrics=metrics,
        prediction_state_hashes=before,after_prediction_state_hashes=after,units=len(rows),independent_units=len(independent),command=command,
        lock=str(args.lock.resolve()) if args.lock else None,lock_sha256=sha(args.lock) if args.lock else None,
        independent_high_predictions=[p for p in predictions if p['source_audited_independent'] and p['high_identity']],
        limitations=['Source-author/reviewed catalytic identity, not per-protein biochemical activity, complex/function precision or physiology.',
            'Sequence-unit Wilson intervals are descriptive; component representatives and genera remain separately reported.',
            'Public Pfam training membership cannot be reconstructed; own prototype/core/genus separation is checked.',
            'No Hmd or split/noncanonical-center high-confidence calibration; unsupported cases remain explicit.',
            'Final identity confirmation must never be used to select the locked stage again.'],
        artifacts={name:sha(out/name) for name in ['predictions.tsv','annotation/run_manifest.json','annotation/catalytic_confirmation.tsv','annotation/confirmation/features.json','annotation/confirmation/feature_manifest.json']})
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.split=='calibration':
        threshold=dict(status='calibration_locked_before_final_confirmation',policy_sha256=sha(policy_path),model_manifest_sha256=report['model_manifest_sha256'],
            feature_stage_fingerprint=fingerprint,families=rules,calibration_summary=str(out/'summary.json'),calibration_summary_sha256=sha(out/'summary.json'),
            source_audit_sha256=sha(audit_path),scope=policy['scope'])
        rule_path.write_text(json.dumps(threshold,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],families=rules,metrics=metrics),indent=2))


if __name__=='__main__':main()
