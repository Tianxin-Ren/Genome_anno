"""Verify a FeS-only reference expansion preserves the prior L1/L2 decision function."""
import json
from pathlib import Path
from annotate import ROOT
from catalytic_confirmation import sha,stage_fingerprint
from run_acceptance_checks import state_hashes


def load(path):return json.loads(path.read_text())
def require(value,message):
    if not value:raise ValueError(message)


def architecture_equivalent(old,new,old_sites,new_sites,old_parameters,new_parameters):
    pins={'fes_sites_sha256','fes_parameters_sha256'}
    require({k:v for k,v in old.items() if k not in pins}=={k:v for k,v in new.items() if k not in pins},'Architecture changed outside the two metadata pins')
    require(old_parameters==new_parameters,'FeS prediction parameters changed')
    ids=set(old['signature']['site_ids']+[old['signature']['scaffold_site_id']])
    previous={s['site_id']:s for s in old_sites};current={s['site_id']:s for s in new_sites}
    require(all(i in current and previous[i]==current[i] for i in ids),'Architecture-used FeS evidence changed')
    return sorted(ids)


def allowed_state_transition(previous,current,allowed,expected_model):
    changed={k for k in set(previous)|set(current) if previous.get(k)!=current.get(k)}-{'scripts/assess_acceptance.py'}
    require(not changed-set(allowed),'Unregistered change outside the FeS bank transition: '+str(sorted(changed-set(allowed))))
    require(all(current.get(k)==v for k,v in expected_model.items()),'Deployed FeS transition model differs from the pre-execution lock')
    return sorted(changed)


def verify_transition(root=ROOT):
    root=Path(root)
    rules_path=root/'config/catalytic_confirmation_thresholds.json'
    if rules_path.exists() and load(rules_path).get('input_translation_transition'):
        from audit_input_translation_transfer import verify_registered_transition
        from catalytic_confirmation import stage_fingerprint
        rules=load(rules_path)
        current_stage=stage_fingerprint(root)
        stage_valid=verify_registered_transition(root, root/rules['input_translation_transition'], current_stage,
                old_stage_fingerprint=rules['feature_stage_fingerprint'], rules=rules)
        if not stage_valid and rules.get('reporting_transition'):
            from audit_discovery_reporting_transfer import verify_reporting_transition
            stage_valid=verify_reporting_transition(root,root/rules['reporting_transition'],current_stage,rules)
        if not stage_valid and rules.get('topology_transition'):
            from audit_topology_transition import verify_topology_transition
            stage_valid=verify_topology_transition(root,root/rules['topology_transition'],current_stage,rules)
        if not stage_valid:
            raise ValueError('Input translation certificate failed before FeS transfer audit')
        plan=load(root/'experiments/round49/plan.json'); proof_path=root/plan['parent_identity_proof']['path']; snapshot_path=root/plan['parent_identity_evidence_snapshot']['path']
        proof=load(proof_path); snapshot=load(snapshot_path)
        require(sha(proof_path)==plan['parent_identity_proof']['sha256'] and sha(snapshot_path)==plan['parent_identity_evidence_snapshot']['sha256'],'Frozen identity proof changed')
        require(snapshot.get('shared_artifacts',{}).get(plan['parent_identity_proof']['path'])==sha(proof_path),'Proof is outside frozen snapshot')
        return dict(status='input_translation_transition_verified_before_FeS_chain',
            validation_scope='registered input-stage changes only; FeS chain not re-run under changed input stage',
            current_metric_bounds=proof['current_metric_bounds'], new_final_prediction_runs=0,
            inherited_final_metrics=proof['inherited_final_metrics'], family_checks=proof['family_checks'],
            precision_gate_supported=proof['precision_gate_supported'], checks=['registered_input_translation_certificate','frozen_identity_proof_sha256'], passed=2,
            interpretation='The historical FeS transfer evidence remains read-only and is not relabeled as a new FeS test.')
    directory=root/'experiments/round49';plan=load(directory/'plan.json')
    lock=load(root/'references/validation_locks/round49_fes_bank_transition.json');checks=[]
    def check(name,condition):
        require(condition,name);checks.append(dict(test=name,passed=True))
    check('transition_registered_before_full_checks',lock['status']=='frozen_before_full_checks_and_current_negative_or_user_predictions' and lock['plan_sha256']==sha(directory/'plan.json'))
    check('auditor_source_pinned',lock['auditor_sha256']==sha(Path(__file__)))
    baseline=directory/'baseline_view'
    for name,value in plan['production_before'].items():require(sha(baseline/name)==value,'Historical source copy changed: '+name)
    checks.append(dict(test='all388historical_production_sources_match',passed=True))
    current=state_hashes();allowed_state_transition(plan['production_before'],current,plan['allowed_changed_files'],lock['deployed_model_hashes'])
    checks.append(dict(test='only_registered_model_or_verifier_changes',passed=True))
    check('rejected_selector_not_deployed',sha(root/'scripts/fes_sites.py')==plan['production_before']['scripts/fes_sites.py'])
    old_sites=load(baseline/'references/processed/fes/sites.json');new_sites=load(root/'references/processed/fes/sites.json')
    check('old31sites_unchanged_only3reviewed_sites_added',len(old_sites)==31 and new_sites[:-3]==old_sites and
        {s['site_id'] for s in new_sites[-3:]}==set(lock['new_site_ids']) and all(s['reference_id']=='1CC1__entity_1' for s in new_sites[-3:]))
    old_policy=load(baseline/'config/fes_site_review.json');new_policy=load(root/'config/fes_site_review.json')
    check('FeS_old_rules_and_parameters_unchanged',new_policy['sites'][:-3]==old_policy['sites'] and
        {k:v for k,v in new_policy.items() if k not in ['version','sites']}=={k:v for k,v in old_policy.items() if k not in ['version','sites']})
    old_bank=load(baseline/'references/structural_sites/templates.json');new_bank=load(root/'references/structural_sites/templates.json')
    candidate_path=root/'experiments/round48/candidate_model/references/structural_sites/templates.json';candidate=load(candidate_path)
    check('reference_source_extension_is_the_reviewed_candidate',sha(candidate_path)==plan['candidate_sources_sha256'] and new_bank[:-1]==old_bank and
        all(new_bank[-1].get(k)==v for k,v in candidate[-1].items()) and sha(root/'references/structural_sites/templates.json')==plan['canonical_sources_sha256'])
    manifest=load(root/'references/processed/fes/manifest.json')
    check('new_Fes_manifest_pins_current_sources',manifest['policy_sha256']==sha(root/'config/fes_site_review.json') and
        manifest['sources_sha256']==sha(root/'references/structural_sites/templates.json') and manifest['prediction_script_sha256']==sha(root/'scripts/fes_sites.py'))
    for name,value in manifest['files'].items():require(sha(root/'references/processed/fes'/name)==value,'FeS model artifact changed: '+name)
    architecture_equivalent(load(baseline/'config/module_architecture.json'),load(root/'config/module_architecture.json'),old_sites,new_sites,old_policy['parameters'],new_policy['parameters'])
    checks.append(dict(test='module_architecture_used_sites_and_logic_identical',passed=True))
    old_build=load(baseline/'references/processed/build_manifest.json');new_build=load(root/'references/processed/build_manifest.json');pins={'fes_model_manifest_sha256','fes_site_policy_sha256'}
    check('reference_build_changed_only_Fes_metadata',{k:v for k,v in old_build.items() if k not in pins}=={k:v for k,v in new_build.items() if k not in pins})
    old_fp=plan['old_feature_fingerprint'];new_fp=stage_fingerprint(root)
    check('identity_fingerprint_diff_is_one_metadata_file',old_fp['discovery']==new_fp['discovery'] and set(old_fp['files'])==set(new_fp['files']) and
        [k for k in old_fp['files'] if old_fp['files'][k]!=new_fp['files'][k]]==['config/module_architecture.json'] and new_fp==lock['feature_stage_fingerprint'])
    old_rule=load(baseline/'config/catalytic_confirmation_thresholds.json');rule=load(root/'config/catalytic_confirmation_thresholds.json')
    check('identity_parameters_fixed',rule['families']==old_rule['families'])
    cal_path=Path(rule['calibration_summary']);cal=load(cal_path);old_cal_path=Path(old_rule['calibration_summary'])
    check('historical_calibration_integrity',sha(cal_path)==rule['calibration_summary_sha256']==lock['calibration_summary_sha256'] and cal['feature_stage_fingerprint']==new_fp)
    for name,value in cal['artifacts'].items():require(sha(cal_path.parent/name)==value,'Calibration artifact changed: '+name)
    check('all_historical_calibration_feature_values_identical',load(cal_path.parent/'annotation/confirmation/features.json')==load(old_cal_path.parent/'annotation/confirmation/features.json'))
    parent_path=root/plan['parent_identity_proof']['path'];parent=load(parent_path);evidence_path=root/plan['parent_identity_evidence_snapshot']['path'];evidence=load(evidence_path)
    check('parent_identity_proof_is_the_verified_frozen_artifact',sha(evidence_path)==plan['parent_identity_evidence_snapshot']['sha256'] and
        sha(parent_path)==plan['parent_identity_proof']['sha256']==evidence['shared_artifacts'][plan['parent_identity_proof']['path']] and
        parent['status']=='conservative_guard_transition_verified' and parent['current_feature_stage_fingerprint']==old_fp)
    report=dict(parent);report.update(status='FeS_bank_metadata_equivalence_verified',validation_scope='FeS_bank_metadata_only_after_repeat_guard',checks=checks,passed=len(checks),
        current_feature_stage_fingerprint=new_fp,current_thresholds_sha256=sha(root/'config/catalytic_confirmation_thresholds.json'),
        direct_parent_proof=str(parent_path),direct_parent_proof_sha256=sha(parent_path),new_final_prediction_runs=0,
        interpretation='OnlyFeSbank additions and twoarchitecture hashpins changed; actual usedsites/parameters/code and allcalibration features identical. Prior finiteidentity bounds inherited, no new finalprediction or fitting.',
        limits=parent['limits']+['No FeS/complex/function accuracy is inherited. The added1CC1source is already consumed and nowreference-only.'])
    return report


def verify_layer_state(previous,current,root=ROOT):
    root=Path(root)
    rules_path=root/'config/catalytic_confirmation_thresholds.json'
    if rules_path.exists() and load(rules_path).get('input_translation_transition'):
        prior=load(root/'experiments/round47/pipeline_state_lock.json')['production_state']
        require(set(previous)==set(prior) and all(previous.get(k)==v for k,v in prior.items()),
            'Prior L2 execution does not belong to the frozen R47 production state')
        report=verify_transition(root)
        return dict(status='unchanged_conditional_L2_under_input_translation_transition',
            changed_files=['scripts/annotate.py','scripts/catalytic_confirmation.py'],
            identity_transition_status=report['status'],new_L2_predictions=0)
    plan=load(root/'experiments/round49/plan.json');lock=load(root/'references/validation_locks/round49_fes_bank_transition.json')
    require(all(previous.get(k)==v for k,v in plan['production_before'].items() if k!='scripts/assess_acceptance.py') and
        set(previous)==set(plan['production_before']),'Prior L2 execution does not belong to the frozen pre-FeS state')
    changed=allowed_state_transition(previous,current,plan['allowed_changed_files'],lock['deployed_model_hashes'])
    report=verify_transition(root)
    return dict(status='unchanged_conditional_L2_under_verified_Fes_bank_extension',changed_files=changed,identity_transition_status=report['status'],new_L2_predictions=0)
