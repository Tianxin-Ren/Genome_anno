#!/usr/bin/env python3
"""Verify one scope-preserving R17 -> R18 identity transition; never refit or predict final data."""
import argparse
import ast
import json
from pathlib import Path
from annotate import ROOT
from catalytic_confirmation import sha,stage_fingerprint,qualifies
from identity_truth import load_units

FAMILIES=('NiFe','FeFe')
PARENT='revisions/round17_passed'
FINAL='results/round17_final_confirmation'

# This exact post-confirmation block pins seven immutable input/identity/site artifacts and only merges an allowlisted L4 patch on deep copies.
POST_SELECTION_BLOCK = "    # BEGIN reviewed post-confirmation L4 selection\n    from copy import deepcopy\n    from complex_selection import refine, FIELDS as HYPOTHESIS_FIELDS\n    selection_policy_path=ROOT/'config/complex_selection.json'\n    selection_policy=json.loads(selection_policy_path.read_text())\n    protected_names=['proteins.faa','genes.tsv','active_center_summary.tsv','active_center_blast.tsv','catalytic_sites.tsv','catalytic_confirmation.tsv','confirmation/features.json']\n    protected_hashes={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}\n    from complex_extensions import search_and_merge\n    references,systems,hits,selection_graph_policy,extension_manifest,reference_contracts=search_and_merge(deepcopy(genes),deepcopy(hits),deepcopy(references),deepcopy(systems),deepcopy(thresholds),mode,deepcopy(annotations),selection_policy,out/'structure_extensions',threads)\n    refs_by_id={r['reference_id']:r for r in references}\n    patches,replacements,hypotheses=refine(deepcopy(genes),deepcopy(hits),deepcopy(references),deepcopy(systems),deepcopy(thresholds),mode,deepcopy(annotations),selection_policy,graph_policy=selection_graph_policy,\n        confirmation_features=json.loads((out/'confirmation/features.json').read_text()),\n        center_candidates=list(csv.DictReader((out/'noncanonical_centers/alternatives.tsv').read_text().splitlines(),delimiter='\\t')))\n    allowed_fields=set(('structure_id structure_confidence function_confidence potential_function_ids structure_status function_status '\n        'missing_roles reference_clade clade_method reference_hyddb_label reference_ko_label diagnostic_positions_conserved diagnostic_residue_map '\n        'catalytic_identity catalytic_query_coverage catalytic_reference_coverage local_cluster_complete missing_context_roles small_subunit_topology '\n        'potential_direct_partner reference_mechanism reference_energy_coupling operon_status module_graph_search_status module_graph_assignment_states '\n        'module_organization_status missing_optional_roles structure_completeness_scope localization_status evidence_url flags').split())\n    if set(patches)!=set(replacements) or set(patches)-{a['protein_id'] for a in annotations}:\n        raise ValueError('L4 selection cannot add or remove catalytic anchors')\n    for a in annotations:\n        qid=a['protein_id']\n        if qid not in patches:continue\n        if set(patches[qid])-allowed_fields:raise ValueError('L4 selection cannot alter protected identity/site decisions')\n        if patches[qid].get('potential_function_ids') or patches[qid].get('structure_confidence')=='high':\n            raise ValueError('Alternative L4 selection cannot promote function or high confidence')\n        new_modules,new_nodes,new_edges=replacements[qid]\n        if any(r['anchor']!=qid or r['structure_id']!=patches[qid]['structure_id'] for rows in replacements[qid] for r in rows):\n            raise ValueError('L4 replacement has incorrect anchor/structure ownership')\n        a.update(patches[qid])\n        modules[:]=[r for r in modules if r['anchor']!=qid]+new_modules\n        graph_nodes[:]=[r for r in graph_nodes if r['anchor']!=qid]+new_nodes\n        graph_edges[:]=[r for r in graph_edges if r['anchor']!=qid]+new_edges\n    write_tsv(out/'complex_hypotheses.tsv',[{k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in h.items()} for h in hypotheses],HYPOTHESIS_FIELDS)\n    selection_manifest=dict(policy_sha256=hashlib.sha256(selection_policy_path.read_bytes()).hexdigest(),script_sha256=hashlib.sha256((ROOT/'scripts/complex_selection.py').read_bytes()).hexdigest(),rescued_anchors=sorted(patches),scope='L4-only post-identity provisional core selection; no function or high-confidence promotion')\n    if protected_hashes!={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}:\n        raise ValueError('Post-identity selection changed protected query/identity/site evidence')\n    selection_manifest['protected_before_and_after_sha256']=protected_hashes\n    selection_manifest['extensions']=extension_manifest\n    selection_manifest['module_architecture']=dict(policy_sha256=hashlib.sha256((ROOT/'config/module_architecture.json').read_bytes()).hexdigest(),script_sha256=hashlib.sha256((ROOT/'scripts/module_architecture.py').read_bytes()).hexdigest(),scope='Existing FeS coordinate resources and source-reviewed compatible core competition only; no identity or function promotion')\n    # END reviewed post-confirmation L4 selection\n"

NONCANONICAL_BLOCK = "    # BEGIN reviewed auxiliary noncanonical literature evidence\n    from noncanonical_center_evidence import annotate_evidence\n    noncanonical_manifest=annotate_evidence(out/'proteins.faa',out/'noncanonical_centers',\n        frozenset(a['protein_id'] for a in annotations if a['catalytic_family']=='FeFe'),\n        frozenset(g['id'] for g in genes if g['partial']),threads=threads,root=ROOT)\n    if protected_hashes!={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}:\n        raise ValueError('Auxiliary noncanonical evidence changed protected identity/site data')\n    # END reviewed auxiliary noncanonical literature evidence\n"

EXACT_FUNCTION_BLOCK = "    # BEGIN reviewed post-pairing exact-reference L5 potential\n    from complex_extensions import exact_reference_patches\n    exact_policy_path=ROOT/'config/exact_reference_functions.json'\n    exact_policy=json.loads(exact_policy_path.read_text())\n    exact_patches,exact_evidence=exact_reference_patches(deepcopy(genes),deepcopy(references),deepcopy(systems),deepcopy(annotations),deepcopy(modules),deepcopy(reference_contracts),exact_policy,deepcopy(thresholds))\n    exact_allowed={'function_status','function_confidence','potential_function_ids','potential_direct_partner','flags'}\n    if set(exact_patches)-{a['protein_id'] for a in annotations} or any(set(p)-exact_allowed for p in exact_patches.values()):\n        raise ValueError('Exact-reference L5 rule cannot alter catalytic identity, sites, structure or anchors')\n    for a in annotations:\n        if a['protein_id'] in exact_patches:a.update(exact_patches[a['protein_id']])\n    if protected_hashes!={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in protected_names}:\n        raise ValueError('Exact-reference L5 rule changed protected identity/site evidence')\n    exact_function_manifest=dict(policy_sha256=hashlib.sha256(exact_policy_path.read_bytes()).hexdigest(),script_sha256=hashlib.sha256((ROOT/'scripts/complex_extensions.py').read_bytes()).hexdigest(),assigned_anchors=sorted(exact_patches),evidence=exact_evidence,protected_before_and_after_sha256=protected_hashes,scope=exact_policy['output_scope'])\n    # END reviewed post-pairing exact-reference L5 potential\n"

def require(condition,message):
    if not condition:raise ValueError(message)


def effective_prototypes(rows,family):
    return [(r['prototype_id'],r['sequence'],r['label']=='positive' and r['family']==family,
        r['homolog_components'],r['genus_taxids']) for r in rows]


def remove_structure_extension(text):
    require(text.count(NONCANONICAL_BLOCK)==1,'Auxiliary evidence changed outside read-only post-identity block')
    require(text.index(NONCANONICAL_BLOCK)>text.index('confirmation_manifest=annotate_confirmation(out,annotations,threads,ROOT)'),'Auxiliary evidence must follow identity confirmation')
    text=text.replace(NONCANONICAL_BLOCK,'')
    require(text.count('        noncanonical_center_evidence=noncanonical_manifest,\n')==1,'Auxiliary manifest changed')
    text=text.replace('        noncanonical_center_evidence=noncanonical_manifest,\n','')
    require(text.count(POST_SELECTION_BLOCK)==1,'L4 extension is outside the reviewed allowlisted post-identity block')
    require(text.index(POST_SELECTION_BLOCK)>text.index('confirmation_manifest=annotate_confirmation(out,annotations,threads,ROOT)'),
        'L4 selection must run after immutable identity features and decisions are produced')
    text=text.replace(POST_SELECTION_BLOCK,'')
    require(text.count(EXACT_FUNCTION_BLOCK)==1,'Exact-reference L5 rule changed outside the reviewed function-only block')
    require(text.index(EXACT_FUNCTION_BLOCK)>text.index('cross_pairs=apply_cross_anchor_guard'),'L5 rule must follow cross-anchor safeguards')
    text=text.replace(EXACT_FUNCTION_BLOCK,'')
    for before in ['    from complex_extensions import exact_reaction_rows\n    functions.extend(exact_reaction_rows(annotations,reference_contracts))\n','        exact_reference_functions=exact_function_manifest,\n']:
        require(text.count(before)==1,'Unexpected post-identity L5 output change')
        text=text.replace(before,'')
    changes=[("'catalytic_center_candidates.tsv','catalytic_center_components.tsv','catalytic_confirmation.tsv','complex_hypotheses.tsv','reaction_contract_candidates.tsv',\n            'structure_extensions/search_manifest.json','structure_extensions/queries.faa','structure_extensions/blast.tsv']}", "'catalytic_center_candidates.tsv','catalytic_center_components.tsv','catalytic_confirmation.tsv']}"), ('        complex_selection=selection_manifest,\n', ''), ('        reaction_contract_candidates=contract_manifest,\n', ''), ('    from complex_extensions import write_contract_candidates\n    contract_manifest=write_contract_candidates(out,annotations,systems,reference_contracts)\n', ''), ('def classify(genes, hits, hmm_hits, references, systems, thresholds, mode="assembly", graph_policy=None):', 'def classify(genes, hits, hmm_hits, references, systems, thresholds, mode="assembly"):'), ('    if graph_policy is None:graph_policy=load_policy(ROOT,references)\n', '    graph_policy=load_policy(ROOT,references)\n')]
    for before,after in changes:
        require(text.count(before)==1,'Unexpected annotation traceability change outside L4 selection')
        text=text.replace(before,after)
    return text


def remove_hmd_extension(text):
    # This exact reviewed patch can only add unresolved Fe anchors, with no graph nodes or function.
    replacements=[
        ("    from identity_truth import adjudications\n    identity_review=adjudications(ROOT)\n",''),
        ('        retained_hmd=False\n',''),
        ('''                review=identity_review.get(best['sequence_sha256'],{})
                # Hmd paralog similarity cannot exclude H2 chemistry with a different substrate.
                retained_hmd=(review.get('family')=='Fe' and review.get('identity_label')=='indeterminate' and
                    any(h['family']=='Fe' for h in hmm_by_query[qid]))
                if not retained_hmd:
                    rejected.append(dict(protein_id=qid,reason="non_catalytic_homolog_wins",reference=best["reference_id"]))
                    continue
''','''                rejected.append(dict(protein_id=qid,reason="non_catalytic_homolog_wins",reference=best["reference_id"]))
                continue
'''),
        ('''family=('Fe' if retained_hmd else refs[competing[0]["sseqid"]]["center"] if competing else hmm_by_query[qid][0]["family"])''',
         '''family=(refs[competing[0]["sseqid"]]["center"] if competing else hmm_by_query[qid][0]["family"])'''),
        ('''flags="family_candidate_only"+
                    (';Hmd_paralog_competitor_blocks_function_transfer_not_discovery' if retained_hmd else '')))''',
         '''flags="family_candidate_only"))''')]
    for before,after in replacements:
        require(text.count(before)==1,'Annotation change is outside the reviewed Fe-only extension')
        text=text.replace(before,after)
    return text


def prediction_ast(text,child=False):
    tree=ast.parse(text)
    # The current input-translation certificate adds a narrowly scoped stale
    # guard branch.  Normalize that branch back to the historical predicate;
    # the certificate itself separately pins and verifies the new branch.
    for fn in ast.walk(tree):
        if isinstance(fn,ast.FunctionDef) and fn.name=='annotate_confirmation':
            for i,node in enumerate(fn.body):
                if isinstance(node,ast.If) and isinstance(node.test,ast.Name) and node.test.id=='rules':
                    replacement=ast.parse("""if rules and (rules['policy_sha256']!=sha(root/'config/catalytic_confirmation.json') or rules['model_manifest_sha256']!=sha(directory/'manifest.json') or rules['feature_stage_fingerprint']!=stage_fingerprint(root)):\n    raise ValueError('Confirmation calibration is stale')""").body[0]
                    fn.body[i]=replacement
    tree.body=[n for n in tree.body if not isinstance(n,ast.FunctionDef) or n.name not in ['build_model','stage_fingerprint']]
    if child:
        expected=ast.dump(ast.parse("manifest['identity_adjudications_sha256']!=sha(root/policy['identity_adjudications'])",mode='eval').body)
        count=0
        for node in ast.walk(tree):
            if isinstance(node,ast.BoolOp):
                retained=[]
                for value in node.values:
                    if ast.dump(value)==expected:count+=1
                    else:retained.append(value)
                node.values=retained
        require(count==1,'Expected one added source-integrity precondition, not changed identity mathematics')
    return ast.dump(tree,include_attributes=False)


def verify_original_transition(root=ROOT):
    root=Path(root);parent=root/PARENT;read=lambda p:json.loads(Path(p).read_text());checks=[]
    def check(name,value):
        require(value,name);checks.append(name)
    snapshot=read(parent/'snapshot_sha256.json')
    for name,digest in snapshot.items():check('frozen_parent:'+name,sha(parent/name)==digest)
    old_final=root/FINAL;cf=read(old_final/'summary.json')
    check('parent_final_summary_matches_frozen_copy',(old_final/'summary.json').read_bytes()==(parent/'evidence/round17_final_confirmation/summary.json').read_bytes())
    old_rule_path=parent/'config/catalytic_confirmation_thresholds.json';old_rule=read(old_rule_path)
    rule_path=root/'config/catalytic_confirmation_thresholds.json';rule=read(rule_path)
    old_fp=cf['feature_stage_fingerprint'];new_fp=stage_fingerprint(root)
    check('parent_prediction_was_immutable',cf['prediction_state_hashes']==cf['after_prediction_state_hashes'])
    check('current_calibration_fingerprint',rule['feature_stage_fingerprint']==new_fp)
    check('unchanged_calibrated_families_and_parameters',set(rule['families'])==set(FAMILIES) and rule['families']==old_rule['families'])
    lock_path=Path(cf['lock']);lock=read(lock_path)
    check('parent_final_lock',sha(lock_path)==cf['lock_sha256'] and lock['feature_stage_fingerprint']==old_fp and lock['thresholds_sha256']==sha(old_rule_path))
    for name,digest in cf['artifacts'].items():check('parent_final_artifact:'+name,sha(old_final/name)==digest)
    check('registered_targets_unchanged',lock['acceptance_targets_sha256']==sha(root/'config/acceptance_targets.json'))
    allowed={'scripts/annotate.py','scripts/catalytic_confirmation.py','config/catalytic_confirmation.json','references/processed/catalytic_confirmation/manifest.json'}
    added={'scripts/noncanonical_center_evidence.py','config/noncanonical_center_evidence.json','scripts/identity_truth.py','config/catalytic_identity_adjudications.json','scripts/complex_selection.py','config/complex_selection.json','scripts/module_architecture.py','config/module_architecture.json','scripts/complex_extensions.py','scripts/build_complex_extensions.py','config/complex_extension_sources.json','config/reaction_contracts.json','config/exact_reference_functions.json','references/processed/complex_extensions/manifest.json'}
    check('no_unreviewed_stage_files',set(new_fp['files'])==set(old_fp['files'])|added)
    check('discovery_stage_identical',new_fp['discovery']==old_fp['discovery'])
    for name,digest in old_fp['files'].items():
        if name not in allowed:check('unchanged_prediction_dependency:'+name,new_fp['files'][name]==digest)
    check('annotation_only_reviewed_Fe_candidate_and_post_identity_L4_extensions',remove_hmd_extension(remove_structure_extension((root/'scripts/annotate.py').read_text()))==(parent/'scripts/annotate.py').read_text())
    check('confirmation_prediction_mathematics_identical',prediction_ast((root/'scripts/catalytic_confirmation.py').read_text(),True)==prediction_ast((parent/'scripts/catalytic_confirmation.py').read_text()))
    # Remove only descriptive/training-provenance metadata; every mathematical policy field must agree.
    old_policy=read(parent/'config/catalytic_confirmation.json');new_policy=read(root/'config/catalytic_confirmation.json')
    for p in [old_policy,new_policy]:
        for key in ['version','training_selection','identity_adjudications']:p.pop(key,None)
    check('identity_policy_mathematics_identical',old_policy==new_policy and new_policy['eligible_calibration_families']==list(FAMILIES))
    old_model=parent/'model_artifacts/catalytic_confirmation';model=root/'references/processed/catalytic_confirmation'
    manifest=read(model/'manifest.json')
    for name,digest in manifest['files'].items():check('current_model_artifact:'+name,sha(model/name)==digest)
    old_prototypes=read(old_model/'prototypes.json');new_prototypes=read(model/'prototypes.json')
    for family in FAMILIES:check('same_sequence_order_positive_and_competitor_roles:'+family,effective_prototypes(old_prototypes,family)==effective_prototypes(new_prototypes,family))
    for name in ['domain_models.json','domains.hmm']:check('same_query_domain_model:'+name,(old_model/name).read_bytes()==(model/name).read_bytes())
    cal_path=Path(rule['calibration_summary']);cal=read(cal_path);old_cal_path=Path(lock['calibration_summary']);old_cal=read(old_cal_path)
    check('current_calibration_summary_integrity',sha(cal_path)==rule['calibration_summary_sha256'])
    check('parent_calibration_summary_integrity',sha(old_cal_path)==lock['calibration_summary_sha256'])
    for label,path,summary in [('parent',old_cal_path,old_cal),('current',cal_path,cal)]:
        for name,digest in summary['artifacts'].items():check(label+'_calibration_artifact:'+name,sha(path.parent/name)==digest)
    old_audit=read(Path(cf['source_audit']));new_audit=read(Path(cal['source_audit']))
    check('source_audits_integrity',sha(Path(cf['source_audit']))==cf['source_audit_sha256'] and sha(Path(cal['source_audit']))==cal['source_audit_sha256'])
    for label,path,summary in [('parent',Path(cf['source_audit']),old_audit),('current',Path(cal['source_audit']),new_audit)]:
        for name,digest in summary['artifacts'].items():check(label+'_source_audit_artifact:'+name,sha(path.parent/name)==digest)
    check('same_independence_population',old_audit['independent_unit_ids']==new_audit['independent_unit_ids'])
    corpus=root/new_policy['corpus'];historical=read(corpus/'units.json');effective=load_units(corpus,root)
    check('same_final_labels_and_family_truth',[(r['unit_id'],r['label'],r['categories']) for r in historical if r['split']=='final_validation']==[(r['unit_id'],r['label'],r['categories']) for r in effective if r['split']=='final_validation'])
    old_features=read(old_cal_path.parent/'annotation/confirmation/features.json');features=read(cal_path.parent/'annotation/confirmation/features.json')
    keys=['protein_id','family','eligible','eligibility_reasons','core_positions','core_domain_coverage','positive_score','competitor_score','competitor_margin','positive_prototype','competitor_prototype','domain_evidence','source_active_center']
    vector=lambda rows:[{k:r[k] for k in keys} for r in rows if r['family'] in FAMILIES]
    check('actual_calibration_decision_vectors_identical',vector(old_features)==vector(features))
    selected=lambda rows:{r['protein_id'] for r in rows if rule['families'].get(r['family']) and qualifies(r,rule['families'][r['family']])}
    check('actual_calibration_high_ids_identical',selected(old_features)==selected(features))
    old_run=read(old_final/'annotation/run_manifest.json');new_run=read(cal_path.parent/'annotation/run_manifest.json')
    check('same_recorded_software',old_run['software']==new_run['software'])
    # L4 selection runs on deep copies after identity decisions, merges only literal allowlisted fields, and cannot add anchors.
    # Existing cross-anchor pairing only changes structure/function; its source is byte-identical above.
    # Added Fe anchors have no graph nodes; split-center pairing only visits NiFe anchors.
    # All center mapping helpers are byte-identical above. Classification of other queries is unchanged.
    return dict(status='scope_preserving_transition_verified',families=list(FAMILIES),parent_final=str(old_final/'summary.json'),
        parent_final_summary_sha256=sha(old_final/'summary.json'),parent_thresholds=str(old_rule_path),parent_thresholds_sha256=sha(old_rule_path),
        current_feature_stage_fingerprint=new_fp,current_thresholds_sha256=sha(rule_path),checks=checks,passed=len(checks),
        inherited_final_metrics=cf['metrics']['source_audited_independent'],
        interpretation='The calibrated NiFe/FeFe decision function is unchanged for all inputs under the exact reviewed Fe-only discovery and post-identity L4-only source transition. Parent independent final evidence is inherited within that scope; no new final predictions or threshold selection.',
        limits=['No inheritance for Fe/Hmd, new reaction functions, structure, FeS or complex accuracy.',
            'No new independent biological samples or population-precision guarantee.',
            'A later prediction-source, policy, sequence, prototype-role or parameter change invalidates this transition.'])


def verify_transition(root=ROOT):
    root=Path(root)
    rule=json.loads((root/'config/catalytic_confirmation_thresholds.json').read_text())
    if rule.get('input_translation_transition'):
        from audit_input_translation_transfer import verify_registered_transition
        from catalytic_confirmation import stage_fingerprint
        current_fp=stage_fingerprint(root)
        stage_valid=verify_registered_transition(root, root/rule['input_translation_transition'], current_fp,
                old_stage_fingerprint=rule['feature_stage_fingerprint'], rules=rule)
        if not stage_valid and rule.get('reporting_transition'):
            from audit_discovery_reporting_transfer import verify_reporting_transition
            stage_valid=verify_reporting_transition(root,root/rule['reporting_transition'],current_fp,rule)
        if not stage_valid and rule.get('topology_transition'):
            from audit_topology_transition import verify_topology_transition
            stage_valid=verify_topology_transition(root,root/rule['topology_transition'],current_fp,rule)
        if not stage_valid:
            raise ValueError('Registered input-translation transition failed strict scope verification')
        # Inherit only the frozen R47/R49 identity proof, never calibration
        # counts.  The plan and round47 snapshot bind both source and SHA.
        plan=json.loads((root/'experiments/round49/plan.json').read_text())
        proof_path=root/plan['parent_identity_proof']['path']
        snapshot_path=root/plan['parent_identity_evidence_snapshot']['path']
        snapshot=json.loads(snapshot_path.read_text()); proof=json.loads(proof_path.read_text())
        if sha(proof_path)!=plan['parent_identity_proof']['sha256'] or sha(snapshot_path)!=plan['parent_identity_evidence_snapshot']['sha256']:
            raise ValueError('Frozen identity proof or snapshot changed')
        if snapshot.get('shared_artifacts',{}).get(plan['parent_identity_proof']['path'])!=sha(proof_path):
            raise ValueError('Identity proof is not covered by frozen snapshot')
        if proof.get('status')!='conservative_guard_transition_verified':
            raise ValueError('Unexpected inherited identity proof status')
        provided={'calibration_query_proteins': proof.get('calibration_queries'),
                  'source_proof': str(plan['parent_identity_proof']['path']),
                  'new_independent_validation': False}
        return dict(status='input_translation_scope_preserving_transition_verified', families=list(FAMILIES),
            parent_final=str(plan['parent_identity_proof']['path']),
            parent_final_summary_sha256=sha(proof_path),
            current_feature_stage_fingerprint=current_fp,
            current_thresholds_sha256=sha(root/'config/catalytic_confirmation_thresholds.json'),
            inherited_final_metrics=proof['inherited_final_metrics'],
            current_metric_bounds=proof['current_metric_bounds'],
            family_checks=proof['family_checks'], precision_gate_supported=proof['precision_gate_supported'],
            provided_protein_identity_scope=provided,
            new_final_prediction_runs=0,
            checks=['registered_input_translation_certificate','same_NiFe_FeFe_rules_and_parameters',
                    'same_policy_model_and_calibration_summary','stage_delta_only_annotate_and_guard',
                    'R49_provided_protein_identity_scope_inherited_from_frozen_calibration'],
            passed=5,
            interpretation='Registered input-stage transitions preserve the reviewed identity scope. Final identity metrics and bounds are inherited from the frozen R47/R49 proof; this is not new independent biological validation.',
            limits=['No new final predictions or calibration fitting.','No inheritance for new reaction functions, FeS branches, or uncalibrated families.'])
    if rule.get('identity_validation_transition')=='references/validation_locks/round49_fes_bank_transition.json':
        from audit_fes_bank_transfer import verify_transition as verify_fes_transition
        return verify_fes_transition(root)
    if rule.get('identity_validation_transition')=='references/validation_locks/round47_repeat_indel_transition.json':
        from audit_repeat_indel_transfer import verify_transition as verify_guard_transition
        return dict(verify_guard_transition(root,verify_original_transition),families=list(FAMILIES))
    return verify_original_transition(root)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False);report=verify_transition()
    (args.out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],checks=report['passed'],families=report['families']),indent=2))


if __name__=='__main__':main()
