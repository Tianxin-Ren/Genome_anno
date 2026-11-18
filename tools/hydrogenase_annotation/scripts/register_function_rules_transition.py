#!/usr/bin/env python3
"""Register only downstream function additions; preserve calibrated identity inputs."""
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'research_workflow')]
from audit_topology_transition import lock_sha,verify_topology_transition
from catalytic_confirmation import stage_fingerprint

HERE=ROOT/'experiments/function_rules_20261112'
CONFIG='config/catalytic_confirmation_thresholds.json'
ACTIVE='experiments/function_rules_20261112/topology_transition_v1.json'
CHANGED=['research_workflow/explain_candidates.py','run_hydrogenase.py']
MODELS=['research_workflow/assay_sequence_identity_model.json','research_workflow/fefe_experimental_backbone_model.json']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def restore_parent_source(name,text):
    if name=='research_workflow/explain_candidates.py':
        for line in [
            "    from function_rules import FIELDS as FUNCTION_RULE_FIELDS, annotate as annotate_function_rules\n",
            "    class_predictions,class_manifest=annotate_function_rules(annotation,out/'function_rules',functional)\n",
            "    for record in functional:record.update(class_predictions[record['candidate_id']])\n",
            "        function_rules_version=class_manifest['version'],\n",
            "        function_rules_manifest_sha256=sha(out/'function_rules/manifest.json'),\n",
            "        function_rules_class_counts=class_manifest['predicted_class_counts'],\n",
        ]:
            assert text.count(line)==1,line
            text=text.replace(line,'')
        assert text.count('*FUNCTION_RULE_FIELDS,')==1
        return text.replace('*FUNCTION_RULE_FIELDS,','')
    for line in [
        "        class_manifest=json.loads((functions/'function_rules/manifest.json').read_text())\n",
        "        manifest.update(function_rules_version=class_manifest['version'],function_class_predictions=class_manifest['predicted_class_counts'])\n",
        "        for name in ['manifest.json','predictions.tsv','evidence.json','member_hits.tsv']:\n",
        "            manifest['output_sha256']['functions/function_rules/'+name]=sha(functions/'function_rules'/name)\n",
    ]:
        assert text.count(line)==1,line
        text=text.replace(line,'')
    return text


def main():
    assert not (ROOT/ACTIVE).exists(),'Successor already exists; do not overwrite a registered record'
    HERE.mkdir(parents=True,exist_ok=True)
    original=(ROOT/CONFIG).read_text();rules=json.loads(original);prior_path=rules['topology_transition']
    prior=json.loads((ROOT/prior_path).read_text())
    assert lock_sha(prior)==prior['lock_sha256']
    fingerprint=stage_fingerprint(ROOT)
    current=hashlib.sha256(json.dumps(fingerprint,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert current==prior['current_stage_fingerprint_sha256'],'Calibrated identity feature stage changed'
    parents={}
    for name in CHANGED:
        restored=restore_parent_source(name,(ROOT/name).read_text())
        assert hashlib.sha256(restored.encode()).hexdigest()==prior['protected_sha256'][name],name
        parents[name]=restored.encode()
    for name,digest in prior['protected_sha256'].items():
        if name not in CHANGED:assert sha(ROOT/name)==digest,name
    target_config=original.replace(prior_path,ACTIVE)
    assert original.count(prior_path)==1
    new_rules=json.loads(target_config)
    assert {k:v for k,v in new_rules.items() if k!='topology_transition'}=={k:v for k,v in rules.items() if k!='topology_transition'}
    edits={CONFIG:target_config};old_sha=sha(ROOT/CONFIG);new_sha=hashlib.sha256(target_config.encode()).hexdigest()
    for name in MODELS:
        text=(ROOT/name).read_text();assert text.count(old_sha)==1,name
        edits[name]=text.replace(old_sha,new_sha)
    originals={name:(ROOT/name).read_bytes() for name in edits}
    parents.update(originals)
    for name,content in parents.items():
        dest=HERE/'parent_baseline'/name;dest.parent.mkdir(parents=True,exist_ok=True)
        with dest.open('xb') as f:f.write(content)
    proof=dict(status='identity_stage_unchanged_downstream_function_additions',parent_lock=prior_path,
        parent_lock_sha256=sha(ROOT/prior_path),identity_stage_fingerprint_sha256=current,
        numerical_thresholds_unchanged=True,identity_models_and_training_unchanged=True,
        original_sources_reconstructed_byte_for_byte={name:hashlib.sha256(parents[name]).hexdigest() for name in CHANGED},
        changed_downstream_sources={name:sha(ROOT/name) for name in CHANGED},
        metadata_only_model_edits=MODELS,
        function_layer_sources={name:sha(ROOT/name) for name in ['research_workflow/function_rules.py','config/function_rules_34.json','references/processed/function_rules/members.faa','research_workflow/candidate_screen.py','run_batch.py','scripts/reorganize_evidence.py','scripts/build_function_rules.py']},
        scope='Adds offline functional class predictions after the unchanged identity stage. Inherited frozen identity replay evidence is not counted as a new accuracy experiment.')
    proof_path=HERE/'identity_preservation_proof.json';proof_path.write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n')
    lock=dict(prior);lock.pop('lock_sha256',None)
    protected=dict(prior['protected_sha256'])
    protected.update({name:sha(ROOT/name) for name in CHANGED})
    protected.update({name:hashlib.sha256(text.encode()).hexdigest() for name,text in edits.items()})
    protected[prior_path]=sha(ROOT/prior_path)
    protected[str(proof_path.relative_to(ROOT))]=sha(proof_path)
    protected[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    lock.update(transition='34-class downstream functional rules; calibrated identity unchanged',protected_sha256=protected,
        parent_function_rules_transition=dict(path=prior_path,sha256=sha(ROOT/prior_path)),
        scope_assertions=[*prior['scope_assertions'],
            'The original identity stage fingerprint, catalytic models and numerical thresholds are unchanged.',
            'Only the downstream function producer and workflow reporting are changed among protected source files.',
            'Functional rules and evidence banks are independently bound by their runtime manifest, outside identity calibration.',
            'Genome/bin association, catalytic identity, architecture, maturation and unmeasured physiological qualifiers remain separate.'])
    lock['lock_sha256']=lock_sha(lock)
    with (ROOT/ACTIVE).open('x') as f:f.write(json.dumps(lock,ensure_ascii=False,indent=2)+'\n')
    try:
        for name,text in edits.items():(ROOT/name).write_text(text)
        assert verify_topology_transition(ROOT,ROOT/ACTIVE,fingerprint,new_rules)
        from assay_sequence_identity import _model
        _model(ROOT)
    except BaseException:
        for name,content in originals.items():(ROOT/name).write_bytes(content)
        (ROOT/ACTIVE).unlink()
        raise
    result=dict(status='verified',active=ACTIVE,identity_stage_unchanged=True,numerical_thresholds_unchanged=True)
    (HERE/'registration_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
