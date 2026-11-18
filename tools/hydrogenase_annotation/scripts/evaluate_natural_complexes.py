#!/usr/bin/env python3
"""Evaluate frozen natural-context complexes and directly assayed reaction contracts."""
import argparse
import json
from pathlib import Path
from annotate import ROOT, annotate, write_tsv
from catalytic_confirmation import sha, table
from run_acceptance_checks import state_hashes


def compare_case(case, anchors, members, complete_statuses, reaction_keys=None):
    truth = {m['query_id'] for m in case['members']}
    relevant = [a for a in anchors if a['protein_id'] in truth]
    complete = [a for a in relevant if a['structure_status'] in complete_statuses]
    exact = []
    predicted_sets = {}
    for a in complete:
        chosen = {m['protein_id'] for m in members if m['anchor'] == a['protein_id'] and str(m['required_for_structure']) == 'True'}
        predicted_sets[a['protein_id']] = sorted(chosen)
        if case['hydrogenase'] and a['protein_id'] == case['query_anchor'] and chosen == truth:
            exact.append(a)
    outputs = [(a, f) for a in relevant for f in a.get('potential_function_ids', '').split(';') if f]
    reaction_keys=reaction_keys or {}
    key=lambda f:reaction_keys.get(f,f)
    compatible={key(f) for f in case['compatible_function_ids']}
    contradicted={key(f) for f in case.get('contradicted_function_ids',[])}
    correct_outputs = [(a, f) for a, f in outputs if a in exact and key(f) in compatible]
    # One assayed forward capability per case; unassayed extra directions are retained as ungradable.
    # The scored object is a member-set/reaction pair; attaching a known reaction to a wrong core is also an error.
    incompatible = [(a, f) for a, f in outputs if not case['hydrogenase'] or key(f) in contradicted or
                    (key(f) in compatible and a not in exact)]
    ungraded = [(a, f) for a, f in outputs if case['hydrogenase'] and key(f) not in compatible and key(f) not in contradicted]
    return dict(case_id=case['id'], source=case['source'], hydrogenase=case['hydrogenase'],
        expected_members=sorted(truth), expected_anchor=case['query_anchor'],
        anchor_found=any(a['protein_id'] == case['query_anchor'] for a in relevant) if case['hydrogenase'] else None,
        relevant_anchors=[{k: a.get(k, '') for k in ['protein_id', 'structure_id', 'catalytic_status', 'catalytic_identity_confidence',
            'structure_status', 'structure_confidence', 'function_status', 'potential_function_ids', 'flags', 'missing_roles']} for a in relevant],
        predicted_complete_member_sets=predicted_sets,
        core_TP=int(bool(exact)), core_FN=int(case['hydrogenase'] and not exact),
        core_FP=sum(a not in exact for a in complete),
        high_core_TP=int(any(a.get('structure_confidence') == 'high' for a in exact)),
        reaction_TP=int(bool(correct_outputs)), reaction_FN=int(case['hydrogenase'] and not correct_outputs),
        reaction_FP=len(incompatible),
        ungraded_reaction_assignments=len(ungraded),
        reaction_output_count=len(outputs),
        reaction_outputs=[dict(anchor=a['protein_id'], function_id=f, matches_assayed_contract=key(f) in compatible) for a, f in outputs],
        contrast_hydrogenase_function_assignments=len(outputs) if not case['hydrogenase'] else 0,
        context='experimental_native_source_development_replay_not_independent_test')


def counts(rows, prefix):
    tp, fp, fn = [sum(r[prefix + '_' + k] for r in rows) for k in ['TP', 'FP', 'FN']]
    from validate_fes_prediction import wilson
    return dict(TP=tp, FP=fp, FN=fn, precision=tp / (tp + fp) if tp + fp else None,
        recall=tp / (tp + fn) if tp + fn else None,
        precision_wilson95=wilson(tp, tp + fp) if tp + fp else None,
        recall_wilson95=wilson(tp, tp + fn) if tp + fn else None)


def verify_corpus(corpus, freeze):
    for name, digest in freeze['artifacts'].items():
        if sha(corpus / name) != digest: raise ValueError('Corpus changed: ' + name)
    for name, digest in freeze['sources'].items():
        if sha(ROOT / name) != digest: raise ValueError('Primary source changed: ' + name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); corpus = args.corpus.resolve()
    freeze = json.loads((corpus / 'freeze.json').read_text())
    if freeze['status'] != 'development_truth_frozen_before_annotation':
        raise ValueError('Frozen natural-context truth required')
    verify_corpus(corpus, freeze)
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    policy = json.loads((ROOT / 'config/complex_function_validation.json').read_text())
    if sha(ROOT / 'config/complex_function_validation.json') != freeze['protocol_sha256']:
        raise ValueError('Protocol changed')
    from complex_extensions import load_contracts
    reaction_keys={c['contract_id']:c['reaction_key'] for c in load_contracts()}
    before = state_hashes(); sources = json.loads((corpus / 'sources.json').read_text())
    truth = json.loads((corpus / 'truth.json').read_text()); rows, runs, ungraded = [], [], []
    execution = dict(status='running', before_hashes=before, corpus_freeze_sha256=sha(corpus / 'freeze.json'), runs=runs)
    def save(): (out / 'execution.json').write_text(json.dumps(execution, indent=2) + '\n')
    save()
    for source, item in sources.items():
        directory = out / source; directory.mkdir()
        annotations = annotate(item['anonymous_genes'], directory, threads=2, mode=item['mode'], input_files=[ROOT / item['input_file']])
        cases = [c for c in truth if c['source'] == source]
        members = table(directory / 'complex_members.tsv')
        manifest=json.loads((directory/'run_manifest.json').read_text())
        replays=set(manifest.get('complex_selection',{}).get('extensions',{}).get('reference_replay_cases',[]))
        for c in cases:
            comparison=compare_case(c,annotations,members,policy['complete_statuses'],reaction_keys)
            comparison['extension_reference_replay']=c['id'] in replays
            rows.append(comparison)
        graded = {m['query_id'] for c in cases for m in c['members']}
        ungraded.extend(dict(source=source, **a) for a in annotations if a['protein_id'] not in graded)
        runs.append(dict(source=source, proteins=item['original_genes'], anchors=len(annotations),
            run_manifest_sha256=sha(directory / 'run_manifest.json')))
        save(); print(source, item['original_genes'], len(annotations), flush=True)
    verify_corpus(corpus, freeze)
    after = state_hashes()
    if before != after: raise ValueError('Code, models or validation policy changed during evaluation')
    metrics = {p: counts(rows, p) for p in ['core', 'reaction']}
    contrast_fp = sum(r['contrast_hydrogenase_function_assignments'] for r in rows)
    target = policy['targets']
    passed = (metrics['core']['precision'] is not None and metrics['core']['precision'] >= target['core_precision'] and
        metrics['core']['recall'] >= target['core_recall'] and metrics['reaction']['precision'] is not None and
        metrics['reaction']['precision'] >= target['reaction_precision'] and metrics['reaction']['recall'] >= target['minimum_reaction_coverage'] and
        contrast_fp <= target['non_H2_complex_false_assignments'])
    (out / 'case_comparisons.json').write_text(json.dumps(rows, indent=2) + '\n')
    write_tsv(out / 'case_comparisons.tsv', [{k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in r.items()} for r in rows])
    write_tsv(out / 'ungraded_anchors.tsv', ungraded)
    summary = dict(status='development_natural_context_evaluated', development_targets_met=passed,
        metrics=metrics, contrast_hydrogenase_function_assignments=contrast_fp, cases=len(rows),
        extension_reference_replay_case_ids=[r['case_id'] for r in rows if r['extension_reference_replay']],
        metrics_by_exposure={name:{p:counts([r for r in rows if r['extension_reference_replay']==is_replay],p) for p in ['core','reaction']}
            for name,is_replay in [('extension_reference_replay',True),('other_exposed_development',False)]},
        ungraded_reaction_assignments=sum(r['ungraded_reaction_assignments'] for r in rows),
        high_exact_core_TP=sum(r['high_core_TP'] for r in rows), runs=runs, ungraded_anchors=len(ungraded),
        targets=target, before_hashes=before, after_hashes=after, corpus_freeze_sha256=sha(corpus / 'freeze.json'),
        limits=freeze['limits'], independent_validation=False)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    execution.update(status='completed', after_hashes=after, summary_sha256=sha(out / 'summary.json')); save()
    print(json.dumps({k: summary[k] for k in ['status', 'development_targets_met', 'metrics', 'contrast_hydrogenase_function_assignments', 'cases', 'ungraded_anchors']}, indent=2))


if __name__ == '__main__': main()
