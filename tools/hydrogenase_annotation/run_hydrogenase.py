#!/usr/bin/env python3
"""Retain protein/DNA candidates and evaluate claim-specific functional evidence."""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_table(path):
    with path.open() as handle: return list(csv.DictReader(handle, delimiter='\t'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def preferred_summary(rows, output):
    """Return auditable preferred-prediction statistics for workflow.json."""
    return dict(
        preferred_function_predictions=dict(Counter(r.get('preferred_function', '') for r in rows)),
        preferred_prediction_confidence=dict(Counter(r.get('prediction_confidence', '') for r in rows)),
        preferred_prediction_output=str(Path(output)),
    )


def record_translation_hypotheses(dna, out, manifest):
    from research_workflow.translation_hypotheses import audit

    output = out/'translation_hypotheses.json'
    result = audit(dna, output)
    split = result.get('split_component_translation', {})
    dna_manifest = json.loads((dna / 'dna_rescue_manifest.json').read_text())
    manifest.update(DNA_rescue_span_count=dna_manifest.get('distinct_native_span_count', 0),
                    translation_hypotheses_status=result['status'],
                    translation_hypotheses_count=len(result['hypotheses']),
                    split_component_translation_count=len(split.get('hypotheses', [])),
                    relay_component_translation_count=len(split.get('upstream_relay_translation', {}).get('hypotheses', [])))
    manifest['output_sha256'][output.name] = sha(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    for flag in ['genbank', 'proteins', 'annotation']: group.add_argument('--'+flag, type=Path)
    parser.add_argument('--genome', type=Path)
    parser.add_argument('--gff', type=Path); parser.add_argument('--bin-map', type=Path)
    parser.add_argument('--mode', choices=['genome','assembly','bin','transcript'], default='assembly')
    parser.add_argument('--out', type=Path, required=True); parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--translation-table', type=int, default=11)
    parser.add_argument('--topology', choices=['auto','linear','circular'], default='auto',
                        help='Explicit replicon topology forwarded to annotation; auto preserves source metadata')
    args = parser.parse_args()
    if not any([args.genome, args.genbank, args.proteins, args.annotation]):
        parser.error('one of --genome, --genbank, --proteins, or --annotation is required')
    if args.genome and (args.genbank or args.annotation):
        parser.error('--genome can only be used standalone or with --proteins')
    if args.threads < 1: parser.error('--threads must be positive')
    if args.annotation and (args.gff or args.bin_map): parser.error('Saved annotation retains its original input scope; omit --gff/--bin-map')
    out = args.out.resolve()
    if out.exists(): parser.error('Use a new output directory; existing results are never overwritten')
    out.mkdir(parents=True)
    local_py = ROOT/'.env/bin/python'
    python = str(local_py) if local_py.is_file() else sys.executable
    env = dict(os.environ)
    env['PATH'] = str(Path(python).parent)+os.pathsep+env.get('PATH','')
    annotation = args.annotation.resolve() if args.annotation else out/'annotation'
    domains = out/'domains'; functions = out/'functions'; fes = out/'fes_refinement'; candidates = out/'candidates'; confirmation = out/'candidate_confirmation'; dna = out/'dna_rescue'
    commands = []
    if not args.annotation:
        flag = next(k for k in ['genbank','proteins','genome'] if getattr(args,k))
        command = [python,str(ROOT/'scripts/annotate.py'),'--'+flag,str(getattr(args,flag).resolve()),'--mode',args.mode,
            '--out',str(annotation),'--threads',str(args.threads),'--translation-table',str(args.translation_table),
            '--topology',args.topology]
        if flag == 'proteins' and args.genome:
            command += ['--genome', str(args.genome.resolve())]
        for option in ['gff','bin_map']:
            if getattr(args,option): command += ['--'+option.replace('_','-'),str(getattr(args,option).resolve())]
        commands.append(('annotation',command))
    commands += [('domains',[python,str(ROOT/'scripts/query_domain_organization.py'),'--runs',str(annotation),'--out',str(domains)]),
        ('candidates',[python,str(ROOT/'research_workflow/candidate_screen.py'),'--annotation',str(annotation),'--domains',str(domains),'--out',str(candidates)]),
        ('dna_rescue',[python,str(ROOT/'research_workflow/dna_rescue.py'),'--annotation',str(annotation),'--screen',str(candidates),'--out',str(dna),'--threads',str(args.threads)]),
        ('candidate_confirmation',[python,str(ROOT/'research_workflow/candidate_confirmation.py'),'--annotation',str(annotation),'--screen',str(candidates),'--out',str(confirmation),'--threads',str(args.threads)]),
        ('functions',[python,str(ROOT/'research_workflow/explain_candidates.py'),'--annotation',str(annotation),'--domain-evidence',str(domains),'--candidate-screen',str(candidates),'--candidate-confirmation',str(confirmation),'--out',str(functions)]),
        ('fes_refinement',[python,str(ROOT/'research_workflow/fes_refinement.py'),str(annotation),str(fes),'--threads',str(args.threads)])]
    manifest = dict(status='running',annotation_reused=bool(args.annotation),annotation=str(annotation),stages=[],
        script_sha256=sha(Path(__file__)),scope='Research predictions with claim-specific uncertainty; no measured query activity or new calibration.')
    path = out/'workflow.json'
    def save(): path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    save()
    try:
        for name, command in commands:
            print(name, flush=True)
            with (out/(name+'.log')).open('w') as log:
                result = subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
            manifest['stages'].append(dict(stage=name,command=command,exit_code=result.returncode)); save()
            if result.returncode:
                from Bio import SeqIO
                predicted = annotation/'prodigal.faa'
                no_predicted_protein = not predicted.exists() or not any(SeqIO.parse(predicted,'fasta'))
                if name != 'annotation' or not args.genome or not no_predicted_protein:
                    raise RuntimeError(f'{name} failed; see {out/(name+".log")}')
                # A failed/empty gene call must not prevent a direct DNA search.
                command=[python,str(ROOT/'research_workflow/dna_rescue.py'),'--genome',str(args.genome.resolve()),'--out',str(dna),'--threads',str(args.threads)]
                with (out/'dna_rescue.log').open('w') as log:
                    rescued=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
                manifest['stages'].append(dict(stage='dna_rescue',command=command,exit_code=rescued.returncode));save()
                if rescued.returncode: raise RuntimeError('DNA rescue failed; see '+str(out/'dna_rescue.log'))
                from research_workflow.candidate_screen import dna_only_catalog
                rows=dna_only_catalog(dna,out/'candidate_function_catalog.tsv',args.genome.resolve())
                manifest.update(status='partial',analysis_scope='DNA_only_unresolved_candidates',
                    protein_annotation_status='unavailable_gene_prediction_failed_or_empty',
                    identity_and_function_status='not_evaluated_without_protein_annotation',
                    candidates=len(rows),core_anchors=0,function_evaluated_candidates=0,
                    protein_candidate_count=0,DNA_candidate_count=len(rows),
                    function_predictions=dict(Counter(r['function_prediction'] for r in rows)),
                    functional_annotations=str(out/'candidate_function_catalog.tsv'),
                    DNA_rescue_status=json.loads((dna/'dna_rescue_manifest.json').read_text())['status'],
                    output_sha256={str(p.relative_to(out)):sha(p) for p in [out/'candidate_function_catalog.tsv',dna/'dna_rescue_manifest.json']})
                record_translation_hypotheses(dna, out, manifest)
                manifest.update(preferred_summary(rows, out/'candidate_function_catalog.tsv'))
                save();print(json.dumps({k:manifest[k] for k in ['status','analysis_scope','candidates','functional_annotations','preferred_function_predictions','preferred_prediction_confidence','preferred_prediction_output']},ensure_ascii=False));return
        from research_workflow.candidate_screen import catalog
        from research_workflow.claim_assertions import VERSION as CLAIM_ASSERTIONS_VERSION, from_results
        from research_workflow.start_hypotheses import audit as audit_member_starts
        from research_workflow.maturation_observation import audit as audit_maturation
        core_rows = read_table(annotation/'catalytic_anchors.tsv')
        evaluated_rows = read_table(functions/'functional_annotations.tsv')
        class_manifest=json.loads((functions/'function_rules/manifest.json').read_text())
        rows=catalog(candidates,functions,out/'candidate_function_catalog.tsv',confirmation,dna_rescue_dir=dna)
        claims = from_results(functions)
        manifest.update(status='complete',candidates=len(rows),core_anchors=len(core_rows),function_evaluated_candidates=len(evaluated_rows),
            protein_candidate_count=sum(r['candidate_kind']=='protein' for r in rows),DNA_candidate_count=sum(r['candidate_kind']=='dna_locus' for r in rows),
            DNA_rescue_status=json.loads((dna/'dna_rescue_manifest.json').read_text())['status'],
            function_predictions=dict(Counter(r['function_prediction'] for r in rows)),
            functional_annotations=str(out/'candidate_function_catalog.tsv'),core_function_annotations=str(functions/'functional_annotations.tsv'),
            claim_assertions_schema_version=CLAIM_ASSERTIONS_VERSION,
            claim_assertions_policy_sha256=sha(ROOT/'research_workflow/claim_assertions.py'),
            fes_refinement=str(fes/'fes_sites.tsv'),
            output_sha256={str(p.relative_to(out)):sha(p) for p in [out/'candidate_function_catalog.tsv',dna/'dna_rescue_manifest.json',candidates/'candidate_screen_manifest.json',confirmation/'candidate_confirmation_manifest.json',functions/'functional_annotations.tsv',functions/'function_attributes.tsv',functions/'manifest.json',claims,domains/'query_domain_manifest.json',fes/'manifest.json']})
        fn_manifest=json.loads((functions/'manifest.json').read_text())
        manifest['fdh_coding_hypotheses']=fn_manifest.get('fdh_coding_hypotheses') or {'status':'not_assessed_mechanistic_mode','candidate_count':0,'supported_conditional_count':0}
        manifest['mechanistic_rule_version']=fn_manifest.get('mechanistic_rule_version','mechanistic_rules_v1')
        manifest.update(function_rules_version=class_manifest['version'],function_class_predictions=class_manifest['predicted_class_counts'])
        for name in ['manifest.json','predictions.tsv','evidence.json','member_hits.tsv']:
            manifest['output_sha256']['functions/function_rules/'+name]=sha(functions/'function_rules'/name)
        if manifest['fdh_coding_hypotheses'].get('path'):
            manifest['fdh_coding_hypotheses']['path']='functions/'+manifest['fdh_coding_hypotheses']['path']
        record_translation_hypotheses(dna, out, manifest)
        member_starts = out/'member_start_hypotheses.json'
        starts = audit_member_starts(dna, claims, member_starts)
        manifest.update(member_start_hypotheses_status=starts['status'],
                        member_start_hypotheses_count=len(starts['records']))
        manifest['output_sha256'][member_starts.name] = sha(member_starts)
        maturation = audit_maturation(annotation, out)
        manifest['maturation_observation_schema_version'] = maturation['schema']
        manifest['output_sha256']['maturation_observation.json'] = sha(out/'maturation_observation.json')
        manifest.update(preferred_summary(rows, out/'candidate_function_catalog.tsv'))
        save(); print(json.dumps({k:manifest[k] for k in ['status','candidates','function_predictions','functional_annotations','preferred_function_predictions','preferred_prediction_confidence','preferred_prediction_output']},ensure_ascii=False))
    except Exception as error:
        manifest.update(status='failed',error=str(error)); save(); raise


if __name__ == '__main__':
    try: main()
    except (OSError, ValueError, RuntimeError) as error: sys.exit(str(error))
