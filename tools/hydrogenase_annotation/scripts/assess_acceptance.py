#!/usr/bin/env python3
"""Assess the complete registered objective; missing independent evidence never passes."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from annotate import ROOT
from run_acceptance_checks import state_hashes
from discovery_evidence import stage_fingerprint
from evaluate_identity_corpus import identity_metrics,family


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle,delimiter="\t"))


def grade_negative_gene(directory,protein_id,expected_sha256):
    """Missing inputs/columns are failed evidence, never successful abstentions."""
    from Bio import SeqIO
    with (directory/'proteins.faa').open() as handle:records=list(SeqIO.parse(handle,'fasta'))
    matching=[r for r in records if r.id==protein_id]
    if len(matching)!=1 or hashlib.sha256(str(matching[0].seq).encode()).hexdigest()!=expected_sha256:
        raise ValueError('Negative benchmark gene missing, duplicated or changed: '+protein_id)
    selected=[]
    for name,key,columns in [('catalytic_anchors.tsv','protein_id',{'protein_id','catalytic_identity_confidence'}),
                             ('reaction_capabilities.tsv','complex_anchor',{'complex_anchor','function_id'})]:
        with (directory/name).open() as handle:
            reader=csv.DictReader(handle,delimiter='\t')
            if not columns<=set(reader.fieldnames or []):raise ValueError('Negative benchmark schema missing: '+name)
            rows=list(reader)
        if any(any(r.get(c) is None for c in columns) for r in rows):raise ValueError('Incomplete benchmark output row')
        selected.append([r for r in rows if r[key]==protein_id])
    anchors,reactions=selected
    return dict(protein_id=protein_id,calibrated_high_identity=any(r['catalytic_identity_confidence']=='calibrated_high' for r in anchors),
                assigned_reactions=sorted({r['function_id'] for r in reactions}),anchor_records=len(anchors))


def hard_negative_panel(cases,grades,contexts,positive_controls,minimum):
    """Recount the preregistered panel; contexts are repeated views, not extra TNs."""
    case_ids={c['case_id'] for c in cases};byunit={c['unit_id']:c for c in cases}
    if len(case_ids)!=len(cases) or len(byunit)!=len(cases) or {g['case_id'] for g in grades}!=case_ids or len(grades)!=len(cases):
        raise ValueError('Negative benchmark cases missing or duplicated')
    support={}
    for name in ['Ehr','NARF']:
        rows=[c for c in cases if c['family']==name]
        support[name]=dict(sequences=len(rows),genera=len({g for r in rows for g in r['genus_taxids']}),components=len({r['homolog_component'] for r in rows}))
    enough=all(support[f][k]>=minimum['new_'+f+'_'+k] for f in support for k in ['sequences','genera','components'])
    if any(r['unit_id'] not in byunit or byunit[r['unit_id']]['family']!='Ehr' for r in contexts):raise ValueError('Unknown native Ehr view')
    context_genera={g for r in contexts for g in byunit[r['unit_id']]['genus_taxids']}
    enough &= len(context_genera)>=minimum['new_Ehr_native_context_genera'] and any(c['family']=='Nuo' for c in cases)
    violations=[r for r in grades+contexts if r['grade']['calibrated_high_identity'] or r['grade']['assigned_reactions']]
    controls_ok=bool(positive_controls) and all(r['identity_passed'] and r['reaction_passed'] for r in positive_controls)
    return dict(panel_passed=bool(enough and not violations and controls_ok),source_support=support,native_context_genera=len(context_genera),
                negative_reference_identity_units=len(cases),native_context_additional_views=len(contexts),violations=violations,positive_controls=positive_controls)


def verify_recorded_inputs(items,historical=None):
    """An old source audit retains its original bytes when a reference bank is versioned."""
    for name,digest in items.items():
        current=hashlib.sha256(Path(name).read_bytes()).hexdigest()
        if current==digest:continue
        previous=(historical or {}).get(name)
        if previous is None or hashlib.sha256(Path(previous).read_bytes()).hexdigest()!=digest:
            raise ValueError('Negative benchmark source/artifact changed: '+str(name))


def verify_hard_negatives(directory):
    """Verify the sealed R45 identity/overclaim challenge against raw outputs."""
    from Bio import SeqIO
    from prepare_fes_calibration import homolog_overlap
    from annotate import extract_genbank
    def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    def load(path):return json.loads(Path(path).read_text())
    def hashes(items,historical=None):verify_recorded_inputs(items,historical)
    base=ROOT/'references/validation/round45_hard_negative_benchmark';manifest_path=base/'benchmark_manifest.json'
    manifest=load(manifest_path);execution=load(directory/'execution.json');plan=load(ROOT/'experiments/round45/plan.json')
    seal=load(ROOT/'experiments/round45/benchmark_completion_seal.json');hashes(seal['artifacts'])
    if seal['preprediction_manifest_sha256']!=sha(manifest_path):raise ValueError('Completed benchmark manifest changed')
    source_state=manifest['production_state'];state=execution['production_state'];current=state_hashes()
    replay_rounds={'round47_repeat_guard_regression_replay':'round47','round49_fes_reference_regression_replay':'round49',
                   'round50_evaluation_semantics_regression_replay':'round50',
                   'round52_state_truth_regression_replay':'round52'}
    replay_round=replay_rounds.get(execution.get('scope'));replay=bool(replay_round)
    runner=ROOT/(f'experiments/{replay_round}/run_hard_negative_replay.py' if replay else 'experiments/round45/run_benchmark.py')
    if replay:
        lock_path=ROOT/'experiments'/replay_round/'pipeline_state_lock.json';locked=load(lock_path)
        if state!=locked['production_state'] or execution['pipeline_lock_sha256']!=sha(lock_path):
            raise ValueError('Negative replay predictor version is not locked')
        allowed_data=set()
        if replay_round in {'round49','round50','round52'}:
            from audit_fes_bank_transfer import verify_transition as verify_fes_transition
            verify_fes_transition(ROOT)
            allowed_data=set(load(ROOT/'experiments/round49/plan.json')['allowed_changed_files'])
        if any(state.get(k)!=v for k,v in source_state.items() if not k.endswith('.py') and k!='config/catalytic_confirmation_thresholds.json' and k not in allowed_data):
            raise ValueError('Negative replay changed source/model data')
        hashes({str(directory/name):value for name,value in execution['artifacts'].items()})
    if (execution['status']!='complete' or execution['production_state']!=state or execution['after_state']!=state or
        source_state!=plan['production_before'] or any(state.get(k)!=current.get(k) for k in set(state)|set(current) if k!='scripts/assess_acceptance.py') or
        execution['manifest_sha256']!=sha(manifest_path) or execution['runner_sha256']!=sha(runner) or
        manifest['plan_sha256']!=sha(ROOT/'experiments/round45/plan.json') or plan['targets_sha256']!=sha(ROOT/'config/acceptance_targets.json')):
        raise ValueError('Negative benchmark execution, model, protocol or target changed')
    hashes(manifest['input_sha256']);hashes(plan['input_sha256'])
    audit=load(ROOT/'results/round45_negative_source_admission/before_source_alignment.json');historical=None
    if replay_round in {'round49','round50','round52'}:
        historical={str(ROOT/name):ROOT/'experiments/round49/baseline_view'/name for name in
            ['references/processed/fes/templates.faa','references/processed/fes/site_windows.faa']}
    hashes(audit['input_sha256'],historical)
    if audit['plan_sha256']!=manifest['plan_sha256'] or audit['policy']!={'identity_percent':70,'shorter_coverage':.7}:raise ValueError('Negative source admission protocol changed')
    units=load(ROOT/'references/validation/round09_identity_corpus/units.json');byid={u['unit_id']:u for u in units}
    sources=load(ROOT/'references/validation/round38_seed_source_audit/resolved_prototype_taxonomy.json')+load(ROOT/'references/validation/round38_seed_source_audit/seed_bindings.json')
    excluded_genera={str(g) for u in sources+[u for u in units if u['split'] in ['development','calibration']] for g in u['genus_taxids']}
    excluded_components={u['homolog_component'] for u in units if u['split'] in ['development','calibration']}
    prechecks={r['unit_id']:r for r in load(ROOT/'results/round45_negative_source_admission/source_prechecks.json')}
    history={hashlib.sha256(r['sequence'].encode()).hexdigest():r for r in load(ROOT/'references/validation/round34_external_function/historical_subjects.json')}
    blocked_roles={'prior_query','prior_exposure','actual_HMM_contiguous_training_spans','deployed_or_reviewed_reference','historical_design_artifact'}
    overlaps=set()
    for path in [ROOT/'results/round45_negative_source_admission/source_alignment.tsv',ROOT/'results/round45_GT1_source_precheck/source_alignment.tsv']:
        for line in path.read_text().splitlines():
            q,s,identity,*coords=line.split('\t')
            if homolog_overlap(float(identity),*map(int,coords[:6]),audit['policy']):overlaps.add(q)
    protected=set()
    for path in [ROOT/'results/round10_final_discovery/annotation/proteins.faa',ROOT/'results/round17_final_confirmation/annotation/proteins.faa']:
        protected.update(hashlib.sha256(str(r.seq).rstrip('*').encode()).hexdigest() for r in SeqIO.parse(path,'fasta'))
    for row in load(ROOT/'references/validation/round35_archaeal_sources/construct_correspondence.json'):
        if row['code'] in ['Mu','Ia','Ps','Na','Fm']:protected.update(row[k] for k in ['native_sequence_sha256','construct_sequence_sha256'])
    cases=manifest['negative_cases'];sequences={r.id:str(r.seq) for r in SeqIO.parse(base/'protein_challenge.faa','fasta')}
    for case in cases:
        if (hashlib.sha256(sequences[case['case_id']].encode()).hexdigest()!=case['sequence_sha256'] or case['sequence_sha256'] in protected or
            not case['genus_taxids'] or set(map(str,case['genus_taxids']))&excluded_genera or case['unit_id'] in overlaps):
            raise ValueError('Negative source sequence/taxon isolation failed: '+case['case_id'])
        if case['unit_id'] in byid:
            u=byid[case['unit_id']];pre=prechecks[case['unit_id']]
            if (u['label']!='negative' or u['split']!=case['original_partition_preserved'] or u['split'] not in ['split_bridge','source_review'] or
                not pre['passes_sequence_taxon_precheck'] or pre['reasons'] or u['homolog_component'] in excluded_components or
                set(history.get(u['sequence_sha256'],{}).get('roles',[]))&blocked_roles or
                any(case[k]!=u[k] for k in ['sequence_sha256','genus_taxids','homolog_component','evidence_tiers']) or
                case['source_accessions']!=[s['accession'] for s in u['sources']]):raise ValueError('Negative catalog source changed or exposed')
        else:
            selection=load(ROOT/'experiments/round45/GT1_source_selection.json');raw=ROOT/'references/raw/round45_hard_negatives'
            record=SeqIO.read(raw/'TGGT1_242580.protein.gp','genbank');dna=SeqIO.read(raw/'GT1_NAR1_native_CDS.fna','fasta')
            if (case['source_accessions']!=[record.id] or case['sequence_sha256']!=selection['protein_sha256'] or str(record.seq)!=sequences[case['case_id']] or
                str(dna.seq.translate(cds=True))!=str(record.seq) or selection['primary_gene'] not in (raw/'PMC12674569.xml').read_text() or
                case['truth_scope']!=selection['truth_scope'] or case['genus_taxids']!=[selection['genus_taxid']]):raise ValueError('Primary Nar1 reference-gene correspondence failed')
    if set(sequences)!={c['case_id'] for c in cases+manifest['positive_controls']}:raise ValueError('Unregistered benchmark protein')
    jobs=manifest['jobs']
    if len(execution['executions'])!=len(jobs):raise ValueError('Negative benchmark jobs missing')
    for job,actual in zip(jobs,execution['executions']):
        run_dir=directory/job['name'];run=load(run_dir/'run_manifest.json')
        expected=[str(ROOT/'.env/bin/python'),str(ROOT/'scripts/annotate.py'),*job['arguments'],'--out',str(run_dir),'--threads','2']
        if (actual['name']!=job['name'] or actual['command']!=expected or actual['exit_code']!=0 or sha(directory/(job['name']+'.log'))!=actual['log_sha256'] or
            run['status']!='complete' or run['annotation_script_sha256']!=state['scripts/annotate.py'] or
            run['normalized_proteins_sha256']!=sha(run_dir/'proteins.faa')):raise ValueError('Negative benchmark run failed or changed')
        hashes({str(run_dir/name):value for name,value in run['output_sha256'].items()})
        hashes({r['path']:r['sha256'] for r in run['input_files']})
        if protected&{hashlib.sha256(str(r.seq).encode()).hexdigest() for r in SeqIO.parse(run_dir/'proteins.faa','fasta')}:raise ValueError('Protected final protein queried')
    grades=[dict(case_id=c['case_id'],family=c['family'],unit_id=c['unit_id'],grade=grade_negative_gene(directory/'protein_identity_challenge',c['case_id'],c['sequence_sha256'])) for c in cases]
    byunit={c['unit_id']:c for c in cases};contexts=[]
    for native in manifest['native_contexts']:
        job=next(j for j in jobs if j.get('unit_id')==native['unit_id']);case=byunit[native['unit_id']]
        if sha(native['context_path'])!=native['context_sha256'] or native['source_sequence_sha256']!=case['sequence_sha256']:raise ValueError('Native context source changed')
        contexts.append(dict(job=job['name'],unit_id=case['unit_id'],family='Ehr',grade=grade_negative_gene(directory/job['name'],job['target_input_gene'],case['sequence_sha256'])))
    controls=[]
    for c in manifest['positive_controls']:
        g=grade_negative_gene(directory/'protein_identity_challenge',c['case_id'],hashlib.sha256(sequences[c['case_id']].encode()).hexdigest())
        controls.append(dict(case_id=c['case_id'],grade=g,identity_passed=g['calibrated_high_identity'],reaction_passed=c['expected_reference_reaction'] in g['assigned_reactions']))
    job=next(j for j in jobs if j['name']=='positive_Hya_native');genes,_=extract_genbank(Path(job['arguments'][1]));gene=next(g for g in genes if g['id']==job['target_input_gene'])
    g=grade_negative_gene(directory/job['name'],gene['id'],hashlib.sha256(gene['sequence'].encode()).hexdigest())
    controls.append(dict(case_id=job['name'],grade=g,identity_passed=g['calibrated_high_identity'],reaction_passed=bool(job['expected_reactions']) and set(job['expected_reactions'])<=set(g['assigned_reactions'])))
    result=hard_negative_panel(cases,grades,contexts,controls,plan['required_panel'])
    legacy=execution['legacy_guard'] if replay else manifest['legacy_guard'];legacy_path=Path(legacy['path']);summary=load(legacy_path/'summary.json')
    if legacy['required_checks']!=manifest['legacy_guard']['required_checks']:raise ValueError('Negative replay omitted legacy scoped guards')
    if (sha(legacy_path/'summary.json')!=legacy['summary_sha256'] or sha(legacy_path/'checks.tsv')!=legacy['checks_sha256'] or
        summary['before_hashes']!=state or summary['after_hashes']!=state):raise ValueError('Legacy scoped guards changed or stale')
    checks={r['test']:r['passed']=='True' for r in table(legacy_path/'checks.tsv')}
    if not all(checks.get(c,False) for c in legacy['required_checks']):result['panel_passed']=False
    for run in summary['runs']:
        run_dir=legacy_path/run['name'];run_manifest=load(run_dir/'run_manifest.json')
        if sha(run_dir/'run_manifest.json')!=run['run_manifest_sha256']:raise ValueError('Legacy guard execution changed')
        hashes({str(run_dir/name):value for name,value in run_manifest['output_sha256'].items()})
    if grades!=load(directory/'negative_results.json') or contexts!=load(directory/'native_context_results.json'):raise ValueError('Reported negative results differ from raw predictions')
    result.update(manifest_sha256=sha(manifest_path),truth_scope=manifest['truth_scope'],legacy_scoped_guards_passed=all(checks.get(c,False) for c in legacy['required_checks']))
    if historical:result['historical_source_input_copies']={name:str(path) for name,path in historical.items()}
    return result


def verify_catalytic_sites(directory):
    """Recount frozen conditional L2 evidence without searching or refitting queries."""
    from Bio import SeqIO
    sys.path.insert(0,str(ROOT/'experiments/round47'))
    from run_natural_validation import grade_sources,summarize
    load=lambda p:json.loads(p.read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    experiment=ROOT/'experiments/round47';base=ROOT/'references/validation/round47_active_center_truth'
    execution=load(directory/'execution.json');locked=load(experiment/'pipeline_state_lock.json')
    protocol=load(experiment/'natural_evaluation_protocol.json');freeze=load(base/'gold_freeze.json')
    current=state_hashes();state=locked['production_state'];layer_transition=None
    if any(current.get(k)!=state.get(k) for k in set(current)|set(state) if k!='scripts/assess_acceptance.py'):
        from audit_fes_bank_transfer import verify_layer_state
        layer_transition=verify_layer_state(state,current)
    if (execution['status']!='complete' or execution['before_state']!=execution['after_state'] or
        execution['before_state']!=state or execution['new_full_pipeline_final_queries']!=0 or
        execution['protocol_sha256']!=sha(experiment/'natural_evaluation_protocol.json') or
        execution['runner_sha256']!=sha(experiment/'run_natural_validation.py') or
        execution['gold_freeze_sha256']!=sha(base/'gold_freeze.json') or protocol['primary_target']!=.98):
        raise ValueError('Conditional L2 predictor, metric or execution changed')
    for name,value in protocol['inputs'].items():
        if sha(Path(name))!=value:raise ValueError('Changed L2 preregistered input: '+name)
    for name,value in freeze['artifacts'].items():
        if sha(base/name)!=value:raise ValueError('Changed L2 gold: '+name)
    for name,value in execution['artifacts'].items():
        if sha(directory/name)!=value:raise ValueError('Changed L2 output: '+name)
    sources=load(base/'sources.json');audit=load(ROOT/'results/round46_L2_source_admission/source_audit.json')
    queries={r.id:str(r.seq) for r in SeqIO.parse(base/'natural_queries.faa','fasta')}
    if set(queries)!={s['protein_id'] for s in sources}:raise ValueError('Missing or extra L2 query')
    for source in sources:
        if (sha(Path(source['source_path']))!=source['source_sha256'] or
            hashlib.sha256(queries[source['protein_id']].encode()).hexdigest()!=source['sequence_sha256'] or
            queries[source['protein_id']]!=source['sequence']):raise ValueError('L2 source sequence changed')
    selected=load(directory/'natural/selected.json');raw=table(directory/'natural/active_center_summary.tsv')
    if sorted((r['protein_id'],r['reference_id'],r['family']) for r in selected)!=sorted((r['protein_id'],r['reference_id'],r['family']) for r in raw):
        raise ValueError('L2 selected outputs disagree')
    grades=grade_sources(sources,load(ROOT/'references/processed/active_centers/models.json'),load(base/'coordinate_correspondence.json'),
        load(ROOT/'experiments/round46/manual_ligand_positions.json')['ligands'],selected,table(directory/'natural/catalytic_sites.tsv'),audit)
    if grades!=load(directory/'site_grades.json'):raise ValueError('L2 grades differ from source coordinates and raw output')
    observed=summarize(grades,.98);reported=load(directory/'summary.json')
    if any(reported.get(k)!=v for k,v in observed.items()):raise ValueError('L2 reported metrics disagree with raw grades')
    boundary={r['protein_id']:r for r in table(directory/'boundary_controls/active_center_summary.tsv')};sites=table(directory/'boundary_controls/catalytic_sites.tsv');controls=[]
    for case in load(ROOT/'results/round46_template_boundary_development/expected.json'):
        selected=boundary.get(case['protein_id'],{});ligands=[r for r in sites if r['protein_id']==case['protein_id'] and r['reference_id']==selected.get('reference_id') and r['residue_role']=='direct_metal_ligand']
        state=ligands[case['ligand_index']]['status'] if len(ligands)>case['ligand_index'] else None
        controls.append(dict(**case,observed_state=state,passed=state==case['expected_state']))
    if controls!=load(directory/'boundary_grades.json'):raise ValueError('L2 boundary controls differ from raw output')
    supported=observed['natural_mapping_target_passed'] and len(controls)==35 and all(c['passed'] for c in controls)
    if reported['site_gate_supported']!=supported:raise ValueError('L2 gate decision disagrees with original metric')
    return dict(**observed,site_gate_supported=supported,boundary_passed=sum(c['passed'] for c in controls),boundary_cases=len(controls),
        primary_source_units=sum(r['L2_conditional_eligible'] for r in audit),new_full_pipeline_final_queries=0,
        source_components=load(experiment/'source_component_report.json')['components'],execution_sha256=sha(directory/'execution.json'),
        layer_transition=layer_transition,
        scope='Conditional L2 coordinates/states only; site intervals are descriptive; no end-to-end, occupancy, activity, architecture or function accuracy claim.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round",type=Path,required=True,help="Completed run_acceptance_checks output")
    parser.add_argument('--out',type=Path,help='Fresh assessment JSON file; preserves the original execution-directory report')
    parser.add_argument("--pilot",type=Path,help="Optional real-genome pilot directory; a small unlabeled pilot cannot satisfy M4")
    parser.add_argument("--identity-development",type=Path,help="Optional frozen identity development diagnostic; cannot pass final blind-test gates")
    parser.add_argument('--identity-final',type=Path,help='Optional locked final candidate-discovery validation; does not validate confirmation/complex/function')
    parser.add_argument('--confirmation-final',type=Path,help='Optional separately locked final catalytic-identity validation; does not validate structure/function')
    parser.add_argument('--localization',type=Path,help='Optional verified auxiliary local-topology run; does not pass independent complex/function accuracy gates')
    parser.add_argument('--hard-negatives',type=Path,help='Optional sealed R45 difficult-homolog identity/reaction-overclaim benchmark; not a population FPR estimate')
    parser.add_argument('--catalytic-sites',type=Path,help='Frozen R47 conditional L2 structure/state validation; not end-to-end accuracy')
    args=parser.parse_args();directory=args.round.resolve()
    if args.out and args.out.exists():parser.error('Use a new assessment output file; existing evidence cannot be overwritten')
    targets=json.loads((ROOT/"config/acceptance_targets.json").read_text())
    execution=json.loads((directory/"execution.json").read_text())
    current_state=state_hashes();evaluated_state=execution['after_hashes']
    changed_since_execution=sorted(k for k in set(current_state)|set(evaluated_state) if current_state.get(k)!=evaluated_state.get(k))
    # Recounting a frozen experiment may repair this verifier, but never its predictor, tests or data.
    immutable=(execution["status"]=="checks_passed" and execution["before_hashes"]==evaluated_state and
        not (set(changed_since_execution)-{'scripts/assess_acceptance.py'}))
    for job in execution["executions"]:
        immutable &= job["exit_code"]==0 and hashlib.sha256((directory/(job["name"]+".log")).read_bytes()).hexdigest()==job["log_sha256"]
    if not immutable:
        raise ValueError("Checks missing, failed, stale, or changed after execution; rerun before assessing")
    references=json.loads((ROOT/"references/processed/references.json").read_text())
    old=table(ROOT/"revisions/v1_20261026/results/validation/checks.tsv")
    current=table(directory/"regression/checks.tsv")
    regression={r["test"] for r in old}<={r["test"] for r in current} and all(r["passed"]=="True" for r in current)
    for run in (directory/"regression").rglob("run_manifest.json"):
        manifest=json.loads(run.read_text())
        regression &= manifest["status"]=="complete" and manifest["annotation_script_sha256"]==execution["after_hashes"]["scripts/annotate.py"]
        regression &= manifest["model_sha256"]==execution["after_hashes"]["references/processed/catalytic_families.hmm"]
    safety=table(directory/"input_safety/checks.tsv")
    safety_ok=len(safety)>=29 and all(r["passed"]=="True" for r in safety)
    units=(directory/"unit.log").read_text()
    mode_ok=safety_ok and "test_bin_map_requires_unique_explicit_ownership" in units and "test_transcript_mode_rejects_eukaryotic_genomic_prediction" in units
    repeat=[r for r in safety if "repeatability" in r["test"]]
    reproducibility=len(repeat)>=11 and all(r["passed"]=="True" for r in repeat)
    from audit_reference_domains import verify_applied_review
    build=json.loads((ROOT/'references/processed/build_manifest.json').read_text())
    applied=verify_applied_review(ROOT,build)
    reference_audit=json.loads((directory/'reference_domains/summary.json').read_text())
    maturation=json.loads((directory/'maturation_domains/summary.json').read_text())
    reference_ok=(reference_audit['status']=='review_complete' and reference_audit['non_maturation_reviewed']==140 and
        reference_audit['structure_systems']==25 and reference_audit['function_records_reviewed']==28 and
        reference_audit['reviewed_references_sha256']==applied['reviewed_references_sha256'] and
        reference_audit['config_sha256']==applied['config_sha256'] and reference_audit['source_sequence_and_reference_kinds_unchanged'] and
        len(references)==214 and all(all(r.get(k) for k in ['protein_accession','sequence_sha256','source_url','evidence_url','evidence_scope']) for r in references) and
        maturation['status']=='passed' and maturation['total']==maturation['passed']==76 and
        sum(r['review_status']=='exclude_non_maturase_homonym' for r in maturation['rows'])==2)
    fes=json.loads((directory/'fes_prediction/summary.json').read_text())
    fes_metrics=fes['metrics']['source_withheld']
    active=json.loads((directory/'active_centers/summary.json').read_text())
    if active['status']!='development_checks_passed' or not active['unchanged'] or not all(r['passed']=='True' for r in table(directory/'active_centers/checks.tsv')):
        raise ValueError('Active-center development checks failed or changed')
    active_metrics=active['metrics']['source_withheld_evaluation']
    split=json.loads((directory/'split_centers/summary.json').read_text())
    if (split['status']!='development_checks_passed' or not split['unchanged'] or split['passed']!=split['checks'] or
        split['policy_sha256']!=hashlib.sha256((ROOT/'config/split_centers.json').read_bytes()).hexdigest()):
        raise ValueError('Split-center development checks failed or changed')
    graph=json.loads((directory/'module_graph/summary.json').read_text())
    if (graph['status']!='development_checks_passed' or not graph['unchanged'] or
        graph['policy_sha256']!=hashlib.sha256((ROOT/'config/module_graph.json').read_bytes()).hexdigest()):
        raise ValueError('Module-graph development checks failed or changed')
    partner=json.loads((directory/'partner_associations/summary.json').read_text())
    if (partner['status']!='development_checks_passed' or not partner['unchanged'] or partner['policy_sha256']!=graph['policy_sha256']):
        raise ValueError('Cross-anchor partner checks failed or changed')
    confirmation_checks=json.loads((directory/'catalytic_confirmation/summary.json').read_text())
    if confirmation_checks['status']!='development_checks_passed' or not confirmation_checks['unchanged']:
        raise ValueError('Catalytic identity source/boundary controls failed or changed')
    output_contract=json.loads((directory/'output_contract/summary.json').read_text())
    contract_checks=table(directory/'output_contract/checks.tsv')
    contract_ok=(output_contract['status']=='passed' and output_contract['runs']>0 and output_contract['empty_candidate_runs']>0 and
        set(output_contract['modes'])=={'genome','bin','assembly','transcript'} and all(r['passed']=='True' for r in contract_checks))
    evidence={
        "M0_REFERENCE":('passed' if reference_ok else 'partial',f"Actual accession/hash/source-pinned role audit of {reference_audit['non_maturation_reviewed']} non-maturation references, 25 systems and 28 scoped function records; evidence tiers={reference_audit['domain_status_counts']}. 74 maturation profiles plus two excluded homonyms rechecked. Partial domains, inferred roles and unresolved references remain explicit; sequences/kinds unchanged. Reference curation is not full-domain/query-function accuracy validation."),
        "M0_REGRESSION":("passed" if regression else "not_met",f"All {len(old)} original cases retained; {len(current)} current cases; actual run script/model hashes checked."),
        "M1_DISCOVERY":("not_met","No frozen independent >=100-positive set or family-wise recall with confidence intervals yet."),
        "M1_GROUP4":("partial","Five source-withheld seed diagnostics plus one external development case pass HMM-only checks. Development case is not a blind test; >=20 independent positives across >=4 branches not yet available."),
        "M1_PRECISION":("not_met","High-confidence catalytic precision has no independent mixed-set denominator or confidence interval."),
        "M1_NEGATIVES":("partial","Five existing decoys remain rejected despite expanded HMM hits. Independent Ehr/NARF and expanded difficult-negative truth set still missing."),
        "M2_CATALYTIC_SITES":("partial",f"{active['template_proteins']} single-chain templates map {active['model_positions']} protein ligand/neighbor positions. One additional physical VhuA/U center maps four ligands across two chains; {split['passed']}/{split['checks']} development checks include two source-reviewed natural operons of the same genus, processing, ambiguity, truncation, ownership and synthetic fusion controls. Paired FCO iron is parsed without loosening the ligand-distance threshold. Required local contexts and whole-template scaffold quality are separate. Actual source-withheld single-chain mapping {active_metrics['correct_ligand_positions']}/{active_metrics['eligible_ligand_positions']}, Wilson95={active_metrics['wilson95']}; {active['passed']}/{active['checks']} checks pass. These development mappings do not establish broad homolog/taxon independent >=98% state mapping, neighbor truth or complete domain-architecture accuracy."),
        "M2_SEC_INPUT":("passed" if safety_ok and regression and "test_group4_n_terminal_sec_candidates_are_not_lost_by_pattern" in units else "not_met","Real FdhF input truncation, N-terminal absent-site, Sec/Cys and ambiguous-letter controls executed; N-terminal Sec variants remain computational candidate boundaries."),
        "M2_FES":("partial",f"Per-query cluster and ligand mapping implemented. Source-withheld development structures: TP={fes_metrics['TP']}, FP={fes_metrics['FP']}, FN={fes_metrics['FN']}, eligible={fes_metrics['eligible_sites']}; precision CI={fes_metrics['precision_wilson95']}. No independent Huc or sealed final calibration; calls remain provisional, occupancy unknown."),
        "M2_MODULE_GRAPH":('partial',f"Coordinate-resource graph for 25 systems/119 nodes; native HoxF plus five boundary controls {graph['passed']}/{graph['checks']}. Seven cross-anchor reference-CDS runs give {partner['passed']}/{partner['checks']} association checks. Fixed-member ambiguity, feasible competing alternatives and diffusible/context sharing are explicit. One VhuA/U catalytic split is represented at L2 with separate precursor and structural windows; complete Vhu complex composition/function is not enabled. Localization remains not predicted. Biological assembly alternatives, broader natural fusion/fission truth, independently validated association decisions and localization/topology remain missing."),
        "M2_COMPLEX":("not_met","Reference retrieval and perturbation controls do not provide independent complex precision/recall or confidence intervals."),
        "P0_INPUT_SCOPE":("passed" if mode_ok else "not_met","Actual genome/bin/assembly/transcript runs plus duplicate/missing-map and unsupported-transcript-DNA guards checked."),
        "P0_CROSS_BIN":("passed" if safety_ok else "not_met","FHL/FdhF and Mvh/Hdr cross-bin/unbinned controls output no combined function and no cross-bin members; catalytic candidates retained."),
        "M3_CALIBRATION":("not_met","Global 70/80 identity transfer thresholds remain uncalibrated. Train/calibration/final-test split and mechanism-specific decision calibration missing."),
        "M3_FUNCTION":("not_met","No independent experimental-function precision, coverage or 95% CI; close-reference success is not sufficient."),
        "M3_SPLIT":("not_met","External development case is excluded from training, but no frozen homolog-cluster/taxonomic final-test split and full leakage audit yet."),
        "OUTPUT_LAYERS":('passed' if contract_ok else 'not_met',f"Actual output-contract audit: {output_contract['runs']} runs, {output_contract['empty_candidate_runs']} empty-candidate runs; all nine tables, stable empty schemas, six evidence layers, separate statuses, software/parameter and artifact hashes, FeS coordinates and source provenance checked. Incomplete scientific models remain explicit in separate gates."),
        "REPRODUCIBILITY":("passed" if reproducibility else "not_met",f"{len(repeat)} repeated fixed-input pairs across supported modes compared annotation/evidence tables byte-for-byte; frozen execution hashes verified."),
        "M4_MAG_PILOT":("not_met","Representative ~1000 real MAG sequences/coordinates, stratified pilot, resource and old/new-result review not run."),
        "M4_FALSE_POSITIVE_BOUND":("not_met","No sufficient independent trustworthy negative set. A <=1e-6 one-sided 95% bound has not been demonstrated."),
        "M4_RELEASE":("not_met","Production release prohibited until every registered required gate is proven; no full-dataset final reannotation claimed.")
    }
    if args.localization:
        from localization import verify_completed
        auxiliary=verify_completed(args.localization)
        status,description=evidence['M2_MODULE_GRAPH']
        description=description.replace('Localization remains not predicted.',
            f"Core graph labels retain their original scope. Separate verified local TMbed evidence covers {auxiliary['selected_protein_owners']} protein owners/{auxiliary['unique_sequences']} unique sequences in {auxiliary['sources']} source runs; residue labels are projected onto actual module resources, without compartment or function promotion.")
        evidence['M2_MODULE_GRAPH']=(status,description)
    if args.pilot:
        pilot=json.loads((args.pilot/"execution.json").read_text())
        evidence["M4_MAG_PILOT"]=("partial" if pilot["status"]=="complete" else "not_met",
            f"Actual supplied-genome pilot: {pilot['input_genomes']} genomes, {len(pilot['executions'])} executions, status {pilot['status']}. Representative ~1000-MAG sampling, taxonomic/quality strata and biological adjudication remain pending.")
    challenge_path=directory/'fes_challenge/summary.json'
    if challenge_path.exists():
        challenge=json.loads(challenge_path.read_text())
        if challenge['before_hashes']!=execution['after_hashes'] or challenge['after_hashes']!=execution['after_hashes']:
            raise ValueError('New structural challenge is stale')
        evidence['M2_FES']=('partial',f"Coordinate-complete sources and exact-sequence corroboration recover FrhB Cys104; template self-test={fes['metrics']['template_self_test']}. Previous source-withheld development={fes_metrics}. New source-selected development challenge={challenge['metrics']}. Predictions outside modeled structural residues or ligand atoms are explicitly ungradable, not FP or TP. No new sealed final test or independent Huc query; query occupancy remains unknown.")
    if args.identity_development:
        diagnostic=json.loads((args.identity_development/'summary.json').read_text())
        if (diagnostic['status']!='development_diagnostic_complete' or diagnostic['split']!='development' or
            diagnostic['model_state_hashes']!=execution['after_hashes'] or diagnostic['after_state_hashes']!=execution['after_hashes']):
            raise ValueError('Identity development diagnostic is missing, changed or stale')
        corpus=Path(diagnostic['corpus'])
        if hashlib.sha256((corpus/'freeze.json').read_bytes()).hexdigest()!=diagnostic['corpus_freeze_sha256']:
            raise ValueError('Identity corpus freeze changed')
        for filename,expected in diagnostic['artifacts'].items():
            path=args.identity_development/filename if filename!='run_manifest.json' else args.identity_development/'annotation'/filename
            if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
                raise ValueError('Identity diagnostic output changed: '+filename)
        total=diagnostic['metrics']['all'];remote=diagnostic['metrics']['remote_Group4']['hmm_discovery']
        hmm=total['hmm_discovery'];supported=total['supported_catalytic']
        evidence['M1_DISCOVERY']=('partial',f"Frozen development identity corpus: HMM-only TP={hmm['TP']}, FN={hmm['FN']}, recall={hmm['recall']}, CI={hmm['recall_wilson95']}. Author-curated/reviewed identity tiers kept separate. Final validation remains reserved; no final recall pass claimed.")
        evidence['M1_GROUP4']=('partial',f"Frozen remote Group4 development cases: TP={remote['TP']}, FN={remote['FN']}, recall={remote['recall']}. Current genus/homolog exposure excluded. Final >=20-case, >=4-branch, >=3-genus acceptance not yet evaluated.")
        evidence['M1_PRECISION']=('partial',f"Current uncalibrated supported-catalytic development status: TP={supported['TP']}, FP={supported['FP']}, precision={supported['precision']}, CI={supported['precision_wilson95']}. Not calibrated high-confidence or final-test precision.")
        context_path=directory/'paralog_contexts/summary.json'
        context=json.loads(context_path.read_text()) if context_path.exists() else None
        if context and (context['status']!='development_context_checks_passed' or context['before_hashes']!=execution['after_hashes'] or context['after_hashes']!=execution['after_hashes']):
            raise ValueError('Native paralog-context checks are failed or stale')
        evidence['M1_NEGATIVES']=('partial',f"Current adjudicated development negatives={diagnostic['labels'].get('negative',0)}, indeterminate={diagnostic['labels'].get('indeterminate',0)}; supported identity FP={supported['FP']}; negative function outputs={len(diagnostic['negative_function_outputs'])}. Native context checks={context['passed'] if context else 0}: exact 13-member MBS source mapping and three complete-genome plus experimental Hmd inputs; MBS is exposed reference replay, and Hmd truth is reaction-scoped. Independent full-context difficult-negative coverage is still incomplete.")
        evidence['M3_SPLIT']=('partial',f"Identity corpus has frozen source hashes, exact-sequence deduplication, genus ownership and homolog-component purging; {diagnostic['units']} development sequence units across {diagnostic['homolog_components']} components/{diagnostic['genera']} genera actually evaluated. Calibration/final identities remain reserved; complex/function-specific split and truth not complete.")
    rule_path=ROOT/'config/catalytic_confirmation_thresholds.json'
    if rule_path.exists() and not args.confirmation_final:
        from catalytic_confirmation import stage_fingerprint as confirmation_fingerprint,sha
        rule=json.loads(rule_path.read_text());cal_path=Path(rule['calibration_summary']);cal=json.loads(cal_path.read_text())
        current_stage=confirmation_fingerprint(ROOT)
        stage_valid=rule['feature_stage_fingerprint']==current_stage
        if not stage_valid and rule.get('input_translation_transition'):
            from audit_input_translation_transfer import verify_registered_transition
            stage_valid=verify_registered_transition(ROOT,ROOT/rule['input_translation_transition'],current_stage,
                old_stage_fingerprint=rule['feature_stage_fingerprint'],rules=rule)
        if not stage_valid and rule.get('reporting_transition'):
            from audit_discovery_reporting_transfer import verify_reporting_transition
            stage_valid=verify_reporting_transition(ROOT,ROOT/rule['reporting_transition'],current_stage,rule)
        if not stage_valid and rule.get('topology_transition'):
            from audit_topology_transition import verify_topology_transition
            stage_valid=verify_topology_transition(ROOT,ROOT/rule['topology_transition'],current_stage,rule)
        if not stage_valid or sha(cal_path)!=rule['calibration_summary_sha256']:
            raise ValueError('Current confirmation calibration is stale')
        evidence['M1_PRECISION']=('partial',f"Current feature/model version has source-audited calibration only: {cal['metrics']['source_audited_independent']}. Previous final-confirmation evidence belongs to the earlier immutable version; no new independent final result is claimed after the primary Hmd truth/discovery correction.")
    if args.identity_final:
        final_path=args.identity_final.resolve();final=json.loads((final_path/'summary.json').read_text())
        if final['status']!='final_discovery_validation_complete' or final['split']!='final_validation' or final['model_state_hashes']!=final['after_state_hashes']:
            raise ValueError('Final discovery validation is missing or changed during execution')
        if final['discovery_stage_fingerprint']!=stage_fingerprint(ROOT):
            raise ValueError('Current discovery stage differs from locked final validation')
        lock_path=Path(final['discovery_lock']);lock=json.loads(lock_path.read_text())
        calibration_path=Path(lock['calibration_summary']);calibration=json.loads(calibration_path.read_text())
        corpus=Path(final['corpus']);audit_path=Path(final['training_audit']);audit=json.loads((audit_path/'summary.json').read_text())
        if (hashlib.sha256(lock_path.read_bytes()).hexdigest()!=final['discovery_lock_sha256'] or lock['stage']!='candidate_discovery_only' or
            lock['model_state_hashes']!=final['model_state_hashes'] or lock['discovery_stage_fingerprint']!=final['discovery_stage_fingerprint'] or
            hashlib.sha256(calibration_path.read_bytes()).hexdigest()!=lock['calibration_summary_sha256'] or calibration['split']!='calibration' or
            calibration['model_state_hashes']!=final['model_state_hashes'] or
            hashlib.sha256((corpus/'freeze.json').read_bytes()).hexdigest()!=final['corpus_freeze_sha256'] or
            hashlib.sha256((audit_path/'summary.json').read_bytes()).hexdigest()!=final['training_audit_summary_sha256'] or
            audit['corpus_freeze_sha256']!=final['corpus_freeze_sha256']):
            raise ValueError('Final discovery lock, calibration or source audit is inconsistent')
        for name,expected in final['artifacts'].items():
            path=final_path/name if name!='run_manifest.json' else final_path/'annotation'/name
            if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Final discovery artifact changed: '+name)
        for name in ['per_unit_audit.tsv']:
            if hashlib.sha256((audit_path/name).read_bytes()).hexdigest()!=audit['artifacts'][name]:raise ValueError('Core independence audit changed')
        frozen=json.loads((corpus/'freeze.json').read_text())
        if hashlib.sha256((corpus/'units.json').read_bytes()).hexdigest()!=frozen['artifacts']['units.json']:raise ValueError('Final truth units changed')
        all_units=json.loads((corpus/'units.json').read_text());independent=set(audit['independent_unit_ids'])
        truth=[r for r in all_units if r['split']=='final_validation' and r['unit_id'] in independent]
        raw_hmm=table(final_path/'annotation/family_candidates.tsv');families={}
        for hit in raw_hmm:families.setdefault(hit['protein_id'],set()).add(hit['family'])
        correct={r['unit_id'] for r in truth if r['unit_id'] in families and (r['label']=='negative' or family(r) in families[r['unit_id']])}
        observed=identity_metrics(truth,correct)
        if observed!=final['metrics']['independent_after_core_audit']['hmm_correct_family_discovery']:raise ValueError('Final discovery counts do not match actual HMM candidates')
        family_metrics={name:identity_metrics([r for r in truth if r['label']=='positive' and family(r)==name],correct) for name in ['NiFe','FeFe']}
        discovery_ok=(observed['TP']+observed['FN']>=100 and all(m['TP']+m['FN']>=20 and m['recall']>=.95 for m in family_metrics.values()))
        evidence['M1_DISCOVERY']=('passed' if discovery_ok else 'not_met',f"Locked final provided-protein HMM discovery, source-audited independent positives={observed['TP']+observed['FN']}: TP={observed['TP']}, FN={observed['FN']}, recall={observed['recall']}, Wilson95={observed['recall_wilson95']}; per-family={family_metrics}. Frozen author-curated/reviewed catalytic identity, not per-protein biochemical activity or complete-genome gene-prediction recall.")
        remote_ids={r['unit_id'] for r in table(audit_path/'per_unit_audit.tsv') if r['remote_Group4_after_core_audit']=='True'}
        remote=[r for r in truth if r['unit_id'] in remote_ids]
        branch_hits={r['protein_id'] for r in raw_hmm if r['profile']=='NiFe_Group4'};rm=identity_metrics(remote,branch_hits)
        branches={c for r in remote for c in r['categories'] if c.startswith('[NiFe] Group 4')};genera={g for r in remote for g in r['genus_taxids']}
        group4_ok=len(remote)>=20 and len(branches)>=4 and len(genera)>=3 and rm['recall']>=.90
        evidence['M1_GROUP4']=('passed' if group4_ok else 'not_met',f"Locked final Group4 profile-only discovery: {len(remote)} source-audited distant positives, {len(branches)} branches/{len(genera)} genera; TP={rm['TP']}, FN={rm['FN']}, recall={rm['recall']}, Wilson95={rm['recall_wilson95']}. Current full/core training comparisons retained; >=70% both-sequence comparable coverage and <=60% max identity required.")
        evidence['M3_CALIBRATION']=('partial','Candidate discovery has a locked final evaluation. Current catalytic identity has source-audited calibration; the earlier identity final belongs to the previous feature/model version. Structure and function thresholds remain uncalibrated; global identity-transfer rules cannot be treated as biochemical boundaries.')
        evidence['M3_SPLIT']=('partial',f"Identity train/development/calibration/final sources frozen; four development-only discovery seeds and actual training spans audited; author FeFe core components add source-only purging before final prediction. Final independent positives={observed['TP']+observed['FN']}. Complex/function-specific truth and split remain incomplete.")
    if args.confirmation_final:
        from catalytic_confirmation import stage_fingerprint as confirmation_fingerprint,qualifies,source_family,sha
        path=args.confirmation_final.resolve();cf=json.loads((path/'summary.json').read_text())
        transfer=None
        if cf['feature_stage_fingerprint']!=confirmation_fingerprint(ROOT):
            from audit_confirmation_transfer import verify_transition
            transfer_path=directory/'confirmation_transfer/summary.json';transfer=verify_transition(ROOT)
            if (not transfer_path.exists() or json.loads(transfer_path.read_text())!=transfer or
                Path(transfer['parent_final'])!=path/'summary.json'):
                raise ValueError('Confirmation transition was not actually checked for this parent/current scope')
        if (cf['status']!='final_confirmation_validation_complete' or
            cf['prediction_state_hashes']!=cf['after_prediction_state_hashes'] or
            (cf['after_prediction_state_hashes']!=execution['after_hashes'] and transfer is None)):
            raise ValueError('Final identity confirmation is missing, stale or changed')
        lock_path=Path(cf['lock']);lock=json.loads(lock_path.read_text());cal_path=Path(lock['calibration_summary']);cal=json.loads(cal_path.read_text())
        audit_path=Path(cf['source_audit']);audit=json.loads(audit_path.read_text());rule_path=Path(transfer['parent_thresholds']) if transfer else ROOT/'config/catalytic_confirmation_thresholds.json';rule=json.loads(rule_path.read_text())
        if (sha(lock_path)!=cf['lock_sha256'] or lock['stage']!='catalytic_identity_confirmation_only' or lock['thresholds_sha256']!=sha(rule_path) or
            lock['calibration_summary_sha256']!=sha(cal_path) or lock['source_audit_sha256']!=sha(audit_path) or
            cf['source_audit_sha256']!=sha(audit_path) or cal['feature_stage_fingerprint']!=cf['feature_stage_fingerprint'] or
            lock['acceptance_targets_sha256']!=sha(ROOT/'config/acceptance_targets.json') or cal['families']!=rule['families']):
            raise ValueError('Confirmation calibration/source/parameter lock changed')
        for filename,digest in cf['artifacts'].items():
            if sha(path/filename)!=digest:raise ValueError('Final confirmation artifact changed: '+filename)
        for filename,digest in audit['artifacts'].items():
            if sha(audit_path.parent/filename)!=digest:raise ValueError('Confirmation source audit artifact changed: '+filename)
        corpus=ROOT/json.loads((ROOT/'config/catalytic_confirmation.json').read_text())['corpus'];freeze=json.loads((corpus/'freeze.json').read_text())
        if sha(corpus/'freeze.json')!=cf['corpus_freeze_sha256'] or sha(corpus/'units.json')!=freeze['artifacts']['units.json']:raise ValueError('Confirmation truth changed')
        truth=[r for r in json.loads((corpus/'units.json').read_text()) if r['split']=='final_validation' and r['unit_id'] in set(audit['independent_unit_ids'])]
        features=json.loads((path/'annotation/confirmation/features.json').read_text());chosen={r['protein_id'] for r in features if rule['families'].get(r['family']) and qualifies(r,rule['families'][r['family']])}
        feature_index={r['protein_id']:r for r in features}
        observed=identity_metrics(truth,chosen);reported=cf['metrics']['source_audited_independent']
        wrong=[r['unit_id'] for r in truth if r['label']=='positive' and r['unit_id'] in chosen and source_family(r)!=feature_index[r['unit_id']]['family']]
        if wrong!=reported['wrong_family_high_outputs']:raise ValueError('Final confirmation family errors do not match truth')
        if observed!={k:v for k,v in reported.items() if k!='wrong_family_high_outputs'}:raise ValueError('Final confirmation metrics do not match predictions')
        raw_high={r['protein_id'] for r in table(path/'annotation/catalytic_confirmation.tsv') if r['identity_confidence']=='calibrated_high'}
        if raw_high!=chosen:raise ValueError('Final pipeline high-confidence labels differ from locked identity rule')
        selection=json.loads((ROOT/'config/catalytic_confirmation.json').read_text())['calibration_selection'];family_checks={}
        for center_family in ['NiFe','FeFe']:
            family_truth=[r for r in truth if r['label']=='negative' or source_family(r)==center_family]
            family_selected={q for q in chosen if feature_index[q]['family']==center_family};m=identity_metrics(family_truth,family_selected)
            family_reported=cf['metrics']['family_'+center_family]
            if m!={k:v for k,v in family_reported.items() if k!='wrong_family_high_outputs'}:raise ValueError('Final family counts changed')
            family_checks[center_family]=bool(rule['families'].get(center_family) and m['TP']>=selection['min_true_positive_outputs_per_family'] and
                m['precision']>=selection['min_precision'] and m['recall']>=selection['min_recall'] and not family_reported['wrong_family_high_outputs'])
        identity_ok=observed['TP']+observed['FN']>=100 and all(family_checks.values()) and observed['precision']>=.99 and not wrong
        evidence['M1_PRECISION']=('passed' if identity_ok else 'partial',f"Separately locked final catalytic identity: source-audited positives={observed['TP']+observed['FN']}, negatives={observed['TN']+observed['FP']}; TP={observed['TP']}, FP={observed['FP']}, FN={observed['FN']}, TN={observed['TN']}; precision={observed['precision']}, Wilson95={observed['precision_wilson95']}, recall={observed['recall']}. NiFe/FeFe family preregistered coverage/output checks={family_checks}. Curated/author identity labels, not per-protein biochemical activity or a population precision guarantee; Hmd/split/noncanonical centers outside new high calibration.")
        if transfer:
            status,reason=evidence['M1_PRECISION']
            if transfer.get('validation_scope') in ['repeat_guard_conservative_subset','FeS_bank_metadata_only_after_repeat_guard']:
                evidence['M1_PRECISION']=('passed' if transfer['precision_gate_supported'] else 'partial',
                    f"Frozen R17 identity evidence inherited through the locked transition ({transfer['status']}): current TP/FP/FN bounds={transfer['current_metric_bounds']}; family coverage checks={transfer['family_checks']}. Potentially affected original high IDs={len(transfer['possibly_changed_high_ids'])}. Fixed prototypes/core scores and thresholds; no new final prediction/fitting runs. Earlier alignment-trace bounds are preserved only within the verified decision-function scope. Finite curated-identity evidence, not new blind activity or population accuracy.")
            else:evidence['M1_PRECISION']=(status,reason+f" R17 evidence inherited only for the unchanged NiFe/FeFe decision function through {transfer['passed']} source/data/code equivalence checks; not a new final test, and no inheritance for the corrected Hmd candidate branch.")
        evidence['M3_CALIBRATION']=('partial',f"Candidate discovery and catalytic identity now have distinct source-separated calibration and preprediction locks. Identity calibration={cal['metrics']['source_audited_independent']}; final identity evaluated without refitting. Structure and function remain uncalibrated; global70/80 identity transfer cannot be treated as mechanism boundaries.")
        if transfer:
            current_rule=json.loads((ROOT/'config/catalytic_confirmation_thresholds.json').read_text())
            current_cal=json.loads(Path(current_rule['calibration_summary']).read_text())
            evidence['M3_CALIBRATION']=('partial',f"Current identity calibration={current_cal['metrics']['source_audited_independent']}; earlier NiFe/FeFe evidence inherited within the verified transition scope ({transfer['status']}), without new final prediction or refitting from final outcomes. Discovery is separately locked; structure and function remain uncalibrated.")
        evidence['M3_SPLIT']=('partial',f"Discovery and confirmation sources have frozen whole-protein/core/genus audits. Confirmation prototypes use development plus exact exposed references; original unresolved function kind remains unchanged for12 source-positive identity correspondences. Final confirmation independent units={len(truth)}. Complex/function truth and their independent splits remain incomplete.")
    natural_path=directory/'natural_complexes/summary.json'
    if natural_path.exists():
        natural=json.loads(natural_path.read_text())
        if (natural['status']!='development_natural_context_evaluated' or
            natural['before_hashes']!=execution['after_hashes'] or natural['after_hashes']!=execution['after_hashes']):
            raise ValueError('Natural-context complex/function evaluation is failed or stale')
        core=natural['metrics']['core'];reaction=natural['metrics']['reaction']
        evidence['M2_COMPLEX']=('not_met',f"Frozen primary-source natural-context development: {natural['cases']} cases, four positive complexes and two non-H2-center contrasts, {sum(r['proteins'] for r in natural['runs'])} proteins from five native genomic/locus sources. Exact assayed core-member sets: {core}; high exact cores={natural['high_exact_core_TP']}. Other {natural['ungraded_anchors']} genomic anchors remain ungradable. No homolog/taxon-independent complete-complex test; independent precision/recall targets are not achieved.")
        evidence['M3_FUNCTION']=('not_met',f"Directly assayed natural-context reaction contracts: {reaction}; non-H2-center contrast hydrogenase-function outputs={natural['contrast_hydrogenase_function_assignments']}; extra unassayed contract outputs={natural['ungraded_reaction_assignments']}. Precision is undefined when no evaluable positive assignments are made. Development cases only; no independent >=95% function precision or representative coverage.")
        if natural.get('extension_reference_replay_case_ids'):
            for gate in ['M2_COMPLEX','M3_FUNCTION']:
                status,description=evidence[gate]
                evidence[gate]=(status,description+f" New post-identity structure templates include these exact development cases: {natural['extension_reference_replay_case_ids']}; exposure-stratified results={natural['metrics_by_exposure']}. Their recovery is reference replay, not independent performance.")
    extension_path=directory/'complex_extensions/summary.json'
    if extension_path.exists():
        extension=json.loads(extension_path.read_text())
        if (extension['status']!='source_review_checks_passed' or extension['passed']!=extension['checks'] or
            extension['model_manifest_sha256']!=hashlib.sha256((ROOT/'references/processed/complex_extensions/manifest.json').read_bytes()).hexdigest() or
            extension['reaction_registry_sha256']!=hashlib.sha256((ROOT/'config/reaction_contracts.json').read_bytes()).hexdigest()):
            raise ValueError('Post-identity reference curation is failed or stale')
        status,description=evidence['M0_REFERENCE']
        evidence['M0_REFERENCE']=(status,description+f" Additional L4-only bank: {extension['additional_structures']} experimentally characterized complex systems/{extension['additional_proteins']} exact native protein correspondences and {extension['typed_contracts']} typed reference capability records programmatically rechecked against their scoped sources. These templates do not enter L1 discovery/confirmation or establish query reaction accuracy.")
        status,description=evidence['M2_MODULE_GRAPH']
        evidence['M2_MODULE_GRAPH']=(status,description+f" Post-identity full-subunit templates add{extension['additional_proteins']}nodes across{extension['additional_structures']}systems; no new domain classifier, FeS topology or independently calibrated function is implied.")
    native_path=directory/'ralbus_context/summary.json'
    if native_path.exists():
        native=json.loads(native_path.read_text())
        if (native['status']!='development_native_assay_challenge_evaluated' or native['before_hashes']!=execution['after_hashes'] or
            native['after_hashes']!=execution['after_hashes']):raise ValueError('New native assay challenge failed or stale')
        for gate,metric in [('M2_COMPLEX','core_metrics'),('M3_FUNCTION','reaction_metrics')]:
            status,description=evidence[gate]
            evidence[gate]=(status,description+f" Additional R.albus native chromosome/plasmid development challenge={native[metric]}; one purified HydABC biological system. HydA2 activity association is not purified-monomer truth; putative HydS is not an absolute negative. No source/genus/homolog independence claimed.")
    function_path=directory/'function_mechanism/summary.json'
    if function_path.exists():
        function=json.loads(function_path.read_text())
        if (function['status']!='development_function_challenges_evaluated' or function['before_hashes']!=execution['after_hashes'] or
            function['after_hashes']!=execution['after_hashes']):raise ValueError('Function mechanism development experiment failed or stale')
        status,description=evidence['M2_COMPLEX']
        evidence['M2_COMPLEX']=(status,description+f" New reaction-scoped native development cases: {function['native_metrics']['core']}; three genomes/{function['native_proteins']} proteins; template replays={function.get('reference_replay_case_ids',[])} include one natural two-polypeptide fusion and one six-subunit non-bifurcating NADP system. Uncharacterized proteins remain ungradable.")
        status,description=evidence['M3_FUNCTION']
        evidence['M3_FUNCTION']=(status,description+f" New native core/reaction-pair metrics={function['native_metrics']['reaction']}. Standalone binding-pattern experiment={function['binding_experiment']['metrics']}, status={function['binding_experiment']['status']}; whole-module nearest-source correct={function['binding_experiment']['nearest_reference_correct']}/10. Failed specificity model remains outside production; no reference-replay/development result passes independent function gates.")
        if any(r.get('noncanonical_literature_evidence') for r in function['native_results']):
            status,description=evidence['M2_CATALYTIC_SITES']
            evidence['M2_CATALYTIC_SITES']=(status,description+' Auxiliary Aw Ser301/3C literature pattern is mapped separately. Three same-source structural models fail the unchanged canonical first-shell geometry screen. Query Ser308 in Syntrophus is a literature-compatible candidate, not experimentally established query ligation; original canonical sites, identity and functional decisions remain unchanged. This is not independent state-mapping validation.')
        if function.get('fusion_boundary_results'):
            status,description=evidence['M2_MODULE_GRAPH']
            controls=function['fusion_boundary_results']
            evidence['M2_MODULE_GRAPH']=(status,description+f" Existing separated Sw role references now test the native two-polypeptide Syntrophus C-before-B fusion as provisional composition, with no new reference protein or identity/function promotion. Full native-context perturbation controls={sum(c['passed'] for c in controls)}/{len(controls)}. Exposed development and synthetic boundaries do not establish independent complex precision/recall.")
    distributed_path=directory/'distributed_function/summary.json'
    if distributed_path.exists():
        distributed=json.loads(distributed_path.read_text())
        distributed_execution=json.loads((distributed_path.parent/'execution.json').read_text())
        if (distributed['status']!='development_distributed_function_evaluated' or not distributed['cross_bin_boundary_passed'] or
            distributed_execution['before_hashes']!=evaluated_state or distributed_execution['after_hashes']!=evaluated_state or
            distributed_execution['summary_sha256']!=hashlib.sha256(distributed_path.read_bytes()).hexdigest()):
            raise ValueError('Distributed native function diagnostic failed or stale')
        native=next(r for r in distributed['rows'] if r['mode']=='genome')
        for gate in ['M2_COMPLEX','M3_FUNCTION','M3_SPLIT']:
            status,description=evidence[gate]
            evidence[gate]=(status,description+f" Source-audited KIST612 diagnostic: observed purified members covered={native['observed_member_coverage']}/4; reaction counts={native['reaction']}; cross-bin boundary passed. Homolog overlap excludes independent final eligibility. HydD membership and untested reactions remain ungraded; this is not an independent precision estimate.")
    sampling_path=directory/'pilot_sampling/freeze.json'
    if sampling_path.exists():
        sampling=json.loads(sampling_path.read_text())
        frozen_pilot=json.loads((ROOT/'references/validation/round30_original_MAG_pilot/freeze.json').read_text())
        selected=sampling_path.parent/'selected_samples.json'
        if (sampling['status']!='metadata_sample_frozen_sequences_unresolved_not_run' or sampling['selected']!=1000 or
            sampling['annotation_runs']!=0 or sampling['source_sha256']!=hashlib.sha256(Path(sampling['source_file']).read_bytes()).hexdigest() or
            sampling['policy_sha256']!=hashlib.sha256((ROOT/'config/mag_pilot_sampling.json').read_bytes()).hexdigest() or
            sampling['artifacts']['selected_samples.json']!=hashlib.sha256(selected.read_bytes()).hexdigest() or
            sampling['artifacts']!=frozen_pilot['artifacts']):
            raise ValueError('Original MAG sampling is changed or not reproducible')
        status,description=evidence['M4_MAG_PILOT']
        evidence['M4_MAG_PILOT']=(status,description+f" Original {sampling['population']}-genome metadata now fixes 1000 bins covering {len(sampling['coverage']['phylum']['selected'])} phylum labels and {len(sampling['coverage']['joint_stratum']['selected'])} joint strata; repeated selection identical. Sequence inputs unresolved and annotation runs=0; manifest preparation does not satisfy the pilot gate.")
    if args.hard_negatives:
        hard_negatives=verify_hard_negatives(args.hard_negatives.resolve())
        evidence['M1_NEGATIVES']=('passed' if hard_negatives['panel_passed'] else 'partial',
            f"Sealed current-model reference-identity/overclaim panel: {hard_negatives['negative_reference_identity_units']} unused negative-identity references, source support={hard_negatives['source_support']}; {hard_negatives['native_context_additional_views']} additional Ehr native-context views across {hard_negatives['native_context_genera']} genera, not extra independent TNs. Wrong calibrated-high identity or assigned-reaction cases={len(hard_negatives['violations'])}. Three exposed positive identity/reaction controls evaluated. Original Nuo/MBS and reaction-scoped Hmd guards retained. This passes only the finite difficult-homolog overclaim check; CIA identity is not universal H2 inactivity, indeterminate Hmd cases are not global activity-negative TNs, and no population1e-6FPR or independent reaction-precision conclusion follows.")
    if args.catalytic_sites:
        catalytic_sites=verify_catalytic_sites(args.catalytic_sites.resolve())
        evidence['M2_CATALYTIC_SITES']=('passed' if catalytic_sites['site_gate_supported'] else 'partial',
            f"Frozen conditional L2 structural validation: {catalytic_sites['primary_source_units']} source-admitted proteins, {len(catalytic_sites['source_components'])}70/70 source components; overall={catalytic_sites['all']}; families={catalytic_sites['by_family']}; ligand/signature/neighbor={catalytic_sites['by_category']}; four-state boundary controls={catalytic_sites['boundary_passed']}/{catalytic_sites['boundary_cases']}. Exact coordinate/residue/state jointly graded against manual ligand positions and independent CE structural correspondence. Ungradable variable context remains excluded from accuracy but included in truth coverage/conservative lower bounds. Most sources overlap upstream identity models; this passes the finite conditional L2 mapping target only, not M3_SPLIT, whole-domain architecture, full H-cluster occupancy, enzyme activity or end-to-end generalization. Source/site intervals are descriptive, not a population98% confidence guarantee.")
    rows=[]
    for requirement in targets["requirements"]:
        status,reason=evidence.get(requirement["id"],("not_met","Missing assessment implementation"))
        rows.append(dict(id=requirement["id"],required=requirement["required"],status=status,target=requirement["target"],evidence=reason))
    complete=all(r["status"]=="passed" for r in rows if r["required"])
    report=dict(status="achieved" if complete else "not_achieved",targets_version=targets["version"],
        target_sha256=hashlib.sha256((ROOT/"config/acceptance_targets.json").read_bytes()).hexdigest(),
        frozen_execution=str(directory/"execution.json"),requirements=rows,
        passed=sum(r["status"]=="passed" for r in rows),required=len(rows),
        next_blocking_work=[r["id"] for r in rows if r["status"]!="passed"])
    report['assessment_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if natural_path.exists():report['natural_complex_development']=str(natural_path)
    if native_path.exists():report['new_native_assay_challenge']=str(native_path)
    report['non_prediction_changes_since_execution']=changed_since_execution
    if args.localization:report['auxiliary_localization']=str(args.localization.resolve()/'localization_manifest.json')
    if args.hard_negatives:report['hard_negative_benchmark']=dict(path=str(args.hard_negatives.resolve()),verification=hard_negatives)
    if args.catalytic_sites:report['conditional_L2_evidence']=dict(path=str(args.catalytic_sites.resolve()),verification=catalytic_sites)
    if args.identity_development:
        report['identity_development']=str(args.identity_development.resolve()/'summary.json')
    if args.identity_final:
        report['identity_final']=str(args.identity_final.resolve()/'summary.json')
    if args.confirmation_final:
        report['confirmation_final']=str(args.confirmation_final.resolve()/'summary.json')
        if transfer:report['confirmation_transition']=str(directory/'confirmation_transfer/summary.json')
    output=args.out.resolve() if args.out else directory/'acceptance.json'
    with output.open('x' if args.out else 'w') as handle:handle.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    if not complete:
        sys.exit(1)


if __name__=="__main__":
    main()
