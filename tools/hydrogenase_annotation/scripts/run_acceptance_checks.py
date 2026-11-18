#!/usr/bin/env python3
"""Record actual per-round executions and source/model hashes without declaring overall success."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from annotate import ROOT


def state_hashes():
    files=[*sorted((ROOT/"scripts").glob("*.py")),*sorted((ROOT/"tests").glob("*.py")),
           *sorted((ROOT/"config").glob("*.json"))]
    files.extend(ROOT/"references/processed"/name for name in ["references.json","references.faa","reference_source_proteins.faa","catalytic_families.hmm",
        "build_manifest.json","profile_manifest.json","maturation_domain_audit/summary.json","maturation_domain_audit/maturation_domains.hmm"])
    files.extend(sorted((ROOT/'references/processed/fes').glob('*')))
    files.extend(sorted((ROOT/'references/processed/active_centers').glob('*')))
    files.extend(sorted((ROOT/'references/processed/split_centers').glob('*')))
    files.extend(sorted((ROOT/'references/processed/catalytic_confirmation').glob('*')))
    if (ROOT/'references/processed/localization/manifest.json').exists():files.append(ROOT/'references/processed/localization/manifest.json')
    fes_policy=json.loads((ROOT/'config/fes_site_review.json').read_text())
    files.extend(ROOT/name for name in fes_policy.get('source_files',{}))
    files.extend(ROOT/c['file'] for s in fes_policy['sites'] for c in s.get('coordination_corroboration',[]))
    files.append(ROOT/'references/raw/structures/download_manifest.json')
    files.extend(ROOT/'references/raw/structures'/(s['pdb_id']+'.cif') for s in json.loads((ROOT/'config/structure_sources.json').read_text())['sources'])
    challenge=ROOT/'references/validation/round21_structural_challenge_development_v4'
    if (challenge/'freeze.json').exists():
        files.append(challenge/'freeze.json')
        files.extend(challenge/name for name in json.loads((challenge/'freeze.json').read_text())['artifacts'])
        files.extend(ROOT/name for name in json.loads((challenge/'freeze.json').read_text()).get('source_files',{}))
    calibration=ROOT/'references/validation/round20_fes_calibration'
    if (calibration/'freeze.json').exists():
        freeze=json.loads((calibration/'freeze.json').read_text())
        files.append(calibration/'freeze.json')
        files.extend(calibration/name for name in freeze['artifacts'])
        files.extend(ROOT/name for name in freeze['source_files'])
        files.extend(ROOT/name for name in freeze['metadata_files'])
        files.append(ROOT/'references/raw/round20_fes_calibration/download_manifest.json')
    split_policy=json.loads((ROOT/'config/split_centers.json').read_text())
    files.extend(ROOT/name for name in split_policy['source_files'])
    files.extend(f for f in sorted((ROOT/'references/processed/reference_domain_audit').glob('*')) if f.is_file() and f.suffix in ['.json','.hmm','.tsv','.faa'])
    files.extend(f for f in sorted((ROOT/'references/curation').glob('*')) if f.is_file())
    files.append(ROOT/'references/raw/catalytic_structures/download_manifest.json')
    files.extend(ROOT/'references/structural_sites'/name for name in ['templates.json','reserved_validation.json'])
    files.extend(sorted((ROOT/'references/discovery_training').glob('*')))
    files.extend(sorted((ROOT/'references/processed/complex_extensions').glob('*')))
    architecture=ROOT/'config/module_architecture.json'
    if architecture.exists():files.extend(ROOT/name for name in json.loads(architecture.read_text())['source_files'])
    literature=ROOT/'config/noncanonical_center_evidence.json'
    if literature.exists():files.extend(ROOT/name for name in json.loads(literature.read_text())['source_files'])
    function=ROOT/'config/round26_function_validation.json'
    if function.exists():files.extend(ROOT/name for name in json.loads(function.read_text())['source_files'])
    distributed=ROOT/'config/round30_distributed_function.json'
    if distributed.exists():files.extend(ROOT/name for name in json.loads(distributed.read_text())['source_files'])
    sampling=ROOT/'config/mag_pilot_sampling.json'
    if sampling.exists():files.append(ROOT/json.loads(sampling.read_text())['population_file'])
    pilot=ROOT/'references/validation/round30_original_MAG_pilot'
    if (pilot/'freeze.json').exists():
        files.append(pilot/'freeze.json')
        files.extend(pilot/name for name in json.loads((pilot/'freeze.json').read_text())['artifacts'])
    files.extend(sorted((ROOT/'experiments/round26').glob('*.py')))
    native_policy=ROOT/'config/round24_ralbus_validation.json'
    if native_policy.exists():files.extend(ROOT/name for name in json.loads(native_policy.read_text())['source_files'])
    natural=ROOT/'references/validation/round22_natural_complex_development'
    if (natural/'freeze.json').exists():
        freeze=json.loads((natural/'freeze.json').read_text())
        files.append(natural/'freeze.json')
        files.extend(natural/name for name in freeze['artifacts'])
        files.extend(ROOT/name for name in freeze['sources'])
    return {str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args();output=args.out.resolve()
    if (output/"execution.json").exists():
        parser.error("Use a new output directory to preserve previous round evidence")
    output.mkdir(parents=True,exist_ok=True)
    jobs=[("unit",["-m","unittest","discover","-s",str(ROOT/"tests"),"-v"]),
          ("structure_extraction",[str(ROOT/"scripts/extract_structure_sites.py")]),
          ("structural_sources",[str(ROOT/"scripts/validate_structure_sources.py")]),
          ('active_centers',[str(ROOT/'scripts/validate_active_centers.py'),'--out',str(output/'active_centers')]),
          ('split_centers',[str(ROOT/'scripts/validate_split_centers.py'),'--out',str(output/'split_centers')]),
          ('catalytic_confirmation',[str(ROOT/'scripts/validate_catalytic_confirmation.py'),'--out',str(output/'catalytic_confirmation')]),
          ('confirmation_transfer',[str(ROOT/'scripts/audit_confirmation_transfer.py'),'--out',str(output/'confirmation_transfer')]),
          ("fes_prediction",[str(ROOT/'scripts/validate_fes_prediction.py'),'--out',str(output/'fes_prediction')]),
          ('fes_challenge',[str(ROOT/'scripts/evaluate_fes_challenge.py'),'--corpus',str(ROOT/'references/validation/round21_structural_challenge_development_v4'),'--out',str(output/'fes_challenge')]),
          ('fes_calibration',[str(ROOT/'scripts/evaluate_fes_challenge.py'),'--corpus',str(ROOT/'references/validation/round20_fes_calibration'),'--out',str(output/'fes_calibration')]),
          ('localization_reference',[str(ROOT/'scripts/validate_localization.py'),'--localization',str(ROOT/'results/round21_localization_reference_verified'),'--out',str(output/'localization_reference')]),
          ("maturation_domains",[str(ROOT/"scripts/audit_maturation_domains.py"),"--out",str(output/"maturation_domains")]),
          ('reference_domains',[str(ROOT/'scripts/audit_reference_domains.py'),'--out',str(output/'reference_domains')]),
          ('module_graph',[str(ROOT/'scripts/validate_module_graph.py'),'--out',str(output/'module_graph')]),
          ('partner_associations',[str(ROOT/'scripts/validate_partner_associations.py'),'--out',str(output/'partner_associations')]),
          ('paralog_contexts',[str(ROOT/'scripts/validate_paralog_contexts.py'),'--out',str(output/'paralog_contexts')]),
          ('complex_extensions',[str(ROOT/'scripts/validate_complex_extensions.py'),'--out',str(output/'complex_extensions')]),
          ('natural_complexes',[str(ROOT/'scripts/evaluate_natural_complexes.py'),'--corpus',str(ROOT/'references/validation/round22_natural_complex_development'),'--out',str(output/'natural_complexes')]),
          ('ralbus_context',[str(ROOT/'scripts/validate_ralbus_context.py'),'--out',str(output/'ralbus_context')]),
          ('function_mechanism',[str(ROOT/'scripts/validate_function_mechanism.py'),'--out',str(output/'function_mechanism')]),
          ('distributed_function',[str(ROOT/'scripts/validate_distributed_function.py'),'--out',str(output/'distributed_function')]),
          ('pilot_sampling',[str(ROOT/'scripts/prepare_mag_pilot.py'),'--out',str(output/'pilot_sampling')]),
          ("regression",[str(ROOT/"scripts/validate_workflow.py"),"--out",str(output/"regression")]),
          ("input_safety",[str(ROOT/"scripts/validate_input_safety.py"),"--out",str(output/"input_safety")]),
          ("discovery",[str(ROOT/"scripts/validate_discovery_model.py"),"--out",str(output/"discovery")]),
          ('output_contract',[str(ROOT/'scripts/validate_output_contract.py'),'--round',str(output)])]
    before=state_hashes();results=[]
    env=dict(os.environ);env["PATH"]=str(ROOT/".env/bin")+os.pathsep+env.get("PATH","")
    report=dict(status="running",before_hashes=before,executions=results)
    record=output/"execution.json";record.write_text(json.dumps(report,indent=2)+"\n")
    for name,command in jobs:
        start=time.monotonic();cmd=[sys.executable,*command]
        with (output/(name+".log")).open("w") as log:
            process=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
        results.append(dict(name=name,command=cmd,exit_code=process.returncode,elapsed_seconds=round(time.monotonic()-start,3),
            log_sha256=hashlib.sha256((output/(name+".log")).read_bytes()).hexdigest()))
        report["after_hashes"]=state_hashes()
        report["status"]="running"
        record.write_text(json.dumps(report,indent=2)+"\n")
        print(name,process.returncode,flush=True)
    report["status"]="checks_passed" if all(r["exit_code"]==0 for r in results) and before==report["after_hashes"] else "failed_or_state_changed"
    report["overall_acceptance"]="not_determined_by_regression_runner"
    record.write_text(json.dumps(report,indent=2)+"\n")
    if report["status"]!="checks_passed":
        sys.exit(1)


if __name__=="__main__":
    main()
