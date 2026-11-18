#!/usr/bin/env python3
"""Measure a frozen new structural challenge without training on its sequences or coordinates."""
import argparse
import json
from pathlib import Path
from annotate import ROOT,write_tsv
from catalytic_confirmation import sha
from fes_sites import annotate_fes
from validate_fes_prediction import compare
from run_acceptance_checks import state_hashes


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--corpus',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    corpus=args.corpus.resolve();freeze=json.loads((corpus/'freeze.json').read_text())
    if freeze['status'] not in ['development_sources_and_truth_frozen_before_prediction','development_truth_coverage_amendment_after_initial_evaluation','development_atom_coverage_amendment_after_initial_evaluation','development_primary_literature_coordination_amendment_after_initial_evaluation','calibration_sources_and_truth_frozen_before_prediction']:
        raise ValueError('Expected frozen challenge sources or explicit development truth amendment')
    for name,digest in freeze['artifacts'].items():
        if sha(corpus/name)!=digest:raise ValueError('Challenge source changed: '+name)
    for name,digest in freeze.get('source_files',{}).items():
        if sha(ROOT/name)!=digest:raise ValueError('Challenge raw coordinates changed: '+name)
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False);before=state_hashes()
    sources=json.loads((corpus/'sources.json').read_text());predictions,_=annotate_fes(corpus/'proteins.faa',out,threads=2)
    rows,metrics=compare(predictions,sources);write_tsv(out/'site_comparisons.tsv',rows)
    per_source={}
    for source in sources:
        _,per_source[source['pdb_id']]=compare([p for p in predictions if p['protein_id'].startswith(source['pdb_id']+'__')],[source])
    after=state_hashes()
    if before!=after:raise ValueError('Predictor changed during structural challenge')
    calibration=freeze['status']=='calibration_sources_and_truth_frozen_before_prediction'
    policy=json.loads((ROOT/'config/round19_structural_challenge.json').read_text());target=policy['planned_development_target']
    if calibration:
        policy=json.loads((ROOT/'config/round20_fes_validation.json').read_text())
        if sha(ROOT/'config/round20_fes_validation.json')!=freeze['protocol_sha256']:raise ValueError('Calibration protocol changed')
        target=dict(precision=policy['structural_precision_target'],min_evaluable_site_recall=policy['min_evaluable_site_recall'])
    target_met=bool(metrics['precision'] is not None and metrics['precision']>=target['precision'] and metrics['exact_site_recall']>=target['min_evaluable_site_recall'])
    summary=dict(status='calibration_challenge_evaluated' if calibration else 'development_challenge_evaluated',development_target_met=target_met,metrics=metrics,per_source=per_source,
        before_hashes=before,after_hashes=after,source_freeze_sha256=sha(corpus/'freeze.json'),gold_source_files=sources,
        artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file()},
        limits=['Source-selected calibration structures, not a sealed final test or population-accuracy estimate.' if calibration else 'Three source-selected development structures, not a sealed final test or population-accuracy estimate.',
            'No independent Huc query; all calibrated high FeS and full-function acceptance remains pending.',
            'Raw first-shell coordinates and asymmetric-unit copies are retained; cofactor occupancy in queries is unknown.'])
    if calibration:
        from Bio import SeqIO
        audit=json.loads((corpus/'source_audit.json').read_text())
        query_ids={r.id for r in SeqIO.parse(corpus/'proteins.faa','fasta')}
        eligible={r['protein_id'] for r in audit if r['declared_split']=='calibration' and r['independent_eligible']} & query_ids
        subset=[]
        for source in sources:
            sites=[s for s in source['sites'] if all(source['pdb_id']+'__entity_'+r['entity_id'] in eligible for r in s['ligand_residues'])]
            subset.append(dict(source,sites=sites))
        _,summary['operationally_independent_calibration_metrics']=compare([p for p in predictions if p['protein_id'] in eligible],subset)
        summary['operationally_independent_protein_ids']=sorted(eligible)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(dict(status=summary['status'],target_met=target_met,metrics=metrics,per_source=per_source),indent=2))


if __name__=='__main__':main()
