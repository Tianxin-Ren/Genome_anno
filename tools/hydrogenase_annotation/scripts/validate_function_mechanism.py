#!/usr/bin/env python3
"""Evaluate scoped native reaction truth and a rejected specificity experiment; never train production."""
import argparse
import copy
import json
import sys
from datetime import datetime,timezone
from pathlib import Path
from annotate import ROOT,annotate,extract_genbank,load_blast,write_tsv
from catalytic_confirmation import sha,table
from complex_extensions import load_contracts,load_bank
from prepare_natural_complexes import anonymous_genes
from run_acceptance_checks import state_hashes
sys.path.insert(0,str(ROOT/'experiments/round26'))
from nucleotide_specificity import load_policy,classify


def verify_sources(policy):
    for name,digest in policy['source_files'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Function challenge source changed: '+name)


def binding_experiment(policy):
    spec,seqs,_=load_policy();raw=ROOT/'references/raw/round26_function_mechanism'
    thresholds=dict(blast_evalue=1e-5,min_identity=0,min_query_coverage=0,min_reference_coverage=0)
    old=json.loads((raw/'assayed_modules.json').read_text());truth={r['id']:('NADP_reference_like' if r['nucleotide']=='NADP' else 'NAD_reference_like') for r in old}
    truth['MOTHE_c19220']='NADP_reference_like'
    hits=load_blast(raw/'binding_alignments.tsv',thresholds)+load_blast(raw/'Moorella_binding_alignments.tsv',thresholds)
    predicted={r['protein_id']:r for r in classify(hits,spec,seqs)}
    rows=[dict(protein_id=q,expected=expected,prediction=predicted.get(q,{}).get('prediction','unknown'),
        source_group='new_source_challenge' if q=='MOTHE_c19220' else 'exposed_development',
        correct=predicted.get(q,{}).get('prediction')==expected,covered=predicted.get(q,{}).get('prediction','unknown')!='unknown') for q,expected in truth.items()]
    covered=sum(r['covered'] for r in rows);correct=sum(r['correct'] for r in rows)
    metrics=dict(correct=correct,wrong=covered-correct,unassigned=len(rows)-covered,samples=len(rows),precision=correct/covered if covered else None,coverage=covered/len(rows),new_source_correct=next(r['correct'] for r in rows if r['source_group']=='new_source_challenge'))
    nearest=[];module_rows=json.loads((raw/'competition_modules.json').read_text());labels={r['id']:r['nucleotide'] for r in module_rows}
    hh=load_blast(raw/'competition_alignments.tsv',thresholds)
    for q,label in labels.items():
        choices=[h for h in hh if h['qseqid']==q and h['sseqid']!=q]
        best=max(choices,key=lambda h:h['bitscore']);nearest.append(dict(protein_id=q,expected=label,prediction=labels[best['sseqid']],best_source=best['sseqid'],bitscore=best['bitscore'],correct=label==labels[best['sseqid']]))
    targets=policy['targets'];passed=bool(metrics['precision'] is not None and metrics['precision']>=targets['nucleotide_precision'] and metrics['coverage']>=targets['nucleotide_coverage'] and metrics['new_source_correct'])
    return dict(status='candidate_preference_model_rejected' if not passed else 'development_only_targets_met',metrics=metrics,rows=rows,
        nearest_reference_rows=nearest,nearest_reference_correct=sum(r['correct'] for r in nearest),
        production_enabled=False,independent_validation=False,rule_before_new_source_sha256=sha(ROOT/'references/model_snapshots/round26_binding_rule_before_new_source/nucleotide_specificity.json'),
        failure='New source betaEFP motif resembles NAD refs despite purifiedNADP specificity. Full-module nearest reference also predictsNAD. Protein-only patterns do not resolve complete-complex specificity.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);args=p.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    protocol=ROOT/'config/round26_function_validation.json';policy=json.loads(protocol.read_text());verify_sources(policy);before=state_hashes()
    contract={c['contract_id']:c for c in load_contracts()};prepared=[]
    replays=set(load_bank()[3]['reference_replay_cases'])
    for case in policy['cases']:
        native=extract_genbank(ROOT/case['genome'])[0];by_id={g['id']:g for g in native}
        if any(n not in by_id for n in case['members']):raise ValueError('Missing exact native member locus')
        genes,mapping=anonymous_genes(native);prepared.append((case,genes,mapping))
    freeze=dict(status='native_development_truth_frozen_before_prediction',created_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=sha(protocol),source_files=policy['source_files'],cases=policy['cases'],
        extension_model_manifest_sha256=sha(ROOT/'references/processed/complex_extensions/manifest.json'),reference_replay_case_ids=sorted(replays&{c['id'] for c in policy['cases']}),
        source_mappings={c['id']:m for c,g,m in prepared},independence='Exposed development only. Experimental new-source binding prediction was frozen separately and failed; never recast as final test.')
    (out/'freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
    execution=dict(status='running',before_hashes=before,runs=[]);save=lambda:(out/'execution.json').write_text(json.dumps(execution,indent=2)+'\n');save();results=[]
    for case,genes,mapping in prepared:
        run=out/case['id'];run.mkdir();annotations=annotate(copy.deepcopy(genes),run,threads=2,mode='genome',input_files=[ROOT/case['genome']])
        q=mapping[case['anchor']];a=next((a for a in annotations if a['protein_id']==q),{});expected={mapping[n] for n in case['members']}
        members=table(run/'complex_members.tsv');chosen={m['protein_id'] for m in members if m['anchor']==q and m['required_for_structure']=='True'}
        complete=a.get('structure_status') in ['complete_supported','core_complete_accessory_missing','composition_supported_reorganized'];exact=bool(complete and chosen==expected)
        assigned={contract[f]['reaction_key'] for f in a.get('potential_function_ids','').split(';') if f};assigned={'H2_NAD_reduction_without_required_Fd' if k=='H2_NAD_reduction' else k for k in assigned}
        positive=set(case['positive']);negative=set(case['negative']);correct=assigned&positive if exact else set()
        fp=(assigned&negative)|((assigned&positive) if not exact else set())
        rows=table(run/'complex_hypotheses.tsv');hypotheses=[r for r in rows if r['anchor']==q]
        results.append(dict(case_id=case['id'],proteins=len(genes),expected_members=sorted(expected),predicted_members=sorted(chosen),
            core=dict(TP=int(exact),FP=int(complete and not exact),FN=int(not exact)),
            reaction=dict(TP=len(correct),FP=len(fp),FN=len(positive-correct),positive_truth=sorted(positive),negative_truth=sorted(negative),assigned=sorted(assigned),ungraded_extra=sorted(assigned-positive-negative)),
            anchor=a or 'not_retained',hypotheses=hypotheses,source=case['primary_source'],extension_reference_replay=case['id'] in replays,independent_validation=False,
            noncanonical_literature_evidence=[r for r in table(run/'noncanonical_centers/alternatives.tsv') if r['protein_id']==q]))
        execution['runs'].append(dict(case_id=case['id'],manifest_sha256=sha(run/'run_manifest.json'),proteins=len(genes)));save();print(case['id'],len(genes),results[-1]['core'],results[-1]['reaction'],flush=True)
    boundaries=[];boundary_protocol=ROOT/'config/round27_six_subunit_validation.json'
    if boundary_protocol.exists():
        boundary_policy=json.loads(boundary_protocol.read_text())
        spec=next(s for s in json.loads((ROOT/'config/complex_extension_sources.json').read_text())['systems'] if s['case_id']==boundary_policy['case_id'])
        case,genes,mapping=next(p for p in prepared if p[0]['id']==boundary_policy['case_id'])
        roles={m['role']:mapping[m['source_locus']] for m in spec['native_case']['members']}
        for name in boundary_policy['boundary_runs']:
            altered=copy.deepcopy(genes)
            if name.startswith('delete_'):altered=[g for g in altered if g['id']!=roles[name.removeprefix('delete_')]]
            elif name=='HydB_one_residue_variant':
                gene=next(g for g in altered if g['id']==roles['HydB'])
                gene['sequence']=('A' if gene['sequence'][0]!='A' else 'G')+gene['sequence'][1:]
            else:raise ValueError('Unreviewed six-subunit perturbation')
            run=out/'six_subunit_boundaries'/name;run.mkdir(parents=True)
            annotations=annotate(altered,run,threads=2,mode='genome',input_files=[ROOT/case['genome']])
            assigned=[a['protein_id'] for a in annotations if boundary_policy['contract_id'] in a['potential_function_ids'].split(';')]
            boundaries.append(dict(case=name,reference_reaction_anchors=assigned,passed=not assigned,run_manifest_sha256=sha(run/'run_manifest.json')))
            print('six_subunit_boundary',name,not assigned,flush=True)
    fusion_boundaries=[];fusion_protocol=ROOT/'config/round29_fusion_validation.json'
    if fusion_protocol.exists():
        fusion_policy=json.loads(fusion_protocol.read_text())
        for name,digest in fusion_policy['source_files'].items():
            if sha(ROOT/name)!=digest:raise ValueError('Fusion boundary source changed: '+name)
        case,genes,mapping=next(p for p in prepared if p[0]['id']==fusion_policy['case_id'])
        q=mapping[fusion_policy['anchor']];partner=mapping[fusion_policy['partner']]
        for name in fusion_policy['boundary_runs']:
            altered=copy.deepcopy(genes);by_id={g['id']:g for g in altered};beta=by_id[partner];mode='genome'
            if name=='delete_partner':altered=[g for g in altered if g['id']!=partner]
            elif name=='delete_gamma_region':beta['sequence']=beta['sequence'][152:]
            elif name=='delete_beta_region':beta['sequence']=beta['sequence'][:176]
            elif name=='reverse_fused_order':beta['sequence']=beta['sequence'][176:585]+'GSGSG'+beta['sequence'][14:152]
            elif name=='partial_partner':beta['partial']=True
            elif name=='cross_bin_partner':
                for gene in altered:gene['bin_id']='bin_one'
                beta['bin_id']='bin_two';beta['contig']='other_bin_contig';mode='bin'
            elif name=='ambiguous_anchor_ligand':
                sequence=by_id[q]['sequence'];by_id[q]['sequence']=sequence[:307]+'X'+sequence[308:]
            else:raise ValueError('Unreviewed natural fusion perturbation')
            run=out/'fusion_boundaries'/name;run.mkdir(parents=True)
            annotations=annotate(altered,run,threads=2,mode=mode,input_files=[ROOT/case['genome']])
            anchor=next((a for a in annotations if a['protein_id']==q),{})
            complete=anchor.get('structure_status') in ['complete_supported','core_complete_accessory_missing','composition_supported_reorganized']
            passed=not complete and not anchor.get('potential_function_ids')
            fusion_boundaries.append(dict(case=name,passed=passed,complete_core=complete,anchor=anchor,
                run_manifest_sha256=sha(run/'run_manifest.json'),scope='Synthetic perturbation of the full native genomic protein/ownership context, not a new natural organism.'))
            print('fusion_boundary',name,passed,flush=True)
    model=binding_experiment(policy);after=state_hashes();verify_sources(policy)
    if before!=after:raise ValueError('Function experiment code/models/sources changed during prediction')
    summary=dict(status='development_function_challenges_evaluated',before_hashes=before,after_hashes=after,protocol_sha256=sha(protocol),freeze_sha256=sha(out/'freeze.json'),
        native_results=results,native_proteins=sum(r['proteins'] for r in results),binding_experiment=model,production_predictor_changed=sha(ROOT/'references/processed/complex_extensions/manifest.json')!=json.loads((ROOT/'revisions/round26_checked/snapshot.json').read_text())['files']['references/processed/complex_extensions/manifest.json'],independent_validation=False,
        reference_replay_case_ids=freeze['reference_replay_case_ids'],six_subunit_boundary_results=boundaries,fusion_boundary_results=fusion_boundaries,
        native_metrics={scope:{k:sum(r[scope][k] for r in results) for k in ['TP','FP','FN']} for scope in ['core','reaction']},limits=[v for v in policy['limits'] if not v.startswith('R25productionpredictorunchanged')]+['Post-identity structure bank can evolve; exposed cases used as templates are explicitly reference replays. Failed nucleotide classifier remains excluded.'],protocol_legacy_limits=policy['limits'],targets=policy['targets'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');write_tsv(out/'binding_comparisons.tsv',model['rows']);execution.update(status='completed',after_hashes=after,summary_sha256=sha(out/'summary.json'));save()
    print(json.dumps(dict(status=summary['status'],native=summary['native_metrics'],binding=model['metrics'],production_enabled=False),indent=2))
    if any(not b['passed'] for b in boundaries):raise SystemExit('Six-subunit boundary emitted a reference reaction')
    if any(not b['passed'] for b in fusion_boundaries):raise SystemExit('Fusion boundary emitted a complete core or reaction')


if __name__=='__main__':main()
