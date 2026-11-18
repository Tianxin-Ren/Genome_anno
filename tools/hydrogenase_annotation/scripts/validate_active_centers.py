#!/usr/bin/env python3
"""Actual source-withheld structural mapping and explicit sequence perturbations."""
import argparse
import json
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, write_tsv
from active_centers import annotate_active_centers, apply_site_guard, structural_record, sha
from validate_fes_prediction import wilson


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();output=args.out.resolve()
    if (output/'summary.json').exists():parser.error('Use a new directory')
    output.mkdir(parents=True,exist_ok=True)
    directory=ROOT/'references/processed/active_centers'
    models=json.loads((directory/'models.json').read_text());truth=json.loads((directory/'structural_truth.json').read_text())
    tracked=[ROOT/'scripts/active_centers.py',Path(__file__),ROOT/'config/active_center_review.json',
        ROOT/'config/active_center_sources.json',directory/'manifest.json',directory/'models.json',directory/'structural_truth.json']
    before={str(p.relative_to(ROOT)):sha(p) for p in tracked}
    # Source split was declared before predictions; it is not a genus/homolog blind split.
    (output/'pre_prediction_plan.json').write_text(json.dumps(dict(source_plan_sha256=sha(ROOT/'config/active_center_sources.json'),
        sources=[dict(pdb_id=r['pdb_id'],entity_id=r['entity_id'],split=r['split']) for r in truth],
        target=.98,gate_status='development_only_requires_independent_taxonomic_homolog_audit',before_hashes=before),indent=2)+'\n')
    checks=[];comparisons=[];metrics={}
    def check(name,passed,detail):checks.append(dict(test=name,passed=bool(passed),detail=json.dumps(detail)))
    for split in ['template','source_withheld_evaluation','development_control']:
        run=output/split;run.mkdir(exist_ok=True);sources=[r for r in truth if r['split']==split]
        fasta=run/'proteins.faa';SeqIO.write([SeqRecord(Seq(r['sequence']),id=r['pdb_id']+'__entity_'+r['entity_id'],description=split) for r in sources],fasta,'fasta')
        predictions,_=annotate_active_centers(fasta,run);by_query={r['protein_id']:r for r in predictions}
        correct=total=0;per_source={}
        for source in sources:
            qid=source['pdb_id']+'__entity_'+source['entity_id'];row=by_query.get(qid,{})
            # Re-extract coordinates from pinned original bytes, not from predictions.
            actual=structural_record(ROOT/'references/raw/catalytic_structures'/(source['pdb_id']+'.cif'),source['entity_id'],source['family'])
            assert actual['source_sha256']==source['source_sha256'] and actual['ligands']==source['ligands']
            ligands=[s for s in row.get('sites',[]) if s['residue_role']=='direct_metal_ligand']
            observed={(s['query_position'],s['observed_residue']) for s in ligands if s['status']=='present'}
            coordinates={(int(p),r) for p,r in actual['ligands'].items()}
            complete=len(coordinates)==(1 if source['family']=='Fe' else 4)
            if complete:
                count=len(coordinates&observed);correct+=count;total+=len(coordinates)
                check('ligand_coordinate_exact:'+qid,coordinates==observed,dict(expected=sorted(coordinates),observed=sorted(observed)))
            else:count=None
            per_source[qid]=dict(eligible_complete_structural_truth=complete,expected=sorted(coordinates),observed=sorted(observed),
                correct=count,active_center_status=row.get('active_center_status','no_mapping'),reference_id=row.get('reference_id',''))
            comparisons.append(dict(protein_id=qid,split=split,**per_source[qid]))
        metrics[split]=dict(correct_ligand_positions=correct,eligible_ligand_positions=total,accuracy=correct/total if total else None,
            wilson95=wilson(correct,total),per_source=per_source)
        if split=='development_control':
            mutant=by_query.get('3F46__entity_1',{});paralog=by_query.get('4YT2__entity_1',{})
            check('real_Hmd_C176A_not_complete_ligand_supported',mutant.get('active_center_status')=='ligands_substituted',mutant.get('active_center_status'))
            check('real_HmdII_keeps_Cys_but_uncovered_signature_stays_uncertain',paralog.get('active_center_status')=='insufficient_alignment' and
                paralog.get('direct_ligand_status')=='present' and paralog.get('reaction_signature_status')=='not_covered',
                dict(status=paralog.get('active_center_status'),ligands=paralog.get('direct_ligand_status'),signature=paralog.get('reaction_signature_status')))
            guarded=[dict(protein_id='4YT2__entity_1',catalytic_status='supported',structure_confidence='high',potential_function_ids='Hmd_function',flags='')]
            apply_site_guard(guarded,predictions)
            check('real_HmdII_full_model_blocks_legacy_function',guarded[0]['catalytic_status']=='candidate' and not guarded[0]['potential_function_ids'],guarded)
    run=output/'perturbations';run.mkdir(exist_ok=True);records=[];expected=[]
    for model in models:
        sequence=model['sequence'];ligands=[s for s in model['sites'] if s['residue_role']=='direct_metal_ligand']
        for site in ligands:
            position=site['reference_position']
            for residue,state in [('A','substituted'),('X','ambiguous')]+([('U','present')] if model['family']=='NiFe' else []):
                qid=model['reference_id']+f'__{position}_{residue}'
                records.append(SeqRecord(Seq(sequence[:position-1]+residue+sequence[position:]),id=qid,description='synthetic_residue_state_control'))
                expected.append(dict(protein_id=qid,position=position,expected_state=state,reference_id=model['reference_id'],
                    expected_query_position=position,expected_residue=residue))
        last=ligands[-1]['reference_position'];qid=model['reference_id']+'__Cterminal_truncation'
        records.append(SeqRecord(Seq(sequence[:last-1]),id=qid,description='synthetic_truncation_control'))
        expected.append(dict(protein_id=qid,position=last,expected_state='not_covered',reference_id=model['reference_id'],expected_query_position=None,expected_residue=None))
    cpi=next(m for m in models if m['pdb_id']=='1FEH');qid='CpI_context_C299A'
    records.append(SeqRecord(Seq(cpi['sequence'][:298]+'A'+cpi['sequence'][299:]),id=qid,description='synthetic_context_not_inactivation_truth'))
    variant=json.loads((directory/'manifest.json').read_text())['exact_active_reference_variants'][0]
    protein=next(r for r in SeqIO.parse(ROOT/variant['source_fasta'],'fasta') if r.id==variant['reference_id'])
    records.append(SeqRecord(protein.seq,id='real_active_Swol_Hyd1A',description='exact_biochemical_reference_noncanonical_site'))
    edited=str(protein.seq);edited=edited[:10]+('A' if edited[10]!='A' else 'G')+edited[11:]
    records.append(SeqRecord(Seq(edited),id='edited_Swol_not_exact',description='variant_exception_must_not_generalize'))
    hmd=next(m for m in models if m['family']=='Fe');sequence=hmd['sequence']
    records.append(SeqRecord(Seq(sequence[:13]+'Y'+sequence[14:]),id='Hmd_H14Y',description='synthetic_canonical_reaction_signature_control'))
    records.append(SeqRecord(Seq(cpi['sequence'][:374]+cpi['sequence'][478:]),id='CpI_scaffold_gap_outside_ligand_contexts',
        description='synthetic_domain_deletion_not_biochemical_truth'))
    fasta=run/'proteins.faa';SeqIO.write(records,fasta,'fasta')
    (run/'expected_states.json').write_text(json.dumps(expected,indent=2)+'\n')
    predictions,_=annotate_active_centers(fasta,run);by_query={r['protein_id']:r for r in predictions}
    for expected_row in expected:
        qid=expected_row['protein_id'];run_row=by_query.get(qid,{})
        row=next((s for s in run_row.get('sites',[]) if s['reference_id']==expected_row['reference_id'] and s['reference_position']==expected_row['position']),{})
        ok=row.get('status')==expected_row['expected_state'] and row.get('query_position')==expected_row['expected_query_position']
        if expected_row['expected_residue'] is not None:ok &= row.get('observed_residue')==expected_row['expected_residue']
        check('perturbation:'+qid,ok,dict(expected=expected_row,observed=row))
    context=by_query.get('CpI_context_C299A',{})
    check('neighbor_C299A_does_not_block_ligand_set',context.get('active_center_status')=='template_ligand_set_supported',context.get('active_center_status'))
    check('active_Swol_exact_reference_retains_experimental_boundary',by_query.get('real_active_Swol_Hyd1A',{}).get('active_center_status')=='experimentally_active_exact_reference_noncanonical_site',by_query.get('real_active_Swol_Hyd1A',{}).get('active_center_status'))
    check('Swol_exception_does_not_generalize_to_unverified_alleles',by_query.get('edited_Swol_not_exact',{}).get('active_center_status')=='ligands_substituted',by_query.get('edited_Swol_not_exact',{}).get('active_center_status'))
    check('covered_Hmd_H14Y_reported_as_changed_signature',by_query.get('Hmd_H14Y',{}).get('active_center_status')=='canonical_reaction_signature_changed',by_query.get('Hmd_H14Y',{}).get('active_center_status'))
    deletion=by_query.get('CpI_scaffold_gap_outside_ligand_contexts',{})
    check('ligand_contexts_and_whole_scaffold_quality_are_separate',deletion.get('active_center_status')=='template_ligand_set_supported' and
        deletion.get('template_span_alignment_status')=='insufficient' and deletion.get('catalytic_domain_integrity')=='not_established_by_active_site_mapping',
        {k:v for k,v in deletion.items() if k!='sites'})
    guarded=[dict(protein_id='CpI_scaffold_gap_outside_ligand_contexts',catalytic_status='supported',structure_confidence='high',potential_function_ids='CpI_function',flags='')]
    apply_site_guard(guarded,predictions)
    check('local_ligands_do_not_rescue_scaffold_deletion_function',guarded[0]['function_status']=='unassigned' and not guarded[0]['potential_function_ids'] and guarded[0]['structure_confidence']!='high',guarded)
    unchanged=all(sha(ROOT/name)==digest for name,digest in before.items());check('frozen_state_unchanged',unchanged,before)
    def serialized(row):return {k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in row.items()}
    write_tsv(output/'checks.tsv',checks)
    write_tsv(output/'structural_comparisons.tsv',[serialized(r) for r in comparisons])
    manifest=json.loads((directory/'manifest.json').read_text())
    summary=dict(status='development_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),
        passed=sum(c['passed'] for c in checks),failed=[c for c in checks if not c['passed']],metrics=metrics,
        template_proteins=manifest['template_proteins'],model_positions=manifest['positions'],
        before_hashes=before,unchanged=unchanged,complete_catalytic_site_gate='not_met',
        limits=['Four source-withheld proteins do not prove broad independent >=98% state-mapping accuracy.',
            'Structural sources can share homologs or genera with earlier reference/development data; no final-test independence claim.',
            'Direct ligand coordinates are experimental truth; neighbor correctness and full catalytic-domain architecture need independent review.',
            'Perturbations test software states, not activity, actual metal assembly, or a biological inactivation rule.',
            'MBH development control includes four protein ligands after admitting the FCO iron paired to nickel; no new independent sample is implied.'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(dict(status=summary['status'],checks=summary['checks'],passed=summary['passed'],metrics={k:{a:b for a,b in v.items() if a!='per_source'} for k,v in metrics.items()},failed=[c['test'] for c in summary['failed']]),indent=2))
    if summary['status']=='failed':raise SystemExit(1)


if __name__=='__main__':main()
