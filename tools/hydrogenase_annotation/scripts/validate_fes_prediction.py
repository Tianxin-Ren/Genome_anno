#!/usr/bin/env python3
"""Evaluate per-site predictions with source-withheld experimental ligand coordinates."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, write_tsv
from fes_sites import annotate_fes, sha

COORDINATE_TRUTH_REVIEW=ROOT/'references/validation/round52_coordinate_truth_review/review.json'
COORDINATE_TRUTH_REVIEW_SHA256='462dd68dfe8e12524798e69bb9bf4c3e417b9dce5aa67838fa78e7b832df9a78'


def truth_digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def source_truth_signature(source):
    return truth_digest({k:source.get(k) for k in ['sequences','canonical_chains','modeled_protein_positions','modeled_ligand_atoms']})


def site_truth_signature(site):
    return truth_digest(dict(cluster_type=site['cluster_type'],contacts=sorted((r['entity_id'],r['position'],r['residue']) for r in site['ligand_residues']),
        validation_coordination_review=site.get('validation_coordination_review')))


def coordinate_reviews(sources):
    if sha(COORDINATE_TRUTH_REVIEW)!=COORDINATE_TRUTH_REVIEW_SHA256:raise ValueError('Coordinate truth review artifact changed')
    manifest=json.loads(COORDINATE_TRUTH_REVIEW.read_text());reviews={}
    for source in sources:
        pdb=source['pdb_id'];digest=source.get('source_sha256',source.get('sha256'));registered=manifest['sources'].get(pdb)
        if not digest and registered is None:continue
        if registered is None or digest!=registered['raw_sha256']:raise ValueError('Unreviewed or changed raw coordinate source: '+pdb)
        if sha(ROOT/registered['raw_path'])!=digest:raise ValueError('Raw coordinate evidence changed: '+pdb)
        if source_truth_signature(source) not in registered['source_truth_signatures']:raise ValueError('Structural sequence or modeled coverage differs from review: '+pdb)
        for site in source['sites']:
            review=registered['sites'].get(site['site_id'])
            if review is None or site_truth_signature(site) not in review['gold_signatures']:raise ValueError('Coordinate truth site version changed: '+site['site_id'])
            evidence=site.get('validation_coordination_review')
            if evidence and sha(ROOT/evidence['evidence_file'])!=evidence['evidence_sha256']:raise ValueError('Noncanonical primary evidence changed')
            reviews[site['site_id']]=review
    return reviews


def ligand_atoms_covered(prediction,source_atoms):
    required={'C':{'SG'},'U':{'SE'},'H':{'ND1','NE2'},'D':{'OD1','OD2'}}
    return all(residue in required and required[residue]<=set(source_atoms.get(str(position),[])) for position,residue in prediction)


def wilson(success,total):
    if not total:return None
    z=1.959963984540054;p=success/total;den=1+z*z/total
    center=(p+z*z/(2*total))/den
    half=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/den
    return [max(0,center-half),min(1,center+half)]


def compare(predictions,sources):
    truth=[];rows=[];used=set();ungradable=[];coordinate_ungradable=[];reviewed=coordinate_reviews(sources)
    for source in sources:
        for site in source['sites']:
            entities={r['entity_id'] for r in site['ligand_residues']}
            if len(entities)!=1:
                raise ValueError('Cross-entity FeS validation needs explicit multi-chain truth handling')
            entity=next(iter(entities))
            canonical={'[4Fe-4S]':4,'[3Fe-4S]':3,'[4Fe-3S]':6,'[2Fe-2S]':4}[site['cluster_type']]
            expected=canonical
            review=site.get('validation_coordination_review')
            if review:
                if (review.get('status')!='primary_literature_reviewed_noncanonical_ligand_set' or
                    review.get('expected_ligand_count')!=5 or site['cluster_type']!='[2Fe-2S]' or
                    sorted(r['residue'] for r in site['ligand_residues'])!=['C','C','C','C','H'] or
                    len({(r['entity_id'],r['position']) for r in site['ligand_residues']})!=5 or not review.get('evidence_url') or
                    len(review.get('evidence_sha256',''))!=64):raise ValueError('Noncanonical FeS truth lacks a valid explicit coordination review')
                expected=5
            observed=len(site['ligand_residues'])
            if observed==expected:
                coordination_status='explicitly_reviewed_noncanonical' if review else 'canonical_count_compatible'
            else:
                coordination_status='unreviewed_additional_contact_residues' if observed>canonical else 'unreviewed_fewer_contact_residues'
            record=dict(site,protein_id=source['pdb_id']+'__entity_'+entity,
                truth_coordination_status=coordination_status,canonical_contact_count=canonical,
                coordinates={(r['position'],r['residue']) for r in site['ligand_residues']})
            geometry=reviewed.get(site['site_id'])
            record['coordinate_truth_status']='reviewed_invariant_local_contact_set' if geometry else 'unversioned_fixture_no_structural_accuracy_claim'
            if geometry and (geometry['contact_residue_set_variable'] or not geometry['all_local_cluster_compositions_complete']):
                record['coordinate_truth_status']='reference_local_contact_set_not_unique' if geometry['contact_residue_set_variable'] else 'reference_cluster_atom_composition_incomplete'
                record['truth_coordination_status']=record['coordinate_truth_status']
                record['local_contact_options']=[{(p,r) for e,p,r in state['possible_contact_residues']} for state in geometry['local_cluster_states']
                    if state['mandatory_contact_residues']==state['possible_contact_residues'] and all(e==entity for e,p,r in state['possible_contact_residues'])]
                coordinate_ungradable.append(record)
            else:(truth if observed==expected else ungradable).append(record)
    tp=fp=ungraded_supported=unmodeled_supported=state_ungraded_supported=0
    coverage={s['pdb_id']+'__entity_'+entity:set(positions) for s in sources for entity,positions in (s.get('modeled_protein_positions') or {}).items()}
    atom_coverage={s['pdb_id']+'__entity_'+entity:positions for s in sources for entity,positions in (s.get('modeled_ligand_atoms') or {}).items()}
    for prediction in predictions:
        if prediction['status']!='template_supported':continue
        coordinates={(r['query_position'],r['observed_residue']) for r in prediction['ligand_residues']}
        exact=[t for t in truth if t['protein_id']==prediction['protein_id'] and t['coordinates']==coordinates and t['cluster_type']==prediction['cluster_type_predicted']]
        coordination_status=''
        if exact and exact[0]['site_id'] not in used:
            tp+=1;used.add(exact[0]['site_id']);status='TP_exact_type_and_ligand_coordinates';matched=exact[0]['site_id']
            coordination_status=exact[0]['truth_coordination_status']
        else:
            unmodeled=prediction['protein_id'] in coverage and any(p not in coverage[prediction['protein_id']] for p,r in coordinates)
            if prediction['protein_id'] in atom_coverage:
                unmodeled |= not ligand_atoms_covered(coordinates,atom_coverage[prediction['protein_id']])
            partial=[t for t in ungradable if t['protein_id']==prediction['protein_id'] and t['coordinates']<=coordinates and t['cluster_type']==prediction['cluster_type_predicted']]
            uncertain=[t for t in coordinate_ungradable if t['protein_id']==prediction['protein_id'] and t['cluster_type']==prediction['cluster_type_predicted'] and coordinates in t['local_contact_options']]
            if unmodeled:
                unmodeled_supported+=1;status='not_scored_predicted_ligands_outside_modeled_structure';matched=''
            elif uncertain:
                state_ungraded_supported+=1;status='not_scored_reference_local_state_not_unique_or_complete';matched=uncertain[0]['site_id']
                coordination_status=uncertain[0]['truth_coordination_status']
            elif partial:
                ungraded_supported+=1;status='not_scored_unreviewed_coordination';matched=partial[0]['site_id']
                coordination_status=partial[0]['truth_coordination_status']
            else:
                fp+=1;status='FP_type_or_ligand_coordinate_mismatch';matched=''
        rows.append(dict(protein_id=prediction['protein_id'],prediction_id=prediction['fes_site_id'],status=status,
            matched_truth=matched,template_site=prediction['template_site_id'],predicted_type=prediction['cluster_type_predicted'],
            predicted_coordinates=json.dumps(sorted(coordinates)),truth_coordination_status=coordination_status))
    for t in truth:
        if t['site_id'] not in used:
            rows.append(dict(protein_id=t['protein_id'],prediction_id='',status='FN_no_exact_supported_prediction',matched_truth=t['site_id'],
                template_site='',predicted_type='',predicted_coordinates=json.dumps(sorted(t['coordinates'])),truth_coordination_status=t['truth_coordination_status']))
    coordination_audit=[dict(site_id=t['site_id'],protein_id=t['protein_id'],cluster_type=t['cluster_type'],
        status=t['truth_coordination_status'],canonical_contact_count=t['canonical_contact_count'],observed_contact_count=len(t['coordinates']),
        observed_contact_residues=sorted(t['coordinates']),coordinate_truth_status=t['coordinate_truth_status'],
        local_contact_options=[sorted(s) for s in t.get('local_contact_options',[])]) for t in truth+ungradable+coordinate_ungradable]
    total_truth=len(truth)+len(ungradable)+len(coordinate_ungradable)
    metrics=dict(TP=tp,FP=fp,FN=len(truth)-tp,eligible_sites=len(truth),incomplete_truth_sites=len(ungradable)+len(coordinate_ungradable),
        supported_predictions_with_incomplete_truth=ungraded_supported+state_ungraded_supported,all_output_sites=len(predictions),
        unreviewed_coordination_sites=len(ungradable),
        additional_observed_contact_sites=sum(t['truth_coordination_status']=='unreviewed_additional_contact_residues' for t in ungradable),
        fewer_observed_contact_sites=sum(t['truth_coordination_status']=='unreviewed_fewer_contact_residues' for t in ungradable),
        supported_predictions_with_unreviewed_coordination=ungraded_supported,truth_coordination_audit=coordination_audit,
        reference_coordinate_state_ungradable_sites=len(coordinate_ungradable),supported_predictions_with_reference_coordinate_uncertainty=state_ungraded_supported,
        ungradable_truth_sites=len(ungradable)+len(coordinate_ungradable),total_truth_sites=total_truth,eligible_truth_fraction=len(truth)/total_truth if total_truth else None,
        coordinate_reviewed_sites=len(reviewed),coordinate_truth_review=dict(path=str(COORDINATE_TRUTH_REVIEW.relative_to(ROOT)),sha256=COORDINATE_TRUTH_REVIEW_SHA256),
        legacy_metric_notes='incomplete_truth_sites and supported_predictions_with_incomplete_truth are deprecated aliases for ungradable truth, including coordinate-state uncertainty; they do not establish missing coordinates.',
        truth_scope='Source-bound local coordinate states and primary-reviewed exceptions determine eligibility. Unversioned software fixtures cannot establish structural accuracy. Contact counts/distances do not prove chemistry, redox state or query occupancy.',
        supported_predictions_outside_modeled_structure=unmodeled_supported,
        scored_supported_predictions=tp+fp,
        supported_prediction_gradability=(tp+fp)/sum(p['status']=='template_supported' for p in predictions) if any(p['status']=='template_supported' for p in predictions) else None,
        supported_output_sites=sum(p['status']=='template_supported' for p in predictions),
        precision=tp/(tp+fp) if tp+fp else None,precision_wilson95=wilson(tp,tp+fp),
        exact_site_recall=tp/len(truth) if truth else None,recall_wilson95=wilson(tp,len(truth)),
        abstained_or_ambiguous=sum(p['status']!='template_supported' for p in predictions),
        confidence_interval_scope='Descriptive site-level interval; sites within proteins/species are correlated, not independent biological replicates')
    return rows,metrics


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();output=args.out.resolve()
    if (output/'summary.json').exists():parser.error('Use a fresh evaluation directory')
    output.mkdir(parents=True,exist_ok=True)
    templates=json.loads((ROOT/'references/structural_sites/templates.json').read_text())
    reserved=json.loads((ROOT/'references/structural_sites/reserved_validation.json').read_text())
    if {s['pdb_id'] for s in templates}&{s['pdb_id'] for s in reserved}:raise ValueError('Structure split leakage')
    policy=json.loads((ROOT/'config/fes_site_review.json').read_text());metrics={};all_rows=[]
    before={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'scripts/fes_sites.py',ROOT/'scripts/validate_fes_prediction.py',
        ROOT/'config/fes_site_review.json',ROOT/'references/processed/fes/manifest.json',
        ROOT/'references/structural_sites/templates.json',ROOT/'references/structural_sites/reserved_validation.json',COORDINATE_TRUTH_REVIEW]}
    inputs=[]
    for split,sources in [('template_self_test',templates),('source_withheld',reserved)]:
        directory=output/split;directory.mkdir(exist_ok=True);fasta=directory/'proteins.faa';proteins=[]
        for source in sources:
            entities={l['entity_id'] for s in source['sites'] for l in s['ligand_residues']}
            for entity in sorted(entities):
                sequence=source['sequences'][entity];pid=source['pdb_id']+'__entity_'+entity
                proteins.append(SeqRecord(Seq(sequence),id=pid,description='source_withheld' if split=='source_withheld' else 'template_self_test'))
                inputs.append(dict(protein_id=pid,split=split,pdb_id=source['pdb_id'],source=source['source'],
                    sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),source_sha256=source['sha256']))
        SeqIO.write(proteins,fasta,'fasta')
        predictions,_=annotate_fes(fasta,directory)
        rows,metrics[split]=compare(predictions,sources)
        all_rows.extend(dict(split=split,**r) for r in rows)
        per_source={}
        for source in sources:
            _,per_source[source['pdb_id']]=compare([p for p in predictions if p['protein_id'].startswith(source['pdb_id']+'__')],[source])
        metrics[split]['per_source']=per_source
    write_tsv(output/'site_comparisons.tsv',all_rows)
    (output/'input_manifest.json').write_text(json.dumps(inputs,indent=2)+'\n')
    corpus=ROOT/policy['donor_gain_validation_corpus']
    freeze=json.loads((corpus/'freeze.json').read_text())
    for name,digest in freeze['artifacts'].items():
        if sha(corpus/name)!=digest:raise ValueError('FeS donor-gain fixture changed: '+name)
    for name,digest in freeze['source_files'].items():
        if sha(ROOT/name)!=digest:raise ValueError('FeS donor-gain fixture source changed: '+name)
    directory=output/'donor_gain_guard';directory.mkdir()
    predictions,_=annotate_fes(corpus/'proteins.faa',directory)
    cases=json.loads((corpus/'cases.json').read_text());guard_cases=[]
    for case in cases['cases']:
        matched=[p for p in predictions if p['protein_id']==case['protein_id'] and p['template_site_id']==cases['template_site_id']]
        if len(matched)!=1:raise ValueError('Expected one reviewed donor-gain site: '+case['protein_id'])
        prediction=matched[0];alternative=prediction['alternative_coordination_states'][0]
        passed=(prediction['status']==case['expected_status'] and alternative['query_position']==case['expected_alternative_position'] and
            (not case['expected_reason'] or case['expected_reason'] in prediction['reasons']) and
            prediction['confidence']=='provisional_structural_homology_not_calibrated' and prediction['cofactor_occupancy']=='unknown')
        guard_cases.append(dict(**case,passed=passed,prediction=prediction))
    (directory/'checks.json').write_text(json.dumps(guard_cases,indent=2)+'\n')
    held=metrics['source_withheld'];guard=policy['validation_guard'];self_test=metrics['template_self_test']
    target_met=bool(held['precision'] is not None and held['precision']>=guard['precision_target'] and
        held['exact_site_recall']>=guard['min_evaluable_site_coverage'])
    unchanged=all(sha(ROOT/name)==digest for name,digest in before.items())
    reviewed_sites=sum(len(source['sites']) for source in templates)
    regression_ok=self_test['TP']==reviewed_sites and self_test['FP']==0 and self_test['FN']==0 and unchanged and all(c['passed'] for c in guard_cases)
    summary=dict(status='execution_and_self_test_passed' if regression_ok else 'failed',
        before_hashes=before,unchanged=unchanged,metrics=metrics,development_numeric_target_met=target_met,
        donor_gain_guard_cases=len(guard_cases),donor_gain_guard_passed=sum(c['passed'] for c in guard_cases),
        complete_FES_gate='not_met',
        limits=['Source-withheld structures have been inspected during development, not a sealed final blind test.',
            'No independent Huc topology query or complete reviewed alternative-ligand coverage yet.',
            'Point precision and a few structures cannot establish broad calibrated high-confidence calls.',
            'Unreviewed extra/fewer observed contact residues and predictions outside modeled donor atoms remain distinct ungradable cases; query occupancy remains unknown.'])
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not regression_ok:raise SystemExit(1)


if __name__=='__main__':main()
