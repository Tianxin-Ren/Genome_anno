#!/usr/bin/env python3
"""Evaluate declared known-reference topology controls; never claim an independent test."""
import argparse
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT,write_tsv
from localization import sha,table,verify_completed,parse_topology
from validate_fes_prediction import wilson


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--localization',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    source=args.localization.resolve();run=verify_completed(source);out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    protocol=ROOT/'config/localization_validation.json';policy=json.loads(protocol.read_text());rows=table(source/'localization.tsv')
    queries={r.id:str(r.seq) for r in SeqIO.parse(source/'queries.faa','fasta')};predictions=parse_topology(source/'raw.pred',queries)
    dataset=ROOT/'tools/TMbed/data/datasets/all_seq.fasta'
    training={str(r.seq).upper() for r in SeqIO.parse(dataset,'fasta')}
    checks=[];metrics={k:dict(TP=0,FP=0,FN=0,TN=0,ungradable=0) for k in ['TM','SP']};source_rows=[]
    for case in policy['cases']:
        owners=[r for r in rows if r['protein_id']==case['locus_tag'] and Path(r['source_run']).name==case['sample']]
        if len(owners)!=1:raise ValueError('Missing/repeated declared reference owner: '+case['reference_id'])
        row=owners[0];sequence=queries[row['query_id']]
        if hashlib.sha256(sequence.encode()).hexdigest()!=case['sequence_sha256']:raise ValueError('Reference sequence correspondence changed: '+case['reference_id'])
        source_rows.append(dict(reference_id=case['reference_id'],query_id=row['query_id'],exact_TMbed_development_sequence_overlap=sequence in training))
        for prop,key in [('TM','transmembrane_segments'),('SP','signal_peptide_predicted')]:
            expected=case['expected_'+prop]
            if expected is None:continue
            observed=int(row[key])>0 if prop=='TM' else row[key]=='True'
            outcome=('TP' if observed else 'FN') if expected else ('FP' if observed else 'TN');metrics[prop][outcome]+=1
            checks.append(dict(reference_id=case['reference_id'],protein_id=row['protein_id'],source_run=row['source_run'],property=prop,
                expected=expected,observed=observed,outcome=outcome,evidence_url=case['evidence_url'],
                training_exact_sequence_overlap=sequence in training,evidence_scope='Known-reference control, not independent accuracy or boundary validation'))
    for m in metrics.values():
        positive=m['TP']+m['FN'];predicted=m['TP']+m['FP']
        m.update(recall=m['TP']/positive if positive else None,precision=m['TP']/predicted if predicted else None,
            recall_wilson95=wilson(m['TP'],positive),precision_wilson95=wilson(m['TP'],predicted))
    write_tsv(out/'reference_checks.tsv',checks);write_tsv(out/'training_exact_overlap.tsv',source_rows)
    summary=dict(status='reference_controls_evaluated_not_independent',declared_proteins=len(policy['cases']),property_checks=len(checks),
        correct=sum(c['outcome'] in ['TP','TN'] for c in checks),disagreements=[c for c in checks if c['outcome'] in ['FP','FN']],metrics=metrics,
        exact_TMbed_development_sequence_overlaps=sum(r['exact_TMbed_development_sequence_overlap'] for r in source_rows),
        generic_training_sequence_count=len(training),protocol_sha256=sha(protocol),localization_manifest_sha256=sha(source/'localization_manifest.json'),
        official_training_dataset_sha256=sha(dataset),source_auxiliary_run=run,
        limits=['Reference controls are exposed and may share generic supervised training or close homologs; exact-overlap counts do not establish homology/genus separation.',
            'No independent hydrogenase-specific calibrated accuracy, exact membrane boundary or cellular-compartment truth is established.',
            'No outputs, thresholds, weights or original hydrogenase labels changed after observing these controls.'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['status','declared_proteins','property_checks','correct','disagreements','metrics','exact_TMbed_development_sequence_overlaps']},indent=2))


if __name__=='__main__':main()
