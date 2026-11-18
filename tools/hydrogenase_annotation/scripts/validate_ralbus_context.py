#!/usr/bin/env python3
"""New native-context development challenge with purified and lower-evidence cases separated."""
import argparse
import copy
import hashlib
import json
from datetime import datetime,timezone
from pathlib import Path
from annotate import ROOT,annotate,extract_genbank,write_tsv
from catalytic_confirmation import sha,table
from complex_extensions import load_contracts
from prepare_natural_complexes import anonymous_genes
from run_acceptance_checks import state_hashes


def verify_sources(policy):
    for name,digest in policy['source_files'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Native challenge source changed: '+name)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    protocol=ROOT/'config/round24_ralbus_validation.json';policy=json.loads(protocol.read_text());verify_sources(policy)
    correspondence=json.loads((ROOT/policy['correspondence']).read_text());native=[]
    for name in policy['genomes']:native.extend(extract_genbank(ROOT/name)[0])
    by_native={g['id']:g for g in native}
    for r in correspondence['members']:
        g=by_native[r['locus']]
        if hashlib.sha256(g['sequence'].encode()).hexdigest()!=r['sequence_sha256'] or g['protein_accession']!=r['protein_accession']:
            raise ValueError('Figure/nativeCDS correspondence mismatch: '+r['locus'])
    genes,mapping=anonymous_genes(native)
    freeze=dict(status='development_source_and_assay_truth_frozen_before_prediction',created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(protocol),source_files=policy['source_files'],source_to_query=mapping,truth=policy,
        protein_count=len(genes),native_replicons=2,independence='not_claimed',
        input_anonymization='Gene,alias,protein_accession cleared; coordinates,partial flags,sequence and replicon ownership retained.')
    (out/'freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
    before=state_hashes();run=out/'native_context';run.mkdir()
    annotations=annotate(copy.deepcopy(genes),run,threads=2,mode='genome',input_files=[ROOT/n for n in policy['genomes']])
    by_query={a['protein_id']:a for a in annotations};c=policy['scored_complex'];qid=mapping[c['anchor']];a=by_query.get(qid,{})
    complete=a.get('structure_status') in ['complete_supported','core_complete_accessory_missing','composition_supported_reorganized']
    predicted={m['protein_id'] for m in table(run/'complex_members.tsv') if m['anchor']==qid and m['required_for_structure']=='True'}
    expected={mapping[n] for n in c['members']};exact=complete and predicted==expected
    contracts={c['contract_id']:c for c in load_contracts()};assigned={contracts[f]['reaction_key'] for f in a.get('potential_function_ids','').split(';') if f}
    positive=set(c['supported_reaction_keys']);negative=set(c['assayed_negative_reaction_keys'])
    associated=by_query.get(mapping[policy['associated_protein']['locus']],{})
    indeterminate=by_query.get(mapping[policy['indeterminate_protein']['locus']],{})
    checks=[dict(test='anonymous_prediction_input',passed=all(not g['gene'] and not g['aliases'] and not g['protein_accession'] for g in genes)),
        dict(test='indeterminate_HydS_physiology_remains_unresolved',passed=not indeterminate or (indeterminate.get('in_vivo_physiology')=='unresolved' and indeterminate.get('net_direction')=='unresolved'))]
    after=state_hashes();verify_sources(policy)
    checks.append(dict(test='predictor_and_sources_unchanged',passed=before==after))
    summary=dict(status='development_native_assay_challenge_evaluated' if all(c['passed'] for c in checks) else 'failed',
        before_hashes=before,after_hashes=after,protocol_sha256=sha(protocol),source_freeze_sha256=sha(out/'freeze.json'),
        proteins=len(genes),replicons=2,complete_strain_genome=False,independent_validation=False,
        core_metrics=dict(TP=int(exact),FP=int(complete and not exact),FN=int(not exact)),
        reaction_metrics=dict(TP=len(assigned&positive),FP=len(assigned&negative),FN=len(positive-assigned),
            ungraded_extra_reactions=sorted(assigned-positive-negative),counting_scope='Two assay directions from one purified biological complex; no independent-replicate claim.'),
        expected_members=sorted(expected),predicted_members=sorted(predicted),HydABC=a or 'not_retained',
        HydA2_activity_association=dict(annotation=associated or 'not_retained',score_scope=policy['associated_protein']['metric_scope']),
        HydS_indeterminate=dict(annotation=indeterminate or 'not_retained',score_scope=policy['indeterminate_protein']['metric_scope']),
        ungraded_anchors=[a['protein_id'] for a in annotations if a['protein_id'] not in {qid,mapping[policy['associated_protein']['locus']],mapping[policy['indeterminate_protein']['locus']]}],
        checks=checks,limits=policy['limits'],run_manifest_sha256=sha(run/'run_manifest.json'))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');write_tsv(out/'checks.tsv',checks)
    print(json.dumps({k:summary[k] for k in ['status','proteins','core_metrics','reaction_metrics','independent_validation']},indent=2))
    if summary['status']=='failed':raise SystemExit(1)


if __name__=='__main__':main()
