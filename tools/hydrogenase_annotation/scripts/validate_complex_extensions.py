#!/usr/bin/env python3
"""Recheck every post-identity reference against native CDSs and scoped reaction metadata."""
import argparse
import json
from pathlib import Path
from annotate import ROOT,extract_genbank,write_tsv
from catalytic_confirmation import sha
from complex_extensions import load_bank
from build_complex_extensions import source_case


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    refs,systems,graphs,manifest,contracts=load_bank();checks=[]
    review=json.loads((ROOT/'config/complex_extension_sources.json').read_text());corpus=ROOT/review['source_corpus']
    sources=json.loads((corpus/'sources.json').read_text());truth={c['id']:c for c in json.loads((corpus/'truth.json').read_text())}
    def check(name,passed,observed):checks.append(dict(test=name,passed=bool(passed),observed=json.dumps(observed)))
    for spec in review['systems']:
        c,source=source_case(spec,truth,sources);genes={g['id']:g for g in extract_genbank(source)[0]}
        for r in [r for r in refs if r['structure_id']==spec['structure_id']]:
            member=next(m for m in c['members'] if m['role']==r['role']);g=genes[r['locus_tag']]
            check(r['reference_id']+':native_translation_and_experimental_member',r['sequence']==g['sequence'] and r['sequence_sha256']==member['sequence_sha256'] and
                r['protein_accession']==member['protein_accession'] and r['kind']==('catalytic' if r['role']==spec['catalytic_role'] else 'module'),
                dict(protein=r['protein_accession'],locus=r['locus_tag'],role=r['role'],kind=r['kind'],source_sha256=sha(source)))
    a=next(r for r in refs if r['role']=='HydA_ET_not_H_cluster')
    check('HydA_electron_transfer_homolog_is_not_a_second_catalyst',a['kind']=='module',a['reference_id'])
    c={c['contract_id']:c for c in contracts}
    check('NAD_and_NADP_bifurcation_not_collapsed',c['Fe04']['electron_acceptors']==['NAD+','ferredoxin_oxidized'] and c['Fe07']['electron_acceptors']==['NADP+','ferredoxin_oxidized'] and c['Fe04']['reaction_key']!=c['Fe07']['reaction_key'],[c['Fe04']['reaction_key'],c['Fe07']['reaction_key']])
    check('NADH_only_and_NADH_Fd_reactions_distinct',c['Fe03']['electron_donors']==['NADH','ferredoxin_reduced'] and c['Fe05']['electron_donors']==['NADH'] and c['Fe05']['ferredoxin_requirement']=='not_required_in_reference_assay', [c['Fe03']['electron_donors'],c['Fe05']['electron_donors']])
    check('Hox_activation_not_NADPH_coupling',c['Ni14']['electron_acceptors']==['NAD+'] and c['Ni12']['electron_donors']==['NADPH'],[c['Ni14']['reaction_key'],c['Ni12']['reaction_key']])
    check('physiology_and_regulation_do_not_invent_direct_acceptors',all(not c[x]['electron_acceptors'] and c[x]['direct_acceptor_status']=='unresolved' for x in ['Ni06','Ni08']) and c['Ni10']['capability_scope']=='regulatory_pathway_potential',[c[x]['capability_scope'] for x in ['Ni06','Ni08','Ni10']])
    check('new_structures_enable_no_legacy_function_ids',all(not s['function_ids'] for s in systems),[s['structure_id'] for s in systems])
    write_tsv(out/'checks.tsv',checks)
    result=dict(status='source_review_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),passed=sum(c['passed'] for c in checks),failed=[c for c in checks if not c['passed']],
        model_manifest_sha256=sha(ROOT/'references/processed/complex_extensions/manifest.json'),reaction_registry_sha256=sha(ROOT/'config/reaction_contracts.json'),
        additional_structures=len(systems),additional_proteins=len(refs),typed_contracts=len(contracts),reference_replay_cases=manifest['reference_replay_cases'],
        scope='Source and reference reaction curation only. Development templates are reference replays, not independent accuracy. No new L1 or FeS model; query potential requires the separate exact whole-core rule.')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','checks','passed','additional_structures','additional_proteins','typed_contracts']},indent=2))
    if result['failed']:raise SystemExit(1)


if __name__=='__main__':main()
