#!/usr/bin/env python3
"""Natural split-center source cases and actual sequence-search boundary controls."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
from annotate import ROOT,annotate,write_tsv
from active_centers import split_structural_record
from split_centers import source_genes


def table(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():parser.error('Use a fresh result directory')
    policy_path=ROOT/'config/split_centers.json';policy=json.loads(policy_path.read_text());before=hashlib.sha256(policy_path.read_bytes()).hexdigest();checks=[]
    for relative,pin in policy['source_files'].items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=pin['sha256']:raise ValueError('Development source changed')
    truth=split_structural_record(ROOT/policy['template']['path'],['14','6'],'HF','420')
    def check(name,passed,detail):checks.append(dict(test=name,passed=bool(passed),observed=json.dumps(detail)))
    check('same_NiFe_pair_four_ligands_on_two_polypeptides',[c['ligands'] for c in truth['components']]==[{'61':'C','64':'C'},{'20':'U','23':'C'}],truth)
    def run(name,genes,anchor,mode='genome'):
        genes=copy.deepcopy(genes)
        for g in genes:g.update(gene='',aliases='',protein_accession='')
        directory=out/name;directory.mkdir();input_path=directory/'input_genes.json';input_path.write_text(json.dumps(genes,indent=2)+'\n')
        results=annotate(genes,directory,mode=mode,threads=1,input_files=[input_path]);a=next(r for r in results if r['protein_id']==anchor)
        sites=[s for s in table(directory/'catalytic_sites.tsv') if s['catalytic_anchor']==anchor and s['residue_role']=='direct_metal_ligand']
        candidates=table(directory/'catalytic_center_candidates.tsv')
        check(name+':no_function_promotion',not a['potential_function_ids'] and a['catalytic_status']=='candidate',a)
        return a,sites,candidates
    mapping=[]
    for case in policy['development_cases']:
        genes=source_genes(ROOT,case);a,sites,_=run(case['case_id'],genes,case['anchor_id'])
        expected={(case['anchor_id'],61,'C'),(case['anchor_id'],64,'C'),(case['partner_id'],20,'U'),(case['partner_id'],23,'C')}
        observed={(s['protein_id'],int(s['query_position']),s['observed_residue']) for s in sites if s['status']=='present'}
        check(case['case_id']+':complete_two_chain_center',a['active_center_status']=='multi_polypeptide_ligand_set_supported' and observed==expected,observed and sorted(observed))
        mapping.append(dict(case_id=case['case_id'],template_source=case['template_source'],expected=sorted(expected),observed=sorted(observed),
            source_scope='Template replay' if case['template_source'] else 'Same-genus source-reviewed development mapping; not independent structural or blind final accuracy'))
    case=policy['development_cases'][1];native=source_genes(ROOT,case);raw=source_genes(ROOT,case,False)
    a,_,candidates=run('historical_X_not_automatically_Sec',raw,case['anchor_id'])
    raw_sites=[s for c in candidates for s in json.loads(c['sites']) if s['center_component_id']=='VhuU' and s['reference_position']==20]
    check('historical_X:ambiguous_not_substituted_or_corrected',a['active_center_status']!='multi_polypeptide_ligand_set_supported' and any(s['status']=='ambiguous' for s in raw_sites),raw_sites)
    for name,change,expected_state in [('U_to_A',lambda s:s[:19]+'A'+s[20:],'substituted'),('U_to_X',lambda s:s[:19]+'X'+s[20:],'ambiguous'),
        ('TGA_truncated_small_CDS',lambda s:s[:19],'not_covered')]:
        genes=copy.deepcopy(native);small=next(g for g in genes if g['id']==case['partner_id']);small['sequence']=change(small['sequence']);small['coordinate_evidence']='synthetic_reference_sequence_perturbation'
        a,_,candidates=run(name,genes,case['anchor_id']);ss=[s for c in candidates for s in json.loads(c['sites']) if s['center_component_id']=='VhuU' and s['reference_position']==20]
        check(name+':residue_state',a['active_center_status']!='multi_polypeptide_ligand_set_supported' and any(s['status']==expected_state for s in ss),ss)
    genes=[g for g in copy.deepcopy(native) if g['id']!=case['partner_id']];a,_,candidates=run('missing_small_component',genes,case['anchor_id'])
    check('missing_component:not_covered_not_mutated',a['split_center_status']=='required_component_not_detected' and all(s['status']=='not_covered' for c in candidates for s in json.loads(c['sites']) if s['center_component_id']=='VhuU'),candidates)
    genes=copy.deepcopy(native);small=next(g for g in genes if g['id']==case['partner_id']);small['sequence']=small['sequence'][1:26];small['coordinate_evidence']='source_reviewed_processed_peptide_on_original_CDS'
    a,sites,_=run('reviewed_mature_25aa_peptide',genes,case['anchor_id'])
    check('mature_peptide:ligand_coordinates_not_precursor_tail',{(int(s['query_position']),s['observed_residue']) for s in sites if s['protein_id']==case['partner_id']}=={(19,'U'),(22,'C')} and a['active_center_status']=='multi_polypeptide_ligand_set_supported',sites)
    genes=copy.deepcopy(native);large=next(g for g in genes if g['id']==case['anchor_id']);small=next(g for g in genes if g['id']==case['partner_id']);large['sequence']+=small['sequence'];genes.remove(small)
    large['coordinate_evidence']='synthetic_protein_concatenation_not_the_published_fusion_construct'
    a,sites,_=run('synthetic_AU_fusion',genes,case['anchor_id'])
    check('fused_sequence:four_sites_on_one_protein',{(int(s['query_position']),s['observed_residue']) for s in sites if s['status']=='present'}=={(61,'C'),(64,'C'),(440,'U'),(443,'C')},sites)
    for name,mode in [('different_bin','bin'),('different_contig','assembly'),('unlinked_transcripts','transcript')]:
        genes=copy.deepcopy(native)
        for g in genes:g.update(bin_id='one')
        small=next(g for g in genes if g['id']==case['partner_id'])
        if mode!='transcript':small.update(bin_id='two' if mode=='bin' else 'one',contig='synthetic_other_contig')
        a,_,_=run(name,genes,case['anchor_id'],mode)
        check(name+':no_cross_ownership_center',a['active_center_status']!='multi_polypeptide_ligand_set_supported' and a['split_center_status']=='required_component_not_detected',a)
    genes=source_genes(ROOT,policy['development_cases'][0]);small=next(g for g in genes if g['id']==policy['development_cases'][0]['partner_id'])
    genes.append(dict(small,id=small['id']+'_tail_variant',sequence=small['sequence'][:27]+'G'*17,start=small['start']+150,end=small['end']+150,
        coordinate_evidence='synthetic_precursor_tail_variant_with_same_structural_window'))
    a,_,candidates=run('two_mature_windows_different_precursor_tails',genes,policy['development_cases'][0]['anchor_id'])
    check('precursor_tail:does_not_choose_the_binding_partner',a['split_center_status']=='ambiguous_component_pairing' and a['active_center_status']!='multi_polypeptide_ligand_set_supported' and len({float(c['partner_window_score']) for c in candidates})==1,candidates)
    write_tsv(out/'checks.tsv',checks)
    summary=dict(status='development_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),passed=sum(c['passed'] for c in checks),failed=[c for c in checks if not c['passed']],
        policy_sha256=before,unchanged=before==hashlib.sha256(policy_path.read_bytes()).hexdigest(),natural_source_mappings=mapping,
        scope='One structural template center, two source-reviewed natural operon cases of the same genus, plus synthetic input/state/processing/ownership controls. No blind independent accuracy or complete-complex/function precision claimed.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='natural_source_mappings'},indent=2))
    if summary['status']!='development_checks_passed' or not summary['unchanged']:raise SystemExit(1)


if __name__=='__main__':main()
