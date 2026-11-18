#!/usr/bin/env python3
"""Actual searches on shared-partner reference-CDS controls; not independent interaction truth."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
from annotate import ROOT,annotate,extract_genbank,write_tsv


def table(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():parser.error('Use a fresh result directory')
    refs=json.loads((ROOT/'references/processed/references.json').read_text());systems=json.loads((ROOT/'config/systems.json').read_text())
    policy=ROOT/'config/module_graph.json';before=hashlib.sha256(policy.read_bytes()).hexdigest();checks=[]
    def check(name,passed,observed):checks.append(dict(test=name,passed=bool(passed),observed=json.dumps(observed)))
    def duplicate(sid,partial=False):
        s=next(s for s in systems if s['structure_id']==sid);ids={r['locus_tag'] for r in refs if r['structure_id']==sid}
        genes=copy.deepcopy([g for g in extract_genbank(ROOT/'references/raw'/(s['sample']+'.gb'))[0] if g['id'] in ids])
        local={r['locus_tag'] for r in refs if r['structure_id']==sid and r['role'] in s['local']}
        catalytic=next(r['locus_tag'] for r in refs if r['structure_id']==sid and r['role']==s['catalytic_role'])
        for g in list(genes):
            if g['id'] in local:
                if partial and g['id']==next(r['locus_tag'] for r in refs if r['reference_id']==sid+'__MvhG'):continue
                genes.append(dict(g,id=g['id']+'_copy',start=g['start']+120000,end=g['end']+120000))
        length=max(g['end'] for g in genes)+1000
        for g in genes:g.update(gene='',aliases='',protein_accession='',contig_length=length,coordinate_evidence='synthetic_reference_CDS_rearrangement')
        return genes,catalytic
    def run(name,sid,genes,mode='genome'):
        result=annotate(genes,out/name,mode=mode,threads=1)
        anchors=[a for a in result if a['structure_id']==sid];pairs=table(out/name/'cross_anchor_associations.tsv')
        check(name+':interaction_not_measured',all(a['physical_interaction_status']=='not_established_by_sequence' for a in anchors),len(anchors))
        return anchors,pairs
    sid='NiFe_MvhHdr';genes,cat=duplicate(sid)
    anchors,pairs=run('two_Mvh_cores_shared_Hdr',sid,genes)
    check('shared_Hdr:two_catalytic_cores_retained',len(anchors)==2 and all(a['catalytic_status']=='supported' for a in anchors),anchors)
    check('shared_Hdr:specific_pairings_withheld',all(not a['potential_function_ids'] and a['structure_confidence']!='high' and a['cross_anchor_pairing_status']=='overlapping_fixed_member_pairing_unresolved' for a in anchors),anchors)
    check('shared_Hdr:three_selected_member_overlaps',len(pairs)==3 and all(p['blocks_specific_pairing']=='True' for p in pairs),pairs)
    genes,cat=duplicate(sid,partial=True);anchors,pairs=run('partial_competing_Mvh_core',sid,genes)
    original=next(a for a in anchors if a['protein_id']==cat)
    check('partial_competitor:complete_core_not_erased',original['potential_function_ids']=='Ni13' and original['structure_confidence']=='high',original)
    check('partial_competitor:overlap_still_reported',bool(pairs) and all(p['pairing_status']=='provisional_member_overlap_unresolved' and p['blocks_specific_pairing']=='False' for p in pairs),pairs)
    genes,_=duplicate('FeFe_TmHydABC');anchors,pairs=run('two_separate_HydABC_clusters','FeFe_TmHydABC',genes)
    check('separate_clusters:both_reference_potentials_retained',len(anchors)==2 and all(set(a['potential_function_ids'].split(';'))=={'Fe03','Fe04'} and a['structure_confidence']=='high' for a in anchors),anchors)
    check('separate_clusters:no_shared_gene_claim',not pairs,pairs)
    genes,_=duplicate('NiFe_periplasmic_Hys');anchors,pairs=run('two_Hys_cores_shared_CycA','NiFe_periplasmic_Hys',genes)
    check('diffusible_partner:nonexclusive_candidates',len(anchors)==2 and all(a['potential_function_ids']=='Ni01' for a in anchors) and bool(pairs) and all(p['pairing_status']=='shared_diffusible_pool_candidate' for p in pairs),pairs)
    genes,_=duplicate('NiFe_HupSL');anchors,pairs=run('two_Hup_cores_shared_Nif_context','NiFe_HupSL',genes)
    check('physiological_context:not_physical_allocation',len(anchors)==2 and all(a['potential_function_ids']=='Ni08' for a in anchors) and bool(pairs) and all(p['pairing_status']=='shared_context_cooccurrence' for p in pairs),pairs)
    genes,cat=duplicate(sid)
    for g in genes:g.update(bin_id='two' if g['id'].endswith('_copy') else 'one',contig='synthetic_bin_two' if g['id'].endswith('_copy') else 'synthetic_bin_one')
    anchors,pairs=run('cross_bin_shared_partner_forbidden',sid,genes,mode='bin')
    check('cross_bin:no_combined_pair',not pairs and next(a for a in anchors if a['protein_id']==cat)['potential_function_ids']=='Ni13' and not next(a for a in anchors if a['protein_id']!=cat)['potential_function_ids'],anchors)
    genes,_=duplicate(sid);hdr_id=next(r['locus_tag'] for r in refs if r['reference_id']==sid+'__HdrA');hdr=next(g for g in genes if g['id']==hdr_id)
    sequence=list(hdr['sequence']);position=next(i for i in range(30,len(sequence)) if sequence[i] not in 'CU');sequence[position]='A' if sequence[position]!='A' else 'G'
    genes.append(dict(hdr,id=hdr_id+'_near_alternative',sequence=''.join(sequence)))
    anchors,pairs=run('unselected_competitive_Hdr_alternative',sid,genes)
    check('competitive_alternative:unselected_claim_included',any(p['selection_a']=='competitive_alternative' or p['selection_b']=='competitive_alternative' for p in pairs),pairs)
    check('competitive_alternative:no_specific_function',all(not a['potential_function_ids'] for a in anchors),anchors)
    write_tsv(out/'checks.tsv',checks)
    summary=dict(status='development_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),passed=sum(c['passed'] for c in checks),
        failed=[c for c in checks if not c['passed']],policy_sha256=before,unchanged=before==hashlib.sha256(policy.read_bytes()).hexdigest(),
        scope='Seven actual annotation runs using synthetic copies/rearrangements of real reference CDSs. Distinguishes fixed-member ambiguity, diffusible/background sharing, partial competitors, alternative realizations and bin ownership. Not independent natural interaction precision, gene-copy stoichiometry or measured activity.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    if summary['status']!='development_checks_passed' or not summary['unchanged']:raise SystemExit(1)


if __name__=='__main__':main()
