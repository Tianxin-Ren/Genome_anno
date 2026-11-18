#!/usr/bin/env python3
"""Actual searches on native HoxF and rearranged reference-CDS controls; not blind accuracy."""
import argparse
import copy
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from annotate import ROOT,annotate,extract_genbank,write_tsv


def table(path):
    with path.open() as handle:return list(csv.DictReader(handle,delimiter='\t'))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():parser.error('Use a fresh result directory')
    refs=json.loads((ROOT/'references/processed/references.json').read_text());systems=json.loads((ROOT/'config/systems.json').read_text())
    policy=ROOT/'config/module_graph.json';before=hashlib.sha256(policy.read_bytes()).hexdigest();checks=[]
    def genes_for(sid):
        s=next(s for s in systems if s['structure_id']==sid);ids={r['locus_tag'] for r in refs if r['structure_id']==sid}
        return copy.deepcopy([g for g in extract_genbank(ROOT/'references/raw'/(s['sample']+'.gb'))[0] if g['id'] in ids])
    def check(name,passed,detail):checks.append(dict(test=name,passed=bool(passed),observed=json.dumps(detail)))
    def run(name,genes,mode='genome'):
        for g in genes:g.update(gene='',aliases='',protein_accession='')
        result=annotate(genes,out/name,mode=mode,threads=1)
        nodes=table(out/name/'module_graph_nodes.tsv');edges=table(out/name/'module_graph_edges.tsv')
        resources=defaultdict(set);duplicate=[]
        for e in edges:
            key=(e['anchor'],e['protein_id']);positions=set(json.loads(e['query_residue_positions']))
            if positions&resources[key]:duplicate.append(key)
            resources[key]|=positions
        check(name+':no_reused_aligned_residues',not duplicate,duplicate)
        check(name+':localization_unknown',all(n['localization_status'].startswith('not_predicted') for n in nodes),len(nodes))
        return result,nodes,edges
    sid='NiFe_HoxFUYH';genes=genes_for(sid)
    result,nodes,edges=run('native_HoxF',genes)
    fids={n['node_id']:n for n in nodes if n['role']=='HoxF'}
    check('native_HoxF:two_supported_modules_on_one_protein',set(fids)=={'HoxF_NuoE_like','HoxF_NuoF_like'} and
        len({n['protein_ids'] for n in fids.values()})==1 and all(n['status']=='homology_supported' for n in fids.values()),fids)
    check('native_HoxF:original_function_retained',any(a['structure_id']==sid and a['potential_function_ids']=='Ni14' for a in result),[a.get('potential_function_ids') for a in result])
    gid=next(r['locus_tag'] for r in refs if r['reference_id']==sid+'__HoxF')
    split=copy.deepcopy(genes);f=next(g for g in split if g['id']==gid);split.remove(f)
    for index,(lo,hi) in enumerate([(0,160),(160,len(f['sequence']))]):
        split.append(dict(f,id=gid+'_synthetic_part_'+str(index),sequence=f['sequence'][lo:hi],start=f['start']+3*lo,end=f['start']+3*hi-1,partial=True))
    result,nodes,edges=run('split_HoxF_CDS',split)
    fids=[n for n in nodes if n['role']=='HoxF']
    check('split_HoxF:two_modules_on_two_CDS',len(fids)==2 and len({n['protein_ids'] for n in fids})==2 and all(n['status']=='homology_supported' for n in fids),fids)
    check('split_HoxF:no_function_or_high_structure',all(not a.get('potential_function_ids') and a['structure_confidence']!='high' for a in result),result)
    sid='FeFe_TmHydABC';genes=genes_for(sid);by_role={r['role']:r['locus_tag'] for r in refs if r['structure_id']==sid}
    b=next(g for g in genes if g['id']==by_role['HydB']);c=next(g for g in genes if g['id']==by_role['HydC'])
    fused=[g for g in copy.deepcopy(genes) if g['id'] not in [b['id'],c['id']]]
    fused.append(dict(b,id='synthetic_HydBC',sequence=b['sequence']+c['sequence'],start=min(b['start'],c['start']),end=max(b['end'],c['end'])))
    result,nodes,edges=run('fused_HydBC',fused)
    bc=[n for n in nodes if n['role'] in ['HydB','HydC']]
    check('fused_HydBC:both_roles_on_nonoverlapping_resources',len(bc)==2 and all(n['protein_ids']=='synthetic_HydBC' and n['status']=='homology_supported' for n in bc),bc)
    check('fused_HydBC:organization_does_not_prove_function',all(not a.get('potential_function_ids') and a['structure_confidence']!='high' for a in result),result)
    fragments=[g for g in copy.deepcopy(genes) if g['id']!=b['id']];cut=len(b['sequence'])//2
    for index,(lo,hi) in enumerate([(0,cut),(cut,len(b['sequence']))]):
        fragments.append(dict(b,id='synthetic_HydB_part_'+str(index),sequence=b['sequence'][lo:hi],start=b['start']+3*lo,end=b['start']+3*hi-1,partial=True))
    result,nodes,edges=run('fragmented_HydB',fragments)
    bn=next(n for n in nodes if n['role']=='HydB')
    check('fragmented_HydB:collinear_union_supported',bn['status']=='homology_supported' and bn['realization']=='collinear_CDS_fragments',bn)
    check('fragmented_HydB:no_function_or_high_structure',all(not a.get('potential_function_ids') and a['structure_confidence']!='high' for a in result),result)
    duplicated=copy.deepcopy(genes);duplicated.append(dict(c,id='synthetic_duplicate_HydC',start=c['start']-400,end=c['end']-400))
    result,_,_=run('duplicate_HydC',duplicated)
    check('duplicate_HydC:ambiguous_pairing_withheld',any('ambiguous_modules' in a['flags'] for a in result) and all(not a.get('potential_function_ids') for a in result),result)
    huc=[g for g in genes_for('NiFe_HucLMS') if g['id']!='MSMEG_RS10955']
    result,nodes,_=run('Huc_without_assembly_accessory',huc)
    check('Huc_without_M:core_separate_from_oligomer_and_function',any(a['structure_status']=='core_complete_accessory_missing' and a['missing_optional_roles']=='HucM' for a in result) and all(not a.get('potential_function_ids') for a in result),result)
    check('Huc_without_M:optional_in_structure_required_for_reference_function',any(n['role']=='HucM' and n['required_for_structure']=='False' and n['required_for_function']=='Ni09' for n in nodes),nodes)
    write_tsv(out/'checks.tsv',checks)
    summary=dict(status='development_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),passed=sum(c['passed'] for c in checks),
        failed=[c for c in checks if not c['passed']],policy_sha256=before,unchanged=before==hashlib.sha256(policy.read_bytes()).hexdigest(),
        scope='Native reference HoxF domain organization plus five reference-CDS rearrangement/ambiguity/absence controls. Not independent complex precision/recall, natural fission proof or validated novel fusion activity.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    if summary['status']!='development_checks_passed' or not summary['unchanged']:raise SystemExit(1)


if __name__=='__main__':main()
