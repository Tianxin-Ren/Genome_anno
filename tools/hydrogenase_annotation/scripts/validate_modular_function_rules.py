#!/usr/bin/env python3
"""Exercise divergent, fused functional modules with fresh BLAST and HMM evidence."""
import argparse
import copy
import json
import sys
from pathlib import Path

from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'research_workflow'))
from function_rules import annotate,load_rules,sha,table,write
from query_domain_organization import annotate_domains


def validate(source,out,split_source=None):
    source,out=Path(source).resolve(),Path(out).resolve()
    if out.exists():raise FileExistsError('Use a fresh validation output')
    out.mkdir(parents=True);annotation=out/'annotation';annotation.mkdir()
    original_genes={r['id']:r for r in table(source/'annotation/genes.tsv')}
    original_sequences={r.id:str(r.seq) for r in SeqIO.parse(source/'annotation/proteins.faa','fasta')}
    by_number={pid.rsplit('_',1)[-1]:pid for pid in original_genes}
    anchors={r['candidate_id']:r for r in table(source/'functions/functional_annotations.tsv')}
    model,_,_=load_rules()
    donors=model['references'][model['HdrA_donor_comparison']['F420_Fd_alternative']]['sequence']
    numbers=['1136','1137','1163','1785','1786','1995','1996','1997','1998','2001']
    cases={
        'renamed_baseline':(True,True),
        'circular_rotation':(True,True),
        'delete_HdrC':(False,True),
        'destroy_HdrC_ligands':(False,True),
        'replace_HdrC_with_NuoI':(False,True),
        'delete_MvhD_fused_region':(False,True),
        'split_remote_HdrA_MvhD':(True,True),
        'unlinked_remote_HdrA_MvhD':(False,True),
        'opposite_strand_remote_HdrA_MvhD':(False,True),
        'cross_bin_remote_MvhD':(False,True),
        'replace_HdrA_with_F420_donor':(False,True),
        'cross_bin_HdrC':(False,True),
        'two_competing_Mvh_cores':(False,True),
        'delete_EchA':(True,False),
        'delete_EchB':(True,False),
        'delete_EchF':(True,False),
        'destroy_EchF_ligands':(True,False),
        'truncate_EchF':(True,False),
        'delete_EchD_fused_region':(True,False),
        'move_EchB_away':(True,False),
        'unconfirmed_Ech_identity':(True,False),
        'destroy_Ech_catalytic_ligands':(True,False),
    }
    genes=[];sequences={};identities=[];expectations=[]
    for case,expected in cases.items():
        local_sequences={n:original_sequences[by_number[n]] for n in numbers}
        deleted={'delete_HdrC':'1137','delete_EchA':'2001','delete_EchB':'1997','delete_EchF':'1998'}
        if case in deleted:del local_sequences[deleted[case]]
        if case=='destroy_HdrC_ligands':
            chars=list(local_sequences['1137'])
            for p in [24,27,30,34,69,72,75,79]:assert chars[p-1]=='C';chars[p-1]='S'
            local_sequences['1137']=''.join(chars)
        if case=='replace_HdrC_with_NuoI':local_sequences['1137']=local_sequences['1998']
        if case=='delete_MvhD_fused_region':local_sequences['1163']=local_sequences['1163'][:622]
        if case in {'split_remote_HdrA_MvhD','unlinked_remote_HdrA_MvhD','opposite_strand_remote_HdrA_MvhD','cross_bin_remote_MvhD'}:
            local_sequences['relay']=local_sequences['1163'][629:]
            local_sequences['1163']=local_sequences['1163'][:622]
        if case=='replace_HdrA_with_F420_donor':local_sequences['1163']=donors+'GGG'+local_sequences['1163'][629:]
        if case=='two_competing_Mvh_cores':local_sequences['second_core']=local_sequences['1786']
        if case=='destroy_EchF_ligands':
            chars=list(local_sequences['1998'])
            for p in [49,52,55,59,82,85,88,92]:assert chars[p-1]=='C';chars[p-1]='S'
            local_sequences['1998']=''.join(chars)
        if case=='truncate_EchF':local_sequences['1998']=local_sequences['1998'][:66]
        if case=='delete_EchD_fused_region':local_sequences['1996']=local_sequences['1996'][148:]
        if case=='destroy_Ech_catalytic_ligands':
            positions={int(s['query_position']) for s in table(source/'annotation/catalytic_sites.tsv')
                       if s['protein_id']==by_number['1996'] and s['residue_role']=='direct_metal_ligand' and s['status']=='present' and s['observed_residue'] in {'C','U'}}
            assert len(positions)==4,positions
            chars=list(local_sequences['1996'])
            for p in positions:chars[p-1]='S'
            local_sequences['1996']=''.join(chars)
        for n,seq in local_sequences.items():
            origin='1163' if n=='relay' else '1786' if n=='second_core' else n
            pid=case+'__gene_'+n;g=copy.deepcopy(original_genes[by_number[origin]])
            g.update(id=pid,contig=case+'__chromosome',bin_id=case)
            if case=='cross_bin_HdrC' and n=='1137':g['bin_id']='foreign_bin'
            if case=='move_EchB_away' and n=='1997':g.update(start=100000,end=100000+3*len(seq)-1)
            if n=='relay':g.update(start=1162000,end=1162000+3*len(seq)-1)
            if n=='relay' and case=='unlinked_remote_HdrA_MvhD':g.update(start=1250000,end=1250000+3*len(seq)-1)
            if n=='relay' and case=='opposite_strand_remote_HdrA_MvhD':g['strand']=-int(g['strand'])
            if n=='relay' and case=='cross_bin_remote_MvhD':g['bin_id']='foreign_bin'
            if n=='second_core':g.update(start=1805000,end=1805000+3*len(seq)-1)
            if case=='circular_rotation':
                length=int(g['contig_length']);g.update(start=(int(g['start'])+900000-1)%length+1,end=(int(g['end'])+900000-1)%length+1)
            genes.append(g);sequences[pid]=seq
            if n in {'1786','1996'}:
                identity=dict(anchors[by_number[n]],candidate_id=pid)
                if case=='unconfirmed_Ech_identity' and n=='1996':identity['identity_status']='evidence_insufficient'
                identities.append(identity)
                expectations.append(dict(case=case,candidate_id=pid,function_id='Ni13' if n=='1786' else 'Ni23',expected_assignment=expected[0 if n=='1786' else 1]))
    split_cases={}
    if split_source is not None:
        split_source=Path(split_source).resolve()
        original={r['id']:r for r in table(split_source/'annotation/genes.tsv')}
        seqs={r.id:str(r.seq) for r in SeqIO.parse(split_source/'annotation/proteins.faa','fasta')}
        numbers={p.rsplit('_',1)[-1]:p for p in original}
        anchor=numbers['928'];identity=next(r for r in table(split_source/'functions/functional_annotations.tsv') if r['candidate_id']==anchor)
        old=json.loads((split_source/'functions/function_rules/evidence.json').read_text())[anchor]
        template=max((r for r in old if r['function_id']=='Ni13'),key=lambda r:r['score'])
        assert not template['missing_required']
        ligand_positions={m['role']:[s['query_position'] for s in m.get('coordination_evidence',{}).get(m['role'],[])] for m in template['members']}
        catalytic_sites={int(s['query_position']) for s in table(split_source/'annotation/catalytic_sites.tsv')
                         if s['protein_id']==anchor and s['residue_role']=='direct_metal_ligand' and s['status']=='present' and s['observed_residue'] in {'C','U'}}
        assert len(catalytic_sites)==4
        split_cases={name:name in {'natural_split','fuse_HdrA_MvhD','renamed_rotation'} for name in [
            'natural_split','fuse_HdrA_MvhD','renamed_rotation','remote_unlinked_MvhD','opposite_strand_MvhD',
            'delete_MvhD','delete_HdrA','delete_HdrC','cross_bin_MvhD','two_Mvh_cores','F420_donor_replacement',
            'HdrC_coordination_loss','HdrB_coordination_loss','catalytic_ligand_loss','unconfirmed_core']}
        for name,expected in split_cases.items():
            scope='split_source_'+name;seq={n:seqs[numbers[n]] for n in ['530','531','563','564','927','928']}
            deletion={'delete_MvhD':'530','delete_HdrA':'531','delete_HdrC':'564'}
            if name in deletion:del seq[deletion[name]]
            if name=='fuse_HdrA_MvhD':seq['531']+='GGG'+seq.pop('530')
            if name=='F420_donor_replacement':seq['531']=donors
            if name=='two_Mvh_cores':seq['second_core']=seq['928']
            for case,role,n in [('HdrC_coordination_loss','HdrC','564'),('HdrB_coordination_loss','HdrB','563'),('catalytic_ligand_loss','MvhA','928')]:
                if name==case:
                    positions=catalytic_sites if role=='MvhA' else ligand_positions[role]
                    assert positions
                    chars=list(seq[n])
                    for p in positions:assert chars[p-1] in 'CU';chars[p-1]='S'
                    seq[n]=''.join(chars)
            for n,sequence in seq.items():
                pid=scope+'__'+n;g=copy.deepcopy(original[numbers['928' if n=='second_core' else n]])
                g.update(id=pid,contig=scope+'__chromosome',bin_id=scope)
                if name=='remote_unlinked_MvhD' and n=='530':g.update(start=400000,end=400000+3*len(sequence)-1)
                if name=='opposite_strand_MvhD' and n=='530':g['strand']=-int(g['strand'])
                if name=='cross_bin_MvhD' and n=='530':g['bin_id']='foreign_bin'
                if n=='second_core':g.update(start=1000000,end=1000000+3*len(sequence)-1)
                if name=='renamed_rotation':
                    length=int(g['contig_length']);g.update(start=(int(g['start'])+800000-1)%length+1,end=(int(g['end'])+800000-1)%length+1)
                genes.append(g);sequences[pid]=sequence
                if n=='928':
                    identities.append(dict(identity,candidate_id=pid,identity_status='evidence_insufficient' if name=='unconfirmed_core' else identity['identity_status']))
                    expectations.append(dict(case=scope,candidate_id=pid,function_id='Ni13',expected_assignment=expected))
    (annotation/'proteins.faa').write_text(''.join('>'+pid+'\n'+seq+'\n' for pid,seq in sequences.items()))
    write(annotation/'genes.tsv',genes,list(genes[0]))
    write(annotation/'catalytic_anchors.tsv',[],['protein_id','contig','bin_id','catalytic_family'])
    manifest=dict(status='complete',input_mode='bin',normalized_proteins_sha256=sha(annotation/'proteins.faa'),
                  output_sha256={name:sha(annotation/name) for name in ['genes.tsv','catalytic_anchors.tsv']},
                  scope='Synthetic perturbations of real sample proteins and coordinates; not an independent genome accuracy benchmark.')
    (annotation/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Scanning fresh query profiles for '+str(len(sequences))+' fixture proteins',flush=True)
    annotate_domains([annotation],out/'domains')
    predictions,_=annotate(annotation,out/'functions/function_rules',identities,2,out/'domains')
    for row in expectations:
        p=predictions[row['candidate_id']];observed=row['function_id'] in p['function_class_ids']
        row.update(observed_assignment=observed,status=p['function_rule_status'],passed=observed==row['expected_assignment'])
        if row['function_id']=='Ni23':assert not set(p['function_class_ids']) & {'Ni19','Ni20'},p
    write(out/'controls.tsv',expectations,list(expectations[0]))
    report=dict(cases=len(cases)+len(split_cases),assertions=len(expectations),passed=sum(r['passed'] for r in expectations),
                failures=[r for r in expectations if not r['passed']],
                policy_sha256=sha(ROOT/'config/function_rules_34.json'),engine_sha256=sha(ROOT/'research_workflow/function_rules.py'),
                scope=manifest['scope'])
    (out/'summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False),flush=True)
    assert not report['failures'],report['failures']
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--split-source-run',type=Path)
    args=p.parse_args();validate(args.source_run,args.out,args.split_source_run)
