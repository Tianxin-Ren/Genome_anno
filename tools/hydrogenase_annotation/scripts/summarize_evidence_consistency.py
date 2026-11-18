#!/usr/bin/env python3
"""Compare reference compositions, without treating copied memberships as replication."""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path('/Volumes/Mac_FN/Hydrogenase')
SOURCE=ROOT/'outputs/01a18a6b-afa7-7331-abea-680ebf4e2278/evidence_six_20261112'
OUT=SOURCE.parent/'function_consistency_20261112'


def read(name):
    with (SOURCE/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def write(name,rows):
    assert rows
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def intersection(sets):
    return set.intersection(*sets) if sets else set()


def helper_modules(names):
    # Name normalization is descriptive; it does not experimentally prove maturation roles.
    result=set()
    for name in names:
        if name.lower()=='hydef':result.update(['hydE','hydF'])
        for family,letters in [('hyp','ABCDEF'),('hyd','EFG'),('hcg','ABCDEFG')]:
            for letter in letters:
                if re.search(r'(?i)(?:^|_)'+family+letter+r'(?:$|_)',name):result.add(family+letter)
    return result


def main():
    assert ' on /Volumes/Mac_FN (smbfs,' in __import__('subprocess').check_output(['/sbin/mount'],text=True)
    assert helper_modules(['HydEF','hydG'])==helper_modules(['hydE','hydF','HydG'])
    evidence=read('complete_evidence.tsv')
    seqs={r['序列ID']:r for r in read('sequences.tsv')}
    systems={r['体系ID']:r for r in read('systems_genomes.tsv')}
    links=defaultdict(list)
    for r in read('evidence_sequence_links.tsv'):links[r['证据ID']].append(r)
    domains=defaultdict(set)
    for r in read('protein_domain_architecture.tsv'):
        if r['supported'].lower()=='true':domains[r['序列ID']].add(r['accession'])
    by_system=defaultdict(list)
    for e in evidence:by_system[e['体系ID']].append(e)
    fingerprints={}
    for sid,es in by_system.items():
        members=links[es[0]['证据ID']]
        core=[r for r in members if r['组件类别'] not in {'成熟辅助因子','调控辅助基因','生理背景基因'}]
        aux=[r for r in members if r['组件类别']=='成熟辅助因子']
        cats=[r for r in members if r['组件类别']=='催化亚基']
        assert core and cats
        signatures={tuple(sorted((r['亚基或基因名'],r['序列ID'],r['组件类别']) for r in links[e['证据ID']])) for e in es}
        fingerprints[sid]={'core':{r['亚基或基因名'] for r in core},'aux_names':{r['亚基或基因名'] for r in aux},
            'aux_modules':helper_modules(r['亚基或基因名'] for r in aux),
            'core_ids':{r['序列ID'] for r in core},'aux_ids':{r['序列ID'] for r in aux},
            'cat_hashes':tuple(sorted({seqs[r['序列ID']]['蛋白SHA256'] for r in cats})),
            'cat_organisms':{seqs[r['序列ID']]['微生物'] for r in cats},
            'cat_domains':{d for r in cats for d in domains[r['序列ID']]},
            'reference_sets_per_system':len(signatures),
            'scope':{r['成员功能范围'] for r in members if r in core}}
    rows=[];details=[]
    for f in read('completion_summary.tsv'):
        fid=f['功能ID'];es=[e for e in evidence if e['功能ID']==fid];sids=sorted({e['体系ID'] for e in es})
        fs=[fingerprints[sid] for sid in sids]
        core_sets={tuple(sorted(x['core'])) for x in fs}
        aux_sets={tuple(sorted(x['aux_modules'])) for x in fs}
        centers=sorted({systems[sid]['活性中心'] for sid in sids})
        category='仅一套参考，不能评价跨体系一致性' if len(sids)==1 else '多套参考，成员名称组合相同' if len(core_sets)==1 else '多套参考，成员名称组合不同，须按模块比较'
        shared_core=intersection([x['core'] for x in fs]);shared_aux=intersection([x['aux_modules'] for x in fs])
        role_issue=fid=='Fe08' and any(x in fingerprints.get('FeFe_Mthermo_HydABCDEF',{}).get('aux_names',set()) for x in ['HydE','HydF'])
        rows.append({'功能ID':fid,'功能类别':f['功能类别'],'齐全记录数':len(es),'去重文献数':len({e['DOI或文献键'].lower() for e in es}),
            '参考体系数':len(sids),'不同催化链组合数_SHA去重':len({x['cat_hashes'] for x in fs}),
            '催化链来源物种标签数':len(set().union(*(x['cat_organisms'] for x in fs))),
            '催化中心类别':';'.join(centers),'参考核心及伙伴名称组合数':len(core_sets),
            '参考核心及伙伴成员数范围':str(min(len(x['core_ids']) for x in fs))+'–'+str(max(len(x['core_ids']) for x in fs)),
            '共同核心及伙伴名称_非必需规则':';'.join(sorted(shared_core)),
            '成熟模块名称组合数':len(aux_sets),'共同成熟模块_名称归一化':';'.join(sorted(shared_aux)),
            '成熟因子参考序列数范围':str(min(len(x['aux_ids']) for x in fs))+'–'+str(max(len(x['aux_ids']) for x in fs)),
            '催化链Pfam组合数':len({tuple(sorted(x['cat_domains'])) for x in fs}),
            '带催化映射风险标记的记录数':sum(bool(e['催化中心风险标记']) for e in es),
            '分类':category,
            '已知组件分类问题':'Moorella HydABCDEF的HydE/HydF被源表归为成熟辅助因子；原研究为六亚基复合体成员。当前4成员统计不代表四亚基酶。' if role_issue else '',
            '范围':'比较的是参考成员清单；名称不等于同源模块，未逐论文独立核定组成、必要性或化学计量。'})
        for sid in sids:
            x=fingerprints[sid]
            details.append({'功能ID':fid,'体系ID':sid,'参考体系':systems[sid]['微生物及体系'],
                '活性中心':systems[sid]['活性中心'],'关联齐全研究数':len({e['DOI或文献键'].lower() for e in es if e['体系ID']==sid}),
                '参考核心及伙伴名称':';'.join(sorted(x['core'])),'成熟因子名称':';'.join(sorted(x['aux_names'])),
                '成熟模块_名称归一化':';'.join(sorted(x['aux_modules'])),'催化链Pfam':';'.join(sorted(x['cat_domains'])),
                '体系序列限制':systems[sid]['范围与限制'],'成员证据范围':';'.join(sorted(x['scope']))})
    assert len(rows)==34 and sum(r['齐全记录数'] for r in rows)==len(evidence)
    assert all(r['去重文献数']>=6 for r in rows)
    shared_helpers=defaultdict(set)
    for sid,x in fingerprints.items():
        for qid in x['aux_ids']:shared_helpers[qid].add(sid)
    report={'complete_records':len(evidence),'classes':len(rows),'reference_systems':len(fingerprints),
        'classification_counts':dict(Counter(r['分类'] for r in rows)),
        'classes_with_multiple_centers':[r['功能ID'] for r in rows if ';' in r['催化中心类别']],
        'classes_with_different_helper_module_sets':[r['功能ID'] for r in rows if r['成熟模块名称组合数']>1],
        'systems_with_multiple_papers':sum(len({e['DOI或文献键'].lower() for e in es})>1 for es in by_system.values()),
        'systems_with_multiple_papers_and_identical_reference_members':sum(len({e['DOI或文献键'].lower() for e in by_system[sid]})>1 and x['reference_sets_per_system']==1 for sid,x in fingerprints.items()),
        'shared_maturation_sequence_count':sum(len(s)>1 for s in shared_helpers.values()),
        'records_with_catalytic_mapping_risk':sum(bool(e['催化中心风险标记']) for e in evidence),
        'known_role_classification_issue':next((r['已知组件分类问题'] for r in rows if r['已知组件分类问题']),None),
        'input_sha256':{name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in ['complete_evidence.tsv','evidence_sequence_links.tsv','sequences.tsv','protein_domain_architecture.tsv']},
        'interpretation':'This is reference-composition diversity, not per-study experimental agreement, stoichiometry, or functional accuracy.'}
    OUT.mkdir(parents=True,exist_ok=True)
    write('function_consistency.tsv',rows);write('system_compositions.tsv',details)
    (OUT/'summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='input_sha256'},ensure_ascii=False),flush=True)
    for r in rows:print('CLASS',r['功能ID'],r['去重文献数'],r['参考体系数'],r['不同催化链组合数_SHA去重'],r['参考核心及伙伴名称组合数'],r['成熟模块名称组合数'],r['催化中心类别'],flush=True)


if __name__=='__main__':main()
