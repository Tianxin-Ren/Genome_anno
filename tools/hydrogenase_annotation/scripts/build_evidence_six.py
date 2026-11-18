#!/usr/bin/env python3
"""Build an accession-linked evidence successor; never modify classifier/model sources."""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import reorganize_evidence as catalog

ROOT = catalog.ROOT
P = catalog.P
CURATION = P/'references/curation/evidence_six_20261112'
OUT = ROOT/'outputs/01a18a6b-afa7-7331-abea-680ebf4e2278/evidence_six_20261112'
WORK = ROOT/'.nas/evidence_six_20261112'


def read(path):
    return json.loads(path.read_text())


def write_tsv(path, rows):
    assert rows
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)


def build_catalog():
    catalog.OUT = OUT; catalog.WORK = WORK
    catalog.EXTRA_PAPERS = read(CURATION/'added_studies.json')
    catalog.EXTRA_SYSTEMS = read(CURATION/'added_systems.json')
    if (CURATION/'genomic_helpers.json').exists():catalog.EXTRA_HELPERS=read(CURATION/'genomic_helpers.json')
    if (CURATION/'study_system_bindings.json').exists():
        catalog.EXTRA_BINDINGS={r['function']+'|'+catalog.citation_key(r['citation']):r['system'] for r in read(CURATION/'study_system_bindings.json')}
    catalog.build()


def annotate_features():
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from audit_maturation_domains import scan_profile_domains
    from active_centers import annotate_active_centers

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'features').mkdir(exist_ok=True)
    payload = read(WORK/'catalog.json')
    rows = payload['sheets']['sequences.tsv']
    query = OUT/'features/search_U_to_X.faa'
    SeqIO.write([SeqRecord(Seq(r['完整CDS蛋白序列'].replace('U','X')),id=r['序列ID'],description='search_only') for r in rows],query,'fasta')
    domains=[]
    for kind,policy_name,model_path in [
        ('reference','reference_domain_policy.json','references/processed/reference_domain_audit/reference_domains.hmm'),
        ('maturation','maturation_domain_policy.json','references/processed/maturation_domain_audit/maturation_domains.hmm'),
    ]:
        policy=read(P/'config'/policy_name)
        models={}
        for m in policy['models']:
            filename=m.get('file',m['accession'].split('.')[0]+'.hmm')
            path=P/'references/raw/domain_review'/filename
            if not path.exists():
                path=P/'references/raw/reference_domain_review'/filename
            assert hashlib.sha256(path.read_bytes()).hexdigest()==m['sha256'],str(path)
            ga=next(line.split()[1:] for line in path.read_text().splitlines() if line.startswith('GA '))
            fields={line.split()[0]:' '.join(line.split()[1:]) for line in path.read_text().split('HMM          ')[0].splitlines()[1:] if line.strip()}
            models[m['profile']]={**m,'accession':fields['ACC'],'model_length':int(fields['LENG']),
                'sequence_ga':float(ga[0].strip(';')),'domain_ga':float(ga[1].strip(';'))}
        found=scan_profile_domains(P/model_path,query,models,policy['min_model_coverage'],OUT/'features'/kind,2)
        domains.extend({'序列ID':qid,'模型集':kind,**d} for qid,ds in found.items() for d in ds)
    write_tsv(OUT/'protein_domains.tsv',domains)
    catalysts={r['序列ID'] for r in rows if '催化亚基' in r['组件类别']}
    selected,_=annotate_active_centers(OUT/'proteins.faa',OUT/'features',threads=2,anchors=catalysts)
    sites=[{'序列ID':r['protein_id'],'催化位点状态':r['active_center_status'],**s} for r in selected for s in r['sites']]
    write_tsv(OUT/'catalytic_site_details.tsv',sites)
    (WORK/'mapped_features.json').write_text(json.dumps({'domains':domains,'centers':selected},ensure_ascii=False))
    print('FEATURES',len(domains),'domain hits;',len(selected),'catalytic chains;',len(sites),'site rows',flush=True)


def grade_record(e):
    """Code the evidence method, separately from whether the full function was measured."""
    text=e['实际测得'];old=e['原证据标签']
    if (e['支持关系'] in {'反例或特异性对照','不支持该类别','更正/背景'} or old in {'EXCLUDED_WRONG_ENZYME','DUPLICATE'}
        or e['功能ID']=='Fe08' and 'membrane_NiFe' in e['体系ID']):
        return 'X','', '对照或更正，不计功能支持证据'
    if e['证据级别'].startswith('D'):
        return 'E4','ECO:0000305','候选或资料不足，不计实验支持'
    catalytic_readout=r'酶活|比活|催化活性|催化性质|活性(?!中心|位点|部位|胶|染色|钯)|反应速率|周转|动力学|ATP合成|ATP生成|还原速率|产氢(?!系统)|产H[₂2]|产生\s*H[₂2]|H[₂2]\s*产生|驱动H[₂2]|测H[₂2]|检测H[₂2]|定量H[₂2]|U/mg|nmol/min/mg|↔|→|共同驱动H[₂2]'
    structural_or_comp=re.search(r'晶体|X射线|冷冻?电镜|光谱|中点电位|N端|测序|Northern|结构\d|活性中心|活性位点|氧化还原模块|共纯化|相互作用',text,re.I)
    activity=re.search(catalytic_readout+r'|消耗|H2生成|氢气生成|转运(?!蛋白)|电化学',text)
    direct_biochem=(
        re.search(r'膜囊泡|膜泡|脂质体|电极|tBLM|重构(?!.*结构)|未处理膜|提取物的.*活性|部分纯化.*活性|吡啶核苷酸.*活性',text)
        or (re.search(r'纯化(?!.*(?:共纯化|晶体|电镜|N端))',text) and re.search(catalytic_readout+r'|电化学|还原(?!型|剂|模块)|氧化(?!态|还原模块)',text) and not structural_or_comp)
    )
    genetic=re.search(r'敲除|缺失|knockout|deletion|突变体|突变株|相关突变|遗传干预|基因干预|失活株|启动子替换|调控改造|调控扰动|回导|基因簇转入|Δ[A-Za-z]',text,re.I)
    if direct_biochem:
        return 'E1','ECO:0000314','直接测定所列端点；不等于该功能所有限定均已测定'
    if genetic:
        return 'E2','ECO:0000315','成员干预/突变表型；直接化学反应的范围另列'
    if re.search(r'表达|转录|蛋白质组|proteom|转录组|活性染色|活性胶',text,re.I) and not activity:
        return 'E3','ECO:0000270','表达/丰度关联；不能独自指定直接电子伙伴'
    if structural_or_comp and not activity:
        return 'E3','ECO:0000269','结构/光谱/成员实验支持；不等于直接完整反应测定'
    if re.search(r'酶活|活性(?!中心|位点|部位|胶|染色|钯)|反应(?!体系)|还原(?!型|剂|模块)|氧化(?!态|还原模块)|产氢(?!系统)|产H[₂2]|H[₂2]\s*产生|动力学|电化学|膜囊泡|膜泡|脂质体|转运(?!蛋白)|ATP合成|ATP生成',text):
        return 'E1','ECO:0000314','直接测定所列端点；不等于该功能所有限定均已测定'
    if re.search(r'表达|转录|蛋白质组|proteom|转录组',text,re.I):
        return 'E3','ECO:0000270','表达/丰度关联；不能独自指定直接电子伙伴'
    return 'E3','ECO:0000269','其他实验支持；采用通用EXP，未猜测更具体实验类型'


def export_construct_fragments(seqs):
    path=CURATION/'construct_fragments.json'
    if not path.exists():return []
    from Bio.Seq import Seq
    specs=read(path)
    by_locus={(r['核苷酸登录号'],r['locus_tag']):r for r in seqs.values()}
    rows=[]
    faa_lines=[]; fna_lines=[]
    for sp in specs:
        parent=by_locus[(sp['parent_accession'],sp['parent_locus_tag'])]
        s1,e1=int(sp['aa_start_1based']),int(sp['aa_end_1based'])
        aa_full=parent['完整CDS蛋白序列']; cds_full=parent['CDS核酸序列']
        assert 1<=s1<=e1<=len(aa_full),(sp['fragment_id'],s1,e1,len(aa_full))
        frag_aa=aa_full[s1-1:e1]
        nt_s0=(s1-1)*3; nt_e0=e1*3
        assert nt_e0<=len(cds_full),(sp['fragment_id'],nt_e0,len(cds_full))
        frag_cds=cds_full[nt_s0:nt_e0]
        assert str(Seq(frag_cds).translate())==frag_aa,sp['fragment_id']
        aa_sha=hashlib.sha256(frag_aa.encode()).hexdigest()
        cds_sha=hashlib.sha256(frag_cds.encode()).hexdigest()
        rows.append({
            '片段ID':sp['fragment_id'],'父序列ID':parent['序列ID'],'体系ID':sp['system_id'],'功能ID':sp['function_id'],
            '亚基或构建名称':sp['role'],'参考文献':sp['citation'],'微生物':parent['微生物'],'参考菌株':parent['参考菌株'],
            '核苷酸登录号':parent['核苷酸登录号'],'locus_tag':parent['locus_tag'],'蛋白登录号':parent['蛋白登录号'],
            'UniProt参考':sp.get('uniprot_reference',''),'母蛋白全长_aa':len(aa_full),
            '片段起点_1based':s1,'片段终点_1based':e1,'片段长度_aa':len(frag_aa),'片段蛋白序列':frag_aa,
            '片段CDS长度_nt':len(frag_cds),'片段CDS序列':frag_cds,'片段蛋白SHA256':aa_sha,'片段CDS_SHA256':cds_sha,
            '边界依据':sp['boundary_source'],'说明':sp['notes'],
        })
        faa_lines.append(f">{sp['fragment_id']} parent={parent['序列ID']} {parent['蛋白登录号']}:{s1}..{e1} {parent['核苷酸登录号']} {parent['locus_tag']}\n"+'\n'.join(frag_aa[i:i+80] for i in range(0,len(frag_aa),80)))
        fna_lines.append(f">{sp['fragment_id']} parent={parent['序列ID']} {parent['locus_tag']}:nt{nt_s0+1}..{nt_e0} {parent['核苷酸登录号']}\n"+'\n'.join(frag_cds[i:i+80] for i in range(0,len(frag_cds),80)))
    if rows:
        write_tsv(OUT/'construct_fragments.tsv',rows)
        (OUT/'construct_fragments.faa').write_text('\n'.join(faa_lines)+'\n')
        (OUT/'construct_fragments_cds.fna').write_text('\n'.join(fna_lines)+'\n')
    return rows


def summarize():
    from evidence_domain_fragments import summarize_domains, structure_correspondences
    payload=read(WORK/'catalog.json'); mapped=read(WORK/'mapped_features.json')
    seqs={r['序列ID']:r for r in payload['sheets']['sequences.tsv']}
    systems={r['体系ID']:r for r in payload['sheets']['systems_genomes.tsv']}
    reference_policy=read(P/'config/reference_domain_policy.json')
    domain=summarize_domains(mapped['domains'],reference_policy)
    architecture=[{'序列ID':qid,**r} for qid,rs in domain.items() for r in rs]
    write_tsv(OUT/'protein_domain_architecture.tsv',architecture)
    structural=structure_correspondences(list(seqs.values()),P)
    centers={r['protein_id']:r for r in mapped['centers']}
    templates={r['reference_id']:r for r in read(P/'references/processed/active_centers/models.json')}
    scaffold={}
    for qid,c in centers.items():
        m=templates[c['reference_id']]
        if (m['window_start']==1 and m['window_end']==len(m['sequence']) and
            c['window_coverage']>=.9 and c['window_identity']>=30 and
            c['active_center_status']=='template_ligand_set_supported' and c['template_span_alignment_status']=='supported'):
            scaffold[qid]={'template':c['reference_id'],'source':m['evidence_url'],
                'template_coverage':c['window_coverage'],'identity_percent':c['window_identity'],
                'scope':'计算同源结构支架支持；完整实验模板跨度及配体均有对应，不是该被测蛋白自身结构的直接实测，也不证明功能。'}
    if scaffold:write_tsv(OUT/'catalytic_scaffold_correspondence.tsv',[{'序列ID':qid,**r} for qid,r in scaffold.items()])
    for qid,r in seqs.items():
        ds=domain.get(qid,[]); c=centers.get(qid); full=[d for d in ds if d['supported']]
        r['结构域检测状态']='单片段/共线片段链达到既有模型覆盖阈值' if full else '模型覆盖不足，保留局部命中；不据此断言蛋白截断' if ds else '已扫描，未检出统计支持的已用模型'
        r['结构域及坐标']=';'.join(f"{d['accession']}:{d['query_segments']}:{d['profile']}:union_coverage={d['model_coverage']:.3f}:{d['status']}" for d in ds)
        r['已有实验结构对应']=json.dumps(structural.get(qid,[]),ensure_ascii=False)
        r['完整催化模板同源支架对应']=json.dumps(scaffold.get(qid,{}),ensure_ascii=False)
        r['催化中心检测状态']=c['active_center_status'] if c else '未检出可用映射' if '催化亚基' in r['组件类别'] else '非氢酶催化链，不适用'
        r['催化模板跨度状态']=c['template_span_alignment_status'] if c else '未评估或不适用'
        r['催化位点_序列映射']=';'.join(f"{s['query_position']}:{s['observed_residue']}:{s['residue_role']}:{s['status']}" for s in c['sites']) if c else ''
        r['特征证据说明']='结构域与位点为本次计算检测，ECO:0000256/序列比对；不证明辅因子实际占据或酶活。'
    qualification=[]
    members=defaultdict(list)
    for link in payload['sheets']['evidence_sequence_links.tsv']:members[link['证据ID']].append(link)
    for e in payload['sheets']['evidence.tsv']:
        level,eco,scope=grade_record(e)
        ids=[x for field in ['催化亚基序列ID','小亚基序列ID','电子伙伴及相关亚基序列ID','成熟辅助基因序列ID'] for x in e[field].split(';') if x]
        cats=[x for x in e['催化亚基序列ID'].split(';') if x]
        aux=[x for x in e['成熟辅助基因序列ID'].split(';') if x]
        e['统一证据等级']=level;e['ECO方法代码']=eco;e['方法证据范围']=scope
        reasons=[]
        if level in {'X','E4'}:reasons.append('候选/对照/更正不作功能支持')
        if not ids or not e['参考基因组或核苷酸记录']:reasons.append('缺基因组与成员序列对应')
        if not cats:reasons.append('缺催化链序列')
        if e['序列对应级别'].startswith('S3'):reasons.append('参考菌株与实验来源尚未对应')
        if cats and not all(any(d['supported'] for d in domain.get(x,[])) or structural.get(x) or scaffold.get(x) for x in cats):reasons.append('催化链结构域覆盖不足且缺完整模板支架支持')
        if cats and not all(x in centers and
            all(s['query_position'] is not None for s in centers[x]['sites'] if s['residue_role']=='direct_metal_ligand')
            for x in cats):reasons.append('催化配体位点尚未完整定位')
        basic=not reasons
        if not aux:reasons.append('缺成熟辅助基因序列')
        system=systems.get(e['体系ID'],{})
        if '仅绑定已知催化亚基' in system.get('范围与限制',''):reasons.append('多亚基组成尚未明确')
        if system.get('本地GenBank文件'):
            source=Path(system['本地GenBank文件'])
            with source.open() as f:header=''.join(next(f,'') for _ in range(25))
            if re.search(r'cDNA|mRNA|partial sequence|complete cds',header,re.I) and not re.search(r'complete genome|chromosome|Assembly:|WGS',header,re.I):
                reasons.append('仅局部或转录本记录，尚缺全基因组原位对应')
        else:reasons.append('缺本地基因组来源')
        e['验收状态']='六项资料齐全的实验支持记录' if not reasons else '待补齐'
        e['待补字段']=';'.join(reasons)
        e['催化链结构域']=';'.join(seqs[x]['结构域及坐标'] for x in cats)
        e['催化位点及状态']=';'.join(seqs[x]['催化位点_序列映射'] for x in cats)
        e['催化中心风险标记']=';'.join(x+':'+centers[x]['active_center_status'] for x in cats if x in centers and centers[x]['active_center_status'] not in
            {'template_ligand_set_supported','experimentally_active_exact_reference_noncanonical_site'})
        e['催化位点解释']='位点为计算映射；非典型或比对不足单独标记，不以计算不确定性否定已存在的实验。带风险标记的映射不得直接作为确定活性位点或自动功能迁移依据。'
        e['辅助基因说明']=';'.join(seqs[x]['亚基或基因名']+'='+seqs[x]['locus_tag'] for x in aux) or '未找到或尚未对应，不视为不存在'
        e['辅助基因范围']='原生基因组成熟候选，可共享；不是全部同一操纵子，也不代表本研究逐个实验验证'
        e['亚基组成说明']=';'.join(r['亚基或基因名']+'='+seqs[r['序列ID']]['locus_tag'] for r in members[e['证据ID']] if '成熟' not in r['组件类别'])
        e['各成员结构域结果']=';'.join(x+'='+seqs[x]['结构域检测状态'] for x in dict.fromkeys(ids))
        e['计数适用范围']='按DOI/PMID去重的实验支持记录，允许间接证据；资料齐全不等于每篇独立测全类别完整命题'
        qualification.append({'证据ID':e['证据ID'],'功能ID':e['功能ID'],'统一证据等级':level,'ECO方法代码':eco,'基础对应通过':basic,'六项资料齐全':not reasons,'待补字段':e['待补字段']})
    write_tsv(OUT/'evidence.tsv',payload['sheets']['evidence.tsv'])
    write_tsv(OUT/'sequences.tsv',list(seqs.values()))
    export_construct_fragments(seqs)
    write_tsv(OUT/'record_completion.tsv',qualification)
    ready=[e for e in payload['sheets']['evidence.tsv'] if not e['待补字段']]
    if ready:write_tsv(OUT/'complete_evidence.tsv',ready)
    summary=[]
    for f in payload['sheets']['function_summary.tsv']:
        es=[e for e in payload['sheets']['evidence.tsv'] if e['功能ID']==f['功能ID']]
        base=[e for e in es if next(q['基础对应通过'] for q in qualification if q['证据ID']==e['证据ID'])]
        n=len({catalog.citation_key(e['DOI或文献键']) for e in base})
        complete=[e for e in es if not e['待补字段']]
        cn=len({catalog.citation_key(e['DOI或文献键']) for e in complete})
        summary.append({'功能ID':f['功能ID'],'功能类别':f['功能类别'],'收录记录数':len(es),'基因组序列结构域位点对应研究数':n,
            '六项资料齐全研究数':cn,'目标研究数':6,'仍需补齐研究数':max(0,6-cn),
            'E1直接端点研究数':len({catalog.citation_key(e['DOI或文献键']) for e in complete if e['统一证据等级']=='E1'}),
            'E2遗传研究数':len({catalog.citation_key(e['DOI或文献键']) for e in complete if e['统一证据等级']=='E2'}),
            'E3结构表达等研究数':len({catalog.citation_key(e['DOI或文献键']) for e in complete if e['统一证据等级']=='E3'}),
            '计数范围':'六项资料齐全的实验支持；非全部完整命题独立验证'})
    write_tsv(OUT/'completion_summary.tsv',summary)
    payload['qualification_summary']=summary
    (WORK/'enriched_catalog.json').write_text(json.dumps(payload,ensure_ascii=False))
    print('COMPLETENESS',json.dumps({'records':len(qualification),'records_with_basic_linkage':sum(r['基础对应通过'] for r in qualification),
        'classes_with_six_basic_linked_studies':sum(r['基因组序列结构域位点对应研究数']>=6 for r in summary),
        'complete_records':len(ready),'classes_with_six_complete_studies':sum(r['六项资料齐全研究数']>=6 for r in summary),
        'remaining_class_study_slots':sum(r['仍需补齐研究数'] for r in summary)},ensure_ascii=False),flush=True)
    assert all(e['成熟辅助基因序列ID'] and e['亚基组成说明'] and e['催化位点及状态'] and e['各成员结构域结果'] for e in ready)
    assert all(r['统一证据等级'] in {'E1','E2','E3'} for r in ready)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('step',choices=['catalog','features','summary','all'])
    args=parser.parse_args()
    if args.step in {'catalog','all'}:build_catalog()
    if args.step in {'features','all'}:annotate_features()
    if args.step in {'summary','all'}:summarize()
