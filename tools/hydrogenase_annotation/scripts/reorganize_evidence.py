#!/usr/bin/env python3
"""Export existing hydrogenase evidence and accession-pinned sequences; no model changes."""
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path('/Volumes/Mac_FN/Hydrogenase')
P = ROOT / 'hydrogenase_annotation'
OUT = ROOT / 'outputs/01a18a6b-afa7-7331-abea-680ebf4e2278/evidence_catalog_20261112'
WORK = ROOT / '.nas/evidence_catalog_20261112'
EXTRA_PAPERS = []
EXTRA_SYSTEMS = []
EXTRA_HELPERS = {}
EXTRA_BINDINGS = {}

# Reaction cores are navigation labels; all original functional qualifiers remain in the class text.
CORE_GROUPS = [
    ('Ni01', 'H2与细胞色素c电子传递'),
    ('Ni02 Ni03 Ni05 Ni09', 'H2氧化与醌相关电子传递'),
    ('Ni04', '醌池驱动产H2'),
    ('Ni06 Ni08', '摄氢的生理用途，直接受体另判'),
    ('Ni07', 'H2与methanophenazine电子传递'),
    ('Ni10', 'H2信号感知'),
    ('Ni11', 'H2与F420氧化还原'),
    ('Ni12 Ni12R Fe08', 'H2与NADP氧化还原，Fd依赖性另列'),
    ('Ni13', 'H2与异二硫键/Fd电子分叉'),
    ('Ni14 Fe05', 'H2与NAD氧化还原，Fd依赖性另列'),
    ('Ni15 Ni16 Ni22 Fe06', 'H2与甲酸/CO2相互转换'),
    ('Ni17', 'CO氧化耦联产H2'),
    ('Ni18 Ni19 Ni20 Ni23 Fe01 Fe02', 'H2与Fd氧化还原'),
    ('Fe03 Fe04 Ni24', 'H2与NAD/Fd电子分叉或汇合'),
    ('Fe07', 'H2与NADP/Fd电子分叉'),
    ('H01', 'H2与methenyl-H4MPT氧化还原'),
    ('Ni21', 'H2经外部蛋白伙伴还原NADP'),
]
CORES = {fid: label for ids, label in CORE_GROUPS for fid in ids.split()}


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(name, rows):
    assert rows, name
    with open(OUT / name, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def citation_key(text):
    text = text.strip().rstrip('/').removeprefix('UNRESOLVED:')
    m = re.search(r'10\.\d{4,9}/[^\s；|]+', text, re.I)
    return m.group().rstrip('.,').lower() if m else text.lower()


def fasta(path):
    result, key, seq = {}, '', []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if key:
                result[key] = ''.join(seq)
            key, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if key:
        result[key] = ''.join(seq)
    return result


def extract_location(dna, loc):
    loc = re.sub(r'\s|[<>]', '', loc)
    if loc.startswith('complement(') and loc.endswith(')'):
        part = extract_location(dna, loc[11:-1])
        return part.translate(str.maketrans('ACGTRYMKBDHVN', 'TGCAYRKMVHDBN'))[::-1]
    if loc.startswith(('join(', 'order(')) and loc.endswith(')'):
        inside = loc[loc.index('(')+1:-1]
        pieces, start, depth = [], 0, 0
        for i, ch in enumerate(inside):
            depth += (ch == '(') - (ch == ')')
            if ch == ',' and depth == 0:
                pieces.append(inside[start:i]); start = i+1
        return ''.join(extract_location(dna, x) for x in pieces + [inside[start:]])
    m = re.fullmatch(r'(\d+)\.\.(\d+)', loc)
    if not m:
        raise ValueError('Unsupported CDS location: ' + loc)
    lo, hi = map(int, m.groups())
    assert 1 <= lo <= hi <= len(dna), (loc, len(dna))
    return dna[lo-1:hi]


GENOMES = {}


def genbank(rel):
    if rel in GENOMES:
        return GENOMES[rel]
    path = P / rel
    result = []
    for rec in path.read_text().split('\n//'):
        if not rec.strip() or '\nORIGIN' not in rec or '\nFEATURES' not in rec:
            continue
        head, tail = rec.split('\nORIGIN', 1)
        dna = ''.join(re.findall('[A-Za-z]+', tail)).upper()
        ver = re.search(r'^VERSION\s+(\S+)', head, re.M)
        acc = ver.group(1) if ver else re.search(r'^ACCESSION\s+(\S+)', head, re.M).group(1)
        org = re.search(r'^  ORGANISM\s+(.+)', head, re.M)
        strain = re.search(r'/strain="([^"]+)"', head)
        definition = re.search(r'^DEFINITION\s+(.+?)(?=^ACCESSION)', head, re.M | re.S)
        desc = re.sub(r'\s+', ' ', definition.group(1)).strip() if definition else ''
        genes = []
        for m in re.finditer(r'^     CDS\s+(.+?)(?=^     \S|\Z)', head, re.M | re.S):
            block = m.group(1)
            loc = re.sub(r'\s+', '', block.split('/')[0])
            qualifiers = [(k,re.sub(r'\s+', ' ', v).strip()) for k,v in re.findall(r'/(\w+)="([^"]*)"',block,re.S)]
            q = dict(qualifiers)
            protein = q.get('translation', '').replace(' ', '')
            if not protein:
                continue
            nums = [int(x) for x in re.findall(r'\d+', loc)]
            aliases = {a for k,v in qualifiers if k=='old_locus_tag' for a in re.split(r'[;, ]+',v)}
            aliases |= {q.get('gene', ''), q.get('locus_tag', ''), q.get('protein_id', '')}
            genes.append(dict(accession=acc, organism=org.group(1) if org else '', strain=strain.group(1) if strain else '',
                locus=q.get('locus_tag', q.get('protein_id', '')), gene=q.get('gene', ''), aliases=aliases,
                protein_accession=q.get('protein_id', ''), product=q.get('product', ''), sequence=protein,
                location=loc, start=min(nums), end=max(nums), strand=-1 if loc.startswith('complement') else 1,
                cds=extract_location(dna, loc), cds_status='GenBank CDS提取；保留原translation及转译例外，未重新预测', source=rel))
        result.append(dict(accession=acc, organism=org.group(1) if org else '', strain=strain.group(1) if strain else '',
                           description=desc, dna=dna, genes=genes, source=rel))
    assert result, path
    GENOMES[rel] = result
    print('read', rel, sum(len(x['genes']) for x in result), flush=True)
    return result


def select(rel, selector):
    matches = [g for record in genbank(rel) for g in record['genes'] if selector in g['aliases']]
    assert len(matches) == 1, (rel, selector, len(matches))
    return matches[0]


def component_type(role, kind, product):
    if kind == 'catalytic': return '催化亚基'
    if kind == 'maturation': return '成熟辅助因子'
    if kind == 'context': return '生理背景基因'
    if role=='HoxU': return '电子传递或其他相关亚基'
    if re.search(r'^(HyaA|HybO|HybS|HucS|HhyS|HysB|HydS|HyhS|HoxK|HoxY|HoxB|HupS|HupU|EchC|MbhJ|VhtG|FrhG|MvhG)$', role) or re.search(r'hydrogenase.*small subunit|small subunit.*hydrogenase', product, re.I): return '小亚基或近端电子传递亚基'
    if re.search(r'Fdh|Hdr|CooS|HydB|HytB|FrhB|HoxF|CAC0764|BMR', role, re.I): return '电子伙伴或底物耦联亚基'
    # Curated complex modules can share names with unrelated maturation proteins.
    if kind != 'module' and re.search(r'hyp[A-F]|hyd[EFG]|peptidase|matur', role+' '+product, re.I): return '成熟辅助因子'
    if role in {'HoxJ','HoxA','HyfR'} or re.search(r'HupT|HupR|HisK|regulat|kinase', role+' '+product, re.I): return '调控辅助基因'
    if re.search(r'membrane|cytochrome b', product, re.I) or role in {'HyaC','HybB','VhtC','HoxZ'}: return '膜或细胞色素亚基'
    return '电子传递或其他相关亚基'


def build():
    mounts = subprocess.check_output(['/sbin/mount'], text=True)
    assert ' on /Volumes/Mac_FN (smbfs,' in mounts
    OUT.mkdir(parents=True, exist_ok=True); WORK.mkdir(parents=True, exist_ok=True)
    funcs = read_tsv(P/'references/processed/function_evidence_grouped_20261111.tsv')
    function_by_id = {r['功能ID']: r for r in funcs}
    assert set(CORES)==set(function_by_id)
    assert extract_location('ACGTAACC', 'complement(join(1..2,5..8))')=='GGTTGT'
    papers = read_tsv(P/'references/raw/user_package_20261111/all_paper_evidence_updated.tsv')
    for r in papers:
        r['_source'] = 'all_paper_evidence_updated.tsv'
    for r in EXTRA_PAPERS:
        p = dict.fromkeys(papers[0], '')
        p.update(r)
        p['_source'] = 'evidence_six_20261112/added_studies.json'
        papers.append(p)
    # Later rows refine the same evidence rather than creating extra independent studies.
    groups = read_tsv(P/'references/processed/function_evidence_groups_20261111.tsv')
    singles = read_tsv(P/'references/processed/function_evidence_single_studies_20261111.tsv')
    for r in read_tsv(P/'references/raw/user_package_20261111/targeted_paper_audit_20261111.tsv'):
        key = citation_key(r['DOI或文献键'])
        matches = [p for p in papers if p['功能ID']==r['功能ID'] and citation_key(p['DOI或文献键'])==key and not p.get('_source','').startswith('evidence_six_20261112/')]
        if matches:
            for p in matches:
                p['实际测得']=r['实验实际测得']; p['尚未证明或限制']=r['限制或未测端点']
                p['原文位置']=r['原文位置或定位说明']; p['获取层级']=r['本轮读取层级']
                p['_source']+='；targeted_paper_audit_20261111.tsv'
    for s in singles:
        key = citation_key(s['研究DOI'])
        matches = [p for p in papers if p['功能ID']==s['功能ID'] and citation_key(p['DOI或文献键'])==key]
        if not matches:
            p = dict.fromkeys(papers[0], '')
            p.update({'功能ID':s['功能ID'],'DOI或文献键':s['研究DOI'],'物种或来源':s['微生物与体系'],
                '酶或体系':s['架构'],'体系去重键':'successor_'+s['证据ID'],'证据分级':'FULL_DIRECT',
                '文献题名或主题_可缩写':s['实测与位置'].split('：')[0]})
            papers.append(p); matches=[p]
        for p in matches:
            p.update({'实际测得':s['实测与位置'],'尚未证明或限制':s['限制'],'原文位置':s['实测与位置'],
                      '获取层级':s['核读层级'],'原始来源URL':s['来源'],'_source':p.get('_source','')+'；function_evidence_single_studies_20261111.tsv'})
            p['_single']=s['证据ID']; p['_publication']=s['发表状态']

    source_proteins = fasta(P/'references/processed/reference_source_proteins.faa')
    refs = json.loads((P/'references/processed/references.json').read_text())
    systems = {}
    seqs = {}
    membership = []
    missing = []

    def add_system(sid, label, source, members, center='', fids=(), note='', binding='S2 同体系参考序列', clade=''):
        if sid in systems: return
        systems[sid]=dict(id=sid,label=label,source=source,center=center,fids=list(fids),note=note,binding=binding,clade=clade)
        seen_members=set()
        for member in list(members)+EXTRA_HELPERS.get(source,[]):
            role, selector, kind=member[:3]
            g=select(member[3] if len(member)>3 else source,selector)
            member_key=(g['accession'],g['locus'],kind)
            if member_key in seen_members:continue
            seen_members.add(member_key)
            qid='SEQ_'+sha(g['accession']+'|'+g['locus']+'|'+sha(g['sequence']))[:14]
            if qid not in seqs:
                seqs[qid]={**g,'id':qid,'systems':set(),'roles':set(),'types':set(),'fids':set(),'ref_ids':set(),'domains':set(),'sites':set()}
            q=seqs[qid]; q['systems'].add(sid); q['roles'].add(role); q['types'].add(component_type(role,kind,g['product']));q['fids'].update(fids)
            membership.append(dict(system=sid,sequence=qid,role=role,kind=kind,scope=binding))

    seed_aliases={
      'NiFe_Ecoli_Hyd1':['Ecoli_K12_HyaAB'], 'NiFe_Ecoli_Hyd2':['Ecoli_K12_HybOCAB'],
      'NiFe_FHL1_Hyc':['Ecoli_FHL1'], 'NiFe_HucLMS':['Msmegmatis_mc2155_Huc'],
      'NiFe_FrhABG':[], 'NiFe_MvhHdr':['M_marburgensis_MvhHdr'], 'NiFe_EchABCDEF':['Mbarkeri_Ech'],
      'NiFe_Pf_Mbh':['Pfur_Mbh'], 'FeFe_TmHydABC':['Thermotoga_maritima_HydABC'],
      'FeFe_SwHyd1ABC':['Syntrophomonas_wolfei_Hyd1ABC'], 'Fe_Hmd':['M_marburgensis_Hmd'],
      'NiFe_periplasmic_Hys':['DvH_HysAB'], 'NiFe_HynABC':['Wsuccinogenes_HydABC'],
      'NiFe_Hhy1h':['M_smegmatis_Hhy'], 'NiFe_VhoABC':['Mbarkeri_Fusaro_Vht'],
      'NiFe_HupSL':[], 'NiFe_H2_sensor_2b':['Rhodobacter capsulatus|HupUV–HupT/HupR'],
      'NiFe_HyhNADPH':['Pyrococcus_furiosus_SHI'], 'NiFe_HoxFUYH':['C_necator_SH'],
      'NiFe_Tonn_Fdh2Mfh2Mnh2':['Tonn_Fdh2Mfh2Mnh2'], 'NiFe_Coo_CODH':['Chydro_Coo'],
      'FeFe_CpI_monomer':['Clostridium_pasteurianum_W5_CpI'], 'FeFe_CrHydA1':['Chlamydomonas_reinhardtii_HydA1'],
      'FeFe_HDCR_Aw':['Acetobacterium_woodii_HDCR'], 'FeFe_HytABCDE':['Clostridium_autoethanogenum_Hyt_FdhA']}
    alias={}
    seed = json.loads((P/'config/systems.json').read_text())
    for s in seed:
        sid=s['structure_id']; rr=[r for r in refs if r['structure_id']==sid]
        members=[(r['role'],r['locus_tag'],r['kind']) for r in rr]
        for r in refs:
            if r['kind']=='maturation' and r['sample']==s['sample']:
                name=r['role'].removeprefix(s['sample']+'_')
                if s['center']=='Fe' or (s['center']=='FeFe') != (name.lower() in {'hyde','hydf','hydg'}):
                    continue
                members.append((r['role'].removeprefix(s['sample']+'_'),r['locus_tag'],'maturation'))
        add_system(sid,s['sample'],'references/raw/'+s['sample']+'.gb',members,s['center'],s['function_ids'],
            '参考基因组CDS，蛋白加工及实验构建形式逐篇另判；成熟基因按基因组归属，不能视为全部属于此操纵子。',clade=s.get('reference_clade',''))
        systems[sid]['small_features']=s.get('reference_small_subunit_topology','未单独核定')
        for k in seed_aliases[sid]:alias[k]=sid
        for r in rr:
            matches=[q for q in seqs.values() if q['accession']==r['source_accession'] and q['locus']==r['locus_tag']]
            assert len(matches)==1,(sid,r['locus_tag'])
            q=matches[0]
            assert sha(q['sequence'])==r['source_sequence_sha256'], r['reference_id']
            assert q['sequence']==source_proteins[r['reference_id']],r['reference_id']
            q['ref_ids'].add(r['reference_id'])
            for e in r.get('reference_role_evidence',{}).get('profile_evidence',[]):
                q['domains'].add(e.get('accession','')+':'+e.get('profile',''))
            if r.get('diagnostic_positions'):
                offset=int(r.get('protein_region_start',1))-1
                q['sites'].add('参考配位位点 '+','.join(str(int(i)+offset) for i in r['diagnostic_positions']))

    def extra(sid, keys, label, source, members, center, fids, **kw):
        add_system(sid,label,'references/raw/'+source,members,center,fids,**kw)
        for k in keys:alias[k]=sid

    # These are explicit locus maps already present in the project, not inferred by gene-name proximity.
    extrefs=json.loads((P/'references/processed/complex_extensions/references.json').read_text())
    ext=[('NiFe_Amobile_HydABCSL','A_mobile_HydABCSL','Acetomicrobium mobile DSM 13181','round22_function_truth/genomes/CP003198.gb',['Ni24']),
         ('FeFe_Sfructo_HndABCD','Solidesulfovibrio_fructosivorans_Hnd','Solidesulfovibrio fructosivorans','round22_function_truth/genomes/U07229.gb',['Fe03','Fe04','Fe08']),
         ('FeFe_Mthermo_HydABCDEF','Moorella_thermoacetica_HydABCDEF','Moorella thermoacetica DSM 521','round26_function_mechanism/CP012369.gb',['Fe08'])]
    for sid,key,label,src,fids in ext:
        rr=[r for r in extrefs if r['structure_id']==sid]
        extra(sid,[key],label,src,[(r['role'],r['locus_tag'],r['kind']) for r in rr],rr[0]['center'],fids,
              binding='S1 已有实验成员至原生CDS对应',note='沿用项目已有逐序列成员对应记录；不代表每篇论文的构建或加工形式均相同。')
        for r in rr:
            assert any(q['sequence']==r['sequence'] and q['protein_accession']==r['protein_accession'] for q in seqs.values()),r['reference_id']
    extra('NiFe_Cnecator_MBH',['Cnecator_H16_HoxKGZ','Cnecator_MBH'],'Cupriavidus necator H16 MBH','cnecator_H16_pHG1.gb',
          [('HoxK','PHG001','module'),('HoxG','PHG002','catalytic'),('HoxZ','PHG003','module')],'NiFe',['Ni02','Ni03','Ni05'])
    extra('NiFe_Pf_SHII',['Pyrococcus_furiosus_SHII'],'Pyrococcus furiosus SHII','pfuriosus_DSM3638.gb',
          [(r,'PF'+str(n),'catalytic' if n==1332 else 'module') for r,n in [('SHII_beta',1329),('SHII_gamma',1330),('SHII_delta',1331),('SHII_alpha',1332)]],'NiFe',['Ni12','Ni12R'])
    extra('NiFe_Tko_Hyh',['T_kodakarensis_Hyh'],'Thermococcus kodakarensis KOD1 Hyh','round60_function_validation/tkodakarensis/NC_006624.gb',
          [(r,'TK'+str(n),'catalytic' if r=='HyhL' else 'module') for r,n in [('HyhB',2072),('HyhG',2071),('HyhS',2070),('HyhL',2069)]],'NiFe',['Ni12','Ni12R'])
    extra('NiFe_Tko_Mbh',['Tko_Mbh'],'Thermococcus kodakarensis KOD1 Mbh','round60_function_validation/tkodakarensis/NC_006624.gb',
          [('Mbh'+chr(65+i),'TK'+str(2080+i),'catalytic' if i==11 else 'module') for i in range(14)]+[('Fd1','TK1694','module'),('Fd2','TK1087','module'),('Fd3','TK2012','module')],'NiFe',['Ni18'])
    extra('NiFe_Synechocystis_Hox',['Synechocystis_Hox'],'Synechocystis sp. PCC 6803 Hox','round66_split_sources/synechocystis/NC_000911.1.gb',
          [(r,l,'catalytic' if r=='HoxH' else 'module') for r,l in [('HoxE','sll1220'),('HoxF','sll1221'),('HoxU','sll1223'),('HoxY','sll1224'),('HoxH','sll1226')]],'NiFe',['Ni14','Ni24'],
          note='5个Hox亚基；转录区插入基因不自动视为催化亚基。预印本的实验构建未逐位核对。')
    extra('NiFe_Savermitilis_Hhy',['S_avermitilis_Hhy'],'Streptomyces avermitilis MA-4680 Hhy','round93_independent_sources/next_screen/S_avermitilis_hhySL/BA000030.4.gb',
          [('HhyL','BAC75078.1','catalytic'),('HhyS','BAC75079.1','module')],'NiFe',['Ni06'])
    extra('NiFe_HpSH',['H_pseudoflava_HpSH'],'Hydrogenophaga pseudoflava DSM 1084 HpSH','round60_function_validation/hydrogenophaga/CP037867.gb',
          [(r,'HPF_'+str(n).zfill(5),'catalytic' if r=='HoxH' else 'module') for r,n in [('HoxF',1495),('HoxU',1500),('HoxY',1505),('HoxH',1510)]],'NiFe',['Ni14','Ni12R'],
          binding='S1 已有实验成员至原生CDS对应',note='作者人工操纵子表达；原生CDS与既有对应记录连接，不把人工排列当原生操纵子。')
    extra('NiFe_Hthermo_Hox',['H_thermoluteolus_SH'],'Hydrogenophilus thermoluteolus TH-1 Hox','round22_function_truth/genomes/MF153396.gb',
          [(r,a,'catalytic' if r=='HoxH' else 'module') for r,a in [('HoxF','AST15075.1'),('HoxU','AST15076.1'),('HoxY','AST15077.1'),('HoxH','AST15078.1')]],'NiFe',['Ni14'])
    extra('NiFe_Mmazei_Ech',['Mmazei_Ech'],'Methanosarcina mazei Gö1 Ech','round71_group4_sources/Mmazei_Ech/NC_003901.1.gb',
          [('Ech'+chr(65+i),'MM_'+str(2320+i),'catalytic' if i==4 else 'module') for i in range(6)],'NiFe',['Ni19','Ni23'])
    extra('NiFe_Tteng_Ech',['Tteng_Ech'],'Caldanaerobacter tengcongensis MB4 Ech','round67_ech_sources/tengcongensis_NC_003869.1.gb',
          [('Ech'+chr(65+i),'TTE0'+str(123+i),'catalytic' if i==4 else 'module') for i in range(6)]+[('Hyp'+letter,'TTE0'+str(129+i),'maturation') for i,letter in enumerate('ABFCDE')],'NiFe',['Ni23'],
          note='Ech及邻近Hyp基因对应来自既有位点记录；参考序列不证明每个基因均在该实验中被测定。')
    extra('FeFe_Tkivui_HDCR',['Thermoanaerobacter_kivui_HDCR'],'Thermoanaerobacter kivui DSM 2030 HDCR','round22_function_truth/genomes/CP009170.gb',
          [('HydA2','TKV_c19960','catalytic'),('HycB4','TKV_c19970','module'),('HycB3','TKV_c19980','module'),('FdhF','TKV_c19990','module'),('FdhD','TKV_c19950','maturation')],'FeFe',['Fe06'])
    extra('FeFe_Tkivui_HydABC',['Thermoanaerobacter_kivui_HydABC'],'Thermoanaerobacter kivui DSM 2030 HydABC','round22_function_truth/genomes/CP009170.gb',
          [('HydA','TKV_c19580','catalytic'),('HydB','TKV_c19590','module'),('HydC','TKV_c19600','module')],'FeFe',['Fe07'],note='受体特异性按每条实测记录；基因名本身不能区分NAD与NADP。')
    extra('NiFe_Tkivui_Ech2',['Tkivui_Ech2'],'Thermoanaerobacter kivui DSM 2030 Ech2','round22_function_truth/genomes/CP009170.gb',
          [(r,l,'catalytic' if r=='Ech2E' else 'module') for r,l in [('Ech2A1','TKV_c19720'),('Ech2A2','TKV_c19710'),('Ech2B','TKV_c19690'),('Ech2C','TKV_c19700'),('Ech2D','TKV_c19750'),('Ech2E','TKV_c19740'),('Ech2F','TKV_c19680'),('HycB2','TKV_c19730')]],'NiFe',['Ni23'],note='纯化Ech2对应；2019总膜Na运输不能单独归因于此Ech2。')
    extra('FeFe_Sacid_HydAB',['Syntrophus_aciditrophicus_HydAB'],'Syntrophus aciditrophicus SB HydAB','round26_function_mechanism/CP000252.gb',
          [('HydA','SYN_01370','catalytic'),('HydB','SYN_01369','module')],'FeFe',['Fe05'],note='原生蛋白序列；N端His实验构建另有既有对应记录。')
    extra('NiFe_Patro_FHL2',['Patrosepticum_FHL2'],'Pectobacterium atrosepticum SCRI1043 FHL-2','round91_Hyf_sources/NC_004547.2_Pa_SCRI1043.gb',
          [(r,'ECA'+str(n),'catalytic' if r=='HyfG' else 'module') for r,n in [('HyfA',1247),('HyfB',1246),('HyfC',1245),('HyfD',1244),('HyfE',1243),('HyfF',1242),('HyfG',1241),('HyfH',1240),('HyfI',1239),('FdhF',1250),('HydN',1249)]]+[('HycI','ECA1237','maturation'),('HycH','ECA1238','maturation')],'NiFe',['Ni22'])
    extra('NiFe_Tguam_FHL2',['Tguamensis_FHL2'],'Trabulsiella guamensis ATCC 49490 FHL-2','round92_FHL_sources/GCF_000734965.1_Trabulsiella_guamensis_ATCC49490.gbff',
          [(r,l,'catalytic' if r=='HyfG' else 'module') for r,l in [('HyfI','GTGU_RS12540'),('HyfH','GTGU_RS12545'),('HyfG','GTGU_RS12550'),('HyfF','GTGU_RS12555'),('HyfE','GTGU_RS12560'),('HyfD','GTGU_RS12565'),('HyfC','GTGU_RS12570'),('HyfB','GTGU_RS12575'),('HyfA','GTGU_RS12580'),('FdhF','GTGU_RS14390')]],'NiFe',['Ni22'],note='HyfA–I原生CDS与FdH伙伴；WGS不同contig保持各自坐标，不人工拼接为单一操纵子。')

    # Additional source sequences are useful as references even when the assayed strain is different.
    extra('NiFe_Tonn_Frh_external',['T_onnurineus_Frh_external'],'Thermococcus onnurineus NA1 Frh及外部伙伴','tonnurineus_NA1.gb',
          [(r,l,'catalytic' if r=='FrhA' else 'module') for r,l in [('FrhA','TON_1559'),('FrhG','TON_1560'),('FrhB','TON_1561'),('Fdh3B','TON_0542'),('TON_0702','TON_0702'),('TON_1376','TON_1376')]],'NiFe',['Ni21'],note='不同论文使用不同外部伙伴；本列表为参考蛋白并集，并非单个复合体。')
    extra('FeFe_Swolfei_Hyd2',['Swolfei_FeFeHyd2'],'Syntrophomonas wolfei Göttingen Hyd2','swolfei_Goettingen.gb',
          [('HydA2','Swol_1925','catalytic'),('HydB2','Swol_1926','module'),('HydC2','Swol_1927','module')],'FeFe',['Ni04'],
          note='保留膜相关FeFe体系作为Ni04跨架构候选证据；不能由邻域或表达直接断言醌驱动产氢。')
    extra('FeFe_Ralbus_HydABC',['Ruminococcus_albus_7_HydABC'],'Ruminococcus albus 7 HydABC','round22_function_truth/genomes/CP002403.gb',
          [('HydA','Rumal_2964','catalytic'),('HydB','Rumal_2965','module'),('HydC','Rumal_2966','module')],'FeFe',['Fe03','Fe04'],
          note='既有功能参考位点；HydB全长序列与assayed_modules.faa已存参考一致。')
    extra('FeFe_Ralbus_HydA2',['Ruminococcus_albus_7_HydA2'],'Ruminococcus albus 7 HydA2','round22_function_truth/genomes/CP002404.gb',
          [('HydA2','Rumal_3407','catalytic'),
           ('HydE','Rumal_0545','maturation','references/raw/round22_function_truth/genomes/CP002403.gb'),
           ('HydF','Rumal_0543','maturation','references/raw/round22_function_truth/genomes/CP002403.gb')],
          'FeFe',['Fe01','Fe05'],
          note='质粒编码催化亚基Rumal_3407，附同菌株染色体HydE/F成熟基因；小亚基及生理电子伙伴未单独核定。')
    extra('FeFe_Awoodii_HydABC',['Acetobacterium_woodii_HydABC'],'Acetobacterium woodii DSM 1030 HydABC','awoodii_DSM1030.gb',
          [('HydA','Awo_c26970','catalytic'),('HydB','Awo_c26980','module'),('HydC','Awo_c27010','module')],'FeFe',['Fe03','Fe04'],
          note='列出已对应HydABC核心；部分来源称HydABCD，本表未强行给额外邻近基因赋予HydD角色。')
    extra('FeFe_Cbei_HDCR',['Clostridium_beijerinckii_NCIMB8052_HDCR'],'Clostridium beijerinckii NCIMB 8052 HDCR','round87_FdhHyd_sources/CP000721.1_genbank.txt',
          [('HycB_Cbei3795','Cbei_3795','module'),('HydA','Cbei_3796','catalytic'),('HycB_Cbei3797','Cbei_3797','module'),('FdhF','Cbei_3801','module')],'FeFe',['Fe06'],
          note='采用既有HDCR成员与GenBank位点；HycB以位点区分，不猜测不同论文的1/2命名次序。')
    for key,sid,role in [('Ra_HydB','FeFe_Ralbus_HydABC','HydB'),('Aw_HydB','FeFe_Awoodii_HydABC','HydB')]:
        assayed=fasta(P/'references/raw/round26_function_mechanism/assayed_modules.faa')[key]
        assert any(seqs[m['sequence']]['sequence']==assayed for m in membership if m['system']==sid and m['role']==role),key
    for spec in EXTRA_SYSTEMS:
        extra(**spec)
    # Mature hydrogenase sources and genome assemblies remain separate from experimental protein constructs.
    genome_rows=read_tsv(P/'references/processed/evidence_genome_accessions_20261110.tsv')
    genome_lookup={r['证据ID']:r['基因组数据库编号'] for r in genome_rows}
    original=read_tsv(P/'references/processed/function_evidence_expansion_20261110.tsv')
    genome_hints=defaultdict(set)
    for f in original:
        for i in range(1,4):
            k=citation_key(f.get('文献'+str(i),'')); acc=genome_lookup.get(f['功能ID']+'-E'+str(i),'')
            if acc:genome_hints[(f['功能ID'],k)].add(acc)

    # Merge repeat entries for a paper/system; distinct functions and distinct measured systems stay separate.
    records={}
    for p in papers:
        fid=p['功能ID']; key=citation_key(p['DOI或文献键']); sk=p['体系去重键'] or 'unresolved'
        sid=alias.get(sk,sk if sk in systems else '')
        sid=EXTRA_BINDINGS.get(fid+'|'+key,sid)
        if not sid and p.get('_single'):
            txt=p['物种或来源']
            if 'Synechocystis' in txt:sid='NiFe_Synechocystis_Hox'
            elif 'Cupriavidus' in txt or 'Ralstonia' in txt:sid='NiFe_Cnecator_MBH'
        k=(fid,key,sid or sk)
        if k not in records:records[k]={**p,'_sid':sid,'_key':key}
        else:
            old=records[k]
            if old.get('证据分级')=='CANDIDATE' and p.get('_source')=='evidence_six_20261112/added_studies.json' and p.get('证据分级') not in {'','CANDIDATE','UNRESOLVED'}:
                old['证据分级']=p['证据分级']
            for field in ['实际测得','尚未证明或限制','原文位置','获取层级','_source']:
                if p.get(field) and p[field] not in old.get(field,''):old[field]=old.get(field,'')+'；'+p[field]
    evidence=[]
    for (fid,key,sk),p in records.items():
        sid=p['_sid']; s=systems.get(sid); old=p['证据分级']; measured=p['实际测得']; limits=p['尚未证明或限制']
        if old in {'DUPLICATE','EXCLUDED_WRONG_ENZYME'}:
            level='X 排除或更正'; stance='更正/背景' if old=='DUPLICATE' else '不支持该类别';basis=old
        elif old=='CANDIDATE' or old=='UNRESOLVED' and not measured:
            level='D 候选或待确认';stance='待确认';basis='原记录未形成可定位的实验结论'
        elif re.search('敲除|缺失|突变|遗传|knockout|deletion',measured,re.I) and not re.search('纯化|重构|动力学|电化学|脂质体',measured):
            level='B 遗传或生理因果';stance='实验支持';basis='现有记录描述成员干预及生理/反应表型'
        elif re.search('纯化|重构|还原活性|酶活|反应速率|周转|动力学|电化学|脂质体|膜囊泡|膜泡|ATP合成|H₂.*还原|H2.*还原',measured):
            level='A 直接反应或生化测定';stance='实验支持';basis='现有记录描述反应、生化或膜体系实测；仅支持实际测得栏所列端点'
        elif old in {'FULL_DIRECT','FULL_LINKED','NONPEER_COMPLETE','CROSS_ARCH_DIRECT'}:
            level='B 遗传或生理因果' if re.search('生长|产氢|存活|表达|培养',measured) else 'C 结构或间接实验';stance='实验支持';basis='沿用已收录实验记录；实验类型见实际测得栏'
        else:
            level='C 结构或间接实验' if measured and old not in {'CANDIDATE','UNRESOLVED'} else 'D 候选或待确认'
            stance='部分支持' if level.startswith('C') else '待确认';basis='现有结构/组分/相关性记录，或尚不能确定实验端点'
        if not measured or p['物种或来源']=='未重新核实':level='D 候选或待确认';stance='待确认'
        if re.search('不是Huc|非氢酶本体|研究对象不是',limits) and old=='EXCLUDED_WRONG_ENZYME':sid='';s=None
        if fid=='Fe08' and ('Hnd' in p['酶或体系'] or 'HydAB'==p['酶或体系']):stance='反例或特异性对照'
        mids=[m for m in membership if m['system']==sid] if sid else []
        cats=[m['sequence'] for m in mids if m['kind']=='catalytic']
        small=[m['sequence'] for m in mids if '小亚基' in seqs[m['sequence']]['types']]
        aux=[m['sequence'] for m in mids if m['kind']=='maturation']
        other=[m['sequence'] for m in mids if m['sequence'] not in cats+small+aux]
        seq_ids=list(dict.fromkeys(cats+small+other+aux))
        links=[]
        for g in groups:
            if g['功能ID']==fid and key in [citation_key(x) for x in g['研究DOI'].split('；')]:links.append(g['证据组ID'])
        eid=fid+'-EV-'+sha(key+'|'+sk)[:10]
        genomic=sorted({seqs[x]['accession'] for x in seq_ids})
        hints=sorted(genome_hints.get((fid,key),set()))
        is_reference_strain = False
        if s:
            if sid in {'NiFe_Ecoli_Hyd1', 'NiFe_Ecoli_Hyd2'}:
                is_reference_strain = bool(re.search(r'MC4100|BW25113|Salmonella', p['物种或来源']))
            elif sid == 'NiFe_FHL1_Hyc':
                is_reference_strain = not bool(re.search(r'Escherichia|E\.?\s*coli', p['物种或来源']))
        binding='S3 近缘或不同菌株参考' if is_reference_strain else s['binding'] if s else 'S0 未绑定序列'
        if not seq_ids: binding='S0 未绑定序列'
        r={'证据ID':eid,'功能ID':fid,'功能类别':function_by_id[fid]['反应级功能'],'证据级别':level,'支持关系':stance,
           '反应核心归组':CORES[fid],
           '实际测得':measured or '现有资料尚未提取实验端点','限制及未测部分':limits,'微生物及实验来源':p['物种或来源'],
           '氢酶及被测体系':p['酶或体系'],'体系ID':sid or sk,'原体系去重键':p['体系去重键'],
           '蛋白参考活性中心':s['center'] if s else '未逐体系绑定','架构分支':s['clade'] if s else '',
           '序列对应级别':binding,'参考基因组或核苷酸记录':';'.join(genomic),'旧表组装线索_菌株待核对':';'.join(hints),
           '催化亚基序列ID':';'.join(cats),'小亚基序列ID':';'.join(small),'电子伙伴及相关亚基序列ID':';'.join(other),
           '成熟辅助基因序列ID':';'.join(aux),'已绑定蛋白数':len(seq_ids),'互补证据组ID':';'.join(links),
           '发表年份':p['发表年份'],'发表状态':p.get('_publication','预印本/学位论文，见来源' if re.search('biorxiv|博士|论文摘要|预印本',p['原始来源URL']+' '+limits) else '按原始来源记录'),
           'DOI或文献键':p['DOI或文献键'],'文献题名或主题':p['文献题名或主题_可缩写'],
           '原文位置':p['原文位置'],'原始来源URL':p['原始来源URL'],'读取层级':p['获取层级'],
           '分级依据':'基于现有实际测得记录整理初分级，未新增逐篇全文复核：'+basis,'原证据标签':old,'整理来源':p.get('_source',''),'体系序列限制':s['note'] if s else '未以物种相近或类名相同自动绑定蛋白。'}
        evidence.append(r)
    # Bare literature links are search leads, not measured evidence; retain them separately.
    candidates=[r for r in evidence if r['微生物及实验来源']=='未重新核实']
    evidence=[r for r in evidence if r['微生物及实验来源']!='未重新核实']
    evidence.sort(key=lambda r:(list(function_by_id).index(r['功能ID']),r['证据级别'],r['DOI或文献键'].lower(),r['体系ID']))
    assert len({r['证据ID'] for r in evidence})==len(evidence)
    links=[]
    for e in evidence:
        for m in membership:
            if m['system']==e['体系ID']:
                q=seqs[m['sequence']]
                links.append({'证据ID':e['证据ID'],'功能ID':e['功能ID'],'体系ID':m['system'],'序列ID':m['sequence'],
                    '亚基或基因名':m['role'],'组件类别':component_type(m['role'],m['kind'],q['product']),
                    '序列对应级别':e['序列对应级别'],'功能证据级别':e['证据级别'],
                    '成员功能范围':'成熟基因支持成熟，不决定直接电子受体' if m['kind']=='maturation' else '体系参考成员；不表示该论文逐个测定所有亚基'})
                q['fids'].add(e['功能ID'])
    seq_rows=[]
    for qid,q in sorted(seqs.items()):
        assert q['sequence'] and re.fullmatch('[A-Z*]+',q['sequence'])
        assert q['cds'] and re.fullmatch('[ACGTRYMKBDHVNSW]+',q['cds'])
        seq_rows.append({'序列ID':qid,'体系ID':';'.join(sorted(q['systems'])),'功能ID':';'.join(sorted(q['fids'])),
            '亚基或基因名':';'.join(sorted(q['roles'])),'组件类别':';'.join(sorted(q['types'])),
            '微生物':q['organism'],'参考菌株':q['strain'],'核苷酸登录号':q['accession'],'locus_tag':q['locus'],
            'gene':q['gene'],'蛋白登录号':q['protein_accession'],'蛋白产品原注释':q['product'],'氨基酸长度_aa':len(q['sequence']),
            '完整CDS蛋白序列':q['sequence'],'CDS长度_nt':len(q['cds']),'CDS核酸序列':q['cds'],
            '起点_1based':q['start'],'终点_1based':q['end'],'链方向':q['strand'],'GenBank位置':q['location'],
            'Cys位置_1based':';'.join(str(i+1) for i,x in enumerate(q['sequence']) if x=='C'),
            'Sec位置_1based':';'.join(str(i+1) for i,x in enumerate(q['sequence']) if x=='U'),
            '已有结构域特征':';'.join(sorted(q['domains'])),'已有参考配位位点':';'.join(sorted(q['sites'])),
            '蛋白SHA256':sha(q['sequence']),'CDS_SHA256':sha(q['cds']),
            '序列来源URL':'https://www.ncbi.nlm.nih.gov/protein/'+q['protein_accession'],
            '基因组来源URL':'https://www.ncbi.nlm.nih.gov/nuccore/'+q['accession'],
            '序列提取说明':q['cds_status'],'本地来源':str(P/q['source'])})
    system_rows=[]
    for sid,s in systems.items():
        mm=[m for m in membership if m['system']==sid]
        q0=seqs[mm[0]['sequence']]
        ee=[e for e in evidence if e['体系ID']==sid]
        system_rows.append({'体系ID':sid,'微生物及体系':s['label'],'参考菌株':q0['strain'],'活性中心':s['center'],
            '功能ID':';'.join(sorted(set(s['fids'])|{e['功能ID'] for e in ee})),
            '核苷酸记录':';'.join(sorted({seqs[m['sequence']]['accession'] for m in mm})),
            '参考成员':'；'.join(m['role']+'='+seqs[m['sequence']]['locus'] for m in mm),
            '蛋白数':len(mm),'关联证据行数':len(ee),'序列对应级别':s['binding'],
            '小亚基FeS组织_已有参考':s.get('small_features','未单独核定；Cys位置不能独自证明簇类型或占据率'),'范围与限制':s['note'],
            '本地GenBank文件':str(P/s['source']),'记录链接':'https://www.ncbi.nlm.nih.gov/nuccore/'+q0['accession']})
    # Retain unbound systems explicitly, with their accession clues and original experimental organism.
    bound=set(systems)
    for key in sorted({e['体系ID'] for e in evidence}-bound):
        ee=[e for e in evidence if e['体系ID']==key]
        system_rows.append({'体系ID':key,'微生物及体系':'；'.join(dict.fromkeys(e['微生物及实验来源']+' / '+e['氢酶及被测体系'] for e in ee)),
            '参考菌株':'未绑定','活性中心':'未逐体系核定','功能ID':';'.join(sorted({e['功能ID'] for e in ee})),
            '核苷酸记录':'；'.join(sorted({e['旧表组装线索_菌株待核对'] for e in ee if e['旧表组装线索_菌株待核对']})),
            '参考成员':'未逐亚基绑定','蛋白数':0,'关联证据行数':len(ee),'序列对应级别':'S0 未绑定序列',
            '小亚基FeS组织_已有参考':'未绑定','范围与限制':'保留证据；组装线索不等于被测成员和序列已经对应。','本地GenBank文件':'','记录链接':''})
    summaries=[]
    for f in funcs:
        fid=f['功能ID']; ee=[e for e in evidence if e['功能ID']==fid]
        qq={l['序列ID'] for l in links if l['功能ID']==fid}
        summaries.append({'功能ID':fid,'功能类别':f['反应级功能'],'反应核心归组':CORES[fid],'原代表体系':f['代表系统'],'类别活性中心':f['活性中心'],
            '证据记录数':len(ee),'A直接测定':sum(e['证据级别'].startswith('A') for e in ee),
            'B遗传生理':sum(e['证据级别'].startswith('B') for e in ee),'C间接实验':sum(e['证据级别'].startswith('C') for e in ee),
            'D候选待确认':sum(e['证据级别'].startswith('D') for e in ee),'X排除更正':sum(e['证据级别'].startswith('X') for e in ee),
            '不同文献键数':len({citation_key(e['DOI或文献键']) for e in ee}),
            '已关联蛋白的证据行数':sum(e['已绑定蛋白数']>0 for e in ee),'蛋白序列记录数':len(qq),
            '候选资料链接数':sum(c['功能ID']==fid for c in candidates),
            '成员与结构特征_类别参考':f['结构及电子传递亚基'],'特异电子模块_类别参考':f['特异模块'],
            '通用成熟基因_类别参考':f['通用成熟基因'],
            '分类使用范围':'按实际测得栏赋予反应结论；类别中的受体、离子和生理限定未测时保留候选，序列相似不自动等于全部功能。'})
    outputs={'evidence.tsv':evidence,'function_summary.tsv':summaries,'systems_genomes.tsv':system_rows,
             'sequences.tsv':seq_rows,'evidence_sequence_links.tsv':links,'candidate_sources.tsv':candidates,
             'complementary_groups.tsv':groups}
    for name,rows in outputs.items():write_tsv(name,rows)
    for filename,field in [('proteins.faa','sequence'),('cds.fna','cds')]:
        with open(OUT/filename,'w') as f:
            for qid,q in sorted(seqs.items()):
                f.write('>'+qid+' '+q['protein_accession']+' '+q['accession']+' '+q['locus']+'\n')
                seq=q[field]
                for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')
    clusters=0
    with open(OUT/'gene_clusters.fna','w') as f:
        for sid,s in systems.items():
            genes=[seqs[m['sequence']] for m in membership if m['system']==sid and m['kind'] not in {'maturation','context'}]
            for acc in sorted({g['accession'] for g in genes}):
                gs=[g for g in genes if g['accession']==acc]; lo=min(g['start'] for g in gs);hi=max(g['end'] for g in gs)
                if hi-lo>100000:continue
                rec=next(r for r in genbank(gs[0]['source']) if r['accession']==acc)
                dna=rec['dna'][max(0,lo-501):min(len(rec['dna']),hi+500)]
                f.write('>'+sid+'__'+acc+' '+acc+':'+str(max(1,lo-500))+'..'+str(min(len(rec['dna']),hi+500))+' genomic_region_not_transcription_proof\n')
                for i in range(0,len(dna),80):f.write(dna[i:i+80]+'\n')
                clusters+=1
    # Referential integrity and whole-source sequence correspondence are the important export checks.
    assert len(summaries)==34
    assert sum(r['证据记录数'] for r in summaries)==len(evidence)
    assert all(l['序列ID'] in seqs for l in links)
    assert len(fasta(OUT/'gene_clusters.fna'))==clusters
    assert all(l['证据ID'] in {e['证据ID'] for e in evidence} for l in links)
    for r in summaries:assert sum(r[k] for k in ['A直接测定','B遗传生理','C间接实验','D候选待确认','X排除更正'])==r['证据记录数']
    payload={'sheets':outputs,'groups':groups,'stats':{'functions':34,'reaction_core_groups':len(CORE_GROUPS),'evidence_rows':len(evidence),
        'unique_citation_keys':len({citation_key(e['DOI或文献键']) for e in evidence}),
        'candidate_source_rows':len(candidates),
        'systems':len(system_rows),'systems_with_sequences':len(systems),'protein_sequences':len(seqs),
        'unique_protein_sha256':len({sha(q['sequence']) for q in seqs.values()}),
        'evidence_with_sequences':sum(e['已绑定蛋白数']>0 for e in evidence),'sequence_links':len(links),
        'gene_clusters':clusters,'grades':dict(Counter(e['证据级别'] for e in evidence)),
        'sequence_binding_counts':dict(Counter(e['序列对应级别'] for e in evidence)),
        'component_counts':dict(Counter(t for q in seqs.values() for t in q['types'])),
        'proteins_with_existing_domain_features':sum(bool(q['domains']) for q in seqs.values()),
        'proteins_with_existing_ligand_positions':sum(bool(q['sites']) for q in seqs.values())}}
    (WORK/'catalog.json').write_text(json.dumps(payload,ensure_ascii=False))
    (OUT/'summary.json').write_text(json.dumps(payload['stats'],ensure_ascii=False,indent=2)+'\n')
    stats=payload['stats']
    guide=[
        ('整理日期','2026-11-12'),
        ('范围',f"34项既有功能、{len(CORE_GROUPS)}个反应核心导航组；{len(evidence)}条证据/对照记录；{len(candidates)}条未核实资料链接另表。保留原功能ID和完整限定，不再用每类3条最强证据作为收录门槛。"),
        ('每行的含义','证据明细每行=一个功能类别×一个文献键×一个被测体系。相同论文支持不同功能可出现多行；不能把行数当作独立论文或独立微生物数。同文同体系的重复条目已合并。'),
        ('A 直接测定','已有记录描述直接反应、生化、纯化或膜体系实测。A只支持“实际测得”栏端点，不表示该类别所有限定都已验证。'),
        ('B 遗传或生理','成员干预、遗传或生理实验支持；生长/表型不单独证明直接电子受体。'),
        ('C 结构或间接','结构、组分、光谱、表达或其他间接实验支持。保留为证据，不强行升级为直接反应。'),
        ('D 候选或待确认','有具体候选体系，但端点或功能尚不明确；不作确定功能正标签。只有链接、没有实验记录的资料移至候选资料表。'),
        ('X 排除或更正','重复来源、错误酶对象等负面记录，保留可追溯性，不计正证据。'),
        ('分级适用范围','本次依据现有实际测得描述整理初分级，保留原证据标签、全文/摘要/检索读取层级和局限；没有重新阅读全文，不能将全部记录解释为已独立复核。'),
        ('序列S1','项目已有实验成员—原生CDS对应；仍须区分实验标签、加工形式及不同论文构建。'),
        ('序列S2','同体系参考菌株的基因组CDS；未逐篇确认实验构建与参考序列完全相同。'),
        ('序列S3','不同菌株或近缘参考；不作为被测蛋白精确序列的证明。'),
        ('序列S0','未绑定成员。空白不代表基因不存在，未把相似物种或架构自动当作被测体系。'),
        ('序列文件',f"{len(seqs)}条蛋白序列记录，按核苷酸登录号+locus_tag+蛋白哈希去重；全局不同AA哈希数={stats['unique_protein_sha256']}。proteins.faa与cds.fna同ID逐条对应，序列表也保留完整字符串。"),
        ('亚基与辅助基因','催化亚基、小亚基/近端电子传递、电子伙伴、膜模块、成熟及调控因子分别标注。成熟因子来自参考基因组，可能被多个氢酶共享，不能视作功能特异标志或全部同一操纵子。'),
        ('特征','保留原产品注释、长度、Cys/Sec位置、已有结构域及参考配位位点、坐标、链、完整AA/CDS和SHA256。Cys数量或间距不能独自证明FeS簇类型/占据率、电子受体或反应方向。'),
        ('基因组与操纵子','体系表给核苷酸记录、菌株、原位点和本地GenBank路径；旧组装线索与实绑记录分开。gene_clusters.fna为参考成员范围两侧各500nt的基因组区域（跨度>100kb不输出），不是转录操纵子证明。'),
        ('互补研究','互补证据组单独列出；每篇论文仍保留各自证据行，不能将互补组再次累加进明细行数。'),
        ('使用方式','先按功能ID/反应核心筛选证据；看实际测得、支持关系、读取层级及限制；再用证据ID→证据序列对应→序列ID取得FASTA。反例/特异性对照不作为该类正标签。'),
        ('来源','主要来自20261111用户包336行、64条定向复核更新及已有后继研究/互补表；序列来自项目已缓存GenBank与既有成员/全长序列映射。本次未改原证据库、分类器或模型。'),
        ('校验','断言校验34类与等级汇总、证据/序列ID唯一、关联无悬空、既有参考全长蛋白/哈希一致、CDS坐标提取及FASTA记录对应；不以导出成功替代功能准确率评测。'),
    ]
    payload['guide']=[{'项目':k,'说明':v} for k,v in guide]
    (WORK/'catalog.json').write_text(json.dumps(payload,ensure_ascii=False))
    (OUT/'README.txt').write_text('\n\n'.join(k+'\n'+v for k,v in guide)+'\n',encoding='utf-8')
    print(json.dumps(payload['stats'],ensure_ascii=False),flush=True)


if __name__=='__main__':
    build()
