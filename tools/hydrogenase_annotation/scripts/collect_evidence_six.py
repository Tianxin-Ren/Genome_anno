#!/usr/bin/env python3
"""Cache accession-pinned public sources for the six-evidence reference successor."""
import hashlib
import json
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from build_evidence_six import CURATION, P

RAW=P/'references/raw/evidence_six_20261112'


def fetch(url,path,expected):
    if path.exists() and expected in path.read_bytes():return path.read_bytes()
    req=urllib.request.Request(url,headers={'User-Agent':'HydrogenaseReferenceCuration/1.0'})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=40) as response:data=response.read()
            if expected not in data:raise ValueError('Unexpected source content: '+url)
            path.with_suffix(path.suffix+'.part').write_bytes(data)
            path.with_suffix(path.suffix+'.part').replace(path)
            return data
        except Exception:
            if attempt==2:raise
            time.sleep(attempt+1)


def collect_study(plan):
    key=plan['source']
    query=('DOI:"'+key+'"' if key.startswith('10.') else 'PMCID:'+key if key.startswith('PMC') else 'EXT_ID:'+key+' AND SRC:MED')
    path=RAW/('study_'+hashlib.sha256(key.encode()).hexdigest()[:14]+'.json')
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':query,'format':'json','resultType':'core','pageSize':3})
    content=json.loads(fetch(url,path,b'resultList'))
    results=content['resultList']['result']
    if not results and key.startswith('10.'):
        cross=RAW/('crossref_'+hashlib.sha256(key.encode()).hexdigest()[:14]+'.json')
        metadata=json.loads(fetch('https://api.crossref.org/works/'+urllib.parse.quote(key,safe=''),cross,b'"message"'))['message']
        if metadata['DOI'].lower()!=key.lower():raise ValueError('Crossref DOI mismatch')
        date=metadata.get('published',metadata.get('issued',{})).get('date-parts',[['']])[0]
        results=[{'doi':metadata['DOI'],'id':'','title':metadata['title'][0],'pubYear':str(date[0]),'isOpenAccess':'N','metadata_source':'Crossref'}]
    if len(results)!=1:raise ValueError((key,'not a unique source',len(results)))
    r=results[0];doi=r.get('doi','');pmid=r.get('id','');pmcid=r.get('pmcid','')
    if key.startswith('10.') and doi.lower()!=key.lower():raise ValueError((key,doi))
    level='已核读出版社摘要/文中片段；Crossref元数据' if r.get('metadata_source')=='Crossref' else '摘要及原始数据库元数据'
    if pmcid and r.get('isOpenAccess')=='Y':
        try:
            fetch(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML',RAW/(pmcid+'.xml'),b'<article')
            level+='；开放全文已缓存，具体核读范围见定位'
        except Exception as e:print('XML unavailable',pmcid,str(e)[:100],flush=True)
    out=[]
    for fid in plan['functions']:
        out.append({'功能ID':fid,'证据分级':plan.get('legacy_grade','PARTIAL'),'发表年份':r.get('pubYear',''),
            '物种或来源':plan['organism'],'酶或体系':plan['enzyme'],'体系去重键':plan['system'],
            'DOI或文献键':doi or 'PMID:'+pmid,'文献题名或主题_可缩写':r['title'],
            '实际测得':plan['measured'],'尚未证明或限制':plan['limits'],
            '独立性及去重':'本次新增来源，按DOI/PMID与同功能去重；不是同一文章的不同入口。',
            '原文位置':plan.get('locator','原始摘要；完整实验范围待正文细化'),
            '获取层级':level,'原始来源URL':'https://doi.org/'+doi if doi else 'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/'})
    print('STUDY',key,doi or pmid,r['title'][:85],flush=True)
    return out


def papers():
    RAW.mkdir(parents=True,exist_ok=True)
    plans=json.loads((CURATION/'study_plan.json').read_text())
    rows=[];errors=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs=[(p,pool.submit(collect_study,p)) for p in plans]
        for p,f in jobs:
            try:rows.extend(f.result())
            except Exception as e:errors.append({'source':p['source'],'error':str(e)});print('ERROR',p['source'],str(e),flush=True)
    (CURATION/'added_studies.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
    (RAW/'study_fetch_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2)+'\n')
    print('COLLECTED',len(rows),'function-study rows;',len(errors),'errors',flush=True)


def genomes():
    RAW.mkdir(parents=True,exist_ok=True)
    plans=json.loads((CURATION/'genome_plan.json').read_text())
    def one(plan):
        acc=plan['accession'];path=RAW/(acc+'.gb')
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'nuccore','id':acc,'rettype':'gbwithparts','retmode':'text'})
        data=fetch(url,path,b'\nORIGIN')
        text=data.decode();version=next(x.split()[1] for x in text.splitlines() if x.startswith('VERSION'))
        if version.split('.')[0]!=acc.split('.')[0]:raise ValueError('Genome accession mismatch')
        if plan['organism_check'].lower() not in text.lower():raise ValueError('Genome organism mismatch')
        row={**plan,'version':version,'path':str(path.relative_to(P)),'url':url,'sha256':hashlib.sha256(data).hexdigest()}
        print('GENOME',acc,version,len(data),flush=True);return row
    rows=[];errors=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs=[(p,pool.submit(one,p)) for p in plans]
        for p,f in jobs:
            try:rows.append(f.result())
            except Exception as e:errors.append({'accession':p['accession'],'error':str(e)});print('ERROR',p['accession'],str(e),flush=True)
    (RAW/'genome_manifest.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
    (RAW/'genome_errors.json').write_text(json.dumps(errors,ensure_ascii=False,indent=2)+'\n')
    print('GENOMES',len(rows),'errors',len(errors),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('mode',nargs='?',default='papers',choices=['papers','genomes'])
    args=parser.parse_args()
    genomes() if args.mode=='genomes' else papers()
