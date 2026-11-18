"""Use existing collinear-domain rules; never promote arbitrary partial hits to full domains."""
import hashlib
import json
import math
from collections import defaultdict

from audit_reference_domains import model_evidence


def summarize_domains(rows, policy):
    groups=defaultdict(list)
    for d in rows:
        if (not 1<=d['model_start']<=d['model_end']<=d['model_length'] or
            not 1<=d['query_start']<=d['query_end']<=d['query_length'] or
            any(not math.isfinite(d[k]) for k in ['sequence_bitscore','domain_bitscore','independent_evalue']) or
            d['independent_evalue']<0):
            raise ValueError('Invalid domain coordinates or scores')
        # Different model versions/libraries are never stitched into one coverage estimate.
        groups[(d['序列ID'],d['模型集'],d['model_sha256'],d['accession'],d['profile'])].append(d)
    result=defaultdict(list)
    for (qid,library,digest,accession,profile),parts in sorted(groups.items()):
        if len({(p['model_length'],p['query_length'],p['sequence_ga'],p['domain_ga']) for p in parts})!=1:
            raise ValueError('Inconsistent lengths or gathering thresholds in one model group')
        unique={tuple(p[k] for k in ['model_start','model_end','query_start','query_end','domain_bitscore']):p for p in parts}
        for e in model_evidence(list(unique.values()),policy):
            chain=e.pop('components')
            models=[[p['model_start'],p['model_end']] for p in chain]
            queries=[[p['query_start'],p['query_end']] for p in chain]
            assert all(a[1]<b[0] for spans in [models,queries] for a,b in zip(spans,spans[1:]))
            assert 0<e['model_coverage']<=1
            result[qid].append({**e,'模型集':library,'模型完整编号':accession,'model_sha256':digest,
                'status':'profile_covered_by_single_segment' if e['supported'] and len(chain)==1 else
                    'profile_covered_by_collinear_segments' if e['supported'] else 'partial_profile_match',
                'model_length':chain[0]['model_length'],'query_length':chain[0]['query_length'],
                'model_segments':json.dumps(models),'query_segments':json.dumps(queries),
                'component_count':len(chain),'raw_component_count':len(parts),
                'largest_model_gap':max([b[0]-a[1]-1 for a,b in zip(models,models[1:])],default=0),
                'largest_query_gap':max([b[0]-a[1]-1 for a,b in zip(queries,queries[1:])],default=0),
                'max_component_independent_evalue':max(p['independent_evalue'] for p in chain),
                'coverage_threshold':policy['min_model_coverage'],
                'scope':'HMM覆盖仅统计选中片段，不填补模型空洞，不证明蛋白完整折叠、催化中心占据或反应功能。'})
    return dict(result)


def structure_correspondences(sequences, root):
    """Only reuse reviewed structure matches pinned to this full reference sequence."""
    refs=json.loads((root/'references/processed/references.json').read_text())
    policy=json.loads((root/'config/reference_domain_policy.json').read_text())
    checked=set()
    by_sha=defaultdict(list)
    for r in refs:
        for s in r.get('reference_role_evidence',{}).get('structural_correspondences',[]):
            if s.get('coverage_scope')=='near_complete_sequence' and s['query_coverage']>=.9:
                by_sha[s['sequence_sha256']].append(s)
    result={}
    for r in sequences:
        digest=hashlib.sha256(r['完整CDS蛋白序列'].encode()).hexdigest()
        matches=[s for s in by_sha[digest] if s['protein_accession']==r['蛋白登录号']]
        for s in matches:
            pdb=s['structure_entity'].split('__entity_')[0]
            source=policy['structural_sources'][pdb]
            if pdb not in checked:
                if hashlib.sha256((root/source['path']).read_bytes()).hexdigest()!=source['sha256']:
                    raise ValueError('Reviewed structural source changed: '+pdb)
                checked.add(pdb)
            if s['source_sha256']!=source['sha256']:
                raise ValueError('Structure correspondence no longer matches source: '+pdb)
        if matches:result[r['序列ID']]=matches
    return result
