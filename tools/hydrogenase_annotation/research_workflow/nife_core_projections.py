"""Uncalibrated NiFe core candidates in role-specific reference coordinates."""
import hashlib, json, math, subprocess, tempfile
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, BLAST_FIELDS, association_allowed, local_cluster, load_blast
from localization import sha
from tool_paths import resolve_tool
def _validate_hit(hit, queries, references):
    if hit['qseqid'] not in queries or hit['sseqid'] not in references: raise ValueError('Unknown BLAST sequence ownership')
    if len(hit['qseq'])!=len(hit['sseq']): raise ValueError('BLAST alignment lengths differ')
    for field, sequence, lo, hi in [('qseq',queries[hit['qseqid']],hit['qstart'],hit['qend']),('sseq',references[hit['sseqid']],hit['sstart'],hit['send'])]:
        observed=hit[field].replace('-','').replace('U','X').replace('O','X')
        if not 1<=lo<=hi<=len(sequence) or observed!=sequence[lo-1:hi].replace('U','X').replace('O','X'):
            raise ValueError('BLAST alignment differs from bound '+field+' sequence')
try:
    from .coding_completeness import supports as _coding_supports
except ImportError:
    from coding_completeness import supports as _coding_supports

VERSION='NiFe-role-class-core-projection-v2'

def _hash(sequence): return hashlib.sha256(sequence.encode()).hexdigest()

def _metric(hit, seqs, refs):
    _validate_hit(hit,seqs,refs)
    if hit['qlen']!=len(seqs[hit['qseqid']]) or hit['slen']!=len(refs[hit['sseqid']]) or not all(math.isfinite(hit[k]) for k in ('evalue','pident','bitscore')):
        raise ValueError('NiFe projection HSP length or score invalid')
    if hit['evalue']<0 or hit['bitscore']<0 or not 0<=hit['pident']<=100: raise ValueError('NiFe projection HSP score out of range')
    qp=hit['qstart']-1; pairs=[]
    for a,b in zip(hit['qseq'],hit['sseq']):
        if a!='-': qp+=1
        if a!='-' and b!='-': pairs.append((qp,a==b))
    return dict(identity=100*sum(v for _,v in pairs)/len(pairs) if pairs else 0,
        reference_coverage=len(pairs)/hit['slen'],query_coverage=len(pairs)/hit['qlen'],query_positions=[p for p,_ in pairs])

def _component(hit, metric, role):
    return dict(role=role,protein_id=hit['qseqid'],query_domain_start=hit['qstart'],query_domain_end=hit['qend'],
        resources={hit['qseqid']:metric['query_positions']},identity=metric['identity'],reference_coverage=metric['reference_coverage'],
        reference_id=hit['sseqid'],raw_target_alignment=hit)

def _review(model, refs):
    root=Path(ROOT);cfg=model['ag_core']['reference_projection']
    if cfg.get('version')!=VERSION: raise ValueError('NiFe reference projection scope changed')
    path=root/cfg['review_path'];review=json.loads(path.read_text());meta={r['reference_id']:r for r in json.loads((root/'references/processed/references.json').read_text())}
    for role in ('large','small'):
        rows=review['roles'][role]; by={r['reference_id']:r for r in rows}; ids=set(model['ag_core'][role+'_core_reference_ids'])
        if len(by)!=len(rows) or set(by)!=ids: raise ValueError('NiFe projection role review differs')
        for rid,row in by.items():
            m=meta[rid]
            if m['center']!='NiFe' or m['kind']!=('catalytic' if role=='large' else 'module') or m.get('reference_role_review_status')!='supported_reference_family' or (role=='small' and m['reference_role_evidence']['role_class']!='NiFe_small_NuoB_related'):
                raise ValueError('Unreviewed NiFe projection reference')
            if row['sequence_sha256']!=_hash(refs[rid]): raise ValueError('NiFe projection reference AA changed')
    bank=root/cfg['counter_bank'];records=list(SeqIO.parse(bank,'fasta'));counter={r.id:str(r.seq).upper() for r in records};declared={r['reference_id']:r for r in review['counters']}
    if len(counter)!=len(records) or len(declared)!=len(review['counters']) or set(counter)!=set(declared) or not counter:
        raise ValueError('NiFe small counterbank ownership changed')
    if any(_hash(s)!=declared[rid]['sequence_sha256'] for rid,s in counter.items()): raise ValueError('NiFe small counter AA changed')
    paths=[path,bank]
    return counter,{str(p):sha(p) for p in paths}

def _shared_large_roles(model, refs):
    """Return only independently certified same-role large-core equivalents."""
    root=Path(ROOT); cfg=model['ag_core']['reference_projection']
    rel=cfg['shared_large_role_review_path']
    path=root/rel
    cert=json.loads(path.read_text())
    if cert.get('version')!='NiFe-large-core-functional-equivalence-v1':raise ValueError('NiFe shared-large certificate version changed')
    allowed=set()
    for row in cert.get('rows',[]):
        allowed.add(row['reference_id'])
    return allowed,{str(path):sha(path)}

def _counter_search(queries, counter, model):
    if not queries:return {},{},{}
    root=Path(ROOT);bank=root/model['ag_core']['reference_projection']['counter_bank'];binary=resolve_tool(root, 'blastp')
    protected=json.loads((root/'references/processed/protected_fefe_sequence_hashes.json').read_text())
    hashes=set(protected.get('sequence_sha256',[]))
    if any(_hash(s.upper().rstrip('*')) in hashes for s in queries.values()): raise ValueError('Protected FeFe sequence cannot enter NiFe core counter search')
    normalized={q:s.upper().replace('U','X').replace('O','X').rstrip('*') for q,s in queries.items()}
    if any(not s or set(s)-set('ACDEFGHIKLMNPQRSTVWYBXZ') for s in normalized.values()): raise ValueError('Invalid NiFe small query alphabet')
    bindings={str(bank):sha(bank),str(binary):sha(binary)};dbsize=sum(map(len,counter.values()))
    cmd=[str(binary),'-query','-','-subject',str(bank),'-dbsize',str(dbsize),'-evalue',str(model['max_evalue']),'-max_target_seqs',str(len(counter)),'-max_hsps','1','-outfmt','6 '+' '.join(BLAST_FIELDS)]
    fasta=''.join('>'+q+'\n'+s+'\n' for q,s in normalized.items())
    proc=subprocess.run(cmd,input=fasta,text=True,capture_output=True,check=True)
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/'counter.tsv';p.write_text(proc.stdout)
        raw=load_blast(p,dict(blast_evalue=model['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
    by={q:[] for q in queries}
    for h in raw:
        metric=_metric(h,normalized,counter);by[h['qseqid']].append(dict(raw_hsp=h,**metric))
    if any(sha(Path(p))!=d for p,d in bindings.items()): raise ValueError('NiFe counter resource changed during search')
    audit=dict(command=cmd,reference_records=len(counter),reference_residues=dbsize,normalized_query_fasta_sha256=_hash(fasta),
        query_AA_sha256={q:_hash(s) for q,s in queries.items()},stdout_sha256=_hash(proc.stdout),normalization='Search only U/O to X; original called proteins unchanged')
    return by,bindings,audit

def prepare_fallbacks(annotation, results, anchors, genes, domains, seqs, raw, refs, identities, model, coding_completeness=None):
    from mvh_hdr_function import _ag_profile, _ag_reference_audit
    mode=json.loads((Path(annotation)/'run_manifest.json').read_text())['input_mode']
    by={g['id']:dict(g,**{k:int(g[k]) for k in ('start','end','strand')}) for g in genes}
    proof = coding_completeness if isinstance(coding_completeness, dict) else {}
    def proof_ok(pid):
        row = proof.get(pid) if isinstance(proof.get(pid), dict) else None
        return bool(row and _coding_supports(row, pid,
                                            expected_protein_sha256=_hash(seqs[pid]),
                                            gene=by[pid]))
    eligible={a['protein_id'] for a in anchors if a.get('catalytic_family')=='NiFe' and results.get(a['protein_id'],{}).get('status')=='unresolved'
        and identities.get(a['protein_id']) in ('calibrated_high','sequence_supported_uncalibrated_core') and (not str(by[a['protein_id']].get('partial')).lower()=='true' or proof_ok(a['protein_id'])) and mode!='transcript'}
    if not eligible:return {},{}
    counter,inputs=_review(model,refs);large=set(model['ag_core']['large_core_reference_ids']);shared,shared_inputs=_shared_large_roles(model,refs);inputs.update(shared_inputs);large_audit_group=large|shared;small=set(model['ag_core']['small_core_reference_ids'])
    large_hits={};small_hits={}
    for h in raw:
        if h['sseqid'] not in large|small or h['evalue']>model['max_evalue']: continue
        metric=_metric(h,seqs,refs)
        if metric['identity']>=model['min_identity'] and metric['reference_coverage']>=model['min_coverage'] and (h['sseqid'] in small or metric['query_coverage']>=model['min_coverage']):
            if h['sseqid'] in large:large_hits.setdefault(h['qseqid'],[]).append((h,metric))
            else:small_hits.setdefault(by[h['qseqid']]['contig'],[]).append((h,metric))
    candidates={q:{} for q in sorted(eligible)};projections={q:[] for q in sorted(eligible)}
    for q in sorted(eligible):
        for ah,am in large_hits.get(q,[]):
            a=_component(ah,am,'MvhA');aa=_ag_reference_audit(a,raw,seqs,refs,large_audit_group,model)
            if aa['status']!='supported':continue
            for gh,gm in small_hits.get(by[q]['contig'],[]):
                pid=gh['qseqid']
                if gh['sseqid'] not in small or pid==q or (str(by[pid].get('partial')).lower()=='true' and not proof_ok(pid)):continue
                if not association_allowed(by[q],by[pid],mode) or not local_cluster([by[q],by[pid]],model['linkage']):continue
                g=_component(gh,gm,'MvhG');profile=_ag_profile(g,domains,seqs,model)
                if not profile:continue
                ga=_ag_reference_audit(g,raw,seqs,refs,small,model)
                if ga['status']!='supported':continue
                witness=dict(large_reference=ah['sseqid'],small_reference=gh['sseqid'],large_identity=am['identity'],small_identity=gm['identity'],small_protein=pid,small_whole_query_coverage=gm['query_coverage'],small_paired_reference_coverage=gm['reference_coverage'],small_coverage_scope='complete_reference_role_and_owned_GA_profile; unrelated fusion tail is not a core role')
                projections[q].append(witness)
                rank=(min(am['identity'],gm['identity']),am['identity']+gm['identity'],ah['sseqid'],gh['sseqid'])
                if pid not in candidates[q] or rank>candidates[q][pid]['rank']:
                    candidates[q][pid]=dict(rank=rank,hypothesis=[a,g],small_profile_evidence=profile,large_reference_audit=aa,small_reference_audit=ga)
    queries={pid:seqs[pid] for pid in sorted({p for cs in candidates.values() for p in cs})};hits,bindings,search=_counter_search(queries,counter,model);inputs.update(bindings)
    output={}
    for q,cs in candidates.items():
        usable=[];excluded=[]
        for pid,c in cs.items():
            g=c['hypothesis'][1];positions=set(g['resources'][pid]);best=max(x['paired_identity'] for x in c['small_reference_audit']['positive_references']);rivals=[]
            for h in hits.get(pid,[]):
                overlap=len(positions&set(h['query_positions']))/len(h['query_positions']) if h['query_positions'] else 0
                if h['reference_coverage']>=model['min_coverage'] and overlap>=.7 and h['identity']>=best*model['ag_core']['competition_identity_ratio']:
                    rivals.append(dict(h,resource_overlap=overlap))
            c['small_counter_audit']=dict(status='unresolved' if rivals else 'supported',positive_best_paired_identity=best,competing_nonhydrogenase_small_references=rivals,all_counter_hits=hits.get(pid,[]))
            (excluded if rivals else usable).append(c)
        e=dict(identity_status=identities[q],search_status='complete',hypotheses=[c['hypothesis'] for c in usable],role_candidates={},small_profile_evidence=[],large_reference_audit={},small_reference_audit={},resource_conflicts=[],coding_completeness_proof_used=[pid for pid in sorted(candidates.get(q,{})) if proof_ok(pid)],
            projection_route=VERSION,interface_status='unmeasured',assembly_status='unresolved',all_reference_projections=projections[q],shared_large_core_reference_ids=sorted(shared),excluded_counter_candidates=excluded,counter_search=search,
            scope='Independent reviewed large/small role classes in their own reference coordinates; local candidate only, no exact lineage, physical interface or electron-partner transfer.')
        if len(usable)==1:e.update({k:usable[0][k] for k in ('small_profile_evidence','large_reference_audit','small_reference_audit','small_counter_audit')})
        output[q]=e
    if any(sha(Path(p))!=d for p,d in inputs.items()): raise ValueError('NiFe projection resource changed during evaluation')
    return output,inputs
