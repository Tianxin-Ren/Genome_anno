"""Uncalibrated FeFe backbone support; homology contexts never replace atomic evidence."""
import csv,hashlib,json,math,subprocess
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT,BLAST_FIELDS,load_blast
from active_centers import map_model
from catalytic_confirmation import core_similarity
from fefe_profile_ligands import _digest,_stockholm,_site_map
def _validate_hit(hit, queries, references):
    if hit['qseqid'] not in queries or hit['sseqid'] not in references: raise ValueError('Unknown BLAST sequence ownership')
    if len(hit['qseq'])!=len(hit['sseq']): raise ValueError('BLAST alignment lengths differ')
    for field, sequence, lo, hi in [('qseq',queries[hit['qseqid']],hit['qstart'],hit['qend']),('sseq',references[hit['sseqid']],hit['sstart'],hit['send'])]:
        observed=hit[field].replace('-','').replace('U','X').replace('O','X')
        if not 1<=lo<=hi<=len(sequence) or observed!=sequence[lo-1:hi].replace('U','X').replace('O','X'):
            raise ValueError('BLAST alignment differs from bound '+field+' sequence')
from tool_paths import resolve_tool

RULE_VERSION='FeFe-active-backbone-homology-pattern-v2'
MODEL_REL='research_workflow/fefe_experimental_backbone_model.json'
ALLOWED_REASONS={'reviewed_single_chain_center_or_span_not_supported','required_domain_architecture_not_colocated_with_center'}

def _sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def _load_model(root=ROOT):
    root=Path(root);p=root/MODEL_REL;m=json.loads(p.read_text())
    if m['version']!=RULE_VERSION or m['reference_positions']!=[138,193,293,297] or m['template_positions']!=[300,355,499,503] or m['profile_match_columns']!=[77,130,272,276] or m['accepted_prior_ligand_states']!=['not_covered','present']:raise ValueError('FeFe backbone source-coordinate scope changed')
    with (root/m['bank_path']).open() as h:records=list(SeqIO.parse(h,'fasta'))
    refs={r.id:str(r.seq) for r in records}
    tam=refs[m['positive_reference']]
    params=json.loads((root/'config/active_center_review.json').read_text())['parameters'];rules=json.loads((root/'config/catalytic_confirmation_thresholds.json').read_text())['families']['FeFe']
    mapper=dict(sequence=tam,window_start=1,window_end=len(tam),reference_id=m['positive_reference'],family='FeFe',source_sha256=m['positive_sequence_sha256'],evidence_url='https://doi.org/10.1021/acs.biochem.5c00784',author_positions={},sites=[dict(reference_position=p,expected_residue='C',residue_role='direct_metal_ligand',requirement='required_for_template_coordination') for p in m['reference_positions']])
    return m,refs,params,rules,mapper

def _queries(annotation):
    with (annotation/'proteins.faa').open() as h:records=list(SeqIO.parse(h,'fasta'))
    seqs={r.id:str(r.seq).upper().rstrip('*') for r in records}
    with (annotation/'genes.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
    genes={r['id']:r for r in rows}
    if len(seqs)!=len(records) or len(genes)!=len(rows) or seqs.keys()!=genes.keys():raise ValueError('FeFe backbone gene/protein ownership mismatch')
    if any(not s or not s.isalpha() for s in seqs.values()):raise ValueError('FeFe backbone input protein alphabet invalid')
    return seqs,genes

def _eligible(features,genes):
    return {q:f for q,f in features.items() if f.get('family')=='FeFe' and q in genes and genes[q].get('partial') in [False,'False'] and f.get('source_active_center',{}).get('direct_ligand_status') in ['not_covered','present']}

def _decision_blocks(r,parameters,context):
    b=[];positive=r.get('positive_comparison',{});negative=r.get('NARF_countercomparison',{})
    numeric=[positive.get('similarity',0),positive.get('core_coverage',0),negative.get('similarity',0),r.get('paired_reference_coverage',0)]+[s.get(k,0) for s in r.get('sites',[]) for k in ['site_context_coverage','site_context_identity']]
    if any(not isinstance(x,(int,float)) or isinstance(x,bool) or not math.isfinite(x) for x in numeric):return ['nonfinite_or_invalid_backbone_evidence']
    if not 0<=positive.get('similarity',0)<=100 or not 0<=negative.get('similarity',0)<=100 or not 0<=positive.get('core_coverage',0)<=1 or not 0<=r.get('paired_reference_coverage',0)<=1:return ['backbone_evidence_out_of_range']
    if positive.get('similarity',0)<parameters['min_core_similarity'] or positive.get('core_coverage',0)<.7:b.append('active_reference_core_comparison_insufficient')
    if r.get('paired_reference_coverage',0)<parameters['min_domain_coverage']:b.append('active_reference_backbone_incomplete')
    if negative.get('search_status') not in ['qualified_hit','no_qualified_hit_under_search_scope'] or positive.get('similarity',0)-negative.get('similarity',0)<parameters['min_competitor_margin']:b.append('NARF_countercomparison_insufficient')
    if r.get('query_partial') is not False:b.append('partial_or_unknown_query')
    if r.get('backbone_span_status')!='supported':b.append('backbone_context_span_insufficient')
    if len(r.get('sites',[]))!=4:b.append('four_homologous_Cys_missing')
    positions=[]
    for site in r.get('sites',[]):
        if site.get('status')!='present' or site.get('observed_residue')!='C':b.append('homologous_Cys_conflict_or_ambiguity')
        if site.get('site_context_coverage',0)<context['min_window_coverage'] or site.get('site_context_identity',0)<context['min_window_identity']:b.append('ligand_homology_context_insufficient')
        p=site.get('query_position')
        if not isinstance(p,int) or isinstance(p,bool) or not 1<=p<=r.get('query_length',0):b.append('homology_site_coordinate_invalid')
        else:positions.append(p)
    if len(set(positions))!=4 or positions!=sorted(positions):b.append('homology_site_coordinates_not_unique_ordered')
    return sorted(set(b))

def _evaluate(q,f,hits,sequence,gene,m,refs,context,parameters,mapper,model_sha):
    positive=hits.get(m['positive_reference']);negative=hits.get(m['negative_reference'])
    if positive is None:return None
    normalized=sequence.replace('U','X').replace('O','X')
    for h in hits.values():
        _validate_hit(h,{q:normalized},refs)
        if h['qlen']!=len(normalized) or h['slen']!=len(refs[h['sseqid']]) or not all(math.isfinite(h[k]) for k in ['evalue','bitscore','pident']) or not 0<=h['evalue']<=m['max_evalue'] or h['bitscore']<0 or not 0<=h['pident']<=100:raise ValueError('FeFe backbone HSP length/score invalid')
    positions=f['core_positions']
    if not isinstance(positions,list):raise ValueError('FeFe query-core positions invalid')
    if not positions:return None
    if len(set(positions))!=len(positions) or any(not isinstance(p,int) or isinstance(p,bool) or not 1<=p<=len(sequence) for p in positions):raise ValueError('FeFe query-core positions invalid')
    mapping=map_model([positive],mapper,context,normalized)
    # map_model supplies sequence-context and repeat/indel validation only. Its
    # structural display labels are replaced: Tam coordinates are homology based.
    sites=[{k:s[k] for k in ['reference_position','query_position','observed_residue','status','coverage_reason','site_context_coverage','site_context_identity']} for s in mapping['sites']]
    pc=core_similarity(positive,positions);nc=core_similarity(negative,positions) if negative else dict(similarity=0.,core_coverage=0.)
    paired=sum(a!='-' and b!='-' for a,b in zip(positive['qseq'],positive['sseq']))
    r=dict(protein_id=q,rule_version=RULE_VERSION,feature_sha256=_digest(f),model_sha256=model_sha,query_sha256=hashlib.sha256(sequence.encode()).hexdigest(),query_length=len(sequence),query_partial=gene.get('partial') not in [False,'False'],positive_reference=m['positive_reference'],positive_reference_sha256=m['positive_sequence_sha256'],negative_reference=m['negative_reference'],negative_reference_sha256=m['negative_sequence_sha256'],positive_comparison=pc,NARF_countercomparison=dict(nc,search_status='qualified_hit' if negative else 'no_qualified_hit_under_search_scope',evalue_max=m['max_evalue']),paired_reference_coverage=paired/len(refs[m['positive_reference']]),backbone_span_status=mapping['template_span_alignment_status'],sites=sites,parameters=parameters,source_mapping_scope=m['source_mapping_scope'],scope=m['scope'])
    r['blockers']=_decision_blocks(r,parameters,context);r['status']='unresolved' if r['blockers'] else 'experimental_backbone_pattern_supported'
    return r

def record_supported(record,feature,parameters,root=ROOT):
    if not record or record.get('status')!='experimental_backbone_pattern_supported':return False
    m,refs,context,rules,_=_load_model(root)
    if parameters!=rules:raise ValueError('FeFe backbone rule parameters differ from saved scope')
    if record.get('feature_sha256')!=_digest(feature) or record.get('protein_id')!=feature.get('protein_id') or record.get('rule_version')!=RULE_VERSION or record.get('model_sha256')!=_sha(Path(root)/MODEL_REL):return False
    if record.get('positive_reference_sha256')!=m['positive_sequence_sha256'] or record.get('negative_reference_sha256')!=m['negative_sequence_sha256'] or record.get('positive_reference')!=m['positive_reference'] or record.get('negative_reference')!=m['negative_reference']:return False
    if record.get('parameters')!=parameters or [s.get('reference_position') for s in record.get('sites',[])]!=m['reference_positions']:return False
    if record.get('NARF_countercomparison',{}).get('evalue_max')!=m['max_evalue']:return False
    return not record.get('blockers') and not _decision_blocks(record,parameters,context)

def _decode(annotation,features,parameters,output,opaque,root):
    m,refs,context,rules,mapper=_load_model(root);queries,genes=_queries(annotation);eligible=_eligible(features,genes)
    if parameters!=rules or set(opaque)!=set(eligible) or len(set(opaque.values()))!=len(opaque):raise ValueError('FeFe backbone input scope/opaque ownership changed')
    reverse={v:k for k,v in opaque.items()};raw=load_blast(output/'blast.tsv',dict(blast_evalue=m['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0));hits={q:{} for q in eligible}
    for h in raw:
        if h['qseqid'] not in reverse or h['sseqid'] not in refs or h['evalue']>m['max_evalue']:raise ValueError('FeFe backbone raw search ownership/scope changed')
        q=reverse[h['qseqid']]
        if h['sseqid'] in hits[q]:raise ValueError('FeFe backbone repeated HSP')
        h=dict(h,qseqid=q);hits[q][h['sseqid']]=h
    results={}
    for q,f in eligible.items():
        r=_evaluate(q,f,hits[q],queries[q],genes[q],m,refs,context,parameters,mapper,_sha(Path(root)/MODEL_REL))
        if r:results[q]=r
    return results

def prepare(annotation,features,output,parameters=None,root=ROOT):
    annotation,output,root=Path(annotation).resolve(),Path(output).resolve(),Path(root).resolve()
    if not any(f.get('family')=='FeFe' and f.get('source_active_center',{}).get('direct_ligand_status') in ['not_covered','present'] for f in features.values()):return {},{}
    queries,genes=_queries(annotation);eligible=_eligible(features,genes)
    if not eligible:return {},{}
    m,refs,context,rules,_=_load_model(root);parameters=parameters or rules
    if parameters!=rules:raise ValueError('FeFe backbone parameters differ from frozen rule scope')
    paths=[annotation/'proteins.faa',annotation/'genes.tsv',root/MODEL_REL,Path(__file__)]
    fp=annotation/'confirmation/features.json'
    if fp.exists():
        saved={f['protein_id']:f for f in json.loads(fp.read_text())}
        if any(saved.get(q)!=f for q,f in eligible.items()):raise ValueError('FeFe backbone feature source mismatch')
        paths.append(fp)
    inputs={str(p.resolve()):_sha(p) for p in paths};output.mkdir(parents=True,exist_ok=False);opaque={q:f'q{i:08d}' for i,q in enumerate(sorted(eligible))}
    (output/'opaque_ids.json').write_text(json.dumps(opaque,sort_keys=True,indent=2)+'\n');(output/'query.faa').write_text(''.join('>'+opaque[q]+'\n'+queries[q].replace('U','X').replace('O','X')+'\n' for q in sorted(eligible)))
    cmd=[str(resolve_tool(root, 'blastp')),'-query',str(output/'query.faa'),'-subject',str(root/m['bank_path']),'-dbsize',str(sum(len(s) for s in refs.values())),'-outfmt','6 '+' '.join(BLAST_FIELDS),'-max_target_seqs',str(len(refs)),'-max_hsps','1','-num_threads','2','-evalue',str(m['max_evalue']),'-out',str(output/'blast.tsv')];binary=_sha(Path(cmd[0]));run=subprocess.run(cmd,check=True,capture_output=True,text=True);(output/'blast.log').write_text(run.stdout+run.stderr)
    records=_decode(annotation,features,parameters,output,opaque,root);(output/'backbone_evidence.json').write_text(json.dumps(records,sort_keys=True,indent=2)+'\n');outputs={p.name:_sha(p) for p in output.iterdir() if p.is_file()}
    manifest=dict(status='complete',rule_version=RULE_VERSION,annotation=str(annotation),parameters=parameters,feature_sha256=_digest(eligible),command=cmd,binary_sha256=binary,bank_records=len(refs),input_sha256=inputs,output_sha256=outputs,scope=m['scope']);(output/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
    if any(_sha(Path(p))!=h for p,h in inputs.items()) or _sha(Path(cmd[0]))!=binary:raise ValueError('FeFe backbone input changed during execution')
    return records,{**inputs,**{str(output/n):h for n,h in outputs.items()},str(output/'manifest.json'):_sha(output/'manifest.json')}

def verify_saved(output,annotation=None,features=None,parameters=None,root=ROOT):
    if annotation is None or features is None or parameters is None:raise ValueError('FeFe backbone verification requires actual input context')
    output=Path(output).resolve();annotation=Path(annotation).resolve();root=Path(root).resolve();manifest=json.loads((output/'manifest.json').read_text())
    if manifest.get('status')!='complete' or manifest.get('rule_version')!=RULE_VERSION or manifest.get('annotation')!=str(annotation) or manifest.get('parameters')!=parameters:raise ValueError('FeFe backbone saved scope changed')
    for p,h in manifest['input_sha256'].items():
        if _sha(Path(p))!=h:raise ValueError('FeFe backbone saved input changed')
    for n,h in manifest['output_sha256'].items():
        if Path(n).name!=n or _sha(output/n)!=h:raise ValueError('FeFe backbone saved output changed')
    if _sha(Path(manifest['command'][0]))!=manifest['binary_sha256']:raise ValueError('FeFe backbone search binary changed')
    queries,genes=_queries(annotation);eligible=_eligible(features,genes)
    if _digest(eligible)!=manifest['feature_sha256']:raise ValueError('FeFe backbone saved feature scope changed')
    opaque=json.loads((output/'opaque_ids.json').read_text());actual=_decode(annotation,features,parameters,output,opaque,root);saved=json.loads((output/'backbone_evidence.json').read_text())
    if saved!=actual:raise ValueError('FeFe backbone saved evidence contradicts actual alignment')
    return saved,manifest
