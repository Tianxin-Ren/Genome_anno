"""Keep complementary FeFe scaffold homology separate from complete core proof."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from annotate import BLAST_FIELDS, load_blast, local_cluster, residue_states, repeat_mapping_ambiguity
from active_centers import collinear_hits
from circular_coordinates import valid_cds_interval
def _validate_hit(hit, queries, references):
    if hit['qseqid'] not in queries or hit['sseqid'] not in references: raise ValueError('Unknown BLAST sequence ownership')
    if len(hit['qseq'])!=len(hit['sseq']): raise ValueError('BLAST alignment lengths differ')
    for field, sequence, lo, hi in [('qseq',queries[hit['qseqid']],hit['qstart'],hit['qend']),('sseq',references[hit['sseqid']],hit['sstart'],hit['send'])]:
        observed=hit[field].replace('-','').replace('U','X').replace('O','X')
        if not 1<=lo<=hi<=len(sequence) or observed!=sequence[lo-1:hi].replace('U','X').replace('O','X'):
            raise ValueError('BLAST alignment differs from bound '+field+' sequence')

ATTRIBUTE = 'FeFe_complementary_scaffold_context'
RULE = 'FeFe-complementary-structural-context-v1'


def sha(p):
    with Path(p).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()


def _rows(p):
    with p.open() as f: return list(csv.DictReader(f,delimiter='\t'))


def _mapping(hits, model, sequence):
    components=collinear_hits(hits); mapping={}; ambiguous=set()
    window=list(range(model['window_start'],model['window_end']+1))
    for hit in components:
        states,piece=residue_states(hit,window)
        if mapping.keys() & piece.keys(): raise ValueError('Overlapping context reference segments')
        mapping.update(piece)
        ambiguous.update(s['reference_position'] for s in states if s['coverage_reason']=='repeat_indel_alignment_ambiguity')
    ambiguous.update(repeat_mapping_ambiguity(mapping,model['sequence'],sequence))
    return {p:dict(protein_id=components[0]['qseqid'],query_position=q,residue=r,
                   reference_residue=model['sequence'][p-1])
            for p,(q,r) in mapping.items() if q is not None and p not in ambiguous}, components


def _local(a,b,thresholds):
    for g in (a,b):
        if (str(g.get('partial','')).lower() in {'true','1','yes'}
                or g.get('coordinate_evidence')!='genomic_CDS' or g.get('topology') not in {'linear','circular'}):
            return False,'incomplete_or_unverified_CDS'
        try:
            if int(g['strand']) not in {-1,1} or not valid_cds_interval(g,allow_unknown_length=False):
                return False,'invalid_CDS_coordinates'
        except (ValueError,TypeError,KeyError): return False,'invalid_CDS_coordinates'
    if any(a.get(k,'')!=b.get(k,'') for k in ('contig','strand','topology','contig_length','bin_id')):
        return False,'coordinate_source_or_orientation_mismatch'
    return (True,'same_source_local_CDS') if local_cluster([a,b],thresholds) else (False,'local_linkage_not_supported')


def _metrics(mapping,model,params):
    contexts=[]
    for site in model['sites']:
        if site['requirement'] not in {'required_for_template_coordination','canonical_reaction_signature'}:continue
        p=site['reference_position']; lo=max(model['window_start'],p-params['site_context_flank']);hi=min(model['window_end'],p+params['site_context_flank'])
        covered=[mapping[n] for n in range(lo,hi+1) if n in mapping]
        contexts.append(dict(reference_position=p,start=lo,end=hi,coverage=len(covered)/(hi-lo+1),
            identity=100*sum(r['residue']==r['reference_residue'] for r in covered)/max(1,len(covered))))
    return dict(window_coverage=len(mapping)/(model['window_end']-model['window_start']+1),
        window_identity=100*sum(r['residue']==r['reference_residue'] for r in mapping.values())/max(1,len(mapping)),
        contexts=contexts,context_alignment_supported=all(c['coverage']>=params['min_window_coverage'] and c['identity']>=params['min_window_identity'] for c in contexts),
        min_required_context_coverage=min(c['coverage'] for c in contexts),min_required_context_identity=min(c['identity'] for c in contexts))


def evaluate(anchors,genes,sequences,models,hits,summaries,params,thresholds):
    """Return hypotheses with explicit residue owners; never concatenate peptides."""
    grouped=defaultdict(list)
    for hit in hits:
        if hit['sseqid'] not in models or hit['evalue']>params['max_evalue']:continue
        _validate_hit(hit,sequences,{rid:m['sequence'] for rid,m in models.items()})
        if hit['qlen']!=len(sequences[hit['qseqid']]) or hit['slen']!=len(models[hit['sseqid']]['sequence']):
            raise ValueError('Context HSP sequence lengths differ')
        grouped[(hit['qseqid'],hit['sseqid'])].append(hit)
    mapped={key:_mapping(value,models[key[1]],sequences[key[0]]) for key,value in grouped.items()}
    results={}
    for anchor in sorted(anchors,key=lambda a:a['protein_id']):
        aid=anchor['protein_id']; summary=summaries.get(aid,{})
        rid=summary.get('reference_id'); model=models.get(rid)
        if anchor.get('catalytic_family')!='FeFe' or not model or summary.get('direct_ligand_status')!='present':continue
        am,ah=mapped.get((aid,rid),({},[])); direct=[s for s in model['sites'] if s['residue_role']=='direct_metal_ligand']
        if len(direct)!=4 or any(s['reference_position'] not in am or am[s['reference_position']]['residue'] not in s['expected_residue'] for s in direct):continue
        original=_metrics(am,model,params)
        if original['context_alignment_supported'] and original['window_coverage']>=params['min_window_coverage']:continue
        candidates=[]
        for pid in sorted(genes):
            if pid==aid:continue
            pm,ph=mapped.get((pid,rid),({},[]))
            if not pm or am.keys() & pm.keys():continue
            union=am|pm; metrics=_metrics(union,model,params)
            local,reason=_local(genes[aid],genes[pid],thresholds)
            blockers=[] if local else [reason]
            if not metrics['context_alignment_supported']:blockers.append('required_site_context_incomplete')
            if metrics['window_coverage']<1: blockers.append('reviewed_window_incomplete')
            competing=[h for h in grouped[(pid,rid)] if h not in ph and h['bitscore']>=max(c['bitscore'] for c in ph)*params['competitive_bitscore_ratio']]
            if competing:blockers.append('alternative_context_alignment')
            candidates.append(dict(members=[dict(role='catalytic_core',protein_id=aid),dict(role='complementary_structural_context',protein_id=pid)],
                association_status='local' if local else 'rejected',association_reason=reason,blockers=blockers,
                new_reference_positions=sorted(pm),reference_position_owners={str(p):r for p,r in sorted(union.items())},
                anchor_alignment=ah,partner_alignment=ph,competing_partner_alignments=competing,metrics=metrics))
        if candidates:results[aid]=dict(rule_version=RULE,anchor_id=aid,reference_id=rid,candidates=candidates,
            original_metrics=original,status='unresolved',selected_members=[],physical_assembly='unknown',query_activity='not_tested',
            scope='Partial complementary scaffold homology and local CDS context only; no complete core, calibrated identity, native assembly or functional partner assignment.')
    owners=defaultdict(set)
    for aid,result in results.items():
        for candidate in result['candidates']:
            if candidate['association_status']=='local':owners[candidate['members'][1]['protein_id']].add(aid)
    for aid,result in results.items():
        local=[c for c in result['candidates'] if c['association_status']=='local']
        for c in local:
            c['alternative_anchors']=sorted(owners[c['members'][1]['protein_id']]-{aid})
        if local:
            ambiguous=len(local)>1 or any(c['alternative_anchors'] or c['competing_partner_alignments'] for c in local)
            result['status']='ambiguous' if ambiguous else 'conditional'
            if not ambiguous:result['selected_members']=local[0]['members']
    return results


def prepare(annotation,output,anchors,genes,sequences,root):
    if not anchors:return {},{}
    annotation,output,root=Path(annotation),Path(output),Path(root)
    required=[annotation/n for n in ('active_center_blast.tsv','active_center_summary.tsv','genes.tsv','proteins.faa','run_manifest.json')]
    if any(not p.is_file() for p in required):return {},{}
    directory=root/'references/processed/active_centers';manifest=json.loads((directory/'manifest.json').read_text())
    before={str(p):sha(p) for p in required}
    run=json.loads((annotation/'run_manifest.json').read_text())
    expected=run.get('output_sha256',{})
    if run.get('status')!='complete' or not run.get('normalized_proteins_sha256') or any(n not in expected for n in ('genes.tsv','active_center_summary.tsv')):
        return {},{}
    if sha(annotation/'proteins.faa')!=run['normalized_proteins_sha256']:
        raise ValueError('Context normalized proteins differ from producer manifest')
    for name in ('genes.tsv','active_center_summary.tsv','active_center_blast.tsv','proteins.faa'):
        if name in expected and sha(annotation/name)!=expected[name]:
            raise ValueError('Context annotation differs from producer manifest: '+name)
    gene_rows=_rows(annotation/'genes.tsv'); supplied=genes if isinstance(genes,dict) else {g['id']:g for g in genes}
    with (annotation/'proteins.faa').open() as handle:
        native={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    if native!=sequences:raise ValueError('Context query sequences differ from annotation')
    for row in gene_rows:
        if row['id'] not in supplied or any(str(supplied[row['id']].get(k,''))!=str(v) for k,v in row.items()):
            raise ValueError('Context gene arguments differ from source table')
    typed={g['id']:{**g,**{k:int(g[k]) for k in ('start','end','strand','contig_length') if g.get(k)}} for g in gene_rows}
    raw=annotation/'active_center_blast.tsv'
    if any(len(line.split('\t'))!=len(BLAST_FIELDS) for line in raw.read_text().splitlines()):raise ValueError('Malformed context HSP row')
    params=manifest['parameters'];thresholds=run['thresholds']
    models={m['reference_id']:m for m in json.loads((directory/'models.json').read_text()) if m['family']=='FeFe'}
    hits=load_blast(raw,dict(blast_evalue=params['max_evalue'],min_identity=0,min_query_coverage=0,min_reference_coverage=0))
    results=evaluate(anchors,typed,sequences,models,hits,{r['protein_id']:r for r in _rows(annotation/'active_center_summary.tsv')},params,thresholds)
    if any(sha(p)!=expected for p,expected in before.items()):raise ValueError('Context evidence changed during analysis')
    output.mkdir(parents=True,exist_ok=False)
    artifact=output/'evidence.json'
    artifact.write_text(json.dumps(dict(rule_version=RULE,results=results,input_sha256=before,
        provenance_scope='Current raw HSPs are sequence-validated and newly pinned here; no historical run-manifest raw-HSP binding is invented.'),ensure_ascii=False,indent=2)+'\n')
    return results,{**before,str(artifact):sha(artifact)}


def attributes(result):
    if not result.get('candidates'):return []
    state={'conditional':'conditional_partial_context','ambiguous':'ambiguous_partial_context'}.get(result['status'],'unresolved')
    return [dict(attribute=ATTRIBUTE,status=state,value=result,evidence_dependencies=['independent FeFe core identity','exact component HSP ownership','separate local CDS linkage'],
        evidence_source='complementary_core_context_evidence.json;complementary_core_context/evidence.json',limitations=result['scope'])]
