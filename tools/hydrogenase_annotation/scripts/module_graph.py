"""Coordinate-based module assignment; linkage is sequence evidence, not interaction proof."""
import hashlib
import json
from collections import defaultdict
from itertools import combinations
from circular_coordinates import circular_overlap, circular_span


def load_policy(root, references):
    path=root/'config/module_graph.json';policy=json.loads(path.read_text())
    if hashlib.sha256((root/'config/systems.json').read_bytes()).hexdigest()!=policy['systems_sha256']:
        raise ValueError('Module graph systems changed')
    refs={r['reference_id']:r for r in references};seen=set()
    for graph in policy['systems']:
        ids=set()
        for node in graph['nodes']:
            r=refs[node['reference_id']]
            if node['module_category'] not in policy['cross_anchor_sharing']['category_modes']:
                raise ValueError('Unreviewed cross-anchor module sharing category')
            if (node['node_id'] in ids or r['structure_id']!=graph['structure_id'] or
                r['sequence_sha256']!=node['reference_sequence_sha256'] or
                not 1<=node['reference_start']<=node['reference_end']<=r['protein_region_end']-r['protein_region_start']+1):
                raise ValueError('Module graph reference/coordinates changed')
            ids.add(node['node_id']);seen.add(node['reference_id'])
    if seen!={r['reference_id'] for r in references if r['kind'] in ['catalytic','module','context']}:
        raise ValueError('Module graph has missing or unexpected references')
    return policy


def project(hit,node):
    """Keep aligned residue resources; gaps cannot become covered domain positions."""
    if len(hit['qseq'])!=len(hit['sseq']):raise ValueError('Unequal BLAST aligned sequence lengths')
    rp=hit['sstart']-1;qp=hit['qstart']-1;pairs=[]
    for q,r in zip(hit['qseq'],hit['sseq']):
        if r!='-':rp+=1
        if q!='-':qp+=1
        if r!='-' and q!='-' and node['reference_start']<=rp<=node['reference_end']:
            pairs.append((rp,qp,q==r))
    if rp!=hit['send'] or qp!=hit['qend'] or not 1<=qp<=hit['qlen'] or not 1<=rp<=hit['slen']:
        raise ValueError('BLAST alignment coordinates do not match aligned residues')
    if not pairs:return None
    return dict(hit=hit,protein_id=hit['qseqid'],reference_positions={p[0] for p in pairs},
        query_positions={p[1] for p in pairs},identity=100*sum(p[2] for p in pairs)/len(pairs),
        reference_start=min(p[0] for p in pairs),reference_end=max(p[0] for p in pairs))


def resource_conflict(parts, selected):
    return any(a['protein_id']==b['protein_id'] and a['query_positions']&b['query_positions'] for a in parts for b in selected)


def _span_within(genes, thresholds, local_cluster):
    from annotate import _circular_contig_length
    if not genes or len({g['contig'] for g in genes})!=1 or len({g['strand'] for g in genes})!=1: return False
    length=_circular_contig_length(genes)
    if length is None:
        return max(g['end'] for g in genes)-min(g['start'] for g in genes)+1 <= thresholds['local_max_span_bp']
    if length <= 0:
        return False
    try:
        return circular_span([(g['start'], g['end']) for g in genes], length) <= thresholds['local_max_span_bp']
    except ValueError:
        return False


def compatible_fragments(a,b,genes,thresholds,local_cluster):
    from annotate import _circular_contig_length
    ga,gb=genes[a['protein_id']],genes[b['protein_id']]
    if (a['protein_id']==b['protein_id'] or a['reference_end']>=b['reference_start'] or
        ga['contig']!=gb['contig'] or ga['strand']!=gb['strand']):return False
    # Preserve reference collinearity on linear replicons.  A verified
    # circular pair may cross the coordinate origin, in which case the
    # canonical locality guard supplies the wrapped geometry.
    length = _circular_contig_length([ga,gb])
    if length == 0: return False
    left,right=(ga,gb) if ga['strand']==1 else (gb,ga)
    gap = right['start']-left['end']-1
    if length is None:
        return 0 <= gap <= thresholds['local_max_intergenic_gap_bp']
    try:
        if circular_overlap((ga['start'],ga['end']), (gb['start'],gb['end']), length): return False
    except ValueError:
        return False
    return gap % length <= thresholds['local_max_intergenic_gap_bp'] and local_cluster([ga,gb], thresholds)


def realizations(parts,node,genes,thresholds,policy,local_cluster):
    length=node['reference_end']-node['reference_start']+1;options=[];truncated=False;visits=0
    parts=sorted(parts,key=lambda p:(p['reference_start'],p['reference_end'],p['protein_id']))
    def extend(chain,index):
        nonlocal truncated,visits
        visits+=1
        if visits>policy['max_assignment_states']:truncated=True;return
        covered=set().union(*(p['reference_positions'] for p in chain));coverage=len(covered)/length
        if coverage>=thresholds['min_reference_coverage']:
            options.append(dict(parts=chain,coverage=coverage,identity=sum(p['identity']*len(p['reference_positions']) for p in chain)/len(covered),
                score=sum(p['identity']*len(p['reference_positions']) for p in chain)/length,fragmented=len(chain)>1))
            return
        if len(options)>policy['max_node_options']:
            truncated=True;return
        for j in range(index,len(parts)):
            if compatible_fragments(chain[-1],parts[j],genes,thresholds,local_cluster):extend(chain+[parts[j]],j+1)
    for i,p in enumerate(parts):extend([p],i+1)
    options.sort(key=lambda o:(-o['score'],o['fragmented'],tuple(p['protein_id'] for p in o['parts'])))
    if len(options)>policy['max_node_options']:truncated=True
    return options[:policy['max_node_options']],truncated


def resolve(anchor,best,system,hits,genes,thresholds,mode,policy,association_allowed,local_cluster):
    graph=next(g for g in policy['systems'] if g['structure_id']==system['structure_id']);nodes=graph['nodes']
    by_reference=defaultdict(list);by_query=defaultdict(list)
    for h in hits:by_reference[h['sseqid']].append(h);by_query[h['qseqid']].append(h)
    pools={};truncated=False
    for node in nodes:
        parts=[]
        for h in ([best] if node['role']==system['catalytic_role'] else by_reference[node['reference_id']]):
            if h['sseqid']!=node['reference_id'] or h['evalue']>thresholds['blast_evalue']:continue
            g=genes[h['qseqid']]
            if not association_allowed(anchor,g,mode):continue
            if node['linkage']=='local' and not _span_within([anchor,g], thresholds,local_cluster):continue
            piece=project(h,node)
            if piece is None or len(piece['reference_positions'])<min(policy['min_fragment_residues'],node['reference_end']-node['reference_start']+1):continue
            if piece['identity']<thresholds['min_identity']:continue
            # Competing references matter only on overlapping residues, allowing unrelated fused modules.
            stronger=False
            for other in by_query[h['qseqid']]:
                if other['qseqid']!=h['qseqid'] or other['sseqid']==h['sseqid'] or other['evalue']>thresholds['blast_evalue']:continue
                if other['bitscore']<=h['bitscore']/thresholds['competitive_score_ratio']:continue
                overlap=sum(other['qstart']<=p<=other['qend'] for p in piece['query_positions'])
                if overlap/len(piece['query_positions'])>=.8:stronger=True;break
            if not stronger:parts.append(piece)
        pools[node['node_id']],limited=realizations(parts,node,genes,thresholds,policy,local_cluster);truncated|=limited
    ordered=sorted(nodes,key=lambda n:(not n['required_for_structure'],len(pools[n['node_id']]),n['node_id']))
    best_solution=None;solutions=[];visits=0
    def signature(selected):
        return tuple((n['node_id'],tuple((p['protein_id'],tuple(sorted(p['query_positions']))) for p in selected[n['node_id']]['parts']))
            for n in nodes if n['required_for_function'] and n['module_category']!='physiological_context' and n['node_id'] in selected)
    def score(selected):
        structural=[selected[n['node_id']] for n in nodes if n['required_for_structure'] and n['node_id'] in selected]
        functional=[selected[n['node_id']] for n in nodes if n['required_for_function'] and n['node_id'] in selected]
        return len(structural),len(functional),sum(o['score'] for o in structural),sum(o['score'] for o in functional)
    def visit(i,selected,resources):
        nonlocal visits,truncated,best_solution,solutions
        visits+=1
        if visits>policy['max_assignment_states']:truncated=True;return
        if best_solution and sum(n['required_for_structure'] for n in ordered[i:])+sum(n['required_for_structure'] and n['node_id'] in selected for n in nodes)<best_solution['score'][0]:return
        if i==len(ordered):
            local_ids={anchor['id']}|{p['protein_id'] for n in nodes if n['required_for_structure'] and n['linkage']=='local' and n['node_id'] in selected for p in selected[n['node_id']]['parts']}
            if not local_cluster([genes[q] for q in local_ids],thresholds):return
            value=dict(selected=dict(selected),score=score(selected),signature=signature(selected))
            if best_solution is None or value['score'][:2]>best_solution['score'][:2]:solutions=[]
            if best_solution is None or value['score'][:2]>=best_solution['score'][:2]:solutions.append(value)
            if best_solution is None or value['score']>best_solution['score']:
                best_solution=value
            return
        node=ordered[i];nid=node['node_id']
        for option in pools[nid]:
            if resource_conflict(option['parts'],resources):continue
            if node['linkage']=='local':
                local_ids={anchor['id']}|{p['protein_id'] for n in nodes if n['linkage']=='local' and n['required_for_structure'] and n['node_id'] in selected for p in selected[n['node_id']]['parts']}
                local_ids.update(p['protein_id'] for p in option['parts'])
                if not _span_within([genes[q] for q in local_ids], thresholds,local_cluster):continue
            selected[nid]=option;visit(i+1,selected,resources+option['parts']);del selected[nid]
        visit(i+1,selected,resources)
    visit(0,{},[])
    if best_solution is None:best_solution=dict(selected={},score=score({}),signature=())
    selected=best_solution['selected'];ambiguous=False;structural_ambiguous=False
    feasible_options=defaultdict(set);competitive_options=defaultdict(set)
    for solution in solutions:
        for nid,option in solution['selected'].items():feasible_options[nid].add(id(option))
    for nid,option in selected.items():competitive_options[nid].add(id(option))
    for alternative in solutions:
        if alternative['signature']==best_solution['signature']:
            for nid,option in alternative['selected'].items():competitive_options[nid].add(id(option))
            continue
        changed=[n for n in nodes if n['module_category']!='physiological_context' and selected.get(n['node_id'])!=alternative['selected'].get(n['node_id'])]
        # Shared catalytic modules must not dilute a discriminating partner-score difference.
        best_score=sum(selected.get(n['node_id'],{}).get('score',0) for n in changed)
        alternative_score=sum(alternative['selected'].get(n['node_id'],{}).get('score',0) for n in changed)
        if changed and alternative_score*thresholds['catalytic_winner_ratio']>=best_score:
            ambiguous=True;structural_ambiguous|=any(n['required_for_structure'] for n in changed)
            for nid,option in alternative['selected'].items():competitive_options[nid].add(id(option))
    missing_structure=[n['node_id'] for n in nodes if n['required_for_structure'] and n['node_id'] not in selected]
    missing_function=[n['node_id'] for n in nodes if n['required_for_function'] and n['node_id'] not in selected]
    chosen=defaultdict(list)
    query_references=defaultdict(set);reference_queries=defaultdict(set)
    for node in nodes:
        for p in selected.get(node['node_id'],{}).get('parts',[]):
            h=p['hit']
            if h not in chosen[node['role']]:chosen[node['role']].append(h)
            query_references[p['protein_id']].add(node['reference_id']);reference_queries[node['reference_id']].add(p['protein_id'])
    fused_queries=sorted(q for q,rr in query_references.items() if len(rr)>1)
    split_references=sorted(r for r,qq in reference_queries.items() if len(qq)>1)
    return dict(nodes=nodes,pools=pools,selected=selected,chosen=dict(chosen),missing_structure=missing_structure,
        missing_function=missing_function,ambiguous=ambiguous,structural_ambiguous=structural_ambiguous,search_truncated=truncated,assignment_states=visits,
        fragmented=any(o['fragmented'] for o in selected.values()) or bool(split_references),fused_queries=fused_queries,
        split_references=split_references,organization_changed=bool(fused_queries or split_references),localization_status=policy['localization_status'],
        feasible_options=feasible_options,competitive_options=competitive_options)


NODE_FIELDS='anchor structure_id node_id role module_category required_for_structure required_for_function status realization protein_ids reference_id reference_start reference_end reference_coverage aligned_pair_identity candidate_count candidate_realizations query_resources source source_scope localization_status graph_search_status ambiguous_assignment cross_anchor_pairing_status cross_anchor_competing_anchors'.split()
EDGE_FIELDS='anchor structure_id node_id protein_id association evidence_scope query_residue_positions reference_residue_positions cross_anchor_pairing_status cross_anchor_competing_anchors'.split()
PAIR_FIELDS='association_id protein_id anchor_a structure_id_a node_a category_a selection_a core_status_a candidate_realizations_a anchor_b structure_id_b node_b category_b selection_b core_status_b candidate_realizations_b overlap_residue_count overlap_query_positions pairing_status blocks_specific_pairing source_a source_b limitation'.split()


def export_graph(anchor,sid,result,genes):
    rows=[];edges=[]
    for node in result['nodes']:
        option=result['selected'].get(node['node_id']);parts=option['parts'] if option else []
        choices=[dict(protein_ids=[p['protein_id'] for p in o['parts']],coverage=o['coverage'],identity=o['identity'],fragmented=o['fragmented'],
            query_resources={p['protein_id']:sorted(p['query_positions']) for p in o['parts']},
            jointly_feasible=id(o) in result['feasible_options'][node['node_id']],
            competitive=id(o) in result['competitive_options'][node['node_id']],selected=o is option) for o in result['pools'][node['node_id']]]
        rows.append(dict(anchor=anchor['id'],structure_id=sid,node_id=node['node_id'],role=node['role'],module_category=node['module_category'],
            required_for_structure=node['required_for_structure'],required_for_function=';'.join(node['required_for_function']),
            status='homology_supported' if option else 'not_detected_or_unassigned',realization='collinear_CDS_fragments' if option and option['fragmented'] else 'single_polypeptide_segment' if option else 'none',
            protein_ids=';'.join(p['protein_id'] for p in parts),reference_id=node['reference_id'],reference_start=node['reference_start'],reference_end=node['reference_end'],
            reference_coverage=option['coverage'] if option else '',aligned_pair_identity=option['identity'] if option else '',candidate_count=len(choices),candidate_realizations=json.dumps(choices),
            query_resources=json.dumps({p['protein_id']:sorted(p['query_positions']) for p in parts}),source=node['source'],source_scope=node['source_scope'],localization_status=result['localization_status']))
        rows[-1].update(graph_search_status='truncated' if result['search_truncated'] else 'complete',ambiguous_assignment=result['ambiguous'])
        for p in parts:
            gene=genes[p['protein_id']]
            association=('same_polypeptide' if p['protein_id']==anchor['id'] else 'local_cluster' if node['linkage']=='local' else
                'same_contig' if gene['contig']==anchor['contig'] else 'same_bin_remote' if gene.get('bin_id') else 'same_genome_remote')
            edges.append(dict(anchor=anchor['id'],structure_id=sid,node_id=node['node_id'],protein_id=p['protein_id'],association=association,
                evidence_scope='Aligned sequence resources and ownership-compatible linkage; transcription/physical interaction not established',
                query_residue_positions=json.dumps(sorted(p['query_positions'])),reference_residue_positions=json.dumps(sorted(p['reference_positions']))))
    return rows,edges


def credible_core(annotation):
    return (annotation['catalytic_status']=='supported' and annotation.get('active_center_span_status')=='supported' and
        annotation.get('local_cluster_complete') in [True,'True'] and annotation.get('module_graph_search_status')=='complete')


def apply_cross_anchor_guard(annotations,nodes,edges,modules,genes,mode,policy,association_allowed):
    """Report alternative physical-association hypotheses without imposing gene copy stoichiometry."""
    anchors={a['protein_id']:a for a in annotations};claims=defaultdict(dict)
    sharing=policy['cross_anchor_sharing']['category_modes']
    default='no_overlapping_fixed_member_claim_detected'
    rank={default:0,'shared_context_cooccurrence':1,'shared_diffusible_pool_candidate':2,
        'provisional_member_overlap_unresolved':3,'overlapping_fixed_member_pairing_unresolved':4}
    node_marks={};anchor_marks={q:dict(status=default,other=set(),proteins=set(),structure=False,functions=set()) for q in anchors}
    for node in nodes:
        qid=node['anchor'];selected=json.loads(node['query_resources'])
        for index,choice in enumerate(json.loads(node['candidate_realizations'])):
            if not choice['jointly_feasible'] or not choice['competitive']:continue
            for protein,positions in choice['query_resources'].items():
                resource=frozenset(positions)
                if (not resource or not all(isinstance(p,int) and not isinstance(p,bool) and 1<=p<=len(genes[protein]['sequence']) for p in resource) or
                    not association_allowed(genes[qid],genes[protein],mode)):
                    raise ValueError('Invalid cross-anchor candidate coordinates or input ownership')
                key=(qid,node['node_id'],tuple(sorted(resource)))
                claim=claims[protein].setdefault(key,dict(node=node,positions=resource,selected=False,realizations=[]))
                claim['selected']|=choice['selected'] and set(selected.get(protein,[]))==resource
                claim['realizations'].append(index)
    rows=[]
    for protein,by_key in sorted(claims.items()):
        for (ka,a),(kb,b) in combinations(sorted(by_key.items()),2):
            na,nb=a['node'],b['node'];qa,qb=na['anchor'],nb['anchor']
            overlap=a['positions']&b['positions']
            if qa==qb or not overlap:continue
            ma,mb=sharing[na['module_category']],sharing[nb['module_category']]
            core_a,core_b=credible_core(anchors[qa]),credible_core(anchors[qb])
            if 'cooccurrence_only' in [ma,mb]:status='shared_context_cooccurrence';blocking=False
            elif ma==mb=='diffusible_pool_candidate':status='shared_diffusible_pool_candidate';blocking=False
            elif core_a and core_b:status='overlapping_fixed_member_pairing_unresolved';blocking=True
            else:status='provisional_member_overlap_unresolved';blocking=False
            digest=hashlib.sha256(json.dumps([protein,ka,kb,sorted(overlap)]).encode()).hexdigest()[:16]
            rows.append(dict(association_id='association_'+digest,protein_id=protein,
                anchor_a=qa,structure_id_a=na['structure_id'],node_a=na['node_id'],category_a=na['module_category'],
                selection_a='selected' if a['selected'] else 'competitive_alternative',core_status_a='supported_local_core' if core_a else 'provisional_or_partial_core',
                candidate_realizations_a=json.dumps(a['realizations']),
                anchor_b=qb,structure_id_b=nb['structure_id'],node_b=nb['node_id'],category_b=nb['module_category'],
                selection_b='selected' if b['selected'] else 'competitive_alternative',core_status_b='supported_local_core' if core_b else 'provisional_or_partial_core',
                candidate_realizations_b=json.dumps(b['realizations']),
                overlap_residue_count=len(overlap),overlap_query_positions=json.dumps(sorted(overlap)),pairing_status=status,
                blocks_specific_pairing=blocking,source_a=na['source'],source_b=nb['source'],
                limitation='Sequence-supported association hypotheses, not physical exclusivity, gene-copy stoichiometry, expression or measured flux. Provisional competing cores do not invalidate a supported core solely by this overlap.'))
            for node,other in [(na,qb),(nb,qa)]:
                qid=node['anchor'];key=(qid,node['node_id']);mark=node_marks.setdefault(key,dict(status=default,other=set()))
                if rank[status]>rank[mark['status']]:mark['status']=status
                mark['other'].add(other);summary=anchor_marks[qid]
                if rank[status]>rank[summary['status']]:summary['status']=status
                summary['other'].add(other);summary['proteins'].add(protein)
                if blocking:
                    summary['structure']|=node['required_for_structure'] in [True,'True']
                    required=node['required_for_function']
                    summary['functions'].update(fid for fid in (required.split(';') if isinstance(required,str) else required) if fid)
    for annotation in annotations:
        qid=annotation['protein_id'];mark=anchor_marks[qid]
        annotation.update(cross_anchor_pairing_status=mark['status'],cross_anchor_competing_anchors=';'.join(sorted(mark['other'])),
            cross_anchor_partner_proteins=';'.join(sorted(mark['proteins'])),cross_anchor_function_ambiguity_ids=';'.join(sorted(mark['functions'])),
            physical_interaction_status='not_established_by_sequence')
        if mark['status'] in ['provisional_member_overlap_unresolved','overlapping_fixed_member_pairing_unresolved']:
            flags=annotation.get('flags','').split(';');flags.append('cross_anchor_'+mark['status']);annotation['flags']=';'.join(dict.fromkeys(f for f in flags if f))
        if mark['structure']:
            if annotation['structure_confidence']=='high':annotation['structure_confidence']='medium'
            if annotation['structure_status'] in ['complete_supported','core_complete_accessory_missing','composition_supported_reorganized']:
                annotation['structure_status']='shared_member_pairing_unresolved'
            annotation['small_subunit_topology']='unresolved'
        if mark['functions']:
            annotation['potential_function_ids']=';'.join(fid for fid in annotation.get('potential_function_ids','').split(';') if fid and fid not in mark['functions'])
            if not annotation['potential_function_ids']:
                annotation.update(function_status='unassigned',function_confidence='unassigned',potential_direct_partner='unresolved')
    for node in nodes:
        mark=node_marks.get((node['anchor'],node['node_id']),dict(status=default,other=set()))
        node.update(cross_anchor_pairing_status=mark['status'],cross_anchor_competing_anchors=';'.join(sorted(mark['other'])))
    for row in [*edges,*modules]:
        key=(row['anchor'],row.get('node_id',row.get('module_node_id')));mark=node_marks.get(key,dict(status=default,other=set()))
        row.update(cross_anchor_pairing_status=mark['status'],cross_anchor_competing_anchors=';'.join(sorted(mark['other'])))
    return rows
