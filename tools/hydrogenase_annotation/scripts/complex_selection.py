"""Compare feasible L4 hypotheses after identity confirmation; never promote a reaction."""
import json
from collections import defaultdict
from annotate import classify, association_allowed, local_cluster
from localization import sha
from module_graph import load_policy, resolve

FIELDS = ('anchor structure_id anchor_bitscore unresolved_competitor_bitscore anchor_outscores_unresolved required_member_ids composition_complete missing_structure '
          'assignment_ambiguous search_truncated organization_changed partial_members decision '
          'reference_function_ids reference_contract_ids composition_complete_before_architecture architecture_status architecture_evidence anchor_competition_guard reviewed_competitor_compatibility previous_structure_id previous_structure_flags evidence_scope').split()
FIELDS+=['fusion_composition_status','fusion_composition_evidence']


def maximal_cores(candidates):
    complete = [c for c in candidates if c['composition_complete']]
    return [c for c in complete if not any(set(c['required_member_ids']) < set(d['required_member_ids']) for d in complete)]


def fusion_composition(graph,sid,anchor,feature,center,genes,policy,review,refs,competitors,score):
    spec=policy.get('fusion_composition',{})
    if sid!=spec.get('structure_id'):return 'outside_reviewed_fusion',{}
    from annotate import ROOT
    from catalytic_confirmation import qualifies
    rule=json.loads((ROOT/'config/catalytic_confirmation_thresholds.json').read_text())['families']['FeFe']
    if (anchor['catalytic_family']!='FeFe' or anchor.get('active_center_span_status')!='supported' or
        feature.get('family')!='FeFe' or feature.get('eligible') is not False or
        feature.get('eligibility_reasons')!=['reviewed_single_chain_center_or_span_not_supported'] or
        not qualifies(dict(feature,eligible=True),rule) or
        center.get('pattern_status')!='literature_ligand_pattern_candidate' or
        center.get('context_alignment_status')!='supported' or center.get('template_span_status')!='supported' or
        center.get('confidence')!='uncalibrated_literature_hypothesis_only' or
        graph['missing_structure'] or graph['structural_ambiguous'] or graph['search_truncated'] or graph['fragmented'] or graph['split_references']):
        return 'fusion_requirements_not_met',{}
    required=[n for n in graph['nodes'] if n['required_for_structure']]
    if {n['role'] for n in required}!=set(spec['reference_role_hashes']):return 'unreviewed_fusion_roles',{}
    parts={};details=[]
    for node in required:
        selected=graph['selected'][node['node_id']]
        if (len(selected['parts'])!=1 or selected['coverage']<spec['minimum_role_reference_coverage'] or
            refs[node['reference_id']]['sequence_sha256']!=spec['reference_role_hashes'][node['role']]):return 'fusion_role_not_complete',{}
        part=selected['parts'][0]
        if genes[part['protein_id']]['partial']:return 'partial_fusion_member',{}
        parts[node['role']]=part
        details.append(dict(role=node['role'],protein_id=part['protein_id'],reference_id=node['reference_id'],
            reference_coverage=selected['coverage'],query_positions=sorted(part['query_positions'])))
    first,second=[parts[r] for r in spec['fused_roles_in_query_order']];catalytic=parts[spec['catalytic_role']]
    if (first['protein_id']!=second['protein_id'] or first['protein_id']==catalytic['protein_id'] or
        graph['fused_queries']!=[first['protein_id']] or max(first['query_positions'])>=min(second['query_positions']) or
        first['query_positions']&second['query_positions'] or
        len(first['query_positions']|second['query_positions'])/len(genes[first['protein_id']]['sequence'])<spec['minimum_fused_query_coverage']):
        return 'fusion_partition_not_supported',{}
    positions=center['query_ligand_positions'];positions=json.loads(positions) if isinstance(positions,str) else positions
    if not set(positions)<=catalytic['query_positions']:return 'center_outside_selected_catalytic_resources',{}
    compatible=[];identity=review[0]['identity_positive_competitors']
    for hit in competitors:
        if hit['bitscore']<score:continue
        ref=refs[hit['sseqid']];source=identity.get(hit['sseqid'])
        if not source or ref['kind']!='unresolved' or source['sequence_sha256']!=ref['sequence_sha256']:
            return 'fusion_blocked_by_unreviewed_or_negative_competitor',{}
        compatible.append(dict(reference_id=ref['reference_id'],sequence_sha256=ref['sequence_sha256'],bitscore=hit['bitscore'],identity_label='positive'))
    return 'reviewed_noncanonical_fusion_composition_candidate',dict(parts=details,center_model_id=center['candidate_model_id'],
        center_query_positions=positions,center_confidence=center['confidence'],sequence_guard=rule,
        primary_source=spec['primary_source'],policy_sha256=sha(ROOT/'config/complex_selection.json'),
        sequence_guard_scope='Necessary sequence support only; canonical-center failure is retained, no high-identity calibration inherited for this exception.',
        compatible_competitors=compatible,scope=spec['scope'])


def refine(genes, hits, references, systems, thresholds, mode, anchors, policy, graph_policy=None, confirmation_features=None, center_candidates=None):
    """Return L4-only patches and replacement graph rows; inputs remain unchanged."""
    if mode not in policy['eligible_input_modes']: return {}, {}, []
    from annotate import ROOT
    if graph_policy is None: graph_policy = load_policy(ROOT, references)
    from module_architecture import load_review,explain,competitor_guard
    architecture_review=load_review()
    features={f['protein_id']:f for f in confirmation_features or []};centers={c['protein_id']:c for c in center_candidates or [] if c.get('pattern_status')=='literature_ligand_pattern_candidate'}
    refs = {r['reference_id']: r for r in references}; by_gene = {g['id']: g for g in genes}
    by_system = {s['structure_id']: s for s in systems}; by_query = defaultdict(list); unresolved_scores = defaultdict(list)
    for h in hits:
        if h['passes_thresholds'] and refs[h['sseqid']]['kind'] == 'catalytic': by_query[h['qseqid']].append(h)
        if h['passes_thresholds'] and refs[h['sseqid']]['kind'] in ['unresolved', 'decoy']:
            unresolved_scores[h['qseqid']].append(h)
    patches, replacements, evidence = {}, {}, []
    for anchor in anchors:
        qid = anchor['protein_id']
        if (anchor['structure_status'] in policy['preserve_statuses'] or
            anchor['catalytic_family'] not in policy['eligible_families'] or
            (anchor.get('active_center_status') != policy['required_active_center_status'] and qid not in centers) or
            anchor.get('active_center_span_status') != policy['required_active_span_status']): continue
        candidates, contexts = [], {}
        for hit in by_query[qid]:
            ref = refs[hit['sseqid']]; sid = ref['structure_id']; system = by_system[sid]
            if system['center'] != anchor['catalytic_family']: continue
            # A generic shared module cannot be vetoed by a different, potentially infeasible system in advance.
            scoped_hits = [h for h in hits if refs[h['sseqid']]['structure_id'] == sid and
                (refs[h['sseqid']]['kind'] != 'catalytic' or h['qseqid'] == qid)]
            graph = resolve(by_gene[qid], hit, system, scoped_hits, by_gene, thresholds, mode, graph_policy, association_allowed, local_cluster)
            members = sorted({p['protein_id'] for n in graph['nodes'] if n['required_for_structure'] and n['node_id'] in graph['selected']
                for p in graph['selected'][n['node_id']]['parts']})
            partial = [q for q in members if by_gene[q]['partial']]
            complete = not graph['missing_structure'] and not graph['structural_ambiguous'] and not graph['search_truncated']
            architecture_status,architecture_evidence=explain(graph,sid,hits,by_gene,architecture_review)
            before_architecture=complete
            if architecture_status in ['unexplained_reviewed_FeS_extension','reviewed_signature_supported_scaffold_ligands_incomplete']:complete=False
            guard,compatible=competitor_guard(anchor,sid,hit['bitscore'],unresolved_scores[qid],refs,architecture_status,architecture_review)
            fusion_status,fusion_evidence=fusion_composition(graph,sid,anchor,features.get(qid,{}),centers.get(qid,{}),by_gene,policy,architecture_review,refs,unresolved_scores[qid],hit['bitscore'])
            if fusion_status=='reviewed_noncanonical_fusion_composition_candidate':
                guard='sequence_supported_noncanonical_fusion_composition';compatible=fusion_evidence['compatible_competitors']
            competitor = max((h['bitscore'] for h in unresolved_scores[qid]), default=0)
            row = dict(anchor=qid, structure_id=sid, anchor_bitscore=hit['bitscore'], unresolved_competitor_bitscore=competitor,
                anchor_outscores_unresolved=hit['bitscore'] > competitor, required_member_ids=members,
                composition_complete=complete,composition_complete_before_architecture=before_architecture,architecture_status=architecture_status,architecture_evidence=architecture_evidence,anchor_competition_guard=guard,reviewed_competitor_compatibility=compatible, missing_structure=graph['missing_structure'], assignment_ambiguous=graph['structural_ambiguous'],
                search_truncated=graph['search_truncated'], organization_changed=graph['organization_changed'], partial_members=partial,
                decision='incomplete_or_ambiguous', reference_function_ids=system['function_ids'],
                reference_contract_ids=system.get('reference_contract_ids',system['function_ids']),
                previous_structure_id=anchor['structure_id'], previous_structure_flags=anchor.get('flags',''),fusion_composition_status=fusion_status,fusion_composition_evidence=fusion_evidence,
                evidence_scope='Alternative reference composition only; not physical assembly, stoichiometry, activity or a query reaction.')
            candidates.append(row); contexts[sid] = scoped_hits
        maximal = maximal_cores(candidates)
        eligible = (len(maximal) == 1 and len(maximal[0]['required_member_ids']) >= policy['minimum_rescue_polypeptides'] and
            maximal[0]['anchor_competition_guard'] in ['strictly_outscores_unresolved','independent_identity_and_identity_positive_competitors','sequence_supported_noncanonical_fusion_composition'] and
            (not maximal[0]['organization_changed'] or maximal[0]['fusion_composition_status']=='reviewed_noncanonical_fusion_composition_candidate') and not maximal[0]['partial_members'] and
            (anchor.get('active_center_status') == policy['required_active_center_status'] or maximal[0]['fusion_composition_status']=='reviewed_noncanonical_fusion_composition_candidate') and
            anchor.get('active_center_span_status') == policy['required_active_span_status'])
        for row in candidates:
            if row['composition_complete']:
                row['decision'] = 'competing_complete_core' if row in maximal else 'complete_subset_not_discriminating'
        if eligible:
            chosen = maximal[0]; sid = chosen['structure_id']
            result, modules, _, nodes, edges = classify(genes, contexts[sid], [], references, systems, thresholds, mode, graph_policy=graph_policy)
            candidate = next((a for a in result if a['protein_id'] == qid), None)
            if candidate and candidate['structure_status'] in policy['preserve_statuses']:
                # The caller enforces the patch key allowlist, so this stage cannot replace L1/L2 decisions.
                keys = ('structure_id structure_confidence function_confidence potential_function_ids structure_status function_status '
                        'missing_roles reference_clade clade_method reference_hyddb_label reference_ko_label '
                        'diagnostic_positions_conserved diagnostic_residue_map catalytic_identity catalytic_query_coverage catalytic_reference_coverage '
                        'local_cluster_complete missing_context_roles small_subunit_topology potential_direct_partner reference_mechanism '
                        'reference_energy_coupling operon_status module_graph_search_status module_graph_assignment_states '
                        'module_organization_status missing_optional_roles structure_completeness_scope localization_status evidence_url flags').split()
                patch = {k: candidate[k] for k in keys}
                patch.update(structure_confidence='medium', function_confidence='unassigned', potential_function_ids='', function_status='unassigned',
                    small_subunit_topology='unresolved', potential_direct_partner='unresolved',
                    clade_method='reference_module_graph_hypothesis; not a phylogenetic placement',
                    flags=';'.join(dict.fromkeys(filter(None, [*[f for f in anchor.get('flags','').split(';') if f.startswith(('identity_confirmation_','confirmation_'))], *candidate.get('flags', '').split(';'),
                        'alternative_complete_core_provisional_function_uncalibrated']))))
                if chosen['fusion_composition_status']=='reviewed_noncanonical_fusion_composition_candidate':
                    patch['module_organization_status']='reviewed_fusion_composition_candidate'
                    patch['flags']+=';noncanonical_fusion_composition_does_not_confirm_catalysis'
                    if 'complete_active_center_requires_review' in anchor.get('flags','').split(';'):
                        patch['flags']+=';complete_active_center_requires_review'
                catalytic_ref=next(r for r in references if r['structure_id']==sid and r['kind']=='catalytic')
                if not catalytic_ref['diagnostic_positions']:
                    patch['diagnostic_positions_conserved']='not_assessed_by_legacy_seed_subset'
                    patch['flags']+=';no_new_active_center_seed_in_structure_extension'
                patches[qid] = patch; replacements[qid] = (modules, nodes, edges)
                chosen['decision'] = 'selected_provisional_core_no_function_transfer'
        evidence.extend(candidates)
    return patches, replacements, evidence
