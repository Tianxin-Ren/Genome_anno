#!/usr/bin/env python3
"""Verify evidence schemas, provenance, six layers and coordinates in actual runs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, COMPLEX_FIELDS, MODULE_FIELDS, REACTION_FIELDS, annotate, extract_genbank, write_tsv
from fes_sites import SITE_FIELDS,site_evidence,full_reference_coordinates
from active_centers import SITE_FIELDS as ACTIVE_SITE_FIELDS, SUMMARY_FIELDS as ACTIVE_SUMMARY_FIELDS
from discovery_evidence import FIELDS as DISCOVERY_FIELDS,chain_compatible
from module_graph import NODE_FIELDS, EDGE_FIELDS, PAIR_FIELDS
from split_centers import CANDIDATE_FIELDS, COMPONENT_FIELDS
from catalytic_confirmation import CONFIRMATION_FIELDS,qualifies,stage_fingerprint as confirmation_fingerprint
from complex_selection import FIELDS as HYPOTHESIS_FIELDS, maximal_cores
from complex_extensions import CONTRACT_FIELDS,load_bank,function_reviews
from module_architecture import load_review as load_architecture_review
from annotate import load_blast,local_cluster,residue_states,repeat_indel_ambiguity,repeat_mapping_ambiguity
from circular_coordinates import valid_cds_interval


def traced_fes_states_match(recorded,hits,model):
    positions=[int(p) for p in model['expected_ligands']]
    return any(residue_states(h,positions,model['expected_ligands'])[0]==recorded for h in hits)


def traced_active_states(hits,model,query):
    mapped={};ambiguous=set()
    for hit in hits:
        _,piece=residue_states(hit,range(model['window_start'],model['window_end']+1))
        if set(mapped)&set(piece):raise ValueError('Overlapping traced active-center HSPs')
        mapped.update(piece);ambiguous.update(repeat_indel_ambiguity(hit))
    ambiguous.update(repeat_mapping_ambiguity(mapped,model['sequence'],query));result={}
    for site in model['sites']:
        p=site['reference_position'];q,residue=mapped.get(p,(None,''))
        state=('ambiguous' if p in ambiguous else 'not_covered' if residue in ['','-'] else 'ambiguous' if residue in 'XBZJO' else 'present' if residue in site['expected_residue'] else 'substituted')
        result[p]=dict(status=state,query_position=None if p in ambiguous else q,observed_residue=residue)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--round',type=Path,required=True)
    args=parser.parse_args();directory=args.round.resolve();checks=[]
    schemas={name:COMPLEX_FIELDS for name in ['complexes.tsv','catalytic_anchors.tsv']}
    schemas.update({name:MODULE_FIELDS for name in ['modules.tsv','complex_members.tsv']})
    schemas.update({name:REACTION_FIELDS for name in ['potential_functions.tsv','reaction_capabilities.tsv']})
    schemas.update({
        'catalytic_sites.tsv':ACTIVE_SITE_FIELDS,
        'active_center_summary.tsv':ACTIVE_SUMMARY_FIELDS,
        'fes_sites.tsv':SITE_FIELDS,
        'maturation.tsv':['protein_id','maturation_gene','reference','source','reference_domain_validation','query_domain_validation','cofactor_assembly'],
        'annotation_evidence.tsv':['protein_id','complex_anchor','layer','criterion','status','reference_id','source','details','limitation']})
    schemas['family_candidates.tsv']=DISCOVERY_FIELDS
    schemas['module_graph_nodes.tsv']=NODE_FIELDS;schemas['module_graph_edges.tsv']=EDGE_FIELDS
    schemas['cross_anchor_associations.tsv']=PAIR_FIELDS
    schemas['catalytic_center_candidates.tsv']=CANDIDATE_FIELDS
    schemas['catalytic_center_components.tsv']=COMPONENT_FIELDS
    schemas['catalytic_confirmation.tsv']=CONFIRMATION_FIELDS
    schemas['complex_hypotheses.tsv']=HYPOTHESIS_FIELDS
    schemas['reaction_contract_candidates.tsv']=CONTRACT_FIELDS
    def check(name,passed,observed):checks.append(dict(test=name,passed=bool(passed),observed=json.dumps(observed)))
    refs=json.loads((ROOT/'references/processed/references.json').read_text())
    function_review={r['function_id']:r for r in json.loads((ROOT/'config/reference_function_review.json').read_text())['records']}
    graph_spec={(g['structure_id'],n['node_id']):n for g in json.loads((ROOT/'config/module_graph.json').read_text())['systems'] for n in g['nodes']}
    extension_refs,extension_systems,extension_graphs,extension_model,contracts=load_bank()
    graph_spec.update({(g['structure_id'],n['node_id']):n for g in extension_graphs for n in g['nodes']})
    function_review=function_reviews(contracts)
    architecture_review,architecture_sites,architecture_parameters=load_architecture_review()
    identity_review={r['reference_id']:r for r in json.loads((ROOT/'config/reference_catalytic_identity_review.json').read_text())['records']}
    contracts={c['contract_id']:c for c in contracts}
    decoy=next(r for r in refs if r['kind']=='decoy')
    source=ROOT/'references/raw'/f"{decoy['sample']}.gb"
    gene=next(g for g in extract_genbank(source)[0] if g['id']==decoy['locus_tag'])
    empty=directory/'output_contract_empty_case'
    result=annotate([gene],empty,input_files=[source],mode='genome')
    check('real_decoy_has_no_catalytic_anchor',len(result)==0,len(result))
    runs=[p for suite in ['regression','input_safety','module_graph','partner_associations','split_centers','catalytic_confirmation','paralog_contexts','natural_complexes','ralbus_context','function_mechanism','distributed_function'] for p in sorted((directory/suite).rglob('run_manifest.json'))]+[empty/'run_manifest.json']
    modes=set();empty_cases=0
    models={r['site_id']:r for r in json.loads((ROOT/'references/processed/fes/sites.json').read_text())}
    window_names={v:k for k,v in json.loads((ROOT/'references/processed/fes/site_window_index.json').read_text()).items()}
    alignment_cache={};active_trace_cache={}
    def actual_hits(run,filename,qid,rid,qs,qe,score):
        key=(str(run),filename)
        if key not in alignment_cache:
            indexed={}
            for h in load_blast(run/filename,dict(blast_evalue=float('inf'),min_identity=0,min_query_coverage=0,min_reference_coverage=0)):
                indexed.setdefault((h['qseqid'],h['sseqid'],h['qstart'],h['qend'],h['bitscore']),[]).append(h)
            alignment_cache[key]=indexed
        return alignment_cache[key].get((qid,rid,int(qs),int(qe),float(score)),[])
    def actual_active_states(run,site,model,query):
        raw=site['alignment_components'];parts=json.loads(raw) if isinstance(raw,str) else raw
        key=(str(run),site['protein_id'],site['reference_id'],json.dumps(parts,sort_keys=True))
        if key not in active_trace_cache:
            filename='active_center_blast.tsv' if site['reference_id'] in canonical_active_models else 'split_center_blast.tsv';hits=[]
            for part in parts:
                choices=[h for h in actual_hits(run,filename,site['protein_id'],site['reference_id'],part['query_start'],part['query_end'],part['bitscore'])
                    if h['sstart']==part['reference_start'] and h['send']==part['reference_end'] and h['evalue']==part['evalue']]
                if len(choices)!=1:return {}
                hits.append(choices[0])
            active_trace_cache[key]=traced_active_states(hits,model,query) if hits else {}
        return active_trace_cache[key]
    canonical_active_models={m['reference_id'] for m in json.loads((ROOT/'references/processed/active_centers/models.json').read_text())}
    for path in runs:
        run=path.parent;label=str(run.relative_to(directory));manifest=json.loads(path.read_text());tables={}
        check(label+':completed',manifest['status']=='complete',manifest['status'])
        if manifest['status']!='complete':continue
        modes.add(manifest['input_mode'])
        for name,fields in schemas.items():
            with (run/name).open() as handle:
                reader=csv.DictReader(handle,delimiter='\t');tables[name]=list(reader)
                missing=set(fields)-set(reader.fieldnames or [])
            check(label+':schema:'+name,not missing,sorted(missing))
        if not tables['catalytic_anchors.tsv']:empty_cases+=1
        check(label+':legacy_aliases',all((run/old).read_bytes()==(run/new).read_bytes() for old,new in [
            ('complexes.tsv','catalytic_anchors.tsv'),('modules.tsv','complex_members.tsv'),('potential_functions.tsv','reaction_capabilities.tsv')]),'byte-identical compatible aliases')
        bad=[]
        for row in tables['reaction_capabilities.tsv']:
            review=function_review[row['function_id']]
            anchor=next(a for a in tables['catalytic_anchors.tsv'] if a['protein_id']==row['complex_anchor'])
            if (row['query_physiology']!='unresolved' or row['query_net_direction']!='unresolved' or
                row['function_id'] not in anchor['potential_function_ids'].split(';') or row['function_id'] in anchor['cross_anchor_function_ambiguity_ids'].split(';') or
                row['physical_interaction_status']!='not_established_by_sequence' or row['association_pairing_status']!=anchor['cross_anchor_pairing_status'] or
                any(row[k]!=review[k] for k in ['molecular_reaction','reference_evidence_scope','capability_scope','reference_context','limitations'])):
                bad.append(row['function_id'])
        check(label+':bounded_reference_function_evidence',not bad,bad)
        hashes=manifest.get('output_sha256',{})
        required_hashes=set(schemas)
        bad=[name for name in required_hashes if hashes.get(name)!=hashlib.sha256((run/name).read_bytes()).hexdigest()]
        check(label+':output_hashes',not bad,bad)
        software=manifest.get('software',{})
        check(label+':software_parameters',all(name in software for name in ['python','biopython','blastp','hmmsearch','prodigal','mafft','iqtree3']) and
            all(k in manifest for k in ['threads','translation_table','thresholds','input_mode','fes_site_policy_sha256']),sorted(software))
        auxiliary=manifest.get('noncanonical_center_evidence',{})
        aux_path=run/'noncanonical_centers/alternatives.tsv'
        with aux_path.open() as handle:
            reader=csv.DictReader(handle,delimiter='\t');aux_rows=list(reader)
            from noncanonical_center_evidence import FIELDS as AUX_FIELDS
            aux_schema=set(AUX_FIELDS)<=set(reader.fieldnames or [])
        ids={r['protein_id'] for r in tables['catalytic_anchors.tsv'] if r['catalytic_family']=='FeFe'}
        check(label+':noncanonical_literature_evidence',aux_schema and auxiliary.get('changes_annotation_decisions') is False and
            auxiliary.get('policy_sha256')==hashlib.sha256((ROOT/'config/noncanonical_center_evidence.json').read_bytes()).hexdigest() and
            auxiliary.get('script_sha256')==hashlib.sha256((ROOT/'scripts/noncanonical_center_evidence.py').read_bytes()).hexdigest() and
            auxiliary.get('output_sha256')==hashlib.sha256(aux_path.read_bytes()).hexdigest() and
            all(r['protein_id'] in ids and r['confidence']=='uncalibrated_literature_hypothesis_only' and
                r['query_activity']=='not_measured' and r['cofactor_occupancy']=='unknown' for r in aux_rows) and
            not any(g['complete_primary_ligand_set'] for g in auxiliary.get('source_geometry',[])),len(aux_rows))
        expected_hashes={'discovery_evidence_policy_sha256':ROOT/'config/discovery_evidence.json',
            'discovery_evidence_script_sha256':ROOT/'scripts/discovery_evidence.py',
            'discovery_profile_manifest_sha256':ROOT/'references/processed/profile_manifest.json',
            'active_center_model_manifest_sha256':ROOT/'references/processed/active_centers/manifest.json',
            'active_center_policy_sha256':ROOT/'config/active_center_review.json',
            'active_center_script_sha256':ROOT/'scripts/active_centers.py',
            'split_center_policy_sha256':ROOT/'config/split_centers.json',
            'split_center_model_manifest_sha256':ROOT/'references/processed/split_centers/manifest.json',
            'split_center_script_sha256':ROOT/'scripts/split_centers.py',
            'module_graph_policy_sha256':ROOT/'config/module_graph.json','module_graph_script_sha256':ROOT/'scripts/module_graph.py'}
        check(label+':discovery_provenance',all(manifest.get(k)==hashlib.sha256(p.read_bytes()).hexdigest() for k,p in expected_hashes.items()),sorted(expected_hashes))
        selection=manifest['complex_selection'];bad=[]
        if (selection['policy_sha256']!=hashlib.sha256((ROOT/'config/complex_selection.json').read_bytes()).hexdigest() or
            selection['script_sha256']!=hashlib.sha256((ROOT/'scripts/complex_selection.py').read_bytes()).hexdigest()):bad.append('selection_source_changed')
        grouped={}
        anchor_index={a['protein_id']:a for a in tables['catalytic_anchors.tsv']}
        for row in tables['complex_hypotheses.tsv']:
            row=dict(row,composition_complete=row['composition_complete']=='True',required_member_ids=json.loads(row['required_member_ids']))
            grouped.setdefault(row['anchor'],[]).append(row)
            if row['anchor'] not in anchor_index or row['structure_id'] not in {s[0] for s in graph_spec}:bad.append('hypothesis_ownership')
        selected=[]
        for q,rows in grouped.items():
            maxima=maximal_cores(rows)
            for row in rows:
                if row['decision']!='selected_provisional_core_no_function_transfer':continue
                a=anchor_index[q];selected.append(q)
                if (maxima!=[row] or len(row['required_member_ids'])<json.loads((ROOT/'config/complex_selection.json').read_text())['minimum_rescue_polypeptides'] or
                    row['anchor_competition_guard'] not in ['strictly_outscores_unresolved','independent_identity_and_identity_positive_competitors','sequence_supported_noncanonical_fusion_composition'] or json.loads(row['partial_members']) or
                    (row['organization_changed']!='False' and row['fusion_composition_status']!='reviewed_noncanonical_fusion_composition_candidate') or a['structure_id']!=row['structure_id'] or
                    a['structure_confidence']=='high' or (a['potential_function_ids'] and q not in manifest['exact_reference_functions']['assigned_anchors']) or
                    a['function_status'] not in ['unassigned','reference_sequence_supported_potential'] or
                    a['clade_method']!='reference_module_graph_hypothesis; not a phylogenetic placement'):
                    bad.append(q+':invalid_L4_rescue_or_unvalidated_function')
        if sorted(selected)!=selection['rescued_anchors']:bad.append('selection_manifest_ownership')
        if (selection['module_architecture']['policy_sha256']!=hashlib.sha256((ROOT/'config/module_architecture.json').read_bytes()).hexdigest() or
            selection['module_architecture']['script_sha256']!=hashlib.sha256((ROOT/'scripts/module_architecture.py').read_bytes()).hexdigest()):bad.append('architecture_provenance')
        raw_hits=load_blast(run/'blast.tsv',manifest['thresholds'],refs)
        reference_index={r['reference_id']:r for r in refs}
        for q,rows in grouped.items():
            for row in rows:
                architecture=json.loads(row['architecture_evidence']);compatibility=json.loads(row['reviewed_competitor_compatibility'])
                if row['composition_complete_before_architecture']=='True' and row['architecture_status'] in ['unexplained_reviewed_FeS_extension','reviewed_signature_supported_scaffold_ligands_incomplete'] and row['composition_complete']:
                    bad.append(q+':unexplained_architecture_counted_complete')
                if row['anchor_competition_guard']=='strictly_outscores_unresolved' and (row['anchor_outscores_unresolved']!='True' or float(row['anchor_bitscore'])<=float(row['unresolved_competitor_bitscore'])):
                    bad.append(q+':false_strict_competitor_win')
                if row['anchor_competition_guard']=='sequence_supported_noncanonical_fusion_composition':
                    proof=json.loads(row['fusion_composition_evidence'])
                    feature=next(f for f in json.loads((run/'confirmation/features.json').read_text()) if f['protein_id']==q)
                    center=next((r for r in aux_rows if r['protein_id']==q),{})
                    spec=json.loads((ROOT/'config/complex_selection.json').read_text())['fusion_composition']
                    parts={r['role']:r for r in proof.get('parts',[])}
                    stronger=[h for h in raw_hits if h['qseqid']==q and h['passes_thresholds'] and h['bitscore']>=float(row['anchor_bitscore']) and reference_index[h['sseqid']]['kind'] in ['unresolved','decoy']]
                    if (row['structure_id']!=spec['structure_id'] or row['fusion_composition_status']!='reviewed_noncanonical_fusion_composition_candidate' or
                        feature['eligible'] or feature['eligibility_reasons']!=['reviewed_single_chain_center_or_span_not_supported'] or
                        not qualifies(dict(feature,eligible=True),spec['sequence_guard']) or
                        center.get('pattern_status')!='literature_ligand_pattern_candidate' or set(parts)!=set(spec['reference_role_hashes']) or
                        proof.get('center_query_positions')!=json.loads(center.get('query_ligand_positions','[]')) or
                        proof.get('center_confidence')!=center.get('confidence') or proof.get('primary_source')!=spec['primary_source'] or
                        {c['reference_id'] for c in compatibility}!={h['sseqid'] for h in stronger} or
                        any(identity_review.get(h['sseqid'],{}).get('identity_label')!='positive' for h in stronger)):
                        bad.append(q+':invalid_noncanonical_fusion_evidence')
                    if row['decision']=='selected_provisional_core_no_function_transfer':
                        a=anchor_index[q]
                        if (a['catalytic_identity_confidence']=='calibrated_high' or a['structure_confidence']!='medium' or a['potential_function_ids'] or
                            a['module_organization_status']!='reviewed_fusion_composition_candidate' or a['active_center_status']=='template_ligand_set_supported'):
                            bad.append(q+':fusion_candidate_promoted_identity_or_function')
                    for role,part in parts.items():
                        node=next(n for n in tables['module_graph_nodes.tsv'] if n['anchor']==q and n['structure_id']==row['structure_id'] and n['role']==role) if row['decision']=='selected_provisional_core_no_function_transfer' else None
                        if node and json.loads(node['query_resources'])!={part['protein_id']:part['query_positions']}:
                            bad.append(q+':fusion_proof_differs_from_selected_resources')
                if row['anchor_competition_guard']=='independent_identity_and_identity_positive_competitors':
                    stronger=[h for h in raw_hits if h['qseqid']==q and h['passes_thresholds'] and h['bitscore']>=float(row['anchor_bitscore']) and reference_index[h['sseqid']]['kind'] in ['unresolved','decoy']]
                    if (row['structure_id'] not in architecture_review['reconciled_structure_ids'] or anchor_index[q]['catalytic_identity_confidence']!='calibrated_high' or
                        row['architecture_status']!='selected_resources_explain_reviewed_signature' or not stronger or
                        {h['sseqid'] for h in stronger}!={c['reference_id'] for c in compatibility} or any(
                        reference_index[h['sseqid']]['kind']!='unresolved' or identity_review.get(h['sseqid'],{}).get('identity_label')!='positive' or
                        identity_review[h['sseqid']]['sequence_sha256']!=reference_index[h['sseqid']]['sequence_sha256'] for h in stronger)):
                        bad.append(q+':unreviewed_identity_competitor_reconciliation')
                for ar in architecture:
                    p=ar['protein_id'];choices=[h for h in raw_hits if h['qseqid']==p and h['sseqid']==architecture_review['signature']['reference_id'] and h['evalue']<=architecture_parameters['max_evalue']]
                    if not choices:bad.append(q+':architecture_missing_actual_alignment');continue
                    hit=max(choices,key=lambda h:h['bitscore']);positions=[];coordination=[]
                    if len(ar['site_evidence'])!=2 or {e['template_site_id'] for e in ar['site_evidence']}!=set(architecture_review['signature']['site_ids']):bad.append(q+':architecture_site_ownership')
                    for evidence in ar['site_evidence']:
                        spec=next(s for s in architecture_sites if s['site_id']==evidence['template_site_id'])
                        raw=site_evidence(hit,spec,architecture_parameters);coordination.append(raw['status'])
                        if any(evidence.get(k)!=raw[k] for k in ['status','reasons','cluster_type_predicted','alternative_coordination_states']):
                            bad.append(q+':architecture_coordination_evidence_changed')
                        if not (raw['status']=='template_supported' or (raw['status']=='ambiguous' and raw['reasons']==['possible_alternative_ligand'])):
                            bad.append(q+':unsupported_architecture_resource_state')
                        states,_=residue_states(hit,[int(p) for p in spec['expected_ligands']],spec['expected_ligands'])
                        lo,hi=spec['domain_start'],spec['domain_end'];seq=spec['reference_sequence']
                        windows,_=residue_states(hit,list(range(lo,hi+1)),{str(i):seq[i-1] for i in range(lo,hi+1)})
                        covered=[w for w in windows if w['query_position'] is not None];coverage=len(covered)/len(windows)
                        identity=100*sum(w['status']=='present' for w in covered)/max(1,len(covered))
                        positions.extend(w['query_position'] for w in states)
                        if (states!=evidence['ligand_residues'] or any(w['status']!='present' for w in states) or
                            coverage<architecture_parameters['min_window_coverage'] or identity<architecture_parameters['min_window_identity'] or
                            abs(coverage-evidence['window_coverage'])>1e-10 or abs(identity-evidence['window_identity'])>1e-10 or evidence['cofactor_occupancy']!='unknown'):
                            bad.append(q+':architecture_ligand_or_context_mapping')
                    scaffold=ar['scaffold_evidence'];spec=architecture_review['scaffold_site']
                    scaffold_states,_=residue_states(hit,[int(p) for p in spec['expected_ligands']],spec['expected_ligands'])
                    scaffold_present=all(w['status']=='present' for w in scaffold_states)
                    if scaffold_states!=scaffold['ligand_residues'] or scaffold['cofactor_occupancy']!='unknown' or scaffold_present!=ar['scaffold_ligands_present']:
                        bad.append(q+':scaffold_ligand_states')
                    missing=sorted(set(positions)-set(ar['selected_query_resources']))
                    if (ar.get('resource_mapping_scope')!='eight_positions_and_windows_only' or
                        ar.get('coordination_status')!=('ambiguous' if 'ambiguous' in coordination else 'template_supported')):
                        bad.append(q+':architecture_coordination_scope_changed')
                    if (len(set(positions))!=8 or sorted(positions)!=ar['signature_query_ligands'] or missing!=ar['unexplained_query_ligands'] or
                        ar['policy_sha256']!=selection['module_architecture']['policy_sha256'] or
                        (ar['status']=='selected_resources_explain_reviewed_signature')!=bool(not missing and scaffold_present)):
                        bad.append(q+':architecture_query_resource_ownership')
        check(label+':post_identity_L4_alternatives_and_guards',not bad,bad)
        bad=[];extension=selection['extensions']
        if (extension!=json.loads((run/'structure_extensions/search_manifest.json').read_text()) or
            extension['model_manifest_sha256']!=hashlib.sha256((ROOT/'references/processed/complex_extensions/manifest.json').read_bytes()).hexdigest() or
            extension['reaction_registry_sha256']!=hashlib.sha256((ROOT/'config/reaction_contracts.json').read_bytes()).hexdigest() or
            extension['script_sha256']!=hashlib.sha256((ROOT/'scripts/complex_extensions.py').read_bytes()).hexdigest()):bad.append('extension_source_changed')
        for name in ['structure_extensions/search_manifest.json','structure_extensions/queries.faa','structure_extensions/blast.tsv']:
            if manifest['output_sha256'].get(name)!=hashlib.sha256((run/name).read_bytes()).hexdigest():bad.append(name+':hash')
        for name,digest in selection['protected_before_and_after_sha256'].items():
            if hashlib.sha256((run/name).read_bytes()).hexdigest()!=digest:bad.append(name+':post_identity_mutation')
        if extension['query_fasta_sha256']!=hashlib.sha256((run/'structure_extensions/queries.faa').read_bytes()).hexdigest() or extension['raw_alignment_sha256']!=hashlib.sha256((run/'structure_extensions/blast.tsv').read_bytes()).hexdigest():bad.append('extension_query_or_alignment_changed')
        check(label+':extension_search_and_protected_identity_artifacts',not bad,bad)
        bad=[]
        for row in tables['reaction_contract_candidates.tsv']:
            c=contracts[row['contract_id']];a=anchor_index[row['complex_anchor']]
            assigned=c['contract_id'] in a['potential_function_ids'].split(';') or bool(c['legacy_function_id'] and c['legacy_function_id'] in a['potential_function_ids'].split(';'))
            expected=('existing_legacy_assignment_not_newly_calibrated' if c['legacy_function_id'] else 'exact_reference_sequence_potential') if assigned else 'unassigned_reference_reaction_candidate'
            if (row['structure_id']!=a['structure_id'] or c['structure_id']!=a['structure_id'] or row['status']!=expected or
                any(row[k]!=c[k] for k in ['reaction_key','hydrogen_role','ferredoxin_requirement','electron_coupling','capability_scope','direct_acceptor_status']) or
                any(json.loads(row[k])!=c[k] for k in ['electron_donors','electron_acceptors']) or
                row['query_net_direction']!='unresolved' or row['query_physiology']!='unresolved' or row['query_activity']!='not_measured'):
                bad.append(row['complex_anchor']+':'+row['contract_id'])
        cm=manifest['reaction_contract_candidates']
        if cm['rows']!=len(tables['reaction_contract_candidates.tsv']) or cm['output_sha256']!=manifest['output_sha256']['reaction_contract_candidates.tsv']:bad.append('contract_candidate_manifest')
        check(label+':typed_reference_candidates_do_not_assign_query_reactions',not bad,bad)
        bad=[];exact=manifest['exact_reference_functions'];ep=json.loads((ROOT/'config/exact_reference_functions.json').read_text())
        if (exact['policy_sha256']!=hashlib.sha256((ROOT/'config/exact_reference_functions.json').read_bytes()).hexdigest() or
            exact['script_sha256']!=hashlib.sha256((ROOT/'scripts/complex_extensions.py').read_bytes()).hexdigest() or
            exact['protected_before_and_after_sha256']!=selection['protected_before_and_after_sha256'] or exact['scope']!=ep['output_scope']):bad.append('exact_reference_provenance')
        sequences={r.id:str(r.seq) for r in SeqIO.parse(run/'proteins.faa','fasta')}
        with (run/'genes.tsv').open() as handle:gene_index={g['id']:g for g in csv.DictReader(handle,delimiter='\t')}
        assigned={q for q,a in anchor_index.items() if a['function_status']=='reference_sequence_supported_potential'}
        if assigned!=set(exact['assigned_anchors']) or assigned!={e['anchor'] for e in exact['evidence']}:bad.append('exact_reference_ownership')
        for e in exact['evidence']:
            a=anchor_index[e['anchor']];sid=a['structure_id'];q=a['protein_id'];ag=gene_index[q]
            rr={r['role']:r for r in extension_refs if r['structure_id']==sid}
            mm=[m for m in tables['complex_members.tsv'] if m['anchor']==q and m['required_for_structure']=='True']
            ids={c['contract_id'] for c in contracts.values() if c['structure_id']==sid and not c['legacy_function_id']}
            if (sid not in ep['eligible_structure_ids'] or manifest['input_mode'] not in ep['eligible_input_modes'] or
                a['structure_status'] not in ep['eligible_structure_statuses'] or a['catalytic_status']!=ep['catalytic_status'] or
                a['active_center_status']!=ep['active_center_status'] or a['active_center_span_status']!=ep['active_span_status'] or
                a['module_organization_status']!='reference_partition_compatible' or a['module_graph_search_status']!='complete' or
                a['cross_anchor_pairing_status'] in ep['forbidden_pairing_statuses'] or a['cross_anchor_function_ambiguity_ids'] or
                a['function_confidence']!='exact_reference_sequence_potential' or a['structure_confidence']=='high' or
                set(a['potential_function_ids'].split(';'))!=ids or set(e['contract_ids'])!=ids or e['structure_id']!=sid or
                len(mm)!=len(rr) or {m['role'] for m in mm}!=set(rr) or len({m['protein_id'] for m in mm})!=len(rr)):
                bad.append(q+':exact_reference_guard')
            expected_members=[];owned=[]
            for m in mm:
                r=rr[m['role']];g=gene_index[m['protein_id']];seq=sequences[m['protein_id']];owned.append(g)
                if (seq!=r['sequence'] or set(seq)&set('XBZJ') or g['partial']!='False' or
                    r['evidence_scope']!=ep['required_reference_scope'] or m['reference_protein']!=r['protein_accession'] or
                    m['structure_id']!=sid or g['contig']!=ag['contig'] or g['strand']!=ag['strand'] or
                    not valid_cds_interval(g, allow_unknown_length=False) or
                    (manifest['input_mode']=='bin' and g.get('bin_id')!=ag.get('bin_id'))):bad.append(q+':member_exactness_or_ownership')
                expected_members.append(dict(role=r['role'],protein_id=m['protein_id'],reference_id=r['reference_id'],native_reference_protein=r['protein_accession'],
                    query_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),reference_sequence_sha256=r['sequence_sha256']))
            if owned and not local_cluster([dict(g,start=int(g['start']),end=int(g['end'])) for g in owned],manifest['thresholds']):
                bad.append(q+':exact_reference_linkage')
            if sorted(e['matched_members'],key=lambda m:m['role'])!=sorted(expected_members,key=lambda m:m['role']):bad.append(q+':exact_reference_evidence')
        check(label+':exact_reference_reaction_potential_and_immutable_identity',not bad,bad)
        cp=manifest['catalytic_confirmation'];feature_manifest_path=run/'confirmation/feature_manifest.json';fm=json.loads(feature_manifest_path.read_text())
        check(label+':confirmation_provenance',cp['policy_sha256']==hashlib.sha256((ROOT/'config/catalytic_confirmation.json').read_bytes()).hexdigest() and
            cp['model_manifest_sha256']==hashlib.sha256((ROOT/'references/processed/catalytic_confirmation/manifest.json').read_bytes()).hexdigest() and
            cp['script_sha256']==hashlib.sha256((ROOT/'scripts/catalytic_confirmation.py').read_bytes()).hexdigest() and
            cp['thresholds_sha256']==hashlib.sha256((ROOT/'config/catalytic_confirmation_thresholds.json').read_bytes()).hexdigest() and
            cp['feature_manifest_sha256']==hashlib.sha256(feature_manifest_path.read_bytes()).hexdigest() and
            all(hashlib.sha256((run/'confirmation'/name).read_bytes()).hexdigest()==digest for name,digest in fm['artifacts'].items()),'Source/model/parameter and raw query evidence hashes')
        bad=[];policy=json.loads((ROOT/'config/discovery_evidence.json').read_text())
        for hit in tables['family_candidates.tsv']:
            parts=json.loads(hit['component_domains']);coverage=sum(d['model_end']-d['model_start']+1 for d in parts)/parts[0]['model_length']
            if abs(coverage-float(hit['profile_coverage']))>1e-10 or json.loads(hit['query_segments'])!=[[d['query_start'],d['query_end']] for d in parts]:bad.append(hit['protein_id']+':coordinates_or_coverage')
            if any(not chain_compatible(a,b,policy) for a,b in zip(parts,parts[1:])):bad.append(hit['protein_id']+':noncollinear')
            if float(hit['independent_evalue'])!=max(d['independent_evalue'] for d in parts):bad.append(hit['protein_id']+':component_evalue')
            if int(hit['evalue_reference_sequences'])!=policy['evalue_reference_sequences']:bad.append(hit['protein_id']+':reference_search_size')
            if len(parts)>1 and hit['evalue_scope']!='max_component_domain_E_not_combined_pvalue':bad.append(hit['protein_id']+':significance_scope')
        check(label+':discovery_segment_coordinates',not bad,bad)
        evidence=tables['annotation_evidence.tsv'];failures=[]
        for anchor in tables['catalytic_anchors.tsv']:
            qid=anchor['protein_id'];layers={e['layer'] for e in evidence if e['complex_anchor']==qid and e['status'] and e['source'] and e['limitation']}
            if layers!={'L1','L2','L3','L4','L5','L6'} or any(not anchor[k] for k in ['catalytic_status','structure_status','function_status','in_vivo_physiology']):failures.append(qid)
        check(label+':six_layers_and_statuses',not failures,failures)
        for row in evidence:json.loads(row['details'])
        proteins={r.id:str(r.seq) for r in SeqIO.parse(run/'proteins.faa','fasta')};bad=[]
        parameters=json.loads((ROOT/'config/catalytic_confirmation_thresholds.json').read_text())
        current_stage=confirmation_fingerprint(ROOT)
        stage_valid=parameters['feature_stage_fingerprint']==current_stage
        if not stage_valid and parameters.get('input_translation_transition'):
            from audit_input_translation_transfer import verify_registered_transition
            stage_valid=verify_registered_transition(ROOT,ROOT/parameters['input_translation_transition'],current_stage,
                old_stage_fingerprint=parameters['feature_stage_fingerprint'],rules=parameters)
        if not stage_valid and parameters.get('reporting_transition'):
            from audit_discovery_reporting_transfer import verify_reporting_transition
            stage_valid=verify_reporting_transition(ROOT,ROOT/parameters['reporting_transition'],current_stage,parameters)
        if not stage_valid and parameters.get('topology_transition'):
            from audit_topology_transition import verify_topology_transition
            stage_valid=verify_topology_transition(ROOT,ROOT/parameters['topology_transition'],current_stage,parameters)
        if not stage_valid:bad.append('stale_confirmation_stage')
        features={r['protein_id']:r for r in json.loads((run/'confirmation/features.json').read_text())}
        confirmation={r['protein_id']:r for r in tables['catalytic_confirmation.tsv']}
        if set(features)!=set(proteins) or set(confirmation)!=set(proteins):bad.append('incomplete_query_ownership')
        high_ids=set()
        for q,f in features.items():
            rule=parameters['families'].get(f['family']);high=bool(rule and qualifies(f,rule))
            if high:high_ids.add(q)
            if (confirmation[q]['identity_confidence']=='calibrated_high')!=high or any(not 1<=p<=len(proteins[q]) for p in f['core_positions']):bad.append(q+':label_or_core_coordinates')
        for anchor in tables['catalytic_anchors.tsv']:
            q=anchor['protein_id']
            if anchor['catalytic_family']!=features[q]['family']:bad.append(q+':identity_family_changed_after_confirmation')
            if (anchor['catalytic_identity_confidence']=='calibrated_high')!=(q in high_ids) or (q in high_ids and anchor['catalytic_status']!='supported'):bad.append(q+':high_identity_not_propagated')
        if cp['high_identity_outputs']!=len(high_ids):bad.append('high_identity_denominator')
        check(label+':confirmation_rule_coordinates_and_labels',not bad,bad)
        bad=[]
        node_map={(n['anchor'],n['node_id']):n for n in tables['module_graph_nodes.tsv']};used={}
        for n in node_map.values():
            spec=graph_spec[n['structure_id'],n['node_id']]
            if (n['reference_id']!=spec['reference_id'] or int(n['reference_start'])!=spec['reference_start'] or
                int(n['reference_end'])!=spec['reference_end'] or not n['localization_status'].startswith('not_predicted')):bad.append(n['node_id']+':reference_or_scope')
        for e in tables['module_graph_edges.tsv']:
            n=node_map[e['anchor'],e['node_id']];q=e['protein_id'];positions=set(json.loads(e['query_residue_positions']));rpositions=set(json.loads(e['reference_residue_positions']))
            key=e['anchor'],q
            if (not positions or min(positions)<1 or max(positions)>len(proteins[q]) or positions&used.get(key,set()) or
                min(rpositions)<int(n['reference_start']) or max(rpositions)>int(n['reference_end']) or
                positions!=set(json.loads(n['query_resources'])[q])):bad.append(q+':reused_or_invalid_resources')
            used.setdefault(key,set()).update(positions)
        with (run/'sequence_evidence.tsv').open() as handle:
            blast_identity={(h['qseqid'],h['sseqid']):float(h['pident']) for h in csv.DictReader(handle,delimiter='\t')}
        for m in tables['modules.tsv']:
            spec=graph_spec[m['structure_id'],m['module_node_id']]
            if float(m['identity'])!=blast_identity[m['protein_id'],spec['reference_id']]:bad.append(m['protein_id']+':legacy_BLAST_identity_changed')
        for a in tables['catalytic_anchors.tsv']:
            if (a['module_organization_status']=='changed_from_reference_unvalidated' or a['module_graph_search_status']=='truncated') and (a['potential_function_ids'] or a['structure_confidence']=='high'):
                bad.append(a['protein_id']+':unvalidated_organization_or_search_limit_transfer')
        check(label+':module_graph_resources_and_guards',not bad,bad)
        bad=[];anchors={a['protein_id']:a for a in tables['catalytic_anchors.tsv']}
        for pair in tables['cross_anchor_associations.tsv']:
            overlap=set(json.loads(pair['overlap_query_positions']));resources=[]
            for side in ['a','b']:
                n=node_map[pair['anchor_'+side],pair['node_'+side]];choices=json.loads(n['candidate_realizations'])
                selected=[choices[i] for i in json.loads(pair['candidate_realizations_'+side])]
                pp=[set(c['query_resources'][pair['protein_id']]) for c in selected]
                if not pp or any(p!=pp[0] for p in pp) or not all(c['jointly_feasible'] and c['competitive'] for c in selected):bad.append(pair['association_id']+':infeasible_claim')
                resources.append(pp[0])
                a=anchors[pair['anchor_'+side]]
                if pair['blocks_specific_pairing']=='True':
                    required=set(n['required_for_function'].split(';'))-{''}
                    if required&set(a['potential_function_ids'].split(';')) or (n['required_for_structure']=='True' and a['structure_confidence']=='high'):bad.append(pair['association_id']+':blocked_pairing_transferred')
            if not overlap or overlap!=resources[0]&resources[1] or len(overlap)!=int(pair['overlap_residue_count']):bad.append(pair['association_id']+':overlap_coordinates')
        check(label+':cross_anchor_claim_provenance_and_guards',not bad,bad)
        bad=[]
        for row in tables['fes_sites.tsv']:
            model=models[row['template_site_id']];sequence=proteins[row['protein_id']]
            if row['cofactor_occupancy']!='unknown' or row['source_sha256']!=model['source_sha256']:bad.append(row['fes_site_id']+':source_or_occupancy')
            recorded=json.loads(row['ligand_residues'])
            if row['alignment_source']=='reviewed_site_window_search':
                choices=[full_reference_coordinates(h,model) for h in actual_hits(run,'fes_site_window_blast.tsv',row['protein_id'],window_names[model['site_id']],row['query_alignment_start'],row['query_alignment_end'],row['bitscore'])]
            elif row['alignment_source']=='full_reference_search':
                choices=actual_hits(run,'fes_blast.tsv',row['protein_id'],model['reference_id'],row['query_alignment_start'],row['query_alignment_end'],row['bitscore'])
            else:choices=[]
            if not traced_fes_states_match(recorded,choices,model):bad.append(row['fes_site_id']+':residue_state_trace')
            for ligand in recorded:
                pos=ligand['query_position'];actual=sequence[pos-1] if pos is not None and 1<=pos<=len(sequence) else ''
                if ligand['expected_residue']!=model['expected_ligands'][str(ligand['reference_position'])]:bad.append(row['fes_site_id']+':reference_residue')
                if pos is not None and actual!=ligand['observed_residue']:bad.append(row['fes_site_id']+':query_coordinate')
        check(label+':FeS_coordinate_provenance',not bad,bad)
        active_models={m['reference_id']:m for name in ['active_centers','split_centers'] for m in json.loads((ROOT/'references/processed'/name/'models.json').read_text())};bad=[]
        for row in tables['catalytic_sites.tsv']:
            model=active_models[row['reference_id']];p=int(row['reference_position']);q=int(row['query_position']) if row['query_position'] else None
            actual=proteins[row['protein_id']][q-1] if q is not None else ''
            site=next(s for s in model['sites'] if s['reference_position']==p)
            traced=actual_active_states(run,row,model,proteins[row['protein_id']]).get(p,{})
            if (row['catalytic_anchor'] not in anchors or not row['source_center_id'] or not row['center_component_id'] or
                row['cofactor_occupancy']!='unknown' or row['source_sha256']!=model['source_sha256'] or row['expected_residue']!=site['expected_residue'] or
                row['status']!=traced.get('status') or q!=traced.get('query_position') or row['observed_residue']!=traced.get('observed_residue') or
                (q is not None and actual!=row['observed_residue'])):bad.append(row['protein_id']+':'+str(p))
        check(label+':complete_active_center_coordinates_and_source',not bad,bad)
        bad=[];components={}
        for row in tables['catalytic_center_components.tsv']:
            components.setdefault(row['center_id'],[]).append(row)
            model=active_models[row['reference_id']]
            if row['reference_accession']!=model['reference_accession'] or row['source_sha256']!=model['source_sha256']:bad.append(row['center_id']+':component_source')
            for h in json.loads(row['alignment_components']):
                if (not row['protein_id'] or not 1<=h['query_start']<=h['query_end']<=len(proteins[row['protein_id']]) or
                    not 1<=h['reference_start']<=h['reference_end']<=len(model['sequence'])):bad.append(row['center_id']+':component_coordinates')
        for row in tables['catalytic_center_candidates.tsv']:
            parts=components.get(row['center_id'],[]);sites=json.loads(row['sites'])
            if row['catalytic_anchor'] not in anchors or len(parts)!=2 or len(sites)!=4 or {p['component_id'] for p in parts}!={'VhuA','VhuU'}:bad.append(row['center_id']+':component_membership')
            for s in sites:
                model=active_models[s['reference_id']];p=s['reference_position'];q=s['query_position']
                actual=proteins[s['protein_id']][q-1] if q is not None else ''
                expected=next(t for t in model['sites'] if t['reference_position']==p)['expected_residue']
                state=actual_active_states(run,s,model,proteins[s['protein_id']]).get(p,{}).get('status') if s['protein_id'] else 'not_covered'
                if (s['catalytic_anchor']!=row['catalytic_anchor'] or s['expected_residue']!=expected or s['status']!=state or
                    s['cofactor_occupancy']!='unknown' or s['source_sha256']!=model['source_sha256'] or actual!=s['observed_residue']):bad.append(row['center_id']+':candidate_site')
        check(label+':split_center_components_sites_and_provenance',not bad,bad)
        summaries={r['protein_id']:r for r in tables['active_center_summary.tsv']};bad=[]
        for anchor in tables['catalytic_anchors.tsv']:
            qid=anchor['protein_id'];status=summaries.get(qid,{}).get('active_center_status','no_reviewed_template_mapping')
            if anchor['active_center_status']!=status:bad.append(qid+':status')
            if status not in ['template_ligand_set_supported','experimentally_active_exact_reference_noncanonical_site','multi_polypeptide_ligand_set_supported','fused_split_template_ligand_set_supported'] and (anchor['catalytic_status']=='supported' or anchor['potential_function_ids']):bad.append(qid+':unguarded_incomplete_center')
            if qid in summaries:
                summary=summaries[qid];contexts=json.loads(summary['site_contexts'])
                if (anchor['active_center_span_status']!=summary['template_span_alignment_status'] or
                    summary['catalytic_domain_integrity']!='not_established_by_active_site_mapping' or
                    float(summary['min_required_context_coverage'])!=min(c['coverage'] for c in contexts) or
                    float(summary['min_required_context_identity'])!=min(c['identity'] for c in contexts)):
                    bad.append(qid+':context_or_integrity_scope')
                if summary['template_span_alignment_status']!='supported' and (anchor['potential_function_ids'] or anchor['structure_confidence']=='high'):
                    bad.append(qid+':local_context_rescued_incomplete_scaffold')
                if status in ['multi_polypeptide_ligand_set_supported','fused_split_template_ligand_set_supported'] and (anchor['potential_function_ids'] or anchor['structure_confidence']=='high'):
                    bad.append(qid+':uncalibrated_split_complex_function')
        check(label+':active_center_confirmation_guard',not bad,bad)
    check('all_supported_input_modes',modes=={'genome','bin','assembly','transcript'},sorted(modes))
    check('empty_outputs_retain_full_schema',empty_cases>0,empty_cases)
    output=directory/'output_contract';output.mkdir(exist_ok=True)
    write_tsv(output/'checks.tsv',checks)
    summary=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',runs=len(runs),empty_candidate_runs=empty_cases,
        modes=sorted(modes),checks=len(checks),passed=sum(c['passed'] for c in checks),failed=[c for c in checks if not c['passed']],
        scope='Output schema, coordinate and traceability contract only; not independent biochemical accuracy')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if summary['failed']:raise SystemExit(1)


if __name__=='__main__':main()
