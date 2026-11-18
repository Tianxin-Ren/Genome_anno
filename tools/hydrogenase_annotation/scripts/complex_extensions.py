"""Source-pinned post-identity structure models and typed reference reaction contracts."""
import copy
import hashlib
import json
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, BLAST_FIELDS, load_blast,association_allowed,local_cluster
from circular_coordinates import valid_cds_interval
from catalytic_confirmation import sha
from module_graph import load_policy
from annotate import write_tsv

CONTRACT_FIELDS=('complex_anchor structure_id contract_id reaction_key electron_donors electron_acceptors hydrogen_role '
    'ferredoxin_requirement electron_coupling capability_scope direct_acceptor_status structure_status status '
    'query_net_direction query_physiology query_activity reference_source limitations').split()


def load_bank(root=ROOT):
    root=Path(root);directory=root/'references/processed/complex_extensions'
    manifest=json.loads((directory/'manifest.json').read_text())
    if manifest.get('status')!='reviewed_development_structure_extension_bank':
        raise ValueError('Structure extension bank is incomplete')
    refs=json.loads((directory/'references.json').read_text());systems=json.loads((directory/'systems.json').read_text());graphs=json.loads((directory/'graphs.json').read_text())
    contracts=load_contracts(root)
    return refs,systems,graphs,manifest,contracts


def load_contracts(root=ROOT):
    root=Path(root);registry=json.loads((root/'config/reaction_contracts.json').read_text())
    contracts=registry['contracts']
    if len({c['contract_id'] for c in contracts})!=len(contracts):raise ValueError('Repeated reaction contract')
    return contracts


def search_and_merge(genes,hits,references,systems,thresholds,mode,anchors,policy,out,threads=2,root=ROOT):
    """Search only a child output directory after identity; do not read or change identity evidence."""
    root=Path(root);out=Path(out)
    if out.name!='structure_extensions':raise ValueError('Expected dedicated post-identity artifact directory')
    out.mkdir(parents=True,exist_ok=False)
    extra_refs,extra_systems,extra_graphs,model,contracts=load_bank(root)
    if {r['reference_id'] for r in references}&{r['reference_id'] for r in extra_refs}:raise ValueError('Reference ID collision')
    if {s['structure_id'] for s in systems}&{s['structure_id'] for s in extra_systems}:raise ValueError('Structure ID collision')
    graph_policy=copy.deepcopy(load_policy(root,references));graph_policy['systems'].extend(extra_graphs)
    eligible=mode in policy['eligible_input_modes'] and any(a['catalytic_family'] in policy['eligible_families'] and
        a.get('active_center_status')==policy['required_active_center_status'] and
        a.get('active_center_span_status')==policy['required_active_span_status'] for a in anchors)
    queries=genes if eligible else []
    faa=out/'queries.faa';SeqIO.write([SeqRecord(Seq(g['sequence']),id=g['id'],description='') for g in queries],faa,'fasta')
    raw=out/'blast.tsv';raw.write_text('')
    if queries:
        with (out/'blast.log').open('w') as log:
            subprocess.run(['blastp','-query',str(faa),'-db',str(root/'references/processed/complex_extensions/seed'),'-evalue','1e-5',
                '-max_target_seqs','1000','-max_hsps','1','-num_threads',str(threads),'-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(raw)],stderr=log,check=True)
    extra_hits=load_blast(raw,thresholds,extra_refs)
    qids={g['id'] for g in queries};rids={r['reference_id'] for r in extra_refs}
    if any(h['qseqid'] not in qids or h['sseqid'] not in rids for h in extra_hits):raise ValueError('Extension hit ownership mismatch')
    manifest=dict(status='post_identity_structure_extension_search_complete',model_manifest_sha256=sha(root/'references/processed/complex_extensions/manifest.json'),
        reaction_registry_sha256=sha(root/'config/reaction_contracts.json'),script_sha256=sha(Path(__file__)),mode=mode,
        eligible_query_proteins=len(queries),source_proteins=len(genes),hits=len(extra_hits),
        reference_replay_cases=model['reference_replay_cases'],additional_structures=len(extra_systems),additional_proteins=len(extra_refs),
        query_fasta_sha256=sha(faa),raw_alignment_sha256=sha(raw),
        scope='L4 templates and reference reaction contracts only; new model references do not enter pre-identity discovery, active-center mapping or confirmation prototypes.')
    (out/'search_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return references+extra_refs,systems+extra_systems,hits+extra_hits,graph_policy,manifest,contracts


def write_contract_candidates(out,annotations,systems,contracts):
    """Expose typed reference hypotheses after pairing guards; never add a query assignment."""
    by_structure={s['structure_id']:s for s in systems};by_contract={c['contract_id']:c for c in contracts};rows=[]
    for a in annotations:
        system=by_structure.get(a['structure_id'])
        if not system:continue
        for cid in system.get('reference_contract_ids',system['function_ids']):
            c=by_contract[cid]
            if c['structure_id']!=a['structure_id']:raise ValueError('Contract/structure mismatch')
            assigned=bool(cid in a['potential_function_ids'].split(';') or c['legacy_function_id'] and c['legacy_function_id'] in a['potential_function_ids'].split(';'))
            rows.append(dict(complex_anchor=a['protein_id'],structure_id=a['structure_id'],contract_id=cid,
                **{k:json.dumps(c[k]) if isinstance(c[k],list) else c[k] for k in ['reaction_key','electron_donors','electron_acceptors','hydrogen_role',
                    'ferredoxin_requirement','electron_coupling','capability_scope','direct_acceptor_status']},
                structure_status=a['structure_status'],status=('exact_reference_sequence_potential' if not c['legacy_function_id'] else 'existing_legacy_assignment_not_newly_calibrated') if assigned else 'unassigned_reference_reaction_candidate',
                query_net_direction='unresolved',query_physiology='unresolved',query_activity='not_measured',reference_source=c['primary_source'],
                limitations=c['limitations']+' This candidate table does not assign a reaction; use reaction_capabilities.tsv for existing query assignments.'))
    path=Path(out)/'reaction_contract_candidates.tsv';write_tsv(path,rows,CONTRACT_FIELDS)
    return dict(registry_sha256=sha(ROOT/'config/reaction_contracts.json'),rows=len(rows),output_sha256=sha(path),
        scope='Typed reference hypotheses and assignment provenance. Exact-reference potential is separately guarded; no query activity, occupancy, physiology or net direction is measured.')


def exact_reference_patches(genes,references,systems,annotations,modules,contracts,policy,thresholds):
    """Replay only a complete experimentally characterized protein set; never infer a variant."""
    by_gene={g['id']:g for g in genes};by_system={s['structure_id']:s for s in systems};by_contract={c['contract_id']:c for c in contracts}
    patches={};evidence=[]
    for a in annotations:
        sid=a['structure_id'];qid=a['protein_id'];mode=a.get('input_mode')
        if (sid not in policy['eligible_structure_ids'] or mode not in policy['eligible_input_modes'] or
            a['structure_status'] not in policy['eligible_structure_statuses'] or a['potential_function_ids'] or a.get('catalytic_status')!=policy['catalytic_status'] or
            a.get('active_center_status')!=policy['active_center_status'] or a.get('active_center_span_status')!=policy['active_span_status'] or
            a.get('module_organization_status')!='reference_partition_compatible' or a.get('module_graph_search_status')!='complete' or
            a.get('cross_anchor_pairing_status') in policy['forbidden_pairing_statuses'] or a.get('cross_anchor_function_ambiguity_ids')):continue
        refs={r['role']:r for r in references if r['structure_id']==sid and r['kind'] in ['catalytic','module']}
        selected=[m for m in modules if m['anchor']==qid and (m['required_for_structure'] is True or m['required_for_structure']=='True')]
        if len(selected)!=len(refs) or {m['role'] for m in selected}!=set(refs):continue
        matched=[];used=set();valid=True
        for m in selected:
            r=refs[m['role']];g=by_gene[m['protein_id']];seq=g['sequence']
            if (m['structure_id']!=sid or m['protein_id'] in used or g['partial'] not in [False,'False'] or
                seq!=r['sequence'] or set(seq)&set('XBZJ') or r.get('evidence_scope')!=policy['required_reference_scope'] or
                m['reference_protein']!=r['protein_accession'] or g['contig']!=by_gene[qid]['contig'] or g['strand']!=by_gene[qid]['strand'] or
                not association_allowed(by_gene[qid],g,mode) or not all(isinstance(g[k],int) and not isinstance(g[k],bool) for k in ['start','end']) or not valid_cds_interval(g, allow_unknown_length=False)):
                valid=False;break
            used.add(m['protein_id']);matched.append(dict(role=r['role'],protein_id=g['id'],reference_id=r['reference_id'],
                native_reference_protein=r['protein_accession'],query_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),reference_sequence_sha256=r['sequence_sha256']))
        if not valid or not local_cluster([by_gene[q] for q in used],thresholds):continue
        ids=by_system[sid]['reference_contract_ids']
        if not ids or any(by_contract[c]['structure_id']!=sid or by_contract[c]['legacy_function_id'] for c in ids):
            raise ValueError('Exact-reference function source is not an enabled new contract')
        patches[qid]=dict(function_status='reference_sequence_supported_potential',function_confidence='exact_reference_sequence_potential',
            potential_function_ids=';'.join(ids),potential_direct_partner=by_system[sid]['direct_partner'],
            flags=';'.join(filter(None,[a.get('flags',''),'exact_reference_complex_sequence_potential_not_measured_query_activity'])))
        evidence.append(dict(anchor=qid,structure_id=sid,contract_ids=ids,matched_members=matched,
            primary_sources=[by_contract[c]['primary_source'] for c in ids],scope=policy['output_scope']))
    return patches,evidence


def function_reviews(contracts,root=ROOT):
    """Join legacy reviews with the separate exact-reference capability records."""
    result={r['function_id']:r for r in json.loads((Path(root)/'config/reference_function_review.json').read_text())['records']}
    for c in contracts:
        if c['legacy_function_id']:continue
        result[c['contract_id']]=dict(function_id=c['contract_id'],structure_id=c['structure_id'],molecular_reaction=c['molecular_reaction'],
            reference_evidence_scope=c['reference_evidence_scope'],capability_scope='exact_reference_sequence_molecular_reaction_potential',reference_context='',
            limitations=c['limitations']+' Exact full core sequence replay supports reference reaction potential; query activity, occupancy, kinetics and in-vivo net flux remain unmeasured.')
    return result


def exact_reaction_rows(annotations,contracts):
    """Serialize new exact-reference assignments without altering legacy function records."""
    by_contract={c['contract_id']:c for c in contracts};reviews=function_reviews(contracts);rows=[]
    for a in annotations:
        for cid in filter(None,a['potential_function_ids'].split(';')):
            c=by_contract.get(cid)
            if c is None:raise ValueError('Missing typed reaction contract: '+cid)
            if c['legacy_function_id']:continue
            if a['function_status']!='reference_sequence_supported_potential' or c['structure_id']!=a['structure_id']:
                raise ValueError('Invalid exact-reference assignment')
            rows.append(dict(complex_anchor=a['protein_id'],function_id=cid,structure_id=a['structure_id'],status='reference_sequence_supported_potential',
                reference_reaction=c['molecular_reaction'],reference_direction=c['hydrogen_role'],query_net_direction='unresolved',source=c['primary_source'],
                **{k:reviews[cid][k] for k in ['molecular_reaction','reference_evidence_scope','capability_scope','reference_context','limitations']},
                query_physiology='unresolved',association_pairing_status=a['cross_anchor_pairing_status'],physical_interaction_status=a['physical_interaction_status']))
    return rows
