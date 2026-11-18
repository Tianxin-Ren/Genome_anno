"""Explain selected module resources with existing reviewed FeS signatures, after identity."""
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT
from catalytic_confirmation import sha
from fes_sites import site_evidence


def load_review(root=ROOT):
    root=Path(root);path=root/'config/module_architecture.json';policy=json.loads(path.read_text())
    spec=policy['signature']
    with (root/'references/processed/references.faa').open() as handle:sequences={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    seq=sequences[spec['reference_id']]
    index={s['site_id']:s for s in json.loads((root/'references/processed/fes/sites.json').read_text())}
    sites=[index[s] for s in spec['site_ids']]
    policy['scaffold_site']=index[spec['scaffold_site_id']]
    identity={r['reference_id']:r for r in json.loads((root/'config/reference_catalytic_identity_review.json').read_text())['records']}
    refs={r['reference_id']:r for r in json.loads((root/'references/processed/references.json').read_text())}
    parameters=json.loads((root/'config/fes_site_review.json').read_text())['parameters']
    policy['identity_positive_competitors']={rid:i for rid,i in identity.items() if i['identity_label']=='positive' and refs[rid]['kind']=='unresolved' and refs[rid]['sequence_sha256']==i['sequence_sha256']}
    return policy,sites,parameters


def signature_evidence(hit,sites,parameters):
    if hit['slen']!=len(sites[0]['reference_sequence']):raise ValueError('Signature hit has incorrect reference coordinates')
    rows=[site_evidence(hit,s,parameters) for s in sites]
    return rows if all(r['status']=='template_supported' for r in rows) and len({p for r in rows for p in r['query_ligand_positions']})==8 else []


def module_resource_evidence(hit,sites,parameters):
    if hit['slen']!=len(sites[0]['reference_sequence']):raise ValueError('Resource hit has incorrect reference coordinates')
    rows=[site_evidence(hit,s,parameters) for s in sites]
    # A potential donor gain changes ligation certainty, not the mapped protein resources.
    supported=lambda r:r['status']=='template_supported' or (r['status']=='ambiguous' and r['reasons']==['possible_alternative_ligand'])
    return rows if all(supported(r) and r['domain_alignment']=='supported' and
        all(s['status']=='present' for s in r['ligand_residues']) for r in rows) and len({p for r in rows for p in r['query_ligand_positions']})==8 else []


def explain(graph,sid,hits,genes,review):
    policy,sites,parameters=review;spec=policy['signature'];role=spec['comparison_roles'].get(sid)
    if role is None:return 'not_assessed_outside_reviewed_comparison',[]
    nodes=[n for n in graph['nodes'] if n['role']==role and n['node_id'] in graph['selected']]
    evidence=[]
    for node in nodes:
        parts=graph['selected'][node['node_id']]['parts']
        if len(parts)!=1:continue
        part=parts[0];q=part['protein_id']
        if genes[q]['partial']:continue
        choices=[h for h in hits if h['qseqid']==q and h['sseqid']==spec['reference_id'] and h['evalue']<=parameters['max_evalue']]
        if not choices:continue
        hit=max(choices,key=lambda h:h['bitscore']);mapped=module_resource_evidence(hit,sites,parameters)
        if not mapped:continue
        scaffold=site_evidence(hit,policy['scaffold_site'],parameters)
        scaffold_present=all(s['status']=='present' for s in scaffold['ligand_residues'])
        positions={p for r in mapped for p in r['query_ligand_positions']};unexplained=sorted(positions-part['query_positions'])
        status=('unexplained_reviewed_FeS_extension' if unexplained else
            'reviewed_signature_supported_scaffold_ligands_incomplete' if not scaffold_present else 'selected_resources_explain_reviewed_signature')
        evidence.append(dict(signature_id=spec['id'],node_id=node['node_id'],role=role,protein_id=q,
            selected_query_resources=sorted(part['query_positions']),signature_query_ligands=sorted(positions),unexplained_query_ligands=unexplained,
            status=status,scaffold_evidence=scaffold,scaffold_ligands_present=scaffold_present,
            resource_mapping_scope='eight_positions_and_windows_only',
            coordination_status='ambiguous' if any(r['status']=='ambiguous' for r in mapped) else 'template_supported',
            site_evidence=mapped,policy_sha256=sha(ROOT/'config/module_architecture.json'),source=spec['source'],scope=spec['site_scope']))
    status=('unexplained_reviewed_FeS_extension' if any(e['unexplained_query_ligands'] for e in evidence) else
        'reviewed_signature_supported_scaffold_ligands_incomplete' if any(not e['scaffold_ligands_present'] for e in evidence) else 'selected_resources_explain_reviewed_signature' if evidence else 'signature_not_supported_not_evidence_of_absence')
    return status,evidence


def competitor_guard(anchor,sid,score,competitors,refs,architecture_status,review):
    stronger=[h for h in competitors if h['bitscore']>=score]
    if not stronger:return 'strictly_outscores_unresolved',[]
    policy,_,_=review;allowed=policy['identity_positive_competitors'] if sid in policy['reconciled_structure_ids'] else {}
    compatible=[]
    for h in stronger:
        r=refs[h['sseqid']];c=allowed.get(h['sseqid'])
        if not c or r['kind']!='unresolved' or r['sequence_sha256']!=c['sequence_sha256']:return 'blocked_by_unreviewed_or_incompatible_competitor',[]
        compatible.append(dict(reference_id=r['reference_id'],sequence_sha256=r['sequence_sha256'],bitscore=h['bitscore'],identity_label=c['identity_label'],source=str(ROOT/'config/reference_catalytic_identity_review.json'),scope=policy['competition_scope']))
    if anchor.get('catalytic_identity_confidence')!='calibrated_high' or architecture_status!='selected_resources_explain_reviewed_signature':
        return 'blocked_without_independent_identity_and_explained_signature',compatible
    return 'independent_identity_and_identity_positive_competitors',compatible
