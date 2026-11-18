#!/usr/bin/env python3
"""Audit accession-pinned reference roles without changing sequence or functional rules."""
import argparse
import copy
import hashlib
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, BLAST_FIELDS, load_blast, write_tsv
from audit_maturation_domains import scan_profile_domains, architecture_supported
from discovery_evidence import best_chain


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_sequence(row):
    if row.get('entryType')=='Inactive' or 'sequence' not in row:return None
    return row['sequence']['value']


def source_id(row):
    sequence=source_sequence(row)
    if sequence is None:return None
    return 'UP_'+row['primaryAccession']+'__'+hashlib.sha256(sequence.encode()).hexdigest()[:12]


def model_evidence(domains, policy):
    groups=defaultdict(list)
    for d in domains:
        if d['sequence_bitscore']>=d['sequence_ga'] and d['domain_bitscore']>=d['domain_ga']:
            groups[d['accession'].split('.')[0]].append(d)
    evidence=[]
    for accession,parts in sorted(groups.items()):
        chain=best_chain([dict(d,protein_id='reference',domain_bitscore=d['domain_bitscore']) for d in parts],policy['chain_policy'])
        coverage=sum(d['model_end']-d['model_start']+1 for d in chain)/chain[0]['model_length']
        evidence.append(dict(accession=accession,profile=chain[0]['profile'],model_coverage=coverage,
            supported=coverage>=policy['min_model_coverage'],evidence_kind='single_domain' if len(chain)==1 else 'non_overlapping_collinear_segments',
            bitscore_scope='component_scores_not_a_combined_significance',components=chain))
    return evidence


def classify_role(kind, rule, evidence, structural=(), curated=(), bounded_primary=False, curated_homology=False):
    accepted={d['accession'] for d in evidence if d['supported']}
    partial={d['accession'] for d in evidence if not d['supported']}
    expected=set(rule.get('all_of',[]))|set(rule.get('any_of',[]))
    architecture=bool(expected) and architecture_supported(accepted,rule)
    full_structure=any(s['coverage_scope']=='near_complete_sequence' for s in structural)
    reviewed_sequence=any('Swiss-Prot' in s['entry_type'] for s in curated)
    if architecture:return 'reviewed_profile_architecture','supported_reference_family'
    if full_structure:
        exact_or_close=any(s['coverage_scope']=='near_complete_sequence' and s['identity']>=90 for s in structural)
        return ('reviewed_experimental_structure_correspondence' if exact_or_close else 'reviewed_bounded_structural_homology'),'supported_reference_scaffold'
    if reviewed_sequence:return 'reviewed_curated_sequence_partial_or_no_domain','supported_reference_sequence_with_domain_limits'
    if bounded_primary and bool((accepted|partial)&expected):
        return 'reviewed_primary_role_with_partial_domain_evidence','supported_reference_role_with_explicit_limits'
    if curated_homology:return 'reviewed_curated_homology_and_operon_inference','supported_reference_role_with_explicit_limits'
    if kind=='unresolved':return 'reviewed_unresolved_reference_not_function_enabled','unresolved_reference_retained'
    return 'pending_reference_role_evidence','requires_more_evidence'


def check_reference(row, review, sequence):
    if (row['reference_id']!=review['reference_id'] or row['protein_accession']!=review['protein_accession'] or
        row['sequence_sha256']!=review['sequence_sha256'] or row['kind']!=review['original_kind'] or
        hashlib.sha256(sequence.encode()).hexdigest()!=row['sequence_sha256']):
        raise ValueError('Reference accession, kind or sequence changed: '+row['reference_id'])


def blast(query, database_fasta, output, name):
    with (output/(name+'.log')).open('w') as log:
        subprocess.run(['makeblastdb','-in',str(database_fasta),'-dbtype','prot','-out',str(output/name)],stdout=log,stderr=log,check=True)
        subprocess.run(['blastp','-query',str(query),'-db',str(output/name),'-evalue','1e-10','-max_target_seqs','1000','-max_hsps','1',
            '-outfmt','6 '+' '.join(BLAST_FIELDS),'-out',str(output/(name+'.tsv'))],stderr=log,check=True)
    thresholds=dict(blast_evalue=1e-10,min_identity=0,min_query_coverage=0,min_reference_coverage=0)
    return load_blast(output/(name+'.tsv'),thresholds)


def audit(records, root, output):
    root=Path(root);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    policy_path=root/'config/reference_domain_policy.json';role_path=root/'config/reference_domain_review.json';structure_path=root/'config/reference_structure_review.json'
    function_path=root/'config/reference_function_review.json'
    policy=json.loads(policy_path.read_text());review=json.loads(role_path.read_text());structures=json.loads(structure_path.read_text())
    roles={r['reference_id']:r for r in review['records']}
    if len(roles)!=len(review['records']) or set(roles)!={r['reference_id'] for r in records if r['kind']!='maturation'}:
        raise ValueError('Reference review is duplicated, missing or contains unused records')
    functions=json.loads(function_path.read_text());function_index={r['function_id']:r for r in functions['records']}
    if len(function_index)!=len(functions['records']):raise ValueError('Duplicate function review')
    reviewed_functions=[];systems=json.loads((root/'config/systems.json').read_text())
    for system in systems:
        for function in system['reference_function_records']:
            row=function_index.pop(function['功能ID'])
            if (row['structure_id']!=system['structure_id'] or row['primary_source']!=system['evidence_url'] or
                row['reference_record_sha256']!=hashlib.sha256(json.dumps(function,sort_keys=True,ensure_ascii=False).encode()).hexdigest() or
                not all(row[k] for k in ['reference_evidence_scope','capability_scope','molecular_reaction','limitations'])):
                raise ValueError('Reference function source/scope changed')
            reviewed_functions.append(row)
    if function_index:raise ValueError('Unused function review')
    for relative,digest in {**policy['source_files'],**policy['curated_source_files']}.items():
        if sha(root/relative)!=digest:raise ValueError('Reference review source changed: '+relative)
    for source in policy['structural_sources'].values():
        if sha(root/source['path'])!=source['sha256']:raise ValueError('Structural source changed')
    for key in ['structural_fasta','curated_fasta','curated_source_audit']:
        if sha(root/policy[key])!=policy[key+'_sha256']:raise ValueError('Prepared reference source changed: '+key)
    sequences={r.id:str(r.seq) for r in SeqIO.parse(root/'references/processed/references.faa','fasta')}
    for row in records:
        if row['kind']!='maturation':check_reference(row,roles[row['reference_id']],sequences[row['reference_id']])
    query=output/'references_U_to_X.faa'
    SeqIO.write([SeqRecord(Seq(sequences[r['reference_id']].replace('U','X')),id=r['reference_id'],description='HMM_only_Sec_normalization') for r in records if r['kind']!='maturation'],query,'fasta')
    original_query=output/'references_original.faa'
    SeqIO.write([SeqRecord(Seq(sequences[r['reference_id']]),id=r['reference_id'],description='original_reference_U_preserved') for r in records if r['kind']!='maturation'],original_query,'fasta')
    models={};texts=[]
    for model in policy['models']:
        path=root/'references/raw/reference_domain_review'/(model['accession']+'.hmm')
        if sha(path)!=model['sha256']:raise ValueError('Curated profile changed')
        fields={line.split()[0]:' '.join(line.split()[1:]) for line in path.read_text().split('HMM          ')[0].splitlines()[1:] if line.strip()}
        ga=fields['GA'].strip(';').split()
        if fields['NAME']!=model['profile'] or fields['ACC']!=model['versioned_accession'] or len(ga)!=2:
            raise ValueError('Profile identity or gathering threshold changed')
        if model['profile'] in models:raise ValueError('Duplicate profile name')
        models[model['profile']]=dict(model,accession=model['versioned_accession'],sequence_ga=float(ga[0]),domain_ga=float(ga[1]));texts.append(path.read_text())
    combined=output/'reference_domains.hmm';combined.write_text(''.join(texts))
    domains=scan_profile_domains(combined,query,models,policy['min_model_coverage'],output)
    structural_hits=blast(original_query,root/policy['structural_fasta'],output,'structures')
    curated_hits=blast(original_query,root/policy['curated_fasta'],output,'curated')
    structural_index=defaultdict(list)
    structural_sequences={r.id:str(r.seq) for r in SeqIO.parse(root/policy['structural_fasta'],'fasta')}
    for row in structures['records']:structural_index[row['reference_id']].append(row)
    source_records={};inactive=[]
    for relative in policy['curated_source_files']:
        payload=json.loads((root/relative).read_text())
        for row in payload.get('results',[payload]):
            key=source_id(row)
            if key is None:
                inactive.append(dict(accession=row.get('primaryAccession'),entry_type=row.get('entryType'),source=relative));continue
            source_records[key]=dict(uniprot=row['primaryAccession'],entry_type=row['entryType'],source=relative,
                source_sha256=sha(root/relative),sequence_sha256=hashlib.sha256(source_sequence(row).encode()).hexdigest(),
                protein_name=row.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value',''))
    results=[];updated=copy.deepcopy(records)
    for row in updated:
        qid=row['reference_id']
        if row['kind']=='maturation':continue
        spec=roles[qid];rule=review['classes'][spec['role_class']];fold=model_evidence(domains[qid],policy)
        structure=[]
        for item in structural_index[qid]:
            if item['protein_accession']!=row['protein_accession'] or item['sequence_sha256']!=row['sequence_sha256']:
                raise ValueError('Reviewed structural correspondence belongs to another reference')
            pdb=item['structure_entity'].split('__entity_')[0]
            if (item['source_sha256']!=policy['structural_sources'][pdb]['sha256'] or
                hashlib.sha256(structural_sequences[item['structure_entity']].encode()).hexdigest()!=item['structure_sequence_sha256']):
                raise ValueError('Reviewed structural source/sequence changed')
            matching=[h for h in structural_hits if h['qseqid']==qid and h['sseqid']==item['structure_entity'] and
                abs(h['pident']-item['identity'])<.001 and (h['qstart'],h['qend'],h['sstart'],h['send'])==
                (item['query_start'],item['query_end'],item['structure_start'],item['structure_end'])]
            if len(matching)!=1:raise ValueError('Reviewed structural alignment changed: '+qid)
            structure.append(item)
        curated=[]
        for h in curated_hits:
            if h['qseqid']==qid and h['pident']==100 and min(h['qcov'],h['scov'])>=.9:
                curated.append(dict(**source_records[h['sseqid']],identity=h['pident'],query_coverage=h['qcov'],source_coverage=h['scov']))
        homolog=spec.get('curated_homology_review');homology=None
        if homolog:
            matching=[h for h in curated_hits if h['qseqid']==qid and h['sseqid']==homolog['source_id'] and
                abs(h['pident']-homolog['identity'])<.001 and (h['qstart'],h['qend'],h['sstart'],h['send'])==
                (homolog['query_start'],homolog['query_end'],homolog['source_start'],homolog['source_end'])]
            if len(matching)!=1 or 'Swiss-Prot' not in source_records[homolog['source_id']]['entry_type']:
                raise ValueError('Reviewed curated homolog changed')
            homology=dict(homolog,source_record=source_records[homolog['source_id']])
        primary=spec.get('primary_role_review')
        if primary and primary['gene_alias'] not in row['aliases']:raise ValueError('Primary gene/CDS mapping changed')
        if primary and any(sequences[qid][int(p)-1]!=residue for p,residue in primary.get('required_positions',{}).items()):
            raise ValueError('Reviewed primary experimental residue changed')
        state,scope=classify_role(row['kind'],rule,fold,structure,curated,bool(primary),bool(homology))
        evidence=dict(role_class=spec['role_class'],role_scope=rule['scope'],profile_evidence=fold,
            structural_correspondences=structure,curated_sequence_correspondences=curated,curated_homology_review=homology,
            primary_role_review=primary,query_function_inference='not_performed_by_reference_audit',cofactor_occupancy='unknown')
        row.update(domain_validation=state,domain_specificity_scope=rule['scope'],reference_role_review_status=scope,
            reference_role_evidence=evidence,reference_role_review_source=spec['primary_source'])
        results.append(dict(reference_id=qid,protein_accession=row['protein_accession'],sequence_sha256=row['sequence_sha256'],kind=row['kind'],
            role_class=spec['role_class'],domain_validation=state,role_review_status=scope,
            domain_evidence=evidence,source=spec['primary_source']))
    if set(structural_index)-set(roles):raise ValueError('Structural review has unused proteins')
    unresolved=[r['reference_id'] for r in results if r['role_review_status']=='requires_more_evidence']
    encoded=(json.dumps(updated,indent=2)+'\n').encode()
    summary=dict(status='review_complete' if not unresolved else 'review_complete_with_open_items',references=len(records),non_maturation_reviewed=len(results),
        profile_models=len(models),rows=results,unresolved_reference_roles=unresolved,role_status_counts=dict(Counter(r['role_review_status'] for r in results)),
        domain_status_counts=dict(Counter(r['domain_validation'] for r in results)),inactive_source_records_not_used=inactive,
        source_sequence_and_reference_kinds_unchanged=all((a['sequence_sha256'],a['kind'])==(b['sequence_sha256'],b['kind']) for a,b in zip(records,updated)),
        reviewed_references_sha256=hashlib.sha256(encoded).hexdigest(),reference_fasta_sha256=sha(root/'references/processed/references.faa'),
        structure_systems=len(systems),function_records_reviewed=len(reviewed_functions),function_source_scopes=reviewed_functions,
        config_sha256={str(p.relative_to(root)):sha(p) for p in [policy_path,role_path,structure_path,function_path,root/'config/systems.json']},
        scope='Reference family, scaffold and role curation with explicit evidence limits. Not independent classifier accuracy, full-domain validation or query biochemical truth.')
    write_tsv(output/'audit.tsv',[{**r,'domain_evidence':json.dumps(r['domain_evidence'],sort_keys=True)} for r in results])
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return updated,summary


def write_reference_tables(records, output):
    write_tsv(output/'reference_proteins.tsv',[{**r,**{k:';'.join(map(str,r[k])) for k in
        ['diagnostic_positions','selenocysteine_positions','cysteine_positions']}} for r in records])
    write_tsv(output/'reference_evidence.tsv',[dict(reference_id=r['reference_id'],kind=r['kind'],
        accession=r['protein_accession'],sequence_sha256=r['sequence_sha256'],source_accession=r['source_accession'],
        source_url=r['source_url'],evidence_url=r['evidence_url'],evidence_scope=r['evidence_scope'],source_product=r['source_product'],
        domain_validation=r['domain_validation'],reference_role_review_status=r.get('reference_role_review_status','maturation_profile_review'),
        reference_role_review_source=r.get('reference_role_review_source',r['evidence_url'])) for r in records])


def verify_applied_review(root, manifest):
    path=root/'references/processed/reference_domain_audit/summary.json'
    if not path.exists() or manifest.get('reference_domain_audit_sha256')!=sha(path):
        raise ValueError('Reference role audit is missing or changed; apply a fresh reference audit')
    summary=json.loads(path.read_text())
    if (summary['status']!='review_complete' or not summary['source_sequence_and_reference_kinds_unchanged'] or
        summary['reviewed_references_sha256']!=sha(root/'references/processed/references.json') or
        summary['reference_fasta_sha256']!=sha(root/'references/processed/references.faa') or
        any(sha(root/p)!=digest for p,digest in summary['config_sha256'].items())):
        raise ValueError('Reference role audit no longer matches references/configuration')
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--apply',action='store_true')
    args=parser.parse_args()
    if (args.out/'summary.json').exists():parser.error('Use a fresh audit directory')
    references_path=ROOT/'references/processed/references.json';records=json.loads(references_path.read_text())
    updated,summary=audit(records,ROOT,args.out)
    if args.apply:
        if args.out.resolve()!=ROOT/'references/processed/reference_domain_audit':
            raise ValueError('Applied review must use the canonical processed reference_domain_audit directory')
        if summary['unresolved_reference_roles']:raise ValueError('Do not apply incomplete reference-role review')
        references_path.write_text(json.dumps(updated,indent=2)+'\n')
        write_reference_tables(updated,references_path.parent)
        manifest_path=references_path.parent/'build_manifest.json';manifest=json.loads(manifest_path.read_text())
        manifest['reference_domain_audit_sha256']=sha(args.out/'summary.json')
        manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['rows','config_sha256','function_source_scopes']},indent=2))


if __name__=='__main__':main()
