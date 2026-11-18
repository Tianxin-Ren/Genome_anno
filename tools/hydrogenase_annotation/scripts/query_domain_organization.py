#!/usr/bin/env python3
"""Scan completed runs for auxiliary query-domain organization, without label transfer."""
import argparse
import hashlib
import json
import math
import shutil
import subprocess
import time
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from annotate import ROOT, write_tsv
from audit_maturation_domains import scan_profile_domains
from discovery_evidence import load_targets
from audit_reference_domains import model_evidence
from localization import sha, table, verify_source_file

OWNER_FIELDS = 'source_run protein_id query_id bin_id contig partial catalytic_family'.split()
DOMAIN_FIELDS = 'accession profile domain_bitscore sequence_bitscore independent_evalue model_coverage model_start model_end model_length query_length query_start query_end passes sequence_ga domain_ga source model_sha256'.split()
TARGET_FIELDS = 'query_id accession profile sequence_evalue sequence_bitscore sequence_bias model_length query_length model_sha256 source'.split()
SCOPE = ('Auxiliary fold organization only. Overlapping profiles are alternative views, not extra subunits. '
         'PF01058 is a shared oxidoreductase fold, not specific HyhS identity. No catalytic identity, '
         'FeS type/occupancy, physical complex, electron partner, mechanism or function is assigned.')


def load_models(root):
    root = Path(root)
    policy_path = root/'config/reference_domain_policy.json'
    policy = json.loads(policy_path.read_text())
    combined = root/'references/processed/reference_domain_audit/reference_domains.hmm'
    headers = {}
    for block in combined.read_text().split('HMMER3/f'):
        if not block.strip():
            continue
        fields = {line.split()[0]: ' '.join(line.split()[1:]) for line in block.split('HMM          ')[0].splitlines()[1:] if line.strip()}
        if 'NAME' in fields:
            headers[fields['NAME']] = fields
    models = {}
    for model in policy['models']:
        fields = headers[model['profile']]
        ga = fields['GA'].strip(';').split()
        if fields['ACC'] != model['versioned_accession'] or len(ga) != 2 or model['profile'] in models:
            raise ValueError('Domain profile identity or gathering threshold invalid')
        models[model['profile']] = dict(model, accession=model['versioned_accession'], sequence_ga=float(ga[0]), domain_ga=float(ga[1]), model_length=int(fields['LENG']))
    resources = {str(policy_path.resolve()): sha(policy_path), str(combined.resolve()): sha(combined)}
    return policy, models, combined, resources


def collect_sources(runs):
    owners = []; queries = {}; inputs = {}
    for run in runs:
        run = Path(run).resolve(); path = run/'run_manifest.json'
        manifest = json.loads(path.read_text())
        if manifest['status'] != 'complete':
            raise ValueError('Source annotation is incomplete')
        inputs[str(path)] = sha(path)
        for name in ['proteins.faa', *manifest['output_sha256']]:
            verify_source_file(run, manifest, name)
            inputs[str(run/name)] = sha(run/name)
        with (run/'proteins.faa').open() as handle: records = list(SeqIO.parse(handle, 'fasta'))
        genes = table(run/'genes.tsv'); anchors = table(run/'catalytic_anchors.tsv')
        by_id = {r.id: str(r.seq).upper() for r in records}; by_gene = {g['id']: g for g in genes}
        by_anchor = {a['protein_id']: a for a in anchors}
        if len(by_id) != len(records) or len(by_gene) != len(genes) or len(by_anchor) != len(anchors) or set(by_id) != set(by_gene) or set(by_anchor)-set(by_id):
            raise ValueError('Repeated or missing source protein/coordinate/anchor ownership')
        for pid, sequence in sorted(by_id.items()):
            if not sequence or set(sequence)-set('ACDEFGHIKLMNPQRSTVWYBJZOUX'):
                raise ValueError('Invalid normalized source protein sequence')
            qid = hashlib.sha256(sequence.encode()).hexdigest(); queries[qid] = sequence
            gene = by_gene[pid]; anchor = by_anchor.get(pid, {})
            if anchor and (anchor.get('bin_id', '') != gene.get('bin_id', '') or anchor['contig'] != gene['contig']):
                raise ValueError('Source anchor ownership differs from gene')
            owners.append(dict(source_run=str(run), protein_id=pid, query_id=qid, bin_id=gene.get('bin_id', ''),
                contig=gene['contig'], partial=gene['partial'], catalytic_family=anchor.get('catalytic_family', 'not_a_catalytic_anchor')))
    return owners, queries, inputs


def validate_domains(domains, queries, models, min_coverage):
    for qid, hits in domains.items():
        if qid not in queries:
            raise ValueError('Unknown domain query ownership')
        for d in hits:
            model = models.get(d['profile'])
            if (not model or d['accession'] != model['accession'] or d['query_length'] != len(queries[qid]) or d['model_length'] != model['model_length'] or
                not 1 <= d['model_start'] <= d['model_end'] <= d['model_length'] or
                not 1 <= d['query_start'] <= d['query_end'] <= d['query_length'] or
                any(not math.isfinite(d[k]) for k in ['independent_evalue', 'domain_bitscore', 'sequence_bitscore']) or d['independent_evalue'] < 0):
                raise ValueError('Invalid query domain coordinates, scores or profile')
            coverage = (d['model_end']-d['model_start']+1)/d['model_length']
            passed = d['sequence_bitscore'] >= model['sequence_ga'] and d['domain_bitscore'] >= model['domain_ga'] and coverage >= min_coverage
            if (d['model_coverage'] != coverage or d['passes'] != passed or
                any(d[k] != model[k] for k in ['sequence_ga', 'domain_ga']) or d['model_sha256'] != model['sha256'] or d['source'] != model['url']):
                raise ValueError('Domain coverage, thresholds or source differs from curated model')


def organization(hits, policy):
    evidence = model_evidence(hits, policy)
    # Keep full single domains as well as the best split chain; repeats are alternatives, never summed.
    candidates = {}
    for e in evidence:
        if e['supported']: candidates.setdefault(e['accession'], []).append(e['components'])
    for d in hits:
        if d['passes']: candidates.setdefault(d['accession'].split('.')[0], []).append([d])
    pairs = {}
    for left in candidates.get('PF01058', []):
        for right in candidates.get('PF02906', []):
            if max(d['query_end'] for d in left) >= min(d['query_start'] for d in right): continue
            pair = dict(PF01058_segments=[[d['query_start'], d['query_end']] for d in left],
                        PF02906_segments=[[d['query_start'], d['query_end']] for d in right])
            pairs[json.dumps(pair, sort_keys=True)] = pair
    overlaps = []
    for i, a in enumerate(hits):
        for j, b in enumerate(hits[i+1:], i+1):
            if max(a['query_start'], b['query_start']) <= min(a['query_end'], b['query_end']):
                overlaps.append([i+1, j+1])
    return dict(raw_domain_count=len(hits), profile_evidence=evidence, overlapping_domain_rows=overlaps,
        ordered_PF01058_PF02906_status='predicted_same_polypeptide_domain_organization' if pairs else 'not_supported_by_this_profile_panel',
        ordered_PF01058_PF02906_candidates=list(pairs.values()), evidence_scope=SCOPE)


def annotate_domains(runs, output, root=ROOT):
    output = Path(output).resolve(); runs = [Path(p).resolve() for p in runs]
    if not runs or len(set(runs)) != len(runs): raise ValueError('Missing or duplicate source run')
    if output.exists(): raise FileExistsError('Use a fresh query-domain output directory')
    if any(output == run or run in output.parents for run in runs): raise ValueError('Output must be outside source annotation runs')
    policy, models, combined, resources = load_models(root)
    owners, queries, inputs = collect_sources(sorted(runs))
    output.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    record = output/'query_domain_manifest.json'
    report = dict(status='running', scope=SCOPE, source_runs=[str(p) for p in sorted(runs)], input_sha256=inputs,
        resource_sha256=resources, protein_owners=len(owners), unique_sequences=len(queries), profile_models=len(models),
        min_model_coverage=policy['min_model_coverage'], evalue_scope='Batch-dependent HMMER independent-domain E values; decisions use curated GA and model coverage, not population FPR.',
        HMM_query_normalization='U_to_X_same_length_only; query_id hashes original uppercase sequence', software=subprocess.check_output(['hmmsearch', '-h'], text=True).splitlines()[1])
    record.write_text(json.dumps(report, indent=2)+'\n')
    try:
        fasta = output/'queries_U_to_X.faa'
        SeqIO.write([SeqRecord(Seq(s.replace('U', 'X')), id=qid, description='query_domain_only') for qid, s in sorted(queries.items())], fasta, 'fasta')
        domains = scan_profile_domains(combined, fasta, models, policy['min_model_coverage'], output) if queries else {}
        if not queries:
            (output/'domains.tblout').write_text('# no query sequences\n')
        raw_targets = load_targets(output/'domains.tblout', {profile: profile for profile in models}) if queries else []
        targets = []
        for target in raw_targets:
            if target['protein_id'] not in queries:
                raise ValueError('Whole-sequence target query is not owned by source proteins')
            model = models[target['profile']]
            target = dict(target, accession=model['accession'], model_length=model['model_length'],
                          query_length=len(queries[target['protein_id']]), model_sha256=model['sha256'], source=model['url'])
            target['query_id'] = target.pop('protein_id')
            targets.append(target)
        if any(t['query_id'] not in queries for t in targets):
            raise ValueError('Whole-sequence target query is not owned by source proteins')
        target_by_query = {}
        for target in targets:
            target_by_query.setdefault(target['query_id'], []).append(target)
        validate_domains(domains, queries, models, policy['min_model_coverage'])
        rows = []; summaries = []
        for owner in owners:
            hits = domains.get(owner['query_id'], [])
            rows.extend(dict(**owner, domain_row=i+1, **d) for i, d in enumerate(hits))
            summaries.append(dict(**owner, **organization(hits, policy)))
        write_tsv(output/'query_domains.tsv', rows, [*OWNER_FIELDS, 'domain_row', *DOMAIN_FIELDS])
        target_rows=[]
        for owner in owners:
            for target in target_by_query.get(owner['query_id'], []):
                target_row=dict(owner); target_row.update(target); target_rows.append(target_row)
        write_tsv(output/'query_target_signals.tsv', target_rows, [*OWNER_FIELDS, *TARGET_FIELDS])
        (output/'target_evidence.json').write_text(json.dumps(target_by_query, indent=2)+'\n')
        summary_fields = list(organization([], policy))
        write_tsv(output/'query_domain_organizations.tsv', [{k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in s.items()} for s in summaries], [*OWNER_FIELDS, *summary_fields])
        (output/'domain_evidence.json').write_text(json.dumps(domains, indent=2)+'\n')
        if any(sha(p) != digest for p, digest in {**inputs, **resources}.items()): raise ValueError('Source or domain model changed during execution')
        report.update(status='complete', inputs_and_models_unchanged=True, elapsed_seconds=time.monotonic()-start,
            domain_rows=len(rows), owners_with_supported_ordered_domains=sum(bool(s['ordered_PF01058_PF02906_candidates']) for s in summaries),
            target_rows=len(target_rows),
            artifacts={p.name: sha(p) for p in output.iterdir() if p.is_file() and p != record})
    except Exception as error:
        report.update(status='failed', error=str(error), elapsed_seconds=time.monotonic()-start)
        record.write_text(json.dumps(report, indent=2)+'\n'); raise
    record.write_text(json.dumps(report, indent=2)+'\n')
    return report


def verify_completed(output, root=ROOT):
    output = Path(output); report = json.loads((output/'query_domain_manifest.json').read_text())
    if report['status'] != 'complete' or not report['inputs_and_models_unchanged']: raise ValueError('Incomplete query-domain execution')
    policy, models, combined, resources = load_models(root)
    if resources != report['resource_sha256']: raise ValueError('Query-domain resources changed after execution')
    owners, queries, inputs = collect_sources([Path(p) for p in report['source_runs']])
    if inputs != report['input_sha256']: raise ValueError('Query-domain source changed after execution')
    if any(sha(output/name) != digest for name, digest in report['artifacts'].items()): raise ValueError('Query-domain output changed after execution')
    domains = json.loads((output/'domain_evidence.json').read_text()); validate_domains(domains, queries, models, policy['min_model_coverage'])
    raw_targets = load_targets(output/'domains.tblout', {profile: profile for profile in models})
    targets = []
    for target in raw_targets:
        if target['protein_id'] not in queries:
            raise ValueError('Whole-sequence target query is not owned by source proteins')
        model = models[target['profile']]
        target = dict(target, accession=model['accession'], model_length=model['model_length'],
                      query_length=len(queries[target['protein_id']]), model_sha256=model['sha256'], source=model['url'])
        target['query_id'] = target.pop('protein_id')
        targets.append(target)
    expected_targets = {}
    for target in targets:
        if target['query_id'] not in queries: raise ValueError('Whole-sequence target query is not owned by source proteins')
        expected_targets.setdefault(target['query_id'], []).append(target)
    if json.loads((output/'target_evidence.json').read_text()) != expected_targets: raise ValueError('Target evidence differs from coordinate evidence')
    target_rows = []
    for owner in owners:
        for target in expected_targets.get(owner['query_id'], []):
            row = dict(owner); row.update(target); target_rows.append(row)
    if table(output/'query_target_signals.tsv') != [{k: str(v) for k, v in r.items() if k in OWNER_FIELDS + TARGET_FIELDS} for r in target_rows]:
        raise ValueError('Query target rows differ from coordinate evidence')
    expected = [dict(**o, **organization(domains.get(o['query_id'], []), policy)) for o in owners]
    rows = table(output/'query_domain_organizations.tsv')
    serialized = [{k: json.dumps(v) if isinstance(v, (list, dict)) else str(v) for k, v in r.items()} for r in expected]
    if rows != serialized: raise ValueError('Query domain organization differs from coordinate evidence')
    raw = [dict(**o, domain_row=i+1, **d) for o in owners for i, d in enumerate(domains.get(o['query_id'], []))]
    if table(output/'query_domains.tsv') != [{k: str(v) for k, v in r.items()} for r in raw]:
        raise ValueError('Query domain rows differ from coordinate evidence')
    if (report['protein_owners'] != len(owners) or report['unique_sequences'] != len(queries) or report['profile_models'] != len(models) or
        report['domain_rows'] != len(raw) or report['target_rows'] != sum(len(expected_targets.get(o['query_id'], [])) for o in owners) or report['owners_with_supported_ordered_domains'] != sum(bool(r['ordered_PF01058_PF02906_candidates']) for r in expected)):
        raise ValueError('Query-domain manifest counts differ from actual evidence')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', nargs='+', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); args = parser.parse_args()
    report = annotate_domains(args.runs, args.out); verify_completed(args.out)
    print(json.dumps({k: report[k] for k in ['status', 'protein_owners', 'unique_sequences', 'profile_models', 'domain_rows', 'owners_with_supported_ordered_domains', 'elapsed_seconds']}, indent=2))


if __name__ == '__main__': main()
