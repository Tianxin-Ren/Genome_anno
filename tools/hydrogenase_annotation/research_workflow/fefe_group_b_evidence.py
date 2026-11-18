"""Resolve an uncertain PF02906 column by two assayed FeFe reference alignments."""
import hashlib
import json
import subprocess
from pathlib import Path

from annotate import BLAST_FIELDS, load_blast
from catalytic_confirmation import core_similarity
from tool_paths import resolve_tool

RULE = 'FeFe-assayed-lineage-concordant-ligands-v1'
MODEL = 'research_workflow/fefe_group_b_model.json'
MAX_EVALUE = 1e-10


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _profile_eligible(record, feature):
    return (record.get('blockers') == ['low_profile_posterior']
            and [s.get('pp') for s in record['sites']] == ['9', '*', '*', '*']
            and all(s['residue_status'] == 'present' and s['within_domain'] for s in record['sites'])
            and feature.get('source_active_center', {}).get('direct_ligand_status') in {'present', 'ambiguous'}
            and feature.get('source_active_center', {}).get('center_organization') == 'single_polypeptide')


def _site_positions(hit, reference_sites):
    query, reference = hit['qstart'] - 1, hit['sstart'] - 1
    result = []
    for q, r in zip(hit['qseq'], hit['sseq']):
        query += q != '-'
        reference += r != '-'
        if r != '-' and reference in reference_sites:
            result.append(query if q == r == 'C' else None)
    return result


def prepare(annotation, output, records, features, queries, parameters, root, *, replay=False):
    from fefe_profile_ligands import _fasta, _stockholm, _site_map, _load_model
    def _validate_hit(hit, queries, references):
        if hit['qseqid'] not in queries or hit['sseqid'] not in references: raise ValueError('Unknown BLAST sequence ownership')
        if len(hit['qseq'])!=len(hit['sseq']): raise ValueError('BLAST alignment lengths differ')
        for field, sequence, lo, hi in [('qseq',queries[hit['qseqid']],hit['qstart'],hit['qend']),('sseq',references[hit['sseqid']],hit['sstart'],hit['send'])]:
            observed=hit[field].replace('-','').replace('U','X').replace('O','X')
            if not 1<=lo<=hi<=len(sequence) or observed!=sequence[lo-1:hi].replace('U','X').replace('O','X'):
                raise ValueError('BLAST alignment differs from bound '+field+' sequence')

    selected = {q: record for q, record in records.items() if _profile_eligible(record, features[q])}
    if not selected:
        return {}
    background = annotation/'confirmation/prototype_blast.tsv'
    if not background.is_file():
        return {}
    background_manifest = background.with_name('feature_manifest.json')
    artifact_hashes = json.loads(background_manifest.read_text())['artifacts']
    if artifact_hashes.get(background.name) != sha(background):
        raise ValueError('FeFe competing prototype search differs from the bound feature run')
    model_path = root/MODEL
    model = json.loads(model_path.read_text())
    if model['version'] != RULE:
        raise ValueError('FeFe assayed ligand rule version mismatch')
    fasta = root/model['sequence_bank']
    refs = _fasta(fasta)
    sites = {ref: row['sites'] for ref, row in model['assayed_references'].items()}
    if set(sites) != set(refs) or len(sites) != 2:
        raise ValueError('Incomplete FeFe assayed reference bank')
    profile, _, template, _ = _load_model(root)
    seqs, posteriors, rf = _stockholm((root/model['source_alignment']).read_text())
    columns = [i for i, letter in enumerate(rf, 1) if letter.lower() == 'x']
    ligand_columns = {columns[i - 1] for i in profile['expected_match_columns']}
    for ref, seq in {**refs, 'ref_1FEH': template['sequence']}.items():
        mapped = list(_site_map(seqs[ref], posteriors[ref], ligand_columns).values())
        expected = sites.get(ref, profile['reference_positions'])
        if (''.join(c for c in seqs[ref] if c not in '-.').upper() != seq
                or [site['query_position'] for site in mapped] != expected
                or any(site['residue'] != 'C' or site['pp'] not in {'9', '*'} for site in mapped)):
            raise ValueError('Assayed reference ligand mapping differs')
    prototype_path = root/'references/processed/catalytic_confirmation/prototypes.json'
    negative = {p['prototype_id']: p['sequence'] for p in json.loads(prototype_path.read_text()) if p['label'] == 'negative'}
    if not negative:
        raise ValueError('FeFe assayed ligand route requires negative homolog competition')
    cutoff = dict(blast_evalue=MAX_EVALUE, min_identity=0, min_query_coverage=0, min_reference_coverage=0)
    query = output/'group_B_queries.faa'
    raw = output/'group_B_assay_blast.tsv'
    binary = resolve_tool(root, 'blastp')
    deps = [model_path, background, background_manifest, prototype_path, Path(__file__), binary,
            *[root/p for p in model['files_sha256']]]
    inputs = {str(p):sha(p) for p in deps}
    normalized = {q: queries[q].replace('U', 'X').replace('O', 'X') for q in selected}
    if not replay:
        query.write_text(''.join(f'>{q}\n{normalized[q]}\n' for q in sorted(selected)))
        subprocess.run([str(binary), '-query', str(query), '-subject', str(fasta),
                        '-evalue', str(MAX_EVALUE), '-max_hsps', '1', '-seg', 'no',
                        '-outfmt', '6 '+' '.join(BLAST_FIELDS), '-out', str(raw)], check=True, capture_output=True)
    domains = {q: range(r['domain_ownership']['query_start'], r['domain_ownership']['query_end'] + 1) for q, r in selected.items()}
    best_negative = {}
    for h in load_blast(background, cutoff):
        q, ref = h['qseqid'], h['sseqid']
        if q in selected and ref in negative and h['evalue'] <= MAX_EVALUE:
            _validate_hit(h, normalized, {ref:negative[ref].replace('U', 'X').replace('O', 'X')})
            value = core_similarity(h, domains[q])['similarity']
            best_negative[q] = max(best_negative.get(q, 0), value)
    positives = {q:{} for q in selected}
    for h in load_blast(raw, cutoff):
        q, ref = h['qseqid'], h['sseqid']
        _validate_hit(h, normalized, refs)
        if ref in positives[q] or h['evalue'] > MAX_EVALUE:
            raise ValueError('Unexpected FeFe assayed reference HSP')
        positions = _site_positions(h, sites[ref])
        if positions == [site['query_position'] for site in selected[q]['sites']]:
            positives[q][ref] = dict(core_similarity(h, domains[q]), query_ligand_positions=positions)
    if any(sha(p) != h for p, h in inputs.items()):
        raise ValueError('FeFe assayed ligand evidence changed during evaluation')
    for q, matches in positives.items():
        negative_score = best_negative.get(q, 0)
        if (matches.keys() == refs.keys() and all(m['core_coverage'] >= parameters['min_domain_coverage']
                and m['similarity'] >= parameters['min_core_similarity']
                and m['similarity'] - negative_score >= parameters['min_competitor_margin'] for m in matches.values())):
            records[q]['assayed_lineage_corroboration'] = dict(rule=RULE, model_sha256=inputs[str(model_path)],
                references=matches, best_known_negative_similarity=negative_score,
                negative_prototype_count=len(negative), scope=model['scope'])
    return inputs


def record_supported(feature, record, parameters):
    proof = record.get('assayed_lineage_corroboration', {})
    if proof.get('rule') != RULE or not _profile_eligible(record, feature):
        return False
    path = Path(__file__).resolve().parents[1]/MODEL
    model = json.loads(path.read_text())
    positions = [site['query_position'] for site in record['sites']]
    matches = proof.get('references', {})
    return bool(proof.get('model_sha256') == sha(path) and proof.get('negative_prototype_count', 0) > 0
        and matches.keys() == model['assayed_references'].keys()
        and all(m['query_ligand_positions'] == positions and m['evalue'] <= MAX_EVALUE
                and m['core_coverage'] >= parameters['min_domain_coverage']
                and m['similarity'] >= parameters['min_core_similarity']
                and m['similarity'] - proof['best_known_negative_similarity'] >= parameters['min_competitor_margin']
                for m in matches.values()))
