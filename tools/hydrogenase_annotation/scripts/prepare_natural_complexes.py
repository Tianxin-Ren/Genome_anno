#!/usr/bin/env python3
"""Freeze experimentally scoped natural complexes before annotation; never train on them."""
import argparse
import copy
import hashlib
import json
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from Bio import SeqIO
from annotate import ROOT, extract_genbank, write_tsv
from catalytic_confirmation import sha


def anonymous_genes(genes):
    anonymous = copy.deepcopy(genes)
    mapping = {}
    contigs = {c: 'contig_' + str(i + 1) for i, c in enumerate(dict.fromkeys(g['contig'] for g in genes))}
    for i, gene in enumerate(anonymous):
        mapping[gene['id']] = 'query_' + str(i + 1)
        gene.update(id=mapping[gene['id']], contig=contigs[gene['contig']], gene='', aliases='', protein_accession='')
    return anonymous, mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    policy_path = ROOT / 'config/complex_function_validation.json'
    policy = json.loads(policy_path.read_text())
    raw = ROOT / 'references/raw/round22_function_truth'
    uniprot = {p.stem: json.loads(p.read_text()) for p in (raw / 'uniprot').glob('*.json') if p.stem not in ['manifest', 'additional_manifest', 'Hnd_genes_search', 'Ralbus_hyd_search']}
    uniprot.update({r['primaryAccession']: r for r in json.loads((raw / 'uniprot/Hnd_genes_search.json').read_text())['results']})
    references = json.loads((ROOT / 'references/processed/references.json').read_text())
    prototypes = json.loads((ROOT / 'references/processed/catalytic_confirmation/prototypes.json').read_text())
    units = json.loads((ROOT / 'references/validation/round09_identity_corpus/units.json').read_text())
    sources, truth, audit, source_hashes = {}, [], [], {str(policy_path.relative_to(ROOT)): sha(policy_path)}
    for case in policy['cases']:
        source = case['source']; path = raw / 'genomes' / (source + '.gb')
        if source not in sources:
            genes = extract_genbank(path)[0]
            anonymous, mapping = anonymous_genes(genes)
            sources[source] = dict(mode=case['mode'], input_file=str(path.relative_to(ROOT)), input_sha256=sha(path),
                                   original_genes=len(genes), anonymous_genes=anonymous, source_to_query=mapping)
        source_hashes[str(path.relative_to(ROOT))] = sha(path)
        paper = ROOT / case['paper']
        if ET.parse(paper).find('.//body') is None:
            raise ValueError('Full experimental paper body required: ' + str(paper))
        source_hashes[case['paper']] = sha(paper)
        genes = {g['id']: g for g in extract_genbank(path)[0]}
        members = []
        for role, (locus, accession) in case['members'].items():
            gene = genes[locus]; entry = uniprot[accession]; sequence = gene['sequence']
            if sequence != entry['sequence']['value']:
                raise ValueError('Native CDS/UniProt mismatch: ' + locus)
            entry_path = raw / 'uniprot' / (accession + '.json')
            if not entry_path.exists(): entry_path = raw / 'uniprot/Hnd_genes_search.json'
            source_hashes[str(entry_path.relative_to(ROOT))] = sha(entry_path)
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            members.append(dict(role=role, source_locus=locus, protein_accession=gene['protein_accession'],
                uniprot=accession, query_id=sources[source]['source_to_query'][locus], sequence_sha256=digest,
                length=len(sequence), start=gene['start'], end=gene['end'], strand=gene['strand']))
            exact_refs = [r['reference_id'] for r in references if r['sequence_sha256'] == digest]
            exact_prototypes = [r['prototype_id'] for r in prototypes if r['sequence_sha256'] == digest]
            genus_units = [r['unit_id'] for r in units if case['genus'] in r.get('genus_names', [])]
            audit.append(dict(case_id=case['id'], role=role, query_id=members[-1]['query_id'],
                source=source, sequence_sha256=digest, exact_functional_references=';'.join(exact_refs),
                exact_confirmation_prototypes=';'.join(exact_prototypes), identity_corpus_genus_units=len(genus_units),
                independence='development_only_not_independent',
                limitation='Exact/genus exposure inventory only; absence of an exact match does not prove homolog/core or taxon independence. No final-test eligibility.'))
        truth.append(dict(case, members=members, query_anchor=sources[source]['source_to_query'][case['anchor']] if case['anchor'] else None))
    (out / 'sources.json').write_text(json.dumps(sources, indent=2) + '\n')
    (out / 'truth.json').write_text(json.dumps(truth, indent=2) + '\n')
    write_tsv(out / 'source_exposure.tsv', audit)
    freeze = dict(status='development_truth_frozen_before_annotation', created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(policy_path), sources=source_hashes,
        artifacts={p.name: sha(p) for p in out.iterdir() if p.is_file()},
        source_independence='Not claimed. Selected after source inspection, with previously exposed genera/structures and possible identity training overlap.',
        targets=policy['targets'], cases=len(truth), positive_complexes=sum(c['hydrogenase'] for c in truth),
        non_H2_center_contrasts=sum(not c['hydrogenase'] for c in truth), distinct_input_sources=len(sources),
        proteins=sum(s['original_genes'] for s in sources.values()),
        members=sum(len(c['members']) for c in truth),
        limits=['Native annotated translations, not independent de novo gene prediction.',
            'Four positive complexes and two mechanism contrasts cannot establish population accuracy.',
            'Other genomic anchors are ungradable; physiological flux and expression remain unknown.',
            '5XFA FeS coordinate truth is not opened or used in L4/L5 evaluation.'])
    (out / 'freeze.json').write_text(json.dumps(freeze, indent=2) + '\n')
    print(json.dumps({k: freeze[k] for k in ['status', 'cases', 'positive_complexes', 'non_H2_center_contrasts', 'proteins', 'members']}, indent=2))


if __name__ == '__main__': main()
