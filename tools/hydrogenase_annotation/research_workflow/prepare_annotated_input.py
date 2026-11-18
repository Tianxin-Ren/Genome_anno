"""Preserve explicitly annotated Sec when legacy GenBank translations encode it as X."""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import Location, SimpleLocation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from annotate import extract_genbank
from input_translation import cds_coordinates, normalized_translation, verified_origin_translation


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def annotated_translation(record, feature):
    return normalized_translation(record, feature)


def prepare(genbank, output):
    genbank = Path(genbank).resolve(); output = Path(output).resolve()
    if output.exists(): raise FileExistsError('Use a fresh annotated-input directory')
    before = sha(genbank)
    with genbank.open() as handle:
        records = list(SeqIO.parse(handle, 'genbank'))
    genes, _ = extract_genbank(genbank)
    features = {}; lengths = {r.id: len(r) for r in records}; topologies = {}
    for record in records:
        topology = record.annotations.get('topology')
        topologies[record.id] = topology if topology in {'linear', 'circular'} else None
        for feature in record.features:
            if feature.type != 'CDS' or 'translation' not in feature.qualifiers: continue
            if 'pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers: continue
            try:
                start, end, strand, _partial = cds_coordinates(record, feature, topologies[record.id])
            except ValueError as exc:
                raise ValueError('Unsupported non-local or multi-part CDS') from exc
            sequence, _edits = annotated_translation(record, feature)
            if end > len(record) and not verified_origin_translation(record, feature, sequence):
                raise ValueError('Origin-spanning CDS translation is not verified')
            key = (record.id, start, end, strand)
            if key in features: raise ValueError('Repeated CDS coordinate ownership')
            features[key] = (record, feature)
    fasta = []; gff = ['##gff-version 3', *[f'##sequence-region {cid} 1 {length}' for cid, length in lengths.items()]]
    gff.extend(f'##topology {cid} {topology}' for cid, topology in topologies.items() if topology)
    changes = []
    for gene in genes:
        if gene.get('coordinate_evidence', 'genomic_CDS') != 'genomic_CDS' or int(gene.get('start', 0)) < 1:
            raise ValueError('CDS lacks verified genomic coordinates')
        record, feature = features[(gene['contig'], gene['start'], gene['end'], gene['strand'])]
        if int(feature.qualifiers.get('codon_start', ['1'])[0]) != 1: raise ValueError('Partial-frame CDS cannot be emitted as phase-zero GFF')
        sequence, edits = annotated_translation(record, feature)
        fasta.append('>'+gene['id']+'\n'+sequence+'\n')
        gff.append('\t'.join([gene['contig'], 'explicit_source_CDS', 'CDS', str(gene['start']), str(gene['end']), '.', '+' if gene['strand']==1 else '-', '0',
            'ID='+gene['id']+';partial='+('01' if gene['partial'] else '00')]))
        changes.extend(dict(protein_id=gene['id'], **e) for e in edits)
    if sha(genbank) != before: raise ValueError('Source changed during normalization')
    output.mkdir(parents=True, exist_ok=False)
    (output/'proteins.faa').write_text(''.join(fasta)); (output/'genes.gff').write_text('\n'.join(gff)+'\n')
    manifest = dict(status='explicit_source_translation_prepared', source_genbank=str(genbank), source_sha256=before,
        protein_count=len(genes), changes=changes, input_script_sha256=sha(Path(__file__)),
        artifacts={name:sha(output/name) for name in ['proteins.faa', 'genes.gff']},
        scope='Annotated contiguous prokaryotic CDS only; explicit Sec exception plus codon and source note required for OTHER. No inference of Sec from unannotated FNA or blanket stop-codon recoding.')
    (output/'input_preparation_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--genbank', type=Path, required=True); parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); print(json.dumps(prepare(args.genbank, args.out), ensure_ascii=False))
