"""Small, evidence-only helpers for recording input translation models.

These functions deliberately do not infer a genetic code from proteins.  In
particular, a protein FASTA supplied without its originating annotation keeps
an unknown effective translation table.
"""
import re
import hashlib
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqFeature import CompoundLocation, ExactPosition, Location

from circular_coordinates import interval_segments


def cds_coordinates(record, feature, topology):
    """Map a local CDS to a one-based interval, including an exact origin join."""
    location = feature.location
    if location is None or any(part.ref or part.ref_db for part in location.parts):
        raise ValueError("CDS coordinates are not local to the sequence record")
    if not isinstance(location, CompoundLocation):
        start, end = int(location.start) + 1, int(location.end)
        if not 1 <= start <= end <= len(record):
            raise ValueError("CDS coordinates are outside the sequence record")
        partial = any(not isinstance(p, ExactPosition) for p in (location.start, location.end))
        return start, end, location.strand or 0, partial

    strand = location.strand
    parts = location.parts
    if (topology != "circular" or location.operator != "join" or len(parts) != 2
            or strand not in (1, -1) or any(not isinstance(p, ExactPosition)
                                           for part in parts for p in (part.start, part.end))):
        raise ValueError("CDS join does not establish an exact circular interval")
    tail, head = parts if strand == 1 else parts[::-1]
    start, end = int(tail.start) + 1, len(record) + int(head.end)
    if (head.start != 0 or tail.end != len(record)
            or interval_segments(start, end, len(record)) != [(start, len(record)), (1, int(head.end))]):
        raise ValueError("CDS join is not contiguous across the circular origin")
    return start, end, strand, False


def verified_origin_translation(record, feature, sequence):
    """Verify the supplied origin-spanning CDS under its explicit source code."""
    qualifiers = feature.qualifiers
    try:
        table_id, = qualifiers["transl_table"]
        table = CodonTable.unambiguous_dna_by_id[int(table_id)]
        if qualifiers.get("codon_start", ["1"]) != ["1"]:
            return False
        dna = feature.extract(record.seq).upper()
        if not qualifiers.get("transl_except"):
            return str(dna.translate(table=int(table_id), cds=True)) == sequence.rstrip("*")
        if len(dna) % 3 or str(dna[:3]) not in table.start_codons or str(dna[-3:]) not in table.stop_codons:
            return False
        normalized, sites = normalized_translation(record, feature, include_unchanged_sec=True)
        positions = {s["protein_position"] for s in sites}
        internal_stops = {i // 3 + 1 for i in range(0, len(dna) - 3, 3)
                          if str(dna[i:i+3]) in table.stop_codons}
        return (normalized == sequence and len(positions) == len(qualifiers["transl_except"])
                and positions == internal_stops == {i for i, aa in enumerate(normalized, 1) if aa == "U"})
    except (KeyError, ValueError, CodonTable.TranslationError):
        return False


def normalized_translation(record, feature, *, include_unchanged_sec=False):
    """Return a GenBank CDS translation, changing only explicitly annotated Sec.

    Legacy records encode Sec as X and use ``OTHER`` plus a note.  No stop
    recoding is attempted unless the source exception, note, frame and TGA
    codon all agree.
    """
    original = feature.qualifiers['translation'][0].upper()
    sequence = list(original)
    if not feature.qualifiers.get('transl_except'):
        return original, []
    offset = int(feature.qualifiers.get('codon_start', ['1'])[0]) - 1
    coordinates = list(feature.location)[offset:]
    dna = str(feature.extract(record.seq))[offset:].upper()
    table = int(feature.qualifiers.get('transl_table', ['11'])[0])
    translated = str(feature.extract(record.seq)[offset:].translate(table=table))
    if translated.endswith('*'):
        translated = translated[:-1]
    if len(translated) != len(sequence):
        raise ValueError('Translation/CDS length conflict for explicit translation exception')
    if translated and sequence[0] == 'M' and dna[:3] in CodonTable.unambiguous_dna_by_id[table].start_codons:
        translated = 'M' + translated[1:]
    notes = feature.qualifiers.get('note', []); edits = []; unresolved = set()
    supported_sec = {i for i, aa in enumerate(sequence) if aa == 'U' and dna[3*i:3*i+3] == 'TGA'}
    for exception in feature.qualifiers.get('transl_except', []):
        match = re.fullmatch(r'\(pos:(.+),(?:aa:)?([^,()]+)\)', exception)
        if not match:
            raise ValueError('Unsupported transl_except syntax')
        location = Location.fromstring(match[1], len(record), circular=record.annotations.get("topology") == "circular")
        genomic = list(location)
        if len(genomic) != 3 or any(part.ref or part.ref_db or
                not isinstance(part.start, ExactPosition) or not isinstance(part.end, ExactPosition)
                for part in location.parts):
            raise ValueError('Translation exception does not identify one exact local codon')
        try:
            indices = [coordinates.index(p) for p in genomic]
        except ValueError as error:
            raise ValueError('Translation exception lies outside CDS') from error
        start = min(indices)
        if start % 3 or sorted(indices) != list(range(start, start+3)):
            raise ValueError('Translation exception is out of frame')
        position = start // 3
        if position >= len(sequence):
            raise ValueError('Translation exception has no protein residue')
        aa = match[2].strip().lower()
        coordinate_text = f'{int(location.start)+1}..{int(location.end)}'
        explicit_note = any('selenocysteine' in n.lower() and coordinate_text in n for n in notes)
        sec = aa in ['sec', 'selenocysteine'] or (aa == 'other' and explicit_note)
        if not sec:
            if sequence[position] == 'X':
                unresolved.add(position)
            continue
        if dna[start:start+3] != 'TGA':
            raise ValueError('Annotated Sec is not encoded by TGA')
        if sequence[position] not in ['X', 'U']:
            raise ValueError('Annotated Sec conflicts with supplied amino acid')
        supported_sec.add(position)
        original_aa = sequence[position]
        if original_aa == 'X':
            sequence[position] = 'U'
        if original_aa == 'X' or include_unchanged_sec:
            edits.append(dict(protein_position=position+1, original=original_aa, normalized='U',
                genomic_positions_1based=[p+1 for p in genomic], codon='TGA', source_exception=exception,
                source_notes=[n for n in notes if 'selenocysteine' in n.lower() and coordinate_text in n],
                scope='Explicit source coding annotation; not evidence of query metal occupancy or activity.'))
    for i, (aa, standard) in enumerate(zip(sequence, translated)):
        if aa == standard or (i in supported_sec and aa == 'U' and standard == '*') or (i in unresolved and aa == 'X'):
            continue
        raise ValueError(f'Unexplained protein/CDS conflict at residue {i+1}')
    return ''.join(sequence), edits


def prodigal_models(path):
    """Return effective Prodigal tables and CDS counts keyed by contig."""
    tables, counts, current = {}, {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('# Sequence Data:'):
            current = None
        elif line.startswith('# Model Data:'):
            match = re.search(r'(?:^|;)transl_table=(\d+)(?:;|$)', line)
            if not match:
                raise ValueError('Prodigal model lacks a translation table')
            current = int(match.group(1))
        elif line and not line.startswith('#'):
            fields = line.split('\t')
            if len(fields) == 9 and fields[2] == 'CDS' and 'Prodigal' in fields[1]:
                if fields[0] in tables and tables[fields[0]] != current:
                    raise ValueError('Conflicting translation models on one contig')
                tables[fields[0]] = current
                counts[fields[0]] = counts.get(fields[0], 0) + 1
    return tables, counts


def _genbank_evidence(path):
    tables, exceptions, features, cds = {}, [], [], 0
    with Path(path).open() as handle:
        for record in SeqIO.parse(handle, 'genbank'):
            for feature in record.features:
                if feature.type != 'CDS':
                    continue
                cds += 1
                q = feature.qualifiers
                value = q.get('transl_table', [None])[0]
                tables.setdefault(record.id, set())
                if value is not None:
                    try:
                        value = int(value)
                    except (TypeError, ValueError):
                        value = str(value)
                    tables.setdefault(record.id, set()).add(value)
                feature_id = q.get('protein_id', q.get('locus_tag', ['']))[0]
                item = {'contig': record.id, 'protein_id': feature_id,
                        'location': str(feature.location), 'translation_table': value}
                if q.get('transl_except'):
                    if 'translation' not in q:
                        exceptions.append(dict(item, transl_except=[], valid_sec=False,
                            exception_status='translation_qualifier_missing'))
                    else:
                        _, edits = normalized_translation(record, feature, include_unchanged_sec=True)
                        for edit in edits:
                            check = dict(edit)
                            check['valid_sec'] = True
                            exceptions.append(dict(item, transl_except=[check]))
                        if not edits:
                            exceptions.append(dict(item, transl_except=[], valid_sec=False,
                                exception_status='not_an_explicit_sec_exception'))
                item['translation_contains_U'] = 'U' in q.get('translation', [''])[0]
                features.append(item)
    used = {contig: (None if not values else sorted(values)[0] if len(values) == 1 else sorted(values))
            for contig, values in tables.items()}
    return used, exceptions, features, cds


def evidence(input_files=(), mode='assembly', requested=None, generated_gff=None):
    """Build manifest-safe translation evidence without changing translation."""
    files = [Path(p) for p in input_files]
    gff = Path(generated_gff) if generated_gff else None
    if gff and gff.exists():
        tables, counts = prodigal_models(gff)
        return dict(used=tables, table_evidence={'status': 'verified_prodigal_model_data',
            'source': str(gff.resolve()), 'source_sha256': hashlib.sha256(gff.read_bytes()).hexdigest(),
            'contigs': tables, 'coding_features': counts},
            exception_evidence={'status': 'not_checked_in_Prodigal_input', 'records': []})
    genbank = next((p for p in files if p.suffix.lower() in ('.gb', '.gbk', '.gbff', '.genbank')), None)
    if genbank and genbank.exists():
        tables, exceptions, features, cds = _genbank_evidence(genbank)
        return dict(used=tables, table_evidence={'status': 'source_qualifiers', 'source': str(genbank.resolve()),
            'source_sha256': hashlib.sha256(genbank.read_bytes()).hexdigest(), 'contigs': tables,
            'features': features, 'cds_count': cds},
            exception_evidence={'status': 'explicit_source_qualifiers', 'source': str(genbank.resolve()),
                'source_sha256': hashlib.sha256(genbank.read_bytes()).hexdigest(), 'records': exceptions})
    # GFF attributes can be evidence, but an accompanying FAA alone cannot.
    source_gff = next((p for p in files if p.suffix.lower() in ('.gff', '.gff3')), None)
    if source_gff and source_gff.exists():
        prep = source_gff.parent / 'input_preparation_manifest.json'
        prepared = None
        if prep.exists():
            import json
            prepared = json.loads(prep.read_text())
            source = Path(prepared.get('source_genbank', ''))
            valid = source.exists() and hashlib.sha256(source.read_bytes()).hexdigest() == prepared.get('source_sha256')
            for name, digest in prepared.get('artifacts', {}).items():
                artifact = prep.parent / name
                valid = valid and artifact.exists() and hashlib.sha256(artifact.read_bytes()).hexdigest() == digest
            if not valid:
                prepared = None
        # A supplied Prodigal GFF retains the effective model in its headers.
        if any(line.startswith('# Model Data:') for line in source_gff.read_text().splitlines()):
            tables, counts = prodigal_models(source_gff)
            return dict(used=tables, table_evidence={'status': 'verified_prodigal_model_data',
                'source': str(source_gff.resolve()), 'source_sha256': hashlib.sha256(source_gff.read_bytes()).hexdigest(),
                'contigs': tables, 'coding_features': counts},
                exception_evidence={'status': 'not_checked_in_Prodigal_input', 'records': [],
                    'preparation_manifest': prepared} if prepared else {'status': 'not_checked_in_Prodigal_input', 'records': []})
        tables = {}
        for line in source_gff.read_text().splitlines():
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) != 9 or fields[2] != 'CDS':
                continue
            match = re.search(r'(?:^|;)transl_table=(\d+)(?:;|$)', fields[8])
            if match:
                table = int(match.group(1)); previous = tables.get(fields[0])
                tables[fields[0]] = table if previous in (None, table) else 'mixed'
        status = 'source_gff_attributes' if tables else 'source_gff_table_unknown'
        return dict(used=tables, table_evidence={'status': status, 'source': str(source_gff.resolve()),
            'source_sha256': hashlib.sha256(source_gff.read_bytes()).hexdigest(), 'contigs': tables},
            exception_evidence={'status': 'explicit_source_preparation_manifest' if prepared else 'not_available',
                'source': str(prep.resolve()) if prepared else None,
                'source_sha256': hashlib.sha256(prep.read_bytes()).hexdigest() if prepared else None,
                'records': (prepared or {}).get('changes', []) if prepared else []})
    return dict(used={}, table_evidence={'status': 'unknown_no_source_evidence', 'source': None, 'contigs': {}},
                exception_evidence={'status': 'not_available', 'records': []})
