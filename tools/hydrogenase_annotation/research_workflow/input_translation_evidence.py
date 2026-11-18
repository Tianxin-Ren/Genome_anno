"""Record effective Prodigal models and translation-sensitive ligand observations."""
import json
import re
from pathlib import Path
from Bio.Seq import Seq
from circular_coordinates import circular_slice, valid_cds_interval
from localization import sha
from nucleotide_input import fasta_records


def prodigal_tables(path):
    tables = {}; counts = {}; current = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('# Sequence Data:'): current = None
        elif line.startswith('# Model Data:'):
            match = re.search(r'(?:^|;)transl_table=(\d+)(?:;|$)', line)
            if not match: raise ValueError('Prodigal model lacks a translation table')
            current = int(match[1])
        elif line and not line.startswith('#'):
            fields = line.split('\t')
            if len(fields) != 9 or fields[2] != 'CDS': continue
            if 'Prodigal' not in fields[1]: continue
            if fields[0] in tables and tables[fields[0]] != current: raise ValueError('Conflicting translation models on one contig')
            tables[fields[0]] = current; counts[fields[0]] = counts.get(fields[0], 0)+1
    return tables, counts


def audit(annotation, genes, sites):
    annotation = Path(annotation); manifest = json.loads((annotation/'run_manifest.json').read_text()); bindings = {}
    gff = annotation/'prodigal.gff'
    if not gff.exists():
        supplied = [r for r in manifest.get('input_files', []) if Path(r['path']).suffix.lower() in ['.gff', '.gff3']]
        if len(supplied) == 1:
            gff = Path(supplied[0]['path'])
            if not gff.exists(): return dict(status='provided_GFF_unavailable', sensitive_sites=[]), bindings
            if sha(gff) != supplied[0]['sha256']: raise ValueError('Input GFF changed after annotation')
    if not gff.exists(): return dict(status='not_a_Prodigal_input', sensitive_sites=[]), bindings
    bindings[str(gff.resolve())] = sha(gff); tables, counts = prodigal_tables(gff)
    if not tables: return dict(status='GFF_has_no_Prodigal_model_header', sensitive_sites=[]), bindings
    requested = manifest.get('translation_table'); raw = {}
    if (annotation/'prodigal.gff').exists() and len(manifest.get('input_files', [])) == 1:
        source = manifest['input_files'][0]; path = Path(source['path'])
        if path.exists():
            if sha(path) != source['sha256']: raise ValueError('Raw genomic input changed after annotation')
            raw = {record.id: record.seq for record in fasta_records(path)}
            bindings[str(path)] = sha(path)
    byid = {g['id']: g for g in genes}; sensitive = []
    for site in sites:
        if site.get('residue_role') != 'direct_metal_ligand': continue
        if site.get('status') != 'substituted' or site.get('observed_residue') != 'W' or 'C' not in site.get('expected_residue', ''): continue
        pid = site['protein_id']; gene = byid.get(pid)
        if not gene or tables.get(gene['contig']) != 4 or gene['contig'] not in raw: continue
        pos = site.get('query_position', '')
        if not pos or not str(pos).isdigit(): continue
        dna = raw[gene['contig']]
        if (int(gene.get('contig_length') or len(dna)) != len(dna)
                or int(gene['strand']) not in (-1, 1)
                or not valid_cds_interval(dict(gene, contig_length=len(dna)), allow_unknown_length=False)):
            continue
        dna = Seq(circular_slice(str(dna), int(gene['start']), int(gene['end'])))
        if int(gene['strand']) == -1: dna = dna.reverse_complement()
        codon = str(dna[3*(int(pos)-1):3*int(pos)]).upper()
        if codon != 'TGA': continue
        sensitive.append(dict(protein_id=pid, query_position=int(pos), reference_position=site['reference_position'],
            observed_called_residue='W', expected_template_residue=site['expected_residue'], genomic_codon='TGA',
            effective_Prodigal_table=4, legacy_requested_table=requested, status='translation_model_dependent_ligand_observation',
            interpretation='W was produced by a table4 gene model. This does not establish a genomic Cys/Sec-to-Trp mutation. DNA alone does not establish Sec; no U is assigned.'))
    known = {cid: code for cid, code in tables.items() if code is not None}
    return dict(status='effective_Prodigal_models_recorded' if len(known) == len(tables) else 'Prodigal_model_headers_incomplete', requested_table_in_legacy_manifest=requested,
        effective_tables_by_contig=tables, coding_features_by_contig=counts,
        contigs_without_effective_model=[cid for cid, code in tables.items() if code is None],
        requested_effective_disagreement=any(v != requested for v in known.values()) if requested is not None and known else None,
        raw_DNA_available_for_codon_checks=bool(raw), sensitive_sites=sensitive,
        scope='Effective codes are parsed from actual Prodigal Model Data headers, including supplied GFF. The legacy table field is a requested value, not proof that meta mode used it. No automatic Sec inference or confidence promotion.'), bindings
