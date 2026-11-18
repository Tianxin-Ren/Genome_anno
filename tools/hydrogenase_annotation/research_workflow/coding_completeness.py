"""Fail-closed, source-bound evidence for protein coding completeness.

This module deliberately proves only that an existing protein is completely
encoded by a bound nucleotide feature.  It never edits a partial flag, makes
an ORF, or infers catalytic function.
"""
import hashlib
import json
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from circular_coordinates import circular_slice
from input_translation import cds_coordinates, normalized_translation, verified_origin_translation
from nucleotide_input import fasta_records

_STATUSES = {"complete_under_bound_translation_model", "incomplete", "unknown", "unsupported"}
RULE_VERSION = "R116_coding_completeness_v2"
_SEC_MATCH = "explicit_genbank_sec"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _seqsha(seq):
    return hashlib.sha256(str(seq).encode()).hexdigest()


def _record_digest(row):
    payload = {k: v for k, v in row.items() if k != "record_digest"}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def supports(evidence, protein_id, expected_protein_sha256=None, feature=None, gene=None):
    """Strict consumer predicate for a validated coding-completeness record."""
    if not isinstance(evidence, dict) or evidence.get("rule_version") != RULE_VERSION:
        return False
    if evidence.get("record_digest") != _record_digest(evidence):
        return False
    if evidence.get("protein_id") != protein_id or evidence.get("status") != "complete_under_bound_translation_model":
        return False
    if evidence.get("coding_completeness_scope") != "protein_coding_completeness_only":
        return False
    if expected_protein_sha256 is not None and evidence.get("protein_sha256") != expected_protein_sha256:
        return False
    if not evidence.get("source_artifact", {}).get("path") or not evidence.get("source_artifact", {}).get("sha256"):
        return False
    if not re.fullmatch(r"[0-9a-f]{64}", str(evidence.get("protein_sha256", ""))):
        return False
    if evidence.get("translated_protein_sha256") != evidence.get("protein_sha256"):
        return False
    mode = evidence.get("translation_match_mode")
    if mode not in {"normalized_initiator_M", "raw_legal_alternate_initiator", _SEC_MATCH}:
        return False
    if not re.fullmatch(r"[0-9a-f]{64}", str(evidence.get("raw_translated_protein_sha256", ""))):
        return False
    if not re.fullmatch(r"[0-9a-f]{64}", str(evidence.get("dna_sha256", ""))):
        return False
    if not re.fullmatch(r"[0-9a-f]{64}", str(evidence.get("translated_protein_sha256", ""))):
        return False
    table_id = evidence.get("translation_table")
    if not isinstance(table_id, int) or table_id not in CodonTable.unambiguous_dna_by_id:
        return False
    table = CodonTable.unambiguous_dna_by_id[table_id]
    if evidence.get("start_codon") not in table.start_codons:
        return False
    if mode == "raw_legal_alternate_initiator":
        if evidence.get("source_partial") is not True or evidence.get("start_codon_normalized") is not False:
            return False
        if evidence.get("translated_protein_sha256") != evidence.get("raw_translated_protein_sha256"):
            return False
    elif mode == "normalized_initiator_M":
        if evidence.get("translated_protein_sha256") != evidence.get("normalized_translated_protein_sha256"):
            return False
        if evidence.get("start_codon_normalized") is not (evidence.get("start_codon") not in {"ATG"}):
            return False
    try:
        start, end, strand = int(evidence.get("start")), int(evidence.get("end")), int(evidence.get("strand"))
    except (TypeError, ValueError):
        return False
    contig_length = int(evidence.get("contig_length", 0) or 0)
    crossing = end > contig_length if contig_length else False
    if start < 1 or end < start or strand not in (-1, 1) or int(evidence.get("dna_length", 0)) != end - start + 1:
        return False
    if crossing:
        if evidence.get("topology") != "circular" or start > contig_length or end > start + contig_length - 1:
            return False
    elif evidence.get("topology") == "circular" and contig_length < end:
        return False
    if evidence.get("terminal_stop_present") is not True or evidence.get("internal_stop_positions") != []:
        return False
    if int(evidence.get("protein_length", 0)) <= 0 or int(evidence.get("translated_length", 0)) != int(evidence.get("protein_length", 0)):
        return False
    artifact = evidence.get("source_artifact", {})
    if not re.fullmatch(r"[0-9a-f]{64}", str(artifact.get("sha256", ""))) or artifact.get("sha256") != evidence.get("input_dna_sha256"):
        return False
    artifact_path = Path(str(artifact.get("path", "")))
    if not artifact_path.is_absolute() or not artifact_path.exists() or _sha(artifact_path) != artifact.get("sha256"):
        return False
    genbank_source = artifact_path.suffix.lower() in {".gb", ".gbk", ".gbff", ".genbank"}
    if mode == _SEC_MATCH and (not genbank_source or evidence.get("source_partial") is not False or
                               gene is not None and "partial" in gene and str(gene["partial"]).lower() != "false"):
        return False
    if crossing or genbank_source:
        if genbank_source:
            try:
                dna, _features, records = _read_genbank(artifact_path)
                record = records.get(evidence.get("contig"))
                owner = dict(evidence)
                if gene is not None:
                    owner.update(gene)
                source_feature = _genbank_feature(record, owner) if record is not None else None
                source_topology = str(record.annotations.get("topology", "")).lower() if record is not None else ""
                source_start, source_end, source_strand, source_partial = cds_coordinates(record, source_feature, source_topology)
                source_table, = source_feature.qualifiers["transl_table"]
                source_protein, = source_feature.qualifiers["translation"]
            except (AttributeError, IndexError, KeyError, OSError, TypeError, ValueError):
                return False
            if (evidence.get("translation_table_source") != "genbank_CDS_qualifier" or
                    not _genbank_topology_allowed(source_topology, mode == _SEC_MATCH) or
                    evidence.get("topology") != source_topology or source_partial or
                    (source_start, source_end, source_strand) != (start, end, strand) or
                    source_table != str(table_id) or
                    bool(source_feature.qualifiers.get("transl_except")) != (mode == _SEC_MATCH)):
                return False
            if mode != _SEC_MATCH and _seqsha(source_protein.rstrip("*").upper()) != evidence["protein_sha256"]:
                return False
        else:
            dna = _read_fasta(artifact_path)
        contig = evidence.get("contig")
        if contig not in dna or len(dna[contig]) != contig_length:
            return False
        nt = circular_slice(dna[contig], start, end) if crossing else dna[contig][start - 1:end]
        if strand == -1:
            nt = str(Seq(nt).reverse_complement())
        if _seqsha(nt) != evidence.get("dna_sha256"):
            return False
        if mode == _SEC_MATCH:
            try:
                model = _source_sec_model(record, source_feature, nt, table_id)
            except (KeyError, ValueError, CodonTable.TranslationError):
                return False
            if any(evidence.get(key) != value for key, value in model.items()):
                return False
    for obj in (feature, gene):
        if obj is None:
            continue
        for key in ("contig", "start", "end", "strand"):
            if key in obj and str(evidence.get(key)) != str(obj.get(key)):
                return False
        if crossing or genbank_source:
            # Confirmation features can carry only the four basic coordinates.
            topology = evidence.get("topology")
            if obj is gene and (str(obj.get("topology", "")).lower() != topology or
                                int(obj.get("contig_length", 0) or 0) != contig_length):
                return False
            if obj is feature and (
                    ("topology" in obj and str(obj.get("topology", "")).lower() != topology) or
                    ("contig_length" in obj and int(obj.get("contig_length", 0) or 0) != contig_length)):
                return False
    return True


def supports_translated_region(evidence, protein_id, query_sequence, feature=None, gene=None):
    """Validate a complete translated region whose terminal stop is absent.

    This is deliberately separate from ``supports``: it never changes the
    source partial/CDS status and can only return the narrowly observed
    ``terminal_stop_missing`` case.
    """
    if not isinstance(evidence, dict) or evidence.get("rule_version") != RULE_VERSION:
        return False
    if evidence.get("record_digest") != _record_digest(evidence) or evidence.get("protein_id") != protein_id:
        return False
    if evidence.get("coding_completeness_scope") != "protein_coding_completeness_only":
        return False
    if evidence.get("status") != "incomplete" or evidence.get("reasons") != ["terminal_stop_missing"]:
        return False
    artifact = evidence.get("source_artifact", {})
    path = Path(str(artifact.get("path", "")))
    if (not path.is_absolute() or not path.is_file() or
            artifact.get("sha256") != evidence.get("input_dna_sha256") or
            _sha(path) != artifact.get("sha256")):
        return False
    try:
        start, end, strand, table_id = (int(evidence[k]) for k in ("start", "end", "strand", "translation_table"))
    except (KeyError, TypeError, ValueError):
        return False
    contig_length = int(evidence.get("contig_length", 0) or 0)
    crossing = end > contig_length if contig_length else False
    if strand not in (-1, 1) or start < 1 or end < start or end - start + 1 != evidence.get("dna_length"):
        return False
    if crossing and (evidence.get("topology") != "circular" or start > contig_length or end > start + contig_length - 1):
        return False
    if not isinstance(table_id, int) or table_id not in CodonTable.unambiguous_dna_by_id:
        return False
    table = CodonTable.unambiguous_dna_by_id[table_id]
    if evidence.get("start_codon") not in table.start_codons or evidence.get("terminal_stop_present") is not False:
        return False
    if (evidence.get("internal_stop_positions") != [] or evidence.get("stop_positions") != [] or
            not query_sequence or evidence.get("protein_length") != len(query_sequence)):
        return False
    if _seqsha(query_sequence) != evidence.get("protein_sha256"):
        return False
    for obj in (feature, gene):
        if obj is not None and any(key in obj and str(evidence.get(key)) != str(obj.get(key))
                                   for key in ("contig", "start", "end", "strand")):
            return False
        if obj is not None and crossing:
            if obj is gene and (str(obj.get("topology", "")).lower() != "circular" or
                                int(obj.get("contig_length", 0) or 0) != contig_length):
                return False
            if obj is feature and (
                    ("topology" in obj and str(obj.get("topology", "")).lower() != "circular") or
                    ("contig_length" in obj and int(obj.get("contig_length", 0) or 0) != contig_length)):
                return False
    dna = _read_fasta(path)
    contig = evidence.get("contig")
    if contig not in dna or (crossing and len(dna[contig]) != contig_length) or (not crossing and end > len(dna[contig])):
        return False
    nt = circular_slice(dna[contig], start, end) if crossing else dna[contig][start - 1:end]
    if strand == -1:
        nt = str(Seq(nt).reverse_complement())
    if len(nt) % 3 or _seqsha(nt) != evidence.get("dna_sha256"):
        return False
    codons = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    if (codons[0] != evidence.get("start_codon") or codons[0] not in table.start_codons or
            codons[-1] in table.stop_codons):
        return False
    if any(c in table.stop_codons for c in codons[:-1]):
        return False
    aas = [_codon_aa(c, table) for c in codons]
    if any(a is None or a == "*" for a in aas):
        return False
    raw = "".join(aas)
    normalized = "M" + raw[1:] if raw else raw
    if normalized != query_sequence and (raw != query_sequence or evidence.get("source_partial") is not True):
        return False
    return True


def _read_fasta(path):
    records = list(fasta_records(path))
    ids = [r.id for r in records]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate FASTA identifiers")
    return {r.id: str(r.seq).upper() for r in records}


def _read_genbank(path):
    """Read GenBank source records and retain CDS features for exact binding."""
    with Path(path).open() as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records or len({r.id for r in records}) != len(records):
        raise ValueError("duplicate or missing GenBank record identifiers")
    return ({r.id: str(r.seq).upper() for r in records},
            {r.id: [f for f in r.features if f.type == "CDS"] for r in records},
            {r.id: r for r in records})


def _genbank_feature(record, gene):
    """Find one translated CDS owned by an annotation row, never by coordinates alone."""
    wanted = {str(gene.get("id", "")), str(gene.get("protein_id", "")),
              str(gene.get("protein_accession", ""))}
    wanted.update(x for x in str(gene.get("aliases", "")).split(";") if x)
    hits = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        values = set()
        for key in ("locus_tag", "protein_id", "gene", "old_locus_tag"):
            values.update(str(x) for x in feature.qualifiers.get(key, []))
        if wanted & values:
            hits.append(feature)
    if len(hits) != 1:
        raise ValueError("GenBank CDS ownership is missing or ambiguous")
    feature = hits[0]
    if (not feature.qualifiers.get("translation") or "pseudo" in feature.qualifiers
            or "pseudogene" in feature.qualifiers):
        raise ValueError("GenBank CDS translation qualifier is missing")
    if int(feature.qualifiers.get("codon_start", ["1"])[0]) != 1:
        raise ValueError("GenBank codon_start is outside coding-completeness scope")
    return feature


def _genbank_topology_allowed(topology, explicit_sec):
    # The ordinary GenBank proof historically requires a circular source.
    return topology == "circular" or (explicit_sec and topology == "linear")


def _genes(path):
    rows = []
    with Path(path).open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) != len(header):
                raise ValueError("genes.tsv column mismatch")
            rows.append(dict(zip(header, values)))
    ids = [r.get("id", "") for r in rows]
    if not all(ids) or len(ids) != len(set(ids)):
        raise ValueError("duplicate or missing gene identifiers")
    return {r["id"]: r for r in rows}


def _effective_tables(gff):
    sys.path.insert(0, str(ROOT / "scripts"))
    from input_translation import prodigal_models
    return prodigal_models(gff)


def _source(manifest):
    files = manifest.get("input_files", [])
    for item in files:
        p = Path(item.get("path", ""))
        if p.suffix.lower() in {".fna", ".fa", ".fasta", ".fas", ".ffn", ".gb", ".gbk", ".gbff", ".genbank"} or p.name.lower().endswith(".fna.gz"):
            return p, item
    return None, None


def _codon_aa(codon, table):
    """Translate a codon only when every IUPAC expansion has one AA."""
    alphabet = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
                "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
                "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
    codon = codon.upper().replace("U", "T")
    if any(c not in alphabet for c in codon):
        return None
    aas = set()
    for a in alphabet[codon[0]]:
        for b in alphabet[codon[1]]:
            for c in alphabet[codon[2]]:
                aas.add(table.forward_table.get(a + b + c, "*" if a + b + c in table.stop_codons else "X"))
    return aas.pop() if len(aas) == 1 else None


def _source_sec_model(record, feature, nt, table_id):
    """Reproduce the exact annotated Sec CDS without extending its bounds."""
    normalized, sites = normalized_translation(record, feature, include_unchanged_sec=True)
    if not sites or not verified_origin_translation(record, feature, normalized):
        raise ValueError("CDS lacks complete, explicit Sec translation evidence")
    coordinates = list(feature.location)
    for site in sites:
        offset = 3 * (site["protein_position"] - 1)
        if site["genomic_positions_1based"] != [p + 1 for p in coordinates[offset:offset + 3]]:
            raise ValueError("Sec qualifier is not in CDS translation order")
    if nt != str(feature.extract(record.seq)).upper():
        raise ValueError("Sec CDS does not match the bound genomic interval")
    table = CodonTable.unambiguous_dna_by_id[table_id]
    codons = [nt[i:i + 3] for i in range(0, len(nt), 3)]
    raw_aas = [_codon_aa(codon, table) for codon in codons[:-1]]
    if any(aa is None or aa == "X" for aa in raw_aas):
        raise ValueError("Sec CDS contains ambiguous translation")
    raw = "".join(raw_aas)
    normalized_sha = _seqsha(normalized)
    return dict(translation_match_mode=_SEC_MATCH, source_sec_sites=sites,
                start_codon=codons[0], start_codon_normalized=raw[0] != "M",
                raw_stop_positions=[i + 1 for i, codon in enumerate(codons) if codon in table.stop_codons],
                stop_positions=[len(codons)], internal_stop_positions=[], terminal_stop_present=True,
                translated_length=len(normalized), raw_translated_protein_sha256=_seqsha(raw),
                normalized_translated_protein_sha256=normalized_sha, translated_protein_sha256=normalized_sha)


def _unknown(pid, gene, partial, reasons, **extra):
    row = dict(protein_id=pid, source_partial=partial, status="unknown", reasons=reasons,
               coding_completeness_scope="protein_coding_completeness_only")
    row.update(extra)
    return row


def _row(pid, gene, partial):
    return dict(protein_id=pid, source_partial=partial, status="unknown", reasons=[],
                coding_completeness_scope="protein_coding_completeness_only",
                contig=gene.get("contig"), start=gene.get("start"), end=gene.get("end"),
                strand=gene.get("strand"))


def prepare(annotation, output, candidate_ids=None, root=ROOT):
    """Create a new coding-completeness evidence directory and return rows/hashes."""
    annotation, output, root = Path(annotation), Path(output), Path(root)
    if output.exists():
        raise FileExistsError("output directory already exists")
    output.mkdir(parents=True)
    manifest_out = output / "coding_completeness_manifest.json"
    def failed(message):
        manifest_out.write_text(json.dumps({"schema": "coding_completeness_manifest_v1",
                                             "status": "failed", "error": str(message)}, indent=2))
    try:
        run_path = annotation / "run_manifest.json"
        run_bytes = run_path.read_bytes()
        run_digest = hashlib.sha256(run_bytes).hexdigest()
        run = json.loads(run_bytes)
    except Exception as exc:
        failed("invalid run manifest: " + str(exc))
        raise
    if run.get("status") not in ("complete", "completed"):
        failed("annotation run is not complete")
        raise ValueError("annotation run is not complete")
    required = [annotation / "genes.tsv", annotation / "proteins.faa"]
    if not all(p.exists() for p in required):
        failed("missing annotation genes.tsv or proteins.faa")
        raise FileNotFoundError("missing annotation genes.tsv or proteins.faa")
    # Bind every artifact before reading any annotation, DNA, or translation table.
    bound = {str(p.resolve()): _sha(p) for p in required if p.exists()}
    bound[str(run_path.resolve())] = run_digest
    missing_inputs = []
    for item in run.get("input_files", []):
        p = Path(item.get("path", ""))
        if not p.exists():
            missing_inputs.append(str(p.resolve()))
            continue
        got = _sha(p)
        if not item.get("sha256"):
            failed("declared input lacks sha256: " + str(p))
            raise ValueError("declared input lacks sha256: " + str(p))
        if item["sha256"] != got:
            failed("declared input hash mismatch: " + str(p))
            raise ValueError("declared input hash mismatch: " + str(p))
        bound[str(p.resolve())] = got
    gff = annotation / "prodigal.gff"
    if gff.exists():
        bound[str(gff.resolve())] = _sha(gff)
    translation_helper = root / "scripts" / "input_translation.py"
    if not translation_helper.exists():
        translation_helper = ROOT / "scripts" / "input_translation.py"
    if not translation_helper.exists():
        failed("translation helper is missing")
        raise FileNotFoundError(str(translation_helper))
    bound[str(translation_helper.resolve())] = _sha(translation_helper)
    bound[str(Path(__file__).resolve())] = _sha(Path(__file__))
    try:
        genes = _genes(required[0])
        proteins = _read_fasta(required[1])
    except Exception as exc:
        failed("invalid annotation records: " + str(exc))
        raise
    declared_proteins = run.get("normalized_proteins_sha256")
    if declared_proteins and declared_proteins != _sha(required[1]):
        failed("normalized proteins hash mismatch")
        raise ValueError("normalized proteins hash mismatch")
    declared_outputs = run.get("output_sha256", {})
    if isinstance(declared_outputs, dict):
        for name, digest in declared_outputs.items():
            p = annotation / name if not Path(name).is_absolute() else Path(name)
            if p.exists() and p in required and _sha(p) != digest:
                failed("declared output hash mismatch: " + str(p))
                raise ValueError("declared output hash mismatch: " + str(p))
    if set(genes) != set(proteins):
        failed("gene/protein identifiers are not one-to-one")
        raise ValueError("gene/protein identifiers are not one-to-one")
    selected_input = list(candidate_ids) if candidate_ids is not None else list(genes)
    if len(selected_input) != len(set(selected_input)):
        failed("candidate IDs contain duplicates")
        raise ValueError("candidate IDs contain duplicates")
    selected = sorted(selected_input)
    missing = sorted(set(selected) - set(genes))
    if missing:
        failed("candidate IDs absent from annotation: " + ",".join(missing))
        raise KeyError("candidate IDs absent from annotation: " + ",".join(missing))
    tables, counts = ({}, {})
    if gff.exists():
        try:
            tables, counts = _effective_tables(gff)
        except Exception as exc:
            failed("effective translation model invalid: " + str(exc))
            raise
    source, source_item = _source(run)
    dna = None
    genbank_records = {}
    genbank_features = {}
    source_is_genbank = bool(source and source.suffix.lower() in {".gb", ".gbk", ".gbff", ".genbank"})
    source_reason = None
    if source is None:
        source_reason = "source_nucleotide_unavailable"
    elif not source.exists():
        source_reason = "source_nucleotide_missing"
    elif source_is_genbank:
        try:
            dna, genbank_features, genbank_records = _read_genbank(source)
        except (OSError, ValueError, KeyError, TypeError):
            source_reason = "genbank_source_validation_failed"
    elif not (source.suffix.lower() in {".fna", ".fa", ".fasta", ".fas", ".ffn"} or source.name.lower().endswith(".fna.gz")):
        source_reason = "input_mode_not_nucleotide"
    else:
        dna = _read_fasta(source)
    rows = []
    for pid in selected:
        gene = genes[pid]
        raw_partial = gene.get("partial")
        partial = str(raw_partial).lower() == "true" if str(raw_partial).lower() in ("true", "false") else None
        r = _row(pid, gene, partial)
        r["raw_source_partial"] = raw_partial
        r["rule_version"] = RULE_VERSION
        r["protein_sha256"] = _seqsha(proteins[pid]); r["protein_length"] = len(proteins[pid])
        if source_reason:
            r["status"] = "unknown"; r["reasons"] = [source_reason]; rows.append(r); continue
        if gene.get("contig") not in dna:
            r["reasons"] = ["contig_not_in_source_fasta"]; rows.append(r); continue
        if (not source_is_genbank and (not gff.exists() or gene.get("contig") not in tables or
                                           not isinstance(tables[gene.get("contig")], int))):
            r["reasons"] = ["effective_translation_model_unavailable"]; rows.append(r); continue
        try:
            contig = gene["contig"]
            contig_length = len(dna[contig])
            source_feature = None
            if source_is_genbank:
                source_feature = _genbank_feature(genbank_records[contig], gene)
                annotated_translation, = source_feature.qualifiers["translation"]
                explicit_sec = bool(source_feature.qualifiers.get("transl_except"))
                if not explicit_sec and annotated_translation.rstrip("*").upper() != proteins[pid]:
                    raise ValueError("GenBank CDS translation conflicts with proteins.faa")
                topology = str(genbank_records[contig].annotations.get("topology", "")).lower()
                start, end, strand, source_partial = cds_coordinates(genbank_records[contig], source_feature, topology)
                called = tuple(int(gene[key]) for key in ("start", "end", "strand", "contig_length"))
                if (not _genbank_topology_allowed(topology, explicit_sec) or gene.get("topology") != topology or
                        called != (start, end, strand, contig_length) or
                        source_partial or partial is not False):
                    raise ValueError("GenBank CDS is outside the bound annotation coordinates")
            else:
                start, end, strand = int(gene["start"]), int(gene["end"]), int(gene["strand"])
                topology = str(gene.get("topology", "")).lower()
            crossing = end > contig_length
            if strand not in (-1, 1) or start < 1 or end < start:
                raise ValueError("invalid_coordinates")
            if crossing:
                if topology != "circular" or (not source_is_genbank and int(gene.get("contig_length", 0) or 0) != contig_length) or end > start + contig_length - 1:
                    raise ValueError("invalid_circular_coordinates")
                nt = circular_slice(dna[contig], start, end)
            else:
                if end > contig_length or (topology == "circular" and not source_is_genbank and int(gene.get("contig_length", 0) or 0) != contig_length):
                    raise ValueError("invalid_coordinates")
                nt = dna[contig][start - 1:end]
            if strand == -1:
                nt = str(Seq(nt).reverse_complement())
            if source_is_genbank:
                source_table, = source_feature.qualifiers["transl_table"]
                table_id = int(source_table)
            else:
                table_id = tables[gene["contig"]]
            r.update(start=start, end=end, strand=strand, translation_table=table_id,
                     translation_table_source="genbank_CDS_qualifier" if source_is_genbank else "prodigal_model_header",
                     dna_length=len(nt), dna_sha256=_seqsha(nt),
                     input_dna_sha256=bound.get(str(source.resolve())),
                     source_artifact={"path": str(source.resolve()), "sha256": bound.get(str(source.resolve()))})
            if crossing or source_is_genbank:
                r["topology"] = topology
                r["contig_length"] = contig_length
            if source_is_genbank and explicit_sec:
                model = _source_sec_model(genbank_records[contig], source_feature, nt, table_id)
                if model["translated_protein_sha256"] != r["protein_sha256"]:
                    raise ValueError("Explicit Sec CDS translation conflicts with proteins.faa")
                r.update(model, status="complete_under_bound_translation_model")
                rows.append(r)
                continue
            if len(nt) % 3:
                r["reasons"] = ["frameshift_cds_length_not_multiple_of_3"]; r["status"] = "incomplete"; rows.append(r); continue
            table = CodonTable.unambiguous_dna_by_id[table_id]
            codons = [nt[i:i + 3] for i in range(0, len(nt), 3)]
            r["start_codon"] = codons[0]; r["stop_positions"] = [i + 1 for i, c in enumerate(codons) if c in table.stop_codons]
            r["terminal_stop_present"] = bool(codons and codons[-1] in table.stop_codons)
            r["internal_stop_positions"] = [i + 1 for i, c in enumerate(codons[:-1]) if c in table.stop_codons]
            if codons[0] not in table.start_codons:
                r["reasons"] = ["canonical_start_missing"]; r["status"] = "incomplete"; rows.append(r); continue
            if r["internal_stop_positions"]:
                r["reasons"] = ["internal_stop"]; r["status"] = "incomplete"; rows.append(r); continue
            if not r["terminal_stop_present"]:
                r["reasons"] = ["terminal_stop_missing"]; r["status"] = "incomplete"; rows.append(r); continue
            aas = [_codon_aa(c, table) for c in codons]
            if any(a is None for a in aas):
                r["reasons"] = ["ambiguous_codon_not_uniquely_translatable"]; rows.append(r); continue
            raw_translated = "".join(aas[:-1])
            normalized = "M" + raw_translated[1:]
            r["translated_length"] = len(raw_translated)
            r["raw_translated_protein_sha256"] = _seqsha(raw_translated)
            r["normalized_translated_protein_sha256"] = _seqsha(normalized)
            if normalized == proteins[pid]:
                r["translation_match_mode"] = "normalized_initiator_M"
                r["start_codon_normalized"] = aas[0] != "M"
                r["translated_protein_sha256"] = _seqsha(normalized)
            elif partial is True and codons[0] in table.start_codons and raw_translated == proteins[pid]:
                r["translation_match_mode"] = "raw_legal_alternate_initiator"
                r["start_codon_normalized"] = False
                r["translated_protein_sha256"] = _seqsha(raw_translated)
            else:
                r["reasons"] = ["protein_translation_mismatch"]; r["status"] = "incomplete"; rows.append(r); continue
            r["status"] = "complete_under_bound_translation_model"
            r["reasons"] = ["source_partial_retained; coding completeness is separate from source partial flag"] if partial else []
        except (KeyError, ValueError, IndexError, CodonTable.TranslationError):
            r["status"] = "unknown"; r["reasons"] = ["coding_feature_validation_failed"]
        r["record_digest"] = _record_digest(r)
        rows.append(r)
    for r in rows:
        if "record_digest" not in r:
            r["record_digest"] = _record_digest(r)
    before = dict(bound)
    after = {p: _sha(p) for p in bound if Path(p).exists()}
    if before != after:
        manifest_out.write_text(json.dumps({"schema": "coding_completeness_v1", "status": "failed", "error": "bound input changed during scan", "input_hashes_before": before, "input_hashes_after": after}, indent=2))
        raise RuntimeError("bound input changed during scan")
    data = {"schema": "coding_completeness_v1", "scope": "protein coding completeness only", "rows": rows}
    data_path = output / "coding_completeness.json"; data_path.write_text(json.dumps(data, indent=2, sort_keys=True))
    manifest = {"schema": "coding_completeness_manifest_v1", "status": "complete", "rule_version": RULE_VERSION,
                "source_annotation": str(annotation.resolve()), "input_hashes": bound,
                "required_files": sorted(bound),
                "missing_inputs": sorted(missing_inputs),
                "output_sha256": _sha(data_path), "effective_tables_by_contig": tables,
                "effective_table_source": str(gff.resolve()) if gff.exists() else None,
                "script_sha256": _sha(Path(__file__)),
                "input_translation_helper_sha256": bound.get(str(translation_helper.resolve()))}
    manifest_out.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return {r["protein_id"]: r for r in rows}, bound


def load_evidence(output, annotation=None, root=ROOT):
    output = Path(output)
    manifest = json.loads((output / "coding_completeness_manifest.json").read_text())
    if manifest.get("status") != "complete":
        raise ValueError("coding completeness manifest is not complete")
    if manifest.get("rule_version") != RULE_VERSION or annotation is None:
        raise ValueError("missing or incompatible evidence owner")
    expected_owner = str(Path(annotation).resolve())
    if manifest.get("source_annotation") != expected_owner:
        raise ValueError("evidence owner mismatch")
    ann = Path(annotation).resolve()
    run_path = ann / "run_manifest.json"
    try:
        run_bytes = run_path.read_bytes(); run = json.loads(run_bytes)
    except Exception as exc:
        raise ValueError("current run manifest unavailable") from exc
    if run.get("status") not in ("complete", "completed"):
        raise ValueError("current annotation run is not complete")
    helper = Path(root).resolve() / "scripts" / "input_translation.py"
    if not helper.exists():
        helper = ROOT / "scripts" / "input_translation.py"
    expected_paths = {run_path, ann / "genes.tsv", ann / "proteins.faa", Path(__file__).resolve(), helper}
    gff_now = ann / "prodigal.gff"
    if gff_now.exists():
        expected_paths.add(gff_now)
    missing_expected = set()
    for item in run.get("input_files", []):
        p = Path(item.get("path", "")).resolve()
        expected_paths.add(p)
        if not p.exists():
            missing_expected.add(str(p))
    actual_existing = {str(p.resolve()) for p in expected_paths if p.exists()}
    declared = manifest.get("input_hashes")
    required = manifest.get("required_files")
    if not isinstance(declared, dict) or not declared or not isinstance(required, list) or set(required) != set(declared):
        raise ValueError("incomplete bound input envelope")
    if set(declared) != actual_existing or sorted(manifest.get("missing_inputs", [])) != sorted(missing_expected):
        raise ValueError("canonical bound input set mismatch")
    for path, digest in declared.items():
        p = Path(path)
        if not p.is_absolute() or not digest or not p.exists() or _sha(p) != digest:
            raise ValueError("bound input missing or changed: " + path)
    script = Path(__file__).resolve()
    if manifest.get("script_sha256") != _sha(script) or manifest.get("input_translation_helper_sha256") != _sha(helper):
        raise ValueError("bound helper changed")
    data_path = output / "coding_completeness.json"
    if _sha(data_path) != manifest.get("output_sha256"):
        raise ValueError("coding completeness output hash mismatch")
    data = json.loads(data_path.read_text()); rows = data.get("rows", [])
    genes = _genes(ann / "genes.tsv")
    proteins = _read_fasta(ann / "proteins.faa")
    if set(genes) != set(proteins):
        raise ValueError("current annotation IDs are inconsistent")
    if len({r.get("protein_id") for r in rows}) != len(rows) or any(r.get("status") not in _STATUSES for r in rows):
        raise ValueError("invalid coding completeness rows")
    for r in rows:
        pid = r.get("protein_id")
        if pid not in genes or r.get("protein_sha256") != _seqsha(proteins[pid]):
            raise ValueError("current protein does not match evidence")
        if r.get("record_digest") != _record_digest(r):
            raise ValueError("record digest mismatch")
    if manifest.get("missing_inputs") and any(r.get("status") == "complete_under_bound_translation_model" for r in rows):
        raise ValueError("complete proof has missing bound inputs")
    return manifest, {r["protein_id"]: r for r in rows}
