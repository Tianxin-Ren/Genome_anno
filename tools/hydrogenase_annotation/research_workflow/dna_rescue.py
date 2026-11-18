"""DNA-level hydrogenase recall audit.

This stage deliberately produces unresolved DNA evidence.  It never repairs a
protein FASTA, invents a CDS, or upgrades identity/function.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import argparse
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "scripts"))
from input_translation import prodigal_models
from discovery_evidence import load_domains, load_targets, search_discovery
from circular_coordinates import circular_slice, interval_segments, valid_cds_interval
from nucleotide_input import fasta_records, materialize_fasta
from tool_paths import resolve_tool
KINDS = {"catalytic", "unresolved", "decoy"}
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
OUTFMT = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qseq sseq sframe"
CANDIDATE_GROUPING = "native_span_frame_dna_and_peptide_v1"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def _candidate_key(row):
    return (row["contig"], int(row["dna_start"]), int(row["dna_end"]), int(row["frame"]),
            row["dna_sha256"], row["translation_sha256"])


def _candidate_id(key, table):
    contig, start, end, frame, dna, _ = key
    # Retain the previous ID of the lowest table supporting the same peptide.
    legacy_key = (contig, start, end, frame, table, dna)
    return "DNA_" + hashlib.sha256("|".join(map(str, legacy_key)).encode()).hexdigest()[:20]


def _span_count(candidates):
    return len({(row["contig"], row["start"], row["end"], row["strand"]) for row in candidates})


def _verify_candidate_groups(candidates, hsps, manifest):
    expected = {}
    for hid, row in hsps.items():
        if not row["represented_by"]:
            expected.setdefault(_candidate_key(row), set()).add(hid)
    for candidate in candidates:
        key = (candidate["contig"], candidate["start"], candidate["end"], candidate["frame"],
               candidate["DNA_sha256"], candidate["translation_sha256"])
        ids = candidate["source_hsp_ids"]
        if not ids or len(set(ids)) != len(ids) or set(ids) != expected.pop(key, None):
            raise ValueError("DNA candidate grouping differs from its unrepresented HSP evidence")
        tables = sorted({int(hsps[hid]["translation_table"]) for hid in ids})
        if (candidate["translation_tables"] != tables
                or candidate["strand"] != (1 if key[3] > 0 else -1)
                or candidate["candidate_id"] != _candidate_id(key, tables[0])):
            raise ValueError("DNA candidate table provenance or identity differs from its HSPs")
    if expected or manifest["distinct_native_span_count"] != _span_count(candidates):
        raise ValueError("DNA candidate grouping or native span count is incomplete")


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def _write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _input_source(annotation, run):
    if run.get("input_mode") == "transcript":
        return None, "input_mode_is_transcript"
    for item in run.get("input_files", []):
        path = Path(item.get("path", ""))
        if path.suffix.lower() in {".fna", ".fa", ".fasta", ".fas", ".gb", ".gbk", ".gbff", ".genbank"} or path.name.lower().endswith(".fna.gz"):
            return path, None
    if not (Path(annotation) / "prodigal.gff").exists():
        return None, "prodigal_gff_missing_for_dna_rescue"
    return None, "input_is_not_nucleotide_fasta"


def _source_format(path):
    name = Path(path).name.lower()
    if name.endswith(".fna.gz") or Path(path).suffix.lower() in {".fna", ".fa", ".fasta", ".fas"}:
        return "fasta"
    if Path(path).suffix.lower() in {".gb", ".gbk", ".gbff", ".genbank"}:
        return "genbank"
    return "other"


def _genbank_subject(source, out):
    with Path(source).open() as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records:
        raise ValueError("GenBank source contains no records")
    ids = [r.id for r in records]
    if len(ids) != len(set(ids)):
        raise ValueError("GenBank source contains duplicate contig IDs")
    tables = set()
    with Path(out).open("w") as handle:
        for record in records:
            seq = str(record.seq).upper()
            if not seq or any(base not in DNA_ALPHABET for base in seq):
                raise ValueError("GenBank source contains invalid nucleotide alphabet")
            for feature in record.features:
                value = feature.qualifiers.get("transl_table", [])
                for item in value:
                    if str(item).isdigit():
                        tables.add(int(item))
            handle.write(f">{record.id}\n{seq}\n")
    return tables


def _load_genes(annotation):
    with (Path(annotation) / "genes.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids = [r.get("id", "") for r in rows]
    if any(not x for x in ids) or len(ids) != len(set(ids)):
        raise ValueError("Annotation genes contain missing or duplicate IDs")
    required = {"id", "contig", "start", "end", "strand"}
    if any(not required <= set(r) for r in rows):
        raise ValueError("Annotation genes lack ownership columns")
    return {r["id"]: r for r in rows}


def _load_proteins(path):
    with Path(path).open() as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len({r.id for r in records}) != len(records):
        raise ValueError("Annotation proteins contain duplicate IDs")
    return {r.id: str(r.seq).rstrip("*") for r in records}


def _reference_query(root):
    meta_path = root / "references/processed/references.json"
    faa_path = root / "references/processed/references.faa"
    refs = [r for r in load(meta_path) if r.get("kind") in KINDS]
    wanted = {r["reference_id"] for r in refs}
    with faa_path.open() as handle:
        seqs = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(handle, "fasta")}
    missing = sorted(wanted - set(seqs))
    if missing:
        raise ValueError("Reference sequence missing from bound FASTA: " + missing[0])
    # BLAST protein alphabets do not encode selenocysteine; retain the
    # normalization explicitly in the manifest rather than silently changing it.
    return refs, {rid: seqs[rid].replace("U", "X") for rid in sorted(wanted)}, meta_path, faa_path


def _prodigal_tables(annotation):
    gff = Path(annotation) / "prodigal.gff"
    if not gff.exists():
        return []
    tables, _counts = prodigal_models(gff)
    return sorted(set(tables.values()))


def _fasta(path):
    records = list(fasta_records(path))
    if not records:
        raise ValueError("DNA input contains no FASTA records")
    contigs = {r.id: str(r.seq).upper() for r in records}
    if len(contigs) != len(records):
        raise ValueError("DNA input contains duplicate FASTA IDs")
    for name, seq in contigs.items():
        if not seq or any(base not in DNA_ALPHABET for base in seq):
            raise ValueError("DNA input contains invalid alphabet: " + name)
    return contigs


def _circular_lengths(genes, contigs):
    """Return only explicitly declared circular contigs with trusted lengths."""
    lengths = {}
    for gene in genes.values():
        if str(gene.get("topology", "")).lower() != "circular":
            continue
        contig = gene.get("contig")
        try:
            declared = int(gene.get("contig_length", 0))
        except (TypeError, ValueError):
            raise ValueError("Circular gene lacks a valid contig length")
        if contig not in contigs or declared != len(contigs[contig]) or declared < 1:
            raise ValueError("Circular contig length differs from native FASTA")
        previous = lengths.setdefault(contig, declared)
        if previous != declared:
            raise ValueError("Conflicting circular contig lengths")
    return lengths


def _least_rotation_offset(sequence):
    """Return the first offset of the lexicographically least rotation in O(n)."""
    doubled = sequence + sequence
    length = len(sequence)
    left, right, matched = 0, 1, 0
    while left < length and right < length and matched < length:
        a, b = doubled[left + matched], doubled[right + matched]
        if a == b:
            matched += 1
            continue
        if a > b:
            left = max(left + matched + 1, right + 1)
        else:
            right = max(right + matched + 1, left + 1)
        matched = 0
    return min(left, right)


def _search_contigs(contigs, offsets):
    return {name: seq[offsets.get(name, 0):] + seq[:offsets.get(name, 0)]
            for name, seq in contigs.items()}


def _seam_fasta(contigs, circular_lengths, out):
    """Write a doubled subject for one-turn seam HSPs; native input is untouched."""
    doubled = {}
    with Path(out).open("w") as handle:
        for contig, length in sorted(circular_lengths.items()):
            sequence = contigs[contig]
            doubled[contig] = sequence + sequence
            handle.write(f">{contig}\n{doubled[contig]}\n")
    return doubled


def _normalize_seam_row(row, native_length):
    """Attach native one-turn coordinates to a HSP from a doubled subject."""
    if row["slen"] != native_length * 2:
        raise ValueError("Circular seam subject length differs from native contig")
    low, high = sorted((row["sstart"], row["send"]))
    span = high - low + 1
    if low < 1 or high > row["slen"] or span < 1:
        raise ValueError("Invalid circular seam HSP coordinates")
    if span > native_length or not (low <= native_length < high):
        return None
    row["seam_sstart"], row["seam_send"], row["seam_slen"] = row["sstart"], row["send"], row["slen"]
    row["native_start"] = ((low - 1) % native_length) + 1
    row["native_end"] = row["native_start"] + span - 1
    row["slen"] = native_length
    row["circular_seam"] = row["native_end"] > native_length
    return row


def _native_hsp_span(row, offset=0):
    length = row["slen"]
    low = row.get("native_start", min(row["sstart"], row["send"]))
    high = row.get("native_end", max(row["sstart"], row["send"]))
    start = (low + offset - 1) % length + 1
    end = start + high - low
    frame = (start - 1) % 3 + 1 if row["sframe"] > 0 else -((length - end) % 3 + 1)
    return start, end, frame


def _parse_float(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("Invalid BLAST " + label) from error
    if not math.isfinite(result):
        raise ValueError("Invalid BLAST " + label)
    return result


def _parse_hsp(line, refs, contigs):
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 15:
        raise ValueError("Invalid tblastn outfmt row")
    names = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qseq sseq sframe".split()
    row = dict(zip(names, fields))
    if row["qseqid"] not in refs or row["sseqid"] not in contigs:
        raise ValueError("BLAST query or subject is not bound to the source")
    integer = ["length", "qstart", "qend", "sstart", "send", "qlen", "slen", "sframe"]
    for key in integer:
        try:
            row[key] = int(row[key])
        except (TypeError, ValueError) as error:
            raise ValueError("Invalid BLAST " + key) from error
    row["pident"] = _parse_float(row["pident"], "pident")
    row["evalue"] = _parse_float(row["evalue"], "evalue")
    row["bitscore"] = _parse_float(row["bitscore"], "bitscore")
    if row["length"] <= 0 or not 0 <= row["pident"] <= 100 or row["evalue"] < 0 or row["bitscore"] < 0:
        raise ValueError("Invalid BLAST score")
    if not 1 <= row["qstart"] <= row["qend"] <= row["qlen"] or row["qlen"] != len(refs[row["qseqid"]]):
        raise ValueError("Invalid BLAST query coordinates or length")
    if row["slen"] != len(contigs[row["sseqid"]]) or not 1 <= abs(row["sframe"]) <= 3:
        raise ValueError("Invalid BLAST subject length or frame")
    if not 1 <= min(row["sstart"], row["send"]) <= max(row["sstart"], row["send"]) <= row["slen"]:
        raise ValueError("Invalid BLAST subject coordinates")
    qseq = row["qseq"].replace("-", "").upper().replace("U", "X")
    sseq = row["sseq"].replace("-", "").upper()
    if qseq != refs[row["qseqid"]][row["qstart"] - 1:row["qend"]]:
        raise ValueError("BLAST query alignment does not match bound reference")
    if len(row["qseq"]) != row["length"] or len(row["sseq"]) != row["length"]:
        raise ValueError("BLAST alignment length differs")
    if len(qseq) != row["qend"] - row["qstart"] + 1:
        raise ValueError("BLAST query alignment coordinates differ")
    if len(sseq) * 3 != abs(row["send"] - row["sstart"]) + 1:
        raise ValueError("BLAST subject coordinates are not whole codons")
    phase = ((min(row["sstart"], row["send"]) - 1) if row["sframe"] > 0
             else (row["slen"] - max(row["sstart"], row["send"])))
    if phase % 3 != abs(row["sframe"]) - 1:
        raise ValueError("BLAST subject frame phase differs")
    row["qseq_ungapped"] = qseq
    row["sseq_ungapped"] = sseq
    return row


def _oriented_translation(dna, start, end, frame, table):
    """Return translated codons covering a BLAST subject span and raw stops."""
    low, high = sorted((start, end))
    if frame > 0:
        if (low - 1) % 3 != frame - 1:
            raise ValueError("BLAST positive frame phase differs")
        nt = dna[low - 1:high]
    else:
        if (len(dna) - high) % 3 != abs(frame) - 1:
            raise ValueError("BLAST negative frame phase differs")
        oriented = str(Seq(dna).reverse_complement())
        nt = oriented[len(dna) - high:len(dna) - low + 1]
    if len(nt) == 0 or len(nt) % 3:
        raise ValueError("BLAST subject span is not a complete codon span")
    aa = str(Seq(nt).translate(table=table, to_stop=False))
    return aa, nt, (low, high)


def _representation(row, genes, proteins, screened, contigs=None, table=None):
    span_start = row.get("native_start", min(row["sstart"], row["send"]))
    span_end = row.get("native_end", max(row["sstart"], row["send"]))
    for gid, gene in genes.items():
        if gene["contig"] != row["sseqid"] or gid not in screened or gid not in proteins:
            continue
        try:
            contig_length = len(contigs[row["sseqid"]]) if contigs is not None else int(gene.get("contig_length", 0) or 0)
            gene_start, gene_end = int(gene["start"]), int(gene["end"])
            crossing = 0 < contig_length < gene_end
            if crossing:
                if not valid_cds_interval(gene, allow_unknown_length=False) or int(gene["contig_length"]) != contig_length:
                    continue
                # Keep the native HSP coordinates, but compare them on the
                # one-turn unwrapped CDS arc (e.g. 950..1050 on a 1000 bp contig).
                interval_segments(gene_start, gene_end, contig_length)
                if span_start >= gene_start:
                    hsp_start, hsp_end = span_start, span_end
                else:
                    hsp_start, hsp_end = span_start + contig_length, span_end + contig_length
                contained = (gene["contig"] == row["sseqid"] and
                             hsp_start >= gene_start and hsp_end <= gene_end and
                             (1 if row["sframe"] > 0 else -1) == int(gene["strand"]))
            else:
                hsp_start, hsp_end = span_start, span_end
                contained = (gene["contig"] == row["sseqid"] and gene["strand"] and
                             gene_start <= span_start and gene_end >= span_end and
                             (1 if row["sframe"] > 0 else -1) == int(gene["strand"]))
        except (KeyError, TypeError, ValueError):
            contained = False
            hsp_start, hsp_end = span_start, span_end
        if not contained:
            continue
        source = proteins[gid]
        span_len = span_end - span_start + 1
        if span_len % 3:
            continue
        if int(gene["strand"]) > 0:
            offset = hsp_start - gene_start
        else:
            offset = gene_end - hsp_end
        if offset < 0 or offset % 3:
            continue
        aa_start = offset // 3
        translated = row.get("translation", "")
        comparable = translated
        normalized = False
        if (aa_start == 0 and source.startswith("M") and translated[:1] and translated[:1] != "M"
                and contigs is not None and table is not None):
            from Bio.Data import CodonTable
            dna = contigs[row["sseqid"]]
            if crossing:
                start_codon = ((gene_start - 1 if int(gene["strand"]) > 0 else gene_end - 3) % contig_length) + 1
                codon = circular_slice(dna, start_codon, start_codon + 2)
                if int(gene["strand"]) < 0:
                    codon = str(Seq(codon).reverse_complement())
            elif int(gene["strand"]) > 0:
                codon = dna[int(gene["start"]) - 1:int(gene["start"]) + 2]
            else:
                codon = str(Seq(dna[int(gene["end"]) - 3:int(gene["end"])]).reverse_complement())
            starts = CodonTable.unambiguous_dna_by_id[int(table)].start_codons
            if codon in starts:
                comparable = "M" + translated[1:]
                normalized = True
        if source[aa_start:aa_start + len(comparable)] == comparable:
            row["represented_start_codon_normalized"] = normalized
            return gid
    return ""


def _run_tblastn(root, query, subject, table, out, threads, max_targets):
    command = [str(resolve_tool(root, "tblastn")), "-query", str(query), "-subject", str(subject),
               "-db_gencode", str(table), "-evalue", "10", "-seg", "no",
               "-max_target_seqs", str(max_targets),
               "-outfmt", "6 " + OUTFMT, "-out", str(out)]
    log = out.with_suffix(".log")
    with log.open("w") as handle:
        result = subprocess.run(command, cwd=root, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("tblastn failed for translation table " + str(table))
    return command, log


def _bound_path_objects(annotation, screen_dir, run_path, screen_path, source, root, raw=False):
    if raw:
        paths = {"raw_nucleotide": source, "dna_rescue_script": Path(__file__),
                 "reference_metadata": root / "references/processed/references.json",
                 "reference_fasta": root / "references/processed/references.faa"}
    else:
        paths = {
            "annotation_run_manifest": run_path,
            "annotation_genes": Path(annotation) / "genes.tsv",
            "annotation_proteins": Path(annotation) / "proteins.faa",
            "screen_manifest": screen_path,
            "screen_table": Path(screen_dir) / "candidate_screen.tsv",
            "dna_rescue_script": Path(__file__),
            "reference_metadata": root / "references/processed/references.json",
            "reference_fasta": root / "references/processed/references.faa",
        }
        gff = Path(annotation) / "prodigal.gff"
        if gff.exists():
            paths["annotation_prodigal_gff"] = gff
        if source is not None and source.exists():
            paths["source_nucleotide"] = source
    binary = resolve_tool(root, "tblastn")
    if binary.exists():
        paths["tblastn_binary"] = binary
    return paths


def _bound_paths(annotation, screen_dir, run_path, screen_path, source, root, raw=False):
    return {name: sha(path) for name, path in _bound_path_objects(annotation, screen_dir, run_path, screen_path, source, root, raw).items() if path is not None and path.exists()}


def _assert_run_reference_bindings(run, root):
    for field, rel in (("references_json_sha256", "references/processed/references.json"),
                       ("reference_fasta_sha256", "references/processed/references.faa")):
        declared = run.get(field)
        if declared and declared != sha(root / rel):
            raise ValueError("Source run reference binding changed: " + field)


def _write_version(root, out):
    binary = resolve_tool(root, "tblastn")
    if not binary.exists():
        raise ValueError("tblastn binary is unavailable")
    version = subprocess.run([str(binary), "-version"], cwd=root, capture_output=True, text=True)
    if version.returncode:
        raise RuntimeError("tblastn version query failed")
    path = out / "tblastn_version.txt"
    path.write_text((version.stdout + version.stderr).strip() + "\n")
    return path


def prepare(annotation=None, screen_dir=None, out=None, threads=2, root=ROOT, genome=None):
    root = Path(root).resolve()
    if out is None:
        raise ValueError("DNA rescue output is required")
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("DNA rescue threads must be positive")
    if genome is not None and (annotation is not None or screen_dir is not None):
        raise ValueError("Raw genome mode cannot receive annotation or screen")
    raw_mode = genome is not None
    annotation = Path(annotation).resolve() if annotation is not None else None
    screen_dir = Path(screen_dir).resolve() if screen_dir is not None else None
    out = Path(out).resolve()
    if raw_mode:
        genome = Path(genome).resolve()
    if out.exists():
        raise FileExistsError("Use a fresh DNA rescue output directory")
    out.mkdir(parents=True)
    manifest_path = out / "dna_rescue_manifest.json"
    base = dict(status="running", source_mode="raw_nucleotide_only" if raw_mode else "annotation_bound",
                annotation=str(annotation) if annotation else None, screen=str(screen_dir) if screen_dir else None,
                profile_confirmation="not_performed",
                annotation_run_manifest_sha256=None, screen_manifest_sha256=None,
                input_sha256={}, output_sha256={}, candidates=0,
                requested_threads=threads, effective_threads=1)
    _write_json(manifest_path, base)
    try:
        run_path = annotation / "run_manifest.json" if annotation else None
        screen_path = screen_dir / "candidate_screen_manifest.json" if screen_dir else None
        if raw_mode:
            if not genome.exists():
                base.update(status="not_evaluated", reason="raw_nucleotide_input_missing", input=str(genome))
                _write_json(manifest_path, base)
                return base
            source = genome
            run, screen = {}, {}
            before_bound = _bound_paths(annotation, screen_dir, None, None, source, root, raw=True)
            if not {"dna_rescue_script", "reference_metadata", "reference_fasta", "tblastn_binary"} <= set(before_bound):
                raise ValueError("Raw DNA rescue binding is incomplete")
            if _source_format(source) == "other":
                base.update(status="not_applicable", reason="raw_input_is_not_nucleotide_fasta")
                base.update(bound_input_sha256=before_bound,
                            bound_input_paths={k: str(v) for k, v in _bound_path_objects(annotation, screen_dir, None, None, source, root, raw=True).items()})
                _write_json(manifest_path, base)
                return base
            source_digest = sha(source)
            if set(before_bound) != {"raw_nucleotide", "dna_rescue_script", "reference_metadata", "reference_fasta", "tblastn_binary"}:
                raise ValueError("Raw DNA rescue binding is incomplete")
            base.update(input_sha256={"source_nucleotide": source_digest}, bound_input_sha256=before_bound,
                        bound_input_paths={k: str(v) for k, v in _bound_path_objects(annotation, screen_dir, None, None, source, root, raw=True).items()})
        elif not run_path.exists() or not screen_path.exists():
            base.update(status="not_evaluated", reason="annotation_or_screen_manifest_missing")
            _write_json(manifest_path, base)
            return base
        else:
            run, screen = load(run_path), load(screen_path)
            base["annotation_run_manifest_sha256"] = sha(run_path)
            base["screen_manifest_sha256"] = sha(screen_path)
            if run.get("status") != "complete":
                raise ValueError("Source annotation is not complete")
            if screen.get("status") != "complete" or Path(screen.get("annotation", "")).resolve() != annotation:
                raise ValueError("Screen manifest ownership or status is invalid")
            _assert_run_reference_bindings(run, root)
            for name, digest in screen.get("output_sha256", {}).items():
                candidate = screen_dir / name
                if not candidate.exists() or sha(candidate) != digest:
                    raise ValueError("Screen output binding changed or missing: " + name)
            for name, digest in run.get("output_sha256", {}).items():
                candidate = annotation / name
                if not candidate.exists() or sha(candidate) != digest:
                    raise ValueError("Annotation output binding changed or missing: " + name)
            before_bound = _bound_paths(annotation, screen_dir, run_path, screen_path, None, root)
            required = {"annotation_run_manifest", "annotation_genes", "annotation_proteins", "screen_manifest", "screen_table", "dna_rescue_script", "reference_metadata", "reference_fasta", "tblastn_binary"}
            if not required <= set(before_bound):
                raise ValueError("Required annotation/screen binding file is missing")
            base["bound_input_sha256"] = before_bound
            base["bound_input_paths"] = {k: str(v) for k, v in _bound_path_objects(annotation, screen_dir, run_path, screen_path, None, root).items()}
        if not raw_mode and run.get("status") != "complete":
            raise ValueError("Source annotation is not complete")
        if not raw_mode:
            source, source_reason = _input_source(annotation, run)
        else:
            source_reason = None
        if source is None:
            base.update(status="not_applicable", reason=source_reason)
            _write_json(manifest_path, base)
            return base
        if not source.exists():
            base.update(status="not_evaluated", reason="source_nucleotide_input_missing", input=str(source))
            _write_json(manifest_path, base)
            return base
        source_digest = sha(source)
        declared = next((x.get("sha256") for x in run.get("input_files", []) if Path(x.get("path", "")).resolve() == source.resolve()), None)
        if not raw_mode and not declared:
            raise ValueError("Source nucleotide input has no declared hash")
        if declared and declared != source_digest:
            raise ValueError("Source nucleotide input hash changed")
        if not raw_mode:
            before_bound = _bound_paths(annotation, screen_dir, run_path, screen_path, source, root)
            base["bound_input_sha256"] = before_bound
            base["bound_input_paths"] = {k: str(v) for k, v in _bound_path_objects(annotation, screen_dir, run_path, screen_path, source, root).items()}
        source_format = _source_format(source)
        subject = materialize_fasta(source, out)
        genbank_tables = set()
        if source_format == "genbank":
            subject = out / "nucleotide_subject.fna"
            genbank_tables = _genbank_subject(source, subject)
        contigs = _fasta(subject)
        refs, refseqs, refs_meta, refs_faa = _reference_query(root)
        query = out / "reference_queries.faa"
        with query.open("w") as handle:
            for rid, seq in refseqs.items():
                handle.write(f">{rid}\n{seq}\n")
        genes = {} if raw_mode else _load_genes(annotation)
        circular_lengths = {} if raw_mode else _circular_lengths(genes, contigs)
        offsets = {name: _least_rotation_offset(contigs[name]) for name in circular_lengths}
        search_subject, search_contigs = subject, contigs
        seam_contigs = {}
        seam_subject = None
        if offsets:
            search_subject = out / "nucleotide_search_subject.fna"
            search_contigs = _search_contigs(contigs, offsets)
            with search_subject.open("w") as handle:
                for name, dna in search_contigs.items():
                    handle.write(f">{name}\n{dna}\n")
            seam_subject = out / "nucleotide_seam_subject.fna"
            seam_contigs = _seam_fasta(search_contigs, circular_lengths, seam_subject)
        screened = set()
        if raw_mode:
            screened, proteins = set(), {}
        else:
            with (screen_dir / "candidate_screen.tsv").open() as handle:
                screened = {r["protein_id"] for r in csv.DictReader(handle, delimiter="\t")}
            proteins = _load_proteins(annotation / "proteins.faa")
        observed_tables = [] if raw_mode else _prodigal_tables(annotation)
        tables = sorted(set(observed_tables) | genbank_tables | {11, 4})
        bound_context_hashes = {"source_nucleotide": source_digest, "references_json": sha(refs_meta), "references_faa": sha(refs_faa), "query_fasta": sha(query)}
        if run_path is not None:
            bound_context_hashes["run_manifest"] = sha(run_path)
        if screen_path is not None:
            bound_context_hashes["screen_manifest"] = sha(screen_path)
        base.update(input=str(source), contig_count=len(contigs), translation_tables=tables,
                    table_evidence={"prodigal_tables": observed_tables, "requested_tables": [11, 4],
                                    "genbank_tables": sorted(genbank_tables),
                                    "table_4_is_hypothesis": True}, reference_count=len(refseqs),
                    input_sha256=bound_context_hashes)
        version_path = _write_version(root, out)
        hsp_rows, commands, seam_commands, logs = [], [], [], [version_path]
        scans = [(False, search_subject, search_contigs)]
        if seam_subject is not None:
            scans.append((True, seam_subject, seam_contigs))
        for table in tables:
            for seam, scan_subject, scan_contigs in scans:
                suffix = "_seam" if seam else ""
                raw = out / f"tblastn_table{table}{suffix}.tsv"
                max_targets = len(scan_contigs) if seam else max(5, len(scan_contigs))
                command, log = _run_tblastn(root, query, scan_subject, table, raw, threads, max_targets)
                (seam_commands if seam else commands).append(command)
                logs.append(log)
                for row_index, line in enumerate(raw.read_text().splitlines(), 1):
                    if not line.strip():
                        continue
                    row = _parse_hsp(line, refseqs, scan_contigs)
                    contig = row["sseqid"]
                    if seam and _normalize_seam_row(row, circular_lengths[contig]) is None:
                        continue
                    aa, nt, _ = _oriented_translation(scan_contigs[contig], row["sstart"], row["send"], row["sframe"], table)
                    if row["sseq_ungapped"] != aa:
                        raise ValueError("BLAST subject alignment does not match translated DNA")
                    start, end, frame = _native_hsp_span(row, offsets.get(contig, 0))
                    row.update(native_start=start, native_end=end, translation=aa)
                    represented = _representation(row, genes, proteins, screened, contigs, table)
                    raw_hash = hashlib.sha256(line.encode()).hexdigest()
                    namespace = "seam|" if seam else ""
                    rid = hashlib.sha256((f"{namespace}{table}|{row_index}|{line}").encode()).hexdigest()[:20]
                    row.update(hsp_id=rid, translation_table=table, table_hypothesis=(table == 4),
                               raw_row_index=row_index, raw_row_sha256=raw_hash,
                               contig=contig, dna_start=start, dna_end=end, frame=frame,
                               translation_stop_positions=[i + 1 for i, x in enumerate(aa) if x == "*"],
                               dna_sha256=hashlib.sha256(nt.encode()).hexdigest(),
                               translation_sha256=hashlib.sha256(aa.encode()).hexdigest(),
                               represented_by=represented,
                               represented_start_codon_normalized=bool(row.get("represented_start_codon_normalized", False)),
                               source_hsp_id=rid)
                    if offsets:
                        row.update(circular_seam=end > len(contigs[contig]),
                                   source_subject="circular_seam" if seam else "canonical_search")
                    hsp_rows.append(row)
        # Keep the primary evidence ahead of any duplicate from the doubled subject.
        unique, seen = [], set()
        for row in hsp_rows:
            key = (row["qseqid"], row["sseqid"], row["dna_start"], row["dna_end"],
                   row["frame"], row["translation_table"], row["translation"])
            if row.get("source_subject") == "circular_seam" and key in seen:
                continue
            seen.add(key)
            unique.append(row)
        ledger_fields = ["hsp_id","raw_row_index","raw_row_sha256","qseqid","sseqid","pident","length","qstart","qend","sstart","send","evalue","bitscore","qlen","slen","qseq","sseq","sframe","translation_table","table_hypothesis","contig","dna_start","dna_end","frame","translation","translation_stop_positions","dna_sha256","translation_sha256","represented_by","represented_start_codon_normalized"]
        if seam_subject is not None:
            ledger_fields += ["seam_sstart","seam_send","seam_slen","circular_seam","source_subject"]
        ledger_fields += ["profile_evidence"]
        ledger = [{k: json.dumps(r[k], ensure_ascii=False) if isinstance(r.get(k), list) else r.get(k, "") for k in ledger_fields if k != "profile_evidence"} for r in unique]
        for row in ledger:
            row["profile_evidence"] = "[]"
        _write_tsv(out / "dna_hsp_ledger.tsv", ledger, ledger_fields)
        candidates_by_span = {}
        for row in unique:
            if row["represented_by"]:
                continue
            key = _candidate_key(row)
            candidate = candidates_by_span.setdefault(key, dict(candidate_id="", contig=row["sseqid"], start=row["dna_start"], end=row["dna_end"], strand=1 if row["frame"] > 0 else -1, frame=row["frame"], translation_tables=[], source_hsp_ids=[], reference_ids=[], families=[], DNA_sha256=row["dna_sha256"], translation_sha256=row["translation_sha256"], identity_status="unconfirmed", function_prediction="potential_homolog_function_unresolved", represented_by="", translation_has_stop=bool(row["translation_stop_positions"]), scope="translated_DNA_interval_not_reconstructed_CDS"))
            candidate["translation_tables"].append(row["translation_table"])
            candidate["source_hsp_ids"].append(row["hsp_id"])
            if row["qseqid"] not in candidate["reference_ids"]: candidate["reference_ids"].append(row["qseqid"])
            family=next((x.get("center", "") for x in refs if x["reference_id"] == row["qseqid"]), "")
            if family and family not in candidate["families"]: candidate["families"].append(family)
        candidates=[]
        for key, candidate in candidates_by_span.items():
            candidate["translation_tables"] = sorted(set(candidate["translation_tables"]))
            candidate["candidate_id"] = _candidate_id(key, candidate["translation_tables"][0])
            candidates.append(candidate)
        if len({x["candidate_id"] for x in candidates}) != len(candidates):
            raise ValueError("DNA candidate ID collision")
        out_faa = out / "dna_rescue_translations.faa"
        with out_faa.open("w") as handle:
            for row in unique:
                handle.write(f">{row['hsp_id']} table={row['translation_table']} frame={row['frame']}\n{row['translation'].replace('*','X')}\n")
        profile_dir = out / "profile_scan"
        profile_evidence = {}
        profile_artifacts = []
        profile_status = "not_performed_empty"
        if unique:
            profile_dir.mkdir()
            profile_paths = {
                "profile_model": root / "references/processed/catalytic_families.hmm",
                "profile_manifest": root / "references/processed/profile_manifest.json",
                "profile_policy": root / "config/discovery_evidence.json",
                "profile_script": root / "scripts/discovery_evidence.py",
                "profile_binary": resolve_tool(root, "hmmsearch"),
                "profile_query": out_faa,
            }
            profile_before = {k: sha(v) for k, v in profile_paths.items()}
            profile_artifacts_dir = profile_dir
            thresholds = load(root / "config/thresholds.json")
            search_discovery(out_faa, profile_artifacts_dir, root, 1, thresholds)
            domains = load_domains(profile_dir / "family.domtbl", {r["profile"]: r["family"] for r in load(profile_paths["profile_manifest"])})
            targets = load_targets(profile_dir / "family.tbl", {r["profile"]: r["family"] for r in load(profile_paths["profile_manifest"])})
            known_hsp = {r["hsp_id"] for r in unique}
            by_hsp = {r["hsp_id"]: r for r in unique}
            for row in domains:
                if row["protein_id"] not in known_hsp:
                    raise ValueError("Profile domain query is not a DNA HSP owner")
                profile_evidence.setdefault(row["protein_id"], []).append(dict(row, evidence_kind="raw_domain"))
            for row in targets:
                if row["protein_id"] not in known_hsp:
                    raise ValueError("Profile target query is not a DNA HSP owner")
                profile_evidence.setdefault(row["protein_id"], []).append(dict(
                    row, evidence_kind="raw_target",
                    query_length=len(by_hsp[row["protein_id"]]["translation"]),
                    query_length_evidence="bound_HSP_translation"))
            profile_after = {k: sha(v) for k, v in profile_paths.items()}
            if profile_after != profile_before:
                raise ValueError("DNA profile scan input changed")
            profile_artifacts = [profile_dir / "family.domtbl", profile_dir / "family.tbl", profile_dir / "hmmsearch.log"]
            profile_status = "reported_family_profile_evidence_only"
            base["profile_scan_input_sha256"] = profile_before
            base["profile_scan_input_paths"] = {k: str(v) for k, v in profile_paths.items()}
            base["profile_scan_output_sha256"] = {str(p.relative_to(out)): sha(p) for p in profile_artifacts}
        for candidate in candidates:
            evidence = []
            for hsp_id in candidate["source_hsp_ids"]:
                evidence.extend(profile_evidence.get(hsp_id, []))
            candidate["profile_evidence"] = evidence
            for item in evidence:
                family = item.get("family", "")
                if family and family not in candidate["families"]:
                    candidate["families"].append(family)
        for row in ledger:
            row["profile_evidence"] = json.dumps(profile_evidence.get(row["hsp_id"], []), ensure_ascii=False)
        _write_tsv(out / "dna_hsp_ledger.tsv", ledger, ledger_fields)
        _write_json(out / "dna_rescue_candidates.json", candidates)
        artifacts=[query, out / "dna_hsp_ledger.tsv", out / "dna_rescue_candidates.json", out_faa, version_path] + logs + profile_artifacts
        if subject != source: artifacts.append(subject)
        if seam_subject is not None:
            artifacts.extend([search_subject, seam_subject])
            base["circular_search_offsets"] = offsets
        artifacts += [out / f"tblastn_table{table}.tsv" for table in tables]
        if seam_subject is not None:
            artifacts += [out / f"tblastn_table{table}_seam.tsv" for table in tables]
        after_bound = _bound_paths(annotation, screen_dir, run_path, screen_path, source, root, raw=raw_mode)
        if after_bound != before_bound:
            raise ValueError("Bound annotation/source inputs changed during DNA rescue")
        base.update(status="complete", source_format=source_format, subject_fasta=str(subject) if subject != source else None,
                    profile_confirmation=profile_status, hsp_count=len(unique), raw_hsp_count=len(hsp_rows), candidate_count=len(candidates), candidates=len(candidates),
                    candidate_grouping=CANDIDATE_GROUPING, distinct_native_span_count=_span_count(candidates),
                    represented_hsp_count=sum(bool(r["represented_by"]) for r in unique), commands=commands,
                    seam_scan=bool(seam_subject), seam_commands=seam_commands,
                    output_sha256={str(p.relative_to(out)): sha(p) for p in artifacts})
        _write_json(manifest_path, base)
        return base
    except Exception as error:
        base.update(status="failed", error=f"{type(error).__name__}: {error}")
        _write_json(manifest_path, base)
        raise


def _verify_circular_hsps(out, ledger, tables, native, search, doubled, offsets):
    with (out / "reference_queries.faa").open() as handle:
        refs = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(handle, "fasta")}
    raw = {(table, route): (out / f"tblastn_table{table}{suffix}.tsv").read_text().splitlines()
           for table in tables for route, suffix in (("canonical_search", ""), ("circular_seam", "_seam"))}
    seen = set()
    for row in ledger.values():
        table, route, index = int(row["translation_table"]), row["source_subject"], int(row["raw_row_index"])
        lines = raw.get((table, route), [])
        owner = (table, route, index)
        if not 1 <= index <= len(lines) or owner in seen:
            raise ValueError("Circular HSP ledger has an invalid raw row owner")
        seen.add(owner)
        line = lines[index - 1]
        seam = route == "circular_seam"
        hit = _parse_hsp(line, refs, doubled if seam else search)
        contig = hit["sseqid"]
        dna = native[contig]
        if seam and _normalize_seam_row(hit, len(dna)) is None:
            raise ValueError("Circular seam ledger row is outside the search seam")
        start, end, frame = _native_hsp_span(hit, offsets.get(contig, 0))
        nt = circular_slice(dna, start, end)
        if frame < 0:
            nt = str(Seq(nt).reverse_complement())
        aa = hit["sseq_ungapped"]
        if (row["raw_row_sha256"] != hashlib.sha256(line.encode()).hexdigest()
                or row["contig"] != contig or row["qseqid"] != hit["qseqid"]
                or any(row[k] != str(hit[k]) for k in ("qstart", "qend", "qseq", "sseq", "sstart", "send", "sframe"))
                or tuple(int(row[k]) for k in ("dna_start", "dna_end", "frame")) != (start, end, frame)
                or row["circular_seam"] != str(end > len(dna))
                or row["dna_sha256"] != hashlib.sha256(nt.encode()).hexdigest()
                or row["translation_sha256"] != hashlib.sha256(aa.encode()).hexdigest()
                or row["translation"] != aa or str(Seq(nt).translate(table=table)) != aa):
            raise ValueError("Circular HSP ledger differs from its native DNA or raw alignment")


def load_candidates(out, annotation=None, screen_dir=None, genome=None, root=ROOT):
    out = Path(out).resolve()
    root = Path(root).resolve()
    manifest = load(out / "dna_rescue_manifest.json")
    status = manifest.get("status")
    if status not in {"complete", "not_applicable", "not_evaluated", "failed"}:
        raise ValueError("DNA rescue manifest has invalid status")
    raw = manifest.get("source_mode") == "raw_nucleotide_only"
    annotation = Path(annotation).resolve() if annotation is not None else (Path(manifest["annotation"]).resolve() if manifest.get("annotation") else None)
    screen_dir = Path(screen_dir).resolve() if screen_dir is not None else (Path(manifest["screen"]).resolve() if manifest.get("screen") else None)
    if annotation is not None and manifest.get("annotation") != str(annotation):
        raise ValueError("DNA rescue annotation ownership mismatch")
    if screen_dir is not None and manifest.get("screen") != str(screen_dir):
        raise ValueError("DNA rescue screen ownership mismatch")
    if raw and (annotation is not None or screen_dir is not None):
        raise ValueError("Raw DNA rescue cannot be consumed with annotation ownership")
    if raw and genome is not None:
        expected = str(Path(genome).resolve())
        if manifest.get("input") != expected:
            raise ValueError("Raw DNA rescue genome ownership mismatch")
    if raw:
        source = Path(genome).resolve() if genome is not None else (Path(manifest.get("input")).resolve() if manifest.get("input") else None)
        run_path = screen_path = None
        canonical = _bound_path_objects(None, None, None, None, source, root, raw=True)
    else:
        run_path = annotation / "run_manifest.json" if annotation else None
        screen_path = screen_dir / "candidate_screen_manifest.json" if screen_dir else None
        run = load(run_path) if run_path and run_path.exists() else {}
        source, _ = _input_source(annotation, run) if annotation else (None, "annotation_missing")
        canonical = _bound_path_objects(annotation, screen_dir, run_path, screen_path, source, root)
    canonical_paths = {k: str(v.resolve()) for k, v in canonical.items()}
    paths = {k: str(Path(v).resolve()) for k, v in manifest.get("bound_input_paths", {}).items()}
    bound = manifest.get("bound_input_sha256", {})
    if paths != canonical_paths:
        raise ValueError("DNA rescue bound input path ownership changed")
    if set(bound) != set(canonical_paths):
        raise ValueError("DNA rescue bound input keys are incomplete")
    for key, path_text in canonical_paths.items():
        path = Path(path_text)
        if not path.exists() or sha(path) != bound[key]:
            raise ValueError("DNA rescue bound input changed or missing: " + key)
    required = set(canonical_paths)
    if status == "complete":
        required.add("tblastn_binary")
        if not raw:
            required.add("source_nucleotide")
            if manifest.get("source_format") != "genbank":
                required.add("annotation_prodigal_gff")
    if not required <= set(bound) or not required <= set(paths):
        raise ValueError("DNA rescue bound input keys are incomplete")
    if status != "complete":
        if not manifest.get("reason") or manifest.get("candidate_count", 0) != 0:
            raise ValueError("Non-scanning DNA rescue status lacks an honest reason")
        if manifest.get("output_sha256"):
            raise ValueError("Non-scanning DNA rescue has output artifacts")
        if (out / "dna_rescue_candidates.json").exists():
            raise ValueError("Non-scanning DNA rescue has candidates")
        return manifest, []
    profile_status = manifest.get("profile_confirmation")
    if profile_status == "reported_family_profile_evidence_only":
        profile_inputs = manifest.get("profile_scan_input_sha256", {})
        expected_profile = {
            "profile_model": root / "references/processed/catalytic_families.hmm",
            "profile_manifest": root / "references/processed/profile_manifest.json",
            "profile_policy": root / "config/discovery_evidence.json",
            "profile_script": root / "scripts/discovery_evidence.py",
            "profile_binary": resolve_tool(root, "hmmsearch"),
            "profile_query": out / "dna_rescue_translations.faa",
        }
        if set(expected_profile) != {"profile_model", "profile_manifest", "profile_policy", "profile_script", "profile_binary", "profile_query"}:
            raise ValueError("DNA profile scan input binding is incomplete")
        declared_profile_paths = {k: str(Path(v).resolve()) for k, v in manifest.get("profile_scan_input_paths", {}).items()}
        if declared_profile_paths != {k: str(v.resolve()) for k, v in expected_profile.items()}:
            raise ValueError("DNA profile scan input path ownership changed")
        if set(profile_inputs) != set(expected_profile) or any(not p.exists() or sha(p) != profile_inputs[k] for k, p in expected_profile.items()):
            raise ValueError("DNA profile scan input binding changed or missing")
        profile_outputs = manifest.get("profile_scan_output_sha256", {})
        required_profile_outputs = {"profile_scan/family.domtbl", "profile_scan/family.tbl", "profile_scan/hmmsearch.log"}
        if set(profile_outputs) != required_profile_outputs or not required_profile_outputs <= set(manifest.get("output_sha256", {})):
            raise ValueError("DNA profile scan artifacts are incomplete")
        for name, digest in profile_outputs.items():
            if not (out / name).exists() or sha(out / name) != digest:
                raise ValueError("DNA profile scan output changed or missing: " + name)
    elif profile_status != "not_performed_empty":
        raise ValueError("DNA profile confirmation status is invalid")
    required_artifacts = {"dna_rescue_candidates.json", "dna_hsp_ledger.tsv", "dna_rescue_translations.faa", "reference_queries.faa", "tblastn_version.txt"}
    raw_keys = set(manifest.get("output_sha256", {}))
    canonical = "circular_search_offsets" in manifest
    if canonical != ("nucleotide_search_subject.fna" in raw_keys) or (canonical and not manifest.get("seam_scan")):
        raise ValueError("Circular canonical search binding is incomplete")
    for table in manifest.get("translation_tables", []):
        required_artifacts |= {f"tblastn_table{table}.tsv", f"tblastn_table{table}.log"}
        if manifest.get("seam_scan"):
            required_artifacts |= {f"tblastn_table{table}_seam.tsv", f"tblastn_table{table}_seam.log"}
    if manifest.get("seam_scan"):
        required_artifacts.add("nucleotide_seam_subject.fna")
    if not required_artifacts <= raw_keys:
        raise ValueError("DNA rescue complete artifacts are incomplete")
    for name, digest in manifest.get("output_sha256", {}).items():
        if not (out / name).exists() or sha(out / name) != digest:
            raise ValueError("DNA rescue output changed or missing: " + name)
    compressed_source = source.name.lower().endswith(".fna.gz")
    expects_subject = manifest.get("source_format") == "genbank" or compressed_source
    if expects_subject:
        expected_name = "nucleotide_subject.fna" if manifest.get("source_format") == "genbank" else ".input_genome.fna"
        subject_path = out / expected_name
        if manifest.get("subject_fasta") != str(subject_path.resolve()):
            raise ValueError("Derived nucleotide subject ownership changed")
        if expected_name not in raw_keys or not subject_path.exists():
            raise ValueError("Derived nucleotide subject artifact is missing")
        if not _fasta(subject_path):
            raise ValueError("Derived nucleotide subject is invalid")
    elif manifest.get("subject_fasta") is not None:
        raise ValueError("Plain FASTA cannot declare a derived subject")
    tables = manifest.get("translation_tables", [])
    if manifest.get("seam_scan"):
        if annotation is None:
            raise ValueError("Circular seam scan requires annotation ownership")
        native_path = Path(manifest["subject_fasta"]) if manifest.get("subject_fasta") else source
        native_contigs = _fasta(native_path)
        genes = _load_genes(annotation)
        circular_lengths = _circular_lengths(genes, native_contigs)
        search_contigs = native_contigs
        if canonical:
            offsets = manifest["circular_search_offsets"]
            if expects_subject:
                if manifest.get("source_format") == "genbank":
                    with source.open() as handle:
                        original = {r.id: str(r.seq).upper() for r in SeqIO.parse(handle, "genbank")}
                else:
                    original = _fasta(source)
                if native_contigs != original:
                    raise ValueError("Circular native subject differs from its source")
            if (not isinstance(offsets, dict) or any(type(offset) is not int for offset in offsets.values())
                    or not circular_lengths or offsets != {name: _least_rotation_offset(native_contigs[name]) for name in circular_lengths}):
                raise ValueError("Circular search offsets differ from native DNA")
            search_contigs = _fasta(out / "nucleotide_search_subject.fna")
            if search_contigs != _search_contigs(native_contigs, offsets):
                raise ValueError("Circular search subject differs from native rotation")
        seam_contigs = _fasta(out / "nucleotide_seam_subject.fna")
        if seam_contigs != {name: search_contigs[name] * 2 for name in circular_lengths}:
            raise ValueError("Circular seam subject does not match search sequence duplication")
    candidates = load(out / "dna_rescue_candidates.json")
    if manifest.get("candidate_count") != len(candidates) or manifest.get("candidates") != len(candidates):
        raise ValueError("DNA rescue candidate count is inconsistent")
    tables = manifest.get("translation_tables")
    if not isinstance(tables, list) or not tables or any(isinstance(t, bool) or not isinstance(t, int) or t <= 0 for t in tables):
        raise ValueError("DNA rescue translation tables are invalid")
    commands = manifest.get("commands", [])
    if len(commands) != len(tables):
        raise ValueError("DNA rescue commands do not match translation tables")
    for table, command in zip(tables, commands):
        if not isinstance(command, list) or not all(isinstance(item, str) for item in command) or "-db_gencode" not in command:
            raise ValueError("DNA rescue command record is invalid")
        pos = command.index("-db_gencode")
        if pos + 1 >= len(command) or command[pos + 1] != str(table):
            raise ValueError("DNA rescue command table does not match manifest")
        if canonical and ("-subject" not in command or command[command.index("-subject") + 1:][:1] != [str(out / "nucleotide_search_subject.fna")]):
            raise ValueError("Circular primary command is not bound to the search subject")
    if manifest.get("seam_scan"):
        seam_commands = manifest.get("seam_commands", [])
        if len(seam_commands) != len(tables) or any(not isinstance(command, list) for command in seam_commands):
            raise ValueError("DNA seam commands do not match translation tables")
        for table, command in zip(tables, seam_commands):
            if "-db_gencode" not in command or command[command.index("-db_gencode") + 1] != str(table):
                raise ValueError("DNA seam command table does not match manifest")
            if canonical and ("-subject" not in command or command[command.index("-subject") + 1:][:1] != [str(out / "nucleotide_seam_subject.fna")]):
                raise ValueError("Circular seam command is not bound to the search subject")
    grouping = manifest.get("candidate_grouping")
    if grouping is not None and grouping != CANDIDATE_GROUPING:
        raise ValueError("Unknown DNA candidate grouping")
    hsp_lengths = {}
    grouped_hsps = {}
    ledger_rows = {}
    with (out / "dna_hsp_ledger.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hid = row["hsp_id"]
            hsp_lengths[hid] = len(row.get("translation", ""))
            if grouping:
                if hid in grouped_hsps:
                    raise ValueError("Duplicate DNA HSP ID in candidate grouping")
                grouped_hsps[hid] = row
            if manifest.get("seam_scan"):
                ledger_rows[row["hsp_id"] if canonical else row["raw_row_sha256"]] = row
    if canonical:
        if len(ledger_rows) != manifest["hsp_count"]:
            raise ValueError("Circular HSP ledger row count differs from manifest")
        _verify_circular_hsps(out, ledger_rows, tables, native_contigs, search_contigs, seam_contigs, offsets)
    elif manifest.get("seam_scan"):
        refs = {record.id: str(record.seq).rstrip("*") for record in SeqIO.parse(out / "reference_queries.faa", "fasta")}
        crossing_hashes = set()
        for table in tables:
            seam_path = out / f"tblastn_table{table}_seam.tsv"
            for line in seam_path.read_text().splitlines():
                if not line.strip():
                    continue
                parsed = _parse_hsp(line, refs, seam_contigs)
                if _normalize_seam_row(parsed, circular_lengths[parsed["sseqid"]]) is None:
                    continue
                digest = hashlib.sha256(line.encode()).hexdigest()
                crossing_hashes.add(digest)
                owner = ledger_rows.get(digest)
                if owner is None or owner.get("sstart") != str(parsed["sstart"]) or owner.get("send") != str(parsed["send"]):
                    raise ValueError("Circular seam ledger is not bound to raw crossing HSP")
        for digest, row in ledger_rows.items():
            if row.get("circular_seam") in {"True", "true", "1"} and digest not in crossing_hashes:
                raise ValueError("Circular seam ledger row lacks a raw crossing HSP")
    if grouping:
        _verify_candidate_groups(candidates, grouped_hsps, manifest)
    for candidate in candidates:
        if candidate.get("identity_status") != "unconfirmed" or candidate.get("function_prediction") != "potential_homolog_function_unresolved":
            raise ValueError("DNA rescue candidate was promoted")
        for hsp_id in candidate.get("source_hsp_ids", []):
            if hsp_id not in hsp_lengths:
                raise ValueError("DNA rescue candidate references unknown HSP")
            if canonical:
                row = ledger_rows[hsp_id]
                fields = {"contig": "contig", "start": "dna_start", "end": "dna_end", "frame": "frame",
                          "DNA_sha256": "dna_sha256", "translation_sha256": "translation_sha256"}
                if row["represented_by"] or any(str(candidate[k]) != row[v] for k, v in fields.items()):
                    raise ValueError("Circular DNA candidate differs from its HSP")
        for evidence in candidate.get("profile_evidence", []):
            hsp_id = evidence.get("protein_id")
            if (hsp_id not in hsp_lengths or (grouping and hsp_id not in candidate["source_hsp_ids"])
                    or int(evidence.get("query_length", -1)) != hsp_lengths[hsp_id]):
                raise ValueError("DNA profile evidence query length is not the HSP peptide length")
    return manifest, candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--annotation", type=Path)
    group.add_argument("--genome", type=Path)
    parser.add_argument("--screen", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    if args.annotation is not None and args.screen is None:
        parser.error("--annotation requires --screen")
    if args.genome is not None and args.screen is not None:
        parser.error("--genome cannot receive --screen")
    if args.threads < 1:
        parser.error("--threads must be positive")
    try:
        print(json.dumps(prepare(args.annotation, args.screen, args.out, args.threads, args.root, genome=args.genome), ensure_ascii=False))
    except (OSError, ValueError, RuntimeError) as error:
        raise SystemExit(str(error))

