"""Evidence-only raw-FNA scanner for a reviewed split-center model."""
from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

if __package__:
    from .dna_rescue import _oriented_translation, _parse_hsp, _run_tblastn
else:
    from dna_rescue import _oriented_translation, _parse_hsp, _run_tblastn
from circular_coordinates import circular_slice, interval_segments


def _sha(value):
    return hashlib.sha256(value.encode()).hexdigest()


def _circular(value, length):
    return (int(value) - 1) % length + 1


def _native_interval(first, nt_length, strand, length):
    if strand > 0:
        start = _circular(first, length)
    else:
        start = _circular(first - nt_length + 1, length)
    return start, start + nt_length - 1


def _native_segments(start, end, length):
    return [list(pair) for pair in interval_segments(start, end, length)]


def _oriented_window(sequence, anchor, strand, width, upstream):
    length = len(sequence)
    topology = str(anchor.get("topology", "")).lower()
    if topology == "circular":
        width = min(width, length)
        def get(start, end):
            span = end - start
            start = _circular(start, length)
            return circular_slice(sequence, start, start + span)
    else:
        def get(start, end):
            if end < 1 or start > length:
                return ""
            return sequence[max(1, start) - 1:min(length, end)]
    if strand > 0:
        first = int(anchor["start"]) - width if upstream else int(anchor["end"]) + 1
        raw = get(first, first + width - 1)
    else:
        first = int(anchor["end"]) + width if upstream else int(anchor["start"]) - 1
        raw = get(first - width + 1, first)
        raw = str(Seq(raw).reverse_complement()) if raw else ""
    return raw, first


def _alignment_sites(row, nt, raw_aa, declared, window_start, window_end):
    offsets = {}
    qpos = row["qstart"] - 1
    subject_aa = covered = matches = 0
    for qchar, schar in zip(row["qseq"], row["sseq"]):
        if qchar != "-":
            qpos += 1
        if schar == "-":
            continue
        if qchar != "-":
            if window_start <= qpos <= window_end:
                covered += 1
                matches += qchar == schar
            if qpos in declared:
                codon = nt[subject_aa * 3:subject_aa * 3 + 3]
                aa = raw_aa[subject_aa]
                if declared[qpos] == "U" and (aa != "*" or codon != "TGA"):
                    return None
                if declared[qpos] == "C" and aa != "C":
                    return None
                offsets[qpos] = subject_aa
        subject_aa += 1
    if covered and offsets.keys() == declared.keys():
        return offsets, covered / (window_end - window_start + 1), 100 * matches / covered


def _candidate_starts(seq, hsp_low, hsp_high, max_nt, internal_stops):
    table = CodonTable.unambiguous_dna_by_id[11]
    starts, stops = set(table.start_codons), set(table.stop_codons)
    result = []
    for start in range(hsp_low, max(1, hsp_low - max_nt) - 1, -3):
        if seq[start - 1:start + 2] not in starts:
            continue
        terminal = None
        for pos in range(start, min(len(seq) - 2, start + max_nt - 1) + 1, 3):
            codon = seq[pos - 1:pos + 2]
            if set(codon) - {"A", "C", "G", "T"}:
                break
            if codon in stops and pos not in internal_stops:
                if pos < hsp_high:
                    break
                terminal = pos
                break
        if terminal is not None and terminal >= hsp_high:
            result.append((start, terminal, terminal + 2))
    return result


def _complete_orf_window(query, conditional, raw_translation, nt, sites, expected_offsets,
                         window_start, window_end):
    """Test a complete coding hypothesis without crediting its conditional X residues."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    mapping = None
    for count, alignment in enumerate(aligner.align(query, conditional)):
        if count >= 128:
            return None
        current = {}
        for (qstart, qend), (sstart, send) in zip(*alignment.aligned):
            current.update(zip(range(qstart, qend), range(sstart, send)))
        if mapping is not None and current != mapping:
            return None
        mapping = current
    if mapping is None or any(mapping.get(pos - 1) != offset for pos, offset in expected_offsets.items()):
        return None
    for position, residue in sites.items():
        offset = mapping.get(position - 1)
        if offset is None:
            return None
        if residue == "U" and (raw_translation[offset] != "*" or nt[offset * 3:offset * 3 + 3] != "TGA"):
            return None
        if residue == "C" and raw_translation[offset] != "C":
            return None
    paired = [(q, s) for q, s in mapping.items() if window_start <= q + 1 <= window_end]
    if not paired:
        return None
    matches = sum(query[q] == conditional[s] and query[q] != "X" for q, s in paired)
    return len(paired) / (window_end - window_start + 1), 100 * matches / len(paired)


def scan(anchors, contigs, model, parameters, root, *, upstream=False):
    """Scan same-orientation raw DNA beside a reviewed coding anchor."""
    sequence = model["sequence"]
    reference_id = model["reference_id"]
    sites = {s["reference_position"]: s["reference_residue"] for s in model["sites"]}
    sec_positions = sorted(position for position, residue in sites.items() if residue == "U")
    if not sec_positions:
        raise ValueError("Conditional component search requires a reviewed U position")
    window_start, window_end = model["window_start"], model["window_end"]
    max_evalue = float(parameters["max_evalue"])
    min_cov = float(parameters["min_window_coverage"])
    min_identity = float(parameters["min_window_identity"])
    max_gap = int(parameters["max_intergenic_gap_bp"])
    max_span = int(parameters["max_span_bp"])
    flank = int(parameters.get("site_context_flank", 15))
    max_precursor = (len(sequence) + flank) * 3
    query_sequence = sequence.replace("U", "X")
    root = Path(root).resolve()
    with tempfile.TemporaryDirectory(prefix="split_component_scan_") as tmp_name:
        tmp = Path(tmp_name)
        query = tmp / "query.faa"
        query.write_text(f">{reference_id}\n{query_sequence}\n")
        subject = tmp / "windows.fna"
        metadata = {}
        with subject.open("w") as handle:
            for index, anchor in enumerate(anchors):
                contig = str(anchor["contig"])
                topology = str(anchor.get("topology", "")).lower()
                if contig not in contigs or topology not in {"linear", "circular"}:
                    continue
                length = len(contigs[contig])
                if int(anchor.get("contig_length") or length) != length:
                    continue
                strand = int(anchor["strand"])
                anchor_len = abs(int(anchor["end"]) - int(anchor["start"])) + 1
                width = max(0, min(max_span - anchor_len, max_gap + max_precursor))
                if topology == "circular":
                    width = min(width, max(0, length - anchor_len))
                else:
                    left = int(anchor["start"]) - 1
                    right = length - int(anchor["end"])
                    remaining = left if (strand > 0) == upstream else right
                    width = min(width, max(0, remaining))
                if width <= 0:
                    continue
                oriented, first = _oriented_window(contigs[contig], anchor, strand, width, upstream)
                if not oriented:
                    continue
                sid = f"split_{index}"
                metadata[sid] = (anchor, contig, strand, first, oriented)
                handle.write(f">{sid}\n{oriented}\n")
        result = {str(a["id"]): [] for a in anchors}
        if not metadata:
            return result
        out = tmp / "tblastn.tsv"
        _run_tblastn(root, query, subject, 11, out, 1, max(5, len(metadata)))
        refs = {reference_id: query_sequence}
        subjects = {sid: data[4] for sid, data in metadata.items()}
        rows = []
        for line_no, line in enumerate(out.read_text().splitlines(), 1):
            if not line.strip():
                continue
            row = _parse_hsp(line, refs, subjects)
            if row["evalue"] > max_evalue or row["sframe"] < 0:
                continue
            anchor, contig, strand, native_first, oriented = metadata[row["sseqid"]]
            aa, nt, _ = _oriented_translation(oriented, row["sstart"], row["send"], row["sframe"], 11)
            if row["sseq_ungapped"] != aa:
                raise ValueError("Split-component alignment does not match its translated window")
            alignment = _alignment_sites(row, nt, aa, sites, window_start, window_end)
            if alignment is None:
                continue
            offsets, coverage, identity = alignment
            if identity < min_identity:
                continue
            hsp_low, hsp_high = sorted((row["sstart"], row["send"]))
            internal_stops = {hsp_low + 3 * offsets[position] for position in sec_positions}
            for start, terminal_start, terminal_end in _candidate_starts(oriented, hsp_low, hsp_high, max_precursor, internal_stops):
                gap = len(oriented) - terminal_end if upstream else start - 1
                if gap > max_gap:
                    continue
                full_nt = oriented[start - 1:terminal_end]
                raw_translation = str(Seq(full_nt).translate(table=11))
                codons = [(position - start) // 3 for position in sorted(internal_stops)]
                observed = [i for i, aa in enumerate(raw_translation[:-1]) if aa == "*"]
                if observed != codons or not raw_translation.endswith("*"):
                    continue
                conditional = "M" + raw_translation[1:-1].replace("*", "X")
                method = "local_HSP"
                window_coverage, window_identity = coverage, identity
                if coverage < min_cov:
                    expected_offsets = {position: (hsp_low + 3 * offset - start) // 3
                                        for position, offset in offsets.items()}
                    complete = _complete_orf_window(query_sequence, conditional, raw_translation, full_nt,
                                                    sites, expected_offsets, window_start, window_end)
                    if complete is None or complete[0] < min_cov or complete[1] < min_identity:
                        continue
                    window_coverage, window_identity = complete
                    method = "complete_ORF_global_uncalibrated_v1"
                oriented_first = native_first + (start - 1) if strand > 0 else native_first - (start - 1)
                nstart, nend = _native_interval(oriented_first, len(full_nt), strand, len(contigs[contig]))
                length = len(contigs[contig])
                if strand > 0:
                    tga_genomic = [_circular(oriented_first + codon * 3, length) for codon in codons]
                    terminal_genomic = _circular(oriented_first + terminal_start - start, length)
                    frame = (nstart - 1) % 3 + 1
                else:
                    tga_genomic = [_circular(oriented_first - codon * 3 - 2, length) for codon in codons]
                    terminal_genomic = _circular(oriented_first - (terminal_start - start) - 2, length)
                    frame = -((length - nend) % 3 + 1)
                if upstream:
                    tga_fields = dict(tga_sites=[
                        dict(reference_position=ref_pos, position=codon + 1, genomic_start=genomic)
                        for ref_pos, codon, genomic in zip(sec_positions, codons, tga_genomic)
                    ], upstream_gap_bp=gap)
                else:
                    tga_fields = dict(tga_position=codons[0] + 1, tga_genomic_start=tga_genomic[0])
                rows.append(dict(anchor_id=str(anchor["id"]), component_id=model["component_id"],
                    reference_id=reference_id, reference_accession=model["reference_accession"],
                    contig=contig, start=nstart, end=nend, strand=strand, frame=frame,
                    nt_length=len(full_nt), nt_sequence=full_nt, DNA_sha256=_sha(full_nt),
                    translation_table=11, translation_raw=raw_translation,
                    translation_conditional=conditional, translation_conditional_sha256=_sha(conditional),
                    **tga_fields,
                    native_segments=_native_segments(nstart, nend, length),
                    terminal_stop=full_nt[-3:], terminal_stop_genomic_start=terminal_genomic,
                    hsp_qstart=row["qstart"], hsp_qend=row["qend"], hsp_sstart=row["sstart"],
                    hsp_send=row["send"], hsp_sframe=row["sframe"], hsp_pident=row["pident"],
                    hsp_length=row["length"], hsp_evalue=row["evalue"], hsp_bitscore=row["bitscore"],
                    hsp_line=line_no, hsp_qseq=row["qseq"], hsp_sseq=row["sseq"],
                    window_coverage=window_coverage, window_identity=window_identity,
                    window_alignment_method=method, hsp_window_coverage=coverage, hsp_window_identity=identity))
    for row in rows:
        bucket = result[row["anchor_id"]]
        key = (row["component_id"], row["contig"], row["start"], row["end"], row["strand"])
        if not any((x["component_id"], x["contig"], x["start"], x["end"], x["strand"]) == key for x in bucket):
            bucket.append(row)
    return result
