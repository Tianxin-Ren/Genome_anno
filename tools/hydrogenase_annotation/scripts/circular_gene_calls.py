"""Recover Prodigal CDSs across a declared circular replicon's origin."""
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from itertools import chain
from pathlib import Path
from urllib.parse import quote

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from circular_coordinates import circular_slice, interval_segments


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def call_prodigal(source, faa, gff, log, table):
    with log.open("w") as handle:
        subprocess.run(["prodigal", "-i", str(source), "-a", str(faa), "-o", str(gff),
                        "-f", "gff", "-p", "meta", "-g", str(table)], stderr=handle, check=True)


def protein_call(record):
    _, start, end, strand, attrs = record.description.split(" # ", 4)
    tags = dict(item.split("=", 1) for item in attrs.split(";") if "=" in item)
    return int(start), int(end), int(strand), tags


def covered_edge_fragment(record, start, end, strand, length):
    low, high, sign, tags = protein_call(record)
    if sign != strand or tags["partial"] == "00":
        return False
    if low <= 3 and high <= end - length:
        low, high = low + length, high + length
    elif high < length - 2 or low < start:
        return False
    return ((low - start) if strand == 1 else (end - high)) % 3 == 0


def bounded_edge_candidate(partial, candidate, length):
    """Match an intact endpoint using the copy next to the central seam."""
    low, high, sign, partial_tags = partial
    start, end, strand, tags = candidate
    if strand != sign or tags["partial"] != "00":
        return None
    if partial_tags["partial"] == "01" and length - 2 <= high <= length and end <= length:
        if low == start <= end <= high and (high - end) % 3 == 0:
            return start, end
    if partial_tags["partial"] == "10" and 1 <= low <= 3 and length < start:
        if low <= start - length <= end - length == high and (start - length - low) % 3 == 0:
            return start - length, end - length
    return None


def model_tables(gff, contigs):
    """Read model headers by sequence ordinal, even for a replicon without CDSs."""
    tables = {}
    for line in gff.read_text().splitlines():
        if line.startswith("# Sequence Data:"):
            number, length = re.match(r"# Sequence Data: seqnum=(\d+);seqlen=(\d+);", line).groups()
            contig, expected_length = contigs[int(number) - 1]
            if int(length) != expected_length:
                raise ValueError("Prodigal sequence length differs from its input")
        elif line.startswith("# Model Data:"):
            table = re.search(r"(?:^|;)transl_table=(\d+)(?:;|$)", line)
            if not table:
                raise ValueError("Prodigal model lacks a translation table")
            tables[contig] = int(table.group(1))
    return tables


def recover_circular_calls(source, records, topologies, output, table):
    """Recover complete seam calls and replace uniquely matched edge fragments."""
    trace = output / "prodigal_circular"
    trace.mkdir()
    native_faa, native_gff = trace / "native.faa", trace / "native.gff"
    seam_fna, seam_faa, seam_gff = (trace / f"seam.{suffix}" for suffix in ("fna", "faa", "gff"))
    source_hash = sha(source)
    call_prodigal(source, native_faa, native_gff, output / "prodigal.log", table)
    circular = {r.id: str(r.seq).upper() for r in records if topologies[r.id] == "circular"}
    SeqIO.write((SeqRecord(Seq(nt + nt), id=contig, description="")
                 for contig, nt in circular.items()), seam_fna, "fasta")
    call_prodigal(seam_fna, seam_faa, seam_gff, trace / "seam.log", table)
    if sha(source) != source_hash:
        raise ValueError("Nucleotide source changed during circular gene calling")

    native_tables = model_tables(native_gff, [(r.id, len(r)) for r in records])
    seam_tables = model_tables(seam_gff, [(contig, 2 * len(nt)) for contig, nt in circular.items()])
    native = list(SeqIO.parse(native_faa, "fasta"))
    by_contig = defaultdict(list)
    for record in native:
        by_contig[record.id.rsplit("_", 1)[0]].append(record)
    seam_rows = {}
    for line in seam_gff.read_text().splitlines():
        if not line.startswith("#") and line.strip():
            fields = line.split("\t")
            key = fields[0], fields[8].split(";", 1)[0].removeprefix("ID=")
            seam_rows[key] = fields

    spanning, edge_calls = [], []
    fragments_by_id = defaultdict(list)
    complete = {}
    edge_fragments = {}
    for contig, dna in circular.items():
        complete[contig] = set()
        edge_fragments[contig] = []
        for record in by_contig[contig]:
            call = protein_call(record)
            if call[3]["partial"] == "00":
                complete[contig].add(call[:3])
            elif ((call[3]["partial"] == "10" and call[0] <= 3) or
                  (call[3]["partial"] == "01" and call[1] >= len(dna) - 2)):
                edge_fragments[contig].append((record, call))
    for record in SeqIO.parse(seam_faa, "fasta"):
        contig = record.id.rsplit("_", 1)[0]
        call = protein_call(record)
        length = len(circular[contig])
        if call[0] <= length < call[1]:
            spanning.append((record, contig, call, call[:2], None))
            continue
        matches = []
        for fragment, fragment_call in edge_fragments[contig]:
            coords = bounded_edge_candidate(fragment_call, call, length)
            if coords and (*coords, call[2]) not in complete[contig]:
                matches.append((fragment, coords))
                fragments_by_id[fragment.id].append(record.id)
        if len(matches) == 1:
            fragment, coords = matches[0]
            edge_calls.append((record, contig, call, coords, fragment))

    retained, removed, rejected = defaultdict(list), set(), []
    replacements = {}
    for record, contig, call, coords, fragment in chain(spanning, edge_calls):
        dna = circular[contig]
        length = len(dna)
        raw_start, raw_end, strand, tags = call
        start, end = coords
        if fragment and (len(fragments_by_id[fragment.id]) != 1 or
                         (contig, protein_call(fragment)[3]["ID"]) in removed):
            continue
        code = native_tables.get(contig)
        if (tags["partial"] != "00" or code is None or code != seam_tables.get(contig)
                or end - start + 1 > length):
            rejected.append(dict(contig=contig, start=start, end=end, strand=strand,
                                 reason="partial_model_conflict_or_multiple_turns"))
            continue
        nt = Seq(circular_slice(dna, start, end))
        if strand == -1:
            nt = nt.reverse_complement()
        try:
            protein = str(nt.translate(table=code, cds=True))
        except (TranslationError, ValueError):
            protein = None
        if not protein or protein != str(record.seq).rstrip("*"):
            rejected.append(dict(contig=contig, start=start, end=end, strand=strand,
                                 reason="translation_does_not_match_source"))
            continue
        if fragment:
            pid = fragment.id
            internal_id = protein_call(fragment)[3]["ID"]
        else:
            ordinal = max((int(r.id.rsplit("_", 1)[1]) for r in by_contig[contig]), default=0)
            ordinal += len(retained[contig]) + 1
            pid = f"{contig}_{ordinal}"
            native_index = next(index for index, r in enumerate(records, 1) if r.id == contig)
            internal_id = f"{native_index}_{ordinal}"
        fields = list(seam_rows[contig, tags["ID"]])
        if fields[2] != "CDS" or fields[3:5] != [str(raw_start), str(raw_end)] or fields[6:8] != [
                "+" if strand == 1 else "-", "0"]:
            raise ValueError("Prodigal seam FAA/GFF ownership mismatch")
        fields[3:5] = [str(start), str(end)]
        marker = "within_boundary_edge_partial" if fragment else "origin_spanning"
        fields[8] = re.sub(r"^ID=[^;]+", f"ID={internal_id}", fields[8]) + f"{marker}=true;"
        record_out = SeqRecord(record.seq, id=pid,
                               description=f"{pid} # {start} # {end} # {strand} # {fields[8]}")
        fragments = [fragment] if fragment else [
            r for r in by_contig[contig] if covered_edge_fragment(r, start, end, strand, length)]
        proof = dict(protein_id=pid, contig=contig, start=start, end=end, strand=strand,
                     segments=interval_segments(start, end, length), effective_table=code,
                     oriented_nucleotide_sha256=hashlib.sha256(str(nt).encode()).hexdigest(),
                     protein_sha256=hashlib.sha256(protein.encode()).hexdigest(),
                     seam_protein_id=record.id, replaced_partial_protein_ids=[r.id for r in fragments])
        row = record_out, "\t".join(fields), proof
        if fragment:
            proof.update(recovery_type=marker, seam_coordinates=[raw_start, raw_end])
            replacements[contig, internal_id] = row
        else:
            removed.update((contig, protein_call(r)[3]["ID"]) for r in fragments)
            retained[contig].append(row)

    faa, gff = output / "prodigal.faa", output / "prodigal.gff"
    with faa.open("w") as handle:
        for record in native:
            key = record.id.rsplit("_", 1)[0], protein_call(record)[3]["ID"]
            if key not in removed:
                SeqIO.write(replacements[key][0] if key in replacements else record, handle, "fasta")
        SeqIO.write((r for hits in retained.values() for r, _, _ in hits), handle, "fasta")
    lines = native_gff.read_text().splitlines()
    with gff.open("w") as handle:
        print(lines[0], file=handle)
        for contig, nt in circular.items():
            print(f"##sequence-region {contig} 1 {len(nt)}", file=handle)
            print(f"##topology {contig} circular", file=handle)
        records_iter = iter(records)
        current = None
        for line in lines[1:]:
            if line.startswith("# Sequence Data:"):
                for _, row, _ in retained.get(current, ()):
                    print(row, file=handle)
                current = next(records_iter).id
            if not line.startswith("#") and line.strip():
                fields = line.split("\t")
                key = fields[0], fields[8].split(";", 1)[0].removeprefix("ID=")
                if key in removed:
                    continue
                if key in replacements:
                    line = replacements[key][1]
            print(line, file=handle)
            if line.startswith("# Model Data:") and current in circular:
                contig = quote(current, safe=".:^*$@!+_?-|")
                print(f"{current}\tworkflow\tregion\t1\t{len(circular[current])}\t.\t+\t.\tID={contig};Is_circular=true",
                      file=handle)
        for _, row, _ in retained.get(current, ()):
            print(row, file=handle)
    proof = dict(schema="prodigal_circular_origin_v1", coordinate_system="one_based_inclusive_unwrapped",
                 source=dict(path=str(source.resolve()), sha256=source_hash),
                 caller_script_sha256=sha(Path(__file__)),
                 circular_contigs={k: len(v) for k, v in circular.items()},
                 recovered=[row for hits in retained.values() for _, _, row in hits] +
                           [row for _, _, row in replacements.values()], rejected=rejected,
                 artifact_sha256={str(p.relative_to(output)): sha(p) for p in
                                  (native_faa, native_gff, seam_fna, seam_faa, seam_gff, faa, gff)})
    (output / "circular_origin.json").write_text(json.dumps(proof, indent=2) + "\n")
