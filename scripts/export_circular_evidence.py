#!/usr/bin/env python3
"""Export measured coverage, clipping and junction alignments for circular contigs."""
import argparse
from array import array
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

try:
    import pysam
except ImportError:
    pysam = None

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.audit_circular_assembly import resolve_bam_reference

BAD_FLAGS = 4 | 256 | 512 | 1024 | 2048
READ_FIELDS = [
    "track", "qname", "read", "query_length", "start_0based", "end_0based",
    "relative_start", "relative_end", "strand", "mapq", "cigar", "flags",
    "mate_ref", "mate_start_0based", "template_length", "passes_filter",
    "spans_junction", "spans_junction_strict", "bridges_junction",
]


def resolve_reference(references, contig):
    matches = [r for r in references if r == contig or
               r.startswith((contig + "_", contig + "ctg_", contig + "."))]
    if contig in matches:
        return contig
    if len(matches) != 1:
        raise ValueError(f"Expected one reference for {contig!r}; found {matches}")
    return resolve_bam_reference(matches, contig)


def read_fasta(path, contig):
    records = {}
    name = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f"Duplicate FASTA reference: {name}")
                records[name] = []
            elif line:
                if name is None:
                    raise ValueError("Sequence before FASTA header")
                records[name].append(line)
    ref = resolve_reference(records, contig)
    sequence = "".join(records[ref]).upper()
    if len(records) != 1 or len(sequence) < 2:
        raise ValueError("Expected a single nonempty circular contig FASTA")
    return ref, sequence


def alignment_blocks(read):
    """Merge reference-contiguous M/=/X blocks, including across insertions."""
    position = read.reference_start
    blocks = []
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8):
            if blocks and blocks[-1][1] == position:
                blocks[-1] = (blocks[-1][0], position + length)
            else:
                blocks.append((position, position + length))
            position += length
        elif op in (2, 3):
            position += length
    return blocks


def passes_filter(read, min_mapq):
    return not read.flag & BAD_FLAGS and read.mapping_quality >= min_mapq


def span_junction(read, boundary, anchor, max_indel=0):
    blocks = alignment_blocks(read)
    left, right = boundary - anchor, boundary + anchor
    if not blocks or blocks[0][0] > left or blocks[-1][1] < right:
        return False
    if max_indel == 0:
        return any(s <= left and e >= right for s, e in blocks)
    # ONT sequencing indels must not turn every otherwise supported junction into a negative.
    position = read.reference_start
    for op, length in read.cigartuples or ():
        if op in (2, 3) and max(position, left) < min(position + length, right):
            if op == 3 or length > max_indel:
                return False
        if op == 1 and left <= position <= right and length > max_indel:
            return False
        if op in (0, 2, 3, 7, 8):
            position += length
    return all(sum(max(0, min(e, high) - max(s, low)) for s, e in blocks) >= 0.70 * anchor
               for low, high in ((left, boundary), (boundary, right)))


def add_blocks(delta, blocks, start, end):
    for left, right in blocks:
        left, right = max(start, left), min(end, right)
        if left < right:
            delta[left - start] += 1
            delta[right - start] -= 1


def accumulate(delta):
    total = 0
    for position in range(len(delta) - 1):
        total += delta[position]
        delta[position] = total
    delta.pop()
    return delta


def whole_coverage(path, ref, length, min_mapq, collect_clips=False):
    filtered = array("q", [0]) * (length + 1)
    all_alignments = array("q", [0]) * (length + 1)
    clips = Counter()
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(ref):
            if read.is_unmapped:
                continue
            blocks = alignment_blocks(read)
            add_blocks(all_alignments, blocks, 0, length)
            if passes_filter(read, min_mapq):
                add_blocks(filtered, blocks, 0, length)
            cigar = read.cigartuples or []
            if collect_clips and blocks and cigar:
                if cigar[0][0] in (4, 5):
                    clips[(blocks[0][0], "left")] += 1
                if cigar[-1][0] in (4, 5):
                    clips[(blocks[-1][1] - 1, "right")] += 1
    return accumulate(filtered), accumulate(all_alignments), clips


def coverage_stats(values):
    return sum(values) / len(values), sum(v > 0 for v in values) / len(values)


def coverage_bins(track, filtered, all_depth, size):
    for start in range(0, len(filtered), size):
        end = min(len(filtered), start + size)
        row = dict(track=track, start_0based=start, end_0based=end)
        for prefix, values in (("", filtered[start:end]), ("all_", all_depth[start:end])):
            mean, fraction = coverage_stats(values)
            row.update({prefix + "mean_depth": mean,
                        prefix + "zero_bases": values.count(0),
                        prefix + "covered_fraction": fraction})
        yield row


def zero_regions(track, policy, values):
    start = None
    for pos in range(len(values) + 1):
        zero = pos < len(values) and values[pos] == 0
        if zero and start is None:
            start = pos
        elif not zero and start is not None:
            yield dict(track=track, policy=policy, start_0based=start,
                       end_0based=pos, length=pos - start)
            start = None


def bridging_pairs(reads, boundary, min_mapq):
    pairs = defaultdict(list)
    for read in reads:
        if (passes_filter(read, min_mapq) and read.is_paired and read.is_proper_pair
                and not read.mate_is_unmapped and read.reference_id == read.next_reference_id):
            pairs[read.query_name].append(read)
    bridges = set()
    for qname, records in pairs.items():
        if len(records) != 2 or sum(r.is_read1 for r in records) != 1 or sum(r.is_read2 for r in records) != 1:
            continue
        left, right = sorted(records, key=lambda r: r.reference_start)
        if (left.is_reverse or not right.is_reverse or
                left.next_reference_start != right.reference_start or
                right.next_reference_start != left.reference_start):
            continue
        fragment_length = max(left.reference_end, right.reference_end) - left.reference_start
        if (50 <= fragment_length <= 2000 and left.template_length == fragment_length
                and right.template_length == -fragment_length
                and left.reference_start < boundary <= right.reference_start):
            bridges.add(qname)
    return bridges


def junction_evidence(path, ref, boundary, window, anchor, min_mapq, track):
    with pysam.AlignmentFile(path, "rb") as bam:
        start = max(0, boundary - window)
        end = min(bam.get_reference_length(ref), boundary + window)
        reads = list(bam.fetch(ref, start, end))
    # shortcut: both mates must occur in the exported window; enlarge it for long insert libraries.
    bridges = bridging_pairs(reads, boundary, min_mapq) if track == "short" else set()
    delta = array("q", [0]) * (end - start + 1)
    spanning = set()
    strict_spanning = set()
    rows = []
    for read in reads:
        passed = passes_filter(read, min_mapq)
        strict = passed and span_junction(read, boundary, anchor)
        spans = passed and span_junction(read, boundary, anchor, 50 if track == "long" else 0)
        if passed:
            add_blocks(delta, alignment_blocks(read), start, end)
        if spans:
            spanning.add(read.query_name)
        if strict:
            strict_spanning.add(read.query_name)
        rows.append(dict(
            track=track, qname=read.query_name,
            read="R1" if read.is_read1 else "R2" if read.is_read2 else "NA",
            query_length=read.query_length or 0, start_0based=read.reference_start,
            end_0based=read.reference_end, relative_start=read.reference_start - boundary,
            relative_end=read.reference_end - boundary, strand="-" if read.is_reverse else "+",
            mapq=read.mapping_quality, cigar=read.cigarstring or "", flags=read.flag,
            mate_ref=read.next_reference_name or "", mate_start_0based=read.next_reference_start,
            template_length=read.template_length, passes_filter=passed, spans_junction=spans,
            spans_junction_strict=strict,
            bridges_junction=passed and read.query_name in bridges,
        ))
    return start, accumulate(delta), rows, spanning, strict_spanning, bridges


def write_table(path, fields, rows):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def validate(args):
    if pysam is None:
        raise RuntimeError("pysam is required; run in the longread_audit environment")
    for key in ("bin_size", "junction_window", "long_anchor", "short_anchor"):
        if getattr(args, key) <= 0:
            raise ValueError(f"{key} must be positive")
    if not 0 <= args.min_mapq <= 255:
        raise ValueError("min_mapq must be between 0 and 255")
    paths = {k: getattr(args, k) for k in ("contig_fa", "long_bam", "long_junction_bam",
             "short_bam", "short_rotated_bam", "rotation_coords")}
    if args.self_blast:
        paths["self_blast"] = args.self_blast
    for path in paths.values():
        if not path or not Path(path).is_file():
            raise FileNotFoundError(path)
    fasta_ref, sequence = read_fasta(args.contig_fa, args.contig)
    length = len(sequence)
    with open(args.rotation_coords) as handle:
        coords = [r for r in csv.DictReader(handle, delimiter="\t") if r.get("contig") == args.contig]
    if len(coords) != 1:
        raise ValueError("Expected exactly one rotation coordinate row for target")
    row = coords[0]
    offset, short_boundary = int(row["rotation_offset"]), int(row["rotated_junction_coord"])
    if int(row["length"]) != length or not 0 < offset < length or short_boundary != length - offset:
        raise ValueError("Rotation coordinates do not match FASTA length / L-offset boundary")
    refs = {}
    lengths = {}
    for key in ("long_bam", "short_bam", "short_rotated_bam", "long_junction_bam"):
        with pysam.AlignmentFile(paths[key], "rb") as bam:
            bam.check_index()
            ref = resolve_reference(bam.references, args.contig)
            refs[key], lengths[key] = ref, bam.get_reference_length(ref)
            if key != "long_junction_bam" and lengths[key] != length:
                raise ValueError(f"BAM reference length mismatch: {key}")
    if row.get("seq_id") and row["seq_id"] != refs["short_rotated_bam"]:
        raise ValueError("Rotation reference differs from rotated BAM reference")
    junction_length = lengths["long_junction_bam"]
    if not 2 <= junction_length <= length or junction_length % 2:
        raise ValueError("Synthetic junction reference must have an even length <= contig length")
    return paths, refs, fasta_ref, sequence, offset, short_boundary, junction_length // 2


def export(args):
    out = Path(args.out_dir)
    # A failed rerun must not leave an old success marker beside partial outputs.
    (out / "manifest.json").unlink(missing_ok=True)
    paths, refs, fasta_ref, sequence, offset, short_boundary, long_boundary = validate(args)
    out.mkdir(parents=True, exist_ok=True)
    length = len(sequence)
    coverage_rows, gaps, depth_rows, read_rows, summary = [], [], [], [], []
    clip_rows = []
    for track, whole_key, junction_key, boundary, anchor in (
        ("long", "long_bam", "long_junction_bam", long_boundary, args.long_anchor),
        ("short", "short_bam", "short_rotated_bam", short_boundary, args.short_anchor),
    ):
        filtered, all_depth, clips = whole_coverage(paths[whole_key], refs[whole_key], length,
                                                   args.min_mapq, track == "long")
        coverage_rows.extend(coverage_bins(track, filtered, all_depth, args.bin_size))
        gaps.extend(zero_regions(track, "filtered_primary", filtered))
        gaps.extend(zero_regions(track, "all_alignments", all_depth))
        for (position, side), count in sorted(clips.items()):
            support = all_depth[position]
            clip_rows.append(dict(position_0based=position, side=side, clipped_alignments=count,
                                  depth=support, clipping_ratio=count / support if support else "NA",
                                  distance_to_end_bp=min(position, length - 1 - position)))
        start, depth, reads, spanning, strict_spanning, bridges = junction_evidence(
            paths[junction_key], refs[junction_key], boundary, args.junction_window,
            anchor, args.min_mapq, track)
        read_rows.extend(reads)
        for i, value in enumerate(depth):
            position = start + i
            original = (position + offset) % length if track == "short" else (position - boundary) % length
            depth_rows.append(dict(track=track, relative_pos=position - boundary,
                                   reference_pos_0based=position, original_pos_0based=original, depth=value))
        genome_mean, genome_fraction = coverage_stats(filtered)
        all_mean, all_fraction = coverage_stats(all_depth)
        junction_mean, junction_fraction = coverage_stats(depth)
        summary.append(dict(track=track, genome_mean_depth=genome_mean,
                            genome_covered_fraction=genome_fraction,
                            all_alignment_mean_depth=all_mean,
                            all_alignment_covered_fraction=all_fraction,
                            junction_window_mean_depth=junction_mean,
                            junction_window_covered_fraction=junction_fraction,
                            spanning_unique_qnames=len(spanning),
                            strict_spanning_unique_qnames=len(strict_spanning), bridge_unique_fragments=len(bridges)))
    tables = {
        "genome_coverage.tsv.gz": (["track", "start_0based", "end_0based", "mean_depth", "zero_bases",
            "covered_fraction", "all_mean_depth", "all_zero_bases", "all_covered_fraction"], coverage_rows),
        "junction_depth.tsv.gz": (["track", "relative_pos", "reference_pos_0based", "original_pos_0based", "depth"], depth_rows),
        "junction_reads.tsv.gz": (READ_FIELDS, read_rows),
        "evidence_summary.tsv": (list(summary[0]), summary),
        "zero_coverage_regions.tsv": (["track", "policy", "start_0based", "end_0based", "length"], gaps),
        "long_clipping.tsv": (["position_0based", "side", "clipped_alignments", "depth", "clipping_ratio", "distance_to_end_bp"], clip_rows),
    }
    for filename, (fields, rows) in tables.items():
        write_table(out / filename, fields, rows)
    manifest = {
        "schema_version": 1, "contig": args.contig, "fasta_reference": fasta_ref, "length": length,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "bam_references": refs,
        "coordinates": {"intervals": "0-based, half-open", "relative_pos": "0 is the first base after L->1",
            "self_blast": "Unmodified BLAST outfmt 6 uses 1-based inclusive coordinates, possibly descending on reverse strand",
            "long_junction_boundary": long_boundary, "short_junction_boundary": short_boundary,
            "short_rotation_offset": offset, "clipping_position": "first/last M,=,X base for left/right clipping"},
        "filters": {"min_mapq": args.min_mapq, "excluded_flags": BAD_FLAGS,
            "all_alignments": "All mapped alignments including secondary/supplementary/MAPQ0; alignment support, not molecule abundance",
            "depth": "M,=,X only; D/N excluded; overlapping paired reads counted separately",
            "spanning": "Long: >=70% M,=,X in each anchor, no N or indel >50 bp. Short/strict: one reference-contiguous M,=,X block, I allowed. Unique QNAME counts; project thresholds.",
            "bridging": "Both mates in exported window; filtered primary proper FR pair, reciprocal positions, consistent TLEN 50..2000 bp"},
        "parameters": vars(args),
        "inputs": {key: {"path": str(Path(path).resolve()), "size": os.stat(path).st_size,
                         "mtime": os.stat(path).st_mtime} for key, path in paths.items()},
        "outputs": list(tables),
        "method": {"reference": "https://doi.org/10.1038/s41587-025-02971-8",
            "workflow": "https://merenlab.org/data/benchmarking-long-read-assemblers/",
            "scope": "Paper-informed raw support and clipping export, not a reproduction of anvi'o. Short reads, MAPQ and junction anchor thresholds are project supplements.",
            "limitations": ["Clipping counts are per alignment and side, not independent molecules or a biological error verdict.",
                "Zero coverage requires raw-read search and inspection; filtered-only gaps may reflect mapping ambiguity.",
                "Query length excludes hard-clipped sequence and is not necessarily the original FASTQ read length.",
                "Synthetic long junction mappings are target-only; high MAPQ there does not establish whole-assembly uniqueness.",
                "Full-genome tracks use original linearized coordinates; use junction tracks to inspect artificial endpoints."]},
        "note": "Measured evidence only; no simulated values or automatic biological circularity decision.",
    }
    with open(out / "manifest.json", "w") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(f"Exported {args.contig} ({length} bp) to {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("contig", "contig_fa", "long_bam", "long_junction_bam", "short_bam",
                 "short_rotated_bam", "rotation_coords", "out_dir"):
        parser.add_argument("--" + name.replace("_", "-"), required=True)
    parser.add_argument("--self-blast", default="", help="Existing BLAST outfmt 6 file, recorded as provenance")
    parser.add_argument("--short-read-context", default="unspecified", help="Sampling relationship to the long reads")
    for name, default in (("bin_size", 1000), ("junction_window", 5000), ("min_mapq", 20),
                          ("long_anchor", 500), ("short_anchor", 30)):
        parser.add_argument("--" + name.replace("_", "-"), type=int, default=default)
    args = parser.parse_args()
    try:
        export(args)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"error: {exc}\n")


if __name__ == "__main__":
    main()
