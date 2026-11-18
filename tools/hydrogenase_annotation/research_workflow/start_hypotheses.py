"""Record conditional adjacent start boundaries for assigned raw-FNA members."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

from .dna_rescue import _fasta, _load_genes, _load_proteins, load, load_candidates, prodigal_models, sha
from circular_coordinates import circular_slice, interval_segments, valid_cds_interval

START_SEARCH_BOUND_NT = 12
START_SEARCH_STEP_NT = 3
START_POLICY = "in_frame_start_offsets_within_plus_or_minus_12nt; table_valid_cds_translation; unchanged_stop"
OVERLAP_SCOPE = "Upstream alternatives only; same-contig valid genomic_CDS neighbors; diagnostics only"


def _sequence_hash(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def adjacent_starts(dna, gene, protein, table):
    """Validate the called CDS and try bounded in-frame, table-valid starts."""
    if (gene["partial"] != "False" or gene["coordinate_evidence"] != "genomic_CDS"
            or table is None or not valid_cds_interval(gene, allow_unknown_length=False)):
        return []
    start, end, strand = (int(gene[key]) for key in ("start", "end", "strand"))
    length = len(dna)
    if strand not in (-1, 1) or int(gene["contig_length"]) != length:
        return []
    circular = gene["topology"] == "circular"

    def coding(left, right):
        if not 6 <= right - left + 1 <= length:
            return None
        if circular:
            span = right - left
            left = (left - 1) % length + 1
            right = left + span
            nt = circular_slice(dna, left, right)
        elif 1 <= left <= right <= length:
            nt = dna[left - 1:right]
        else:
            return None
        if strand == -1:
            nt = str(Seq(nt).reverse_complement())
        if set(nt) - set("ACGT"):
            return None
        try:
            translated = str(Seq(nt).translate(table=table, cds=True))
        except (TranslationError, KeyError, ValueError):
            return None
        return dict(start=left, end=right, native_segments=interval_segments(left, right, length),
                    start_codon=nt[:3], stop_codon=nt[-3:], dna_sha256=_sequence_hash(nt),
                    protein_length=len(translated), protein_sha256=_sequence_hash(translated), translation=translated)

    called = coding(start, end)
    if called is None or called["translation"] != protein:
        return []
    alternatives = []
    for offset in range(-START_SEARCH_BOUND_NT, START_SEARCH_BOUND_NT + 1, START_SEARCH_STEP_NT):
        if offset == 0:
            continue
        variant = coding(start + offset, end) if strand == 1 else coding(start, end - offset)
        if variant and variant["stop_codon"] == called["stop_codon"]:
            alternatives.append(dict(variant, start_offset_nt=offset))
    return [called, *alternatives] if alternatives else []


def _neighbor_overlaps(genes_by_contig, gene, called, alternative, length):
    """Report only newly added upstream bases overlapping valid same-contig CDSs."""
    if alternative.get("start_offset_nt", 0) >= 0:
        return []

    def intersection(left, right):
        return sum(max(0, min(a2, b2) - max(a1, b1) + 1)
                   for a1, a2 in left for b1, b2 in right)

    strand = int(gene["strand"])
    added_nt = -alternative["start_offset_nt"]
    first_added = called["start"] - added_nt if strand == 1 else called["end"] + 1
    first_added = (first_added - 1) % length + 1
    added_segments = interval_segments(first_added, first_added + added_nt - 1, length)
    called_segments = called["native_segments"]
    alternative_segments = alternative["native_segments"]
    hits = []
    for neighbor in genes_by_contig.get(gene["contig"], ()):
        if neighbor.get("id") == gene.get("id") or neighbor.get("coordinate_evidence") != "genomic_CDS":
            continue
        try:
            if not valid_cds_interval(neighbor, allow_unknown_length=False) or int(neighbor["contig_length"]) != length:
                continue
            segments = interval_segments(int(neighbor["start"]), int(neighbor["end"]), length)
            neighbor_strand = int(neighbor["strand"])
        except (KeyError, TypeError, ValueError):
            continue
        if neighbor_strand not in (-1, 1):
            continue
        added_overlap = intersection(added_segments, segments)
        if not added_overlap:
            continue
        hits.append(dict(protein_id=neighbor["id"],
                         strand_relation="same" if neighbor_strand == strand else "opposite",
                         added_overlap_nt=added_overlap,
                         called_overlap_nt=intersection(called_segments, segments),
                         alternative_overlap_nt=intersection(alternative_segments, segments),
                         neighbor_native_segments=segments,
                         partial=neighbor.get("partial") == "True"))
    return sorted(hits, key=lambda hit: (hit["protein_id"], hit["strand_relation"], hit["neighbor_native_segments"]))


def audit(dna_rescue, claims_path, output):
    """Publish a diagnostic sidecar without changing the called members."""
    dna_rescue, claims_path, output = map(Path, (dna_rescue, claims_path, output))
    manifest_path = dna_rescue / "dna_rescue_manifest.json"
    manifest = load(manifest_path)
    result = dict(schema="member_start_hypotheses_v1", status="not_applicable", records=[],
                  start_search_bound_nt=START_SEARCH_BOUND_NT, start_search_step_nt=START_SEARCH_STEP_NT,
                  start_policy=START_POLICY,
                  overlap_scope=OVERLAP_SCOPE,
                  scope="Conditional start boundaries only; physical start selection is unresolved.")
    if (manifest["status"] == "complete" and manifest.get("source_format") == "fasta"
            and manifest.get("source_mode") != "raw_nucleotide_only"):
        manifest, _ = load_candidates(dna_rescue)
        paths = manifest["bound_input_paths"]
        bound = manifest["bound_input_sha256"]
        result.update(status="complete", script_sha256=sha(__file__),
                      dna_rescue_manifest_sha256=sha(manifest_path), claim_assertions_sha256=sha(claims_path),
                      source_inputs={key: {"path": paths[key], "sha256": bound[key]} for key in (
                          "source_nucleotide", "annotation_genes", "annotation_proteins", "annotation_prodigal_gff")})
        claims = load(claims_path)
        if claims["source_binding"]["status"] != "producer_manifest_verified":
            raise ValueError("Selected members lack producer-bound claims")
        owners = defaultdict(list)
        for cid, claim in sorted(claims["assertions"].items()):
            if claim["association_status"] == "supported":
                for member in claim["member_evidence"]:
                    owners[member["protein_id"]].append(dict(candidate_id=cid, role=member["role"]))
        genes = _load_genes(manifest["annotation"])
        genes_by_contig = defaultdict(list)
        for item in genes.values():
            genes_by_contig[item.get("contig")].append(item)
        proteins = _load_proteins(paths["annotation_proteins"])
        contigs = _fasta(paths["source_nucleotide"])
        tables, _ = prodigal_models(paths["annotation_prodigal_gff"])
        result["selected_member_count"] = len(owners)
        for pid, assignments in sorted(owners.items()):
            gene = genes[pid]
            contig = gene["contig"]
            table = tables.get(contig)
            candidates = adjacent_starts(contigs[contig], gene, proteins[pid], table)
            if candidates:
                called, *alternatives = candidates
                for alternative in alternatives:
                    alternative["added_cds_overlaps"] = _neighbor_overlaps(
                        genes_by_contig, gene, called, alternative, len(contigs[contig]))
                result["records"].append(dict(protein_id=pid, assignments=assignments, contig=contig,
                    contig_length=len(contigs[contig]), strand=int(gene["strand"]), topology=gene["topology"],
                    translation_table=table, translation_table_source="prodigal_model_header",
                    start_search_bound_nt=START_SEARCH_BOUND_NT, start_search_step_nt=START_SEARCH_STEP_NT,
                    start_policy=START_POLICY,
                    overlap_scope=OVERLAP_SCOPE,
                    status="conditional_boundary", called=called, alternatives=alternatives))
        if (any(sha(item["path"]) != item["sha256"] for item in result["source_inputs"].values())
                or sha(claims_path) != result["claim_assertions_sha256"]
                or sha(manifest_path) != result["dna_rescue_manifest_sha256"]):
            raise ValueError("Start hypothesis source changed during audit")
    with output.open("x") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return result
