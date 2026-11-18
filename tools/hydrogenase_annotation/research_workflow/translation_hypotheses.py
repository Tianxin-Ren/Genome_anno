"""Audit internal UGA at NiFe ligand sites using the bound DNA search."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

if __package__:
    from .dna_rescue import _fasta, _load_genes, _load_proteins, load, load_candidates, prodigal_models, sha
else:
    from dna_rescue import _fasta, _load_genes, _load_proteins, load, load_candidates, prodigal_models, sha
from circular_coordinates import circular_overlap, circular_slice, interval_segments

RULE = "NiFe_ligand_UGA_v3"
MAX_DOWNSTREAM_NT = 300
STOP_CODONS = frozenset(CodonTable.unambiguous_dna_by_id[11].stop_codons)


def ligand_uga(row, dna, ligand_sites):
    """Require four aligned ligands and continuous stop-free flanks at each UGA."""
    if (len(ligand_sites) != 4 or float(row["evalue"]) > 1e-10
            or (row.get("translation_table") == "11" and "*" not in row["sseq"])):
        return None
    ref_pos, aa_pos = int(row["qstart"]) - 1, -1
    mapped = {}
    for ref, aa in zip(row["qseq"], row["sseq"]):
        ref_pos += ref != "-"
        aa_pos += aa != "-"
        if ref != "-" and aa != "-" and ref_pos in ligand_sites:
            mapped[ref_pos] = aa_pos
    if len(mapped) != 4:
        return None
    start, end, frame = (int(row[key]) for key in ("dna_start", "dna_end", "frame"))
    crossing = end > len(dna)
    if crossing and str(row.get("circular_seam")).lower() != "true":
        return None
    segment = Seq(circular_slice(dna, start, end) if crossing else dna[start - 1:end])
    nt = str(segment if frame > 0 else segment.reverse_complement())
    standard = str(Seq(nt).translate(table=11))
    sites = [mapped[p] for p in sorted(ligand_sites)]
    stops = [i for i, aa in enumerate(standard) if aa == "*"]
    if not 1 <= len(stops) <= 2 or any(i not in sites for i in stops):
        return None
    if any(standard[i] != "C" for i in sites if i not in stops):
        return None
    for i in stops:
        if nt[3 * i:3 * i + 3] != "TGA" or i < 15 or i + 15 >= len(standard):
            return None
    starts = [start + 3 * i if frame > 0 else end - 3 * i - 2 for i in stops]
    sec = "".join("U" if aa == "*" else aa for aa in standard)
    a, b, c, d = (sec[i] for i in sites)
    proposal = dict(start=start, end=end, frame=frame, raw_translation_table_11=standard,
                conditional_sec_translation=sec, uga_codon_starts_1based=starts,
                ligand_pattern=f"{a}XX{b}/{c}XX{d}")
    if crossing:
        proposal["native_segments"] = interval_segments(start, end, len(dna))
    return proposal


def coding_extent(proposal, dna, genes, proteins):
    """Extend a verified ordinary CDS across only the aligned ligand UGAs."""
    forward = proposal["frame"] > 0
    uga_starts = proposal["uga_codon_starts_1based"]
    anchors = []
    for gene in genes:
        start, end = int(gene["start"]), int(gene["end"])
        if (gene.get("partial") != "False" or gene.get("coordinate_evidence") != "genomic_CDS"
                or not 1 <= start <= end <= len(dna)
                or (end - 2 if forward else start) != uga_starts[0]):
            continue
        offset = proposal["start"] - start if forward else end - proposal["end"]
        if offset < 0 or offset % 3:
            continue
        nt = Seq(dna[start - 1:end])
        if not forward:
            nt = nt.reverse_complement()
        try:
            ordinary = str(nt.translate(table=11, cds=True))
        except (ValueError, CodonTable.TranslationError):
            continue
        if ordinary == proteins.get(gene["id"]) and "X" not in ordinary:
            anchors.append((gene["id"], start, end, ordinary))
    if len(anchors) != 1:
        return {"status": "ambiguous_ordinary_cds_anchor" if anchors else "no_verified_ordinary_cds_anchor"}

    protein_id, start, end, ordinary = anchors[0]
    ordinary_cds = dict(protein_id=protein_id, start=start, end=end, protein_length=len(ordinary),
                        protein_sha256=hashlib.sha256(ordinary.encode()).hexdigest())
    # shortcut: bound the search beyond the aligned span; longer tails need additional alignment evidence.
    left = start - 1 if forward else max(0, proposal["start"] - 1 - MAX_DOWNSTREAM_NT)
    right = min(len(dna), proposal["end"] + MAX_DOWNSTREAM_NT) if forward else end
    region = dna[left:right]
    nt = region if forward else str(Seq(region).reverse_complement())
    required_nt = proposal["end"] - start + 1 if forward else end - proposal["start"] + 1
    uga_offsets = {pos - start if forward else end - pos - 2 for pos in uga_starts}
    observed = set()
    for i in range(end - start - 2, len(nt) - 2, 3):
        codon = nt[i:i + 3]
        if set(codon) - set("ACGT"):
            return dict(status="ambiguous_sequence_before_terminal_stop", ordinary_cds=ordinary_cds)
        if i in uga_offsets and codon == "TGA":
            observed.add(i)
        elif codon in STOP_CODONS:
            if i < required_nt or observed != uga_offsets:
                return dict(status="stop_conflicts_with_aligned_span", ordinary_cds=ordinary_cds)
            coding = nt[:i + 3]
            raw = str(Seq(coding).translate(table=11))
            sec = "M" + raw[1:-1].replace("*", "U")
            return dict(status="conditional_start_to_stop", ordinary_cds=ordinary_cds,
                        start=start if forward else end - i - 2, end=start + i + 2 if forward else end,
                        terminal_stop=codon, nucleotide_sha256=hashlib.sha256(coding.encode()).hexdigest(),
                        raw_translation_table_11=raw, conditional_sec_translation=sec)
    at_edge = right == len(dna) if forward else left == 0
    return dict(status="contig_edge_before_terminal_stop" if at_edge else "terminal_stop_outside_search_limit",
                ordinary_cds=ordinary_cds)


def circular_coding_extent(proposal, dna, genes, proteins):
    """Lift a verified seam span into a contiguous window for the same CDS audit."""
    length = len(dna)
    lifted = dict(proposal, start=proposal["start"] + length, end=proposal["end"] + length,
                  uga_codon_starts_1based=[p + length for p in proposal["uga_codon_starts_1based"]])
    copies = [dict(g, start=int(g["start"]) + offset, end=int(g["end"]) + offset)
              for g in genes if g.get("topology") == "circular" and int(g.get("contig_length") or 0) == length
              and g.get("coordinate_evidence") == "genomic_CDS"
              for offset in (0, length, 2 * length)]
    result = coding_extent(lifted, dna * 3, copies, proteins)
    if "start" in result and result["end"] - result["start"] + 1 > length:
        result = dict(status="conditional_extent_exceeds_one_turn", ordinary_cds=result["ordinary_cds"])
    for region in (result, result.get("ordinary_cds", {})):
        if "start" in region:
            span = region["end"] - region["start"]
            start = (region["start"] - 1) % length + 1
            region.update(start=start, end=start + span, native_segments=interval_segments(start, start + span, length))
    return result


def audit(dna_rescue, output):
    dna_rescue, output = Path(dna_rescue), Path(output)
    manifest, candidates = load_candidates(dna_rescue)
    manifest_path = dna_rescue / "dna_rescue_manifest.json"
    source = manifest.get("input")
    result = dict(rule=RULE, status=manifest["status"], source=source,
                  dna_rescue_manifest_sha256=sha(manifest_path), script_sha256=sha(__file__), hypotheses=[])
    if manifest["status"] != "complete":
        result["reason"] = manifest["reason"]
    else:
        metadata = Path(manifest["bound_input_paths"]["reference_metadata"])
        refs = {r["reference_id"]: r["diagnostic_positions"] for r in load(metadata)
                if r.get("kind") == "catalytic" and r.get("center") == "NiFe"}
        contigs = _fasta(manifest.get("subject_fasta") or source)
        candidates_by_hsp = {hsp: c["candidate_id"] for c in candidates for hsp in c["source_hsp_ids"]}
        machinery = defaultdict(set)
        if manifest["source_format"] == "genbank":
            for record in SeqIO.parse(source, "genbank"):
                machinery[record.id].update(gene.lower() for f in record.features
                                            for gene in f.qualifiers.get("gene", [])
                                            if gene.lower() in {"sela", "selb", "selc", "seld"})
        genes = defaultdict(list)
        proteins, tables = {}, {}
        if manifest.get("annotation"):
            for gene in _load_genes(manifest["annotation"]).values():
                genes[gene["contig"], gene["strand"]].append(gene)
            gff = manifest["bound_input_paths"].get("annotation_prodigal_gff")
            if gff:
                tables, _ = prodigal_models(gff)
                proteins = _load_proteins(manifest["bound_input_paths"]["annotation_proteins"])
        groups = defaultdict(list)
        with (dna_rescue / "dna_hsp_ledger.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sites = refs.get(row["qseqid"])
                if not sites:
                    continue
                contig = row["contig"]
                proposal = ligand_uga(row, contigs[contig], sites)
                if proposal:
                    key = (contig, 1 if proposal["frame"] > 0 else -1, tuple(proposal["uga_codon_starts_1based"]))
                    groups[key].append((row, proposal))
        for (contig, strand, _), matches in sorted(groups.items()):
            matches.sort(key=lambda hit: (-float(hit[0]["bitscore"]), hit[0]["qseqid"], hit[0]["hsp_id"]))
            row, proposal = matches[0]
            dna = contigs[contig]
            circular = "native_segments" in proposal
            local_genes = genes[contig, str(strand)]
            overlap = [dict(id=g["id"], start=int(g["start"]), end=int(g["end"]))
                       for g in local_genes if int(g["start"]) > 0
                       and (circular_overlap((int(g["start"]), int(g["end"])), (proposal["start"], proposal["end"]), len(dna))
                            if circular else int(g["start"]) <= proposal["end"] and int(g["end"]) >= proposal["start"])]
            extent = {"status": "table11_cds_model_not_available"}
            if tables.get(contig) == 11:
                extent = (circular_coding_extent if circular else coding_extent)(proposal, dna, local_genes, proteins)
            result["hypotheses"].append(dict(
                **proposal, contig=contig, strand=strand, representative_hsp=row["hsp_id"],
                supporting_hsps=[hit["hsp_id"] for hit, _ in matches],
                mapped_catalytic_references=sorted({hit["qseqid"] for hit, _ in matches}),
                source_dna_candidate_ids=sorted({candidates_by_hsp[hit["hsp_id"]] for hit, _ in matches
                                                 if hit["hsp_id"] in candidates_by_hsp}),
                overlapping_cds=sorted(overlap, key=lambda g: (g["start"], g["end"], g["id"])),
                coding_extent=extent,
                annotated_same_replicon_sec_genes=sorted(machinery[contig]), SECIS="not_assessed",
                translation_status="conditional_aligned_span", identity_status="unconfirmed",
                competing_explanations=["termination_or_pseudogene", "sequence_or_assembly_error",
                                        "alternative_genetic_code", "non_Sec_recoding"],
            ))
        if __package__:
            from .split_component_translation import audit as audit_split_components
        else:
            from split_component_translation import audit as audit_split_components
        result["split_component_translation"] = audit_split_components(manifest, contigs)
    with output.open("x") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dna-rescue", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    audit(args.dna_rescue, args.out)
