"""Bind local split-center DNA hypotheses to reviewed annotation evidence."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq

if __package__:
    from .dna_rescue import ROOT, _load_genes, _load_proteins, load, sha
    from .split_component_scan import scan
else:
    from dna_rescue import ROOT, _load_genes, _load_proteins, load, sha
    from split_component_scan import scan
from circular_coordinates import circular_overlap, circular_slice, valid_cds_interval
from tool_paths import resolve_tool


RULE = "split_component_internal_UGA_v2"
CANDIDATES = "catalytic_center_candidates.tsv"
COMPONENTS = "catalytic_center_components.tsv"
RELAY_MODEL = "research_workflow/upstream_relay_model.json"


def _rows(directory, name, expected):
    path = directory / name
    if sha(path) != expected[name]:
        raise ValueError(f"Split translation annotation output changed: {name}")
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _bound_anchor(row, components, genes, proteins, models, contigs, tables, *, bind_short=True):
    pid = row["catalytic_anchor"]
    a, partner = models
    gene = genes.get(pid)
    if (not gene or gene.get("partial") != "False" or gene.get("coordinate_evidence") != "genomic_CDS"
            or tables.get(gene["contig"]) != 11 or gene["strand"] not in ("1", "-1")):
        return None
    if len(components) != 2 or row["source_sha256"] != a["source_sha256"]:
        return None
    by_ref = {r["reference_id"]: r for r in components}
    large = by_ref.get(a["reference_id"], {})
    small = by_ref.get(partner["reference_id"], {})
    if (large.get("protein_id") != pid or large.get("template_span_status") != "supported"
            or large.get("component_status") != "template_ligand_set_supported"):
        return None
    partner_id = small.get("protein_id")
    fragment = None
    if bind_short:
        if row["center_status"] == "required_component_not_detected":
            if (row["organization"] != "incomplete" or json.loads(row["component_proteins"]) != [pid]
                    or partner_id != "" or small.get("component_status") != "not_detected_or_unpredicted"):
                return None
        else:
            fragment = genes.get(partner_id)
            if (row["organization"] != "multiple_polypeptides" or json.loads(row["component_proteins"]) != [pid, partner_id]
                    or not fragment or fragment.get("partial") != "True" or fragment["contig"] != gene["contig"]
                    or fragment["strand"] != gene["strand"] or fragment.get("coordinate_evidence") != "genomic_CDS"
                    or small.get("component_status") == "template_ligand_set_supported"):
                return None
    if any(c["source_sha256"] != a["source_sha256"] for c in components):
        return None
    seq = proteins.get(pid, "")
    sites = [s for s in json.loads(row["sites"]) if s.get("reference_id") == a["reference_id"]]
    expected = {s["reference_position"]: s["reference_residue"] for s in a["sites"]}
    if len(sites) != len(expected) or {s["reference_position"] for s in sites} != expected.keys():
        return None
    for site in sites:
        pos = site.get("query_position")
        if (site.get("status") != "present" or site.get("protein_id") != pid
                or not isinstance(pos, int) or not 1 <= pos <= len(seq)
                or site.get("observed_residue") != seq[pos - 1]
                or seq[pos - 1] != expected[site["reference_position"]]):
            return None
    if len({s["query_position"] for s in sites}) != len(sites):
        return None
    dna = contigs.get(gene["contig"])
    if not dna or not valid_cds_interval(gene, allow_unknown_length=False) or int(gene["contig_length"]) != len(dna):
        return None
    start, end, strand = (int(gene[k]) for k in ("start", "end", "strand"))
    nt = Seq(circular_slice(dna, start, end))
    if strand < 0:
        nt = nt.reverse_complement()
    try:
        observed = str(nt.translate(table=11, cds=True))
    except (ValueError, CodonTable.TranslationError):
        return None
    if observed != seq:
        return None
    fragment_evidence = None
    if fragment:
        if not valid_cds_interval(fragment, allow_unknown_length=False) or int(fragment["contig_length"]) != len(dna):
            return None
        fstart, fend = int(fragment["start"]), int(fragment["end"])
        fragment_dna = Seq(circular_slice(dna, fstart, fend))
        if strand < 0:
            fragment_dna = fragment_dna.reverse_complement()
        if len(fragment_dna) % 3 or str(fragment_dna[-3:]) != "TGA":
            return None
        peptide = str(fragment_dna[:-3].translate(table=11))
        if "*" in peptide or "X" in peptide or peptide != proteins.get(partner_id):
            return None
        stop_start = fend - 2 if strand > 0 else fstart
        fragment_evidence = dict(protein_id=partner_id, start=fstart, end=fend,
                                 terminal_TGA_start=(stop_start - 1) % len(dna) + 1,
                                 protein_sha256=hashlib.sha256(peptide.encode()).hexdigest())
    return dict(gene, start=start, end=end, strand=strand, center_id=row["center_id"],
                protein_sha256=hashlib.sha256(seq.encode()).hexdigest(), mapped_sites=sites,
                called_component_fragment=fragment_evidence)


def _overlaps(genes, contig, proposals, length):
    spans = [(r["start"], r["end"]) for r in proposals]
    overlap = []
    for gene in genes.values():
        if gene["contig"] != contig or not valid_cds_interval(gene, allow_unknown_length=False):
            continue
        interval = (int(gene["start"]), int(gene["end"]))
        if any(circular_overlap(span, interval, length) for span in spans):
            overlap.append(dict(protein_id=gene["id"], start=interval[0], end=interval[1], strand=int(gene["strand"])))
    return overlap


def _upstream_relay(anchors, contigs, genes, root):
    model_path = root / RELAY_MODEL
    model_digest = sha(model_path)
    model = load(model_path)
    sites = {site["reference_position"]: site["reference_residue"] for site in model["sites"]}
    sequence = model["sequence"]
    if (hashlib.sha256(sequence.encode()).hexdigest() != model["sequence_sha256"]
            or sorted(sites.values()) != ["C", "C", "C", "U", "U"]
            or any(sequence[pos - 1] != residue for pos, residue in sites.items())):
        raise ValueError("Invalid reviewed upstream relay model")
    results = scan(anchors, contigs, model, model["parameters"], root, upstream=True)
    if sha(model_path) != model_digest:
        raise ValueError("Upstream relay model changed during search")

    result = dict(status="evaluated_conditional_upstream_relay", model_sha256=model_digest, hypotheses=[], searches=[])
    groups = defaultdict(lambda: defaultdict(list))
    owners = defaultdict(set)
    for anchor in anchors:
        pid = anchor["id"]
        for proposal in results[pid]:
            key = (anchor["contig"], anchor["strand"],
                   tuple(site["genomic_start"] for site in proposal["tga_sites"]),
                   proposal["terminal_stop_genomic_start"])
            groups[pid][key].append(proposal)
            owners[key].add(pid)
    for anchor in anchors:
        pid, contig = anchor["id"], anchor["contig"]
        loci = groups[pid]
        outcome = "no_qualified_conditional_ORF"
        if len(loci) > 1 or any(len(owners[key]) > 1 for key in loci):
            outcome = "ambiguous_conditional_locus_or_anchor"
        elif loci:
            outcome = "unique_conditional_locus"
            (_, _, stops, terminal), alternatives = next(iter(loci.items()))
            alternatives.sort(key=lambda r: (r["start"], r["end"]))
            result["hypotheses"].append(dict(
                contig=contig, strand=anchor["strand"], catalytic_anchor=pid,
                tga_genomic_starts=list(stops), terminal_stop_genomic_start=terminal,
                coding_start_status="multiple_compatible_starts" if len(alternatives) > 1 else "one_compatible_start",
                coding_alternatives=alternatives, anchor_protein_sha256=anchor["protein_sha256"],
                anchor_coordinates=[anchor["start"], anchor["end"]], anchor_mapped_sites=anchor["mapped_sites"],
                component_reference_id=model["reference_id"], component_reference_accession=model["reference_accession"],
                translation_status="conditional_internal_TGA_no_U_assigned", SECIS="not_assessed",
                Sec_machinery="not_assessed", Sec_incorporation="unresolved", component_membership="unresolved",
                catalytic_identity="unconfirmed", query_activity="not_tested",
                overlapping_called_CDS=_overlaps(genes, contig, alternatives, len(contigs[contig])),
                competing_explanations=["ordinary_termination", "sequence_or_assembly_error", "other_coding_frame"]))
        result["searches"].append(dict(catalytic_anchor=pid, status=outcome, locus_count=len(loci),
                                       candidate_intervals=sorted([r["start"], r["end"]] for rs in loci.values() for r in rs)))
    return result


def audit(manifest, contigs, root=ROOT):
    """Search native sequence only when the reviewed split anchor has no called partner."""
    result = dict(rule=RULE, status="not_applicable", hypotheses=[], searches=[])
    if manifest.get("status") != "complete" or manifest.get("source_format") != "fasta" or not manifest.get("annotation"):
        return result
    root = Path(root)
    annotation = Path(manifest["annotation"])
    run_path = annotation / "run_manifest.json"
    run = load(run_path)
    expected = run.get("output_sha256", {})
    if (not {CANDIDATES, COMPONENTS} <= expected.keys()
            or not run.get("split_center_model_manifest_sha256") or not run.get("normalized_proteins_sha256")):
        return result
    eligible = {("required_component_not_detected", "unresolved"),
                ("component_ligand_evidence_incomplete", "competitive")}
    all_rows = _rows(annotation, CANDIDATES, expected)
    rows = [r for r in all_rows if (r["center_status"], r["pairing_status"]) in eligible]
    if not all_rows:
        result["status"] = "no_missing_or_incomplete_split_component"
        return result
    model_dir = root / "references/processed/split_centers"
    model_manifest_path = model_dir / "manifest.json"
    models_path = model_dir / "models.json"
    model_manifest = load(model_manifest_path)
    policy_path = root / "config/split_centers.json"
    if model_manifest["status"] != "complete":
        raise ValueError("Split translation model incomplete")
    models = load(models_path)
    if len(models) != 2 or sum(s["reference_residue"] == "U" for s in models[1]["sites"]) != 1:
        raise ValueError("Split translation requires one reviewed U site in the short component")
    component_rows = defaultdict(list)
    for r in _rows(annotation, COMPONENTS, expected):
        component_rows[r["catalytic_anchor"], r["center_id"]].append(r)
    genes_path = annotation / "genes.tsv"
    if sha(genes_path) != expected.get("genes.tsv"):
        raise ValueError("Split translation coordinates are not bound to the source annotation")
    genes = _load_genes(annotation)
    proteins_path = annotation / "proteins.faa"
    if sha(proteins_path) != run.get("normalized_proteins_sha256"):
        raise ValueError("Split translation protein FASTA is not bound to the source annotation")
    proteins = _load_proteins(proteins_path)
    tables = run.get("translation_table_used") or {}
    anchored = []
    counts = Counter(r["catalytic_anchor"] for r in rows)
    for row in rows:
        pid = row["catalytic_anchor"]
        if counts[pid] != 1:
            continue
        gene = _bound_anchor(row, component_rows[pid, row["center_id"]], genes, proteins, models, contigs, tables)
        if gene:
            anchored.append(gene)
    relay_anchors = {}
    for row in all_rows:
        pid = row["catalytic_anchor"]
        if pid not in relay_anchors:
            gene = _bound_anchor(row, component_rows[pid, row["center_id"]], genes, proteins, models, contigs, tables,
                                 bind_short=False)
            if gene:
                relay_anchors[pid] = gene
    if not anchored and not relay_anchors:
        result["status"] = "no_complete_bound_split_anchor" if rows else "no_missing_or_incomplete_split_component"
        return result
    files = {"annotation_manifest": run_path, "split_candidates": annotation / CANDIDATES,
             "split_components": annotation / COMPONENTS, "annotation_genes": genes_path,
             "annotation_proteins": proteins_path, "split_model_manifest": model_manifest_path,
             "split_models": models_path, "split_policy": policy_path, "component_audit_script": Path(__file__),
             "component_scan_script": Path(__file__).with_name("split_component_scan.py"),
             "tblastn": resolve_tool(root, "tblastn")}
    pins = {name: sha(path) for name, path in files.items()}
    if relay_anchors:
        result["upstream_relay_translation"] = _upstream_relay(list(relay_anchors.values()), contigs, genes, root)
    results = scan(anchored, contigs, models[1], model_manifest["parameters"], root) if anchored else {}
    if pins != {name: sha(path) for name, path in files.items()}:
        raise ValueError("Bound split translation inputs changed during search")
    if not anchored:
        result.update(status="no_complete_bound_split_anchor" if rows else "no_missing_or_incomplete_split_component",
                      evidence_sha256=pins)
        return result

    def locus(proposal):
        return (proposal["tga_genomic_start"], proposal["terminal_stop_genomic_start"])

    groups = defaultdict(lambda: defaultdict(list))
    owners = defaultdict(set)
    for gene in anchored:
        pid, contig = gene["id"], gene["contig"]
        fragment = gene["called_component_fragment"]
        for proposal in results.get(pid, []):
            if fragment and proposal["tga_genomic_start"] != fragment["terminal_TGA_start"]:
                continue
            key = locus(proposal)
            groups[pid][key].append(proposal)
            owners[contig, gene["strand"], key].add(pid)

    result.update(status="evaluated_conditional_split_component", evidence_sha256=pins)
    for gene in anchored:
        pid, contig = gene["id"], gene["contig"]
        loci = groups[pid]
        outcome = "no_qualified_conditional_ORF"
        if len(loci) > 1 or any(len(owners[contig, gene["strand"], key]) > 1 for key in loci):
            outcome = "ambiguous_conditional_locus_or_anchor"
        elif loci:
            outcome = "unique_conditional_locus"
            (tga, terminal), alternatives = next(iter(loci.items()))
            alternatives.sort(key=lambda r: (r["start"], r["end"]))
            overlap = _overlaps(genes, contig, alternatives, len(contigs[contig]))
            result["hypotheses"].append(dict(contig=contig, strand=gene["strand"], catalytic_anchor=pid,
                tga_genomic_start=tga, terminal_stop_genomic_start=terminal,
                coding_start_status="multiple_compatible_starts" if len(alternatives) > 1 else "one_compatible_start",
                coding_alternatives=alternatives, called_component_fragment=gene["called_component_fragment"],
                anchor_protein_sha256=gene["protein_sha256"],
                anchor_coordinates=[gene["start"], gene["end"]], anchor_mapped_sites=gene["mapped_sites"],
                component_reference_id=models[1]["reference_id"], component_reference_accession=models[1]["reference_accession"],
                translation_status="conditional_internal_TGA_no_U_assigned", SECIS="not_assessed",
                Sec_machinery="not_assessed", Sec_incorporation="unresolved", component_membership="unresolved",
                catalytic_identity="unconfirmed", query_activity="not_tested", overlapping_called_CDS=overlap,
                competing_explanations=["ordinary_termination", "sequence_or_assembly_error", "other_coding_frame"]))
        result["searches"].append(dict(catalytic_anchor=pid, status=outcome, locus_count=len(loci),
                                      candidate_intervals=sorted([r["start"], r["end"]] for rows in loci.values() for r in rows)))
    return result
