#!/usr/bin/env python3
"""Create focused, comparable GenBank windows for Clinker visualization.

The analysis FASTA/GBK files remain unchanged. This helper only crops the
visualization inputs to module/anchor genes plus a bounded amount of context.
"""

import argparse
import csv
import glob
import os
import re
from copy import deepcopy

try:
    from Bio import SeqIO, SeqFeature
except ImportError:  # pragma: no cover - exercised only in a missing dependency env
    SeqIO = None
    SeqFeature = None


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-gbk", required=True)
    parser.add_argument("--topology-tsv", required=True)
    parser.add_argument("--refs-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--module-id", required=True)
    parser.add_argument("--target-kos", default="")
    parser.add_argument("--context-bp", type=int, default=10000)
    parser.add_argument("--flank-genes", type=int, default=25,
                        help="Number of CDS features retained on each side of the focus region")
    return parser.parse_args()


def read_topology_anchor_ids(topology_tsv):
    """Return target protein IDs used as left/right conserved-context anchors."""
    anchors = set()
    if not topology_tsv or not os.path.isfile(topology_tsv):
        return anchors
    with open(topology_tsv, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for column in ("left_matched_genes", "right_matched_genes"):
                anchors.update(pid.strip() for pid in row.get(column, "").split(",") if pid.strip())
            if anchors:
                break
    return anchors


def read_topology_target_ids(topology_tsv):
    """Return target CDS IDs between the conserved left and right anchors."""
    target_ids = set()
    if not topology_tsv or not os.path.isfile(topology_tsv):
        return target_ids
    with open(topology_tsv, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target_ids.update(
                pid.strip() for pid in row.get("intervening_genes", "").split(",")
                if pid.strip() and pid.strip() != "CROSS_CONTIG_SPLIT"
            )
            if target_ids:
                break
    return target_ids


def feature_ids(feature):
    ids = set()
    for key in ("locus_tag", "protein_id", "orig_protein_id", "gene_id"):
        for value in feature.qualifiers.get(key, []):
            ids.add(str(value).strip())
    return ids


def _id_keys(values):
    keys = set()
    for value in values:
        text = str(value).strip()
        if not text:
            continue
        keys.add(text)
        keys.add(text.split("___")[-1])
        keys.add(text.rsplit("|", 1)[-1])
    return keys


def feature_kos(feature):
    kos = set()
    for values in feature.qualifiers.values():
        for value in values:
            kos.update(re.findall(r"K\d{5}", str(value)))
    return kos


def feature_is_operon_core(feature):
    values = []
    for key in ("gene", "label", "note", "product"):
        values.extend(str(value) for value in feature.qualifiers.get(key, []))
    return any(re.search(r"\bOperon_", value) for value in values)


def _feature_bounds(feature):
    return int(feature.location.start), int(feature.location.end)


def _compact_cluster(features, max_gap=20000):
    """Keep the densest contiguous CDS cluster when markers are scattered."""
    ordered = sorted(features, key=lambda f: _feature_bounds(f)[0])
    if len(ordered) <= 1:
        return ordered
    clusters = [[ordered[0]]]
    for feature in ordered[1:]:
        prev_end = _feature_bounds(clusters[-1][-1])[1]
        if _feature_bounds(feature)[0] - prev_end <= max_gap:
            clusters[-1].append(feature)
        else:
            clusters.append([feature])
    # Prefer the cluster with the most genes, then the shortest span.
    return min(clusters, key=lambda c: (-len(c), _feature_bounds(c[-1])[1] - _feature_bounds(c[0])[0]))


def focus_feature_ids(record, module_kos, anchor_ids=(), target_ids=(), include_module_core=True):
    """Select CDS features marking the target anchors or reference module core."""
    selected_anchors = []
    selected_target = []
    selected_core = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        ids = feature_ids(feature)
        id_keys = _id_keys(ids)
        if id_keys.intersection(_id_keys(anchor_ids)):
            selected_anchors.append(feature)
        if id_keys.intersection(_id_keys(target_ids)):
            selected_target.append(feature)
        # For the target, module KOs are the precise target-gene definition.
        # Reference GBKs use Operon_ labels to retain the corresponding core
        # genes even when their KO qualifiers are absent.
        is_core = bool(feature_kos(feature).intersection(module_kos))
        if not target_ids and not anchor_ids:
            is_core = is_core or feature_is_operon_core(feature)
        if include_module_core and is_core:
            selected_core.append(feature)
    # Older target GBKs may lack KO qualifiers entirely. In that case retain
    # their explicit Operon_ core labels as a safe fallback.
    if not selected_core and (target_ids or anchor_ids):
        selected_core = [
            feature for feature in record.features
            if feature.type == "CDS" and feature_is_operon_core(feature)
        ]
    # The module core defines the visual focus when present. Including every
    # nearby topology anchor can make the target window span most of a long
    # locus, while reference windows remain centered on the compact core.
    if selected_core:
        return _compact_cluster(selected_core)
    if selected_target:
        return _compact_cluster(selected_target)
    elif selected_anchors:
        # For KO-free samples, use the outermost matched anchors to retain the
        # complete interval containing the target genes, then add 25 genes on
        # each side. A single median anchor would discard that interval.
        ordered = sorted(selected_anchors, key=lambda f: sum(_feature_bounds(f)) / 2)
        return [ordered[0], ordered[-1]] if len(ordered) > 1 else ordered
    return []


def crop_record(record, focus_features, context_bp, flank_genes=25):
    """Crop one GenBank record around focus genes without truncating CDSs.

    The full interval between the outer focus genes is retained. Thus a true
    target insertion or long region between conserved anchors remains visible.
    """
    if not focus_features:
        return record, 0, len(record.seq), False

    seq_len = len(record.seq)
    cds_features = [feature for feature in record.features if feature.type == "CDS"]
    cds_features.sort(key=lambda feature: _feature_bounds(feature)[0])
    focus_ids = {id(feature) for feature in focus_features}
    focus_indices = [i for i, feature in enumerate(cds_features) if id(feature) in focus_ids]
    if focus_indices:
        first = max(0, min(focus_indices) - max(0, int(flank_genes)))
        last = min(len(cds_features) - 1, max(focus_indices) + max(0, int(flank_genes)))
        focus_start = _feature_bounds(cds_features[first])[0]
        focus_end = _feature_bounds(cds_features[last])[1]
    else:
        focus_start = min(_feature_bounds(feature)[0] for feature in focus_features)
        focus_end = max(_feature_bounds(feature)[1] for feature in focus_features)
    # Keep the configured CDS flank count, while retaining a bp minimum for
    # sparse regions. The complete focus interval is always preserved.
    window_start = max(0, focus_start - max(0, int(context_bp)))
    window_end = min(seq_len, focus_end + max(0, int(context_bp)))

    # Extend to complete any CDS touched by a window edge, then keep the bounds
    # stable. This avoids emitting partial CDS annotations with full translations.
    changed = True
    while changed:
        changed = False
        for feature in cds_features:
            start, end = _feature_bounds(feature)
            if start < window_start < end:
                window_start = start
                changed = True
            if start < window_end < end:
                window_end = end
                changed = True

    if window_start <= 0 and window_end >= seq_len:
        return record, 0, seq_len, False

    cropped = deepcopy(record)
    cropped.seq = record.seq[window_start:window_end]
    new_features = []
    for feature in record.features:
        if feature.type != "CDS":
            continue
        start, end = _feature_bounds(feature)
        if start < window_start or end > window_end:
            continue
        updated = deepcopy(feature)
        updated.location = SeqFeature.FeatureLocation(
            start - window_start,
            end - window_start,
            strand=feature.location.strand,
        )
        new_features.append(updated)

    new_features.insert(
        0,
        SeqFeature.SeqFeature(
            SeqFeature.FeatureLocation(0, len(cropped.seq)),
            type="source",
            qualifiers={},
        ),
    )
    cropped.features = new_features
    cropped.description = f"{record.description} [focused {window_start}-{window_end}]"
    return cropped, window_start, window_end, True


def _write_focused(source_path, out_path, module_kos, anchor_ids, target_ids, context_bp, flank_genes):
    with open(source_path, encoding="utf-8") as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in {source_path}")
    focused_records = []
    reports = []
    for record in records:
        selected = focus_feature_ids(record, module_kos, anchor_ids, target_ids)
        focused, start, end, cropped = crop_record(record, selected, context_bp, flank_genes)
        focused_records.append(focused)
        reports.append((record.id, len(record.seq), end - start, len(selected), cropped))
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        SeqIO.write(focused_records, handle, "genbank")
    return reports


def prepare_inputs(target_gbk, topology_tsv, refs_dir, out_dir, module_id, target_kos, context_bp, flank_genes=25):
    if SeqIO is None or SeqFeature is None:
        raise RuntimeError("Biopython is required to prepare focused Clinker inputs")
    os.makedirs(out_dir, exist_ok=True)
    module_kos = {ko.strip() for ko in target_kos.split(",") if ko.strip()}
    anchor_ids = read_topology_anchor_ids(topology_tsv)
    target_ids = read_topology_target_ids(topology_tsv)
    reports = []

    # Use a stable target filename; the source GenBank basename can vary
    # between samples and previously caused the target track to be omitted.
    target_out = os.path.join(out_dir, "target_locus.focused.gbk")
    reports.extend(("target",) + report for report in _write_focused(
        target_gbk, target_out, module_kos, anchor_ids, target_ids, context_bp, flank_genes
    ))

    # Focus each reference on its ATPase/module core. Other GenBank files in
    # this directory are pair-specific target records and use the matching
    # row's anchor IDs where possible.
    reference_files = [
        path for path in sorted(glob.glob(os.path.join(refs_dir, "*_locus.gbk")))
        if not path.endswith("_target_locus.gbk")
    ]
    for source_path in reference_files:
        out_path = os.path.join(out_dir, os.path.basename(source_path).replace(".gbk", ".focused.gbk"))
        reports.extend(("reference",) + report for report in _write_focused(
            source_path, out_path, module_kos, (), (), context_bp, flank_genes
        ))

    pair_target_files = sorted(glob.glob(os.path.join(refs_dir, "*_target_locus.gbk")))
    for source_path in pair_target_files:
        out_path = os.path.join(out_dir, os.path.basename(source_path).replace(".gbk", ".focused.gbk"))
        reports.extend(("pair_target",) + report for report in _write_focused(
            source_path, out_path, module_kos, anchor_ids, target_ids, context_bp, flank_genes
        ))

    original_order = os.path.join(refs_dir, "ordered_clinker_files.txt")
    order_out = os.path.join(out_dir, "ordered_clinker_files.txt")
    ordered_focus = []
    if os.path.isfile(original_order):
        with open(original_order, encoding="utf-8") as handle:
            for line in handle:
                source = line.strip()
                if not source:
                    continue
                focused_path = os.path.join(out_dir, os.path.basename(source).replace(".gbk", ".focused.gbk"))
                if os.path.isfile(focused_path):
                    ordered_focus.append(focused_path)
    if not ordered_focus:
        ordered_focus = [
            os.path.join(out_dir, os.path.basename(path).replace(".gbk", ".focused.gbk"))
            for path in reference_files
        ]
    with open(order_out, "w", encoding="utf-8") as handle:
        handle.write("".join(f"{path}\n" for path in ordered_focus))

    for kind, record_id, original_bp, focused_bp, n_focus, cropped in reports:
        action = "cropped" if cropped else "kept full (no marker or window covers record)"
        print(
            f"{kind} {record_id}: {original_bp} -> {focused_bp} bp; "
            f"{n_focus} marker genes; {action}",
            flush=True,
        )
    print(
        f"Prepared Clinker windows for {module_id}: up to {flank_genes} CDS flanks per side "
        f"plus {context_bp} bp minimum context; "
        f"target genes={len(target_ids)}, target anchors={len(anchor_ids)}, module KOs={len(module_kos)}",
        flush=True,
    )
    return target_out


def main():
    args = parse_args()
    target_kos = args.target_kos
    target_out = prepare_inputs(
        args.target_gbk,
        args.topology_tsv,
        args.refs_dir,
        args.out_dir,
        args.module_id,
        target_kos,
        args.context_bp,
        args.flank_genes,
    )
    print(f"Focused target GenBank: {target_out}", flush=True)


if __name__ == "__main__":
    main()
