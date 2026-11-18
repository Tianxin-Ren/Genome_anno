#!/usr/bin/env python3
"""Lightweight integration helpers for the unified BioB annotation route."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
from collections import defaultdict


ALLOWED_STATES = {
    "supported_specific_function",
    "supported_broad_function",
    "ambiguous_between_candidates",
    "partial_or_truncated_evidence",
    "conflicting_evidence",
    "unresolved",
}


def read_fasta_ids(path):
    ids = []
    if not path or not os.path.exists(path):
        return ids
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_fasta_records(path):
    records = {}
    current = None
    chunks = []
    if not path or not os.path.exists(path):
        return records
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if current is not None:
        records[current] = "".join(chunks)
    return records


def read_contig_lengths(path):
    return {name: len(seq) for name, seq in read_fasta_records(path).items()}


def attrs_to_dict(text):
    attrs = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def parse_cds_gff(path):
    features = {}
    contig_order = defaultdict(int)
    if not path or not os.path.exists(path):
        return features
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = attrs_to_dict(parts[8])
            pid = attrs.get("ID") or attrs.get("protein_id") or attrs.get("Name")
            if not pid:
                continue
            contig_order[parts[0]] += 1
            partial = "partial=true" in parts[8].lower() or attrs.get("partial", "").lower() == "true"
            features[pid] = {
                "feature_id": pid,
                "contig": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "score": parts[5],
                "gene_index_on_contig": contig_order[parts[0]],
                "partial": partial,
            }
    return features


def add_evidence(evidence, protein_id, source, evidence_class, value, detail=None, score=None):
    if not protein_id or not value:
        return
    evidence[protein_id].append({
        "source": source,
        "evidence_class": evidence_class,
        "value": value,
        "detail": detail or "",
        "score": score or "",
    })


def parse_mmseqs(path, evidence):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            detail = f"identity={parts[2]};aln_len={parts[3]};evalue={parts[10]}"
            add_evidence(evidence, parts[0], "mmseqs2_swissprot", "full_length_sequence_or_ortholog", parts[1], detail, parts[11])


def parse_eggnog(path, evidence):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            pid = parts[0]
            desc = parts[7] if len(parts) > 7 else ""
            name = parts[8] if len(parts) > 8 else ""
            kos = parts[11] if len(parts) > 11 else ""
            pfams = parts[20] if len(parts) > 20 else ""
            for ko in re.split(r"[,\s]+", kos.replace("ko:", "")):
                if ko.startswith("K"):
                    add_evidence(evidence, pid, "eggnog_mapper", "full_length_sequence_or_ortholog", ko, desc)
            if name and name != "-":
                add_evidence(evidence, pid, "eggnog_mapper", "full_length_sequence_or_ortholog", name, desc)
            for pfam in re.split(r"[,;\s]+", pfams):
                if pfam and pfam != "-":
                    add_evidence(evidence, pid, "eggnog_mapper_pfam", "domain_and_active_site", pfam, desc)


def parse_kofam(path, evidence):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if not parts:
                continue
            above = parts[0] == "*"
            offset = 1 if above else 0
            if len(parts) <= offset + 1:
                continue
            pid = parts[offset]
            ko = parts[offset + 1]
            detail = "above_threshold=true" if above else "above_threshold=false"
            add_evidence(evidence, pid, "kofamscan", "profile_hmm", ko, detail)


def parse_interpro(path, evidence):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            pid = parts[0]
            analysis = parts[3]
            accession = parts[4]
            desc = parts[5] if len(parts) > 5 else ""
            value = f"{analysis}:{accession}"
            add_evidence(evidence, pid, "interproscan", "domain_and_active_site", value, desc)


def parse_hmmer_tbl(path, evidence, source):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            target = parts[0]
            query = parts[2]
            add_evidence(evidence, target, source, "profile_hmm", query, f"evalue={parts[4]}", parts[5])


def parse_simple_tsv(path, evidence, source, evidence_class):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pid = row.get("feature_id") or row.get("protein_id") or row.get("Protein_ID") or row.get("query")
            value = row.get("prediction") or row.get("value") or row.get("signal_peptide") or row.get("topology")
            if pid and value:
                add_evidence(evidence, pid, source, evidence_class, value, json.dumps(row, sort_keys=True))


def canonical_function(items, gene_id):
    """Normalize heterogeneous hits into a namespace-preserving primary ID."""
    xrefs = defaultdict(set)
    products = []
    for item in items:
        value = str(item.get("value", ""))
        source = item.get("source", "")
        detail = str(item.get("detail", ""))
        for ec in re.findall(r"\bEC[: ]?(\d+\.\d+\.\d+\.[\d-]+)\b", value + " " + detail):
            xrefs["EC"].add(ec)
        for ko in re.findall(r"\bK\d{5}\b", value):
            xrefs["KO"].add(ko)
        if source == "interproscan" and ":" in value:
            analysis, accession = value.split(":", 1)
            analysis = analysis.upper()
            if analysis == "INTERPRO":
                ns = "InterPro"
            elif analysis in {"PFAM", "PFAM-A"}:
                ns = "Pfam"
            else:
                ns = analysis
            xrefs[ns].add(accession)
        if source == "eggnog_mapper_pfam":
            xrefs["Pfam"].add(value)
        if source == "eggnog_mapper" and value not in xrefs["KO"] and not value.startswith("K"):
            if value and value not in {"-", "unknown"}:
                products.append(value)
        if source == "mmseqs2_swissprot":
            match = re.search(r"(?:sp|tr)\|([^|]+)\|", value)
            accession = match.group(1) if match else value
            if accession:
                xrefs["UniProtKB"].add(accession)
            name = value.split("|")[-1] if "|" in value else ""
            if name:
                products.append(name.replace("_", " "))
        if source == "biob_hmm":
            xrefs["BioB_HMM"].add(value)

    # The primary namespace is a deliberate policy, not an attempted accession conversion.
    if len(xrefs["KO"]) == 1:
        namespace, identifier, level, reason = "KEGG_KO", f"KO:{next(iter(xrefs['KO']))}", "specific_function", "single KOfam/orthology KO"
    elif len(xrefs["KO"]) > 1:
        identifier = "AMBIGUOUS:KO:" + ",".join(sorted(xrefs["KO"]))
        namespace, level, reason = "AMBIGUOUS", "ambiguous", "multiple incompatible KO candidates"
    elif xrefs["EC"] and len(xrefs["EC"]) == 1:
        identifier = f"EC:{next(iter(xrefs['EC']))}"
        namespace, level, reason = "EC", "enzyme_reaction", "EC recovered from convergent evidence"
    elif xrefs["InterPro"]:
        identifier = f"InterPro:{sorted(xrefs['InterPro'])[0]}"
        namespace, level, reason = "InterPro", "protein_family", "InterPro family/domain evidence"
    elif xrefs["Pfam"]:
        identifier = f"Pfam:{sorted(xrefs['Pfam'])[0]}"
        namespace, level, reason = "Pfam", "domain_family", "Pfam domain evidence"
    elif xrefs["BioB_HMM"]:
        identifier = f"BioB_HMM:{sorted(xrefs['BioB_HMM'])[0]}"
        namespace, level, reason = "BioB_HMM", "protein_family", "curated profile-HMM evidence"
    elif xrefs["UniProtKB"]:
        identifier = f"UniProtKB:{sorted(xrefs['UniProtKB'])[0]}"
        namespace, level, reason = "UniProtKB", "homology_only", "Swiss-Prot homology retained without KO/domain convergence"
    else:
        identifier = f"UNKNOWN:{gene_id}"
        namespace, level, reason = "UNKNOWN", "unknown", "no normalized functional accession"
    return identifier, namespace, level, (products[0] if products else "hypothetical protein"), reason, {
        key: sorted(values) for key, values in sorted(xrefs.items()) if values
    }


def stable_gene_id(sample, pid, feature, fallback_index):
    contig = feature.get("contig", "unknown")
    index = feature.get("gene_index_on_contig") or fallback_index
    return f"{sample}|{contig}|gene_{int(index):06d}"


def consensus_state(items, feature):
    if feature and feature.get("partial"):
        return "partial_or_truncated_evidence"
    if not items:
        return "unresolved"
    values_by_class = defaultdict(set)
    for item in items:
        values_by_class[item["evidence_class"]].add(item["value"])
    ko_values = {item["value"] for item in items if re.match(r"^K\d{5}$", item["value"])}
    if len(ko_values) > 1:
        return "ambiguous_between_candidates"
    if any(item["source"] in {"kofamscan", "biob_hmm"} for item in items):
        return "supported_specific_function"
    if values_by_class:
        return "supported_broad_function"
    return "unresolved"


def confidence_levels(items, feature):
    """Return independent structural (S0-S3) and functional (F0-F4) levels."""
    if not feature:
        structural = "S0"
    elif feature.get("partial"):
        structural = "S1"
    elif feature.get("contig_edge"):
        structural = "S1"
    elif feature.get("protein_length", 0) < 30:
        structural = "S1"
    else:
        structural = "S3"
    sources = {x["source"] for x in items}
    classes = {x["evidence_class"] for x in items}
    if not items:
        functional = "F0"
    elif "uniprot_swissprot" in sources or "biob_hmm" in sources:
        functional = "F4"
    elif "kofamscan" in sources and ("eggnog_mapper" in sources or "interproscan" in sources):
        functional = "F3"
    elif "kofamscan" in sources or "eggnog_mapper" in sources:
        functional = "F2"
    elif classes:
        functional = "F1"
    else:
        functional = "F0"
    return structural, functional


def write_json(path, payload):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_unresolved_faa(records, unresolved_ids, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        for pid in unresolved_ids:
            seq = records.get(pid)
            if seq:
                handle.write(f">{pid}\n")
                for i in range(0, len(seq), 80):
                    handle.write(seq[i:i + 80] + "\n")


def run_consensus(args):
    protein_ids = read_fasta_ids(args.faa)
    records = read_fasta_records(args.faa)
    features = parse_cds_gff(args.gff)
    contig_lengths = read_contig_lengths(args.fasta) if getattr(args, "fasta", "") else {}
    for feature in features.values():
        feature["protein_length"] = 0
        feature["contig_length"] = contig_lengths.get(feature["contig"], 0)
        feature["contig_edge"] = bool(feature["contig_length"] and
                                       (feature["start"] <= 5 or feature["end"] >= feature["contig_length"] - 4))
    for pid, seq in records.items():
        if pid in features:
            features[pid]["protein_length"] = len(seq)
    evidence = defaultdict(list)
    parse_mmseqs(args.mmseqs, evidence)
    parse_eggnog(args.eggnog, evidence)
    parse_kofam(args.kofam, evidence)
    parse_interpro(args.interpro, evidence)
    parse_simple_tsv(args.signalp, evidence, "signalp", "localization_and_topology")
    parse_simple_tsv(args.topology, evidence, "deeptmhmm", "localization_and_topology")
    parse_hmmer_tbl(args.biob_hmm, evidence, "biob_hmm")

    os.makedirs(os.path.dirname(args.candidates), exist_ok=True)
    os.makedirs(os.path.dirname(args.inventory), exist_ok=True)

    unresolved = []
    state_counts = defaultdict(int)
    source_counts = defaultdict(int)
    ko_inventory = defaultdict(list)
    domain_inventory = defaultdict(list)
    conflicts = []

    with open(args.candidates, "w", newline="") as cand_handle, open(args.consensus, "w") as con_handle:
        cand_writer = csv.writer(cand_handle, delimiter="\t")
        cand_writer.writerow(["feature_id", "evidence_source", "evidence_class", "candidate", "detail", "score"])
        for pid in protein_ids:
            items = evidence.get(pid, [])
            gene_id = stable_gene_id(args.sample, pid, features.get(pid, {}), protein_ids.index(pid) + 1)
            canonical_id, canonical_ns, function_level, product, selection_reason, xrefs = canonical_function(items, gene_id)
            for item in items:
                cand_writer.writerow([pid, item["source"], item["evidence_class"], item["value"], item["detail"], item["score"]])
                source_counts[item["source"]] += 1
                if re.match(r"^K\d{5}$", item["value"]):
                    ko_inventory[item["value"]].append(pid)
                if item["evidence_class"] == "domain_and_active_site":
                    domain_inventory[item["value"]].append(pid)
            state = consensus_state(items, features.get(pid))
            structural_confidence, functional_confidence = confidence_levels(items, features.get(pid))
            state_counts[state] += 1
            if state == "unresolved":
                unresolved.append(pid)
            if state in {"ambiguous_between_candidates", "conflicting_evidence"}:
                conflicts.append((pid, state, "multiple incompatible candidates require lineage/site/neighborhood review"))
            payload = {
                "feature_id": pid,
                "gene_id": gene_id,
                "state": state,
                "evidence_classes": sorted({item["evidence_class"] for item in items}),
                "candidate_count": len(items),
                "candidates": items[:25],
                "location": features.get(pid, {}),
                "structural_confidence": structural_confidence,
                "functional_confidence": functional_confidence,
                "canonical_function_id": canonical_id,
                "canonical_function_namespace": canonical_ns,
                "function_level": function_level,
                "selection_reason": selection_reason,
                "all_xrefs": xrefs,
            }
            con_handle.write(json.dumps(payload, sort_keys=True) + "\n")

    with open(args.conflicts, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["feature_id", "conflict_type", "description"])
        for row in conflicts:
            writer.writerow(row)

    write_unresolved_faa(records, unresolved, args.unresolved)

    # Stable, analysis-facing tables. Raw source files remain untouched above.
    os.makedirs(os.path.dirname(args.genes), exist_ok=True)
    with open(args.genes, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", "gene_id", "feature_id", "contig", "start", "end", "strand",
                         "gene_type", "protein_length", "partial", "contig_edge", "product",
                         "canonical_function_id", "canonical_function_namespace", "function_level",
                         "KO", "EC", "InterPro", "Pfam", "UniProtKB", "all_xrefs",
                         "structural_confidence", "functional_confidence", "conflict_flag",
                         "selection_reason", "evidence_summary"])
        for fallback_index, pid in enumerate(protein_ids, 1):
            feature = features.get(pid, {})
            items = evidence.get(pid, [])
            gene_id = stable_gene_id(args.sample, pid, feature, fallback_index)
            gene_state = consensus_state(items, feature)
            s_conf, f_conf = confidence_levels(items, feature)
            canonical_id, canonical_ns, function_level, product, selection_reason, xrefs = canonical_function(items, gene_id)
            writer.writerow([args.sample, gene_id, pid, feature.get("contig", ""), feature.get("start", ""),
                             feature.get("end", ""), feature.get("strand", ""), "CDS",
                             len(records.get(pid, "")), str(bool(feature.get("partial"))).lower(),
                             str(bool(feature.get("contig_edge"))).lower(), product,
                             canonical_id, canonical_ns, function_level,
                             ";".join(xrefs.get("KO", [])), ";".join(xrefs.get("EC", [])),
                             ";".join(xrefs.get("InterPro", [])), ";".join(xrefs.get("Pfam", [])),
                             ";".join(xrefs.get("UniProtKB", [])), json.dumps(xrefs, sort_keys=True), s_conf, f_conf,
                             str(gene_state == "ambiguous_between_candidates").lower(),
                             selection_reason,
                             ";".join(sorted({x["source"] for x in items})) or "none"])
    with open(args.evidence_table, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", "gene_id", "feature_id", "source", "evidence_type", "hit_id", "detail", "score", "evidence_rank"])
        for fallback_index, pid in enumerate(protein_ids, 1):
            gene_id = stable_gene_id(args.sample, pid, features.get(pid, {}), fallback_index)
            for item in evidence.get(pid, []):
                writer.writerow([args.sample, gene_id, pid, item["source"], item["evidence_class"], item["value"], item["detail"], item["score"], 1])
    for path, header in [(args.modules, ["genome_id", "module_id", "module_name", "required_components", "detected_components", "missing_components", "module_completeness", "module_confidence"]),
                         (args.unknown_clusters, ["cluster_id", "member_count", "genome_count", "representative_gene", "domain_evidence", "priority"])]:
        with open(path, "w", newline="") as handle:
            csv.writer(handle, delimiter="\t").writerow(header)
    os.makedirs(os.path.dirname(args.report), exist_ok=True)
    with open(args.report, "w") as handle:
        handle.write(f"<html><head><meta charset='utf-8'><title>{args.sample} annotation report</title></head><body>")
        handle.write(f"<h1>{args.sample}</h1><p>Proteins: {len(protein_ids)}; unresolved: {len(unresolved)}</p>")
        handle.write("<p>Confidence is reported independently as structural S0-S3 and functional F0-F4.</p></body></html>\n")

    inventory = {
        "sample": args.sample,
        "protein_count": len(protein_ids),
        "state_counts": dict(sorted(state_counts.items())),
        "evidence_source_counts": dict(sorted(source_counts.items())),
        "ko_count": len(ko_inventory),
        "domain_signature_count": len(domain_inventory),
        "unresolved_count": len(unresolved),
    }
    write_json(args.inventory, inventory)


def trnascan_rows(path):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("-") or line.lower().startswith("sequence"):
                continue
            parts = line.split()
            if len(parts) < 6 or not parts[2].isdigit() or not parts[3].isdigit():
                continue
            contig, idx, start, end, trna_type, codon = parts[:6]
            yield contig, idx, int(start), int(end), trna_type, codon


def rfam_rows(path):
    if not path or not os.path.exists(path):
        return
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            target_name = parts[0]
            accession = parts[1]
            contig = parts[2]
            try:
                start = int(parts[7])
                end = int(parts[8])
            except ValueError:
                continue
            yield contig, accession, target_name, start, end


def aragorn_rows(path):
    if not path or not os.path.exists(path):
        return
    current = None
    idx = 0
    coord_re = re.compile(r"([tc]RNA|tmRNA)-?([A-Za-z]*)\s+c?\[?(\d+),(\d+)\]?")
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                idx = 0
                continue
            if not current:
                continue
            match = coord_re.search(line)
            if not match:
                continue
            idx += 1
            kind, name, start, end = match.groups()
            yield current, idx, kind, name or kind, int(start), int(end)


def run_rna_gff(args):
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write("##gff-version 3\n")
        for contig, idx, start, end, trna_type, codon in trnascan_rows(args.trna) or []:
            strand = "+" if start <= end else "-"
            lo, hi = sorted((start, end))
            attrs = f"ID={contig}:tRNA_{idx};product=tRNA-{trna_type};anticodon={codon};source=tRNAscan-SE"
            handle.write(f"{contig}\ttRNAscan-SE\ttRNA\t{lo}\t{hi}\t.\t{strand}\t.\t{attrs}\n")
        for contig, accession, query, start, end in rfam_rows(args.rfam) or []:
            strand = "+" if start <= end else "-"
            lo, hi = sorted((start, end))
            attrs = f"ID={contig}:Rfam_{accession}_{lo}_{hi};Name={query};Dbxref=Rfam:{accession};source=Infernal_Rfam"
            handle.write(f"{contig}\tInfernal_Rfam\tncRNA\t{lo}\t{hi}\t.\t{strand}\t.\t{attrs}\n")
        for contig, idx, kind, name, start, end in aragorn_rows(args.aragorn) or []:
            strand = "+" if start <= end else "-"
            lo, hi = sorted((start, end))
            feature_type = "tmRNA" if kind == "tmRNA" else "tRNA"
            attrs = f"ID={contig}:aragorn_{feature_type}_{idx};product={name};source=ARAGORN"
            handle.write(f"{contig}\tARAGORN\t{feature_type}\t{lo}\t{hi}\t.\t{strand}\t.\t{attrs}\n")


def run_cross_genome(args):
    """Create a compact multi-genome inventory without requiring a database."""
    import glob
    rows = []
    unknown = defaultdict(list)
    for path in sorted(glob.glob(os.path.join(args.root, "*", "tables", "genes.tsv"))):
        sample = os.path.basename(os.path.dirname(os.path.dirname(path)))
        faa_path = os.path.join(args.root, sample, "structural", "proteins.faa")
        protein_records = read_fasta_records(faa_path)
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            counts = defaultdict(int)
            total = 0
            for row in reader:
                total += 1
                counts[row.get("functional_confidence", "F0")] += 1
                if row.get("product") == "hypothetical protein":
                    seq = protein_records.get(row.get("feature_id", ""), "")
                    # Exact-sequence grouping is deliberately conservative and
                    # reproducible; MMseqs2 clustering can replace this later.
                    unknown[seq or row.get("gene_id", "")].append((sample, row.get("gene_id", "")))
            rows.append((sample, total, counts))
    os.makedirs(os.path.dirname(args.genomes), exist_ok=True)
    with open(args.genomes, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_id", "protein_count", "F4_count", "F3_count", "F2_count", "F1_count", "F0_count", "specific_function_rate", "family_or_better_rate"])
        for sample, total, counts in rows:
            denom = total or 1
            writer.writerow([sample, total, counts["F4"], counts["F3"], counts["F2"], counts["F1"], counts["F0"],
                             f"{(counts['F4'] + counts['F3']) / denom:.6f}", f"{(total - counts['F0']) / denom:.6f}"])
    with open(args.unknown, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cluster_id", "member_count", "genome_count", "representative_gene", "domain_evidence", "priority"])
        for i, (cluster_key, members) in enumerate(sorted(unknown.items()), 1):
            samples = {sample for sample, _ in members}
            representative = members[0][1]
            writer.writerow([f"UC{i:06d}", len(members), len(samples), representative, "none",
                             "high" if len(samples) > 1 else "normal"])
    with open(args.report, "w") as handle:
        handle.write("<html><head><meta charset='utf-8'><title>Unified annotation report</title></head><body>")
        handle.write(f"<h1>Unified annotation report</h1><p>Genomes: {len(rows)}</p>")
        handle.write("<p>Rates distinguish specific function (F3/F4) from family-or-better coverage (F1-F4).</p></body></html>\n")


def run_gene_evidence(args):
    features = parse_cds_gff(args.gff)
    ids = read_fasta_ids(args.faa)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as handle:
        if not ids:
            ids = sorted(features)
        for pid in ids:
            feature = features.get(pid, {})
            payload = {
                "sample": args.sample,
                "feature_id": pid,
                "tool": "pyrodigal",
                "mode": args.mode,
                "translation_table": args.translation_table,
                "confidence": "not_reported_by_cli_wrapper",
                "prediction_score": feature.get("score", ""),
                "partial": bool(feature.get("partial", False)),
                "location": feature,
            }
            handle.write(json.dumps(payload, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    rna = sub.add_parser("rna-gff")
    rna.add_argument("--trna", required=True)
    rna.add_argument("--rfam", required=True)
    rna.add_argument("--aragorn", required=True)
    rna.add_argument("--output", required=True)
    rna.set_defaults(func=run_rna_gff)

    gene = sub.add_parser("gene-evidence")
    gene.add_argument("--sample", required=True)
    gene.add_argument("--gff", required=True)
    gene.add_argument("--faa", required=True)
    gene.add_argument("--mode", required=True)
    gene.add_argument("--translation-table", required=True)
    gene.add_argument("--output", required=True)
    gene.set_defaults(func=run_gene_evidence)

    con = sub.add_parser("consensus")
    for name in ["sample", "fasta", "faa", "gff", "mmseqs", "eggnog", "kofam", "interpro", "signalp", "topology", "biob-hmm"]:
        con.add_argument(f"--{name}", required=True)
    for name in ["candidates", "consensus", "conflicts", "unresolved", "inventory",
                 "genes", "evidence-table", "modules", "unknown-clusters", "report"]:
        con.add_argument(f"--{name}", required=True)
    con.set_defaults(func=run_consensus)

    summary = sub.add_parser("cross-genome")
    summary.add_argument("--root", required=True)
    summary.add_argument("--genomes", required=True)
    summary.add_argument("--unknown", required=True)
    summary.add_argument("--report", required=True)
    summary.set_defaults(func=run_cross_genome)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
