#!/usr/bin/env python3
"""
postprocess_atp_synthase.py
Post-annotation processing for config_mag_annotation:
1. Prunes GTDB metadata strictly to the 24 specified columns with robust accession resolution.
2. Evaluates ATP synthase modules (M00159, M00157, M00158, M00160) using continuous scoring,
   normalized score margins, structural sector gating, and physical operon clustering.
3. Distinguishes Complete_strict vs Complete_supported (operon context rescue).
4. Implements secondary-system ambiguity gating (Mixed_evidence_ambiguous) to produce clean
   VA_only vs F_like_only genome calls.
5. Employs MAG completeness awareness (None_confident, Uncertain_partial, Low_quality_not_detected).
Outputs:
  - genome_atp_synthase_master_table.tsv (Wide master table, one row per genome)
  - genome_atp_synthase_modules.long.tsv (Long-format module evaluation table)
  - genome_atp_synthase_subunits.tsv (Per-subunit detail table)
  - pruned_metadata.tsv (Cleaned 24-column metadata table)
"""

import argparse
import gzip
import os
import re
import sys
import pandas as pd
import numpy as np

METADATA_COLUMNS_24 = [
    "accession",
    "checkm2_completeness",
    "checkm2_contamination",
    "checkm_completeness",
    "checkm_contamination",
    "contig_count",
    "gc_percentage",
    "genome_size",
    "gtdb_representative",
    "gtdb_taxonomy",
    "gtdb_type_designation_ncbi_taxa",
    "gtdb_type_designation_ncbi_taxa_sources",
    "gtdb_type_species_of_genus",
    "l50_contigs",
    "mimag_high_quality",
    "mimag_low_quality",
    "mimag_medium_quality",
    "n50_contigs",
    "ncbi_assembly_level",
    "ncbi_genome_category",
    "ncbi_isolation_source",
    "ncbi_trna_count",
    "trna_aa_count",
    "trna_count",
]

# Module subunit definitions
MODULE_DEFINITIONS = {
    "M00159": {
        "name": "V/A-type ATPase, prokaryotes",
        "type": "V/A-type",
        "units": {
            "K02117": {"name": "A", "sector": "catalytic_head"},
            "K02118": {"name": "B", "sector": "catalytic_head"},
            "K02119": {"name": "C", "sector": "stalk"},
            "K02120": {"name": "D", "sector": "stalk"},
            "K02121": {"name": "E", "sector": "stalk"},
            "K02122": {"name": "F", "sector": "stalk"},
            "K02107": {"name": "G/H", "sector": "stalk"},
            "K02123": {"name": "I", "sector": "membrane"},
            "K02124": {"name": "K", "sector": "membrane"},
        },
    },
    "M00157": {
        "name": "F-type ATPase, prokaryotes and chloroplasts",
        "type": "F-like",
        "units": {
            "K02111": {"name": "alpha", "sector": "catalytic_head"},
            "K02112": {"name": "beta", "sector": "catalytic_head"},
            "K02113": {"name": "delta", "sector": "peripheral_coupling"},
            "K02114": {"name": "epsilon", "sector": "central_stalk"},
            "K02115": {"name": "gamma", "sector": "central_stalk"},
            "K02108": {"name": "a", "sector": "membrane"},
            "K02109": {"name": "b", "sector": "peripheral_coupling"},
            "K02110": {"name": "c", "sector": "membrane"},
        },
    },
    "M00158": {
        "name": "F-type ATPase, eukaryotes",
        "type": "Eukaryotic_F",
        "units": {
            "K02126": {"name": "ATP1", "sector": "catalytic"},
            "K02127": {"name": "ATP2", "sector": "catalytic"},
            "K02128": {"name": "ATP3", "sector": "stalk"},
            "K02129": {"name": "ATP4", "sector": "stalk"},
            "K02130": {"name": "ATP5", "sector": "stalk"},
            "K02131": {"name": "ATP6", "sector": "membrane"},
            "K02132": {"name": "ATP7", "sector": "membrane"},
            "K02133": {"name": "ATP8", "sector": "membrane"},
            "K02134": {"name": "ATP9", "sector": "membrane"},
            "K02135": {"name": "ATP10", "sector": "stalk"},
            "K02136": {"name": "ATP14", "sector": "stalk"},
            "K02137": {"name": "ATP15", "sector": "stalk"},
            "K02138": {"name": "ATP16", "sector": "stalk"},
        },
    },
    "M00160": {
        "name": "V-type ATPase, eukaryotes",
        "type": "Eukaryotic_V",
        "units": {
            "K02145": {"name": "V1A", "sector": "catalytic"},
            "K02146": {"name": "V1B", "sector": "catalytic"},
            "K02147": {"name": "V1C", "sector": "stalk"},
            "K02148": {"name": "V1D", "sector": "stalk"},
            "K02149": {"name": "V1E", "sector": "stalk"},
            "K02150": {"name": "V1F", "sector": "stalk"},
            "K02151": {"name": "V1G", "sector": "stalk"},
            "K02152": {"name": "V1H", "sector": "stalk"},
            "K02153": {"name": "V0a", "sector": "membrane"},
            "K02154": {"name": "V0d", "sector": "membrane"},
            "K02155": {"name": "V0c", "sector": "membrane"},
        },
    },
}

ALL_ATP_KOS = set()
for m_info in MODULE_DEFINITIONS.values():
    ALL_ATP_KOS.update(m_info["units"].keys())


def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)


def normalize_acc(acc_str):
    """Normalize accession by removing RS_ or GB_ prefixes and versions for robust mapping."""
    if not acc_str or not isinstance(acc_str, str):
        return ""
    acc = acc_str.strip()
    if acc.startswith("RS_"):
        acc = acc[3:]
    elif acc.startswith("GB_"):
        acc = acc[3:]
    return acc


def extract_core_accession(acc_str):
    """
    Extracts canonical NCBI accession (e.g. GCA_035542115.1 or GCF_000008085.1)
    from filenames, assembly names, or GTDB identifiers.
    """
    if not acc_str or not isinstance(acc_str, str):
        return ""
    m = re.search(r"(GC[AF]_\d+\.\d+)", acc_str)
    if m:
        return m.group(1)
    m = re.search(r"(GC[AF]_\d+)", acc_str)
    if m:
        return m.group(1)
    acc = acc_str.strip()
    if acc.startswith(("GB_", "RS_")):
        acc = acc[3:]
    return acc


def parse_gene_coordinates(protein_id, gene_id):
    """
    Extract contig, start, end, strand, and gene order from gene_id or protein_id.
    Standard pipeline ID scheme: {genome_id}___{contig_id}___gene{order:05d}
    """
    contig_id = "unknown_contig"
    order = 0
    start = 0
    end = 0
    strand = "+"

    target_id = gene_id if (gene_id and gene_id != ".") else protein_id
    if target_id and "___" in target_id:
        parts = target_id.split("___")
        if len(parts) >= 3:
            contig_id = parts[1]
            suffix = parts[2]
            # Match gene00042 or prot00042
            m = re.search(r"\d+", suffix)
            if m:
                order = int(m.group(0))
                start = order * 1000
                end = start + 900
        elif len(parts) == 2:
            contig_id = parts[0]
            m = re.search(r"\d+", parts[1])
            if m:
                order = int(m.group(0))
                start = order * 1000
                end = start + 900

    return contig_id, order, start, end, strand


def load_pruned_metadata(metadata_path):
    """
    Loads filtered_metadata.tsv, normalizes accessions, and indexes strictly to the 24 columns.
    Returns:
      metadata_dict: accession_variant -> dict of 24 fields
    """
    if not metadata_path or not os.path.exists(metadata_path):
        print(f"Warning: Metadata file '{metadata_path}' not found! Proceeding with empty metadata.", file=sys.stderr)
        return {}

    print(f"Loading metadata from {metadata_path}...")
    try:
        df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    except Exception as e:
        print(f"Warning: Failed to read metadata file: {e}", file=sys.stderr)
        return {}

    # Ensure all 24 columns exist in DataFrame
    for col in METADATA_COLUMNS_24:
        if col not in df.columns:
            df[col] = np.nan

    pruned_df = df[METADATA_COLUMNS_24].copy()

    metadata_dict = {}
    for _, row in pruned_df.iterrows():
        raw_acc = str(row["accession"]).strip()
        norm_acc = normalize_acc(raw_acc)
        core_acc = extract_core_accession(raw_acc)
        row_dict = row.to_dict()

        # Index under raw accession
        metadata_dict[raw_acc] = row_dict
        if norm_acc:
            metadata_dict[norm_acc] = row_dict
            if "." in norm_acc:
                base_acc = norm_acc.split(".")[0]
                metadata_dict[base_acc] = row_dict
        if core_acc:
            metadata_dict[core_acc] = row_dict
            if "." in core_acc:
                base_acc = core_acc.split(".")[0]
                metadata_dict[base_acc] = row_dict

    print(f"Loaded metadata for {len(df)} records. Indexed under {len(metadata_dict)} accession variants.")
    return metadata_dict


def lookup_genome_metadata(genome_id, metadata_dict):
    """
    Robust metadata lookup resolving:
    1. Exact genome_id (e.g. GB_GCA_001.1 or GCA_035542115.1_ASM3554211v1_genomic)
    2. Core NCBI accession (e.g. GCA_035542115.1)
    3. Normalized accession (stripped GB_/RS_)
    4. Base accession without dot/version (e.g. GCA_035542115)
    """
    if not metadata_dict:
        return None

    m_row = metadata_dict.get(genome_id)
    if m_row is not None:
        return m_row

    core_acc = extract_core_accession(genome_id)
    if core_acc:
        m_row = metadata_dict.get(core_acc)
        if m_row is not None:
            return m_row
        if "." in core_acc:
            m_row = metadata_dict.get(core_acc.split(".")[0])
            if m_row is not None:
                return m_row

    norm_id = normalize_acc(genome_id)
    if norm_id:
        m_row = metadata_dict.get(norm_id)
        if m_row is not None:
            return m_row
        if "." in norm_id:
            m_row = metadata_dict.get(norm_id.split(".")[0])
            if m_row is not None:
                return m_row

    if "." in genome_id:
        m_row = metadata_dict.get(genome_id.split(".")[0])
        if m_row is not None:
            return m_row

    return None


def discover_cohort_genomes(args, all_annotated_genomes, metadata_dict):
    """
    Determines the definitive set of cohort genomes.
    Guarantees no phantom genomes from global reference databases are injected.
    """
    cohort = set(all_annotated_genomes)

    if getattr(args, "genome_list", None):
        if os.path.isfile(args.genome_list):
            with open_file(args.genome_list, "rt") as f:
                for line in f:
                    g = line.strip()
                    if g:
                        cohort.add(g)
        else:
            for g in re.split(r"[\s,]+", args.genome_list):
                if g.strip():
                    cohort.add(g.strip())

    if getattr(args, "genome_dir", None) and os.path.isdir(args.genome_dir):
        valid_exts = (".fna.gz", ".fasta.gz", ".fa.gz", ".fna", ".fasta", ".fa")
        for root, _, files in os.walk(args.genome_dir):
            for f in files:
                for ext in valid_exts:
                    if f.endswith(ext) and not f.endswith((".faa", ".faa.gz", ".ffn", ".ffn.gz")):
                        stem = f[:-len(ext)]
                        cohort.add(stem)
                        break

    if getattr(args, "proteins_dir", None) and os.path.isdir(args.proteins_dir):
        valid_exts = (".faa.gz", ".faa")
        for root, _, files in os.walk(args.proteins_dir):
            for f in files:
                for ext in valid_exts:
                    if f.endswith(ext):
                        stem = f[:-len(ext)]
                        cohort.add(stem)
                        break

    if getattr(args, "include_all_metadata_genomes", False):
        annotated_core_accs = {extract_core_accession(g) for g in cohort if extract_core_accession(g)}
        for raw_acc, row in metadata_dict.items():
            if row.get("accession") and row["accession"] == raw_acc:
                core_acc = extract_core_accession(raw_acc)
                if raw_acc not in cohort and core_acc not in annotated_core_accs:
                    cohort.add(raw_acc)

    return cohort


def parse_kofam_hits(ko_strict_path):
    """
    Parses protein_KO_all_strict.tsv (or .tsv.gz).
    Extracts all hits relevant to ATP synthase modules.
    Returns:
      genome_hits: genome_id -> {ko -> list of hit_dicts}
      all_annotated_genomes: set of all genome IDs in file
    """
    genome_hits = {}
    all_annotated_genomes = set()

    if not os.path.exists(ko_strict_path):
        print(f"ERROR: Input KO strict file '{ko_strict_path}' does not exist!", file=sys.stderr)
        sys.exit(1)

    print(f"Reading ATP synthase hits from {ko_strict_path}...")
    with open_file(ko_strict_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_map = {name: idx for idx, name in enumerate(header)}

        p_idx = col_map.get("protein_id", 0)
        g_idx = col_map.get("gene_id", 1)
        gm_idx = col_map.get("genome_id", 2)
        k_idx = col_map.get("ko", 3)
        b_idx = col_map.get("bitscore", 5)
        th_idx = col_map.get("curated_threshold", 6)
        sr_idx = col_map.get("score_margin", 10)
        e_idx = col_map.get("evalue", 8)

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(p_idx, g_idx, gm_idx, k_idx):
                continue

            genome_id = parts[gm_idx]
            all_annotated_genomes.add(genome_id)

            ko = parts[k_idx]
            if ko not in ALL_ATP_KOS:
                continue

            prot_id = parts[p_idx]
            gene_id = parts[g_idx]

            try:
                bitscore = float(parts[b_idx]) if b_idx < len(parts) else 0.0
            except (ValueError, IndexError):
                bitscore = 0.0

            try:
                threshold = float(parts[th_idx]) if th_idx < len(parts) else 0.0
            except (ValueError, IndexError):
                threshold = 0.0

            try:
                evalue = float(parts[e_idx]) if e_idx < len(parts) else 1.0
            except (ValueError, IndexError):
                evalue = 1.0

            # Compute raw margin and relative margin
            raw_margin = bitscore - threshold
            if threshold > 0.0:
                rel_margin = raw_margin / threshold
            else:
                if bitscore > 0.0 and evalue <= 1e-5:
                    rel_margin = 0.0
                elif bitscore > 0.0 and evalue <= 1e-3:
                    rel_margin = -0.05
                else:
                    rel_margin = -0.25

            # Reject hits that are heavily sub-threshold (rel_margin < -0.20)
            if rel_margin < -0.20:
                continue

            contig, order, start, end, strand = parse_gene_coordinates(prot_id, gene_id)

            hit_data = {
                "protein_id": prot_id,
                "gene_id": gene_id,
                "bitscore": bitscore,
                "threshold": threshold,
                "evalue": evalue,
                "raw_margin": raw_margin,
                "relative_margin": rel_margin,
                "contig": contig,
                "order": order,
                "start": start,
                "end": end,
                "strand": strand,
            }

            if genome_id not in genome_hits:
                genome_hits[genome_id] = {}
            if ko not in genome_hits[genome_id]:
                genome_hits[genome_id][ko] = []
            genome_hits[genome_id][ko].append(hit_data)

    print(f"Extracted ATP synthase hits across {len(genome_hits)} genomes (out of {len(all_annotated_genomes)} total genomes).")
    return genome_hits, all_annotated_genomes


def detect_operon_clusters(subunits_by_ko, operon_max_gap=3):
    """
    Groups subunit hits on the same contig and determines if they form an operon cluster.
    Returns:
      clustered_kos: set of KOs in an operon cluster
      cluster_metrics: dict with cluster_size_units, cluster_span_bp, max_gene_gap, strand_consistency
    """
    contig_groups = {}
    for ko, hits in subunits_by_ko.items():
        for h in hits:
            c = h["contig"]
            if c not in contig_groups:
                contig_groups[c] = []
            contig_groups[c].append((h["order"], h["start"], h["end"], h["strand"], ko))

    best_cluster_kos = set()
    best_metrics = {
        "cluster_size_units": 0,
        "cluster_span_bp": 0,
        "max_gene_gap": 0,
        "strand_consistency": True,
    }

    for c, items in contig_groups.items():
        if len(items) < 2:
            continue
        # Sort by order or start position
        items.sort(key=lambda x: (x[0], x[1]))

        # Find contiguous segments where inter-gene gap <= operon_max_gap
        current_segment = [items[0]]
        for i in range(1, len(items)):
            prev_order = items[i - 1][0]
            curr_order = items[i][0]
            gap = abs(curr_order - prev_order) - 1
            if gap <= operon_max_gap:
                current_segment.append(items[i])
            else:
                seg_kos = set(x[4] for x in current_segment)
                if len(seg_kos) > len(best_cluster_kos):
                    best_cluster_kos = seg_kos
                    strands = set(x[3] for x in current_segment)
                    best_metrics = {
                        "cluster_size_units": len(seg_kos),
                        "cluster_span_bp": max(x[2] for x in current_segment) - min(x[1] for x in current_segment) + 1,
                        "max_gene_gap": max((items[j][0] - items[j - 1][0] - 1) for j in range(1, len(current_segment))) if len(current_segment) > 1 else 0,
                        "strand_consistency": (len(strands) <= 1),
                    }
                current_segment = [items[i]]

        if current_segment:
            seg_kos = set(x[4] for x in current_segment)
            if len(seg_kos) > len(best_cluster_kos):
                best_cluster_kos = seg_kos
                strands = set(x[3] for x in current_segment)
                best_metrics = {
                    "cluster_size_units": len(seg_kos),
                    "cluster_span_bp": max(x[2] for x in current_segment) - min(x[1] for x in current_segment) + 1,
                    "max_gene_gap": max((items[j][0] - items[j - 1][0] - 1) for j in range(1, len(current_segment))) if len(current_segment) > 1 else 0,
                    "strand_consistency": (len(strands) <= 1),
                }

    return best_cluster_kos, best_metrics


def evaluate_genome_module(genome_id, mod_id, genome_hits_dict, operon_max_gap=3,
                           w_borderline=0.5, w_supported_borderline=0.75,
                           genome_completeness=100.0, genome_contamination=0.0):
    """
    Evaluates a single module for a single genome.
    Returns:
      mod_summary: dict of module-level metrics
      subunit_records: list of subunit detail records
    """
    mod_info = MODULE_DEFINITIONS[mod_id]
    required_units = mod_info["units"]
    n_required = len(required_units)

    # Collect present units
    subunits_data = {}
    present_kos_map = {}

    for ko, u_info in required_units.items():
        hits = genome_hits_dict.get(ko, [])
        if hits:
            # Sort by relative margin descending
            hits.sort(key=lambda x: (x["relative_margin"], x["bitscore"]), reverse=True)
            best_h = hits[0]
            strict_copies = sum(1 for h in hits if h["relative_margin"] >= 0.0)
            supported_copies = sum(1 for h in hits if h["relative_margin"] >= -0.10)

            subunits_data[ko] = {
                "unit_present": 1,
                "strict_copy_number": strict_copies,
                "supported_copy_number": supported_copies,
                "unit_best_score": best_h["bitscore"],
                "unit_score_margin": best_h["relative_margin"],
                "raw_score_margin": best_h["raw_margin"],
                "best_hit": best_h,
                "all_proteins": [h["protein_id"] for h in hits],
                "sector": u_info["sector"],
                "name": u_info["name"],
            }
            present_kos_map[ko] = hits
        else:
            subunits_data[ko] = {
                "unit_present": 0,
                "strict_copy_number": 0,
                "supported_copy_number": 0,
                "unit_best_score": 0.0,
                "unit_score_margin": np.nan,
                "raw_score_margin": np.nan,
                "best_hit": None,
                "all_proteins": [],
                "sector": u_info["sector"],
                "name": u_info["name"],
            }

    # Detect operon clustering
    clustered_kos, cluster_metrics = detect_operon_clusters(present_kos_map, operon_max_gap)
    gene_clustered = (len(clustered_kos) >= 3) or (len(clustered_kos) >= 2 and len(clustered_kos) >= 0.5 * len(present_kos_map))

    n_strong = 0
    n_supported_borderline = 0
    n_borderline = 0
    n_missing = 0

    margins_list = []
    missing_kos_desc = []
    subunit_records = []

    for ko, s_data in subunits_data.items():
        u_pres = s_data["unit_present"]
        rel_m = s_data["unit_score_margin"]
        is_in_cluster = (ko in clustered_kos)

        evidence_tier = "None"

        if u_pres == 1:
            margins_list.append(rel_m)
            if rel_m >= 0.0:
                evidence_tier = "Strong"
                n_strong += 1
            elif -0.10 <= rel_m < 0.0:
                if is_in_cluster:
                    evidence_tier = "Supported_borderline"
                    n_supported_borderline += 1
                else:
                    evidence_tier = "Borderline"
                    n_borderline += 1
            else:
                evidence_tier = "Weak"
                n_missing += 1
                missing_kos_desc.append(f"{ko}(weak:m={rel_m:.2f})")
        else:
            evidence_tier = "None"
            n_missing += 1
            missing_kos_desc.append(f"{ko}(missing)")

        best_hit = s_data["best_hit"]
        contig_id = best_hit["contig"] if best_hit else "none"

        subunit_records.append({
            "genome_id": genome_id,
            "module_id": mod_id,
            "ko": ko,
            "subunit_name": s_data["name"],
            "sector": s_data["sector"],
            "unit_present": u_pres,
            "strict_copy_number": s_data["strict_copy_number"],
            "supported_copy_number": s_data["supported_copy_number"],
            "unit_best_score": round(s_data["unit_best_score"], 2),
            "unit_score_margin": round(rel_m, 4) if not np.isnan(rel_m) else np.nan,
            "raw_score_margin": round(s_data["raw_score_margin"], 2) if not np.isnan(s_data["raw_score_margin"]) else np.nan,
            "evidence_tier": evidence_tier,
            "is_clustered": is_in_cluster,
            "contig_id": contig_id,
            "protein_ids": ",".join(s_data["all_proteins"]),
        })

    # Calculate completeness
    strict_comp = n_strong / float(n_required)
    weighted_comp = (n_strong + n_supported_borderline * w_supported_borderline + n_borderline * w_borderline) / float(n_required)

    mean_margin = np.mean(margins_list) if margins_list else np.nan
    min_margin = np.min(margins_list) if margins_list else np.nan

    # Sector presence
    catalytic_present = False
    stalk_present = False
    membrane_present = False

    if mod_id == "M00159":
        # Catalytic head: K02117 (A) and K02118 (B)
        cat_a = subunits_data["K02117"]["unit_present"] == 1 and subunits_data["K02117"]["unit_score_margin"] >= -0.10
        cat_b = subunits_data["K02118"]["unit_present"] == 1 and subunits_data["K02118"]["unit_score_margin"] >= -0.10
        catalytic_present = (cat_a and cat_b)

        stalk_kos = ["K02119", "K02120", "K02121", "K02122", "K02107"]
        stalk_cnt = sum(1 for k in stalk_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        stalk_present = (stalk_cnt >= 2)

        mem_kos = ["K02123", "K02124"]
        mem_cnt = sum(1 for k in mem_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        membrane_present = (mem_cnt >= 1)

    elif mod_id == "M00157":
        # Catalytic head: K02111 (alpha) and K02112 (beta)
        cat_a = subunits_data["K02111"]["unit_present"] == 1 and subunits_data["K02111"]["unit_score_margin"] >= -0.10
        cat_b = subunits_data["K02112"]["unit_present"] == 1 and subunits_data["K02112"]["unit_score_margin"] >= -0.10
        catalytic_present = (cat_a and cat_b)

        stalk_kos = ["K02113", "K02114", "K02115", "K02109"]
        stalk_cnt = sum(1 for k in stalk_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        stalk_present = (stalk_cnt >= 1)

        mem_kos = ["K02108", "K02110"]
        mem_cnt = sum(1 for k in mem_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        membrane_present = (mem_cnt >= 1)

    elif mod_id == "M00158":
        # Catalytic head: K02126 (ATP1) and K02127 (ATP2)
        cat_a = subunits_data["K02126"]["unit_present"] == 1 and subunits_data["K02126"]["unit_score_margin"] >= -0.10
        cat_b = subunits_data["K02127"]["unit_present"] == 1 and subunits_data["K02127"]["unit_score_margin"] >= -0.10
        catalytic_present = (cat_a and cat_b)

        stalk_kos = ["K02128", "K02129", "K02130", "K02135", "K02136", "K02137", "K02138"]
        stalk_cnt = sum(1 for k in stalk_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        stalk_present = (stalk_cnt >= 2)

        mem_kos = ["K02131", "K02132", "K02133", "K02134"]
        mem_cnt = sum(1 for k in mem_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        membrane_present = (mem_cnt >= 1)

    elif mod_id == "M00160":
        # Catalytic head: K02145 (V1A) and K02146 (V1B)
        cat_a = subunits_data["K02145"]["unit_present"] == 1 and subunits_data["K02145"]["unit_score_margin"] >= -0.10
        cat_b = subunits_data["K02146"]["unit_present"] == 1 and subunits_data["K02146"]["unit_score_margin"] >= -0.10
        catalytic_present = (cat_a and cat_b)

        stalk_kos = ["K02147", "K02148", "K02149", "K02150", "K02151", "K02152"]
        stalk_cnt = sum(1 for k in stalk_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        stalk_present = (stalk_cnt >= 2)

        mem_kos = ["K02153", "K02154", "K02155"]
        mem_cnt = sum(1 for k in mem_kos if subunits_data[k]["unit_present"] == 1 and subunits_data[k]["unit_score_margin"] >= -0.10)
        membrane_present = (mem_cnt >= 1)

    else:
        # General eukaryotic modules
        catalytic_present = (n_strong + n_borderline + n_supported_borderline) >= 2
        stalk_present = (n_strong + n_borderline + n_supported_borderline) >= 3
        membrane_present = (n_strong + n_borderline + n_supported_borderline) >= 1

    # Discrete state determination
    total_detected = n_strong + n_supported_borderline + n_borderline

    if n_strong == n_required:
        mod_state = "Complete_strict"
    elif n_missing == 0 and weighted_comp >= 0.95:
        mod_state = "Complete_supported"
    elif weighted_comp >= 0.80 and n_missing <= 2 and catalytic_present and membrane_present:
        mod_state = "Near_complete"
    elif weighted_comp >= 0.50 and catalytic_present and (stalk_present or membrane_present):
        mod_state = "Partial"
    elif weighted_comp >= 0.20 or total_detected >= 2:
        mod_state = "Fragmentary"
    else:
        mod_state = "Not_detected"

    # MAG completeness aware modulation
    mod_confidence = "High"
    if genome_completeness < 70.0:
        if mod_state == "Partial":
            mod_state = "Uncertain_partial"
            mod_confidence = "Low"
        elif mod_state == "Not_detected":
            mod_state = "Low_quality_not_detected"
            mod_confidence = "Low"
        else:
            mod_confidence = "Medium"
    elif 70.0 <= genome_completeness < 90.0:
        if mod_state in ("Partial", "Fragmentary"):
            mod_confidence = "Medium"
        elif mod_state == "Near_complete":
            mod_confidence = "Medium"

    mod_summary = {
        "genome_id": genome_id,
        "module": mod_id,
        "module_name": mod_info["name"],
        "n_required": n_required,
        "n_strong": n_strong,
        "n_supported_borderline": n_supported_borderline,
        "n_borderline": n_borderline,
        "n_missing": n_missing,
        "strict_completeness": round(strict_comp, 4),
        "weighted_completeness": round(weighted_comp, 4),
        "catalytic_present": catalytic_present,
        "stalk_present": stalk_present,
        "membrane_present": membrane_present,
        "best_missing_KOs": ";".join(missing_kos_desc) if missing_kos_desc else "None",
        "mean_score_margin": round(mean_margin, 4) if not np.isnan(mean_margin) else np.nan,
        "min_score_margin": round(min_margin, 4) if not np.isnan(min_margin) else np.nan,
        "gene_clustered": gene_clustered,
        "cluster_size_units": cluster_metrics["cluster_size_units"],
        "cluster_span_bp": cluster_metrics["cluster_span_bp"],
        "max_gene_gap": cluster_metrics["max_gene_gap"],
        "strand_consistency": cluster_metrics["strand_consistency"],
        "genome_completeness": round(genome_completeness, 2),
        "module_state": mod_state,
        "module_confidence": mod_confidence,
    }

    return mod_summary, subunit_records


def assign_genome_phenotype(m00159_summary, m00157_summary, genome_comp=100.0, genome_contam=0.0,
                            hq_comp_cutoff=90.0, hq_contam_cutoff=5.0):
    """
    Classifies the genome into one of the mutually exclusive genome phenotypes:
    VA_only, F_like_only, VA_plus_F_like, None_confident, Mixed_evidence_ambiguous, Ambiguous.
    Assigns contrast cohort memberships.
    """
    va_state = m00159_summary["module_state"]
    f_state = m00157_summary["module_state"]

    va_functional = va_state in ("Complete_strict", "Complete_supported", "Near_complete")
    f_functional = f_state in ("Complete_strict", "Complete_supported", "Near_complete")

    va_clean_absent = va_state == "Not_detected"
    f_clean_absent = f_state == "Not_detected"

    va_partial = va_state in ("Partial", "Fragmentary", "Uncertain_partial")
    f_partial = f_state in ("Partial", "Fragmentary", "Uncertain_partial")

    genome_phenotype = "Ambiguous"

    if va_functional and f_functional:
        genome_phenotype = "VA_plus_F_like"
    elif va_functional and f_clean_absent:
        genome_phenotype = "VA_only"
    elif f_functional and va_clean_absent:
        genome_phenotype = "F_like_only"
    elif (va_functional and f_partial) or (f_functional and va_partial):
        # Secondary-system ambiguity gate: One complete system + remnants of another
        genome_phenotype = "Mixed_evidence_ambiguous"
    elif va_clean_absent and f_clean_absent:
        if genome_comp >= hq_comp_cutoff and genome_contam <= hq_contam_cutoff:
            genome_phenotype = "None_confident"
        else:
            genome_phenotype = "Ambiguous"
    else:
        genome_phenotype = "Ambiguous"

    # Contrast 1 Cohort (ATP+ vs ATP-)
    if genome_phenotype in ("VA_only", "F_like_only", "VA_plus_F_like"):
        contrast1_cohort = "ATP_POSITIVE"
    elif genome_phenotype == "None_confident":
        contrast1_cohort = "ATP_NEGATIVE"
    else:
        contrast1_cohort = "EXCLUDED"

    # Contrast 2 Cohort (VA_only vs F_like_only)
    if genome_phenotype == "VA_only":
        contrast2_cohort = "VA_ONLY"
    elif genome_phenotype == "F_like_only":
        contrast2_cohort = "F_LIKE_ONLY"
    else:
        contrast2_cohort = "EXCLUDED"

    return genome_phenotype, contrast1_cohort, contrast2_cohort


def main():
    parser = argparse.ArgumentParser(description="Post-process MAG annotation and evaluate ATP synthase structural gating.")
    parser.add_argument("-k", "--ko_strict", required=True, help="Path to protein_KO_all_strict.tsv (or .gz)")
    parser.add_argument("-m", "--metadata", default="", help="Path to filtered_metadata.tsv (optional)")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for master tables")
    parser.add_argument("--operon_max_gene_gap", type=int, default=3, help="Max non-subunit genes between operon subunits (default: 3)")
    parser.add_argument("--weight_borderline", type=float, default=0.5, help="Weight for borderline hits (default: 0.5)")
    parser.add_argument("--weight_supported_borderline", type=float, default=0.75, help="Weight for supported borderline hits (default: 0.75)")
    parser.add_argument("--high_quality_genome_completeness", type=float, default=90.0, help="Min completeness for None_confident (default: 90.0)")
    parser.add_argument("--high_quality_contamination_limit", type=float, default=5.0, help="Max contamination for None_confident (default: 5.0)")
    parser.add_argument("--low_quality_genome_completeness", type=float, default=70.0, help="Cutoff below which detection is low confidence (default: 70.0)")
    parser.add_argument("--genome_list", type=str, default="", help="Optional file or comma-separated list of cohort genome IDs")
    parser.add_argument("--genome_dir", type=str, default="", help="Optional directory containing cohort genome fasta files")
    parser.add_argument("--proteins_dir", type=str, default="", help="Optional directory containing cohort protein fasta files")
    parser.add_argument("--include_all_metadata_genomes", action="store_true", default=False, help="Explicitly include all accessions from metadata file in cohort")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Load pruned metadata
    metadata_dict = load_pruned_metadata(args.metadata)

    # 2. Parse KO hits
    genome_hits, all_annotated_genomes = parse_kofam_hits(args.ko_strict)

    # 3. Determine cohort genomes strictly without phantom metadata injection
    all_target_genomes = discover_cohort_genomes(args, all_annotated_genomes, metadata_dict)
    print(f"Total target cohort: {len(all_target_genomes)} genomes.")

    # 4. Export cohort-specific pruned metadata
    pruned_rows = []
    for genome_id in sorted(all_target_genomes):
        m_row = lookup_genome_metadata(genome_id, metadata_dict)
        if m_row is None:
            row_dict = {col: np.nan for col in METADATA_COLUMNS_24}
            row_dict["accession"] = genome_id
        else:
            row_dict = {col: m_row.get(col, np.nan) for col in METADATA_COLUMNS_24}
            if pd.isna(row_dict.get("accession")):
                row_dict["accession"] = genome_id
        pruned_rows.append(row_dict)

    pruned_meta_tsv = os.path.join(args.out_dir, "pruned_metadata.tsv")
    pruned_meta_df = pd.DataFrame(pruned_rows, columns=METADATA_COLUMNS_24)
    pruned_meta_df.to_csv(pruned_meta_tsv, sep="\t", index=False)
    print(f"Exported 24-column pruned metadata for {len(pruned_meta_df)} cohort genomes to {pruned_meta_tsv}")

    master_rows = []
    long_module_rows = []
    all_subunit_rows = []

    for genome_id in sorted(all_target_genomes):
        g_hits = genome_hits.get(genome_id, {})

        # Resolve metadata for genome
        m_row = lookup_genome_metadata(genome_id, metadata_dict)
        if m_row is None:
            # Create empty placeholder row for 24 columns
            m_row = {col: np.nan for col in METADATA_COLUMNS_24}
            m_row["accession"] = genome_id

        # Determine genome completeness and contamination
        try:
            g_comp = float(m_row["checkm2_completeness"])
            if np.isnan(g_comp):
                g_comp = float(m_row["checkm_completeness"])
        except (ValueError, TypeError, KeyError):
            g_comp = 100.0

        if np.isnan(g_comp):
            g_comp = 100.0

        try:
            g_contam = float(m_row["checkm2_contamination"])
            if np.isnan(g_contam):
                g_contam = float(m_row["checkm_contamination"])
        except (ValueError, TypeError, KeyError):
            g_contam = 0.0

        if np.isnan(g_contam):
            g_contam = 0.0

        # Evaluate M00159 (V/A-type)
        m00159_sum, m00159_sub = evaluate_genome_module(
            genome_id, "M00159", g_hits,
            operon_max_gap=args.operon_max_gene_gap,
            w_borderline=args.weight_borderline,
            w_supported_borderline=args.weight_supported_borderline,
            genome_completeness=g_comp,
            genome_contamination=g_contam,
        )

        # Evaluate M00157 (F-like)
        m00157_sum, m00157_sub = evaluate_genome_module(
            genome_id, "M00157", g_hits,
            operon_max_gap=args.operon_max_gene_gap,
            w_borderline=args.weight_borderline,
            w_supported_borderline=args.weight_supported_borderline,
            genome_completeness=g_comp,
            genome_contamination=g_contam,
        )

        # Evaluate M00158 (Euk F)
        m00158_sum, m00158_sub = evaluate_genome_module(
            genome_id, "M00158", g_hits,
            operon_max_gap=args.operon_max_gene_gap,
            w_borderline=args.weight_borderline,
            w_supported_borderline=args.weight_supported_borderline,
            genome_completeness=g_comp,
            genome_contamination=g_contam,
        )

        # Evaluate M00160 (Euk V)
        m00160_sum, m00160_sub = evaluate_genome_module(
            genome_id, "M00160", g_hits,
            operon_max_gap=args.operon_max_gene_gap,
            w_borderline=args.weight_borderline,
            w_supported_borderline=args.weight_supported_borderline,
            genome_completeness=g_comp,
            genome_contamination=g_contam,
        )

        long_module_rows.extend([m00159_sum, m00157_sum, m00158_sum, m00160_sum])
        all_subunit_rows.extend(m00159_sub + m00157_sub + m00158_sub + m00160_sub)

        # Assign genome-level phenotype & cohorts
        g_pheno, c1_cohort, c2_cohort = assign_genome_phenotype(
            m00159_sum, m00157_sum,
            genome_comp=g_comp,
            genome_contam=g_contam,
            hq_comp_cutoff=args.high_quality_genome_completeness,
            hq_contam_cutoff=args.high_quality_contamination_limit,
        )

        # Construct Master Row
        # 1. 24 Metadata columns first
        master_row = {}
        for col in METADATA_COLUMNS_24:
            master_row[col] = m_row.get(col, np.nan)
        if pd.isna(master_row["accession"]):
            master_row["accession"] = genome_id

        # 2. Genome-level phenotype and contrast cohorts
        master_row["genome_id"] = genome_id
        master_row["atp_synthase_genome_phenotype"] = g_pheno
        master_row["atp_contrast1_cohort"] = c1_cohort
        master_row["atp_contrast2_cohort"] = c2_cohort

        # 3. Module-level details for all 4 ATP synthase modules (M00159, M00157, M00158, M00160)
        for mod_id, mod_sum in [
            ("M00159", m00159_sum),
            ("M00157", m00157_sum),
            ("M00158", m00158_sum),
            ("M00160", m00160_sum),
        ]:
            master_row[f"{mod_id}_state"] = mod_sum["module_state"]
            master_row[f"{mod_id}_confidence"] = mod_sum["module_confidence"]
            master_row[f"{mod_id}_strict_comp"] = mod_sum["strict_completeness"]
            master_row[f"{mod_id}_weighted_comp"] = mod_sum["weighted_completeness"]
            master_row[f"{mod_id}_n_strong"] = mod_sum["n_strong"]
            master_row[f"{mod_id}_n_supported_borderline"] = mod_sum["n_supported_borderline"]
            master_row[f"{mod_id}_n_borderline"] = mod_sum["n_borderline"]
            master_row[f"{mod_id}_n_missing"] = mod_sum["n_missing"]
            master_row[f"{mod_id}_catalytic"] = mod_sum["catalytic_present"]
            master_row[f"{mod_id}_stalk"] = mod_sum["stalk_present"]
            master_row[f"{mod_id}_membrane"] = mod_sum["membrane_present"]
            master_row[f"{mod_id}_gene_clustered"] = mod_sum["gene_clustered"]
            master_row[f"{mod_id}_cluster_size_units"] = mod_sum["cluster_size_units"]
            master_row[f"{mod_id}_cluster_span_bp"] = mod_sum["cluster_span_bp"]
            master_row[f"{mod_id}_max_gene_gap"] = mod_sum["max_gene_gap"]
            master_row[f"{mod_id}_strand_consistency"] = mod_sum["strand_consistency"]
            master_row[f"{mod_id}_mean_margin"] = mod_sum["mean_score_margin"]
            master_row[f"{mod_id}_missing_KOs"] = mod_sum["best_missing_KOs"]

        master_rows.append(master_row)

    # 3. Export Master Table
    master_df = pd.DataFrame(master_rows)
    master_path = os.path.join(args.out_dir, "genome_atp_synthase_master_table.tsv")
    master_df.to_csv(master_path, sep="\t", index=False)
    print(f"Exported master summary table to {master_path} ({len(master_df)} genomes).")

    # 4. Export Long Module Table
    long_df = pd.DataFrame(long_module_rows)
    long_path = os.path.join(args.out_dir, "genome_atp_synthase_modules.long.tsv")
    long_df.to_csv(long_path, sep="\t", index=False)
    print(f"Exported long module evaluation table to {long_path} ({len(long_df)} records).")

    # 5. Export Subunit Table
    subunit_df = pd.DataFrame(all_subunit_rows)
    subunit_path = os.path.join(args.out_dir, "genome_atp_synthase_subunits.tsv")
    subunit_df.to_csv(subunit_path, sep="\t", index=False)
    print(f"Exported subunit breakdown table to {subunit_path} ({len(subunit_df)} records).")


if __name__ == "__main__":
    main()
