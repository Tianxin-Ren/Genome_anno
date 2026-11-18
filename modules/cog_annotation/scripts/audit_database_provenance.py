#!/usr/bin/env python3
"""
audit_database_provenance.py
Gate 0: Database Provenance & Profile Compatibility Audit for COG Annotation.
Audits compatibility between NCBI CDD Cog_LE PSSMs (via cddid.tbl) and
official COG2024 definitions (cog-24.def.tab, cog-24.fun.tab, bitscore_specific.txt).

Key functions:
1. Calculates SHA256 checksums of database assets:
   - Mandatory: cddid.tbl, Cog_LE files, cog-24.def.tab, cog-24.fun.tab
   - Recommended: bitscore_specific.txt
2. Performs set intersection analysis between COG2024 official definitions and Cog_LE PSSM registry.
3. Ingests profile-specific bit score thresholds from bitscore_specific.txt if provided.
4. Exports clean mapping table (cdd_pssm_mapping.tsv.gz) mapping pssm_id -> cog_id, pssm_length, description, specific_bitscore_threshold.
5. Generates official category definitions table (cog_category_definition.tsv) from cog-24.fun.tab.
6. Generates audit report with canonical method designation:
   "CDD Cog_LE-based COG annotation mapped to COG2024 definitions".
"""

import argparse
import gzip
import hashlib
import json
import os
import sys
import glob

# Fallback NCBI categories if cog-24.fun.tab is missing or unreadable
FALLBACK_CATEGORIES = [
    ("J", "Translation, ribosomal structure and biogenesis", "INFORMATION STORAGE AND PROCESSING"),
    ("A", "RNA processing and modification", "INFORMATION STORAGE AND PROCESSING"),
    ("K", "Transcription", "INFORMATION STORAGE AND PROCESSING"),
    ("L", "Replication, recombination and repair", "INFORMATION STORAGE AND PROCESSING"),
    ("B", "Chromatin structure and dynamics", "INFORMATION STORAGE AND PROCESSING"),
    ("D", "Cell cycle control, cell division, chromosome partitioning", "CELLULAR PROCESSES AND SIGNALING"),
    ("Y", "Nuclear structure", "CELLULAR PROCESSES AND SIGNALING"),
    ("V", "Defense mechanisms", "CELLULAR PROCESSES AND SIGNALING"),
    ("T", "Signal transduction mechanisms", "CELLULAR PROCESSES AND SIGNALING"),
    ("M", "Cell wall/membrane/envelope biogenesis", "CELLULAR PROCESSES AND SIGNALING"),
    ("N", "Cell motility", "CELLULAR PROCESSES AND SIGNALING"),
    ("Z", "Cytoskeleton", "CELLULAR PROCESSES AND SIGNALING"),
    ("W", "Extracellular structures", "CELLULAR PROCESSES AND SIGNALING"),
    ("U", "Intracellular trafficking, secretion, and vesicular transport", "CELLULAR PROCESSES AND SIGNALING"),
    ("O", "Posttranslational modification, protein turnover, chaperones", "CELLULAR PROCESSES AND SIGNALING"),
    ("X", "Mobilome: prophages, transposons", "CELLULAR PROCESSES AND SIGNALING"),
    ("C", "Energy production and conversion", "METABOLISM"),
    ("G", "Carbohydrate transport and metabolism", "METABOLISM"),
    ("E", "Amino acid transport and metabolism", "METABOLISM"),
    ("F", "Nucleotide transport and metabolism", "METABOLISM"),
    ("H", "Coenzyme transport and metabolism", "METABOLISM"),
    ("I", "Lipid transport and metabolism", "METABOLISM"),
    ("P", "Inorganic ion transport and metabolism", "METABOLISM"),
    ("Q", "Secondary metabolites biosynthesis, transport and catabolism", "METABOLISM"),
    ("R", "General function prediction only", "POORLY CHARACTERIZED"),
    ("S", "Function unknown", "POORLY CHARACTERIZED"),
]

def compute_sha256(filepath):
    """Compute SHA256 checksum of a file."""
    if not filepath or not os.path.isfile(filepath):
        return None
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def parse_cddid(cddid_path):
    """
    Parse cddid.tbl from NCBI CDD.
    Standard columns: PSSM-Id, Accession, Short Name, Description, PSSM-Length.
    Returns:
        pssm_dict: {pssm_id: {"cog_id": cog_id, "length": int, "name": str, "desc": str}}
        cog_to_pssm: {cog_id: [pssm_id, ...]}
    """
    pssm_dict = {}
    cog_to_pssm = {}
    
    if not os.path.isfile(cddid_path):
        raise FileNotFoundError(f"cddid.tbl not found at: {cddid_path}")

    opener = gzip.open if cddid_path.endswith(".gz") else open
    with opener(cddid_path, "rt", encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            pssm_id = parts[0].strip()
            acc = parts[1].strip()
            if not acc.upper().startswith("COG"):
                continue
            cog_id = acc.upper()
            short_name = parts[2].strip() if len(parts) > 2 else ""
            desc = parts[3].strip() if len(parts) > 3 else ""
            length = 0
            if len(parts) > 4 and parts[4].strip().isdigit():
                length = int(parts[4].strip())

            pssm_dict[pssm_id] = {
                "cog_id": cog_id,
                "length": length,
                "name": short_name,
                "desc": desc
            }
            cog_to_pssm.setdefault(cog_id, []).append(pssm_id)

    return pssm_dict, cog_to_pssm

def parse_cog_defs(cog_def_path):
    """
    Parse cog-24.def.tab or cog-20.def.tab.
    Columns: COG_id, functional_category, COG_name, gene_name, pathway, pubmed_id, pdb_id.
    """
    cog_defs = {}
    if not os.path.isfile(cog_def_path):
        raise FileNotFoundError(f"COG definition table not found at: {cog_def_path}")

    opener = gzip.open if cog_def_path.endswith(".gz") else open
    with opener(cog_def_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if not parts:
                continue
            cog_id = parts[0].strip().upper()
            if not cog_id.startswith("COG"):
                continue
            category = parts[1].strip() if len(parts) > 1 else "S"
            name = parts[2].strip() if len(parts) > 2 else ""
            gene = parts[3].strip() if len(parts) > 3 else ""
            pathway = parts[4].strip() if len(parts) > 4 else ""
            cog_defs[cog_id] = {
                "category": category,
                "name": name,
                "gene": gene,
                "pathway": pathway
            }
    return cog_defs

def parse_cog_fun(cog_fun_path):
    """
    Parse cog-24.fun.tab or fun-20.tab.
    Extracts category letter -> {description, group, color}.
    """
    categories = []
    if not cog_fun_path or not os.path.isfile(cog_fun_path):
        return FALLBACK_CATEGORIES

    opener = gzip.open if cog_fun_path.endswith(".gz") else open
    with opener(cog_fun_path, "rt", encoding="utf-8", errors="replace") as f:
        current_group = "UNASSIGNED"
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Look for functional group headers or tab lines
            parts = line.split("\t")
            if len(parts) >= 2:
                code_cand = parts[0].strip().replace("[", "").replace("]", "")
                if len(code_cand) == 1 and code_cand.isalpha():
                    code = code_cand.upper()
                    # Check if color is part 1 or desc is part 1
                    desc = parts[-1].strip()
                    group = parts[1].strip() if len(parts) > 2 and parts[1] != desc else current_group
                    categories.append((code, desc, group))
            elif line.isupper() and len(line) > 3 and not "\t" in line:
                current_group = line

    if not categories:
        return FALLBACK_CATEGORIES
    return categories

def parse_bitscore_specific(bitscore_path):
    """
    Parse bitscore_specific.txt from NCBI CDD.
    Format: PSSM-Id / Accession \t threshold_bitscore.
    Returns: dict {pssm_id_or_cog_id: float_threshold}
    """
    thresholds = {}
    if not bitscore_path or not os.path.isfile(bitscore_path):
        return thresholds

    opener = gzip.open if bitscore_path.endswith(".gz") else open
    with opener(bitscore_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) >= 2:
                model_id = parts[0].strip()
                try:
                    thresh = float(parts[1].strip())
                    thresholds[model_id] = thresh
                except ValueError:
                    continue
    return thresholds

def main():
    parser = argparse.ArgumentParser(
        description="Gate 0: Audit Database Provenance and Cog_LE vs COG2024 Compatibility."
    )
    parser.add_argument("--cddid", required=True, help="Path to cddid.tbl (or cddid.tbl.gz)")
    parser.add_argument("--cog_def", required=True, help="Path to cog-24.def.tab")
    parser.add_argument("--cog_fun", default="", help="Path to cog-24.fun.tab (mandatory for official category defs)")
    parser.add_argument("--bitscore_specific", default="", help="Path to bitscore_specific.txt (recommended)")
    parser.add_argument("--cogle_db", default="", help="Path or prefix to Cog_LE PSSM database files")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for audit report and mapping")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"INFO: Auditing cddid={args.cddid}, cog_def={args.cog_def}, cog_fun={args.cog_fun}")

    # 1. Compute checksums
    cddid_sha = compute_sha256(args.cddid)
    cog_def_sha = compute_sha256(args.cog_def)
    cog_fun_sha = compute_sha256(args.cog_fun) if args.cog_fun else None
    bitscore_sha = compute_sha256(args.bitscore_specific) if args.bitscore_specific else None
    
    cogle_checksums = {}
    if args.cogle_db:
        db_dir = os.path.dirname(args.cogle_db) if os.path.dirname(args.cogle_db) else "."
        db_prefix = os.path.basename(args.cogle_db)
        for fpath in glob.glob(os.path.join(db_dir, f"{db_prefix}*")):
            if os.path.isfile(fpath):
                cogle_checksums[os.path.basename(fpath)] = compute_sha256(fpath)

    # 2. Parse assets
    pssm_dict, cog_to_pssm = parse_cddid(args.cddid)
    cog_defs = parse_cog_defs(args.cog_def)
    categories = parse_cog_fun(args.cog_fun)
    specific_thresholds = parse_bitscore_specific(args.bitscore_specific)

    cogle_cog_ids = set(cog_to_pssm.keys())
    cog2024_ids = set(cog_defs.keys())

    intersection = cogle_cog_ids & cog2024_ids
    only_cog2024 = sorted(cog2024_ids - cogle_cog_ids)
    only_cogle = sorted(cogle_cog_ids - cog2024_ids)

    n_cog2024 = len(cog2024_ids)
    n_cogle = len(cogle_cog_ids)
    n_intersection = len(intersection)
    coverage_fraction = round(n_intersection / max(1, n_cog2024), 6)

    # 3. Methodological designation
    method_designation = "CDD Cog_LE-based COG annotation mapped to COG2024 definitions"

    print(f"INFO: COG2024 definitions: {n_cog2024} COGs")
    print(f"INFO: Cog_LE PSSMs: {n_cogle} COGs ({len(pssm_dict)} total PSSM models)")
    print(f"INFO: Intersection: {n_intersection} ({coverage_fraction*100:.2f}% of COG2024)")
    print(f"INFO: Specific Bitscore Thresholds Loaded: {len(specific_thresholds)} models")

    # 4. Write audit reports
    report_data = {
        "method_designation": method_designation,
        "n_cog2024": n_cog2024,
        "n_cogle_cogs": n_cogle,
        "n_cogle_pssms": len(pssm_dict),
        "n_intersection": n_intersection,
        "n_only_cog2024": len(only_cog2024),
        "n_only_cogle": len(only_cogle),
        "coverage_fraction": coverage_fraction,
        "n_specific_bitscore_thresholds": len(specific_thresholds),
        "database_checksums": {
            "cddid_sha256": cddid_sha,
            "cog_def_sha256": cog_def_sha,
            "cog_fun_sha256": cog_fun_sha,
            "bitscore_specific_sha256": bitscore_sha,
            "cogle_files": cogle_checksums
        },
        "source_paths": {
            "cddid": os.path.abspath(args.cddid),
            "cog_def": os.path.abspath(args.cog_def),
            "cog_fun": os.path.abspath(args.cog_fun) if args.cog_fun else "",
            "bitscore_specific": os.path.abspath(args.bitscore_specific) if args.bitscore_specific else "",
            "cogle_db": os.path.abspath(args.cogle_db) if args.cogle_db else ""
        }
    }

    report_path = os.path.join(args.out_dir, "cog_database_provenance_report.json")
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=2)

    # Write only_cog2024.tsv
    only_cog2024_path = os.path.join(args.out_dir, "only_cog2024.tsv")
    with open(only_cog2024_path, "w", encoding="utf-8") as f:
        f.write("cog_id\tcategory\tname\tgene\n")
        for cid in only_cog2024:
            info = cog_defs.get(cid, {})
            f.write(f"{cid}\t{info.get('category','')}\t{info.get('name','')}\t{info.get('gene','')}\n")

    # Write only_cogle.tsv
    only_cogle_path = os.path.join(args.out_dir, "only_cogle.tsv")
    with open(only_cogle_path, "w", encoding="utf-8") as f:
        f.write("cog_id\tpssm_ids\n")
        for cid in only_cogle:
            pssms = ",".join(cog_to_pssm.get(cid, []))
            f.write(f"{cid}\t{pssms}\n")

    # Write official category definitions (cog_category_definition.tsv)
    cat_def_path = os.path.join(args.out_dir, "cog_category_definition.tsv")
    with open(cat_def_path, "w", encoding="utf-8") as f:
        f.write("category_code\tcategory_name\tfunctional_class\n")
        for code, name, fclass in categories:
            f.write(f"{code}\t{name}\t{fclass}\n")

    # Write cdd_pssm_mapping.tsv with specific_bitscore_threshold
    mapping_path = os.path.join(args.out_dir, "cdd_pssm_mapping.tsv")
    with open(mapping_path, "w", encoding="utf-8") as f:
        f.write("pssm_id\tcog_id\tpssm_length\tcdd_name\tcdd_desc\tcog_category\tcog_name\tspecific_bitscore_threshold\n")
        for pssm_id, p_info in sorted(pssm_dict.items(), key=lambda x: x[0]):
            cid = p_info["cog_id"]
            c_def = cog_defs.get(cid, {})
            # Lookup specific bitscore threshold by pssm_id or cog_id
            spec_thresh = specific_thresholds.get(pssm_id) or specific_thresholds.get(cid)
            spec_str = f"{spec_thresh:.1f}" if spec_thresh is not None else "NA"
            f.write(
                f"{pssm_id}\t{cid}\t{p_info['length']}\t{p_info['name']}\t"
                f"{p_info['desc']}\t{c_def.get('category','S')}\t{c_def.get('name','')}\t{spec_str}\n"
            )

    print(f"SUCCESS: Gate 0 Audit complete. Report saved to {report_path}")

if __name__ == "__main__":
    main()
