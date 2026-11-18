#!/usr/bin/env python3
"""
Aggregate dbCAN, CAMPER, and Phydon results for all MAGs into a unified trait matrix.

Output:
  all_mag_traits_matrix.tsv
"""
import sys
import os
import glob
import re
import argparse
from collections import defaultdict


def parse_dbcan_overview(overview_path):
    """
    Parses dbCAN overview.txt file to extract CAZyme counts and class distributions.
    Standard dbCAN overview.txt format:
    Gene ID   HMMER   Diamond   dbCAN_sub   #ofTools
    """
    cazyme_counts = defaultdict(int)
    cazyme_classes = defaultdict(int)
    total_cazymes = 0

    if not os.path.exists(overview_path) or os.path.getsize(overview_path) == 0:
        return total_cazymes, cazyme_classes, ""

    try:
        with open(overview_path, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                # The last column is usually #ofTools (integer). A gene is confident if tools >= 2 (or >=1)
                gene_id = parts[0]
                annotations = " ".join(parts[1:])
                # Match CAZyme classes: GH, GT, PL, CE, AA, CBM
                hits = re.findall(r'\b(GH|GT|PL|CE|AA|CBM)\d+', annotations)
                if hits:
                    total_cazymes += 1
                    for fam in set(hits):
                        cazyme_counts[fam] += 1
                        cls_name = re.match(r'^[A-Z]+', fam).group()
                        cazyme_classes[cls_name] += 1
    except Exception as e:
        print(f"WARNING: Error parsing dbCAN overview {overview_path}: {e}", file=sys.stderr)

    # Top 5 most abundant families
    sorted_fams = sorted(cazyme_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    top_fams_str = ";".join(f"{k}({v})" for k, v in sorted_fams)
    return total_cazymes, cazyme_classes, top_fams_str


def parse_camper_hits(camper_path):
    """
    Parses CAMPER hits file (TSV or HMM tblout) to summarize polyphenol enzyme counts.
    """
    total_hits = 0
    unique_enzymes = set()

    if not os.path.exists(camper_path) or os.path.getsize(camper_path) == 0:
        return total_hits, ""

    try:
        with open(camper_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                total_hits += 1
                if len(parts) >= 2:
                    # Model name or enzyme name
                    enzyme = parts[1] if not parts[0].startswith(">") else parts[0]
                    unique_enzymes.add(enzyme.split(".")[0])
                else:
                    # Generic HMM tblout
                    tokens = line.split()
                    if len(tokens) >= 3:
                        unique_enzymes.add(tokens[2])
    except Exception as e:
        print(f"WARNING: Error parsing CAMPER file {camper_path}: {e}", file=sys.stderr)

    enzymes_str = ";".join(sorted(unique_enzymes)[:8])
    return total_hits, enzymes_str


def parse_growth_prediction(growth_path):
    """
    Parses gRodon2 (or legacy Phydon) output TSV.
    Columns: MAG, DoublingTime_hours, MaxGrowthRate_per_hour, CUB_HEG, GrowthCategory, Method, Mode, Temperature, Status
    """
    data = {
        "DoublingTime_hours": "NA",
        "MaxGrowthRate_per_hour": "NA",
        "CUB_HEG": "NA",
        "GrowthCategory": "NA",
        "Method": "gRodon2",
        "Status": "Missing"
    }

    if not os.path.exists(growth_path) or os.path.getsize(growth_path) == 0:
        return data

    try:
        with open(growth_path, "r", encoding="utf-8", errors="ignore") as f:
            headers = [h.strip() for h in f.readline().split("\t") if h.strip()]
            line = f.readline().strip()
            if line:
                row = [c.strip() for c in line.split("\t")]
                if len(row) >= len(headers):
                    for idx, h in enumerate(headers[:-1]):
                        if h in data:
                            data[h] = re.sub(r'[\r\n\t]+', ' ', row[idx]).strip()
                    last_h = headers[-1]
                    if last_h in data:
                        data[last_h] = re.sub(r'[\r\n\t]+', ' ', " ".join(row[len(headers)-1:])).strip()
                elif len(row) >= 2:
                    for i in range(min(5, len(row))):
                        if i < len(headers) and headers[i] in data:
                            data[headers[i]] = row[i]
                    if len(row) >= 6:
                        data["Method"] = row[5]
                    if len(row) >= 7:
                        data["Status"] = re.sub(r'[\r\n\t]+', ' ', " ".join(row[6:])).strip()
                    else:
                        data["Status"] = "Error"
    except Exception as e:
        print(f"WARNING: Error parsing growth rate file {growth_path}: {e}", file=sys.stderr)

    return data


def parse_drep_cdb(cdb_path):
    """
    Parses dRep Cdb.csv to map clusters and their constituent MAGs.
    Cdb.csv headers: genome,secondary_cluster,primary_cluster
    """
    cluster_members = defaultdict(list)
    genome_to_cluster = {}

    if not cdb_path or not os.path.exists(cdb_path) or os.path.getsize(cdb_path) == 0:
        return cluster_members, genome_to_cluster

    try:
        with open(cdb_path, "r", encoding="utf-8", errors="ignore") as f:
            header_line = f.readline()
            headers = [h.strip() for h in header_line.split(",")]
            g_idx = headers.index("genome") if "genome" in headers else 0
            c_idx = headers.index("secondary_cluster") if "secondary_cluster" in headers else 1

            for line in f:
                parts = [p.strip() for p in line.split(",")]
                if len(parts) > max(g_idx, c_idx):
                    # Clean extension
                    genome = re.sub(r'\.(fa|fasta|fna)(\.gz)?$', '', parts[g_idx])
                    cluster = parts[c_idx]
                    cluster_members[cluster].append(genome)
                    genome_to_cluster[genome] = cluster
    except Exception as e:
        print(f"WARNING: Error parsing dRep Cdb.csv {cdb_path}: {e}", file=sys.stderr)

    return cluster_members, genome_to_cluster


def parse_drep_clusters(cdb_r1, cdb_r2, single_cdb):
    """
    Parses single or two-round dRep Cdb.csv files.
    Returns a dict: mag -> {'cluster': str, 'size': int, 'members': list[str]}
    """
    if cdb_r1 and cdb_r2 and os.path.exists(cdb_r1) and os.path.exists(cdb_r2):
        # Two-round dRep
        r1_members, r1_g2c = parse_drep_cdb(cdb_r1)
        r2_members, r2_g2c = parse_drep_cdb(cdb_r2)
        
        mapping = {}
        for r2_cluster, r1_winners in r2_members.items():
            all_raw_mags = []
            for w in r1_winners:
                r1_c = r1_g2c.get(w, "")
                all_raw_mags.extend(r1_members.get(r1_c, [w]))
            
            # De-duplicate raw mags while preserving order
            seen = set()
            dedup_mags = [m for m in all_raw_mags if not (m in seen or seen.add(m))]
            
            # Each r1_winner that survived to r2
            for w in r1_winners:
                mapping[w] = {
                    "cluster": f"R2_{r2_cluster}",
                    "size": len(dedup_mags),
                    "members": dedup_mags
                }
        return mapping

    if single_cdb and os.path.exists(single_cdb):
        cluster_members, g2c = parse_drep_cdb(single_cdb)
        mapping = {}
        for g, c in g2c.items():
            mems = cluster_members.get(c, [g])
            mapping[g] = {
                "cluster": c,
                "size": len(mems),
                "members": mems
            }
        return mapping

    return {}


def parse_all_dbcan_sources(results_dir):
    """
    Pre-scans both batch dbcan overviews and individual dbcan overviews.
    Returns: dict[mag -> (total_cazymes, cazyme_classes_dict, top_fams_str)]
    """
    data = {}
    # 1. Batch overview files (supports both overview.txt and dbCAN 5.2+ overview.tsv)
    batch_pattern_txt = os.path.join(results_dir, "batches", "*", "dbcan", "overview.txt")
    batch_pattern_tsv = os.path.join(results_dir, "batches", "*", "dbcan", "overview.tsv")
    b_files = sorted(set(glob.glob(batch_pattern_txt) + glob.glob(batch_pattern_tsv)))
    for b_file in b_files:
        if not os.path.exists(b_file) or os.path.getsize(b_file) == 0:
            continue
        try:
            with open(b_file, "r", encoding="utf-8", errors="ignore") as f:
                f.readline()  # header
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue
                    gene_id = parts[0]
                    # Check if gene_id is MAG___gene_id
                    if "___" in gene_id:
                        mag, _ = gene_id.split("___", 1)
                    else:
                        continue
                    
                    if mag not in data:
                        data[mag] = (defaultdict(int), defaultdict(int), 0)
                    
                    caz_counts, caz_classes, tot = data[mag]
                    annotations = " ".join(parts[1:])
                    hits = re.findall(r'\b(GH|GT|PL|CE|AA|CBM)\d+', annotations)
                    if hits:
                        tot += 1
                        for fam in set(hits):
                            caz_counts[fam] += 1
                            cls_name = re.match(r'^[A-Z]+', fam).group()
                            caz_classes[cls_name] += 1
                        data[mag] = (caz_counts, caz_classes, tot)
        except Exception as e:
            print(f"WARNING: Error parsing batch dbcan file {b_file}: {e}", file=sys.stderr)

    # Format pre-scanned batch data
    formatted = {}
    for mag, (caz_counts, caz_classes, tot) in data.items():
        sorted_fams = sorted(caz_counts.items(), key=lambda x: x[1], reverse=True)[:5]
        top_fams_str = ";".join(f"{k}({v})" for k, v in sorted_fams)
        formatted[mag] = (tot, caz_classes, top_fams_str)
    return formatted


def parse_all_grodon2_sources(results_dir):
    """
    Pre-scans batch grodon2 tables.
    Returns: dict[mag -> growth_dict]
    """
    data = {}
    batch_pattern = os.path.join(results_dir, "batches", "*", "grodon2", "*grodon2*.tsv")
    for b_file in glob.glob(batch_pattern):
        if not os.path.exists(b_file) or os.path.getsize(b_file) == 0:
            continue
        try:
            with open(b_file, "r", encoding="utf-8", errors="ignore") as f:
                header_line = f.readline()
                if not header_line:
                    continue
                headers = [h.strip() for h in header_line.split("\t")]
                if "MAG" not in headers:
                    continue
                mag_idx = headers.index("MAG")
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    row = [c.strip() for c in line.split("\t")]
                    if len(row) >= 2:
                        m_id = row[mag_idx]
                        item = {}
                        if len(row) >= len(headers):
                            for idx, h in enumerate(headers[:-1]):
                                item[h] = re.sub(r'[\r\n\t]+', ' ', row[idx]).strip()
                            last_h = headers[-1]
                            item[last_h] = re.sub(r'[\r\n\t]+', ' ', " ".join(row[len(headers)-1:])).strip()
                        else:
                            for i in range(min(5, len(row))):
                                if i < len(headers):
                                    item[headers[i]] = row[i]
                            if len(row) >= 6:
                                item["Method"] = row[5]
                            if len(row) >= 7:
                                item["Status"] = re.sub(r'[\r\n\t]+', ' ', " ".join(row[6:])).strip()
                            else:
                                item["Status"] = "Error"
                        data[m_id] = item
        except Exception as e:
            print(f"WARNING: Error parsing batch grodon2 file {b_file}: {e}", file=sys.stderr)
    return data


def parse_all_camper_sources(results_dir):
    """
    Pre-scans batch camper tables.
    Returns: dict[mag -> (total_hits, enzymes_str)]
    """
    data = {}
    batch_pattern = os.path.join(results_dir, "batches", "*", "camper", "*camper*.tsv")
    for b_file in glob.glob(batch_pattern):
        if not os.path.exists(b_file) or os.path.getsize(b_file) == 0:
            continue
        try:
            with open(b_file, "r", encoding="utf-8", errors="ignore") as f:
                header_line = f.readline()
                headers = [h.strip() for h in header_line.split("\t")]
                has_mag_col = "MAG" in headers
                mag_col_idx = headers.index("MAG") if has_mag_col else -1
                target_col_idx = headers.index("Target_HMM") if "Target_HMM" in headers else (2 if (has_mag_col and len(headers) > 2) else 1)
                
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if has_mag_col and len(parts) > mag_col_idx:
                        mag = parts[mag_col_idx]
                        enzyme = parts[target_col_idx] if len(parts) > target_col_idx else (parts[1] if len(parts) >= 2 else "Enzyme")
                    elif "___" in parts[0]:
                        mag, _ = parts[0].split("___", 1)
                        enzyme = parts[target_col_idx] if len(parts) > target_col_idx else (parts[1] if len(parts) >= 2 else "Enzyme")
                    else:
                        continue

                    if mag not in data:
                        data[mag] = (0, set())
                    tot, enz_set = data[mag]
                    tot += 1
                    enz_set.add(enzyme.split(".")[0])
                    data[mag] = (tot, enz_set)
        except Exception as e:
            print(f"WARNING: Error parsing batch camper file {b_file}: {e}", file=sys.stderr)

    formatted = {}
    for mag, (tot, enz_set) in data.items():
        enzymes_str = ";".join(sorted(enz_set)[:8])
        formatted[mag] = (tot, enzymes_str)
    return formatted


def parse_all_genomespot_sources(results_dir):
    """
    Pre-scans batch genomespot tables and individual prediction files.
    Returns: dict[mag -> dict of traits]
    """
    data = {}
    # 1. Scan batch summary tables
    batch_pattern = os.path.join(results_dir, "batches", "*", "genomespot", "*summary*.tsv")
    for b_file in glob.glob(batch_pattern):
        if not os.path.exists(b_file) or os.path.getsize(b_file) == 0:
            continue
        try:
            with open(b_file, "r", encoding="utf-8", errors="ignore") as f:
                header_line = f.readline()
                if not header_line:
                    continue
                headers = [h.strip() for h in header_line.split("\t")]
                if "MAG" not in headers:
                    continue
                mag_idx = headers.index("MAG")
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = [p.strip() for p in line.split("\t")]
                    if len(parts) > mag_idx:
                        mag = parts[mag_idx]
                        item = {}
                        for idx, h in enumerate(headers):
                            if idx < len(parts):
                                item[h] = parts[idx]
                            else:
                                item[h] = "NA"
                        data[mag] = item
        except Exception as e:
            print(f"WARNING: Error parsing batch genomespot file {b_file}: {e}", file=sys.stderr)

    # 2. Fallback scan for individual predictions if not already found
    indiv_patterns = [
        os.path.join(results_dir, "batches", "*", "genomespot", "*.predictions.tsv"),
        os.path.join(results_dir, "genomespot", "*", "*.predictions.tsv"),
        os.path.join(results_dir, "genomespot", "*.predictions.tsv")
    ]
    for pattern in indiv_patterns:
        for p_file in glob.glob(pattern):
            if not os.path.exists(p_file) or os.path.getsize(p_file) == 0:
                continue
            bname = os.path.basename(p_file)
            mag = bname.split(".")[0]
            if mag in data:
                continue
            item = {
                "OGT": "NA", "Tmin": "NA", "Tmax": "NA",
                "pH_opt": "NA", "pH_min": "NA", "pH_max": "NA",
                "Salinity_opt": "NA", "Salinity_min": "NA", "Salinity_max": "NA",
                "Oxygen_tolerance": "NA", "Oxygen_prob": "NA", "Status": "OK"
            }
            try:
                with open(p_file, "r", encoding="utf-8", errors="ignore") as f:
                    f.readline()
                    for line in f:
                        parts = [p.strip() for p in line.split("\t")]
                        if len(parts) < 2:
                            continue
                        target = parts[0].lower()
                        val = parts[1]
                        if "temperature_optimum" in target:
                            item["OGT"] = val
                        elif "temperature_min" in target:
                            item["Tmin"] = val
                        elif "temperature_max" in target:
                            item["Tmax"] = val
                        elif "ph_optimum" in target:
                            item["pH_opt"] = val
                        elif "ph_min" in target:
                            item["pH_min"] = val
                        elif "ph_max" in target:
                            item["pH_max"] = val
                        elif "salinity_optimum" in target:
                            item["Salinity_opt"] = val
                        elif "salinity_min" in target:
                            item["Salinity_min"] = val
                        elif "salinity_max" in target:
                            item["Salinity_max"] = val
                        elif "oxygen" in target:
                            item["Oxygen_tolerance"] = val
                            if len(parts) >= 3 and parts[2]:
                                item["Oxygen_prob"] = parts[2]
                data[mag] = item
            except Exception:
                pass
    return data


def main():
    parser = argparse.ArgumentParser(description="Merge dbCAN, CAMPER, and gRodon2 MAG traits.")
    parser.add_argument("--mags", nargs="*", default=[], help="List of MAG identifiers")
    parser.add_argument("--mags-file", default=None, help="File containing list of MAG identifiers (one per line)")
    parser.add_argument("--results-dir", required=True, help="Path to results/mag_traits directory")
    parser.add_argument("--output", required=True, help="Path to output all_mag_traits_matrix.tsv")
    parser.add_argument("--cdb", default=None, help="Optional path to single dRep Cdb.csv")
    parser.add_argument("--cdb-r1", default=None, help="Optional path to dRep Round 1 Cdb.csv (ANI 0.99)")
    parser.add_argument("--cdb-r2", default=None, help="Optional path to dRep Round 2 Cdb.csv (ANI 0.95)")
    args = parser.parse_args()

    # Collect MAGs
    all_mags = set(args.mags)
    if args.mags_file and os.path.exists(args.mags_file):
        with open(args.mags_file, "r", encoding="utf-8") as f:
            for line in f:
                val = line.strip()
                if val:
                    all_mags.add(val)

    results_dir = args.results_dir
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    # Pre-parse batch sources
    batch_dbcan = parse_all_dbcan_sources(results_dir)
    batch_camper = parse_all_camper_sources(results_dir)
    batch_grodon2 = parse_all_grodon2_sources(results_dir)
    batch_genomespot = parse_all_genomespot_sources(results_dir)

    # Parse dRep mapping
    drep_mapping = parse_drep_clusters(args.cdb_r1, args.cdb_r2, args.cdb)
    has_drep = bool(drep_mapping)

    # If no mags specified directly, discover from parsed data or drep
    if not all_mags:
        batch_mags_pattern = os.path.join(results_dir, "batches", "*", "mag_list.tsv")
        for ml_file in glob.glob(batch_mags_pattern):
            if os.path.exists(ml_file):
                try:
                    with open(ml_file, "r", encoding="utf-8") as f:
                        for line in f:
                            line = line.strip()
                            if line:
                                parts = line.split("\t")
                                if parts:
                                    all_mags.add(parts[0].strip())
                except Exception:
                    pass
        all_mags.update(drep_mapping.keys())
        all_mags.update(batch_grodon2.keys())
        all_mags.update(batch_dbcan.keys())
        all_mags.update(batch_camper.keys())
        all_mags.update(batch_genomespot.keys())

    header = ["MAG"]
    if has_drep:
        header.extend(["dRep_Cluster", "Cluster_Size", "Cluster_Members"])

    header.extend([
        "DoublingTime_hours",
        "MaxGrowthRate_per_hour",
        "CUB_HEG",
        "GrowthCategory",
        "Total_CAZymes",
        "GH_GlycosideHydrolases",
        "GT_GlycosylTransferases",
        "PL_PolysaccharideLyases",
        "CE_CarbohydrateEsterases",
        "AA_AuxiliaryActivities",
        "CBM_BindingModules",
        "Top_CAZyme_Families",
        "Total_CAMPER_Hits",
        "CAMPER_Key_Enzymes",
        "GenomeSPOT_OGT_C",
        "GenomeSPOT_Tmin_C",
        "GenomeSPOT_Tmax_C",
        "GenomeSPOT_pH_Opt",
        "GenomeSPOT_pH_Min",
        "GenomeSPOT_pH_Max",
        "GenomeSPOT_Salinity_Opt_pct",
        "GenomeSPOT_Salinity_Min_pct",
        "GenomeSPOT_Salinity_Max_pct",
        "GenomeSPOT_Oxygen_Tolerance",
        "GenomeSPOT_Oxygen_Prob",
        "GenomeSPOT_Status",
        "Growth_Predict_Method",
        "Growth_Predict_Status"
    ])

    rows = []
    for mag in sorted(all_mags):
        # 1. dbCAN
        if mag in batch_dbcan:
            total_caz, caz_classes, top_caz = batch_dbcan[mag]
        else:
            dbcan_overview = os.path.join(results_dir, "dbcan", mag, "overview.txt")
            total_caz, caz_classes, top_caz = parse_dbcan_overview(dbcan_overview)

        # 2. CAMPER
        if mag in batch_camper:
            total_camper, camper_enzymes = batch_camper[mag]
        else:
            camper_file = os.path.join(results_dir, "camper", mag, f"{mag}_camper_hits.tsv")
            total_camper, camper_enzymes = parse_camper_hits(camper_file)

        # 3. gRodon2
        if mag in batch_grodon2:
            growth_data = batch_grodon2[mag]
        else:
            growth_file = os.path.join(results_dir, "grodon2", mag, f"{mag}_grodon2_growth.tsv")
            if not os.path.exists(growth_file):
                legacy_file = os.path.join(results_dir, "phydon", mag, f"{mag}_phydon_growth.tsv")
                if os.path.exists(legacy_file):
                    growth_file = legacy_file
            growth_data = parse_growth_prediction(growth_file)

        # 4. GenomeSPOT
        gspot_data = batch_genomespot.get(mag, {})

        row = [mag]
        if has_drep:
            cluster_info = drep_mapping.get(mag, {"cluster": "Unclustered", "size": 1, "members": [mag]})
            row.extend([
                cluster_info["cluster"],
                str(cluster_info["size"]),
                ";".join(cluster_info["members"])
            ])

        row.extend([
            str(growth_data.get("DoublingTime_hours", "NA")),
            str(growth_data.get("MaxGrowthRate_per_hour", "NA")),
            str(growth_data.get("CUB_HEG", "NA")),
            str(growth_data.get("GrowthCategory", "NA")),
            str(total_caz),
            str(caz_classes.get("GH", 0)),
            str(caz_classes.get("GT", 0)),
            str(caz_classes.get("PL", 0)),
            str(caz_classes.get("CE", 0)),
            str(caz_classes.get("AA", 0)),
            str(caz_classes.get("CBM", 0)),
            top_caz if top_caz else "None",
            str(total_camper),
            camper_enzymes if camper_enzymes else "None",
            str(gspot_data.get("OGT", "NA")),
            str(gspot_data.get("Tmin", "NA")),
            str(gspot_data.get("Tmax", "NA")),
            str(gspot_data.get("pH_opt", "NA")),
            str(gspot_data.get("pH_min", "NA")),
            str(gspot_data.get("pH_max", "NA")),
            str(gspot_data.get("Salinity_opt", "NA")),
            str(gspot_data.get("Salinity_min", "NA")),
            str(gspot_data.get("Salinity_max", "NA")),
            str(gspot_data.get("Oxygen_tolerance", "NA")),
            str(gspot_data.get("Oxygen_prob", "NA")),
            str(gspot_data.get("Status", "Missing")),
            str(growth_data.get("Method", "gRodon2")),
            str(growth_data.get("Status", "Missing"))
        ])
        cleaned_row = [re.sub(r'[\r\n\t]+', ' ', str(c)).strip() for c in row]
        if len(cleaned_row) == len(header):
            rows.append(cleaned_row)
        else:
            print(f"WARNING: Skipping malformed row for {mag} with {len(cleaned_row)} fields (expected {len(header)})", file=sys.stderr)

    with open(args.output, "w", encoding="utf-8") as out_f:
        out_f.write("\t".join(header) + "\n")
        for r in rows:
            out_f.write("\t".join(r) + "\n")

    print(f"INFO: Successfully generated comprehensive MAG traits matrix for {len(rows)} MAGs at {args.output}")


if __name__ == "__main__":
    main()
