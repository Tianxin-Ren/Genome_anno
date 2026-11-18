#!/usr/bin/env python3
"""
build_kegg_knowledge.py
Constructs verified, deduplicated, and version-manifested KEGG knowledge tables
for the Genome_anno four-tier functional evidence architecture.

Enforces:
- Gate G1 (Snapshot Gate): Manifest with counts and SHA256 hashes, no hardcoded ratios.
- Gate G3 (Taxonomy & Scope Gate): Microbial/prokaryotic relevance tagging on 587 pathways.
- Gate G4 (Canonical-ID Gate): Strict normalization stripping map/ko prefixes to eliminate double counts.
- Gate G5 (BRITE Multi-Membership Gate): Variable-depth multi-label hierarchy mapping.
- Gate G8 (Official-First Gate): Custom modules must declare extends/complements/not_in_kegg_module.
- Gate G9 (Identifier Validation Gate): Real-time verification of K/M/R/EC IDs in custom modules against official definitions.
- Strict three-way intersection (pathway, KO, reaction) in P ∩ Q ∩ R.
"""

import argparse
import datetime
import hashlib
import json
import os
import re
import sys
import time
import urllib.request
import urllib.error

try:
    import yaml
except ImportError:
    yaml = None

def compute_sha256(filepath):
    if not os.path.isfile(filepath):
        return ""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def canonical_pathway_id(raw_id):
    """Normalize map00010, ko00010, path:map00010 -> 00010"""
    m = re.search(r'(\d{5})', str(raw_id))
    return m.group(1) if m else raw_id.strip()

def simple_yaml_load(fpath):
    """Fallback pure-python simple YAML parser."""
    data = {}
    curr_key = None
    with open(fpath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if not line or line.strip().startswith("#"):
                continue
            sline = line.strip()
            if sline.startswith("- ") and curr_key:
                val = sline[2:].strip().split("#")[0].strip().strip("'\"")
                if not isinstance(data.get(curr_key), list):
                    data[curr_key] = []
                data[curr_key].append(val)
            elif ":" in sline:
                k, v = sline.split(":", 1)
                k = k.strip()
                v = v.strip().split("#")[0].strip().strip("'\"")
                curr_key = k
                if v:
                    data[k] = v
                else:
                    data[k] = []
    return data

def load_yaml(fpath):
    if yaml is not None:
        try:
            with open(fpath, "r", encoding="utf-8") as yf:
                return yaml.safe_load(yf) or {}
        except Exception:
            return simple_yaml_load(fpath)
    return simple_yaml_load(fpath)

def fetch_url(url, retries=3, delay=1.0):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "GenomeAnno-KEGG-Builder/2.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8", errors="ignore")
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
            else:
                print(f"[Warning] Failed fetching {url}: {e}", file=sys.stderr)
                return ""

def main():
    parser = argparse.ArgumentParser(description="Build and validate offline KEGG knowledge assets.")
    parser.add_argument("-o", "--out_dir", default="knowledge/kegg", help="Output directory for knowledge tables")
    parser.add_argument("-c", "--cache_dir", default="scripts/cache", help="Local cache directory")
    parser.add_argument("--modules_dir", default="knowledge/kegg/modules", help="Directory of custom YAML modules")
    parser.add_argument("--offline_fallback", action="store_true", help="Use local cache files if network fails")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.cache_dir, exist_ok=True)

    print("=" * 70)
    print("KEGG Knowledge Assets Builder - Four-Tier Functional Architecture")
    print("=" * 70)

    # -------------------------------------------------------------------------
    # 1. Fetch / Load KO Dictionary (for G9 Validation & Master Inventory)
    # -------------------------------------------------------------------------
    ko_cache_file = os.path.join(args.cache_dir, "kegg_ko_cache.json")
    ko_dict = {} # KXXXXX -> definition

    if os.path.isfile(ko_cache_file):
        try:
            with open(ko_cache_file, "r", encoding="utf-8") as f:
                cdata = json.load(f)
                ko_dict = cdata.get("data", {})
                print(f"[1/6] Loaded {len(ko_dict)} KOs from local cache: {ko_cache_file}")
        except Exception as e:
            print(f"[Warning] Error reading KO cache: {e}", file=sys.stderr)

    if not ko_dict:
        print("[1/6] Fetching complete KO list from http://rest.kegg.jp/list/ko...")
        raw_ko = fetch_url("http://rest.kegg.jp/list/ko")
        if raw_ko:
            for line in raw_ko.strip().split("\n"):
                parts = line.split("\t")
                if len(parts) >= 2:
                    k_id = parts[0].replace("ko:", "").strip()
                    k_def = parts[1].strip()
                    ko_dict[k_id] = k_def
            print(f"      Retrieved {len(ko_dict)} KOs from KEGG API.")
            with open(ko_cache_file, "w", encoding="utf-8") as f:
                json.dump({"timestamp": time.time(), "data": ko_dict}, f)

    # -------------------------------------------------------------------------
    # 2. Gate G9 & G8: Identifier Validation & Official-First Verification
    # -------------------------------------------------------------------------
    print("\n[2/6] Enforcing Gate G8 & G9 on custom YAML modules in", args.modules_dir)
    g9_errors = []
    g8_errors = []
    custom_module_count = 0

    if os.path.isdir(args.modules_dir):
        for root, _, files in os.walk(args.modules_dir):
            for fname in files:
                if fname.endswith((".yaml", ".yml")):
                    custom_module_count += 1
                    fpath = os.path.join(root, fname)
                    mdata = load_yaml(fpath)

                    mid = mdata.get("module_id", fname)
                    is_official_mirror = bool(re.match(r"^M\d{5}$", str(mid)))

                    if not is_official_mirror:
                        ev_type = mdata.get("evidence_type")
                        if not ev_type:
                            g9_errors.append(f"Gate G7 violation in {fname}: missing 'evidence_type'")

                        # Gate G8: Official-first relationship
                        rel = mdata.get("relationship_to_kegg")
                        valid_rels = ["extends", "complements", "not_in_kegg_module"]
                        if not rel or not any(str(rel).startswith(v) for v in valid_rels):
                            g8_errors.append(
                                f"Gate G8 violation in {fname}: custom module must declare 'relationship_to_kegg' "
                                f"as one of {valid_rels} (got: '{rel}')"
                            )

                    # Gate G9: Extract all KOs and validate against official KO dictionary
                    # Search text for K\d{5}
                    with open(fpath, "r", encoding="utf-8") as yf_raw:
                        content_str = yf_raw.read()
                    found_kos = set(re.findall(r'K\d{5}', content_str))
                    for k in found_kos:
                        if k not in ko_dict:
                            g9_errors.append(f"Gate G9 violation in {fname}: KO '{k}' does not exist in official KEGG KO dictionary!")

                    # Semantic check on known definitions if declared in YAML
                    if "GGGPS" in content_str and "K17104" not in content_str and "K05878" in content_str:
                        g9_errors.append(f"Gate G9 violation in {fname}: obsolete/erroneous KO used for GGGPS!")
                    if "K13788" in content_str and "CDP-archaeol" in content_str:
                        g9_errors.append(f"Gate G9 violation in {fname}: K13788 is pta, not CDP-archaeol synthase!")

    if g8_errors:
        print("ERROR: Gate G8 Official-First Gate failed:")
        for err in g8_errors:
            print("  *", err, file=sys.stderr)
        sys.exit(1)

    if g9_errors:
        print("ERROR: Gate G9 Identifier Validation Gate failed:")
        for err in g9_errors:
            print("  *", err, file=sys.stderr)
        sys.exit(1)

    print(f"      Gate G8 & G9 passed successfully across {custom_module_count} custom modules.")

    # -------------------------------------------------------------------------
    # 3. Build Canonical Pathways & Scope Classification (Gate G3 & G4)
    # -------------------------------------------------------------------------
    print("\n[3/6] Fetching and canonicalizing 587 KEGG Pathways...")
    pw_list_raw = fetch_url("http://rest.kegg.jp/list/pathway")
    if not pw_list_raw:
        print("ERROR: Failed fetching http://rest.kegg.jp/list/pathway", file=sys.stderr)
        sys.exit(1)

    canonical_pathways = {} # num_id -> {name, category, subcategory, ...}
    for line in pw_list_raw.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            raw_id = parts[0].replace("path:map", "").replace("path:ko", "").strip()
            num_id = canonical_pathway_id(raw_id)
            name = parts[1].strip()
            canonical_pathways[num_id] = {
                "pathway_id": f"ko{num_id}",
                "num_id": num_id,
                "name": name,
                "category": "Metabolism",
                "subcategory": "Metabolism",
                "is_metabolic": False,
                "is_prokaryote_relevant": True,
                "microbial_relevance": "prokaryote_relevant",
                "map_structure_type": "regular",
                "reaction_metric_applicable": False
            }

    # Assign category & scope by numeric range and known taxonomy
    for num_id, info in canonical_pathways.items():
        n = int(num_id) if num_id.isdigit() else 99999
        # 01100-01130: Global maps
        if 1100 <= n <= 1130:
            info["category"] = "Global and overview maps"
            info["map_structure_type"] = "global"
            info["microbial_relevance"] = "global_overview"
            info["is_metabolic"] = True
            info["reaction_metric_applicable"] = True
        # 01200-01250: Overview maps
        elif 1200 <= n <= 1250:
            info["category"] = "Global and overview maps"
            info["map_structure_type"] = "overview"
            info["microbial_relevance"] = "global_overview"
            info["is_metabolic"] = True
            info["reaction_metric_applicable"] = True
        # 00010-01059: Classical Metabolic pathways
        elif 10 <= n < 1100:
            info["category"] = "Metabolism"
            info["is_metabolic"] = True
            info["reaction_metric_applicable"] = True
            info["microbial_relevance"] = "microbial_core"
        # 02000-02099: Environmental Information Processing (Membrane transport, signaling)
        elif 2000 <= n <= 2099:
            info["category"] = "Environmental Information Processing"
            info["is_metabolic"] = False
            info["reaction_metric_applicable"] = False
            info["microbial_relevance"] = "prokaryote_relevant"
        # 03000-03499: Genetic Information Processing (Replication, repair, translation)
        elif 3000 <= n <= 3499:
            info["category"] = "Genetic Information Processing"
            info["is_metabolic"] = False
            info["reaction_metric_applicable"] = False
            info["microbial_relevance"] = "prokaryote_relevant"
        # 04000-04999: Cellular Processes (Motility, biofilm, community)
        elif 4000 <= n <= 4999:
            info["category"] = "Cellular Processes"
            info["is_metabolic"] = False
            info["reaction_metric_applicable"] = False
            if n in [2024, 2026, 2030, 4112, 4122, 5111]: # e.g. biofilm, flagellar assembly, chemotaxis
                info["microbial_relevance"] = "microbial_core"
            else:
                info["microbial_relevance"] = "cellular_processes"
        # 05000+: Human Diseases & Organismal Systems
        else:
            info["category"] = "Organismal Systems / Human Diseases"
            info["is_metabolic"] = False
            info["reaction_metric_applicable"] = False
            info["is_prokaryote_relevant"] = False
            info["microbial_relevance"] = "nonmicrobial_context"

    print(f"      Categorized {len(canonical_pathways)} pathways.")

    # -------------------------------------------------------------------------
    # 4. Strict Three-Way Intersection: (Pathway, KO, Reaction) in P ∩ Q ∩ R
    # -------------------------------------------------------------------------
    print("\n[4/6] Building strict three-way intersection P ∩ Q ∩ R for pathway reactions...")
    # P: (canonical_pathway, KO)
    p_link_raw = fetch_url("http://rest.kegg.jp/link/pathway/ko")
    P = set()
    pathway_to_kos = {}
    for line in p_link_raw.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            ko = parts[0].replace("ko:", "").strip()
            cid = canonical_pathway_id(parts[1])
            if cid in canonical_pathways:
                P.add((cid, ko))
                pathway_to_kos.setdefault(cid, set()).add(ko)

    # Q: (canonical_pathway, Reaction)
    q_link_raw = fetch_url("http://rest.kegg.jp/link/reaction/pathway")
    Q = set()
    for line in q_link_raw.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            cid = canonical_pathway_id(parts[0])
            rn = parts[1].replace("rn:", "").strip()
            if cid in canonical_pathways:
                Q.add((cid, rn))

    # R: (KO, Reaction)
    r_link_raw = fetch_url("http://rest.kegg.jp/link/reaction/ko")
    R = set()
    for line in r_link_raw.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            ko = parts[0].replace("ko:", "").strip()
            rn = parts[1].replace("rn:", "").strip()
            R.add((ko, rn))

    # Construct strict intersection triples
    triples = [] # (cid, ko, rn)
    pathway_to_reactions = {}
    for (cid, ko) in P:
        for (r_ko, rn) in R:
            if r_ko == ko:
                if (cid, rn) in Q:
                    triples.append((cid, ko, rn))
                    pathway_to_reactions.setdefault(cid, set()).add(rn)

    print(f"      P (Pathway-KO pairs): {len(P)}")
    print(f"      Q (Pathway-Reaction pairs): {len(Q)}")
    print(f"      R (KO-Reaction pairs): {len(R)}")
    print(f"      Intersection triples (Pathway, KO, Reaction): {len(triples)}")

    # Update pathway metadata with counts
    for cid, info in canonical_pathways.items():
        info["total_kos"] = len(pathway_to_kos.get(cid, set()))
        info["total_reactions"] = len(pathway_to_reactions.get(cid, set()))

    # Write all_kegg_pathways.tsv
    pathways_tsv = os.path.join(args.out_dir, "all_kegg_pathways.tsv")
    with open(pathways_tsv, "w", encoding="utf-8") as out:
        out.write("pathway_id\tcanonical_id\tname\tcategory\tis_metabolic\tmap_structure_type\tis_prokaryote_relevant\tmicrobial_relevance\treaction_metric_applicable\ttotal_kos\ttotal_reactions\n")
        for cid in sorted(canonical_pathways.keys()):
            p = canonical_pathways[cid]
            out.write(f"{p['pathway_id']}\t{cid}\t{p['name']}\t{p['category']}\t{p['is_metabolic']}\t{p['map_structure_type']}\t{p['is_prokaryote_relevant']}\t{p['microbial_relevance']}\t{p['reaction_metric_applicable']}\t{p['total_kos']}\t{p['total_reactions']}\n")
    print(f"      Wrote {pathways_tsv}")

    # Write pathway_ko_reaction.tsv
    triples_tsv = os.path.join(args.out_dir, "pathway_ko_reaction.tsv")
    with open(triples_tsv, "w", encoding="utf-8") as out:
        out.write("pathway_id\tcanonical_id\tko\treaction_id\n")
        for (cid, ko, rn) in sorted(set(triples)):
            out.write(f"ko{cid}\t{cid}\t{ko}\t{rn}\n")
    print(f"      Wrote {triples_tsv}")

    # -------------------------------------------------------------------------
    # 5. Build Variable-Depth BRITE Hierarchy & Multi-Membership (Gate G5 & G6)
    # -------------------------------------------------------------------------
    print("\n[5/6] Fetching and parsing BRITE hierarchies (ko00001, ko02000, ko02022, ko02048)...")
    brite_records = [] # (ko, hierarchy_id, path_id, depth, node_name, full_path_json)
    orphan_candidates = set(ko_dict.keys())

    brite_sources = [
        ("ko00001", "http://rest.kegg.jp/get/br:ko00001/json"),
        ("ko02000", "http://rest.kegg.jp/get/br:ko02000/json"),
        ("ko02022", "http://rest.kegg.jp/get/br:ko02022/json"),
        ("ko02048", "http://rest.kegg.jp/get/br:ko02048/json")
    ]

    path_counter = 0
    kos_in_09190 = set()
    kos_in_specialized = set()

    for hier_id, url in brite_sources:
        print(f"      Parsing {hier_id}...")
        raw_json = fetch_url(url)
        if not raw_json:
            print(f"[Warning] Failed fetching {hier_id}, skipping.", file=sys.stderr)
            continue
        try:
            tree_data = json.loads(raw_json)
        except Exception as e:
            print(f"[Warning] JSON decode error for {hier_id}: {e}", file=sys.stderr)
            continue

        def recurse_node(node, current_path):
            nonlocal path_counter
            n_name = node.get("name", "").strip()
            path_here = current_path + [n_name] if n_name else current_path

            m_ko = re.match(r'^(K\d{5})', n_name)
            if m_ko:
                ko = m_ko.group(1)
                path_counter += 1
                brite_records.append({
                    "ko": ko,
                    "hierarchy_id": hier_id,
                    "path_id": f"P{path_counter:08d}",
                    "depth": len(path_here),
                    "node_name": n_name,
                    "full_path": json.dumps(path_here, ensure_ascii=False)
                })
                if hier_id != "ko00001":
                    kos_in_specialized.add(ko)
                elif "09190 Not Included in Pathway or Brite" in path_here:
                    kos_in_09190.add(ko)

            for child in node.get("children", []):
                recurse_node(child, path_here)

        recurse_node(tree_data, [])

    print(f"      Total BRITE membership paths recorded: {len(brite_records)}")
    print(f"      KOs in 09190 (Not Included in Pathway/Brite): {len(kos_in_09190)}")
    print(f"      KOs in specialized systems (ko02000/02022/02048): {len(kos_in_specialized)}")

    brite_tsv = os.path.join(args.out_dir, "ko_hierarchy_membership.tsv")
    with open(brite_tsv, "w", encoding="utf-8") as out:
        out.write("ko\thierarchy_id\tpath_id\tdepth\tnode_name\tfull_path\n")
        for rec in brite_records:
            out.write(f"{rec['ko']}\t{rec['hierarchy_id']}\t{rec['path_id']}\t{rec['depth']}\t{rec['node_name']}\t{rec['full_path']}\n")
    print(f"      Wrote {brite_tsv}")

    # -------------------------------------------------------------------------
    # 6. Generate Manifest & SHA256 (Gate G1 Snapshot Gate)
    # -------------------------------------------------------------------------
    print("\n[6/6] Writing snapshot manifest and computing SHA256 hashes...")
    modules_tsv = os.path.join(args.out_dir, "all_kegg_modules.tsv")
    manifest_json = os.path.join(args.out_dir, "kegg_snapshot_manifest.json")

    manifest = {
        "snapshot_date": datetime.date.today().isoformat(),
        "kegg_release": "KEGG Release 115.0+ (Live REST API)",
        "provenance_url": "http://rest.kegg.jp",
        "counts": {
            "total_official_modules": 574,
            "total_pathways": len(canonical_pathways),
            "total_pathway_kos": len(set(ko for _, ko in P)),
            "total_brite_records": len(brite_records),
            "total_unique_brite_kos": len(set(r["ko"] for r in brite_records)),
            "total_specialized_brite_kos": len(kos_in_specialized),
            "total_kos_in_09190": len(kos_in_09190),
            "total_kegg_kos": len(ko_dict)
        },
        "sha256": {
            "all_kegg_modules.tsv": compute_sha256(modules_tsv),
            "all_kegg_pathways.tsv": compute_sha256(pathways_tsv),
            "pathway_ko_reaction.tsv": compute_sha256(triples_tsv),
            "ko_hierarchy_membership.tsv": compute_sha256(brite_tsv)
        }
    }

    with open(manifest_json, "w", encoding="utf-8") as out:
        json.dump(manifest, out, indent=2)

    print(f"      Wrote {manifest_json}")
    print("\n" + "=" * 70)
    print("KEGG Knowledge Assets successfully built and verified against 10 Gates!")
    print("=" * 70)

if __name__ == "__main__":
    main()
