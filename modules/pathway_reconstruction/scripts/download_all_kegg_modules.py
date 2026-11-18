#!/usr/bin/env python3
"""
download_all_kegg_modules.py
Fetches the complete catalog of all KEGG Modules (~585 modules) from the official KEGG REST API:
1. Calls https://rest.kegg.jp/list/module to retrieve all module IDs and names.
2. Batch-fetches module entries via https://rest.kegg.jp/get/M00001+M00002+... (10 per batch).
3. Parses DEFINITION, CLASS, and NAME fields.
4. Outputs knowledge/kegg/all_kegg_modules.tsv.
"""

import argparse
import os
import re
import sys
import time
import urllib.request

def fetch_url(url, retries=3, delay=1.0):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "GenomeAnno-KEGG-Fetcher/1.0"})
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8", errors="ignore")
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
            else:
                print(f"Warning: Failed fetching {url}: {e}", file=sys.stderr)
                return ""

def main():
    parser = argparse.ArgumentParser(description="Download all official KEGG Module definitions.")
    parser.add_argument("-o", "--output", default="knowledge/kegg/all_kegg_modules.tsv", help="Output TSV path")
    parser.add_argument("--batch_size", type=int, default=10, help="Modules per GET request (max 10 for KEGG API)")

    args = parser.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    print("Step 1: Fetching list of all KEGG modules from https://rest.kegg.jp/list/module...")
    list_content = fetch_url("https://rest.kegg.jp/list/module")
    if not list_content:
        print("ERROR: Could not retrieve module list from KEGG REST API.", file=sys.stderr)
        sys.exit(1)

    all_modules = []
    for line in list_content.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            m_id = parts[0].replace("md:", "").strip()
            m_name = parts[1].strip()
            if m_id.startswith("M"):
                all_modules.append((m_id, m_name))

    print(f"Found {len(all_modules)} official KEGG modules.")

    # Batch fetch module definitions
    print(f"Step 2: Batch-fetching definitions in chunks of {args.batch_size}...")
    module_definitions = {} # m_id -> {name, definition, category, subcategory}

    for i in range(0, len(all_modules), args.batch_size):
        chunk = all_modules[i:i + args.batch_size]
        chunk_ids = [m[0] for m in chunk]
        query_url = "https://rest.kegg.jp/get/" + "+".join(chunk_ids)

        data = fetch_url(query_url)
        if not data:
            continue

        # Split multi-entry text (delimited by ///)
        entries = data.split("///")
        for entry in entries:
            entry = entry.strip()
            if not entry:
                continue

            # Parse entry
            m_id = None
            name = None
            definition = None
            category = "Metabolism"
            subcategory = "Pathway"

            for line in entry.split("\n"):
                if line.startswith("ENTRY"):
                    parts = line.split()
                    if len(parts) >= 2:
                        m_id = parts[1].strip()
                elif line.startswith("NAME"):
                    name = line[12:].strip()
                elif line.startswith("DEFINITION"):
                    definition = line[12:].strip()
                elif line.startswith("CLASS"):
                    cls_parts = [c.strip() for c in line[12:].split(";")]
                    if len(cls_parts) >= 2:
                        category = cls_parts[1]
                    if len(cls_parts) >= 3:
                        subcategory = cls_parts[2]

            if m_id and definition:
                module_definitions[m_id] = {
                    "name": name or m_id,
                    "definition": definition,
                    "category": category,
                    "subcategory": subcategory
                }

        print(f"  Processed {min(i + args.batch_size, len(all_modules))}/{len(all_modules)} modules...")
        time.sleep(0.2) # Friendly rate limiting

    print(f"Successfully retrieved {len(module_definitions)} module definitions.")

    # Write unified TSV
    with open(args.output, "w", encoding="utf-8") as out:
        out.write("module_id\tname\tcategory\tsubcategory\tdefinition\n")
        for m_id, info in sorted(module_definitions.items()):
            out.write(f"{m_id}\t{info['name']}\t{info['category']}\t{info['subcategory']}\t{info['definition']}\n")

    print(f"Complete KEGG module catalog written to {args.output}")

if __name__ == "__main__":
    main()
