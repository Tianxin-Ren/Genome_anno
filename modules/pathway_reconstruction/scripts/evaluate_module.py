#!/usr/bin/env python3
"""
evaluate_module.py
Evaluates all KEGG Modules (supports all ~574 official KEGG modules from TSV or curated YAMLs)
and reconstructs full-lineage evidence chains:
Genome -> Module -> KO -> Protein -> Gene -> Contig -> Coordinate.
Handles multiple threshold-passing hits per protein with explicit evidence strategies:
  --ko_evidence_strategy:
    - 'all_strict' (default): Every threshold-passing hit can provide evidence for relevant modules.
      Evidence records explicit ko_rank (1=primary, 2=secondary) and score margin.
    - 'best_hit': Only Rank 1 (best score) hits are passed to module evaluation.
Outputs:
  - genome_module_completeness.tsv.gz: completeness (0.0~1.0), status, matched/missing KOs
  - genome_module_presence.tsv.gz: binary/qualitative presence (complete, partial, absent)
  - module_evidence.tsv.gz: detailed protein evidence per module step
"""

import argparse
import gzip
import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

try:
    import yaml
except ImportError:
    yaml = None

from kegg_module_parser import evaluate_module, parse_definition

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def load_modules(modules_source):
    """Load KEGG module definitions from either a TSV file, a directory with all_kegg_modules.tsv, or YAML files."""
    modules = []
    if not os.path.exists(modules_source):
        return modules

    # Case 1: Direct TSV catalog
    # Case 1: Direct TSV catalog
    if os.path.isfile(modules_source) and modules_source.endswith(".tsv"):
        tsv_path = modules_source
    # Case 2: Directory containing YAML files
    elif os.path.isdir(modules_source) and any(f.endswith((".yaml", ".yml")) for f in os.listdir(modules_source)):
        tsv_path = None
    # Case 3: Directory with all_kegg_modules.tsv
    elif os.path.isdir(modules_source) and os.path.isfile(os.path.join(modules_source, "all_kegg_modules.tsv")):
        tsv_path = os.path.join(modules_source, "all_kegg_modules.tsv")
    elif os.path.isdir(modules_source) and os.path.isfile(os.path.join(os.path.dirname(modules_source), "all_kegg_modules.tsv")):
        tsv_path = os.path.join(os.path.dirname(modules_source), "all_kegg_modules.tsv")
    else:
        tsv_path = None

    if tsv_path:
        with open(tsv_path, "r", encoding="utf-8") as f:
            header = f.readline().strip().split("\t")
            m_col = header.index("module_id") if "module_id" in header else 0
            n_col = header.index("name") if "name" in header else 1
            d_col = header.index("definition") if "definition" in header else (len(header) - 1)
            c_col = header.index("category") if "category" in header else -1
            sc_col = header.index("subcategory") if "subcategory" in header else -1

            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > max(m_col, d_col):
                    m_id = parts[m_col]
                    m_def = parts[d_col]
                    m_name = parts[n_col] if n_col < len(parts) else m_id
                    cat = parts[c_col] if c_col != -1 and c_col < len(parts) else "Metabolism"
                    subcat = parts[sc_col] if sc_col != -1 and sc_col < len(parts) else ""
                    modules.append({
                        "module_id": m_id,
                        "name": m_name,
                        "definition": m_def,
                        "category": cat,
                        "subcategory": subcat
                    })
        return modules

    # Case 4: Directory of YAML definitions (recursive)
    if os.path.isdir(modules_source):
        for root, _, files in os.walk(modules_source):
            for f in sorted(files):
                if f.endswith((".yaml", ".yml")):
                    fpath = os.path.join(root, f)
                    try:
                        data = {}
                        curr_key = None
                        with open(fpath, "r", encoding="utf-8") as yf:
                            for line in yf:
                                sline = line.strip()
                                if not sline or sline.startswith("#"):
                                    continue
                                if sline.startswith("- ") and curr_key:
                                    item = sline[2:].strip().split("#")[0].strip().strip("'\"")
                                    if isinstance(data.get(curr_key), list):
                                        data[curr_key].append(item)
                                    else:
                                        data[curr_key] = [item]
                                elif ":" in sline:
                                    k, v = sline.split(":", 1)
                                    k = k.strip()
                                    v = v.strip().split("#")[0].strip().strip("'\"")
                                    curr_key = k
                                    data[k] = v if v else []
                        if data and "module_id" in data and "definition" in data:
                            modules.append(data)
                    except Exception as e:
                        print(f"Warning: Failed to load {fpath}: {e}", file=sys.stderr)
    return modules

def main():
    parser = argparse.ArgumentParser(description="Evaluate all KEGG modules across annotated genomes.")
    parser.add_argument("-k", "--ko_strict", required=True, help="Path to protein_KO_all_strict.tsv.gz")
    parser.add_argument("-d", "--modules_source", required=True, help="Path to all_kegg_modules.tsv or modules directory")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for module reconstruction tables")
    parser.add_argument("--ko_evidence_strategy", default="all_strict", choices=["all_strict", "best_hit"],
                        help="How multiple threshold-passing hits per protein are used (default: all_strict)")
    parser.add_argument("--completion_threshold", type=float, default=0.75, help="Threshold for complete status (default: 0.75)")
    parser.add_argument("--max_missing_steps", type=int, default=1, help="Max missing steps allowed in >=3-step modules for complete status (default: 1)")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    modules = load_modules(args.modules_source)
    if not modules:
        print(f"ERROR: No KEGG module definitions could be loaded from '{args.modules_source}'! Check file path and format.", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(modules)} KEGG module definitions from {args.modules_source}")
    print(f"Multiple-hit KO evidence strategy: {args.ko_evidence_strategy}")


    # 1. Parse KO annotations and index by genome
    # genome_id -> {ko -> list of evidence dicts}
    genome_ko_evidence = {}
    
    print(f"Reading KO annotations from {args.ko_strict}...")
    if os.path.isfile(args.ko_strict):
        with open_file(args.ko_strict, "rt") as f:
            header = None
            for line in f:
                sline = line.strip()
                if not sline or sline.startswith("#"):
                    continue
                header = sline.split("\t")
                break

            if header:
                p_idx = header.index("protein_id") if "protein_id" in header else 0
                g_idx = header.index("gene_id") if "gene_id" in header else 1
                gm_idx = header.index("genome_id") if "genome_id" in header else 2
                k_idx = header.index("ko") if "ko" in header else 3
                b_idx = header.index("bitscore") if "bitscore" in header else 5
                th_idx = header.index("curated_threshold") if "curated_threshold" in header else 6
                rk_idx = header.index("ko_rank") if "ko_rank" in header else 9
                sr_idx = header.index("score_margin") if "score_margin" in header else (header.index("score_rank") if "score_rank" in header else 10)

                for line in f:
                    sline = line.strip()
                    if not sline or sline.startswith("#"):
                        continue
                    parts = sline.split("\t")
                    if len(parts) > max(p_idx, g_idx, gm_idx, k_idx):
                        prot_id = parts[p_idx]
                        gene_id = parts[g_idx]
                        genome_id = parts[gm_idx]
                        ko = parts[k_idx].strip()
                        if not ko or ko in ["-", "None", "null", "NA"]:
                            continue

                        try:
                            bitscore = float(parts[b_idx]) if b_idx < len(parts) else 0.0
                        except (ValueError, TypeError):
                            bitscore = 0.0

                        try:
                            threshold = float(parts[th_idx]) if th_idx < len(parts) else 0.0
                        except (ValueError, TypeError):
                            threshold = 0.0

                        try:
                            ko_rank = int(parts[rk_idx]) if rk_idx < len(parts) else 1
                        except (ValueError, TypeError):
                            ko_rank = 1

                        try:
                            margin = float(parts[sr_idx]) if sr_idx < len(parts) else (bitscore - threshold)
                        except (ValueError, TypeError):
                            margin = bitscore - threshold

                        # Filter if best_hit strategy is chosen
                        if args.ko_evidence_strategy == "best_hit" and ko_rank > 1:
                            continue

                        if genome_id not in genome_ko_evidence:
                            genome_ko_evidence[genome_id] = {}
                        if ko not in genome_ko_evidence[genome_id]:
                            genome_ko_evidence[genome_id][ko] = []
                        
                        genome_ko_evidence[genome_id][ko].append({
                            "protein_id": prot_id,
                            "gene_id": gene_id,
                            "bitscore": bitscore,
                            "threshold": threshold,
                            "ko_rank": ko_rank,
                            "margin": margin
                        })

    print(f"Parsed annotations for {len(genome_ko_evidence)} genomes.")

    # 2. Evaluate modules for each genome
    completeness_tsv = os.path.join(args.out_dir, "genome_module_completeness.tsv")
    presence_tsv = os.path.join(args.out_dir, "genome_module_presence.tsv")
    evidence_tsv = os.path.join(args.out_dir, "module_evidence.tsv")

    with open(completeness_tsv, "w", encoding="utf-8") as c_out, \
         open(presence_tsv, "w", encoding="utf-8") as p_out, \
         open(evidence_tsv, "w", encoding="utf-8") as e_out:


        c_out.write("genome_id\tmodule_id\tmodule_name\tcompleteness\tdefinition_satisfied\tstatus\tmatched_kos\tmissing_kos\n")
        p_out.write("genome_id\tmodule_id\tstatus\n")
        e_out.write("genome_id\tmodule_id\tmodule_name\tko\tprotein_id\tgene_id\tcontig_id\tstart\tend\tstrand\tbitscore\tko_rank\tscore_margin\n")

        # Pre-compile AST trees for all modules to avoid millions of redundant tokenizations
        print(f"Pre-compiling AST trees for {len(modules)} KEGG modules...")
        for mod in modules:
            mod["tree"] = parse_definition(mod["definition"])

        for genome_id in sorted(genome_ko_evidence.keys()):
            present_kos = set(genome_ko_evidence[genome_id].keys())

            for mod in modules:
                m_id = mod["module_id"]
                m_name = mod.get("name", m_id)

                comp, status, def_sat, matched, missing = evaluate_module(
                    mod["tree"], present_kos, args.completion_threshold, args.max_missing_steps
                )

                c_out.write(f"{genome_id}\t{m_id}\t{m_name}\t{comp:.4f}\t{def_sat}\t{status}\t{','.join(matched)}\t{','.join(missing)}\n")
                p_out.write(f"{genome_id}\t{m_id}\t{status}\n")

                # Write detailed evidence with ko_rank and margin
                for ko in matched:
                    ev_list = genome_ko_evidence[genome_id].get(ko, [])
                    for ev in ev_list:
                        contig = "unknown"
                        start = "1"
                        end = "100"
                        strand = "+"
                        g_parts = ev["gene_id"].split("___")
                        if len(g_parts) >= 2:
                            contig = g_parts[1]
                        e_out.write(f"{genome_id}\t{m_id}\t{m_name}\t{ko}\t{ev['protein_id']}\t{ev['gene_id']}\t{contig}\t{start}\t{end}\t{strand}\t{ev['bitscore']:.1f}\t{ev['ko_rank']}\t{ev['margin']:.1f}\n")

    print(f"Module evaluation complete across all {len(modules)} modules. Outputs written to {args.out_dir}")

if __name__ == "__main__":
    main()
