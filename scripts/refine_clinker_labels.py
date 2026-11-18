#!/usr/bin/env python3
import sys
import os
import json
import re

def clean_id_spaces(val):
    """Removes GenBank 80-char line-wrapping whitespace artifacts ONLY from structured identifiers."""
    if not val or not isinstance(val, str):
        return val
    if ":" in val or ";" in val or "(" in val or "[EC:" in val:
        return val
    if not re.search(r"\s", val):
        return val
    no_space = re.sub(r"\s+", "", val)
    if re.match(r"^[a-zA-Z0-9_.\-]+$", no_space):
        if ("___" in val or
            re.search(r"GC[AF]_", no_space) or
            re.search(r"\bprot\d+", no_space, re.I) or
            re.search(r"\bgene\d+", no_space, re.I) or
            re.match(r"^[A-Z]{4,6}\d+(?:\.\d+)?$", no_space) or
            re.match(r"^[a-zA-Z0-9]+_\d+$", no_space)):
            return no_space
    return val

def extract_metric_str_from_note(note_str):
    """Extracts score, threshold, evalue from note and formats as '(sc: score/threshold, e: evalue)'."""
    if not note_str:
        return ""
    m_sc = re.search(r"score=([^\s|]+)", note_str)
    m_th = re.search(r"threshold=([^\s|]+)", note_str)
    m_ev = re.search(r"evalue=([^\s|]+)", note_str)
    sc = m_sc.group(1) if m_sc else ""
    th = m_th.group(1) if m_th else ""
    ev = m_ev.group(1) if m_ev else ""

    parts = []
    sc_s, th_s = "", ""
    if sc:
        try:
            v = float(sc)
            sc_s = f"{v:.1f}" if v % 1 != 0 else f"{int(v)}"
        except Exception:
            sc_s = sc
    if th:
        try:
            v = float(th)
            th_s = f"{v:.1f}" if v % 1 != 0 else f"{int(v)}"
        except Exception:
            th_s = th
    if sc_s and th_s:
        parts.append(f"sc: {sc_s}/{th_s}")
    elif sc_s:
        parts.append(f"sc: {sc_s}")

    if ev:
        try:
            v = float(ev)
            if v == 0.0:
                parts.append("e: 0.0")
            elif v < 0.001 or v >= 10000:
                parts.append(f"e: {v:.1e}")
            else:
                parts.append(f"e: {v:.3f}")
        except Exception:
            parts.append(f"e: {ev}")

    joined = ", ".join(parts)
    return f" ({joined})" if parts else ""

def extract_contig_from_genes(genes):
    """Extracts genomic contig ID from qualifiers or structured IDs of genes in locus."""
    for g in genes:
        names = g.get("names", {})
        orig_c = names.get("orig_coord", "")
        if orig_c and ":" in orig_c:
            return orig_c.split(":")[0].strip()
        note = names.get("note", "")
        m_c = re.search(r"orig_coord=([^:]+):", note)
        if m_c:
            return m_c.group(1).strip()
        for k in ["locus_tag", "protein_id", "orig_protein_id"]:
            val = names.get(k, "")
            if "___" in val:
                parts = val.split("___")
                if len(parts) >= 2 and parts[1]:
                    return parts[1].strip()
    return ""

def refine_clinker_html(html_path):
    """
    Post-processes Clinker HTML file:
    1. Annotates track/cluster headers with the matched contig ID (e.g. 'GCA_964582085.1 (CAZOPV010000010.1)').
    2. Eliminates GenBank 80-column line-wrapping space corruption in identifiers (e.g. 'CA ZOPV' -> 'CAZOPV').
    3. Sets each gene's primary label to its informative KO / gene / product description with explicit
       quality metrics (sc: score/threshold, e: evalue) for both target and reference genomes.
    """
    if not html_path or not os.path.exists(html_path):
        return
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()

    match = re.search(r"const data\s*=\s*(\{.*?\});", content, re.DOTALL)
    if not match:
        return

    try:
        data = json.loads(match.group(1))
    except Exception as e:
        print(f"Notice: could not parse clinker data JSON ({e})", file=sys.stderr)
        return

    updated_count = 0
    for cl in data.get("clusters", []):
        cl_name = cl.get("name", "")
        all_cluster_genes = []
        for lc in cl.get("loci", []):
            all_cluster_genes.extend(lc.get("genes", []))

        # 1. Track / Cluster Contig ID Annotation
        ctg_id = extract_contig_from_genes(all_cluster_genes)
        is_target = ("target" in cl_name.lower()) or (cl_name.endswith("_locus") and not re.search(r"GC[AF]_", cl_name))

        if is_target:
            clean_target = re.sub(r'(_target)?_locus$', '', cl_name)
            base_title = clean_target if clean_target else "Target"
        else:
            acc_match = re.search(r"(GC[AF]_\d+(?:[._]\d+)?)", cl_name)
            if not acc_match and all_cluster_genes:
                first_lt = all_cluster_genes[0].get("names", {}).get("locus_tag", "")
                acc_match = re.search(r"(GC[AF]_\d+(?:[._]\d+)?)", first_lt)
            if not acc_match:
                acc_match = re.search(r"(GC[AF]_\d+(?:[._]\d+)?)", os.path.basename(html_path))

            base_title = acc_match.group(1) if acc_match else cl_name
            base_title = re.sub(r'(_target)?_locus$', '', base_title)
            if re.match(r"^GC[AF]_\d+_\d+$", base_title):
                base_title = base_title.rsplit('_', 1)[0] + '.' + base_title.rsplit('_', 1)[1]

        # Check for module ID in path or cluster name
        mod_m = re.search(r"\b(M0015[789]|M00160)\b", html_path) or re.search(r"\b(M0015[789]|M00160)\b", cl_name)
        mod_str = f" [{mod_m.group(1)}]" if mod_m and mod_m.group(1) not in base_title else ""

        if ctg_id and ctg_id not in cl.get("name", ""):
            cl["name"] = f"{base_title}{mod_str} ({ctg_id})"
            for lc in cl.get("loci", []):
                s = lc.get("start", 0)
                e = lc.get("end", 0)
                lc["name"] = f"{ctg_id}:{s}-{e}"
        elif base_title:
            cl["name"] = f"{base_title}{mod_str}"

        # 2. Gene Label Refinement
        for lc in cl.get("loci", []):
            for g in lc.get("genes", []):
                names = g.get("names", {})
                # Clean identifier spaces in names dictionary
                for k in ["locus_tag", "protein_id", "orig_protein_id"]:
                    if k in names and names[k]:
                        names[k] = clean_id_spaces(names[k])
                if "gene" in names and names["gene"] and "___" in names["gene"]:
                    names["gene"] = clean_id_spaces(names["gene"])

                lbl = names.get("label", "")
                lt = names.get("locus_tag", "")
                gn = names.get("gene", "")
                prod = names.get("product", "")
                note = names.get("note", "")

                clean_lt = lt.split("___")[-1] if "___" in lt else lt
                m_metric = extract_metric_str_from_note(note)

                # Priority 1: Rich label in /gene qualifier (from format_gene_qualifiers)
                if gn and (gn.startswith("K") or ":" in gn) and gn != clean_lt and gn != lt:
                    if "(sc:" not in gn and m_metric:
                        g["label"] = f"{gn}{m_metric}"
                    else:
                        g["label"] = gn
                    updated_count += 1
                # Priority 2: Structured KO annotation present in label qualifier
                elif lbl and (lbl != lt and lbl != clean_lt):
                    if "(sc:" not in lbl and m_metric:
                        g["label"] = f"{lbl}{m_metric}"
                    else:
                        g["label"] = lbl
                    updated_count += 1
                # Priority 3: Reconstruct from note KO + gene / product if label wasn't refined
                else:
                    m_ko = re.search(r"KO:\s*(K\d{5})", note)
                    if m_ko:
                        ko_id = m_ko.group(1)
                        if gn and prod and prod != "hypothetical protein" and gn != clean_lt and gn != lt and gn != ko_id:
                            base_lbl = f"{ko_id}: {gn}; {prod}"
                        elif prod and prod != "hypothetical protein":
                            base_lbl = f"{ko_id}: {prod}"
                        elif gn and gn != clean_lt and gn != lt and gn != ko_id:
                            base_lbl = f"{ko_id}: {gn}"
                        else:
                            base_lbl = f"{ko_id}: {clean_lt}"
                        g["label"] = f"{base_lbl}{m_metric}"
                        updated_count += 1
                    # Priority 4: Existing g["label"] has KO prefix or colon
                    elif g.get("label") and (g["label"].startswith("K") or ":" in g["label"]):
                        if "(sc:" not in g["label"] and m_metric:
                            g["label"] = f"{g['label']}{m_metric}"
                        updated_count += 1
                    # Priority 5: Informative gene symbol + product
                    elif gn and gn != lt and gn != clean_lt:
                        base_lbl = f"{gn}: {prod}" if (prod and prod != "hypothetical protein") else gn
                        g["label"] = f"{base_lbl}{m_metric}"
                        updated_count += 1
                    # Priority 6: Product alone
                    elif prod and prod != "hypothetical protein":
                        g["label"] = f"{prod}{m_metric}"
                        updated_count += 1
                    # Priority 7: Fallback to clean locus tag without contig prefix
                    elif clean_lt and g.get("label") != clean_lt:
                        g["label"] = f"{clean_lt}{m_metric}"
                        updated_count += 1

                # Clean any lingering wrapping spaces in final label ONLY if purely an identifier
                if "label" in g and isinstance(g["label"], str):
                    if ":" not in g["label"] and ";" not in g["label"]:
                        g["label"] = clean_id_spaces(g["label"])

    new_json = json.dumps(data)
    new_content = content[:match.start(1)] + new_json + content[match.end(1):]
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(new_content)
    print(f"Successfully refined Clinker track headers & {updated_count} gene labels in: {html_path}", flush=True)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        refine_clinker_html(sys.argv[1])
