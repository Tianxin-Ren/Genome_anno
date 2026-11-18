#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import gzip
import json
import re
import glob
import math
import pandas as pd
try:
    from Bio import SeqIO, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    SeqIO = None
    SeqFeature = None
    SeqRecord = None
    Seq = None

DEFAULT_KO_DEFINITIONS = {
    # M00159 (Prokaryotic V/A-type)
    "K02117": "atpA; V/A-type H+/Na+-transporting ATPase subunit A",
    "K02118": "atpB; V/A-type H+/Na+-transporting ATPase subunit B",
    "K02119": "atpC; V/A-type H+/Na+-transporting ATPase subunit C",
    "K02120": "atpD; V/A-type H+/Na+-transporting ATPase subunit D",
    "K02121": "atpE; V/A-type H+/Na+-transporting ATPase subunit E",
    "K02122": "atpF; V/A-type H+/Na+-transporting ATPase subunit F",
    "K02107": "atpH; V/A-type H+/Na+-transporting ATPase subunit G/H",
    "K02123": "atpI; V/A-type H+/Na+-transporting ATPase subunit I",
    "K02124": "atpK; V/A-type H+/Na+-transporting ATPase subunit K",
    # M00157 (Prokaryotic F-type)
    "K02111": "atpA; F-type H+/Na+-transporting ATPase subunit alpha",
    "K02112": "atpB; F-type H+/Na+-transporting ATPase subunit a",
    "K02113": "atpC; F-type H+/Na+-transporting ATPase subunit gamma",
    "K02114": "atpD; F-type H+/Na+-transporting ATPase subunit delta",
    "K02115": "atpE; F-type H+/Na+-transporting ATPase subunit c",
    "K02108": "atpF; F-type H+/Na+-transporting ATPase subunit b",
    "K02109": "atpG; F-type H+/Na+-transporting ATPase subunit beta",
    "K02110": "atpH; F-type H+/Na+-transporting ATPase subunit epsilon",
    # M00158 (Eukaryotic F-type)
    "K02126": "ATP5F1A; ATP synthase F1 subunit alpha",
    "K02127": "ATP5PB; ATP synthase F(0) complex subunit b",
    "K02128": "ATP5F1C; ATP synthase F1 subunit gamma",
    "K02129": "ATP5F1D; ATP synthase F1 subunit delta",
    "K02130": "ATP5MC; ATP synthase F(0) complex subunit c",
    "K02131": "ATP5PD; ATP synthase F(0) complex subunit d",
    "K02132": "ATP5PF; ATP synthase F(0) complex subunit F6",
    "K02133": "ATP5ME; ATP synthase F(0) complex subunit e",
    "K02134": "ATP5F1B; ATP synthase F1 subunit beta",
    "K02135": "ATP5F1E; ATP synthase F1 subunit epsilon",
    "K02136": "ATP5PO; ATP synthase F1 subunit OSCP",
    "K02137": "ATP6; ATP synthase F(0) complex subunit a",
    "K02138": "ATP8; ATP synthase F(0) complex subunit 8",
    # M00160 (Eukaryotic V-type)
    "K02145": "ATP6V1A; V-type H+-transporting ATPase subunit A",
    "K02146": "ATP6V0A; V-type H+-transporting ATPase subunit a",
    "K02147": "ATP6V1B; V-type H+-transporting ATPase subunit B",
    "K02148": "ATP6V0B; V-type H+-transporting ATPase subunit b",
    "K02149": "ATP6V1C; V-type H+-transporting ATPase subunit C",
    "K02150": "ATP6V0C; V-type H+-transporting ATPase subunit c",
    "K02151": "ATP6V1D; V-type H+-transporting ATPase subunit D",
    "K02152": "ATP6V1E; V-type H+-transporting ATPase subunit E",
    "K02153": "ATP6V1F; V-type H+-transporting ATPase subunit F",
    "K02154": "ATP6V1G; V-type H+-transporting ATPase subunit G",
    "K02155": "ATP6V1H; V-type H+-transporting ATPase subunit H",
    # M00595 (SOX system, thiosulfate oxidation)
    "K17224": "soxB; S-sulfosulfanyl-L-cysteine sulfohydrolase",
    "K17222": "soxA; L-cysteine S-thiosulfotransferase",
    "K17223": "soxX; L-cysteine S-thiosulfotransferase",
    "K17227": "soxZ; sulfur-oxidizing protein SoxZ",
    "K17226": "soxY; sulfur-oxidizing protein SoxY",
    "K17225": "soxC; sulfane dehydrogenase subunit SoxC",
    "K22622": "soxD; S-disulfanyl-L-cysteine oxidoreductase SoxD",
}

MODULE_TARGET_KOS = {
    "M00595": {"K17224", "K17222", "K17223", "K17227", "K17226", "K17225", "K22622"},
    "M00159": {"K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02107", "K02123", "K02124"},
    "M00157": {"K02111", "K02112", "K02113", "K02114", "K02115", "K02108", "K02109", "K02110"},
    "M00158": {"K02126", "K02127", "K02128", "K02129", "K02130", "K02131", "K02132", "K02133", "K02134", "K02135", "K02136", "K02137", "K02138"},
    "M00160": {"K02145", "K02146", "K02147", "K02148", "K02149", "K02150", "K02151", "K02152", "K02153", "K02154", "K02155"},
}

def load_fasta_records(path):
    """Safely loads FASTA records whether uncompressed or gzipped."""
    if not path or not os.path.exists(path):
        return []
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as h:
            return list(SeqIO.parse(h, "fasta"))
    return list(SeqIO.parse(path, "fasta"))

CODON_MAP = {
    'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT',
    'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA', 'R': 'CGT',
    'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT',
    '*': 'TAA'
}

def amino_acid_to_dna(aa_seq):
    """Converts amino acid sequence to pseudo-coding DNA."""
    return "".join(CODON_MAP.get(aa.upper(), 'NNN') for aa in str(aa_seq))

def load_kegg_ko_definitions():
    """Loads full KEGG KO definitions dictionary with fallback."""
    defs = DEFAULT_KO_DEFINITIONS.copy()
    cand_paths = [
        os.path.join(os.path.dirname(__file__), "cache/kegg_ko_cache.json"),
        os.path.join(os.path.dirname(__file__), "../scripts/cache/kegg_ko_cache.json"),
        "scripts/cache/kegg_ko_cache.json"
    ]
    for p in cand_paths:
        if os.path.exists(p):
            try:
                with open(p) as f:
                    k_data = json.load(f).get("data", {})
                    defs.update(k_data)
                break
            except Exception:
                pass
    return defs

def format_metric_str(score, threshold, evalue):
    """Formats (sc: score/threshold, e: evalue) or (sc: score, e: evalue) concisely."""
    parts = []
    sc_str = ""
    if score not in [None, "", "None"]:
        try:
            sc_val = float(score)
            sc_str = f"{sc_val:.1f}" if sc_val % 1 != 0 else f"{int(sc_val)}"
        except Exception:
            sc_str = str(score).strip()
    th_str = ""
    if threshold not in [None, "", "None"]:
        try:
            th_val = float(threshold)
            th_str = f"{th_val:.1f}" if th_val % 1 != 0 else f"{int(th_val)}"
        except Exception:
            th_str = str(threshold).strip()

    if sc_str and th_str:
        parts.append(f"sc: {sc_str}/{th_str}")
    elif sc_str:
        parts.append(f"sc: {sc_str}")

    if evalue not in [None, "", "None"]:
        try:
            ev_val = float(evalue)
            if ev_val == 0.0:
                parts.append("e: 0.0")
            elif ev_val < 0.001 or ev_val >= 10000:
                parts.append(f"e: {ev_val:.1e}")
            else:
                parts.append(f"e: {ev_val:.3f}")
        except Exception:
            parts.append(f"e: {str(evalue).strip()}")

    return f" ({', '.join(parts)})" if parts else ""

def fix_gbk_wrapped_qualifiers(gbk_path):
    """
    Biopython SeqIO.write wraps GenBank lines at 79 columns.
    When long identifiers like /locus_tag="..." or /protein_id="..." are wrapped across multiple lines,
    standard GenBank parsers (including Clinker) join continuation lines with a space (' '),
    corrupting identifiers (e.g. '...___CA ZOPV...').
    This helper post-processes the GBK file to unwrap multi-line quoted identifiers into single lines.
    """
    if not gbk_path or not os.path.exists(gbk_path):
        return
    try:
        with open(gbk_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        fixed_lines = []
        target_quals = {"locus_tag", "protein_id", "orig_protein_id", "gene", "label"}
        in_unwrappable = False

        for line in lines:
            line_str = line.rstrip("\r\n")
            stripped = line_str.strip()
            m_qual = re.match(r'^\s+/([a-zA-Z0-9_]+)="(.*)$', line_str)
            if m_qual:
                q_name = m_qual.group(1)
                q_val = m_qual.group(2)
                is_closed = q_val.endswith('"') and not q_val.endswith('\\"')
                if q_name in target_quals and not is_closed:
                    fixed_lines.append(line_str)
                    in_unwrappable = True
                else:
                    fixed_lines.append(line_str)
                    in_unwrappable = False
            elif in_unwrappable:
                prev = fixed_lines.pop()
                fixed_lines.append(prev + stripped)
                if stripped.endswith('"') and not stripped.endswith('\\"'):
                    in_unwrappable = False
            else:
                fixed_lines.append(line_str)

        with open(gbk_path, "w", encoding="utf-8") as f:
            f.write("\n".join(fixed_lines) + "\n")
    except Exception as e:
        print(f"Notice: fix_gbk_wrapped_qualifiers encountered error on {gbk_path}: {e}", file=sys.stderr)

def format_gene_qualifiers(p, is_op, ko, ko_definitions, default_def="", c_inf=None, hit=None):
    """Constructs clean, non-redundant label, gene, product, locus_tag and note qualifiers for GenBank and Clinker."""
    coord_str = ""
    if c_inf:
        ctg = c_inf.get("contig", "")
        s = c_inf.get("start")
        e = c_inf.get("end")
        st = c_inf.get("strand", "+")
        if s is not None and e is not None:
            coord_str = f"{ctg}:{s}..{e}({st})"

    if ko:
        ko_def = ko_definitions.get(ko, default_def) or default_def
        if ko_def:
            parts = ko_def.split(";", 1)
            sym = parts[0].split(",")[-1].strip()
            desc = parts[1].strip() if len(parts) > 1 else sym
            clean_desc = re.sub(r"\[EC:[^\]]+\]", "", desc).strip()
            if clean_desc and clean_desc != sym and sym != ko:
                base_label = f"{ko}: {sym}; {clean_desc}"
            else:
                base_label = f"{ko}: {clean_desc or sym}"
            gene_val = sym
            prod_val = desc
        else:
            base_label = f"{ko}: Target operon protein {ko}" if is_op else f"{ko}: Ref protein {ko}"
            gene_val = ko
            prod_val = f"Target operon protein {ko}" if is_op else f"Reference protein {ko}"

        sc = hit.get("score", "") if hit else ""
        th = hit.get("threshold", "") if hit else ""
        ev = hit.get("evalue", "") if hit else ""
        m_str = format_metric_str(sc, th, ev)
        label_val = f"{base_label}{m_str}"

        note_items = [f"KO: {ko}", f"orig_protein_id={p.id}"]
        if coord_str:
            note_items.append(f"orig_coord={coord_str}")
        if sc or th or ev:
            sc_part = f"score={sc}" if sc else ""
            th_part = f"threshold={th}" if th else ""
            ev_part = f"evalue={ev}" if ev else ""
            note_items.append(" ".join(filter(None, [sc_part, th_part, ev_part])))

        qualifiers = {
            "locus_tag": [p.id],
            "protein_id": [p.id],
            "gene": [label_val],
            "orig_gene": [gene_val],
            "product": [prod_val],
            "label": [label_val],
            "note": [" | ".join(note_items)],
            "orig_protein_id": [p.id],
            "translation": [str(p.seq)]
        }
        if coord_str:
            qualifiers["orig_coord"] = [coord_str]
    elif is_op:
        clean_id = p.id.split("___")[-1] if "___" in p.id else p.id
        label_val = f"Operon_{clean_id}"
        note_items = [f"Operon protein {clean_id}", f"orig_protein_id={p.id}"]
        if coord_str:
            note_items.append(f"orig_coord={coord_str}")
        qualifiers = {
            "locus_tag": [p.id],
            "protein_id": [p.id],
            "gene": [label_val],
            "orig_gene": [clean_id],
            "product": [f"Operon protein {clean_id}"],
            "label": [label_val],
            "note": [" | ".join(note_items)],
            "orig_protein_id": [p.id],
            "translation": [str(p.seq)]
        }
        if coord_str:
            qualifiers["orig_coord"] = [coord_str]
    else:
        clean_id = p.id.split("___")[-1] if "___" in p.id else p.id
        label_val = clean_id
        note_items = ["unannotated flanking protein", f"orig_protein_id={p.id}"]
        if coord_str:
            note_items.append(f"orig_coord={coord_str}")
        qualifiers = {
            "locus_tag": [p.id],
            "protein_id": [p.id],
            "gene": [label_val],
            "orig_gene": [clean_id],
            "product": ["hypothetical protein"],
            "label": [label_val],
            "note": [" | ".join(note_items)],
            "orig_protein_id": [p.id],
            "translation": [str(p.seq)]
        }
        if coord_str:
            qualifiers["orig_coord"] = [coord_str]
    return qualifiers

GENERIC_ACC_BLACKLIST = frozenset({
    "GCA", "GCF", "GB", "RS", "NZ", "CP", "ASM", "GENOMIC", "PROTEIN", "PROT", "GENE", "CONTIG", "SCAFFOLD"
})

def is_better_hit(candidate, current):
    """
    Returns True if candidate has a strictly higher bitscore than current.
    If bitscores are equal, uses evalue (lower is better) as tiebreaker.
    If evalues are equal, uses ko_rank (lower is better) as tiebreaker.
    If current is None or empty, returns True.
    """
    if not current:
        return True
    if not candidate:
        return False
    if not isinstance(current, dict):
        return True
    if not isinstance(candidate, dict):
        return False

    def _parse_score(val):
        try:
            return float(val)
        except (ValueError, TypeError):
            return -float("inf")

    def _parse_evalue(val):
        try:
            return float(val)
        except (ValueError, TypeError):
            return float("inf")

    def _parse_rank(val):
        try:
            return int(val)
        except (ValueError, TypeError):
            return 999999

    c_score = _parse_score(candidate.get("score"))
    cur_score = _parse_score(current.get("score"))
    if c_score > cur_score:
        return True
    elif c_score < cur_score:
        return False

    # Bitscores are equal, compare evalue (lower is better)
    c_eval = _parse_evalue(candidate.get("evalue"))
    cur_eval = _parse_evalue(current.get("evalue"))
    if c_eval < cur_eval:
        return True
    elif c_eval > cur_eval:
        return False

    # Evalues are equal, compare rank (lower is better)
    c_rank = _parse_rank(candidate.get("rank"))
    cur_rank = _parse_rank(current.get("rank"))
    return c_rank < cur_rank

def register_genome_alias_keys(raw_id, hit_info, g_map):
    if not raw_id:
        return
    cands = {raw_id, raw_id.replace("GB_", "").replace("RS_", "")}
    for c in list(cands):
        if "___prot" in c:
            cands.add(c.replace("___prot", "___gene"))
        elif "___gene" in c:
            cands.add(c.replace("___gene", "___prot"))
        if "prot" in c:
            cands.add(re.sub(r"prot(\d+)", r"gene\1", c))
        if "gene" in c:
            cands.add(re.sub(r"gene(\d+)", r"prot\1", c))
        if "___" in c:
            parts = c.split("___")
            sub = c.split("___", 1)[-1]
            cands.add(sub)
            cands.add(sub.replace("gene", "prot"))
            cands.add(sub.replace("prot", "gene"))
            last = c.rsplit("___", 1)[-1]
            cands.add(last)
            cands.add(last.replace("gene", "prot"))
            cands.add(last.replace("prot", "gene"))
            if len(parts) >= 3:
                ctg = parts[1]
                tail = parts[2]
                m_tail = re.search(r"(\d+)$", tail)
                if m_tail:
                    n = int(m_tail.group(1))
                    cands.add(f"{ctg}_{n}")
                    cands.add(f"{ctg}_{n:05d}")
                    cands.add(f"{ctg}___prot{n:05d}")
                    cands.add(f"{ctg}___gene{n:05d}")
                    cands.add(f"prot{n:05d}")
                    cands.add(f"gene{n:05d}")
                    cands.add(f"prot{n}")
                    cands.add(f"gene{n}")
                    cands.add(str(n))
        if "_" in c:
            last_ = c.rsplit("_", 1)[-1]
            cands.add(last_)
            m_num = re.search(r"(\d+)$", last_)
            if m_num:
                n = int(m_num.group(1))
                cands.add(f"prot{n:05d}")
                cands.add(f"gene{n:05d}")
                cands.add(f"prot{n}")
                cands.add(f"gene{n}")
                cands.add(str(n))
        m_any = re.search(r"(?:prot|gene|_)?(\d+)$", c)
        if m_any:
            n = int(m_any.group(1))
            cands.add(f"prot{n:05d}")
            cands.add(f"gene{n:05d}")
            cands.add(f"prot{n}")
            cands.add(f"gene{n}")
            cands.add(str(n))
    for k in cands:
        if k and (k not in g_map or is_better_hit(hit_info, g_map[k])):
            g_map[k] = hit_info

def load_ref_flanking_annotations_from_table(gtdb_ko_table, target_accessions, target_pids):
    """
    Fast streaming extraction of KO annotations from large GTDB KO table (e.g. protein_KO_all_strict.tsv)
    specifically for the reference genomes and target genomes/contigs in the visualized loci.
    """
    if not gtdb_ko_table or not os.path.exists(gtdb_ko_table):
        return {}

    clean_accs = set()
    for acc in target_accessions:
        if not acc:
            continue
        acc_str = str(acc).strip()
        clean = acc_str.replace("GB_", "").replace("RS_", "")
        if len(clean) >= 6 and clean.upper() not in GENERIC_ACC_BLACKLIST:
            clean_accs.add(clean)
        if "." in clean:
            prefix = clean.split(".")[0]
            if len(prefix) >= 6 and prefix.upper() not in GENERIC_ACC_BLACKLIST:
                clean_accs.add(prefix)
        # Only add accession numeric digits for genuine GCA/GCF accessions to avoid matching generic contig numbers
        if (clean.startswith("GCA_") or clean.startswith("GCF_") or clean.startswith("GB_") or clean.startswith("RS_")) and "_" in clean:
            num = clean.split("_")[-1].split(".")[0]
            if len(num) >= 7 and num.isdigit():
                clean_accs.add(num)
        if "___" in clean:
            for part in clean.split("___"):
                part_clean = part.split(".")[0]
                if len(part_clean) >= 6 and part_clean.upper() not in GENERIC_ACC_BLACKLIST:
                    clean_accs.add(part_clean)

    clean_accs_tuple = tuple(k for k in clean_accs if len(k) >= 6 and k.upper() not in GENERIC_ACC_BLACKLIST)

    # Build target protein lookup set for precise reference gene filtering
    target_pids_set = set()
    for pid in target_pids:
        if not pid:
            continue
        p_str = str(pid).strip().replace("GB_", "").replace("RS_", "")
        target_pids_set.add(pid)
        target_pids_set.add(p_str)
        if "___" in p_str:
            target_pids_set.add(p_str.split("___", 1)[-1])
            target_pids_set.add(p_str.rsplit("___", 1)[-1])
        if "_" in p_str:
            target_pids_set.add(p_str.rsplit("_", 1)[-1])

    target_accs_raw = tuple(a for a in target_accessions if a and len(a) >= 6 and a.upper() not in GENERIC_ACC_BLACKLIST)

    anno_map = {}
    by_genome = {}
    opener = gzip.open if gtdb_ko_table.endswith(".gz") else open
    try:
        with opener(gtdb_ko_table, "rt", encoding="utf-8", errors="ignore") as f:
            header_line = f.readline()
            if not header_line:
                return {}
            header = [c.strip().lower() for c in header_line.strip().split("\t")]
            prot_idx = header.index("protein_id") if "protein_id" in header else 0
            gene_idx = header.index("gene_id") if "gene_id" in header else -1
            genome_idx = header.index("genome_id") if "genome_id" in header else -1
            ko_idx = header.index("ko") if "ko" in header else -1
            def_idx = header.index("ko_definition") if "ko_definition" in header else -1
            score_idx = header.index("bitscore") if "bitscore" in header else -1
            eval_idx = header.index("evalue") if "evalue" in header else -1
            rank_idx = header.index("ko_rank") if "ko_rank" in header else -1
            thresh_idx = -1
            for cand_th in ["curated_threshold", "threshold", "score_threshold", "thrshld", "cut_off"]:
                if cand_th in header:
                    thresh_idx = header.index(cand_th)
                    break

            for line in f:
                # Fast line-level prefix filter across cols 0, 1, 2 (within first 350 chars)
                if clean_accs_tuple and not any(acc in line[:350] for acc in clean_accs_tuple):
                    continue

                parts = line.strip().split("\t")
                if len(parts) <= max(prot_idx, ko_idx if ko_idx >= 0 else 0):
                    continue

                pid = parts[prot_idx] if prot_idx < len(parts) else ""
                gid = parts[gene_idx] if (gene_idx >= 0 and gene_idx < len(parts)) else ""
                genome_id = parts[genome_idx] if (genome_idx >= 0 and genome_idx < len(parts)) else ""
                if not genome_id:
                    m_g = re.search(r"(GC[AF]_\d+\.\d+)", pid)
                    genome_id = m_g.group(1) if m_g else ""
                clean_g_acc = genome_id.replace("GB_", "").replace("RS_", "").strip() if genome_id else ""

                # If target_pids_set is provided: only keep genes that match target_pids_set,
                # OR belong to the target genome accessions
                if target_pids_set:
                    clean_pid = pid.replace("GB_", "").replace("RS_", "")
                    clean_gid = gid.replace("GB_", "").replace("RS_", "")
                    pid_suf = pid.rsplit("___", 1)[-1] if "___" in pid else ""
                    gid_suf = gid.rsplit("___", 1)[-1] if "___" in gid else ""
                    is_in_pids = (
                        pid in target_pids_set or
                        gid in target_pids_set or
                        clean_pid in target_pids_set or
                        clean_gid in target_pids_set or
                        pid_suf in target_pids_set or
                        gid_suf in target_pids_set
                    )
                    if not is_in_pids and not any(acc in line[:350] for acc in target_accs_raw):
                        continue

                ko = parts[ko_idx] if (ko_idx >= 0 and ko_idx < len(parts)) else ""
                definition = parts[def_idx] if (def_idx >= 0 and def_idx < len(parts)) else ""
                score = parts[score_idx] if (score_idx >= 0 and score_idx < len(parts)) else ""
                threshold = parts[thresh_idx] if (thresh_idx >= 0 and thresh_idx < len(parts)) else ""
                evalue = parts[eval_idx] if (eval_idx >= 0 and eval_idx < len(parts)) else ""
                rank = parts[rank_idx] if (rank_idx >= 0 and rank_idx < len(parts)) else ""

                hit_info = {
                    "ko": ko,
                    "definition": definition,
                    "score": score,
                    "threshold": threshold,
                    "evalue": evalue,
                    "rank": rank
                }

                # Register globally by full unique ID (keeping highest scoring hit)
                for raw_id in [pid, gid]:
                    if not raw_id:
                        continue
                    clean_id = raw_id.replace("GB_", "").replace("RS_", "")
                    for rid in (raw_id, clean_id):
                        if rid not in anno_map or is_better_hit(hit_info, anno_map[rid]):
                            anno_map[rid] = hit_info

                # Register shorthand aliases ONLY within this genome's dedicated sub-dictionary (keeping highest scoring hit)
                if clean_g_acc:
                    target_sub_maps = []
                    acc_match = re.search(r"(GC[AF]_\d+(?:\.\d+)?)", clean_g_acc or pid)
                    pure_acc = acc_match.group(1) if acc_match else clean_g_acc
                    g_keys = {clean_g_acc, pure_acc, pure_acc.split(".")[0]}
                    for gk in g_keys:
                        if gk:
                            if gk not in by_genome:
                                by_genome[gk] = {}
                            target_sub_maps.append(by_genome[gk])
                    for g_map in target_sub_maps:
                        register_genome_alias_keys(pid, hit_info, g_map)
                        if gid:
                            register_genome_alias_keys(gid, hit_info, g_map)

        # Attach by_genome namespace to anno_map
        anno_map["__by_genome__"] = by_genome
        for g_acc, g_dict in by_genome.items():
            if g_acc not in anno_map:
                anno_map[g_acc] = g_dict

        print(f"Loaded {len(anno_map)} protein/gene annotations across {len(by_genome)} genomes from {gtdb_ko_table}", flush=True)
    except Exception as e:
        print(f"Notice: could not parse gtdb_ko_table ({e})", file=sys.stderr)

    return anno_map

def lookup_ref_flanking_anno(pid, anno_map, genome_acc=None):
    """Looks up protein annotation in reference flanking annotation map with per-genome namespace isolation."""
    if not anno_map or not pid:
        return None

    # Determine genome-isolated sub-map
    sub_map = None
    target_acc = genome_acc
    if not target_acc:
        m_g = re.search(r"(GC[AF]_\d+\.\d+)", pid)
        if m_g:
            target_acc = m_g.group(1)

    if target_acc:
        clean_g = target_acc.replace("GB_", "").replace("RS_", "").strip()
        cand_keys = [
            clean_g,
            clean_g.split(".")[0],
            f"GB_{clean_g}",
            f"RS_{clean_g}",
            target_acc
        ]
        for ck in cand_keys:
            if ck in anno_map and isinstance(anno_map[ck], dict):
                sub_map = anno_map[ck]
                break
            if "__by_genome__" in anno_map and ck in anno_map["__by_genome__"]:
                sub_map = anno_map["__by_genome__"][ck]
                break
        if not sub_map and "__by_genome__" in anno_map:
            for bg_k, bg_v in anno_map["__by_genome__"].items():
                if clean_g in bg_k or bg_k in clean_g or clean_g.split(".")[0] in bg_k:
                    sub_map = bg_v
                    break

    # If sub_map is resolved: query strictly within this genome
    if sub_map:
        if pid in sub_map:
            return sub_map[pid]
        clean_p = pid.replace("GB_", "").replace("RS_", "")
        if clean_p in sub_map:
            return sub_map[clean_p]
        for cand in [clean_p.replace("gene", "prot"), clean_p.replace("prot", "gene")]:
            if cand in sub_map:
                return sub_map[cand]
        if "___" in pid:
            sub = pid.split("___", 1)[-1]
            for s in [sub, sub.replace("gene", "prot"), sub.replace("prot", "gene")]:
                if s in sub_map:
                    return sub_map[s]
            suf = pid.rsplit("___", 1)[-1]
            for s in [suf, suf.replace("gene", "prot"), suf.replace("prot", "gene")]:
                if s in sub_map:
                    return sub_map[s]
        if "_" in pid:
            suf = pid.rsplit("_", 1)[-1]
            for s in [suf, suf.replace("gene", "prot"), suf.replace("prot", "gene")]:
                if s in sub_map:
                    return sub_map[s]
        m = re.search(r"(?:gene|prot|_)?(\d+)$", pid)
        if m:
            n = int(m.group(1))
            for cand in [f"prot{n:05d}", f"gene{n:05d}", f"prot{n}", f"gene{n}", str(n)]:
                if cand in sub_map:
                    return sub_map[cand]
        return None

    # Fallback: exact full ID lookup only (never use stripped suffix across undefined genomes)
    if pid in anno_map and not isinstance(anno_map[pid], dict):
        return anno_map[pid]
    clean_p = pid.replace("GB_", "").replace("RS_", "")
    if clean_p in anno_map and not isinstance(anno_map[clean_p], dict):
        return anno_map[clean_p]

    return None

def find_ko_for_pid(pid, pid_to_ko_map):
    """Robustly looks up KO for a protein ID without false positive substring matches."""
    if not pid_to_ko_map:
        return ""
    if pid in pid_to_ko_map:
        return pid_to_ko_map[pid]
    clean_p = pid.replace("GB_", "").replace("RS_", "")
    for p_key, ko in pid_to_ko_map.items():
        clean_k = p_key.replace("GB_", "").replace("RS_", "")
        if clean_p == clean_k:
            return ko
        if clean_p.endswith(f"_{clean_k}") or clean_k.endswith(f"_{clean_p}"):
            return ko
    return ""

def parse_kofam_tsv(kofam_path):
    """
    Parses KofamScan / pykofamsearch tabular output.
    Header format:
    # gene name           KO     thrshld  score   E-value KO definition
      u499975_00001       K17993  555.00  535.5  7.7e-162 sulfhydrogenase subunit alpha [EC:1.12.1.3 1.12.1.5]
    """
    if not kofam_path or not os.path.exists(kofam_path):
        return {}

    gene_to_hits = {}
    with open(kofam_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            if tokens[0] == "*":
                tokens = tokens[1:]
            if len(tokens) < 5:
                continue
            gene_name = tokens[0]
            ko = tokens[1]
            threshold = tokens[2]
            try:
                score = float(tokens[3])
                evalue = float(tokens[4])
            except ValueError:
                score = 0.0
                evalue = 1.0
            definition = " ".join(tokens[5:]) if len(tokens) > 5 else ""

            gene_to_hits.setdefault(gene_name, []).append({
                "ko": ko,
                "score": score,
                "threshold": threshold,
                "evalue": evalue,
                "definition": definition
            })

    gene_to_best_anno = {}
    for g, hits in gene_to_hits.items():
        sorted_hits = sorted(hits, key=lambda x: (-x["score"], x["evalue"]))
        best = sorted_hits[0]
        gene_to_best_anno[g] = best
        if "_" in g:
            suf = g.rsplit("_", 1)[-1]
            if suf not in gene_to_best_anno:
                gene_to_best_anno[suf] = best

    return gene_to_best_anno

def lookup_target_anno(pid, anno_dict):
    """Looks up target gene annotation with cross-gene/prot translation and suffix tolerance."""
    if not anno_dict or not pid:
        return None
    if pid in anno_dict:
        return anno_dict[pid]

    # 1. Swap gene <-> prot directly
    swaps = [
        pid.replace("gene", "prot"),
        pid.replace("prot", "gene"),
        re.sub(r"gene(\d+)", r"prot\1", pid),
        re.sub(r"prot(\d+)", r"gene\1", pid),
    ]
    for s in swaps:
        if s in anno_dict:
            return anno_dict[s]

    # 2. Clean prefix
    clean_p = pid.replace("GB_", "").replace("RS_", "")
    if clean_p in anno_dict:
        return anno_dict[clean_p]
    for s in [clean_p.replace("gene", "prot"), clean_p.replace("prot", "gene")]:
        if s in anno_dict:
            return anno_dict[s]

    # 3. Splitting by ___
    if "___" in pid:
        sub = pid.split("___", 1)[-1]
        for s in [sub, sub.replace("gene", "prot"), sub.replace("prot", "gene")]:
            if s in anno_dict:
                return anno_dict[s]
        sub_last = pid.rsplit("___", 1)[-1]
        for s in [sub_last, sub_last.replace("gene", "prot"), sub_last.replace("prot", "gene")]:
            if s in anno_dict:
                return anno_dict[s]

    # 4. Splitting by _
    if "_" in pid:
        suf = pid.rsplit("_", 1)[-1]
        for s in [suf, suf.replace("gene", "prot"), suf.replace("prot", "gene")]:
            if s in anno_dict:
                return anno_dict[s]

    # 5. Regex numeric match (e.g. gene00004 -> prot00004 or 00004)
    m = re.search(r"(?:gene|prot)?(\d{4,})$", pid)
    if m:
        num = m.group(1)
        for cand in [f"gene{num}", f"prot{num}", num]:
            if cand in anno_dict:
                return anno_dict[cand]

    # 6. Suffix matching against all valid annotation entries
    for k, v in anno_dict.items():
        if k.startswith("__") or not isinstance(v, dict) or "ko" not in v:
            continue
        if pid == k or pid.endswith(f"_{k}") or k.endswith(f"_{pid}"):
            return v
        k_clean = k.replace("gene", "prot")
        p_clean = pid.replace("gene", "prot")
        if p_clean == k_clean or p_clean.endswith(f"_{k_clean}") or k_clean.endswith(f"_{p_clean}"):
            return v

    return None

def format_target_gene_qualifiers(pid, hit, trans_seq, ko_definitions, c_info=None):
    """Formats non-redundant qualifiers for target genome genes using KofamScan annotations."""
    coord_str = ""
    if c_info:
        ctg = c_info.get("contig", "")
        s = c_info.get("start")
        e = c_info.get("end")
        st = c_info.get("strand", "+")
        if s is not None and e is not None:
            coord_str = f"{ctg}:{s}..{e}({st})"

    if hit:
        ko = hit.get("ko", "")
        def_from_hit = hit.get("definition", "")
        ko_def = ko_definitions.get(ko, def_from_hit)

        sym = ""
        desc = def_from_hit
        if ko_def and ";" in ko_def:
            parts = ko_def.split(";", 1)
            sym = parts[0].split(",")[-1].strip()
            desc = parts[1].strip()
        elif def_from_hit:
            desc = def_from_hit

        clean_desc = re.sub(r"\[EC:[^\]]+\]", "", desc).strip()
        sc = hit.get("score", "")
        th = hit.get("threshold", "")
        ev = hit.get("evalue", "")
        m_str = format_metric_str(sc, th, ev)

        if sym and clean_desc and sym != clean_desc and sym != ko:
            base_label = f"{ko}: {sym}; {clean_desc}"
            gene_val = sym
        else:
            base_label = f"{ko}: {clean_desc or desc}"
            gene_val = sym or ko
        label_val = f"{base_label}{m_str}"

        note_parts = [f"KO: {ko}", f"orig_protein_id={pid}"]
        if coord_str:
            note_parts.append(f"orig_coord={coord_str}")
        if sc or th or ev:
            sc_part = f"score={sc}" if sc else ""
            th_part = f"threshold={th}" if th else ""
            ev_part = f"evalue={ev}" if ev else ""
            note_parts.append(" ".join(filter(None, [sc_part, th_part, ev_part])))
        note_val = " | ".join(note_parts)
        prod_val = desc or def_from_hit
    else:
        label_val = pid
        gene_val = pid
        prod_val = "hypothetical protein"
        note_parts = ["unannotated", f"orig_protein_id={pid}"]
        if coord_str:
            note_parts.append(f"orig_coord={coord_str}")
        note_val = " | ".join(note_parts)

    quals = {
        "locus_tag": [pid],
        "protein_id": [pid],
        "gene": [label_val],
        "orig_gene": [gene_val],
        "product": [prod_val],
        "label": [label_val],
        "note": [note_val],
        "orig_protein_id": [pid],
        "translation": [trans_seq]
    }
    if coord_str:
        quals["orig_coord"] = [coord_str]
    return quals


def parse_args():
    parser = argparse.ArgumentParser(description="Stage 1: Locate consensus anchors in target, evaluate edge safety margin, and extract locus slices.")
    parser.add_argument("--target_name", required=True, help="Target genome sample ID")
    parser.add_argument("--target_fa", required=True, help="Target genome assembly FASTA")
    parser.add_argument("--target_faa", required=True, help="Target genome predicted proteome FAA")
    parser.add_argument("--target_gff", default="", help="Target genome GFF3 annotation")
    parser.add_argument("--anchors_tsv", required=True, help="TSV of consensus anchors from Stage 0")
    parser.add_argument("--refs_tsv", required=True, help="TSV of validated references from Stage 0")
    parser.add_argument("--seed_operon_faa", default="", help="Extracted seed query operon FAA from Stage 0")
    parser.add_argument("--expected_locus_length", type=int, default=15000, help="Expected locus span in positive references (bp)")
    parser.add_argument("--min_edge_bp", type=int, default=10000, help="Minimum configured contig edge margin (bp)")
    parser.add_argument("--flank_bp", type=int, default=20000, help="Flanking sequence window size around anchors (bp)")
    parser.add_argument("--n_flank_genes", type=int, default=25, help="Number of flanking genes to extract on each side (default: 25)")
    parser.add_argument("--max_visualized_refs", type=int, default=10, help="Maximum number of confirmed syntenic references to include in Clinker plot (default: 10)")
    parser.add_argument("--consensus_contig_threshold", type=float, default=0.75, help="Minimum fraction of positive references converging on the same contig (legacy/fixed_ratio mode, default: 0.75)")
    parser.add_argument("--significance_mode", default="empirical_max", choices=["empirical_max", "fast_screen", "fixed_ratio"], help="Localization gate mode: 'empirical_max' (recommended default), 'fast_screen', or 'fixed_ratio'")
    parser.add_argument("--significance_alpha", type=float, default=0.01, help="Alpha threshold for statistical significance (default: 0.01)")
    parser.add_argument("--n_permutations", type=int, default=1000, help="Number of permutations for empirical MAX-contig null (default: 1000)")
    parser.add_argument("--allow_insufficient_evidence_locus", action="store_true", help="Extract 100kb locus even if localization has insufficient statistical evidence")
    parser.add_argument("--extended_window_bp", type=int, default=100000, help="Extended locus window size around center (bp, default: 100000)")
    parser.add_argument("--target_kos", default="", help="Comma-separated target KOs for the module")
    parser.add_argument("--diamond_env", default="", help="Conda environment for diamond")
    parser.add_argument("--target_annotation_tsv", default="", help="Optional KofamScan TSV annotation file for the target genome")
    parser.add_argument("--gtdb_ko_table", default="", help="Optional GTDB r232 protein KO annotation table for reference flanking genes")
    parser.add_argument("--ko_index_cache", default="", help="Optional pre-filtered target KO index cache TSV")
    parser.add_argument("--module_id", default="M00159", help="ATP synthase module ID (e.g. M00159, M00157, M00158, M00160)")
    parser.add_argument("--module_name", default="", help="ATP synthase module human-readable name")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPU threads for DIAMOND searches")
    parser.add_argument("--out_topology_tsv", required=True, help="Output TSV with anchor topology and edge safety metrics")
    parser.add_argument("--max_candidate_loci", type=int, default=5, help="Maximum number of candidate loci to rank and track (default: 5)")
    parser.add_argument("--out_target_locus_fa", required=True, help="Output target locus slice FASTA")
    parser.add_argument("--out_target_locus_gbk", required=True, help="Output target locus slice GBK")
    parser.add_argument("--out_refs_locus_dir", required=True, help="Directory to save reference locus slices")
    parser.add_argument("--out_extended_proteins_tsv", default="", help="Output TSV with ±100kb locus proteins and annotations across target and positive references")
    return parser.parse_args()


def run_diamond_search(query_faa, db_faa, out_tsv, diamond_env, threads=4):
    """Runs DIAMOND blastp with fallback for local environments."""
    cmd_prefix = f"conda run -p {diamond_env} " if diamond_env and os.path.isabs(diamond_env) else (f"conda run -n {diamond_env} " if diamond_env else "")
    db_prefix = out_tsv + ".dmnd"
    t_arg = f"--threads {threads}" if threads and threads > 1 else ""

    def execute_diamond(prefix):
        subprocess.run(
            f"{prefix}diamond makedb --in '{db_faa}' -d '{db_prefix}' {t_arg} --quiet",
            shell=True, check=True, stderr=subprocess.PIPE
        )
        subprocess.run(
            f"{prefix}diamond blastp -q '{query_faa}' -d '{db_prefix}' -o '{out_tsv}' -e 1e-4 --id 25.0 --query-cover 25.0 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen {t_arg} --quiet",
            shell=True, check=True, stderr=subprocess.PIPE
        )
        if os.path.exists(db_prefix + ".dmnd"):
            os.remove(db_prefix + ".dmnd")

    try:
        execute_diamond(cmd_prefix)
    except Exception as e1:
        try:
            execute_diamond("")
        except Exception as e2:
            raise RuntimeError(
                f"[extract_anchor_loci] DIAMOND blastp execution failed with conda env ({e1}) and system path ({e2}). "
                f"Fail-fast enforced; mock alignment rows will NOT be generated."
            )

def parse_gff_coordinates(gff_file):
    """Extracts protein_id -> (contig, start, end, strand) mapping from GFF3 (plain or gzipped)."""
    coords = {}
    if not gff_file or not os.path.exists(gff_file):
        return coords
    opener = gzip.open if gff_file.endswith(".gz") else open
    try:
        with opener(gff_file, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 9 and parts[2] in ["CDS", "gene", "polypeptide"]:
                    contig = parts[0]
                    full_ctg = contig
                    if "___" in contig:
                        c_short = contig.split("___")[1]
                    else:
                        c_short = contig
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attrs = parts[8]
                    info = {"contig": c_short, "full_contig": full_ctg, "start": start, "end": end, "strand": strand}

                    m_pid = re.search(r"protein_id=([^;]+)", attrs)
                    if m_pid:
                        pid = m_pid.group(1).strip()
                        coords[pid] = info
                        coords[pid.split()[0]] = info
                        if "___" in pid:
                            coords[pid.split("___")[-1]] = info
                            coords[pid.split("___", 1)[-1]] = info

                    m_id = re.search(r"ID=([^;]+)", attrs)
                    if m_id:
                        gid = m_id.group(1).strip()
                        coords[gid] = info
                        coords[gid.split()[0]] = info
                        if "___" in gid:
                            coords[gid.split("___")[-1]] = info
                            coords[gid.split("___", 1)[-1]] = info

                    m_name = re.search(r"Name=([^;]+)", attrs)
                    if m_name:
                        nid = m_name.group(1).strip()
                        coords[nid] = info
                        coords[nid.split()[0]] = info
    except Exception as e:
        print(f"Notice: GFF parsing encountered error ({e}), continuing with fallbacks.", file=sys.stderr)
    return coords

def discover_companion_gff(faa_file, target_name):
    """Auto-discovers companion genes/*.gff.gz or metadata/*.tsv.gz relative to the proteome file."""
    if not faa_file:
        return ""
    faa_dir = os.path.dirname(os.path.abspath(faa_file))
    clean_name = target_name.replace("GB_", "").replace("RS_", "")
    cand_dirs = [
        faa_dir,
        os.path.join(os.path.dirname(faa_dir), "genes"),
        os.path.join(os.path.dirname(os.path.dirname(faa_dir)), "genes"),
        os.path.join(os.path.dirname(faa_dir), "metadata"),
        os.path.join(os.path.dirname(os.path.dirname(faa_dir)), "metadata"),
    ]
    for d in cand_dirs:
        if os.path.isdir(d):
            matches = glob.glob(os.path.join(d, f"*{clean_name}*.gff*")) or glob.glob(os.path.join(d, f"*{target_name}*.gff*"))
            if matches:
                return matches[0]
    return ""

def parse_fasta_header_coordinates(faa_file):
    """
    Extracts coordinates and contig mapping from FAA headers or protein IDs.
    Supports:
      1) Prodigal standard format: # start # end # strand # ...
      2) Pyrodigal bracket format: [start - end] [strand=...]
      3) Key-value format: start=... end=... strand=...
      Records without observed coordinates are omitted.
    """
    coords = {}
    if not faa_file or not os.path.exists(faa_file):
        return coords
    for r in load_fasta_records(faa_file):
        desc = r.description
        rid = r.id

        # Determine contig
        if "___" in rid:
            subparts = rid.split("___")
            contig = subparts[1]
            full_ctg = rid.rsplit("___", 1)[0]
        elif "_" in rid:
            contig = rid.rsplit("_", 1)[0]
            full_ctg = contig
        else:
            contig = rid
            full_ctg = rid

        start, end, strand = None, None, "+"

        # 1. Prodigal # format
        parts = desc.split("#")
        if len(parts) >= 4:
            try:
                start = int(parts[1].strip())
                end = int(parts[2].strip())
                strand = "+" if int(parts[3].strip()) >= 0 else "-"
            except Exception:
                pass

        # 2. Pyrodigal [start - end] [strand=+] format
        if start is None:
            m = re.search(r"\[\s*(\d+)\s*-\s*(\d+)\s*\]", desc)
            if m:
                start = int(m.group(1))
                end = int(m.group(2))
                sm = re.search(r"\[\s*strand\s*=\s*([+-])\s*\]", desc)
                strand = sm.group(1) if sm else "+"

        # 3. start=... end=...
        if start is None:
            m = re.search(r"start=(\d+)\s+end=(\d+)", desc)
            if m:
                start = int(m.group(1))
                end = int(m.group(2))
                sm = re.search(r"strand=([+-])", desc)
                strand = sm.group(1) if sm else "+"

        # A protein identifier's ordinal is not a genomic coordinate.
        if start is None:
            continue

        info = {"contig": contig, "full_contig": full_ctg, "start": start, "end": end, "strand": strand}
        coords[rid] = info
        r_clean = rid.split()[0]
        coords[r_clean] = info
        if "___" in rid:
            coords[rid.split("___")[-1]] = info
            coords[rid.split("___", 1)[-1]] = info
    return coords


def group_reference_proteins_by_coordinate(ref_records, ref_coords):
    """Group reference CDSs using observed contig and start fields."""
    grouped = {}
    for record in ref_records:
        info = ref_coords.get(record.id)
        contig = info.get("full_contig") or info.get("contig") if isinstance(info, dict) else None
        if not isinstance(info, dict) or not contig or info.get("start") is None:
            continue
        grouped.setdefault(contig, []).append(record)
    return {contig: list(enumerate(sorted(records, key=lambda record: int(ref_coords[record.id]["start"]))))
            for contig, records in grouped.items()}


def is_same_contig(c1, c2):
    """Determines whether two contig identifiers refer to the same biological contig."""
    if not c1 or not c2:
        return False
    if c1 == c2:
        return True
    c1_sub = c1.split("___")[-1].strip("_") if "___" in c1 else c1.strip("_")
    c2_sub = c2.split("___")[-1].strip("_") if "___" in c2 else c2.strip("_")
    c1_clean = c1_sub.split()[0]
    c2_clean = c2_sub.split()[0]
    if c1_clean == c2_clean:
        return True
    return (c1_sub == c2_sub or c1.endswith(c2) or c2.endswith(c1) or
            c1_sub in c2_sub or c2_sub in c1_sub)


def find_matching_record(ctg_name, records_dict):
    """Finds matching SeqRecord or entry in records_dict using robust contig ID matching."""
    if not ctg_name or not records_dict:
        return None
    if ctg_name in records_dict:
        return records_dict[ctg_name]
    c_clean = ctg_name.split("___")[-1].strip("_").split()[0] if "___" in ctg_name else ctg_name.strip("_").split()[0]
    for k, rec in records_dict.items():
        k_id = getattr(rec, "id", str(k))
        k_clean = k_id.split("___")[-1].strip("_").split()[0] if "___" in k_id else k_id.strip("_").split()[0]
        desc = getattr(rec, "description", "")
        if (k == ctg_name or k == c_clean or k_clean == c_clean or
            k.startswith(c_clean + " ") or ctg_name.startswith(k_clean + " ") or
            c_clean in desc or is_same_contig(k_id, ctg_name)):
            return rec
    return None


def dereplicate_reference_cohort(pos_refs):
    """
    Phylogenetic dereplication of positive references:
    - Merges identical assembly accessions (e.g. GCF_009778815.1 and GCA_009783705.1).
    - Clusters by GTDB species (s__...) or ANI >= 99.0% into K independent evolutionary clusters.
    - Assigns cluster weight w_r = 1 / |C_k| and computes N_eff = K.
    Returns: (clusters, ref_to_cluster, weights, n_eff)
    """
    clusters = {}
    ref_to_cluster = {}
    weights = {}

    if pos_refs is None or (isinstance(pos_refs, pd.DataFrame) and pos_refs.empty):
        return clusters, ref_to_cluster, weights, 0

    records = []
    if isinstance(pos_refs, pd.DataFrame):
        for _, row in pos_refs.iterrows():
            records.append({
                "accession": str(row.get("accession", "")).strip(),
                "taxonomy": str(row.get("taxonomy", "")).strip(),
                "ani": float(row.get("ANI", 0.0)) if pd.notnull(row.get("ANI")) else 0.0
            })
    elif isinstance(pos_refs, list):
        for item in pos_refs:
            if isinstance(item, dict):
                records.append(item)
            else:
                records.append({"accession": str(item), "taxonomy": "", "ani": 0.0})

    # Step 1: Group by assembly accession number (GCF/GCA duplicate assemblies)
    assembly_groups = {}
    for rec in records:
        acc = rec["accession"]
        if not acc or acc == "NONE":
            continue
        clean_acc = acc.replace("GB_", "").replace("RS_", "").strip()
        m = re.search(r"(?:GC[AF]_)?(\d+\.\d+)", clean_acc)
        acc_num = m.group(1) if m else clean_acc
        assembly_groups.setdefault(acc_num, []).append(rec)

    # Step 2: Group by GTDB species or ANI >= 99.0
    cluster_idx = 1
    processed_acc_nums = set()

    for acc_num, group in assembly_groups.items():
        if acc_num in processed_acc_nums:
            continue
        c_id = f"Cluster_{cluster_idx:02d}"
        cluster_recs = list(group)
        processed_acc_nums.add(acc_num)

        sp_set = set()
        for r in group:
            m_sp = re.search(r"s__([^;]+)", r.get("taxonomy", ""))
            if m_sp and m_sp.group(1).strip() and m_sp.group(1).strip() != "":
                sp_set.add(m_sp.group(1).strip())

        for other_num, other_group in assembly_groups.items():
            if other_num in processed_acc_nums:
                continue
            other_sp_set = set()
            max_ani = 0.0
            for r in other_group:
                m_sp = re.search(r"s__([^;]+)", r.get("taxonomy", ""))
                if m_sp and m_sp.group(1).strip() and m_sp.group(1).strip() != "":
                    other_sp_set.add(m_sp.group(1).strip())
                if r.get("ani", 0.0) > max_ani:
                    max_ani = r.get("ani", 0.0)

            if (sp_set and other_sp_set and sp_set.intersection(other_sp_set)) or max_ani >= 99.0:
                cluster_recs.extend(other_group)
                processed_acc_nums.add(other_num)

        c_accs = [r["accession"].replace("GB_", "").replace("RS_", "").strip() for r in cluster_recs]
        clusters[c_id] = c_accs
        for a in c_accs:
            ref_to_cluster[a] = c_id
            weights[a] = 1.0 / len(c_accs)
        cluster_idx += 1

    return clusters, ref_to_cluster, weights, len(clusters)


def compute_empirical_max_contig_null(
    target_records,
    gff_coords,
    ctg_scores,
    winning_ctg_id,
    cluster_scores,
    n_permutations=1000,
    alpha=0.01,
    random_seed=42
):
    """
    Computes the Empirical MAX-Contig Permutation Null (S_max) distribution:
    - S_obs: observed score on winning contig.
    - S_max^(b): maximum score achieved across all target contigs under random pseudo-locus sampling.
    - P_emp = (1 + sum(S_max^(b) >= S_obs)) / (B + 1).
    - Dominance ratio: D = S_obs / (S_runner_up + 1e-5).
    """
    import numpy as np
    np.random.seed(random_seed)

    s_obs = float(ctg_scores.get(winning_ctg_id, 0.0))
    runner_ups = [v for k, v in ctg_scores.items() if not is_same_contig(k, winning_ctg_id)]
    s_runner_up = float(max(runner_ups)) if runner_ups else 0.0
    dominance_ratio = float(s_obs) / float(s_runner_up + 1e-5) if s_runner_up > 0 else 999.0

    ctg_cds_counts = {}
    for pid, inf in gff_coords.items():
        c = inf.get("contig", "")
        if c in target_records:
            ctg_cds_counts[c] = ctg_cds_counts.get(c, 0) + 1
    for c in target_records:
        if c not in ctg_cds_counts:
            ctg_cds_counts[c] = 1

    total_cds = sum(ctg_cds_counts.values())
    ctg_names = list(ctg_cds_counts.keys())
    ctg_probs = np.array([ctg_cds_counts[c] / total_cds for c in ctg_names])
    cluster_weights = list(cluster_scores.values()) if cluster_scores else [s_obs]
    n_weights = len(cluster_weights)
    if n_weights == 0 or len(ctg_names) == 0:
        null_s_max = np.zeros(n_permutations)
    elif len(ctg_names) == 1:
        null_s_max = np.full(n_permutations, sum(cluster_weights))
    else:
        # Vectorized sampling: generate all (n_permutations, n_weights) choices in C level
        choice_indices = np.random.choice(len(ctg_names), size=(n_permutations, n_weights), p=ctg_probs)
        weights_arr = np.array(cluster_weights, dtype=float)
        null_s_max = np.zeros(n_permutations, dtype=float)
        num_ctgs = len(ctg_names)
        for b in range(n_permutations):
            null_s_max[b] = np.max(np.bincount(choice_indices[b], weights=weights_arr, minlength=num_ctgs))

    p_emp = float(1.0 + np.sum(null_s_max >= s_obs)) / float(n_permutations + 1.0)
    s_null_mean = float(np.mean(null_s_max))
    s_null_95 = float(np.percentile(null_s_max, 95))
    stat_support = "SIGNIFICANT" if p_emp < alpha else "NOT_SIGNIFICANT"

    return p_emp, s_obs, s_null_mean, s_null_95, dominance_ratio, stat_support


def extract_extended_100kb_locus_proteins(
    args,
    pos_refs,
    ref_candidates,
    ref_evaluations,
    target_records,
    gff_coords,
    target_faa_dict,
    target_anno_dict,
    ref_flank_annotations,
    ko_definitions,
    out_tsv_path,
    target_contig_id=None,
    winning_geometry=None
):
    """
    Extracts all proteins in a ±100kb window (up to contig boundaries) centered on the
    weighted anchor midpoint of the target locus contig.
    Outputs a standardized TSV table with topology labels, reference support counts, and KO annotations.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv_path)), exist_ok=True)
    header = [
        "target_id", "module_id", "genome_id", "genome_type", "contig_id",
        "protein_id", "gene_symbol", "start", "end", "strand", "length_aa",
        "ko", "ko_definition", "label", "mapped_ref_count", "mapped_ref_ratio"
    ]

    def write_empty(reason):
        print(f"[{args.target_name}] {reason}", flush=True)
        with open(out_tsv_path, "w", encoding="utf-8") as f:
            f.write("\t".join(header) + "\n")

    n_pos = len(pos_refs) if pos_refs is not None else 0
    if n_pos == 0:
        n_pos = len(ref_candidates)
    if n_pos == 0:
        write_empty("No confirmed positive references available. Emitting header-only 100kb locus TSV.")
        return

    def get_canonical_target_contig(ctg_name):
        if not ctg_name or not target_records:
            return None
        if ctg_name in target_records:
            return target_records[ctg_name]
        c_strip = ctg_name.split("___")[-1].strip("_") if "___" in ctg_name else ctg_name.strip("_")
        for k, rec in target_records.items():
            k_strip = k.split("___")[-1].strip("_") if "___" in k else k.strip("_")
            if k == ctg_name or k.endswith(ctg_name) or ctg_name.endswith(k) or k_strip == c_strip or k_strip in c_strip or c_strip in k_strip:
                return rec
        return None

    # 1. Evaluate Contig Consensus across positive references
    ctg_ref_map = {}
    for r_acc, cand in ref_candidates.items():
        ev = ref_evaluations.get(r_acc, {})
        r_best_c = ev.get("best_contig", "")
        r_score = ev.get("score", 0.0)
        r_cluster = ev.get("best_cluster", [])
        if r_best_c and r_score > 0 and r_cluster:
            canon_rec = get_canonical_target_contig(r_best_c)
            canon_id = canon_rec.id if canon_rec else r_best_c
            ctg_ref_map.setdefault(canon_id, []).append(r_acc)

    if not ctg_ref_map and not target_contig_id:
        write_empty("No positive references localized with positive synteny scores. Emitting header-only 100kb locus TSV.")
        return

    # Identify consensus target contig
    mapped_refs = set(r for r_list in ctg_ref_map.values() for r in r_list)
    n_mapped = len(mapped_refs)
    if target_contig_id:
        canon_t_rec = get_canonical_target_contig(target_contig_id)
        winning_ctg_id = canon_t_rec.id if canon_t_rec else target_contig_id
        converged_refs = ctg_ref_map.get(winning_ctg_id, [])
        if not converged_refs:
            converged_refs = [
                r for r, ev in ref_evaluations.items()
                if (ev.get("best_contig") and is_same_contig(ev.get("best_contig"), winning_ctg_id)) or
                   any(is_same_contig(h.get("contig"), winning_ctg_id) for h in ev.get("best_cluster", []))
            ]
        if not converged_refs and ctg_ref_map:
            # Fallback to top contig from ctg_ref_map if chosen contig had no direct hits
            winning_ctg_id = max(ctg_ref_map.keys(), key=lambda c: (len(ctg_ref_map[c]), sum(ref_evaluations[r]["score"] for r in ctg_ref_map[c])))
            converged_refs = ctg_ref_map[winning_ctg_id]
    else:
        winning_ctg_id = max(ctg_ref_map.keys(), key=lambda c: (len(ctg_ref_map[c]), sum(ref_evaluations[r]["score"] for r in ctg_ref_map[c])))
        converged_refs = ctg_ref_map[winning_ctg_id]

    consensus_ratio = float(len(converged_refs)) / float(n_mapped) if n_mapped > 0 else 0.0

    sig_mode = getattr(args, "significance_mode", "empirical_max")
    threshold = float(getattr(args, "consensus_contig_threshold", 0.75))

    if sig_mode == "fixed_ratio":
        if consensus_ratio < threshold:
            write_empty(
                f"Consensus contig ratio ({len(converged_refs)}/{n_mapped} mapped refs = {consensus_ratio:.1%}, out of {n_pos} total pos refs) "
                f"on contig '{winning_ctg_id}' is below threshold ({threshold:.1%}). Emitting header-only 100kb locus TSV."
            )
            return

    geom_label = f" ({winning_geometry})" if winning_geometry else ""
    print(
        f"[{args.target_name}] Candidate locus contig identified: '{winning_ctg_id}'{geom_label} "
        f"supported by {len(converged_refs)}/{n_mapped} mapped positive references ({consensus_ratio:.1%}, total pos refs={n_pos}). Extracting ±100kb locus proteins.",
        flush=True
    )

    winning_rec = get_canonical_target_contig(winning_ctg_id)
    target_ctg_len = len(winning_rec.seq) if winning_rec else 0

    # 2. Calculate Weighted Midpoint across converged references
    # Weight = span length of anchor hits on winning contig
    weights = []
    midpoints = []
    ref_anchor_clusters = {}
    ref_intervals = {}

    for r_acc in converged_refs:
        ev = ref_evaluations[r_acc]
        cl = ev.get("best_cluster", [])
        ref_anchor_clusters[r_acc] = cl
        if cl:
            r_start = min(h["start"] for h in cl)
            r_end = max(h["end"] for h in cl)
            span = max(1, r_end - r_start)
            mid = (r_start + r_end) / 2.0
            weights.append(span)
            midpoints.append(mid)

            # Compute UP-DOWN breakpoint interval for this reference
            u_hits = [h for h in cl if h["side"] == "UP"]
            d_hits = [h for h in cl if h["side"] == "DOWN"]
            if u_hits and d_hits:
                med_u = sum((h["start"] + h["end"]) / 2.0 for h in u_hits) / len(u_hits)
                med_d = sum((h["start"] + h["end"]) / 2.0 for h in d_hits) / len(d_hits)
                if med_u <= med_d:
                    u_bound = max(h["end"] for h in u_hits)
                    d_bound = min(h["start"] for h in d_hits)
                else:
                    u_bound = max(h["end"] for h in d_hits)
                    d_bound = min(h["start"] for h in u_hits)
                int_s = min(u_bound, d_bound)
                int_e = max(u_bound, d_bound)
                ref_intervals[r_acc] = (int_s, int_e)

    if weights:
        weighted_center = sum(w * m for w, m in zip(weights, midpoints)) / sum(weights)
    else:
        weighted_center = target_ctg_len / 2.0 if target_ctg_len > 0 else 0.0

    center_int = int(round(weighted_center))
    win_bp = int(getattr(args, "extended_window_bp", 100000))
    t_win_start = max(0, center_int - win_bp)
    t_win_end = min(target_ctg_len, center_int + win_bp) if target_ctg_len > 0 else center_int + win_bp

    print(
        f"[{args.target_name}] Calculated weighted anchor center on {winning_ctg_id}: {center_int} bp "
        f"(window: {t_win_start}-{t_win_end} bp, length={t_win_end - t_win_start} bp, total_weight={sum(weights):.0f}).",
        flush=True
    )

    # Resolve module target KOs
    mod_target_kos = set()
    if getattr(args, "target_kos", ""):
        mod_target_kos = {k.strip() for k in args.target_kos.split(",") if k.strip()}
    if not mod_target_kos and args.module_id in MODULE_TARGET_KOS:
        mod_target_kos = MODULE_TARGET_KOS[args.module_id]

    # 3. Extract and Classify Target Proteins in Window
    target_prots_in_win = []
    seen_target_pids = set()

    for pid, inf in gff_coords.items():
        c_name = inf.get("contig", "")
        c_full = inf.get("full_contig", "")
        if not (is_same_contig(c_name, winning_ctg_id) or is_same_contig(c_full, winning_ctg_id)):
            continue
        p_s = inf["start"]
        p_e = inf["end"]
        if p_e > t_win_start and p_s < t_win_end:
            clean_p = pid.split("___")[-1] if "___" in pid else pid
            if clean_p not in seen_target_pids:
                seen_target_pids.add(clean_p)
                seen_target_pids.add(pid)
                target_prots_in_win.append((pid, inf))

    target_prots_in_win.sort(key=lambda x: x[1]["start"])

    # Prepare mapped target pids per reference
    ref_mapped_target_pids = {}
    for r_acc in converged_refs:
        m_set = set()
        for h in ref_anchor_clusters.get(r_acc, []):
            tp = h.get("t_pid", "")
            m_set.add(tp)
            m_set.add(tp.split("___")[-1])
            m_set.add(tp.replace("gene", "prot"))
            m_set.add(tp.replace("prot", "gene"))
        ref_mapped_target_pids[r_acc] = m_set

    rows = []

    for pid, inf in target_prots_in_win:
        p_s = inf["start"]
        p_e = inf["end"]
        st = inf.get("strand", "+")
        clean_p = pid.split("___")[-1] if "___" in pid else pid

        # Evaluate votes across converged references
        k_mapped = 0
        k_interval = 0
        k_flank = 0

        for r_acc in converged_refs:
            m_set = ref_mapped_target_pids[r_acc]
            if pid in m_set or clean_p in m_set:
                k_mapped += 1
            elif r_acc in ref_intervals:
                int_s, int_e = ref_intervals[r_acc]
                if p_e > int_s and p_s < int_e:
                    k_interval += 1
                else:
                    k_flank += 1
            else:
                k_flank += 1

        if k_mapped > 0:
            label = "anchor_mapped"
            mapped_cnt = k_mapped
        elif k_interval > 0:
            label = "interval"
            mapped_cnt = k_interval
        else:
            label = "flanking"
            mapped_cnt = k_flank

        ratio_denom = len(converged_refs) if len(converged_refs) > 0 else (n_mapped if n_mapped > 0 else n_pos)
        ratio = round(float(mapped_cnt) / float(ratio_denom), 4)

        # Look up KO & Definition
        t_hit = lookup_target_anno(pid, target_anno_dict)
        if not t_hit and clean_p != pid:
            t_hit = lookup_target_anno(clean_p, target_anno_dict)
        if not t_hit:
            for r_acc in converged_refs:
                for h in ref_anchor_clusters.get(r_acc, []):
                    if h.get("t_pid") in (pid, clean_p):
                        orig_p = h.get("orig_pid", "")
                        ref_hit = lookup_ref_flanking_anno(orig_p, ref_flank_annotations, genome_acc=r_acc)
                        if ref_hit and ref_hit.get("ko"):
                            t_hit = ref_hit
                            break
                if t_hit:
                    break

        ko = t_hit.get("ko", "") if t_hit else ""
        ko_def = t_hit.get("definition", "") if t_hit else ""
        if ko and not ko_def and ko in ko_definitions:
            ko_def = ko_definitions[ko]
        gene_sym = ko_def.split(";")[0].strip() if (ko_def and ";" in ko_def) else ""

        trans_seq = target_faa_dict.get(pid, "") or target_faa_dict.get(clean_p, "")
        if not trans_seq:
            for alt_id in [pid.replace("gene", "prot"), pid.replace("prot", "gene"), f"{winning_ctg_id}___{pid}"]:
                if alt_id in target_faa_dict:
                    trans_seq = target_faa_dict[alt_id]
                    break
        len_aa = len(trans_seq) if trans_seq else max(1, (p_e - p_s + 1) // 3)

        rows.append({
            "target_id": args.target_name,
            "module_id": args.module_id,
            "genome_id": args.target_name,
            "genome_type": "target",
            "contig_id": winning_ctg_id,
            "protein_id": clean_p,
            "gene_symbol": gene_sym,
            "start": p_s,
            "end": p_e,
            "strand": st,
            "length_aa": len_aa,
            "ko": ko,
            "ko_definition": ko_def,
            "label": label,
            "mapped_ref_count": mapped_cnt,
            "mapped_ref_ratio": ratio
        })

    # 4. Extract and Classify Reference Proteins in Window for each Converged Reference
    for r_acc in converged_refs:
        cand = ref_candidates.get(r_acc)
        if not cand:
            continue
        r_ctg = cand.get("best_ctg", "")
        op_pids = cand.get("op_pids", set())
        pid_to_ko = cand.get("pid_to_ko", {})
        r_coords = cand.get("ref_coords", {})
        r_faa_dict = cand.get("ref_faa_dict", {})

        # Compute operon midpoint on reference contig
        op_cds = [r_coords[p] for p in op_pids if p in r_coords]
        if op_cds:
            r_op_s = min(c["start"] for c in op_cds)
            r_op_e = max(c["end"] for c in op_cds)
            r_op_mid = (r_op_s + r_op_e) / 2.0
        else:
            genes_in_ord = cand.get("u_prots", []) + cand.get("core_prots", []) + cand.get("d_prots", [])
            coords_found = [r_coords[p.id] for p in genes_in_ord if p.id in r_coords]
            r_op_mid = (min(c["start"] for c in coords_found) + max(c["end"] for c in coords_found)) / 2.0 if coords_found else 0.0

        r_win_s = max(0, int(r_op_mid - win_bp))
        r_win_e = int(r_op_mid + win_bp)

        # Set of reference proteins that successfully mapped to target
        ref_mapped_pids = {h["orig_pid"] for h in ref_anchor_clusters.get(r_acc, [])}

        ref_prots_in_win = []
        seen_ref_pids = set()
        for p, inf in r_coords.items():
            if inf.get("contig") != r_ctg and inf.get("full_contig") != r_ctg:
                continue
            p_s = inf["start"]
            p_e = inf["end"]
            if p_e > r_win_s and p_s < r_win_e:
                clean_ref_p = p.split("___")[-1] if "___" in p else p
                if clean_ref_p not in seen_ref_pids:
                    seen_ref_pids.add(clean_ref_p)
                    seen_ref_pids.add(p)
                    ref_prots_in_win.append((p, inf))

        ref_prots_in_win.sort(key=lambda x: x[1]["start"])

        for pid, inf in ref_prots_in_win:
            p_s = inf["start"]
            p_e = inf["end"]
            st = inf.get("strand", "+")
            clean_p = pid.split("___")[-1] if "___" in pid else pid

            # Classify label
            is_op = (pid in op_pids or clean_p in op_pids or
                     pid_to_ko.get(pid) in mod_target_kos or pid_to_ko.get(clean_p) in mod_target_kos)
            is_mapped = (pid in ref_mapped_pids or clean_p in ref_mapped_pids)

            if is_op:
                label = "ref_target"
            elif is_mapped:
                label = "ref_anchor_mapped"
            else:
                label = "ref_flanking"

            # Look up KO & Definition
            ko = pid_to_ko.get(pid) or pid_to_ko.get(clean_p) or ""
            ko_def = ""
            if not ko:
                ref_hit = lookup_ref_flanking_anno(pid, ref_flank_annotations, genome_acc=r_acc)
                if not ref_hit and clean_p != pid:
                    ref_hit = lookup_ref_flanking_anno(clean_p, ref_flank_annotations, genome_acc=r_acc)
                if ref_hit:
                    ko = ref_hit.get("ko", "")
                    ko_def = ref_hit.get("definition", "")

            if ko and not ko_def and ko in ko_definitions:
                ko_def = ko_definitions[ko]
            gene_sym = ko_def.split(";")[0].strip() if (ko_def and ";" in ko_def) else ""

            seq_str = r_faa_dict.get(pid, "") or r_faa_dict.get(clean_p, "")
            len_aa = len(seq_str) if seq_str else max(1, (p_e - p_s + 1) // 3)

            rows.append({
                "target_id": args.target_name,
                "module_id": args.module_id,
                "genome_id": cand.get("clean_acc", r_acc),
                "genome_type": "positive_ref",
                "contig_id": r_ctg,
                "protein_id": clean_p,
                "gene_symbol": gene_sym,
                "start": p_s,
                "end": p_e,
                "strand": st,
                "length_aa": len_aa,
                "ko": ko,
                "ko_definition": ko_def,
                "label": label,
                "mapped_ref_count": 1,
                "mapped_ref_ratio": 1.0
            })

    df_out = pd.DataFrame(rows, columns=header)
    df_out.to_csv(out_tsv_path, sep="\t", index=False)
    print(
        f"[{args.target_name}] Successfully wrote {len(rows)} proteins to ±100kb locus TSV: {out_tsv_path} "
        f"({len(target_prots_in_win)} target, {len(rows) - len(target_prots_in_win)} ref across {len(converged_refs)} refs).",
        flush=True
    )

def main():
    args = parse_args()
    os.makedirs(args.out_refs_locus_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_topology_tsv), exist_ok=True)

    # Clean out any stale outputs from previous runs in out_refs_locus_dir
    if os.path.exists(args.out_refs_locus_dir):
        for stale_f in os.listdir(args.out_refs_locus_dir):
            if stale_f.endswith((".gbk", ".fa", ".faa", ".txt", ".tsv")):
                try:
                    os.remove(os.path.join(args.out_refs_locus_dir, stale_f))
                except Exception:
                    pass

    ko_definitions = load_kegg_ko_definitions()
    target_anno_dict = parse_kofam_tsv(args.target_annotation_tsv) if args.target_annotation_tsv else {}
    if not target_anno_dict:
        target_clean_acc = args.target_name.replace("GB_", "").replace("RS_", "")
        t_accs = [args.target_name, target_clean_acc]

        cache_file = args.ko_index_cache
        if not cache_file and args.out_topology_tsv:
            potential_cache = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.out_topology_tsv)), "..", "..", "cache", "target_kos_gtdb_index.tsv"))
            if os.path.exists(potential_cache) and os.path.getsize(potential_cache) > 0:
                cache_file = potential_cache

        if cache_file and os.path.exists(cache_file):
            target_anno_dict = load_ref_flanking_annotations_from_table(cache_file, t_accs, set())
            if target_anno_dict:
                print(f"[{args.target_name}] Loaded {len(target_anno_dict)} gene annotations for target genome from cache: {cache_file}", flush=True)

        if not target_anno_dict and args.gtdb_ko_table and os.path.exists(args.gtdb_ko_table):
            target_anno_dict = load_ref_flanking_annotations_from_table(args.gtdb_ko_table, t_accs, set())
            if target_anno_dict:
                print(f"[{args.target_name}] Loaded {len(target_anno_dict)} gene annotations for target genome from GTDB KO table: {args.gtdb_ko_table}", flush=True)

        if target_anno_dict:
            t_sub = target_anno_dict.get(target_clean_acc) or target_anno_dict.get("__by_genome__", {}).get(target_clean_acc)
            if t_sub and isinstance(t_sub, dict):
                for k, v in t_sub.items():
                    if k not in target_anno_dict or is_better_hit(v, target_anno_dict[k]):
                        target_anno_dict[k] = v
    elif target_anno_dict:
        print(f"[{args.target_name}] Loaded {len(target_anno_dict)} gene annotations from target annotation TSV: {args.target_annotation_tsv}", flush=True)

    # 1. Load consensus anchors
    anchors_df = pd.read_csv(args.anchors_tsv, sep="\t")
    anchor_dict = {r["anchor_id"]: r for _, r in anchors_df.iterrows()}

    refs_df = pd.read_csv(args.refs_tsv, sep="\t")
    pos_refs = refs_df[refs_df["validation_status"] == "CONFIRMED_POSITIVE"]
    if pos_refs.empty:
        print(f"[{args.target_name}] No validated positive reference locus; reporting unresolved localization.", flush=True)
        pd.DataFrame([{
            "target_id": args.target_name, "module_id": args.module_id,
            "module_name": args.module_name or args.module_id,
            "nearest_positive_ref": "NONE", "same_contig": "NO",
            "syntenic_geometry": "UNRESOLVED", "localization_status": "NO_VALIDATED_POSITIVE_REF",
            "gate1_status": "INSUFFICIENT_EVIDENCE", "locator_score": 0.0,
            "bilateral_context_score": 0.0, "target_contig": "",
            "target_slice_start": 0, "target_slice_end": 0,
            "ancestral_interval_start": 0, "ancestral_interval_end": 0,
            "observed_anchor_gap": -1, "observed_interanchor_span": -1,
            "expected_interanchor_span": args.expected_locus_length,
            "expected_locus_length": args.expected_locus_length,
            "intervening_genes": "", "intervening_genes_count": 0,
        }]).to_csv(args.out_topology_tsv, sep="\t", index=False)
        for output_path in (args.out_target_locus_fa, args.out_target_locus_gbk):
            os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
            open(output_path, "w").close()
        os.makedirs(os.path.dirname(os.path.abspath(args.out_extended_proteins_tsv)), exist_ok=True)
        pd.DataFrame(columns=["genome_id", "protein_id", "start", "end", "strand"]).to_csv(
            args.out_extended_proteins_tsv, sep="\t", index=False)
        return

    # 2. Extract target coordinates from GFF and FAA
    if not args.target_faa or not os.path.exists(args.target_faa):
        raise FileNotFoundError(f"[{args.target_name}] Target proteome file '{args.target_faa}' does not exist! Please check configuration paths.")

    target_faa_dict = {r.id: str(r.seq) for r in load_fasta_records(args.target_faa)}

    target_gff = args.target_gff
    if not target_gff or not os.path.exists(target_gff):
        target_gff = discover_companion_gff(args.target_faa, args.target_name)

    gff_coords = parse_gff_coordinates(target_gff) if target_gff else {}
    faa_coords = parse_fasta_header_coordinates(args.target_faa)
    for k, v in faa_coords.items():
        if k not in gff_coords:
            gff_coords[k] = v

    if not gff_coords:
        raise ValueError(f"[{args.target_name}] Could not extract gene coordinates or contig associations from '{args.target_faa}' or GFF '{target_gff}'!")

    print(f"[{args.target_name}] Successfully indexed {len(gff_coords)} CDS coordinates/contigs for target proteome.", flush=True)

    target_records = {r.id: r for r in load_fasta_records(args.target_fa)}
    if not target_records:
        raise ValueError(f"Target genome FASTA {args.target_fa} is empty!")

    # 3. Process Positive References and Extract Flanking Genes
    ref_clusters, ref_to_cluster, ref_weights, n_eff = dereplicate_reference_cohort(pos_refs)
    print(f"[{args.target_name}] Reference cohort dereplication: {len(pos_refs)} positive references grouped into {n_eff} independent evolutionary clusters: {ref_clusters}", flush=True)

    operon_records = load_fasta_records(args.seed_operon_faa) if args.seed_operon_faa and os.path.exists(args.seed_operon_faa) else []

    # 3.1 Direct search of seed operon against target proteome to locate target operon genes (robust to missing kofam.tsv)
    target_op_pids = set()
    target_pid_to_ko = {}
    if args.seed_operon_faa and os.path.exists(args.seed_operon_faa) and args.target_faa and os.path.exists(args.target_faa):
        target_op_blast = os.path.join(args.out_refs_locus_dir, f"{args.target_name}_vs_seed_operon.tsv")
        try:
            run_diamond_search(args.seed_operon_faa, args.target_faa, target_op_blast, args.diamond_env, threads=args.threads)
            if os.path.exists(target_op_blast) and os.path.getsize(target_op_blast) > 0:
                b_top = pd.read_csv(target_op_blast, sep="\t", header=None)
                valid_top = b_top[(b_top[2] >= 28.0) & (b_top[10] <= 1e-4)]
                for _, b_row in valid_top.iterrows():
                    q_id = str(b_row[0])
                    s_id = str(b_row[1])
                    target_op_pids.add(s_id)
                    ko_m = re.search(r"(K\d{5})", q_id)
                    if ko_m and s_id not in target_pid_to_ko:
                        target_pid_to_ko[s_id] = ko_m.group(1)
            print(f"[{args.target_name}] Direct seed operon search against target proteome detected {len(target_op_pids)} candidate operon genes.", flush=True)
        except Exception as e:
            print(f"Warning during target operon search: {e}", file=sys.stderr)

    if not target_op_pids and operon_records and target_faa_dict:
        for op_r in operon_records:
            ko_m = re.search(r"(K\d{5})", op_r.id) or re.search(r"(K\d{5})", op_r.description)
            r_ko = ko_m.group(1) if ko_m else ""
            for t_id, t_seq in target_faa_dict.items():
                if op_r.id in t_id or t_id in op_r.id or (op_r.seq and t_seq and str(op_r.seq) == str(t_seq)):
                    target_op_pids.add(t_id)
                    if r_ko and t_id not in target_pid_to_ko:
                        target_pid_to_ko[t_id] = r_ko

    for p_id, ko_val in target_pid_to_ko.items():
        if p_id not in target_anno_dict:
            target_anno_dict[p_id] = {
                "ko": ko_val,
                "score": 500.0,
                "threshold": "100.0",
                "evalue": 1e-50,
                "definition": ko_definitions.get(ko_val, {}).get("definition", "")
            }

    ref_candidates = {}
    flanks_query_records = []
    seen_qids = set()

    for ref_idx, (_, r) in enumerate(pos_refs.iterrows()):
        acc = str(r["accession"])
        clean_acc = acc.replace("GB_", "").replace("RS_", "")
        ref_faa = str(r.get("proteome_faa", ""))
        ref_genome_fa = str(r.get("genome_fasta", ""))
        if not ref_faa or not os.path.exists(ref_faa):
            continue

        ref_records = load_fasta_records(ref_faa)
        if not ref_records:
            continue

        ref_coords = parse_fasta_header_coordinates(ref_faa)
        original_faa = str(r.get("original_proteome_faa", ""))
        ref_gff = discover_companion_gff(original_faa or ref_faa, acc)
        if ref_gff:
            ref_coords.update(parse_gff_coordinates(ref_gff))

        # Locate operon genes in this reference and map to KO
        op_pids = set()
        ref_pid_to_ko = {}

        # 1. From Stage 0 ko_proteins if recorded in refs_df
        ko_prots_raw = r.get("ko_proteins", "")
        if pd.notna(ko_prots_raw) and ko_prots_raw:
            try:
                ko_dict = json.loads(str(ko_prots_raw)) if isinstance(ko_prots_raw, str) else ko_prots_raw
                if isinstance(ko_dict, dict):
                    for ko_k, pid_v in ko_dict.items():
                        op_pids.add(str(pid_v))
                        ref_pid_to_ko[str(pid_v)] = str(ko_k)
            except Exception as e:
                print(f"Notice: could not parse ko_proteins for {acc}: {e}", file=sys.stderr)

        # 2. DIAMOND search against seed_operon_faa (only if not already resolved via Stage 0 KO map)
        tmp_op_blast = os.path.join(args.out_refs_locus_dir, f"{clean_acc}_vs_seed_operon.tsv")
        if not op_pids and args.seed_operon_faa and os.path.exists(args.seed_operon_faa):
            try:
                run_diamond_search(args.seed_operon_faa, ref_faa, tmp_op_blast, args.diamond_env, threads=args.threads)
                if os.path.exists(tmp_op_blast) and os.path.getsize(tmp_op_blast) > 0:
                    b_op = pd.read_csv(tmp_op_blast, sep="\t", header=None)
                    valid_op = b_op[(b_op[2] >= 28.0) & (b_op[10] <= 1e-4)]
                    for _, b_row in valid_op.iterrows():
                        q_id = str(b_row[0])
                        s_id = str(b_row[1])
                        op_pids.add(s_id)
                        ko_m = re.search(r"(K\d{5})", q_id)
                        if ko_m and s_id not in ref_pid_to_ko:
                            ref_pid_to_ko[s_id] = ko_m.group(1)
            except Exception as e:
                print(f"Warning during operon search in {acc}: {e}", file=sys.stderr)

        # 3. Fallback: match by sequence or ID against seed_operon records
        if not op_pids:
            for op_r in operon_records:
                ko_m = re.search(r"(K\d{5})", op_r.id) or re.search(r"(K\d{5})", op_r.description)
                r_ko = ko_m.group(1) if ko_m else ""
                for rec in ref_records:
                    if op_r.id in rec.id or rec.id in op_r.id or (op_r.seq and rec.seq and str(op_r.seq) == str(rec.seq)):
                        op_pids.add(rec.id)
                        if r_ko and rec.id not in ref_pid_to_ko:
                            ref_pid_to_ko[rec.id] = r_ko

        # Ensure op_pids can match ref_records even if prefixes or whitespaces differ
        ref_pids_set = {rec.id for rec in ref_records}
        if op_pids and not any(pid in ref_pids_set for pid in op_pids):
            for rec in ref_records:
                for op_p in list(op_pids):
                    if op_p == rec.id or op_p in rec.id or rec.id in op_p:
                        op_pids.add(rec.id)
                        if op_p in ref_pid_to_ko and rec.id not in ref_pid_to_ko:
                            ref_pid_to_ko[rec.id] = ref_pid_to_ko[op_p]

        # Group by contig and find primary operon contig
        ctg_to_prots = group_reference_proteins_by_coordinate(ref_records, ref_coords)

        ctg_op_counts = {c: sum(1 for _, rec in prots if rec.id in op_pids) for c, prots in ctg_to_prots.items()}
        if not ctg_op_counts or max(ctg_op_counts.values()) == 0:
            if ctg_to_prots:
                best_op_ctg = max(ctg_to_prots.keys(), key=lambda c: len(ctg_to_prots[c]))
                prots_on_ctg = ctg_to_prots[best_op_ctg]
                mid_idx = len(prots_on_ctg) // 2
                min_idx = max(0, mid_idx - 2)
                max_idx = min(len(prots_on_ctg) - 1, mid_idx + 2)
            else:
                continue
        else:
            best_op_ctg = max(ctg_op_counts.keys(), key=lambda c: ctg_op_counts[c])
            prots_on_ctg = ctg_to_prots[best_op_ctg]
            op_indices_on_ctg = [idx for idx, (_, rec) in enumerate(prots_on_ctg) if rec.id in op_pids]
            min_idx = min(op_indices_on_ctg)
            max_idx = max(op_indices_on_ctg)

        n_flank = args.n_flank_genes
        u_prots = [rec for _, rec in prots_on_ctg[max(0, min_idx - n_flank) : min_idx]]
        core_prots = [rec for _, rec in prots_on_ctg[min_idx : max_idx + 1]]
        d_prots = [rec for _, rec in prots_on_ctg[max_idx + 1 : min(len(prots_on_ctg), max_idx + 1 + n_flank)]]

        raw_ani = r.get("ANI")
        if raw_ani is not None and pd.notna(raw_ani) and str(raw_ani).strip() not in ["NA", "nan", ""]:
            try:
                r_ani = float(raw_ani)
                ani_weight = max(0.5, r_ani / 100.0)
            except (ValueError, TypeError):
                r_ani = 0.0
                ani_weight = 1.0
        else:
            r_ani = 0.0
            ani_weight = 1.0

        r_tier = str(r.get("tier", "Tier_B"))
        is_single_ctg_op = bool(ctg_op_counts and max(ctg_op_counts.values()) >= len(op_pids))
        tier_weight = 1.5 if "Tier_A" in r_tier else (1.2 if "Tier_B" in r_tier else 1.0)
        single_contig_mult = 1.3 if is_single_ctg_op else 1.0
        quality_mult = tier_weight * ani_weight * single_contig_mult

        cand_data = {
            "acc": acc, "clean_acc": clean_acc, "best_ctg": best_op_ctg,
            "ref_faa": ref_faa, "ref_genome_fa": ref_genome_fa, "ref_coords": ref_coords,
            "u_prots": u_prots, "core_prots": core_prots, "d_prots": d_prots,
            "op_pids": op_pids,
            "pid_to_ko": ref_pid_to_ko,
            "ref_records": ref_records,
            "ref_faa_dict": {rec.id: str(rec.seq) for rec in ref_records},
            "ani": r_ani,
            "tier": r_tier,
            "single_contig_operon": is_single_ctg_op,
            "quality_mult": quality_mult
        }
        ref_candidates[clean_acc] = cand_data

        for u in u_prots:
            qid = f"UP__{clean_acc}__{u.id}"
            if qid not in seen_qids:
                flanks_query_records.append(SeqRecord(u.seq, id=qid, description=f"side=UP|ref={clean_acc}"))
                seen_qids.add(qid)
        for d in d_prots:
            qid = f"DOWN__{clean_acc}__{d.id}"
            if qid not in seen_qids:
                flanks_query_records.append(SeqRecord(d.seq, id=qid, description=f"side=DOWN|ref={clean_acc}"))
                seen_qids.add(qid)

    # Calculate reference cohort median homologous interanchor span (P1-1)
    ref_cohort_interanchor_spans = []
    for clean_acc, cand_data in ref_candidates.items():
        u_p = cand_data.get("u_prots", [])
        d_p = cand_data.get("d_prots", [])
        r_coords = cand_data.get("ref_coords", {})
        if u_p and d_p and r_coords:
            u_innermost = u_p[-1].id
            d_innermost = d_p[0].id
            if u_innermost in r_coords and d_innermost in r_coords:
                u_inf = r_coords[u_innermost]
                d_inf = r_coords[d_innermost]
                if u_inf.get("contig") == d_inf.get("contig"):
                    if u_inf["start"] <= d_inf["start"]:
                        span = max(0, d_inf["start"] - u_inf["end"])
                    else:
                        span = max(0, u_inf["start"] - d_inf["end"])
                    if span > 0:
                        ref_cohort_interanchor_spans.append(span)

    expected_interanchor_span = int(pd.Series(ref_cohort_interanchor_spans).median()) if ref_cohort_interanchor_spans else int(args.expected_locus_length)
    if expected_interanchor_span <= 0:
        expected_interanchor_span = int(args.expected_locus_length)

    # Fallback to consensus anchor symbols if no flank records
    if not flanks_query_records:
        for sid, akey in [("UP__REF__Anchor_L1", "L1"), ("UP__REF__Anchor_L2", "L2"), ("DOWN__REF__Anchor_R1", "R1"), ("DOWN__REF__Anchor_R2", "R2")]:
            sym = anchor_dict.get(akey, {}).get("gene_symbol", "")
            if sym:
                flanks_query_records.append(SeqRecord(Seq("MKVFISLLLLASV"), id=sid, description=f"symbol={sym}"))

    # 4. Search Reference Flanks against Target Proteome
    flanks_query_fa = os.path.join(args.out_refs_locus_dir, "all_refs_flanks_query.faa")
    with open(flanks_query_fa, "w") as fq:
        SeqIO.write(flanks_query_records, fq, "fasta")

    tmp_blast = os.path.join(args.out_refs_locus_dir, "all_flanks_vs_target.tsv")
    if os.path.exists(args.target_faa) and os.path.exists(flanks_query_fa):
        run_diamond_search(flanks_query_fa, args.target_faa, tmp_blast, args.diamond_env, threads=args.threads)

    b_df = pd.read_csv(tmp_blast, sep="\t", header=None) if (os.path.exists(tmp_blast) and os.path.getsize(tmp_blast) > 0) else pd.DataFrame()

    parsed_hits = []
    if not b_df.empty:
        for _, mrow in b_df.iterrows():
            qid = str(mrow[0])
            t_pid = str(mrow[1])
            pident = float(mrow[2])
            bitscore = float(mrow[11])
            if pident < 25.0 or t_pid not in gff_coords:
                continue

            c_inf = gff_coords[t_pid]
            parts = qid.split("__", 2)
            side = parts[0]
            ref_acc = parts[1] if len(parts) > 1 else "unknown"
            orig_pid = parts[2] if len(parts) > 2 else t_pid
            parsed_hits.append({
                "qid": qid, "side": side, "ref_acc": ref_acc, "orig_pid": orig_pid,
                "t_pid": t_pid, "contig": c_inf["contig"], "start": c_inf["start"],
                "end": c_inf["end"], "strand": c_inf["strand"], "identity": pident,
                "bitscore": bitscore
            })

    # 5. Collinear Synteny Clustering per Reference and Contig
    ref_evaluations = {}
    ref_ctg_evaluations = {}
    max_syntenic_gap = 35000

    for r_acc in ref_candidates.keys():
        r_hits = [h for h in parsed_hits if h["ref_acc"] == r_acc]
        ref_ctg_evaluations[r_acc] = {}
        if not r_hits:
            ref_evaluations[r_acc] = {"score": 0.0, "best_cluster": [], "best_contig": "", "n_u": 0, "n_d": 0, "has_bilateral": False}
            continue

        ctg_hits = {}
        for h in r_hits:
            ctg_hits.setdefault(h["contig"], []).append(h)

        best_ctg_score = -1e9
        best_ctg_cluster = []
        best_ctg_id = ""

        for ctg_id, c_list in ctg_hits.items():
            s_hits = sorted(c_list, key=lambda x: x["start"])
            sub_clusters = []
            curr = []
            for h in s_hits:
                if not curr:
                    curr.append(h)
                else:
                    if (h["start"] - curr[-1]["end"]) <= max_syntenic_gap:
                        curr.append(h)
                    else:
                        sub_clusters.append(curr)
                        curr = [h]
            if curr:
                sub_clusters.append(curr)

            # Evaluate individual sub-clusters as well as bilateral bridging pairs across missing operon gaps
            candidate_clusters = list(sub_clusters)
            max_breakpoint_gap = max(450000, 2 * int(getattr(args, "flank_bp", 30000)) + int(getattr(args, "expected_locus_length", 15000)))
            for i in range(len(sub_clusters)):
                for j in range(i + 1, len(sub_clusters)):
                    c1 = sub_clusters[i]
                    c2 = sub_clusters[j]
                    has_u = any(h["side"] == "UP" for h in c1) or any(h["side"] == "UP" for h in c2)
                    has_d = any(h["side"] == "DOWN" for h in c1) or any(h["side"] == "DOWN" for h in c2)
                    if has_u and has_d:
                        gap = c2[0]["start"] - c1[-1]["end"]
                        if gap <= max_breakpoint_gap:
                            candidate_clusters.append(c1 + c2)

            ctg_best_score = -1e9
            ctg_best_cl = []
            ctg_best_u = 0
            ctg_best_d = 0

            cand_ref = ref_candidates.get(r_acc, {})
            q_mult = cand_ref.get("quality_mult", 1.0)

            for cl in candidate_clusters:
                u_hits = [h for h in cl if h["side"] == "UP"]
                d_hits = [h for h in cl if h["side"] == "DOWN"]
                n_u = len(set(h["orig_pid"] for h in u_hits))
                n_d = len(set(h["orig_pid"] for h in d_hits))
                bit_sum = sum(h["bitscore"] for h in cl)
                has_bilateral = (n_u >= 2 and n_d >= 2)
                bilateral_bonus = 10000.0 if (n_u >= 3 and n_d >= 3) else (5000.0 if has_bilateral else 0.0)
                score = (bit_sum + bilateral_bonus + 100.0 * (n_u + n_d)) * q_mult

                if score > ctg_best_score:
                    ctg_best_score = score
                    ctg_best_cl = cl
                    ctg_best_u = n_u
                    ctg_best_d = n_d

            ref_ctg_evaluations[r_acc][ctg_id] = {
                "score": ctg_best_score,
                "cluster": ctg_best_cl,
                "n_u": ctg_best_u,
                "n_d": ctg_best_d,
                "has_bilateral": (ctg_best_u >= 2 and ctg_best_d >= 2)
            }

            if ctg_best_score > best_ctg_score:
                best_ctg_score = ctg_best_score
                best_ctg_cluster = ctg_best_cl
                best_ctg_id = ctg_id

        u_cnt = len(set(h["orig_pid"] for h in best_ctg_cluster if h["side"] == "UP"))
        d_cnt = len(set(h["orig_pid"] for h in best_ctg_cluster if h["side"] == "DOWN"))
        ref_evaluations[r_acc] = {
            "score": best_ctg_score,
            "best_cluster": best_ctg_cluster,
            "best_contig": best_ctg_id,
            "n_u": u_cnt,
            "n_d": d_cnt,
            "has_bilateral": (u_cnt >= 2 and d_cnt >= 2)
        }

    # 6. Hierarchical Lexicographic Ranking of Target Contigs (Localization-First & Top-N Tracking)
    def find_target_record(ctg_name, records_dict):
        return find_matching_record(ctg_name, records_dict)

    cand_contigs = list(set(h["contig"] for h in parsed_hits if h.get("contig")))
    contig_geometry = {}
    contig_geom_rank = {}
    contig_conserved_genes = {}
    contig_supporting_clusters = {}
    contig_converged_clusters = {}
    contig_syntenic_clusters = {}
    contig_order_score = {}
    contig_total_score = {}
    contig_cluster_scores = {}
    contig_loc_score = {}
    contig_bilat_score = {}
    contig_has_operon = {}
    contig_boundary_anchors = {}
    contig_ancestral_intervals = {}
    contig_observed_spans = {}
    contig_locus_hits = {}
    contig_locus_refs = {}
    contig_locus_type = {}

    target_ko_list = [k.strip() for k in getattr(args, "target_kos", "").split(",") if k.strip()]
    target_ko_set = set(target_ko_list)

    for ctg_id in cand_contigs:
        # References with valid syntenic cluster on this contig
        syn_refs = [
            r for r in ref_candidates.keys()
            if ref_ctg_evaluations.get(r, {}).get(ctg_id, {}).get("score", 0.0) > 0 and
               ref_ctg_evaluations.get(r, {}).get(ctg_id, {}).get("cluster")
        ]
        # References whose global best contig is this contig
        conv_refs = [
            r for r, ev in ref_evaluations.items()
            if is_same_contig(ev.get("best_contig"), ctg_id) and ev.get("score", 0.0) > 0
        ]

        # Independent evolutionary clusters converging on this contig
        conv_clusters = set(ref_to_cluster.get(r.replace("GB_", "").replace("RS_", "").strip(), r) for r in conv_refs)
        syn_clusters = set(ref_to_cluster.get(r.replace("GB_", "").replace("RS_", "").strip(), r) for r in syn_refs)

        contig_converged_clusters[ctg_id] = len(conv_clusters)
        contig_syntenic_clusters[ctg_id] = len(syn_clusters)
        contig_supporting_clusters[ctg_id] = len(conv_clusters)

        # Conserved genes strictly within syntenic clusters (filtering unclustered spurious hits)
        syn_pids = set()
        for r in syn_refs:
            cl = ref_ctg_evaluations[r][ctg_id].get("cluster", [])
            for h in cl:
                if h.get("t_pid"):
                    syn_pids.add(h["t_pid"])
        contig_conserved_genes[ctg_id] = len(syn_pids)

        # Check for target operon genes directly present on this contig
        target_op_cds = []
        for p, inf in gff_coords.items():
            if inf.get("contig") == ctg_id or is_same_contig(inf.get("contig"), ctg_id):
                clean_p = p.split("___")[-1] if "___" in p else p
                if (p in target_op_pids or clean_p in target_op_pids or
                    p in target_pid_to_ko or clean_p in target_pid_to_ko):
                    target_op_cds.append(inf)
                    continue
                t_anno = lookup_target_anno(p, target_anno_dict)
                ko_val = t_anno.get("ko") if t_anno else (inf.get("ko") or target_pid_to_ko.get(p) or target_pid_to_ko.get(clean_p))
                if ko_val and (not target_ko_set or ko_val in target_ko_set):
                    target_op_cds.append(inf)
        has_op = len(target_op_cds) > 0
        contig_has_operon[ctg_id] = has_op

        # Collect all syntenic hits on this contig across the reference cohort
        syn_ctg_hits = [h for r in syn_refs for h in ref_ctg_evaluations[r][ctg_id].get("cluster", [])]
        syn_u = [h for h in syn_ctg_hits if h.get("side") == "UP"]
        syn_d = [h for h in syn_ctg_hits if h.get("side") == "DOWN"]

        def group_syntenic_hits(hit_list, max_gap):
            if not hit_list:
                return []
            s_hits = sorted(hit_list, key=lambda x: x["start"])
            grps = []
            cur = [s_hits[0]]
            for h in s_hits[1:]:
                if (h["start"] - cur[-1]["end"]) <= max_gap:
                    cur.append(h)
                else:
                    grps.append(cur)
                    cur = [h]
            if cur:
                grps.append(cur)
            return grps

        u_groups = group_syntenic_hits(syn_u, max_syntenic_gap)
        d_groups = group_syntenic_hits(syn_d, max_syntenic_gap)

        eval_rec = find_matching_record(ctg_id, target_records)
        ctg_len = len(eval_rec.seq) if eval_rec else 100000

        max_cohort_gap = max(450000, 2 * int(getattr(args, "flank_bp", 30000)) + int(getattr(args, "expected_locus_length", 15000)))

        candidate_loci = []

        # 1. Evaluate Bilateral Candidate Pairs (Across Cohort Flank Fusion or Single-Ref Bilateral)
        for u_grp in u_groups:
            u_s = min(h["start"] for h in u_grp)
            u_e = max(h["end"] for h in u_grp)
            n_u_pids = len(set(h["orig_pid"] for h in u_grp))
            u_refs = set(h["ref_acc"] for h in u_grp)

            for d_grp in d_groups:
                d_s = min(h["start"] for h in d_grp)
                d_e = max(h["end"] for h in d_grp)
                n_d_pids = len(set(h["orig_pid"] for h in d_grp))
                d_refs = set(h["ref_acc"] for h in d_grp)

                if u_e <= d_s:
                    inner_gap = d_s - u_e
                    int_s, int_e = u_e, d_s
                    left_h = max(u_grp, key=lambda h: h["end"])
                    right_h = min(d_grp, key=lambda h: h["start"])
                    collinear = True
                elif d_e <= u_s:
                    inner_gap = u_s - d_e
                    int_s, int_e = d_e, u_s
                    left_h = max(d_grp, key=lambda h: h["end"])
                    right_h = min(u_grp, key=lambda h: h["start"])
                    collinear = False
                else:
                    inner_gap = 0
                    int_s = min(u_s, d_s)
                    int_e = max(u_e, d_e)
                    left_h = min(u_grp + d_grp, key=lambda h: h["start"])
                    right_h = max(u_grp + d_grp, key=lambda h: h["end"])
                    collinear = True

                op_in_int = any(c["end"] >= int_s and c["start"] <= int_e for c in target_op_cds)
                allowed_gap = max_cohort_gap * 2 if op_in_int else max_cohort_gap

                if inner_gap <= allowed_gap:
                    all_refs = u_refs | d_refs
                    ref_cl_set = set(ref_to_cluster.get(r.replace("GB_", "").replace("RS_", "").strip(), r) for r in all_refs)
                    tot_bit = sum(h["bitscore"] for h in u_grp + d_grp)
                    q_weights = [ref_candidates.get(r, {}).get("quality_mult", 1.0) for r in all_refs]
                    avg_q = sum(q_weights) / max(1.0, len(q_weights))

                    bilat_sc = 10000.0 if (n_u_pids >= 2 and n_d_pids >= 2) else 5000.0
                    op_bonus = 25000.0 if op_in_int else 0.0
                    conv_bonus = 4000.0 * len(ref_cl_set)
                    coll_bonus = 1000.0 if collinear else 0.0

                    sc = (tot_bit * avg_q) + bilat_sc + op_bonus + conv_bonus + coll_bonus + 100.0 * (n_u_pids + n_d_pids) - (0.02 * inner_gap)

                    candidate_loci.append({
                        "type": "BILATERAL",
                        "score": sc,
                        "hits": u_grp + d_grp,
                        "refs": all_refs,
                        "left_anchor": left_h.get("t_pid", "NONE"),
                        "right_anchor": right_h.get("t_pid", "NONE"),
                        "left_hit": left_h,
                        "right_hit": right_h,
                        "interval_start": int_s,
                        "interval_end": int_e,
                        "observed_span": inner_gap,
                        "n_u": n_u_pids,
                        "n_d": n_d_pids,
                        "collinear": collinear,
                        "has_op_in_interval": op_in_int
                    })

        # 2. Evaluate Unilateral Candidates (If No Bilateral Exists or as Lower-Tier Options)
        for u_grp in u_groups:
            u_s = min(h["start"] for h in u_grp)
            u_e = max(h["end"] for h in u_grp)
            n_u_pids = len(set(h["orig_pid"] for h in u_grp))
            u_refs = set(h["ref_acc"] for h in u_grp)
            q_weights = [ref_candidates.get(r, {}).get("quality_mult", 1.0) for r in u_refs]
            avg_q = sum(q_weights) / max(1.0, len(q_weights))
            tot_bit = sum(h["bitscore"] for h in u_grp)
            u_bound_h = max(u_grp, key=lambda h: h["end"])
            max_ani = max((float(ref_candidates.get(r, {}).get("ani", 80.0)) for r in u_refs), default=80.0)
            has_tier_a = any("Tier_A" in ref_candidates.get(r, {}).get("tier", "") for r in u_refs)
            tier_bonus = 3000.0 if has_tier_a else 0.0
            candidate_loci.append({
                "type": "UNILATERAL_LEFT",
                "score": (tot_bit * avg_q) + 100.0 * n_u_pids + tier_bonus + (max_ani * 30.0),
                "hits": u_grp,
                "refs": u_refs,
                "left_anchor": u_bound_h.get("t_pid", "NONE"),
                "right_anchor": "NONE",
                "left_hit": u_bound_h,
                "right_hit": None,
                "interval_start": u_e,
                "interval_end": min(ctg_len, u_e + int(getattr(args, "expected_locus_length", 15000))),
                "observed_span": 0,
                "n_u": n_u_pids,
                "n_d": 0,
                "collinear": True,
                "has_op_in_interval": False
            })

        for d_grp in d_groups:
            d_s = min(h["start"] for h in d_grp)
            d_e = max(h["end"] for h in d_grp)
            n_d_pids = len(set(h["orig_pid"] for h in d_grp))
            d_refs = set(h["ref_acc"] for h in d_grp)
            q_weights = [ref_candidates.get(r, {}).get("quality_mult", 1.0) for r in d_refs]
            avg_q = sum(q_weights) / max(1.0, len(q_weights))
            tot_bit = sum(h["bitscore"] for h in d_grp)
            d_bound_h = min(d_grp, key=lambda h: h["start"])
            max_ani = max((float(ref_candidates.get(r, {}).get("ani", 80.0)) for r in d_refs), default=80.0)
            has_tier_a = any("Tier_A" in ref_candidates.get(r, {}).get("tier", "") for r in d_refs)
            tier_bonus = 3000.0 if has_tier_a else 0.0
            candidate_loci.append({
                "type": "UNILATERAL_RIGHT",
                "score": (tot_bit * avg_q) + 100.0 * n_d_pids + tier_bonus + (max_ani * 30.0),
                "hits": d_grp,
                "refs": d_refs,
                "left_anchor": "NONE",
                "right_anchor": d_bound_h.get("t_pid", "NONE"),
                "left_hit": None,
                "right_hit": d_bound_h,
                "interval_start": max(0, d_s - int(getattr(args, "expected_locus_length", 15000))),
                "interval_end": d_s,
                "observed_span": 0,
                "n_u": 0,
                "n_d": n_d_pids,
                "collinear": True,
                "has_op_in_interval": False
            })

        if candidate_loci:
            best_cand_locus = max(candidate_loci, key=lambda x: (
                1 if x["type"] == "BILATERAL" else 0,
                1 if x.get("has_op_in_interval") else 0,
                x["score"]
            ))
        else:
            best_cand_locus = {
                "type": "NONE", "score": 0.0, "hits": [], "refs": set(),
                "left_anchor": "NONE", "right_anchor": "NONE",
                "left_hit": None, "right_hit": None,
                "interval_start": 0, "interval_end": 0, "observed_span": 0,
                "n_u": 0, "n_d": 0, "collinear": True, "has_op_in_interval": False
            }

        is_bilat = (best_cand_locus["type"] == "BILATERAL")
        is_unilat = best_cand_locus["type"].startswith("UNILATERAL")

        if is_bilat:
            geom = "BILATERAL_SAME_CONTIG"
            geom_rank = 3
        elif is_unilat:
            geom = "UNILATERAL_SAME_CONTIG"
            geom_rank = 2
        elif syn_refs:
            geom = "UNILATERAL_SAME_CONTIG"
            geom_rank = 1
        else:
            geom = "SCATTERED_HITS"
            geom_rank = 0

        contig_geometry[ctg_id] = geom
        contig_geom_rank[ctg_id] = geom_rank

        o_score = 50 if best_cand_locus["collinear"] and (best_cand_locus["n_u"] > 0 and best_cand_locus["n_d"] > 0) else 0
        contig_order_score[ctg_id] = o_score

        c_score = 0.0
        c_cl_scores = {}
        for c_k, r_list in ref_clusters.items():
            k_best = 0.0
            for r_a in r_list:
                r_sc = ref_ctg_evaluations.get(r_a, {}).get(ctg_id, {}).get("score", 0.0)
                if r_sc > k_best:
                    k_best = r_sc
            c_score += k_best
            c_cl_scores[c_k] = k_best
        contig_total_score[ctg_id] = c_score
        contig_cluster_scores[ctg_id] = c_cl_scores

        # Calculate locator_score and bilateral_context_score (P0-2, P1-1)
        n_u = best_cand_locus["n_u"]
        n_d = best_cand_locus["n_d"]
        flank_norm = min(5.0, float(getattr(args, "n_flank_genes", 25)))
        L_ratio = min(1.0, n_u / max(1.0, flank_norm))
        R_ratio = min(1.0, n_d / max(1.0, flank_norm))
        O_order = 1.0 if best_cand_locus["collinear"] else 0.5
        S_cohort = len(set(ref_to_cluster.get(r.replace("GB_", "").replace("RS_", "").strip(), r) for r in best_cand_locus["refs"])) / max(1.0, float(n_eff))

        if geom == "BILATERAL_SAME_CONTIG":
            loc_score = round(math.sqrt(max(1e-5, L_ratio * R_ratio)) * O_order * S_cohort, 4)
            bilat_score = loc_score
        elif geom == "UNILATERAL_SAME_CONTIG":
            loc_score = round(0.6 * max(L_ratio, R_ratio) * O_order * S_cohort, 4)
            bilat_score = 0.0
        else:
            loc_score = round(0.5 * min(L_ratio, R_ratio) * O_order * S_cohort, 4)
            bilat_score = 0.0

        contig_loc_score[ctg_id] = loc_score
        contig_bilat_score[ctg_id] = bilat_score

        # Identify boundary anchors and ancestral interval (P0-1)
        contig_boundary_anchors[ctg_id] = (best_cand_locus["left_anchor"], best_cand_locus["right_anchor"])
        contig_ancestral_intervals[ctg_id] = (best_cand_locus["interval_start"], best_cand_locus["interval_end"])
        contig_observed_spans[ctg_id] = best_cand_locus["observed_span"]
        contig_locus_hits[ctg_id] = best_cand_locus["hits"]
        contig_locus_refs[ctg_id] = best_cand_locus["refs"]
        contig_locus_type[ctg_id] = best_cand_locus["type"]

    def compute_expected_interanchor_span_for_pair(u_rank, d_rank, default_len=15000):
        cluster_spans = []
        for c_id, r_list in ref_clusters.items():
            c_vals = []
            for r_acc in r_list:
                clean_acc = r_acc.replace("GB_", "").replace("RS_", "")
                cand = ref_candidates.get(clean_acc)
                if not cand:
                    continue
                u_prots = cand.get("u_prots", [])
                d_prots = cand.get("d_prots", [])
                r_coords = cand.get("ref_coords", {})
                if not u_prots or not d_prots or not r_coords:
                    continue
                u_idx = len(u_prots) - 1 - u_rank
                d_idx = d_rank
                if 0 <= u_idx < len(u_prots) and 0 <= d_idx < len(d_prots):
                    u_p = u_prots[u_idx]
                    d_p = d_prots[d_idx]
                    if u_p.id in r_coords and d_p.id in r_coords:
                        u_c = r_coords[u_p.id]
                        d_c = r_coords[d_p.id]
                        if u_c.get("contig") == d_c.get("contig"):
                            if u_c["start"] <= d_c["start"]:
                                sp = max(0, d_c["start"] - u_c["end"])
                            else:
                                sp = max(0, u_c["start"] - d_c["end"])
                            if sp > 0:
                                c_vals.append(sp)
            if c_vals:
                cluster_spans.append(float(pd.Series(c_vals).median()))
        if cluster_spans:
            return int(pd.Series(cluster_spans).median())
        return int(default_len)

    def get_contig_ranking_key(c):
        n_conv = contig_converged_clusters.get(c, 0)
        n_syn = contig_syntenic_clusters.get(c, 0)
        g_rank = contig_geom_rank.get(c, 0)
        tot_sc = contig_total_score.get(c, 0.0)
        n_genes = contig_conserved_genes.get(c, 0)
        o_sc = contig_order_score.get(c, 0)
        bilat_sc = contig_bilat_score.get(c, 0.0)
        loc_sc = contig_loc_score.get(c, 0.0)

        # Multi-cluster tier (Pure target-blind syntenic ranking):
        # Tier 3: Strong multi-cluster convergence (>= 2 independent evolutionary clusters converge)
        # Tier 2: Single cluster with bilateral conservation (1 cluster, bilateral)
        # Tier 1: Single cluster with unilateral conservation OR >= 2 secondary syntenic clusters
        # Tier 0: Weak / scattered hits
        if n_conv >= 2:
            tier = 3
        elif n_conv == 1 and g_rank >= 3:
            tier = 2
        elif n_conv == 1 or n_syn >= 2:
            tier = 1
        else:
            tier = 0

        return (
            tier,
            n_conv,
            g_rank,
            bilat_sc,
            loc_sc,
            tot_sc,
            n_syn,
            n_genes,
            o_sc
        )

    ranked_contigs = []
    if cand_contigs:
        ranked_contigs = sorted(
            cand_contigs,
            key=get_contig_ranking_key,
            reverse=True
        )
        chosen_contig_id = ranked_contigs[0]
        other_scores = [contig_total_score[c] for c in cand_contigs if c != chosen_contig_id]
        max_other_score = max(other_scores) if other_scores else 0.0
        chosen_geometry = contig_geometry[chosen_contig_id]
        dominance_ratio = contig_total_score[chosen_contig_id] / (max_other_score + 1e-5) if max_other_score > 0 else 999.0
    else:
        chosen_contig_id = list(target_records.keys())[0] if target_records else ""
        chosen_geometry = "NONE"
        dominance_ratio = 1.0

    # Export Top-N Ranked Candidates Table (P0-2, Target-Blind)
    max_cands = int(getattr(args, "max_candidate_loci", 5))
    top_cands = ranked_contigs[:max_cands] if cand_contigs else [chosen_contig_id]
    cand_records = []
    contig_matched_expected_spans = {}

    for rank_idx, c_id in enumerate(top_cands, 1):
        if not c_id:
            continue
        c_geom = contig_geometry.get(c_id, "NONE")
        c_loc = contig_loc_score.get(c_id, 0.0)
        c_bil = contig_bilat_score.get(c_id, 0.0)
        c_l_anc, c_r_anc = contig_boundary_anchors.get(c_id, ("NONE", "NONE"))
        c_int_s, c_int_e = contig_ancestral_intervals.get(c_id, (0, 0))
        c_obs_sp = contig_observed_spans.get(c_id, 0)
        
        # Calculate pair-specific expected interanchor span if bilateral
        c_syn_refs = [
            r for r in ref_candidates.keys()
            if ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("score", 0.0) > 0 and
               ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("cluster")
        ]
        c_u_h = [h for r in c_syn_refs for h in ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("cluster", []) if h.get("side") == "UP"]
        c_d_h = [h for r in c_syn_refs for h in ref_ctg_evaluations.get(r, {}).get(c_id, {}).get("cluster", []) if h.get("side") == "DOWN"]
        c_u_rank = 0
        c_d_rank = 0
        if c_u_h:
            u_last = max(c_u_h, key=lambda h: h["end"])
            u_ref = ref_candidates.get(u_last.get("ref_acc", ""), {})
            u_plist = [p.id for p in u_ref.get("u_prots", [])]
            if u_last.get("orig_pid") in u_plist:
                c_u_rank = max(0, len(u_plist) - 1 - u_plist.index(u_last.get("orig_pid")))
        if c_d_h:
            d_first = min(c_d_h, key=lambda h: h["start"])
            d_ref = ref_candidates.get(d_first.get("ref_acc", ""), {})
            d_plist = [p.id for p in d_ref.get("d_prots", [])]
            if d_first.get("orig_pid") in d_plist:
                c_d_rank = max(0, d_plist.index(d_first.get("orig_pid")))

        c_exp_sp = compute_expected_interanchor_span_for_pair(c_u_rank, c_d_rank, default_len=args.expected_locus_length)
        contig_matched_expected_spans[c_id] = c_exp_sp

        cand_records.append({
            "rank": rank_idx,
            "target_id": args.target_name,
            "module_id": args.module_id,
            "contig_id": c_id,
            "geometry": c_geom,
            "locator_score": c_loc,
            "bilateral_context_score": c_bil,
            "supporting_clusters": contig_supporting_clusters.get(c_id, 0),
            "converged_clusters": contig_converged_clusters.get(c_id, 0),
            "conserved_genes": contig_conserved_genes.get(c_id, 0),
            "left_boundary_anchor": c_l_anc,
            "right_boundary_anchor": c_r_anc,
            "ancestral_interval_start": c_int_s,
            "ancestral_interval_end": c_int_e,
            "observed_interanchor_span": c_obs_sp,
            "expected_interanchor_span": c_exp_sp,
            "is_primary": "YES" if rank_idx == 1 else "NO"
        })
    df_cands = pd.DataFrame(cand_records)
    cands_tsv_1 = os.path.join(os.path.dirname(args.out_topology_tsv), f"{args.target_name}_{args.module_id}_candidate_loci_ranked.tsv")
    cands_tsv_2 = os.path.join(os.path.dirname(args.out_topology_tsv), "candidate_loci_ranked.tsv")
    df_cands.to_csv(cands_tsv_1, sep="\t", index=False)
    df_cands.to_csv(cands_tsv_2, sep="\t", index=False)

    expected_interanchor_span = contig_matched_expected_spans.get(chosen_contig_id, expected_interanchor_span)

    if chosen_geometry == "BILATERAL_SAME_CONTIG":
        localization_status = "STRONG_CANDIDATE"
    elif chosen_geometry in ["UNILATERAL_SAME_CONTIG", "SPLIT_CONTIG"]:
        localization_status = "ASSEMBLY_UNRESOLVED_CANDIDATE"
    else:
        localization_status = "NO_LOCALIZABLE_LOCUS"

    refs_on_chosen = [r for r in ref_candidates.keys() if ref_ctg_evaluations.get(r, {}).get(chosen_contig_id, {}).get("score", 0) > 0]
    locus_refs = contig_locus_refs.get(chosen_contig_id, set())
    locus_hits = contig_locus_hits.get(chosen_contig_id, [])

    def ref_selection_key(r):
        cand_r = ref_candidates.get(r, {})
        in_loc = 1 if r in locus_refs else 0
        tier = cand_r.get("tier", "")
        tier_val = 3000.0 if "Tier_A" in tier else (1500.0 if "Tier_B" in tier else 0.0)
        ani_val = float(cand_r.get("ani", 80.0)) * 50.0
        single_ctg = 2000.0 if cand_r.get("single_contig_operon", False) else 0.0
        r_loc_bits = sum(h["bitscore"] for h in locus_hits if h.get("ref_acc") == r)
        r_tot_sc = ref_ctg_evaluations.get(r, {}).get(chosen_contig_id, {}).get("score", 0.0)
        return (in_loc, tier_val + ani_val + single_ctg, r_loc_bits, r_tot_sc)

    if refs_on_chosen:
        best_ref_acc = max(refs_on_chosen, key=ref_selection_key)
        best_eval = ref_ctg_evaluations[best_ref_acc][chosen_contig_id]
        best_cluster = [h for h in locus_hits if h.get("ref_acc") == best_ref_acc] or best_eval.get("cluster", [])
    elif ref_evaluations:
        ranked_refs_all = sorted(ref_evaluations.keys(), key=lambda k: ref_evaluations[k]["score"], reverse=True)
        best_ref_acc = ranked_refs_all[0]
        best_eval = ref_evaluations[best_ref_acc]
        best_cluster = best_eval.get("best_cluster", [])
    else:
        best_ref_acc = "NONE"
        best_eval = {"score": 0.0, "best_cluster": [], "best_contig": chosen_contig_id, "n_u": 0, "n_d": 0, "has_bilateral": False}
        best_cluster = []

    ranked_refs = sorted(
        ref_candidates.keys(),
        key=ref_selection_key if refs_on_chosen else (lambda k: ref_evaluations.get(k, {}).get("score", 0)),
        reverse=True
    )

    # Auxiliary Empirical MAX-Contig Permutation Null
    n_perms = int(getattr(args, "n_permutations", 1000))
    alpha_val = float(getattr(args, "significance_alpha", 0.01))
    p_emp, s_obs, s_null_mean, s_null_95, _, stat_support = compute_empirical_max_contig_null(
        target_records=target_records,
        gff_coords=gff_coords,
        ctg_scores=contig_total_score,
        winning_ctg_id=chosen_contig_id,
        cluster_scores=contig_cluster_scores.get(chosen_contig_id, {}),
        n_permutations=n_perms,
        alpha=alpha_val
    )
    if n_eff <= 2 and stat_support == "NOT_SIGNIFICANT":
        stat_support = "NOT_SIGNIFICANT (LOW_STATISTICAL_POWER)"

    print(
        f"[{args.target_name}] Hierarchical Lexicographic Ranking selected candidate contig: '{chosen_contig_id}' "
        f"({chosen_geometry}, status={localization_status}, conserved_genes={contig_conserved_genes.get(chosen_contig_id, 0)}, "
        f"supporting_clusters={contig_supporting_clusters.get(chosen_contig_id, 0)}/{n_eff}, "
        f"score={contig_total_score.get(chosen_contig_id, 0.0):.1f}, dominance_ratio={dominance_ratio:.2f}, "
        f"P_emp={p_emp:.4f}, support={stat_support})",
        flush=True
    )

    # 7. Determine Target Locus Boundaries
    chosen_contig = find_target_record(chosen_contig_id, target_records)
    if not chosen_contig:
        chosen_contig = list(target_records.values())[0]
        chosen_contig_id = chosen_contig.id
    else:
        chosen_contig_id = chosen_contig.id
    chosen_contig_len = len(chosen_contig.seq)

    active_hits = locus_hits if locus_hits else best_cluster
    c_min = min(h["start"] for h in active_hits) if active_hits else 0
    c_max = max(h["end"] for h in active_hits) if active_hits else chosen_contig_len

    # Select references for visualization (up to args.max_visualized_refs):
    # Strictly deduplicate by clean accession so RS_/GB_ variants never appear twice.
    visualized_refs = []
    seen_visualized = set()

    def add_visualized_ref(acc_raw):
        if not acc_raw:
            return
        c_acc = acc_raw.replace("GB_", "").replace("RS_", "")
        if c_acc in ref_candidates and c_acc not in seen_visualized:
            visualized_refs.append(c_acc)
            seen_visualized.add(c_acc)

    # Priority 1: References that directly support the chosen consensus syntenic locus
    for r_acc in ranked_refs:
        if r_acc in locus_refs:
            add_visualized_ref(r_acc)

    # Priority 2: References with syntenic hits matching the chosen target contig in the locus neighborhood
    for r_acc in ranked_refs:
        if r_acc not in ref_candidates:
            continue
        ev = ref_evaluations.get(r_acc, {})
        if is_same_contig(ev.get("best_contig"), chosen_contig_id) and ev.get("score", 0) > 0:
            r_cl = ev.get("best_cluster", [])
            if r_cl and active_hits:
                r_min = min(h["start"] for h in r_cl)
                r_max = max(h["end"] for h in r_cl)
                if r_max >= (c_min - 50000) and r_min <= (c_max + 50000):
                    add_visualized_ref(r_acc)
            else:
                add_visualized_ref(r_acc)

    # Priority 2: Other references with syntenic hits on the chosen contig
    for r_acc in ranked_refs:
        ev = ref_evaluations.get(r_acc, {})
        if is_same_contig(ev.get("best_contig"), chosen_contig_id) and ev.get("score", 0) > 0:
            add_visualized_ref(r_acc)

    # Priority 3: Other references in ranked_refs with syntenic hits on any contig
    for r_acc in ranked_refs:
        if ref_evaluations.get(r_acc, {}).get("score", 0) > 0:
            add_visualized_ref(r_acc)

    # Priority 4: Any remaining references in ranked_refs
    for r_acc in ranked_refs:
        add_visualized_ref(r_acc)

    # Priority 5: Any remaining candidates from ref_candidates
    for r_acc in ref_candidates.keys():
        add_visualized_ref(r_acc)

    # Priority 6: Any remaining confirmed positive references from Stage 0 (pos_refs)
    for _, r in pos_refs.iterrows():
        add_visualized_ref(str(r["accession"]))

    # Fallback if somehow still empty
    if not visualized_refs and best_ref_acc != "NONE":
        add_visualized_ref(best_ref_acc)
    if not visualized_refs and not pos_refs.empty:
        add_visualized_ref(str(pos_refs.iloc[0]["accession"]))
    if not visualized_refs and not pos_refs.empty:
        visualized_refs = [str(pos_refs.iloc[0]["accession"]).replace("GB_", "").replace("RS_", "")]

    # Cap by max_visualized_refs (default 5)
    visualized_refs = visualized_refs[:args.max_visualized_refs]

    if best_ref_acc == "NONE" and visualized_refs:
        best_ref_acc = visualized_refs[0]

    print(f"[{args.target_name}] References selected for Clinker visualization (max={args.max_visualized_refs}, contig={chosen_contig_id}): {visualized_refs}", flush=True)

    required_edge_margin = max(args.min_edge_bp, int(args.expected_locus_length * 1.5))

    # Pre-calculate reference locus spans and extents to guarantee target slice length >= ref locus length
    ref_locus_spans = []
    ref_u_spans = []
    ref_d_spans = []
    all_ref_accs = set(visualized_refs) | set(ref_candidates.keys())
    all_ref_pids = set()

    # Index proteins across all positive reference candidates
    for r_acc in ref_candidates.keys():
        cand = ref_candidates.get(r_acc)
        if not cand:
            continue
        all_ref_accs.add(r_acc)
        if cand.get("acc"):
            all_ref_accs.add(cand["acc"])
        if cand.get("clean_acc"):
            all_ref_accs.add(cand["clean_acc"])
        if "." in r_acc:
            all_ref_accs.add(r_acc.split(".")[0])

        r_coords = cand.get("ref_coords", {})
        r_ctg = cand.get("best_ctg", "")
        op_cds = [r_coords[p] for p in cand.get("op_pids", set()) if p in r_coords]
        if op_cds:
            r_op_mid = (min(c["start"] for c in op_cds) + max(c["end"] for c in op_cds)) / 2.0
        else:
            genes_in_ord = cand.get("u_prots", []) + cand.get("core_prots", []) + cand.get("d_prots", [])
            coords_found = [r_coords[p.id] for p in genes_in_ord if p.id in r_coords]
            r_op_mid = (min(c["start"] for c in coords_found) + max(c["end"] for c in coords_found)) / 2.0 if coords_found else 0.0

        for pid, inf in r_coords.items():
            if inf.get("contig") == r_ctg or inf.get("full_contig") == r_ctg:
                if abs((inf["start"] + inf["end"]) / 2.0 - r_op_mid) <= getattr(args, "extended_window_bp", 100000):
                    all_ref_pids.add(pid)

    for r_acc in visualized_refs:
        cand = ref_candidates.get(r_acc)
        if not cand:
            continue
        genes_in_order = cand["u_prots"] + cand["core_prots"] + cand["d_prots"]
        for p in genes_in_order:
            all_ref_pids.add(p.id)

        ref_coords = cand.get("ref_coords", {})
        has_coords = all(p.id in ref_coords for p in genes_in_order)
        if has_coords and genes_in_order:
            r_min = min(ref_coords[p.id]["start"] for p in genes_in_order)
            r_max = max(ref_coords[p.id]["end"] for p in genes_in_order)
            span = (r_max - r_min) + 6000
            ref_locus_spans.append(span)

            op_cds = [ref_coords[p.id] for p in cand["core_prots"] if p.id in ref_coords]
            if op_cds:
                op_s = min(c["start"] for c in op_cds)
                op_e = max(c["end"] for c in op_cds)
                ref_u_spans.append(max(0, op_s - r_min + 3000))
                ref_d_spans.append(max(0, r_max - op_e + 3000))
            else:
                ref_u_spans.append(span // 2)
                ref_d_spans.append(span // 2)
        else:
            span = sum(len(p.seq) * 3 + 120 for p in genes_in_order) if genes_in_order else (args.expected_locus_length + 2 * args.flank_bp)
            ref_locus_spans.append(span)
            u_s = sum(len(p.seq) * 3 + 120 for p in cand["u_prots"])
            d_s = sum(len(p.seq) * 3 + 120 for p in cand["d_prots"])
            ref_u_spans.append(u_s)
            ref_d_spans.append(d_s)

    max_ref_span = max(ref_locus_spans) if ref_locus_spans else (args.expected_locus_length + 2 * args.flank_bp)
    max_u_span = max(ref_u_spans) if ref_u_spans else args.flank_bp
    max_d_span = max(ref_d_spans) if ref_d_spans else args.flank_bp

    # Load GTDB protein KO annotations for reference flanking genes AND target genes
    ref_flank_annotations = {}
    if args.ko_index_cache or (args.gtdb_ko_table and os.path.exists(args.gtdb_ko_table)):
        target_clean_acc = args.target_name.replace("GB_", "").replace("RS_", "")
        combined_accs = set(all_ref_accs) | {args.target_name, target_clean_acc, chosen_contig.id}
        if "." in chosen_contig.id:
            combined_accs.add(chosen_contig.id.split(".")[0])
        if "." in target_clean_acc:
            combined_accs.add(target_clean_acc.split(".")[0])

        contig_pids = set()
        for p, inf in gff_coords.items():
            c_ctg = inf.get("contig", "")
            c_full = inf.get("full_contig", "")
            if c_ctg == chosen_contig.id or c_full == chosen_contig.id or chosen_contig.id.endswith(f"___{c_ctg}"):
                contig_pids.add(p)

        combined_pids = set(all_ref_pids) | contig_pids

        # 1. Try pre-filtered target KO micro-index cache first (instant for target operon genes)
        cache_file = args.ko_index_cache
        if not cache_file and args.out_topology_tsv:
            potential_cache = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.out_topology_tsv)), "..", "..", "cache", "target_kos_gtdb_index.tsv"))
            if os.path.exists(potential_cache) and os.path.getsize(potential_cache) > 0:
                cache_file = potential_cache

        def merge_genome_isolated_annos(src_annos, dst_target_dict, dst_ref_dict, t_clean):
            if not src_annos:
                return

            def update_dict_better(dst, src):
                for k, v in src.items():
                    if k not in dst or is_better_hit(v, dst[k]):
                        dst[k] = v

            # 1. Update target annotations ONLY with target's own annotations
            t_map = src_annos.get(t_clean) or (src_annos.get("__by_genome__", {}).get(t_clean) if isinstance(src_annos, dict) else None)
            if t_map and isinstance(t_map, dict):
                update_dict_better(dst_target_dict, t_map)
            # Also if any full unique ID matches target genome accession prefix, include it
            for k, v in src_annos.items():
                if k.startswith("__"):
                    continue
                if isinstance(v, dict) and "ko" in v:
                    if t_clean in k:
                        if k not in dst_target_dict or is_better_hit(v, dst_target_dict[k]):
                            dst_target_dict[k] = v

            # 2. Update reference annotations by genome
            if "__by_genome__" in src_annos:
                if "__by_genome__" not in dst_ref_dict:
                    dst_ref_dict["__by_genome__"] = {}
                for g_k, g_v in src_annos["__by_genome__"].items():
                    if g_k != t_clean:  # do not pollute ref dict with target
                        sub_dict = dst_ref_dict["__by_genome__"].setdefault(g_k, {})
                        update_dict_better(sub_dict, g_v)
                        top_dict = dst_ref_dict.setdefault(g_k, {})
                        update_dict_better(top_dict, g_v)
            for k, v in src_annos.items():
                if k.startswith("__"):
                    continue
                if isinstance(v, dict) and "ko" in v:
                    # Global full IDs
                    if k not in dst_ref_dict or is_better_hit(v, dst_ref_dict[k]):
                        dst_ref_dict[k] = v

        if cache_file and os.path.exists(cache_file):
            cache_annos = load_ref_flanking_annotations_from_table(cache_file, combined_accs, combined_pids)
            if cache_annos:
                merge_genome_isolated_annos(cache_annos, target_anno_dict, ref_flank_annotations, target_clean_acc)

        # 2. Check dedicated flanking KO cache before falling back to full GTDB KO table
        flank_cache_file = ""
        if args.ko_index_cache:
            flank_cache_file = os.path.join(os.path.dirname(os.path.abspath(args.ko_index_cache)), "flanking_ko_cache.tsv")
        elif args.out_topology_tsv:
            flank_cache_file = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.out_topology_tsv)), "..", "..", "cache", "flanking_ko_cache.tsv"))

        if flank_cache_file and os.path.exists(flank_cache_file) and os.path.getsize(flank_cache_file) > 10:
            flank_annos = load_ref_flanking_annotations_from_table(flank_cache_file, combined_accs, combined_pids)
            if flank_annos:
                merge_genome_isolated_annos(flank_annos, target_anno_dict, ref_flank_annotations, target_clean_acc)
                print(f"[{args.target_name}] Fast-loaded {len(flank_annos)} flanking annotations from cache: {flank_cache_file}", flush=True)

        # 3. Reference flanking genes (and any unannotated target genes) require querying full GTDB KO table
        # Check which proteins in reference loci or target contig are missing annotation
        missing_ref_pids = set(all_ref_pids) - set(ref_flank_annotations.keys())
        missing_target_pids = set(contig_pids) - set(target_anno_dict.keys())
        missing_pids = missing_ref_pids | missing_target_pids

        if missing_pids and args.gtdb_ko_table and os.path.exists(args.gtdb_ko_table):
            print(f"[{args.target_name}] Querying full GTDB KO table for {len(missing_ref_pids)} ref proteins and {len(missing_target_pids)} target proteins across {len(combined_accs)} genomes...", flush=True)
            all_table_annotations = load_ref_flanking_annotations_from_table(args.gtdb_ko_table, combined_accs, missing_pids)
            if all_table_annotations:
                merge_genome_isolated_annos(all_table_annotations, target_anno_dict, ref_flank_annotations, target_clean_acc)
                print(f"[{args.target_name}] Loaded {len(all_table_annotations)} flanking annotations from full GTDB KO table (isolated per-genome).", flush=True)

                # Persist newly discovered flanking annotations to flank_cache_file for subsequent tasks
                if flank_cache_file:
                    try:
                        os.makedirs(os.path.dirname(os.path.abspath(flank_cache_file)), exist_ok=True)
                        write_header = not os.path.exists(flank_cache_file) or os.path.getsize(flank_cache_file) == 0
                        with open(flank_cache_file, "a", encoding="utf-8") as f_flank:
                            if write_header:
                                f_flank.write("protein_id\tko\tgenome_id\tko_definition\tbitscore\tevalue\tko_rank\n")
                            for pid, info in all_table_annotations.items():
                                if pid.startswith("__"): continue
                                if isinstance(info, dict) and "ko" in info:
                                    gid = info.get("genome_id", "")
                                    ko = info.get("ko", "")
                                    k_def = info.get("def", "")
                                    sc = info.get("score", "")
                                    ev = info.get("evalue", "")
                                    rk = info.get("rank", "")
                                    f_flank.write(f"{pid}\t{ko}\t{gid}\t{k_def}\t{sc}\t{ev}\t{rk}\n")
                        print(f"[{args.target_name}] Appended {len(all_table_annotations)} annotations to flanking cache: {flank_cache_file}", flush=True)
                    except Exception as ce:
                        print(f"Notice: Failed to update flanking cache {flank_cache_file}: {ce}", file=sys.stderr)

    # Identify whether target operon genes exist on chosen_contig
    chosen_target_op_cds = []
    for p, inf in gff_coords.items():
        if inf.get("contig") == chosen_contig.id or is_same_contig(inf.get("contig"), chosen_contig.id):
            clean_p = p.split("___")[-1] if "___" in p else p
            if (p in target_op_pids or clean_p in target_op_pids or
                p in target_pid_to_ko or clean_p in target_pid_to_ko):
                chosen_target_op_cds.append(inf)
                continue
            t_anno = lookup_target_anno(p, target_anno_dict)
            ko_val = t_anno.get("ko") if t_anno else (inf.get("ko") or target_pid_to_ko.get(p) or target_pid_to_ko.get(clean_p))
            if ko_val and (not target_ko_set or ko_val in target_ko_set):
                chosen_target_op_cds.append(inf)

    active_hits = locus_hits if locus_hits else best_cluster
    flank_expansion = max(args.flank_bp, 20000)

    if chosen_target_op_cds:
        # Case A: Operon is present in the target on chosen_contig.
        # Focus slice tightly around the operon and immediate proximal flanking synteny.
        op_min = min(c["start"] for c in chosen_target_op_cds)
        op_max = max(c["end"] for c in chosen_target_op_cds)
        op_mid = (op_min + op_max) // 2

        # Flank reach: only consider syntenic hits within proximity to avoid stretching over huge unrelated gene insertions (>45 kb)
        max_flank_reach = max(int(args.flank_bp * 1.5), 45000)
        prox_hits = [h for h in active_hits if (h.get("end", 0) >= op_min - max_flank_reach) and (h.get("start", 0) <= op_max + max_flank_reach)]

        if prox_hits:
            p_min = min(h["start"] for h in prox_hits)
            p_max = max(h["end"] for h in prox_hits)
            raw_start = min(op_min - flank_expansion, p_min - flank_expansion // 2)
            raw_end = max(op_max + flank_expansion, p_max + flank_expansion // 2)
        else:
            raw_start = op_min - flank_expansion
            raw_end = op_max + flank_expansion

        u_hits = [h for h in active_hits if h["side"] == "UP"]
        d_hits = [h for h in active_hits if h["side"] == "DOWN"]
        if u_hits and d_hits:
            u_max = max(h["end"] for h in u_hits)
            d_min = min(h["start"] for h in d_hits)
            if u_max <= d_min:
                observed_anchor_gap = max(0, d_min - u_max)
            else:
                d_max = max(h["end"] for h in d_hits)
                u_min = min(h["start"] for h in u_hits)
                observed_anchor_gap = max(0, u_min - d_max)
        else:
            observed_anchor_gap = 0

    elif active_hits:
        # Case B: Operon is absent (Gene Loss / Deletion candidate).
        u_hits = [h for h in active_hits if h["side"] == "UP"]
        d_hits = [h for h in active_hits if h["side"] == "DOWN"]
        if u_hits and d_hits:
            med_u = float(pd.Series([(h["start"] + h["end"]) / 2.0 for h in u_hits]).median())
            med_d = float(pd.Series([(h["start"] + h["end"]) / 2.0 for h in d_hits]).median())
            if med_u <= med_d:
                u_max = max(h["end"] for h in u_hits)
                d_min = min(h["start"] for h in d_hits)
                observed_anchor_gap = max(0, d_min - u_max)
                junction_mid = (u_max + d_min) // 2
            else:
                d_max = max(h["end"] for h in d_hits)
                u_min = min(h["start"] for h in u_hits)
                observed_anchor_gap = max(0, u_min - d_max)
                junction_mid = (d_max + u_min) // 2
        elif u_hits:
            u_max = max(h["end"] for h in u_hits)
            observed_anchor_gap = 0
            junction_mid = u_max
        elif d_hits:
            d_min = min(h["start"] for h in d_hits)
            observed_anchor_gap = 0
            junction_mid = d_min
        else:
            observed_anchor_gap = 0
            junction_mid = (c_min + c_max) // 2

        u_ext = max(max_u_span, flank_expansion)
        d_ext = max(max_d_span, flank_expansion)

        if observed_anchor_gap <= max(args.flank_bp * 2, 60000):
            raw_start = min(junction_mid - u_ext, c_min - flank_expansion)
            raw_end = max(junction_mid + d_ext, c_max + flank_expansion)
        else:
            # Massive insertion between anchors: center compactly around breakpoint
            max_del_span = max(int(max_ref_span * 2.0), 90000)
            raw_start = junction_mid - max_del_span // 2
            raw_end = junction_mid + max_del_span // 2

    else:
        # Case C: No hits and no operon
        mid = chosen_contig_len // 2
        half_span = max(max_ref_span // 2, 25000)
        raw_start = mid - half_span
        raw_end = mid + half_span
        observed_anchor_gap = -1

    slice_start = max(0, raw_start)
    slice_end = min(chosen_contig_len, raw_end)

    # Guarantee target slice is at least as long as the reference locus (avoiding truncated target issue)
    curr_len = slice_end - slice_start
    desired_len = max(max_ref_span, curr_len)
    if curr_len < desired_len and chosen_contig_len > curr_len:
        deficit = min(desired_len - curr_len, chosen_contig_len - curr_len)
        add_left = min(slice_start, deficit // 2)
        slice_start -= add_left
        remaining = deficit - add_left
        add_right = min(chosen_contig_len - slice_end, remaining)
        slice_end += add_right
        leftover = remaining - add_right
        if leftover > 0 and slice_start > 0:
            slice_start = max(0, slice_start - leftover)

    # Bound maximum target slice length to prevent unreadable over-extended plots (>100 kb)
    max_allowed_slice = max(int(max_ref_span * 2.5), 100000)
    if (slice_end - slice_start) > max_allowed_slice:
        mid_pt = op_mid if chosen_target_op_cds else ((slice_start + slice_end) // 2)
        half_lim = max_allowed_slice // 2
        slice_start = max(0, mid_pt - half_lim)
        slice_end = min(chosen_contig_len, slice_start + max_allowed_slice)

    dist_left_edge = slice_start
    dist_right_edge = chosen_contig_len - slice_end

    # Index all CDS features by contig from gff_coords for accurate per-contig evaluation
    cds_by_contig = {}
    seen_all_cds_locs = set()
    for pid, c_info in sorted(gff_coords.items(), key=lambda kv: (kv[1].get("start", 0), kv[1].get("end", 0))):
        c_ctg = c_info.get("contig", "")
        c_full = c_info.get("full_contig", "")
        g_start = int(c_info.get("start", 0))
        g_end = int(c_info.get("end", 0))
        strand_val = 1 if c_info.get("strand", "+") == "+" else -1
        loc_key = (c_ctg, g_start, g_end, strand_val)
        if loc_key in seen_all_cds_locs:
            continue
        seen_all_cds_locs.add(loc_key)
        item = {
            "pid": pid,
            "start": g_start,
            "end": g_end,
            "strand": strand_val
        }
        cds_by_contig.setdefault(c_ctg, []).append(item)
        if c_full and c_full != c_ctg:
            cds_by_contig.setdefault(c_full, []).append(item)

    def get_cds_for_contig(ctg_name):
        if not ctg_name:
            return []
        if ctg_name in cds_by_contig:
            return cds_by_contig[ctg_name]
        c_strip = ctg_name.split("___")[-1].strip("_") if "___" in ctg_name else ctg_name.strip("_")
        for k, cds_list in cds_by_contig.items():
            k_strip = k.split("___")[-1].strip("_") if "___" in k else k.strip("_")
            if k == ctg_name or k.endswith(ctg_name) or ctg_name.endswith(k) or k_strip == c_strip:
                return cds_list
        return []

    chosen_contig_cds = get_cds_for_contig(chosen_contig.id)

    # Write topology TSV with a summary row for EACH visualized reference
    topology_rows = []
    for r_acc in visualized_refs:
        ev = ref_evaluations.get(r_acc, {})
        r_cl = ev.get("best_cluster", [])
        u_h = [h for h in r_cl if h["side"] == "UP"]
        d_h = [h for h in r_cl if h["side"] == "DOWN"]

        r_l_ident = float(pd.Series([h["identity"] for h in u_h]).mean()) if u_h else 0.0
        r_r_ident = float(pd.Series([h["identity"] for h in d_h]).mean()) if d_h else 0.0
        r_l_found = len(u_h) > 0
        r_r_found = len(d_h) > 0

        best_c = ev.get("best_contig", "")
        same_contig = "YES" if (best_c and (is_same_contig(best_c, chosen_contig.id) or is_same_contig(best_c, chosen_contig_id))) else ("YES" if not best_c else "NO")

        # Determine CDS on the specific contig where this reference matched
        eval_ctg = chosen_contig.id if same_contig == "YES" else (best_c if best_c else chosen_contig.id)
        eval_cds = get_cds_for_contig(eval_ctg)

        # Detailed anchor matched genes and other genes analysis
        u_target_pids = sorted(list(set(h["t_pid"] for h in u_h)))
        d_target_pids = sorted(list(set(h["t_pid"] for h in d_h)))
        u_hit_locs = set((h["start"], h["end"]) for h in u_h)
        d_hit_locs = set((h["start"], h["end"]) for h in d_h)
        all_matched_pids = set(u_target_pids).union(set(d_target_pids))
        all_hit_locs = u_hit_locs.union(d_hit_locs)

        left_matched_cnt = len(u_target_pids)
        right_matched_cnt = len(d_target_pids)
        left_matched_str = ",".join(u_target_pids)
        right_matched_str = ",".join(d_target_pids)

        # Determine anchor spatial arrangement and robust boundary partitioning
        u_ctgs = set(h.get("contig", "") for h in u_h)
        d_ctgs = set(h.get("contig", "") for h in d_h)
        is_same_ctg_anchors = bool(u_h and d_h and u_ctgs == d_ctgs and len(u_ctgs) == 1)

        anchor_order = "YES"
        obs_gap = 0
        intervening_genes_list = []
        intervening_genes_cnt = 0
        intervening_genes_str = ""

        if u_h and d_h:
            if is_same_ctg_anchors:
                # Calculate median spatial position of UP hits vs DOWN hits on the target contig
                med_u = float(pd.Series([(h["start"] + h["end"]) / 2.0 for h in u_h]).median())
                med_d = float(pd.Series([(h["start"] + h["end"]) / 2.0 for h in d_h]).median())

                if med_u <= med_d:
                    anchor_order = "YES"  # Collinear: UP is on left, DOWN is on right
                    split_mid = (med_u + med_d) / 2.0
                    core_u = [h for h in u_h if (h["start"] + h["end"]) / 2.0 <= split_mid] or u_h
                    core_d = [h for h in d_h if (h["start"] + h["end"]) / 2.0 >= split_mid] or d_h
                    u_reg_s = min(h["start"] for h in core_u)
                    u_reg_e = max(h["end"] for h in core_u)
                    d_reg_s = min(h["start"] for h in core_d)
                    d_reg_e = max(h["end"] for h in core_d)
                    gap_s, gap_e = (u_reg_e, d_reg_s) if d_reg_s >= u_reg_e else (None, None)
                    obs_gap = max(0, d_reg_s - u_reg_e)
                else:
                    anchor_order = "INVERTED"  # Inverted: DOWN is on left, UP is on right
                    split_mid = (med_u + med_d) / 2.0
                    core_d = [h for h in d_h if (h["start"] + h["end"]) / 2.0 <= split_mid] or d_h
                    core_u = [h for h in u_h if (h["start"] + h["end"]) / 2.0 >= split_mid] or u_h
                    d_reg_s = min(h["start"] for h in core_d)
                    d_reg_e = max(h["end"] for h in core_d)
                    u_reg_s = min(h["start"] for h in core_u)
                    u_reg_e = max(h["end"] for h in core_u)
                    gap_s, gap_e = (d_reg_e, u_reg_s) if u_reg_s >= d_reg_e else (None, None)
                    obs_gap = max(0, u_reg_s - d_reg_e)

                if gap_s is not None and gap_e > gap_s:
                    intervening_genes_list = [c["pid"] for c in eval_cds if c["end"] > gap_s and c["start"] < gap_e and c["pid"] not in all_matched_pids and (c["start"], c["end"]) not in all_hit_locs]
                intervening_genes_cnt = len(intervening_genes_list)
                intervening_genes_str = ",".join(intervening_genes_list)
            else:
                anchor_order = "CROSS_CONTIG"
                obs_gap = 0
                intervening_genes_cnt = 0
                intervening_genes_str = "CROSS_CONTIG_SPLIT"
                u_reg_s = min(h["start"] for h in u_h) if u_h else 0
                u_reg_e = max(h["end"] for h in u_h) if u_h else 0
                d_reg_s = min(h["start"] for h in d_h) if d_h else 0
                d_reg_e = max(h["end"] for h in d_h) if d_h else 0
        elif u_h:
            anchor_order = "UNILATERAL_LEFT"
            obs_gap = 0
            u_reg_s = min(h["start"] for h in u_h)
            u_reg_e = max(h["end"] for h in u_h)
            d_reg_s, d_reg_e = 0, 0
        elif d_h:
            anchor_order = "UNILATERAL_RIGHT"
            obs_gap = 0
            d_reg_s = min(h["start"] for h in d_h)
            d_reg_e = max(h["end"] for h in d_h)
            u_reg_s, u_reg_e = 0, 0
        else:
            anchor_order = "NONE"
            obs_gap = 0
            u_reg_s, u_reg_e, d_reg_s, d_reg_e = 0, 0, 0, 0

        # Other genes inside left anchor region [u_reg_s, u_reg_e] strictly on eval_ctg
        if u_h and u_reg_e > u_reg_s:
            left_other = [c["pid"] for c in eval_cds if c["end"] >= u_reg_s and c["start"] <= u_reg_e and c["pid"] not in u_target_pids and (c["start"], c["end"]) not in u_hit_locs]
            left_other_cnt = len(left_other)
        else:
            left_other_cnt = 0

        # Other genes inside right anchor region [d_reg_s, d_reg_e] strictly on eval_ctg
        if d_h and d_reg_e > d_reg_s:
            right_other = [c["pid"] for c in eval_cds if c["end"] >= d_reg_s and c["start"] <= d_reg_e and c["pid"] not in d_target_pids and (c["start"], c["end"]) not in d_hit_locs]
            right_other_cnt = len(right_other)
        else:
            right_other_cnt = 0

        # Locus slice coordinates and gene metrics for eval_ctg
        if same_contig == "YES":
            r_slice_s = slice_start
            r_slice_e = slice_end
            r_dist_l = dist_left_edge
            r_dist_r = dist_right_edge
            r_ctg_len = chosen_contig_len
            r_ctg_name = chosen_contig.id
        else:
            eval_rec = find_target_record(eval_ctg, target_records)
            r_ctg_len = len(eval_rec.seq) if eval_rec else chosen_contig_len
            r_ctg_name = eval_rec.id if eval_rec else eval_ctg
            if r_cl:
                c_min_r = min(h["start"] for h in r_cl)
                c_max_r = max(h["end"] for h in r_cl)
            elif eval_cds:
                c_min_r = min(c["start"] for c in eval_cds)
                c_max_r = max(c["end"] for c in eval_cds)
            else:
                c_min_r, c_max_r = 0, r_ctg_len
            flank_exp = max(args.flank_bp, 20000)
            r_slice_s = max(0, c_min_r - flank_exp)
            r_slice_e = min(r_ctg_len, c_max_r + flank_exp)
            curr_span = r_slice_e - r_slice_s
            desired_span = max(max_ref_span, curr_span)
            if curr_span < desired_span and r_ctg_len > curr_span:
                deficit = min(desired_span - curr_span, r_ctg_len - curr_span)
                add_l = min(r_slice_s, deficit // 2)
                r_slice_s -= add_l
                rem = deficit - add_l
                add_r = min(r_ctg_len - r_slice_e, rem)
                r_slice_e += add_r
                leftover = rem - add_r
                if leftover > 0 and r_slice_s > 0:
                    r_slice_s = max(0, r_slice_s - leftover)
            max_allowed_alt = max(int(max_ref_span * 2.5), 100000)
            if (r_slice_e - r_slice_s) > max_allowed_alt:
                mid_alt = (r_slice_s + r_slice_e) // 2
                half_alt = max_allowed_alt // 2
                r_slice_s = max(0, mid_alt - half_alt)
                r_slice_e = min(r_ctg_len, r_slice_s + max_allowed_alt)
            r_dist_l = r_slice_s
            r_dist_r = r_ctg_len - r_slice_e

        # Genes across the evaluated target slice on eval_ctg
        slice_cds = [c["pid"] for c in eval_cds if c["end"] > r_slice_s and c["start"] < r_slice_e]
        target_locus_tot_genes = len(slice_cds)
        target_locus_oth_genes = len([c["pid"] for c in eval_cds if c["end"] > r_slice_s and c["start"] < r_slice_e and c["pid"] not in all_matched_pids and (c["start"], c["end"]) not in all_hit_locs])

        is_edge_safe = (r_dist_l >= required_edge_margin and r_dist_r >= required_edge_margin)
        if chosen_geometry == "BILATERAL_SAME_CONTIG":
            gate1_status = "PASS" if is_edge_safe else "ASSEMBLY_UNRESOLVED"
        elif chosen_geometry in ["UNILATERAL_SAME_CONTIG", "SPLIT_CONTIG"]:
            gate1_status = "ASSEMBLY_UNRESOLVED"
        else:
            gate1_status = "INSUFFICIENT_EVIDENCE"

        chosen_left_anchor, chosen_right_anchor = contig_boundary_anchors.get(chosen_contig_id, ("NONE", "NONE"))
        chosen_interval_start, chosen_interval_end = contig_ancestral_intervals.get(chosen_contig_id, (0, chosen_contig_len))
        chosen_obs_span = contig_observed_spans.get(chosen_contig_id, 0)
        chosen_loc_score = contig_loc_score.get(chosen_contig_id, 0.0)
        chosen_bilat_score = contig_bilat_score.get(chosen_contig_id, 0.0)

        topology_rows.append({
            "target_id": args.target_name,
            "module_id": args.module_id,
            "module_name": args.module_name or args.module_id,
            "nearest_positive_ref": r_acc,
            "left_anchor_1": r_l_found,
            "left_anchor_2": r_l_found,
            "right_anchor_1": r_r_found,
            "right_anchor_2": r_r_found,
            "left_anchors_matched_count": left_matched_cnt,
            "right_anchors_matched_count": right_matched_cnt,
            "left_matched_genes": left_matched_str,
            "right_matched_genes": right_matched_str,
            "left_anchor_region_other_genes": left_other_cnt,
            "right_anchor_region_other_genes": right_other_cnt,
            "intervening_genes_count": intervening_genes_cnt,
            "intervening_genes": intervening_genes_str,
            "target_locus_total_genes": target_locus_tot_genes,
            "target_locus_other_genes": target_locus_oth_genes,
            "anchor_left_identity": round(r_l_ident, 2),
            "anchor_right_identity": round(r_r_ident, 2),
            "anchor_order_consistent": anchor_order,
            "anchor_orientation_consistent": "YES",
            "same_contig": same_contig,
            "distance_to_left_edge": r_dist_l,
            "distance_to_right_edge": r_dist_r,
            "required_edge_margin": required_edge_margin,
            "expected_locus_length": args.expected_locus_length,
            "observed_anchor_gap": obs_gap,
            "target_contig_length": r_ctg_len,
            "target_slice_start": r_slice_s,
            "target_slice_end": r_slice_e,
            "target_contig": r_ctg_name,
            "localization_status": localization_status,
            "syntenic_geometry": chosen_geometry,
            "locator_score": chosen_loc_score,
            "bilateral_context_score": chosen_bilat_score,
            "left_boundary_anchor": chosen_left_anchor,
            "right_boundary_anchor": chosen_right_anchor,
            "ancestral_interval_start": chosen_interval_start,
            "ancestral_interval_end": chosen_interval_end,
            "observed_interanchor_span": chosen_obs_span,
            "expected_interanchor_span": expected_interanchor_span,
            "n_effective_refs": n_eff,
            "k_effective_refs": contig_supporting_clusters.get(chosen_contig_id, 0),
            "localization_p_value": round(p_emp, 6),
            "statistical_support": stat_support,
            "empirical_s_obs": round(s_obs, 2),
            "empirical_s_null_mean": round(s_null_mean, 2),
            "empirical_s_null_95": round(s_null_95, 2),
            "dominance_ratio": round(dominance_ratio, 2),
            "gate1_status": gate1_status
        })

    if not topology_rows:
        fallback_edge_safe = (dist_left_edge >= required_edge_margin and dist_right_edge >= required_edge_margin)
        if chosen_geometry == "BILATERAL_SAME_CONTIG":
            fallback_gate1 = "PASS" if fallback_edge_safe else "ASSEMBLY_UNRESOLVED"
        elif chosen_geometry in ["UNILATERAL_SAME_CONTIG", "SPLIT_CONTIG"]:
            fallback_gate1 = "ASSEMBLY_UNRESOLVED"
        else:
            fallback_gate1 = "INSUFFICIENT_EVIDENCE"

        chosen_left_anchor, chosen_right_anchor = contig_boundary_anchors.get(chosen_contig_id, ("NONE", "NONE"))
        chosen_interval_start, chosen_interval_end = contig_ancestral_intervals.get(chosen_contig_id, (0, chosen_contig_len))
        chosen_obs_span = contig_observed_spans.get(chosen_contig_id, 0)
        chosen_loc_score = contig_loc_score.get(chosen_contig_id, 0.0)
        chosen_bilat_score = contig_bilat_score.get(chosen_contig_id, 0.0)

        topology_rows.append({
            "target_id": args.target_name,
            "module_id": args.module_id,
            "module_name": args.module_name or args.module_id,
            "nearest_positive_ref": best_ref_acc,
            "left_anchor_1": False,
            "left_anchor_2": False,
            "right_anchor_1": False,
            "right_anchor_2": False,
            "left_anchors_matched_count": 0,
            "right_anchors_matched_count": 0,
            "left_matched_genes": "",
            "right_matched_genes": "",
            "left_anchor_region_other_genes": 0,
            "right_anchor_region_other_genes": 0,
            "intervening_genes_count": 0,
            "intervening_genes": "",
            "target_locus_total_genes": 0,
            "target_locus_other_genes": 0,
            "anchor_left_identity": 0.0,
            "anchor_right_identity": 0.0,
            "anchor_order_consistent": "YES",
            "anchor_orientation_consistent": "YES",
            "same_contig": "YES",
            "distance_to_left_edge": dist_left_edge,
            "distance_to_right_edge": dist_right_edge,
            "required_edge_margin": required_edge_margin,
            "expected_locus_length": args.expected_locus_length,
            "observed_anchor_gap": observed_anchor_gap,
            "target_contig_length": chosen_contig_len,
            "target_slice_start": slice_start,
            "target_slice_end": slice_end,
            "target_contig": chosen_contig.id,
            "localization_status": localization_status,
            "syntenic_geometry": chosen_geometry,
            "locator_score": chosen_loc_score,
            "bilateral_context_score": chosen_bilat_score,
            "left_boundary_anchor": chosen_left_anchor,
            "right_boundary_anchor": chosen_right_anchor,
            "ancestral_interval_start": chosen_interval_start,
            "ancestral_interval_end": chosen_interval_end,
            "observed_interanchor_span": chosen_obs_span,
            "expected_interanchor_span": expected_interanchor_span,
            "n_effective_refs": n_eff,
            "k_effective_refs": contig_supporting_clusters.get(chosen_contig_id, 0),
            "localization_p_value": round(p_emp, 6),
            "statistical_support": stat_support,
            "empirical_s_obs": round(s_obs, 2),
            "empirical_s_null_mean": round(s_null_mean, 2),
            "empirical_s_null_95": round(s_null_95, 2),
            "dominance_ratio": round(dominance_ratio, 2),
            "gate1_status": fallback_gate1
        })

    pd.DataFrame(topology_rows).to_csv(args.out_topology_tsv, sep="\t", index=False)
    print(f"Gate 1 Topology evaluated for {len(topology_rows)} references: target slice=[{slice_start}..{slice_end}] ({slice_end - slice_start} bp, max_ref_span={max_ref_span} bp)", flush=True)

    # 8. Extract Target Locus Slice and GenBank File
    target_faa_dict = {r.id: str(r.seq) for r in load_fasta_records(args.target_faa)}
    locus_seq = chosen_contig.seq[slice_start:slice_end]
    target_name_clean = re.sub(r'[^a-zA-Z0-9_]', '_', args.target_name)[:16]
    target_subrec = SeqRecord(
        locus_seq,
        id=target_name_clean,
        name=target_name_clean,
        description=f"{args.target_name} [{args.module_id}] ({chosen_contig.id}:{slice_start}-{slice_end})"
    )
    target_subrec.annotations["molecule_type"] = "DNA"

    target_features = []
    seen_feature_locs = set()
    sorted_cds = sorted(
        gff_coords.items(),
        key=lambda kv: (kv[1].get("contig", ""), kv[1].get("start", 0), kv[1].get("end", 0))
    )

    for pid, c_info in sorted_cds:
        c_ctg = c_info.get("contig", "")
        c_full = c_info.get("full_contig", "")
        matched_ctg = (
            c_ctg == chosen_contig.id or
            c_full == chosen_contig.id or
            chosen_contig.id.endswith(f"___{c_ctg}") or
            chosen_contig.id.endswith(f"_{c_ctg}") or
            (c_ctg and chosen_contig.id.startswith(c_ctg))
        )
        if not matched_ctg:
            continue

        g_start = c_info["start"]
        g_end = c_info["end"]
        strand_val = 1 if c_info.get("strand", "+") == "+" else -1
        loc_key = (g_start, g_end, strand_val)
        if loc_key in seen_feature_locs:
            continue

        if g_end > slice_start and g_start < slice_end:
            seen_feature_locs.add(loc_key)
            rel_start = max(0, min(len(locus_seq) - 1, g_start - slice_start))
            rel_end = max(rel_start + 1, min(len(locus_seq), g_end - slice_start))

            trans_seq = target_faa_dict.get(pid, "")
            if not trans_seq:
                for alt_id in [pid.replace("gene", "prot"), pid.replace("prot", "gene"), f"{chosen_contig.id}___{pid}"]:
                    if alt_id in target_faa_dict:
                        trans_seq = target_faa_dict[alt_id]
                        break
            if not trans_seq and (rel_end - rel_start) >= 30:
                try:
                    sub_dna = locus_seq[rel_start:rel_end] if strand_val == 1 else locus_seq[rel_start:rel_end].reverse_complement()
                    trans_seq = str(sub_dna.translate(to_stop=True))
                except Exception:
                    trans_seq = ""

            if trans_seq:
                cand_pids = [
                    pid,
                    pid.replace("gene", "prot"),
                    pid.replace("prot", "gene"),
                    c_info.get("protein_id", ""),
                    c_info.get("gene_id", "")
                ]
                t_hit = None
                for c_pid in cand_pids:
                    if c_pid:
                        t_hit = lookup_target_anno(c_pid, target_anno_dict)
                        if t_hit:
                            break

                # Homology fallback from active_hits / reference flanking hits
                if not t_hit and active_hits:
                    for h in active_hits:
                        h_tid = h.get("target_id", "")
                        if any(c_pid and (c_pid == h_tid or c_pid in h_tid or h_tid in c_pid) for c_pid in cand_pids):
                            ref_qid = h.get("query_id", "")
                            ref_hit = lookup_ref_flanking_anno(ref_qid, ref_flank_annotations, genome_acc=h.get("ref_acc")) if ref_flank_annotations else None
                            if ref_hit and ref_hit.get("ko"):
                                t_hit = {
                                    "ko": ref_hit["ko"],
                                    "definition": f"{ref_hit.get('definition', '')} (homology {h['identity']:.0f}% to {h['ref_acc']})",
                                    "score": str(h.get("score", "")),
                                    "threshold": str(ref_hit.get("threshold", "")),
                                    "evalue": str(h.get("eval", ""))
                                }
                                break

                quals = format_target_gene_qualifiers(pid, t_hit, trans_seq, ko_definitions, c_info=c_info)
                feat = SeqFeature.SeqFeature(
                    SeqFeature.FeatureLocation(rel_start, rel_end, strand=strand_val),
                    type="CDS",
                    qualifiers=quals
                )
                target_features.append(feat)

    target_features.sort(key=lambda x: int(x.location.start))
    target_subrec.features = target_features

    with open(args.out_target_locus_fa, "w") as f_fa:
        SeqIO.write(target_subrec, f_fa, "fasta")

    with open(args.out_target_locus_gbk, "w") as f_gbk:
        SeqIO.write(target_subrec, f_gbk, "genbank")
    fix_gbk_wrapped_qualifiers(args.out_target_locus_gbk)

    print(f"Stage 1 Target Locus Slice: {args.out_target_locus_fa} ({len(locus_seq)} bp, {len(target_features)} CDS)", flush=True)

    # 9. Extract Reference Loci for Visualized References
    # Clean previous gbk files in out_refs_locus_dir to prevent obsolete files from appearing
    for fname in os.listdir(args.out_refs_locus_dir):
        if fname.endswith(".gbk") or fname.endswith(".fa"):
            try:
                os.remove(os.path.join(args.out_refs_locus_dir, fname))
            except Exception:
                pass

    for r_acc in visualized_refs:
        cand = ref_candidates.get(r_acc)
        if not cand:
            continue

        genes_in_order = cand["u_prots"] + cand["core_prots"] + cand["d_prots"]
        if not genes_in_order:
            continue

        ref_coords = cand["ref_coords"]
        ref_genome_fa = cand["ref_genome_fa"]
        ref_ctg_id = cand["best_ctg"]
        has_coords = all(p.id in ref_coords for p in genes_in_order)

        ref_subrec = None

        if has_coords and ref_genome_fa and os.path.exists(ref_genome_fa):
            try:
                ref_genome_records = {rec.id: rec for rec in load_fasta_records(ref_genome_fa)}
                full_ref_ctg = find_matching_record(ref_ctg_id, ref_genome_records)
                if full_ref_ctg:
                    r_min = min(ref_coords[p.id]["start"] for p in genes_in_order)
                    r_max = max(ref_coords[p.id]["end"] for p in genes_in_order)
                    pad_ref = 3000
                    r_slice_s = max(0, r_min - pad_ref)
                    r_slice_e = min(len(full_ref_ctg.seq), r_max + pad_ref)

                    ref_slice_seq = full_ref_ctg.seq[r_slice_s:r_slice_e]
                    ref_name_clean = re.sub(r'[^a-zA-Z0-9_]', '_', r_acc)[:16]
                    ref_subrec = SeqRecord(
                        ref_slice_seq,
                        id=ref_name_clean,
                        name=ref_name_clean,
                        description=f"{r_acc} ({ref_ctg_id}:{r_slice_s}-{r_slice_e})"
                    )
                    ref_subrec.annotations["molecule_type"] = "DNA"

                    r_feats = []
                    pid_to_ko_map = cand.get("pid_to_ko", {})
                    for p in genes_in_order:
                        c_inf = ref_coords[p.id]
                        rel_s = max(0, c_inf["start"] - r_slice_s)
                        rel_e = min(len(ref_slice_seq), c_inf["end"] - r_slice_s)
                        st_val = 1 if c_inf["strand"] == "+" else -1
                        is_op = (p.id in cand["op_pids"]) or any(op in p.id or p.id in op for op in cand["op_pids"])
                        ko = find_ko_for_pid(p.id, pid_to_ko_map)
                        def_val = ""
                        fl_hit = None
                        if ref_flank_annotations:
                            clean_r_acc = r_acc.replace("GB_", "").replace("RS_", "")
                            fl_hit = lookup_ref_flanking_anno(p.id, ref_flank_annotations, genome_acc=clean_r_acc)
                            if fl_hit:
                                if not ko:
                                    ko = fl_hit.get("ko", "")
                                if not def_val:
                                    def_val = fl_hit.get("definition", "")
                        quals = format_gene_qualifiers(p, is_op, ko, ko_definitions, default_def=def_val, c_inf=c_inf, hit=fl_hit)

                        feat = SeqFeature.SeqFeature(
                            SeqFeature.FeatureLocation(rel_s, rel_e, strand=st_val),
                            type="CDS",
                            qualifiers=quals
                        )
                        r_feats.append(feat)

                    r_feats.sort(key=lambda x: int(x.location.start))
                    ref_subrec.features = r_feats
            except Exception as e:
                print(f"Notice: genomic slice extraction failed for {r_acc} ({e}), falling back to stitched CDS.", file=sys.stderr)
                ref_subrec = None

        if ref_subrec is None:
            # Fallback: construct stitched record with real strands
            dna_parts = []
            r_feats = []
            pid_to_ko_map = cand.get("pid_to_ko", {})
            curr_pos = 0
            spacer = "A" * 120

            for p in genes_in_order:
                gene_dna = amino_acid_to_dna(p.seq)
                g_len = len(gene_dna)
                c_inf = ref_coords.get(p.id, {})
                st_val = 1 if c_inf.get("strand", "+") == "+" else -1
                is_op = (p.id in cand["op_pids"]) or any(op in p.id or p.id in op for op in cand["op_pids"])
                ko = find_ko_for_pid(p.id, pid_to_ko_map)
                def_val = ""
                fl_hit = None
                if ref_flank_annotations:
                    clean_r_acc = r_acc.replace("GB_", "").replace("RS_", "")
                    fl_hit = lookup_ref_flanking_anno(p.id, ref_flank_annotations, genome_acc=clean_r_acc)
                    if fl_hit:
                        if not ko:
                            ko = fl_hit.get("ko", "")
                        if not def_val:
                            def_val = fl_hit.get("definition", "")
                quals = format_gene_qualifiers(p, is_op, ko, ko_definitions, default_def=def_val, c_inf=c_inf, hit=fl_hit)

                feat = SeqFeature.SeqFeature(
                    SeqFeature.FeatureLocation(curr_pos, curr_pos + g_len, strand=st_val),
                    type="CDS",
                    qualifiers=quals
                )
                r_feats.append(feat)
                dna_parts.append(gene_dna)
                dna_parts.append(spacer)
                curr_pos += g_len + len(spacer)

            ref_full_dna = Seq("".join(dna_parts)) if dna_parts else Seq("A" * 1000)
            ref_name_clean = re.sub(r'[^a-zA-Z0-9_]', '_', r_acc)[:16]
            ref_subrec = SeqRecord(
                ref_full_dna,
                id=ref_name_clean,
                name=ref_name_clean,
                description=f"Positive reference locus from {r_acc}"
            )
            ref_subrec.annotations["molecule_type"] = "DNA"
            ref_subrec.features = r_feats

        # Write reference locus files
        ref_fa_out = os.path.join(args.out_refs_locus_dir, f"{r_acc}_locus.fa")
        with open(ref_fa_out, "w") as rf:
            SeqIO.write(ref_subrec, rf, "fasta")

        ref_gbk_out = os.path.join(args.out_refs_locus_dir, f"{r_acc}_locus.gbk")
        with open(ref_gbk_out, "w") as rg:
            SeqIO.write(ref_subrec, rg, "genbank")
        fix_gbk_wrapped_qualifiers(ref_gbk_out)

        print(f"Stage 1 Reference Locus Slice written: {ref_gbk_out} ({len(ref_subrec.features)} CDS)", flush=True)

        # Write dedicated pairwise target locus GBK for this reference
        # If the reference matched a different contig from chosen_contig, extract the locus on that specific contig!
        ev = ref_evaluations.get(r_acc, {})
        r_best_c = ev.get("best_contig", "")
        pair_target_gbk = os.path.join(args.out_refs_locus_dir, f"{args.target_name}_vs_{r_acc}_target_locus.gbk")
        is_alt = bool(r_best_c and not (is_same_contig(r_best_c, chosen_contig.id) or is_same_contig(r_best_c, chosen_contig_id)))
        if not is_alt:
            # Same contig as chosen_contig: clean up any old cross-contig target gbk if present
            if os.path.exists(pair_target_gbk):
                try:
                    os.remove(pair_target_gbk)
                except Exception:
                    pass
        else:
            alt_rec = find_target_record(r_best_c, target_records)
            if alt_rec:
                r_cl = ev.get("best_cluster", [])
                alt_cds = get_cds_for_contig(r_best_c)
                if r_cl:
                    c_min_r = min(h["start"] for h in r_cl)
                    c_max_r = max(h["end"] for h in r_cl)
                elif alt_cds:
                    c_min_r = min(c["start"] for c in alt_cds)
                    c_max_r = max(c["end"] for c in alt_cds)
                else:
                    c_min_r, c_max_r = 0, len(alt_rec.seq)

                flank_exp = max(args.flank_bp, 20000)
                alt_s = max(0, c_min_r - flank_exp)
                alt_e = min(len(alt_rec.seq), c_max_r + flank_exp)
                curr_alt_len = alt_e - alt_s
                desired_alt_len = max(max_ref_span, curr_alt_len)
                if curr_alt_len < desired_alt_len and len(alt_rec.seq) > curr_alt_len:
                    deficit = min(desired_alt_len - curr_alt_len, len(alt_rec.seq) - curr_alt_len)
                    add_l = min(alt_s, deficit // 2)
                    alt_s -= add_l
                    rem = deficit - add_l
                    add_r = min(len(alt_rec.seq) - alt_e, rem)
                    alt_e += add_r
                    leftover = rem - add_r
                    if leftover > 0 and alt_s > 0:
                        alt_s = max(0, alt_s - leftover)

                max_allowed_alt = max(int(max_ref_span * 2.5), 100000)
                if (alt_e - alt_s) > max_allowed_alt:
                    mid_alt = (alt_s + alt_e) // 2
                    half_alt = max_allowed_alt // 2
                    alt_s = max(0, mid_alt - half_alt)
                    alt_e = min(len(alt_rec.seq), alt_s + max_allowed_alt)

                alt_seq = alt_rec.seq[alt_s:alt_e]
                clean_id = re.sub(r'[^a-zA-Z0-9_]', '_', f"{args.target_name}_{r_best_c}")[:16]
                alt_subrec = SeqRecord(alt_seq, id=clean_id, name=clean_id, description=f"{args.target_name} ({r_best_c}:{alt_s}-{alt_e})")
                alt_subrec.annotations["molecule_type"] = "DNA"

                alt_features = []
                for c in alt_cds:
                    if c["end"] > alt_s and c["start"] < alt_e:
                        rel_s = max(0, c["start"] - alt_s)
                        rel_e = min(len(alt_seq), c["end"] - alt_s)
                        st = c["strand"]
                        t_pid = c["pid"]
                        trans_seq = target_faa_dict.get(t_pid, "")
                        if not trans_seq:
                            for alt_id in [t_pid.replace("gene", "prot"), t_pid.replace("prot", "gene"), f"{r_best_c}___{t_pid}"]:
                                if alt_id in target_faa_dict:
                                    trans_seq = target_faa_dict[alt_id]
                                    break
                        if not trans_seq and (rel_e - rel_s) >= 30:
                            try:
                                sub_dna = alt_seq[rel_s:rel_e] if st == 1 else alt_seq[rel_s:rel_e].reverse_complement()
                                trans_seq = str(sub_dna.translate(to_stop=True))
                            except Exception:
                                trans_seq = ""

                        # KO annotation lookup
                        t_hit = None
                        if target_anno_dict:
                            for cand_p in [t_pid, t_pid.replace("gene", "prot"), t_pid.replace("prot", "gene"), t_pid.split("___")[-1]]:
                                t_hit = lookup_target_anno(cand_p, target_anno_dict)
                                if t_hit:
                                    break
                        if not t_hit and r_cl:
                            for h in r_cl:
                                if h.get("t_pid") == t_pid:
                                    ref_qid = h.get("qid", "")
                                    orig_p = h.get("orig_pid", "")
                                    ref_hit = lookup_ref_flanking_anno(orig_p, ref_flank_annotations, genome_acc=r_acc) if ref_flank_annotations else None
                                    if ref_hit and ref_hit.get("ko"):
                                        t_hit = {
                                            "ko": ref_hit["ko"],
                                            "definition": f"{ref_hit.get('definition', '')} (homology {h['identity']:.0f}% to {r_acc})",
                                            "score": str(h.get("bitscore", "")),
                                            "threshold": str(ref_hit.get("threshold", "")),
                                            "evalue": ""
                                        }
                                        break

                        c_inf_alt = gff_coords.get(t_pid, {})
                        quals = format_target_gene_qualifiers(t_pid, t_hit, trans_seq, ko_definitions, c_info=c_inf_alt)
                        feat = SeqFeature.SeqFeature(
                            SeqFeature.FeatureLocation(rel_s, rel_e, strand=st),
                            type="CDS",
                            qualifiers=quals
                        )
                        alt_features.append(feat)
                alt_features.sort(key=lambda x: int(x.location.start))
                alt_subrec.features = alt_features
                with open(pair_target_gbk, "w") as fg:
                    SeqIO.write(alt_subrec, fg, "genbank")
                fix_gbk_wrapped_qualifiers(pair_target_gbk)
                print(f"Stage 1 Pairwise Target Locus GBK written for cross-contig ref {r_acc}: {pair_target_gbk} (contig {r_best_c}, {len(alt_features)} CDS)", flush=True)

    # Write ordered list of GBK files for Clinker visualization
    # References are ordered from most basal to closest relative (Rank N -> ... -> Rank 1)
    # so that Rank 1 is positioned directly adjacent to the target track in Clinker
    clinker_order_txt = os.path.join(args.out_refs_locus_dir, "ordered_clinker_files.txt")
    with open(clinker_order_txt, "w") as cof:
        for r_acc in reversed(visualized_refs):
            gbk_path = os.path.join(args.out_refs_locus_dir, f"{r_acc}_locus.gbk")
            if os.path.exists(gbk_path):
                cof.write(gbk_path + "\n")

    # 10. Extract extended ±100kb locus proteins table across target and positive references
    out_ext_tsv = args.out_extended_proteins_tsv
    if not out_ext_tsv:
        out_ext_tsv = os.path.join(os.path.dirname(args.out_topology_tsv), f"{args.target_name}_100kb_locus_proteins.tsv")
    extract_extended_100kb_locus_proteins(
        args=args,
        pos_refs=pos_refs,
        ref_candidates=ref_candidates,
        ref_evaluations=ref_evaluations,
        target_records=target_records,
        gff_coords=gff_coords,
        target_faa_dict=target_faa_dict,
        target_anno_dict=target_anno_dict,
        ref_flank_annotations=ref_flank_annotations,
        ko_definitions=ko_definitions,
        out_tsv_path=out_ext_tsv,
        target_contig_id=chosen_contig_id,
        winning_geometry=chosen_geometry
    )

    print("Stage 1 Extraction and Collinear Synteny Filtering Completed Successfully!", flush=True)

if __name__ == "__main__":
    main()
