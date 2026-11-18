#!/usr/bin/env python3
"""
Ultra-fast Bakta merge script.
Optimized for ID consistency and supports FAA/FFN/GFF/TSV/GBK merging.
"""
import sys, os, re, argparse
from collections import defaultdict

def simplify_contig_id(full_name):
    # 1. Megahit: k141_12345
    m = re.match(r'^(k\d+_\d+)', full_name)
    if m: return m.group(1)
    
    # 2. SPAdes: NODE_1_length_...
    m = re.match(r'^(NODE_\d+)', full_name)
    if m: return m.group(1)
    
    # 3. Traditional scaffold/ctg naming
    m = re.match(r'^([a-zA-Z0-9]+)ctg', full_name)
    if m: return m.group(1)
    
    # 4. Fallback: Keep the first word of the ID (up to first space or end of string)
    return full_name.split()[0] if full_name.strip() else full_name

def merge_bakta(gff_list, faa_list, ffn_list, gbk_list, out_gff, out_faa, out_ffn, out_tsv, out_gbk):
    gff_list = [os.path.abspath(p) for p in gff_list]
    faa_list = [os.path.abspath(p) for p in faa_list]
    ffn_list = [os.path.abspath(p) for p in ffn_list]
    gbk_list = [os.path.abspath(p) for p in gbk_list]

    global_maps = defaultdict(dict)
    counters = defaultdict(int)
    
    # --- Phase 1: Build Mapping ---
    print(f"[merge] Phase 1: Building mapping...", file=sys.stderr)
    for gff_path in gff_list:
        if not os.path.exists(gff_path): continue
        with open(gff_path) as f:
            for line in f:
                if line.startswith('#') or not line.strip() or line.startswith('>'): continue
                parts = line.split('\t')
                if len(parts) < 9: continue
                
                contig_simple = simplify_contig_id(parts[0])
                attr = parts[8]
                
                m_lt = re.search(r'locus_tag=([^; \n\t]+)', attr)
                m_id = re.search(r'ID=([^; \n\t]+)', attr)
                
                if m_lt and m_id:
                    old_id = m_id.group(1)
                    old_lt = m_lt.group(1)
                    if old_id not in global_maps[gff_path]:
                        counters[contig_simple] += 1
                        new_id = f"{contig_simple}_{counters[contig_simple]:05d}"
                        global_maps[gff_path][old_id] = new_id
                        global_maps[gff_path][old_lt] = new_id
                elif m_id:
                    old_id = m_id.group(1)
                    global_maps[gff_path][old_id] = f"{contig_simple}_ext_{old_id}"

    word_pattern = re.compile(r'[a-zA-Z0-9_.-]+')
    def fast_replace(line, mapping):
        if not mapping: return line
        def subst(match):
            word = match.group(0)
            return mapping.get(word, word)
        return word_pattern.sub(subst, line)

    # --- Phase 2: GFF ---
    fasta_lines = []
    with open(out_gff, 'w') as out:
        out.write("##gff-version 3\n")
        for gff_path in gff_list:
            mapping = global_maps[gff_path]
            in_fasta_section = False
            with open(gff_path) as f:
                for line in f:
                    if line.startswith('##gff-version') or not line.strip(): continue
                    if line.startswith('##FASTA'):
                        in_fasta_section = True
                        continue
                    if in_fasta_section:
                        fasta_lines.append(line)
                        continue
                    if line.startswith('#'):
                        out.write(line); continue
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        parts[8] = fast_replace(parts[8], mapping)
                        out.write('\t'.join(parts))
                    else:
                        out.write(line)
        if fasta_lines:
            out.write("##FASTA\n")
            for line in fasta_lines:
                out.write(line)

    # --- Phase 3: FAA ---
    with open(out_faa, 'w') as out:
        for faa_path in faa_list:
            gff_ref = faa_path.replace('.faa', '.gff3')
            mapping = global_maps.get(gff_ref, {})
            with open(faa_path) as f:
                for line in f:
                    out.write(fast_replace(line, mapping) if line.startswith('>') else line)

    # --- Phase 4: FFN (Predicted genes, including non-coding) ---
    print(f"[merge] Phase 4: Writing FFN...", file=sys.stderr)
    with open(out_ffn, 'w') as out:
        for ffn_path in ffn_list:
            gff_ref = ffn_path.replace('.ffn', '.gff3')
            mapping = global_maps.get(gff_ref, {})
            with open(ffn_path) as f:
                for line in f:
                    out.write(fast_replace(line, mapping) if line.startswith('>') else line)

    # --- Phase 5: TSV ---
    with open(out_tsv, 'w') as out:
        hdr_written = False
        for gff_path in gff_list:
            tsv = gff_path.replace('.gff3', '.tsv')
            if not os.path.exists(tsv): continue
            mapping = global_maps[gff_path]
            with open(tsv) as f:
                for line in f:
                    if line.startswith('#'):
                        if not hdr_written: out.write(line); hdr_written = True
                        continue
                    out.write(fast_replace(line, mapping))

    # --- Phase 6: GBK ---
    with open(out_gbk, 'w') as out:
        for gbk_path in gbk_list:
            gff_ref = gbk_path.replace('.gbk', '.gff3')
            mapping = global_maps.get(gff_ref, {})
            with open(gbk_path) as f:
                for line in f:
                    out.write(fast_replace(line, mapping) if '/' in line else line)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--gffs", nargs="+", required=True)
    p.add_argument("--faas", nargs="+", required=True)
    p.add_argument("--ffns", nargs="+", required=True)
    p.add_argument("--gbks", nargs="+", default=[])
    p.add_argument("--out_gff", required=True)
    p.add_argument("--out_faa", required=True)
    p.add_argument("--out_ffn", required=True)
    p.add_argument("--out_tsv", required=True)
    p.add_argument("--out_gbk", required=True)
    a = p.parse_args()
    merge_bakta(a.gffs, a.faas, a.ffns, a.gbks, a.out_gff, a.out_faa, a.out_ffn, a.out_tsv, a.out_gbk)
