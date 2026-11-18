#!/usr/bin/env python3
import argparse
import numpy as np
import sys
import os
import re
import gzip

# =====================================================================
# 策略库：涵盖细菌、高GC菌、古菌及质粒
# 约定：第一个模式通常为正向，第二个为反向（针对 DnaA box 等）
# =====================================================================
MOTIF_STRATEGIES = {
    "Bac_DnaA": [r'TTAT[AC]CA[AC]A', r'T[GT]TG[GT]ATAA'],
    "Bac_Actino": [r'TT[GT]TCCACA', r'TGTGG[AC]AA'],
    "Arc_ORB": [r'GTT[CT]C[ACGT]{1,4}G[AG]AAC', r'C[TC]GTT[CT]C'],
    "Plasmid_Rep": [r'CTTG[ACGT]{1,2}CAAG', r'TG[ACGT]{2}CA[ACGT]{2}TG']
}

def get_at_content(seq_str):
    if not seq_str: return 0.0
    s = seq_str.upper()
    at_count = s.count('A') + s.count('T')
    gc_count = s.count('G') + s.count('C')
    total = at_count + gc_count
    return (at_count / total) * 100 if total > 0 else 0.0

def find_overlap(seq, min_overlap=50, max_overlap=20000):
    """检测基因组首尾是否存在物理重叠"""
    seq_len = len(seq)
    if seq_len < min_overlap * 2: return 0, ""
    seed = seq[:min_overlap]
    search_start = max(min_overlap, seq_len - max_overlap - min_overlap)
    search_area = seq[search_start:]
    idx = search_area.find(seed)
    while idx != -1:
        actual_idx = search_start + idx
        overlap_len = seq_len - actual_idx
        if seq[actual_idx:] == seq[:overlap_len]:
            return overlap_len, seq[:overlap_len]
        idx = search_area.find(seed, idx + 1)
    return 0, ""

def calculate_multi_curves(seq):
    # 使用numpy加速大型基因组计算
    s_arr = np.frombuffer(seq.encode('ascii'), dtype=np.int8)
    a_c = (s_arr == ord('A')).cumsum()
    c_c = (s_arr == ord('C')).cumsum()
    g_c = (s_arr == ord('G')).cumsum()
    t_c = (s_arr == ord('T')).cumsum()
    ry = (a_c + g_c) - (c_c + t_c)
    mk = (a_c + c_c) - (g_c + t_c)
    # GC = g-c, AT = a-t
    return {"GC": (ry - mk) / 2.0, "AT": (ry + mk) / 2.0, "RY": ry, "MK": mk}

def evaluate_strategy(seq, approx_oric, bg_freq, patterns, radius):
    seq_len = len(seq)
    start_idx = int((approx_oric - radius) % seq_len)
    end_idx = int((approx_oric + radius) % seq_len)
    region = seq[start_idx:end_idx] if start_idx < end_idx else seq[start_idx:] + seq[:end_idx]
    
    # 记录所有 motif 的位置和方向
    all_hits = []
    for p_idx, p in enumerate(patterns):
        strand = "+" if p_idx == 0 else "-"
        for m in re.finditer(p, region):
            all_hits.append({"pos": m.start(), "strand": strand, "end": m.end()})
    
    all_hits.sort(key=lambda x: x["pos"])
    if len(all_hits) < 3: return None
    
    # 寻找 1000bp 内最密集的区域
    best_count, best_reg, best_strands = 0, None, ""
    for i in range(len(all_hits)):
        for j in range(i + 2, len(all_hits)): # 至少 3 个
            if all_hits[j]["pos"] - all_hits[i]["pos"] <= 1000:
                count = j - i + 1
                if count >= best_count:
                    best_count = count
                    best_reg = (max(0, all_hits[i]["pos"]-50), all_hits[j]["end"]+50)
                    best_strands = "".join([h["strand"] for h in all_hits[i:j+1]])

    if best_count >= 3:
        ex_s = int((start_idx + best_reg[0]) % seq_len)
        ex_e = int((start_idx + best_reg[1]) % seq_len)
        o_seq = seq[ex_s:ex_e] if ex_s < ex_e else seq[ex_s:] + seq[:ex_e]
        
        # 格式化序列：大写 motif
        fmt_chars = list(o_seq.lower())
        up_o = o_seq.upper()
        for p in patterns:
            for m in re.finditer(p, up_o):
                for k in range(m.start(), m.end()):
                    if k < len(fmt_chars): fmt_chars[k] = fmt_chars[k].upper()
        
        enrich = best_count / (len(o_seq) * bg_freq) if bg_freq > 0 else 0.0
        return (ex_s, ex_e, best_count, "".join(fmt_chars), enrich, best_strands)
    return None

def main():
    parser = argparse.ArgumentParser(description="Analyze Z-curve and DnaA box to find oriC and determine genome topology.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file (supports .gz)")
    parser.add_argument("--out_oric", help="Path to save predicted oriC sequences")
    parser.add_argument("-r", "--radius", type=int, default=50000, help="Search radius around skew extrema (default: 50000)")
    parser.add_argument("--header", action="store_true", help="Print TSV header")
    parser.add_argument("--min_overlap", type=int, default=50, help="Min overlap for circularity check (default: 50)")
    parser.add_argument("--min_ratio", type=float, default=0.7, help="Min length ratio to be classified as Primary_Chr (default: 0.7)")
    args = parser.parse_args()

    # 读取序列
    all_seqs = []
    header, current_seq = "", []
    open_func = gzip.open if args.input.endswith(".gz") else open
    with open_func(args.input, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header: all_seqs.append((header, "".join(current_seq).upper()))
                header = line[1:]; current_seq = []
            else: current_seq.append(line)
        if header: all_seqs.append((header, "".join(current_seq).upper()))

    if not all_seqs: return
    max_l = max(len(s[1]) for s in all_seqs)

    if args.header:
        columns = [
            "File", "ID", "Class", "Topology", "Confidence", 
            "Phys_Over", "Sym", "Obs_Clo", "Exp_Clo", "Clo_Err", 
            "Gen_AT%", "Ori_AT%", "AT_Diff", "Enrich", "Anchor", 
            "Strategy", "Start", "End", "Boxes_Strand"
        ]
        print("\t".join(columns))

    for h, seq in all_seqs:
        s_len = len(seq)
        if s_len < 5000: continue
        
        over_len, over_seq = find_overlap(seq, min_overlap=args.min_overlap)
        
        b_class = "Primary_Chr" if s_len >= max_l * args.min_ratio else ("Secondary_Chr" if s_len > 500000 else "Plasmid")
        gen_at = get_at_content(seq)
        curves = calculate_multi_curves(seq)
        
        # 3. 矩阵竞标 (寻找 oriC)
        best = {"enrich": -1.0, "anchor": "Min_GC", "motif": "None", "start": "NA", "end": "NA", "count": 0, "seq": "", "strands": ""}
        anchors = {}
        for k, v in curves.items():
            anchors[f"Min_{k}"], anchors[f"Max_{k}"] = int(np.argmin(v)), int(np.argmax(v))

        for a_n, a_p in anchors.items():
            for m_n, pats in MOTIF_STRATEGIES.items():
                bg_count = len(set([m.start() for p in pats for m in re.finditer(p, seq)]))
                bg_f = bg_count / s_len if s_len > 0 else 0
                res = evaluate_strategy(seq, a_p, bg_f, pats, args.radius)
                if res and res[4] > best["enrich"]:
                    best.update({
                        "enrich": res[4], "anchor": a_n, "motif": m_n, 
                        "start": res[0], "end": res[1], "count": res[2], 
                        "seq": res[3], "strands": res[5]
                    })

        # 4. 拓扑校验
        c_k = best["anchor"].split("_")[1] if "_" in best["anchor"] else "GC"
        v_c = curves[c_k]
        a_p = anchors.get(best["anchor"], 0)
        t_p = int(np.argmax(v_c)) if "Min" in best["anchor"] else int(np.argmin(v_c))
        
        sym = abs(a_p - t_p) / s_len
        amp = abs(v_c[t_p] - v_c[a_p])
        obs_c = abs(v_c[-1] - v_c[0]) / amp if amp > 0 else 1.0
        exp_c = abs(2.0 - (1.0/sym)) if sym > 0 else 1.0
        err = abs(obs_c - exp_c)
        
        is_phys = over_len > 0
        is_genomic = err <= 0.20 and 0.40 <= sym <= 0.60
        
        if is_phys and is_genomic: topo = "Circular (Confirmed)"
        elif is_phys: topo = "Circular (Physical)"
        elif is_genomic: topo = "Circular (Genomic)"
        elif err > 0.30: topo = "Linear"
        else: topo = "Uncertain"
        
        ori_at = get_at_content(best["seq"])
        at_diff = ori_at - gen_at if best["seq"] else 0.0
        
        conf = "Low"
        if best["count"] >= 3:
            if topo.startswith("Circular") and (at_diff > 4.0 or best["enrich"] > 30):
                conf = "HIGH"
            elif best["enrich"] > 10 or at_diff > 2.0:
                conf = "Med"

        # 输出行列信息
        boxes_info = f"{best['count']} ({best['strands']})" if best['strands'] else str(best['count'])
        row = [
            os.path.basename(args.input), h.split()[0], b_class, topo, conf,
            str(over_len), f"{sym:.3f}", f"{obs_c:.3f}", f"{exp_c:.3f}", f"{err:.3f}",
            f"{gen_at:.1f}", f"{ori_at:.1f}", f"{at_diff:.1f}", f"{best['enrich']:.1f}", 
            best['anchor'], best['motif'], str(best['start']), str(best['end']), boxes_info
        ]
        print("\t".join(row))

        if args.out_oric and best["count"] >= 3:
            with open(args.out_oric, 'a') as f:
                f.write(f">{h.split()[0]}|{b_class}|Conf:{conf}|Topo:{topo}|Boxes:{boxes_info}\n{best['seq']}\n")

if __name__ == "__main__":
    main()