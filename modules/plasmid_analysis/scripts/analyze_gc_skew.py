#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re

def read_fasta_all(file_path):
    sequences = []
    header, seq = "", []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header: sequences.append((header, "".join(seq).upper()))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header: sequences.append((header, "".join(seq).upper()))
    return sequences

def calculate_cumulative_gc_skew(seq, window_size=1000, step_size=100):
    seq_len = len(seq)
    cgc_skew, positions = [], []
    current_cgc = 0.0
    for i in range(0, seq_len - window_size + 1, step_size):
        window = seq[i:i+window_size]
        g_count, c_count = window.count('G'), window.count('C')
        skew = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0.0
        current_cgc += skew
        cgc_skew.append(current_cgc)
        positions.append(i + window_size // 2)
    return np.array(positions), np.array(cgc_skew)

def search_dnaa_boxes(search_seq):
    pattern_fwd = re.compile(r'TTAT[AC]CA[AC]A')
    pattern_rev = re.compile(r'T[GT]TG[GT]ATAA')
    boxes = []
    for match in pattern_fwd.finditer(search_seq): boxes.append(match.start())
    for match in pattern_rev.finditer(search_seq): boxes.append(match.start())
    return sorted(boxes)

def format_oric_sequence(oric_seq):
    pattern_fwd = re.compile(r'TTAT[AC]CA[AC]A')
    pattern_rev = re.compile(r'T[GT]TG[GT]ATAA')
    seq_chars = list(oric_seq.lower())
    upper_seq = oric_seq.upper()
    matches = list(pattern_fwd.finditer(upper_seq)) + list(pattern_rev.finditer(upper_seq))
    for match in matches:
        start, end = match.span()
        for i in range(start, end):
            seq_chars[i] = seq_chars[i].upper()
    return "".join(seq_chars)

def refine_oric_region(seq, approx_oric, bg_freq, window_radius=50000):
    """
    【升级核心】支持动态半径与全序列扫描
    """
    seq_len = len(seq)
    
    # 如果搜索窗口覆盖了全基因组（常见于质粒或设定了超大半径）
    if window_radius * 2 >= seq_len:
        search_region = seq
        start_idx = 0
    else:
        start_idx = int((approx_oric - window_radius) % seq_len)
        end_idx = int((approx_oric + window_radius) % seq_len)
        
        if start_idx < end_idx:
            search_region = seq[start_idx:end_idx]
        else:
            search_region = seq[start_idx:] + seq[:end_idx]
            
    boxes = search_dnaa_boxes(search_region)
    if not boxes:
        return None, None, 0, "", 0.0
        
    best_density = 0
    best_region = None
    
    # 聚类算法：寻找 1000bp 内最密集的 box 簇
    for i in range(len(boxes)):
        for j in range(i, len(boxes)):
            if boxes[j] - boxes[i] <= 1000:
                density = j - i + 1
                if density > best_density:
                    best_density = density
                    best_region = (max(0, boxes[i] - 50), boxes[j] + 59) 

    if best_density >= 3:
        exact_start = int((start_idx + best_region[0]) % seq_len)
        exact_end = int((start_idx + best_region[1]) % seq_len)
        
        if exact_start < exact_end:
            oric_seq = seq[exact_start:exact_end]
        else:
            oric_seq = seq[exact_start:] + seq[:exact_end]
            
        formatted_oric_seq = format_oric_sequence(oric_seq)
        
        region_len = len(oric_seq)
        expected_boxes = region_len * bg_freq
        enrichment = best_density / expected_boxes if expected_boxes > 0 else 0.0
            
        return exact_start, exact_end, best_density, formatted_oric_seq, enrichment
        
    return None, None, 0, "", 0.0

def main():
    parser = argparse.ArgumentParser(description="终极拓扑分析脚本 (支持宏观漂移校验与微观聚类扫描)")
    parser.add_argument("-i", "--input", required=True, help="输入的 FASTA 文件")
    parser.add_argument("-p", "--plot_prefix", help="输出图谱的文件名前缀 (可选)")
    parser.add_argument("--out_oric", help="输出提取到的精准 oriC 序列的 FASTA 文件 (可选)")
    parser.add_argument("-r", "--radius", type=int, default=50000, help="在 GC-skew 极值点两侧搜索 DnaA 簇的半径 (默认: 50000 bp)")
    parser.add_argument("--header", action="store_true", help="是否在 stdout 输出 TSV 表头")
    args = parser.parse_args()

    file_basename = os.path.basename(args.input)
    genome_id = os.path.splitext(file_basename)[0] 

    all_sequences = read_fasta_all(args.input)
    if not all_sequences:
        sys.exit(0)

    max_len = max([len(seq) for _, seq in all_sequences])

    if args.header:
        print("File\tContig_ID\tClass\tLength(bp)\tSkew_oriC\tSkew_ter\tSymmetry\tClosure\tTopology\tExact_oriC_Start\tExact_oriC_End\tDnaA_Boxes\tEnrichment")

    oric_fasta_buffer = []

    for header, seq in all_sequences:
        seq_len = len(seq)
        contig_id = header.split()[0]
        
        if seq_len < 10000: continue
            
        if seq_len == max_len:
            bio_class, strict_mode = "Primary_Chr", True
        elif seq_len >= 500000 or seq_len >= (max_len * 0.15):
            bio_class, strict_mode = "Secondary_Chr", True
        else:
            bio_class, strict_mode = "Plasmid", False

        total_boxes = len(search_dnaa_boxes(seq))
        bg_freq = total_boxes / seq_len if seq_len > 0 else 0

        positions, cgc_skew = calculate_cumulative_gc_skew(seq)
        if len(cgc_skew) < 10: continue
            
        min_idx, max_idx = np.argmin(cgc_skew), np.argmax(cgc_skew)
        approx_oric, approx_ter = positions[min_idx], positions[max_idx]
        
        symmetry_ratio = abs(approx_oric - approx_ter) / seq_len
        amplitude = cgc_skew[max_idx] - cgc_skew[min_idx]
        closure_ratio = abs(cgc_skew[-1] - cgc_skew[0]) / amplitude if amplitude > 0 else 1.0

        if strict_mode:
            is_symmetric, is_closed = 0.35 <= symmetry_ratio <= 0.65, closure_ratio < 0.20
        else:
            is_symmetric, is_closed = 0.10 <= symmetry_ratio <= 0.90, closure_ratio < 0.25 

        # 调用带自定义半径的 refine_oric_region
        exact_start, exact_end, box_count, formatted_oric_seq, enrichment = refine_oric_region(seq, approx_oric, bg_freq, args.radius)

        topology = "Linear_or_Frag"
        if is_symmetric and is_closed: 
            topology = "Circular"
        elif not is_symmetric and is_closed: 
            topology = "Circular_Asym"

        # 生物学豁免规则：跨界且富集度高，则强行判定为环状
        if topology == "Linear_or_Frag" and exact_start is not None and exact_end is not None:
            if exact_start > exact_end and enrichment >= 50.0:
                topology = "Circular_Rescued"

        ex_start_str = str(exact_start) if exact_start is not None else "NA"
        ex_end_str = str(exact_end) if exact_end is not None else "NA"
        
        print(f"{file_basename}\t{contig_id}\t{bio_class}\t{seq_len}\t{approx_oric}\t{approx_ter}\t{symmetry_ratio:.3f}\t{closure_ratio:.3f}\t{topology}\t{ex_start_str}\t{ex_end_str}\t{box_count}\t{enrichment:.1f}")

        if box_count >= 3 and args.out_oric:
            header_name = f">{genome_id}|{contig_id}|{bio_class}_oriC_{box_count}boxes_Enrich{enrichment:.1f}x_{exact_start}-{exact_end}\n"
            oric_fasta_buffer.append(header_name)
            for i in range(0, len(formatted_oric_seq), 80):
                oric_fasta_buffer.append(formatted_oric_seq[i:i+80] + "\n")

        if args.plot_prefix:
            plt.figure(figsize=(8, 4))
            line_color = '#ff7f0e' if topology == "Circular_Rescued" else '#1f77b4'
            plt.plot(positions, cgc_skew, color=line_color)
            
            plt.axvline(x=approx_oric, color='green', linestyle='--', alpha=0.4, label='GC-skew min (approx)')
            plt.axvline(x=approx_ter, color='red', linestyle='--', alpha=0.4, label='GC-skew max')
            
            if exact_start is not None:
                # 绘制一条红线，代表实际找到的 oriC
                plt.axvline(x=exact_start, color='darkgreen', linestyle='-', linewidth=2, label=f'Exact oriC ({box_count} boxes, {enrichment:.1f}x)')

            plt.title(f"{contig_id} ({bio_class} | {topology} | Len: {seq_len//1000}kb)")
            plt.legend()
            plt.tight_layout()
            
            clean_contig_id = contig_id.replace('|', '_')
            plt.savefig(f"{args.plot_prefix}_{bio_class}_{clean_contig_id}.png", dpi=150)
            plt.close()

    if args.out_oric and oric_fasta_buffer:
        with open(args.out_oric, 'a') as f:
            f.writelines(oric_fasta_buffer)

if __name__ == "__main__":
    main()