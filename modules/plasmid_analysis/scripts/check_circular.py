#!/usr/bin/env python3
import argparse
import sys

def read_fasta(file_path):
    """极其轻量级的 FASTA 解析器，生成 (header, sequence)"""
    header = ""
    seq = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

def find_overlap(seq, min_overlap=50, max_overlap=20000):
    """
    寻找序列首尾的完美重叠区域。
    返回: overlap_length, overlap_sequence
    """
    seq_len = len(seq)
    if seq_len < min_overlap * 2:
        return 0, ""

    # 提取序列开头的 'min_overlap' 长度作为搜索种子
    seed = seq[:min_overlap]
    
    # 为了提速，只在序列尾部 (max_overlap + min_overlap) 的范围内搜索
    search_start = max(min_overlap, seq_len - max_overlap - min_overlap)
    search_area = seq[search_start:]
    
    # 在尾部区域寻找种子的位置
    idx = search_area.find(seed)
    
    while idx != -1:
        # 计算在完整序列中的真实物理索引
        actual_idx = search_start + idx
        # 假设这就是重叠的起点，计算重叠总长度
        overlap_len = seq_len - actual_idx
        
        # 严格验证：尾部序列是否与头部序列完全一致
        if seq[actual_idx:] == seq[:overlap_len]:
            return overlap_len, seq[:overlap_len]
        
        # 如果不完全匹配，继续寻找下一个可能的种子位置
        idx = search_area.find(seed, idx + 1)
        
    return 0, ""

def main():
    parser = argparse.ArgumentParser(description="检测基因组序列是否成环，并输出重叠区域。")
    parser.add_argument("-i", "--input", required=True, help="输入的 FASTA 文件路径")
    parser.add_argument("-o", "--out_fasta", help="输出重叠区域序列到指定的 FASTA 文件 (可选)")
    parser.add_argument("--min", type=int, default=50, help="最小重叠长度阈值 (默认: 50 bp)")
    parser.add_argument("--max", type=int, default=20000, help="最大重叠搜索范围 (默认: 20000 bp)")
    
    args = parser.parse_args()

    found_circular = False

    for header, seq in read_fasta(args.input):
        overlap_len, overlap_seq = find_overlap(seq, args.min, args.max)
        
        if overlap_len > 0:
            found_circular = True
            print(f"✅ [成环] 序列 '{header}': 发现首尾重叠！")
            print(f"   -> 基因组总长: {len(seq)} bp")
            print(f"   -> 重叠区长度: {overlap_len} bp")
            
            if args.out_fasta:
                with open(args.out_fasta, 'a') as out_f:
                    out_f.write(f">{header}_overlap_{overlap_len}bp\n")
                    # 每 80 个字符换行，保持 FASTA 格式规范
                    for i in range(0, len(overlap_seq), 80):
                        out_f.write(overlap_seq[i:i+80] + "\n")
                print(f"   -> 重叠序列已保存至: {args.out_fasta}")
            else:
                print(f"   -> 重叠序列预览: {overlap_seq[:30]}...{overlap_seq[-30:]}")
        else:
            print(f"❌ [未成环] 序列 '{header}': 未在指定范围内检测到完全重叠。")

    if not found_circular:
        sys.exit(1) # 如果需要集成到流程中，可以通过退出码判断状态

if __name__ == "__main__":
    main()