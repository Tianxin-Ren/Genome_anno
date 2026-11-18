import pandas as pd
import numpy as np
from Bio import Align
import sys
import os

def get_features(row, prefix):
    k_ko = str(row.get("Kofam_KO", "")).strip()
    e_ko = str(row.get("EggNOG_KO", "")).strip()
    og = str(row.get("EggNOG_OGs", "")).strip()
    pid = str(row.get("Protein_ID", "")).strip()
    
    f = ""
    if k_ko and k_ko != "nan" and k_ko != "": f = k_ko.split(";")[0].strip()
    elif e_ko and e_ko != "nan" and e_ko != "": f = e_ko.split(",")[0].strip()
    elif og and og != "nan" and og != "": f = og.split("@")[0].strip()
    else: f = "[UNANNOTATED]"
    return f, pid

def main(path1, path2, out_file):
    df1 = pd.read_csv(path1, sep="\t")
    df2 = pd.read_csv(path2, sep="\t")
    name1 = os.path.basename(path1).replace("_functional_summary.tsv", "")
    name2 = os.path.basename(path2).replace("_functional_summary.tsv", "")
    
    # 1. Prepare Features
    f1 = [get_features(r, name1)[0] for _, r in df1.iterrows()]
    f2 = [get_features(r, name2)[0] for _, r in df2.iterrows()]

    # 2. Setup Aligner
    aligner = Align.PairwiseAligner()
    all_f = sorted(list(set(f1 + f2)))
    f_to_char = {f: chr(i % 5000 + 33) for i, f in enumerate(all_f)}
    matrix = {}
    for fa in all_f:
        for fb in all_f:
            if fa == fb:
                matrix[(f_to_char[fa], f_to_char[fb])] = 1 if fa == "[UNANNOTATED]" else 20
            else:
                matrix[(f_to_char[fa], f_to_char[fb])] = -10
    
    aligner.substitution_matrix = Align.substitution_matrices.Array(data=matrix)
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    # 3. Automatic Circular Detection (Optimized to avoid Overflow)
    best_total_score = -float('inf')
    best_rev = False
    best_idx2_start = 0

    for rev in [False, True]:
        curr_f2 = f2[::-1] if rev else f2
        target_str = "".join([f_to_char[f] for f in curr_f2 + curr_f2])
        query_str = "".join([f_to_char[f] for f in f1])
        
        # Step A: Find best score first
        aligner.mode = 'local'
        try:
            # We don't call align() here to avoid the iterator explosion
            score = aligner.score(target_str, query_str)
            if score > best_total_score:
                # Step B: Only get one alignment if score is better
                als = aligner.align(target_str, query_str)
                # Use next() to avoid len() which causes the OverflowError
                it = iter(als)
                best_al = next(it, None)
                if best_al:
                    best_total_score = score
                    best_rev = rev
                    best_idx2_start = best_al.aligned[0][0][0] % len(f2)
        except Exception as e:
            print(f"Warning during orientation search: {e}")

    # 4. Final Re-alignment at the optimal rotation
    f2_proc = f2[::-1] if best_rev else f2
    f2_rotated = f2_proc[best_idx2_start:] + f2_proc[:best_idx2_start]
    
    aligner.mode = 'global'
    str1 = "".join([f_to_char[f] for f in f1])
    str2 = "".join([f_to_char[f] for f in f2_rotated])
    
    final_als = aligner.align(str1, str2)
    best_final = next(iter(final_als)) # Fixed to avoid OverflowError
    
    a1_str, a2_str = format(best_final).split("\n")[0], format(best_final).split("\n")[2]

    # 5. Build Result
    res = []
    p1, p2 = 0, 0
    df2_proc = (df2.iloc[::-1] if best_rev else df2).reset_index(drop=True)
    
    # Mapping back char to feature
    c_to_f = {v: k for k, v in f_to_char.items()}
    
    for c1, c2 in zip(a1_str, a2_str):
        row = {}
        feat1 = c_to_f.get(c1, "-") if c1 != "-" else "-"
        feat2 = c_to_f.get(c2, "-") if c2 != "-" else "-"
        
        if c1 != "-":
            orig = df1.iloc[p1]
            row["GenA_ID"], row["GenA_KO"] = orig["Protein_ID"], feat1
            p1 += 1
        else: row["GenA_ID"] = row["GenA_KO"] = "-"
            
        if c2 != "-":
            orig = df2_proc.iloc[(best_idx2_start + p2) % len(df2_proc)]
            row["GenB_ID"], row["GenB_KO"] = orig["Protein_ID"], feat2
            p2 += 1
        else: row["GenB_ID"] = row["GenB_KO"] = "-"
            
        if feat1 == feat2 and feat1 != "-":
            row["Status"] = "==" if feat1 != "[UNANNOTATED]" else ".."
        elif feat1 != "-" and feat2 != "-": row["Status"] = "!="
        elif feat1 != "-": row["Status"] = "->"
        else: row["Status"] = "<-"
        res.append(row)

    with open(out_file, "w") as f:
        f.write(f"# --- Smart Circular Synteny (Robust Mode) ---\n")
        f.write(f"# Orientation: {'REVERSED' if best_rev else 'FORWARD'}\n")
        f.write(f"# Optimal Rotation Index: {best_idx2_start}\n")
        f.write(f"# Score: {best_total_score}\n")
        pd.DataFrame(res).to_csv(f, sep="\t", index_label="Line")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
