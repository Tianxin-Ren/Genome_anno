#!/usr/bin/env python3
"""
merge_kofam_results.py
Assembles and maps shard-level PyKofamSearch results back to the original proteome space.
Primary Artifact (Canonical Asset):
  - protein_KO_all_strict.tsv: All hits passing KOfam curated adaptive score thresholds (or relaxed threshold_factor),
    with fallback to best E-value when no hit passes score. Retains multiple KO candidates with ko_rank,
    score_margin, bitscore, evalue, is_strict, and pass_type.
    Matches schemas/v1/kofam_hit.schema.json.
Secondary Derivatives:
  - protein_KO_best.tsv: Derivative retaining only the single best-scoring/fallback KO per protein.
  - genome_KO_counts.long.tsv: Canonical long-format sparse table (genome_id, KO, count).
    Matches schemas/v1/genome_annotation.schema.json.
  - Optional: genome_KO_counts.dense.tsv & genome_KO_presence.dense.tsv
"""

import argparse
import gzip
import os
import sys

def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)

def load_ko_score_types(ko_list_path):
    """
    Loads mapping of id_ko -> score_type ('domain' or 'full') from:
    1. KOfam ko_list TSV file (header contains: koid, score_type, ...)
    2. Serialized database.pkl (pickle containing ko_to_data dict)
    """
    ko_types = {}
    if not ko_list_path or not os.path.exists(ko_list_path):
        return ko_types

    # Case 1: pickle format
    if ko_list_path.endswith((".pkl", ".pickle", ".pgz")):
        try:
            import pickle
            open_fn = gzip.open if ko_list_path.endswith((".gz", ".pgz")) else open
            with open_fn(ko_list_path, "rb") as f:
                data = pickle.load(f)
                if isinstance(data, tuple) and len(data) >= 1 and isinstance(data[0], dict):
                    ko_to_data = data[0]
                    for id_ko, k_info in ko_to_data.items():
                        if isinstance(k_info, dict) and "score_type" in k_info and k_info["score_type"]:
                            ko_types[id_ko] = str(k_info["score_type"]).strip().lower()
            print(f"Loaded {len(ko_types)} KO score types from serialized database {ko_list_path}")
            return ko_types
        except Exception as e:
            print(f"Warning: Failed to load score_type from pickle {ko_list_path}: {e}", file=sys.stderr)

    # Case 2: TSV ko_list format
    try:
        with open_file(ko_list_path, "rt") as f:
            header = None
            k_idx = 0
            st_idx = 2
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if header is None:
                    header = [p.lower() for p in parts]
                    k_idx = header.index("koid") if "koid" in header else 0
                    st_idx = header.index("score_type") if "score_type" in header else (2 if len(header) > 2 else 2)
                    continue
                if len(parts) > max(k_idx, st_idx):
                    koid = parts[k_idx]
                    st = parts[st_idx].strip().lower()
                    if st in ["domain", "full"]:
                        ko_types[koid] = st
        print(f"Loaded {len(ko_types)} KO score types from ko_list {ko_list_path}")
    except Exception as e:
        print(f"Warning: Failed to load score_type from {ko_list_path}: {e}", file=sys.stderr)

    return ko_types

def process_single_shard(args_tuple):
    """
    Parses a single PyKofamSearch shard TSV output and immediately arbitrates
    candidate hits for all sequences in this shard using in-worker early pruning.

    Memory Optimization:
      - Maintains at most 1 best fallback candidate per sequence (O(1) space).
      - Drops sub-threshold non-qualifying hits on the fly.
      - Drastically cuts peak RAM from dozens of gigabytes to megabytes.

    Returns:
      (arbitrated_hits_dict, raw_hit_count, score_pass_seq_count, fallback_seq_count)
    """
    if len(args_tuple) == 4:
        tsv_file, threshold_factor, evalue_fallback_cutoff, ko_score_types = args_tuple
    else:
        tsv_file, threshold_factor, evalue_fallback_cutoff = args_tuple
        ko_score_types = {}

    if not os.path.exists(tsv_file) or os.path.getsize(tsv_file) == 0:
        return {}, 0, 0, 0

    candidates = {} # u_id -> [list_of_score_passing_hits, best_fallback_hit_or_None]
    raw_hit_count = 0

    with open_file(tsv_file, "rt") as f:
        header = None
        q_col = 0
        t_col = 1
        s_col = 3
        th_col = 2
        st_col = None
        e_col = 4
        d_col = 5

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if header is None:
                header = [p.lower() for p in parts]

                def find_idx(candidates_names):
                    for c in candidates_names:
                        if c in header:
                            return header.index(c)
                    return None

                q_idx = find_idx(["id_protein", "query", "protein_id", "gene_id"])
                q_col = q_idx if q_idx is not None else 0

                t_idx = find_idx(["id_ko", "target", "ko"])
                t_col = t_idx if t_idx is not None else 1

                s_idx = find_idx(["score", "bitscore"])
                s_col = s_idx if s_idx is not None else (3 if len(header) > 3 else 2)

                th_idx = find_idx(["threshold", "curated_threshold"])
                th_col = th_idx if th_idx is not None else (2 if len(header) > 2 else 3)

                st_col = find_idx(["score_type", "type"])

                e_idx = find_idx(["e-value", "evalue", "e_value"])
                e_col = e_idx if e_idx is not None else (4 if len(header) > 4 else 5)

                d_idx = find_idx(["definition", "description", "ko_definition"])
                d_col = d_idx if d_idx is not None else (5 if len(header) > 5 else (len(parts) - 1))
                continue

            if len(parts) > max(q_col, t_col, s_col):
                raw_hit_count += 1
                u_id = parts[q_col].split()[0]
                ko = parts[t_col]
                try:
                    score = float(parts[s_col])
                except (ValueError, IndexError):
                    score = 0.0
                try:
                    threshold = float(parts[th_col]) if th_col is not None and th_col < len(parts) else 0.0
                except (ValueError, IndexError):
                    threshold = 0.0

                # Resolve score_type: check TSV column first, then ko_score_types mapping, else default "full"
                if st_col is not None and st_col < len(parts) and parts[st_col].strip().lower() in ["full", "domain"]:
                    score_type = parts[st_col].strip().lower()
                elif ko_score_types and ko in ko_score_types:
                    score_type = ko_score_types[ko]
                else:
                    score_type = "full"

                try:
                    evalue = float(parts[e_col]) if (e_col is not None and e_col < len(parts)) else 0.0
                except (ValueError, IndexError):
                    evalue = 0.0
                desc = parts[d_col] if (d_col is not None and d_col < len(parts)) else ""

                # Evaluate score threshold pass
                if threshold > 0.0:
                    effective_cutoff = threshold * threshold_factor
                    is_score_pass = (score >= effective_cutoff)
                    is_strict = (score >= threshold)
                else:
                    is_score_pass = (score > 0.0 and evalue <= 1e-5)
                    is_strict = is_score_pass

                if u_id not in candidates:
                    candidates[u_id] = [[], None]

                if is_score_pass:
                    candidates[u_id][0].append({
                        "ko": ko,
                        "bitscore": score,
                        "threshold": threshold,
                        "score_type": score_type,
                        "evalue": evalue,
                        "description": desc,
                        "margin": score - threshold,
                        "is_score_pass": True,
                        "is_strict": is_strict
                    })
                else:
                    # Early Pruning: only retain the single candidate with lowest E-value <= evalue_fallback_cutoff
                    if evalue <= evalue_fallback_cutoff:
                        curr_best_fb = candidates[u_id][1]
                        if curr_best_fb is None or (evalue, -score) < (curr_best_fb["evalue"], -curr_best_fb["bitscore"]):
                            candidates[u_id][1] = {
                                "ko": ko,
                                "bitscore": score,
                                "threshold": threshold,
                                "score_type": score_type,
                                "evalue": evalue,
                                "description": desc,
                                "margin": score - threshold,
                                "is_score_pass": False,
                                "is_strict": False
                            }

    # Arbitrate final surviving hits for this shard
    arbitrated = {}
    score_pass_seqs = 0
    fallback_seqs = 0

    for u_id, (score_passing, best_fb) in candidates.items():
        if score_passing:
            score_passing.sort(key=lambda x: (-x["bitscore"], x["evalue"]))
            for idx, h in enumerate(score_passing):
                h["ko_rank"] = idx + 1
                h["score_margin"] = round(h["margin"], 2)
                h["pass_type"] = "score_pass"
            arbitrated[u_id] = score_passing
            score_pass_seqs += 1
        elif best_fb is not None:
            best_fb["ko_rank"] = 1
            best_fb["score_margin"] = round(best_fb["margin"], 2)
            best_fb["is_strict"] = False
            best_fb["pass_type"] = "evalue_fallback"
            arbitrated[u_id] = [best_fb]
            fallback_seqs += 1

    return arbitrated, raw_hit_count, score_pass_seqs, fallback_seqs

def parse_pykofam_tsv(tsv_file, threshold_factor=0.80, evalue_fallback_cutoff=1e-3, ko_score_types=None):
    """
    Parses a PyKofamSearch TSV output file and returns arbitrated hits dict.
    Preserved for backward-compatibility.
    """
    if ko_score_types is None:
        ko_score_types = {}
    arbitrated, _, _, _ = process_single_shard((tsv_file, threshold_factor, evalue_fallback_cutoff, ko_score_types))
    return arbitrated

def main():
    parser = argparse.ArgumentParser(description="Merge and map PyKofamSearch shard hits.")
    parser.add_argument("-s", "--shards_dir", required=True, help="Directory containing shard TSV outputs")
    parser.add_argument("-m", "--mapping", required=True, help="protein_to_unique.tsv.gz mapping table")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for merged artifacts")
    parser.add_argument("-p", "--threads", type=int, default=32,
                        help="Number of parallel worker processes for shard parsing (default: 32)")
    parser.add_argument("--ko_count_strategy", default="best_hit", choices=["best_hit", "all_strict", "weighted"],
                        help="Strategy for accumulating genome KO counts when a protein has multiple qualifying hits (default: best_hit)")
    parser.add_argument("--threshold_factor", type=float, default=0.80,
                        help="Relative score threshold cutoff factor (default: 0.80, meaning score >= 0.80 * threshold)")
    parser.add_argument("--evalue_fallback_cutoff", type=float, default=1e-3,
                        help="E-value threshold for fallback hit selection when no hit passes score threshold (default: 1e-3)")
    parser.add_argument("--ko_list", default="",
                        help="Path to KOfam ko_list TSV or database.pkl to populate accurate 'domain' vs 'full' score_type")
    parser.add_argument("--export_dense", action="store_true", help="Also export optional dense genome x KO matrices")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # 0. Load KO score types if provided
    ko_score_types = load_ko_score_types(args.ko_list)

    # 1. Discover all shard TSV files
    shard_files = []
    for root, _, files in os.walk(args.shards_dir):
        for f in files:
            if f.endswith(".tsv") or f.endswith(".tsv.gz"):
                shard_files.append(os.path.join(root, f))
    shard_files.sort()

    if not shard_files:
        print(f"ERROR: No shard TSV files found in {args.shards_dir}!", file=sys.stderr)
        sys.exit(1)

    workers = min(max(1, args.threads), len(shard_files))
    print(f"Processing {len(shard_files)} shard TSV files using {workers} parallel worker processes "
          f"(threshold_factor={args.threshold_factor}, evalue_fallback_cutoff={args.evalue_fallback_cutoff})...")

    tasks = [(fpath, args.threshold_factor, args.evalue_fallback_cutoff, ko_score_types) for fpath in shard_files]

    unique_hits = {}
    total_raw_hits = 0
    score_pass_count = 0
    fallback_count = 0

    if workers > 1:
        import multiprocessing
        with multiprocessing.Pool(processes=workers) as pool:
            for s_hits, raw_cnt, sp_cnt, fb_cnt in pool.imap_unordered(process_single_shard, tasks):
                for u_id, h_list in s_hits.items():
                    if u_id not in unique_hits:
                        unique_hits[u_id] = h_list
                    else:
                        unique_hits[u_id].extend(h_list)
                total_raw_hits += raw_cnt
                score_pass_count += sp_cnt
                fallback_count += fb_cnt
    else:
        for task in tasks:
            s_hits, raw_cnt, sp_cnt, fb_cnt = process_single_shard(task)
            for u_id, h_list in s_hits.items():
                if u_id not in unique_hits:
                    unique_hits[u_id] = h_list
                else:
                    unique_hits[u_id].extend(h_list)
            total_raw_hits += raw_cnt
            score_pass_count += sp_cnt
            fallback_count += fb_cnt

    print(f"Loaded and arbitrated {total_raw_hits} raw KO hits across {len(unique_hits)} unique sequences from {len(shard_files)} shard files.")
    print(f"Arbitration complete: {score_pass_count} sequences with score-passing hits, {fallback_count} sequences rescued via E-value fallback (total annotated: {len(unique_hits)}).")

    # 2. Map back to original proteins using mapping table with buffered batch writing
    primary_strict_path = os.path.join(args.out_dir, "protein_KO_all_strict.tsv")
    best_path = os.path.join(args.out_dir, "protein_KO_best.tsv")
    long_matrix_path = os.path.join(args.out_dir, "genome_KO_counts.long.tsv")

    genome_ko_counts = {} # (genome_id, ko) -> count
    all_genomes = set()
    all_kos = set()

    BUFFER_SIZE = 50000
    ps_buffer = []
    b_buffer = []

    print(f"Streaming mapping table {args.mapping} to generate primary artifacts (buffer_size={BUFFER_SIZE})...")
    with open(primary_strict_path, "w", encoding="utf-8") as ps_out, \
         open(best_path, "w", encoding="utf-8") as b_out:

        # Write headers
        ps_out.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tko_rank\tscore_margin\tis_strict\tpass_type\n")
        b_out.write("protein_id\tgene_id\tgenome_id\tko\tko_definition\tbitscore\tcurated_threshold\tscore_type\tevalue\tscore_margin\tis_strict\tpass_type\n")

        with open_file(args.mapping, "rt") as mf:
            m_header = mf.readline().strip().split("\t")
            # Expected columns: protein_id, unique_id, aa_length, genome_id
            for line in mf:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    p_id = parts[0]
                    u_id = parts[1]
                    g_id = parts[3] if len(parts) > 3 else (p_id.split("___")[0] if "___" in p_id else "unknown")
                    # Construct gene_id if structured
                    if "___prot" in p_id:
                        gene_id = p_id.replace("___prot", "___gene")
                    else:
                        gene_id = p_id

                    all_genomes.add(g_id)

                    hits = unique_hits.get(u_id, [])
                    if hits:
                        # Write best hit (rank 1)
                        best_h = hits[0]
                        b_buffer.append(f"{p_id}\t{gene_id}\t{g_id}\t{best_h['ko']}\t{best_h['description']}\t{best_h['bitscore']}\t{best_h['threshold']}\t{best_h['score_type']}\t{best_h['evalue']}\t{best_h['score_margin']}\t{best_h['is_strict']}\t{best_h['pass_type']}\n")
                        if len(b_buffer) >= BUFFER_SIZE:
                            b_out.write("".join(b_buffer))
                            b_buffer.clear()

                        # Write all qualifying hits (all score-passing or the 1 fallback hit)
                        for h in hits:
                            ps_buffer.append(f"{p_id}\t{gene_id}\t{g_id}\t{h['ko']}\t{h['description']}\t{h['bitscore']}\t{h['threshold']}\t{h['score_type']}\t{h['evalue']}\t{h['ko_rank']}\t{h['score_margin']}\t{h['is_strict']}\t{h['pass_type']}\n")
                            if len(ps_buffer) >= BUFFER_SIZE:
                                ps_out.write("".join(ps_buffer))
                                ps_buffer.clear()

                        # Accumulate count according to specified strategy
                        if args.ko_count_strategy == "best_hit":
                            key = (g_id, best_h["ko"])
                            genome_ko_counts[key] = genome_ko_counts.get(key, 0) + 1
                            all_kos.add(best_h["ko"])
                        elif args.ko_count_strategy == "all_strict":
                            for h in hits:
                                key = (g_id, h["ko"])
                                genome_ko_counts[key] = genome_ko_counts.get(key, 0) + 1
                                all_kos.add(h["ko"])
                        elif args.ko_count_strategy == "weighted":
                            weight = 1.0 / len(hits)
                            for h in hits:
                                key = (g_id, h["ko"])
                                genome_ko_counts[key] = genome_ko_counts.get(key, 0.0) + weight
                                all_kos.add(h["ko"])

        # Flush remaining buffers
        if b_buffer:
            b_out.write("".join(b_buffer))
            b_buffer.clear()
        if ps_buffer:
            ps_out.write("".join(ps_buffer))
            ps_buffer.clear()


    # 3. Write secondary derivative: genome_KO_counts.long.tsv
    print(f"Writing long-format sparse matrix {long_matrix_path}...")
    with open(long_matrix_path, "w", encoding="utf-8") as lf_out:
        lf_out.write("genome_id\tko\tcount\n")
        for (g_id, ko), cnt in sorted(genome_ko_counts.items()):
            lf_out.write(f"{g_id}\t{ko}\t{cnt}\n")

    # 4. Optional Dense Export
    if args.export_dense:
        dense_count_path = os.path.join(args.out_dir, "genome_KO_counts.dense.tsv")
        dense_pres_path = os.path.join(args.out_dir, "genome_KO_presence.dense.tsv")
        print("Writing optional dense matrices...")
        sorted_kos = sorted(all_kos)
        sorted_genomes = sorted(all_genomes)

        with open(dense_count_path, "w", encoding="utf-8") as c_out, \
             open(dense_pres_path, "w", encoding="utf-8") as p_out:

            
            c_out.write("genome_id\t" + "\t".join(sorted_kos) + "\n")
            p_out.write("genome_id\t" + "\t".join(sorted_kos) + "\n")

            for g in sorted_genomes:
                counts = [str(genome_ko_counts.get((g, ko), 0)) for ko in sorted_kos]
                pres = ["1" if genome_ko_counts.get((g, ko), 0) > 0 else "0" for ko in sorted_kos]
                c_out.write(f"{g}\t" + "\t".join(counts) + "\n")
                p_out.write(f"{g}\t" + "\t".join(pres) + "\n")

    print(f"Merging complete. Primary asset: {primary_strict_path}")

if __name__ == "__main__":
    main()
