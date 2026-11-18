#!/usr/bin/env python3
import sys, os, re, math, csv, json
import pandas as pd
import argparse

def parse_hmmer_tbl(file_path):
    records = []
    if not os.path.exists(file_path): return pd.DataFrame()
    is_hmmscan = False
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"):
                if "hmmscan" in line:
                    is_hmmscan = True
                continue
            parts = line.split()
            if len(parts) < 19: continue
            
            # hmmscan target is HMM (parts[0]) and query is Protein (parts[2])
            # hmmsearch target is Protein (parts[0]) and query is HMM (parts[2])
            if is_hmmscan:
                target_name = parts[2]
                query_name = parts[0]
            else:
                target_name = parts[0]
                query_name = parts[2]
                
            full_hit_name = f"{query_name}|{parts[3]}" if parts[3] != "-" else query_name
            evalue = float(parts[4])
            full_score = float(parts[5])
            # tblout columns 7/8 are the best-domain E-value/score.  Keeping
            # this separate is important for KOfam models whose GA threshold
            # is defined on the domain rather than the full sequence.
            domain_evalue = float(parts[7])
            domain_score = float(parts[8])
            records.append({
                'Protein': target_name,
                'HMM_Match': full_hit_name,
                'Raw_Query': query_name,
                'E_value': evalue,
                'Score': full_score,  # backwards-compatible alias
                'Full_Score': full_score,
                'Domain_Score': domain_score,
                'Domain_E_value': domain_evalue,
            })
    return pd.DataFrame(records)

def parse_hmmer_domtbl(file_path):
    records = []
    if not os.path.exists(file_path): return pd.DataFrame()
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 23: continue
            target_name = parts[0]
            query_name = parts[3]
            query_acc = parts[4]
            full_hit_name = f"{query_acc}|{query_name}" if query_acc != "-" else query_name
            i_evalue = float(parts[12])
            score = float(parts[13])
            records.append({
                'Protein': target_name,
                'HMM_Match': full_hit_name,
                'E_value': i_evalue,
                'Score': score,
                'Full_Score': None,
                'Domain_Score': score,
                'Domain_E_value': i_evalue,
            })
    return pd.DataFrame(records)

def load_ko_thresholds(ko_list_path):
    ko_thresholds = {}
    if ko_list_path and os.path.exists(ko_list_path):
        with open(ko_list_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                ko = (row.get('knum') or '').strip().upper()
                if not ko.startswith('K'):
                    continue
                raw = (row.get('threshold') or '').strip()
                try:
                    threshold = float(raw)
                except (TypeError, ValueError):
                    threshold = None
                ko_thresholds[ko] = {
                    'threshold': threshold,
                    'score_type': (row.get('score_type') or '').strip().lower(),
                    'mlen': float(row['mlen']) if row.get('mlen') else None,
                }
    return ko_thresholds


def load_target_negatives(path):
    """Load target-specific negatives; negatives are never global."""
    result = {}
    if not path or not os.path.exists(path):
        return result
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            target = (row.get('target_ko') or '').strip().upper()
            negative = (row.get('negative_ko') or '').strip().upper()
            if target and negative:
                result.setdefault(target, []).append({
                    'ko': negative,
                    'type': (row.get('candidate_type') or '').strip(),
                    'reason': (row.get('reason') or '').strip(),
                })
    return result


def load_fasta_stats(faa_path):
    stats = {}
    if not faa_path or not os.path.exists(faa_path):
        return stats
    current = None
    seq = []
    def save():
        if current is not None:
            s = ''.join(seq).replace('*', '')
            stats[current] = {'length': len(s), 'internal_stop': '*' in ''.join(seq[:-1])}
    with open(faa_path) as f:
        for line in f:
            if line.startswith('>'):
                save()
                current = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        save()
    return stats


def model_type(raw_query, target_set):
    q = str(raw_query).upper()
    if re.fullmatch(r'K\d{5}', q):
        return 'KOfam_KO'
    if target_set and (raw_query in target_set or q in {x.upper() for x in target_set}):
        return 'HYdDB_subtype'
    if any(x in q.lower() for x in ('hyd', 'nife', 'fefe', 'fe_')):
        return 'HYdDB_subtype'
    return 'other'

def resolve_single_mag(unified_path, pfam_path, ko_thresholds, faa_path, output_path,
                       evalue_cutoff=1e-5, log_factor=0.2, target_set=set(),
                       target_negatives=None, competition_fraction=0.10):
    # 2. 解析文件
    print(f"Parsing HMMer results...", file=sys.stderr)
    df_unified = parse_hmmer_tbl(unified_path)
    df_pfam = parse_hmmer_domtbl(pfam_path)

    # 3. 预处理 Pfam (分组)
    pfam_dict = {}
    if not df_pfam.empty:
        df_pfam_valid = df_pfam[df_pfam['E_value'] <= evalue_cutoff]
        for protein, group in df_pfam_valid.groupby('Protein'):
            pfam_dict[protein] = ",".join(group['HMM_Match'].unique())

    # 4. 预处理 Unified (分组，避免在循环中切片)
    unified_groups = {}
    if not df_unified.empty:
        for protein, group in df_unified.groupby('Protein'):
            unified_groups[protein] = group.sort_values('E_value')
    target_negatives = target_negatives or {}
    fasta_stats = load_fasta_stats(faa_path)

    # 5. 核心裁决 (高性能流式迭代)
    print(f"Starting arbitration...", file=sys.stderr)
    
    # 确保输出目录存在
    out_dir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(out_dir, exist_ok=True)
    
    processed_proteins = set()
    
    with open(output_path, "w") as out:
        out.write("Protein_ID\tScore_Annotation\tEvalue_Annotation\tBest_Evalue\tTarget_Annotation\tPfam_Domains\tAll_Similar_Hits\tDecision_Logic\tModel_Type\tEvidence_Level\tTarget_KO\tTarget_Score\tTarget_Score_Type\tStrongest_Competitor\tCompetitor_Score\tDelta_Target\tStrongest_Negative\tNegative_Type\tNegative_Score\tDelta_Negative\tThreshold_Status\tLength_Status\tCatalytic_Motif_Status\tContext_Status\n")
        
        def process_and_write(prot):
            score_ann = "Hypothetical protein"
            evalue_ann = "Hypothetical protein"
            best_evalue = "-"
            target_ann = "-"
            status = "No_Hit"
            similar = []
            model = 'none'
            evidence = 'candidate'
            target_ko = '-'
            target_score = None
            target_score_type = '-'
            competitor = '-'
            competitor_score = None
            strongest_negative = '-'
            negative_score = None
            delta_target = None
            delta_negative = None
            negative_kind = '-'
            threshold_status = 'not_applicable'
            length_status = 'not_checked'
            motif_status = 'not_checked'
            context_status = 'not_checked'

            # 从预分组字典中直接获取该蛋白的所有 Hits (O(1) 复杂度)
            hits = unified_groups.get(prot, pd.DataFrame())
            
            score_hits = []
            evalue_hits = []
            target_hits = []
            
            if not hits.empty:
                for _, row in hits.iterrows():
                    q = row['Raw_Query']
                    
                    # 1. 检查是否满足 Score 官方限制阈值
                    is_score_pass = False
                    meta = ko_thresholds.get(str(q).upper(), {})
                    if meta.get('threshold') is not None:
                        score_value = row['Domain_Score'] if meta.get('score_type') == 'domain' else row['Full_Score']
                        if score_value is not None and score_value >= meta['threshold']:
                            is_score_pass = True
                            
                    # 2. 检查是否满足 E-value 阈值
                    is_evalue_pass = (row['E_value'] <= evalue_cutoff)
                    
                    # 3. 检查是否为指定的目标 Profile 且 E-value 达标
                    is_target_pass = (target_set and q in target_set and row['E_value'] <= evalue_cutoff)
                    
                    row_dict = row.to_dict()
                    if is_score_pass:
                        row_dict['Pass_Type'] = 'threshold'
                        score_hits.append(row_dict)
                    if is_evalue_pass:
                        row_dict['Pass_Type'] = 'evalue'
                        evalue_hits.append(row_dict)
                    if is_target_pass:
                        target_hits.append(row_dict)
            
            # 裁决结果
            has_score = len(score_hits) > 0
            has_evalue = len(evalue_hits) > 0
            
            if has_score:
                score_ann = score_hits[0]['HMM_Match']
            if has_evalue:
                evalue_ann = evalue_hits[0]['HMM_Match']
                best_evalue = str(evalue_hits[0]['E_value'])
            if target_hits:
                # 因为 hits 本身是按 E-value 从小到大排序的，所以 target_hits 中的第一个即为 E-value 最小的
                target_ann = target_hits[0]['HMM_Match']

            # New evidence arbitration. KOfam and HYdDB are deliberately
            # scored in separate namespaces; a raw maximum across databases
            # is never used as a functional call.
            if not hits.empty:
                typed = []
                for _, row in hits.iterrows():
                    q = row['Raw_Query']
                    typed.append((q, row, model_type(q, target_set)))
                kofam = [(q, r) for q, r, t in typed if t == 'KOfam_KO']
                hyd = [(q, r) for q, r, t in typed if t == 'HYdDB_subtype']
                if kofam:
                    model = 'KOfam_KO'
                    passing = []
                    for q, r in kofam:
                        meta = ko_thresholds.get(str(q).upper(), {})
                        stype = meta.get('score_type', '')
                        score = r['Domain_Score'] if stype == 'domain' else r['Full_Score']
                        if meta.get('threshold') is not None and score is not None and score >= meta['threshold']:
                            passing.append((q, r, score, stype))
                    if passing:
                        passing.sort(key=lambda x: x[2], reverse=True)
                        target_ko, r, target_score, target_score_type = passing[0]
                        threshold_status = 'pass'
                        target_ann = r['HMM_Match']
                        score_ann = r['HMM_Match']
                        competitors = []
                        for q, rr in kofam:
                            if q == target_ko:
                                continue
                            cmeta = ko_thresholds.get(str(q).upper(), {})
                            cscore = rr['Domain_Score'] if target_score_type == 'domain' else rr['Full_Score']
                            if cscore is not None:
                                competitors.append((q, rr, cscore))
                        if competitors:
                            competitor, _, competitor_score = max(competitors, key=lambda x: x[2])
                            delta_target = target_score - competitor_score
                        neg_defs = target_negatives.get(str(target_ko).upper(), [])
                        neg_hits = [(q, rr, next(n for n in neg_defs if q == n['ko'])) for q, rr in kofam if any(q == n['ko'] for n in neg_defs)]
                        if neg_hits:
                            strongest_negative, nr, ndef = max(
                                neg_hits,
                                key=lambda x: (x[1]['Domain_Score'] if target_score_type == 'domain' else x[1]['Full_Score']) or -math.inf)
                            negative_kind = ndef['type'] or 'unspecified'
                            negative_score = nr['Domain_Score'] if target_score_type == 'domain' else nr['Full_Score']
                            delta_negative = target_score - negative_score
                        hard_negative = negative_kind == 'non_target_hard_negative' or 'hard' in negative_kind
                        soft_negative = negative_kind == 'non_target_soft_neighbor' or 'soft' in negative_kind
                        if hard_negative and delta_negative is not None and delta_negative <= 0:
                            evidence = 'rejected'
                            status = 'Rejected_Negative'
                        elif delta_target is not None and delta_target <= 0:
                            evidence = 'rejected'
                            status = 'Rejected_Competitor'
                        elif competitor_score is not None and delta_target <= abs(target_score) * competition_fraction:
                            evidence = 'medium'
                            status = 'Ambiguous_Competitor'
                        elif negative_score is not None and delta_negative <= abs(target_score) * competition_fraction:
                            evidence = 'medium'
                            status = 'Candidate_Soft_Negative' if soft_negative else 'Ambiguous_Negative'
                        else:
                            # Motif/context checks are not available in the
                            # legacy tblout interface, so never overstate a
                            # call as high confidence here.
                            evidence = 'medium'
                    else:
                        threshold_status = 'fail_or_missing'
                        evidence = 'candidate'
                        status = 'Below_KOfam_Threshold'
                elif hyd:
                    model = 'HYdDB_subtype'
                    ranked = sorted(hyd, key=lambda x: x[1]['Full_Score'] or -math.inf, reverse=True)
                    target_ko, r = ranked[0]
                    target_score = r['Full_Score']
                    target_score_type = 'full'
                    target_ann = r['HMM_Match']
                    score_ann = r['HMM_Match']
                    if len(ranked) > 1:
                        competitor, cr = ranked[1]
                        competitor_score = cr['Full_Score']
                        delta_target = target_score - competitor_score
                        evidence = 'medium' if delta_target <= abs(target_score) * competition_fraction else 'candidate'
                        status = 'Subtype_Ambiguous' if evidence == 'medium' else 'Subtype_Candidate'
                    else:
                        evidence = 'candidate'
                        status = 'Subtype_Candidate'
                
            # Preserve the type-specific arbitration status above. The legacy
            # status labels are only used for non-KOfam/non-HYdDB profiles.
            if model == 'none':
                if has_score and has_evalue:
                    status = "Both_Hits"
                elif has_score:
                    status = "Score_Only"
                elif has_evalue:
                    status = "Evalue_Only"
                
            # 处理相似比对 (以最广的 evalue_hits 为基准计算类似 hits)
            if has_evalue:
                best_e = max(evalue_hits[0]['E_value'], 1e-300)
                log_e = abs(math.log10(best_e))
                dynamic_ratio = 10 ** (log_e * log_factor)
                
                similar = []
                for h in evalue_hits[1:]:
                    if h['Pass_Type'] == 'threshold':
                        similar.append(f"{h['HMM_Match']}(threshold)")
                    else:
                        if h['E_value'] <= evalue_hits[0]['E_value'] * dynamic_ratio:
                            similar.append(f"{h['HMM_Match']}(evalue)")

            pfam_info = pfam_dict.get(prot, "-")
            if prot in fasta_stats and target_score is not None:
                expected = ko_thresholds.get(str(target_ko).upper(), {}).get('mlen')
                length_status = 'ok' if expected and fasta_stats[prot]['length'] >= 0.5 * expected else 'short_or_unknown'
                if fasta_stats[prot]['internal_stop']:
                    length_status = 'internal_stop'
                    evidence = 'candidate' if evidence == 'high' else evidence
            if pfam_info != "-":
                if score_ann == "Hypothetical protein":
                    score_ann = f"Hypothetical protein (contains {pfam_info} Pfam domain)"
                if evalue_ann == "Hypothetical protein":
                    evalue_ann = f"Hypothetical protein (contains {pfam_info} Pfam domain)"
                if status == "No_Hit":
                    status = "Domain_Assisted"

            sim_str = ",".join(similar) if similar else "-"
            def fmt(v): return '-' if v is None else str(v)
            out.write("\t".join(map(str, [prot, score_ann, evalue_ann, best_evalue, target_ann, pfam_info, sim_str, status,
                model, evidence, target_ko, fmt(target_score), target_score_type, competitor, fmt(competitor_score),
                fmt(delta_target), strongest_negative, negative_kind, fmt(negative_score), fmt(delta_negative), threshold_status,
                length_status, motif_status, context_status])) + "\n")

        # 第一阶段：流式读取 FASTA 文件
        if faa_path and os.path.exists(faa_path):
            with open(faa_path) as f:
                for line in f:
                    if line.startswith(">"):
                        prot = line[1:].split()[0]
                        process_and_write(prot)
                        processed_proteins.add(prot)
            print(f"Processed {len(processed_proteins)} protein IDs streamed from FASTA.", file=sys.stderr)

        # 第二阶段：处理遗留蛋白
        all_hits = set(unified_groups.keys()) | set(pfam_dict.keys())
        leftovers = all_hits - processed_proteins
        if leftovers:
            print(f"Processing {len(leftovers)} leftover proteins not in FASTA...", file=sys.stderr)
            for prot in sorted(list(leftovers)):
                process_and_write(prot)

    print(f"Done.", file=sys.stderr)

def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--unified", required=True)
    p.add_argument("--pfam", required=True)
    p.add_argument("--ko_list", required=True)
    p.add_argument("--faa", required=False, help="Protein FASTA file to ensure all proteins are included")
    p.add_argument("--output", required=True)
    p.add_argument("--evalue_cutoff", type=float, default=1e-5)
    p.add_argument("--log_factor", type=float, default=0.2)
    p.add_argument("--target_profiles", default="", help="Comma-separated list of target profiles to find the minimum E-value for")
    p.add_argument("--negative_map", default="", help="Target-specific target_to_negative.tsv")
    p.add_argument("--competition_fraction", type=float, default=0.10,
                   help="Near-tie fraction; calibrate with validation data")
    args = p.parse_args()

    # 解析目标 Profile 集合
    import re
    target_set = set()
    if args.target_profiles:
        items = re.split(r'[\s,]+', args.target_profiles)
        target_set = {x.strip() for x in items if x.strip()}

    # 1. 加载 KO 阈值
    ko_thresholds = load_ko_thresholds(args.ko_list)
    target_negatives = load_target_negatives(args.negative_map)

    # 2. 运行裁决
    resolve_single_mag(
        unified_path=args.unified,
        pfam_path=args.pfam,
        ko_thresholds=ko_thresholds,
        faa_path=args.faa,
        output_path=args.output,
        evalue_cutoff=args.evalue_cutoff,
        log_factor=args.log_factor,
        target_set=target_set
        ,target_negatives=target_negatives,
        competition_fraction=args.competition_fraction
    )

if __name__ == "__main__":
    main()
    import os
    os._exit(0)
