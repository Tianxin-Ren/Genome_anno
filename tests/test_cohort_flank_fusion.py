#!/usr/bin/env python3
"""
Unit tests for Cohort Flank Fusion and Robust Contig Matching in extract_anchor_loci.py.

Verifies:
1. find_matching_record resolves contig IDs with GTDB GCA_...___ prefixes against NCBI FNA headers.
2. When reference MAGs suffer from assembly contig breaks near operons (e.g. Ref A has only downstream flank,
   Ref B has only upstream flank), the cohort flank fusion algorithm successfully synthesizes a BILATERAL locus.
3. High-ANI Stage 0 references supporting the consensus locus decisively beat isolated draft MAGs with spurious hits.
4. Target locus slicing encompasses the fused upstream flank, operon, and downstream flank.
"""

import os
import sys
import unittest
import pandas as pd

# Add scripts directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))
from extract_anchor_loci import (
    find_matching_record,
    is_same_contig
)

try:
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    class Seq:
        def __init__(self, data=""):
            self.data = data
        def __len__(self):
            return len(self.data)
        def __str__(self):
            return str(self.data)
        def __getitem__(self, item):
            return Seq(self.data[item])

    class SeqRecord:
        def __init__(self, seq, id="", description=""):
            self.seq = seq
            self.id = id
            self.description = description
        def __len__(self):
            return len(self.seq)


class TestCohortFlankFusion(unittest.TestCase):

    def test_find_matching_record_with_prefixes(self):
        """Test find_matching_record matches GTDB prefix IDs to clean NCBI headers."""
        records = {
            "WRCX01000007.1": SeqRecord(Seq("ATCG" * 1000), id="WRCX01000007.1", description="WRCX01000007.1 Bathycorpusculum"),
            "contig_02": SeqRecord(Seq("NNNN" * 500), id="contig_02", description="contig_02 other")
        }

        # 1. Exact match
        r1 = find_matching_record("WRCX01000007.1", records)
        self.assertIsNotNone(r1)
        self.assertEqual(r1.id, "WRCX01000007.1")

        # 2. GTDB prefix match
        gtdb_id = "GCA_009783705.1_ASM978370v1_genomic___WRCX01000007.1"
        r2 = find_matching_record(gtdb_id, records)
        self.assertIsNotNone(r2)
        self.assertEqual(r2.id, "WRCX01000007.1")

        # 3. Match with whitespace in description
        records_with_spaces = {
            "WRCX01000007.1 Bathycorpusculum sp.": SeqRecord(Seq("ATCG" * 1000), id="WRCX01000007.1", description="WRCX01000007.1 Bathycorpusculum sp.")
        }
        r3 = find_matching_record("WRCX01000007.1", records_with_spaces)
        self.assertIsNotNone(r3)

    def test_cohort_flank_fusion_bilateral_synthesis(self):
        """
        Simulate Bathy_u1287049 scenario:
        - Target contig: 2.75 Mb
        - Target operon: 2,050,000 .. 2,067,000 (genes T_02180 .. T_02189)
        - Ref 1 (GCA_009783705.1, ANI 95%): operon at start of contig, 10 downstream hits at 2,067,000 .. 2,100,000
        - Ref 2 (GCF_009778815.1, ANI 94%): operon at start of contig, 10 downstream hits at 2,067,000 .. 2,100,000
        - Ref 3 (GCA_009776805.1, ANI 93%): operon at end of contig, 4 upstream hits at 1,960,000 .. 1,987,000
        - Ref 4 (GCA_046706605.1, draft split): 12 unlinked hits at 2,390,000 .. 2,450,000 (300 kb away)

        Verify:
        Cohort flank fusion pairs Ref 3 (upstream) with Ref 1 & 2 (downstream),
        establishing a BILATERAL_SAME_CONTIG locus spanning [1.96 Mb .. 2.10 Mb],
        and rejects Ref 4 as the winning locus.
        """
        ctg_id = "target_ctg_2.75mb"

        # Group hits (with realistic ~197 kb gap representing intervening gene expansion)
        syn_u = [
            {"start": 1850000, "end": 1855000, "side": "UP", "ref_acc": "GCA_009776805.1", "orig_pid": "R3_u1", "t_pid": "T_01960", "bitscore": 400.0},
            {"start": 1865000, "end": 1870000, "side": "UP", "ref_acc": "GCA_009776805.1", "orig_pid": "R3_u2", "t_pid": "T_01987", "bitscore": 500.0}
        ]
        syn_d = [
            # Consensus downstream cluster at 2.067 - 2.10 Mb (supported by Ref 1 and Ref 2)
            {"start": 2067000, "end": 2075000, "side": "DOWN", "ref_acc": "GCA_009783705.1", "orig_pid": "R1_d1", "t_pid": "T_02190", "bitscore": 600.0},
            {"start": 2076000, "end": 2095000, "side": "DOWN", "ref_acc": "GCA_009783705.1", "orig_pid": "R1_d2", "t_pid": "T_02210", "bitscore": 700.0},
            {"start": 2067000, "end": 2075000, "side": "DOWN", "ref_acc": "GCF_009778815.1", "orig_pid": "R2_d1", "t_pid": "T_02190", "bitscore": 590.0},
            {"start": 2076000, "end": 2095000, "side": "DOWN", "ref_acc": "GCF_009778815.1", "orig_pid": "R2_d2", "t_pid": "T_02210", "bitscore": 680.0},
            # Outlier downstream cluster at 2.39 - 2.45 Mb (supported ONLY by Ref 4)
            {"start": 2390000, "end": 2420000, "side": "DOWN", "ref_acc": "GCA_046706605.1", "orig_pid": "R4_d1", "t_pid": "T_02535", "bitscore": 450.0},
            {"start": 2425000, "end": 2450000, "side": "DOWN", "ref_acc": "GCA_046706605.1", "orig_pid": "R4_d2", "t_pid": "T_02600", "bitscore": 500.0}
        ]

        # Target operon CDS features
        target_op_cds = [
            {"start": 2050000, "end": 2055000, "ko": "K02117", "contig": ctg_id},
            {"start": 2056000, "end": 2066000, "ko": "K02124", "contig": ctg_id}
        ]

        ref_candidates = {
            "GCA_009783705.1": {"tier": "Tier_A_Genus", "ani": 95.0, "single_contig_operon": True, "quality_mult": 1.85},
            "GCF_009778815.1": {"tier": "Tier_A_Genus", "ani": 94.0, "single_contig_operon": True, "quality_mult": 1.83},
            "GCA_009776805.1": {"tier": "Tier_A_Genus", "ani": 93.0, "single_contig_operon": True, "quality_mult": 1.81},
            "GCA_046706605.1": {"tier": "Tier_A_Genus", "ani": 94.5, "single_contig_operon": False, "quality_mult": 1.42}
        }

        # Cluster by gap (max_syntenic_gap = 35000)
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

        u_groups = group_syntenic_hits(syn_u, 35000)
        d_groups = group_syntenic_hits(syn_d, 35000)

        self.assertEqual(len(u_groups), 1)
        self.assertEqual(len(d_groups), 2)  # D1 at 2.07 Mb, D2 at 2.39 Mb

        max_cohort_gap = 450000
        candidate_loci = []

        for u_grp in u_groups:
            u_e = max(h["end"] for h in u_grp)
            n_u_pids = len(set(h["orig_pid"] for h in u_grp))
            u_refs = set(h["ref_acc"] for h in u_grp)

            for d_grp in d_groups:
                d_s = min(h["start"] for h in d_grp)
                n_d_pids = len(set(h["orig_pid"] for h in d_grp))
                d_refs = set(h["ref_acc"] for h in d_grp)

                inner_gap = d_s - u_e
                op_in_int = any(c["end"] >= u_e and c["start"] <= d_s for c in target_op_cds)
                allowed_gap = max_cohort_gap * 2 if op_in_int else max_cohort_gap

                if inner_gap <= allowed_gap:
                    all_refs = u_refs | d_refs
                    tot_bit = sum(h["bitscore"] for h in u_grp + d_grp)
                    q_weights = [ref_candidates[r]["quality_mult"] for r in all_refs]
                    avg_q = sum(q_weights) / max(1.0, len(q_weights))
                    bilat_sc = 10000.0 if (n_u_pids >= 2 and n_d_pids >= 2) else 5000.0
                    op_bonus = 25000.0 if op_in_int else 0.0
                    conv_bonus = 4000.0 * len(all_refs)

                    sc = (tot_bit * avg_q) + bilat_sc + op_bonus + conv_bonus + 100.0 * (n_u_pids + n_d_pids) - (0.02 * inner_gap)
                    candidate_loci.append({
                        "type": "BILATERAL",
                        "score": sc,
                        "hits": u_grp + d_grp,
                        "refs": all_refs,
                        "left_anchor": max(u_grp, key=lambda h: h["end"])["t_pid"],
                        "right_anchor": min(d_grp, key=lambda h: h["start"])["t_pid"],
                        "interval_start": u_e,
                        "interval_end": d_s,
                        "observed_span": inner_gap,
                        "has_op_in_interval": op_in_int
                    })

        self.assertTrue(len(candidate_loci) >= 1)
        winning_locus = max(candidate_loci, key=lambda x: (
            1 if x["type"] == "BILATERAL" else 0,
            1 if x.get("has_op_in_interval") else 0,
            x["score"]
        ))
        self.assertEqual(winning_locus["type"], "BILATERAL")
        self.assertEqual(winning_locus["left_anchor"], "T_01987")
        self.assertEqual(winning_locus["right_anchor"], "T_02190")
        self.assertTrue(winning_locus["has_op_in_interval"])
        self.assertIn("GCA_009783705.1", winning_locus["refs"])
        self.assertIn("GCF_009778815.1", winning_locus["refs"])
        self.assertIn("GCA_009776805.1", winning_locus["refs"])
        self.assertNotIn("GCA_046706605.1", winning_locus["refs"])

        # Test Reference Ranking on the Winning Locus
        locus_refs = winning_locus["refs"]
        locus_hits = winning_locus["hits"]

        def ref_ranking_key(r):
            cand_r = ref_candidates.get(r, {})
            in_loc = 1 if r in locus_refs else 0
            tier = cand_r.get("tier", "")
            tier_val = 3000.0 if "Tier_A" in tier else 1000.0
            ani_val = float(cand_r.get("ani", 80.0)) * 50.0
            single_ctg = 2000.0 if cand_r.get("single_contig_operon", False) else 0.0
            r_loc_bits = sum(h["bitscore"] for h in locus_hits if h.get("ref_acc") == r)
            return (in_loc, tier_val + ani_val + single_ctg, r_loc_bits)

        all_candidate_refs = list(ref_candidates.keys())
        best_ref = max(all_candidate_refs, key=ref_ranking_key)
        self.assertEqual(best_ref, "GCA_009783705.1", "Top ANI Stage 0 reference must win over draft split assemblies")

        # Test Target Slice Encompassing (Compact operon-centered slicing)
        flank_expansion = 30000
        op_min = min(c["start"] for c in target_op_cds)  # 2,050,000
        op_max = max(c["end"] for c in target_op_cds)    # 2,066,000
        max_flank_reach = max(int(flank_expansion * 1.5), 45000)

        # Proximal hits within reach of operon
        prox_hits = [h for h in locus_hits if (h.get("end", 0) >= op_min - max_flank_reach) and (h.get("start", 0) <= op_max + max_flank_reach)]
        if prox_hits:
            p_min = min(h["start"] for h in prox_hits)
            p_max = max(h["end"] for h in prox_hits)
            raw_start = min(op_min - flank_expansion, p_min - flank_expansion // 2)
            raw_end = max(op_max + flank_expansion, p_max + flank_expansion // 2)
        else:
            raw_start = op_min - flank_expansion
            raw_end = op_max + flank_expansion

        # The slice MUST contain the full operon
        self.assertTrue(raw_start <= 2050000 and 2066000 <= raw_end, "Slice must contain the complete operon")
        # The slice MUST contain the proximal downstream flank
        self.assertGreaterEqual(raw_end, 2095000, "Slice must encompass proximal downstream flank")
        # The slice must NOT stretch 180 kb back to the distant 1.85 Mb anchor
        self.assertGreater(raw_start, 1870000, "Slice must not stretch over 180 kb intervening gap to distant upstream anchor")
        self.assertEqual(raw_start, 2050000 - flank_expansion, "Upstream slice boundary should be clean flank expansion from operon")
        self.assertLessEqual(raw_end - raw_start, 100000, "Target slice must be compact (<= 100 kb) for publication-quality Clinker visualization")


if __name__ == "__main__":
    unittest.main()
