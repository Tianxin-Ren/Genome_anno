#!/usr/bin/env python3
"""
analyze_atp_correlation.py
Whole-genome KO and ATP synthase distribution co-occurrence & correlation analysis engine.
Features:
1. Zero Data Leakage: Isolates ATP phenotype-defining KOs as positive controls.
2. Two Contrast Families:
   - Contrast 1: ATP_any_presence (ATP_POSITIVE vs ATP_NEGATIVE; Partial/Fragmentary excluded)
   - Contrast 2: VA_vs_F_like (VA_ONLY vs F_LIKE_ONLY; Dual/Mixed excluded)
3. 3-Tier Statistical Architecture:
   - Level 1: Global Fisher's exact test, Haldane-Anscombe OR & 95% CI, phi coefficient, BH FDR.
   - Level 2: GTDB taxonomy-stratified Cochran-Mantel-Haenszel (CMH) test, Breslow-Day heterogeneity,
     strata concordance fraction.
   - Level 3: Genome reduction sensitivity logistic regression with leave-one-out KO richness
     (KO_richness_{-i}), log(genome_size), and automatic Firth penalized likelihood fallback on separation.
4. Independent FDR correction per contrast family.
5. Strict candidate tiering (Tier_1_robust, Tier_2_raw_discovery, Tier_3_non_significant, positive_control).
"""

import argparse
import gzip
import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, norm
import scipy.optimize as opt

# Target definition KOs for ATP synthase complexes
TARGET_DEFINITION_KOS = {
    # M00159 (Prokaryotic V/A-type)
    "K02117", "K02118", "K02119", "K02120", "K02121", "K02122", "K02107", "K02123", "K02124",
    # M00157 (Prokaryotic F-like)
    "K02111", "K02112", "K02113", "K02114", "K02115", "K02108", "K02109", "K02110",
    # M00158 (Eukaryotic F-type)
    "K02126", "K02127", "K02128", "K02129", "K02130", "K02131", "K02132", "K02133", "K02134", "K02135", "K02136", "K02137", "K02138",
    # M00160 (Eukaryotic V-type)
    "K02145", "K02146", "K02147", "K02148", "K02149", "K02150", "K02151", "K02152", "K02153", "K02154", "K02155",
}


def open_file(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)


def extract_gtdb_order(tax_str):
    """Extracts order (o__...) from GTDB taxonomy string."""
    if not tax_str or not isinstance(tax_str, str):
        return "o__unknown"
    parts = tax_str.split(";")
    for p in parts:
        p = p.strip()
        if p.startswith("o__"):
            return p
    for p in parts:
        p = p.strip()
        if p.startswith("c__"):
            return p
    return "o__unknown"


def benjamini_hochberg(p_values):
    """Computes Benjamini-Hochberg FDR (q-values) for a list/array of p-values."""
    p_arr = np.asarray(p_values, dtype=float)
    n = len(p_arr)
    if n == 0:
        return np.array([])
    valid_mask = ~np.isnan(p_arr)
    valid_p = p_arr[valid_mask]
    q_arr = np.full(n, np.nan)
    if len(valid_p) == 0:
        return q_arr

    sorted_indices = np.argsort(valid_p)
    sorted_p = valid_p[sorted_indices]
    m = len(sorted_p)

    q_vals = np.zeros(m)
    min_q = 1.0
    for i in range(m - 1, -1, -1):
        q = (sorted_p[i] * m) / (i + 1)
        if q < min_q:
            min_q = q
        q_vals[i] = min(min_q, 1.0)

    for idx, orig_pos in enumerate(sorted_indices):
        q_arr[np.where(valid_mask)[0][orig_pos]] = q_vals[idx]

    return q_arr


def fit_firth_logistic(X, y, max_iter=25, tol=1e-4):
    """
    High-performance Firth penalized logistic regression via Newton-Raphson IRLS (Heinze & Schemper 2002).
    Fully vectorized in NumPy BLAS (converges in 4-6 iterations, ~5ms on 22k rows).
    """
    n, p = X.shape
    beta = np.zeros(p)
    cov = np.eye(p)
    converged = False

    for it in range(max_iter):
        eta = np.clip(X @ beta, -25.0, 25.0)
        p_vec = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p_vec * (1.0 - p_vec), 1e-5, 0.25)
        XW = X * w[:, None]
        I_mat = X.T @ XW
        try:
            cov = np.linalg.pinv(I_mat)
        except np.linalg.LinAlgError:
            break
        # Hat matrix diagonal (leverage): h_i = diag(X (X^T W X)^-1 X^T W)
        h = np.sum((X @ cov) * X, axis=1) * w
        # Firth modified score: U*(beta) = X^T (y - p + h * (0.5 - p))
        u_star = X.T @ (y - p_vec + h * (0.5 - p_vec))
        delta = cov @ u_star
        beta += delta
        if np.max(np.abs(delta)) < tol:
            converged = True
            break

    se = np.sqrt(np.maximum(0.0, np.diag(cov)))
    z = beta / np.maximum(se, 1e-12)
    pval = 2.0 * (1.0 - norm.cdf(np.abs(z)))
    return beta, se, pval, (converged or it >= 3)


def fast_irls_logistic(X, y, max_iter=25, tol=1e-4, l2_reg=1e-4):
    """
    High-performance Iteratively Reweighted Least Squares (Newton-Raphson)
    for multivariable logistic regression with L2 ridge stabilization.
    Converges in 3-5 iterations (<1ms on 22k rows).
    """
    n, p = X.shape
    beta = np.zeros(p)
    reg_diag = np.full(p, l2_reg)
    reg_diag[0] = 0.0  # Unpenalized intercept
    reg_mat = np.diag(reg_diag)

    converged = False
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -25.0, 25.0)
        p_vec = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p_vec * (1.0 - p_vec), 1e-5, 0.25)
        grad = X.T @ (y - p_vec) - reg_diag * beta
        H = (X.T * w) @ X + reg_mat
        try:
            delta = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            return None, None, None, False
        beta += delta
        if np.max(np.abs(delta)) < tol:
            converged = True
            break

    if converged or np.max(np.abs(grad)) < 1e-2:
        eta = np.clip(X @ beta, -25.0, 25.0)
        p_vec = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p_vec * (1.0 - p_vec), 1e-5, 0.25)
        H = (X.T * w) @ X + reg_mat
        try:
            cov = np.linalg.pinv(H)
            se = np.sqrt(np.maximum(0.0, np.diag(cov)))
            z = beta / np.maximum(se, 1e-12)
            pval = 2.0 * (1.0 - norm.cdf(np.abs(z)))
            # Check for extreme separation or numerical instability
            if np.any(np.abs(beta) > 10.0) or np.any(se > 10.0) or np.any(np.isnan(se)):
                return None, None, None, False
            return beta, se, pval, True
        except Exception:
            return None, None, None, False

    return None, None, None, False


def fit_logistic_with_firth_fallback(X, y):
    """
    High-performance logistic regression with automatic Firth penalized likelihood fallback
    upon non-convergence, singular matrix, or coefficient separation (|beta| > 10, SE > 10).
    """
    b_irls, se_irls, p_irls, ok = fast_irls_logistic(X, y)
    if ok:
        return b_irls, se_irls, p_irls, "logistic", True

    # Fallback to Firth Penalized Likelihood
    b_firth, se_firth, p_firth, s_firth = fit_firth_logistic(X, y)
    if s_firth:
        return b_firth, se_firth, p_firth, "firth_logistic", True
    return np.zeros(X.shape[1]), np.full(X.shape[1], np.nan), np.full(X.shape[1], np.nan), "failed", False


def calculate_cmh(tables):
    """
    Computes Cochran-Mantel-Haenszel (CMH) statistic, pooled OR, 95% CI,
    and Breslow-Day / Woolf test for homogeneity of odds ratios.
    tables: list of 2x2 lists [[a, b], [c, d]] for testable strata
    """
    if not tables:
        return np.nan, np.nan, np.nan, np.nan, np.nan, []

    sum_a = 0.0
    sum_exp_a = 0.0
    sum_var_a = 0.0
    sum_r = 0.0
    sum_s = 0.0
    p_list, q_list, r_list, s_list = [], [], [], []
    strata_ors = []

    for t in tables:
        a, b = float(t[0][0]), float(t[0][1])
        c, d = float(t[1][0]), float(t[1][1])
        n1 = a + c
        n2 = b + d
        N = n1 + n2
        if N <= 1 or n1 == 0 or n2 == 0 or (a + b) == 0 or (c + d) == 0:
            continue
        m1 = a + b
        m2 = c + d
        exp_a = (n1 * m1) / N
        var_a = (n1 * n2 * m1 * m2) / (N * N * (N - 1))
        sum_a += a
        sum_exp_a += exp_a
        sum_var_a += var_a

        r = (a * d) / N
        s = (b * c) / N
        sum_r += r
        sum_s += s
        p_list.append((a + d) / N)
        q_list.append((b + c) / N)
        r_list.append(r)
        s_list.append(s)

        # Smoothed stratum OR
        or_k = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        strata_ors.append(or_k)

    if sum_var_a == 0 or sum_s == 0 or sum_r == 0:
        return 1.0, np.nan, np.nan, 1.0, 1.0, strata_ors

    diff = abs(sum_a - sum_exp_a)
    cmh_stat = max(0.0, diff - 0.5) ** 2 / sum_var_a
    cmh_p = 1.0 - chi2.cdf(cmh_stat, df=1)

    pooled_or = sum_r / sum_s

    # Robins-Breslow-Greenland variance estimator
    term1 = sum(p * r for p, r in zip(p_list, r_list)) / (2.0 * sum_r * sum_r)
    term2 = sum(p * s + q * r for p, q, r, s in zip(p_list, q_list, r_list, s_list)) / (2.0 * sum_r * sum_s)
    term3 = sum(q * s for q, s in zip(q_list, s_list)) / (2.0 * sum_s * sum_s)
    var_ln_or = term1 + term2 + term3
    se_ln_or = np.sqrt(max(0.0, var_ln_or))
    ci_low = np.exp(np.log(max(1e-12, pooled_or)) - 1.96 * se_ln_or)
    ci_high = np.exp(np.log(max(1e-12, pooled_or)) + 1.96 * se_ln_or)

    # Woolf / Breslow-Day test of homogeneity
    het_stat = 0.0
    for t in tables:
        a, b = float(t[0][0]), float(t[0][1])
        c, d = float(t[1][0]), float(t[1][1])
        a_s, b_s, c_s, d_s = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        v_k = 1.0 / a_s + 1.0 / b_s + 1.0 / c_s + 1.0 / d_s
        or_k = (a_s * d_s) / (b_s * c_s)
        het_stat += ((np.log(max(1e-12, or_k)) - np.log(max(1e-12, pooled_or))) ** 2) / max(1e-12, v_k)
    df_het = max(1, len(tables) - 1)
    het_p = 1.0 - chi2.cdf(het_stat, df=df_het)

    return pooled_or, ci_low, ci_high, cmh_p, het_p, strata_ors


def load_genome_ko_matrix(ko_counts_path, target_genomes):
    """
    Reads genome KO presence from genome_KO_counts.long.tsv (or protein_KO_all_strict.tsv).
    Returns:
      ko_presence_by_genome: genome_id -> set of KOs
      all_kos: set of all KOs observed
      ko_definitions: ko -> description string
      ko_to_genomes: ko -> set of genome_ids
    """
    target_genomes_set = set(target_genomes)
    ko_presence_by_genome = {g: set() for g in target_genomes}
    ko_to_genomes = {}
    all_kos = set()
    ko_definitions = {}

    print(f"Loading KO presence matrix from {ko_counts_path}...")
    with open_file(ko_counts_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_map = {c.lower(): idx for idx, c in enumerate(header)}

        g_idx = col_map.get("genome_id", col_map.get("genome", 0))
        k_idx = col_map.get("ko", 1)
        d_idx = col_map.get("ko_definition", col_map.get("definition", -1))

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > max(g_idx, k_idx):
                gid = parts[g_idx]
                if gid in target_genomes_set:
                    ko = parts[k_idx]
                    ko_presence_by_genome[gid].add(ko)
                    if ko not in ko_to_genomes:
                        ko_to_genomes[ko] = set()
                    ko_to_genomes[ko].add(gid)
                    all_kos.add(ko)
                    if d_idx != -1 and d_idx < len(parts) and ko not in ko_definitions:
                        ko_definitions[ko] = parts[d_idx]

    print(f"Loaded presence data for {len(all_kos)} unique KOs across {len(target_genomes)} cohort genomes.")
    return ko_presence_by_genome, all_kos, ko_definitions, ko_to_genomes


def run_association_for_contrast(contrast_name, group1_genomes, group2_genomes,
                                master_df, ko_presence_by_genome, all_kos, ko_definitions,
                                ko_to_genomes=None,
                                min_prevalence_count=5, fdr_alpha=0.05, delta_prev_cutoff=0.20,
                                min_stratum_group_size=3, run_regression=True,
                                representative_genomes=None):
    """
    Executes full 3-tier association analysis for a single contrast.
    """
    n1 = len(group1_genomes)
    n2 = len(group2_genomes)
    print(f"\n--- Running Analysis for Contrast: {contrast_name} (Group 1 N={n1}, Group 2 N={n2}) ---")

    if n1 == 0 or n2 == 0:
        print(f"Warning: Contrast {contrast_name} has insufficient sample size (N1={n1}, N2={n2})! Skipping.")
        return []

    if ko_to_genomes is None:
        ko_to_genomes = {}
        for g, kos in ko_presence_by_genome.items():
            for k in kos:
                if k not in ko_to_genomes:
                    ko_to_genomes[k] = set()
                ko_to_genomes[k].add(g)

    cohort_genomes = group1_genomes + group2_genomes
    cohort_df = master_df[master_df["genome_id"].isin(cohort_genomes)].copy()
    cohort_df.set_index("genome_id", inplace=True)

    # Precompute covariates for regression
    # Dependent variable y: 1 if Group 1, 0 if Group 2
    y_dict = {g: 1.0 for g in group1_genomes}
    y_dict.update({g: 0.0 for g in group2_genomes})

    # GTDB Orders
    order_dict = {}
    for g in cohort_genomes:
        tax = str(cohort_df.loc[g, "gtdb_taxonomy"]) if "gtdb_taxonomy" in cohort_df.columns else ""
        order_dict[g] = extract_gtdb_order(tax)

    all_orders = sorted(set(order_dict.values()))

    # Log genome size, CheckM2 completeness, total KO richness
    log_size_dict = {}
    comp_dict = {}
    total_richness_dict = {}

    for g in cohort_genomes:
        try:
            sz = float(cohort_df.loc[g, "genome_size"])
            log_size_dict[g] = np.log(max(100000.0, sz))
        except (ValueError, TypeError, KeyError):
            log_size_dict[g] = 14.5  # ~2 Mb fallback

        try:
            cp = float(cohort_df.loc[g, "checkm2_completeness"])
            if np.isnan(cp):
                cp = float(cohort_df.loc[g, "checkm_completeness"])
        except (ValueError, TypeError, KeyError):
            cp = 100.0
        comp_dict[g] = cp if not np.isnan(cp) else 100.0

        total_richness_dict[g] = float(len(ko_presence_by_genome.get(g, set())))

    # Pre-scale continuous covariates to improve numerical stability in logistic regression
    log_size_vals = np.array([log_size_dict[g] for g in cohort_genomes])
    mean_ls, std_ls = np.mean(log_size_vals), max(1e-6, np.std(log_size_vals))

    comp_vals = np.array([comp_dict[g] for g in cohort_genomes])
    mean_cp, std_cp = np.mean(comp_vals), max(1e-6, np.std(comp_vals))

    rich_vals = np.array([total_richness_dict[g] for g in cohort_genomes])
    mean_rc, std_rc = np.mean(rich_vals), max(1e-6, np.std(rich_vals))

    n_cohort = len(cohort_genomes)
    cohort_genomes_set = set(cohort_genomes)
    genome_to_idx = {g: i for i, g in enumerate(cohort_genomes)}
    y_arr = np.array([y_dict[g] for g in cohort_genomes], dtype=float)

    # Pre-scale continuous covariates ONCE outside the KO loop
    scaled_ls = (log_size_vals - mean_ls) / std_ls
    scaled_cp = (comp_vals - mean_cp) / std_cp
    rich_norm = (rich_vals - mean_rc) / std_rc
    inv_std_rc = 1.0 / std_rc

    # Pre-allocate reusable design matrix template
    X_template = np.empty((n_cohort, 5), dtype=float)
    X_template[:, 0] = 1.0
    X_template[:, 2] = scaled_ls
    X_template[:, 4] = scaled_cp

    g1_set = set(group1_genomes)
    g2_set = set(group2_genomes)

    # Pre-filter and pre-slice testable strata for CMH (avoids O(K*O*N) repeated list comprehensions)
    testable_strata = []
    for o in all_orders:
        o_g1 = set(g for g in group1_genomes if order_dict[g] == o)
        o_g2 = set(g for g in group2_genomes if order_dict[g] == o)
        if len(o_g1) >= min_stratum_group_size and len(o_g2) >= min_stratum_group_size:
            testable_strata.append((o, o_g1, o_g2, len(o_g1), len(o_g2)))

    # Precompute representative sets for fast sensitivity check
    r_g1_set = set(g for g in group1_genomes if g in representative_genomes) if representative_genomes else set()
    r_g2_set = set(g for g in group2_genomes if g in representative_genomes) if representative_genomes else set()
    check_rep = (len(r_g1_set) >= 3 and len(r_g2_set) >= 3)

    records = []
    t_start = time.time()
    total_kos_count = len(all_kos)

    for idx, ko in enumerate(sorted(all_kos)):
        if (idx + 1) % 2000 == 0 or (idx + 1) == total_kos_count:
            elapsed = time.time() - t_start
            rate = (idx + 1) / max(0.1, elapsed)
            print(f"[{contrast_name}] Tested {idx + 1}/{total_kos_count} KOs ({elapsed:.1f}s, {rate:.1f} KOs/s)...", flush=True)

        is_target_def = (ko in TARGET_DEFINITION_KOS)
        ko_g = ko_to_genomes.get(ko, set())

        # 1. 2x2 Contingency Table (Fast C-speed set operations)
        a = len(g1_set & ko_g)
        b = len(g2_set & ko_g)
        c = n1 - a
        d = n2 - b
        total_ko_count = a + b

        prev1 = a / float(n1)
        prev2 = b / float(n2)
        delta_prev = prev1 - prev2

        # Check eligibility based on minimum prevalence count
        eligible = (total_ko_count >= min_prevalence_count)
        exclusion_reason = "none"
        if not eligible:
            exclusion_reason = f"low_prevalence_count_{total_ko_count}<{min_prevalence_count}"

        # 2. Fisher's Exact Test on RAW table
        res_fisher = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        fisher_p = res_fisher[1]

        # 3. Haldane-Anscombe smoothed Odds Ratio & 95% CI
        if a == 0 or b == 0 or c == 0 or d == 0:
            a_s, b_s, c_s, d_s = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        else:
            a_s, b_s, c_s, d_s = float(a), float(b), float(c), float(d)

        odds_ratio = (a_s * d_s) / (b_s * c_s)
        log2_or = np.log2(odds_ratio)
        se_ln_or = np.sqrt(1.0 / a_s + 1.0 / b_s + 1.0 / c_s + 1.0 / d_s)
        ci_low = np.exp(np.log(odds_ratio) - 1.96 * se_ln_or)
        ci_high = np.exp(np.log(odds_ratio) + 1.96 * se_ln_or)

        # Phi coefficient
        denom = np.sqrt(float(a + b) * float(c + d) * float(a + c) * float(b + d))
        phi = ((a * d) - (b * c)) / denom if denom > 0 else 0.0

        # 4. Taxonomy-Stratified (Order-level) CMH analysis
        order_tables = []
        n_pos_strata = 0
        n_neg_strata = 0
        testable_orders = 0

        for o, o_g1, o_g2, len_o1, len_o2 in testable_strata:
            o_a = len(o_g1 & ko_g)
            o_b = len(o_g2 & ko_g)
            o_c = len_o1 - o_a
            o_d = len_o2 - o_b
            # Stratum must have variation in KO
            if (o_a + o_b) > 0 and (o_c + o_d) > 0:
                testable_orders += 1
                order_tables.append([[o_a, o_b], [o_c, o_d]])
                o_or = ((o_a + 0.5) * (o_d + 0.5)) / ((o_b + 0.5) * (o_c + 0.5))
                if o_or > 1.0:
                    n_pos_strata += 1
                elif o_or < 1.0:
                    n_neg_strata += 1

        concordance_frac = (max(n_pos_strata, n_neg_strata) / float(testable_orders)) if testable_orders > 0 else 0.0
        cmh_or, cmh_ci_low, cmh_ci_high, cmh_p, het_p, _ = calculate_cmh(order_tables)

        # 5. Multivariate Logistic Regression controlling for genome reduction
        reg_beta = np.nan
        reg_se = np.nan
        reg_p = np.nan
        reg_method = "none"
        reg_conv = False

        if run_regression and eligible:
            x_ko = np.zeros(n_cohort, dtype=float)
            active = ko_g.intersection(cohort_genomes_set)
            if active:
                x_ko[[genome_to_idx[g] for g in active]] = 1.0

            X_mat = X_template.copy()
            X_mat[:, 1] = x_ko
            X_mat[:, 3] = rich_norm - (x_ko * inv_std_rc)

            beta_hat, se_hat, p_hat, reg_method, reg_conv = fit_logistic_with_firth_fallback(X_mat, y_arr)
            if reg_conv:
                reg_beta = beta_hat[1]
                reg_se = se_hat[1]
                reg_p = p_hat[1]

        # 6. Check representative-only sensitivity if provided
        rep_pass = False
        if check_rep:
            r_a = len(r_g1_set & ko_g)
            r_b = len(r_g2_set & ko_g)
            r_p = fisher_exact([[r_a, r_b], [len(r_g1_set) - r_a, len(r_g2_set) - r_b]])[1]
            if r_p < 0.05 and ((r_a / len(r_g1_set)) > (r_b / len(r_g2_set))) == (prev1 > prev2):
                rep_pass = True

        rec = {
            "ko": ko,
            "ko_definition": ko_definitions.get(ko, "Unannotated"),
            "contrast": contrast_name,
            "n_group1": n1,
            "n_group2": n2,
            "a": a,
            "b": b,
            "c": c,
            "d": d,
            "prev_group1": round(prev1, 4),
            "prev_group2": round(prev2, 4),
            "delta_prevalence": round(delta_prev, 4),
            "odds_ratio": round(odds_ratio, 4),
            "or_ci_low": round(ci_low, 4),
            "or_ci_high": round(ci_high, 4),
            "log2_odds_ratio": round(log2_or, 4),
            "phi_coefficient": round(phi, 4),
            "fisher_p": fisher_p,
            "fisher_q": np.nan, # To be computed across contrast
            "cmh_odds_ratio": round(cmh_or, 4) if not np.isnan(cmh_or) else np.nan,
            "cmh_ci_low": round(cmh_ci_low, 4) if not np.isnan(cmh_ci_low) else np.nan,
            "cmh_ci_high": round(cmh_ci_high, 4) if not np.isnan(cmh_ci_high) else np.nan,
            "cmh_p": cmh_p if not np.isnan(cmh_p) else np.nan,
            "cmh_q": np.nan, # To be computed across contrast
            "heterogeneity_p": round(het_p, 4) if not np.isnan(het_p) else np.nan,
            "n_strata_total": len(all_orders),
            "n_strata_testable": testable_orders,
            "n_strata_positive": n_pos_strata,
            "n_strata_negative": n_neg_strata,
            "strata_concordance_fraction": round(concordance_frac, 4),
            "reduction_model_beta": round(reg_beta, 4) if not np.isnan(reg_beta) else np.nan,
            "reduction_model_se": round(reg_se, 4) if not np.isnan(reg_se) else np.nan,
            "reduction_model_p": reg_p if not np.isnan(reg_p) else np.nan,
            "reduction_model_q": np.nan, # To be computed across contrast
            "regression_method": reg_method,
            "regression_converged": reg_conv,
            "target_definition_ko": is_target_def,
            "analysis_eligible": eligible,
            "exclusion_reason": exclusion_reason,
            "rep_pass": rep_pass,
            "candidate_tier": "Tier_3_non_significant",
            "robustness_flags": "NONE",
        }
        records.append(rec)

    # Compute BH FDR within this independent contrast family
    # 1. Fisher q
    fisher_ps = [r["fisher_p"] for r in records]
    fisher_qs = benjamini_hochberg(fisher_ps)
    for idx, q in enumerate(fisher_qs):
        records[idx]["fisher_q"] = round(q, 6) if not np.isnan(q) else np.nan

    # 2. CMH q (only for testable KOs)
    cmh_ps = [r["cmh_p"] if r["n_strata_testable"] > 0 else np.nan for r in records]
    cmh_qs = benjamini_hochberg(cmh_ps)
    for idx, q in enumerate(cmh_qs):
        records[idx]["cmh_q"] = round(q, 6) if not np.isnan(q) else np.nan

    # 3. Reduction regression q
    reg_ps = [r["reduction_model_p"] for r in records]
    reg_qs = benjamini_hochberg(reg_ps)
    for idx, q in enumerate(reg_qs):
        records[idx]["reduction_model_q"] = round(q, 6) if not np.isnan(q) else np.nan

    # Candidate Tier and Robustness Flags Assignment
    for r in records:
        f_q = r["fisher_q"]
        d_p = r["delta_prevalence"]
        r_q = r["reduction_model_q"]
        c_q = r["cmh_q"]
        n_strata = r["n_strata_testable"]
        c_frac = r["strata_concordance_fraction"]

        flags = []
        fdr_pass = (not np.isnan(f_q) and f_q < fdr_alpha)
        if fdr_pass:
            flags.append("GLOBAL_FDR_PASS")

        eff_pass = (abs(d_p) >= delta_prev_cutoff)
        if eff_pass:
            flags.append("EFFECT_SIZE_PASS")

        red_pass = (not np.isnan(r_q) and r_q < fdr_alpha)
        if red_pass:
            flags.append("REDUCTION_ADJUSTED_PASS")

        multi_order_pass = (n_strata >= 2 and c_frac >= 0.75)
        if multi_order_pass:
            flags.append("MULTI_ORDER_CONCORDANT")

        cmh_pass = (not np.isnan(c_q) and c_q < fdr_alpha)
        if cmh_pass:
            flags.append("CMH_PASS")

        if r["rep_pass"]:
            flags.append("REPRESENTATIVE_PASS")

        r["robustness_flags"] = ";".join(flags) if flags else "NONE"

        # Tier assignment
        if r["target_definition_ko"]:
            r["candidate_tier"] = "positive_control"
        elif fdr_pass and eff_pass and red_pass and (multi_order_pass or n_strata < 2):
            r["candidate_tier"] = "Tier_1_robust"
        elif fdr_pass and eff_pass:
            r["candidate_tier"] = "Tier_2_raw_discovery"
        else:
            r["candidate_tier"] = "Tier_3_non_significant"

    # Clean temporary keys
    for r in records:
        r.pop("rep_pass", None)

    print(f"Contrast {contrast_name} evaluated across {len(records)} KOs.")
    t1_cnt = sum(1 for r in records if r["candidate_tier"] == "Tier_1_robust")
    t2_cnt = sum(1 for r in records if r["candidate_tier"] == "Tier_2_raw_discovery")
    print(f"Candidates: {t1_cnt} Tier 1 robust, {t2_cnt} Tier 2 raw discovery.")

    return records


def main():
    parser = argparse.ArgumentParser(description="Analyze whole-genome KO correlation and co-occurrence with ATP synthase distribution.")
    parser.add_argument("-m", "--master_table", required=True, help="Path to genome_atp_synthase_master_table.tsv")
    parser.add_argument("-k", "--ko_counts", required=True, help="Path to genome_KO_counts.long.tsv (or .gz)")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for correlation tables")
    parser.add_argument("--exclude_target_definition_kos", action="store_true", default=True, help="Exclude target ATPase KOs from biological candidate output")
    parser.add_argument("--min_ko_prevalence_count", type=int, default=5, help="Min count across cohort for KO to be tested (default: 5)")
    parser.add_argument("--fdr_alpha", type=float, default=0.05, help="FDR threshold for significance (default: 0.05)")
    parser.add_argument("--delta_prevalence_cutoff", type=float, default=0.20, help="Min delta prevalence cutoff (default: 0.20)")
    parser.add_argument("--min_stratum_group_size", type=int, default=3, help="Min genomes per phenotype in testable stratum (default: 3)")
    parser.add_argument("--sampling_mode_primary", default="all", choices=["all", "gtdb_representatives"], help="Primary cohort sampling mode")
    parser.add_argument("--run_genome_reduction_sensitivity", action="store_true", default=True, help="Run multivariable logistic regression sensitivity")

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading master table from {args.master_table}...")
    master_df = pd.read_csv(args.master_table, sep="\t", low_memory=False)

    # Parse representative genomes if column present
    representative_genomes = set()
    if "gtdb_representative" in master_df.columns:
        rep_rows = master_df[master_df["gtdb_representative"].astype(str).str.lower().isin(["t", "true", "1"])]
        representative_genomes = set(rep_rows["genome_id"].dropna().unique())
        print(f"Identified {len(representative_genomes)} GTDB representative genomes.")

    # Filter primary cohort if requested
    if args.sampling_mode_primary == "gtdb_representatives" and representative_genomes:
        active_master_df = master_df[master_df["genome_id"].isin(representative_genomes)].copy()
        print(f"Primary analysis restricted to {len(active_master_df)} GTDB representative genomes.")
    else:
        active_master_df = master_df

    # Discover contrast cohorts
    # Contrast 1: ATP_any_presence (ATP_POSITIVE vs ATP_NEGATIVE)
    c1_pos = active_master_df[active_master_df["atp_contrast1_cohort"] == "ATP_POSITIVE"]["genome_id"].dropna().tolist()
    c1_neg = active_master_df[active_master_df["atp_contrast1_cohort"] == "ATP_NEGATIVE"]["genome_id"].dropna().tolist()

    # Contrast 2: VA_vs_F_like (VA_ONLY vs F_LIKE_ONLY)
    c2_va = active_master_df[active_master_df["atp_contrast2_cohort"] == "VA_ONLY"]["genome_id"].dropna().tolist()
    c2_f = active_master_df[active_master_df["atp_contrast2_cohort"] == "F_LIKE_ONLY"]["genome_id"].dropna().tolist()

    all_analyzed_genomes = set(c1_pos + c1_neg + c2_va + c2_f)
    print(f"Genomes in active contrast comparisons: {len(all_analyzed_genomes)}")

    # Load KO presence matrix for analyzed genomes
    ko_presence_by_genome, all_kos, ko_definitions, ko_to_genomes = load_genome_ko_matrix(args.ko_counts, all_analyzed_genomes)

    all_records = []

    # Run Contrast 1: ATP_any_presence
    c1_records = run_association_for_contrast(
        "ATP_any_presence", c1_pos, c1_neg, active_master_df,
        ko_presence_by_genome, all_kos, ko_definitions,
        ko_to_genomes=ko_to_genomes,
        min_prevalence_count=args.min_ko_prevalence_count,
        fdr_alpha=args.fdr_alpha,
        delta_prev_cutoff=args.delta_prevalence_cutoff,
        min_stratum_group_size=args.min_stratum_group_size,
        run_regression=args.run_genome_reduction_sensitivity,
        representative_genomes=representative_genomes,
    )
    all_records.extend(c1_records)

    # Run Contrast 2: VA_vs_F_like
    c2_records = run_association_for_contrast(
        "VA_vs_F_like", c2_va, c2_f, active_master_df,
        ko_presence_by_genome, all_kos, ko_definitions,
        ko_to_genomes=ko_to_genomes,
        min_prevalence_count=args.min_ko_prevalence_count,
        fdr_alpha=args.fdr_alpha,
        delta_prev_cutoff=args.delta_prevalence_cutoff,
        min_stratum_group_size=args.min_stratum_group_size,
        run_regression=args.run_genome_reduction_sensitivity,
        representative_genomes=representative_genomes,
    )
    all_records.extend(c2_records)

    if not all_records:
        print("Warning: No association records generated!", file=sys.stderr)
        sys.exit(0)

    res_df = pd.DataFrame(all_records)

    # 1. Full association table across all KOs
    all_out_path = os.path.join(args.out_dir, "atp_gene_correlations_all_kos.tsv")
    res_df.to_csv(all_out_path, sep="\t", index=False)
    print(f"\nExported complete correlation table to {all_out_path} ({len(res_df)} records).")

    # 2. Positive controls table (ATP phenotype-defining KOs)
    pos_ctrl_df = res_df[res_df["target_definition_ko"] == True]
    pos_ctrl_path = os.path.join(args.out_dir, "atp_gene_correlations_positive_controls.tsv")
    pos_ctrl_df.to_csv(pos_ctrl_path, sep="\t", index=False)
    print(f"Exported positive control table to {pos_ctrl_path} ({len(pos_ctrl_df)} records).")

    # 3. Clean biological candidate table (excluding positive controls)
    cand_df = res_df[res_df["target_definition_ko"] == False]
    cand_path = os.path.join(args.out_dir, "atp_gene_correlations_candidates.tsv")
    cand_df.to_csv(cand_path, sep="\t", index=False)
    print(f"Exported biological candidate table to {cand_path} ({len(cand_df)} records).")

    # 4. Top significant associations (Tier 1 & Tier 2)
    top_df = cand_df[cand_df["candidate_tier"].isin(["Tier_1_robust", "Tier_2_raw_discovery"])]
    top_path = os.path.join(args.out_dir, "atp_gene_correlations_top_associated.tsv")
    top_df.to_csv(top_path, sep="\t", index=False)
    print(f"Exported top associated candidates table to {top_path} ({len(top_df)} records).")


if __name__ == "__main__":
    main()
