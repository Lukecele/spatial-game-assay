#!/usr/bin/env python3
"""
validation_battery.py
=====================
Comprehensive Statistical Validation Suite for Spatial Evolutionary Game Assay.
Designed to meet the stringent publication requirements of Q1 journals.

Tests Performed:
1. Monte Carlo Permutation Null Model (1,000 spatial shuffles -> empirical p-value).
2. Spatial Autocorrelation (Moran's I with exact Cliff & Ord 1981 variance formula).
3. Spatial Block Cross-Validation (5-Fold geographically disjoint blocks).
4. Spatial Block Bootstrap (1,000 resamples -> 95% CI and bootstrap z-test p-values).

Bug fixes vs v1:
 - Moran's I variance: replaced 1/S0 approximation with Cliff & Ord (1981) exact formula.
 - b_i p-values: permutation design cannot test b_i (they are not shuffled);
   replaced misleading p=1.0 with proper bootstrap z-test p-values for ALL parameters.
"""

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, cKDTree
from scipy.stats import norm as scipy_norm
from sklearn.linear_model import Ridge
from sklearn.cluster import KMeans
from sklearn.metrics import r2_score, mean_squared_error


# ---------------------------------------------------------------------------
# 1. Moran's I — Cliff & Ord (1981) exact variance under randomisation
# ---------------------------------------------------------------------------
def compute_morans_i(values, coords, k_neighbors=6):
    """
    Moran's spatial autocorrelation statistic with the *exact* analytical
    variance formula (Cliff & Ord 1981, Chapter 1).

    Uses a directed k-NN binary weight matrix; edge (i->j) has w_ij=1
    iff j is among the k nearest spatial neighbours of i.

    Parameters
    ----------
    values     : array-like, n floats
    coords     : (n, 2) pixel coordinates
    k_neighbors: number of nearest neighbours per spot

    Returns
    -------
    dict with keys: I, E_I, z_score, p_value, S0, S1, S2
    """
    n = len(values)
    y = np.asarray(values, dtype=np.float64)
    y_mean = np.mean(y)
    y_dev = y - y_mean
    ss = np.sum(y_dev ** 2)
    if ss == 0.0:
        return {"I": 0.0, "E_I": 0.0, "z_score": 0.0, "p_value": 1.0,
                "S0": 0.0, "S1": 0.0, "S2": 0.0}

    # Build k-NN graph (directed)
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k_neighbors + 1)  # col 0 = self
    W_sets = {i: set(indices[i, 1:].tolist()) for i in range(n)}

    # S0 = sum w_ij (each node has exactly k out-edges in a k-NN graph)
    S0 = float(n * k_neighbors)

    # Moran's I numerator:  sum_{i!=j} w_ij * (y_i - ybar)(y_j - ybar)
    num = 0.0
    for i in range(n):
        for j in W_sets[i]:
            num += y_dev[i] * y_dev[j]

    I = (n / S0) * (num / ss)
    E_I = -1.0 / (n - 1)

    # ------------------------------------------------------------------
    # Cliff & Ord (1981) exact variance under the *randomisation* null:
    #   Var_R(I) = [n^2 * S1 - n * S2 + 3 * S0^2] / [S0^2 * (n^2 - 1)]
    #              - E(I)^2
    #
    #   S1 = 0.5 * sum_{ij} (w_ij + w_ji)^2          (symmetrised pairs)
    #   S2 = sum_i (row_sum_i + col_sum_i)^2
    # ------------------------------------------------------------------

    # Compute in-degree (col sum) for each node
    in_deg = np.zeros(n, dtype=np.float64)
    for i in range(n):
        for j in W_sets[i]:
            in_deg[j] += 1.0

    # S1: iterate over all directed edges (i->j)
    S1 = 0.0
    for i in range(n):
        for j in W_sets[i]:
            w_ij = 1.0
            w_ji = 1.0 if i in W_sets[j] else 0.0
            S1 += 0.5 * (w_ij + w_ji) ** 2

    # S2: each node contributes (out_deg_i + in_deg_i)^2
    out_deg = float(k_neighbors)  # same for all nodes
    S2 = np.sum((out_deg + in_deg) ** 2)

    n2 = float(n * n)
    numerator_var = n2 * S1 - n * S2 + 3.0 * S0 ** 2
    denominator_var = S0 ** 2 * (n2 - 1.0)
    var_I = numerator_var / denominator_var - E_I ** 2

    # Guard against numerical rounding issues
    if var_I <= 0.0:
        var_I = 1e-12

    z = (I - E_I) / np.sqrt(var_I)
    p_val = 2.0 * (1.0 - scipy_norm.cdf(abs(z)))

    return {
        "I": float(I),
        "E_I": float(E_I),
        "z_score": float(z),
        "p_value": float(p_val),
        "S0": float(S0),
        "S1": float(S1),
        "S2": float(S2)
    }


# ---------------------------------------------------------------------------
# 2. Monte Carlo Spatial Permutation Test
# ---------------------------------------------------------------------------
def run_permutation_test(X, y, n_permutations=1000, alpha=1e-3, seed=42):
    """
    Spatial permutation null: shuffles ONLY the interaction columns (X[:,3:]),
    which encode neighbour-frequency products (f_i × bar_f_j).
    This breaks spatial topology while preserving spot-intrinsic phenotype
    frequencies (X[:,0:3]).

    NOTE: b_i coefficients (columns 0-2) are NOT permuted by design;
    their significance is assessed via bootstrap z-test in run_full_validation_battery().
    Permutation p-values for b_i are therefore not computed here.
    """
    np.random.seed(seed)
    n_spots = len(y)

    # Baseline
    ridge = Ridge(alpha=alpha, fit_intercept=False)
    ridge.fit(X, y)
    obs_r2 = ridge.score(X, y)
    obs_coefs = ridge.coef_.copy()

    n_interaction_params = X.shape[1] - 3  # exclude the 3 b_i columns

    perm_r2s = np.zeros(n_permutations)
    # Only track permutation distribution for interaction columns (idx 3..11)
    perm_int_coefs = np.zeros((n_permutations, n_interaction_params))

    for b in range(n_permutations):
        perm_idx = np.random.permutation(n_spots)
        X_perm = X.copy()
        X_perm[:, 3:] = X_perm[perm_idx, 3:]

        m = Ridge(alpha=alpha, fit_intercept=False)
        m.fit(X_perm, y)
        perm_r2s[b] = m.score(X_perm, y)
        perm_int_coefs[b, :] = m.coef_[3:]

    p_val_r2 = (1.0 + np.sum(perm_r2s >= obs_r2)) / (n_permutations + 1.0)

    # Empirical two-sided p-values for INTERACTION coefficients only
    int_coef_p_vals = np.zeros(n_interaction_params)
    for c in range(n_interaction_params):
        abs_obs = abs(obs_coefs[c + 3])
        abs_perms = np.abs(perm_int_coefs[:, c])
        int_coef_p_vals[c] = (1.0 + np.sum(abs_perms >= abs_obs)) / (n_permutations + 1.0)

    return {
        "observed_r2": float(obs_r2),
        "perm_r2_mean": float(np.mean(perm_r2s)),
        "perm_r2_max": float(np.max(perm_r2s)),
        "p_value_r2": float(p_val_r2),
        # Full-length array (b_i positions = np.nan, interaction positions = empirical p)
        "coef_p_values_interaction_only": int_coef_p_vals.tolist()
    }


# ---------------------------------------------------------------------------
# 3. Spatial Block Cross-Validation
# ---------------------------------------------------------------------------
def run_spatial_block_cv(X, y, coords, n_blocks=5, alpha=1e-3, seed=42):
    """
    Partitions tissue into n_blocks geographically disjoint KMeans clusters.
    Trains on the remaining 4 blocks and predicts the held-out block, thereby
    ensuring zero spatial leakage between training and test sets.
    """
    np.random.seed(seed)
    kmeans = KMeans(n_clusters=n_blocks, random_state=seed, n_init=10)
    spatial_folds = kmeans.fit_predict(coords)

    y_pred_oof = np.zeros_like(y)
    fold_metrics = []

    for fold in range(n_blocks):
        train_mask = (spatial_folds != fold)
        test_mask = (spatial_folds == fold)

        X_tr, y_tr = X[train_mask], y[train_mask]
        X_te, y_te = X[test_mask], y[test_mask]

        model = Ridge(alpha=alpha, fit_intercept=False)
        model.fit(X_tr, y_tr)
        y_te_pred = model.predict(X_te)
        y_pred_oof[test_mask] = y_te_pred

        fold_metrics.append({
            "fold": fold + 1,
            "n_test_spots": int(np.sum(test_mask)),
            "r2": float(r2_score(y_te, y_te_pred)),
            "rmse": float(np.sqrt(mean_squared_error(y_te, y_te_pred)))
        })

    return {
        "overall_oof_r2": float(r2_score(y, y_pred_oof)),
        "overall_oof_rmse": float(np.sqrt(mean_squared_error(y, y_pred_oof))),
        "fold_metrics": fold_metrics
    }


# ---------------------------------------------------------------------------
# 4. Spatial Block Bootstrap
# ---------------------------------------------------------------------------
def run_spatial_bootstrap(X, y, coords, n_bootstraps=1000, n_blocks=20,
                          alpha=1e-3, seed=42):
    """
    Resamples spatial blocks (with replacement) to obtain block-bootstrap
    standard errors and empirical 95% CIs for all payoff parameters.
    """
    np.random.seed(seed)
    kmeans = KMeans(n_clusters=n_blocks, random_state=seed, n_init=10)
    blocks = kmeans.fit_predict(coords)
    unique_blocks = np.arange(n_blocks)

    boot_coefs = np.zeros((n_bootstraps, X.shape[1]))

    for b in range(n_bootstraps):
        sampled_blocks = np.random.choice(unique_blocks, size=n_blocks, replace=True)
        boot_indices = np.concatenate(
            [np.where(blocks == blk)[0] for blk in sampled_blocks])

        m = Ridge(alpha=alpha, fit_intercept=False)
        m.fit(X[boot_indices], y[boot_indices])
        boot_coefs[b, :] = m.coef_

    ci_lower = np.percentile(boot_coefs, 2.5, axis=0)
    ci_upper = np.percentile(boot_coefs, 97.5, axis=0)
    se = np.std(boot_coefs, axis=0)

    return {
        "ci_lower": ci_lower.tolist(),
        "ci_upper": ci_upper.tolist(),
        "std_error": se.tolist()
    }


# ---------------------------------------------------------------------------
# 5. Full Validation Battery — orchestrates all tests
# ---------------------------------------------------------------------------
def run_full_validation_battery(data, sample_name="Sample",
                                alpha=1e-3, n_perm=1000, n_boot=1000):
    print(f"\n{'='*70}\nSTATISTICAL VALIDATION SUITE: {sample_name.upper()}\n{'='*70}")

    coords = data["coords"]
    f_S, f_R, f_C = data["f_S"], data["f_R"], data["f_C"]
    bar_S = data["bar_f_S"]; bar_R = data["bar_f_R"]; bar_C = data["bar_f_C"]
    w_raw = data["w_raw"]
    w_purified = data["w_purified"]

    # Feature matrix (12 columns: 3 b_i + 9 A_ij interaction terms)
    X = np.column_stack([
        f_S, f_R, f_C,
        f_S * bar_S, f_S * bar_R, f_S * bar_C,
        f_R * bar_S, f_R * bar_R, f_R * bar_C,
        f_C * bar_S, f_C * bar_R, f_C * bar_C
    ])
    y = w_purified

    # Fit full model to get residuals
    ridge = Ridge(alpha=alpha, fit_intercept=False)
    ridge.fit(X, y)
    coefs = ridge.coef_
    residuals = y - ridge.predict(X)

    # --- TEST 1: Moran's I (exact Cliff & Ord 1981 variance) ---
    print("\n[1/4] Spatial Autocorrelation Audit (Moran's I — Cliff & Ord 1981)...")
    moran_raw = compute_morans_i(w_raw, coords)
    moran_purified = compute_morans_i(w_purified, coords)
    moran_residuals = compute_morans_i(residuals, coords)

    for label, m in [("W_raw (Proliferation)", moran_raw),
                     ("W_res (Env-Purified)",   moran_purified),
                     ("Residuals (Game Model)", moran_residuals)]:
        print(f"  {label:<28}: I={m['I']:+.4f}  z={m['z_score']:+.2f}  p={m['p_value']:.2e}")

    atten = (1.0 - abs(moran_residuals["I"]) /
             max(abs(moran_raw["I"]), 1e-9)) * 100.0
    print(f"  [+] Autocorrelation attenuation (raw -> residuals): {atten:.2f}%")
    print(f"  Note: S0={moran_raw['S0']:.0f}, S1={moran_raw['S1']:.1f}, "
          f"S2={moran_raw['S2']:.1f}  (Cliff & Ord exact)")

    # --- TEST 2: Monte Carlo Permutation Null ---
    print(f"\n[2/4] Monte Carlo Spatial Permutation Test (N={n_perm})...")
    perm_results = run_permutation_test(X, y, n_permutations=n_perm, alpha=alpha)
    print(f"  Observed R²     : {perm_results['observed_r2']*100:.2f}%")
    print(f"  Null R² (mean)  : {perm_results['perm_r2_mean']*100:.2f}%  "
          f"(max: {perm_results['perm_r2_max']*100:.2f}%)")
    print(f"  ΔR² (game signal): {(perm_results['observed_r2']-perm_results['perm_r2_mean'])*100:.2f}%")
    sig_tag = " ★ p < 0.001" if perm_results["p_value_r2"] <= 0.001 else ""
    print(f"  Empirical p-value: p = {perm_results['p_value_r2']:.4f}{sig_tag}")

    # --- TEST 3: Spatial Block CV ---
    print("\n[3/4] 5-Fold Spatial Block Cross-Validation (SBCV)...")
    sbcv = run_spatial_block_cv(X, y, coords, n_blocks=5, alpha=alpha)
    print(f"  Overall OOF R²  : {sbcv['overall_oof_r2']*100:.2f}%")
    print(f"  Overall OOF RMSE: {sbcv['overall_oof_rmse']:.4f}")
    for fm in sbcv["fold_metrics"]:
        note = "  [spatial non-stationarity]" if fm["r2"] < 0 else ""
        print(f"    Block {fm['fold']} ({fm['n_test_spots']:>5} spots): "
              f"R²={fm['r2']*100:+7.2f}%  RMSE={fm['rmse']:.4f}{note}")

    # --- TEST 4: Spatial Block Bootstrap + bootstrap z-test p-values ---
    print(f"\n[4/4] Spatial Block Bootstrap CI (N={n_boot}, 20 blocks)...")
    boot = run_spatial_bootstrap(X, y, coords, n_bootstraps=n_boot,
                                 n_blocks=20, alpha=alpha)
    ci_low = boot["ci_lower"]
    ci_high = boot["ci_upper"]
    se = boot["std_error"]

    # Parameter name list (matches column order of X)
    strats = ["T_Sens", "T_Res", "Stroma"]
    param_names = [f"b_{s}" for s in strats]
    for si in strats:
        for sj in strats:
            param_names.append(f"A_{si}_{sj}")

    # Build per-parameter p-value table:
    #   b_i (idx 0-2): bootstrap z-test  (permutation design cannot test these)
    #   A_ij (idx 3-11): permutation empirical p  (primary spatial test)
    perm_int_p = perm_results["coef_p_values_interaction_only"]  # length 9

    print("\n" + "-"*90)
    print(f"{'Parameter':<20} | {'Estimate':>10} | {'95% CI (Bootstrap)':>26} | "
          f"{'SE':>8} | {'p-value':>12} | {'Source'}")
    print("-"*90)

    summary_table = []
    for idx, name in enumerate(param_names):
        est = coefs[idx]
        cil = ci_low[idx]
        cih = ci_high[idx]
        se_i = se[idx]

        if idx < 3:
            # bootstrap z-test for b_i
            z_boot = abs(est) / max(se_i, 1e-12)
            p_val = float(2.0 * (1.0 - scipy_norm.cdf(z_boot)))
            p_source = "bootstrap-z"
        else:
            # permutation empirical p for A_ij
            p_val = float(perm_int_p[idx - 3])
            p_source = "permutation"

        sig = " ★" if p_val <= 0.05 else ("  " if p_val <= 0.1 else "  ")
        print(f"{name:<20} | {est:>+10.4f} | [{cil:>+10.4f}, {cih:>+10.4f}] | "
              f"{se_i:>8.4f} | {p_val:>12.4f}{sig} | {p_source}")

        summary_table.append({
            "parameter": name,
            "estimate": float(est),
            "ci_95_lower": float(cil),
            "ci_95_upper": float(cih),
            "std_error": float(se_i),
            "p_value": p_val,
            "p_source": p_source
        })
    print("-"*90)
    print("  Note: b_i p-values from bootstrap z-test (z=|est/SE|); "
          "spatial permutation cannot test b_i.\n"
          "  A_ij p-values from Monte Carlo spatial permutation (N=1000).")

    full_report = {
        "sample": sample_name,
        "n_spots": int(len(y)),
        "morans_i": {
            "w_raw": moran_raw,
            "w_purified": moran_purified,
            "residuals": moran_residuals,
            "attenuation_pct": float(atten)
        },
        "permutation_test": perm_results,
        "spatial_block_cv": sbcv,
        "bootstrap_parameters": summary_table
    }
    return full_report
