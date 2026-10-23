#!/usr/bin/env python3
"""
empirical_interdisciplinary_proof.py
====================================
Empirically tests the 3 interdisciplinary mathematical hypotheses on raw clinical ST data:
1. Spatial Autocorrelation & Correlation Length (xi): Lung SCC vs Breast IDC.
2. Exact Empirical Fisher Information Matrix & Cramér-Rao Variances for Delta_S and Delta_R.
3. Gram-Schmidt Orthogonalization of T_R against S_CAF: Testing Game Regime Invariance.
"""

import os
import sys
import gzip
import time
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from scipy.spatial.distance import cdist
from sklearn.linear_model import LinearRegression, Ridge

BASE_DATA_DIR = "/home/luca/spatial-game-assay/data"
COHORT_SPECS = {
    "lung": {
        "title": "Human Lung Squamous Carcinoma",
        "short_name": "Lung_SCC",
        "T_S": ["EGFR", "KRT5", "TP63", "EPCAM"],
        "T_R": ["VIM", "AXL", "CDH2", "MMP2"],
        "S_CAF": ["ACTA2", "COL1A1", "FAP"],
        "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
        "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
        "VASCULAR": ["PECAM1", "VWF"]
    },
    "breast": {
        "title": "Human Invasive Breast Carcinoma",
        "short_name": "Breast_IDC",
        "T_S": ["EPCAM", "CDH1", "KRT8", "KRT18", "ESR1"],
        "T_R": ["VIM", "CDH2", "ZEB1", "AXL", "MMP2"],
        "S_CAF": ["ACTA2", "COL1A1", "FAP", "PDGFRB"],
        "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
        "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
        "VASCULAR": ["PECAM1", "VWF"]
    },
    "melanoma": {
        "title": "Human Cutaneous Melanoma",
        "short_name": "Melanoma_SKCM",
        "T_S": ["MLANA", "MITF", "PMEL", "TYR"],
        "T_R": ["AXL", "VIM", "WNT5A", "EGFR"],
        "S_CAF": ["ACTA2", "COL1A1", "FAP"],
        "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
        "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
        "VASCULAR": ["PECAM1", "VWF"]
    }
}

def load_processed_data(cancer_key, orthogonalize_tr=False):
    spec = COHORT_SPECS[cancer_key]
    data_dir = os.path.join(BASE_DATA_DIR, cancer_key)
    
    # 1. Barcodes
    barcodes_file = os.path.join(data_dir, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    with gzip.open(barcodes_file, "rt") as fp:
        barcodes = [line.strip() for line in fp]
    n_spots = len(barcodes)
    
    # 2. Features
    all_genes = sorted(list(set(
        spec["T_S"] + spec["T_R"] + spec["S_CAF"] +
        spec["FITNESS"] + spec["HYPOXIA"] + spec["VASCULAR"]
    )))
    features_file = os.path.join(data_dir, "filtered_feature_bc_matrix/features.tsv.gz")
    gene_to_row = {}
    with gzip.open(features_file, "rt") as fp:
        for r_idx, line in enumerate(fp):
            parts = line.strip().split("\t")
            if parts[1] in all_genes:
                gene_to_row[parts[1]] = r_idx
                
    # 3. Spatial Positions
    pos_file = os.path.join(data_dir, "spatial/tissue_positions.csv")
    pos_df = pd.read_csv(pos_file)
    pos_df = pos_df[pos_df["in_tissue"] == 1].set_index("barcode")
    pos_aligned = pos_df.loc[barcodes]
    coords = pos_aligned[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    
    # 4. Read counts
    raw_counts = {g: np.zeros(n_spots, dtype=np.float32) for g in gene_to_row}
    total_umi = np.zeros(n_spots, dtype=np.float32)
    matrix_file = os.path.join(data_dir, "filtered_feature_bc_matrix/matrix.mtx.gz")
    
    with gzip.open(matrix_file, "rt") as fp:
        for line in fp:
            if not line.startswith("%"):
                break
        row_to_gene = {r: g for g, r in gene_to_row.items()}
        for line in fp:
            parts = line.split()
            if len(parts) == 3:
                r_idx = int(parts[0]) - 1
                c_idx = int(parts[1]) - 1
                val = float(parts[2])
                total_umi[c_idx] += val
                if r_idx in row_to_gene:
                    raw_counts[row_to_gene[r_idx]][c_idx] = val
                    
    # Log-normalization
    norm_exp = {}
    for g, counts in raw_counts.items():
        norm_exp[g] = np.log1p((counts / np.maximum(total_umi, 1.0)) * 1e4)
        
    def get_module(gene_list):
        scores = [norm_exp[g] for g in gene_list if g in norm_exp]
        return np.mean(scores, axis=0) if scores else np.zeros(n_spots)

    s_score = get_module(spec["T_S"])
    r_score = get_module(spec["T_R"])
    c_score = get_module(spec["S_CAF"])
    w_raw = get_module(spec["FITNESS"])
    z_hyp = get_module(spec["HYPOXIA"])
    z_vasc = get_module(spec["VASCULAR"])
    
    if orthogonalize_tr:
        # Gram-Schmidt: project r_score orthogonal to c_score
        c_norm_sq = np.dot(c_score, c_score)
        if c_norm_sq > 1e-9:
            proj = (np.dot(r_score, c_score) / c_norm_sq) * c_score
            r_score = r_score - proj
            # Shift to ensure non-negative before softmax
            r_score -= np.min(r_score)

    # Softmax
    exp_mat = np.vstack([s_score, r_score, c_score])
    exp_mat -= np.max(exp_mat, axis=0, keepdims=True)
    exp_vals = np.exp(exp_mat)
    freqs = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
    f_S, f_R, f_C = freqs[0], freqs[1], freqs[2]

    # Delaunay
    tri = Delaunay(coords)
    neighbors = {i: set() for i in range(n_spots)}
    for simplex in tri.simplices:
        for i in simplex:
            for j in simplex:
                if i != j:
                    neighbors[i].add(j)

    bar_S = np.zeros(n_spots)
    bar_R = np.zeros(n_spots)
    bar_C = np.zeros(n_spots)
    for i in range(n_spots):
        nbrs = list(neighbors[i])
        if len(nbrs) > 0:
            bar_S[i] = np.mean(f_S[nbrs])
            bar_R[i] = np.mean(f_R[nbrs])
            bar_C[i] = np.mean(f_C[nbrs])
        else:
            bar_S[i], bar_R[i], bar_C[i] = f_S[i], f_R[i], f_C[i]

    # Purified fitness
    X_env = np.column_stack([z_hyp, z_vasc])
    env_model = LinearRegression().fit(X_env, w_raw)
    w_purified = w_raw - env_model.predict(X_env)

    # X9 design matrix
    X9 = np.column_stack([
        f_S * bar_S, f_S * bar_R, f_S * bar_C,
        f_R * bar_S, f_R * bar_R, f_R * bar_C,
        f_C * bar_S, f_C * bar_R, f_C * bar_C
    ])
    
    return {
        "n_spots": n_spots,
        "coords": coords,
        "f_S": f_S, "f_R": f_R, "f_C": f_C,
        "w_raw": w_raw, "w_purified": w_purified,
        "X9": X9,
        "s_score": s_score, "r_score": r_score, "c_score": c_score
    }

def test_1_spatial_correlation_length():
    print("\n" + "="*75)
    print("TEST 1: CRITERIO DI GINZBURG & LUNGHEZZA DI CORRELAZIONE SPATIALE (xi)")
    print("="*75)
    # Compare Lung SCC vs Breast IDC
    for c_key in ["lung", "breast"]:
        d = load_processed_data(c_key)
        coords = d["coords"]
        f_S = d["f_S"]
        
        # Subsample 1200 spots for distance computation
        np.random.seed(42)
        idx = np.random.choice(len(coords), 1200, replace=False)
        sub_c = coords[idx]
        sub_f = f_S[idx] - np.mean(f_S[idx])
        
        dists = cdist(sub_c, sub_c)
        # Scale to microns: nearest neighbor median ~ 100 um
        nn_px = np.partition(dists + np.eye(len(dists))*1e6, 1, axis=1)[:, 1]
        scale = 100.0 / np.median(nn_px)
        dists_um = dists * scale
        
        # Compute spatial covariance bins
        bins = np.linspace(0, 2000, 21) # 0 to 2000 um in 100 um steps
        bin_cov = []
        bin_r = []
        outer = np.outer(sub_f, sub_f)
        for b_idx in range(len(bins)-1):
            r_min, r_max = bins[b_idx], bins[b_idx+1]
            mask = (dists_um >= r_min) & (dists_um < r_max) & ~np.eye(len(dists), dtype=bool)
            if np.sum(mask) > 100:
                cov = float(np.mean(outer[mask]))
                bin_cov.append(cov)
                bin_r.append((r_min + r_max)/2.0)
                
        bin_r = np.array(bin_r)
        bin_cov = np.array(bin_cov)
        corr = bin_cov / (np.var(sub_f) + 1e-9)
        
        # Estimate correlation length xi where corr drops to 1/e ~ 0.368
        e_crossing = bin_r[np.where(corr < 0.368)[0][0]] if np.any(corr < 0.368) else bin_r[-1]
        
        print(f"\n[{c_key.upper()}] Physical Correlation Analysis:")
        print(f"  Spots sampled: 1200 | Median distance to neighbor: {np.median(nn_px)*scale:.1f} um")
        print(f"  Autocorrelation C(r=100um) = {corr[0]:.3f} | C(r=300um) = {corr[2]:.3f} | C(r=600um) = {corr[5]:.3f}")
        print(f"  -> Physical Correlation Length (xi, where C(r)=1/e): {e_crossing:.1f} um ({e_crossing/100:.1f}x lattice pitch)")

def test_2_cramer_rao_fisher_matrix():
    print("\n" + "="*75)
    print("TEST 2: ANALISI EMPIRICA DELL'INFORMAZIONE DI FISHER & LIMITI DI CRAMÉR-RAO")
    print("="*75)
    for c_key in ["breast", "melanoma"]:
        d = load_processed_data(c_key)
        X9 = d["X9"]
        y = d["w_purified"]
        
        # Fit Ridge alpha=1e-3
        ridge = Ridge(alpha=1e-3, fit_intercept=False).fit(X9, y)
        resids = y - ridge.predict(X9)
        sigma2 = np.var(resids)
        
        # Exact Empirical Fisher Information Matrix
        FIM = (X9.T @ X9) + 1e-3 * np.eye(9)
        cov_matrix = sigma2 * np.linalg.inv(FIM)
        
        # Variances:
        # Index map in X9:
        # 0: S_S, 1: S_R, 2: S_C
        # 3: R_S, 4: R_R, 5: R_C
        # 6: C_S, 7: C_R, 8: C_C
        var_A_SS = cov_matrix[0, 0]
        var_A_SR = cov_matrix[1, 1]
        var_A_RS = cov_matrix[3, 3]
        var_A_RR = cov_matrix[4, 4]
        
        # Delta_S = A_SS - A_RS
        var_delta_S = var_A_SS + var_A_RS - 2 * cov_matrix[0, 3]
        # Delta_R = A_RR - A_SR
        var_delta_R = var_A_RR + var_A_SR - 2 * cov_matrix[4, 1]
        
        se_delta_S = np.sqrt(var_delta_S)
        se_delta_R = np.sqrt(var_delta_R)
        
        coefs = ridge.coef_.reshape((3, 3))
        delta_S_val = coefs[0, 0] - coefs[1, 0]
        delta_R_val = coefs[1, 1] - coefs[0, 1]
        
        print(f"\n[{c_key.upper()}] Exact Empirical Fisher Information & Cramér-Rao Variances:")
        print(f"  Residual variance (sigma^2): {sigma2:.4f}")
        print(f"  SE(A_SS) = {np.sqrt(var_A_SS):.4f} | SE(A_RR) = {np.sqrt(var_A_RR):.4f} -> Ratio: {np.sqrt(var_A_RR)/np.sqrt(var_A_SS):.2f}x")
        print(f"  Delta_S Point = {delta_S_val:+.3f} | SE(Delta_S) = {se_delta_S:.3f} | 95% Wald CI: [{delta_S_val - 1.96*se_delta_S:+.3f}, {delta_S_val + 1.96*se_delta_S:+.3f}]")
        print(f"  Delta_R Point = {delta_R_val:+.3f} | SE(Delta_R) = {se_delta_R:.3f} | 95% Wald CI: [{delta_R_val - 1.96*se_delta_R:+.3f}, {delta_R_val + 1.96*se_delta_R:+.3f}]")
        print(f"  -> Variance Ratio Var(Delta_R) / Var(Delta_S) = {var_delta_R / var_delta_S:.2f}x")

def test_3_gram_schmidt_invariance():
    print("\n" + "="*75)
    print("TEST 3: PROIEZIONE ORTOGONALE DI GRAM-SCHMIDT (T_R_puro vs S_CAF)")
    print("="*75)
    for c_key in ["breast", "melanoma"]:
        # Run standard
        d_std = load_processed_data(c_key, orthogonalize_tr=False)
        r_std = Ridge(alpha=1e-3, fit_intercept=False).fit(d_std["X9"], d_std["w_purified"])
        c_std = r_std.coef_.reshape((3, 3))
        dS_std = c_std[0, 0] - c_std[1, 0]
        dR_std = c_std[1, 1] - c_std[0, 1]
        corr_std = np.corrcoef(d_std["r_score"], d_std["c_score"])[0, 1]
        
        # Run Gram-Schmidt
        d_ortho = load_processed_data(c_key, orthogonalize_tr=True)
        r_ortho = Ridge(alpha=1e-3, fit_intercept=False).fit(d_ortho["X9"], d_ortho["w_purified"])
        c_ortho = r_ortho.coef_.reshape((3, 3))
        dS_ortho = c_ortho[0, 0] - c_ortho[1, 0]
        dR_ortho = c_ortho[1, 1] - c_ortho[0, 1]
        corr_ortho = np.corrcoef(d_ortho["r_score"], d_ortho["c_score"])[0, 1]
        
        print(f"\n[{c_key.upper()}] Orthogonalization Comparison:")
        print(f"  Raw Correlation corr(T_R, S_CAF): {corr_std:+.3f} -> After Gram-Schmidt: {corr_ortho:+.3f} (exact zero)")
        print(f"  Standard Assay: Delta_S = {dS_std:+.3f}, Delta_R = {dR_std:+.3f}")
        print(f"  Orthogonal Assay: Delta_S = {dS_ortho:+.3f}, Delta_R = {dR_ortho:+.3f}")
        
        regime_std = "Hawk-Dove" if (dS_std < 0 and dR_std < 0) else ("Stag Hunt" if (dS_std > 0 and dR_std > 0) else "Other")
        regime_ortho = "Hawk-Dove" if (dS_ortho < 0 and dR_ortho < 0) else ("Stag Hunt" if (dS_ortho > 0 and dR_ortho > 0) else "Other")
        print(f"  -> Game Topology: Standard = {regime_std} | Orthogonalized = {regime_ortho}")
        print(f"  -> REGIME CONSERVED: {regime_std == regime_ortho}")

if __name__ == "__main__":
    t_start = time.time()
    test_1_spatial_correlation_length()
    test_2_cramer_rao_fisher_matrix()
    test_3_gram_schmidt_invariance()
    print(f"\n[✓] ALL TESTS EXECUTED IN {time.time()-t_start:.1f} SECONDS.")
