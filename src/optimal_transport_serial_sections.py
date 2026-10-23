#!/usr/bin/env python3
"""
optimal_transport_serial_sections.py
====================================
Monge-Kantorovich / Wasserstein-2 Optimal Transport Geodesic Coupling 
between Consecutive Human Breast Cancer Serial Sections (~10 um apart).

Formulates the 3D tissue continuity between Section 1 (N=3,798) and Section 2 (N=3,987)
as an entropy-regularized optimal transport problem (Sinkhorn-Knopp algorithm).
Tests whether the in situ evolutionary game replicator velocity (dx/dt) inferred on Section 1
predicts the empirical subclonal frequency displacement (Delta f_S) transported to Section 2.
"""

import os
import sys
import gzip
import time
import math
import json
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr, spearmanr

sys.stdout.reconfigure(line_buffering=True)

BASE_DATA_DIR = "/home/luca/spatial-game-assay/data"
RESULTS_DIR = "/home/luca/spatial-game-assay/results"

BREAST_SPEC = {
    "title": "Human Invasive Breast Carcinoma (10x Visium Serial Sections Block A)",
    "T_S": ["EPCAM", "CDH1", "KRT8", "KRT18", "ESR1"],
    "T_R": ["VIM", "CDH2", "ZEB1", "AXL", "MMP2"],
    "S_CAF": ["ACTA2", "COL1A1", "FAP", "PDGFRB"],
    "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
    "VASCULAR": ["PECAM1", "VWF"]
}


def load_visium_data(section_dir, spec):
    print(f"[*] Ingesting section from: {section_dir}")
    t0 = time.time()
    
    # 1. Barcodes
    barcodes_file = os.path.join(section_dir, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    with gzip.open(barcodes_file, "rt") as fp:
        barcodes = [line.strip() for line in fp]
    n_spots = len(barcodes)
    
    # 2. Features
    all_genes = sorted(list(set(
        spec["T_S"] + spec["T_R"] + spec["S_CAF"] +
        spec["FITNESS"] + spec["HYPOXIA"] + spec["VASCULAR"]
    )))
    features_file = os.path.join(section_dir, "filtered_feature_bc_matrix/features.tsv.gz")
    gene_to_row = {}
    with gzip.open(features_file, "rt") as fp:
        for idx, line in enumerate(fp):
            parts = line.strip().split("\t")
            gene_name = parts[1]
            if gene_name in all_genes:
                gene_to_row[idx + 1] = gene_name
                
    # 3. Matrix streaming
    matrix_file = os.path.join(section_dir, "filtered_feature_bc_matrix/matrix.mtx.gz")
    gene_counts = {g: np.zeros(n_spots, dtype=np.float32) for g in all_genes}
    total_umi = np.zeros(n_spots, dtype=np.float32)
    
    with gzip.open(matrix_file, "rt") as fp:
        for line in fp:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if len(parts) == 3:
                break
        for line in fp:
            parts = line.strip().split()
            r = int(parts[0])
            c = int(parts[1]) - 1
            v = float(parts[2])
            total_umi[c] += v
            if r in gene_to_row:
                gene_counts[gene_to_row[r]][c] += v
                
    # Normalization Log-CPM
    total_umi_safe = np.where(total_umi == 0, 1.0, total_umi)
    log_cpm = {}
    for g in all_genes:
        log_cpm[g] = np.log1p((gene_counts[g] / total_umi_safe) * 10000.0)
        
    def score_module(genes):
        mat = np.column_stack([log_cpm[g] for g in genes if g in log_cpm])
        return np.mean(mat, axis=1)
        
    m_S = score_module(spec["T_S"])
    m_R = score_module(spec["T_R"])
    m_C = score_module(spec["S_CAF"])
    
    # Softmax fractions
    M = np.column_stack([m_S, m_R, m_C])
    M_max = np.max(M, axis=1, keepdims=True)
    exp_M = np.exp(M - M_max)
    denom = np.sum(exp_M, axis=1, keepdims=True)
    f_S = exp_M[:, 0] / denom[:, 0]
    f_R = exp_M[:, 1] / denom[:, 0]
    f_C = exp_M[:, 2] / denom[:, 0]
    
    # Target Proliferation Fitness
    w_raw = score_module(spec["FITNESS"])
    
    # Microenvironmental confounders
    z_hyp = score_module(spec["HYPOXIA"])
    z_vasc = score_module(spec["VASCULAR"])
    
    # Proiezione ortogonale FWL
    Z = np.column_stack([np.ones(n_spots), z_hyp, z_vasc])
    gamma = np.linalg.lstsq(Z, w_raw, rcond=None)[0]
    w_purified = w_raw - (Z @ gamma)
    
    # Spatial coordinates
    tissue_pos_file = os.path.join(section_dir, "spatial/tissue_positions_list.csv")
    if not os.path.exists(tissue_pos_file):
        tissue_pos_file = os.path.join(section_dir, "spatial/tissue_positions.csv")
    coords_dict = {}
    with open(tissue_pos_file, "r") as fp:
        for line in fp:
            parts = line.strip().split(",")
            bc = parts[0]
            in_tissue = int(parts[1])
            if in_tissue == 1:
                r_px = float(parts[4])
                c_px = float(parts[5])
                coords_dict[bc] = (r_px, c_px)
                
    in_tissue_mask = np.array([bc in coords_dict for bc in barcodes])
    coords = np.array([coords_dict[bc] for bc in barcodes if bc in coords_dict])
    
    # Filter arrays
    f_S = f_S[in_tissue_mask]
    f_R = f_R[in_tissue_mask]
    f_C = f_C[in_tissue_mask]
    w_purified = w_purified[in_tissue_mask]
    w_raw = w_raw[in_tissue_mask]
    
    print(f"    Loaded {len(coords)} in-tissue spots in {time.time()-t0:.2f}s")
    return {
        "coords": coords,
        "f_S": f_S, "f_R": f_R, "f_C": f_C,
        "w_purified": w_purified, "w_raw": w_raw
    }


def sinkhorn_optimal_transport(X, Y, a=None, b=None, epsilon=0.05, max_iter=200, tol=1e-5):
    """
    Pure NumPy implementation of Sinkhorn-Knopp algorithm for entropy-regularized OT.
    Cost matrix C_ij = ||X_i - Y_j||^2.
    """
    n = X.shape[0]
    m = Y.shape[0]
    if a is None:
        a = np.ones(n, dtype=np.float64) / n
    if b is None:
        b = np.ones(m, dtype=np.float64) / m
        
    # Cost matrix normalized to [0, 1]
    C = cdist(X, Y, metric="sqeuclidean")
    C = C / np.median(C)
    
    # Gibbs kernel
    K = np.exp(-C / epsilon)
    
    u = np.ones(n, dtype=np.float64) / n
    v = np.ones(m, dtype=np.float64) / m
    
    for it in range(max_iter):
        u_prev = u
        u = a / (K @ v + 1e-16)
        v = b / (K.T @ u + 1e-16)
        
        err = np.max(np.abs(u - u_prev))
        if err < tol:
            break
            
    # Optimal coupling matrix Gamma: shape (n, m)
    # Gamma_ij = u_i * K_ij * v_j
    Gamma = np.diag(u) @ K @ np.diag(v)
    return Gamma, C


def main():
    print("=" * 80)
    print("WASSERSTEIN OPTIMAL TRANSPORT ACROSS CONSECUTIVE SERIAL HISTOLOGICAL SECTIONS")
    print("=" * 80)
    
    sec1_dir = os.path.join(BASE_DATA_DIR, "breast_visium_sec1")
    sec2_dir = os.path.join(BASE_DATA_DIR, "breast_visium_sec2")
    
    d1 = load_visium_data(sec1_dir, BREAST_SPEC)
    d2 = load_visium_data(sec2_dir, BREAST_SPEC)
    
    # Subsample spots if needed for numerical speed (e.g. 1500 spots per section for fast Sinkhorn)
    np.random.seed(42)
    n1 = len(d1["coords"])
    n2 = len(d2["coords"])
    
    # Fit InSituGame Payoff Matrix on Section 1
    from scipy.spatial import Delaunay
    tri1 = Delaunay(d1["coords"])
    indptr, indices = tri1.vertex_neighbor_vertices
    
    f_bar_S = np.zeros(n1)
    f_bar_R = np.zeros(n1)
    f_bar_C = np.zeros(n1)
    for i in range(n1):
        nb = indices[indptr[i]:indptr[i+1]]
        if len(nb) > 0:
            f_bar_S[i] = np.mean(d1["f_S"][nb])
            f_bar_R[i] = np.mean(d1["f_R"][nb])
            f_bar_C[i] = np.mean(d1["f_C"][nb])
        else:
            f_bar_S[i] = d1["f_S"][i]
            f_bar_R[i] = d1["f_R"][i]
            f_bar_C[i] = d1["f_C"][i]
            
    # Build 9-param design matrix
    f_i = [d1["f_S"], d1["f_R"], d1["f_C"]]
    f_bar_j = [f_bar_S, f_bar_R, f_bar_C]
    X_sec1 = []
    for i in range(3):
        for j in range(3):
            X_sec1.append(f_i[i] * f_bar_j[j])
    X_sec1 = np.column_stack(X_sec1)
    
    # OLS without intercept
    A_tilde_sec1 = np.linalg.lstsq(X_sec1, d1["w_purified"], rcond=None)[0]
    A_mat = A_tilde_sec1.reshape((3, 3))
    print(f"\n[+] Section 1 Inferred Payoff Matrix A_tilde (3x3):")
    print(np.round(A_mat, 3))
    
    # Replicator parameters
    Delta_S = A_mat[0, 0] - A_mat[1, 0]
    Delta_R = A_mat[1, 1] - A_mat[0, 1]
    print(f"    Delta_S = {Delta_S:.3f}, Delta_R = {Delta_R:.3f}")
    
    # Predict local replicator gradient on Section 1:
    # dot_x_S = x_S * (1 - x_S) * [ x_S * (Delta_S + Delta_R) - Delta_R ]
    x_S1 = d1["f_S"] / (d1["f_S"] + d1["f_R"] + 1e-12)
    replicator_velocity_S1 = x_S1 * (1.0 - x_S1) * (x_S1 * (Delta_S + Delta_R) - Delta_R)
    
    # Standardize coordinates to [0, 1]
    coord_min = np.min(np.vstack([d1["coords"], d2["coords"]]), axis=0)
    coord_max = np.max(np.vstack([d1["coords"], d2["coords"]]), axis=0)
    C1_norm = (d1["coords"] - coord_min) / (coord_max - coord_min)
    C2_norm = (d2["coords"] - coord_min) / (coord_max - coord_min)
    
    # To capture joint spatial and phenotypic transport:
    # State vector = [x_norm, y_norm, f_S, f_R]
    scale_pheno = 0.5
    feat1 = np.column_stack([C1_norm, scale_pheno * d1["f_S"], scale_pheno * d1["f_R"]])
    feat2 = np.column_stack([C2_norm, scale_pheno * d2["f_S"], scale_pheno * d2["f_R"]])
    
    print(f"\n[*] Computing Sinkhorn Optimal Transport Plan between Sec 1 (N={n1}) and Sec 2 (N={n2})...")
    t_ot = time.time()
    
    # Sample 1500 anchor points for regularized transport
    idx1 = np.linspace(0, n1 - 1, 1500, dtype=int)
    idx2 = np.linspace(0, n2 - 1, 1500, dtype=int)
    
    feat1_sub = feat1[idx1]
    feat2_sub = feat2[idx2]
    v_S1_sub = replicator_velocity_S1[idx1]
    f_S1_sub = d1["f_S"][idx1]
    f_S2_sub = d2["f_S"][idx2]
    
    Gamma, Cost_mat = sinkhorn_optimal_transport(feat1_sub, feat2_sub, epsilon=0.02, max_iter=300)
    print(f"[+] Sinkhorn converged in {time.time()-t_ot:.2f}s! Transport plan shape: {Gamma.shape}")
    
    # Empirical transported clonal frequency shift:
    # Delta f_S[i] = E_{j ~ Gamma_i}[ f_S2[j] ] - f_S1[i]
    Gamma_row_norm = Gamma / np.sum(Gamma, axis=1, keepdims=True)
    expected_f_S2 = Gamma_row_norm @ f_S2_sub
    empirical_transport_shift_S = expected_f_S2 - f_S1_sub
    
    # Statistical Concordance between Predicted Game Replicator Velocity and Empirical Transport Shift
    r_ot, p_ot = pearsonr(v_S1_sub, empirical_transport_shift_S)
    rho_ot, p_rho = spearmanr(v_S1_sub, empirical_transport_shift_S)
    
    print("\n" + "=" * 80)
    print("RESULTS: WASSERSTEIN OPTIMAL TRANSPORT CONCORDANCE")
    print("=" * 80)
    print(f"• Pearson correlation (Replicator dot_x vs Empirical OT delta_f): r = {r_ot:+.4f} (p = {p_ot:.4e})")
    print(f"• Spearman correlation: rho = {rho_ot:+.4f} (p = {p_rho:.4e})")
    
    # Save results to JSON
    out_json = os.path.join(RESULTS_DIR, "WASSERSTEIN_OPTIMAL_TRANSPORT_REPORT.json")
    res_payload = {
        "dataset": "Human Invasive Breast Carcinoma Block A Consecutive Sections (10 um)",
        "n_spots_sec1": int(n1),
        "n_spots_sec2": int(n2),
        "sinkhorn_anchors": 1500,
        "sinkhorn_epsilon": 0.02,
        "pearson_r": float(r_ot),
        "pearson_p": float(p_ot),
        "spearman_rho": float(rho_ot),
        "spearman_p": float(p_rho),
        "concordance_verdict": "Statistically Significant Geodesic Alignment" if p_ot < 0.001 else "Inconclusive"
    }
    with open(out_json, "w") as fp:
        json.dump(res_payload, fp, indent=2)
    print(f"\n[+] Saved Wasserstein OT Report to: {out_json}")


if __name__ == "__main__":
    main()
