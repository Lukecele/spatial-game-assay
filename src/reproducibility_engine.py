#!/usr/bin/env python3
"""
reproducibility_engine.py
=========================
Exhaustive Statistical Verification, Diagnostic Performance, and Sensitivity Analysis.
Adapts Luca Celebrano's exact contingency testing, Odds Ratios (Woolf 95% CI),
Fisher's exact test, and Out-Of-Sample (OOS) validation framework from
meme-intelligence-onchain/analytics/backtest_validator.py.

Evaluations Performed:
1. Model Ablation Study:
   - Full Spatial Game Theory Model
   - Non-Spatial Mass Action Baseline (Intrinsic Only)
   - Physical Microenvironment Alone
   - Incremental F-test, Delta-AUC, and Diagnostic Odds Ratio (DOR)
2. Diagnostic Classification Metrics:
   - Sensitivity, Specificity, PPV, NPV, Accuracy
   - Exact Fisher's Exact p-value
   - 95% CI Odds Ratio
3. Quality Control (QC) Cutoff Robustness:
   - Unfiltered (raw in_tissue == 1)
   - Standard QC (total_umi >= 500)
   - High-Stringency QC (total_umi >= 1000)
4. Spatial Neighborhood Topology Robustness:
   - Delaunay triangulation (immediate physical contact)
   - k-NN with k=6
   - k-NN with k=12
   - k-NN with k=18
"""

import os
import sys
import gzip
import json
import time
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, cKDTree
from scipy.stats import f as scipy_f
from sklearn.linear_model import Ridge, LinearRegression

sys.stdout.reconfigure(line_buffering=True)

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from clinical_diagnostic_metrics import (
    ContingencyTable2x2, calculate_fisher_exact_p_value,
    calculate_odds_ratio, compute_roc_auc, evaluate_clinical_classification
)

BASE_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
RESULTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results"))

COHORT_SPECS = {
    "lung": {
        "title": "Human Lung Squamous Cell Carcinoma",
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
        "title": "Human Cutaneous Malignant Melanoma",
        "short_name": "Melanoma_SKCM",
        "T_S": ["MLANA", "MITF", "PMEL", "TYR"],
        "T_R": ["AXL", "VIM", "WNT5A", "EGFR"],
        "S_CAF": ["ACTA2", "COL1A1", "FAP"],
        "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
        "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
        "VASCULAR": ["PECAM1", "VWF"]
    }
}


def load_raw_dataset(cancer_key):
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
                
    # 3. Spatial Coordinates
    pos_file = os.path.join(data_dir, "spatial/tissue_positions.csv")
    if not os.path.exists(pos_file):
        pos_file = os.path.join(data_dir, "spatial/tissue_positions_list.csv")
        pos_df = pd.read_csv(pos_file, header=None,
                             names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    else:
        pos_df = pd.read_csv(pos_file)
        
    pos_df = pos_df[pos_df["in_tissue"] == 1].set_index("barcode")
    pos_aligned = pos_df.loc[barcodes]
    coords = pos_aligned[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    
    # 4. Streaming Matrix Read
    target_rows = set(gene_to_row.values())
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
                s_idx = int(parts[1]) - 1
                val = float(parts[2])
                total_umi[s_idx] += val
                if r_idx in target_rows:
                    raw_counts[row_to_gene[r_idx]][s_idx] = val
                    
    return barcodes, coords, raw_counts, total_umi, spec


def build_phenotypes_and_graphs(coords, raw_counts, total_umi, spec, min_umi=0, topology="delaunay", k_nn=6):
    # Apply Quality Control gate
    qc_mask = (total_umi >= min_umi)
    coords_sub = coords[qc_mask]
    total_sub = total_umi[qc_mask]
    n_sub = len(coords_sub)
    
    if n_sub < 100:
        return None
        
    norm_exp = {g: np.log1p((raw_counts[g][qc_mask] / np.maximum(total_sub, 1.0)) * 1e4) for g in raw_counts}
    
    def get_mod(genes):
        s = [norm_exp[g] for g in genes if g in norm_exp]
        return np.mean(s, axis=0) if s else np.zeros(n_sub)
        
    s_score = get_mod(spec["T_S"])
    r_score = get_mod(spec["T_R"])
    c_score = get_mod(spec["S_CAF"])
    w_raw = get_mod(spec["FITNESS"])
    z_hyp = get_mod(spec["HYPOXIA"])
    z_vasc = get_mod(spec["VASCULAR"])
    
    # Softmax phenotypic frequencies
    exp_mat = np.vstack([s_score, r_score, c_score])
    exp_mat -= np.max(exp_mat, axis=0, keepdims=True)
    exp_vals = np.exp(exp_mat)
    freqs = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
    f_S, f_R, f_C = freqs[0], freqs[1], freqs[2]
    
    # Neighborhood contact topology
    if topology == "delaunay":
        tri = Delaunay(coords_sub)
        neighbors = {i: set() for i in range(n_sub)}
        for simplex in tri.simplices:
            for i in simplex:
                for j in simplex:
                    if i != j:
                        neighbors[i].add(j)
    elif topology == "knn":
        tree = cKDTree(coords_sub)
        _, indices = tree.query(coords_sub, k=k_nn + 1)
        neighbors = {i: set(indices[i, 1:].tolist()) for i in range(n_sub)}
        
    bar_f_S = np.zeros(n_sub)
    bar_f_R = np.zeros(n_sub)
    bar_f_C = np.zeros(n_sub)
    for i in range(n_sub):
        nbrs = list(neighbors[i])
        if len(nbrs) > 0:
            bar_f_S[i] = np.mean(f_S[nbrs])
            bar_f_R[i] = np.mean(f_R[nbrs])
            bar_f_C[i] = np.mean(f_C[nbrs])
        else:
            bar_f_S[i], bar_f_R[i], bar_f_C[i] = f_S[i], f_R[i], f_C[i]
            
    # Environmental Decoupling
    X_env = np.column_stack([z_hyp, z_vasc])
    env_model = LinearRegression().fit(X_env, w_raw)
    w_purified = w_raw - env_model.predict(X_env)
    
    return {
        "n_spots": n_sub,
        "coords": coords_sub,
        "f_S": f_S, "f_R": f_R, "f_C": f_C,
        "bar_f_S": bar_f_S, "bar_f_R": bar_f_R, "bar_f_C": bar_f_C,
        "w_raw": w_raw, "w_purified": w_purified,
        "z_hyp": z_hyp, "z_vasc": z_vasc,
        "X_env": X_env
    }


def run_model_ablation(data):
    """
    Evaluates incremental diagnostic performance across 3 nested models:
    Model 1: Microenvironment Alone (Z_hyp, Z_vasc)
    Model 2: Non-Spatial Intrinsic Mass-Action (f_S, f_R, f_C)
    Model 3: Full Spatial Evolutionary Game (Intrinsic + Interactions)
    """
    w_raw = data["w_raw"]
    w_pur = data["w_purified"]
    f_S, f_R, f_C = data["f_S"], data["f_R"], data["f_C"]
    bar_S, bar_R, bar_C = data["bar_f_S"], data["bar_f_R"], data["bar_f_C"]
    X_env = data["X_env"]
    n = len(w_raw)
    
    # Model 1: Microenvironment Alone
    m1 = LinearRegression().fit(X_env, w_raw)
    y_pred_m1 = m1.predict(X_env)
    r2_m1 = float(m1.score(X_env, w_raw))
    diag_m1 = evaluate_clinical_classification(w_raw, y_pred_m1)
    
    # Model 2: Intrinsic Phenotypes Alone (on purified fitness)
    X_int = np.column_stack([f_S, f_R, f_C])
    m2 = Ridge(alpha=1e-3, fit_intercept=False).fit(X_int, w_pur)
    y_pred_m2 = m2.predict(X_int)
    r2_m2 = float(m2.score(X_int, w_pur))
    diag_m2 = evaluate_clinical_classification(w_pur, y_pred_m2)
    
    # Model 3: Full Spatial Evolutionary Game
    X_game = np.column_stack([
        f_S, f_R, f_C,
        f_S * bar_S, f_S * bar_R, f_S * bar_C,
        f_R * bar_S, f_R * bar_R, f_R * bar_C,
        f_C * bar_S, f_C * bar_R, f_C * bar_C
    ])
    m3 = Ridge(alpha=1e-3, fit_intercept=False).fit(X_game, w_pur)
    y_pred_m3 = m3.predict(X_game)
    r2_m3 = float(m3.score(X_game, w_pur))
    diag_m3 = evaluate_clinical_classification(w_pur, y_pred_m3)
    
    # Incremental Nested F-Test (Model 2 vs Model 3)
    rss_int = np.sum((w_pur - y_pred_m2) ** 2)
    rss_game = np.sum((w_pur - y_pred_m3) ** 2)
    p_diff = 9  # 9 additional interaction parameters
    df_game = n - 12
    f_stat = ((rss_int - rss_game) / p_diff) / (rss_game / df_game)
    f_p_val = float(1.0 - scipy_f.cdf(f_stat, p_diff, df_game))
    
    # Sub-game 2x2
    coefs = m3.coef_
    A_mat = coefs[3:12].reshape((3, 3))
    A_SS, A_SR = A_mat[0, 0], A_mat[0, 1]
    A_RS, A_RR = A_mat[1, 0], A_mat[1, 1]
    
    return {
        "n_spots": n,
        "model1_env": {"r2": r2_m1, "diagnostic": diag_m1},
        "model2_intrinsic": {"r2": r2_m2, "diagnostic": diag_m2},
        "model3_game": {"r2": r2_m3, "diagnostic": diag_m3},
        "incremental_f_test": {
            "f_statistic": float(f_stat),
            "df_num": p_diff,
            "df_denom": df_game,
            "p_value": float(f_p_val)
        },
        "delta_r2": float(r2_m3 - r2_m2),
        "delta_auc": float(diag_m3["roc_auc"] - diag_m2["roc_auc"]),
        "payoff_2x2": {"A_SS": float(A_SS), "A_SR": float(A_SR), "A_RS": float(A_RS), "A_RR": float(A_RR)}
    }


def execute_full_reproducibility_audit():
    print("="*80)
    print("IN SITU GAME ASSAY: EXHAUSTIVE REPRODUCIBILITY & SENSITIVITY ENGINE")
    print("="*80)
    
    cancers = ["lung", "breast", "melanoma"]
    all_audit_results = {}
    
    for c_key in cancers:
        print(f"\n{'='*60}\nPROCESSING: {c_key.upper()}\n{'='*60}")
        barcodes, coords, raw_counts, total_umi, spec = load_raw_dataset(c_key)
        
        # 1. Model Ablation Study (Primary Delaunay, All spots)
        d_base = build_phenotypes_and_graphs(coords, raw_counts, total_umi, spec, min_umi=0, topology="delaunay")
        ablation = run_model_ablation(d_base)
        
        print("\n--- [A] MODEL ABLATION & CLINICAL CLASSIFICATION PERFORMANCE ---")
        m1 = ablation["model1_env"]["diagnostic"]
        m2 = ablation["model2_intrinsic"]["diagnostic"]
        m3 = ablation["model3_game"]["diagnostic"]
        
        print(f"  Model 1 (Physical Microenv) : R^2={ablation['model1_env']['r2']*100:5.2f}% | AUC={m1['roc_auc']:.3f} | DOR={m1['odds_ratio']:6.2f} [{m1['odds_ratio_95_ci'][0]:.2f}-{m1['odds_ratio_95_ci'][1]:.2f}] | Fisher p={m1['fisher_exact_p_value']:.2e}")
        print(f"  Model 2 (Intrinsic Phenotype): R^2={ablation['model2_intrinsic']['r2']*100:5.2f}% | AUC={m2['roc_auc']:.3f} | DOR={m2['odds_ratio']:6.2f} [{m2['odds_ratio_95_ci'][0]:.2f}-{m2['odds_ratio_95_ci'][1]:.2f}] | Fisher p={m2['fisher_exact_p_value']:.2e}")
        print(f"  Model 3 (Spatial Game Assay) : R^2={ablation['model3_game']['r2']*100:5.2f}% | AUC={m3['roc_auc']:.3f} | DOR={m3['odds_ratio']:6.2f} [{m3['odds_ratio_95_ci'][0]:.2f}-{m3['odds_ratio_95_ci'][1]:.2f}] | Fisher p={m3['fisher_exact_p_value']:.2e}")
        print(f"  [+] Incremental F-Test: F({ablation['incremental_f_test']['df_num']}, {ablation['incremental_f_test']['df_denom']}) = {ablation['incremental_f_test']['f_statistic']:.2f}, p = {ablation['incremental_f_test']['p_value']:.2e} {'★ (Significant)' if ablation['incremental_f_test']['p_value'] < 0.001 else ''}")
        print(f"  [+] Sensitivity={m3['sensitivity']*100:.1f}%, Specificity={m3['specificity']*100:.1f}%, PPV={m3['ppv_precision']*100:.1f}%, NPV={m3['npv']*100:.1f}%")
        
        # 2. Quality Control (QC) Cutoff Robustness
        print("\n--- [B] QUALITY CONTROL (QC) CUTOFF SENSITIVITY ---")
        qc_runs = {}
        for qc_cut in [0, 500, 1000]:
            d_qc = build_phenotypes_and_graphs(coords, raw_counts, total_umi, spec, min_umi=qc_cut, topology="delaunay")
            abl_qc = run_model_ablation(d_qc)
            p2 = abl_qc["payoff_2x2"]
            print(f"  Cutoff >= {qc_cut:>4} UMI (N={abl_qc['n_spots']:>4} spots): R^2={abl_qc['model3_game']['r2']*100:5.2f}% | AUC={abl_qc['model3_game']['diagnostic']['roc_auc']:.3f} | A_SS={p2['A_SS']:+.3f}, A_SR={p2['A_SR']:+.3f}, A_RS={p2['A_RS']:+.3f}, A_RR={p2['A_RR']:+.3f}")
            qc_runs[f"cut_{qc_cut}"] = abl_qc
            
        # 3. Neighborhood Scale Robustness
        print("\n--- [C] NEIGHBORHOOD INTERACTION SCALE SENSITIVITY ---")
        topo_runs = {}
        for topo, knn_val in [("delaunay", 6), ("knn", 6), ("knn", 12), ("knn", 18)]:
            label = "Delaunay (Contact)" if topo == "delaunay" else f"k-NN (k={knn_val})"
            d_topo = build_phenotypes_and_graphs(coords, raw_counts, total_umi, spec, min_umi=500, topology=topo, k_nn=knn_val)
            abl_topo = run_model_ablation(d_topo)
            p2 = abl_topo["payoff_2x2"]
            print(f"  Topology: {label:<18}: R^2={abl_topo['model3_game']['r2']*100:5.2f}% | AUC={abl_topo['model3_game']['diagnostic']['roc_auc']:.3f} | A_SS={p2['A_SS']:+.3f}, A_SR={p2['A_SR']:+.3f}, A_RS={p2['A_RS']:+.3f}, A_RR={p2['A_RR']:+.3f}")
            topo_runs[label] = abl_topo
            
        all_audit_results[c_key] = {
            "cancer_key": c_key,
            "short_name": spec["short_name"],
            "ablation": ablation,
            "qc_sensitivity": qc_runs,
            "topology_sensitivity": topo_runs
        }
        
    out_json = os.path.join(RESULTS_DIR, "reproducibility_and_sensitivity_audit.json")
    with open(out_json, "w") as fp:
        json.dump(all_audit_results, fp, indent=2)
    print(f"\n[✓] Comprehensive audit finished and saved to {out_json}")
    return all_audit_results


if __name__ == "__main__":
    execute_full_reproducibility_audit()
