#!/usr/bin/env python3
"""
run_replicate_blind_oos_test.py
===============================
Blind Out-Of-Sample (OOS) Replicate Validation on Human Breast Cancer Serial Sections.
Evaluates model reproducibility, parameter concordance, and cross-section predictive validity
between two consecutive histological slices (~10 um apart) of human invasive ductal carcinoma:
- 10x Visium Block A Section 1 (data/breast_visium_sec1)
- 10x Visium Block A Section 2 (data/breast_visium_sec2)

Integrates Luca Celebrano's rigorous statistical methodology from meme-intelligence-onchain:
- Blind Out-Of-Sample generalization (trained exclusively on Sec1, evaluated blind on Sec2)
- Woolf logit 95% Confidence Interval Odds Ratios
- Exact Fisher's test via math.lgamma
- Full 2x2 contingency tables (Sensitivity, Specificity, PPV, NPV, ROC-AUC)
- Matrix parameter concordance (Pearson r, Spearman rho, MAE)
- Evolutionary Game Nash Equilibrium preservation test
"""

import os
import sys
import gzip
import json
import time
import math
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge, LinearRegression

sys.stdout.reconfigure(line_buffering=True)
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from clinical_diagnostic_metrics import (
    ContingencyTable2x2, calculate_fisher_exact_p_value,
    calculate_odds_ratio, compute_roc_auc, evaluate_clinical_classification
)

BASE_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
RESULTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results"))

BREAST_SPEC = {
    "title": "Human Invasive Breast Carcinoma (10x Visium Consecutive Serial Sections)",
    "T_S": ["EPCAM", "CDH1", "KRT8", "KRT18", "ESR1"],
    "T_R": ["VIM", "CDH2", "ZEB1", "AXL", "MMP2"],
    "S_CAF": ["ACTA2", "COL1A1", "FAP", "PDGFRB"],
    "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
    "VASCULAR": ["PECAM1", "VWF"]
}


def load_visium_section(section_dir, spec):
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
        for r_idx, line in enumerate(fp):
            parts = line.strip().split("\t")
            if parts[1] in all_genes:
                gene_to_row[parts[1]] = r_idx
                
    # 3. Spatial Coordinates
    pos_file = os.path.join(section_dir, "spatial/tissue_positions_list.csv")
    if not os.path.exists(pos_file):
        pos_file = os.path.join(section_dir, "spatial/tissue_positions.csv")
        pos_df = pd.read_csv(pos_file)
    else:
        pos_df = pd.read_csv(pos_file, header=None,
                             names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
                             
    pos_df = pos_df[pos_df["in_tissue"] == 1].set_index("barcode")
    # Match with filtered barcodes
    valid_barcodes = [bc for bc in barcodes if bc in pos_df.index]
    pos_aligned = pos_df.loc[valid_barcodes]
    coords = pos_aligned[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    
    bc_to_idx = {bc: i for i, bc in enumerate(valid_barcodes)}
    
    # 4. Matrix streaming
    target_rows = set(gene_to_row.values())
    raw_counts = {g: np.zeros(len(valid_barcodes), dtype=np.float32) for g in gene_to_row}
    total_umi = np.zeros(len(valid_barcodes), dtype=np.float32)
    
    matrix_file = os.path.join(section_dir, "filtered_feature_bc_matrix/matrix.mtx.gz")
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
                bc = barcodes[s_idx]
                if bc in bc_to_idx:
                    sub_idx = bc_to_idx[bc]
                    val = float(parts[2])
                    total_umi[sub_idx] += val
                    if r_idx in target_rows:
                        raw_counts[row_to_gene[r_idx]][sub_idx] = val
                        
    # 5. Normalization & Signatures
    norm_exp = {g: np.log1p((raw_counts[g] / np.maximum(total_umi, 1.0)) * 1e4) for g in raw_counts}
    
    def get_mod(genes):
        s = [norm_exp[g] for g in genes if g in norm_exp]
        return np.mean(s, axis=0) if s else np.zeros(len(valid_barcodes))
        
    s_score = get_mod(spec["T_S"])
    r_score = get_mod(spec["T_R"])
    c_score = get_mod(spec["S_CAF"])
    w_raw = get_mod(spec["FITNESS"])
    z_hyp = get_mod(spec["HYPOXIA"])
    z_vasc = get_mod(spec["VASCULAR"])
    
    # Softmax fractions
    exp_mat = np.vstack([s_score, r_score, c_score])
    exp_mat -= np.max(exp_mat, axis=0, keepdims=True)
    exp_vals = np.exp(exp_mat)
    freqs = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
    f_S, f_R, f_C = freqs[0], freqs[1], freqs[2]
    
    # Spatial Delaunay contact graph
    tri = Delaunay(coords)
    neighbors = {i: set() for i in range(len(valid_barcodes))}
    for simplex in tri.simplices:
        for i in simplex:
            for j in simplex:
                if i != j:
                    neighbors[i].add(j)
                    
    bar_f_S = np.zeros(len(valid_barcodes))
    bar_f_R = np.zeros(len(valid_barcodes))
    bar_f_C = np.zeros(len(valid_barcodes))
    for i in range(len(valid_barcodes)):
        nbrs = list(neighbors[i])
        if len(nbrs) > 0:
            bar_f_S[i] = np.mean(f_S[nbrs])
            bar_f_R[i] = np.mean(f_R[nbrs])
            bar_f_C[i] = np.mean(f_C[nbrs])
        else:
            bar_f_S[i], bar_f_R[i], bar_f_C[i] = f_S[i], f_R[i], f_C[i]
            
    # Microenvironmental decoupling
    X_env = np.column_stack([z_hyp, z_vasc])
    env_model = LinearRegression().fit(X_env, w_raw)
    w_purified = w_raw - env_model.predict(X_env)
    
    elapsed = time.time() - t0
    print(f"  [+] Loaded {len(valid_barcodes)} spots in {elapsed:.1f}s | Median UMI = {np.median(total_umi):.0f}")
    
    return {
        "barcodes": valid_barcodes,
        "n_spots": len(valid_barcodes),
        "coords": coords,
        "total_umi": total_umi,
        "f_S": f_S, "f_R": f_R, "f_C": f_C,
        "bar_f_S": bar_f_S, "bar_f_R": bar_f_R, "bar_f_C": bar_f_C,
        "w_raw": w_raw, "w_purified": w_purified,
        "X_env": X_env, "env_model": env_model
    }


def construct_design_matrix(sec_data):
    f_S, f_R, f_C = sec_data["f_S"], sec_data["f_R"], sec_data["f_C"]
    bar_S, bar_R, bar_C = sec_data["bar_f_S"], sec_data["bar_f_R"], sec_data["bar_f_C"]
    
    # Intrinsic design (3 cols)
    X_int = np.column_stack([f_S, f_R, f_C])
    
    # Full Spatial Game design (12 cols: 3 intrinsic + 9 interaction)
    X_game = np.column_stack([
        f_S, f_R, f_C,
        f_S * bar_S, f_S * bar_R, f_S * bar_C,
        f_R * bar_S, f_R * bar_R, f_R * bar_C,
        f_C * bar_S, f_C * bar_R, f_C * bar_C
    ])
    return X_int, X_game


def solve_subgame_nash(A_SS, A_SR, A_RS, A_RR):
    """Solves exact 2x2 sub-game Nash equilibria between S and R."""
    delta_S = A_SS - A_RS
    delta_R = A_RR - A_SR
    
    # Check dominance
    if delta_S > 0 and delta_R < 0:
        return "Sensitive Dominance", {"x_pure": 1.0}
    elif delta_S < 0 and delta_R > 0:
        return "Resistant Dominance", {"x_pure": 0.0}
    elif delta_S > 0 and delta_R > 0:
        x_sep = delta_R / (delta_S + delta_R)
        return "Stag Hunt (Bistable)", {"x_sep": float(x_sep), "x_equilibria": [0.0, 1.0]}
    else:  # delta_S < 0 and delta_R < 0
        denom = (A_SS - A_SR - A_RS + A_RR)
        x_star = (A_RR - A_SR) / denom if abs(denom) > 1e-9 else 0.5
        return "Hawk-Dove / Snowdrift (Coexistence)", {"x_star": float(x_star)}


def execute_blind_oos_validation():
    print("="*80)
    print("BLIND OUT-OF-SAMPLE (OOS) VALIDATION: SERIAL HISTOLOGICAL SECTIONS")
    print("Dataset: 10x Visium Block A Section 1 & Section 2 (~10 um consecutive slices)")
    print("="*80)
    
    sec1_dir = os.path.join(BASE_DATA_DIR, "breast_visium_sec1")
    sec2_dir = os.path.join(BASE_DATA_DIR, "breast_visium_sec2")
    
    sec1 = load_visium_section(sec1_dir, BREAST_SPEC)
    sec2 = load_visium_section(sec2_dir, BREAST_SPEC)
    
    X_int1, X_game1 = construct_design_matrix(sec1)
    X_int2, X_game2 = construct_design_matrix(sec2)
    
    # =========================================================================
    # STEP 1: IN-SAMPLE TRAINING ON SECTION 1
    # =========================================================================
    print("\n[+] STEP 1: In-Sample Model Training on Section 1 (N = 3,798 spots)...")
    m_int1 = Ridge(alpha=1e-3, fit_intercept=False).fit(X_int1, sec1["w_purified"])
    m_game1 = Ridge(alpha=1e-3, fit_intercept=False).fit(X_game1, sec1["w_purified"])
    
    pred_int1 = m_int1.predict(X_int1)
    pred_game1 = m_game1.predict(X_game1)
    
    r2_int1 = float(m_int1.score(X_int1, sec1["w_purified"]))
    r2_game1 = float(m_game1.score(X_game1, sec1["w_purified"]))
    
    diag_int1 = evaluate_clinical_classification(sec1["w_purified"], pred_int1)
    diag_game1 = evaluate_clinical_classification(sec1["w_purified"], pred_game1)
    
    b_vec1 = m_game1.coef_[:3]
    A_mat1 = m_game1.coef_[3:].reshape((3, 3))
    game_type1, eq1 = solve_subgame_nash(A_mat1[0,0], A_mat1[0,1], A_mat1[1,0], A_mat1[1,1])
    
    print(f"  Sec1 Intrinsic Baseline: R^2 = {r2_int1*100:5.2f}% | AUC = {diag_int1['roc_auc']:.3f} | DOR = {diag_int1['odds_ratio']:.2f}")
    print(f"  Sec1 Spatial Game Assay: R^2 = {r2_game1*100:5.2f}% | AUC = {diag_game1['roc_auc']:.3f} | DOR = {diag_game1['odds_ratio']:.2f}")
    print(f"  Sec1 Game Type: {game_type1} | Equil = {eq1}")
    print(f"  Sec1 Payoff Matrix A:\n{np.round(A_mat1, 3)}")
    
    # =========================================================================
    # STEP 2: INDEPENDENT TRAINING ON SECTION 2 (Internal Ground Truth)
    # =========================================================================
    print("\n[+] STEP 2: Independent Model Training on Section 2 (N = 3,987 spots)...")
    m_int2 = Ridge(alpha=1e-3, fit_intercept=False).fit(X_int2, sec2["w_purified"])
    m_game2 = Ridge(alpha=1e-3, fit_intercept=False).fit(X_game2, sec2["w_purified"])
    
    pred_int2 = m_int2.predict(X_int2)
    pred_game2 = m_game2.predict(X_game2)
    
    r2_int2 = float(m_int2.score(X_int2, sec2["w_purified"]))
    r2_game2 = float(m_game2.score(X_game2, sec2["w_purified"]))
    
    diag_int2 = evaluate_clinical_classification(sec2["w_purified"], pred_int2)
    diag_game2 = evaluate_clinical_classification(sec2["w_purified"], pred_game2)
    
    b_vec2 = m_game2.coef_[:3]
    A_mat2 = m_game2.coef_[3:].reshape((3, 3))
    game_type2, eq2 = solve_subgame_nash(A_mat2[0,0], A_mat2[0,1], A_mat2[1,0], A_mat2[1,1])
    
    print(f"  Sec2 Intrinsic Baseline: R^2 = {r2_int2*100:5.2f}% | AUC = {diag_int2['roc_auc']:.3f} | DOR = {diag_int2['odds_ratio']:.2f}")
    print(f"  Sec2 Spatial Game Assay: R^2 = {r2_game2*100:5.2f}% | AUC = {diag_game2['roc_auc']:.3f} | DOR = {diag_game2['odds_ratio']:.2f}")
    print(f"  Sec2 Game Type: {game_type2} | Equil = {eq2}")
    print(f"  Sec2 Payoff Matrix A:\n{np.round(A_mat2, 3)}")
    
    # =========================================================================
    # STEP 3: BLIND OUT-OF-SAMPLE (OOS) TESTING (Sec1 Model -> Sec2 Data)
    # =========================================================================
    print("\n" + "="*80)
    print("STEP 3: BLIND OUT-OF-SAMPLE (OOS) GENERALIZATION TEST")
    print("Train: Section 1 | Test: Section 2 (Unseen, Held-Out Consecutive Slice)")
    print("="*80)
    
    # Blind predictions
    blind_pred_int2 = m_int1.predict(X_int2)
    blind_pred_game2 = m_game1.predict(X_game2)
    
    # Blind OOS R^2
    ss_tot2 = np.sum((sec2["w_purified"] - np.mean(sec2["w_purified"])) ** 2)
    ss_res_int2 = np.sum((sec2["w_purified"] - blind_pred_int2) ** 2)
    ss_res_game2 = np.sum((sec2["w_purified"] - blind_pred_game2) ** 2)
    
    oos_r2_int = float(1.0 - (ss_res_int2 / ss_tot2))
    oos_r2_game = float(1.0 - (ss_res_game2 / ss_tot2))
    
    # Blind Diagnostic Performance (Luca Celebrano contingency framework)
    blind_diag_int = evaluate_clinical_classification(sec2["w_purified"], blind_pred_int2)
    blind_diag_game = evaluate_clinical_classification(sec2["w_purified"], blind_pred_game2)
    
    print(f"\n[BLIND OOS PREDICTIVE PERFORMANCE ON SECTION 2]")
    print(f"  1. Non-Spatial Intrinsic Baseline (OOS):")
    print(f"     OOS R^2 = {oos_r2_int*100:5.2f}% | ROC-AUC = {blind_diag_int['roc_auc']:.3f}")
    print(f"     DOR     = {blind_diag_int['odds_ratio']:.2f} [95% CI: {blind_diag_int['odds_ratio_95_ci'][0]:.2f} - {blind_diag_int['odds_ratio_95_ci'][1]:.2f}]")
    print(f"     Fisher Exact p = {blind_diag_int['fisher_exact_p_value']:.2e}")
    print(f"     Sensitivity = {blind_diag_int['sensitivity']*100:.1f}% | Specificity = {blind_diag_int['specificity']*100:.1f}%")
    print(f"     PPV = {blind_diag_int['ppv_precision']*100:.1f}% | NPV = {blind_diag_int['npv']*100:.1f}%")
    
    print(f"\n  2. In Situ Spatial Evolutionary Game Assay (OOS):")
    print(f"     OOS R^2 = {oos_r2_game*100:5.2f}% | ROC-AUC = {blind_diag_game['roc_auc']:.3f}")
    print(f"     DOR     = {blind_diag_game['odds_ratio']:.2f} [95% CI: {blind_diag_game['odds_ratio_95_ci'][0]:.2f} - {blind_diag_game['odds_ratio_95_ci'][1]:.2f}]")
    print(f"     Fisher Exact p = {blind_diag_game['fisher_exact_p_value']:.2e}")
    print(f"     Sensitivity = {blind_diag_game['sensitivity']*100:.1f}% | Specificity = {blind_diag_game['specificity']*100:.1f}%")
    print(f"     PPV = {blind_diag_game['ppv_precision']*100:.1f}% | NPV = {blind_diag_game['npv']*100:.1f}%")
    
    delta_oos_r2 = oos_r2_game - oos_r2_int
    delta_oos_auc = blind_diag_game['roc_auc'] - blind_diag_int['roc_auc']
    print(f"\n  [★ OOS ADVANTAGE OVER INTRINSIC BASELINE]")
    print(f"     Delta OOS R^2 = +{delta_oos_r2*100:.2f}% explained variance gain")
    print(f"     Delta OOS AUC = +{delta_oos_auc:.3f} discriminative ability gain")
    
    # =========================================================================
    # STEP 4: PARAMETER CONCORDANCE & REPRODUCIBILITY
    # =========================================================================
    flat_A1 = A_mat1.flatten()
    flat_A2 = A_mat2.flatten()
    
    pearson_r_A, pearson_p_A = pearsonr(flat_A1, flat_A2)
    spearman_rho_A, spearman_p_A = spearmanr(flat_A1, flat_A2)
    mae_A = float(np.mean(np.abs(flat_A1 - flat_A2)))
    
    pearson_r_b, _ = pearsonr(b_vec1, b_vec2)
    
    print("\n" + "="*80)
    print("STEP 4: PARAMETER CONCORDANCE ACROSS SERIAL HISTOLOGICAL CUTS")
    print("="*80)
    print(f"  Payoff Matrix Pearson r      = {pearson_r_A:+.3f} (p = {pearson_p_A:.3e})")
    print(f"  Payoff Matrix Spearman rho   = {spearman_rho_A:+.3f} (p = {spearman_p_A:.3e})")
    print(f"  Payoff Matrix Mean Abs Error = {mae_A:.3f}")
    print(f"  Intrinsic Vector Pearson r   = {pearson_r_b:+.3f}")
    print(f"  Game Type Replicate Status   : Sec1 = {game_type1} | Sec2 = {game_type2}")
    
    # Save comprehensive results
    summary_results = {
        "dataset": "10x Visium Human Breast Cancer Block A Serial Sections",
        "cohorts": {
            "section_1": {"n_spots": sec1["n_spots"], "median_umi": float(np.median(sec1["total_umi"]))},
            "section_2": {"n_spots": sec2["n_spots"], "median_umi": float(np.median(sec2["total_umi"]))}
        },
        "section_1_in_sample": {
            "r2_intrinsic": r2_int1,
            "r2_game": r2_game1,
            "diagnostic_intrinsic": diag_int1,
            "diagnostic_game": diag_game1,
            "payoff_matrix_3x3": A_mat1.tolist(),
            "game_type": game_type1,
            "equilibria": eq1
        },
        "section_2_in_sample": {
            "r2_intrinsic": r2_int2,
            "r2_game": r2_game2,
            "diagnostic_intrinsic": diag_int2,
            "diagnostic_game": diag_game2,
            "payoff_matrix_3x3": A_mat2.tolist(),
            "game_type": game_type2,
            "equilibria": eq2
        },
        "blind_oos_evaluation": {
            "oos_r2_intrinsic": oos_r2_int,
            "oos_r2_game": oos_r2_game,
            "delta_oos_r2": delta_oos_r2,
            "delta_oos_auc": delta_oos_auc,
            "diagnostic_intrinsic_oos": blind_diag_int,
            "diagnostic_game_oos": blind_diag_game
        },
        "parameter_concordance": {
            "payoff_matrix_pearson_r": float(pearson_r_A),
            "payoff_matrix_pearson_p": float(pearson_p_A),
            "payoff_matrix_spearman_rho": float(spearman_rho_A),
            "payoff_matrix_spearman_p": float(spearman_p_A),
            "payoff_matrix_mae": mae_A,
            "intrinsic_vector_pearson_r": float(pearson_r_b)
        }
    }
    
    out_file = os.path.join(RESULTS_DIR, "serial_replicate_blind_oos_validation.json")
    with open(out_file, "w") as fp:
        json.dump(summary_results, fp, indent=2)
    print(f"\n[✓] Blind OOS validation successfully saved to: {out_file}")
    return summary_results


if __name__ == "__main__":
    execute_blind_oos_validation()
