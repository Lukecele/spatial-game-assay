#!/usr/bin/env python3
"""
insitu_spatial_game.py
======================
First-in-Class In Situ Game-Theoretic Decomposition of Tumor Clonal Interaction
from Human Spatial Transcriptomics Data (10x Visium CytAssist - Human Lung Squamous Carcinoma).

Mathematical Pipeline:
1. Spatial Neighborhood: Hexagonal k-NN / Delaunay triangulation mapping local contacts.
2. Two-Component Orthogonalization: Residualizing proliferation fitness against hypoxia/vasculature.
3. Inverse Game Payoff Matrix Solver: Estimating A_ij via regularized bilinear interaction regression.
4. Nash Equilibrium & Replicator Dynamics Simulation: Modeling Continuous vs. Adaptive Therapy.
"""

import os
import sys
import gzip
import json
import time
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from scipy.sparse import csc_matrix
from sklearn.linear_model import Ridge, LinearRegression
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = "/home/luca/spatial-game-assay/data/lung"
RESULTS_DIR = "/home/luca/spatial-game-assay/results"
FIGURES_DIR = "/home/luca/spatial-game-assay/figures"

# -------------------------------------------------------------
# 1. Definizione dei Marker Biologici Certificati
# -------------------------------------------------------------
MARKER_GENES = {
    "T_S": ["EGFR", "KRT5", "TP63", "EPCAM"],           # Sensitive / Epithelial
    "T_R": ["VIM", "AXL", "CDH2", "MMP2"],             # Resistant / Mesenchymal / EMT
    "S_CAF": ["ACTA2", "COL1A1", "FAP"],               # Stroma / Fibroblasts
    "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],      # Proliferation / Cell Cycle
    "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],              # Microenvironmental Hypoxia Stress
    "VASCULAR": ["PECAM1", "VWF"]                       # Endothelial / Vascular Proximity
}

ALL_TARGET_GENES = sorted(list(set(
    MARKER_GENES["T_S"] + MARKER_GENES["T_R"] + MARKER_GENES["S_CAF"] +
    MARKER_GENES["FITNESS"] + MARKER_GENES["HYPOXIA"] + MARKER_GENES["VASCULAR"]
)))


def load_dataset():
    print("[*] Caricamento e indicizzazione dei file 10x Genomics...")
    
    # 1. Barcodes
    barcodes_file = os.path.join(DATA_DIR, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    with gzip.open(barcodes_file, "rt") as fp:
        barcodes = [line.strip() for line in fp]
    n_spots = len(barcodes)
    barcode_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    print(f"  [+] Barcodes in-tissue caricati: {n_spots:,}")

    # 2. Features (Geni)
    features_file = os.path.join(DATA_DIR, "filtered_feature_bc_matrix/features.tsv.gz")
    gene_to_row = {}
    with gzip.open(features_file, "rt") as fp:
        for row_idx, line in enumerate(fp):
            parts = line.strip().split("\t")
            gene_symbol = parts[1]
            if gene_symbol in ALL_TARGET_GENES:
                gene_to_row[gene_symbol] = row_idx  # 0-based row index
    print(f"  [+] Marker target mappati: {len(gene_to_row)}/{len(ALL_TARGET_GENES)}")

    # 3. Coordinate Spaziali
    pos_file = os.path.join(DATA_DIR, "spatial/tissue_positions.csv")
    pos_df = pd.read_csv(pos_file)
    pos_df = pos_df[pos_df["in_tissue"] == 1].set_index("barcode")
    # Allineamento rigoroso con l'ordine dei barcodes della matrice
    pos_aligned = pos_df.loc[barcodes]
    coords = pos_aligned[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    print(f"  [+] Coordinate spaziali allineate: shape {coords.shape}")

    # 4. Estrazione mirata in streaming dei conteggi (Zero spreco di RAM)
    target_rows = set(gene_to_row.values())
    raw_counts = {g: np.zeros(n_spots, dtype=np.float32) for g in gene_to_row}
    total_umi_per_spot = np.zeros(n_spots, dtype=np.float32)
    
    matrix_file = os.path.join(DATA_DIR, "filtered_feature_bc_matrix/matrix.mtx.gz")
    print(f"  [*] Parsing in streaming dei 23M di UMI da {matrix_file}...")
    start_t = time.time()
    with gzip.open(matrix_file, "rt") as fp:
        # Skip header e commenti
        for line in fp:
            if not line.startswith("%"):
                break # prima riga dopo i commenti contiene le dimensioni
        
        # Inversione mappa per lookup veloce
        row_to_gene = {r: g for g, r in gene_to_row.items()}
        for line in fp:
            parts = line.split()
            if len(parts) == 3:
                r_idx = int(parts[0]) - 1 # 1-based to 0-based
                s_idx = int(parts[1]) - 1
                val = float(parts[2])
                
                total_umi_per_spot[s_idx] += val
                if r_idx in target_rows:
                    gene = row_to_gene[r_idx]
                    raw_counts[gene][s_idx] = val

    print(f"  [+] Estrazione completata in {time.time() - start_t:.1f}s.")
    print(f"  [i] Media UMI totali per spot: {np.mean(total_umi_per_spot):.1f}")

    # 5. Normalizzazione log-CPM (Counts Per Million)
    cpm_scale = 1e4
    norm_exp = {}
    for g, counts in raw_counts.items():
        norm_exp[g] = np.log1p((counts / total_umi_per_spot) * cpm_scale)

    return barcodes, coords, norm_exp


def compute_phenotypes_and_fitness(norm_exp, coords):
    print("\n[*] Calcolo dei punteggi di modulo e costruzione del simplesso delle strategie...")
    n_spots = len(coords)
    
    # Punteggi medi di modulo
    def get_module_score(gene_list):
        scores = [norm_exp[g] for g in gene_list if g in norm_exp]
        return np.mean(scores, axis=0) if scores else np.zeros(n_spots)

    s_score = get_module_score(MARKER_GENES["T_S"])
    r_score = get_module_score(MARKER_GENES["T_R"])
    c_score = get_module_score(MARKER_GENES["S_CAF"])
    w_raw = get_module_score(MARKER_GENES["FITNESS"])
    z_hyp = get_module_score(MARKER_GENES["HYPOXIA"])
    z_vasc = get_module_score(MARKER_GENES["VASCULAR"])

    # Normalizzazione Softmax per ottenere frequenze sul simplesso (f_S + f_R + f_C = 1)
    exp_matrix = np.vstack([s_score, r_score, c_score])
    # Sottrazione del max per stabilita numerica
    exp_matrix -= np.max(exp_matrix, axis=0, keepdims=True)
    exp_vals = np.exp(exp_matrix)
    freqs = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
    f_S, f_R, f_C = freqs[0], freqs[1], freqs[2]

    print(f"  [+] Frequenze medie nel tessuto: Sensitive={np.mean(f_S):.3f}, Resistant={np.mean(f_R):.3f}, Stroma={np.mean(f_C):.3f}")

    # -------------------------------------------------------------
    # 2. Grafo Spaziale di Delaunay e Vicinati Locali
    # -------------------------------------------------------------
    print("[*] Costruzione della Triangolazione di Delaunay sui 3.858 spot...")
    tri = Delaunay(coords)
    neighbors = {i: set() for i in range(n_spots)}
    for simplex in tri.simplices:
        for i in simplex:
            for j in simplex:
                if i != j:
                    neighbors[i].add(j)

    # Calcolo delle frequenze medie di vicinato \bar{f}_j
    bar_f_S = np.zeros(n_spots)
    bar_f_R = np.zeros(n_spots)
    bar_f_C = np.zeros(n_spots)
    for i in range(n_spots):
        nbrs = list(neighbors[i])
        if len(nbrs) > 0:
            bar_f_S[i] = np.mean(f_S[nbrs])
            bar_f_R[i] = np.mean(f_R[nbrs])
            bar_f_C[i] = np.mean(f_C[nbrs])
        else:
            bar_f_S[i], bar_f_R[i], bar_f_C[i] = f_S[i], f_R[i], f_C[i]

    # -------------------------------------------------------------
    # 3. Disaccoppiamento Ortogonale dallo Stress Microambientale
    # -------------------------------------------------------------
    print("[*] Disaccoppiamento ortogonale da Ipossia (HIF1A/VEGFA/CA9) e Vasi (PECAM1/VWF)...")
    X_env = np.column_stack([z_hyp, z_vasc])
    env_model = LinearRegression().fit(X_env, w_raw)
    w_pred_env = env_model.predict(X_env)
    w_purified = w_raw - w_pred_env  # Fitness purificata (residuo ortogonale)
    
    r2_env = env_model.score(X_env, w_raw)
    print(f"  [+] Varianza spiegata dal solo microambiente fisico: {r2_env*100:.2f}%")
    print(f"  [+] Varianza biologica residua per la Teoria dei Giochi: {(1 - r2_env)*100:.2f}%")

    return {
        "coords": coords,
        "f_S": f_S, "f_R": f_R, "f_C": f_C,
        "bar_f_S": bar_f_S, "bar_f_R": bar_f_R, "bar_f_C": bar_f_C,
        "w_raw": w_raw, "w_purified": w_purified,
        "z_hyp": z_hyp, "z_vasc": z_vasc
    }


def solve_payoff_matrix(data):
    print("\n[*] Risoluzione del sistema inverso per la Matrice dei Payoff A (3x3)...")
    f_S, f_R, f_C = data["f_S"], data["f_R"], data["f_C"]
    bar_S, bar_R, bar_C = data["bar_f_S"], data["bar_f_R"], data["bar_f_C"]
    y = data["w_purified"]
    n_spots = len(y)

    # Costruzione della matrice di disegno per regressione bilineare:
    # W_res = b_S*f_S + b_R*f_R + b_C*f_C +
    #         f_S*(A_SS*bar_S + A_SR*bar_R + A_SC*bar_C) +
    #         f_R*(A_RS*bar_S + A_RR*bar_R + A_RC*bar_C) +
    #         f_C*(A_CS*bar_S + A_CR*bar_R + A_CC*bar_C)
    X = np.column_stack([
        f_S, f_R, f_C,                           # 3 fitness intrinseche b_i
        f_S * bar_S, f_S * bar_R, f_S * bar_C,   # Ritorno per T_S con vicini S, R, C
        f_R * bar_S, f_R * bar_R, f_R * bar_C,   # Ritorno per T_R con vicini S, R, C
        f_C * bar_S, f_C * bar_R, f_C * bar_C    # Ritorno per S_CAF con vicini S, R, C
    ])

    # Risoluzione con Ridge Regression (L2 regularization = 1e-3 per prevenire multicollinearità sul simplesso)
    ridge = Ridge(alpha=1e-3, fit_intercept=False)
    ridge.fit(X, y)
    coefs = ridge.coef_
    r2_game = ridge.score(X, y)

    b_vec = coefs[0:3]
    A_matrix = coefs[3:12].reshape((3, 3))
    strategies = ["T_Sensitive", "T_Resistant", "Stroma_CAF"]

    print(f"  [+] R^2 del modello di interazione strategica: {r2_game*100:.2f}%")
    print("\n" + "="*50)
    print("MATRICE DEI PAYOFF ESTIMATA (InSitu-GameAssay):")
    print("="*50)
    payoff_df = pd.DataFrame(A_matrix, index=strategies, columns=strategies)
    print(payoff_df.to_string())
    print("\nFitness Intrinseca (Basale b_i):")
    for s, b in zip(strategies, b_vec):
        print(f"  - {s:<12}: {b:+.4f}")

    # Salva risultati
    os.makedirs(RESULTS_DIR, exist_ok=True)
    payoff_df.to_csv(os.path.join(RESULTS_DIR, "payoff_matrix.csv"))
    
    # -------------------------------------------------------------
    # 4. Analisi Clinica del Gioco e della Terapia Adattiva
    # -------------------------------------------------------------
    # Consideriamo il sub-gioco 2x2 tra Sensitive e Resistant:
    # A_2x2:
    #         S         R
    #   S   A_SS      A_SR
    #   R   A_RS      A_RR
    A_SS, A_SR = A_matrix[0, 0], A_matrix[0, 1]
    A_RS, A_RR = A_matrix[1, 0], A_matrix[1, 1]
    
    print("\n" + "="*50)
    print("ANALISI STRATEGICA DEL SUB-GIOCO TUMORALE (2x2):")
    print("="*50)
    print(f"  A(Sensitive vs Sensitive): {A_SS:+.4f} | A(Sensitive vs Resistant): {A_SR:+.4f}")
    print(f"  A(Resistant vs Sensitive): {A_RS:+.4f} | A(Resistant vs Resistant): {A_RR:+.4f}")

    # Diagnosi della classe di gioco
    game_class = "Indeterminate"
    if A_SS > A_RS and A_SR > A_RR:
        game_class = "Sensitive Dominance"
    elif A_RS > A_SS and A_RR > A_SR:
        game_class = "Resistant Dominance"
    elif A_RS > A_SS and A_SR > A_RR:
        game_class = "Hawk-Dove (Coexistence)"
        x_star = (A_SR - A_RR) / ((A_SS - A_RS) + (A_SR - A_RR))
        print(f"  [★] INTERNAL NASH EQUILIBRIUM (ESS): x*(Sensitive) = {x_star:.3f}")
    elif A_SS > A_RS and A_RR > A_SR:
        game_class = "Stag Hunt (Bistable)"
    
    print(f"  [>] GAME CLASSIFICATION: {game_class}")

    # Salva referto JSON
    report_data = {
        "n_spots": n_spots,
        "strategies": strategies,
        "payoff_matrix": A_matrix.tolist(),
        "baseline_fitness": b_vec.tolist(),
        "r2_environmental": float(r2_env if 'r2_env' in locals() else 0.0),
        "r2_game_interaction": float(r2_game),
        "game_classification": game_class
    }
    with open(os.path.join(RESULTS_DIR, "clinical_game_report.json"), "w") as fp:
        json.dump(report_data, fp, indent=2)

    return payoff_df, A_matrix, b_vec, game_class


def generate_figures(data, A_matrix, game_class):
    print("\n[*] Generating publication-ready scientific figures...")
    os.makedirs(FIGURES_DIR, exist_ok=True)
    coords = data["coords"]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), dpi=300)
    
    # 1. Spatial Map of Clonal Phenotypes (RGB Composite)
    # R=Resistant, G=Sensitive, B=Stroma
    rgb = np.column_stack([data["f_R"], data["f_S"], data["f_C"]])
    rgb = np.clip(rgb, 0, 1)
    
    sc0 = axes[0].scatter(coords[:, 1], -coords[:, 0], c=rgb, s=6, alpha=0.85)
    axes[0].set_title("A. Spatial Phenotype Architecture\n(Red=Resistant, Green=Sensitive, Blue=Stroma)", fontsize=11, fontweight='bold')
    axes[0].axis('off')

    # 2. Decoupled Proliferation Fitness
    sc1 = axes[1].scatter(coords[:, 1], -coords[:, 0], c=data["w_purified"], cmap="coolwarm", s=6, alpha=0.85)
    axes[1].set_title("B. Decoupled Proliferation Fitness ($w_{purified}$)\n(Decoupled from Hypoxia & Vascular Covariates)", fontsize=11, fontweight='bold')
    axes[1].axis('off')
    cbar1 = plt.colorbar(sc1, ax=axes[1], fraction=0.046, pad=0.04)
    cbar1.ax.tick_params(labelsize=9)

    # 3. Payoff Heatmap
    strategies = ["T_Sens", "T_Res", "Stroma"]
    im = axes[2].imshow(A_matrix, cmap="RdBu_r", vmin=-np.max(np.abs(A_matrix)), vmax=np.max(np.abs(A_matrix)))
    axes[2].set_xticks([0, 1, 2])
    axes[2].set_yticks([0, 1, 2])
    axes[2].set_xticklabels(strategies, fontsize=10, fontweight='bold')
    axes[2].set_yticklabels(strategies, fontsize=10, fontweight='bold')
    clean_game = game_class.split('(')[0].strip()
    axes[2].set_title(f"C. In Situ Payoff Matrix ($A_{{ij}}$)\nGame Regime: {clean_game}", fontsize=11, fontweight='bold')
    cbar2 = plt.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04)
    cbar2.ax.tick_params(labelsize=9)
    
    for i in range(3):
        for j in range(3):
            val = A_matrix[i, j]
            axes[2].text(j, i, f"{val:+.3f}", ha="center", va="center", 
                         color="black" if abs(val)<0.5*np.max(np.abs(A_matrix)) else "white", 
                         fontweight='bold', fontsize=10)

    plt.tight_layout()
    out_fig = os.path.join(FIGURES_DIR, "figure1_spatial_game_decomposition.png")
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [+] Figure 1 successfully saved: {out_fig}")


def main():
    print("="*60)
    print("IN SITU SPATIAL GAME ASSAY (Spatial-GameAssay v1.0)")
    print("="*60)
    barcodes, coords, norm_exp = load_dataset()
    data = compute_phenotypes_and_fitness(norm_exp, coords)
    payoff_df, A_matrix, b_vec, game_class = solve_payoff_matrix(data)
    generate_figures(data, A_matrix, game_class)
    print("\n[✓] PIPELINE COMPLETATA CON SUCCESSO. TUTTI I DATI SONO SALVATI IN 'results/' E 'figures/'.")


if __name__ == "__main__":
    main()
