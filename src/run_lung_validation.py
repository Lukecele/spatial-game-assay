#!/usr/bin/env python3
"""
run_lung_validation.py
======================
Executes the statistical validation suite on the Human Lung SCC dataset.
"""

import os
import sys
import gzip
import json
import time
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from sklearn.linear_model import LinearRegression

# Import validation functions
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from validation_battery import run_full_validation_battery

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "lung"))
RESULTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results"))

MARKER_GENES = {
    "T_S": ["EGFR", "KRT5", "TP63", "EPCAM"],
    "T_R": ["VIM", "AXL", "CDH2", "MMP2"],
    "S_CAF": ["ACTA2", "COL1A1", "FAP"],
    "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
    "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
    "VASCULAR": ["PECAM1", "VWF"]
}

ALL_TARGET_GENES = sorted(list(set(
    MARKER_GENES["T_S"] + MARKER_GENES["T_R"] + MARKER_GENES["S_CAF"] +
    MARKER_GENES["FITNESS"] + MARKER_GENES["HYPOXIA"] + MARKER_GENES["VASCULAR"]
)))


def load_lung_data():
    barcodes_file = os.path.join(DATA_DIR, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    with gzip.open(barcodes_file, "rt") as fp:
        barcodes = [line.strip() for line in fp]
    n_spots = len(barcodes)

    features_file = os.path.join(DATA_DIR, "filtered_feature_bc_matrix/features.tsv.gz")
    gene_to_row = {}
    with gzip.open(features_file, "rt") as fp:
        for row_idx, line in enumerate(fp):
            parts = line.strip().split("\t")
            if parts[1] in ALL_TARGET_GENES:
                gene_to_row[parts[1]] = row_idx

    pos_file = os.path.join(DATA_DIR, "spatial/tissue_positions.csv")
    pos_df = pd.read_csv(pos_file)
    pos_df = pos_df[pos_df["in_tissue"] == 1].set_index("barcode")
    pos_aligned = pos_df.loc[barcodes]
    coords = pos_aligned[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values

    target_rows = set(gene_to_row.values())
    raw_counts = {g: np.zeros(n_spots, dtype=np.float32) for g in gene_to_row}
    total_umi = np.zeros(n_spots, dtype=np.float32)

    matrix_file = os.path.join(DATA_DIR, "filtered_feature_bc_matrix/matrix.mtx.gz")
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

    cpm_scale = 1e4
    norm_exp = {g: np.log1p((c / total_umi) * cpm_scale) for g, c in raw_counts.items()}

    def get_module(genes):
        s = [norm_exp[g] for g in genes if g in norm_exp]
        return np.mean(s, axis=0) if s else np.zeros(n_spots)

    s_score = get_module(MARKER_GENES["T_S"])
    r_score = get_module(MARKER_GENES["T_R"])
    c_score = get_module(MARKER_GENES["S_CAF"])
    w_raw = get_module(MARKER_GENES["FITNESS"])
    z_hyp = get_module(MARKER_GENES["HYPOXIA"])
    z_vasc = get_module(MARKER_GENES["VASCULAR"])

    exp_matrix = np.vstack([s_score, r_score, c_score])
    exp_matrix -= np.max(exp_matrix, axis=0, keepdims=True)
    exp_vals = np.exp(exp_matrix)
    freqs = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
    f_S, f_R, f_C = freqs[0], freqs[1], freqs[2]

    tri = Delaunay(coords)
    neighbors = {i: set() for i in range(n_spots)}
    for simplex in tri.simplices:
        for i in simplex:
            for j in simplex:
                if i != j:
                    neighbors[i].add(j)

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

    X_env = np.column_stack([z_hyp, z_vasc])
    env_model = LinearRegression().fit(X_env, w_raw)
    w_purified = w_raw - env_model.predict(X_env)

    return {
        "coords": coords,
        "f_S": f_S, "f_R": f_R, "f_C": f_C,
        "bar_f_S": bar_f_S, "bar_f_R": bar_f_R, "bar_f_C": bar_f_C,
        "w_raw": w_raw, "w_purified": w_purified,
        "z_hyp": z_hyp, "z_vasc": z_vasc
    }


if __name__ == "__main__":
    t0 = time.time()
    data = load_lung_data()
    print(f"Data loaded in {time.time()-t0:.1f}s. Running validation battery...")
    report = run_full_validation_battery(data, sample_name="Human_Lung_SCC", n_perm=1000, n_boot=1000)
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_json = os.path.join(RESULTS_DIR, "validation_report_lung.json")
    with open(out_json, "w") as fp:
        json.dump(report, fp, indent=2)
    print(f"\n[✓] Validation battery completed and saved to {out_json}")
