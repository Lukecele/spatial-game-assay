#!/usr/bin/env python3
"""
calculate_spatial_correlation_length.py
=======================================
Calculates the two-point spatial correlation function C(r) for phenotype distributions
in Lung SCC vs Breast IDC to test the Ginzburg Criterion / Mean-Field hypothesis.
"""

import os
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit

def exp_decay(r, A, xi, c):
    return A * np.exp(-r / xi) + c

def compute_correlation_length(cancer_key, n_samples=1500):
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", cancer_key))
    # Read tissue positions
    pos_file = os.path.join(data_dir, "spatial", "tissue_positions_list.csv")
    if not os.path.exists(pos_file):
        pos_file = os.path.join(data_dir, "spatial", "tissue_positions.csv")
    
    df = pd.read_csv(pos_file, header=None if "list" in pos_file else 0)
    if "list" in pos_file:
        df = df[df[1] == 1] # in tissue
        coords = df[[4, 5]].values
    else:
        df = df[df["in_tissue"] == 1]
        coords = df[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values

    # Subsample for fast distance calculation on PC
    np.random.seed(42)
    idx = np.random.choice(len(coords), min(n_samples, len(coords)), replace=False)
    sub_coords = coords[idx]

    # Convert coordinates to microns (approx 100 um pitch corresponds to ~120-150 pixels in standard Visium)
    dists = cdist(sub_coords, sub_coords)
    dists_flat = dists[np.triu_indices_from(dists, k=1)]
    
    # Scale: median distance to nearest neighbor is ~100 um
    nn_dist = np.partition(dists + np.eye(len(dists))*1e6, 1, axis=1)[:, 1]
    pixel_to_um = 100.0 / np.median(nn_dist)
    dists_um = dists_flat * pixel_to_um

    return np.median(nn_dist), pixel_to_um

print("Computing spatial scales...")
med_lung, scale_lung = compute_correlation_length("lung")
med_breast, scale_breast = compute_correlation_length("breast")
print(f"Lung: scale = {scale_lung:.4f} um/px | Nearest neighbor = {med_lung:.1f} px")
print(f"Breast: scale = {scale_breast:.4f} um/px | Nearest neighbor = {med_breast:.1f} px")
