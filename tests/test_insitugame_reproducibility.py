#!/usr/bin/env python3
"""
test_insitugame_reproducibility.py
==================================
Comprehensive Unit, Mathematical, and Empirical Reproducibility Test Suite
for InSituGame (In Situ Spatial Evolutionary Game Assay).
Author: Luca Celebrano (School of Medicine and Surgery, Università degli Studi di Napoli Federico II)

Tests Covered:
1. Exact Mathematical Contingency & Fisher Exact Test (via math.lgamma)
2. Woolf Logit Odds Ratio and 95% Confidence Interval Validation
3. Microenvironmental Decoupling Residual Orthogonality Proof
4. Analytical 2x2 Nash Equilibrium Classification Boundaries
5. Replicator Dynamics and Adaptive Dosing Trajectory Invariants
6. End-to-End Assay Deterministic Reproducibility on Real Clinical Data
"""

import os
import sys
import math
import unittest
import numpy as np

# Ensure src is on path
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(BASE_DIR, "src"))

from insitugame import (
    InSituGameAssay, EnvironmentalDecoupler,
    solve_subgame_nash, simulate_adaptive_therapy,
    ContingencyTable2x2, calculate_odds_ratio,
    calculate_fisher_exact_p_value, compute_roc_auc,
    evaluate_clinical_classification
)


class TestClinicalMetrics(unittest.TestCase):
    def test_contingency_table_structure(self):
        table = ContingencyTable2x2(TP=100, FP=20, FN=30, TN=150)
        self.assertEqual(table.total, 300)
        self.assertEqual(table.to_dict(), {"TP": 100, "FP": 20, "FN": 30, "TN": 150})

    def test_woolf_odds_ratio_ci(self):
        # Table: a=100, b=20, c=30, d=150
        # Haldane-Anscombe correction adds 0.5 to each cell:
        # a=100.5, b=20.5, c=30.5, d=150.5 -> OR = (100.5 * 150.5) / (20.5 * 30.5) = 15125.25 / 625.25 = 24.19
        table = ContingencyTable2x2(TP=100, FP=20, FN=30, TN=150)
        odds_ratio, (ci_low, ci_high) = calculate_odds_ratio(table)
        self.assertAlmostEqual(odds_ratio, 24.19, delta=0.1)
        self.assertTrue(ci_low < odds_ratio < ci_high)
        self.assertTrue(ci_low > 1.0)  # Significantly increased odds

    def test_fisher_exact_lgamma(self):
        # Extreme table should have very small p-value
        table = ContingencyTable2x2(TP=500, FP=50, FN=50, TN=500)
        p_val = calculate_fisher_exact_p_value(table)
        self.assertTrue(p_val < 1e-10)
        self.assertTrue(p_val >= 0.0)

        # Null table (no association) should have p-value close to 1.0
        null_table = ContingencyTable2x2(TP=50, FP=50, FN=50, TN=50)
        null_p = calculate_fisher_exact_p_value(null_table)
        self.assertAlmostEqual(null_p, 1.0, delta=0.05)

    def test_roc_auc_properties(self):
        # Perfect discrimination
        y_true = np.array([1, 1, 1, 0, 0, 0])
        y_score_perf = np.array([0.9, 0.8, 0.7, 0.3, 0.2, 0.1])
        self.assertEqual(compute_roc_auc(y_true, y_score_perf), 1.0)

        # Inverted discrimination
        y_score_inv = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
        self.assertEqual(compute_roc_auc(y_true, y_score_inv), 0.0)


class TestEnvironmentalDecoupling(unittest.TestCase):
    def test_residual_orthogonality(self):
        """Proof of mathematical orthogonality between purified fitness and environmental covariates."""
        np.random.seed(42)
        n = 500
        z_hyp = np.random.randn(n)
        z_vasc = np.random.randn(n)
        # Raw fitness has both microenvironmental and independent signal
        w_raw = 2.5 * z_hyp - 1.8 * z_vasc + np.random.randn(n)

        decoupler = EnvironmentalDecoupler()
        w_purified = decoupler.fit_transform(w_raw, z_hyp, z_vasc)

        # Residuals must have correlation = 0 with z_hyp and z_vasc (within machine precision)
        corr_hyp = np.corrcoef(w_purified, z_hyp)[0, 1]
        corr_vasc = np.corrcoef(w_purified, z_vasc)[0, 1]
        self.assertAlmostEqual(corr_hyp, 0.0, places=7)
        self.assertAlmostEqual(corr_vasc, 0.0, places=7)


class TestNashEquilibriumSolver(unittest.TestCase):
    def test_hawk_dove_snowdrift(self):
        # Negative diagonals, positive off-diagonals -> Hawk-Dove
        game_type, details = solve_subgame_nash(A_SS=-1.0, A_SR=2.5, A_RS=2.0, A_RR=-2.0)
        self.assertEqual(game_type, "Hawk-Dove / Snowdrift (Coexistence)")
        # x* = (A_RR - A_SR) / (A_SS - A_SR - A_RS + A_RR) = (-2.0 - 2.5) / (-1.0 - 2.5 - 2.0 - 2.0) = -4.5 / -7.5 = 0.60
        self.assertAlmostEqual(details["x_star"], 0.60, delta=0.01)

    def test_stag_hunt_bistable(self):
        # A_SS > A_RS and A_RR > A_SR -> Stag Hunt
        game_type, details = solve_subgame_nash(A_SS=1.0, A_SR=0.0, A_RS=0.0, A_RR=1.0)
        self.assertEqual(game_type, "Stag Hunt (Bistable)")
        # delta_S = 1, delta_R = 1 -> x_sep = 1 / (1 + 1) = 0.5
        self.assertAlmostEqual(details["x_sep"], 0.50, delta=0.01)

    def test_sensitive_dominance(self):
        game_type, details = solve_subgame_nash(A_SS=2.0, A_SR=1.0, A_RS=1.0, A_RR=0.0)
        self.assertEqual(game_type, "Sensitive Dominance")
        self.assertEqual(details["x_pure"], 1.0)

    def test_resistant_dominance(self):
        game_type, details = solve_subgame_nash(A_SS=0.0, A_SR=1.0, A_RS=2.0, A_RR=3.0)
        self.assertEqual(game_type, "Resistant Dominance")
        self.assertEqual(details["x_pure"], 0.0)


class TestAdaptiveTherapyDynamics(unittest.TestCase):
    def test_adaptive_therapy_superiority(self):
        """Verifies that adaptive dosing maintains sensitive population better than continuous MTD."""
        # Classic Stag hunt matrix where continuous MTD eradicates sensitive cells causing resistant takeover
        A_2x2 = np.array([
            [0.547, -0.077],
            [-0.044, 0.070]
        ])
        
        sim = simulate_adaptive_therapy(
            A_2x2=A_2x2,
            x0=0.80,
            T_max=15.0,
            dt=0.01,
            threshold=0.30,
            d_S_mtd=1.2,
            d_S_adt=0.7
        )
        
        # MTD rapidly eradicates sensitive cells, causing resistant release
        self.assertTrue(sim["final_sensitive_mtd"] < 0.05)
        # Adaptive therapy sustains sensitive cells comfortably above the separatrix
        self.assertTrue(sim["final_sensitive_adaptive"] > 0.25)


class TestEndToEndClinicalDataReproducibility(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Loads authentic clinical Lung SCC dataset."""
        import gzip
        import pandas as pd
        
        cls.data_dir = os.path.join(BASE_DIR, "data/lung")
        cls.cache_path = os.path.join(BASE_DIR, "data/cache/clinical_pan_cancer_cache.npz")
        
        has_raw = os.path.exists(cls.data_dir)
        has_cache = os.path.exists(cls.cache_path)
        
        if not has_raw and not has_cache:
            raise unittest.SkipTest("Neither clinical data/lung directory nor clinical cache found.")
            
        spec = {
            "T_S": ["EGFR", "KRT5", "TP63", "EPCAM"],
            "T_R": ["VIM", "AXL", "CDH2", "MMP2"],
            "S_CAF": ["ACTA2", "COL1A1", "FAP"],
            "FITNESS": ["MKI67", "PCNA", "TOP2A", "CDK1"],
            "HYPOXIA": ["HIF1A", "VEGFA", "CA9"],
            "VASCULAR": ["PECAM1", "VWF"]
        }
        cls.spec = spec
        all_genes = sorted(list(set(spec["T_S"] + spec["T_R"] + spec["S_CAF"] + spec["FITNESS"] + spec["HYPOXIA"] + spec["VASCULAR"])))
        
        if not has_raw and has_cache:
            cache = np.load(cls.cache_path)
            cls.coords = cache["lung_coords"]
            cls.total_umi = cache["lung_total_umi"]
            cls.barcodes = [f"spot_{i:05d}" for i in range(len(cls.coords))]
            cls.raw_counts = {
                g: (np.expm1(cache[f"lung_gene_{g}"]) * np.maximum(cls.total_umi, 1.0)) / 1e4
                for g in all_genes if f"lung_gene_{g}" in cache
            }
            return
            
        # 1. Barcodes
        with gzip.open(os.path.join(cls.data_dir, "filtered_feature_bc_matrix/barcodes.tsv.gz"), "rt") as fp:
            cls.barcodes = [line.strip() for line in fp]
        n_spots = len(cls.barcodes)
        
        # 2. Features
        gene_to_row = {}
        with gzip.open(os.path.join(cls.data_dir, "filtered_feature_bc_matrix/features.tsv.gz"), "rt") as fp:
            for r_idx, line in enumerate(fp):
                parts = line.strip().split("\t")
                if parts[1] in all_genes:
                    gene_to_row[parts[1]] = r_idx
                    
        # 3. Spatial Coords
        pos_file = os.path.join(cls.data_dir, "spatial/tissue_positions.csv")
        if not os.path.exists(pos_file):
            pos_file = os.path.join(cls.data_dir, "spatial/tissue_positions_list.csv")
            pos_df = pd.read_csv(pos_file, header=None, names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
        else:
            pos_df = pd.read_csv(pos_file)
            
        pos_df = pos_df[pos_df["in_tissue"] == 1].set_index("barcode")
        pos_aligned = pos_df.loc[cls.barcodes]
        cls.coords = pos_aligned[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
        
        # 4. Matrix Stream
        target_rows = set(gene_to_row.values())
        cls.raw_counts = {g: np.zeros(n_spots, dtype=np.float32) for g in gene_to_row}
        cls.total_umi = np.zeros(n_spots, dtype=np.float32)
        
        row_to_gene = {r: g for g, r in gene_to_row.items()}
        with gzip.open(os.path.join(cls.data_dir, "filtered_feature_bc_matrix/matrix.mtx.gz"), "rt") as fp:
            for line in fp:
                if not line.startswith("%"):
                    break
            for line in fp:
                parts = line.split()
                if len(parts) == 3:
                    r_idx = int(parts[0]) - 1
                    s_idx = int(parts[1]) - 1
                    val = float(parts[2])
                    cls.total_umi[s_idx] += val
                    if r_idx in target_rows:
                        cls.raw_counts[row_to_gene[r_idx]][s_idx] = val

    def test_lung_assay_deterministic_reproducibility(self):
        """Verifies exact deterministic reproducibility of InSituGameAssay on clinical Lung SCC data."""
        assay = InSituGameAssay(self.spec)
        assay.fit(self.coords, self.raw_counts, self.total_umi)
        
        summary = assay.summary()
        
        # 1. Check Nash Classification
        self.assertEqual(summary["nash_game_type"], "Stag Hunt (Bistable)")
        
        # 2. Check Payoff Matrix Values (A_SS = +0.547, A_SR = -0.077, A_RS = -0.044, A_RR = +0.070)
        A = np.array(summary["payoff_matrix_3x3"])
        self.assertAlmostEqual(A[0, 0], 0.547, delta=0.01)
        self.assertAlmostEqual(A[0, 1], -0.077, delta=0.01)
        self.assertAlmostEqual(A[1, 0], -0.044, delta=0.01)
        self.assertAlmostEqual(A[1, 1], 0.070, delta=0.01)
        
        # 3. Check Separatrix Tipping Point (x_sep = 0.200)
        self.assertAlmostEqual(summary["nash_equilibria"]["x_sep"], 0.200, delta=0.01)
        
        # 4. Check Explained Variance R^2 (~30.9%)
        self.assertAlmostEqual(summary["explained_variance_r2"], 0.3086, delta=0.01)
        
        # 5. Check Clinical Classification & Diagnostic Odds Ratio
        diag = summary["diagnostic_metrics"]
        self.assertTrue(diag["is_statistically_significant"])
        self.assertTrue(diag["odds_ratio"] > 6.0)
        self.assertTrue(diag["fisher_exact_p_value"] < 1e-10)


if __name__ == "__main__":
    unittest.main()
