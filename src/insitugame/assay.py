"""
assay.py
========
Core InSituGameAssay class for end-to-end spatial evolutionary game quantification.
"""

import os
import gzip
from typing import Dict, List, Any, Optional
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from sklearn.linear_model import Ridge, LinearRegression

from .decoupling import EnvironmentalDecoupler
from .dynamics import solve_subgame_nash, simulate_adaptive_therapy
from .metrics import evaluate_clinical_classification


class InSituGameAssay:
    def __init__(self, gene_spec: Dict[str, List[str]]):
        self.spec = gene_spec
        self.decoupler = EnvironmentalDecoupler()
        self.fitted = False
        self.payoff_matrix = None
        self.intrinsic_vector = None
        self.nash_classification = None
        self.nash_details = None
        
    def fit(self, coords: np.ndarray, raw_counts: Dict[str, np.ndarray], total_umi: np.ndarray):
        """
        Fit the InSituGame assay on spatial transcriptomics data.

        Implements the full pipeline:
        1. Log-CPM normalization and softmax phenotypic frequency projection onto the frequency simplex.
        2. Delaunay spatial neighborhood graph construction.
        3. Microenvironmental decoupling via Frisch-Waugh-Lovell orthogonalization (hypoxia + vasculature).
        4. Estimation of the 9-parameter absorbed payoff matrix A_tilde (b_i + A_ij) via Ridge regression
           on the 12-column design matrix [f_i, f_i * bar_f_j]. The 3 intrinsic regressors f_i are
           algebraically collinear with the 9 bilinear interaction terms under the simplex constraint
           (sum bar_f_j = 1), so the absorbed parameterization A_tilde_ij = b_i + A_ij is the
           identifiable quantity. Ridge(alpha=1e-3) provides numerical stability without materially
           altering the estimated regime classification.
        5. Nash equilibrium classification of the 2x2 (T_S vs T_R) subgame.

        Parameters
        ----------
        coords : np.ndarray, shape (n_spots, 2)
            Pixel coordinates (row, col) for each spatial spot.
        raw_counts : dict of str -> np.ndarray, shape (n_spots,)
            Raw UMI counts per gene per spot.
        total_umi : np.ndarray, shape (n_spots,)
            Total UMI count per spot for library-size normalization.

        Returns
        -------
        self : InSituGameAssay
            Fitted assay instance.
        """
        n_spots = len(coords)
        norm_exp = {g: np.log1p((raw_counts[g] / np.maximum(total_umi, 1.0)) * 1e4) for g in raw_counts}
        
        def get_mod(genes):
            s = [norm_exp[g] for g in genes if g in norm_exp]
            return np.mean(s, axis=0) if s else np.zeros(n_spots)
            
        s_score = get_mod(self.spec["T_S"])
        r_score = get_mod(self.spec["T_R"])
        c_score = get_mod(self.spec["S_CAF"])
        w_raw = get_mod(self.spec["FITNESS"])
        z_hyp = get_mod(self.spec["HYPOXIA"])
        z_vasc = get_mod(self.spec["VASCULAR"])
        
        # Softmax phenotypic frequencies
        exp_mat = np.vstack([s_score, r_score, c_score])
        exp_mat -= np.max(exp_mat, axis=0, keepdims=True)
        exp_vals = np.exp(exp_mat)
        freqs = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
        f_S, f_R, f_C = freqs[0], freqs[1], freqs[2]
        
        # Spatial Delaunay contact graph
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
                
        # Environmental decoupling
        w_purified = self.decoupler.fit_transform(w_raw, z_hyp, z_vasc)
        
        # Game design matrix
        X_game = np.column_stack([
            f_S, f_R, f_C,
            f_S * bar_f_S, f_S * bar_f_R, f_S * bar_f_C,
            f_R * bar_f_S, f_R * bar_f_R, f_R * bar_f_C,
            f_C * bar_f_S, f_C * bar_f_R, f_C * bar_f_C
        ])
        
        # Ridge regression
        model = Ridge(alpha=1e-3, fit_intercept=False).fit(X_game, w_purified)
        coefs = model.coef_
        
        self.intrinsic_vector = coefs[:3]
        self.payoff_matrix = coefs[3:].reshape((3, 3))
        
        A = self.payoff_matrix
        self.nash_classification, self.nash_details = solve_subgame_nash(A[0, 0], A[0, 1], A[1, 0], A[1, 1])
        self.fitted = True
        
        y_pred = model.predict(X_game)
        self.diagnostic_metrics = evaluate_clinical_classification(w_purified, y_pred)
        self.r2 = float(model.score(X_game, w_purified))
        
        return self
        
    def summary(self) -> Dict[str, Any]:
        """
        Return a structured summary of the fitted assay results.

        Returns
        -------
        dict with keys:
            nash_game_type : str
                Classified game regime (e.g. 'Hawk-Dove / Coexistence', 'Stag Hunt / Bistable',
                'Dominant Strategy', 'Indeterminate').
            nash_equilibria : dict
                Interior equilibrium x* (Hawk-Dove) or separatrix x_sep (Stag Hunt), where x
                denotes the fraction of sensitive clones T_S within the T_S + T_R subspace.
            payoff_matrix_3x3 : list of list
                Inferred 3x3 absorbed payoff matrix A_tilde (rows = focal clone, cols = neighbor clone).
                Note: these are absorbed parameters A_tilde_ij = b_i + A_ij, not raw A_ij.
            intrinsic_fitness_vector : list
                Intrinsic fitness coefficients b_i for each phenotype [T_S, T_R, S_CAF].
            explained_variance_r2 : float
                OLS R^2 of the full game model on w_purified.
            diagnostic_metrics : dict
                Clinical classification metrics (AUC, sensitivity, specificity, DOR, Fisher p).

        Raises
        ------
        RuntimeError
            If called before fit().
        """
        if not self.fitted:
            raise RuntimeError("Assay has not been fitted yet.")
        return {
            "nash_game_type": self.nash_classification,
            "nash_equilibria": self.nash_details,
            "payoff_matrix_3x3": self.payoff_matrix.tolist(),
            "intrinsic_fitness_vector": self.intrinsic_vector.tolist(),
            "explained_variance_r2": self.r2,
            "diagnostic_metrics": self.diagnostic_metrics
        }
