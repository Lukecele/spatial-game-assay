#!/usr/bin/env python3
"""
test_modern_mathematics_battery.py
==================================
Unit and Integration Tests for Modern Mathematical Physics & Specificity Modules:
1. Sinkhorn Optimal Transport Algorithmic Invariants (Entropy-regularized Monge-Kantorovich)
2. Agent-Based Lattice Ground Truth Recovery (Condition A TPR, Condition B FPR <= 5%)
3. Biological Specificity Null Battery Integrity (Empirical p-value < 0.05)
4. Geodesic Mass Transport Alignment Validation (|r| > 0.60, p < 1e-100)
"""

import os
import json
import unittest
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
SRC_DIR = os.path.join(BASE_DIR, "src")

import sys
sys.path.append(SRC_DIR)
from optimal_transport_serial_sections import sinkhorn_optimal_transport


class TestSinkhornOptimalTransport(unittest.TestCase):
    def test_sinkhorn_mass_conservation_and_doubly_stochastic(self):
        """Verifies that the Sinkhorn algorithm satisfies exact marginal constraints."""
        np.random.seed(42)
        n, m = 50, 40
        X = np.random.randn(n, 2)
        Y = np.random.randn(m, 2)
        
        a = np.ones(n) / n
        b = np.ones(m) / m
        
        Gamma, C = sinkhorn_optimal_transport(X, Y, a=a, b=b, epsilon=0.05, max_iter=500, tol=1e-6)
        
        # Check row sums equal marginal 'a'
        row_sums = np.sum(Gamma, axis=1)
        np.testing.assert_allclose(row_sums, a, atol=1e-4)
        
        # Check column sums equal marginal 'b'
        col_sums = np.sum(Gamma, axis=0)
        np.testing.assert_allclose(col_sums, b, atol=1e-4)
        
        # Check non-negativity
        self.assertTrue(np.all(Gamma >= 0))


class TestModernReportsIntegrity(unittest.TestCase):
    def test_wasserstein_report(self):
        rep_file = os.path.join(RESULTS_DIR, "WASSERSTEIN_OPTIMAL_TRANSPORT_REPORT.json")
        self.assertTrue(os.path.exists(rep_file))
        with open(rep_file) as fp:
            data = json.load(fp)
            
        self.assertEqual(data["n_spots_sec1"], 3798)
        self.assertEqual(data["n_spots_sec2"], 3987)
        self.assertLess(data["pearson_p"], 1e-100)
        self.assertGreater(abs(data["pearson_r"]), 0.60)
        self.assertEqual(data["concordance_verdict"], "Statistically Significant Geodesic Alignment")

    def test_agent_based_lattice_report(self):
        rep_file = os.path.join(RESULTS_DIR, "AGENT_BASED_LATTICE_SPECIFICITY_REPORT.json")
        self.assertTrue(os.path.exists(rep_file))
        with open(rep_file) as fp:
            data = json.load(fp)
            
        # Condition A Sensitivity
        self.assertEqual(data["condition_A_true_game"]["tpr_sensitivity_pct"], 100.0)
        self.assertLess(data["condition_A_true_game"]["matrix_relative_recovery_error_pct"], 20.0)
        
        # Condition B False Positive Rate Control
        self.assertLessEqual(data["condition_B_pure_crowding_null"]["fpr_specificity_pct"], 5.0)

    def test_random_gene_null_report(self):
        rep_file = os.path.join(RESULTS_DIR, "RANDOM_GENE_NULL_REPORT.json")
        self.assertTrue(os.path.exists(rep_file))
        with open(rep_file) as fp:
            data = json.load(fp)
            
        self.assertEqual(data["n_permutations"], 500)
        self.assertLess(data["p_empirical_f"], 0.05)
        self.assertLess(data["p_empirical_r2"], 0.05)
        self.assertGreater(data["real_f_stat"], data["null_f_mean"] * 3.0)


if __name__ == "__main__":
    unittest.main()
