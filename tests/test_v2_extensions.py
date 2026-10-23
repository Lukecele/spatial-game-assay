#!/usr/bin/env python3
"""
test_v2_extensions.py
=====================
Comprehensive Unit & Integration Test Suite for InSituGame v2 Extensions:
1. SubCellularMeshAdapter (Xenium / Single-Cell Delaunay Mesh)
2. SpatialSomaticTracer (RNA BAM VAF Tracking & Correlation)
3. PharmacogenomicReplicatorEngine (GDSC2 Hill Kinetics & PK/PD Replicator Dynamics)
"""

import unittest
import numpy as np
from insitugame import (
    SubCellularMeshAdapter,
    SpatialSomaticTracer,
    LocusPileup,
    PharmacogenomicReplicatorEngine,
    DrugPKPDProfile,
    PRESET_PROFILES
)


class TestSubCellularMeshAdapter(unittest.TestCase):
    def test_xenium_mesh_topology(self):
        np.random.seed(42)
        n_cells = 120
        centroids = np.random.uniform(0, 500, size=(n_cells, 2))  # 500 um tissue patch

        adapter = SubCellularMeshAdapter(platform_name="10x_Xenium").fit_mesh(centroids)
        summary = adapter.get_summary()

        self.assertEqual(summary["total_cells"], 120)
        self.assertGreater(summary["mean_contact_neighbors"], 3.0)
        self.assertLess(summary["mean_contact_neighbors"], 10.0)

        # Phenotype composition: 3 classes
        cell_types = np.random.choice([0, 1, 2], size=n_cells)
        f_cell, bar_f = adapter.compute_local_composition(cell_types, n_types=3)

        self.assertEqual(f_cell.shape, (120, 3))
        self.assertEqual(bar_f.shape, (120, 3))
        np.testing.assert_allclose(np.sum(f_cell, axis=1), 1.0)
        np.testing.assert_allclose(np.sum(bar_f, axis=1), 1.0)


class TestSpatialSomaticTracer(unittest.TestCase):
    def test_vaf_calculation_and_concordance(self):
        tracer = SpatialSomaticTracer()
        barcodes = [f"spot_{i}" for i in range(100)]

        # Simulate true biological gradient of resistant clone f_R
        np.random.seed(1337)
        f_R_true = np.linspace(0.05, 0.95, 100)

        vaf_dict = tracer.ingest_mock_pileup(
            barcodes=barcodes,
            f_R_true=f_R_true,
            variant_key="PIK3CA_H1047R",
            mean_depth=30,
            noise_std=0.04
        )

        self.assertEqual(len(vaf_dict), 100)
        res = tracer.correlate_with_phenotype(barcodes, f_R_true, "PIK3CA_H1047R")

        self.assertEqual(res["concordance_verdict"], "High Genomic-Phenotypic Concordance")
        self.assertGreater(res["pearson_r"], 0.85)
        self.assertLess(res["pearson_p"], 1e-15)


class TestPharmacogenomicReplicatorEngine(unittest.TestCase):
    def test_hill_kill_kinetics(self):
        profile = PRESET_PROFILES["breast_paclitaxel"]
        engine = PharmacogenomicReplicatorEngine(profile)

        # Zero concentration -> Zero kill rate
        d_S_0, d_R_0 = engine.compute_hill_kill_rate(0.0)
        self.assertEqual(d_S_0, 0.0)
        self.assertEqual(d_R_0, 0.0)

        # High concentration -> Approaches d_max
        d_S_inf, d_R_inf = engine.compute_hill_kill_rate(100.0)
        self.assertAlmostEqual(d_S_inf, profile.d_max_per_day, delta=0.05)
        self.assertAlmostEqual(d_R_inf, profile.d_max_per_day, delta=0.05)

        # Therapeutic concentration -> Selective differential d_S > d_R
        d_S_med, d_R_med = engine.compute_hill_kill_rate(profile.ic50_S_uM * 2.0)
        self.assertGreater(d_S_med, d_R_med)

    def test_adaptive_vs_continuous_pkpd_simulation(self):
        profile = PRESET_PROFILES["melanoma_dabrafenib"]
        engine = PharmacogenomicReplicatorEngine(profile)

        # Stag Hunt parameters in Melanoma
        delta_S = 1.633
        delta_R = 1.515

        res_mtd = engine.simulate_schedule(
            delta_S=delta_S,
            delta_R=delta_R,
            x0_S=0.65,
            mode="continuous_mtd",
            t_max_days=90.0
        )

        res_adapt = engine.simulate_schedule(
            delta_S=delta_S,
            delta_R=delta_R,
            x0_S=0.65,
            mode="adaptive",
            t_max_days=90.0,
            adaptive_threshold_R=0.45
        )

        # Continuous MTD causes competitive release -> Fast resistance progression
        # Adaptive therapy halts drug when resistance grows -> Superior TTP
        self.assertGreaterEqual(res_adapt["time_to_progression_days"], res_mtd["time_to_progression_days"])
        self.assertLess(res_adapt["total_doses_administered"], res_mtd["total_doses_administered"])


if __name__ == "__main__":
    unittest.main()
