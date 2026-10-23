#!/usr/bin/env python3
"""
pharmacogenomic_pkpd.py
=======================
Pharmacokinetic / Pharmacodynamic (PK/PD) Hill-Kinetics Simulator for InSituGame.
Calibrates continuous anti-cancer drug exposure and subclonal kill rates directly from
Sanger GDSC2 and Broad DepMap pharmacogenomic dose-response benchmarks.

Features:
1. One-compartment open pharmacokinetic model with first-order elimination (k_e = ln(2)/t_half).
2. Hill-Emax pharmacodynamic cell kill kinetics for Sensitive (T_S) and Resistant (T_R) clones.
3. Coupled Replicator Dynamics ODE solver under realistic drug decay and multi-dose pulsing.
4. Quantitative Time-to-Progression (TTP) benchmark comparing MTD vs Adaptive Therapy.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Any
import numpy as np


@dataclass
class DrugPKPDProfile:
    drug_name: str
    cancer_type: str
    ic50_S_uM: float
    ic50_R_uM: float
    hill_coefficient: float = 1.0
    d_max_per_day: float = 1.5
    half_life_hours: float = 24.0
    c_peak_uM: float = 0.10


# Curated Sanger GDSC2 / Clinical PK Profiles
PRESET_PROFILES = {
    "breast_paclitaxel": DrugPKPDProfile(
        drug_name="Paclitaxel",
        cancer_type="Breast Invasive Carcinoma",
        ic50_S_uM=0.0084,
        ic50_R_uM=0.0420,
        hill_coefficient=1.2,
        d_max_per_day=1.4,
        half_life_hours=20.0,
        c_peak_uM=0.035
    ),
    "melanoma_dabrafenib": DrugPKPDProfile(
        drug_name="Dabrafenib",
        cancer_type="Cutaneous Melanoma",
        ic50_S_uM=0.0025,
        ic50_R_uM=0.0850,
        hill_coefficient=1.5,
        d_max_per_day=1.8,
        half_life_hours=8.0,
        c_peak_uM=0.040
    ),
    "lung_cisplatin": DrugPKPDProfile(
        drug_name="Cisplatin",
        cancer_type="Lung Squamous Cell Carcinoma",
        ic50_S_uM=1.85,
        ic50_R_uM=5.20,
        hill_coefficient=1.0,
        d_max_per_day=1.2,
        half_life_hours=36.0,
        c_peak_uM=4.00
    )
}


class PharmacogenomicReplicatorEngine:
    """
    Couples empirical pharmacokinetics and pharmacodynamics with in situ replicator dynamics.
    """

    def __init__(self, profile: DrugPKPDProfile):
        self.profile = profile
        self.k_e = np.log(2.0) / (profile.half_life_hours / 24.0)  # Elimination rate in day^-1

    def compute_hill_kill_rate(self, conc_uM: float) -> Tuple[float, float]:
        """Computes instantaneous kill rates d_S(C) and d_R(C) via Hill kinetics."""
        n = self.profile.hill_coefficient
        d_max = self.profile.d_max_per_day

        c_n = conc_uM ** n
        ic50_S_n = self.profile.ic50_S_uM ** n
        ic50_R_n = self.profile.ic50_R_uM ** n

        d_S = float(d_max * (c_n / (c_n + ic50_S_n + 1e-18)))
        d_R = float(d_max * (c_n / (c_n + ic50_R_n + 1e-18)))
        return d_S, d_R

    def simulate_schedule(
        self,
        delta_S: float,
        delta_R: float,
        x0_S: float = 0.60,
        mode: str = "adaptive",
        t_max_days: float = 120.0,
        dt_days: float = 0.1,
        dose_interval_days: float = 7.0,
        adaptive_threshold_R: float = 0.50
    ) -> Dict[str, Any]:
        """
        Simulates forward replicator dynamics under PK/PD drug administration.
        
        mode: 'continuous_mtd' or 'adaptive'
        """
        n_steps = int(t_max_days / dt_days)
        time_grid = np.linspace(0, t_max_days, n_steps)

        x_S = np.zeros(n_steps)
        x_R = np.zeros(n_steps)
        conc = np.zeros(n_steps)
        drug_on = np.zeros(n_steps, dtype=bool)

        x_S[0] = x0_S
        x_R[0] = 1.0 - x0_S
        conc[0] = self.profile.c_peak_uM
        drug_on[0] = True

        last_dose_time = 0.0

        for t_idx in range(1, n_steps):
            t_curr = time_grid[t_idx]

            # PK Clearance decay: dC/dt = -k_e * C
            conc[t_idx] = conc[t_idx - 1] * np.exp(-self.k_e * dt_days)

            # Dosing logic
            if (t_curr - last_dose_time) >= dose_interval_days:
                curr_R = x_R[t_idx - 1]
                if mode == "continuous_mtd":
                    conc[t_idx] += self.profile.c_peak_uM
                    last_dose_time = t_curr
                    drug_on[t_idx] = True
                elif mode == "adaptive":
                    # Administer only if resistant fraction is below safety threshold
                    if curr_R < adaptive_threshold_R:
                        conc[t_idx] += self.profile.c_peak_uM
                        last_dose_time = t_curr
                        drug_on[t_idx] = True

            # PD Hill Kill Rates
            d_S, d_R = self.compute_hill_kill_rate(conc[t_idx])
            delta_kill = d_S - d_R

            # In Situ Replicator Dynamics:
            # dx_S/dt = x_S * (1 - x_S) * [ x_S * (Delta_S + Delta_R) - Delta_R - delta_kill ]
            xs_prev = x_S[t_idx - 1]
            endogenous_grad = xs_prev * (delta_S + delta_R) - delta_R
            dx_dt = xs_prev * (1.0 - xs_prev) * (endogenous_grad - delta_kill)

            # Euler integration step
            xs_next = float(np.clip(xs_prev + dx_dt * dt_days, 1e-4, 1.0 - 1e-4))
            x_S[t_idx] = xs_next
            x_R[t_idx] = 1.0 - xs_next

        # Time to Treatment Progression (TTP): time until resistance >= 90%
        progression_indices = np.where(x_R >= 0.90)[0]
        ttp_days = float(time_grid[progression_indices[0]]) if len(progression_indices) > 0 else float(t_max_days)

        return {
            "mode": mode,
            "drug": self.profile.drug_name,
            "cancer": self.profile.cancer_type,
            "t_max_days": t_max_days,
            "time_to_progression_days": ttp_days,
            "final_sensitive_pct": float(x_S[-1] * 100.0),
            "final_resistant_pct": float(x_R[-1] * 100.0),
            "total_doses_administered": int(np.sum(drug_on))
        }
