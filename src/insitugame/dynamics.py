"""
dynamics.py
===========
Analytical Nash equilibrium solver for 2x2 sub-games and
continuous non-linear replicator dynamical ODE simulations with adaptive therapy.
"""

from typing import Dict, Any, Tuple
import numpy as np


def solve_subgame_nash(A_SS: float, A_SR: float, A_RS: float, A_RR: float) -> Tuple[str, Dict[str, Any]]:
    """
    Solves the exact Nash equilibrium classification for the 2x2 Sensitive vs Resistant sub-game.
    """
    delta_S = A_SS - A_RS
    delta_R = A_RR - A_SR
    
    if delta_S > 0 and delta_R < 0:
        return "Sensitive Dominance", {"x_pure": 1.0, "status": "Sensitive clone unconditionally dominates"}
    elif delta_S < 0 and delta_R > 0:
        return "Resistant Dominance", {"x_pure": 0.0, "status": "Resistant clone unconditionally dominates"}
    elif delta_S > 0 and delta_R > 0:
        x_sep = delta_R / (delta_S + delta_R)
        return "Stag Hunt (Bistable)", {
            "x_sep": float(x_sep),
            "x_equilibria": [0.0, 1.0],
            "status": f"Bistable separatrix tipping point at x = {x_sep:.3f}"
        }
    else:  # delta_S < 0 and delta_R < 0
        denom = (A_SS - A_SR - A_RS + A_RR)
        x_star = (A_RR - A_SR) / denom if abs(denom) > 1e-9 else 0.5
        return "Hawk-Dove / Snowdrift (Coexistence)", {
            "x_star": float(x_star),
            "status": f"Stable internal coexistence at x* = {x_star:.3f}"
        }


def simulate_adaptive_therapy(
    A_2x2: np.ndarray,
    x0: float = 0.8,
    T_max: float = 15.0,
    dt: float = 0.01,
    threshold: float = 0.4,
    d_S_mtd: float = 1.2,
    d_S_adt: float = 0.7
) -> Dict[str, Any]:
    """
    Simulates canonical Replicator Dynamics:
    dx_i / dt = x_i * [ (Ax)_i - x^T A x - drug_i(t) ]
    comparing continuous Maximum Tolerated Dose (MTD) vs Adaptive Therapy.
    """
    n_steps = int(T_max / dt)
    t_span = np.linspace(0, T_max, n_steps)
    
    # Initial vector [x_S, x_R]
    x_init = np.array([x0, 1.0 - x0])
    
    # 1. Continuous MTD Chemotherapy
    x_mtd = np.zeros((n_steps, 2))
    x_mtd[0] = x_init
    drug_mtd = np.array([d_S_mtd, 0.0])
    
    for k in range(n_steps - 1):
        fit_mtd = A_2x2 @ x_mtd[k] - drug_mtd
        avg_fit_mtd = np.dot(x_mtd[k], fit_mtd)
        dx_mtd = x_mtd[k] * (fit_mtd - avg_fit_mtd)
        x_mtd[k + 1] = np.clip(x_mtd[k] + dt * dx_mtd, 1e-4, 1.0)
        x_mtd[k + 1] /= np.sum(x_mtd[k + 1])
        
    # 2. Adaptive Therapy
    x_adt = np.zeros((n_steps, 2))
    x_adt[0] = x_init
    drug_on = np.zeros(n_steps, dtype=bool)
    
    for k in range(n_steps - 1):
        is_drug = (x_adt[k, 0] > threshold)
        drug_on[k] = is_drug
        drug_adt = np.array([d_S_adt, 0.0]) if is_drug else np.array([0.0, 0.0])
        
        fit_adt = A_2x2 @ x_adt[k] - drug_adt
        avg_fit_adt = np.dot(x_adt[k], fit_adt)
        dx_adt = x_adt[k] * (fit_adt - avg_fit_adt)
        x_adt[k + 1] = np.clip(x_adt[k] + dt * dx_adt, 1e-4, 1.0)
        x_adt[k + 1] /= np.sum(x_adt[k + 1])
        
    drug_on[-1] = (x_adt[-1, 0] > threshold)
    
    return {
        "time": t_span.tolist(),
        "x_mtd_sensitive": x_mtd[:, 0].tolist(),
        "x_adaptive_sensitive": x_adt[:, 0].tolist(),
        "drug_active": drug_on.tolist(),
        "final_sensitive_mtd": float(x_mtd[-1, 0]),
        "final_sensitive_adaptive": float(x_adt[-1, 0])
    }
