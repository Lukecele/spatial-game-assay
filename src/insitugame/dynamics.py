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
    A_mat: np.ndarray,
    b_vec: np.ndarray,
    x0: float = 0.8,
    T_max: float = 120.0,
    dt: float = 0.05,
    threshold: float = 0.5,
    d_S: float = 1.2,
    d_R: float = 0.7
) -> Dict[str, Any]:
    """
    Simulates Replicator Dynamics comparing:
    1. Maximum Tolerated Dose (MTD, continuous drug)
    2. Adaptive Therapy (pulsed drug when sensitive fraction > threshold)
    """
    time_steps = int(T_max / dt)
    t_span = np.linspace(0, T_max, time_steps)
    
    # 1. Continuous MTD
    x_mtd = np.zeros(time_steps)
    x_mtd[0] = x0
    for t in range(time_steps - 1):
        x = x_mtd[t]
        y = 1.0 - x
        f_S = b_vec[0] + A_mat[0, 0] * x + A_mat[0, 1] * y - d_S
        f_R = b_vec[1] + A_mat[1, 0] * x + A_mat[1, 1] * y - d_R
        phi = x * f_S + y * f_R
        dx = x * (f_S - phi) * dt
        x_mtd[t + 1] = np.clip(x + dx, 1e-4, 1.0 - 1e-4)
        
    # 2. Adaptive Therapy
    x_adp = np.zeros(time_steps)
    x_adp[0] = x0
    drug_on = np.zeros(time_steps, dtype=bool)
    
    for t in range(time_steps - 1):
        x = x_adp[t]
        y = 1.0 - x
        is_drug = (x > threshold)
        drug_on[t] = is_drug
        
        kill_S = d_S if is_drug else 0.0
        kill_R = d_R if is_drug else 0.0
        
        f_S = b_vec[0] + A_mat[0, 0] * x + A_mat[0, 1] * y - kill_S
        f_R = b_vec[1] + A_mat[1, 0] * x + A_mat[1, 1] * y - kill_R
        phi = x * f_S + y * f_R
        dx = x * (f_S - phi) * dt
        x_adp[t + 1] = np.clip(x + dx, 1e-4, 1.0 - 1e-4)
        
    drug_on[-1] = (x_adp[-1] > threshold)
    
    return {
        "time": t_span.tolist(),
        "x_mtd": x_mtd.tolist(),
        "x_adaptive": x_adp.tolist(),
        "drug_active": drug_on.tolist(),
        "final_sensitive_mtd": float(x_mtd[-1]),
        "final_sensitive_adaptive": float(x_adp[-1])
    }
