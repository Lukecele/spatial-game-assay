"""
insitugame: In Situ Spatial Evolutionary Game Assay
===================================================
Author: Luca Celebrano (School of Medicine and Surgery, Università degli Studi di Napoli Federico II)
License: GNU General Public License v3.0 (GPL-3.0)

A computational oncology assay to infer subclonal payoff matrices,
solve Nash equilibria, decouple microenvironmental confounders,
and simulate adaptive therapy regimens directly from intact spatial transcriptomic biopsies.
"""

__version__ = "1.0.0"
__author__ = "Luca Celebrano"

from .assay import InSituGameAssay
from .decoupling import EnvironmentalDecoupler
from .dynamics import solve_subgame_nash, simulate_adaptive_therapy
from .metrics import (
    ContingencyTable2x2,
    calculate_odds_ratio,
    calculate_fisher_exact_p_value,
    compute_roc_auc,
    evaluate_clinical_classification
)

__all__ = [
    "InSituGameAssay",
    "EnvironmentalDecoupler",
    "solve_subgame_nash",
    "simulate_adaptive_therapy",
    "ContingencyTable2x2",
    "calculate_odds_ratio",
    "calculate_fisher_exact_p_value",
    "compute_roc_auc",
    "evaluate_clinical_classification"
]
