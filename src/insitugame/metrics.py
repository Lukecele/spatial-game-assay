"""
metrics.py
==========
Diagnostic and statistical evaluation battery adapted from Luca Celebrano's
framework in meme-intelligence-onchain/analytics/backtest_validator.py.
"""

import math
from typing import Dict, Any, Tuple, List
from dataclasses import dataclass
import numpy as np


@dataclass
class ContingencyTable2x2:
    TP: int
    FP: int
    FN: int
    TN: int
    
    @property
    def total(self) -> int:
        return self.TP + self.FP + self.FN + self.TN
        
    def to_dict(self) -> Dict[str, int]:
        return {"TP": self.TP, "FP": self.FP, "FN": self.FN, "TN": self.TN}


def calculate_odds_ratio(table: ContingencyTable2x2) -> Tuple[float, Tuple[float, float]]:
    a = table.TP + 0.5
    b = table.FP + 0.5
    c = table.FN + 0.5
    d = table.TN + 0.5
    
    odds_ratio = (a * d) / (b * c)
    se_ln_or = math.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))
    ci_lower = math.exp(math.log(odds_ratio) - 1.96 * se_ln_or)
    ci_upper = math.exp(math.log(odds_ratio) + 1.96 * se_ln_or)
    return float(odds_ratio), (float(ci_lower), float(ci_upper))


def _log_hypergeometric_prob(a: int, b: int, c: int, d: int) -> float:
    n = a + b + c + d
    return (
        math.lgamma(a + b + 1) +
        math.lgamma(c + d + 1) +
        math.lgamma(a + c + 1) +
        math.lgamma(b + d + 1) -
        (math.lgamma(a + 1) +
         math.lgamma(b + 1) +
         math.lgamma(c + 1) +
         math.lgamma(d + 1) +
         math.lgamma(n + 1))
    )


def calculate_fisher_exact_p_value(table: ContingencyTable2x2) -> float:
    a, b = table.TP, table.FP
    c, d = table.FN, table.TN
    
    log_p_observed = _log_hypergeometric_prob(a, b, c, d)
    r1 = a + b
    r2 = c + d
    c1 = a + c
    
    min_x = max(0, c1 - r2)
    max_x = min(r1, c1)
    
    total_p = 0.0
    for x in range(min_x, max_x + 1):
        cur_a = x
        cur_b = r1 - x
        cur_c = c1 - x
        cur_d = r2 - cur_c
        log_p = _log_hypergeometric_prob(cur_a, cur_b, cur_c, cur_d)
        if log_p <= log_p_observed + 1e-9:
            total_p += math.exp(log_p)
            
    return float(min(1.0, total_p))


def compute_roc_auc(y_true_binary: np.ndarray, y_score: np.ndarray) -> float:
    pos_mask = (y_true_binary == 1)
    neg_mask = (y_true_binary == 0)
    n_pos = int(np.sum(pos_mask))
    n_neg = int(np.sum(neg_mask))
    
    if n_pos == 0 or n_neg == 0:
        return 0.5
        
    pos_scores = y_score[pos_mask]
    neg_scores = y_score[neg_mask]
    
    n_concordant = 0.0
    for ps in pos_scores:
        n_concordant += np.sum(ps > neg_scores) + 0.5 * np.sum(ps == neg_scores)
        
    return float(n_concordant / (n_pos * n_neg))


def evaluate_clinical_classification(y_true_continuous: np.ndarray, y_pred_continuous: np.ndarray) -> Dict[str, Any]:
    median_true = float(np.median(y_true_continuous))
    median_pred = float(np.median(y_pred_continuous))
    
    y_true_bin = (y_true_continuous >= median_true).astype(int)
    y_pred_bin = (y_pred_continuous >= median_pred).astype(int)
    
    tp = int(np.sum((y_true_bin == 1) & (y_pred_bin == 1)))
    fp = int(np.sum((y_true_bin == 0) & (y_pred_bin == 1)))
    fn = int(np.sum((y_true_bin == 1) & (y_pred_bin == 0)))
    tn = int(np.sum((y_true_bin == 0) & (y_pred_bin == 0)))
    
    table = ContingencyTable2x2(TP=tp, FP=fp, FN=fn, TN=tn)
    sens = float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0
    spec = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
    ppv = float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0
    npv = float(tn / (tn + fn)) if (tn + fn) > 0 else 0.0
    acc = float((tp + tn) / table.total) if table.total > 0 else 0.0
    
    odds_ratio, ci = calculate_odds_ratio(table)
    fisher_p = calculate_fisher_exact_p_value(table)
    auc = compute_roc_auc(y_true_bin, y_pred_continuous)
    
    return {
        "contingency_table": table.to_dict(),
        "sensitivity": sens,
        "specificity": spec,
        "ppv_precision": ppv,
        "npv": npv,
        "accuracy": acc,
        "odds_ratio": odds_ratio,
        "odds_ratio_95_ci": ci,
        "fisher_exact_p_value": fisher_p,
        "roc_auc": auc,
        "is_statistically_significant": bool(fisher_p < 0.05 and ci[0] > 1.0)
    }
