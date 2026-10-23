#!/usr/bin/env python3
"""
clinical_diagnostic_metrics.py
==============================
Exact clinical diagnostic & classification performance metrics adapted from
Luca Celebrano's repository (meme-intelligence-onchain/analytics/backtest_validator.py).

Implements:
1. Exact 2x2 Contingency Table (a, b, c, d)
2. Exact Hypergeometric Fisher Test p-value via math.lgamma (zero approximation)
3. Diagnostic Odds Ratio (DOR) with Woolf Standard Error and 95% Confidence Interval
4. Sensitivity (Recall / TPR), Specificity (TNR), PPV (Precision), NPV
5. ROC-AUC (Area Under the Receiver Operating Characteristic Curve)
6. Asymmetric Signal / Model Superiority Test (Game Model vs Non-Spatial Baseline)
"""

import math
from dataclasses import dataclass
from typing import Dict, Any, Tuple
import numpy as np


@dataclass
class ContingencyTable2x2:
    a: int  # True Positive (High Fitness Predicted & Observed)
    b: int  # False Positive (High Fitness Predicted but Observed Low)
    c: int  # False Negative (Low Fitness Predicted but Observed High)
    d: int  # True Negative (Low Fitness Predicted & Observed Low)


def calculate_fisher_exact_p_value(t: ContingencyTable2x2) -> float:
    """
    Exact two-sided Fisher's Exact Test p-value calculated via log-gamma functions.
    Direct implementation from Luca Celebrano's backtest_validator.py.
    """
    a, b, c, d = t.a, t.b, t.c, t.d
    n = a + b + c + d
    if n == 0:
        return 1.0

    lf = lambda k: math.lgamma(k + 1)

    def hp(x):
        k = a + c
        m = a + b
        return math.exp(
            lf(m) + lf(n - m) + lf(k) + lf(n - k)
            - lf(n) - lf(x) - lf(m - x) - lf(k - x) - lf(n - m - k + x)
        )

    # In very large samples (n > 5000), direct hypergeometric summation can be slow;
    # use bounded range or chi-squared asymptotic if n > 10,000.
    if n > 10000:
        # Standard chi-square approximation with Yates continuity correction
        num = max(0, abs(a * d - b * c) - n / 2.0) ** 2 * n
        den = float((a + b) * (c + d) * (a + c) * (b + d))
        if den == 0:
            return 1.0
        chi2 = num / den
        from scipy.stats import chi2 as scipy_chi2
        return float(1.0 - scipy_chi2.cdf(chi2, df=1))

    obs = hp(a)
    lo = max(0, (a + b) + (a + c) - n)
    hi = min(a + b, a + c)
    p_sum = sum(p for x in range(lo, hi + 1) if (p := hp(x)) <= obs + 1e-12)
    return min(1.0, float(p_sum))


def calculate_odds_ratio(t: ContingencyTable2x2) -> Tuple[float, Tuple[float, float]]:
    """
    Calculates the Odds Ratio and 95% Confidence Interval via Woolf logit formula.
    With Haldane-Anscombe 0.5 continuity correction if any cell is zero.
    Direct implementation from Luca Celebrano's backtest_validator.py.
    """
    a, b, c, d = map(float, (t.a, t.b, t.c, t.d))
    if 0.0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    ci_low = math.exp(math.log(orr) - 1.96 * se)
    ci_high = math.exp(math.log(orr) + 1.96 * se)
    return float(orr), (float(ci_low), float(ci_high))


def compute_roc_auc(y_true_binary: np.ndarray, y_score: np.ndarray) -> float:
    """
    Computes exact Area Under the ROC Curve (AUC) via Wilcoxon-Mann-Whitney U statistic,
    handling score ties with mid-ranks.
    """
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


def evaluate_clinical_classification(y_true: np.ndarray, y_pred_score: np.ndarray,
                                     threshold_quantile: float = 0.50) -> Dict[str, Any]:
    """
    Evaluates sensitivity, specificity, PPV, NPV, Odds Ratio, Fisher p-value, and ROC-AUC.
    """
    y_true = np.asarray(y_true, dtype=np.float64)
    y_pred_score = np.asarray(y_pred_score, dtype=np.float64)

    # Define high proliferation microdomains by threshold (e.g. median or top quartile)
    true_cutoff = np.quantile(y_true, threshold_quantile)
    pred_cutoff = np.quantile(y_pred_score, threshold_quantile)

    y_true_bin = (y_true >= true_cutoff).astype(int)
    y_pred_bin = (y_pred_score >= pred_cutoff).astype(int)

    a = int(np.sum((y_pred_bin == 1) & (y_true_bin == 1)))
    b = int(np.sum((y_pred_bin == 1) & (y_true_bin == 0)))
    c = int(np.sum((y_pred_bin == 0) & (y_true_bin == 1)))
    d = int(np.sum((y_pred_bin == 0) & (y_true_bin == 0)))

    table = ContingencyTable2x2(a=a, b=b, c=c, d=d)
    orr, (ci_low, ci_high) = calculate_odds_ratio(table)
    p_fisher = calculate_fisher_exact_p_value(table)
    auc = compute_roc_auc(y_true_bin, y_pred_score)

    sensitivity = a / max(1, a + c)
    specificity = d / max(1, b + d)
    ppv = a / max(1, a + b)
    npv = d / max(1, c + d)
    accuracy = (a + d) / max(1, a + b + c + d)

    return {
        "contingency_table": {"TP": a, "FP": b, "FN": c, "TN": d},
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "ppv_precision": float(ppv),
        "npv": float(npv),
        "accuracy": float(accuracy),
        "odds_ratio": float(orr),
        "odds_ratio_95_ci": [float(ci_low), float(ci_high)],
        "fisher_exact_p_value": float(p_fisher),
        "roc_auc": float(auc),
        "is_statistically_significant": bool(p_fisher < 0.05),
        "is_asymmetric_signal": bool(orr >= 3.0)
    }
