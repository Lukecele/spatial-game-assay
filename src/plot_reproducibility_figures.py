#!/usr/bin/env python3
"""
plot_reproducibility_figures.py
===============================
Generates high-resolution multi-panel publication figures displaying:
- Panel A: Model Ablation (Incremental R^2 and AUC across Lung, Breast, Melanoma)
- Panel B: Quality Control & Spatial Scale Sensitivity of Payoffs
- Panel C: Serial Histological Section Payoff Concordance (Sec1 vs Sec2, r = +0.811)
- Panel D: Blind Out-Of-Sample ROC Curves and Diagnostic Odds Ratios
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results"))
FIGURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "figures"))
os.makedirs(FIGURES_DIR, exist_ok=True)

# 1. Load results
with open(os.path.join(RESULTS_DIR, "reproducibility_and_sensitivity_audit.json")) as f:
    audit_data = json.load(f)

with open(os.path.join(RESULTS_DIR, "serial_replicate_blind_oos_validation.json")) as f:
    oos_data = json.load(f)

# Set high-quality style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
plt.subplots_adjust(wspace=0.28, hspace=0.32)

# =============================================================================
# PANEL A: MODEL ABLATION ON IDENTICAL TARGET (w_purified) ACROSS COHORTS
# =============================================================================
ax_a = axes[0, 0]
cancers = ["lung", "breast", "melanoma"]
labels = ["Lung SCC\n(N=3,858)", "Breast IDC\n(N=4,992)", "Melanoma SKCM\n(N=3,458)"]
x = np.arange(len(cancers))
width = 0.15

r2_m2 = [27.84, 23.49, 32.48]
r2_m2b = [28.26, 28.38, 33.44]
r2_m2c = [30.77, 26.58, 35.52]
r2_m2d = [31.21, 33.32, 36.58]
r2_m3 = [30.86, 32.29, 36.78]

rects2 = ax_a.bar(x - 2*width, r2_m2, width, label="$M_2$: Mass Action ($f_i$)", color="#FFE082", edgecolor="#FF8F00")
rects2b = ax_a.bar(x - width, r2_m2b, width, label="$M_{2b}$: Co-localization ($f_i f_j$)", color="#B3E5FC", edgecolor="#0288D1")
rects2c = ax_a.bar(x, r2_m2c, width, label="$M_{2c}$: Regional ($f_i + \\bar{f}_i$)", color="#FFAB91", edgecolor="#D84315")
rects2d = ax_a.bar(x + width, r2_m2d, width, label="$M_{2d}$: Quadratics ($f_i f_j + \\bar{f}_i \\bar{f}_j$)", color="#CE93D8", edgecolor="#6A1B9A")
rects3 = ax_a.bar(x + 2*width, r2_m3, width, label="$M_3$: Spatial Game ($f_i \\bar{f}_j$)", color="#A5D6A7", edgecolor="#2E7D32")

ax_a.set_ylabel("Explained Variance on $w_{purified}$ ($R^2$, %)", fontsize=11, fontweight="bold")
ax_a.set_title("A. Model Ablation on Decoupled Target ($w_{purified}$)\n(Breast: p < 1e-16; Melanoma: p = 6e-14; Lung: p = 0.27 N.S.)", fontsize=10, fontweight="bold", loc="left")
ax_a.set_xticks(x)
ax_a.set_xticklabels(labels, fontsize=10)
ax_a.legend(frameon=True, fontsize=7.5, loc="upper left", ncol=2)
ax_a.set_ylim(0, 45)

for rect in rects3:
    height = rect.get_height()
    ax_a.annotate(f"{height:.1f}%",
                  xy=(rect.get_x() + rect.get_width() / 2, height),
                  xytext=(0, 3), textcoords="offset points",
                  ha='center', va='bottom', fontsize=8, fontweight="bold", color="#1B5E20")

# =============================================================================
# PANEL B: QC CUTOFF ROBUSTNESS (BREAST CANCER ABSORBED PAYOFFS)
# =============================================================================
ax_b = axes[0, 1]
cutoffs = ["Unfiltered\n(0 UMI)", "Standard QC\n(>=500 UMI)", "High Stringency\n(>=1000 UMI)"]
spots_cnt = [audit_data["breast"]["qc_sensitivity"]["cut_0"]["n_spots"],
             audit_data["breast"]["qc_sensitivity"]["cut_500"]["n_spots"],
             audit_data["breast"]["qc_sensitivity"]["cut_1000"]["n_spots"]]

a_ss = [audit_data["breast"]["qc_sensitivity"][k]["payoff_2x2"]["A_SS"] for k in ["cut_0", "cut_500", "cut_1000"]]
a_sr = [audit_data["breast"]["qc_sensitivity"][k]["payoff_2x2"]["A_SR"] for k in ["cut_0", "cut_500", "cut_1000"]]
a_rs = [audit_data["breast"]["qc_sensitivity"][k]["payoff_2x2"]["A_RS"] for k in ["cut_0", "cut_500", "cut_1000"]]
a_rr = [audit_data["breast"]["qc_sensitivity"][k]["payoff_2x2"]["A_RR"] for k in ["cut_0", "cut_500", "cut_1000"]]

x_b = np.arange(len(cutoffs))
ax_b.plot(x_b, a_sr, 'o-', color="#2E7D32", lw=2, ms=7, label=r"$\tilde{A}_{SR}$ (Sens vs Res)")
ax_b.plot(x_b, a_rs, 's-', color="#1565C0", lw=2, ms=7, label=r"$\tilde{A}_{RS}$ (Res vs Sens)")
ax_b.plot(x_b, a_ss, '^-', color="#C62828", lw=2, ms=7, label=r"$\tilde{A}_{SS}$ (Sens vs Sens)")
ax_b.plot(x_b, a_rr, 'd-', color="#6A1B9A", lw=2, ms=7, label=r"$\tilde{A}_{RR}$ (Res vs Res)")
ax_b.axhline(0, color="grey", linestyle="--", alpha=0.7)

ax_b.set_ylabel(r"Absorbed Payoff Magnitude ($\tilde{A}_{ij}$)", fontsize=11, fontweight="bold")
ax_b.set_title("B. Quality Control Stringency Sensitivity (Breast IDC)\n(Stable Hawk-Dove Absorbed Sign Pattern Across Thresholds)", fontsize=10, fontweight="bold", loc="left")
ax_b.set_xticks(x_b)
ax_b.set_xticklabels([f"{c}\nN={s}" for c, s in zip(cutoffs, spots_cnt)], fontsize=9)
ax_b.legend(frameon=True, fontsize=8.5, loc="center right")
ax_b.set_ylim(-2.5, 3.5)

# =============================================================================
# PANEL C: SERIAL SECTION PAYOFF CONCORDANCE (SECTION 1 vs SECTION 2)
# =============================================================================
ax_c = axes[1, 0]
A1 = np.array(oos_data["section_1_in_sample"]["payoff_matrix_3x3"]).flatten()
A2 = np.array(oos_data["section_2_in_sample"]["payoff_matrix_3x3"]).flatten()

param_labels = [
    "$A_{SS}$", "$A_{SR}$", "$A_{SC}$",
    "$A_{RS}$", "$A_{RR}$", "$A_{RC}$",
    "$A_{CS}$", "$A_{CR}$", "$A_{CC}$"
]

colors_c = ["#C62828", "#2E7D32", "#1565C0", "#E65100", "#6A1B9A", "#00838F", "#4E342E", "#37474F", "#000000"]

ax_c.scatter(A1, A2, s=80, c=colors_c, zorder=5, edgecolors="black", linewidth=1.2)
for i, txt in enumerate(param_labels):
    ax_c.annotate(f" {txt}", (A1[i], A2[i]), fontsize=10, fontweight="bold")

# Linear fit
slope, intercept = np.polyfit(A1, A2, 1)
x_line = np.linspace(-6, 7, 50)
ax_c.plot(x_line, slope * x_line + intercept, color="#D32F2F", linestyle="--", lw=1.8,
          label=f"Fit: y = {slope:.2f}x + {intercept:.2f}")

ax_c.axhline(0, color="gray", lw=0.8, linestyle=":")
ax_c.axvline(0, color="gray", lw=0.8, linestyle=":")

r_pearson = oos_data["parameter_concordance"]["payoff_matrix_pearson_r"]
p_pearson = oos_data["parameter_concordance"]["payoff_matrix_pearson_p"]
ax_c.set_xlabel("Section 1 Payoff Parameters ($A^{(1)}$)", fontsize=11, fontweight="bold")
ax_c.set_ylabel("Section 2 Payoff Parameters ($A^{(2)}$)", fontsize=11, fontweight="bold")
ax_c.set_title(f"C. Replicate Concordance Across Consecutive Slices (~10 um)\n(Pearson r = +{r_pearson:.3f}, p = {p_pearson:.4f})", fontsize=11, fontweight="bold", loc="left")
ax_c.legend(frameon=True, fontsize=9, loc="upper left")

# =============================================================================
# PANEL D: BLIND OUT-OF-SAMPLE (OOS) GENERALIZATION
# =============================================================================
ax_d = axes[1, 1]

models_d = ["Non-Spatial Intrinsic", "In Situ Game Assay"]
oos_r2_vals = [
    oos_data["blind_oos_evaluation"]["oos_r2_intrinsic"] * 100,
    oos_data["blind_oos_evaluation"]["oos_r2_game"] * 100
]
dors = [
    oos_data["blind_oos_evaluation"]["diagnostic_intrinsic_oos"]["odds_ratio"],
    oos_data["blind_oos_evaluation"]["diagnostic_game_oos"]["odds_ratio"]
]
ci_low = [
    oos_data["blind_oos_evaluation"]["diagnostic_intrinsic_oos"]["odds_ratio_95_ci"][0],
    oos_data["blind_oos_evaluation"]["diagnostic_game_oos"]["odds_ratio_95_ci"][0]
]
ci_high = [
    oos_data["blind_oos_evaluation"]["diagnostic_intrinsic_oos"]["odds_ratio_95_ci"][1],
    oos_data["blind_oos_evaluation"]["diagnostic_game_oos"]["odds_ratio_95_ci"][1]
]

x_d = np.arange(len(models_d))
y_err = [
    [dors[0] - ci_low[0], dors[1] - ci_low[1]],
    [ci_high[0] - dors[0], ci_high[1] - dors[1]]
]

bar_d = ax_d.bar(x_d, dors, yerr=y_err, capsize=6, color=["#FFCC80", "#81C784"], edgecolor=["#E65100", "#2E7D32"], width=0.45)
ax_d.set_ylabel("Diagnostic Odds Ratio (DOR, 95% CI)", fontsize=11, fontweight="bold")
ax_d.set_title("D. Blind Out-of-Sample Performance on Held-Out Section 2\n(Trained on Sec1, Tested Blind on Sec2; Fisher p = 4.8e-12)", fontsize=11, fontweight="bold", loc="left")
ax_d.set_xticks(x_d)
ax_d.set_xticklabels(models_d, fontsize=10, fontweight="bold")
ax_d.set_ylim(0, 3.5)

for i, rect in enumerate(bar_d):
    h = rect.get_height()
    ax_d.text(rect.get_x() + rect.get_width() / 2, h + 0.35,
              f"DOR = {dors[i]:.2f}\nOOS R² = {oos_r2_vals[i]:.2f}%\nAUC = {oos_data['blind_oos_evaluation']['diagnostic_' + ('intrinsic' if i==0 else 'game') + '_oos']['roc_auc']:.3f}",
              ha='center', va='bottom', fontsize=9, fontweight="bold")

plt.tight_layout()
out_png = os.path.join(FIGURES_DIR, "figure_reproducibility_battery.png")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"[✓] Multi-panel validation figure saved to: {out_png}")
