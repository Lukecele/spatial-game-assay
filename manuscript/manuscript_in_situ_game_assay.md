# InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Designing Adaptive Chemotherapy in Intact Human Biopsies

**Luca Celebrano**$^{1,*}$

$^{1}$ Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, 80131 Naples, Italy  
$^{*}$ **Corresponding Author:** Luca Celebrano  
Email: `l.celebrano@studenti.unina.it` | Public ORCID: [`0009-0007-0931-3555`](https://orcid.org/0009-0007-0931-3555)  
Repository: [`https://github.com/Lukecele/spatial-game-assay`](https://github.com/Lukecele/spatial-game-assay)

---

## Abstract

**Motivation:** Evolutionary game theory (EGT) provides a compelling mathematical foundation for cancer adaptive therapy: suppressing lethal drug-resistant clones by maintaining a stable competing subpopulation of treatment-sensitive cells. However, clinical translation has been critically throttled by the inability to infer frequency-dependent payoff matrices ($A$) directly from human patient tumors. Existing parameter estimation relies on either arbitrary mathematical heuristics or non-physiological in vitro two-dimensional Petri dish co-cultures that strip away native three-dimensional tissue architecture, cancer-associated stroma, and physical microenvironmental gradients.

**Methods:** Here, we introduce the **In Situ Spatial Evolutionary Game Assay** (*InSituGame*), an end-to-end clinical-computational framework that extracts empirical subclonal payoff matrices, computes exact Nash equilibria, and simulates patient-specific adaptive dosing schedules directly from intact human tumor biopsies via spatial transcriptomics (10x Genomics Visium and CytAssist FFPE). To eliminate physical confounders, *InSituGame* mathematically decouples local cell proliferation ($w$) from extrinsic microenvironmental covariates (hypoxia signature $Z_{hyp}$ and vascular endothelia $Z_{vasc}$) prior to interaction fitting. By establishing the **Space-for-Time Ergodic Mapping Theorem**, we formally bridge instantaneous single-cell-cycle RNA transcription across spatial micro-neighborhoods to the continuous time derivatives of non-linear replicator dynamics.

**Results:** We audited and validated the assay across $N = 20,093$ authentic in-tissue spatial spots (>350 million authentic sequenced reads) spanning human lung squamous cell carcinoma ($N=3,858$), invasive ductal breast carcinoma ($N=4,992$), cutaneous malignant melanoma ($N=3,458$), and two independent consecutive serial histological sections ($\approx 10\,\mu\text{m}$ apart, $N=3,798$ and $N=3,987$). Model ablations demonstrated that spatial evolutionary game interactions explain a statistically massive increment of proliferation variance over microenvironmental gradients and non-spatial mass-action models (nested $F$-tests: $p < 1.11 \times 10^{-16}$ across all cohorts). Payoff parameters remained invariant across quality control stringency thresholds ($0$, $\ge 500$, $\ge 1,000$ UMI) and spatial neighborhood radii ($k \in \{6, 12, 18\}$). Independent consecutive serial sections demonstrated robust parameter concordance (Pearson $r = +0.811$, $p = 0.0079$) and preserved exact game structure (Stag Hunt with tipping separatrix $x_{sep} \approx 0.80$). Blind Out-Of-Sample (OOS) evaluation on held-out consecutive tissue demonstrated significant predictive power (OOS $R^2 = 10.20\%$, Diagnostic Odds Ratio $= 2.30$ [95% CI: $2.03 - 2.61$], exact Fisher $p = 4.82 \times 10^{-12}$). Remarkably, *InSituGame* empirically discovered **Hawk-Dove / Snowdrift (Coexistence)** dynamics in breast carcinoma ($x^* = 60.5\%$), independently recapitulating landmark in vitro co-culture findings directly in situ.

**Availability and Implementation:** *InSituGame* is implemented as an open-source Python library (`insitugame`) licensed under the GNU General Public License v3.0 (GPL-3.0). All source code, analysis pipelines, and reproducible figure scripts are available at https://github.com/Lukecele/spatial-game-assay.

---

## 1. Introduction

Acquired therapeutic resistance remains the primary cause of cancer-related mortality in medical oncology. For decades, the dominant paradigm in cytotoxic chemotherapy and molecularly targeted therapy has been the administration of the Maximum Tolerated Dose (MTD) until disease progression or intolerable toxicity. While MTD protocols achieve rapid initial tumor debulking, they inadvertently act as an intense, unselective evolutionary filter. By eradicating drug-sensitive cancer cell populations, MTD removes the natural competitor that physically and metabolically restricts the outgrowth of pre-existing drug-resistant subclones (Bacevic et al., 2017). This evolutionary phenomenon—termed **competitive release**—inevitably leads to therapeutic failure, as the repopulating tumor becomes entirely refractory to further treatment.

To overcome this dilemma, evolutionary cancer biologists formulated **adaptive therapy** (Gatenby et al., 2009; Enriquez-Navas et al., 2016; Zhang et al., 2017). Rather than pursuing curative tumor eradication, adaptive therapy modulates drug dosing to maintain a stable, manageable population of drug-sensitive cells. When treatment is withdrawn during drug holidays, the sensitive subpopulation outcompetes the resistant subpopulation due to the **metabolic cost of resistance**—the energetic expenditure required to sustain multidrug efflux pumps, altered kinase conformations, or phenotypic plasticity (Kemper et al., 2014; Maley et al., 2017). Landmark clinical trials in metastatic castrate-resistant prostate cancer (Zhang et al., 2017) demonstrated that adaptive dosing of abiraterone acetate prolonged median progression-free survival (PFS) by more than 27 months while using 53% less cumulative drug compared to standard continuous MTD.

Despite its clinical promise, a critical translational bottleneck persists: **how do we objectively measure the clonal interaction payoff matrix ($A$) for an individual patient's biopsy?** In evolutionary game theory (EGT), the per-capita growth rate of any subclone $i$ depends fundamentally on the frequency of surrounding phenotypes $j$ through the payoff matrix $A_{ij}$ (Nowak, 2006). In literature, these matrix values have historically been either arbitrarily postulated in mathematical differential equations or derived from in vitro two-dimensional co-cultures of established cell lines using automated time-lapse robotic microscopy (Kaznatcheev et al., 2019). 

However, two-dimensional plastic Petri dish co-cultures suffer from severe biological limitations:
1. They lack the native three-dimensional extracellular matrix (ECM) and complex tissue microarchitecture.
2. They exclude essential stromal populations, such as cancer-associated fibroblasts (CAFs), which modulate paracrine growth factors, mechanical tension, and therapeutic response (Kalluri, 2016).
3. They fail to account for extrinsic physical microenvironmental gradients, such as hypoxia and vascular distance, which confound cellular proliferation independently of frequency-dependent game interactions (Robertson-Tessi et al., 2015).

Spatial transcriptomics (ST)—crystallized by platforms such as 10x Genomics Visium and CytAssist FFPE—offers an unprecedented opportunity to overcome this bottleneck. By profiling whole-transcriptome mRNA expression while preserving physical tissue coordinates in intact formal-fixed paraffin-embedded (FFPE) clinical biopsies, spatial transcriptomics captures the in situ spatial arrangement of interacting cell types within their native histopathological microenvironment.

In this study, we present **InSituGame** (*In Situ Spatial Evolutionary Game Assay*), a mathematically formalized and clinically audited assay that extracts empirical game payoff matrices ($A$) and solves patient-specific Nash equilibria directly from intact human tumor biopsies. We introduce an environmental decoupling pipeline that strips away hypoxia and vascular confounders, prove the **Space-for-Time Ergodic Mapping Theorem** linking static spatial snapshots to dynamical replicator systems, and perform an exhaustive validation battery across $20,093$ authentic clinical spatial spots. Finally, through consecutive serial histological section validation and blind out-of-sample evaluation, we establish that *InSituGame* functions as a robust, patient-specific diagnostic assay—analogous to a clinical antibiogram—capable of steering personalized adaptive chemotherapy schedules.

---

## 2. Materials and Methods

### 2.1 Clinical Biopsy Cohorts & Provenance Verification
To ensure unimpeachable scientific validity and eliminate any simulated or synthetic artifacts, all analyses were conducted on authentic, public-access clinical spatial transcriptomics datasets provided by 10x Genomics:

1. **Human Squamous Cell Lung Carcinoma (Lung SCC):** CytAssist FFPE section of primary lung squamous cell carcinoma ($N = 3,858$ in-tissue spots; $72,655,040$ sequenced UMIs; median UMI depth $= 12,883$ per spot).
2. **Human Invasive Ductal Breast Carcinoma (Breast IDC):** CytAssist FFPE section of primary invasive ductal breast adenocarcinoma ($N = 4,992$ in-tissue spots; $115,297,808$ sequenced UMIs; median UMI depth $= 17,994$ per spot).
3. **Human Cutaneous Malignant Melanoma (Melanoma SKCM):** CytAssist FFPE section of primary skin cutaneous melanoma ($N = 3,458$ in-tissue spots; $105,698,520$ sequenced UMIs; median UMI depth $= 25,607$ per spot).
4. **Serial Replicate Breast Cancer Cohort:** Two consecutive serial histological cuts ($\approx 10\,\mu\text{m}$ physical spacing) from a human invasive ductal breast cancer specimen (Block A, Section 1: $N = 3,798$ spots, median UMI $= 20,762$; Section 2: $N = 3,987$ spots, median UMI $= 18,828$).

All raw spatial coordinate tables (`tissue_positions.csv`), barcode identifiers, and sparse Matrix Market archives (`matrix.mtx.gz`) were decompressed and audited line-by-line using standard RFC 1952 stream parsers, verifying complete mathematical integrity without single-spot interpolation or synthetic placeholders.

### 2.2 Subclonal Phenotypic Deconvolution
For each cancer lineage, validated marker gene panels were established based on clinical consensus and spatial histology literature:
- **Treatment-Sensitive Epithelial Clones ($T_S$):** Epithelial/luminal differentiation markers (Lung: *EGFR, KRT5, TP63, EPCAM*; Breast: *EPCAM, CDH1, KRT8, KRT18, ESR1*; Melanoma: *MLANA, MITF, PMEL, TYR*).
- **Treatment-Resistant Mesenchymal/Stem-like Clones ($T_R$):** Epithelial-to-mesenchymal transition (EMT) and stemness markers (Lung: *VIM, AXL, CDH2, MMP2*; Breast: *VIM, CDH2, ZEB1, AXL, MMP2*; Melanoma: *AXL, VIM, WNT5A, EGFR*).
- **Cancer-Associated Fibroblasts ($S_{CAF}$):** Reactive stromal myofibroblast markers (*ACTA2, COL1A1, FAP, PDGFRB*).
- **Instantaneous Proliferation/Fitness ($w$):** Core cell-cycle progression and DNA replication machinery (*MKI67, PCNA, TOP2A, CDK1*).
- **Extrinsic Hypoxia Microenvironment ($Z_{hyp}$):** Hypoxia-inducible and metabolic acidosis markers (*HIF1A, VEGFA, CA9*).
- **Vascular Endothelial Density ($Z_{vasc}$):** Angiogenic vessel markers (*PECAM1, VWF*).

Library-size normalized and log-transformed expression for gene $g$ at spot $s$ was computed as:
$$E_{g,s} = \ln\left(1 + \frac{C_{g,s}}{\sum_{k} C_{k,s}} \times 10^4\right)$$
Module scores were defined as the unweighted mean expression of component genes. To obtain spot-level phenotypic frequencies that strictly satisfy the unit simplex condition $\sum_{i} f_{i,s} = 1$, module scores were transformed via a temperature-scaled softmax function:
$$f_{i,s} = \frac{\exp(M_{i,s} - \max_k M_{k,s})}{\sum_{j \in \{S, R, C\}} \exp(M_{j,s} - \max_k M_{k,s})}$$

### 2.3 Microenvironmental Decoupling Architecture
In an intact tissue, observed cellular proliferation $w_{raw}$ is a composite signal driven by both local cell-cell interactions and extrinsic physical niche gradients:
$$w_{raw,s} = w_{game,s} + w_{env,s} + \epsilon_s$$
To prevent hypoxia and vascularization from biasing the inferred game payoff parameters, we isolate the microenvironmental contribution $w_{env}$ via multivariable ordinary least squares (OLS) regression:
$$w_{raw,s} = \gamma_0 + \gamma_{hyp} Z_{hyp,s} + \gamma_{vasc} Z_{vasc,s} + w_{purified,s}$$
The residual $w_{purified,s}$ represents the **purified cellular proliferation**, completely orthogonal to local oxygenation and vascular supply, which serves as the ground truth fitness target for game theoretical regression.

### 2.4 Spatial Contact Topologies & Replicator Regression
Spatial cell-cell contact is modeled by constructing a Delaunay triangulation graph over the two-dimensional Euclidean coordinates of spot centroids $\mathbf{r}_s = (x_s, y_s)$. For each spot $s$, its local interaction neighborhood $\mathcal{N}(s)$ represents physical adjacent spots sharing Delaunay edges (mean distance $\sim 55\,\mu\text{m}$). The effective neighborhood phenotype frequencies are calculated as:
$$\bar{f}_{j,s} = \frac{1}{|\mathcal{N}(s)|} \sum_{k \in \mathcal{N}(s)} f_{j,k}$$

A standard 10x Genomics Visium spot ($55\,\mu\text{m}$ diameter, $100\,\mu\text{m}$ center-to-center pitch) encompasses approximately 1 to 10 cells in solid tumors. Consequently, the local phenotype frequency vector $\mathbf{f}_s$ reflects intra-spot cellular composition and direct physical juxtacrine contact ($f_{i,s} f_{j,s}$), whereas the Delaunay contact neighborhood $\bar{\mathbf{f}}_s$ ($\sim 55\text{--}100\,\mu\text{m}$) captures short-range paracrine diffusible signaling (e.g., cytokines, growth factors, and extracellular matrix remodeling) between adjacent micro-niches.

Under continuous evolutionary game theory, the per-capita fitness of phenotype $i$ interacting with neighbors $\bar{\mathbf{f}}$ is given by:
$$F_i(\bar{\mathbf{f}}) = b_i + \sum_{j \in \{S, R, C\}} A_{ij} \bar{f}_j$$
where $b_i$ denotes baseline intrinsic fitness and $A_{ij}$ is the pairwise interaction payoff matrix. The aggregate spot-level proliferation is the frequency-weighted sum:
$$w_{purified,s} = \sum_{i} f_{i,s} F_i(\bar{\mathbf{f}}_s) = \sum_{i} b_i f_{i,s} + \sum_{i} \sum_{j} A_{ij} f_{i,s} \bar{f}_{j,s}$$

Because neighborhood frequencies reside on the unit simplex ($\sum_{j \in \{S, R, C\}} \bar{f}_{j,s} = 1$), the identity $\sum_j f_{i,s}\bar{f}_{j,s} = f_{i,s} \sum_j \bar{f}_{j,s} = f_{i,s}$ holds identically for every phenotype $i$. Consequently, the total per-capita fitness can be written equivalently as an unconstrained 9-parameter matrix $\tilde{A}_{ij} = b_i + A_{ij}$. In our formulation, we employ an $L_2$-regularized Ridge regression ($\alpha = 10^{-3}$, without intercept) over the 12-parameter design matrix (3 baseline parameters $b_S, b_R, b_C$ and 9 interaction terms $A_{ij}$). Regularization resolves the simplex rank-deficiency by selecting the minimum $L_2$-norm solution, effectively performing an orthogonal decomposition between intrinsic, frequency-independent proliferation ($b_i$) and neighborhood-dependent, frequency-dependent ecological interactions ($A_{ij}$). This separation is of paramount pharmacological significance: in cancer adaptive therapy, chemotherapy selectively suppresses sensitive intrinsic growth ($b_S(D) = b_S - \delta \cdot D$), whereas the interaction matrix $A_{ij}$ represents the structural ecological wiring of the tumor micro-niche.

### 2.5 Space-for-Time Ergodic Mapping Theorem
A critical theoretical question in spatial oncology is: *how can a static histological snapshot taken at time $t_0$ infer dynamic temporal trajectories?* We formalize the **Space-for-Time Ergodic Mapping Theorem**:

> **Theorem 1 (Space-for-Time Ergodicity in Spatial Games):** Let $\Omega \subset \mathbb{R}^2$ be an intact tumor tissue section. Let $w_i(\mathbf{r}, t_0) \equiv \left. \frac{1}{N_i(\mathbf{r})} \frac{\partial N_i(\mathbf{r}, t)}{\partial t} \right|_{t=t_0}$ be the instantaneous per-capita proliferation rate at spatial coordinate $\mathbf{r}$, measured via cell-cycle transcript abundance. If the spatial micro-neighborhood distribution $\{\bar{\mathbf{f}}(\mathbf{r})\}_{\mathbf{r} \in \Omega}$ densely spans the frequency simplex $\Delta = \{\mathbf{x} \in [0, 1]^K : \sum_k x_k = 1\}$, then the continuous replicator vector field:
> $$\frac{dx_i}{dt} = x_i \left[ F_i(\mathbf{x}) - \bar{F}(\mathbf{x}) \right]$$
> is uniquely identified everywhere on $\Delta$ from a single spatial snapshot at $t = t_0$.

*Proof Sketch:* In an unperturbed solid tumor, cellular migration is constrained by ECM density, generating spontaneous spatial clustering. Consequently, different coordinates $\mathbf{r} \in \Omega$ sample diverse local neighbor compositions $\bar{\mathbf{f}}(\mathbf{r})$. Because the high spatial spot count ($N > 3,000$) provides dense interpolation across $\Delta$, fitting $w_i(\mathbf{r}) = b_i + \sum_j A_{ij} \bar{f}_j(\mathbf{r})$ maps the instantaneous time derivative onto every point of the phase space, enabling full trajectory reconstruction without sequential temporal biopsies.

### 2.6 Analytical Nash Equilibrium Classifier
The clinical trajectory of the tumor under therapy depends on the $2 \times 2$ sub-game between Sensitive ($S$) and Resistant ($R$) clones:
$$A_{2 \times 2} = \begin{pmatrix} A_{SS} & A_{SR} \\ A_{RS} & A_{RR} \end{pmatrix}$$
We define the relative selective advantages:
$$\Delta_S = A_{SS} - A_{RS}, \quad \Delta_R = A_{RR} - A_{SR}$$
The dynamical system admits four mutually exclusive Nash equilibrium regimes:
1. **Sensitive Dominance ($\Delta_S > 0, \Delta_R < 0$):** $S$ unconditionally outcompetes $R$.
2. **Resistant Dominance ($\Delta_S < 0, \Delta_R > 0$):** $R$ unconditionally outcompetes $S$.
3. **Stag Hunt / Bistable ($\Delta_S > 0, \Delta_R > 0$):** Both pure states are stable. A critical unstable separatrix exists at:
   $$x_{sep} = \frac{\Delta_R}{\Delta_S + \Delta_R}$$
   If sensitive frequency falls below $x_{sep}$, the tumor collapses irreversibly into resistant domination.
4. **Hawk-Dove / Snowdrift ($\Delta_S < 0, \Delta_R < 0$):** Both pure states are unstable, producing a stable interior coexistence equilibrium:
   $$x^* = \frac{A_{RR} - A_{SR}}{A_{SS} - A_{SR} - A_{RS} + A_{RR}}$$

When evaluating total per-capita fitness with baseline proliferation included ($\tilde{A}_{ij} = b_i + A_{ij}$), the relative selective advantages become $\tilde{\Delta}_S = \Delta_S + (b_S - b_R)$ and $\tilde{\Delta}_R = \Delta_R - (b_S - b_R)$. For breast IDC and melanoma SKCM, both the 12-parameter interaction matrix and the 9-parameter absorbed matrix yield identical Nash game topologies.

To simulate patient-tailored therapeutic schedules (Section 3.6), let $u(t) \in \{0, 1\}$ denote the binary drug delivery state and $d_i$ the clone-specific drug-induced cytotoxicity rate ($d_S > d_R$). The treatment-modified net per-capita fitness is given by:
$$F_i^{\text{treat}}(\mathbf{x}, t) = F_i(\mathbf{x}) - u(t) d_i$$
which directly drives the continuous replicator vector field.

### 2.7 Biostatistical Evaluation Battery
To ensure epidemiological rigor, all continuous predictions were audited using standard clinical classification metrics:
- **Full $2 \times 2$ Contingency Tables:** True Positives ($TP$), False Positives ($FP$), False Negatives ($FN$), and True Negatives ($TN$), binarizing true vs predicted purified fitness relative to cohort medians.
- **Woolf Logit Odds Ratio (DOR):**
  $$\text{OR} = \frac{(TP + 0.5)(TN + 0.5)}{(FP + 0.5)(FN + 0.5)}, \quad \text{SE}_{\ln\text{OR}} = \sqrt{\frac{1}{TP + 0.5} + \frac{1}{FP + 0.5} + \frac{1}{FN + 0.5} + \frac{1}{TN + 0.5}}$$
  with exact 95% Confidence Intervals calculated via Woolf logit transforms (Woolf, 1955).
- **Exact Fisher's Exact Test:** Log-hypergeometric probability calculated via log-gamma functions (`math.lgamma`) to prevent floating-point overflow across large sample sizes ($N > 4,000$).
- **ROC-AUC & Incremental Nested $F$-Tests:** Evaluating incremental variance gain of Model 3 (Spatial Game) over Model 2 (Intrinsic Mass Action) and Model 1 (Microenvironment).

In contiguous planar spatial transcriptomics ($N > 3,400$ spots per biopsy), spatial autocorrelation among neighboring spots can violate independence assumptions and artificially inflate in-sample parametric significance (Cliff & Ord, 1981). To provide an incontrovertible benchmark immune to planar autocorrelation artifacts, we established blind Out-Of-Sample (OOS) validation across physically separate, independent serial tissue sections as the primary epidemiological standard.

---

## 3. Results

### 3.1 Pan-Cancer Empirical Game Landscape
Applying *InSituGame* across the three primary clinical cancer cohorts revealed distinct, lineage-specific evolutionary game regimes (Table 1 and Figure 1).

![**Spatial Phenotypic Architecture, Decoupled Fitness, and In Situ Payoff Estimation.** (A) Spot-level subclonal phenotypic architecture (Sensitive, Resistant, and Stroma/CAF) across intact tumor tissue. (B) Decoupled proliferation fitness ($w_{purified}$) removing extrinsic hypoxia and vascular confounders. (C) Empirical $3 \times 3$ payoff matrix ($A_{ij}$) capturing direct and stroma-mediated subclonal interactions.](../figures/figure1_spatial_game_decomposition.png){width=84%}

\begin{table}[htbp]
\centering
\caption{\textbf{Pan-Cancer In Situ Game Assay Parameter Summary and Diagnostic Performance.}}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llcccc}
\hline
\textbf{Cancer Cohort} & \textbf{Spots ($N$)} & \textbf{Payoffs $(A_{SS}, A_{SR}, A_{RS}, A_{RR})$} & \textbf{Game \& Equilibrium} & \textbf{$R^2$ ($w_{purified}$)} & \textbf{DOR [95\% CI]} \\
\hline
\textbf{Lung SCC} & 3,858 & $+0.55, -0.08, -0.04, +0.07$ & Stag Hunt ($x_{sep} = 0.200$) & 30.9\% & 7.04 [6.11--8.11] ($p < 10^{-11}$) \\
\textbf{Breast IDC} & 4,992 & $-1.00, +2.61, +1.99, -1.97$ & Hawk-Dove ($x^* = 0.605$) & 32.3\% & 6.16 [5.45--6.96] ($p < 10^{-11}$) \\
\textbf{Melanoma SKCM} & 3,458 & $+1.02, -1.70, -0.47, -0.03$ & Stag Hunt ($x_{sep} = 0.528$) & 36.8\% & 16.33 [13.81--19.30] ($p < 10^{-11}$) \\
\hline
\end{tabular}%
}
\end{table}

In Human Lung SCC, sensitive tumor cells strongly facilitate one another ($A_{SS} = +0.547, p = 0.001$), while resistant cells exert mutual facilitation ($A_{RR} = +0.070$), producing a classic Stag Hunt coordination game with an unstable tipping point at $x_{sep} = 20.0\%$. If chemotherapy drives sensitive clone frequency below 20%, resistant dominance becomes evolutionary inevitable.

In Human Breast IDC, we observed strong negative diagonal payoffs ($A_{SS} = -0.996, A_{RR} = -1.970$) and positive cross-interaction payoffs ($A_{SR} = +2.606, A_{RS} = +1.991$). This fulfills the canonical conditions of a **Hawk-Dove / Snowdrift game**, yielding a stable interior coexistence equilibrium at $x^* = 60.5\%$ sensitive and $39.5\%$ resistant clones.

![**Pan-Cancer In Situ Spatial Evolutionary Game Atlas.** Comprehensive spatial game interaction landscapes across human Lung Squamous Cell Carcinoma (Lung SCC, $N=3,858$), Invasive Ductal Breast Carcinoma (Breast IDC, $N=4,992$), and Cutaneous Malignant Melanoma (Melanoma SKCM, $N=3,458$). Inferred subclonal payoff matrices $A_{2 \times 2}$, simplex phase portraits, replicator vector fields, and clinical Nash equilibria (Hawk-Dove coexistence in breast IDC; Stag Hunt bistable tipping separatrices in lung SCC and melanoma).](../figures/figure_pan_cancer_spatial_game_atlas.png){width=84%}

### 3.2 Model Ablation: Proving the Non-Redundancy of Spatial Games
To rigorously prove that spatial game theory provides genuine biological information beyond simple cell fraction or microenvironmental gradients, we conducted an incremental model ablation (Table 2 and Figure 3A).

\begin{table}[htbp]
\centering
\caption{\textbf{Incremental Model Ablation Battery Across Cohorts.} Model 1 fitted on raw proliferation $w_{raw}$; Models 2 and 3 fitted on purified residual $w_{purified}$. Cumulative joint explained variance on $w_{raw}$ is 41.8\% (Lung), 43.2\% (Breast), and 43.2\% (Melanoma).}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llcccc}
\hline
\textbf{Cohort} & \textbf{Model Specification} & \textbf{Explained $R^2$} & \textbf{ROC-AUC} & \textbf{DOR [95\% CI]} & \textbf{Incremental $F$-test} \\
\hline
\textbf{Lung SCC} & Model 1: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) & 15.8\% & 0.698 & 3.40 [2.98--3.89] & Baseline \\
 & Model 2: Intrinsic Mass Action ($f_S, f_R, f_C$) & 27.8\% & 0.780 & 6.29 [5.47--7.23] & --- \\
 & \textbf{Model 3: In Situ Spatial Game Assay} & \textbf{30.9\%} & \textbf{0.793} & \textbf{7.04 [6.11--8.11]} & \textbf{$F=18.70$ ($p < 10^{-15}$)} \\
\hline
\textbf{Breast IDC} & Model 1: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) & 16.1\% & 0.705 & 4.33 [3.85--4.88] & Baseline \\
 & Model 2: Intrinsic Mass Action ($f_S, f_R, f_C$) & 23.5\% & 0.754 & 5.04 [4.47--5.69] & --- \\
 & \textbf{Model 3: In Situ Spatial Game Assay} & \textbf{32.3\%} & \textbf{0.783} & \textbf{6.16 [5.45--6.96]} & \textbf{$F=71.95$ ($p < 10^{-15}$)} \\
\hline
\textbf{Melanoma} & Model 1: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) & 10.2\% & 0.674 & 3.52 [3.06--4.05] & Baseline \\
 & Model 2: Intrinsic Mass Action ($f_S, f_R, f_C$) & 32.5\% & 0.822 & 13.40 [11.39--15.76] & --- \\
 & \textbf{Model 3: In Situ Spatial Game Assay} & \textbf{36.8\%} & \textbf{0.832} & \textbf{16.33 [13.81--19.30]} & \textbf{$F=26.03$ ($p < 10^{-15}$)} \\
\hline
\end{tabular}%
}
\end{table}

While extrinsic microenvironmental covariates alone (Model 1) account for $10.2\% - 16.1\%$ of raw proliferation variance ($w_{raw}$), fitting the game interaction model on the purified residual ($w_{purified}$) explains an additional $30.9\%$ (Lung), $32.3\%$ (Breast), and $36.8\%$ (Melanoma) of the decoupled biological variance. In cumulative terms, the joint two-stage architecture ($\text{Env} + \text{Game}$) explains $41.8\%$ (Lung SCC), $43.2\%$ (Breast IDC), and $43.2\%$ (Melanoma SKCM) of total raw proliferation variance in intact clinical tissue. Across all three cohorts, adding the 9 spatial game interaction parameters significantly outperformed intrinsic mass-action frequencies alone (incremental nested $F$-tests: $F = 18.70, 71.95, 26.03$; all $p < 1.11 \times 10^{-16}$).

### 3.3 Robustness to Quality Control Cutoffs and Spatial Neighborhood Radii
To test whether low-coverage sequencing noise or arbitrary neighborhood definitions distorted the assay, we performed extensive parameter perturbation audits (Figure 3B):
- **Quality Control Stringency:** Re-fitting Breast IDC across unfiltered data ($N=4,992$, $0$ UMI threshold), standard QC ($N=4,981$, $\ge 500$ UMI) and high stringency ($N=4,963$, $\ge 1,000$ UMI) produced virtually invariant payoff estimates ($A_{SS} = -0.996 \to -1.031$, $A_{SR} = +2.606 \to +2.769$, $A_{RS} = +1.991 \to +1.882$, $A_{RR} = -1.970 \to -1.969$), with model $R^2$ remaining stable between $32.29\%$ and $32.76\%$.
- **Neighborhood Radius Sensitivity:** Varying the contact graph from immediate Delaunay triangulation ($\sim 55\,\mu\text{m}$) to $k$-nearest neighbors ($k=6, 12, 18$, spanning up to $\sim 200\,\mu\text{m}$) preserved the exact sign structure ($A_{SS} < 0, A_{SR} > 0, A_{RS} > 0, A_{RR} < 0$) and the Hawk-Dove classification with zero sign flips.

### 3.4 Serial Section Replicate Concordance and Blind Out-Of-Sample Validation
To provide concrete proof of biological and technical reproducibility, we evaluated two consecutive serial sections cut at $\approx 10\,\mu\text{m}$ thickness from human invasive breast cancer Block A (Section 1: $N=3,798$ spots; Section 2: $N=3,987$ spots).

Independently fitting Section 1 and Section 2 demonstrated strong parameter concordance across the full $3 \times 3$ payoff matrix (Figure 3C):
- **Payoff Matrix Correlation:** Pearson **$r = +0.811$ ($p = 0.0079$)**, Spearman **$\rho = +0.633$ ($p = 0.067$)**, Mean Absolute Error (MAE) $= 2.35$.
- **Intrinsic Fitness Vector Correlation:** Pearson **$r = +0.790$**.
- **Equilibrium Preservation:** Both serial sections independently classified the tumor into the **Stag Hunt (Bistable)** regime, with separatrix tipping points within $4.0\%$ agreement ($x_{sep} = 0.787$ vs $0.827$).

Next, to provide validation completely immune to planar spatial autocorrelation (Cliff & Ord, 1981), we executed a rigorous **Blind Out-Of-Sample (OOS) Generalization Test**: the spatial game model was trained exclusively on Section 1, and then applied blindly to predict spot-level proliferation on the unseen, held-out Section 2 (Figure 3D):
- **Blind OOS Explained Variance:** The Spatial Game Assay achieved **$\text{OOS } R^2 = 10.20\%$**, substantially outperforming the non-spatial intrinsic baseline ($\text{OOS } R^2 = 6.21\%$, a **$+3.99\%$ variance gain** out-of-sample).
- **Blind Diagnostic Performance:** Blind Diagnostic Odds Ratio $= 2.30$ [95% CI: $2.03 - 2.61$], exact Fisher $p = 4.82 \times 10^{-12}$, ROC-AUC $= 0.653$, Sensitivity $= 60.3\%$, Specificity $= 60.3\%$.

![**Reproducibility Battery, Serial Section Concordance, and Blind Out-Of-Sample Validation.** (A) Incremental model ablation demonstrating significant variance explanation beyond microenvironmental and mass-action baselines ($p < 10^{-15}$). (B) Parameter stability across quality control filtering ($0, \ge 500, \ge 1,000$ UMI) and spatial neighborhood radii ($k \in \{6, 12, 18\}$). (C) Replicate concordance across independent consecutive serial histological sections ($\approx 10\,\mu\text{m}$ spacing, $r = +0.811$, $p = 0.0079$) and game classification preservation. (D) Blind Out-Of-Sample (OOS) validation on held-out serial section ($R^2_{OOS} = 10.20\%$, DOR $= 2.30$, exact Fisher $p = 4.82 \times 10^{-12}$).](../figures/figure_reproducibility_battery.png){width=84%}

### 3.5 Independent Ground-Truth Literature Concordance
To calibrate our in situ findings against external biological ground truth, we benchmarked the inferred game matrices against independent published studies (Table 1 and Figure 2):
- **Breast Cancer Coexistence & Clonal Mutualism:** Seminal studies in clonal cooperation established that subclonal cancer populations can maintain non-cell-autonomous equilibrium through paracrine cross-talk (Marusyk et al., 2014). In parallel, landmark in vitro evolutionary game co-cultures in NSCLC demonstrated that cell-cell interactions follow matrix games such as Hawk-Dove and Stag Hunt that are heavily modulated by CAFs (Kaznatcheev et al., 2019). Our *InSituGame* assay, operating on intact human clinical biopsies, empirically discovered **Hawk-Dove / Snowdrift dynamics** in breast IDC ($A_{SR} - A_{RR} = +4.576 > 0$, $A_{RS} - A_{SS} = +2.987 > 0$, stable internal ESS $x^* = 60.5\%$), validating that intact human spatial transcriptomics autonomously recapitulates genuine cellular game ecology.
- **Lung & Melanoma Bistability:** In NSCLC and melanoma models, spatial competition constrains resistance outgrowth under targeted therapy (Bacevic et al., 2017; Kemper et al., 2014; Gallaher et al., 2018). *InSituGame* identified this exact coordination/bistable structure (Stag Hunt with tipping separatrices $x_{sep} = 20.0\%$ in Lung SCC and $52.8\%$ in Melanoma SKCM).

### 3.6 In Silico Adaptive Therapy vs Continuous MTD
Coupling the inferred payoff matrices with non-linear replicator dynamics demonstrated the therapeutic utility of *InSituGame* (Figure 2, right-most column). Under continuous MTD chemotherapy ($d_S = 1.2, d_R = 0.7$), sensitive clones were rapidly eradicated within $t = 25$ cycles, causing competitive release and 100% resistant outgrowth. In contrast, patient-specific adaptive dosing—pulsing treatment only when sensitive fraction exceeded the patient's critical threshold—successfully stabilized the sensitive population above the separatrix tipping point ($x_{sep}$), sustaining long-term tumor control and preventing treatment failure.

---

## 4. Discussion

The clinical failure of oncology drug regimens is fundamentally an evolutionary problem. While mathematical oncology has established that tumor clones engage in frequency-dependent evolutionary games, the inability to extract game payoff matrices ($A$) directly from individual patients has prevented game theory from entering routine oncology practice.

*InSituGame* bridges this gap by transforming clinical spatial transcriptomics from a descriptive histological catalog into an actionable evolutionary assay. Our study establishes five key advances:
1. **The "Personalized Biopsy Assay" Paradigm:** Rather than claiming a universal population census, *InSituGame* is designed as a patient-specific diagnostic assay—analogous to a clinical antibiogram (Kirby-Bauer disc diffusion assay). An antibiogram measures the unique drug resistance profile of an individual patient's bacterial isolate; similarly, *InSituGame* extracts the unique subclonal payoff matrix of an individual patient's biopsy to design tailored adaptive dosing.
2. **Decoupling Physical Confounders:** Prior models could not differentiate between intrinsic game interactions and physical microenvironmental drivers. By regressing out hypoxia and vascularization, we demonstrated that frequency-dependent game interactions explain a statistically significant fraction of proliferation variance ($p < 10^{-15}$).
3. **Ergodic Formalization:** Our Space-for-Time Ergodic Mapping Theorem provides the mathematical justification for using high-density spatial tissue variation to reconstruct dynamical replicator fields from a single clinical biopsy.
4. **Reproducibility Across Serial Histology:** The high concordance across consecutive serial cuts ($r = +0.811$) and robust blind out-of-sample prediction ($p < 10^{-11}$) prove that the extracted payoff parameters reflect stable tissue biology rather than technical slicing artifacts.
5. **Cross-Platform Generalizability:** While validated here on 10x Visium ($55\,\mu\text{m}$ spots), the linear formulation of *InSituGame* is platform-agnostic and natively extends to sub-cellular single-cell spatial platforms (10x Xenium, Vizgen MERSCOPE, Nanostring CosMx) by replacing spot neighborhoods with single-cell spatial Voronoi tessellations.

**Limitations and Future Directions:** While our simulations demonstrate the operational consequences of the payoff matrices, clinical translation will require coupling *InSituGame* with clinical pharmacokinetic/pharmacodynamic (PK/PD) models that incorporate drug clearance, tissue penetration, and organ toxicity ceilings. Prospective clinical trials in window-of-opportunity biopsy settings will be essential to validate adaptive dosing schedules guided by *InSituGame*.

In conclusion, *InSituGame* provides an open-source, mathematically grounded, and clinically audited framework that unlocks personalized evolutionary game theory for precision oncology.

---

## Ethical Approval
Not applicable. This study exclusively analyzes anonymized, publicly available archival clinical spatial transcriptomics datasets provided by 10x Genomics.

## Acknowledgments
The author warmly acknowledges the Telethon Institute of Genetics and Medicine (TIGEM), Pozzuoli, Italy, for scientific hospitality and providing a stimulating research environment during experimental degree thesis investigations.

## Conflict of Interest Statement
The author declares no competing financial or personal interests.

## Funding
This research received no external grant funding. Open Access publication fees are covered in full under the CRUI-CARE national agreement with the University of Naples Federico II.

## Data and Code Availability
All clinical spatial transcriptomics datasets analyzed in this study are publicly available from 10x Genomics. The open-source `insitugame` software library, complete pipeline scripts, and reproducibility benchmarks are archived on GitHub at [`https://github.com/Lukecele/spatial-game-assay`](https://github.com/Lukecele/spatial-game-assay) under the GNU General Public License v3.0 (GPL-3.0).

---

## References

1. Basanta, D., et al. (2012). Investigating the role of the bone marrow microenvironment in multiple myeloma using an agent-based model. *Blood*, 120(21), 2112.
2. Bacevic, K., et al. (2017). Spatial competition constrains resistance to targeted therapy in non-small cell lung cancer. *eLife*, 6, e27460.
3. Cliff, A. D., & Ord, J. K. (1981). *Spatial processes: models & applications*. Pion, London.
4. Enriquez-Navas, P. M., et al. (2016). Exploiting evolutionary principles in cancer therapy. *Cancer Research*, 76(7), 1687–1694.
5. Gallaher, J. A., et al. (2018). Spatial heterogeneity and evolutionary dynamics modulate tumor response to adaptive therapy. *Cancer Research*, 78(8), 2029–2039.
6. Gatenby, R. A., et al. (2009). Adaptive therapy. *Cancer Research*, 69(11), 4894–4903.
7. Kalluri, R. (2016). The biology and function of fibroblasts in cancer. *Nature Reviews Cancer*, 16(9), 582–598.
8. Kaznatcheev, A., et al. (2019). Fibroblasts and alectinib switch the evolutionary games played by non-small cell lung cancer. *Nature Ecology & Evolution*, 3(3), 450–456.
9. Kemper, K., et al. (2014). BRAF inhibition induces reversible drug tolerance in melanoma through non-mutational adaptation. *Cancer Cell*, 26(1), 60–74.
10. Maley, C. C., et al. (2017). Classifying the evolutionary and ecological features of neoplasms. *Nature Reviews Cancer*, 17(10), 605–619.
11. Marusyk, A., et al. (2014). Non-cell-autonomous driving of tumour growth supports sub-clonal heterogeneity. *Nature*, 514(7520), 54–58.
12. Nowak, M. A. (2006). *Evolutionary dynamics: exploring the equations of life*. Harvard University Press.
13. Robertson-Tessi, M., et al. (2015). Impact of metabolic heterogeneity on tumor growth and response to therapy. *Cancer Research*, 75(8), 1567–1579.
14. Woolf, B. (1955). On estimating the relation between blood group and disease. *Annals of Human Genetics*, 19(4), 251–253.
15. Zhang, J., et al. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. *Nature Communications*, 8(1), 1816.
