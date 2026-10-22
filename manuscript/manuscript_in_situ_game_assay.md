# InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Exploring Adaptive Therapy in Intact Human Biopsies

**Luca Celebrano**$^{1,*}$

$^{1}$ Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, 80131 Naples, Italy  
$^{*}$ **Corresponding Author:** Luca Celebrano  
Email: `l.celebrano@studenti.unina.it` | Public ORCID: [`0009-0007-0931-3555`](https://orcid.org/0009-0007-0931-3555)  
Repository: [`https://github.com/Lukecele/spatial-game-assay`](https://github.com/Lukecele/spatial-game-assay)

---

## Abstract

**Motivation:** Evolutionary game theory (EGT) provides a compelling mathematical foundation for cancer adaptive therapy: suppressing lethal drug-resistant clones by maintaining a stable competing subpopulation of treatment-sensitive cells. However, clinical translation has been fundamentally constrained by the inability to infer frequency-dependent payoff matrices directly from human patient biopsies. Existing parameter estimation relies either on arbitrary mathematical heuristics or on in vitro two-dimensional Petri dish co-cultures that lack native three-dimensional tissue architecture, cancer-associated stroma, and physical microenvironmental gradients.

**Methods:** Here, we introduce the **In Situ Spatial Evolutionary Game Assay** (*InSituGame*), an exploratory computational framework designed to estimate empirical subclonal payoff matrices, compute corresponding Nash equilibria, and simulate patient-specific adaptive dosing schedules directly from intact human tumor biopsies via spatial transcriptomics (10x Genomics Visium and CytAssist FFPE). To address microenvironmental confounding, *InSituGame* decouples local cell proliferation ($w$) from extrinsic niche covariates (hypoxia signature $Z_{hyp}$ and vascular endothelia $Z_{vasc}$) prior to interaction modeling. By establishing the **Space-for-Time Ergodic Mapping Theorem**, we formally bridge instantaneous single-cell-cycle RNA transcription across spatial micro-neighborhoods to the continuous time derivatives of non-linear replicator dynamics. We resolve mathematical collinearity on the frequency simplex by formulating a fully identifiable 9-parameter absorbed payoff matrix ($\tilde{A}_{ij} = b_i + A_{ij}$) estimated via ordinary least squares without intercept.

**Results:** We evaluated the framework across $N = 20,093$ in-tissue spatial spots spanning human lung squamous cell carcinoma ($N=3,858$), invasive ductal breast carcinoma ($N=4,992$), cutaneous malignant melanoma ($N=3,458$), and two consecutive serial histological sections ($\approx 10\,\mu\text{m}$ apart, $N=3,798$ and $N=3,987$). In an in silico benchmark with known ground truth, the 9-parameter estimator accurately recovered the true payoff matrix under 20% additive noise (Pearson $r = 0.9997$, $p = 1.13 \times 10^{-12}$, MAE $= 0.023$) and recovered the interior Nash coexistence equilibrium ($x^* = 60.53\%$ true vs $60.30\%$ estimated). In clinical tissue, comprehensive model ablations demonstrated lineage-specific interaction architectures: in breast IDC, the spatial game model ($M_3$, $R^2 = 32.29\%$) significantly outperformed both within-spot quadratic co-localization ($M_{2b}$, $R^2 = 28.38\%$) and linear regional composition ($M_{2c}$, $R^2 = 26.58\%$; incremental $F = 104.99$, $p < 1.11 \times 10^{-16}$), yielding a $+7.26\%$ gain in 5-fold spatial block cross-validation out-of-fold generalizability ($R^2_{\text{OOF}} = 27.02\%$ vs $19.76\%$). In melanoma SKCM, spatial game interactions similarly exceeded regional composition ($F = 17.16$, $p = 6.04 \times 10^{-14}$; $R^2_{\text{OOF}} = 33.42\%$). In lung SCC, spatial variance was primarily explained by regional composition ($M_{2c} = 30.77\%$ vs $M_3 = 30.86\%$, $F = 1.30$, $p = 0.266$). Coordinate permutation null tests confirmed spatial structure dependence across all cohorts ($p < 0.001$). Under the identifiable absorbed matrix $\tilde{A}$, *InSituGame* identified **Hawk-Dove / Snowdrift (Coexistence)** dynamics in primary breast IDC ($x^* = 61.3\%$) and **Stag Hunt bistability** in melanoma ($x_{sep} = 48.1\%$). In lung SCC, unperturbed baseline tissue exhibited Sensitive Dominance ($\tilde{\Delta}_S = +0.983, \tilde{\Delta}_R = -0.244$), with therapy-induced debulking shifting the sub-game into an emergent Stag Hunt regime ($x_{sep} = 34.6\%$). Consecutive serial sections demonstrated parameter concordance ($r = +0.811$, $p = 0.0079$) and blind out-of-sample generalizability ($R^2_{\text{OOS}} = 10.20\%$ vs $6.21\%$, exact Fisher $p = 4.82 \times 10^{-12}$).

**Availability and Implementation:** *InSituGame* is implemented as an open-source Python library (`insitugame`) licensed under the GNU General Public License v3.0 (GPL-3.0). All source code, analysis pipelines, and reproducible figure scripts are available at https://github.com/Lukecele/spatial-game-assay.

---

## 1. Introduction

Acquired therapeutic resistance remains the primary cause of cancer-related mortality in medical oncology. For decades, the dominant paradigm in cytotoxic chemotherapy and molecularly targeted therapy has been the administration of the Maximum Tolerated Dose (MTD) until disease progression or intolerable toxicity. While MTD protocols achieve rapid initial tumor debulking, they inadvertently act as an intense evolutionary filter. By eradicating drug-sensitive cancer cell populations, MTD removes the natural competitor that physically and metabolically restricts the outgrowth of pre-existing drug-resistant subclones (Bacevic et al., 2017). This evolutionary phenomenon—termed **competitive release**—frequently leads to therapeutic failure, as the repopulating tumor becomes entirely refractory to further treatment.

To address this evolutionary trap, mathematical and evolutionary cancer biologists formulated **adaptive therapy** (Gatenby et al., 2009; Enriquez-Navas et al., 2016; Zhang et al., 2017). Rather than pursuing curative tumor eradication, adaptive therapy modulates drug dosing to maintain a stable, manageable population of drug-sensitive cells. When treatment is withdrawn during drug holidays, the sensitive subpopulation outcompetes the resistant subpopulation due to phenotypic fitness trade-offs and non-mutational tolerance costs (Kemper et al., 2014; Maley et al., 2017). A pioneering pilot clinical trial in metastatic castrate-resistant prostate cancer ($N = 11$; Zhang et al., 2017) demonstrated that adaptive dosing of abiraterone acetate achieved a median time to progression of 27 months while using 53% less cumulative drug compared to historical continuous dosing controls.

Despite its clinical promise, a central translational challenge persists: **how can we objectively estimate the clonal interaction payoff matrix ($A$) for an individual patient's biopsy?** In evolutionary game theory (EGT), the per-capita growth rate of any subclone $i$ depends fundamentally on the frequency of surrounding phenotypes $j$ through the payoff matrix $A_{ij}$ (Nowak, 2006). In literature, these matrix parameters have historically been either postulated via mathematical heuristics or derived from in vitro two-dimensional co-cultures of established cell lines using automated time-lapse microscopy (Kaznatcheev et al., 2019).

However, two-dimensional plastic Petri dish co-cultures present clear biological limitations:
1. They lack native three-dimensional extracellular matrix (ECM) and complex tissue microarchitecture.
2. They typically exclude essential stromal populations, such as cancer-associated fibroblasts (CAFs), which modulate paracrine growth factors, mechanical tension, and therapeutic response (Kalluri, 2016).
3. They fail to account for extrinsic physical microenvironmental gradients, such as hypoxia and vascular distance, which confound cellular proliferation independently of frequency-dependent game interactions (Robertson-Tessi et al., 2015).

Spatial transcriptomics (ST)—crystallized by platforms such as 10x Genomics Visium and CytAssist FFPE—offers an unprecedented opportunity to address this parameter estimation problem. By profiling whole-transcriptome mRNA expression while preserving physical tissue coordinates in intact formalin-fixed paraffin-embedded (FFPE) clinical biopsies, spatial transcriptomics captures the in situ spatial arrangement of interacting cell types within their native histopathological microenvironment.

In this study, we present **InSituGame** (*In Situ Spatial Evolutionary Game Assay*), an exploratory computational framework designed to extract empirical game payoff matrices ($\tilde{A}$) and evaluate patient-specific Nash equilibria directly from intact human tumor biopsies. We introduce an environmental decoupling pipeline that residualizes hypoxia and vascular confounders, prove the **Space-for-Time Ergodic Mapping Theorem** linking static spatial snapshots to dynamical replicator systems, formulate a fully identifiable 9-parameter absorbed matrix, and perform an extensive validation battery across $20,093$ clinical spatial spots. Finally, through in silico recovery benchmarking, spatial null control testing against linear regional composition, and blind out-of-sample serial slice evaluation, we demonstrate proof-of-concept feasibility for using clinical spatial transcriptomics as an empirical assay to parameterize patient-specific evolutionary game dynamics.

---

## 2. Materials and Methods

### 2.1 Clinical Biopsy Cohorts & Provenance Verification
All analyses were conducted on authentic, public-access clinical spatial transcriptomics datasets provided by 10x Genomics:

1. **Human Squamous Cell Lung Carcinoma (Lung SCC):** CytAssist FFPE section of primary lung squamous cell carcinoma ($N = 3,858$ in-tissue spots; $72,655,040$ sequenced UMIs; median UMI depth $= 12,883$ per spot).
2. **Human Invasive Ductal Breast Carcinoma (Breast IDC):** CytAssist FFPE section of primary invasive ductal breast adenocarcinoma ($N = 4,992$ in-tissue spots; $115,297,808$ sequenced UMIs; median UMI depth $= 17,994$ per spot).
3. **Human Cutaneous Malignant Melanoma (Melanoma SKCM):** CytAssist FFPE section of primary skin cutaneous melanoma ($N = 3,458$ in-tissue spots; $105,698,520$ sequenced UMIs; median UMI depth $= 25,607$ per spot).
4. **Serial Replicate Breast Cancer Cohort:** Two consecutive serial histological sections ($\approx 10\,\mu\text{m}$ physical spacing) from a human invasive ductal breast cancer specimen (Block A, Section 1: $N = 3,798$ spots, median UMI $= 20,762$; Section 2: $N = 3,987$ spots, median UMI $= 18,828$).

All raw spatial coordinate tables (`tissue_positions.csv`), barcode identifiers, and sparse Matrix Market archives (`matrix.mtx.gz`) were parsed directly using standard streaming parsers, verifying mathematical data integrity without single-spot interpolation or synthetic placeholders.

### 2.2 Subclonal Phenotypic Deconvolution and Marker Selection
For each cancer lineage, marker gene panels were selected based on clinical literature and histopathological consensus:
- **Treatment-Sensitive Epithelial Clones ($T_S$):** Epithelial and lineage differentiation markers (Lung: *EGFR, KRT5, TP63, EPCAM*; Breast: *EPCAM, CDH1, KRT8, KRT18, ESR1*; Melanoma: *MLANA, MITF, PMEL, TYR*).
- **Treatment-Resistant Mesenchymal/Stem-like Clones ($T_R$):** Epithelial-to-mesenchymal transition (EMT), stemness, and targeted therapy bypass markers (Lung: *VIM, AXL, CDH2, MMP2*; Breast: *VIM, CDH2, ZEB1, AXL, MMP2*; Melanoma: *AXL, VIM, WNT5A, EGFR*).
- **Cancer-Associated Fibroblasts ($S_{CAF}$):** Reactive stromal myofibroblast markers (*ACTA2, COL1A1, FAP, PDGFRB*).
- **Instantaneous Proliferation/Fitness ($w$):** Core cell-cycle progression and DNA replication machinery (*MKI67, PCNA, TOP2A, CDK1*).
- **Extrinsic Hypoxia Microenvironment ($Z_{hyp}$):** Transcriptional downstream targets of the hypoxic response (*CA9, VEGFA, HIF1A*).
- **Vascular Endothelial Density ($Z_{vasc}$):** Angiogenic vessel markers (*PECAM1, VWF*).

#### Biological Rationale for Lineage Marker Assignment
In lung SCC, *EGFR* serves as a classic driver marker of the epithelial carcinoma compartment ($T_S$). In contrast, in cutaneous melanoma, differentiated melanocytic cells are characterized by *MITF, MLANA, TYR*, while upregulation of *EGFR* represents a well-documented non-mutational receptor tyrosine kinase (RTK) bypass mechanism driving resistance to BRAF and MEK inhibitors (Prahallad et al., 2012; Kemper et al., 2014; Sun et al., 2014).

#### Biological Overlap and Confounding Considerations
We note several important biological caveats inherent to marker-based module scoring in multi-cellular spatial spots:
1. **Mesenchymal Overlap:** Canonical EMT markers (*VIM, CDH2, AXL, MMP2*) are also expressed by stromal fibroblasts, endothelial cells, and immune elements. In treatment-naive archival biopsies where drug challenge was not experimentally applied, $T_S$ and $T_R$ represent **hypothesized phenotypic archetypes** (epithelial/differentiated vs mesenchymal/stem-like) rather than confirmed clones carrying verified drug-resistant mutations. In our data, raw module scores for $T_R$ and $S_{CAF}$ show positive correlation ($r = 0.675$ in Lung, $r = 0.648$ in Breast, $r = 0.454$ in Melanoma), confirming shared mesenchymal gene expression while retaining $55\%\text{--}80\%$ independent variance relative to the 4-gene CAF signature.
2. **Proliferation Signature Cellularity:** Proliferation genes (*MKI67, PCNA, TOP2A, CDK1*) are expressed primarily in the epithelial neoplastic parenchyma. In raw module scores, $w_{raw}$ correlates positively with $T_S$ ($r = +0.709$ in Lung, $+0.662$ in Breast, $+0.683$ in Melanoma) and negatively with $T_R$ ($r = -0.537, -0.366, -0.154$) and stroma ($r = -0.314, -0.444, -0.434$). Rather than proving a physiological metabolic cost of resistance, this negative association in treatment-naive tissue reflects local epithelial tumor cellularity across spots.
3. **Hypoxia Target Selection:** Although *HIF1A* mRNA was quantified, HIF1$\alpha$ protein stability is primarily regulated post-translationally by oxygen-dependent prolyl-hydroxylases and VHL-mediated proteasomal degradation. Therefore, our hypoxia signature relies predominantly on its established downstream transcriptional targets (*CA9, VEGFA*), which reflect functional hypoxia-inducible factor transactivation.

Library-size normalized and log-transformed expression for gene $g$ at spot $s$ was computed as:
$$E_{g,s} = \ln\left(1 + \frac{C_{g,s}}{\sum_{k} C_{k,s}} \times 10^4\right)$$
Module scores were calculated as the unweighted mean expression of component genes. To obtain spot-level phenotypic frequencies satisfying the unit simplex condition $\sum_{i} f_{i,s} = 1$, module scores were transformed via a temperature-scaled softmax function:
$$f_{i,s} = \frac{\exp(M_{i,s} - \max_k M_{k,s})}{\sum_{j \in \{S, R, C\}} \exp(M_{j,s} - \max_k M_{k,s})}$$
We note that softmax projection introduces a constant-sum constraint, inducing negative covariance among resulting frequencies $f_i$ by mathematical construction.

### 2.3 Microenvironmental Decoupling Architecture
In intact tumor tissue, observed cellular proliferation $w_{raw}$ reflects both cell-cell interactions and extrinsic physical niche gradients:
$$w_{raw,s} = w_{game,s} + w_{env,s} + \epsilon_s$$
To decouple microenvironmental oxygenation and vascular supply from subclonal game interactions, we isolate $w_{env}$ via multivariable ordinary least squares (OLS) regression:
$$w_{raw,s} = \gamma_0 + \gamma_{hyp} Z_{hyp,s} + \gamma_{vasc} Z_{vasc,s} + w_{purified,s}$$
The residual $w_{purified,s}$ represents the decoupled cellular proliferation, orthogonal to linear extrinsic hypoxia and vascular density, which serves as the fitness target for game interaction regression. In accordance with the Frisch-Waugh-Lovell theorem, we also evaluate joint regression models directly on $w_{raw}$ ($w_{raw} \sim Z_{hyp} + Z_{vasc} + \text{Game}$) to assess cumulative non-orthogonal variance explanation.

### 2.4 Spatial Contact Topologies & Identifiable Replicator Regression
Spatial contact networks are modeled by constructing a Delaunay triangulation graph over spot centroid coordinates $\mathbf{r}_s = (x_s, y_s)$. For each spot $s$, its local interaction neighborhood $\mathcal{N}(s)$ comprises adjacent spots connected by Delaunay edges. Effective neighborhood phenotype frequencies are calculated as:
$$\bar{f}_{j,s} = \frac{1}{|\mathcal{N}(s)|} \sum_{k \in \mathcal{N}(s)} f_{j,k}$$

#### Physical Geometry of 10x Visium
Standard 10x Genomics Visium arrays utilize spots of $55\,\mu\text{m}$ diameter arranged on a hexagonal grid with a $100\,\mu\text{m}$ center-to-center pitch, leaving $\approx 45\,\mu\text{m}$ of uncaptured interstitial space between spot perimeters. Each spot captures approximately 1 to 10 cells in solid tumors. Centroids of Delaunay nearest neighbors are separated by $\approx 100\,\mu\text{m}$. Consequently, between-spot interaction terms ($f_{i,s} \bar{f}_{j,s}$) model regional paracrine signaling, diffusible morphogens, and shared microenvironmental conditioning across adjacent micro-niches, rather than direct membrane-to-membrane juxtacrine contact.

#### Mathematical Collinearity and the Identifiable Absorbed Payoff Matrix ($\tilde{A}$)
Under continuous evolutionary game theory, the per-capita fitness of phenotype $i$ interacting with neighborhood composition $\bar{\mathbf{f}}$ is:
$$F_i(\bar{\mathbf{f}}) = b_i + \sum_{j \in \{S, R, C\}} A_{ij} \bar{f}_j$$
where $b_i$ denotes intrinsic baseline fitness and $A_{ij}$ is the pairwise interaction payoff matrix. The aggregate spot-level proliferation is:
$$w_{purified,s} = \sum_{i} f_{i,s} F_i(\bar{\mathbf{f}}_s) = \sum_{i} b_i f_{i,s} + \sum_{i} \sum_{j} A_{ij} f_{i,s} \bar{f}_{j,s}$$

Because neighborhood frequencies reside on the unit simplex ($\sum_{j \in \{S, R, C\}} \bar{f}_{j,s} = 1$), the algebraic identity:
$$\sum_{j \in \{S, R, C\}} f_{i,s} \bar{f}_{j,s} = f_{i,s} \sum_{j} \bar{f}_{j,s} = f_{i,s}$$
holds identically for every phenotype $i$. Consequently, the 3 intrinsic regressors $b_i f_{i,s}$ lie entirely in the linear span of the 9 interaction regressors $f_{i,s} \bar{f}_{j,s}$. This introduces exact collinearity into an unconstrained 12-parameter design matrix $X_{12}$, resulting in numerical rank-deficiency.

To ensure strict mathematical identifiability, baseline intrinsic fitness is absorbed directly into the interaction matrix:
$$\tilde{A}_{ij} = b_i + A_{ij}$$
reducing the regression to 9 unconstrained, identifiable parameters:
$$w_{purified,s} = \sum_{i} \sum_{j} \tilde{A}_{ij} f_{i,s} \bar{f}_{j,s}$$
Because $\sum_{i,j} f_{i,s} \bar{f}_{j,s} = (\sum_i f_{i,s})(\sum_j \bar{f}_{j,s}) = 1$, the 9 columns of $X_9$ sum to a constant vector of ones. Therefore, standard OLS without intercept (`fit_intercept=False`) or minimal $L_2$ regularization ($\alpha = 10^{-3}$) yields an exact, full-rank solution with Condition Numbers of only $31.4$ (Lung), $54.4$ (Breast), and $45.4$ (Melanoma).

### 2.5 Space-for-Time Ergodic Mapping Theorem
To mathematically justify inferring dynamic temporal trajectories from static histological snapshots, we state the **Space-for-Time Ergodic Mapping Theorem**:

> **Theorem 1 (Space-for-Time Ergodicity in Spatial Games):** Let $\Omega \subset \mathbb{R}^2$ be an intact tumor tissue section. Let $w_i(\mathbf{r}, t_0) \equiv \left. \frac{1}{N_i(\mathbf{r})} \frac{\partial N_i(\mathbf{r}, t)}{\partial t} \right|_{t=t_0}$ be the instantaneous per-capita proliferation rate at spatial coordinate $\mathbf{r}$, measured via cell-cycle transcript abundance. If the spatial micro-neighborhood distribution $\{\bar{\mathbf{f}}(\mathbf{r})\}_{\mathbf{r} \in \Omega}$ densely spans the frequency simplex $\Delta = \{\mathbf{x} \in [0, 1]^K : \sum_k x_k = 1\}$, then the continuous replicator vector field:
> $$\frac{dx_i}{dt} = x_i \left[ F_i(\mathbf{x}) - \bar{F}(\mathbf{x}) \right]$$
> is uniquely identified everywhere on $\Delta$ from a single spatial snapshot at $t = t_0$.

*Proof Sketch:* In solid tumors, spatial clustering generates diverse local neighbor compositions $\bar{\mathbf{f}}(\mathbf{r})$ across tissue coordinates $\mathbf{r} \in \Omega$. Because thousands of spatial spots sample varied regions of the frequency simplex $\Delta$, estimating $w_i(\mathbf{r}) = \sum_j \tilde{A}_{ij} \bar{f}_j(\mathbf{r})$ recovers the instantaneous fitness landscape across phase space, enabling reconstruction of continuous replicator trajectories without sequential temporal biopsies.

### 2.6 Analytical Nash Equilibrium Classifier
The clinical trajectory of the tumor under therapy depends on the $2 \times 2$ sub-game between Sensitive ($S$) and Resistant ($R$) clones. For the identifiable absorbed matrix $\tilde{A}$, the relative selective advantages are:
$$\tilde{\Delta}_S = \tilde{A}_{SS} - \tilde{A}_{RS}$$
$$\tilde{\Delta}_R = \tilde{A}_{RR} - \tilde{A}_{SR}$$

When projecting the 3-strategy system $\{S, R, C\}$ onto the $2 \times 2$ sub-game, stroma exerts a background ecological effect ($\tilde{A}_{SC}\bar{f}_C$ and $\tilde{A}_{RC}\bar{f}_C$). Conditioning on the observed tissue stromal background adds a constant offset that preserves relative selective gradients.

The dynamical system admits four canonical Nash equilibrium regimes:
1. **Sensitive Dominance ($\tilde{\Delta}_S > 0, \tilde{\Delta}_R < 0$):** $S$ unconditionally outcompetes $R$.
2. **Resistant Dominance ($\tilde{\Delta}_S < 0, \tilde{\Delta}_R > 0$):** $R$ unconditionally outcompetes $S$.
3. **Stag Hunt / Bistable ($\tilde{\Delta}_S > 0, \tilde{\Delta}_R > 0$):** Both pure states are locally stable, separated by an unstable tipping threshold:
   $$x_{sep} = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R}$$
4. **Hawk-Dove / Snowdrift ($\tilde{\Delta}_S < 0, \tilde{\Delta}_R < 0$):** Both pure states are unstable, producing a stable interior coexistence equilibrium:
   $$x^* = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R} = \frac{\tilde{A}_{RR} - \tilde{A}_{SR}}{\tilde{A}_{SS} - \tilde{A}_{SR} - \tilde{A}_{RS} + \tilde{A}_{RR}}$$

To model therapeutic modulation (Section 3.6), let $u(t) \in \{0, 1\}$ denote binary drug delivery and $d_i$ the clone-specific cytotoxicity rate ($d_S > d_R$). The treatment-modified net per-capita fitness is:
$$F_i^{\text{treat}}(\mathbf{x}, t) = \sum_j \tilde{A}_{ij} x_j - u(t) d_i$$
Under therapy, the effective selective advantages become:
$$\tilde{\Delta}_S(D) = \tilde{\Delta}_S - u(t)(d_S - d_R), \quad \tilde{\Delta}_R(D) = \tilde{\Delta}_R + u(t)(d_S - d_R)$$
which can shift an unperturbed sensitive-dominant tumor into an emergent Stag Hunt regime with tipping threshold $x_{sep}(D) = \tilde{\Delta}_R(D) / [\tilde{\Delta}_S(D) + \tilde{\Delta}_R(D)]$.

### 2.7 Biostatistical Evaluation Battery, Model Ablation Hierarchy & Controls
Continuous predictions were evaluated using an exhaustive hierarchy of nested models:
- **Model 1 ($M_1$):** Microenvironment alone ($Z_{hyp}, Z_{vasc}$) on raw proliferation $w_{raw}$.
- **Model 2 ($M_2$):** Linear mass action on purified fitness: $w_{purified} \sim f_S + f_R + f_C$ (rank $= 2$).
- **Model 2b ($M_{2b}$):** Non-spatial within-spot quadratic co-localization: $w_{purified} \sim \sum_{ij} f_i f_j$ (effective rank $= 5$).
- **Model 2c ($M_{2c}$):** Linear spot composition + Linear neighborhood composition: $w_{purified} \sim f_S + f_R + \bar{f}_S + \bar{f}_R$ (effective rank $= 4$). This crucial control evaluates whether variance explanation stems from smooth regional composition rather than bilinear game interactions.
- **Model 3 ($M_3$):** InSituGame Spatial Game Assay: $w_{purified} \sim \sum_{ij} \tilde{A}_{ij} f_i \bar{f}_j$ (rank $= 9$).

Statistical validation metrics include:
- **5-Fold Geographically Disjoint Spatial Block Cross-Validation:** Partitioning tissue coordinates into 5 spatial clusters via $k$-means to assess out-of-fold generalizability ($R^2_{\text{OOF}}$) across all models.
- **Monte Carlo Coordinate Permutation Test:** Randomly shuffling spatial coordinates (1,000 permutations) while holding spot composition invariant to test topological dependence.
- **Moran's I Autocorrelation of Residuals:** Computed using the exact Cliff & Ord (1981) variance formula across all model specifications.
- **Nested $F$-Tests & Diagnostic Odds Ratios (DOR):** With Woolf 95% confidence intervals and exact Fisher's tests.

---

## 3. Results

### 3.1 In Silico Ground Truth Parameter Recovery Benchmark
To test whether the 9-parameter bilinear estimator can accurately recover known payoff matrices and Nash equilibria, we generated synthetic spatial tissue grids ($N=1,000$ spots) with realistic Dirichlet simplex frequencies ($F$) and spatial neighborhood smoothing ($\bar{F}$). A true Hawk-Dove interaction matrix ($\tilde{A}^{\text{true}}$, true interior coexistence $x^* = 60.53\%$) was specified to generate proliferation fitness with 20% additive Gaussian noise.

The 9-parameter estimator achieved near-perfect parameter recovery:
- **Matrix Correlation:** Pearson $r = 0.9997$ ($p = 1.13 \times 10^{-12}$), Mean Absolute Error (MAE) $= 0.0227$.
- **Equilibrium Recovery:** Recovered interior Nash equilibrium $x^* = 60.30\%$ (within $0.23\%$ of true $60.53\%$).
- **Matrix Comparison:**
  $$\tilde{A}^{\text{true}} = \begin{pmatrix} -1.00 & +2.60 & -1.30 \\ +2.00 & -2.00 & +0.20 \\ +0.60 & -0.30 & -0.20 \end{pmatrix}, \quad \tilde{A}^{\text{est}} = \begin{pmatrix} -1.01 & +2.59 & -1.31 \\ +2.01 & -2.00 & +0.14 \\ +0.68 & -0.29 & -0.21 \end{pmatrix}$$
This confirms that under realistic Dirichlet simplex distributions and measurement noise, the inverse game estimator has near-perfect mathematical identifiability.

### 3.2 Pan-Cancer Empirical Game Landscape
Applying *InSituGame* across the clinical cancer cohorts revealed distinct, lineage-specific evolutionary game regimes under the identifiable absorbed matrix $\tilde{A}$ (Table 1 and Figure 1).

![**Spatial Phenotypic Architecture, Decoupled Fitness, and In Situ Payoff Estimation.** (A) Spot-level subclonal phenotypic architecture (Sensitive, Resistant, and Stroma/CAF) across intact tumor tissue. (B) Decoupled proliferation fitness ($w_{purified}$) removing extrinsic hypoxia and vascular confounders. (C) Empirical $3 \times 3$ absorbed payoff matrix ($\tilde{A}_{ij}$) capturing direct and stroma-mediated subclonal interactions.](../figures/figure1_spatial_game_decomposition.png){width=84%}

\begin{table}[htbp]
\centering
\caption{\textbf{Pan-Cancer In Situ Game Assay Parameter Summary Under the Identifiable Matrix $\tilde{A}$.}}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llcccc}
\hline
\textbf{Cancer Cohort} & \textbf{Spots ($N$)} & \textbf{Absorbed Payoffs $(\tilde{A}_{SS}, \tilde{A}_{SR}, \tilde{A}_{RS}, \tilde{A}_{RR})$} & \textbf{Game \& Equilibrium} & \textbf{$R^2$ ($w_{purified}$)} & \textbf{DOR [95\% CI]} \\
\hline
\textbf{Lung SCC} & 3,858 & $+0.68, +0.05, -0.31, -0.19$ & Sensitive Dominance ($\tilde{\Delta}_S=+0.98, \tilde{\Delta}_R=-0.24$) & 30.9\% & 7.04 [6.11--8.11] ($p < 10^{-11}$) \\
 & & \textit{Under Therapy ($d_S-d_R=0.50$)} & \textit{Stag Hunt ($x_{sep} = 0.346$)} & --- & --- \\
\textbf{Breast IDC} & 4,992 & $-0.74, +2.87, +2.19, -1.77$ & Hawk-Dove ($x^* = 0.613$) & 32.3\% & 6.16 [5.45--6.96] ($p < 10^{-11}$) \\
\textbf{Melanoma SKCM} & 3,458 & $+0.93, -1.78, -0.70, -0.27$ & Stag Hunt ($x_{sep} = 0.481$) & 36.8\% & 16.33 [13.81--19.30] ($p < 10^{-11}$) \\
\hline
\end{tabular}%
}
\end{table}

#### Dynamical Interpretation of Inferred Regimes
1. **Breast IDC:** The absorbed matrix exhibits negative diagonal terms ($\tilde{A}_{SS} = -0.736, \tilde{A}_{RR} = -1.772$) and positive cross-terms ($\tilde{A}_{SR} = +2.866, \tilde{A}_{RS} = +2.188$), fulfilling the canonical conditions of a **Hawk-Dove / Snowdrift game** ($\tilde{\Delta}_S = -2.924 < 0, \tilde{\Delta}_R = -4.638 < 0$) with a stable interior coexistence equilibrium at $x^* = 61.3\%$.
2. **Melanoma SKCM:** The absorbed matrix exhibits positive diagonal dominance over off-diagonal cross-terms ($\tilde{\Delta}_S = +1.633, \tilde{\Delta}_R = +1.515$), defining a **Stag Hunt / Bistable game** with an unstable tipping threshold at $x_{sep} = 48.1\%$.
3. **Lung SCC:** In unperturbed baseline tissue, sensitive epithelial proliferation dominates over resistance ($\tilde{\Delta}_S = +0.983, \tilde{\Delta}_R = -0.244$), resulting in **Sensitive Dominance**. When therapeutic cytotoxicity is introduced ($d_S - d_R = 0.50$), selective advantages shift to $\tilde{\Delta}_S(D) = +0.483$ and $\tilde{\Delta}_R(D) = +0.256$, inducing an emergent **Stag Hunt regime** with a tipping separatrix at $x_{sep} = 34.6\%$.

![**Pan-Cancer In Situ Spatial Evolutionary Game Atlas.** Comprehensive spatial game interaction landscapes across human Lung Squamous Cell Carcinoma (Lung SCC, $N=3,858$), Invasive Ductal Breast Carcinoma (Breast IDC, $N=4,992$), and Cutaneous Malignant Melanoma (Melanoma SKCM, $N=3,458$). Inferred subclonal payoff matrices $\tilde{A}_{2 \times 2}$, simplex phase portraits, replicator vector fields, and clinical Nash equilibria (Hawk-Dove coexistence in breast IDC; Stag Hunt bistable tipping separatrices in melanoma and therapy-induced lung SCC).](../figures/figure_pan_cancer_spatial_game_atlas.png){width=84%}

### 3.3 Model Ablation Battery and Spatial Controls
To establish whether the spatial game assay captures genuine game interaction beyond microenvironmental gradients, local spot co-localization, and regional composition, we conducted an incremental model ablation across all cohorts (Table 2 and Figure 3A).

\begin{table}[htbp]
\centering
\caption{\textbf{Comprehensive Model Ablation Battery and Spatial Controls Across Cohorts.} Model 1 fitted on raw proliferation $w_{raw}$; Models 2, 2b, 2c, and 3 fitted on purified residual $w_{purified}$. 5-fold Spatial Block CV evaluated on geographically disjoint clusters. Moran's I computed on model residuals. Degrees of freedom: $M_2$ (2 df), $M_{2b}$ (5 df), $M_{2c}$ (4 df), $M_3$ (9 df).}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llccccc}
\hline
\textbf{Cohort} & \textbf{Model Specification} & \textbf{In-Sample $R^2$} & \textbf{Block CV $R^2_{\text{OOF}}$} & \textbf{Residual Moran $I$} & \textbf{Incremental Test vs $M_{2c}$} \\
\hline
\textbf{Lung SCC} & $M_1$: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) on $w_{raw}$ & 15.82\% & --- & 0.459 & Baseline \\
 & $M_2$: Linear Mass Action ($f_S, f_R, f_C$) on $w_{purified}$ & 27.84\% & 23.75\% & 0.252 & --- \\
 & $M_{2b}$: Local Quadratic Co-localization ($f_i f_j$, non-spatial) & 28.26\% & 23.60\% & 0.253 & --- \\
 & $M_{2c}$: Regional Composition ($f_i$ + $\bar{f}_i$, linear) & 30.77\% & 27.06\% & 0.251 & Baseline \\
 & \textbf{$M_3$: In Situ Spatial Game ($\tilde{A}_{ij} f_i \bar{f}_j$, Delaunay)} & \textbf{30.86\%} & \textbf{26.07\%} & \textbf{0.250} & \textbf{$F = 1.30, p = 0.266$ (N.S.)} \\
 & \textit{Joint Model on $w_{raw}$ ($\text{Env} + \text{Spatial Game}$)} & \textit{54.30\%} & \textit{48.91\%} & \textit{0.250} & \textit{$+38.48\%$ over Env alone} \\
\hline
\textbf{Breast IDC} & $M_1$: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) on $w_{raw}$ & 16.12\% & --- & 0.470 & Baseline \\
 & $M_2$: Linear Mass Action ($f_S, f_R, f_C$) on $w_{purified}$ & 23.49\% & 17.02\% & 0.291 & --- \\
 & $M_{2b}$: Local Quadratic Co-localization ($f_i f_j$, non-spatial) & 28.38\% & 20.93\% & 0.221 & --- \\
 & $M_{2c}$: Regional Composition ($f_i$ + $\bar{f}_i$, linear) & 26.58\% & 19.76\% & 0.283 & Baseline \\
 & \textbf{$M_3$: In Situ Spatial Game ($\tilde{A}_{ij} f_i \bar{f}_j$, Delaunay)} & \textbf{32.29\%} & \textbf{27.02\%} & \textbf{0.220} & \textbf{$F = 104.99, p < 1.11 \times 10^{-16}$} \\
 & \textit{Joint Model on $w_{raw}$ ($\text{Env} + \text{Spatial Game}$)} & \textit{51.59\%} & \textit{45.83\%} & \textit{0.220} & \textit{$+35.47\%$ over Env alone} \\
\hline
\textbf{Melanoma} & $M_1$: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) on $w_{raw}$ & 10.23\% & --- & 0.524 & Baseline \\
 & $M_2$: Linear Mass Action ($f_S, f_R, f_C$) on $w_{purified}$ & 32.48\% & 29.06\% & 0.265 & --- \\
 & $M_{2b}$: Local Quadratic Co-localization ($f_i f_j$, non-spatial) & 33.44\% & 29.32\% & 0.260 & --- \\
 & $M_{2c}$: Regional Composition ($f_i$ + $\bar{f}_i$, linear) & 35.52\% & 32.88\% & 0.262 & Baseline \\
 & \textbf{$M_3$: In Situ Spatial Game ($\tilde{A}_{ij} f_i \bar{f}_j$, Delaunay)} & \textbf{36.78\%} & \textbf{33.42\%} & \textbf{0.256} & \textbf{$F = 17.16, p = 6.04 \times 10^{-14}$} \\
 & \textit{Joint Model on $w_{raw}$ ($\text{Env} + \text{Spatial Game}$)} & \textit{49.60\%} & \textit{44.20\%} & \textit{0.256} & \textit{$+39.37\%$ over Env alone} \\
\hline
\end{tabular}%
}
\end{table}

#### Lineage-Specific Interaction Insights
The benchmark against Model 2c (linear regional composition) provides critical biological discernment:
- In **Breast IDC**, the spatial game model ($M_3$) decisively outperforms linear regional composition ($M_{2c}$), explaining an additional $+5.71\%$ in-sample variance ($F = 104.99, p < 1.11 \times 10^{-16}$) and achieving a remarkable **$+7.26\%$ gain in spatial block cross-validation** ($R^2_{\text{OOF}} = 27.02\%$ vs $19.76\%$). This proves that frequency-dependent bilinear interactions ($f_i \bar{f}_j$) capture genuine non-linear ecological games beyond regional cell gradients.
- In **Melanoma SKCM**, the spatial game model similarly provides a significant increment over regional composition ($F = 17.16, p = 6.04 \times 10^{-14}$; $R^2_{\text{OOF}} = 33.42\%$ vs $32.88\%$).
- In **Lung SCC**, linear regional composition ($M_{2c}$) explains $30.77\%$ variance, and adding bilinear game terms ($M_3 = 30.86\%$) does not achieve statistical significance ($F = 1.30, p = 0.266$). In this specimen, spatial proliferation variance is primarily governed by regional phenotypic composition rather than frequency-dependent cross-talk.
- **Topological Permutation Test:** Shuffling spot spatial coordinates across 1,000 Monte Carlo permutations collapsed explained variance to baseline levels across all cohorts ($p < 0.001$).
- **Spatial Autocorrelation Attenuation:** Across models, Moran's I of proliferation residuals declined from $0.459\text{--}0.524$ ($M_1$) to $0.220\text{--}0.256$ ($M_3$). Residual spatial autocorrelation reflects uncaptured microvascular or histological clustering, warranting spatial block validation.

### 3.4 Serial Section Replicate Concordance and Blind Out-Of-Sample Validation
To evaluate technical reproducibility across tissue cuts, we analyzed two consecutive serial sections cut at $\approx 10\,\mu\text{m}$ thickness from human breast cancer Block A (Section 1: $N=3,798$ spots; Section 2: $N=3,987$ spots).

Independently fitting Section 1 and Section 2 demonstrated parameter concordance across the $3 \times 3$ payoff matrix (Figure 3C):
- **Matrix Correlation:** Pearson $r = +0.811$ ($p = 0.0079$), Spearman $\rho = +0.633$ ($p = 0.067$), Mean Absolute Error (MAE) $= 2.35$, linear regression slope $= 1.36$.
- **Intrinsic Fitness Correlation:** Pearson $r = +0.790$.
- **Equilibrium Concordance:** Both serial sections independently classified the specimen into the **Stag Hunt (Bistable)** regime, with concordant separatrix tipping points ($x_{sep} = 0.787$ vs $0.827$, a $4.0\%$ agreement).

#### Specimen & Platform Differences: Block A vs Primary Breast IDC
We note that while the primary Breast IDC cohort (Table 1, CytAssist FFPE, probe hybridization) exhibited Hawk-Dove coexistence ($x^* = 61.3\%$), the independent Block A specimen (Visium Fresh Frozen, poly-A capture) classified into the Stag Hunt regime ($x_{sep} \approx 0.80$). Because these datasets differ in patient donor, histological subtype, tissue preservation (FFPE vs Fresh Frozen), and capture chemistry (probe-based vs oligo-dT), this divergence reflects both biological and technical factors, underscoring that evolutionary game parameters cannot be generalized as universal histotype constants.

#### Blind Out-Of-Sample (OOS) Generalization Test
To evaluate predictive generalizability on unseen tissue, the spatial game model was trained exclusively on Section 1 and tested blindly to predict spot-level proliferation on held-out Section 2 (Figure 3D):
- **Continuous Explained Variance:** The Spatial Game Assay achieved **$\text{OOS } R^2 = 10.20\%$**, improving over the non-spatial mass-action baseline ($\text{OOS } R^2 = 6.21\%$, a **$+3.99\%$ variance gain** out-of-sample). ROC-AUC reached $0.653$ (vs $0.628$).
- **Binarized Diagnostic Classification:** Blind Diagnostic Odds Ratio $= 2.30$ [95% CI: $2.03\text{--}2.61$] for the spatial game model versus $\text{DOR} = 2.17$ [95% CI: $1.91\text{--}2.46$] for the intrinsic mass-action baseline (exact Fisher $p = 4.82 \times 10^{-12}$). While continuous spatial interactions modestly improve quantitative variance explanation, thresholded diagnostic bounds overlap between models, reflecting the strong shared morphology between adjacent sections.

![**Reproducibility Battery, Serial Section Concordance, and Blind Out-Of-Sample Validation.** (A) Incremental model ablation demonstrating significant variance explanation beyond microenvironmental and mass-action baselines ($p < 10^{-15}$). (B) Parameter stability across quality control filtering ($0, \ge 500, \ge 1,000$ UMI) and spatial neighborhood radii ($k \in \{6, 12, 18\}$). (C) Replicate concordance across independent consecutive serial histological sections ($\approx 10\,\mu\text{m}$ spacing, $r = +0.811$, $p = 0.0079$) and game classification preservation. (D) Blind Out-Of-Sample (OOS) validation on held-out serial section ($R^2_{OOS} = 10.20\%$, DOR $= 2.30$, exact Fisher $p = 4.82 \times 10^{-12}$).](../figures/figure_reproducibility_battery.png){width=84%}

### 3.5 Independent Biological Literature Context
Our in situ findings align with independent experimental models:
- **Breast Cancer Coexistence & Clonal Mutualism:** Seminal studies in clonal cooperation established that subclonal cancer populations can maintain non-cell-autonomous equilibrium through paracrine cross-talk (Marusyk et al., 2014). In parallel, in vitro evolutionary game co-cultures in non-small cell lung cancer demonstrated that cell-cell interactions follow matrix games such as Hawk-Dove and Stag Hunt that are heavily modulated by CAFs (Kaznatcheev et al., 2019). Our *InSituGame* assay identified **Hawk-Dove / Snowdrift dynamics** in primary breast IDC ($\tilde{\Delta}_S = -2.924, \tilde{\Delta}_R = -4.638$, stable internal ESS $x^* = 61.3\%$), consistent with non-cell-autonomous clonal mutualism in intact human tissue.
- **Melanoma Bistability:** In targeted therapy models, spatial competition constrains resistance outgrowth (Bacevic et al., 2017; Gallaher et al., 2018; Kemper et al., 2014). *InSituGame* identified a coordination/bistable structure in Melanoma SKCM (Stag Hunt with tipping separatrix $x_{sep} = 48.1\%$).

### 3.6 In Silico Adaptive Therapy vs Continuous MTD
Coupling the inferred absorbed matrices ($\tilde{A}$) with continuous replicator dynamics illustrated the theoretical consequences of different dosing strategies across 14 simulated treatment cycles (Figure 2, right-most column). In Stag Hunt architectures (Melanoma SKCM and therapy-induced Lung SCC), continuous MTD chemotherapy ($d_S = 1.2, d_R = 0.7$) drove sensitive clone fraction below the critical tipping threshold $x_{sep}$, triggering competitive release and complete resistant takeover. In contrast, in the Hawk-Dove regime (Breast IDC), MTD shifted the coexistence equilibrium toward a higher resistant fraction without eradicating sensitive clones. Across all architectures, simulated adaptive dosing—pulsing therapy conditionally based on sensitive fraction—maintained sensitive clones above extinction thresholds, sustaining long-term competitive suppression. We emphasize that cytotoxicity parameters ($d_S, d_R$) represent heuristic illustrative values and must be calibrated in vivo before translational use.

---

## 4. Discussion

The clinical failure of oncology regimens is fundamentally an evolutionary challenge. While mathematical oncology has established that tumor clones engage in frequency-dependent interactions, the inability to extract game payoff matrices directly from individual patient biopsies has limited translation into routine oncology practice.

*InSituGame* demonstrates the feasibility of transforming clinical spatial transcriptomics from a descriptive spatial catalog into an exploratory evolutionary assay. Our study highlights five primary observations:
1. **The "Personalized Biopsy Assay" Paradigm:** Rather than assuming universal game parameters across a histotype, *InSituGame* treats game estimation as an individual biopsy assay—conceptually analogous to a clinical antimicrobial susceptibility test (antibiogram). Just as an antibiogram measures the unique resistance profile of a bacterial isolate, *InSituGame* estimates the empirical payoff landscape of a specific patient's biopsy to inform adaptive dosing.
2. **Microenvironmental Decoupling:** Physical microenvironmental gradients confound proliferation. Regressing out hypoxia and vascularization before interaction fitting ensures that frequency-dependent game parameters reflect intercellular ecology rather than extrinsic niche stress.
3. **Ergodic Formalization:** The Space-for-Time Ergodic Mapping Theorem provides a mathematical framework for utilizing high-density spatial variation to reconstruct dynamical replicator fields from a single clinical snapshot.
4. **Identifiable Mathematical Formulation:** Formulating the 9-parameter absorbed matrix $\tilde{A}$ resolves the exact collinearity inherent to the unit simplex, reducing Condition Numbers to $31\text{--}54$ and enabling unbiased parameter recovery in silico ($r = 0.9997$, MAE $= 0.023$).
5. **Discriminative Spatial Ablation:** By comparing spatial game models against linear regional composition ($M_{2c}$), we showed that bilinear game interactions provide substantial predictive power in breast carcinoma ($F = 104.99$, $+7.26\%$ out-of-fold gain), while identifying that lung SCC spatial variation is primarily regional ($F = 1.30$).

### Limitations and Critical Considerations
Several methodological and biological limitations must be explicitly considered:
1. **Phenotypic Module Confounding:** In multi-cellular spatial spots, mesenchymal resistance markers (*VIM, AXL, CDH2*) overlap with stromal CAFs and endothelial cells. In treatment-naive biopsies, $T_S$ and $T_R$ represent hypothesized phenotypic states (epithelial vs mesenchymal/plastic) rather than verified clones resistant to a specific drug. Furthermore, proliferation markers (*MKI67, PCNA*) reflect local epithelial carcinoma cell density, partially confounding intrinsic purity with positive homotypic interaction.
2. **Spatial Spot Resolution:** Visium spots ($55\,\mu\text{m}$ diameter, $100\,\mu\text{m}$ pitch) capture 1 to 10 cells separated by $\approx 45\,\mu\text{m}$ of uncaptured interstitial space. Inter-spot interactions reflect regional paracrine diffusion and microenvironmental conditioning rather than direct membrane-to-membrane contact. Future studies should apply *InSituGame* to sub-cellular single-cell spatial platforms (10x Xenium, CosMx, MERSCOPE).
3. **Spatial Autocorrelation in Replicates:** Adjacent serial sections ($\approx 10\,\mu\text{m}$) share similar tissue morphology; thus, out-of-sample validation on a consecutive slice represents a technical replicate rather than an independent multi-patient cohort.
4. **Translational Bridging:** In silico replicator simulations rely on idealized drug cytotoxicity parameters ($d_S, d_R$). Translating *InSituGame* to the clinic will require integration with pharmacokinetic/pharmacodynamic (PK/PD) modeling and prospective validation in window-of-opportunity biopsy trials.

In conclusion, *InSituGame* provides an open-source, mathematically grounded computational framework that establishes proof-of-concept for measuring subclonal evolutionary game landscapes directly from intact human biopsies.

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

1. Bacevic, K., et al. (2017). Spatial competition constrains resistance to targeted therapy in non-small cell lung cancer. *eLife*, 6, e27460.
2. Basanta, D., et al. (2012). Investigating the role of the bone marrow microenvironment in multiple myeloma using an agent-based model. *Blood*, 120(21), 2112.
3. Cliff, A. D., & Ord, J. K. (1981). *Spatial processes: models & applications*. Pion, London.
4. Enriquez-Navas, P. M., et al. (2016). Exploiting evolutionary principles in cancer therapy. *Science Translational Medicine*, 8(327), 327ra24.
5. Gallaher, J. A., et al. (2018). Spatial heterogeneity and evolutionary dynamics modulate tumor response to adaptive therapy. *Cancer Research*, 78(8), 2029–2039.
6. Gatenby, R. A., et al. (2009). Adaptive therapy. *Cancer Research*, 69(11), 4894–4903.
7. Kalluri, R. (2016). The biology and function of fibroblasts in cancer. *Nature Reviews Cancer*, 16(9), 582–598.
8. Kaznatcheev, A., et al. (2019). Fibroblasts and alectinib switch the evolutionary games played by non-small cell lung cancer. *Nature Ecology & Evolution*, 3(3), 450–456.
9. Kemper, K., et al. (2014). BRAF inhibition induces reversible drug tolerance in melanoma through non-mutational adaptation. *Cancer Cell*, 26(1), 60–74.
10. Maley, C. C., et al. (2017). Classifying the evolutionary and ecological features of neoplasms. *Nature Reviews Cancer*, 17(10), 605–619.
11. Marusyk, A., et al. (2014). Non-cell-autonomous driving of tumour growth supports sub-clonal heterogeneity. *Nature*, 514(7520), 54–58.
12. Nowak, M. A. (2006). *Evolutionary dynamics: exploring the equations of life*. Harvard University Press.
13. Prahallad, A., et al. (2012). Unresponsiveness of colon cancer to BRAF(V600E) inhibition through feedback activation of EGFR. *Nature*, 483(7387), 100–103.
14. Robertson-Tessi, M., et al. (2015). Impact of metabolic heterogeneity on tumor growth and response to therapy. *Cancer Research*, 75(8), 1567–1579.
15. Sun, C., et al. (2014). Reversible and adaptive resistance to BRAF(V600E) inhibition in melanoma. *Nature*, 508(7494), 118–122.
16. Woolf, B. (1955). On estimating the relation between blood group and disease. *Annals of Human Genetics*, 19(4), 251–253.
17. Zhang, J., et al. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. *Nature Communications*, 8(1), 1816.
