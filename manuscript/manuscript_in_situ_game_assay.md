# InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Exploring Adaptive Therapy in Intact Human Biopsies

**Luca Celebrano**$^{1,*}$

$^{1}$ Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, 80131 Naples, Italy  
$^{*}$ **Corresponding Author:** Luca Celebrano  
Email: `l.celebrano@studenti.unina.it` | Public ORCID: [`0009-0007-0931-3555`](https://orcid.org/0009-0007-0931-3555)  
Repository: [`https://github.com/Lukecele/spatial-game-assay`](https://github.com/Lukecele/spatial-game-assay)

---

## Abstract

**Motivation:** Evolutionary game theory (EGT) provides a compelling mathematical foundation for cancer adaptive therapy: suppressing lethal drug-resistant clones by maintaining a stable competing subpopulation of treatment-sensitive cells. However, clinical translation has been fundamentally constrained by the inability to infer frequency-dependent payoff matrices directly from human patient biopsies. Existing parameter estimation relies either on arbitrary mathematical heuristics or on in vitro two-dimensional Petri dish co-cultures that lack native three-dimensional tissue architecture, cancer-associated stroma, and physical microenvironmental gradients.

**Methods:** Here, we introduce the **In Situ Spatial Evolutionary Game Assay** (*InSituGame*), an exploratory computational framework designed to estimate empirical subclonal payoff matrices, compute corresponding Nash equilibria, and simulate patient-specific adaptive dosing schedules directly from intact human tumor biopsies via spatial transcriptomics (10x Genomics Visium and CytAssist FFPE). To address microenvironmental confounding, *InSituGame* decouples local cell proliferation ($w$) from extrinsic niche covariates (hypoxia signature $Z_{hyp}$ and vascular endothelia $Z_{vasc}$) prior to interaction modeling. By establishing the **Space-for-Time Ergodic Mapping Framework (Proposition 1)**, we formally bridge instantaneous single-cell-cycle RNA transcription across spatial micro-neighborhoods to the continuous time derivatives of non-linear replicator dynamics. We resolve mathematical collinearity on the frequency simplex by formulating a fully identifiable 9-parameter absorbed payoff matrix ($\tilde{A}_{ij} = b_i + A_{ij}$) estimated via ordinary least squares without intercept.

**Results:** We evaluated the framework across $N = 20,093$ in-tissue spatial spots spanning human lung squamous cell carcinoma ($N=3,858$), invasive ductal breast carcinoma ($N=4,992$), cutaneous malignant melanoma ($N=3,458$), and two consecutive serial histological sections ($\approx 10\,\mu\text{m}$ apart, $N=3,798$ and $N=3,987$). In an in silico numerical consistency benchmark, the 9-parameter bilinear estimator accurately recovered known payoff parameters under 20% additive noise (Pearson $r = 0.9997$, $p = 1.13 \times 10^{-12}$, MAE $= 0.023$) and recovered the interior Nash coexistence equilibrium ($x^* = 60.53\%$ true vs $60.30\%$ estimated), verifying algebraic code correctness and numerical identifiability without matrix singularities. In clinical tissue, comprehensive model ablations on the decoupled target $w_{purified}$ demonstrated lineage-specific interaction architectures: in breast IDC, the spatial game model ($M_3$, $R^2 = 32.29\%$) significantly improved over linear regional composition ($M_{2c}$, $R^2 = 26.58\%$; incremental $F = 104.99$, $p < 1.11 \times 10^{-16}$), yielding a $+7.26\%$ gain in 5-fold spatial block cross-validation out-of-fold generalizability ($R^2_{\text{OOF}} = 27.02\%$ vs $19.76\%$), while generic quadratics ($M_{2d}$, $R^2_{\text{OOF}} = 28.46\%$) showed comparable predictive capacity with higher parameterization (10 vs 8 parameters). In melanoma SKCM, spatial game interactions provided a modest increment over regional composition ($F = 17.16$, $p = 6.04 \times 10^{-14}$; $R^2_{\text{OOF}} = 33.42\%$ vs $32.88\%$). In lung SCC, spatial proliferation was predominantly governed by regional composition ($M_{2c} = 30.77\%$ vs $M_3 = 30.86\%$, $F = 1.30$, $p = 0.266$, $R^2_{\text{OOF}} = 26.07\%$ vs $27.06\%$), classifying this specimen as indeterminate / dominated by regional composition rather than robust frequency-dependent game interactions. Coordinate permutation null tests confirmed spatial structure dependence across all cohorts ($p < 0.001$). Under the identifiable absorbed matrix $\tilde{A}$, spatial block bootstrap 95% confidence intervals confirmed robust **Hawk-Dove / Snowdrift (Coexistence)** dynamics in primary breast IDC ($x^* \in [0.526, 0.725]$, point estimate $61.3\%$; $\tilde{\Delta}_S \in [-3.93, -1.55]$, $\tilde{\Delta}_R \in [-6.37, -2.58]$) and suggestive **Stag Hunt bistability** in melanoma (point estimate $x_{sep} = 48.1\%$, though the 95% bootstrap confidence interval for $\tilde{\Delta}_R$ spans zero: $[-0.97, +3.57]$, reflecting selective parameter uncertainty). In lung SCC, bootstrap confidence intervals for $\tilde{\Delta}_R$ likewise spanned zero ($[-1.685, +1.168]$), verifying that this specimen cannot be reliably classified into a discrete game. Consecutive serial sections demonstrated parameter concordance ($r = +0.811$, $p = 0.0079$) and blind out-of-sample generalizability ($R^2_{\text{OOS}} = 10.20\%$ vs $6.21\%$, exact Fisher $p = 4.82 \times 10^{-12}$).

**Availability and Implementation:** *InSituGame* is implemented as an open-source Python library (`insitugame`) licensed under the GNU General Public License v3.0 (GPL-3.0). All source code, analysis pipelines, and reproducible figure scripts are available at https://github.com/Lukecele/spatial-game-assay.

---

## 1. Introduction

Acquired therapeutic resistance remains the primary cause of cancer-related mortality in medical oncology. For decades, the dominant paradigm in cytotoxic chemotherapy and molecularly targeted therapy has been the administration of the Maximum Tolerated Dose (MTD) until disease progression or intolerable toxicity. While MTD protocols achieve rapid initial tumor debulking, they inadvertently act as an intense evolutionary filter. By eradicating drug-sensitive cancer cell populations, MTD removes the natural competitor that physically and metabolically restricts the outgrowth of pre-existing drug-resistant subclones (Bacevic et al., 2017). This evolutionary phenomenon—termed **competitive release**—frequently leads to therapeutic failure, as the repopulating tumor becomes entirely refractory to further treatment.

To address this evolutionary trap, mathematical and evolutionary cancer biologists formulated **adaptive therapy** (Gatenby et al., 2009; Enriquez-Navas et al., 2016; Zhang et al., 2017). Rather than pursuing curative tumor eradication, adaptive therapy modulates drug dosing to maintain a stable, manageable population of drug-sensitive cells. When treatment is withdrawn during drug holidays, the sensitive subpopulation outcompetes the resistant subpopulation due to phenotypic fitness trade-offs and non-mutational tolerance costs (Kemper et al., 2014; Maley et al., 2017). A pioneering pilot clinical trial in metastatic castrate-resistant prostate cancer ($N = 11$; Zhang et al., 2017) demonstrated that adaptive dosing of abiraterone acetate achieved a median time to progression of 27 months while using 53% less cumulative drug compared to historical continuous dosing controls.

Despite its clinical promise, a central translational challenge persists: **how can we objectively estimate the clonal interaction payoff matrix ($A$) for an individual patient's biopsy?** In evolutionary game theory (EGT), mathematical modeling of subclonal competition, clonal cooperation, and stroma-mediated niche ecology has provided foundational insights into cancer progression (Basanta et al., 2012; Nowak, 2006), where the per-capita growth rate of any subclone $i$ depends fundamentally on the frequency of surrounding phenotypes $j$ through the payoff matrix $A_{ij}$. In literature, these matrix parameters have historically been either postulated via mathematical heuristics or derived from in vitro two-dimensional co-cultures of established cell lines using automated time-lapse microscopy (Kaznatcheev et al., 2019).

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
Module scores were calculated as the unweighted mean expression of component genes. To obtain spot-level phenotypic frequencies satisfying the unit simplex condition $\sum_{i} f_{i,s} = 1$, module scores were transformed via a softmax function with explicit unit temperature $T = 1.0$:
$$f_{i,s} = \frac{\exp((M_{i,s} - \max_k M_{k,s})/T)}{\sum_{j \in \{S, R, C\}} \exp((M_{j,s} - \max_k M_{k,s})/T)}$$
The unit temperature $T = 1.0$ preserves the intrinsic natural variance and relative scale of log-normalized expression modules. Choosing $T \gg 1$ would smooth frequencies toward a degenerate uniform distribution ($1/3, 1/3, 1/3$), inflating multicollinearity in $X_9$, whereas setting $T \ll 1$ would force binary winner-take-all assignments, artificially extinguishing multi-cellular spot mixtures. We note that the constant-sum simplex constraint naturally induces negative covariance among resulting frequencies $f_i$ by mathematical construction.

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

### 2.5 Space-for-Time Ergodic Mapping Framework

To mathematically justify inferring dynamic temporal trajectories from static histological snapshots, we formulate the **Space-for-Time Ergodic Mapping Framework**:

> **Proposition 1 (Space-for-Time Ergodic Mapping in Spatial Games):** Let $\Omega \subset \mathbb{R}^2$ be an intact tumor tissue section. Let $w_i(\mathbf{r}, t_0)$ denote the instantaneous net per-capita proliferation rate of subclone $i$ at coordinate $\mathbf{r}$, approximated via steady-state cell-cycle transcript abundance ($w$). If the spatial micro-neighborhood distribution $\{\bar{\mathbf{f}}(\mathbf{r})\}_{\mathbf{r} \in \Omega}$ densely samples the frequency simplex $\Delta = \{\mathbf{x} \in [0, 1]^K : \sum_k x_k = 1\}$, then the continuous replicator vector field:
> $$\frac{dx_i}{dt} = x_i \left[ F_i(\mathbf{x}) - \bar{F}(\mathbf{x}) \right]$$
> is uniquely reconstructed across $\Delta$ up to strategy-independent additive constants from a single spatial snapshot at $t = t_0$.

*Mathematical and Biological Rationale:* In solid tumors, spatial clustering generates diverse local neighbor compositions $\bar{\mathbf{f}}(\mathbf{r})$ across tissue coordinates $\mathbf{r} \in \Omega$. Because thousands of spatial spots sample varied regions of the frequency simplex $\Delta$, estimating the empirical fitness response $w_i(\mathbf{r}) = \sum_j \tilde{A}_{ij} \bar{f}_j(\mathbf{r})$ reconstructs the instantaneous fitness landscape across phase space. Two essential nuances govern this mapping:
1. **Phenotypic Steady-State Proxy:** Single-spot transcriptomic cell-cycle scoring measures instantaneous mRNA machinery commitment (*MKI67, PCNA, TOP2A, CDK1*), which serves as a steady-state phenotypic proxy for the per-capita mitotic rate $\frac{1}{N_i}\frac{\partial N_i}{\partial t}$, assuming that transcriptional cell-cycle engagement translates linearly into cell division.
2. **Replicator Gauge Invariance:** In evolutionary dynamics, replicator trajectories are invariant under row-wise additive shifts ($F_i(\mathbf{x}) \to F_i(\mathbf{x}) + c(\mathbf{x})$ cancels identically in $F_i - \bar{F}$). Because absorbing intrinsic fitness into $\tilde{A}_{ij} = b_i + A_{ij}$ represents a strategy-specific constant that acts identically across simplex mixtures ($\sum_j \bar{f}_j = 1$), the relative selective gradients ($\tilde{\Delta}_S, \tilde{\Delta}_R$), Nash equilibria, and replicator phase orbits are strictly preserved.

### 2.6 Analytical Nash Equilibrium Classifier and Stromal Conditioning

The clinical trajectory of the tumor under therapy depends on the sub-game between Sensitive ($S$) and Resistant ($R$) clones. For the identifiable absorbed matrix $\tilde{A}$, the relative selective advantages are:
$$\tilde{\Delta}_S = \tilde{A}_{SS} - \tilde{A}_{RS}$$
$$\tilde{\Delta}_R = \tilde{A}_{RR} - \tilde{A}_{SR}$$

#### Stroma in the $2 \times 2$ Sub-Game Projection
When projecting the 3-strategy system $\{S, R, C\}$ onto the 2-clone tumor sub-game, stroma ($C$) is modeled as a conditioning microenvironmental background. The spot-level fitness equations for $S$ and $R$ evaluate to:
$$w_S = \tilde{A}_{SS} \bar{f}_S + \tilde{A}_{SR} \bar{f}_R + \tilde{A}_{SC} \bar{f}_C$$
$$w_R = \tilde{A}_{RS} \bar{f}_S + \tilde{A}_{RR} \bar{f}_R + \tilde{A}_{RC} \bar{f}_C$$
Evaluating the relative fitness differential between clones yields:
$$w_S - w_R = (\tilde{A}_{SS} - \tilde{A}_{RS})\bar{f}_S + (\tilde{A}_{SR} - \tilde{A}_{RR})\bar{f}_R + (\tilde{A}_{SC} - \tilde{A}_{RC})\bar{f}_C$$
Conditional on the observed mean tissue stromal fraction $\langle \bar{f}_C \rangle$, the term $(\tilde{A}_{SC} - \tilde{A}_{RC})\langle \bar{f}_C \rangle$ acts as an additive constant offset in relative fitness. Because replicator dynamics depend strictly on fitness differentials, this microenvironmental term shifts baseline net growth rates without altering the bilinear interaction curvature $(\tilde{\Delta}_S, \tilde{\Delta}_R)$, enabling rigorous projection onto the $2 \times 2$ clonal game.

In the 2-clone tumor sub-game, we denote $x_S$ as the sensitive clone frequency and $x_R = 1 - x_S$ as the resistant clone frequency. The dynamical system admits four canonical Nash equilibrium regimes:
1. **Sensitive Dominance ($\tilde{\Delta}_S > 0, \tilde{\Delta}_R < 0$):** $S$ unconditionally outcompetes $R$.
2. **Resistant Dominance ($\tilde{\Delta}_S < 0, \tilde{\Delta}_R > 0$):** $R$ unconditionally outcompetes $S$.
3. **Stag Hunt / Bistable ($\tilde{\Delta}_S > 0, \tilde{\Delta}_R > 0$):** Both pure states are locally stable, separated by an unstable tipping threshold:
   $$x_{S,sep} = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R}, \quad x_{R,sep} = 1 - x_{S,sep} = \frac{\tilde{\Delta}_S}{\tilde{\Delta}_S + \tilde{\Delta}_R}$$
4. **Hawk-Dove / Snowdrift ($\tilde{\Delta}_S < 0, \tilde{\Delta}_R < 0$):** Both pure states are unstable, producing a stable interior coexistence equilibrium:
   $$x_S^* = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R} = \frac{\tilde{A}_{RR} - \tilde{A}_{SR}}{\tilde{A}_{SS} - \tilde{A}_{SR} - \tilde{A}_{RS} + \tilde{A}_{RR}}, \quad x_R^* = 1 - x_S^*$$

#### Therapeutic Modulation and Adaptive Dosing Rule
To model therapeutic modulation (Section 3.6), let $u(t) \in \{0, 1\}$ denote binary drug delivery and $d_i$ the clone-specific cytotoxicity rate ($d_S > d_R$). The treatment-modified per-capita fitness is:
$$F_i^{\text{treat}}(\mathbf{x}, t) = \sum_j \tilde{A}_{ij} x_j - u(t) d_i$$
Under therapy, the effective selective advantages become:
$$\tilde{\Delta}_S(D) = \tilde{\Delta}_S - u(t)(d_S - d_R), \quad \tilde{\Delta}_R(D) = \tilde{\Delta}_R + u(t)(d_S - d_R)$$
The operational objective of cancer adaptive therapy is to maintain a target sensitive population $x_S \ge x_S^{\text{target}}$ to competitively suppress resistance. On the resistant frequency axis ($x_R$), this corresponds to withholding or reducing therapy whenever resistance remains below the adaptive dosing threshold:
$$x_R^{\text{thresh}} = 1 - x_S^{\text{target}}$$
For instance, targeting the stable coexistence Nash equilibrium in Breast IDC ($x_S^* = 0.613$) establishes an adaptive dosing threshold of $x_R^{\text{thresh}} = 1 - 0.613 = 0.387$.

### 2.7 Biostatistical Evaluation Battery, Model Ablation Hierarchy & Controls
Continuous predictions were evaluated using an exhaustive hierarchy of nested models:
- **Model 1 ($M_1$):** Microenvironment alone ($Z_{hyp}, Z_{vasc}$) on raw proliferation $w_{raw}$.
- **Model 2 ($M_2$):** Linear mass action on purified fitness: $w_{purified} \sim f_S + f_R + f_C$ (rank $= 2$).
- **Model 2b ($M_{2b}$):** Non-spatial within-spot quadratic co-localization: $w_{purified} \sim \sum_{ij} f_i f_j$ (effective rank $= 5$).
- **Model 2c ($M_{2c}$):** Linear spot composition + Linear neighborhood composition: $w_{purified} \sim f_S + f_R + \bar{f}_S + \bar{f}_R$ (effective rank $= 4$). This crucial control evaluates whether variance explanation stems from smooth regional composition rather than bilinear game interactions.
- **Model 2d ($M_{2d}$):** Generic quadratics (within-spot and neighbor squares/cross-terms): $w_{purified} \sim \sum_{i \le j} f_i f_j + \sum_{i \le j} \bar{f}_i \bar{f}_j$ (effective rank $= 10$). This control tests whether empirical predictive gains reflect generic non-linear spatial smoothing rather than frequency-dependent subclonal cross-talk.
- **Model 3 ($M_3$):** InSituGame Spatial Game Assay: $w_{purified} \sim \sum_{ij} \tilde{A}_{ij} f_i \bar{f}_j$ (effective rank $= 8$ under simplex constraint; 9 design columns). The incremental $F$-test against $M_{2c}$ evaluates whether bilinear interaction terms provide significant explanatory power with $df_1 = 4$ ($8 - 4 = 4$) and $df_2 = N - 9$.

Statistical validation metrics include:
- **5-Fold Geographically Disjoint Spatial Block Cross-Validation:** Partitioning tissue coordinates into 5 spatial clusters via $k$-means to assess out-of-fold generalizability ($R^2_{\text{OOF}}$) across all models.
- **Spatial Block Bootstrap 95% Confidence Intervals:** 200 spatial cluster resamples to obtain non-parametric confidence bounds for absorbed payoff parameters $\tilde{A}$, selective advantages ($\tilde{\Delta}_S, \tilde{\Delta}_R$), and Nash equilibria ($x^*, x_{sep}$) robust to residual spatial autocorrelation.
- **Monte Carlo Coordinate Permutation Test:** Randomly shuffling spatial coordinates (1,000 permutations) while holding spot composition invariant to test topological dependence.
- **Moran's I Autocorrelation of Residuals:** Computed using the exact Cliff & Ord (1981) variance formula across all model specifications.
- **Nested $F$-Tests & Diagnostic Odds Ratios (DOR):** With Woolf 95% confidence intervals and exact Fisher's tests.

---

## 3. Results

### 3.1 In Silico Numerical Code Correctness and Identifiability Benchmark
To verify that the 9-parameter bilinear estimator exhibits numerical stability and parameter identifiability under simplex constraints, we conducted an in silico benchmark on synthetic spatial tissue grids ($N=1,000$ spots) with Dirichlet simplex frequencies ($F$) and spatial neighborhood smoothing ($\bar{F}$). A ground-truth Hawk-Dove interaction matrix ($\tilde{A}^{\text{true}}$, true interior coexistence $x^* = 60.53\%$) was specified to generate proliferation fitness with 20% additive Gaussian noise.

The 9-parameter estimator achieved near-perfect numerical recovery:
- **Matrix Correlation:** Pearson $r = 0.9997$ ($p = 1.13 \times 10^{-12}$), Mean Absolute Error (MAE) $= 0.0227$.
- **Equilibrium Recovery:** Recovered interior Nash equilibrium $x^* = 60.30\%$ (within $0.23\%$ of true $60.53\%$).
- **Matrix Comparison:**
  $$\tilde{A}^{\text{true}} = \begin{pmatrix} -1.00 & +2.60 & -1.30 \\ +2.00 & -2.00 & +0.20 \\ +0.60 & -0.30 & -0.20 \end{pmatrix}, \quad \tilde{A}^{\text{est}} = \begin{pmatrix} -1.01 & +2.59 & -1.31 \\ +2.01 & -2.00 & +0.14 \\ +0.68 & -0.29 & -0.21 \end{pmatrix}$$
This benchmark confirms that our mathematical formulation resolves the linear dependency inherent to the frequency simplex, verifying that the inverse game estimator accurately recovers parameters without matrix singularities or algebraic distortion under known model-matched conditions.

### 3.2 Pan-Cancer Empirical Game Landscape
Applying *InSituGame* across the clinical cancer cohorts revealed distinct, lineage-specific evolutionary game regimes under the identifiable absorbed matrix $\tilde{A}$ (Table 1 and Figure 1).

![**Spatial Phenotypic Architecture, Decoupled Fitness, and In Situ Payoff Estimation.** (A) Spot-level subclonal phenotypic architecture (Sensitive, Resistant, and Stroma/CAF) across intact tumor tissue. (B) Decoupled proliferation fitness ($w_{purified}$) removing extrinsic hypoxia and vascular confounders. (C) Empirical $3 \times 3$ absorbed payoff matrix ($\tilde{A}_{ij}$) capturing direct and stroma-mediated subclonal interactions.](../figures/figure1_spatial_game_decomposition.png){width=84%}

\begin{table}[htbp]
\centering
\caption{\textbf{Pan-Cancer In Situ Game Assay Parameter Summary Under the Identifiable Absorbed Matrix $\tilde{A}$.} Inferred selective advantages ($\tilde{\Delta}_S, \tilde{\Delta}_R$) and Nash equilibria ($x^*, x_{sep}$) reported with 95\% spatial block bootstrap confidence intervals (200 spatial cluster resamples). Lung SCC is classified as indeterminate because the incremental test against regional composition is not significant ($p = 0.266$) and its bootstrap confidence interval spans zero; its high diagnostic odds ratio (DOR $= 7.04$) is driven by linear regional composition rather than bilinear game terms. In Melanoma SKCM, while point estimates indicate Stag Hunt bistability ($x_{sep} = 0.481$), the bootstrap confidence interval for $\tilde{\Delta}_R$ spans zero ($[-0.97, +3.57]$), indicating suggestive bistability with parameter uncertainty.}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llcccc}
\hline
\textbf{Cancer Cohort} & \textbf{Spots ($N$)} & \textbf{Absorbed Payoffs $(\tilde{A}_{SS}, \tilde{A}_{SR}, \tilde{A}_{RS}, \tilde{A}_{RR})$} & \textbf{Inferred Game \& Equilibrium [95\% Bootstrap CI]} & \textbf{$R^2$ ($w_{purified}$)} & \textbf{DOR [95\% CI]} \\
\hline
\textbf{Lung SCC} & 3,858 & $+0.68, +0.05, -0.31, -0.19$ & Indeterminate / Regional ($p = 0.266$) [$\tilde{\Delta}_S: +0.46\text{ to }+1.56, \tilde{\Delta}_R: -1.69\text{ to }+1.17$] & 30.9\% & 7.04 [6.11--8.11] ($p < 10^{-11}$) \\
\textbf{Breast IDC} & 4,992 & $-0.74, +2.87, +2.19, -1.77$ & Hawk-Dove, $x^* = 0.613$ [0.526--0.725] [$\tilde{\Delta}_S: -3.93\text{ to }-1.55, \tilde{\Delta}_R: -6.37\text{ to }-2.58$] & 32.3\% & 6.16 [5.45--6.96] ($p < 10^{-11}$) \\
\textbf{Melanoma SKCM} & 3,458 & $+0.93, -1.78, -0.70, -0.27$ & Suggestive Stag Hunt, $x_{sep} = 0.481$ [0.168--0.637] [$\tilde{\Delta}_S: +0.47\text{ to }+2.60, \tilde{\Delta}_R: -0.97\text{ to }+3.57$] & 36.8\% & 16.33 [13.81--19.30] ($p < 10^{-11}$) \\
\hline
\end{tabular}%
}
\end{table}

#### Dynamical Interpretation of Inferred Regimes
1. **Breast IDC:** The absorbed matrix exhibits negative diagonal terms ($\tilde{A}_{SS} = -0.736, \tilde{A}_{RR} = -1.772$) and positive cross-terms ($\tilde{A}_{SR} = +2.866, \tilde{A}_{RS} = +2.188$), fulfilling the canonical conditions of a **Hawk-Dove / Snowdrift game** ($\tilde{\Delta}_S = -2.924 < 0, \tilde{\Delta}_R = -4.638 < 0$) with a stable interior coexistence equilibrium at $x^* = 61.3\%$. Spatial block bootstrap resamples confirm that this coexistence regime is 100% robust, with 95% bootstrap confidence intervals for both selective parameters strictly negative ($\tilde{\Delta}_S \in [-3.929, -1.551], \tilde{\Delta}_R \in [-6.366, -2.582]$) and the interior Nash equilibrium strictly confined within the simplex ($x^* \in [0.526, 0.725]$).
2. **Melanoma SKCM (Suggestive Stag Hunt / Bistability):** Point estimates of the absorbed matrix yield $\tilde{\Delta}_S = +1.633$ and $\tilde{\Delta}_R = +1.515$, suggesting a **Stag Hunt / Bistable game** with an unstable tipping threshold at $x_{sep} = 48.1\%$. However, non-parametric spatial block bootstrapping demonstrates that while $\tilde{\Delta}_S$ is robustly positive ($[+0.465, +2.601]$), the 95% confidence interval for $\tilde{\Delta}_R$ spans zero ($[-0.973, +3.571]$). While the upper portion of the resamples supports bistability, this parameter uncertainty indicates that discrete classification in melanoma remains suggestive rather than definitive.
3. **Lung SCC (Indeterminate / Regional Composition):** In unperturbed baseline tissue, point estimates show $\tilde{\Delta}_S = +0.983$ and $\tilde{\Delta}_R = -0.244$. However, spatial block bootstrap reveals that the 95% confidence interval for $\tilde{\Delta}_R$ spans zero ($[-1.685, +1.168]$), and Model 3 does not significantly improve over linear regional composition ($F = 1.30, p = 0.266$; $R^2_{\text{OOF}} = 26.07\%$ vs $27.06\%$). Consequently, this specimen cannot be reliably classified into a discrete game regime, and its spatial proliferation variance is predominantly governed by regional composition rather than frequency-dependent subclonal cross-talk. Importantly, although Model 3 does not provide an incremental game benefit, the spatial feature set as a whole exhibits high binarized diagnostic discrimination for proliferative spots (DOR $= 7.04$ [6.11--8.11], $p < 10^{-11}$), driven by the strong baseline explanatory power of regional composition ($M_{2c} = 30.77\%$).

![**Pan-Cancer In Situ Spatial Evolutionary Game Atlas.** Comprehensive spatial game interaction landscapes across human Lung Squamous Cell Carcinoma (Lung SCC, $N=3,858$), Invasive Ductal Breast Carcinoma (Breast IDC, $N=4,992$), and Cutaneous Malignant Melanoma (Melanoma SKCM, $N=3,458$). Inferred subclonal absorbed payoff matrices $\tilde{A}_{2 \times 2}$, simplex phase portraits, replicator vector fields, and clinical Nash equilibria (Hawk-Dove coexistence in breast IDC; Stag Hunt bistable tipping separatrix in melanoma; indeterminate / regional composition in lung SCC).](../figures/figure_pan_cancer_spatial_game_atlas.png){width=84%}

### 3.3 Model Ablation Battery and Spatial Controls
To establish whether the spatial game assay captures genuine game interaction beyond microenvironmental gradients, local spot co-localization, regional composition, and generic non-linearities, we conducted an incremental model ablation across all cohorts (Table 2 and Figure 3A).

\begin{table}[htbp]
\centering
\caption{\textbf{Comprehensive Model Ablation Battery and Spatial Controls Across Cohorts.} Model 1 fitted on raw proliferation $w_{raw}$; Models 2, 2b, 2c, 2d, and 3 fitted on purified residual $w_{purified}$. 5-fold Spatial Block CV evaluated on geographically disjoint clusters. Moran's I computed on model residuals. Degrees of freedom for incremental test vs $M_{2c}$: $df_1 = 4$ (effective parameter difference under simplex constraint: rank 8 for $M_3$ minus rank 4 for $M_{2c}$), $df_2 = N - 9$.}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llccccc}
\hline
\textbf{Cohort} & \textbf{Model Specification} & \textbf{In-Sample $R^2$} & \textbf{Block CV $R^2_{\text{OOF}}$} & \textbf{Residual Moran $I$} & \textbf{Incremental Test vs $M_{2c}$} \\
\hline
\textbf{Lung SCC} & $M_1$: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) on $w_{raw}$ & 15.82\% & --- & 0.459 & Baseline \\
 & $M_2$: Linear Mass Action ($f_S, f_R, f_C$) on $w_{purified}$ & 27.84\% & 23.75\% & 0.252 & --- \\
 & $M_{2b}$: Local Quadratic Co-localization ($f_i f_j$, non-spatial) & 28.26\% & 23.60\% & 0.253 & --- \\
 & $M_{2c}$: Regional Composition ($f_i$ + $\bar{f}_i$, linear) & 30.77\% & 27.06\% & 0.251 & Baseline \\
 & $M_{2d}$: Generic Quadratics ($f_i f_j + \bar{f}_i \bar{f}_j$, non-game) & 31.21\% & 26.40\% & 0.250 & --- \\
 & \textbf{$M_3$: In Situ Spatial Game ($\tilde{A}_{ij} f_i \bar{f}_j$, Delaunay)} & \textbf{30.86\%} & \textbf{26.07\%} & \textbf{0.250} & \textbf{$F = 1.30, p = 0.266$ (N.S.)} \\
 & \textit{Joint Model on $w_{raw}$ ($\text{Env} + \text{Spatial Game}$)} & \textit{54.30\%} & \textit{48.91\%} & \textit{0.250} & \textit{$+38.48\%$ over Env alone} \\
\hline
\textbf{Breast IDC} & $M_1$: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) on $w_{raw}$ & 16.12\% & --- & 0.470 & Baseline \\
 & $M_2$: Linear Mass Action ($f_S, f_R, f_C$) on $w_{purified}$ & 23.49\% & 17.02\% & 0.291 & --- \\
 & $M_{2b}$: Local Quadratic Co-localization ($f_i f_j$, non-spatial) & 28.38\% & 20.93\% & 0.221 & --- \\
 & $M_{2c}$: Regional Composition ($f_i$ + $\bar{f}_i$, linear) & 26.58\% & 19.76\% & 0.283 & Baseline \\
 & $M_{2d}$: Generic Quadratics ($f_i f_j + \bar{f}_i \bar{f}_j$, non-game) & 33.32\% & 28.46\% & 0.220 & --- \\
 & \textbf{$M_3$: In Situ Spatial Game ($\tilde{A}_{ij} f_i \bar{f}_j$, Delaunay)} & \textbf{32.29\%} & \textbf{27.02\%} & \textbf{0.220} & \textbf{$F = 104.99, p < 1.11 \times 10^{-16}$} \\
 & \textit{Joint Model on $w_{raw}$ ($\text{Env} + \text{Spatial Game}$)} & \textit{51.59\%} & \textit{45.83\%} & \textit{0.220} & \textit{$+35.47\%$ over Env alone} \\
\hline
\textbf{Melanoma} & $M_1$: Microenvironment Alone ($Z_{hyp}, Z_{vasc}$) on $w_{raw}$ & 10.23\% & --- & 0.524 & Baseline \\
 & $M_2$: Linear Mass Action ($f_S, f_R, f_C$) on $w_{purified}$ & 32.48\% & 29.06\% & 0.265 & --- \\
 & $M_{2b}$: Local Quadratic Co-localization ($f_i f_j$, non-spatial) & 33.44\% & 29.32\% & 0.260 & --- \\
 & $M_{2c}$: Regional Composition ($f_i$ + $\bar{f}_i$, linear) & 35.52\% & 32.88\% & 0.262 & Baseline \\
 & $M_{2d}$: Generic Quadratics ($f_i f_j + \bar{f}_i \bar{f}_j$, non-game) & 36.58\% & 32.99\% & 0.258 & --- \\
 & \textbf{$M_3$: In Situ Spatial Game ($\tilde{A}_{ij} f_i \bar{f}_j$, Delaunay)} & \textbf{36.78\%} & \textbf{33.42\%} & \textbf{0.256} & \textbf{$F = 17.16, p = 6.04 \times 10^{-14}$} \\
 & \textit{Joint Model on $w_{raw}$ ($\text{Env} + \text{Spatial Game}$)} & \textit{49.60\%} & \textit{44.20\%} & \textit{0.256} & \textit{$+39.37\%$ over Env alone} \\
\hline
\end{tabular}%
}
\end{table}

#### Lineage-Specific Interaction Insights and Control Comparisons
The benchmark across Models 2c, 2d, and 3 provides critical biological and statistical discernment:
- In **Breast IDC**, the spatial game model ($M_3$) significantly improves over linear regional composition ($M_{2c}$), explaining an additional $+5.71\%$ in-sample variance ($F = 104.99, p < 1.11 \times 10^{-16}$) and achieving a notable **$+7.26\%$ gain in spatial block cross-validation** ($R^2_{\text{OOF}} = 27.02\%$ vs $19.76\%$). When compared against generic quadratic terms ($M_{2d}$, $R^2_{\text{OOF}} = 28.46\%$), both models explain substantial non-linear variance, indicating that while bilinear cross-products ($f_i \bar{f}_j$) capture frequency-dependent ecological interactions grounded in evolutionary game theory, generic non-linear regional smoothing also accounts for spatial heterogeneity. Crucially, $M_3$ achieves comparable out-of-fold generalization to $M_{2d}$ while utilizing a more parsimonious (8 effective parameters vs 10 in $M_{2d}$) and mechanistically interpretable parameterization directly mapping to replicator dynamics and clinical frequency dependence.
- In **Melanoma SKCM**, the spatial game model provides a modest increment over regional composition ($F = 17.16, p = 6.04 \times 10^{-14}$; $R^2_{\text{OOF}} = 33.42\%$ vs $32.88\%$) and generic quadratics ($M_{2d} = 32.99\%$).
- In **Lung SCC**, linear regional composition ($M_{2c}$) explains $30.77\%$ in-sample and $27.06\%$ out-of-fold variance. Adding bilinear game terms ($M_3 = 30.86\%$, $R^2_{\text{OOF}} = 26.07\%$) does not achieve statistical significance ($F = 1.30, p = 0.266$). In this specimen, spatial proliferation is predominantly governed by regional phenotypic composition rather than frequency-dependent cross-talk.
- **Topological Permutation Test:** Shuffling spot spatial coordinates across 1,000 Monte Carlo permutations collapsed explained variance to baseline levels across all cohorts ($p < 0.001$).
- **Spatial Autocorrelation and Residual Dependence:** Across models, Moran's I of proliferation residuals declined from $0.459\text{--}0.524$ ($M_1$) to $0.252\text{--}0.291$ ($M_2$) due to environmental residualization and phenotypic deconvolution, with minimal further drop in $M_3$ ($0.220\text{--}0.256$). Because residual spatial autocorrelation persists ($I \approx 0.22\text{--}0.26$), standard OLS $F$-tests are subject to variance inflation due to sample non-independence; consequently, our primary evidentiary foundation rests on 5-fold Spatial Block Cross-Validation and Spatial Block Bootstrap intervals rather than uncorrected parametric $p$-values.

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

![**Reproducibility Battery, Serial Section Concordance, and Blind Out-Of-Sample Validation.** (A) Model ablation on the decoupled target $w_{purified}$ across cohorts (Breast IDC: $p < 10^{-16}$, Melanoma SKCM: $p = 6.0 \times 10^{-14}$, Lung SCC: $p = 0.27$ N.S.). (B) Absorbed payoff parameter stability across quality control filtering stringency ($0, \ge 500, \ge 1,000$ UMI in Breast IDC), demonstrating invariant Hawk-Dove sign topology. (C) Replicate concordance across independent consecutive serial histological sections ($\approx 10\,\mu\text{m}$ spacing, $r = +0.811$, $p = 0.0079$) and game classification preservation. (D) Blind Out-Of-Sample (OOS) validation on held-out serial section ($R^2_{OOS} = 10.20\%$, DOR $= 2.30$, exact Fisher $p = 4.82 \times 10^{-12}$).](../figures/figure_reproducibility_battery.png){width=84%}

### 3.5 Independent Biological Literature Context
Our in situ findings align with independent experimental models:
- **Breast Cancer Coexistence & Clonal Mutualism:** Seminal studies in clonal cooperation established that subclonal cancer populations can maintain non-cell-autonomous equilibrium through paracrine cross-talk (Marusyk et al., 2014) and microenvironmental niche remodeling (Basanta et al., 2012). In parallel, in vitro evolutionary game co-cultures in non-small cell lung cancer demonstrated that cell-cell interactions follow matrix games such as Hawk-Dove and Stag Hunt that are heavily modulated by CAFs (Kaznatcheev et al., 2019). Our *InSituGame* assay identified **Hawk-Dove / Snowdrift dynamics** in primary breast IDC ($\tilde{\Delta}_S = -2.924, \tilde{\Delta}_R = -4.638$, stable internal ESS $x^* = 61.3\%$), consistent with non-cell-autonomous clonal mutualism in intact human tissue.
- **Melanoma Bistability:** In targeted therapy models, spatial competition constrains resistance outgrowth (Bacevic et al., 2017; Gallaher et al., 2018; Kemper et al., 2014). *InSituGame* identified a suggestive coordination/bistable structure in Melanoma SKCM (Stag Hunt with tipping separatrix $x_{sep} = 48.1\%$).

### 3.6 In Silico Adaptive Therapy vs Continuous MTD

Coupling the inferred absorbed matrices ($\tilde{A}$) with continuous replicator dynamics illustrated the theoretical consequences of different dosing strategies across simulated treatment cycles (Figure 2, right-most column):

1. **Hawk-Dove Coexistence (Breast IDC):** The patient biopsy exhibits an initial resistant fraction of $x_{R,0} = 62.4\%$ ($x_{S,0} = 37.6\%$). The unperturbed natural coexistence equilibrium resides at $x_S^* = 61.3\%$ ($x_R^* = 38.7\%$). Under continuous MTD chemotherapy ($d_S = 1.2, d_R = 0.0$), the stable interior equilibrium shifts toward higher resistance: $x_S^*(D) = 45.5\%$ ($x_R^*(D) = 54.5\%$). In Figure 2 (Col 4), the continuous MTD trajectory relaxes from $x_{R,0} = 62.4\%$ down toward $54.5\%$, rather than eradicating the tumor. In contrast, simulated adaptive therapy with a threshold of $x_R^{\text{thresh}} = 1 - 0.613 = 0.387$ conditionally halts treatment when resistance drops, stabilizing the subclonal mixture near its natural lower-resistance state ($x_R^* = 38.7\%$).
2. **Suggestive Stag Hunt Bistability (Melanoma SKCM):** The patient biopsy exhibits $x_{S,0} = 64.7\%$ and $x_{R,0} = 35.3\%$. The unperturbed tipping separatrix is $x_{S,sep} = 48.1\%$ ($x_{R,sep} = 51.9\%$). Because $x_{S,0} > 48.1\%$, the untreated tumor naturally resides in the sensitive-dominated attractor basin. Under continuous MTD chemotherapy ($d_S = 1.2, d_R = 0.0$), the effective separatrix shifts dramatically upward:
   $$x_{S,sep}(D) = \frac{\tilde{\Delta}_R + 1.2}{\tilde{\Delta}_S + \tilde{\Delta}_R} = \frac{1.515 + 1.2}{3.148} \approx 86.2\%$$
   which corresponds on the resistant axis to a collapse of the separatrix down to $x_{R,sep}(D) \approx 13.8\%$. Because the patient's initial sensitive frequency ($x_{S,0} = 64.7\%$) is now well below the treatment-elevated threshold ($86.2\%$), continuous MTD inadvertently forces the tumor across the separatrix into the resistant attractor basin, driving $x_R \to 100\%$ (competitive release). Conversely, adaptive therapy pulses treatment conditionally whenever $x_R$ approaches the safety threshold ($x_R^{\text{thresh}} = 0.42$), maintaining $x_S$ safely above the tipping point and averting resistance fixation.

#### Sensitivity of Induced Separatrices to Assumed Cytotoxicity Differentials
We note that when exploring therapy-induced game shifts (e.g. converting sensitive dominance to an induced Stag Hunt regime), effective selective advantages follow $\tilde{\Delta}_S(D) = \tilde{\Delta}_S - \delta$ and $\tilde{\Delta}_R(D) = \tilde{\Delta}_R + \delta$, where $\delta = d_S - d_R$. The induced tipping separatrix satisfies:
$$x_{sep}(\delta) = \frac{\tilde{\Delta}_R + \delta}{\tilde{\Delta}_S + \tilde{\Delta}_R}$$
In the Lung SCC specimen ($\tilde{\Delta}_S = +0.983, \tilde{\Delta}_R = -0.244$), setting $\delta = 0.30$ yields $x_{sep} = 0.076$, $\delta = 0.50$ yields $x_{sep} = 0.346$, and $\delta = 0.80$ yields $x_{sep} = 0.752$. Because this threshold varies widely as a function of the user-specified drug differential $\delta$, and because the baseline specimen is indeterminate ($p=0.266$), therapeutic modulation must be interpreted strictly as an illustrative parametric exploration rather than an established clinical regimen. We emphasize that cytotoxicity parameters ($d_S, d_R$) represent heuristic illustrative values and must be calibrated in vivo before translational use.

---

## 4. Discussion

The clinical failure of oncology regimens is fundamentally an evolutionary challenge. While mathematical oncology has established that tumor clones engage in frequency-dependent interactions, the inability to extract game payoff matrices directly from individual patient biopsies has limited translation into routine oncology practice.

*InSituGame* demonstrates the feasibility of transforming clinical spatial transcriptomics from a descriptive spatial catalog into an exploratory evolutionary assay. Our study highlights five primary observations:
1. **The "Personalized Biopsy Assay" Paradigm:** Rather than assuming universal game parameters across a histotype, *InSituGame* treats game estimation as an individual biopsy assay—conceptually analogous to a clinical antimicrobial susceptibility test (antibiogram). Just as an antibiogram measures the unique resistance profile of a bacterial isolate, *InSituGame* estimates the empirical payoff landscape of a specific patient's biopsy to inform adaptive dosing.
2. **Microenvironmental Decoupling:** Physical microenvironmental gradients confound proliferation. Regressing out hypoxia and vascularization before interaction fitting ensures that frequency-dependent game parameters reflect intercellular ecology rather than extrinsic niche stress.
3. **Ergodic Formalization:** The Space-for-Time Ergodic Mapping Framework (Proposition 1) provides a formal mathematical bridge for utilizing high-density spatial variation to reconstruct dynamical replicator fields from a single clinical snapshot up to gauge-invariant additive constants.
4. **Identifiable Mathematical Formulation:** Formulating the 9-parameter absorbed matrix $\tilde{A}$ resolves the exact collinearity inherent to the unit simplex, reducing Condition Numbers to $31\text{--}54$ and enabling unbiased parameter recovery in silico ($r = 0.9997$, MAE $= 0.023$).
5. **Discriminative Spatial Ablation and Non-Linear Controls:** By comparing spatial game models against linear regional composition ($M_{2c}$) and generic quadratics ($M_{2d}$), we demonstrated that bilinear game interactions provide substantial predictive power in breast carcinoma ($F = 104.99$, $+7.26\%$ out-of-fold gain over $M_{2c}$), while identifying that lung SCC spatial variation is primarily governed by regional composition ($F = 1.30, p = 0.266$). Furthermore, evaluating generic non-linearities showed that regional quadratic smoothing accounts for a substantial fraction of spatial variance, underscoring that bilinear cross-terms should be interpreted as a parsimonious, game-theoretically grounded parameterization rather than an exclusive mechanism of non-linear variance.

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

## Declaration of Generative AI and AI-Assisted Technologies in the Writing Process
During the preparation of this work, the author utilized LLM-based agentic tools (Claude and Gemini via Antigravity) for code development, data analysis pipelines, and manuscript formatting assistance. The author reviewed, edited, and validated all analyses and text, and takes full responsibility for the scientific integrity and conclusions of the manuscript.

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
