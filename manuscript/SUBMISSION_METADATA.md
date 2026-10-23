# Official Manuscript Submission Metadata Sheet

This metadata sheet contains all standardized fields required for submitting the manuscript to **Bioinformatics (Oxford Academic)** and uploading the preprint to **bioRxiv**.

---

### 1. Basic Article Information
- **Full Title:**  
  `InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Exploring Adaptive Therapy in Intact Human Biopsies`
- **Short Title / Running Head:**  
  `In Situ Spatial Evolutionary Game Assay`
- **Manuscript Type:**  
  `Original Paper / Original Research Article`
- **Primary Journal:**  
  `Bioinformatics (Oxford Academic)`
- **Subject Category / Track:**  
  - For Bioinformatics (OUP): `Original Paper - Systems Biology & Evolutionary Dynamics`  
  - For bioRxiv: `Bioinformatics` (or `Cancer Biology`)
- **Word Count (Main Text):** ~4,400 words  
- **Number of Display Items:** 3 Figures, 2 Tables, 17 References.

---

### 2. Author & Affiliation Details
- **Author 1 (Sole & Corresponding Author):**
  - **Given Name:** `Luca`
  - **Family Name:** `Celebrano`
  - **Academic Degree:** `MD Candidate / Medical Student`
  - **Institutional Affiliation:**  
    `Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Via Sergio Pansini 5, 80131 Naples, Italy`
  - **Institutional Email:** `l.celebrano@studenti.unina.it`
  - **ORCID iD:** [`0009-0007-0931-3555`](https://orcid.org/0009-0007-0931-3555)
  - **Role:** Conceptualization, Methodology, Mathematical Proofs, Software Development, Data Curation, Formal Analysis, Validation, Visualization, Writing - Original Draft, Writing - Review & Editing.

---

### 3. Abstract
```text
Motivation: Evolutionary game theory (EGT) provides a compelling mathematical foundation for cancer adaptive therapy: suppressing lethal drug-resistant clones by maintaining a stable competing subpopulation of treatment-sensitive cells. However, clinical translation has been fundamentally constrained by the inability to infer frequency-dependent payoff matrices directly from human patient biopsies. Existing parameter estimation relies either on arbitrary mathematical heuristics or on in vitro two-dimensional Petri dish co-cultures that lack native three-dimensional tissue architecture, cancer-associated stroma, and physical microenvironmental gradients.

Methods: Here, we introduce the In Situ Spatial Evolutionary Game Assay (InSituGame), an exploratory computational framework designed to estimate empirical subclonal payoff matrices, compute corresponding Nash equilibria, and simulate patient-specific adaptive dosing schedules directly from intact human tumor biopsies via spatial transcriptomics (10x Genomics Visium and CytAssist FFPE). To address microenvironmental confounding, InSituGame decouples local cell proliferation (w) from extrinsic niche covariates (hypoxia signature Z_hyp and vascular endothelia Z_vasc) prior to interaction modeling. By establishing the Space-for-Time Ergodic Mapping Framework (Proposition 1), we formally bridge instantaneous single-cell-cycle RNA transcription across spatial micro-neighborhoods to the continuous time derivatives of non-linear replicator dynamics. We resolve mathematical collinearity on the frequency simplex by formulating a fully identifiable 9-parameter absorbed payoff matrix (A_tilde_ij = b_i + A_ij) estimated via ordinary least squares without intercept.

Results: We evaluated the framework across N = 20,093 in-tissue spatial spots spanning human lung squamous cell carcinoma (N=3,858), invasive ductal breast carcinoma (N=4,992), cutaneous malignant melanoma (N=3,458), and two consecutive serial histological sections (~10 μm apart, N=3,798 and N=3,987). In an in silico numerical consistency benchmark, the 9-parameter bilinear estimator accurately recovered known payoff parameters under 20% additive noise (Pearson r = 0.9997, p = 1.13 x 10^-12, MAE = 0.023) and recovered the interior Nash coexistence equilibrium (x* = 60.53% true vs 60.30% estimated), verifying algebraic code correctness and numerical identifiability without matrix singularities. In clinical tissue, comprehensive model ablations on the decoupled target w_purified demonstrated lineage-specific interaction architectures: in breast IDC, the spatial game model (M_3, R^2 = 32.29%) significantly improved over linear regional composition (M_2c, R^2 = 26.58%; incremental F = 104.99, p < 1.11 x 10^-16), yielding a +7.26% gain in 5-fold spatial block cross-validation out-of-fold generalizability (R^2_OOF = 27.02% vs 19.76%), while generic quadratics (M_2d, R^2_OOF = 28.46%) showed comparable predictive capacity with higher parameterization (10 vs 8 parameters). In melanoma SKCM, spatial game interactions provided a modest increment over regional composition (F = 17.16, p = 6.04 x 10^-14; R^2_OOF = 33.42% vs 32.88%). In lung SCC, spatial proliferation was predominantly governed by regional composition (M_2c = 30.77% vs M_3 = 30.86%, F = 1.30, p = 0.266, R^2_OOF = 26.07% vs 27.06%), classifying this specimen as indeterminate / dominated by regional composition rather than robust frequency-dependent game interactions. Coordinate permutation null tests confirmed spatial structure dependence across all cohorts (p < 0.001). Under the identifiable absorbed matrix A_tilde, spatial block bootstrap 95% confidence intervals confirmed robust Hawk-Dove / Snowdrift (Coexistence) dynamics in primary breast IDC (x* in [0.526, 0.725], point estimate 61.3%; Delta_S in [-3.93, -1.55], Delta_R in [-6.37, -2.58]) and suggestive Stag Hunt bistability in melanoma (point estimate x_sep = 48.1%, though the 95% bootstrap confidence interval for Delta_R spans zero: [-0.97, +3.57], reflecting selective parameter uncertainty). In lung SCC, bootstrap confidence intervals for Delta_R likewise spanned zero ([-1.685, +1.168]), verifying that this specimen cannot be reliably classified into a discrete game. Consecutive serial sections demonstrated parameter concordance (r = +0.811, p = 0.0079) and blind out-of-sample generalizability (R^2_OOS = 10.20% vs 6.21%, exact Fisher p = 4.82 x 10^-12).

Availability and Implementation: InSituGame is implemented as an open-source Python library (insitugame) licensed under the GNU General Public License v3.0 (GPL-3.0). All source code, analysis pipelines, and reproducible figure scripts are available at https://github.com/Lukecele/spatial-game-assay.
```

---

### 4. Keywords
`Spatial Transcriptomics`, `Evolutionary Game Theory`, `Adaptive Therapy`, `Mathematical Oncology`, `Payoff Matrix`, `Replicator Dynamics`, `Nash Equilibrium`, `In Situ Assay`, `10x Genomics Visium`

---

### 5. Statements and Declarations
- **Conflict of Interest / Competing Interests:**  
  `The author declares no competing financial, personal, or commercial interests.`
- **Funding Statement:**  
  `This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors. Open Access publication fees are eligible for 100% coverage (0€ APC) under the national CRUI-CARE Transformative Agreement with the University of Naples Federico II.`
- **Data Availability Statement:**  
  `All clinical spatial transcriptomics datasets analyzed in this study are publicly available from 10x Genomics (Human Lung Squamous Cell Carcinoma CytAssist FFPE, Human Breast Carcinoma CytAssist FFPE, Human Cutaneous Melanoma CytAssist FFPE, and Invasive Ductal Carcinoma Visium Serial Sections 1 and 2). Complete raw matrix loaders and parsed spatial coordinate tables are bundled in the accompanying repository.`
- **Code Availability Statement:**  
  `All software and analysis routines are implemented in the open-source Python package 'insitugame', publicly available under the GNU General Public License v3.0 (GPL-3.0) at https://github.com/Lukecele/spatial-game-assay.`
- **Ethical Approval:**  
  `Not applicable. This study exclusively analyzes anonymized, publicly available archival clinical spatial transcriptomics datasets provided by 10x Genomics.`
- **Acknowledgments:**  
  `The author warmly acknowledges the Telethon Institute of Genetics and Medicine (TIGEM), Pozzuoli, Italy, for scientific hospitality and providing a stimulating research environment during experimental degree thesis investigations.`

---

### 6. Suggested Independent Expert Reviewers
1. **Prof. Alexander R. A. Anderson, Ph.D.**  
   - Affiliation: Chair, Department of Integrated Mathematical Oncology (IMO), H. Lee Moffitt Cancer Center & Research Institute, Tampa, FL, USA  
   - Email: `Alexander.Anderson@moffitt.org`  
   - Expertise: Spatial evolutionary mathematical modeling, adaptive therapy clinical trials, replicator dynamics in tumor progression.
2. **Prof. David Basanta, Ph.D.**  
   - Affiliation: Associate Member, Department of Integrated Mathematical Oncology (IMO), H. Lee Moffitt Cancer Center & Research Institute, Tampa, FL, USA  
   - Email: `David.Basanta@moffitt.org`  
   - Expertise: Evolutionary game theory in cancer, stroma-tumor interactions, mathematical oncology.
3. **Artem Kaznatcheev, Ph.D.**  
   - Affiliation: Oxford University / Moffitt Cancer Center  
   - Email: `artem.kaznatcheev@magd.ox.ac.uk`  
   - Expertise: Empirical cancer game assay development, frequency-dependent fitness measurement, theoretical computer science and game theory.
4. **Prof. Carlo C. Maley, Ph.D.**  
   - Affiliation: Director, Arizona Cancer Evolution Center (ACE), Arizona State University, Tempe, AZ, USA  
   - Email: `maley@asu.edu`  
   - Expertise: Clonal competition, evolutionary ecology of neoplasms, somatic evolution.
