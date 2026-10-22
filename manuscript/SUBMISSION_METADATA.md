# Official Manuscript Submission Metadata Sheet

This metadata sheet contains all standardized fields required for submitting the manuscript to **bioRxiv**, **Bioinformatics (Oxford Academic)**, or **Computers in Biology and Medicine (Elsevier)**.

---

### 1. Basic Article Information
- **Full Title:**  
  `InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Designing Adaptive Chemotherapy in Intact Human Biopsies`
- **Short Title / Running Head:**  
  `In Situ Spatial Evolutionary Game Assay`
- **Manuscript Type:**  
  `Original Paper / Original Research Article`
- **Subject Category / Track:**  
  - For bioRxiv: `Bioinformatics` (or `Cancer Biology`)  
  - For Bioinformatics (OUP): `Original Paper - Systems Biology & Evolutionary Dynamics`  
  - For Computers in Biology and Medicine: `Research Paper - Computational Oncology & AI`
- **Word Count (Main Text):** ~4,250 words  
- **Number of Display Items:** 3 Figures, 0 Tables (all data in embedded vectors / figures), 15 References.

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
Motivation: Evolutionary game theory (EGT) provides a compelling mathematical foundation for cancer adaptive therapy: suppressing lethal drug-resistant clones by maintaining a stable competing subpopulation of treatment-sensitive cells. However, clinical translation has been critically throttled by the inability to infer frequency-dependent payoff matrices (A) directly from human patient tumors. Existing parameter estimation relies on either arbitrary mathematical heuristics or non-physiological in vitro two-dimensional Petri dish co-cultures that strip away native three-dimensional tissue architecture, cancer-associated stroma, and physical microenvironmental gradients.

Methods: Here, we introduce the In Situ Spatial Evolutionary Game Assay (InSituGame), an end-to-end clinical-computational framework that extracts empirical subclonal payoff matrices, computes exact Nash equilibria, and simulates patient-specific adaptive dosing schedules directly from intact human tumor biopsies via spatial transcriptomics (10x Genomics Visium and CytAssist FFPE). To eliminate physical confounders, InSituGame mathematically decouples local cell proliferation (w) from extrinsic microenvironmental covariates (hypoxia signature Z_hyp and vascular endothelia Z_vasc) prior to interaction fitting. By establishing the Space-for-Time Ergodic Mapping Theorem, we formally bridge instantaneous single-cell-cycle RNA transcription across spatial micro-neighborhoods to the continuous time derivatives of non-linear replicator dynamics.

Results: We audited and validated the assay across N = 20,093 authentic in-tissue spatial spots (>350 million authentic sequenced reads) spanning human lung squamous cell carcinoma (N=3,858), invasive ductal breast carcinoma (N=4,992), cutaneous malignant melanoma (N=3,458), and two independent consecutive serial histological sections (~10 μm apart, N=3,798 and N=3,987). Model ablations demonstrated that spatial evolutionary game interactions explain a statistically massive increment of proliferation variance over microenvironmental gradients and non-spatial mass-action models (nested F-tests: p < 1.11 x 10^-16 across all cohorts). Payoff parameters remained invariant across quality control stringency thresholds (0, >= 500, >= 1,000 UMI) and spatial neighborhood radii (k in {6, 12, 18}). Independent consecutive serial sections demonstrated robust parameter concordance (Pearson r = +0.811, p = 0.0079) and preserved exact game structure (Stag Hunt with tipping separatrix x_sep ~ 0.80). Blind Out-Of-Sample (OOS) evaluation on held-out consecutive tissue demonstrated significant predictive power (OOS R^2 = 10.20%, Diagnostic Odds Ratio = 2.30 [95% CI: 2.03 - 2.61], exact Fisher p = 4.82 x 10^-12). Remarkably, InSituGame empirically discovered Hawk-Dove / Snowdrift (Coexistence) dynamics in breast carcinoma (x* = 60.5%), independently recapitulating landmark in vitro co-culture findings directly in situ.

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
