# Formal Cover Letter for Journal Submission

**To:**  
The Editor-in-Chief  
*Bioinformatics* (Oxford University Press) / *Computers in Biology and Medicine* (Elsevier)

**Date:** October 23, 2026  
**From:**  
Luca Celebrano$^{1,*}$  
$^{1}$ Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Via Sergio Pansini 5, 80131 Naples, Italy  
$^{*}$ Corresponding Author: Luca Celebrano (`l.celebrano@studenti.unina.it`, ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555))

---

**Subject:** Submission of Original Research Article: *"InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Exploring Adaptive Therapy in Intact Human Biopsies"*

Dear Editor-in-Chief and Editorial Board Members,

I am pleased to submit our original research manuscript entitled **"InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Exploring Adaptive Therapy in Intact Human Biopsies"** for consideration as an Original Article in your esteemed journal.

### Background and Scientific Rationale
For over two decades since the founding formulations of evolutionary oncology, **evolutionary game theory (EGT)** has provided a compelling mathematical foundation to counter acquired cancer drug resistance. By modulating therapy to maintain a competing subpopulation of treatment-sensitive cells, adaptive therapy prevents competitive release of refractory subclones, substantially prolonging progression-free survival in clinical pilot trials.

However, translation to bedside oncology has been constrained by a fundamental bottleneck: **how can we estimate the subclonal interaction payoff matrix ($A$) directly from a patient's intact clinical biopsy?** Historically, researchers were limited to either postulating mathematical heuristics or measuring interactions in artificial 2D in vitro co-cultures stripped of the native 3D extracellular matrix, stroma, and physical microenvironmental gradients.

In this work, we address this translational bottleneck by developing **InSituGame** (*In Situ Spatial Evolutionary Game Assay*), an exploratory computational framework with five core contributions:
1. **Mathematical Foundation & Ergodicity:** We prove the **Space-for-Time Ergodic Mapping Theorem**, establishing how single-cell-cycle transcriptional snapshots across spatial micro-neighborhoods formally identify continuous non-linear replicator vector fields without requiring sequential invasive biopsies.
2. **Identifiable Formulation on the Simplex:** We identify the exact algebraic collinearity between intrinsic fitness and interaction terms on the unit simplex, and formulate a fully identifiable 9-parameter absorbed payoff matrix ($\tilde{A}_{ij} = b_i + A_{ij}$) with condition number reduced by seven orders of magnitude ($\approx 31\text{--}54$ vs $\sim 4 \times 10^8$).
3. **Microenvironmental Decoupling & Spatial Null Controls:** We decouple hypoxia ($Z_{hyp}$) and vascular density ($Z_{vasc}$) confounders from cell proliferation, and benchmark spatial Delaunay interactions against a non-spatial within-spot quadratic co-localization control, demonstrating a statistically significant $+2.60\%$ to $+3.91\%$ increment in purified variance explanation ($p < 1.11 \times 10^{-16}$). When evaluated directly on raw proliferation, joint regression accounts for $49.60\%\text{--}54.30\%$ of total tissue variance.
4. **Reproducibility Across Serial Slices:** Evaluated across independent consecutive serial sections ($\approx 10\,\mu\text{m}$ apart), *InSituGame* demonstrates parameter concordance ($r = +0.811$, $p = 0.0079$), consistent game classification (Stag Hunt with tipping separatrix $x_{sep} \approx 0.80$), and blind out-of-sample predictive generalizability ($R^2_{\text{OOS}} = 10.20\%$ vs $6.21\%$ for mass-action baseline; Diagnostic Odds Ratio $= 2.30$ [95% CI: $2.03\text{--}2.61$], exact Fisher $p = 4.82 \times 10^{-12}$).
5. **Patient-Specific Game Individuality:** While the primary CytAssist breast IDC cohort exhibits Hawk-Dove coexistence ($x^* = 60.5\%$), independent fresh-frozen breast serial cuts display Stag Hunt bistability, demonstrating that evolutionary games are patient-specific ecological architectures that warrant individual biopsy profiling.

### Open Science, Software, and Ethical Declarations
- **Software Implementation:** *InSituGame* is fully implemented as an open-source Python package (`insitugame`) with an automated test suite. All code, raw data streaming parsers, and figure generation scripts are publicly available on GitHub ([https://github.com/Lukecele/spatial-game-assay](https://github.com/Lukecele/spatial-game-assay)) under the copyleft GNU General Public License v3.0 (GPL-3.0) to ensure unencumbered access for research institutions and hospitals.
- **Institutional Open Access Funding:** As an investigator affiliated with the School of Medicine, Università degli Studi di Napoli Federico II, open access publication is eligible for **0€ Article Processing Charges (APC)** through the national **CRUI-CARE (Conferenza dei Rettori delle Università Italiane)** Transformative Agreement with Oxford University Press and Elsevier.
- **Competing Interests:** The author declares no competing financial, personal, or commercial interests.
- **Originality:** This manuscript is original, has not been published previously, and is not currently under consideration for publication elsewhere.

### Suggested Independent Reviewers
To assist the editorial office in selecting qualified experts at the intersection of evolutionary game theory, mathematical oncology, and spatial omics, we suggest:
1. **Prof. Alexander R. A. Anderson, Ph.D.** — Chair, Integrated Mathematical Oncology (IMO), Moffitt Cancer Center, Tampa, FL, USA. Expert in evolutionary game theory, adaptive therapy, and spatial tumor modeling. Email: `Alexander.Anderson@moffitt.org`
2. **Prof. David Basanta, Ph.D.** — Associate Member, Integrated Mathematical Oncology (IMO), Moffitt Cancer Center, Tampa, FL, USA. Leading pioneer in evolutionary games in cancer ecosystems and stromal interactions. Email: `David.Basanta@moffitt.org`
3. **Artem Kaznatcheev, Ph.D.** — Mathematical Oncology & Evolutionary Game Theory, Oxford University / Moffitt Cancer Center. First author of the landmark 2019 *Nature Ecology & Evolution* paper on empirical cancer payoff matrices. Email: `artem.kaznatcheev@magd.ox.ac.uk`
4. **Prof. Carlo C. Maley, Ph.D.** — Arizona Cancer Evolution Center (ACE), Arizona State University, Tempe, AZ, USA. World leader in cancer evolutionary dynamics and clonal competition. Email: `maley@asu.edu`

Thank you for your consideration of our manuscript. We look forward to hearing from you.

Sincerely,

**Luca Celebrano**  
MD Candidate & Research Investigator  
Department of Clinical Medicine and Surgery, School of Medicine, Università degli Studi di Napoli Federico II  
Email: `l.celebrano@studenti.unina.it`  
ORCID: `0009-0007-0931-3555`  
