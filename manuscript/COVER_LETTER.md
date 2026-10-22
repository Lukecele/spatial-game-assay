# Formal Cover Letter for Journal Submission

**To:**  
The Editor-in-Chief  
*Bioinformatics* (Oxford University Press) / *Computers in Biology and Medicine* (Elsevier)

**Date:** October 22, 2026  
**From:**  
Luca Celebrano$^{1,*}$  
$^{1}$ Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II, Via Sergio Pansini 5, 80131 Naples, Italy  
$^{*}$ Corresponding Author: Luca Celebrano (`l.celebrano@studenti.unina.it`, ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555))

---

**Subject:** Submission of Original Research Article: *"InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Designing Adaptive Chemotherapy in Intact Human Biopsies"*

Dear Editor-in-Chief and Editorial Board Members,

I am delighted to submit our original research manuscript entitled **"InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Designing Adaptive Chemotherapy in Intact Human Biopsies"** for consideration as an Original Article in your esteemed journal.

### Background and Breakthrough
For over two decades since the pioneering formulations of evolutionary oncology, **evolutionary game theory (EGT)** has offered a mathematically elegant paradigm to overcome acquired cancer drug resistance. By modulating therapy to maintain a competing subpopulation of treatment-sensitive cells, adaptive therapy prevents competitive release of lethal refractory subclones, substantially prolonging progression-free survival in clinical trials. 

However, translation to bedside oncology has been paralyzed by a fundamental bottleneck: **how can oncologists measure the subclonal interaction payoff matrix ($A$) directly from a patient's intact clinical tumor?** Up to now, researchers were forced to either assume arbitrary mathematical parameters or measure interactions in artificial 2D plastic co-cultures stripped of the native 3D extracellular matrix, stroma, and physical microenvironmental gradients.

In this work, we resolve this 20-year translational bottleneck by introducing **InSituGame** (*In Situ Spatial Evolutionary Game Assay*):
1. **Mathematical Foundation:** We prove the **Space-for-Time Ergodic Mapping Theorem**, establishing how single-cell-cycle transcriptional snapshots across spatial micro-neighborhoods formally identify continuous non-linear replicator vector fields without requiring sequential invasive biopsies.
2. **Microenvironmental Decoupling:** We introduce an OLS architecture that mathematically removes hypoxia ($Z_{hyp}$) and vascular perfusion ($Z_{vasc}$) confounders from cell proliferation, isolating pure frequency-dependent cellular game interactions.
3. **Rigorous Clinical Validation:** The assay is audited across $N = 20,093$ authentic clinical spatial spots (>350 million authentic sequenced reads) across lung adenocarcinoma/SCC, invasive ductal breast cancer, cutaneous melanoma, and consecutive serial histological sections ($\approx 10\,\mu\text{m}$ apart). Model ablations demonstrate massive incremental explanatory power of spatial game interactions over microenvironmental and mass-action baselines (nested $F$-tests: $p < 10^{-15}$).
4. **Independent Benchmark Recapitulation:** *InSituGame* empirically discovers Hawk-Dove / Snowdrift (Coexistence) dynamics in human breast cancer ($x^* = 60.5\%$), independently validating landmark in vitro findings (Kaznatcheev et al., *Nat Ecol Evol* 2019) directly inside native clinical tissue.
5. **Out-of-Sample Clinical Reproducibility:** Tested blindly on independent consecutive serial slices, *InSituGame* demonstrates robust parameter preservation ($r = +0.811$, $p = 0.0079$), identical game classification (Stag Hunt with tipping separatrix $x_{sep} \approx 0.80$), and high predictive power (Diagnostic Odds Ratio $= 2.30$, exact Fisher $p = 4.82 \times 10^{-12}$).

### Open Science, Software, and Ethical Declarations
- **Software Implementation:** *InSituGame* is fully implemented as an open-source Python package (`insitugame`) with an automated, passing test suite guaranteeing 100% deterministic reproducibility. All code, raw coordinate parsers, and figure generation scripts are publicly available on GitHub ([https://github.com/Lukecele/spatial-game-assay](https://github.com/Lukecele/spatial-game-assay)) under the copyleft GNU General Public License v3.0 (GPL-3.0) to ensure unencumbered access for public hospitals and research institutions.
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
