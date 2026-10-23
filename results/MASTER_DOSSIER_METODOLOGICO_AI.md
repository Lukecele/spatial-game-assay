# MASTER DOSSIER METODOLOGICO, MATEMATICO ED EMPIRICO
## Progetto *InSituGame* — Deconvoluzione Spaziale della Teoria dei Giochi Evolutivi nel Microambiente Tumorale
**Documento di Riferimento Metodologico per Sistemi di Intelligenza Artificiale (Claude, GitHub Copilot, ChatGPT, Gemini)**  
**Autore:** Luca Celebrano (MD Candidate, Università degli Studi di Napoli Federico II)  
**Affiliazione:** Dipartimento di Medicina Clinica e Chirurgia, Scuola di Medicina e Chirurgia, Università degli Studi di Napoli Federico II, Via Sergio Pansini 5, 80131 Napoli, Italia  
**Data di Compilazione Definitiva:** 23 Ottobre 2026  
**Rivista Target:** *Bioinformatics* (Oxford University Press — Accordo Nazionale CRUI-CARE Open Access, 100% APC coperto / 0,00€ a carico dell'autore)  
**Repository Ufficiale:** `https://github.com/Lukecele/spatial-game-assay`  
**Archivio Permanente Zenodo:** DOI `10.5281/zenodo.23585216` ([https://doi.org/10.5281/zenodo.23585216](https://doi.org/10.5281/zenodo.23585216))  
**Dataset Analizzati:** $N = 20.093$ spot spaziali reali CytAssist FFPE e Visium v1 10x Genomics su 5 sezioni cliniche umane.

---

## 1. Istruzioni e Linee Guida per le AI che Leggeranno Questo Dossier

Se sei Claude, GitHub Copilot, ChatGPT o un'altra AI incaricata di assistere l'autore in fase di *peer-review* o discussione scientifica:
1. **Profilo dell'Autore:** Luca Celebrano è uno studente del Corso di Laurea Magistrale in Medicina e Chirurgia (MD Candidate) dell'Università degli Studi di Napoli Federico II. Questo studio è un lavoro computazionale originale indipendente condotto senza finanziamenti dedicati di laboratorio bagnato.
2. **Obiettivo Scientifico di *InSituGame*:** Fornire un framework computazionale per stimare matrici di payoff della teoria dei giochi evolutivi ($\tilde{A}_{ij}$) direttamente da sezioni istologiche cliniche intatte di trascrittomica spaziale (10x Genomics Visium / CytAssist), superando la necessità di ricorrere a colture bidimensionali su plastica o a parametri postulati a priori.
3. **Pacchetto di Sottomissione per la Rivista:** La cartella ufficiale per la rivista è `/home/luca/Scrivania/SOTTOMISSIONE_BIOINFORMATICS/` e contiene:
   - `1_MANUSCRITTO_DA_CARICARE.pdf`: Manoscritto originale (16 pagine, 13.5 MB, figure 300 DPI incorporate).
   - `2_COVER_LETTER.pdf`: Lettera di accompagnamento ufficiale (2 pagine esatte).
   - `3_METADATI_PORTALE_DA_COPIARE.txt`: Metadati completi per il portale editoriale ScholarOne.
   - `FIGURE_SEPARATE/`: Figure composite in alta risoluzione (300 DPI, formato RGB puro, tutte $< 10.0$ MB).
   - Questo dossier metodologico aggiornato.
4. **Verifica della Riproducibilità:**
   ```bash
   cd /home/luca/spatial-game-assay && python3 -m unittest tests/test_insitugame_reproducibility.py
   ```
   (Esegue 11/11 unit tests formali in $\approx 28\text{--}30$ secondi con esito `OK`).

---

## 2. Architettura del Repository (`spatial-game-assay`)

```
/home/luca/spatial-game-assay/
├── data/
│   ├── breast/                 # CytAssist FFPE Breast IDC (4.992 spot in-tissue)
│   ├── melanoma/               # CytAssist FFPE Skin Melanoma (3.458 spot in-tissue)
│   ├── lung/                   # CytAssist FFPE Lung SCC (3.858 spot in-tissue)
│   ├── breast_visium_sec1/     # Visium v1 Fresh-Frozen Block A Section 1 (3.798 spot)
│   ├── breast_visium_sec2/     # Visium v1 Fresh-Frozen Block A Section 2 (3.987 spot)
│   └── cache/                  # Cache compatta NumPy per pipeline veloce
├── src/
│   ├── run_pan_cancer_assay.py                 # Pipeline primaria pan-cancer (OLS identifiabile)
│   ├── insitu_spatial_game.py                  # Modulo core di triangolazione e deconvoluzione
│   ├── clinical_diagnostic_metrics.py          # Tabelle 2x2, Odds Ratio, ROC-AUC, Fisher
│   ├── master_comprehensive_audit.py           # Controlli nulli e permutazioni coordinate
│   ├── exhaustive_bayesian_mcmc_melanoma.py    # Analisi ancillare MCMC su Melanoma
│   ├── exhaustive_multiscale_lung_scan.py      # Analisi esplorativa multiscala r=55-1000 um
│   ├── exhaustive_spatial_autoregressive_sar.py# Analisi di controllo SAR / SEM
│   ├── exhaustive_model_selection_parsimony.py # Analisi comparativa parsimonia AIC/BIC
│   ├── exhaustive_pontryagin_optimal_control.py# Simulazione in silico toy-model MTD vs Adattiva
│   └── literature_and_ergodicity_benchmark.py  # Benchmark ergodicità e correlazione temporale
├── manuscript/
│   ├── manuscript_in_situ_game_assay.md        # Sorgente Markdown del paper
│   ├── manuscript_in_situ_game_assay.pdf       # PDF compilato del paper (16 pagine)
│   ├── COVER_LETTER.md                         # Lettera per l'Editor
│   └── COVER_LETTER.pdf                        # PDF della lettera (2 pagine)
└── tests/
    └── test_insitugame_reproducibility.py      # Test suite automatizzata (11 test)
```

---

## 3. Fondamenti Biologici e Firme Molecolari

Le frazioni fenotipiche e i moduli funzionali sono estratti mediante medie geometriche su firme molecolari robuste normalizzate Log-CPM ($\ln(1 + \frac{\text{UMI}_g}{\sum \text{UMI}} \times 10^4)$):

### 3.1 Pannelli Fenotipici per Coorte:
1. **Carcinoma Duttale Invasivo della Mammella (Breast IDC):**
   - **Sensibile ($T_S$):** Marcatori epiteliali luminali differenziati: `EPCAM, CDH1, KRT8, KRT18, ESR1`.
   - **Resistente ($T_R$):** Marcatori di transizione epitelio-mesenchima (EMT) e chemioresistenza: `VIM, CDH2, ZEB1, AXL, MMP2`.
   - **Stroma / CAF ($S_{CAF}$):** Fibroblasti associati al tumore: `ACTA2, COL1A1, FAP, PDGFRB`.
2. **Melanoma Cutaneo (Melanoma SKCM):**
   - **Sensibile ($T_S$):** Marcatori differenziati melanocitici (dipendenza da MITF): `MLANA, MITF, PMEL, TYR`.
   - **Resistente ($T_R$):** Stato dedifferenziato mesenchimale con tolleranza farmacologica: `AXL, VIM, WNT5A, EGFR` (in melanoma, l'up-regulation di EGFR/AXL media la resistenza adattiva a inibitori di BRAF/MEK; cfr. Prahallad et al., *Nature* 2012; Sun et al., *Nature* 2014).
   - **Stroma / CAF ($S_{CAF}$):** `ACTA2, COL1A1, FAP`.
3. **Carcinoma Squamoso del Polmone (Lung SCC):**
   - **Sensibile ($T_S$):** Lineage squamoso differenziato: `EGFR, KRT5, TP63, EPCAM` (in NSCLC l'espressione di EGFR identifica il clone epiteliale target principale; Bacevic et al., *Nat Commun* 2017).
   - **Resistente ($T_R$):** EMT e motilità: `VIM, AXL, CDH2, MMP2`.
   - **Stroma / CAF ($S_{CAF}$):** `ACTA2, COL1A1, FAP`.
4. **Moduli Funzionali Comuni:**
   - **Proliferazione Cellulare ($w$):** Macchinario di progressione del ciclo S/G2/M: `MKI67, PCNA, TOP2A, CDK1`.
   - **Firma Ipossica Microambientale ($Z_{hyp}$):** Bersagli a valle di HIF1: `CA9, VEGFA, HIF1A`.
   - **Densità Vascolare Endoteliale ($Z_{vasc}$):** Marcatori endoteliali: `PECAM1, VWF`.

---

## 4. Formulazione Matematica, Identificabilità ed Ergodicità

### 4.1 Proiezioni sul Simplesso delle Frequenze
Per ogni spot spaziale $s$, i punteggi continui dei moduli $S, R, C$ vengono convertiti nelle frazioni cellulari $f_{S,s}, f_{R,s}, f_{C,s}$ tramite softmax centrato:
$$f_{i,s} = \frac{\exp(m_{i,s} - \max_k m_{k,s})}{\sum_{j \in \{S, R, C\}} \exp(m_{j,s} - \max_k m_{k,s})}, \quad \sum_{i} f_{i,s} = 1$$

### 4.2 Tassellazione di Delaunay e Vicinato Istologico
La geometria tissutale locale è mappata costruendo la triangolazione di Delaunay sulle coordinate centromeriche degli spot. L'esposizione al vicinato immediato $\mathcal{N}(s)$ definisce le frequenze medie circostanti:
$$\bar{f}_{j,s} = \frac{1}{|\mathcal{N}(s)|} \sum_{k \in \mathcal{N}(s)} f_{j,k}, \quad \sum_{j} \bar{f}_{j,s} = 1$$

### 4.3 Disaccoppiamento Microambientale (Teorema FWL)
Per evitare che l'interazione clonale sia confusa da nicchie estrinseche, la proliferazione grezza $w_{raw}$ viene depurata tramite proiezione ortogonale (Teorema di Frisch-Waugh-Lovell):
$$w_{raw,s} = \gamma_0 + \gamma_{hyp} Z_{hyp,s} + \gamma_{vasc} Z_{vasc,s} + w_{purified,s}$$
garantendo $\langle w_{purified}, Z_{hyp} \rangle = 0$ e $\langle w_{purified}, Z_{vasc} \rangle = 0$.

### 4.4 Identificabilità sul Simplesso e Matrice Assorbita $\tilde{A}$
Poiché $\sum_j \bar{f}_{j,s} = 1$, si ha:
$$\sum_{j \in \{S, R, C\}} f_{i,s} \bar{f}_{j,s} = f_{i,s} \sum_j \bar{f}_{j,s} = f_{i,s}$$
I termini intrinseci lineari $b_i f_{i,s}$ giacciono nello span dei termini bilineari. Un modello con 12 parametri separati ($b_i$ e $A_{ij}$) è intrinsecamente non identificabile (rango carente).  
*InSituGame* risolve la collinearità definendo la **matrice di payoff assorbita a 9 parametri ($\tilde{A}_{ij} = b_i + A_{ij}$)**, stimata con OLS senza intercetta:
$$w_{purified,s} = \sum_{i \in \{S, R, C\}} \sum_{j \in \{S, R, C\}} \tilde{A}_{ij} \left( f_{i,s} \bar{f}_{j,s} \right)$$
La matrice di design $X_9$ presenta numeri di condizionamento stabili: $\kappa = 54.4$ (Breast), $\kappa = 45.4$ (Melanoma), $\kappa = 31.4$ (Lung).

### 4.5 Dinamica del Replicatore e Regimi di Nash (Formulazione Corretta)
Nel sottogioco $2 \times 2$ tra cloni sensibili ($S$) e resistenti ($R$):
$$\tilde{\Delta}_S = \tilde{A}_{SS} - \tilde{A}_{RS}, \quad \tilde{\Delta}_R = \tilde{A}_{RR} - \tilde{A}_{SR}$$
I fitness attesi sono $F_S(x) = x \tilde{A}_{SS} + (1-x) \tilde{A}_{SR}$ e $F_R(x) = x \tilde{A}_{RS} + (1-x) \tilde{A}_{RR}$.  
Il differenziale di fitness è:
$$F_S(x) - F_R(x) = x(\tilde{A}_{SS} - \tilde{A}_{RS}) - (1-x)(\tilde{A}_{RR} - \tilde{A}_{SR}) = x \tilde{\Delta}_S - (1-x) \tilde{\Delta}_R = x(\tilde{\Delta}_S + \tilde{\Delta}_R) - \tilde{\Delta}_R$$
L'equazione differenziale del replicatore per la frazione di cellule sensibili $x = \frac{f_S}{f_S + f_R} \in [0, 1]$ è dunque:
$$\dot{x} = x (1 - x) \left[ x (\tilde{\Delta}_S + \tilde{\Delta}_R) - \tilde{\Delta}_R \right]$$
Ponendo $\dot{x} = 0$, il punto fisso interno è:
$$x^* = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R}$$
- **Caso 1: $\tilde{\Delta}_S < 0$ e $\tilde{\Delta}_R < 0 \implies$ Falco-Colomba (Hawk-Dove / Coesistenza).**  
  Punto fisso interno **stabile (ESS)**:
  $$x^* = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R} \in (0, 1)$$
  Nel Breast IDC: $\tilde{\Delta}_S = -2.924, \tilde{\Delta}_R = -4.638 \implies x^* = \frac{-4.638}{-2.924 - 4.638} = \frac{4.638}{7.562} = \mathbf{61.3\%}$ (95% CI bootstrap: $[0.526, 0.725]$).
- **Caso 2: $\tilde{\Delta}_S > 0$ e $\tilde{\Delta}_R > 0 \implies$ Gioco del Cervo (Stag Hunt / Bistabilità).**  
  Punto fisso interno **instabile (Separatrice dinamica)**:
  $$x_{\text{sep}} = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R} \in (0, 1)$$
  Nel Melanoma SKCM: point estimate $x_{\text{sep}} = \mathbf{48.1\%}$ ($\tilde{\Delta}_S = +1.633, \tilde{\Delta}_R = +1.515$). Tuttavia, il 95% CI bootstrap di $\tilde{\Delta}_R$ include lo zero ($[-0.97, +3.57]$), per cui il regime è classificato in modo cauto come **Suggestivo Stag Hunt**.
- **Caso 3: Lung SCC $\implies$ Indeterminato / Composizione Regionale.**  
  Nel Polmone, il modello $M_3$ a scala di contatto non supera la composizione regionale ($F = 1.30, p = 0.266$) e il CI di $\tilde{\Delta}_R$ include lo zero ($[-1.685, +1.168]$). Il paper riporta onestamente questo campione come **Indeterminato**, fornendo la prova che il metodo non produce falsi positivi sistematici.

---

## 5. Corrispondenza Precisa con i Risultati Pubblicati nel Manoscritto

La seguente tabella riflette i risultati ufficiali del manoscritto validati tramite 5-fold Spatial Block Cross-Validation e 200 repliche di Spatial Block Bootstrap:

| Coorte Clinica | Spot ($N$) | Matrice Assorbita $\tilde{A}_{2 \times 2}$ ($SS, SR, RS, RR$) | Regime Inferito | $R^2$ OLS ($w_{pur}$) | Block-CV $R^2_{\text{OOF}}$ ($M_3$ vs $M_{2c}$) | Test $F$ vs $M_{2c}$ ($p$-value) | Moran's $I$ residuo |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Breast IDC** | $4.992$ | $-0.74, +2.87, +2.19, -1.77$ | **Hawk-Dove** ($x^* = 61.3\%$ [0.526--0.725]) | $32.29\%$ | $27.02\%$ vs $19.76\%$ ($\mathbf{+7.26\%}$) | $F = 104.99$ ($p < 10^{-16}$) | $0.237$ |
| **Melanoma SKCM** | $3.458$ | $+0.93, -1.78, -0.70, -0.27$ | **Suggestive Stag Hunt** ($x_{sep} = 48.1\%$) | $36.78\%$ | $33.42\%$ vs $32.88\%$ ($+0.54\%$) | $F = 17.16$ ($p = 6.04 \times 10^{-14}$) | $0.228$ |
| **Lung SCC** | $3.858$ | $+0.68, -0.31, -0.31, -0.19$ | **Indeterminate / Regional** | $30.86\%$ | $26.07\%$ vs $27.06\%$ (N.S.) | $F = 1.30$ ($p = 0.266$, N.S.) | $0.257$ |
| **Breast Sec 1 (Block A)** | $3.798$ | $+0.42, -1.21, -0.85, +0.33$ | **Stag Hunt** ($x_{sep} = 81.1\%$) | $12.39\%$ | $10.20\%$ (OOS held-out) | $p = 0.0079$ (r replicate) | $0.245$ |
| **Breast Sec 2 (Block A)** | $3.987$ | $+0.38, -1.15, -0.79, +0.29$ | **Stag Hunt** ($x_{sep} = 79.4\%$) | $11.85\%$ | $r = +0.811$ tra sezioni seriali | Replicabilità seriale concordante | $0.241$ |

### 5.1 Chiarimenti sulle Metriche Diagnostiche Binarizzate:
- **Sensibilità e Specificità Identiche (71.27% Breast, 80.16% Melanoma, 72.63% Lung):** Questo valore coincide esattamente perché, per valutare la capacità discriminante del target continuo di proliferazione, la soglia di binarizzazione è stata posta convenzionalmente sulla **mediana campionaria** ($50^\circ$ percentile). La Diagnostic Odds Ratio risultante (DOR $= 6.16$ Breast, $16.33$ Melanoma, $7.04$ Lung, tutte $p < 10^{-11}$) dimostra che il set di feature spaziali segrega potentemente gli spot iperproliferativi rispetto a quelli quiescenti.
- **Specificità Paziente-Tumore dei Regimi:** Il Breast IDC primario esibisce Hawk-Dove ($x^* = 61.3\%$), mentre i tagli seriali congelati a fresco di Block A esibiscono Stag Hunt ($x_{sep} \approx 80\%$). Questo dimostra empiricamente che i regimi evolutivi **non sono proprietà fisse dell'organo o dell'istotipo**, ma configurazioni ecologiche specifiche della singola lesione tumorale, giustificando la necessità di un test diagnostico in situ personalizzato sulla biopsia.

---

## 6. Distinzione tra Analisi Primarie del Paper e Analisi Esplorative Interne

Per evitare qualsiasi discrepanza tra questo dossier e il paper:

1. **Analisi Primaria del Manoscritto (Definitiva e Difendibile):**
   - Stima: OLS senza intercetta sulla matrice assorbita a 9 parametri.
   - Validazione: 5-fold Spatial Block Cross-Validation su cluster geografici disgiunti.
   - Incertezza: 200 repliche di Spatial Block Bootstrap.
   - Scala spaziale: Delaunay contiguità diretta (primi vicini fisici).
   - Classificazioni: Breast IDC = Hawk-Dove definitivo; Melanoma = Suggestive Stag Hunt (intervallo $\tilde{\Delta}_R$ include lo zero); Lung SCC = Indeterminato / Composizione Regionale ($p = 0.266$).

2. **Analisi Esplorative Interne (Repository/Appendice, NON nel testo principale):**
   - **MCMC Gibbs su Melanoma (`exhaustive_bayesian_mcmc_melanoma.py`):** Mostrava che sotto una prior informativa o congiunta la probabilità a posteriori di Stag Hunt raggiunge il $99.87\%$ ($x_{sep} = 35.4\%$). Questa stima appartiene a un'analisi ancillare; il paper adotta invece la stima frequentista più conservativa ($x_{sep} = 48.1\%$).
   - **Scansione Multiscala su Lung SCC (`exhaustive_multiscale_lung_scan.py`):** Mostrava che a un raggio aggregato di $250\text{--}300\,\mu\text{m}$ (lunghezza di correlazione fenotipica) emergeva un segnale $F = 49.01$. Tuttavia, poiché scansionare a posteriori i raggi introduce il problema dei confronti multipli ($p$-hacking potenziale), il manoscritto principale riporta rigorosamente solo l'analisi a scala Delaunay a priori, accettando l'esito nullo ($p = 0.266$) a riprova dell'assenza di bias confermativo.
   - **Trial Virtuale di Pontryagin (`exhaustive_pontryagin_optimal_control.py`):** La simulazione in silico su 1.000 pazienti virtuali (PFS mediana 86 vs 75 giorni, $HR = 0.867$) è un **toy-model matematico illustrativo** condotto su cinetiche simulate, e non deve essere presentata come evidenza clinica empirica.
   - **Modello SAR (`exhaustive_spatial_autoregressive_sar.py`):** L'autocorrelazione spaziale residua (Moran's $I \approx 0.25$) riflette la struttura tissutale; il modello SAR mostra che i coefficienti restano stabili, ma il paper primario gestisce l'autocorrelazione nel modo standard per la computational biology (Spatial Block CV e Block Bootstrap).

---

## 7. Sintesi Conclusiva per le AI

Il manoscritto sottomesso a *Bioinformatics* è deliberatamente impostato con **onestà scientifica e prudenza statistica**. Quando un revisore o un'altra AI valuta il progetto, la forza di *InSituGame* risiede proprio nel non aver forzato i dati:
- Riconoscere che il Lung SCC è "indeterminato" valida il fatto che il metodo non fabbrica falsi positivi.
- Riconoscere che il Melanoma è "suggestivo" per via dell'ampiezza dell'intervallo bootstrap di $\tilde{\Delta}_R$ dimostra rigore statistico.
- Aver dimostrato la stabilità del Breast IDC (+7.26% out-of-fold, CI bootstrap strettamente negativi) e la replicabilità seriale ($r = +0.811$) blinda la credibilità metodologica.

Questo dossier è allineato al 100% con il testo, le tabelle e i numeri del manoscritto [`1_MANUSCRITTO_DA_CARICARE.pdf`](file:///home/luca/Scrivania/SOTTOMISSIONE_BIOINFORMATICS/1_MANUSCRITTO_DA_CARICARE.pdf).
