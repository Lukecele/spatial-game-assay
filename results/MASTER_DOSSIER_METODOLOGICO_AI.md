# MASTER DOSSIER METODOLOGICO, MATEMATICO ED EMPIRICO (v2, corretto)
## Progetto *InSituGame* — Deconvoluzione spaziale della teoria dei giochi evolutivi nel microambiente tumorale

**Autore:** Luca Celebrano (MD Candidate, Università degli Studi di Napoli Federico II)
**Rivista target:** *Bioinformatics* (OUP)
**Repository:** https://github.com/Lukecele/spatial-game-assay
**Zenodo:** https://doi.org/10.5281/zenodo.23585216
**Dataset:** N = 20.093 spot in-tissue (10x Genomics CytAssist FFPE e Visium v1), 4 sezioni indipendenti dal punto di vista del campione biologico (3 tumori + 1 blocco con due sezioni seriali).

> **Stato di questa versione.** Questo dossier è stato riallineato al testo del manoscritto (Tabelle 1-2, Sezioni 3.2-3.6). Ogni numero marcato **[MS]** è presente nel manoscritto. Ogni numero marcato **[DA VERIFICARE]** NON compare nel manoscritto e va ripreso dall'output degli script prima di essere usato. Vedi Sezione 8 (Errata rispetto alla v1).

---

## 1. Istruzioni per le AI che leggono questo dossier

1. Luca Celebrano è MD candidate; il lavoro è computazionale, indipendente, senza finanziamenti di laboratorio bagnato.
2. Obiettivo: stimare matrici di payoff assorbite (Ã_ij) da trascrittomica spaziale, in modo esplorativo (proof-of-concept).
3. **Non affermare** che il lavoro dimostri payoff evolutivi causali: le Ã_ij sono coefficienti di una regressione bilineare su moduli di espressione e possono riflettere cellularità, purezza e rapporto epitelio/stroma (vedi Sezione 7).
4. Pacchetto di sottomissione: `/home/luca/Scrivania/SOTTOMISSIONE_BIOINFORMATICS/` (manoscritto PDF, cover letter PDF, metadati ScholarOne, `FIGURE_SEPARATE/`, questo dossier).
5. Riproducibilità:
   ```bash
   cd /home/luca/spatial-game-assay && python3 -m unittest tests/test_insitugame_reproducibility.py
   ```
   (11 test, circa 28-30 s.) I test verificano la correttezza numerica del codice, non la validità biologica dell'inferenza.

---

## 2. Architettura del repository

```
spatial-game-assay/
├── data/ (breast, melanoma, lung, breast_visium_sec1, breast_visium_sec2, cache)
├── src/
│   ├── run_pan_cancer_assay.py                  # pipeline primaria (OLS identificabile)
│   ├── insitu_spatial_game.py                   # core: Delaunay + deconvoluzione
│   ├── clinical_diagnostic_metrics.py           # 2x2, OR, ROC-AUC, Fisher
│   ├── master_comprehensive_audit.py            # null e permutazioni
│   ├── exhaustive_bayesian_mcmc_melanoma.py     # ancillare
│   ├── exhaustive_multiscale_lung_scan.py       # esplorativo
│   ├── exhaustive_spatial_autoregressive_sar.py # controllo SAR/SEM
│   ├── exhaustive_model_selection_parsimony.py  # AIC/BIC
│   ├── exhaustive_pontryagin_optimal_control.py # toy-model MTD vs adattiva
│   └── literature_and_ergodicity_benchmark.py
├── manuscript/ (md, pdf, cover letter)
└── tests/test_insitugame_reproducibility.py
```

---

## 3. Fondamenti biologici e firme molecolari

Normalizzazione: E_{g,s} = ln(1 + C_{g,s}/ΣC_{k,s} × 10⁴). Punteggio di modulo = media non pesata dei geni.

**Breast IDC:** T_S = EPCAM, CDH1, KRT8, KRT18, ESR1; T_R = VIM, CDH2, ZEB1, AXL, MMP2; CAF = ACTA2, COL1A1, FAP, PDGFRB.
**Melanoma SKCM:** T_S = MLANA, MITF, PMEL, TYR; T_R = AXL, VIM, WNT5A, EGFR; CAF = ACTA2, COL1A1, FAP, PDGFRB.
**Lung SCC:** T_S = EGFR, KRT5, TP63, EPCAM; T_R = VIM, AXL, CDH2, MMP2; CAF = ACTA2, COL1A1, FAP, PDGFRB.
**Comuni:** proliferazione w = MKI67, PCNA, TOP2A, CDK1; ipossia Z_hyp = CA9, VEGFA, HIF1A; vascolare Z_vasc = PECAM1, VWF.

> Nota: il pannello CAF nel manoscritto include PDGFRB per tutte le coorti; nella v1 del dossier era omesso per melanoma e polmone. Verificare nello script quale versione è effettivamente usata.

**Correlazioni tra moduli grezzi [MS]:** T_R vs CAF r = 0,675 (Lung), 0,648 (Breast), 0,454 (Melanoma). w vs T_S r = +0,709 / +0,662 / +0,683; w vs T_R r = −0,537 / −0,366 / −0,154; w vs stroma r = −0,314 / −0,444 / −0,434. Queste correlazioni sono la ragione principale per cui le Ã vanno interpretate con cautela (Sezione 7).

---

## 4. Formulazione matematica

### 4.1 Simplesso
f_{i,s} = softmax((M_{i,s} − max_k M_{k,s})/T), con T = 1, i ∈ {S, R, C}, Σ_i f_{i,s} = 1.

### 4.2 Delaunay
f̄_{j,s} = media di f_{j,k} sui vicini Delaunay k ∈ N(s). Distanza tipica tra centroidi ≈ 100 µm (spot da 55 µm): i termini f_i f̄_j descrivono segnalazione regionale/paracrina, non contatto diretto.

### 4.3 Decoupling microambientale
w_raw = γ0 + γ_hyp Z_hyp + γ_vasc Z_vasc + w_purified. Il residuo è ortogonale a Z_hyp e Z_vasc **per costruzione dell'OLS**. Attenzione: il teorema FWL garantisce l'equivalenza con la regressione congiunta solo se i regressori del gioco sono anch'essi residualizzati su Z. Il manoscritto riporta anche il modello congiunto su w_raw (Tabella 2) come controllo.

### 4.4 Identificabilità
Σ_j f_i f̄_j = f_i, quindi i termini b_i f_i stanno nello span dei 9 termini bilineari: il modello a 12 parametri è rango-deficiente. Si usa la matrice assorbita Ã_ij = b_i + A_ij (9 parametri) con OLS senza intercetta. Numeri di condizione [MS]: 31,4 (Lung), 54,4 (Breast), 45,4 (Melanoma).

**Conteggio dei parametri (corretto).** Le 9 colonne di X9 sono linearmente indipendenti (somma = vettore di uni): M3 ha rango **9**, non 8. M2c (f_S, f_R, f̄_S, f̄_R + costante, che è nello span di X9) ha **5** parametri. Il test F annidato ha quindi df1 = 9 − 5 = 4 e df2 = N − 9: il df1 del manoscritto è corretto, ma le etichette "rank 8" (M3) e "rank 4" (M2c) vanno sostituite con 9 e 5. M2d ha 6 + 6 termini quadratici con una dipendenza condivisa: rango **11**, non 10.

### 4.5 Dinamica del replicatore (sotto-gioco S vs R, senza stroma)
Δ̃_S = Ã_SS − Ã_RS, Δ̃_R = Ã_RR − Ã_SR.
ẋ = x(1−x)[x(Δ̃_S + Δ̃_R) − Δ̃_R], x = f_S/(f_S+f_R).
Punto fisso interno: x* = Δ̃_R/(Δ̃_S + Δ̃_R).
- Δ̃_S < 0 e Δ̃_R < 0: Hawk-Dove, x* stabile.
- Δ̃_S > 0 e Δ̃_R > 0: Stag Hunt, x* = separatrice instabile.

### 4.6 Correzione per lo stroma (NUOVO — non presente nel manoscritto)
Con f̄_S = (1−c) y, f̄_R = (1−c)(1−y), f̄_C = c, dove c è la frazione stromale del vicinato e y la frazione di sensibili tra i cloni tumorali:

w_S − w_R = (1−c)[ y(Δ̃_S + Δ̃_R) − Δ̃_R ] + κ c, con κ = Ã_SC − Ã_RC.

Il termine κc **non** è un semplice offset irrilevante: sposta il punto fisso a

y* = [ Δ̃_R − κ c/(1−c) ] / (Δ̃_S + Δ̃_R).

Valori illustrativi (c = 0,3; Ã_SC e Ã_RC letti dalla Figura 2, non riportati nel testo del manoscritto):
- Breast: κ = −1,09 − 0,37 = −1,46 → y* ≈ **53%** (vs 61,3% senza correzione).
- Melanoma: κ = +0,50 − 0,03 = +0,47 → y_sep ≈ **42%** (vs 48,1%).

Il valore di c va sostituito con la frazione stromale media osservata ⟨f̄_C⟩ per coorte, oppure si analizza il gioco a 3 strategie. Il regime qualitativo (Hawk-Dove / Stag Hunt) non cambia, ma le soglie numeriche sì, e quindi anche le soglie di dosaggio adattivo (es. 0,387).

### 4.7 Proposizione 1 (space-for-time)
Va presentata come **assunzione/proposizione condizionale**, non come teorema dimostrato: richiede che (i) l'espressione di cicline sia proporzionale al tasso di divisione per-capita, (ii) i vicinati campionino densamente il simplesso, (iii) non ci siano confondenti non modellati. Il testo "uniquely reconstructed" non è supportato da una dimostrazione.

---

## 5. Risultati del manoscritto [MS]

| Coorte | Spot | Ã (SS, SR, RS, RR) | Regime | R² M3 | R² OOF M3 vs M2c | Test F vs M2c | Moran I (M3) |
|---|---|---|---|---|---|---|---|
| Breast IDC | 4.992 | −0,74, +2,87, +2,19, −1,77 | Hawk-Dove, x* = 61,3% [0,526–0,725] | 32,29% | 27,02% vs 19,76% (+7,26%) | F = 104,99, p < 1,11e−16 | 0,220 |
| Melanoma SKCM | 3.458 | +0,93, −1,78, −0,70, −0,27 | Suggestive Stag Hunt, x_sep = 48,1% [0,168–0,637] | 36,78% | 33,42% vs 32,88% (+0,54%) | F = 17,16, p = 6,04e−14 | 0,256 |
| Lung SCC | 3.858 | +0,68, +0,05, −0,31, −0,19 | Indeterminato / regionale | 30,86% | 26,07% vs 27,06% | F = 1,30, p = 0,266 | 0,250 |

Selettivi: Breast Δ̃_S = −2,924 [−3,93; −1,55], Δ̃_R = −4,638 [−6,37; −2,58]. Melanoma Δ̃_S = +1,633 [+0,47; +2,60], Δ̃_R = +1,515 [−0,97; +3,57]. Lung Δ̃_S = +0,983, Δ̃_R = −0,244 [−1,685; +1,168].
Confronto con M2d (generico quadratico) OOF: Breast 28,46% (**superiore** a M3, 27,02%), Melanoma 32,99%, Lung 26,40%. In Breast quindi M3 non batte il controllo non-game M2d: il vantaggio di M3 è parsimonia e interpretabilità, non potere predittivo.

**Sezioni seriali Block A (Visium fresh-frozen) [MS]:**
- Concordanza parametri 3×3: Pearson r = +0,811 (p = 0,0079), Spearman ρ = +0,633 (p = 0,067), MAE = 2,35, pendenza 1,36. I 9 parametri non sono indipendenti e il risultato è in parte trainato da A_RR.
- Entrambe le sezioni: Stag Hunt, x_sep = **0,787 vs 0,827** (ordine sezione 1 vs 2 da confermare).
- OOS (training sez. 1, test sez. 2): R² = 10,20% vs 6,21% (baseline mass-action); AUC 0,653 vs 0,628; DOR 2,30 [2,03–2,61] vs 2,17 [1,91–2,46] (intervalli sovrapposti); Fisher p = 4,82e−12.
- I valori per-sezione di Ã, R² e Moran I per Block A: **[DA VERIFICARE]** dall'output degli script (la v1 del dossier riportava matrici che non coincidono con x_sep ≈ 0,80).

**Benchmark in silico [MS]:** N = 1.000 spot sintetici, 20% rumore, r = 0,9997, MAE = 0,023, x* vero 60,53% vs stimato 60,30%. È un test **model-matched**: verifica la correttezza algebrica, non la robustezza a modello mal specificato.

**Metriche binarizzate [MS]:** DOR = 7,04 [6,11–8,11] (Lung), 6,16 [5,45–6,96] (Breast), 16,33 [13,81–19,30] (Melanoma). La soglia è la mediana, quindi sensibilità e specificità sono quasi uguali per costruzione; il DOR del Lung è spiegato dalla composizione regionale (M2c), non dai termini di gioco. I valori esatti di sensibilità/specificità citati nella v1 **[DA VERIFICARE]**.

**Simulazione terapeutica [MS, illustrativa]:** MTD con d_S = 1,2, d_R = 0: in Breast x*_R(D) = 54,5%; in Melanoma x_sep(D) = (1,515 + 1,2)/3,148 ≈ 86,2%. Parametri citotossici euristici, non calibrati.

---

## 6. Analisi primarie vs esplorative

**Primarie (nel manoscritto):** OLS senza intercetta su 9 parametri; 5-fold spatial block CV (k-means); 200 repliche di block bootstrap; Delaunay a priori; Moran I residuo (Cliff & Ord); permutazione coordinate (1.000).
**Esplorative (solo repository):** MCMC melanoma (P(Stag Hunt) = 99,87%, x_sep = 35,4%, con prior informativa); scansione multiscala Lung r = 250-300 µm (F = 49,01, rischio di p-hacking, non riportata); trial virtuale Pontryagin (PFS 86 vs 75 giorni, HR 0,867: toy-model, non evidenza clinica); SAR/SEM.

---

## 7. Limiti noti e interpretazione prudente

1. **Numerosità effettiva.** Gli spot sono pseudo-replicati. I campioni biologici indipendenti sono 4 (1 per istotipo + 1 blocco Block A con due sezioni tecniche). "Pan-cancer", "lineage-specific" e "patient-specific" sono ipotesi, non risultati.
2. **Confondimento piattaforma.** Breast FFPE/probe vs Block A fresh-frozen/poly-A differiscono per donatore, istotipo, conservazione e chimica: la divergenza di regime non è attribuibile alla biologia.
3. **Interpretazione delle Ã.** T_R correla con lo stroma e w con T_S; le Ã possono catturare cellularità, purezza tumorale e mix epitelio/stroma. Controlli mancanti: null con moduli genici casuali, deconvoluzione (cell2location/RCTD), purezza come covariata, sensibilità a T del softmax e ai pannelli genici (leave-one-gene-out).
4. **Inferenza.** Moran I ≈ 0,22-0,26: i p-value parametrici (es. 1e−16) sono anti-conservativi. Le prove principali sono block-CV e block bootstrap; indicare il numero di blocchi e se DOR/AUC/Fisher sono out-of-fold.
5. **Stroma.** Vedi 4.6: la proiezione 2×2 va corretta o sostituita con il gioco a 3 strategie.
6. **Proposizione 1** è un'assunzione (4.7).
7. **Terapia.** I parametri d_S, d_R sono euristici; nessuna evidenza clinica.

---

## 8. Errata rispetto alla v1 del dossier

| Voce | v1 (errata) | Corretto |
|---|---|---|
| Lung, Ã (SS, SR, RS, RR) | +0,68, **−0,31**, −0,31, −0,19 | +0,68, **+0,05**, −0,31, −0,19 |
| Moran I M3 Breast / Melanoma / Lung | 0,237 / 0,228 / 0,257 | 0,220 / 0,256 / 0,250 |
| Block A x_sep | 81,1% e 79,4% | 0,787 e 0,827 |
| Block A matrici 2×2 della v1 | +0,42, −1,21, −0,85, +0,33 ecc. | non presenti nel manoscritto; implicano x_sep ≈ 55%, incoerente con ≈ 0,80: **[DA VERIFICARE]** |
| Rango M3 / M2c / M2d | 8 / 4 / 10 | 9 / 5 / 11 (df1 = 4 invariato) |
| "Allineato al 100% col manoscritto" | presente | rimosso: la v1 non lo era |
| Soglia adattiva e x*, x_sep | senza correzione stroma | aggiunta Sezione 4.6 |
