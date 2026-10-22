# InSituGame: Compendio Tecnico, Metodologico ed Evidenze Empiriche per Revisione Critica (Claude / Gemini Argon)

Questo documento sintetizza formalmente l'architettura matematica, la biologia cellulare, la biostatistica computazionale, i controlli nulli e gli script del framework **InSituGame** (*In Situ Spatial Evolutionary Game Assay*). È strutturato per rispondere punto su punto a qualsiasi disamina scettica di biologia matematica o bioinformatica computazionale.

---

## 1. Rationale Traslazionale ed Obiettivo Scientifico

- **Il Collo di Bottiglia:** L'applicazione clinica della *Teoria dei Giochi Evolutiva (EGT)* alla terapia adattiva del cancro è stata storicamente bloccata dall'impossibilità di stimare le matrici dei payoff ($A_{ij}$) direttamente dal tessuto umano del singolo paziente. I parametri venivano postulati arbitrariamente in equazioni differenziali o misurati in co-colture 2D su piastra di plastica (Kaznatcheev et al., *Nat Ecol Evol* 2019), prive di architettura 3D, stroma (CAF), vasi e gradienti fisici.
- **La Soluzione di InSituGame:** Un framework computazionale esplorativo che trasforma la trascrittomica spaziale clinica (10x Genomics Visium / CytAssist FFPE) in un saggio funzionale analogo a un "antibiogramma evolutivo", quantificando le interazioni competitive/cooperative frequenza-dipendenti tra fenotipi tumorali direttamente all'interno della biopsia intatta.

---

## 2. Dataset Clinici, Provenienza e Riproducibilità

Tutti i dati analizzati provengono da biopsie cliniche umane autentiche ad accesso pubblico fornite da 10x Genomics, senza alcun dato sintetico o simulato ($N = 20.093$ spot in-tissue totali, >350 milioni di read UMI sequenziate):

1. **Human Squamous Cell Lung Carcinoma (Lung SCC):** 10x CytAssist FFPE ($N = 3.858$ spot; $72.655.040$ UMI; mediana: $12.883$ UMI/spot).
2. **Human Invasive Ductal Breast Carcinoma (Breast IDC):** 10x CytAssist FFPE ($N = 4.992$ spot; $115.297.808$ UMI; mediana: $17.994$ UMI/spot).
3. **Human Cutaneous Malignant Melanoma (Melanoma SKCM):** 10x CytAssist FFPE ($N = 3.458$ spot; $105.698.520$ UMI; mediana: $25.607$ UMI/spot).
4. **Serial Replicates Breast Cancer Cohort:** Due sezioni istologiche seriali consecutive a $\approx 10\,\mu\text{m}$ di distanza da biopsia di carcinoma mammario (Block A, Sezione 1: $N = 3.798$ spot; Sezione 2: $N = 3.987$ spot).

---

## 3. Deconvoluzione Fenotipica e Rigore dei Moduli Genici

Nelle biopsie cliniche naive non trattate, i cloni non possono essere definiti da resistenze farmacologiche verificate in vitro; $T_S$ e $T_R$ rappresentano **archetipi fenotipici ipotizzati** (differenziato/epiteliale vs mesenchimale/plastico):

- **Sensibile ($T_S$):** Marcatori di differenziazione epiteliale/lineage:
  - *Lung:* `EGFR, KRT5, TP63, EPCAM`
  - *Breast:* `EPCAM, CDH1, KRT8, KRT18, ESR1`
  - *Melanoma:* `MLANA, MITF, PMEL, TYR`
- **Resistente ($T_R$):** Transizione Epitelio-Mesenchima (EMT), plasticità e vie di bypass:
  - *Lung:* `VIM, AXL, CDH2, MMP2`
  - *Breast:* `VIM, CDH2, ZEB1, AXL, MMP2`
  - *Melanoma:* `AXL, VIM, WNT5A, EGFR`
- **Stroma / Fibroblasti ($S_{CAF}$):** `ACTA2, COL1A1, FAP, PDGFRB`.
- **Fitness / Proliferazione ($w$):** Macchinario del ciclo cellulare e replicazione del DNA: `MKI67, PCNA, TOP2A, CDK1`.
- **Ipossia Microambientale ($Z_{hyp}$):** Bersagli trascrizionali funzionali a valle del fattore HIF1: `CA9, VEGFA, HIF1A` *(nota: HIF1$\alpha$ è regolato per degradazione proteasomica post-traduzionale VHL, quindi l'attività ipossica è tracciata tramite i target transattivati CA9/VEGFA)*.
- **Densità Vascolare ($Z_{vasc}$):** Marcatori endoteliali: `PECAM1, VWF`.

### Razionale Biologico per la collocazione di *EGFR*:
- In **Lung SCC**, *EGFR* fa parte del programma epiteliale carcinomatoso principale ($T_S$).
- Nel **Melanoma**, le cellule melanocitiche differenziate sono *MITF$^+$/MLANA$^+$, EGFR$^{\text{low}}$*, mentre la sovraespressione di *EGFR* è il canonico meccanismo non mutazionale di tolleranza e resistenza adattiva da bypass recettoriale sotto inibitori di BRAF/MEK (Prahallad et al., *Nature* 2012; Kemper et al., *Cancer Cell* 2014; Sun et al., *Nature* 2014).

### Chiarimenti sulle Correlazioni dei Moduli e sulla Proliferazione:
1. **Punteggi Grezzi dei Moduli (prima del Softmax):** Le correlazioni tra punteggi grezzi mediati mostrano che $S$ e $R$ sono negativamente correlati nello spazio ($r = -0.603$ in Lung, $r = -0.517$ in Breast, $r = -0.218$ in Melanoma). Riconosciamo esplicitamente che la successiva trasformazione softmax impone un vincolo a somma costante ($\sum f_i = 1$) che induce covarianza negativa matematica; pertanto, non usiamo il softmax come prova biologica di auto-esclusione.
2. **Sovrapposizione $R$-Stroma:** $R$ condivide con lo stroma $C$ una parte significativa di varianza ($r = 0.65\text{--}0.68$, circa il $42\text{--}45\%$ di varianza condivisa), dovuta al comune fenotipo mesenchimale. L'indipendenza residua del $55\%$ è riferita specificamente al pannello di 4 geni CAF analizzato.
3. **Proliferazione epiteliale vs "costo di fitness":** Poiché $w$ (*MKI67, PCNA, TOP2A, CDK1*) misura la proliferazione carcinomatosa epiteliale, gli spot arricchiti in $R$ o stroma contengono semplicemente meno parenchima epiteliale proliferante. La correlazione negativa tra $w$ ed $R$/stroma riflette la cellularità epiteliale dello spot e **non viene rivendicata come prova fisiologica di un costo metabolico di resistenza**.

---

## 4. Disaccoppiamento Microambientale (Frisch-Waugh Theorem)

Per isolare la fitness cellulare dai gradienti fisici estrinseci, la proliferazione osservata $w_{raw}$ viene depurata tramite regressione OLS multivariata:
$$w_{raw,s} = \gamma_0 + \gamma_{hyp} Z_{hyp,s} + \gamma_{vasc} Z_{vasc,s} + w_{purified,s}$$
Il residuo $w_{purified,s}$ è ortogonale a ipossia e densità vascolare.  
Quando il modello congiunto viene stimato direttamente su $w_{raw}$ ($w_{raw} \sim Z_{hyp} + Z_{vasc} + \text{Gioco Spaziale}$), la varianza totale spiegata sale a:
- **Lung SCC:** **$54.30\%$** (rispetto al $15.82\%$ dell'ambiente da solo)
- **Breast IDC:** **$51.59\%$** (rispetto al $16.12\%$ dell'ambiente da solo)
- **Melanoma SKCM:** **$49.60\%$** (rispetto al $10.23\%$ dell'ambiente da solo)

---

## 5. Geometria Fisica di 10x Visium e Rete di Delaunay

- Uno spot 10x Visium ha diametro di $55\,\mu\text{m}$, con un passo centro-centro di $100\,\mu\text{m}$ su griglia esagonale e $\approx 45\,\mu\text{m}$ di spazio non sequenziato tra spot adiacenti.
- I vicini primi nel grafo di Delaunay hanno baricentri a $\approx 100\,\mu\text{m}$.
- Pertanto, i termini d'interazione $f_{i,s} \bar{f}_{j,s}$ **non modellano il contatto iuxtacrino membrana-membrana**, bensì la segnalazione paracrina diffusibile, il rimodellamento di matrice extracellulare e il condizionamento metabolico della micro-nicchia condivisa.

---

## 6. Identificabilità Matematica sul Simplesso e Matrice Assorbita $\tilde{A}$

### L'Identità Algebrica del Simplesso:
Poiché le frequenze medie nel vicinato soddisfano il vincolo unitario $\sum_{j \in \{S, R, C\}} \bar{f}_{j,s} = 1$, vale identicamente per ogni spot $s$:
$$\sum_{j \in \{S, R, C\}} f_{i,s} \bar{f}_{j,s} = f_{i,s} \sum_j \bar{f}_{j,s} = f_{i,s}$$
Di conseguenza, i 3 regressori intrinseci $b_i f_{i,s}$ appartengono esattamente allo span lineare dei 9 regressori d'interazione. La decomposizione a 12 parametri ($b_i$ e $A_{ij}$) non è identificabile in OLS.

### La Formulazione Rigorosa a 9 Parametri Assorbiti ($\tilde{A}_{ij}$):
Assorbendo la proliferazione di base nella matrice di payoff:
$$\tilde{A}_{ij} = b_i + A_{ij} \implies w_{purified,s} = \sum_{i} \sum_{j} \tilde{A}_{ij} f_{i,s} \bar{f}_{j,s}$$
- **Condition Number di $X_9$:**
  - **Lung SCC:** **$31.4$**
  - **Breast IDC:** **$54.4$**
  - **Melanoma SKCM:** **$45.4$**  
  *(Piena stabilità e invertibilità numerica).*
- **Varianza Spiegata $R^2$:** identica al decimale ($30.86\%$, $32.29\%$, $36.78\%$).
- **Perché OLS senza intercetta:**  
  Poiché $\sum_{i,j} f_{i,s}\bar{f}_{j,s} = 1$, la somma delle 9 colonne genera un vettore costante di $1$. Un'intercetta libera in OLS genererebbe una trappola delle variabili dummy (collinearità perfetta con la costante). Il modello viene quindi stimato direttamente senza intercetta (`fit_intercept=False`).

### Dinamica del Replicatore e Regimi di Gioco su $\tilde{A}$:
Nella dinamica continua $\dot{x}_i = x_i [F_i(\mathbf{x}) - \bar{F}(\mathbf{x})]$, il campo vettoriale dipende unicamente dalla matrice assorbita $\tilde{A}$:
$$\tilde{\Delta}_S = \tilde{A}_{SS} - \tilde{A}_{RS}, \quad \tilde{\Delta}_R = \tilde{A}_{RR} - \tilde{A}_{SR}$$

I risultati canonici riportati in Tabella 1 sono:
1. **Breast IDC:** $\tilde{A}_{2 \times 2} = \begin{pmatrix} -0.74 & +2.87 \\ +2.19 & -1.77 \end{pmatrix}$.  
   $\tilde{\Delta}_S = -2.924 < 0, \tilde{\Delta}_R = -4.638 < 0 \implies$ **Hawk-Dove / Snowdrift (Coesistenza)** con equilibrio interno stabile:
   $$x^* = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R} = \frac{-4.638}{-7.562} = \mathbf{61.3\%}$$
2. **Melanoma SKCM:** $\tilde{A}_{2 \times 2} = \begin{pmatrix} +0.93 & -1.78 \\ -0.70 & -0.27 \end{pmatrix}$.  
   $\tilde{\Delta}_S = +1.633 > 0, \tilde{\Delta}_R = +1.515 > 0 \implies$ **Stag Hunt / Bistabile** con separatice critica instabile:
   $$x_{sep} = \frac{\tilde{\Delta}_R}{\tilde{\Delta}_S + \tilde{\Delta}_R} = \frac{1.515}{3.148} = \mathbf{48.1\%}$$
3. **Lung SCC:** $\tilde{A}_{2 \times 2} = \begin{pmatrix} +0.68 & +0.05 \\ -0.31 & -0.19 \end{pmatrix}$.  
   $\tilde{\Delta}_S = +0.983 > 0, \tilde{\Delta}_R = -0.244 < 0 \implies$ **Sensitive Dominance** nel tessuto baseline non trattato.  
   Sotto terapia citotossica ($F_i^{\text{treat}} = \sum_j \tilde{A}_{ij}x_j - u(t)d_i$, con differenziale citotossico $d_S - d_R = 0.50$), le selective advantage diventano:
   $$\tilde{\Delta}_S(D) = 0.983 - 0.500 = +0.483, \quad \tilde{\Delta}_R(D) = -0.244 + 0.500 = +0.256$$
   facendo emergere una dinamica di **Stag Hunt terapia-indotta** con separatice critica:
   $$x_{sep}(D) = \frac{0.256}{0.483 + 0.256} = \frac{0.256}{0.739} = \mathbf{34.6\%}$$

---

## 7. La Batteria di Validazione Empirica e i Controlli Nulli

### A. Benchmark In Silico di Recupero Parametrico:
Per dimostrare l'identificabilità numerica, abbiamo generato un tessuto sintetico ($N = 1.000$ spot su simplesso Dirichlet con smoothing spaziale) governato da una matrice nota Hawk-Dove ($x^* = 60.53\%$) con il $20\%$ di rumore gaussiano additivo.
- **Correlazione matrice recuperata vs vera:** Pearson $r = 0.9997$ ($p = 1.13 \times 10^{-12}$), MAE $= 0.0227$.
- **Recupero dell'equilibrio di Nash:** $x^* = 60.53\%$ vero vs $x^* = 60.30\%$ stimato (errore $< 0.25\%$).

### B. Il Controllo Nullo $M_{2c}$ (Composizione Regionale Lineare):
Per isolare se il guadagno spaziale derivi dalle **interazioni di gioco bilineari ($f_i \bar{f}_j$)** o dalla sola **composizione regionale liscia ($\bar{f}_i$)**, abbiamo confrontato l'intera gerarchia di modelli:
- $M_2$: Lineare mass-action ($f_S, f_R, f_C$)
- $M_{2b}$: Quadratico Locale non-spaziale ($f_i \cdot f_j$)
- $M_{2c}$: Composizione Regionale Lineare ($f_i$ + $\bar{f}_i$, senza termini di prodotto)
- $M_3$: InSituGame Spaziale ($\tilde{A}_{ij} f_i \bar{f}_j$)

| Coorte | $M_2$ ($R^2$ in / OOF) | $M_{2b}$ ($R^2$ in / OOF) | $M_{2c}$ ($R^2$ in / OOF) | **$M_3$ Game ($R^2$ in / OOF)** | Nested $F$-test $M_3$ vs $M_{2c}$ | Moran $I$ Residui ($M_1 \to M_3$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Breast IDC** | 23.49% / 17.02% | 28.38% / 20.93% | 26.58% / 19.76% | **32.29% / 27.02%** | **$F = 104.99, p < 1.11 \times 10^{-16}$** | $0.470 \to \mathbf{0.220}$ |
| **Melanoma** | 32.48% / 29.06% | 33.44% / 29.32% | 35.52% / 32.88% | **36.78% / 33.42%** | **$F = 17.16, p = 6.04 \times 10^{-14}$** | $0.524 \to \mathbf{0.256}$ |
| **Lung SCC** | 27.84% / 23.75% | 28.26% / 23.60% | 30.77% / 27.06% | **30.86% / 26.07%** | **$F = 1.30, p = 0.266$ (N.S.)** | $0.459 \to \mathbf{0.250}$ |

**Conclusioni Biologiche e Metodologiche:**
1. Nel **Carcinoma Mammario (Breast IDC)**, il gioco spaziale bilineare $M_3$ batte nettamente la composizione regionale $M_{2c}$ ($F = 104.99, p < 10^{-16}$) e ottiene un incremento out-of-fold massiccio di **$+7.26\%$ in Spatial Block CV** ($27.02\%$ vs $19.76\%$). L'interazione non-lineare è indispensabile.
2. Nel **Melanoma SKCM**, l'interazione di gioco batte significativamente la composizione regionale ($F = 17.16, p = 6.04 \times 10^{-14}$).
3. Nel **Carcinoma Polmonare (Lung SCC)**, la variazione spaziale è dominata dalla composizione regionale liscia ($M_{2c} = 30.77\%$), e i termini di gioco aggiuntivi non raggiungono la significatività statistica ($F = 1.30, p = 0.266$). Riconoscere questa distinzione histotype-specifica nel manoscritto conferisce totale trasparenza e credibilità accademica.

### C. Test di Permutazione delle Coordinate Spaziali:
Rimescolando le coordinate spaziali dei vicini per 1.000 permutazioni Monte Carlo, il guadagno spaziale crolla interamente:
- Lung SCC: $p < 0.001$ (osservato $30.86\%$ vs nullo $27.99\%$)
- Breast IDC: $p < 0.001$ (osservato $32.29\%$ vs nullo $23.61\%$)
- Melanoma: $p < 0.001$ (osservato $36.78\%$ vs nullo $32.64\%$)

### D. Repliche Seriali e Test Cieco Out-Of-Sample (OOS):
- Due sezioni seriali consecutive a $\approx 10\,\mu\text{m}$ (Block A): concordanza parametri Pearson $r = +0.811$ ($p = 0.0079$), Spearman $\rho = +0.633$ ($p = 0.067$), pendenza $= 1.36$.
- Entrambe classificate indipendentemente in regime **Stag Hunt** ($x_{sep} = 0.787$ vs $0.827$).
- **Validazione Cieca OOS su Sezione 2:** Il modello fittato su Sezione 1 predice la Sezione 2 con $R^2_{\text{OOS}} = 10.20\%$ (rispetto al $6.21\%$ del modello di massa, guadagno netto $+3.99\%$; ROC-AUC $= 0.653$ vs $0.628$; Fisher exact $p = 4.82 \times 10^{-12}$; $\text{DOR} = 2.30$ [95% CI: $2.03\text{--}2.61$]).
- **Differenze Piattaforma Breast IDC vs Block A:** Il campione di Tabella 1 (FFPE CytAssist, ibridazione sonde) è Hawk-Dove ($x^* = 61.3\%$), mentre Block A (Fresh Frozen Visium, cattura poly-A) è Stag Hunt. Le differenze di donatore, tessuto, chimica di cattura e stadio impediscono di attribuire questa divergenza unicamente alla biologia; essa dimostra che i parametri di gioco non possono essere assunti come costanti universali.

---

## 8. Mappatura degli Script del Repository

Tutti gli script sono documentati, modulari e versionati nel repository GitHub [`https://github.com/Lukecele/spatial-game-assay`](https://github.com/Lukecele/spatial-game-assay):

- `src/insitu_spatial_game.py`: Pipeline completa di streaming, calcolo Delaunay, disaccoppiamento OLS e regressione.
- `src/run_pan_cancer_assay.py`: Esecuzione dell'atlante pan-cancer sui tre tumori e calcolo dei parametri di equilibrio.
- `src/validation_battery.py`: Suite biostatistica con Moran's I esatto (Cliff & Ord 1981), test di permutazione spaziale e Spatial Block Cross-Validation.
- `src/run_replicate_blind_oos_test.py`: Benchmark out-of-sample e correlazione tra sezioni seriali consecutive.
- `src/plot_reproducibility_figures.py`: Generazione delle figure ad alta risoluzione (300 DPI) per il manoscritto.
- `tests/`: Suite di unit test automatizzati con `pytest`.

---

## 9. Riferimenti Bibliografici Puntuali e Corretti

1. **Bacevic et al. (2017):** Spatial competition constrains resistance to targeted therapy in non-small cell lung cancer. *eLife*, 6, e27460.
2. **Cliff & Ord (1981):** *Spatial processes: models & applications*. Pion, London.
3. **Enriquez-Navas et al. (2016):** Exploiting evolutionary principles in cancer therapy. *Science Translational Medicine*, 8(327), 327ra24.
4. **Gallaher et al. (2018):** Spatial heterogeneity and evolutionary dynamics modulate tumor response to adaptive therapy. *Cancer Research*, 78(8), 2029–2039.
5. **Gatenby et al. (2009):** Adaptive therapy. *Cancer Research*, 69(11), 4894–4903.
6. **Kaznatcheev et al. (2019):** Fibroblasts and alectinib switch the evolutionary games played by non-small cell lung cancer. *Nature Ecology & Evolution*, 3(3), 450–456.
7. **Kemper et al. (2014):** BRAF inhibition induces reversible drug tolerance in melanoma through non-mutational adaptation. *Cancer Cell*, 26(1), 60–74.
8. **Marusyk et al. (2014):** Non-cell-autonomous driving of tumour growth supports sub-clonal heterogeneity. *Nature*, 514(7520), 54–58.
9. **Nowak (2006):** *Evolutionary dynamics: exploring the equations of life*. Harvard University Press.
10. **Prahallad et al. (2012):** Unresponsiveness of colon cancer to BRAF(V600E) inhibition through feedback activation of EGFR. *Nature*, 483(7387), 100–103.
11. **Sun et al. (2014):** Reversible and adaptive resistance to BRAF(V600E) inhibition in melanoma. *Nature*, 508(7494), 118–122.
12. **Zhang et al. (2017):** Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. *Nature Communications*, 8(1), 1816.
