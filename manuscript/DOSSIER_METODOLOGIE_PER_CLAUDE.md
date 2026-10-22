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

### Matrici di Correlazione Empirica tra Moduli (Dati Reali):
Contro l'obiezione che $R$ sia indistinguibile dallo stroma o che $w$ sia un mero artefatto:
1. **$S$ e $R$ si auto-escludono nello spazio:** $r = -0.603$ (Lung), $r = -0.517$ (Breast), $r = -0.218$ (Melanoma).
2. **$R$ e Stroma $C$ sono correlati ma distinti:** $r = +0.675$ (Lung, $r^2 = 0.45$), $r = +0.648$ (Breast, $r^2 = 0.42$), $r = +0.454$ (Melanoma, $r^2 = 0.21$). Oltre il **55%–80%** della varianza di $R$ è completamente indipendente dai fibroblasti.
3. **Costo Metabolico della Resistenza:** La proliferazione $w$ è fortemente positiva con l'epitelio sensibile $S$ ($r = +0.66\text{--}0.71$), ma è **negativa con $R$** ($r = -0.15\text{--}-0.54$) e con lo stroma $C$ ($r = -0.31\text{--}-0.44$), confermando lo stato quiescente/costoso dei cloni resistenti in assenza di farmaco.

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
Di conseguenza, i 3 regressori intrinseci $b_i f_{i,s}$ appartengono esattamente allo span lineare dei 9 regressori d'interazione. La matrice non vincolata a 12 parametri ($X_{12}$) è quindi intrinsecamente singolare, con un Condition Number elevato ($\approx 3.6\text{--}4.1 \times 10^8$).

### La Formulazione Rigorosa a 9 Parametri Assorbiti ($\tilde{A}_{ij}$):
Assorbendo la proliferazione di base nella matrice di payoff:
$$\tilde{A}_{ij} = b_i + A_{ij} \implies w_{purified,s} = \sum_{i} \sum_{j} \tilde{A}_{ij} f_{i,s} \bar{f}_{j,s}$$
- **Condition Number di $X_9$:**
  - **Lung SCC:** **$31.4$**
  - **Breast IDC:** **$54.4$**
  - **Melanoma SKCM:** **$45.4$**  
  *(Crollo di 7 ordini di grandezza, piena stabilità e invertibilità numerica).*
- **Varianza Spiegata $R^2$:** identica al decimale ($30.86\%$, $32.29\%$, $36.78\%$).
- **Perché Ridge ($L_2$, $\alpha = 10^{-3}$) senza intercetta:**  
  Poiché $\sum_{i,j} f_{i,s}\bar{f}_{j,s} = 1$, la somma delle 9 colonne genera un vettore costante di $1$. L'aggiunta di un'intercetta libera in un OLS genererebbe una trappola delle variabili dummy (collinearità perfetta con la costante). La regressione Ridge senza intercetta è la formulazione algebricamente corretta per dati composizionali (Aitchison, 1986).

### Dinamica del Replicatore e Risoluzione Lung SCC:
Nella dinamica continua $\dot{x}_i = x_i [F_i(\mathbf{x}) - \bar{F}(\mathbf{x})]$, il campo vettoriale dipende da $\tilde{A}$:
$$\tilde{\Delta}_S = \tilde{A}_{SS} - \tilde{A}_{RS} = \Delta_S + (b_S - b_R)$$
$$\tilde{\Delta}_R = \tilde{A}_{RR} - \tilde{A}_{SR} = \Delta_R - (b_S - b_R)$$

- **Nel Lung SCC:** Nel tessuto non perturbato, la forte proliferazione intrinseca epiteliale ($b_S - b_R = +0.392$) produce $\tilde{\Delta}_S = +0.983$ e $\tilde{\Delta}_R = -0.244$, determinando **Sensitive Dominance** nel tumore vergine.
- **Sotto Terapia ($d_S > d_R$):** Il debulking chemioterapico abbatte selettivamente il tasso di crescita epiteliale. Non appena il danno supera il delta basale, $\tilde{\Delta}_R$ transisce in territorio positivo, facendo emergere la dinamica di **Stag Hunt (Bistabile)** con separatice critica $x_{sep} = 20.0\%$.
- **Breast e Melanoma:** Sia la matrice a 12 parametri che la matrice a 9 parametri assorbiti classificano negli stessi identici regimi (Hawk-Dove per Breast con $x^* = 60.5\%$ vs $61.3\%$; Stag Hunt per Melanoma con $x_{sep} = 0.528$ vs $0.481$).

---

## 7. La Batteria di Validazione Empirica e i Controlli Nulli

### A. Il Controllo Nullo Spaziale (Quadratico Locale Non-Spaziale):
Per verificare se il vicinato Delaunay apporti informazione oltre la co-localizzazione nello stesso spot, abbiamo testato:
- **Modello 2:** Lineare ($f_S, f_R, f_C$)
- **Modello 2b:** Quadratico Locale non-spaziale ($f_i \cdot f_j$ sullo stesso spot)
- **Modello 3:** InSituGame Spaziale ($f_i \cdot \bar{f}_j$ con vicini Delaunay)

| Coorte | M2 (Lineare) | M2b (Quadratico Locale) | **M3 (InSituGame Spaziale)** | Incremento Spaziale Netto | Nested $F$-test vs M2b |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Lung SCC** | 27.84% | 28.26% | **30.86%** | **+2.60%** | **$F = 48.25, p < 1.11 \times 10^{-16}$** |
| **Breast IDC** | 23.49% | 28.38% | **32.29%** | **+3.91%** | **$F = 95.92, p < 1.11 \times 10^{-16}$** |
| **Melanoma** | 32.48% | 33.44% | **36.78%** | **+3.34%** | **$F = 60.74, p < 1.11 \times 10^{-16}$** |

Il vicinato spaziale offre un guadagno significativo di varianza purificata che non può essere replicato dalla co-localizzazione quadratica locale.

### B. Test di Permutazione Monte Carlo (1.000 Shuffle delle Coordinate):
- Lung SCC: $R^2$ reale $= 30.86\%$ vs permutato $= 27.99\%$ ($p = 0.001$)
- Breast IDC: $R^2$ reale $= 32.29\%$ vs permutato $= 23.61\%$ ($p = 0.001$)
- Melanoma: $R^2$ reale $= 36.78\%$ vs permutato $= 32.64\%$ ($p = 0.001$)
Distruggere la topologia spaziale azzera il vantaggio del modello, dimostrando che l'effetto dipende dalla reale geometria cellulare.

### C. Spatial Block Cross-Validation (5 Folds Geograficamente Disgiunti):
- Lung SCC: $R^2_{\text{OOF}} = 26.07\%$
- Breast IDC: $R^2_{\text{OOF}} = 27.02\%$
- Melanoma SKCM: $R^2_{\text{OOF}} = 33.42\%$
Il modello generalizza su blocchi tissutali fisicamente non contigui.

### D. Attenuazione dell'Autocorrelazione Spaziale (Moran's I):
Il modello assorbe gran parte dell'autocorrelazione spaziale tissutale:
- Lung: Moran's I passa da $0.575$ a $0.250$ (**56.48% di attenuazione**)
- Breast: Moran's I passa da $0.626$ a $0.220$ (**64.86% di attenuazione**)
- Melanoma: Moran's I passa da $0.640$ a $0.256$ (**59.90% di attenuazione**)

### E. Repliche Seriali e Test Cieco Out-Of-Sample (OOS):
- Due fette seriali a $\approx 10\,\mu\text{m}$ (Block A): concordanza parametri Pearson $r = +0.811$ ($p = 0.0079$), Spearman $\rho = +0.633$ ($p = 0.067$), pendenza $= 1.36$.
- Stessa classificazione di gioco: entrambe classificate come **Stag Hunt** ($x_{sep} = 0.787$ vs $0.827$, concordanza entro il 4%).
- **Validazione Cieca OOS su Sezione 2:** Il modello fittato su Sezione 1 predice la proliferazione di Sezione 2 con $R^2_{\text{OOS}} = 10.20\%$ (rispetto al $6.21\%$ del modello di massa, guadagno netto $+3.99\%$; ROC-AUC $= 0.653$ vs $0.628$; Fisher exact $p = 4.82 \times 10^{-12}$; $\text{DOR} = 2.30$ [95% CI: $2.03\text{--}2.61$]).
- **Divergenza Biologica Breast IDC vs Block A:** Il campione di Tabella 1 (FFPE CytAssist) è un gioco Hawk-Dove ($x^* = 60.5\%$), mentre il campione seriale Block A (Fresh Frozen) è uno Stag Hunt. Questo dimostra empiricamente che i giochi evolutivi **non sono costanti universali di istotipo, ma architetture eterogenee paziente-specifiche**, motivando la necessità di un test diagnostico su biopsia.

---

## 8. Mappatura degli Script del Repository

Tutti gli script sono documentati, modulari e versionati nel repository GitHub [`https://github.com/Lukecele/spatial-game-assay`](https://github.com/Lukecele/spatial-game-assay):

- `src/insitu_spatial_game.py`: Pipeline completa di streaming, calcolo Delaunay, disaccoppiamento OLS e regressione Ridge.
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
