# Guida Operativa alla Sottomissione su arXiv per lo Studio 1 (InSituGame)

Questa cartella contiene tutti i file e i metadati pronti all'uso per depositare lo **Studio 1 (`InSituGame`)** su **arXiv** nelle sezioni di Biologia Quantitativa e Sistemi Dinamici (**`q-bio.QM`**, **`q-bio.CB`**, **`math.DS`**).

---

## 🎯 Perché arXiv è la Sede Ideale per InSituGame

1. **La casa della Mathematical Oncology:**  
   Tutti i principali pionieri e ricercatori dell'oncologia matematica ed evolutiva (come David Basanta, Jeffrey West, Artem Kaznatcheev, Robert Gatenby, Alexander Anderson dell'Integrated Mathematical Oncology del Moffitt Cancer Center) depositano regolarmente i loro preprint su **arXiv `q-bio.QM`** e **`q-bio.CB`**.
2. **Nessun blocco su "Adaptive Therapy":**  
   Mentre *bioRxiv* ha un filtro rigido che blocca qualsiasi manoscritto che citi terapie adattative o biopsie umane per paura di dare "raccomandazioni cliniche non validate da trial", *arXiv* comprende che si tratta di **metodi quantitativi computazionali e biomatematica pura**.
3. **Indicizzazione e Citabilità Immediata:**  
   Una volta approvato, il paper riceve un identificatore permanente `arXiv:2610.xxxxx`, viene indicizzato automaticamente su Google Scholar, Semantic Scholar, NASA ADS, e può essere citato liberamente dalla comunità scientifica.

---

## 📁 Contenuto della Cartella

| File | Descrizione |
| :--- | :--- |
| **`InSituGame_Manuscript_Full.pdf`** | Il PDF finale completo di 21 pagine con testo, formule matematiche, tabelle e le 3 figure ad alta risoluzione integrate. |
| **`METADATA_ARXIV_COPIA_INCOLLA.txt`** | File di testo puro con tutti i campi (Titolo, Autore, Categorie, Abstract condensato a 1,886 caratteri) pronto per il copia-incolla rapido. |
| **`insitugame_arxiv_tex_source.tar.gz`** | Archivio con sorgente LaTeX e figure per sottomissione avanzata (opzionale se si preferisce TeX rispetto al PDF). |

---

## 🚀 Istruzioni Passo-Passo per la Sottomissione

1. **Accedi a arXiv:**
   - Vai su **[https://arxiv.org/submit](https://arxiv.org/submit)**.
   - Accedi con il tuo account arXiv (o registrati gratuitamente in 30 secondi se non ne hai ancora uno).

2. **Step 1: Verify & Agreement**
   - Conferma i tuoi dati personali e accetta i Termini di Servizio di arXiv.

3. **Step 2: License**
   - Seleziona **`arXiv.org perpetual, non-exclusive license to distribute this article`** (scelta standard consigliata) oppure **`Creative Commons Attribution 4.0 International (CC BY 4.0)`**.

4. **Step 3: Select Primary Category**
   - **Primary Subject:** Seleziona **`Quantitative Biology (q-bio)`** -> **`Quantitative Methods (q-bio.QM)`**.
   - Clicca su *Continue*.

5. **Step 4: Upload File**
   - Seleziona il file **`InSituGame_Manuscript_Full.pdf`** (5.5 MB) da questa cartella e caricalo.
   - *Nota:* Se arXiv chiede se il file è stato prodotto con TeX/LaTeX, seleziona "No / External Tool" per la sottomissione diretta del PDF, oppure carica l'archivio `insitugame_arxiv_tex_source.tar.gz`.

6. **Step 5: Metadata Entry (Copia e Incolla da `METADATA_ARXIV_COPIA_INCOLLA.txt`)**
   - **Title:**  
     `InSituGame: In Situ Spatial Evolutionary Game Assay for Inferring Tumor Subclonal Payoff Matrices and Adaptive Therapy Regimens in Spatial Transcriptomics`
   - **Authors:**  
     `Luca Celebrano`
   - **Abstract:**  
     Incolla il testo dell'abstract condensato (1,886 caratteri, rigorosamente sotto il tetto massimo di 1,920 caratteri imposto dal modulo web di arXiv).
   - **Comments:**  
     `21 pages, 3 figures, 2 tables. Open-source Python package 'insitugame' available at https://github.com/Lukecele/spatial-game-assay`
   - **Report-no / DOI:**  
     `10.5281/zenodo.23585216`
   - **Cross-lists (Categorie Secondarie):**  
     Aggiungi:
     - `q-bio.CB` (Cell Behavior)
     - `math.DS` (Dynamical Systems)

7. **Step 6: Preview & Final Submission**
   - Clicca su **Preview** per verificare l'impaginazione e i metadati.
   - Clicca su **Submit**.

---

## ⏱️ Tempi di Pubblicazione su arXiv

- Le sottomissioni inviate prima delle 14:00 EST (20:00 ora italiana) vengono elaborate quella sera stessa e annunciate la notte successiva (alle 00:00 UTC).
- Riceverai una mail di conferma con l'identificativo arXiv (`arXiv:2610.xxxxx`).
- Il preprint sarà immediatamente visibile online e collegabile a Zenodo e GitHub.
