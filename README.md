# Explainable ML for Heart Failure Gene Signatures

BSc Data Science dissertation (University of Sheffield, 2026), supervised by Dr Zeyneb Kurt. A related paper was accepted at UKCI 2026.

The project classifies **failing vs non-failing human left-ventricle tissue** from microarray gene expression. It derives a compact **25-gene signature** and checks it biologically with pathway enrichment. The model is trained on one cohort and tested on an independent cohort run on a different array platform.

| Cohort | Role | Platform | Samples | Non-failing | Ischemic CM | Dilated CM |
|---|---|---|---|---|---|---|
| [GSE57338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57338) | training | Affymetrix HuGene 1.1 ST | 313 | 136 | 95 | 82 |
| [GSE5406](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5406) | external test | Affymetrix HG-U133A | 210 | 16 | 108 | 86 |

Ischemic and idiopathic dilated cardiomyopathy are pooled into one **heart failure** class (label 1), and non-failing donor hearts are label 0.

---

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 0 | `R/00_install_packages.R` | Installs the CRAN and Bioconductor packages |
| 0 | `python/00_extract_cel_files.py` | Unpacks `GSE57338_RAW.tar` into `.CEL` files |
| 0 | `R/00_download_metadata_GSE57338.R` | Sample labels for GSE57338 (a copy is already in `data/`) |
| 1 | `R/01_rma_normalise_GSE57338.R` | `oligo` RMA (background correction, quantile normalisation, log2), probes → gene symbols, duplicate probes averaged |
| 1 | `R/01_download_GSE5406.R` | GEO series matrix via `GEOquery`, `hgu133a.db` probe → symbol mapping, duplicates averaged |
| 2 | `python/02_label_GSE57338.py`, `02_label_GSE5406.py` | Transposes to samples × genes, attaches binary labels |
| 3 | `python/03_check_missing_values.py` | Checks for missing values |
| 3 | `python/03_outlier_detection.py` | Ward dendrogram plus inter-sample-correlation z-scores (z < −2.5) |
| 4 | `python/04_top_variance_genes.py` | Keeps the 5,000 most variable genes, then PCA / t-SNE |
| 5 | `python/05_anova_top50.py` | ANOVA F-test keeps the top 50 discriminating genes, then PCA / t-SNE |
| 6 | `python/06_align_external_test_set.py` | Aligns GSE5406 to the 5,000 training genes |
| 7 | `R/07_pathway_enrichment.R` | Reactome + KEGG over-representation analysis of the 25-gene signature (`clusterProfiler`, `ReactomePA`) |
| – | `R/probe_to_symbol_map.R`, `python/04b_probe_ids_to_symbols.py` | Alternative probe-ID → symbol mapping for matrices that still have probe IDs |

## Results: pathway enrichment of the 25-gene signature

Signature: HMGN2, SERPINA3, LAD1, FREM1, S1PR3, HMOX2, SH2D1B, S100A8, PGLYRP1, CXCL8, PROK2, VCAN, THBS1, LCN2, CASTOR2, NUB1, ANXA3, TNFAIP6, FCGR3B, ADGRE2, VNN1, PADI4, HCK, FCN1, GPR84.

| Database | Pathway | Genes | Adjusted p |
|---|---|---|---|
| Reactome | Neutrophil degranulation | 10 / 22 | 6.6 × 10⁻⁷ |
| Reactome | Antimicrobial peptides | 3 / 22 | 0.044 |
| KEGG | IL-17 signalling pathway | 3 / 15 | 0.019 |
| KEGG | Phagocytosis | 4 / 15 | 0.019 |
| KEGG | Bladder cancer | 2 / 15 | 0.039 |
| KEGG | Malaria | 2 / 15 | 0.043 |

Full tables are in [`results/pathways/`](results/pathways/).

---

## Getting started

### 1. Environment

```bash
git clone https://github.com/dunciu/heart-failure-gene-signature.git
cd heart-failure-gene-signature
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
Rscript R/00_install_packages.R                        # R ≥ 4.3
```

### 2. Data

Expression data is public on GEO but too large for git (the GSE57338 raw archive alone is 1.4 GB), so it isn't included:

- **GSE57338:** download `GSE57338_RAW.tar` from the [GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57338) (Supplementary file) into `data/`.
- **GSE5406:** downloaded automatically by `R/01_download_GSE5406.R`.

The two metadata files (sample labels) are already in `data/`.

### 3. Run (from the repository root)

```bash
python python/00_extract_cel_files.py
Rscript R/01_rma_normalise_GSE57338.R
Rscript R/01_download_GSE5406.R
python python/02_label_GSE57338.py
python python/02_label_GSE5406.py
python python/03_check_missing_values.py
python python/03_outlier_detection.py GSE57338
python python/04_top_variance_genes.py
python python/05_anova_top50.py
python python/06_align_external_test_set.py
Rscript R/07_pathway_enrichment.R
```

Every script reads from and writes to `data/` (ignored by git), except the pathway results, which go to `results/pathways/`.

## Notes and caveats

- **Cross-platform test set:** GSE57338 is RMA-normalised from raw CEL files, while GSE5406 uses the GEO series-matrix values as provided. Check that both are on a comparable log2 scale before applying a trained model.
- **Missing genes in the test set:** `06_align_external_test_set.py` fills training genes that are absent from GSE5406 with **0**. On a log2 expression scale that's an extreme value, not a neutral one. Restricting both cohorts to their shared genes, or imputing with the training mean, avoids this.
- **Imbalanced test cohort:** GSE5406 has only 16 non-failing samples out of 210, so report balanced accuracy, AUC and per-class recall, not plain accuracy.
- **Feature selection and evaluation:** the ANOVA top-50 selection uses labels. Any cross-validated estimate on GSE57338 must run it inside each fold. The external GSE5406 test isn't affected.

## Acknowledgements

Supervised by Dr Zeyneb Kurt (University of Sheffield). Data from the NCBI Gene Expression Omnibus: GSE57338 (Liu et al., 2015) and GSE5406 (Hannenhalli et al., 2006).

## Author

**Duncan Alexander**, BSc Data Science, University of Sheffield · [github.com/dunciu](https://github.com/dunciu)
