# One-off setup: CRAN + Bioconductor packages used by the R scripts.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages(c("dplyr", "tibble"))
BiocManager::install(c(
  "GEOquery", "Biobase", "AnnotationDbi",
  "oligo", "pd.hugene.1.1.st.v1", "hugene11sttranscriptcluster.db",   # GSE57338 (HuGene 1.1 ST)
  "hgu133a.db",                                                         # GSE5406 (HG-U133A)
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot"         # pathway enrichment
))
