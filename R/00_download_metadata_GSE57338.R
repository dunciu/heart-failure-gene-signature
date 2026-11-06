# Downloads sample labels for GSE57338 (the CEL files themselves come from GSE57338_RAW.tar).
# Writes data/GSE57338_Metadata.csv with geo_accession, title, source_name_ch1, characteristics_ch1.
# NOTE: reconstructed for this repo (the original download step wasn't saved); a copy of the
# metadata file is already in data/.
suppressPackageStartupMessages({ library(GEOquery); library(Biobase); library(dplyr) })
gset <- getGEO("GSE57338", GSEMatrix = TRUE, AnnotGPL = FALSE)
meta <- pData(gset[[1]]) %>% select(geo_accession, title, source_name_ch1, characteristics_ch1)
write.csv(meta, "data/GSE57338_Metadata.csv", row.names = FALSE)
print(table(meta$source_name_ch1))
