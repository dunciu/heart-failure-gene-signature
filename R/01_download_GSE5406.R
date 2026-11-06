# --- INSTALLATION (Run this once in the console if you haven't yet) ---
# BiocManager::install("hgu133a.db")

# --- LIBRARIES ---
suppressPackageStartupMessages({
  library(GEOquery)
  library(hgu133a.db) 
  library(dplyr)
  library(tibble)
  library(Biobase)
  library(AnnotationDbi)
})

print("Libraries loaded successfully.")

# --- STEP 1: DOWNLOAD DATA ---
gse_id <- "GSE5406"
print(paste("Downloading", gse_id, "..."))
gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

data_matrix <- gset[[1]]

# --- STEP 2: EXTRACT RAW MATRICES ---
# Expression data (The numbers)
expression_data <- exprs(data_matrix)

# Metadata (The sample info: Heart Failure vs Non-Failing)
pheno_data <- pData(data_matrix)

# --- STEP 3: MAP PROBES TO GENES ---
probe_ids <- rownames(expression_data)

# Aggressively clean the IDs: force to character, remove spaces, and strip quotes
clean_probe_ids <- as.character(probe_ids)
clean_probe_ids <- trimws(clean_probe_ids)
clean_probe_ids <- gsub('\"', '', clean_probe_ids)

print("Diagnostic: First 5 cleaned probe IDs...")
print(head(clean_probe_ids))

# Map using the correct database for Affymetrix U133A Array
gene_mapping <- AnnotationDbi::select(hgu133a.db,
                                      keys = clean_probe_ids,
                                      columns = c("SYMBOL"),
                                      keytype = "PROBEID")

# Remove lines where no gene symbol was found (NA)
gene_mapping <- na.omit(gene_mapping)

# --- STEP 4: MERGE & CLEAN ---
# Convert expression matrix to a dataframe
expr_df <- as.data.frame(expression_data) %>%
  rownames_to_column(var = "PROBEID")

# Clean the PROBEID column in the dataframe so it perfectly matches the mapping
expr_df$PROBEID <- trimws(gsub('\"', '', as.character(expr_df$PROBEID)))

# Join the gene symbols with the data
clean_data <- inner_join(gene_mapping, expr_df, by = "PROBEID")

# --- STEP 5: HANDLE DUPLICATES ---
# Multiple probes can map to one gene. We average them to get one value per gene.
final_dataset <- clean_data %>%
  select(-PROBEID) %>%               
  group_by(SYMBOL) %>%               
  summarise(across(everything(), mean)) %>% 
  column_to_rownames(var = "SYMBOL") 

print(paste("Cleaning complete. Final dataset has", nrow(final_dataset), "unique genes."))

# --- STEP 6: EXPORT TO CSV ---
# 1. Expression Matrix 
write.csv(final_dataset, "data/GSE5406_Expression_Matrix.csv", row.names = TRUE)

# 2. Metadata (The Labels)
clean_metadata <- pheno_data %>%
  select(geo_accession, title, source_name_ch1, characteristics_ch1)

write.csv(clean_metadata, "data/GSE5406_Metadata.csv", row.names = FALSE)

print("SUCCESS: Files 'GSE5406_Expression_Matrix.csv' and 'GSE5406_Metadata.csv' created.")