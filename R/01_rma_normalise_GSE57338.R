# --- STEP 1: LOAD LIBRARIES ---
library(oligo)
library(Biobase)
# Swapped to 1.1 to correctly match GSE57338!
library(pd.hugene.1.1.st.v1)
library(hugene11sttranscriptcluster.db)

# --- STEP 2: READ FILES ---
# Run from the repository root; CEL files come from python/00_extract_cel_files.py
setwd("data/unzipped_cel_files")
cel_files <- list.celfiles(listGzipped = FALSE) 
raw_data <- read.celfiles(cel_files)

# --- STEP 3: BACKGROUND CORRECT, NORMALIZE, AND LOG2 ---
print("Running full RMA (Background Correction, Normalization, Log2)...")

# Standard rma() automatically handles the big three preprocessing steps!
eset_norm <- rma(raw_data)

# Extract the fully normalized, log2-transformed matrix
expression_matrix <- exprs(eset_norm) 

# --- STEP 4: MAP PROBES TO GENES ---
print("Mapping Probes to Genes...")
gene_symbols <- mapIds(hugene11sttranscriptcluster.db, 
                       keys = rownames(expression_matrix), 
                       column = "SYMBOL", 
                       keytype = "PROBEID", 
                       multiVals = "first")

# --- STEP 5: CLEAN AND SAVE ---
final_df <- data.frame(GeneSymbol = gene_symbols, expression_matrix)

# Drop any probes that don't match to a known gene
final_df <- final_df[!is.na(final_df$GeneSymbol), ]

# Average the expressions if multiple probes map to the exact same gene
final_df <- aggregate(. ~ GeneSymbol, data = final_df, FUN = mean)

print("Saving NORMALIZED and MAPPED data...")
write.csv(final_df, "../GSE57338_Cleaned_Normalized.csv", row.names = FALSE)  # -> data/
setwd("../..")
print("Done! You now have clean, machine-learning ready data.")