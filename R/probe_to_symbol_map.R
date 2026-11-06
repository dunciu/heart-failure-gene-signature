# We can comment out the install line since you already have it!
# BiocManager::install("hugene11sttranscriptcluster.db")

suppressPackageStartupMessages({
  library(hugene11sttranscriptcluster.db)
  library(AnnotationDbi)
})

# Explicitly tell R to use AnnotationDbi for both keys() and select()
all_probes <- AnnotationDbi::keys(hugene11sttranscriptcluster.db, keytype="PROBEID")
mapping <- AnnotationDbi::select(hugene11sttranscriptcluster.db, 
                                 keys=all_probes, 
                                 columns="SYMBOL", 
                                 keytype="PROBEID")

# Drop any probes that don't have a known gene name
mapping <- na.omit(mapping)

# Save as a simple CSV dictionary
write.csv(mapping, "data/Probe_to_Symbol_Map.csv", row.names=FALSE)
print("Success! Mapping dictionary saved as 'Probe_to_Symbol_Map.csv'.")