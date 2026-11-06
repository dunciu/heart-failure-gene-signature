# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)    # Human genome database
  library(ReactomePA)      # Reactome pathway database
  library(enrichplot)      # For plotting
  library(dplyr)           # For data manipulation
})

print("Libraries loaded. Ready for pathway analysis.")

# 2. DEFINE YOUR 25-GENE MASTER SIGNATURE
my_genes <- c(
  "HMGN2", "SERPINA3", "LAD1", "FREM1", "S1PR3", "HMOX2", 
  "SH2D1B", "S100A8", "PGLYRP1", "CXCL8", "PROK2", "VCAN", 
  "THBS1", "LCN2", "CASTOR2", "NUB1", "ANXA3", "TNFAIP6", 
  "FCGR3B", "ADGRE2", "VNN1", "PADI4", "HCK", "FCN1", "GPR84"
)

# 3. CONVERT SYMBOLS TO ENTREZ IDs
# Databases like KEGG and Reactome strictly require Entrez IDs (numbers) instead of text symbols.
gene_conversion <- bitr(my_genes, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)

entrez_ids <- gene_conversion$ENTREZID
print(paste("Successfully mapped", length(entrez_ids), "out of 25 genes."))

# ==========================================
# 4. PERFORM OVER-REPRESENTATION ANALYSIS (ORA)
# ==========================================
# This function automatically runs the Fisher's Exact Test and calculates the Adjusted p-value!

# A. Reactome Database (Excellent for metabolism & lipids)
reactome_res <- enrichPathway(gene = entrez_ids, 
                              organism = "human", 
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.05, 
                              readable = TRUE) # This cleanly turns Entrez IDs back to Gene Symbols in the output

# B. KEGG Database (Excellent for signaling cascades)
kegg_res <- enrichKEGG(gene = entrez_ids, 
                       organism = "hsa", 
                       pvalueCutoff = 0.05)
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# ==========================================
# 5. TARGETED HYPOTHESIS TESTING (LIPIDS)
# ==========================================
# Convert results to a standard dataframe to search for lipid/fatty acid keywords
reactome_df <- as.data.frame(reactome_res)

lipid_pathways <- reactome_df %>%
  filter(grepl("lipid|fatty acid|sphingolipid|cholesterol", Description, ignore.case = TRUE))

print("\n=== SPECIFIC LIPID/FATTY ACID PATHWAYS FOUND ===")
if(nrow(lipid_pathways) > 0) {
  print(lipid_pathways[, c("ID", "Description", "p.adjust", "geneID")])
} else {
  print("No strictly lipid-named pathways passed the p.adjust < 0.05 threshold. Check broader pathways!")
}

# ==========================================
# 6. EXPORT AND VISUALIZE
# ==========================================
# Save the full results to a CSV for your dissertation appendix
dir.create("results/pathways", recursive = TRUE, showWarnings = FALSE)
write.csv(reactome_df, "results/pathways/Reactome_Pathway_Analysis_Results.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_res), "results/pathways/KEGG_Pathway_Analysis_Results.csv", row.names = FALSE)

# Generate a publishable Dot Plot for the top 15 Reactome pathways
dotplot(reactome_res, showCategory = 15, 
        title = "Top Enriched Reactome Pathways\n(25-Gene Consensus Signature)")