# Run after R/07_pathway_enrichment.R in the same R session (uses reactome_df and kegg_res).
# Show the top 10 most significant Reactome pathways
print("=== TOP 10 REACTOME PATHWAYS ===")
head(reactome_df[, c("Description", "p.adjust", "Count")], 10)

# Show the top 10 most significant KEGG pathways
print("=== TOP 10 KEGG PATHWAYS ===")
head(as.data.frame(kegg_res)[, c("Description", "p.adjust", "Count")], 10)