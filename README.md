# PWA_FDR
Basic R function to calculate FDR for Pathways analysis.
Just used the gene list as an example, that gene list which are generted from LLMs using prompts. 
#need to Load necessary libraries
library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
#######################################
# Define function
pathway_analysis <- function(gene_list, plot_bar = TRUE, plot_dot = TRUE) {
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL")
  
  # Remove any potential NAs from entrez_ids
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # Perform pathway analysis
  enrich_result <- enrichKEGG(gene = entrez_ids, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05)
  
  # Calculate FDR-adjusted p-values
  enrich_result@result$padj <- p.adjust(enrich_result@result$pvalue, method = "BH")
  
  # Plot bar plot if specified
  if (plot_bar) {
    bar_plot <- barplot(enrich_result, showCategory = 15, main = "KEGG PWA", ylab = "FDR")
    print(bar_plot)
  }
  
  # Plot dot plot if specified
  if (plot_dot) {
    dot_plot <- dotplot(enrich_result, showCategory = 15)
    dot_plot <- dot_plot + ggtitle("KEGG PWA") + ylab("FDR")
    print(dot_plot)
  }
  
  # Return results
  return(enrich_result@result)
}

# Example usage
gene_list <- c("HLA-A01:01", "CD40LG", "IL1B", "TNF", "IFNG", "CXCL9", "CXCL10", "CXCR3", "CD27", "CD154", "TNFAIP3", "CXCL10", "IL6", "IL1B", "NFKB1", "JUN", "TNF", "FASLG", "TRAF1")
results <- pathway_analysis(gene_list)
print(results)
####################################################################################################################################

