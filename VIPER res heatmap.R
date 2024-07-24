getwd()
head(rownames(result_all_new))
#VIPER analysis results plotting heatmaps
# Map Entrez IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(result_all_new),
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")

# Replace rownames with gene symbols
rownames(result_all_new) <- gene_symbols

# Display the updated data frame
head(result_all_new)


#The significance of correlation vs the strength of correlation

# Install and load pheatmap if not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

merged_results<-read.csv('VIPER5CANCERPR.csv',sep=',', stringsAsFactors = F)

View(merged_results)
str(merged_results)
# Select only the r_value columns for the heatmap
hmsig_data <- merged_results[, c("p_blca", "p_brca", "p_tall","p_hrbc","p_tnbc")]
View()
# Drop rows with any NA values
hmsig_data_clean <- na.omit(hmsig_data)
dim(hmsig_data_clean) #188 5


# Define color breaks
color_breaks <- c(0, 0.001, 0.01, 0.05, 1)

# Define color palette (needs one less color than the number of breaks)
heatmap_colors <- c("#86ebc9", "#869ceb", "grey", "white")


# Create the heatmap
pheatmap(as.matrix(hmsig_data_clean ), 
         color = heatmap_colors, 
         main = "Heatmap of Pearson p values co-regulator vs GR activity",
         breaks = color_breaks
)

viperp <- pheatmap(
  hmsig_data_clean, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  main = "VIPER p Heatmap",
  color = heatmap_colors,
  breaks = color_breaks,
  show_rownames = FALSE,  # Remove row names
  border_color = "black",  # Remove borders around cells
  cellwidth = 20,  
  cellheight = NA,  # Automatic height for cells
  treeheight_row = 100,  # Adjust the dendrogram height for rows
  treeheight_col = 50,  # Adjust the dendrogram height for columns
  fontsize_row = 10,  # Adjust row font size (if row names were to be displayed)
  fontsize_col = 10,  # Adjust column font size
  fontsize = 12,  # Adjust overall font size
  clustering_distance_rows = "euclidean",  # Distance metric for rows
  clustering_distance_cols = "euclidean",  # Distance metric for columns
  clustering_method = "complete"  # Clustering method
)









sigcor_blca <- merged_results %>%
  filter(abs(merged_results$r_blca) > 0.5, merged_results$p_blca < 0.05)
View(sigcor_blca)


sigcor_pos <- merged_results %>%
  filter(merged_results$r_blca > 0.5, merged_results$p_blca < 0.05,
         merged_results$r_brca > 0.5, merged_results$p_brca < 0.05,
         merged_results$r_tall > 0.5, merged_results$p_tall < 0.05,
         merged_results$r_hrbc > 0.5, merged_results$p_hrbc < 0.05,
         merged_results$r_tnbc > 0.5, merged_results$p_tnbc < 0.05)
View(sigcor_pos)

sigcor_neg <- merged_results %>%
  filter(merged_results$r_blca < -0.5, merged_results$p_blca < 0.05,
         merged_results$r_brca < -0.5, merged_results$p_brca < 0.05,
         merged_results$r_tall < -0.5, merged_results$p_tall < 0.05,
         merged_results$r_hrbc < -0.5, merged_results$p_hrbc < 0.05,
         merged_results$r_tnbc < -0.5, merged_results$p_tnbc < 0.05)
View(sigcor_neg)

sigcor_pos_brca <- merged_results %>%
  filter(
         merged_results$r_brca > 0.5, merged_results$p_brca < 0.05,
         
         merged_results$r_hrbc > 0.5, merged_results$p_hrbc < 0.05,
         merged_results$r_tnbc > 0.5, merged_results$p_tnbc < 0.05)
View(sigcor_pos_brca)
print(rownames(sigcor_pos_brca))

sigcor_neg_brca <- merged_results %>%
  filter(
         merged_results$r_brca < -0.5, merged_results$p_brca < 0.05,
         
         merged_results$r_hrbc < -0.5, merged_results$p_hrbc < 0.05,
         merged_results$r_tnbc < -0.5, merged_results$p_tnbc < 0.05)
View(sigcor_neg_brca)





