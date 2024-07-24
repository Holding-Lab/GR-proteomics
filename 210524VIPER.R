#21/05/2024 Reproduce VIPER analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("aracne.networks")
BiocManager::install("mixtools")
library(viper)
library(aracne.networks)

#Bladder cancer blca_msk_tcga_2020
setwd("~/Downloads/blca_msk_tcga_2020")
getwd()
Blca <- read.table(
  "data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt", 
  header = TRUE)
View(Blca)
str(Blca) #data.frame':	20490 obs. of  298 variables:
duplicates <- duplicated(Blca[, 2])
if (any(duplicates)) {
  cat("Duplicates found in the first column:\n")
  print(Blca[duplicates, 2])
} else {
  cat("No duplicates found in the first column.\n")
}

Blca_unique <- Blca[!duplicated(Blca[, c('Entrez_Gene_Id')]), ]
dim(Blca_unique) #[1] 20443   298
View(Blca_unique)
rownames(Blca_unique) <- Blca_unique[, 2]
Blca_unique <- Blca_unique[, -c(1, 2)]
Blca_unique <- Blca_unique[!duplicated(Blca_unique), ]
Blca_unique <- Blca_unique[-c(1),]
str(Blca_unique) #'data.frame':	19847 obs. of  296 variables:
View(Blca_unique)

any_missing <- any(is.na(Blca_unique))
if (any_missing) {
  cat("There are missing values or NAs in the dataframe.")
} else {
  cat("There are no missing values or NAs in the dataframe.")
} #There are no missing values or NAs in the dataframe.

nBlca <- as.matrix(Blca_unique)
str(nBlca)
class(nBlca)
dim(nBlca)
BlcaRes<-viper(nBlca,regulonblca, verbose =TRUE)
dim(BlcaRes)
str(BlcaRes)

vBlcaRes <- t(BlcaRes)
df_BlcaRes<-as.data.frame(vBlcaRes, row.names = NULL, optional = FALSE,
                        make.names = TRUE,
                        stringsAsFactors = default.stringsAsFactors())
View(df_BlcaRes) #5531 genes
str(df_BlcaRes)

# Initialize an empty dataframe to store the results
results <- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)

# Loop through each variable (column) in the dataframe
for (var in colnames(df_BlcaRes)) {
  if (var != "2908") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_BlcaRes[[var]], df_BlcaRes[["2908"]])
    
    # Store the results
    results <- rbind(results, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}

# Print the results
print(results)
View(results)
dim(resultss)
#Breast cancers 
#positive control sample datasets

print(regulonbrca)

print(dset)
View(dset)
class(dset)
str(dset)
vpres2 <- viper(dset, regulonbrca, verbose = TRUE) #worked well
vBrcaRes <- t(vpres2)

df_BrcaRes<-as.data.frame(vBrcaRes, row.names = NULL, optional = FALSE,
                        make.names = TRUE,
                        stringsAsFactors = default.stringsAsFactors())
dim(df_BrcaRes)
# Initialize an empty dataframe to store the results
resultss <- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)

# Loop through each variable (column) in the dataframe
for (var in colnames(df_BrcaRes)) {
  if (var != "2908") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_BrcaRes[[var]], df_BrcaRes[["2908"]])
    
    # Store the results
    resultss <- rbind(resultss, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}
# Print the results
dim(resultss)
View(resultss)

#T CELL Leukemia

alleset <- read.table("tStJude_fpkm.txt", header = TRUE, sep = "\t") #data.frame
View(alleset)
row.names(alleset) <-alleset$Gene
alleset <- alleset[, -1]
alleset<-as.matrix(alleset)
regulonall <- aracne2regulon(afile = "ALLnetwork.txt", 
                            alleset, format = "3col")

str(regulonall)
class(regulonall)
print(regulonall) #5334 regulators, 17851 targets and 300641 interactions

Allviper <- viper(alleset, regulonall, verbose = TRUE)
class(Allviper)
tAllviper<-t(Allviper)
df_TallRes<-as.data.frame(tAllviper)
View(df_TallRes)
# Initialize an empty dataframe to store the results
resultsss <- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)

# Loop through each variable (column) in the dataframe
for (var in colnames(df_TallRes)) {
  if (var != "NR3C1") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_TallRes[[var]], df_TallRes[["NR3C1"]])
    
    # Store the results
    resultsss <- rbind(resultsss, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}
# Print the results
dim(resultsss)
View(resultsss)

BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library(AnnotationDbi)
library(org.Hs.eg.db)

gene_names <- rownames(resultsss)

# Map gene symbols to Entrez IDs
gene_ids <- mapIds(org.Hs.eg.db, keys = gene_names, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Debug: Check the mapping
cat("First few gene ID mappings:\n")
print(head(gene_ids))

# Debug: Check valid indices
cat("Number of valid indices:\n")
print(sum(valid_indices))

valid_indices <- !is.na(gene_ids)
if (all(!valid_indices)) {
  stop("None of the gene names could be mapped to Entrez IDs.")
}

resultsss_valid <- resultsss[valid_indices, ]

# Update row names to Entrez IDs
rownames(resultsss_valid) <- gene_ids[valid_indices]

# Debug: Check the final row names of resultsss_valid
cat("First few row names of resultsss after mapping:\n")
print(head(rownames(resultsss_valid)))

View(resultsss_valid)

print(resultsss_valid)
# Find the shared row names across results, resultss, and resultsss_valid
shared_row_names <- Reduce(intersect, list(rownames(results), rownames(resultss), rownames(resultsss_valid),
                                           rownames(result_hrbc),
                                           rownames(result_tnbc)))
length(shared_row_names) #2771

# Convert row names to a column named "ID"
results$ID <- rownames(results)
resultss$ID <- rownames(resultss)
resultsss_valid$ID <- rownames(resultsss_valid)
result_hrbc_new$ID <- rownames(result_hrbc_new)
result_tnbc_new$ID <- rownames(result_tnbc_new)
result_all_new$ID <- rownames(result_all_new)
head(results)

# Merge the dataframes using full outer join
merged_df_1_2 <- merge(results, resultss, by = "ID", all = TRUE)
merged_df_all <- merge(merged_df_1_2, resultsss_valid, by = "ID", all = TRUE)
merged_df_all_tnbc <- merge(merged_df_all, result_tnbc, by = "ID", all = TRUE)
View(merged_df_all)

# Install and load purrr package
install.packages("purrr")
library(purrr)


dfs <- list(results, resultss, result_all_new, result_hrbc_new, result_tnbc_new)

# Merge all data frames by retaining intersecting rows

merged_results <- reduce(dfs, function(x, y) merge(x, y, by = "ID"))

# Print the first few rows of the merged data frame
View(merged_results)
str(merged_results)





rownames(merged_results) <- merged_results$ID
merged_results$ID <- NULL
str(merged_results)
colnames(merged_results) <- c("r_blca", "p_blca", "r_brca", "p_brca", "r_tall", "p_tall",
                              "r_hrbc", "p_hrbc", "r_tnbc", "p_tnbc")
View(merged_results)
write.csv(merged_results, file = "VIPER5CANCERPR.csv", row.names = TRUE)
#subset genes which are non-significant across all sample types 
notsig_results <- merged_results %>% filter(merged_results$p_blca > 0.05,
         merged_results$p_brca > 0.05,
         merged_results$p_tall > 0.05,
         merged_results$p_hrbc > 0.05,
         merged_results$p_tnbc > 0.05)
#no such genes, all genes are significant in at least one sample

write.csv(merged_df_all,file = "viper3cancers.csv", row.names = TRUE)
# Install and load pheatmap if not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)
View(merged_results)
# Select only the r_value columns for the heatmap
heatmap_data <- merged_results[, c("r_blca", "r_brca", "r_tall","r_hrbc","r_tnbc")]
str()
# Drop rows with any NA values
heatmap_data_clean <- na.omit(heatmap_data)

# Define color palette with a specific color for NA values
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_colors <- c(heatmap_colors, "grey") # Add grey for NA


pheatmap(as.matrix(heatmap_data_clean), 
        color = heatmap_colors, 
        main = "Heatmap of Pearson r values co-regulator vs GR activity"
        )


viperr <- pheatmap(
  as.matrix(heatmap_data_clean), 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  main = "VIPER r Heatmap",
  color = heatmap_colors,
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







# Install VennDiagram package if not already installed
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)
row_names_results <- rownames(results)
row_names_resultss <- rownames(resultss)
row_names_resultsss_valid <- rownames(resultsss_valid)
# Create a list of row names
venn_list <- list(
  Blca = row_names_results,
  Brca = row_names_resultss,
  T_ALL = row_names_resultsss_valid
)


venn.diagram(
  x = venn_list,
  category.names = c("Blca", "Brca", "T-ALL"),
  filename = "venn_diagram.png",
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw"
)

cat("Venn diagram saved as venn_diagram.png\n")
getwd()
E239_trim<-read.csv('E239DIAtrim.csv',sep=',', stringsAsFactors = F)

gene_239 <- rownames(E239_trim)

# Map gene symbols to Entrez IDs
ids_239 <- mapIds(org.Hs.eg.db, keys = gene_239, column = "ENTREZID", keytype = "UNIPORT", multiVals = "first")

#HR+ HER- BC
HRBCeset <- read.table("Brca_hr_clean_exp.txt", header = TRUE, sep = "\t") #data.frame
View(HRBCeset)
row.names(HRBCeset) <-HRBCeset$Gene
HRBCeset <- HRBCeset[, -1]
HRBCeset<-as.matrix(HRBCeset)

#HRBCnetwork.txt is worong due to misuse of TF list
regulonhrbc <- aracne2regulon(afile = "HRBCnetwork.txt", 
                              HRBCeset, format = "3col")
str(regulonhrbc)
print(regulonhrbc) #14674 regulators, 14674 targets and 759102 interactions

HRBCviper <- viper(HRBCeset, regulonhrbc, verbose = TRUE)
class(HRBCviper)
dim(HRBCviper) #[1] 13606   332
tHRBCviper<-t(HRBCviper)
df_HRBCRes<-as.data.frame(tHRBCviper)
View(df_HRBCRes)

regulonhrbc_new <- aracne2regulon(afile = "HR_ARACNE_network.txt", 
                              HRBCeset, format = "3col")
print(regulonhrbc_new) # 273 regulators, 14674 targets and 134563 interactions

HRBCviper_new <- viper(HRBCeset, regulonhrbc_new, verbose = TRUE)
dim(HRBCviper_new) #[1] 273 332
tHRBCviper_new<-t(HRBCviper_new)
df_HRBCRes_new<-as.data.frame(tHRBCviper_new)

result_hrbc_new<- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)
# Loop through each variable (column) in the dataframe
for (var in colnames(df_HRBCRes_new)) {
  if (var != "2908") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_HRBCRes_new[[var]], df_HRBCRes_new[["2908"]])
    
    # Store the results
    result_hrbc_new <- rbind(result_hrbc_new, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}

# Print the results
View(result_hrbc_new)

result_hrbc<- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)
# Loop through each variable (column) in the dataframe
for (var in colnames(df_HRBCRes)) {
  if (var != "2908") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_HRBCRes[[var]], df_HRBCRes[["2908"]])
    
    # Store the results
    result_hrbc <- rbind(result_hrbc, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}

# Print the results
View(result_hrbc)



#TNBC
TNBCeset <- read.table("Brca_TNBC_clean_exp.txt", header = TRUE, sep = "\t") #data.frame
View(TNBCeset)
row.names(TNBCeset) <-TNBCeset$Gene
TNBCeset <- TNBCeset[, -1]
TNBCeset<-as.matrix(TNBCeset)

#TNBCnetwork.txt is worong due to misuse of TF list
regulontnbc <- aracne2regulon(afile = "TNBCnetwork.txt", 
                              TNBCeset, format = "3col")
str(regulontnbc)
print(regulontnbc) #14674 regulators, 14674 targets and 1008889 interactions

TNBCviper <- viper(TNBCeset, regulontnbc, verbose = TRUE)
class(TNBCviper)
dim(TNBCviper) #[1] 14635    89
tTNBCviper<-t(TNBCviper)
df_TNBCRes<-as.data.frame(tTNBCviper)
View(df_TNBCRes)

#Use correct TNBC ARACNE networks
getwd() #"/Users/weiye/Downloads/blca_msk_tcga_2020"
regulontnbc_new <- aracne2regulon(afile = "TNBC_ARACNE_network.txt", 
                              TNBCeset, format = "3col")

str(regulontnbc_new )
print(regulontnbc_new) # 273 regulators, 14468 targets and 75054 interactions

TNBCviper_new <- viper(TNBCeset, regulontnbc_new, verbose = TRUE)
class(TNBCviper_new)
dim(TNBCviper_new) #[1] 273  89
tTNBCviper_new<-t(TNBCviper_new)
df_TNBCRes_new<-as.data.frame(tTNBCviper_new)
View(df_TNBCRes_new)

# Initialize an empty dataframe to store the results
result_tnbc_new<- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)

# Loop through each variable (column) in the dataframe
for (var in colnames(df_TNBCRes_new)) {
  if (var != "2908") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_TNBCRes_new[[var]], df_TNBCRes_new[["2908"]])
    
    # Store the results
    result_tnbc_new <- rbind(result_tnbc_new, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}

# Print the results
View(result_tnbc_new)

# Initialize an empty dataframe to store the results
result_tnbc<- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)

# Loop through each variable (column) in the dataframe
for (var in colnames(df_TNBCRes)) {
  if (var != "2908") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_TNBCRes[[var]], df_TNBCRes[["2908"]])
    
    # Store the results
    result_tnbc <- rbind(result_tnbc, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}

# Print the results
View(result_tnbc)
dim(result_tnbc)


#re-do T-ALL 23/06/2024
#Use correct ALL ARACNE networks

regulonall_new <- aracne2regulon(afile = "ALL_ARCNE_network.txt", 
                                 alleset, format = "3col")

str(regulonall_new)
print(regulonall_new) # 271 regulators, 16855 targets and 123339 interactions

ALLviper_new <- viper(alleset, regulonall_new, verbose = TRUE)
dim(ALLviper_new) #[1] 243  264
tALLviper_new<-t(ALLviper_new)
df_ALLRes_new<-as.data.frame(tALLviper_new)
View(df_ALLRes_new)

# Initialize an empty dataframe to store the results
result_all_new<- data.frame(r_value = numeric(), p_value = numeric(), row.names = NULL)

# Loop through each variable (column) in the dataframe
for (var in colnames(df_ALLRes_new)) {
  if (var != "NR3C1") {  # Skip the target variable itself
    # Perform Pearson correlation test
    test <- cor.test(df_ALLRes_new[[var]], df_ALLRes_new[["NR3C1"]])
    
    # Store the results
    result_all_new <- rbind(result_all_new, data.frame(r_value = test$estimate, p_value = test$p.value, row.names = var))
  }
}

# Print the results
View(result_all_new)

install.packages('dendextend')
install.packages('gplots')
install.packages('viridis')
install.packages("devtools")
library('dendextend')
library('gplots')
library('viridis')
library("devtools")


# Convert the relevant columns to a numeric matrix
heatmap_datafc2 <- as.matrix(diff001fc2[, adj.P.Val.columns001fc2])
View(heatmap_datafc2)
str(heatmap_datafc2)
# Perform hierarchical clustering (rows and columns)
row_dendrogramfc2 <- hclust(dist(heatmap_datafc2))
col_dendrogramfc2 <- hclust(dist(t(heatmap_datafc2)))
str(row_dendrogramfc2)
plot(row_dendrogramfc2, main = "Row Dendrogram")


# df [heatmap_data] is the r values df. 
heatmap_data_matrix <- as.matrix(heatmap_data)
View(heatmap_data_matrix)
row_heatmap_data<- hclust(dist(heatmap_data_matrix))
col_den_heatmap_data <- hclust(dist(t(heatmap_data_matrix)))

plot(row_heatmap_data, main = "Row Dendrogram r values")


den<-as.dendrogram(row_heatmap_data)
cut_tree <- dendextend::cutree(
  den,
  h = 2,
  order_clusters_as_data = FALSE
)
table(cut_tree)


desired_branch <- dendextend::cutree(
  den,
  k = 8,
  order_clusters_as_data = FALSE)

print(desired_branch)
str(desired_branch) #Named int

# Extract the names and values from desired_branch
cluster_ids <- as.integer(desired_branch)
gene_ids <- names(desired_branch)

# Create a data frame
desired_branch_df <- data.frame(GeneID = gene_ids, ClusterID = cluster_ids, stringsAsFactors = FALSE)

# Display the data frame
head(desired_branch_df)
write.csv(desired_branch_df, "5cancervipercluster.csv", row.names = FALSE)

