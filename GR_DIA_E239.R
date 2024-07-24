R.version
# Set the number of significant digits
options(digits = 20)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qPLEXanalyzer")
library(qPLEXanalyzer)
library(gridExtra)
data(exp2_Xlink)
View(exp2_Xlink$metadata)
setwd('~/Downloads')
x<-read.delim('E239_DIA_Report_MinImputation.tsv',sep='\t', stringsAsFactors = F)
View(x)
ncol(x)
xfiltered <- subset(x, !grepl(";", Genes)) 

duplicated_genes <- duplicated(xfiltered$Genes)
View(duplicated_genes)
duplicated_rows <- xfiltered[duplicated_genes, ]
print(duplicated_rows)
str(duplicated_rows)
View(duplicated_rows)
unique_xfiltered <- xfiltered[!duplicated(xfiltered$Genes), ]
View(unique_xfiltered)
unique_xfiltereded<- subset(unique_xfiltered, !grepl(";", Protein.Group))
View(unique_xfiltereded)
names(unique_xfiltereded)[names(unique_xfiltereded) == "Protein.Group"] <- "Accessions"
View(unique_xfiltereded)
write.csv(unique_xfiltereded, file = "E239DIAtrim.csv", row.names = FALSE)


#generate metadata
meta<-read.delim('E239_FragPipe_experimental_annotation.tsv',sep='\t', stringsAsFactors = F)
View(meta)
colnames(meta)[colnames(meta) == "file"] <- "SampleName"
colnames(meta)[colnames(meta) == "condition"] <- "SampleGroup"
colnames(meta)[colnames(meta) == "replicate"] <- "BioRep"
meta$TechRep <- NA

# Check the metadata sample names
print(meta$SampleName)


# Check the experimental data columns
print(unique_xfiltereded[, 9:52])

# Get the column names of unique_xfiltereded[, 9:52]
desired_column_order <- colnames(unique_xfiltereded[, 9:52])
print(desired_column_order)

row_order <- order(match(meta$SampleName, desired_column_order))
metadata_reordered <- meta[row_order, ]
metadata_reordered <- metadata_reordered[, -(2:3)]
View(metadata_reordered)
write.csv(metadata_reordered, file = "E239meta.csv", row.names = FALSE)
E239meta<-read.csv('E239meta.csv',sep=',', stringsAsFactors = F)
View(E239meta)
# Replace elements in "SampleGroup" containing "GR" with "GR"
E239meta$SampleGroup <- gsub(".*GR.*", "GR_IP", E239meta$SampleGroup)
E239meta$SampleGroup <- gsub(".*IgG.*", "IgG_IP", E239meta$SampleGroup)
E239meta$SampleGroup <- gsub("_", ".", E239meta$SampleGroup)
E239meta$SampleGroup <- gsub(" ", "", E239meta$SampleGroup)
# View the updated dataframe
head(E239meta)
View(E239meta)

E239meta.original<-read.csv('E239meta.csv',sep=',', stringsAsFactors = F)
E239meta.original$SampleGroup <- gsub("_", ".", E239meta.original$SampleGroup)
E239meta.original$SampleGroup <- gsub(" ", "", E239meta.original$SampleGroup)
View(E239meta.original)

write.csv(E239meta.original, file = "E239meta.original.csv", row.names = FALSE)

msnXo<-convertToMSnset(unique_xfiltereded, 
                      metadata= E239meta.original,
                      indExpData=c(9:52), type='protein',Accessions= 1) #use original metadata

P11<-intensityBoxplot(msnXo, title = "No normalization")

MSnset_norm_gso <- groupScaling(msnXo, 
                               scalingFunction = median, 
                               groupingColumn = "SampleGroup") #in total 44 groups 

P12<-intensityBoxplot(MSnset_norm_gso, title = "Within Group Scaling")

hierarchicalPlot(MSnset_norm_gso, 
                 colourBy='SampleGroup',title = "Hierarchical Clustering, Groups Scaled")

contrasto1 <- c('PrimaryBreastEpis.GR vs PrimaryBreastEpis.IgG' 
                = "PrimaryBreastEpis.GR - PrimaryBreastEpis.IgG")
PBEpiso <- computeDiffStats(MSnset_norm_gso, contrasts = contrasto1)
diffPBEpiso <- getContrastResults(PBEpiso, 
                               contrast = contrasto1)

contrasto2 <- c('Primary.NHU.GR vs Primary.NHU.IgG' 
                = "Primary.NHU.GR - Primary.NHU.IgG")
NHUo <- computeDiffStats(MSnset_norm_gso, contrasts = contrasto2)
diffNHUo <- getContrastResults(NHUo, 
                                  contrast = contrasto2)

contrasto3 <- c('PrimaryCD4andT.GR vs PrimaryCD4andT.IgG' 
                = "PrimaryCD4andT.GR - PrimaryCD4andT.IgG")

CD4o <- computeDiffStats(MSnset_norm_gso, contrasts = contrasto3)
diffCD4o <- getContrastResults(CD4o, 
                               contrast = contrasto3)

contrasto4 <- c('MCF7.GR vs MCF7.IgG' 
                = "MCF7.GR - MCF7.IgG")

MCF7o <- computeDiffStats(MSnset_norm_gso, contrasts = contrasto4)
diffMCF7o <- getContrastResults(MCF7o, 
                               contrast = contrasto4)

contrasto5 <- c('KMBC2.GR vs KMBC2.IgG' 
                = "KMBC2.GR - KMBC2.IgG")

KMBC2o <- computeDiffStats(MSnset_norm_gso, contrasts = contrasto5)
diffKMBC2o <- getContrastResults(KMBC2o, 
                                contrast = contrasto5)

contrasto6 <- c('Jurkat.GR vs Jurkat.IgG' 
                = "Jurkat.GR - Jurkat.IgG")

Jurkato <- computeDiffStats(MSnset_norm_gso, contrasts = contrasto6)
diffJurkato <- getContrastResults(Jurkato, 
                                 contrast = contrasto6)

contrasto7 <- c('X231.GR vs X231.IgG' 
                = "X231.GR - X231.IgG")

X231o <- computeDiffStats(MSnset_norm_gso, contrasts = contrasto7)
diffX231o <- getContrastResults(X231o, 
                                  contrast = contrasto7)

olist<-list(diffPBEpiso, diffNHUo, diffCD4o,
            diffMCF7o, diffKMBC2o, diffJurkato, diffX231o)

merged.diffo<-data.frame()

for (i in 1:length(olist)) {
  # Extract the "a" and "b" columns from the current data frame
  ab_columns <- olist[[i]][, c("log2FC", "adj.P.Val")]
  
  # If merged_df is empty, assign ab_columns to it
  if (nrow(merged.diffo) == 0) {
    merged.diffo <- ab_columns
  } else {
    # Otherwise, append ab_columns to the existing merged_df
    merged.diffo <- cbind(merged.diffo, ab_columns)
  }
}

head(merged.diffo)
# Rename columns with sample names
colnames(merged.diffo) <- c("log2FC_Epis", "adj.P.Val_Epis", 
                            "log2FC_NHU", "adj.P.Val_NHU", 
                            "log2FC_CD4", "adj.P.Val_CD4",
                            "log2FC_MCF7", "adj.P.Val_MCF7",
                            "log2FC_KMBC2", "adj.P.Val_KMBC2",
                            "log2FC_Jurkat", "adj.P.Val_Jurkat",
                            "log2FC_231", "adj.P.Val_231")
head(merged.diffo)

# Extract the relevant columns containing adj.P.Val values
adj.P.Val.columns <- colnames(merged.diffo)[grep("adj.P.Val", colnames(merged.diffo))]

# Convert the relevant columns to a numeric matrix
heatmap_data <- as.matrix(merged.diffo[, adj.P.Val.columns])

# Perform hierarchical clustering (rows and columns)
row_dendrogram <- hclust(dist(heatmap_data))
col_dendrogram <- hclust(dist(t(heatmap_data)))

# Create a heatmap with clustering
heatmap(
  heatmap_data,
  Colv = as.dendrogram(col_dendrogram),
  Rowv = as.dendrogram(row_dendrogram),
  scale = "none",  # You can change scaling options as needed
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Choose your color palette
  main = "Heatmap of adj.P.Val",
  xlab = "Samples",
  ylab = "Genes"
)

# Extract the relevant columns containing adj.P.Val values
FC_columns <- colnames(merged.diffo)[grep("log2FC", colnames(merged.diffo))]

# Convert the relevant columns to a numeric matrix
heatmap_dataFC <- as.matrix(merged.diffo[, FC_columns])

# Perform hierarchical clustering (rows and columns)
row_dendrogramFC <- hclust(dist(heatmap_dataFC))
col_dendrogramFC <- hclust(dist(t(heatmap_dataFC)))

# Create a heatmap with clustering
heatmap(
  heatmap_dataFC,
  Colv = as.dendrogram(col_dendrogramFC),
  Rowv = as.dendrogram(row_dendrogramFC),
  scale = "none",  # You can change scaling options as needed
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Choose your color palette
  main = "Heatmap of log2FC",
  xlab = "Samples",
  ylab = "Genes"
)

pheatmap(
  heatmap_data,
  Colv = as.dendrogram(col_dendrogram),
  Rowv = as.dendrogram(row_dendrogram),
  scale = "none",  # You can change scaling options as needed
  col = c("#86ebc9", "#869ceb",
          "#b986eb","white"),# Choose your color palette
  breaks = color_breaks,
  main = "Heatmap of adj.P.Val all sample",
  xlab = "Samples"
)


#drop CD4+ T 001 sample
View(unique_xfiltereded)
drop001 <- unique_xfiltereded[, -c(21, 25)]
# Find the column numbers containing "_GR_IP"
GRIPcolumn_numbers <- grep("_GR_IP", colnames(drop001))


str(GRIPcolumn_numbers)
E239meta.original<-read.csv('E239meta.original.csv',sep=',', stringsAsFactors = F)
E239meta001 <- E239meta.original[-c(13, 17), ]
View(E239meta.original)
View(E239meta001)

msn001<-convertToMSnset(drop001, 
                       metadata= E239meta001,
                       indExpData=c(9:50), type='protein',Accessions= 1)
str(msn001)
convertToMSnset()
intensityBoxplot(msn001, title = "No normalization")

MSnset_gsnorm_001 <- groupScaling(msn001, 
                                scalingFunction = median, 
                                groupingColumn = "SampleGroup") #in total 44 groups 
str(MSnset_gsnorm_001)

intensityBoxplot(MSnset_gsnorm_001, title = "Within Group Scaling")

hierarchicalPlot(MSnset_gsnorm_001, 
                 colourBy='SampleGroup',title = "Hierarchical Clustering, Groups Scaled")

contrasto1 <- c('PrimaryBreastEpis.GR vs PrimaryBreastEpis.IgG' 
                = "PrimaryBreastEpis.GR - PrimaryBreastEpis.IgG")
Epis001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto1)

str(Epis001) #List of 3
diffEpis001 <- getContrastResults(Epis001, 
                                  contrast = contrasto1,
                                  )
getContrastResults()
str(diffEpis001) #data.frame
head(diffEpis001)
print(diffEpis001$adj.P.Val) #there are a lot of intrinsically tied p.adj values


print(correlationepis) #0.997

contrasto2 <- c('Primary.NHU.GR vs Primary.NHU.IgG' 
                = "Primary.NHU.GR - Primary.NHU.IgG")
NHU001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto2)
diffNHU001 <- getContrastResults(NHU001, 
                               contrast = contrasto2)

contrasto3 <- c('PrimaryCD4andT.GR vs PrimaryCD4andT.IgG' 
                = "PrimaryCD4andT.GR - PrimaryCD4andT.IgG")

CD4001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto3)
diffCD4001 <- getContrastResults(CD4001, 
                               contrast = contrasto3)

correlationCD4 <- cor(diffCD4001$adj.P.Val , diffCD4o$adj.P.Val) #0.998

contrasto4 <- c('MCF7.GR vs MCF7.IgG' 
                = "MCF7.GR - MCF7.IgG")

MCF7001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto4)
diffMCF7001 <- getContrastResults(MCF7001, 
                                contrast = contrasto4)

contrasto5 <- c('KMBC2.GR vs KMBC2.IgG' 
                = "KMBC2.GR - KMBC2.IgG")

KMBC2001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto5)
diffKMBC2001 <- getContrastResults(KMBC2001, 
                                 contrast = contrasto5)

contrasto6 <- c('Jurkat.GR vs Jurkat.IgG' 
                = "Jurkat.GR - Jurkat.IgG")

Jurkat001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto6)
diffJurkat001 <- getContrastResults(Jurkat001, 
                                  contrast = contrasto6)

contrasto7 <- c('X231.GR vs X231.IgG' 
                = "X231.GR - X231.IgG")

X231001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto7)
diffX231001 <- getContrastResults(X231001, 
                                contrast = contrasto7)


list001<-list(diffEpis001, diffNHU001, diffCD4001,
            diffMCF7001, diffKMBC2001, diffJurkat001, diffX231001)

merged.diff001<-data.frame()

for (i in 1:length(list001)) {
  # Extract the "a" and "b" columns from the current data frame
  ab_columns <- list001[[i]][, c("log2FC", "adj.P.Val")]
  
  # If merged_df is empty, assign ab_columns to it
  if (nrow(merged.diff001) == 0) {
    merged.diff001 <- ab_columns
  } else {
    # Otherwise, append ab_columns to the existing merged_df
    merged.diff001 <- cbind(merged.diff001, ab_columns)
  }
}

head(merged.diff001)
View(merged.diff001)

colnames(merged.diff001) <- c("log2FC_Epis", "adj.P.Val_Epis", 
                            "log2FC_NHU", "adj.P.Val_NHU", 
                            "log2FC_CD4", "adj.P.Val_CD4",
                            "log2FC_MCF7", "adj.P.Val_MCF7",
                            "log2FC_KMBC2", "adj.P.Val_KMBC2",
                            "log2FC_Jurkat", "adj.P.Val_Jurkat",
                            "log2FC_231", "adj.P.Val_231")
write.csv(merged.diff001, file = "merged.diff001.csv", row.names = TRUE)
# Extract the relevant columns containing adj.P.Val values
adj.P.Val.columns001 <- colnames(merged.diff001)[grep("adj.P.Val", colnames(merged.diff001))]

#start with comparison within GR IP
contrasto8 <- c('X231.GR vs PrimaryBreastEpis.GR' 
                = "X231.GR - PrimaryBreastEpis.GR")

breast001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto8)
diffbreast001 <- getContrastResults(breast001, 
                                  contrast = contrasto8)
View(diffbreast001)
write.csv(diffbreast001, file = "diffbreast001.csv", row.names = TRUE)
head(diffbreast001)

diffbreast001 <- diffbreast001[c("Accessions", "log2FC", "adj.P.Val")]
str(diffbreast001)
View(diffbreast001)
diffbreast001$ranking <- -log10(diffbreast001$adj.P.Val) * diffbreast001$log2FC

filtered_diffbreast001 <- diffbreast001[diffbreast001$adj.P.Val < 0.05, ]
View(filtered_diffbreast001)
write.csv(filtered_diffbreast001, file = "Sigbreast001.csv", row.names = TRUE)

#start with comparison within GR IP Epis vs blood (EpisVSBlood)
contrasto9 <- c('PrimaryBreastEpis.GR vs PrimaryCD4andT.GR'
                = "PrimaryBreastEpis.GR - PrimaryCD4andT.GR")

EpisVSBlood001 <- computeDiffStats(MSnset_gsnorm_001, contrasts = contrasto9)
diffEpisVSBlood001 <- getContrastResults(EpisVSBlood001, 
                                    contrast = contrasto9)
View(diffEpisVSBlood001)
write.csv(diffbreast001, file = "diffbreast001.csv", row.names = TRUE)
head(diffbreast001)

diffEpisVSBlood001 <- diffEpisVSBlood001[c("Accessions", "log2FC", "adj.P.Val")]
#seemingly to many significantly different proteins, 
#due to the intrinsic intensity level
str(diffbreast001)
View(diffbreast001)
diffbreast001$ranking <- -log10(diffbreast001$adj.P.Val) * diffbreast001$log2FC

#drop CD4 001 sample, save GR IPs only 
# Find the column numbers containing "_GR_IP"
GRIPcolumn_numbers <- grep("_GR_IP", colnames(drop001))
print(GRIPcolumn_numbers)
str(GRIPcolumn_numbers)
# Subset the dataframe to retain rows with "GR" in SampleName
E239metaGR001<- E239meta001[grep("GR", E239meta001$SampleName), ]
View(E239metaGR001)
msnGR001<-convertToMSnset(drop001, 
                        metadata= E239metaGR001,
                        indExpData=c(GRIPcolumn_numbers), type='protein',Accessions= 1)

MSnset_normalGR001 <- normalizeScaling(msnGR001, 
                                  scalingFunction = median) 
intensityPlot(MSnset_normalGR001, title = "Scaling across all tissues")
intensityBoxplot(MSnset_normalGR001, title = "Scaling across all tissues")

contrasto8 <- c('X231.GR vs PrimaryBreastEpis.GR' 
                = "X231.GR - PrimaryBreastEpis.GR")

breastGR001 <- computeDiffStats(MSnset_normalGR001, contrasts = contrasto8)
diffbreastGR001 <- getContrastResults(breastGR001, 
                                    contrast = contrasto8)
View(diffbreast001)
write.csv(diffbreast001, file = "diffbreast001.csv", row.names = TRUE)
head(diffbreast001)

diffbreastGR001 <- diffbreastGR001[c("Accessions", "log2FC", "adj.P.Val")]
View(diffbreastGR001)

mergedbreast <- merge(diffbreast001, diffbreastGR001, 
                      by = "row.names", all = TRUE)

correlationP <- cor(-log10(mergedbreast$adj.P.Val.x), -log10(mergedbreast$adj.P.Val.y))
print(correlationP) #0.85084735864282645679
plot(-log10(mergedbreast$adj.P.Val.x), -log10(mergedbreast$adj.P.Val.y), 
     xlab = "-log10 In group norm",  # Label for the x-axis
     ylab = "-log10 Norm across tissues",  # Label for the y-axis
     main = "231 vs breast epis",  # Title for the plot
     col = "blue"  # Point color
)
correlationFC<-cor(mergedbreast$log2FC.x, mergedbreast$log2FC.y)
plot(mergedbreast$log2FC.x, mergedbreast$log2FC.y)
print(correlationFC) #r=0.99999987871169759845

diffbreastGR001$ranking <- -log10(diffbreastGR001$adj.P.Val) * diffbreastGR001$log2FC

filtered_diffbreastGR001 <- diffbreastGR001[diffbreastGR001$adj.P.Val < 0.05, ]
View(filtered_diffbreast001)
write.csv(filtered_diffbreastGR001, 
          file = "SigbreastGR001.csv", row.names = TRUE)
#normalized across GR
contrasto9 <- c('PrimaryBreastEpis.GR vs PrimaryCD4andT.GR'
                = "PrimaryBreastEpis.GR - PrimaryCD4andT.GR")

EpisVSBloodGR001 <- computeDiffStats(MSnset_normalGR001, contrasts = contrasto9)
diffEpisVSBloodGR001 <- getContrastResults(EpisVSBloodGR001, 
                                         contrast = contrasto9)
View(diffEpisVSBloodGR001)


diffEpisVSBloodGR001 <- diffEpisVSBloodGR001[c("Accessions", "log2FC", "adj.P.Val")]
filtered_diffEBGR001 <- diffEpisVSBloodGR001[diffEpisVSBloodGR001$adj.P.Val < 0.05, ]
filtered_diffEBGR001$ranking <- -log10(filtered_diffEBGR001$adj.P.Val) * filtered_diffEBGR001$log2FC

write.csv(filtered_diffEBGR001, 
          file = "SigEBGR001.csv", row.names = TRUE)


mergedEB <- merge(diffEpisVSBlood001, diffEpisVSBloodGR001, 
                      by = "row.names", all = TRUE)
cor(-log10(mergedEB$adj.P.Val.x), -log10(mergedEB$adj.P.Val.y)) #0.59805210818626830527


# Convert the relevant columns to a numeric matrix
getwd()

merged.diff001<-read.csv('merged.diff001.csv',
                     sep=',', stringsAsFactors = F)
heatmap_data001 <- as.matrix(merged.diff001[, adj.P.Val.columns001])
View(heatmap_data001)
# Perform hierarchical clustering (rows and columns)
row_dendrogram001 <- hclust(dist(heatmap_data001))
col_dendrogram001 <- hclust(dist(t(heatmap_data001)))

install.packages("pheatmap")
library("pheatmap")

color_breaks = c(0,0.001,0.01,0.05, 1)

pph<-pheatmap(
  heatmap_data001,
  Colv = as.dendrogram(col_dendrogram001),
  Rowv = as.dendrogram(row_dendrogram001),
  scale = "none",  # You can change scaling options as needed
  col = c("#86ebc9", "#869ceb",
                           "#b986eb","white"),# Choose your color palette
  breaks = color_breaks,
  main = "Heatmap of adj.P.Val (Drop CD4 sample 001",
  xlab = "Samples"
)

heatmap_data001_df <- as.data.frame(heatmap_data001)
conserved <- heatmap_data001_df[apply(heatmap_data001_df < 0.05, 1, all), ]
View(conserved)

# and df_a is the DataFrame to be subset

# Extract the UniProtKB accession numbers from the row names of conserved
accession_numbers <- rownames(conserved)

# Subset df_a based on the accession column
ConservedGR <- unique_xfiltereded[unique_xfiltereded$Accessions %in% accession_numbers, ]


write.csv(ConservedGR, file = "ConservedGR.csv", row.names = FALSE)

# Extract the relevant columns containing log2FC
log2FC.columns001 <- colnames(merged.diff001)[grep("log2FC", colnames(merged.diff001))]

dfFC001<-merged.diff001[, log2FC.columns001]
View(dfFC001)

columns_from_df_a <- unique_xfiltereded[, c(1, 4)]

# Get the row names of df_b
row_names_from_df_b <- rownames(dfFC001)

# Create a new data frame by combining columns from df_a with all columns from df_b
combined_df <- data.frame(Names = row_names_from_df_b, columns_from_df_a, dfFC001)
View(combined_df)
ConservedGRFC <- combined_df[combined_df$Accessions %in% accession_numbers, ]
View(ConservedGRFC)

ConservedGRFC0 <- ConservedGRFC[apply(ConservedGRFC <= 0.00, 1, all), ]

View(ConservedGRFC0)

#pre-filter log2FC less than 2 proteins, the union of FC>4 from each sample type
head(merged.diff001)
boxplot<-boxplot(merged.diff001[, c(2,4,6,8,10,12,14)],
        main="In group Normalized", ylab="log2FC")

# Define the condition column names
condition_columns <- c("log2FC_Epis", "log2FC_NHU", 
                       "log2FC_CD4", "log2FC_MCF7", 
                       "log2FC_KMBC2", "log2FC_Jurkat", "log2FC_231")

# Create an empty set to store unique rownames
unique_rownames <- character(0)

# Loop through each condition and collect unique rownames
for (col in condition_columns) {
  # Filter rows where log2FC is greater than or equal to 2
  condition <- merged.diff001[merged.diff001[[col]] >= 2, ]
  # Get the unique rownames from this condition
  unique_rownames_condition <- rownames(condition)
  # Update the set of unique rownames
  unique_rownames <- union(unique_rownames, unique_rownames_condition)
}

# Subset the original dataframe based on the unique rownames
diff001fc2 <- merged.diff001[unique_rownames, ]
View(diff001fc2)
rownames(diff001fc2) <- diff001fc2[, 1]
str(diff001fc2)
boxplot(diff001fc2[, c(2,4,6,8,10,12,14)],
        main="In group Normalize", ylab="log2FC")

# Extract the relevant columns containing adj.P.Val values
adj.P.Val.columns001fc2 <- colnames(diff001fc2)[grep("adj.P.Val", 
                                                     colnames(diff001fc2))]

# Convert the relevant columns to a numeric matrix
heatmap_datafc2 <- as.matrix(diff001fc2[, adj.P.Val.columns001fc2])
View(heatmap_datafc2)
str(heatmap_datafc2)
# Perform hierarchical clustering (rows and columns)
row_dendrogramfc2 <- hclust(dist(heatmap_datafc2))
col_dendrogramfc2 <- hclust(dist(t(heatmap_datafc2)))
str(row_dendrogramfc2)
plot(row_dendrogramfc2, main = "Row Dendrogram")

install.packages("pheatmap")
library("pheatmap")

color_breaks = c(0,0.001,0.01,0.05, 1)

phmfc2<-pheatmap(
  heatmap_datafc2,
  Colv = as.dendrogram(col_dendrogramfc2),
  Rowv = as.dendrogram(row_dendrogramfc2),
  scale = "none",  # You can change scaling options as needed
  col = c("#86ebc9", "#869ceb",
          "grey","white"),# Choose your color palette
  breaks = color_breaks,
  main = "Heatmap of adj.P.Val log2fC>2",
  xlab = "Samples",
  treeheight_row = 200,
  show_rownames = FALSE,
  legend = FALSE)

#070724 reproduce
pheatmap(heatmap_datafc2,
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  color = c("#86ebc9", "#869ceb",
            "grey","white"),
  breaks = color_breaks,
  main = "DIA-NN GR adj.P log2fC>2",
  display_numbers = FALSE,  # Remove row names
  border_color = NA,  # Remove borders around cells
  cellwidth = 40,  # Automatic width for cells
  cellheight = NA,  # Automatic height for cells
  treeheight_row = 100,  # Adjust the dendrogram height for rows
  treeheight_col = 50,  # Adjust the dendrogram height for columns
  fontsize_row = 10,  # Adjust row font size (if row names were to be displayed)
  fontsize_col = 10,  # Adjust column font size
  fontsize = 12,  # Adjust overall font size
  clustering_distance_rows = "euclidean",  # Distance metric for rows
  clustering_distance_cols = "euclidean",  # Distance metric for columns
  clustering_method = "complete",# Clustering method
  show_rownames = FALSE,
)

color = heatmap_colors,
breaks = color_breaks,

pheatmap(
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



install.packages('dendextend')
install.packages('gplots')
install.packages('viridis')
install.packages("devtools")
library('dendextend')
library('gplots')
library('viridis')
library("devtools")
den<-as.dendrogram(row_dendrogramfc2)
dend_list <- get_subdendrograms(den, 2)
# Plotting the result
par(mfrow = c(2,3))
par(mar = c(1, 1, 1, 1)) 
plot(den, main = "Original dendrogram")
sapply(dend_list, plot)


cut_tree <- dendextend::cutree(
  den,
  h = 1,
  order_clusters_as_data = FALSE
)
table(cut_tree)
# Then, we plot the tree with group labels
install.packages("pals")
library("pals")
den %>%
  dendextend::color_branches(.,
                             h = 0.1,
                             col = pals::alphabet2(),
                             groupLabels = FALSE) %>%
  plot

cut_tree0.1 <- dendextend::cutree(
  den,
  h = 0.1,
  order_clusters_as_data = FALSE
)

table(cut_tree0.1)

desired_branch <- dendextend::cutree(
  den,
  k = 35,
  order_clusters_as_data = FALSE)

table(desired_branch)
str(desired_branch)
desired_branch[desired_branch == 19]

# Extract names that match the condition
block11 <- names(desired_branch[desired_branch == 11])

# Convert the matching names to a character vector
block11vector <- as.character(block11)
print(block11vector)
block11vector <- gsub('"', '', block11vector)
block11df <- data.frame(Column1 = block11vector)
View(block11df)
# Export the data frame to a CSV file
write.csv(block11df, "block11.csv", row.names = FALSE)

block12 <- names(desired_branch[desired_branch == 12])
block12vector <- as.character(block12)
block12df <- data.frame(Column1 = block12vector)
View(block12df)
write.csv(block12df, "block12.csv", row.names = FALSE)

block12 <- names(desired_branch[desired_branch == 12])
block12vector <- as.character(block12)
block12df <- data.frame(Column1 = block12vector)
View(block12df)
write.csv(block12df, "block12.csv", row.names = FALSE)

# Loop through values 12 to 16
for (value in 8:10) {
  # Create a vector with elements equal to 'value'
  block_vector <- names(desired_branch[desired_branch == value])
  block_vector <- as.character(block_vector)
  
  # Create a data frame with one column
  block_df <- data.frame(Column1 = block_vector)
  
  # Create a CSV file with a filename that includes the value
  filename <- paste("block", value, ".csv", sep = "")
  write.csv(block_df, filename, row.names = FALSE)
  
  cat("CSV file", filename, "created.\n")
}

desired_branch4 <- dendextend::cutree(
  den,
  k = 4,
  order_clusters_as_data = FALSE)

table(desired_branch4)

blockA <- names(desired_branch4[desired_branch4 == 1])
blockAvector <- as.character(blockA)
blockAdf <- data.frame(Column1 = blockAvector)
View(blockAdf)
write.csv(blockAdf, "blockA.csv", row.names = FALSE)

#play with log2fc
boxplot(diff001fc2[, c(2,4,6,8,10,12,14)],
        main="In group Normalize", ylab="log2FC")

# Convert the relevant columns to a numeric matrix
hmfc2 <- as.matrix(diff001fc2[, c(2,4,6,8,10,12,14)])
View(hmfc2)

# Perform hierarchical clustering (rows and columns)
row_dendrogramhmfc2.1 <- hclust(dist(hmfc2^2))
col_dendrogramhmfc2.1 <- hclust(dist(t(hmfc2^2)))
install.packages("pheatmap")
library("pheatmap")

color_b = c(0,2,10,20,40,60,80,100,250,10000)

# Reset mfrow to 1x1 (one plot per device)
par(mfrow = c(1, 1))



hmfc3<-pheatmap(
  hmfc2,
  Colv = as.dendrogram(col_dendrogramhmfc2),
  Rowv = as.dendrogram(row_dendrogramhmfc2),
  scale = "none",  # You can change scaling options as needed
  col = c("black", "grey", "#F4A582", "red", 
          "green", "#92C5DE", "#4393C3", "#2166AC"),# Choose your color palette
  breaks = color_b,
  main = "Heatmap of log2 FC",
  xlab = "Samples",
  treeheight_row =400,
  show_rownames = FALSE,
  legend = TRUE)

hmfc3.1<-pheatmap(
  hmfc2^2,
  Colv = as.dendrogram(col_dendrogramhmfc2.1),
  Rowv = as.dendrogram(row_dendrogramhmfc2.1),
  scale = "none",  # You can change scaling options as needed
  col = c("black", "grey", "#F4A582", "red", 
          "green", "#92C5DE", "#4393C3", "#2166AC","blue"),# Choose your color palette
  breaks = color_b,
  main = "Heatmap of liner FC",
  xlab = "Samples",
  treeheight_row =400,
  show_rownames = FALSE,
  legend = FALSE)


denfc<-as.dendrogram(row_dendrogramhmfc2.1)
cut_treeFC <- dendextend::cutree(
  denfc,
  h = 100,
  order_clusters_as_data = FALSE
)
table(cut_treeFC)

desired_branchFC <- dendextend::cutree(
  denfc,
  k = 16, #k is the number of total clusters
  order_clusters_as_data = FALSE)

table(desired_branchFC)


blockFCB <- names(desired_branchFC[desired_branchFC == 2])
blockFCBvector <- as.character(blockFCB)
blockFCBdf <- data.frame(Column1 = blockFCBvector)
View(blockFCBdf)
write.csv(blockFCBdf, "blockFCB.csv", row.names = FALSE)

# Loop through values 1 to 4
for (value in 1:4) {
  # Create a vector with elements equal to 'value'
  block_vector <- names(desired_branchFC[desired_branchFC == value])
  block_vector <- as.character(block_vector)
  
  # Create a data frame with one column
  block_df <- data.frame(Column1 = block_vector)
  
  # Create a CSV file with a filename that includes the value
  filename <- paste("blockFC", value, ".csv", sep = "")
  write.csv(block_df, filename, row.names = FALSE)
  
  cat("CSV file", filename, "created.\n")
}

E239_FragPipe_DE<-read.csv('E239_FragPipe_DE.csv',sep=',', stringsAsFactors = F)
head(E239_FragPipe_DE)

# Get the column names that match your criteria
E239_col<- grep("^([^_]+)_GR_vs_\\1_IgG_", colnames(E239_FragPipe_DE))
str(E239_col)
# Subset the dataframe using the selected column names
E239_coll <- c(1, 2, 54, 145, 236, 327, E239_col)
print(E239_coll)
str(E239_coll)
subset_E239_FragPipe_DE <- E239_FragPipe_DE[, E239_coll]
head(subset_E239_FragPipe_DE)
View(subset_E239_FragPipe_DE)


# Get the column names that contain "log2.fold.change"
log2_fold_change_cols <- grep("log2.fold.change", colnames(subset_E239_FragPipe_DE), value = TRUE)
print(log2_fold_change_cols)
# Subset the dataframe
E239_FragPipe_FC <- subset_E239_FragPipe_DE[, c("Protein.ID", "Gene.Name", log2_fold_change_cols)]
View(E239_FragPipe_FC)

# Convert the relevant columns to a numeric matrix
heatmap_dataFC_FragPipe <- as.matrix(E239_FragPipe_FC[, log2_fold_change_cols])

# Perform hierarchical clustering (rows and columns)
row_dendrogramFC_FragPipe <- hclust(dist(heatmap_dataFC_FragPipe))
col_dendrogramFC_FragPipe <- hclust(dist(t(heatmap_dataFC_FragPipe)))

heatmap(
  heatmap_dataFC_FragPipe,
  Colv = as.dendrogram(col_dendrogramFC_FragPipe),
  Rowv = as.dendrogram(row_dendrogramFC_FragPipe),
  scale = "none",  # You can change scaling options as needed
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Choose your color palette
  main = "Heatmap of log2FC FragPipe",
  xlab = "Samples",
  ylab = "Genes",
  margins = c(5, 5)
)

# Get the column names that contain "log2.fold.change"
p.adjcols <- grep("p.adj", colnames(subset_E239_FragPipe_DE), value = TRUE)
print(p.adjcols)
# Subset the dataframe
E239_FragPipe_p <- subset_E239_FragPipe_DE[, c("Protein.ID", "Gene.Name", p.adjcols)]
View(E239_FragPipe_p)

# Convert the relevant columns to a numeric matrix
heatmap_datap_FragPipep <- as.matrix(E239_FragPipe_p[, p.adjcols])

# Perform hierarchical clustering (rows and columns)
row_dendrogramp_FragPipe <- hclust(dist(heatmap_datap_FragPipep))
col_dendrogramp_FragPipe <- hclust(dist(t(heatmap_datap_FragPipep)))

heatmap(
  heatmap_datap_FragPipep,
  Colv = as.dendrogram(col_dendrogramp_FragPipe),
  Rowv = as.dendrogram(row_dendrogramp_FragPipe),
  scale = "none",  # You can change scaling options as needed
  col = colorRampPalette(c("blue", "white", "red"))(100),  # Choose your color palette
  main = "Heatmap of p.adj FragPipe",
  xlab = "Samples",
  ylab = "Genes",
  margins = c(5, 5)
)

msnX<-convertToMSnset(unique_xfiltereded, 
                      metadata= E239meta,
                      indExpData=c(9:52), type='protein',Accessions= 1)

intensityPlot(msnX, title = "No normalization")
p1
intensityBoxplot(msnX, title = "Peptide intensity distribution")
rliPlot(msnX, title = "Relative Peptide intensity")
hierarchicalPlot(msnX)

#groupScaling() Performs scaling normalization on the intensities within group (median or mean)
MSnset_norm_gs <- groupScaling(msnX, 
                               scalingFunction = median, 
                               groupingColumn = "SampleGroup") #IgG+ each tissue GR-IP
View(MSnset_norm_gs)

intensityPlot(MSnset_norm_gs, title = "Within Group Scaling")
intensityBoxplot(MSnset_norm_gs, title = "Within Group Scaling")
hierarchicalPlot(MSnset_norm_gs, 
                 colourBy='SampleGroup',title = "Hierarchical Clustering, Groups Scaled")
?computeDiffStats()

contrast <- c('PrimaryBreastEpis.GR vs MCF7.GR' = "PrimaryBreastEpis.GR - MCF7.GR")
contrast2 <- c('PrimaryBreastEpis.GR vs IgG.IP' =
                "PrimaryBreastEpis.GR - IgG.IP")

install.packages("statmod")
library("statmod")

diffstats2 <- computeDiffStats(MSnset_norm_gs, contrasts = contrast2)
diffexp2 <- getContrastResults(diffstats2, 
                              contrast = contrast2
                              )
maVolPlot(diffstats, contrast = contrast, plotType = "MA"
          )

maVolPlot(diffstats, contrast = contrast, plotType = "Volcano")


write.csv(diffexp, file = "diffexp.csv", row.names = FALSE)
