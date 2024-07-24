#ADD CD4+T cell RIME data in (18/01/2023)
getwd()
setwd("/Users/weiye/Downloads")
E033<-read.csv('Protein Report for E033_GR_RIME.csv',sep=',', stringsAsFactors = F,
               na.strings=c(""," ","NA"))
#remove headers and rename column names
E033 <- E033[-(1:43),]
colnames(E033) <- c(E033[1,])
E033 <- E033[-(1),]
View(E033)
#drop NA Database sources(exclude non-human proteins)
E033<-E033[!is.na(E033$`Database sources`), ]
#Split clustered proteins into indiviuals
IDs<-strsplit(as.character(E033$`Alternate IDs`),',',
              fixed = TRUE)

E033protein<-data.frame(protein=unlist(IDs),
                        ProteinCount=rep(E033$`Total spectrum count`,
                                     sapply(IDs,FUN = length)),
                        Biological=rep(E033$`Biological sample name`,
                                       sapply(IDs,FUN = length)))

E033protein <- E033protein[!is.na(E033protein$protein), ]

View(E033protein)
str(E033protein)
tidy_E033protein <- pivot_wider(E033protein, 
                                names_from = Biological,
                                values_from = ProteinCount)
View(tidy_E033protein)
str(tidy_E033protein)
tidy_E033protein <- as.data.frame(tidy_E033protein)
View(tidy_E033protein)
str(tidy_E033protein)

tidy_E033protein <- tidy_E033protein %>%
  mutate_at(vars(GR1:IgG3), as.numeric)

columns_to_replace <- c("GR1", "GR2", "GR3", "IgG1", "IgG2", "IgG3")

tidy_E033protein <- tidy_E033protein %>%
  mutate_at(vars(columns_to_replace), ~replace(as.numeric(.), is.na(.), 0))
View(tidy_E033protein)
str(tidy_E033protein)

E033peptide<-data.frame(protein=unlist(IDs),
                  Pepcount=rep(E033$`Exclusive unique peptide count`,
                          sapply(IDs,FUN = length)),
                  Biological=rep(E033$`Biological sample name`,
                                 sapply(IDs,FUN = length)))
View(E033peptide)
E033peptide<-E033peptide[!is.na(E033peptide$protein), ]
head(E033peptide)
biological_cou <- table(E033peptide$Biological)
print(biological_cou)

library(tidyr)
tidy_E033peptide <- pivot_wider(E033peptide, 
                                names_from = Biological,
                                values_from = Pepcount)
str(tidy_E033peptide)
unique(E033peptide$Biological)
tidy_E033peptide <- tidy_E033peptide[, c("protein", "GR1", "GR2", "GR3",
                       "IgG1", "IgG2", "IgG3")]


tidy_E033peptide<-as.data.frame(tidy_E033peptide)
tidy_E033peptide[is.na(tidy_E033peptide)] <- 0
# Assuming your dataframe is named 'tidy_E033peptide'
tidy_E033peptide[, 2:7] <- lapply(tidy_E033peptide[, 2:7], as.numeric)
View(tidy_E033peptide)
str(tidy_E033peptide)

tidy_E033peptide$Peptide <- pmax(tidy_E033peptide$GR1, 
                                 tidy_E033peptide$GR2, 
                                 tidy_E033peptide$GR3)

E033DIApeptide<-read.csv('E033_DIA_Report_Peptide_Count.csv',sep=',', stringsAsFactors = F,
               na.strings=c(""," ","NA"))
View(E033DIApeptide)
DIAIDs<-strsplit(as.character(E033DIApeptide$Genes),';',
              fixed = TRUE)
tidy_E033DIApeptide<-data.frame(protein=unlist(DIAIDs),
                        Pepcount=rep(E033DIApeptide$Total_Peptides,
                                     sapply(DIAIDs,FUN = length))
                        )
View(tidy_E033DIApeptide)
str(tidy_E033DIApeptide)

#count number of variables
sum(!is.na(E033DIApeptide$ER_1_DIA.Unique.count.Peptides))
sum(!is.na(E033DIApeptide$ER_2_DIA.Unique.count.Peptides.))
sum(!is.na(E033DIApeptide$ER_3_DIA.Unique.count.Peptides.))

biological_counts <- table(E033DIApeptide$ER_1_DIA.Unique.count.Peptides.)
print(biological_counts)
library(ggplot2)

mergedE033peptide <- merge(tidy_E033peptide, tidy_E033DIApeptide, by = "protein", all = TRUE)
mergedE033peptide[is.na(mergedE033peptide)] <- 0
# Directly assign new names to the last two columns
colnames(mergedE033peptide)[c(8, 9)] <- c("DDA", "DIA")
# Calculate the rank for the "DDA" and "DIA" columns
mergedE033peptide$DDA_Rank <- rank(mergedE033peptide$DDA)
mergedE033peptide$DIA_Rank <- rank(mergedE033peptide$DIA)
str(mergedE033peptide)
View(mergedE033peptide)
plot(mergedE033peptide$DDA,
     mergedE033peptide$DIA, xlab = "DDA peptide count", ylab = "DIA peptide count")

plot(mergedE033peptide$DDA_Rank,
     mergedE033peptide$DIA_Rank, xlab = "DDA Ranks", ylab = "DIA Ranks")
library(dplyr)
CommonE033peptide <- inner_join(tidy_E033peptide, tidy_E033DIApeptide, by = "protein")
CommonE033peptide[is.na(CommonE033peptide)] <- 0
colnames(CommonE033peptide)[c(8, 9)] <- c("DDA", "DIA")
CommonE033peptide$DDA_Rank <- rank(CommonE033peptide$DDA)
CommonE033peptide$DIA_Rank <- rank(CommonE033peptide$DIA)

View(CommonE033peptide)



install.packages("ggpubr")
library("ggpubr")
ggpaired(mergedE033peptide, cond1 = "DDA_Rank", cond2 = "DIA_Rank",
         fill = "condition", palette = "jco",
         title = "peptide count (Union proteins)")+ ylab("Rank")

ggpaired(CommonE033peptide, cond1 = "DDA_Rank", cond2 = "DIA_Rank",
         fill = "condition", palette = "jco",
         title = "peptide count (Common proteins)") + ylab("Rank")

biological_counts <- table(E0333$Biological)
print(biological_counts)

# Counts of detected proteins
dda_counts <- c(236, 278, 244)
dia_counts <- c(331, 331, 331)

# Create a 2x3 contingency table
data <- matrix(c(dda_counts, dia_counts), nrow = 2, byrow = TRUE)
colnames(data) <- c("Replicate 1", "Replicate 2", "Replicate 3")
rownames(data) <- c("DDA", "DIA")

# Perform a fisher exact test
resultfisher <- fisher.test(data)
print(resultfisher) #0.336

#peptide counts 
dda_counts <- c(236, 278, 244)
dia_counts <- c(420, 466, 434)
# p-value = 0.8234 Fisher's Exact Test for Count Data

# Filter the dataframe to select rows with Biological category starting with "GR"
subset_df <- E0333[grep("^GR", E0333$Biological), ]
# Count the unique proteins in the subset
unique_proteins <- unique(subset_df$protein)
print(unique_proteins)
# Count the number of unique proteins
num_unique_proteins <- length(unique_proteins)
# Print the number of unique proteins
print(num_unique_proteins) #295

E033DIA<-read.delim('E033_DIA_NN_Stripped.pg_matrix.tsv',sep='\t', 
                    stringsAsFactors = F)
intersection <- intersect(unique_proteins, E033DIA$Genes)
length(intersection)
unions<- union(unique_proteins, E033DIA$Genes)
length(unions)
install.packages("VennDiagram")   # Install & load VennDiagram package
library("VennDiagram")

venn.diagram(list(B = 0:331, A = 141:436), 
             fill = c("#FF6666", "lightblue"), 
             alpha = c(0.5, 0.5), lwd =0, 
             filename ="venn_diagram4.jpeg")

View(E033DIA)
MetaE033DIA<-read.delim('E033_experiment_annotation (2).tsv',sep='\t', 
                    stringsAsFactors = F)
View(MetaE033DIA)
head(MetaE033DIA)
# Rename the columns
colnames(MetaE033DIA) <- c("TechRep", "SampleName","sample_name",
                           "SampleGroup", "BioRep")

#perform DE analysis with qPLEX analyzer on DIA E033
# Verify the new column names
colnames(MetaE033DIA)

head(E033DIA)
colnames(E033DIA)[6:11] <- MetaE033DIA$SampleName
View(E033DIA)
head()

msnE033DIA<-convertToMSnset(E033DIA, 
                       metadata= MetaE033DIA,
                       indExpData=c(6:11), type='protein',Accessions= 4)
intensityBoxplot(msnE033DIA, title = "No normalization")

msnE033DIA_norm <- groupScaling(msnE033DIA, 
                                scalingFunction = median, 
                                groupingColumn = "SampleGroup")

intensityBoxplot(msnE033DIA_norm, title = "In group")

conE033DIA <- c('GR vs IgG' 
                = "GR - IgG")
E033q <- computeDiffStats(msnE033DIA_norm, contrasts = conE033DIA)
diffE033 <- getContrastResults(E033q, 
                                  contrast = conE033DIA)
View(diffE033)#DE results for DIA
# Check for duplicated values in the specified column
tied_values <-diffE033$adj.P.Val[duplicated(diffE033$adj.P.Val) | duplicated(diffE033$adj.P.Val, fromLast = TRUE)]
length(unique(tied_values)) #67
length(unique(diffE033$adj.P.Val)) #211
# Print the result
print(tied_values)


#working on DDA total spectrum counts qPLEX analyzer, df(tidy_E033protein)
View(tidy_E033protein)
str(tidy_E033protein)

colnames(tidy_E033protein)[colnames(tidy_E033protein) == 'protein'] <- 'Accessions'
# Add 0.0001 to columns 2 through 7, getting rid of zero
tidy_E033protein[, 2:7] <- tidy_E033protein[, 2:7] + 0.0001

msnE033DDA<-convertToMSnset(tidy_E033protein, 
                            metadata= DDAmeta,
                            indExpData=c(2:7), type='protein',Accessions= 1)

DDAmeta <- MetaE033DIA
View(DDAmeta)
msnE033DDA_norm <- groupScaling(msnE033DDA, 
                                scalingFunction = median, 
                                groupingColumn = "SampleGroup")
conE033DIA <- c('GR vs IgG' 
                = "GR - IgG")
E033DDAq <- computeDiffStats(msnE033DDA_norm, contrasts = conE033DIA)
diffE033DDA <- getContrastResults(E033DDAq, 
                               contrast = conE033DIA)
View(diffE033DDA)
str(diffE033DDA)
tied_values <-diffE033$adj.P.Val[duplicated(diffE033$adj.P.Val) | duplicated(diffE033$adj.P.Val, fromLast = TRUE)]
length(unique(diffE033DDA$adj.P.Val))
options(digits = 20)

CommonDEres <- inner_join(diffE033DDA, diffE033, by = "Accessions")
View(CommonDEres)

cor(-log10(CommonDEres$adj.P.Val.x), -log10(CommonDEres$adj.P.Val.y))
plot(-log10(CommonDEres$adj.P.Val.x), 
     -log10(CommonDEres$adj.P.Val.y), 
     main = "Common proteins DDA vs DIA",
     xlab = "-log10 p.adj DDA", 
     ylab = "-log10 p.adj DIA", pch = 16, col = "blue")


diffE033$"-log10P.adj"<--log10(diffE033$adj.P.Val)
str(diffE033) #data.frame
install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
EnhancedVolcano(diffE033,
                lab = rownames(diffE033),
                x = 'log2FC',
                y = 'adj.P.Val',
                pCutoff = 0.01,
                FCcutoff = 2,
                ylim = c(0, 25),
                xlim = c(0, 8),
                pointSize = 0.5,
                labSize = 5,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                boxedLabels = TRUE,
                selectLab = c("C3",
                              "HMGB1",
                              "MCM5",
                              "LRPPRC",
                              "STAT5B",
                              "HDAC2",
                              "SMARCC1",
                              "NCOA3")
) + geom_text_repel(aes(label = c("C3",
                                  "HMGB1",
                                  "MCM5",
                                  "LRPPRC",
                                  "STAT5B",
                                  "HDAC2",
                                  "SMARCC1",
                                  "NCOA3")))
library(ggplot2)
library(ggrepel)

# Your code for creating the EnhancedVolcano plot
EnhancedVolcano(
  diffE033,
  lab = rownames(diffE033),
  x = 'log2FC',
  y = 'adj.P.Val',
  pCutoff = 0.01,
  FCcutoff = 2,
  ylim = c(0, 25),
  xlim = c(0, 8),
  pointSize = 0.5,
  labSize = 5,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  selectLab = c("C3", "HMGB1", "MCM5", "LRPPRC", "STAT5B", "HDAC2", "SMARCC1", "NCOA3")) +
  geom_text_repel(aes(label = ifelse(rownames(diffE033) %in% 
      c("C3", "HMGB1", "MCM5", "LRPPRC", "STAT5B", "HDAC2", "SMARCC1", "NCOA3"), 
                                     rownames(diffE033), "")))
                  
dev.new(width = 25, height = 25)
EnhancedVolcano()