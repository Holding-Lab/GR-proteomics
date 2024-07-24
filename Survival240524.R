#240524 reproduce 
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(stringr)
library(survival)
setwd("/Users/weiye/Downloads/all_phase2_target_2018_pub")
#clinical metadata data
TargetClinical<-read.csv('~/Downloads/all_phase2_target_2018_pub/data_clinical_sample.txt',
                         sep='\t', stringsAsFactors = F)
View(TargetClinical)
#focus on "T cell ALL"
TALLClinical<-TargetClinical[TargetClinical$Cell.of.tumor.origin == 'T Cell ALL',]
View(TALLClinical)

clinical_sample<-read_delim("data_clinical_sample.txt")
clinical_sample <- as.data.frame(clinical_sample)
clinical_patient<-read_delim("data_clinical_patient.txt")
clinical_patient <- as.data.frame(clinical_patient)

ALL_clinical <- merge(clinical_patient, clinical_sample, by = "#Patient Identifier", 
                      all.y = TRUE)
TALL_clinical<-ALL_clinical[ALL_clinical$`Cell of tumor origin`== 
                              "T Cell ALL",]
unique(TALL_clinical$`Cell of tumor origin`)
TALL_clinical<-
  TALL_clinical[!is.na(TALL_clinical$`Cell of tumor origin`), ]



#St Jude ALL RNA-seq data
getwd()
StJude_dir <- "StJude (All patient ALL)"
txt_files<-fs::dir_ls(StJude_dir)
txt_files
readr::read_delim(txt_files[1]) #read the first file in StJude_dir (directory)

StJude_expression <- txt_files %>% 
  map_dfr(read_delim, .id = "source") #works well
StJude_expression<-as.data.frame(StJude_expression)
View(StJude_expression)
StJude_fpkm<-StJude_expression[, c(1,2,3)]
StJude_norm<-StJude_expression[, c(1,2,4)]

StJude_fpkm <- StJude_fpkm %>%
  pivot_wider(names_from = gene, values_from = fpkm)
View(StJude_fpkm)
str(StJude_fpkm) #tibble [264 × 19,465] (S3: tbl_df/tbl/data.frame)

StJude_norm <- StJude_norm %>%
  pivot_wider(names_from = gene, values_from = fpkm_normalized)

IDDD<-str_extract(StJude_norm$source, "TARGET-10-......") #worked! 

StJude_norm$IDDD <- IDDD
StJude_norm<-as.data.frame(StJude_norm)

norm <- merge(StJude_norm, TALL_clinical, by.x = "IDDD", 
              by.y = "#Patient Identifier", all.x = TRUE, all.y = FALSE)

norm[norm$`Overall Survival Status` == "0:LIVING" & !is.na(norm$`Overall Survival Status`), ]$`Overall Survival Status` <- "0"
norm[norm$`Overall Survival Status` == "1:DECEASED" & !is.na(norm$`Overall Survival Status`), ]$`Overall Survival Status` <- "1"

norm$`Overall Survival Status` <- as.numeric(norm$`Overall Survival Status`)
norm <- norm[!duplicated(norm$IDDD), ]
rownames(norm) <- norm$IDDD
norm$`Overall Survival Days` <- as.numeric(norm$`Overall Survival Days`)
str(norm$`Overall Survival Days`)
str(norm$`Overall Survival Status`)
View(norm)
dim(norm) #264 19516
tail(colnames(norm), 55)
head(colnames(norm), 5)
head(norm$"Time To Event (days)")
head(norm$"Overall Survival Days")

summary(norm$`Overall Survival Days`)
summary(norm$`Overall Survival Status`)

# Remove rows with NA values in 'Overall Survival Days' or 'Overall Survival Status'
norm <- norm[complete.cases(norm$`Overall Survival Days`, norm$`Overall Survival Status`), ]

# Verify that the rows with NAs have been removed
summary(norm$`Overall Survival Days`)
summary(norm$`Overall Survival Status`)

# List of genes (from column 3 to 19465)
genes <- colnames(norm)[3:19465]
View(norm)
# Create a list to store results
results <- list()

# Loop through each gene
for (gene in genes) {
  # Create a formula for the Cox model
  formula <- as.formula(paste("Surv(`Overall Survival Days`, `Overall Survival Status`) ~ `", gene, "`", sep = ""))
  
  # Fit the Cox proportional hazards model
  cox_model <- coxph(formula, data = norm)
  
  # Summarize the model
  summary_cox <- summary(cox_model)
  
  # Store the results in the list
  results[[gene]] <- summary_cox
}

# Function to interpret the coefficient and hazard ratio
interpret_hazard <- function(coef, exp_coef) {
  if (coef > 0) {
    interpretation <- "High expression is associated with increased hazard (worse prognosis)."
  } else if (coef < 0) {
    interpretation <- "High expression is associated with decreased hazard (better prognosis)."
  } else {
    interpretation <- "Expression level is not associated with hazard."
  }
  return(interpretation)
}

# Print the summary and interpretation of the Cox model for each gene
for (gene in genes) {
  cat("Gene:", gene, "\n")
  print(results[[gene]])
  
  coef <- results[[gene]]$coefficients[1, "coef"]
  exp_coef <- results[[gene]]$coefficients[1, "exp(coef)"]
  interpretation <- interpret_hazard(coef, exp_coef)
  
  cat("Interpretation:", interpretation, "\n\n")
}


#26/06/2024 reproduce 5 cancer types: T-ALL, brca, pan-BC, ER+, TNBC 
setwd("~/Downloads/blca_msk_tcga_2020")

getwd() #"/Users/weiye/Downloads/blca_msk_tcga_2020"

library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(stringr)
library(survival)

#blca

blca_clinical<-read.csv('data_clinical_patient.txt',
                        sep='\t', stringsAsFactors = F)

# Remove the first 4 rows from the data frame blca_clinical
blca_clinical <- blca_clinical[-c(1:4), ]

View(blca_clinical)

blca_clinical$Overall.Survival..Months. <- as.numeric(blca_clinical$Overall.Survival..Months.)

blca_clinical <- blca_clinical[complete.cases(blca_clinical$Overall.Survival..Months., 
                                              blca_clinical$Overall.Survival.Status), ]

blca_expression<-read.csv('data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt',
                        sep='\t', stringsAsFactors = F)

View(blca_expression)
blca_expression$Hugo_Symbol
# Check for duplicated values in the Hugo_Symbol column
duplicates <- blca_expression$Hugo_Symbol[duplicated(blca_expression$Hugo_Symbol)]
print(duplicates)
# Remove duplicated rows based on the Hugo_Symbol column
blca_expression_unique <- blca_expression[!duplicated(blca_expression$Hugo_Symbol), ]
View(blca_expression_unique)
# Display the data frame without duplicated rows
print(blca_expression_unique)

rownames(blca_expression_unique) <- blca_expression_unique$Hugo_Symbol
View(blca_expression_unique)

blca_expression_unique <- blca_expression_unique[ ,-c(1,2)]

tblca_expression_unique <- t(blca_expression_unique)
View(tblca_expression_unique)

transform_id <- function(id) {
  # Replace periods with hyphens
  id <- gsub("\\.", "-", id)
  # Remove the last section after the third hyphen
  id <- sub("-[0-9]+$", "", id)
  return(id)
}

rownames(tblca_expression_unique) <- sapply(rownames(tblca_expression_unique), transform_id)
tblca_expression_unique <- as.data.frame(tblca_expression_unique)
tblca_expression_unique$TCGA_Id <- rownames(tblca_expression_unique)
dim(tblca_expression_unique) #296 20460
View(tblca_expression_unique)
blca_exp_meta <- merge(tblca_expression_unique, blca_clinical, by.x = "TCGA_Id", by.y = "X.Patient.Identifier",
                       all = FALSE)

 
View(blca_exp_meta)

blca_exp_meta[blca_exp_meta$Overall.Survival.Status == "0:LIVING" & !is.na(blca_exp_meta$Overall.Survival.Status), ]$Overall.Survival.Status<- "0"
blca_exp_meta[blca_exp_meta$Overall.Survival.Status == "1:DECEASED" & !is.na(blca_exp_meta$Overall.Survival.Status), ]$Overall.Survival.Status<- "1"
blca_exp_meta$Overall.Survival.Status <- as.numeric(blca_exp_meta$Overall.Survival.Status)
print(blca_exp_meta$Overall.Survival.Status)

View(blca_exp_meta)



# Define the genes of interest
blca_sur_genes <- colnames(blca_exp_meta)[2:20461]

# Initialize list to store results
results <- list()
blca_results_df <- data.frame(Gene = character(), HR = numeric(), CI_Lower = numeric(), CI_Upper = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene and perform survival analysis
for (gene in blca_sur_genes) {
  # Ensure there are no missing values for the current gene
  gene_data <- blca_exp_meta[[gene]]
  if (all(!is.na(gene_data))) {
    # Create the survival formula dynamically
    # Use backticks to handle gene names with special characters
    formula <- as.formula(paste("Surv(Overall.Survival..Months., Overall.Survival.Status) ~ `", gene, "`", sep = ""))
    
    # Fit the Cox proportional hazards model
    cox_model <- tryCatch({
      coxph(formula, data = blca_exp_meta)
    }, error = function(e) {
      cat("Error in model fitting for gene:", gene, "\n")
      NULL
    })
    
    if (!is.null(cox_model)) {
      summary_cox <- summary(cox_model)
      
      # Extract the hazard ratio, confidence intervals, and p-value
      hr <- summary_cox$coefficients[1, "exp(coef)"]
      ci_lower <- summary_cox$conf.int[1, 1]
      ci_upper <- summary_cox$conf.int[1, 3]
      p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
      
      # Store the results in a data frame
      blca_results_df <- rbind(blca_results_df, data.frame(Gene = gene, HR = hr, CI_Lower = ci_lower, CI_Upper = ci_upper, p_value = p_value))
      
      # Also store the summary in the results list if needed
      results[[gene]] <- summary_cox
    }
  } else {
    cat("Skipping gene", gene, "due to missing values\n")
  }
}

# Display the results data frame
print(blca_results_df)
str(blca_results_df)
View(blca_results_df)


#T-ALL

all_expression<-read.csv('StJude_fpkm.txt',
                        sep='\t', stringsAsFactors = F)

dim(all_expression)
View(all_expression)

#264 19464
tail(colnames(all_expression), 55)

#T Allclinical metadata data
setwd('~/Downloads/all_phase2_target_2018_pub')
TargetClinical<-read.csv('~/Downloads/all_phase2_target_2018_pub/data_clinical_sample.txt',
                         sep='\t', stringsAsFactors = F)
View(TargetClinical)
#focus on "T cell ALL"
TALLClinical<-TargetClinical[TargetClinical$Cell.of.tumor.origin == 'T Cell ALL',]
View(TALLClinical)

clinical_sample<-read_delim("data_clinical_sample.txt")
clinical_sample <- as.data.frame(clinical_sample)
clinical_patient<-read_delim("data_clinical_patient.txt")
clinical_patient <- as.data.frame(clinical_patient)

ALL_clinical <- merge(clinical_patient, clinical_sample, by = "#Patient Identifier", 
                      all.y = TRUE)
TALL_clinical<-ALL_clinical[ALL_clinical$`Cell of tumor origin`== 
                              "T Cell ALL",]
unique(TALL_clinical$`Cell of tumor origin`)
TALL_clinical<-
  TALL_clinical[!is.na(TALL_clinical$`Cell of tumor origin`), ]

View(TALL_clinical)

#norm is the working df for T-ALL


# Define the genes of interest
tall_sur_genes <- colnames(norm)[3:19465]

# Initialize list to store results
tall_res <- list()
tall_res_df <- data.frame(Gene = character(), HR = numeric(), CI_Lower = numeric(), CI_Upper = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene and perform survival analysis
for (gene in tall_sur_genes) {
  # Ensure there are no missing values for the current gene
  gene_data <- norm[[gene]]
  if (all(!is.na(gene_data))) {
    # Create the survival formula dynamically
    # Use backticks to handle gene names with special characters
    formula <- as.formula(paste("Surv(`Overall Survival Days`, `Overall Survival Status`) ~ `", gene, "`", sep = ""))
 
    # Fit the Cox proportional hazards model
    cox_model <- tryCatch({
      coxph(formula, data = norm)
    }, error = function(e) {
      cat("Error in model fitting for gene:", gene, "\n")
      NULL
    })
    
    if (!is.null(cox_model)) {
      summary_cox <- summary(cox_model)
      
      # Extract the hazard ratio, confidence intervals, and p-value
      hr <- summary_cox$coefficients[1, "exp(coef)"]
      ci_lower <- summary_cox$conf.int[1, 1]
      ci_upper <- summary_cox$conf.int[1, 3]
      p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
      
      # Store the results in a data frame
      tall_res_df <- rbind(tall_res_df, data.frame(Gene = gene, HR = hr, CI_Lower = ci_lower, CI_Upper = ci_upper, p_value = p_value))
      
      # Also store the summary in the results list if needed
      results[[gene]] <- summary_cox
    }
  } else {
    cat("Skipping gene", gene, "due to missing values\n")
  }
}

# Display the results data frame
print(tall_res_df)
str(tall_res_df)
View(tall_res_df)


#Pan-BC, TNBC and ER+ BC
setwd("~/Downloads/brca_tcga_pub")
getwd()
Brca <- read.table(
  "data_mrna_agilent_microarray_zscores_ref_diploid_samples.txt", 
  header = TRUE)
View(Brca)
str(Brca) #data.frame':	17267 genes of  528 samples


Brca_unique <- Brca[!duplicated(Brca[, c('Hugo_Symbol')]), ]
View(Brca_unique)
rownames(Brca_unique) <- Brca_unique[, 1]

Brca_unique <- Brca_unique[, -c(1, 2)]
View(Brca_unique)
dim(Brca_unique)#[1] 17232   526

# Check for rows with any NA values
na_rows <- apply(Brca_unique, 1, function(row) any(is.na(row)))

# Print the indices of rows that contain NA values (if any)
na_row_indices <- which(na_rows)
if (length(na_row_indices) > 0) {
  cat("Rows with NA values: ", paste(na_row_indices, collapse = ", "), "\n")
} else {
  cat("No rows with NA values found.\n")
}

# Remove rows that contain NA values
Brca_clean <- Brca_unique[!na_rows, ]

# Confirm the rows with NA values have been removed
dim(Brca_clean) #[1] 14674   526
View(Brca_clean)
Brca_clean <- Brca_clean[, -c(1, 2)]


colnames(Brca_clean) <- gsub("\\.01$", "", colnames(Brca_clean))
colnames(Brca_clean) <- gsub("\\.", "-", colnames(Brca_clean))

tBrca_clean<-t(Brca_clean)

Brca_clean_exp<-as.data.frame(tBrca_clean)
View(Brca_clean_exp)

# Ensure the row names are in character format
rownames(Brca_clean) <- as.character(rownames(Brca_clean))

# Write the row names to a text file, one row name per line 
write.table(rownames(Brca_clean), file = "Brca_clean_gene.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



getwd()

TCGA_BRCA_meta<-read.csv('TCGABRCAmetadataNature2012.csv',sep=',', stringsAsFactors = F)
View(TCGA_BRCA_meta)
head(TCGA_BRCA_meta)

colnames(TCGA_BRCA_meta) <- as.character(unlist(TCGA_BRCA_meta[1, ]))

TCGA_BRCA_meta <- TCGA_BRCA_meta[-1, ]
View(TCGA_BRCA_meta)



# Convert row names of Brca_clean_exp to a column
Brca_clean_exp$TCGA_ID <- rownames(Brca_clean_exp)

# Convert row names of Brca_clean_exp to a column if it doesn't already exist
if (!"Complete_TCGA_ID" %in% colnames(Brca_clean_exp)) {
  Brca_clean_exp$Complete_TCGA_ID <- rownames(Brca_clean_exp)
}

# Rename the column in TCGA_BRCA_meta for consistency
colnames(TCGA_BRCA_meta)[which(colnames(TCGA_BRCA_meta) == "Complete TCGA ID")] <- "Complete_TCGA_ID"

# Merge data frames by Complete_TCGA_ID
brca_exp_meta <- merge(Brca_clean_exp, TCGA_BRCA_meta, by = "Complete_TCGA_ID")

# Save the common entries to a CSV file
write.csv(brca_exp_meta, "merged_common_entries.csv", row.names = FALSE)

# Print the first few rows of the merged data frame to verify
dim(brca_exp_meta)
brca_exp_meta$`OS event` <- as.numeric(brca_exp_meta$`OS event`)
brca_exp_meta$`OS Time` <- as.numeric(brca_exp_meta$`OS Time`)
# Define the genes of interest
brca_sur_genes <- colnames(brca_exp_meta)[2:14690]

# Initialize list to store results
brca_res <- list()
brca_res_df <- data.frame(Gene = character(), HR = numeric(), CI_Lower = numeric(), CI_Upper = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene and perform survival analysis
for (gene in brca_sur_genes) {
  # Ensure there are no missing values for the current gene
  gene_data <- brca_exp_meta[[gene]]
  if (all(!is.na(gene_data))) {
    # Create the survival formula dynamically
    # Use backticks to handle gene names with special characters
    formula <- as.formula(paste("Surv(`OS Time`, `OS event`) ~ `", gene, "`", sep = ""))
    
    # Fit the Cox proportional hazards model
    cox_model <- tryCatch({
      coxph(formula, data = brca_exp_meta)
    }, error = function(e) {
      cat("Error in model fitting for gene:", gene, "\n")
      NULL
    })
    
    if (!is.null(cox_model)) {
      summary_cox <- summary(cox_model)
      
      # Extract the hazard ratio, confidence intervals, and p-value
      hr <- summary_cox$coefficients[1, "exp(coef)"]
      ci_lower <- summary_cox$conf.int[1, 1]
      ci_upper <- summary_cox$conf.int[1, 3]
      p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
      
      # Store the results in a data frame
      brca_res_df <- rbind(brca_res_df, data.frame(Gene = gene, HR = hr, CI_Lower = ci_lower, CI_Upper = ci_upper, p_value = p_value))
      
      # Also store the summary in the results list if needed
      brca_res[[gene]] <- summary_cox
    }
  } else {
    cat("Skipping gene", gene, "due to missing values\n")
  }
}

View(brca_res_df)

str(brca_exp_meta$`OS event`)
str(brca_exp_meta$`OS Time`)

#TNBC


TNBC <- read.table("Brca_TNBC_clean_exp.txt", 
  header = TRUE)

TNBC<-as.data.frame(TNBC)
View(TNBC )


rownames(TNBC) <- TNBC[, 1]

TNBC<-TNBC[, -(1)]
colnames(TNBC) <- gsub("\\.01$", "", colnames(TNBC))
colnames(TNBC) <- gsub("\\.", "-", colnames(TNBC))

tnbc_column_names <- colnames(TNBC)
print(tnbc_column_names)
common_columns <- intersect(tnbc_column_names, brca_exp_meta$Complete_TCGA_ID)
print(common_columns)
tnbc_exp_meta <- brca_exp_meta[brca_exp_meta$Complete_TCGA_ID %in% common_columns, ]

View(tnbc_exp_meta)

# Print the first few rows of the merged data frame to verify
dim(tnbc_exp_meta)
tnbc_exp_meta$`OS event` <- as.numeric(tnbc_exp_meta$`OS event`)
tnbc_exp_meta$`OS Time` <- as.numeric(tnbc_exp_meta$`OS Time`)
# Define the genes of interest
tnbc_sur_genes <- colnames(tnbc_exp_meta)[2:14690]

# Initialize list to store results
tnbc_res <- list()
tnbc_res_df <- data.frame(Gene = character(), HR = numeric(), CI_Lower = numeric(), CI_Upper = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene and perform survival analysis
for (gene in tnbc_sur_genes) {
  # Ensure there are no missing values for the current gene
  gene_data <- tnbc_exp_meta[[gene]]
  if (all(!is.na(gene_data))) {
    # Create the survival formula dynamically
    # Use backticks to handle gene names with special characters
    formula <- as.formula(paste("Surv(`OS Time`, `OS event`) ~ `", gene, "`", sep = ""))
    
    # Fit the Cox proportional hazards model
    cox_model <- tryCatch({
      coxph(formula, data = tnbc_exp_meta)
    }, error = function(e) {
      cat("Error in model fitting for gene:", gene, "\n")
      NULL
    })
    
    if (!is.null(cox_model)) {
      summary_cox <- summary(cox_model)
      
      # Extract the hazard ratio, confidence intervals, and p-value
      hr <- summary_cox$coefficients[1, "exp(coef)"]
      ci_lower <- summary_cox$conf.int[1, 1]
      ci_upper <- summary_cox$conf.int[1, 3]
      p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
      
      # Store the results in a data frame
      tnbc_res_df <- rbind(tnbc_res_df, data.frame(Gene = gene, HR = hr, CI_Lower = ci_lower, CI_Upper = ci_upper, p_value = p_value))
      
      # Also store the summary in the results list if needed
      tnbc_res[[gene]] <- summary_cox
    }
  } else {
    cat("Skipping gene", gene, "due to missing values\n")
  }
}

View(tnbc_res_df)
print(tnbc_res)
#ER+ BC

ERBC <- read.table("Brca_hr_clean_exp.txt", 
                   header = TRUE)

ERBC<-as.data.frame(ERBC)
View(ERBC )

rownames(ERBC) <- ERBC[, 1]

ERBC<-ERBC[, -(1)]
colnames(ERBC) <- gsub("\\.01$", "", colnames(ERBC))
colnames(ERBC) <- gsub("\\.", "-", colnames(ERBC))

erbc_column_names <- colnames(ERBC)
print(tnbc_column_names)
common_columns <- intersect(erbc_column_names, brca_exp_meta$Complete_TCGA_ID)
print(common_columns)
erbc_exp_meta <- brca_exp_meta[brca_exp_meta$Complete_TCGA_ID %in% common_columns, ]

dim(erbc_exp_meta)
View(erbc_exp_meta)

erbc_sur_genes <- colnames(erbc_exp_meta)[2:14690]

# Initialize list to store results
erbc_res <- list()
erbc_res_df <- data.frame(Gene = character(), HR = numeric(), CI_Lower = numeric(), CI_Upper = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene and perform survival analysis
for (gene in erbc_sur_genes) {
  # Ensure there are no missing values for the current gene
  gene_data <- erbc_exp_meta[[gene]]
  if (all(!is.na(gene_data))) {
    # Create the survival formula dynamically
    # Use backticks to handle gene names with special characters
    formula <- as.formula(paste("Surv(`OS Time`, `OS event`) ~ `", gene, "`", sep = ""))
    
    # Fit the Cox proportional hazards model
    cox_model <- tryCatch({
      coxph(formula, data = erbc_exp_meta)
    }, error = function(e) {
      cat("Error in model fitting for gene:", gene, "\n")
      NULL
    })
    
    if (!is.null(cox_model)) {
      summary_cox <- summary(cox_model)
      
      # Extract the hazard ratio, confidence intervals, and p-value
      hr <- summary_cox$coefficients[1, "exp(coef)"]
      ci_lower <- summary_cox$conf.int[1, 1]
      ci_upper <- summary_cox$conf.int[1, 3]
      p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
      
      # Store the results in a data frame
      erbc_res_df <- rbind(erbc_res_df, data.frame(Gene = gene, HR = hr, CI_Lower = ci_lower, CI_Upper = ci_upper, p_value = p_value))
      
      # Also store the summary in the results list if needed
      erbc_res[[gene]] <- summary_cox
    }
  } else {
    cat("Skipping gene", gene, "due to missing values\n")
  }
}

View(erbc_res_df)


# Find the intersection of Gene values
common_genes <- Reduce(intersect, list(blca_results_df$Gene, tall_res_df$Gene, brca_res_df$Gene, tnbc_res_df$Gene, erbc_res_df$Gene))

# Subset each data frame by the common Gene values
blca_subset <- blca_results_df[blca_results_df$Gene %in% common_genes, ]
tall_subset <- tall_res_df[tall_res_df$Gene %in% common_genes, ]
brca_subset <- brca_res_df[brca_res_df$Gene %in% common_genes, ]
tnbc_subset <- tnbc_res_df[tnbc_res_df$Gene %in% common_genes, ]
erbc_subset <- erbc_res_df[erbc_res_df$Gene %in% common_genes, ]

names(blca_subset)[2:5] <- paste0("blca_", names(blca_subset)[2:5])
names(tall_subset)[2:5] <- paste0("tall_", names(tall_subset)[2:5])
names(brca_subset)[2:5] <- paste0("brca_", names(brca_subset)[2:5])
names(tnbc_subset)[2:5] <- paste0("tnbc_", names(tnbc_subset)[2:5])
names(erbc_subset)[2:5] <- paste0("erbc_", names(erbc_subset)[2:5])


# Merge the data frames on the Gene column
merged_sur_res <- Reduce(function(x, y) merge(x, y, by = "Gene"), list(blca_subset, tall_subset, brca_subset, tnbc_subset, erbc_subset))

# Print the first few rows of the merged data frame to verify
View(merged_sur_res)
dim(merged_sur_res)


# Assuming the merged data frame is named merged_df

# Extract p-values from the merged data frame
p_values_df <- merged_sur_res %>%
  dplyr::select(Gene, starts_with("blca_p_value"), starts_with("tall_p_value"), starts_with("brca_p_value"), starts_with("tnbc_p_value"), starts_with("erbc_p_value"))

# Rename columns for clarity
colnames(p_values_df) <- c("Gene", "BLCA", "TALL", "BRCA", "TNBC", "ERBC")

# Set the gene names as row names
rownames(p_values_df) <- p_values_df$Gene

# Remove the Gene column for matrix conversion
p_values_matrix <- as.matrix(p_values_df[,-1])
View(p_values_matrix)
# Check for the first few rows of the matrix to verify
head(p_values_matrix)

# Define color breaks
color_breaks <- c(0, 0.001, 0.01, 0.05, 1)

# Define color palette (needs one less color than the number of breaks)
heatmap_colors <- c("#86ebc9", "#869ceb", "grey", "white")

# Generate the heatmap
library(pheatmap)

#full transcriptomic survival p values heatmaps

pheatmap(p_values_matrix, 
                  cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                color = heatmap_colors,
                  breaks = color_breaks,
                  main = "P-values Heatmap")

GRcore<-read.csv('ConservedGR.csv',
                 sep=',', stringsAsFactors = F)
View(GRcore)

p_values_core <- p_values_df[p_values_df$Gene %in% GRcore$Genes, ]

View(p_values_core)
head(p_values_core)
write_csv(p_values_core,file = "survival GR core.csv")
write_csv(p_values_df,file = "survival full transcript.csv")














filtered_genes <- subset(p_values_core, ERBC < 0.05)
print(filtered_genes)
str(filtered_genes)
write_csv(filtered_genes,file = "ER possurvival GR core.csv")

tnbc_filtered_genes <- subset(p_values_core, TNBC < 0.05)
write_csv(tnbc_filtered_genes,file = "TNBCsurvival GR core.csv")

p_values_core_matrix <- as.matrix(p_values_core[,-1])

dim(p_values_core_matrix)

# Define color breaks
color_breaks <- c(0, 0.001, 0.01, 0.05, 1)

# Define color palette (needs one less color than the number of breaks)
heatmap_colorss <- c("red", "blue", "grey", "white")

# Generate the heatmap
library(pheatmap)



pheatmap(p_values_core_matrix, 
                cluster_rows = TRUE, 
                cluster_cols = TRUE, 
                color = heatmap_colorss,
                breaks = color_breaks,
                 main = "GR conserved P-values Heatmap")


# Extract p-values from the merged data frame
hr_values_df <- merged_sur_res %>%
  dplyr::select(Gene, starts_with("blca_HR"), starts_with("tall_HR"), starts_with("brca_HR"), starts_with("tnbc_HR"), starts_with("erbc_HR"))

# Rename columns for clarity
colnames(hr_values_df) <- c("Gene", "BLCA", "TALL", "BRCA", "TNBC", "ERBC")

# Set the gene names as row names
rownames(hr_values_df) <- hr_values_df$Gene

# Remove the Gene column for matrix conversion
hr_values_matrix <- as.matrix(hr_values_df[,-1])
View(hr_values_matrix)
# Check for the first few rows of the matrix to verify
head(hr_values_matrix)


pheatmap(-log(hr_values_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         color = colorRampPalette(c("white", "blue"))(50),
         main = "GR conserved HR Heatmap")


hr_values_core <- hr_values_df[hr_values_df$Gene %in% GRcore$Genes, ]
hr_values_core <- hr_values_core[rownames(hr_values_core) != "BNC1", ]
hr_values_core_matrix <- as.matrix(hr_values_core[,-1])
View(hr_values_core_matrix)

# Check for the first few rows of the matrix to verify
head(hr_values_matrix)


pheatmap(log10(hr_values_core_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "GR conserved HRs Heatmap")

library(pheatmap)
library(grid)



pheatmap(
  -log(hr_values_matrix), 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  color = colorRampPalette(c("white", "blue"))(50),
  main = "GR Conserved HR Heatmap",
  display_numbers = FALSE,  # Remove row names
  border_color = NA,  # Remove borders around cells
  cellwidth = NA,  # Automatic width for cells
  cellheight = NA,  # Automatic height for cells
  treeheight_row = 50,  # Adjust the dendrogram height for rows
  treeheight_col = 50,  # Adjust the dendrogram height for columns
  fontsize_row = 10,  # Adjust row font size (if row names were to be displayed)
  fontsize_col = 10,  # Adjust column font size
  fontsize = 12,  # Adjust overall font size
  clustering_distance_rows = "euclidean",  # Distance metric for rows
  clustering_distance_cols = "euclidean",  # Distance metric for columns
  clustering_method = "complete"  # Clustering method
)

# Custom function to make dendrogram lines bolder
bold_dendrogram <- function(hm, line_width = 1.5) {
  grid::grid.gedit("dendrogram_segment", gp = grid::gpar(lwd = line_width))
}

# Adjust the dendrogram line width
bold_dendrogram(pheatmap(
  log(hr_values_matrix), 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  main = "GR Conserved HR Heatmap",
  display_numbers = FALSE,  # Remove row names
  border_color = NA,  # Remove borders around cells
  cellwidth = NA,  # Automatic width for cells
  cellheight = NA,  # Automatic height for cells
  treeheight_row = 50,  # Adjust the dendrogram height for rows
  treeheight_col = 50,  # Adjust the dendrogram height for columns
  fontsize_row = 10,  # Adjust row font size (if row names were to be displayed)
  fontsize_col = 10,  # Adjust column font size
  fontsize = 12,  # Adjust overall font size
  clustering_distance_rows = "euclidean",  # Distance metric for rows
  clustering_distance_cols = "euclidean",  # Distance metric for columns
  clustering_method = "complete"  # Clustering method
))


# Load necessary library
library(pheatmap)
library(grid)




#REVISE SURVIVAL P values 

# Define color breaks
color_breaks <- c(0, 0.001, 0.01, 0.05, 1)

# Define color palette (needs one less color than the number of breaks)
heatmap_colorss <- c("red", "blue", "grey", "white")

# Create the heatmap without displaying row names
hmsurp <- pheatmap(
  p_values_core_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  main = "GR Conserved p Heatmap",
  color = heatmap_colorss,
  breaks = color_breaks,
  show_rownames = FALSE,  # Remove row names
  border_color = "black",  # Remove borders around cells
  cellwidth = 20,  
  cellheight = NA,  # Automatic height for cells
  treeheight_row = 75,  # Adjust the dendrogram height for rows
  treeheight_col = 50,  # Adjust the dendrogram height for columns
  fontsize_row = 10,  # Adjust row font size (if row names were to be displayed)
  fontsize_col = 10,  # Adjust column font size
  fontsize = 12,  # Adjust overall font size
  clustering_distance_rows = "euclidean",  # Distance metric for rows
  clustering_distance_cols = "euclidean",  # Distance metric for columns
  clustering_method = "complete"  # Clustering method
)

# Load necessary library
library(pheatmap)

# Example data (replace this with your actual data)
set.seed(123)  # For reproducibility
p_values_core_matrix <- matrix(runif(100, min = 0, max = 1), nrow = 10)

# Define color breaks
color_breaks <- c(0, 0.001, 0.01, 0.05, 1)

# Define color palette (needs one less color than the number of breaks)
heatmap_colorss <- c("red", "blue", "grey", "white")

# Create the heatmap without displaying row names
hmsurp <- pheatmap(
  p_values_core_matrix, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  main = "GR Conserved p Heatmap",
  color = heatmap_colorss,
  breaks = color_breaks,
  show_rownames = FALSE,  # Remove row names
  border_color = "black",  # Borders around cells
  cellwidth = 30,  
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









# Function to make dendrogram lines bolder
bold_dendrogram_lines <- function(hm, line_width = 1.5) {
  # Modify row dendrogram
  if (!is.null(hm$gtable$grobs$grob[[4]]$children[[1]]$children)) {
    row_dend <- hm$gtable$grobs$grob[[4]]$children[[1]]$children[[1]]
    row_dend$gp$lwd <- line_width
    hm$gtable$grobs$grob[[4]]$children[[1]]$children[[1]] <- row_dend
  }
  
  # Modify column dendrogram
  if (!is.null(hm$gtable$grobs$grob[[6]]$children[[1]]$children)) {
    col_dend <- hm$gtable$grobs$grob[[6]]$children[[1]]$children[[1]]
    col_dend$gp$lwd <- line_width
    hm$gtable$grobs$grob[[6]]$children[[1]]$children[[1]] <- col_dend
  }
  
  grid.newpage()
  grid.draw(hm$gtable)
}

# Apply the function to make dendrogram lines bolder
bold_dendrogram_lines(hmsurp, line_width = 10)






