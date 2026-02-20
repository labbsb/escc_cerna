# Task: Survival Analysis

library(TCGAbiolinks)            # Downloads and processes TCGA data
library(survminer)               # Visualizes survival curves
library(survival)                # Fits Cox models and Kaplan-Meier survival curves
library(SummarizedExperiment)    # Stores large omics datasets with metadata
library(tidyverse)               # Loads core tidy data packages
library(DESeq2)                  # Performs normalization and differential expression analysis on count data
library(tidyr)                   # Manipulates data frames
library(dplyr)                   # Reshapes data
library(miRBaseConverter)        # Converts between different miRNA naming formats

# Get clinical data for TCGA-ESCA cohort -------------------
clinical_esca <- GDCquery_clinic("TCGA-ESCA")

# Check if the clinical dataset contains columns for survival analysis 
# (vital status, last follow-up, and days to death)
any(colnames(clinical_esca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))

# Find the positions (column indices) of the above clinical variables
which(colnames(clinical_esca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))

# Extract those columns from the clinical dataset (columns 57, 66, 67 correspond 
# to survival-related variables)
clinical_esca[,c(57,66,67)]

# Filter clinical data to keep only ESCC cases:
#   - Include patients with primary diagnosis = "Squamous cell carcinoma, NOS"
#   - Exclude cases where tissue/organ of origin = "Cardia, NOS"
clinical_escc <- clinical_esca[clinical_esca$primary_diagnosis == "Squamous cell carcinoma, NOS"
                               & clinical_esca$tissue_or_organ_of_origin != "Cardia, NOS",]

# Change certain values the way they are encoded
clinical_escc$deceased <- ifelse(clinical_escc$vital_status == "Alive", FALSE, TRUE)

# Create an "overall survival" variable:
#   - Use days_to_death for patients who died
#   - Use days_to_last_follow_up for patients still alive
clinical_escc$overall_survival <- ifelse(clinical_escc$vital_status == "Alive",
                                         clinical_escc$days_to_last_follow_up,
                                         clinical_escc$days_to_death)

# Identify patients with overall survival of 30 days or less
clinical_escc$submitter_id[clinical_escc$overall_survival <= 30]

# Exclude patients with overall survival < 30 days
# submitter_id: TCGA-L5-A43H, TCGA-L5-A43J, TCGA-IG-A3Y9, TCGA-IG-A50L, TCGA-LN-A9FO
clinical_escc <- clinical_escc[clinical_escc$overall_survival > 30, ]







# ------------Get mRNA expression data for TCGA-ESCA cohort ------------------



# Query RNA-Seq gene expression data (STAR counts) 
# for primary tumor samples from the TCGA-ESCA project
query_esca_mrna <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor"),
  access = "open")

# Download queried RNA-Seq gene expression data (STAR counts) for TCGA-ESCA
GDCdownload(query_esca_mrna)

# Prepare and extract RNA-Seq count matrix for TCGA-ESCA
tcga_esca_data_mrna <- GDCprepare(query_esca_mrna, summarizedExperiment = TRUE) # create SummarizedExperiment object
esca_matrix_mrna <- assay(tcga_esca_data_mrna, "unstranded")                    # extract unstranded counts
esca_matrix_mrna[1:10,1:10]                                                     # preview first 10x10 entries
dim(esca_matrix_mrna)                                                           # check dimensions of the matrix

# Extract gene-level (row) and sample-level (column) metadata 
# from the SummarizedExperiment object
gene_metadata_mrna <- as.data.frame(rowData(tcga_esca_data_mrna))               # gene annotations
coldata_esca_mrna <- as.data.frame(colData(tcga_esca_data_mrna))                # sample annotations

# Count unique patients who have more than one RNA-seq sample
sum(table(coldata_esca_mrna$patient) > 1)                                       # answer: 0 patients

# Filter expression matrix and sample metadata to include only ESCC patients
# Match sample IDs (submitter_id) in clinical_escc with those in colData
samples_to_keep <- coldata_esca_mrna$submitter_id %in% clinical_escc$submitter_id

# Subset expression matrix and sample metadata to keep only ESCC patients
escc_matrix_mrna <- esca_matrix_mrna[, samples_to_keep]
coldata_escc_mrna <- coldata_esca_mrna[samples_to_keep, ]

# Prepare DESeq2 dataset from ESCC count matrix and sample metadata 
# (for downstream variance-stabilizing transformation and survival analysis) 
dds_mrna <- DESeqDataSetFromMatrix(countData = escc_matrix_mrna,
                                   colData = coldata_escc_mrna,
                                   design = ~ 1)


# Filter out low-count genes: keep only genes with at least 10 reads across all samples
keep <- rowSums(counts(dds_mrna)) >= 10
dds_mrna <- dds_mrna[keep,]

# Apply variance-stabilizing transformation (VST) to normalized counts
vsd_mrna <- vst(dds_mrna, blind=FALSE)

# Extract transformed expression matrix
escc_matrix_vst_mrna <- assay(vsd_mrna)

# Read gene list
gene_list <- readLines("mrna.txt")

# Display the loaded gene list
gene_list

# Create output folder for KM plots of ESCC genes
dir.create("survival_plots_mrna", showWarnings = FALSE)

# Initialize storage
logrank_results <- list()

# --------------- LOOP: Log Rank test & Kaplan-Meier estimator ----------------
for (gene in gene_list) {
  message("Processing gene: ", gene)
  
  df <- escc_matrix_vst_mrna %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_id') %>%
    gather(key = 'case_id', value = 'counts', -gene_id) %>%
    left_join(gene_metadata_mrna, by = "gene_id") %>%
    filter(gene_name == gene)
  
  if (nrow(df) == 0) {
    message("Gene not found: ", gene)
    next
  }
  
  # Stratify based on median
  med <- median(df$counts, na.rm = TRUE)
  df$strata <- ifelse(df$counts >= med, "HIGH", "LOW")
  
  # Skip if only one group
  if (length(unique(df$strata)) < 2) {
    message("Only one group for: ", gene)
    next
  }
  
  
  # Prepare clinical IDs
  df$case_id <- substr(df$case_id, 1, 12)
  df <- merge(df, clinical_escc, by.x = "case_id", by.y = "submitter_id")
  
  
  
  # Check for necessary columns
  if (!all(c("overall_survival", "deceased", "strata") %in% colnames(df))) {
    message("Missing clinical data for: ", gene)
    next
  }
  
  # ------------- Log-Rank Test -------------
  fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = df)
  logrank_results[[gene]] <- data.frame(
    gene = gene,
    strata = names(fit2$n),
    n = fit2$n,
    observed = fit2$obs,
    expected = fit2$exp,
    chisq = fit2$chisq,
    pvalue_logrank = 1 - pchisq(fit2$chisq, df = 1)
  )
  
  
  
  # ------------- Plot -------------
  fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = df)
  p <- ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE,
    title = paste("Survival for", gene),
    palette = c("blue", "red")
  )
  
  ggsave(
    filename = paste0("survival_plots_mrna/", gene, "_survival.png"),
    plot = p$plot, width = 6, height = 5, dpi = 300
  )
}

# ------------------ RESULT: Log-Rank test  ------------------

# Combine data frames
logrank_df <- do.call(rbind, logrank_results)


# Save to CSV
write.csv(logrank_df, "mrna_logrank_results.csv", row.names = FALSE)




# ------------Get lncRNA expression data for TCGA-ESCA cohort ------------------

# Read gene list
lncRNA_list <- readLines("lncRNA.txt")

# Create output folder for KM plots of ESCC genes
dir.create("survival_plots_lncRNA", showWarnings = FALSE)

# Initialize storage
logrank_results <- list()

# --------------- LOOP: Log Rank test & Kaplan-Meier estimator ----------------
for (lncRNA in lncRNA_list) {
  message("Processing gene: ", lncRNA)
  
  df <- escc_matrix_vst_mrna %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_id') %>%
    gather(key = 'case_id', value = 'counts', -gene_id) %>%
    left_join(gene_metadata_mrna, by = "gene_id") %>%
    filter(gene_name == lncRNA)
  
  if (nrow(df) == 0) {
    message("lncRNA not found: ", lncRNA)
    next
  }
  
  # Stratify based on median
  med <- median(df$counts, na.rm = TRUE)
  df$strata <- ifelse(df$counts >= med, "HIGH", "LOW")
  
  # Skip if only one group
  if (length(unique(df$strata)) < 2) {
    message("Only one group for: ", lncRNA)
    next
  }
  
  
  # Prepare clinical IDs
  df$case_id <- substr(df$case_id, 1, 12)
  df <- merge(df, clinical_escc, by.x = "case_id", by.y = "submitter_id")
  
  
  
  # Check for necessary columns
  if (!all(c("overall_survival", "deceased", "strata") %in% colnames(df))) {
    message("Missing clinical data for: ", lncRNA)
    next
  }
  
  # ------------- Log-Rank Test -------------
  fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = df)
  logrank_results[[lncRNA]] <- data.frame(
    lncRNA = lncRNA,
    strata = names(fit2$n),
    n = fit2$n,
    observed = fit2$obs,
    expected = fit2$exp,
    chisq = fit2$chisq,
    pvalue_logrank = 1 - pchisq(fit2$chisq, df = 1)
  )
  
  
  
  # ------------- Plot -------------
  fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = df)
  p <- ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE,
    title = paste("Survival for", lncRNA),
    palette = c("blue", "red")
  )
  
  ggsave(
    filename = paste0("survival_plots_lncRNA/", lncRNA, "_survival.png"),
    plot = p$plot, width = 6, height = 5, dpi = 300
  )
}

# ------------------ RESULT: Log-Rank test  ------------------

# Combine data frames
logrank_df <- do.call(rbind, logrank_results)


# Save to CSV
write.csv(logrank_df, "lncRNA_logrank_results.csv", row.names = FALSE)









# ------------Get miRNA expression data for TCGA-ESCA cohort ------------------


# Build query for miRNA-Seq isoform-level expression (BCGSC workflow) in primary tumors
query_esca_mirna <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "miRNA-Seq",
  workflow.type = "BCGSC miRNA Profiling",
  data.type = "Isoform Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open")

# Download queried miRNA-Seq isoform level expression data (BCGSC) for TCGA-ESCA
GDCdownload(query_esca_mirna)

# Prepare and extract miRNA-Seq count matrix for TCGA-ESCA
tcga_esca_data_mirna <- GDCprepare(query_esca_mirna, summarizedExperiment = TRUE)

# Show column names (sample identifiers and metadata fields)
colnames(tcga_esca_data_mirna)

# Remove non-essential metadata columns from the miRNA expression data
tcga_esca_data_mirna <- tcga_esca_data_mirna[, !names(tcga_esca_data_mirna) %in% c("miRNA_ID", "isoform_coords", "reads_per_million_miRNA_mapped", "cross-mapped")]

# Aggregate isoform-level counts to miRNA-level:
# For each (miRNA_region, sample barcode) pair, sum the read counts across isoforms
df_grouped <- tcga_esca_data_mirna %>%
  group_by(miRNA_region, barcode) %>%
  summarise(read_count = sum(read_count, na.rm = TRUE), .groups = "drop")

# Reshape data into wide matrix format (miRNAs x samples)
# - `miRNA_region` stays as row identifier
# - Each column becomes a sample (barcode)
# - Cell values are the summed read counts from Step 1
# - Missing values are filled with 0
df_wide <- df_grouped %>%
  pivot_wider(
    id_cols = miRNA_region,
    names_from = barcode,
    values_from = read_count,
    values_fill = list(read_count = 0)
  )

# Convert wide tibble to data frame with miRNA names as row identifiers
# - Transform tibble into a standard data.frame
# - Set `miRNA_region` as rownames (so miRNAs act like gene IDs in the mRNA matrix)
# - Remove the redundant `miRNA_region` column
df_wide_df <- as.data.frame(df_wide)
rownames(df_wide_df) <- df_wide_df$miRNA_region
df_wide_df$miRNA_region <- NULL

# Filter out unwanted rows (non-mature miRNA categories)
# - Remove rows labeled as "precursor", "stemloop", or "unannotated"
# - Keeps only mature/annotated miRNAs for downstream analysis
df_filtered <- df_wide_df[!rownames(df_wide_df) %in% c("precursor", "stemloop", "unannotated"), ]

# Clean row names by removing the "mature," prefix
rownames(df_filtered) <- sub("^mature,", "", rownames(df_filtered))

# Save filtered data as mirna_counts
mirna_counts <- df_filtered

# mirna_metadata: rows = rownames of mirna_counts (miRNA IDs)
mirna_metadata <- data.frame(accession_number = rownames(mirna_counts))

# Extract mimat IDs from your metadata
accession_list <- mirna_metadata$accession_number

# Convert to miRNA names (miRBase version 22)
converted <- miRNA_AccessionToName(accession_list, targetVersion = "v21")

# Display number of NAs
sum(is.na(converted$TargetName))       # number of NAs: 0

# Identify duplicated mimat IDs (SourceID)
duplicated_ids <- converted$TargetName[duplicated(converted$TargetName)]
duplicated_ids <- unique(duplicated_ids)  # just unique IDs
duplicated_ids  # no duplicates


# Keep all rows from mirna_metadata (left join)
mirna_metadata <- merge(
  mirna_metadata,
  converted[, c("Accession", "TargetName")],  # only keep relevant columns
  by.x = "accession_number", 
  by.y = "Accession",
  all.x = TRUE
)

# Rename TargetName column to miRNA_id
colnames(mirna_metadata)[colnames(mirna_metadata) == "TargetName"] <- "miRNA_id"

# Count NAs
sum(is.na(mirna_metadata$miRNA_id))

# Set rownames to miRNA_id
rownames(mirna_metadata) <- mirna_metadata$miRNA_id

# Assign miRNA_id as rownames for counts matrix
rownames(mirna_counts) <- mirna_metadata$miRNA_id

# coldata_mirna: rows = colnames of mirna_counts (TCGA sample barcodes)
coldata_mirna <- data.frame(barcode = colnames(mirna_counts))
rownames(coldata_mirna) <- coldata_mirna$barcode

# Add submitter_id (first 12 characters of barcode)
coldata_mirna$submitter_id <- substr(coldata_mirna$barcode, 1, 12)

# Ensure values are numeric
mirna_matrix <- as.matrix(sapply(mirna_counts, as.numeric))

# Set rownames
rownames(mirna_matrix) <- rownames(mirna_counts)

# Check result
class(mirna_matrix)
dim(mirna_matrix)

# Filter expression matrix and sample metadata to include only ESCC patients
# Match sample IDs (submitter_id) in clinical_escc with those in colData
samples_to_keep_mirna <- coldata_mirna$submitter_id %in% clinical_escc$submitter_id

# Subset expression matrix and sample metadata to keep only ESCC patients
escc_mirna_matrix <- mirna_matrix[, samples_to_keep_mirna]                      # "TCGA-LN-A49V" is missing from miRNA-Seq
coldata_escc_mirna <- coldata_mirna[samples_to_keep_mirna, ]

# Count unique patients who have more than one RNA-seq sample
sum(table(coldata_escc_mirna$submitter_id) > 1)  

# Show submitter_id with more than one sample
dup_ids <- names(table(coldata_escc_mirna$submitter_id)[table(coldata_escc_mirna$submitter_id) > 1])
dup_ids                                                                         # "TCGA-IC-A6RF" & "TCGA-IG-A3YB"


# Make sure submitter_id matches column order
sample_ids <- coldata_escc_mirna$submitter_id
names(sample_ids) <- colnames(escc_mirna_matrix)

sample_ids

# Collapse duplicates and round to integers in one step
escc_matrix_avg <- sapply(unique(sample_ids), function(pid) {
  samples <- names(sample_ids[sample_ids == pid])
  # Average if multiple samples, otherwise take the single column
  mat <- escc_mirna_matrix[, samples, drop = FALSE]
  as.integer(round(rowMeans2(mat)))
})


# Create coldata_mirna_avg
coldata_mirna_avg <- data.frame(
  submitter_id = colnames(escc_matrix_avg),
  row.names = colnames(escc_matrix_avg)
)

# Set rownames
rownames(escc_matrix_avg) <- rownames(mirna_counts)


# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds_mirna <- DESeqDataSetFromMatrix(countData = escc_matrix_avg,
                                    colData = coldata_mirna_avg,
                                    design = ~ 1)


# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds_mirna)) >= 10
dds_mirna <- dds_mirna[keep,]

# vst 
vst_mirna <- varianceStabilizingTransformation(dds_mirna, blind=FALSE)
escc_matrix_vst_mirna <- assay(vst_mirna)
escc_matrix_vst_mirna[1:10,1:10]


# Read miRNA list
mirna_list <- mirna_metadata$miRNA_id

mirna_list

# Create output folder
dir.create("survival_plots_mirna", showWarnings = FALSE)



# Initialize storage
logrank_results_mirna <- list()

# ------------------ LOOP ------------------
for (mirna in mirna_list) {
  message("Processing gene: ", mirna)
  
  df_mirna <- escc_matrix_vst_mirna %>%
    as.data.frame() %>%
    rownames_to_column(var = 'miRNA_id') %>%
    gather(key = 'case_id', value = 'counts', -miRNA_id) %>%
    left_join(mirna_metadata, by = "miRNA_id") %>%
    filter(miRNA_id == mirna)
  
  if (nrow(df_mirna) == 0) {
    message("miRNA not found: ", mirna)
    next
  }
  
  # Stratify based on median
  med_mirna <- median(df_mirna$counts, na.rm = TRUE)
  df_mirna$strata <- ifelse(df_mirna$counts >= med_mirna, "HIGH", "LOW")
  
  # Skip if only one group
  if (length(unique(df_mirna$strata)) < 2) {
    message("Only one group for: ", mirna)
    next
  }
  
  # Prepare clinical IDs
  df_mirna <- merge(df_mirna, clinical_escc, by.x = "case_id", by.y = "submitter_id")
  
  # Check for necessary columns
  if (!all(c("overall_survival", "deceased", "strata") %in% colnames(df_mirna))) {
    message("Missing clinical data for: ", mirna)
    next
  }
  
  # ------------- Log-Rank Test -------------
  fit2_mirna <- survdiff(Surv(overall_survival, deceased) ~ strata, data = df_mirna)
  logrank_results_mirna[[mirna]] <- data.frame(
    mirna = mirna,
    strata = names(fit2_mirna$n),
    n = fit2_mirna$n,
    observed = fit2_mirna$obs,
    expected = fit2_mirna$exp,
    chisq = fit2_mirna$chisq,
    pvalue_logrank = 1 - pchisq(fit2_mirna$chisq, df = 1)
  )
  
  
  # ------------- Plot -------------
  fit_mirna <- survfit(Surv(overall_survival, deceased) ~ strata, data = df_mirna)
  p <- ggsurvplot(
    fit_mirna, data = df_mirna, pval = TRUE, risk.table = TRUE,
    title = paste("Survival for", mirna),
    palette = c("blue", "red")
  )
  
  ggsave(
    filename = paste0("survival_plots_mirna/", mirna, "_survival.png"),
    plot = p$plot, width = 6, height = 5, dpi = 300
  )
}


# ------------------ RESULTS ------------------

# Combine data frames
logrank_df_mirna <- do.call(rbind, logrank_results_mirna)

# Save to CSV
write.csv(logrank_df_mirna, "logrank_mirna.csv", row.names = FALSE)






