#############################
## 1. Load libraries
#############################
library(DESeq2)
library(dplyr)
library(ggplot2)
setwd("/mnt/sdb2/")
#############################
## 2. Read count matrix
#############################

data <- read.table("ksarkytbayev/nextflow_24_04_25/mirna_quant/edger_qc/mature_counts.csv", header = TRUE, sep = ",")

rownames(data)<- data$X
data$X <- NULL
data <- t(data)


data <- data +1

# Identify control and sample columns
control_cols <- grep("^control_", colnames(data), value = TRUE)
sample_cols <- setdiff(colnames(data), control_cols)

# Reorder the columns
data <- data[, c(control_cols, sample_cols)]

# Get column names
all_cols <- colnames(data)

# Separate control and sample columns
control_cols <- grep("^control_", all_cols, value = TRUE)
sample_cols <- setdiff(all_cols, control_cols)

# Sort each group
control_cols <- sort(control_cols)
sample_cols <- sort(sample_cols)

# Reorder data
data <- data[, c(control_cols, sample_cols)]
#############################
## 4. Sample vector (patients)
#############################
sample_vector <- c(
  "control","control","control","control",
  "X001T","X001T","X001T",
  "X002T","X002T","X002T",
  "X003T","X003T","X003T",
  "X004T","X004T","X004T",
  "X005T","X005T","X005T",
  "X007T","X007T","X007T",
  "X008T","X008T","X008T",
  "X011T","X011T","X011T",
  "X012T","X012T","X012T",
  "X015T","X015T","X015T",
  "X016T","X016T","X016T",
  "X017T","X017T","X017T",
  "X023T","X023T","X023T",
  "X025T","X025T",
  "X029T","X029T","X029T",
  "X033T","X033T","X033T"
)

#############################
## 5. Stage information
#############################
stage_map <- list(
  II  = c("X012T", "X017T"),
  III = c("X001T","X002T","X003T","X004T","X005T",
          "X007T","X008T","X015T","X016T",
          "X023T","X029T","X033T"),
  IV  = c("X011T","X025T")
)

#############################
## 6. Build sample metadata
#############################
sample_ids <- colnames(data)

sample_info <- data.frame(
  sample_id = sample_ids,
  group = ifelse(grepl("^control_", sample_ids), "control", "tumor"),
  patient = sample_vector,
  stage = NA,
  row.names = sample_ids
)

# Assign stages
for (st in names(stage_map)) {
  sample_info$stage[sample_info$patient %in% stage_map[[st]]] <- st
}

# Controls
sample_info$stage[sample_info$group == "control"] <- "control"

# Factors
sample_info$group   <- factor(sample_info$group, levels = c("control","tumor"))
sample_info$stage   <- factor(sample_info$stage, levels = c("control","II","III","IV"))
sample_info$patient <- factor(sample_info$patient)

#############################
## 7. DESeq A: POOLED
#############################
dds_pooled <- DESeqDataSetFromMatrix(
  countData = data,
  colData   = sample_info,
  design    = ~ group
)

dds_pooled <- DESeq(dds_pooled)

res_pooled <- results(dds_pooled, contrast = c("group","tumor","control"))

write.csv(
  as.data.frame(res_pooled),
  file = "DESeq_pooled_tumor_vs_control.csv"
)

#############################
## 8. DESeq B: BY STAGE
#############################
dds_stage <- DESeqDataSetFromMatrix(
  countData = data,
  colData   = sample_info,
  design    = ~ stage
)

dds_stage <- DESeq(dds_stage)

res_stage_2  <- results(dds_stage, contrast = c("stage","II","control"))
res_stage_3 <- results(dds_stage, contrast = c("stage","III","control"))
res_stage_4  <- results(dds_stage, contrast = c("stage","IV","control"))

write.csv(as.data.frame(res_stage_2),  "DESeq_stage_II_vs_control.csv")
write.csv(as.data.frame(res_stage_3), "DESeq_stage_III_vs_control.csv")
write.csv(as.data.frame(res_stage_4),  "DESeq_stage_IV_vs_control.csv")

#############################
## 9. DESeq C: BY INDIVIDUAL
#############################
dds_individual <- DESeqDataSetFromMatrix(
  countData = data,
  colData   = sample_info,
  design    = ~ patient
)

dds_individual <- DESeq(dds_individual)

patients <- setdiff(levels(sample_info$patient), "control")

# Store individual patient results
patient_results <- list()

for (p in patients) {
  res <- results(dds_individual, contrast = c("patient", p, "control"))
  write.csv(
    as.data.frame(res),
    file = paste0("DESeq_patient_", p, "_vs_control.csv")
  )
  # Store for later use
  patient_results[[p]] <- res
}

#############################
## 10. Optional: VST for PCA
#############################
vsd <- vst(dds_pooled, blind = FALSE)
saveRDS(vsd, "VST_pooled.rds")

#############################
## 11. Filter significant results
#############################

replace_dot_with_dash_rownames <- function(df) {
  stopifnot(is.data.frame(df))
  rownames(df) <- gsub("\\.", "-", rownames(df))
  return(df)
}

## POOLED
res_pooled_sig <- res_pooled %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(padj),
                padj < 0.05, 
                abs(log2FoldChange) > 1)
res_pooled_sig <- replace_dot_with_dash_rownames(res_pooled_sig)

## Stage II vs control
res_stage_2_sig <- res_stage_2 %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(padj),
                padj < 0.05,
                abs(log2FoldChange) > 1)
res_stage_2_sig <- replace_dot_with_dash_rownames(res_stage_2_sig)

## Stage III vs control
res_stage_3_sig <- res_stage_3 %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(padj),
                padj < 0.05,
                abs(log2FoldChange) > 1)
res_stage_3_sig <- replace_dot_with_dash_rownames(res_stage_3_sig)

## Stage IV vs control
res_stage_4_sig <- res_stage_4 %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(padj),
                padj < 0.05,
                abs(log2FoldChange) > 1)
res_stage_4_sig <- replace_dot_with_dash_rownames(res_stage_4_sig)

#############################
## 12. ALL STAGES ANALYSIS
## Find miRNAs significant in ALL three stages with consistent direction
#############################

# Get miRNAs present in all three stages
mirnas_stage_2 <- rownames(res_stage_2_sig)
mirnas_stage_3 <- rownames(res_stage_3_sig)
mirnas_stage_4 <- rownames(res_stage_4_sig)

# Intersection: present in ALL stages
mirnas_all_stages <- Reduce(intersect, list(mirnas_stage_2, mirnas_stage_3, mirnas_stage_4))

# Check consistent direction
all_stages_up <- character()
all_stages_down <- character()

for (mirna in mirnas_all_stages) {
  fc2 <- res_stage_2_sig[mirna, "log2FoldChange"]
  fc3 <- res_stage_3_sig[mirna, "log2FoldChange"]
  fc4 <- res_stage_4_sig[mirna, "log2FoldChange"]
  
  # All positive
  if (fc2 > 0 && fc3 > 0 && fc4 > 0) {
    all_stages_up <- c(all_stages_up, mirna)
  }
  # All negative
  if (fc2 < 0 && fc3 < 0 && fc4 < 0) {
    all_stages_down <- c(all_stages_down, mirna)
  }
}

n_all_stages_up <- length(all_stages_up)
n_all_stages_down <- length(all_stages_down)

cat("ALL STAGES analysis (sig in II, III, IV with consistent direction)\n")
cat("Upregulated in all stages:", n_all_stages_up, "\n")
cat("Downregulated in all stages:", n_all_stages_down, "\n")

#############################
## 13. ALL SAMPLES ANALYSIS
## Find miRNAs significant in ALL individual patients with consistent direction
#############################

# Process each patient's results
patient_sig_list <- list()

for (p in names(patient_results)) {
  res_p_sig <- patient_results[[p]] %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(padj),
                  padj < 0.05,
                  abs(log2FoldChange) > 1)
  res_p_sig <- replace_dot_with_dash_rownames(res_p_sig)
  patient_sig_list[[p]] <- res_p_sig
}

# Get miRNAs present in ALL patients
all_patient_mirnas <- lapply(patient_sig_list, rownames)
mirnas_all_samples <- Reduce(intersect, all_patient_mirnas)

# Check consistent direction across ALL samples
all_samples_up <- character()
all_samples_down <- character()

for (mirna in mirnas_all_samples) {
  fcs <- sapply(patient_sig_list, function(df) {
    if (mirna %in% rownames(df)) {
      df[mirna, "log2FoldChange"]
    } else {
      NA
    }
  })
  
  # All positive
  if (all(fcs > 0, na.rm = TRUE)) {
    all_samples_up <- c(all_samples_up, mirna)
  }
  # All negative
  if (all(fcs < 0, na.rm = TRUE)) {
    all_samples_down <- c(all_samples_down, mirna)
  }
}

n_all_samples_up <- length(all_samples_up)
n_all_samples_down <- length(all_samples_down)

cat("\nALL SAMPLES analysis (sig in all 16 patients with consistent direction)\n")
cat("Upregulated in all samples:", n_all_samples_up, "\n")
cat("Downregulated in all samples:", n_all_samples_down, "\n")

#############################
## 14. INDIVIDUAL SAMPLE COUNTS (for dots)
#############################

patient_counts <- data.frame()

for (p in names(patient_sig_list)) {
  res_p_sig <- patient_sig_list[[p]]
  
  n_up <- sum(res_p_sig$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res_p_sig$log2FoldChange < 0, na.rm = TRUE)
  
  patient_counts <- rbind(patient_counts, data.frame(
    patient = p,
    upregulated = n_up,
    downregulated = n_down
  ))
}

# Save individual counts
write.csv(patient_counts, "individual_patient_miRNA_counts.csv", row.names = FALSE)

cat("\nIndividual patient counts saved to: individual_patient_miRNA_counts.csv\n")

#############################
## 15. POOLED COUNTS
#############################

n_pos_pooled <- sum(res_pooled_sig$log2FoldChange > 0, na.rm = TRUE)
n_neg_pooled <- sum(res_pooled_sig$log2FoldChange < 0, na.rm = TRUE)

cat("\nPOOLED analysis\n")
cat("Upregulated:", n_pos_pooled, "\n")
cat("Downregulated:", n_neg_pooled, "\n")

#############################
## 16. CREATE BARPLOT DATA
#############################

bar_data <- data.frame(
  dataset = c("Pooled", "Pooled", 
              "All Stages", "All Stages",
              "All Samples", "All Samples"),
  direction = c("Upregulated", "Downregulated",
                "Upregulated", "Downregulated",
                "Upregulated", "Downregulated"),
  count = c(n_pos_pooled, n_neg_pooled,
            n_all_stages_up, n_all_stages_down,
            n_all_samples_up, n_all_samples_down)
)

# Set factor levels for proper ordering
bar_data$dataset <- factor(bar_data$dataset, 
                           levels = c("Pooled", "All Stages", "All Samples"))

#############################
## 17. CREATE DOT DATA (individual samples)
#############################

# Prepare data for dots overlay on "All Samples" category
dots_data_up <- data.frame(
  dataset = "All Samples",
  direction = "Upregulated",
  count = patient_counts$upregulated,
  patient = patient_counts$patient
)

dots_data_down <- data.frame(
  dataset = "All Samples",
  direction = "Downregulated",
  count = patient_counts$downregulated,
  patient = patient_counts$patient
)

dots_data <- rbind(dots_data_up, dots_data_down)
dots_data$dataset <- factor(dots_data$dataset, levels = c("Pooled", "All Stages", "All Samples"))

#############################
## 18. CREATE THE BARPLOT WITH DOTS
#############################

geom_point(data = dots_data, 
           aes(x = dataset, y = count, color = direction),
           position = position_dodge(width = 0.9),
           size = 2.5, alpha = 0.6, shape = 16) +
  # Add figure label 'a' in upper right
  annotate("text", x = Inf, y = Inf, label = "a", 
           hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold") +
  theme_minimal() +
  labs(
    title = "Number of Significantly Expressed miRNAs Across Three Differential Expression Strategies",
    x = "Analysis Strategy",
    y = "Number of miRNAs"
  ) +
  scale_fill_manual(values = c("Upregulated" = "coral3", "Downregulated" = "deepskyblue4"),
                    name = "Regulation") +
  scale_color_manual(values = c("Upregulated" = "coral3", "Downregulated" = "deepskyblue4"),
                     guide = "none") +  # Hide color legend since it duplicates fill
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
print(p)

ggsave(
  filename = "anelorda/plots/miRNA_barplot_with_dots.jpeg",
  plot = p + theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  ),
  device = "jpeg",
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nPlot saved to: anelorda/plots/miRNA_barplot_with_dots.jpeg\n")
