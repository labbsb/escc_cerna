library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)
library(dplyr)      
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(stringr)
library(gridExtra)
library(grid)
library(SummarizedExperiment)
library(matrixStats)
setwd("/mnt/sdb2/")

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

sample_vector <- c("control", "control", "control", "control", "X001T", "X001T", "X001T", "X002T", "X002T", 
                   "X002T", "X003T", "X003T", "X003T", "X004T", "X004T", "X004T", "X005T", "X005T", 
                   "X005T", "X007T", "X007T", "X007T", "X008T", "X008T", "X008T", "X011T", "X011T", 
                   "X011T", "X012T", "X012T", "X012T", "X015T", "X015T", "X015T", "X016T", "X016T", 
                   "X016T", "X017T", "X017T", "X017T", "X023T", "X023T", "X023T", "X025T", "X025T", 
                   "X029T", "X029T", "X029T", "X033T", "X033T", "X033T"
)


row_names <- colnames(data)
sample_info <- data.frame(sample = sample_vector, row.names = row_names)

# Convert sample_info$sample into a factor with all levels
sample_info$sample <- factor(sample_info$sample, levels = unique(sample_info$sample))

# Set control as the reference level
sample_info$sample <- relevel(sample_info$sample, ref = "control")

# Now, create the DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = data, colData = sample_info, design = ~ sample)
dds <- DESeq(dds)



normalized_counts_mirna <- counts(dds, normalized = TRUE)
write.csv(normalized_counts_mirna, "deseq2_mirna_13/normalised_counts_miRNA.csv", row.names = TRUE)

stopifnot(length(sample_vector) == ncol(data))

sample_info <- data.frame(
  sample = sample_vector,
  row.names = colnames(data)
)

# =========================
# Create DESeqDataSet (RAW)
# =========================
dds_raw <- DESeqDataSetFromMatrix(
  countData = data,
  colData   = sample_info,
  design    = ~ sample
)

# =========================
# Collapse replicates
#   - tumors collapsed
#   - controls kept separate
# =========================
control_idx <- which(sample_info$sample == "control")

group_for_collapse <- as.character(sample_info$sample)
group_for_collapse[control_idx] <- paste0(
  "control_", seq_along(control_idx)
)
group_for_collapse <- factor(group_for_collapse)

# Sanity check
table(group_for_collapse)

dds_collapsed <- collapseReplicates(
  dds_raw,
  groupby = group_for_collapse,
  run     = colnames(dds_raw)
)

# =========================
# Clean design variable
# =========================
colData(dds_collapsed)$condition <- ifelse(
  grepl("^control", colnames(dds_collapsed)),
  "control",
  "tumor"
)

colData(dds_collapsed)$condition <- factor(
  colData(dds_collapsed)$condition,
  levels = c("control", "tumor")
)

design(dds_collapsed) <- ~ condition

# =========================
# Run DESeq
# =========================
dds_collapsed <- DESeq(dds_collapsed)

# =========================
# Normalized counts
# =========================
normalized_counts_mirna <- counts(dds_collapsed, normalized = TRUE)

write.csv(
  normalized_counts_mirna,
  "deseq2_mirna_13/normalised_counts_miRNA.csv",
  row.names = TRUE
)

# =========================
# VST + PCA
# =========================
vsd_mirna <- varianceStabilizingTransformation(
  dds_collapsed,
  blind = TRUE
)

# Extract PCA components
pca_data <- plotPCA(vsd_mirna, intgroup = "sample", returnData = TRUE)




# Calculate percentage variance explained by each PC
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Generate PCA plot with ggplot2
# Generate PCA plot with ggplot2
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size = 1.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(expression(bold("f  ") ~ "    miRNA")) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )


ggsave(
  "anelorda/new_controls_deseq/pca_collapsed.tiff",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300,
  compression = "lzw",
  bg = "white"
)
ggsave(
  "anelorda/new_controls_deseq/pca_collapsed.jpeg",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300,
  quality = 100,
  bg = "white"
)



# Select the top 500 most variable genes
topVarGenes_mirna <- head(order(rowVars(assay(vsd_mirna)), decreasing = TRUE), 300)

# Subset VSD data
mat_mirna <- assay(vsd_mirna)[topVarGenes_mirna, ]

# Scale rows (z-score normalization)
mat_mirna <- t(scale(t(mat_mirna)))


heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

save_heatmap <- function(filename, device = c("tiff", "jpeg")) {
  device <- match.arg(device)
  
  if (device == "tiff") {
    tiff(filename, width = 12, height = 9, units = "in",
         res = 400, compression = "lzw", bg = "white")
  } else {
    jpeg(filename, width = 7, height = 5, units = "in",
         res = 400, quality = 100, bg = "white")
  }
  
  hm <- pheatmap(
    mat_mirna,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    color = heatmap_colors,
    main = "",
    fontsize = 7,
    show_rownames = FALSE
  )
  
  grid::grid.draw(hm$gtable)
  
  grid.text(
    expression(bold("i") * "          miRNA"),
    x = unit(0, "npc") + unit(3, "mm"),
    y = unit(1, "npc") - unit(6, "mm"),  # move DOWN
    just = "left",
    gp = gpar(fontsize = 13),
  )
  
  
  dev.off()
}
save_heatmap("new_controls_deseq/Heatmap_miRNA.tiff", "tiff")
save_heatmap("new_controls_deseq/Heatmap_miRNA.jpeg", "jpeg")

setwd('anelorda/new_controls_deseq/')
# Define the function to extract and save results
extract_results <- function(dds, sample_name, data_type) {
  # Get the results for the comparison: sample vs control
  res <- results(dds, contrast = c("sample", sample_name, "control"))
  
  # Convert to a data frame
  res_df <- as.data.frame(res)
  
  # Save unfiltered results
  unfiltered_filename <- paste0("res_", data_type, "_", sample_name, "_vs_control_unfiltered.csv")
  write.csv(res_df, file = unfiltered_filename, row.names = TRUE)
  
  # Assign unfiltered results to environment
  assign(paste0("res_", data_type, "_", sample_name, "_vs_control_unfiltered"), res_df, envir = .GlobalEnv)
  
  # Filter by log2FoldChange and padj
  res_filtered <- res_df[which(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05), ]
  
  # Save filtered results
  filtered_filename <- paste0("res_", data_type, "_", sample_name, "_vs_control_filtered.csv")
  write.csv(res_filtered, file = filtered_filename, row.names = TRUE)
  
  # Assign filtered results to environment
  assign(paste0("res_", data_type, "_", sample_name, "_vs_control_filtered"), res_filtered, envir = .GlobalEnv)
  
  return(list(filtered = res_filtered, unfiltered = res_df))
}

# Extract miRNA results for each sample (excluding control)
for (sample in setdiff(unique(sample_info$sample), "control")) {
  extract_results(dds, sample, "mirna")
}


# Function to generate a volcano plot
plot_volcano <- function(res_df, sample_name) {
  res_df$Significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    theme_minimal() +
    labs(title = sample_name,
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted p-value") +
    theme(legend.position = "none", plot.title = element_text(size = 8))
}

# Create JPEG (instead of PDF)
jpeg("volcano_plots_mirna.jpeg", width = 12, height = 12, units = "in", res = 400)

# Store plots
mirna_plots <- list()

# Generate plots for each sample
for (sample in setdiff(unique(sample_info$sample), "control")) {
  res_mirna <- results(dds, contrast = c("sample", sample, "control"))
  mirna_plots[[length(mirna_plots) + 1]] <- plot_volcano(as.data.frame(res_mirna), sample)
}

# Print 16 miRNA volcano plots (single image)
grid.arrange(grobs = mirna_plots, ncol = 4, nrow = 4)

dev.off()  # Close the JPEG file


sample_names <- unique(sample_vector[sample_vector != "control"])

# Read all data frames into a named list
deg_mirna <- lapply(sample_names, function(sample) {
  df_name <- paste0("res_mirna_", sample, "_vs_control_filtered")
  get(df_name)  # Retrieve the data frame
})
names(deg_mirna) <- sample_names

# Find common genes across all data frames
common_mirna <- Reduce(intersect, lapply(deg_mirna, rownames))

# Create a new data frame with log2FoldChange values
logfc_mirna <- do.call(cbind, lapply(deg_mirna, function(df) {
  df[common_mirna, "log2FoldChange", drop = FALSE]
}))

# Rename columns with sample names
colnames(logfc_mirna) <- sample_names

# Set row names as gene IDs
rownames(logfc_mirna) <- common_mirna

# View final dataframe
head(logfc_mirna)
logfc_mirna <- logfc_mirna[
  apply(logfc_mirna, 1, function(x) all(x >= 0) | all(x <= 0)), 
]
logfc_mirna$average <- rowMeans(logfc_mirna, na.rm = TRUE)


rownames(logfc_mirna) <- gsub("\\.", "-", rownames(logfc_mirna))


head(logfc_mirna)

logfc_mrna <- read.csv("../deseq_results/mrna/logfc_mrna.csv", row.names = 1, header = TRUE)

common_mrna <- rownames(logfc_mrna)
common_mirna <- rownames(logfc_mirna)

common_mirna <- gsub("\\.", "-", common_mirna)
head(common_mirna)



###MIRTARBASE STRONG EVIDENCE
mirtarbase_strong_evidence <- read.table("../mirtarbase/miRTarBase_human.csv", sep = ",", header = TRUE, fill = TRUE, quote = "")
head(mirtarbase_strong_evidence)



filtered_mirtarbase <- mirtarbase_strong_evidence[
  mirtarbase_strong_evidence$miRNA %in% common_mirna & 
    mirtarbase_strong_evidence$Target.Gene %in% common_mrna, 
]

head(filtered_mirtarbase)


filtered_unique_pairs <- unique(filtered_mirtarbase[, c("miRNA", "Target.Gene")])

head(filtered_unique_pairs)

###MIRTARBASE ALL
mirtarbase <- read.table("../mirtarbase/hsa_MTI.csv", sep =",", header = TRUE, fill=TRUE, quote ="")
head(mirtarbase)

filtered_mirtarbase_all <- mirtarbase[
  mirtarbase$miRNA %in% common_mirna & 
    mirtarbase$Target.Gene %in% common_mrna, 
]


head(filtered_mirtarbase_all)
filtered_mirtarbase_all <- unique(filtered_mirtarbase_all[, c ("miRNA", "Target.Gene")])


# Ensure column names match for merging
logfc_mirna$miRNA <- rownames(logfc_mirna)
logfc_mrna$Target.Gene <- rownames(logfc_mrna)

# Select relevant columns
logfc_mirna_avg <- logfc_mirna[, c("miRNA", "average")]
logfc_mrna_avg <- logfc_mrna[, c("Target.Gene", "average")]

# Merge with filtered_unique_pairs
filtered_unique_pairs <- merge(filtered_unique_pairs, logfc_mirna_avg, by = "miRNA", all.x = TRUE)
filtered_unique_pairs <- merge(filtered_unique_pairs, logfc_mrna_avg, by = "Target.Gene", all.x = TRUE)

# Rename columns for clarity
colnames(filtered_unique_pairs)[3:4] <- c("logFC_miRNA", "logFC_mRNA")

# Merge with filtered_mirtarbase_all
filtered_mirtarbase_all <- merge(filtered_mirtarbase_all, logfc_mirna_avg, by = "miRNA", all.x = TRUE)
filtered_mirtarbase_all <- merge(filtered_mirtarbase_all, logfc_mrna_avg, by = "Target.Gene", all.x = TRUE)

# Rename columns for clarity
colnames(filtered_mirtarbase_all)[3:4] <- c("logFC_miRNA", "logFC_mRNA")
# View results
head(filtered_mirtarbase_all)

# Merge filtered_unique_pairs with filtered_npinter based on miRNA
filtered_unique_pairs_lncRNA <- merge(
  filtered_unique_pairs, 
  filtered_npinter[, c("Target_miRNA", "lncRNA_Name")], 
  by.x = "miRNA", 
  by.y = "Target_miRNA", 
  all.x = TRUE
)

# View the updated dataset
head(filtered_unique_pairs_lncRNA)

filtered_unique_pairs_lncRNA <- filtered_unique_pairs_lncRNA[!is.na(filtered_unique_pairs_lncRNA$lncRNA_Name), ]
filtered_unique_pairs_pos_mirna <- filtered_unique_pairs[filtered_unique_pairs$logFC_miRNA > 0 & filtered_unique_pairs$logFC_mRNA < 0, ]
head(filtered_unique_pairs_pos_mirna)
filtered_unique_pairs_neg_mirna <- filtered_unique_pairs[filtered_unique_pairs$logFC_miRNA < 0 & filtered_unique_pairs$logFC_mRNA > 0, ]

# View the first few rows
head(filtered_unique_pairs_lncRNA)

# Create a new data frame where lncRNAs are also treated as target genes
filtered_df_rearranged <- filtered_df[, c("miRNA", "Target.Gene", "logFC_miRNA", "logFC_mRNA")]
colnames(filtered_df_rearranged)[2] <- "Target"

lncrna_df <- filtered_df[, c("miRNA", "lncRNA_Name", "logFC_miRNA", "logFC_lncRNA")]
colnames(lncrna_df) <- c("miRNA", "Target", "logFC_miRNA", "logFC_mRNA")

# Combine both data frames
final_filtered_df <- rbind(filtered_df_rearranged, lncrna_df)

# Remove rows where lncRNA_Name was NA
final_filtered_df <- final_filtered_df[!is.na(final_filtered_df$Target), ]
final_filtered_df <- final_filtered_df %>% distinct()

# Display the head of the new data frame
head(final_filtered_df)
