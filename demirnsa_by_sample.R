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


setwd("~/../../media/dell/DATA1/anelorda/")

data <- read.table("deseq2_mirna_13/deseq2_input_13control.csv", header = TRUE, sep = ",")

rownames(data) <- data$miRNA

data <- data[, -1]

data <- data +1
sample_vector <- c("control", "control", "control", "control", "control", "control", "control", "control", "control", 
                   "control", "control", "control", "control", "X001T", "X001T", "X001T", "X002T", "X002T", 
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


# Perform variance stabilizing transformation (force direct transformation)
vsd_mirna <- varianceStabilizingTransformation(dds, blind = TRUE)
#dds_collapsed <- collapseReplicates(dds, groupby = sample_info$sample, run = rownames(sample_info))


# Extract PCA components
pca_data <- plotPCA(vsd_mirna, intgroup = "sample", returnData = TRUE)




# Calculate percentage variance explained by each PC
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Generate PCA plot with ggplot2
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size = 1.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle(" miRNA") +
  theme(legend.position = "right")

ggsave("deseq2_mirna_13/PCA_miRNA.pdf", plot = pca_plot, width = 7, height = 5, dpi = 300)

# Select the top 500 most variable genes
topVarGenes_mirna <- head(order(rowVars(assay(vsd_mirna)), decreasing = TRUE), 300)

# Subset VSD data
mat_mirna <- assay(vsd_mirna)[topVarGenes_mirna, ]

# Scale rows (z-score normalization)
mat_mirna <- t(scale(t(mat_mirna)))


heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)


# Open PDF device
pdf("Heatmap_miRNA.pdf", width = 10, height = 8)

# Generate heatmap without default title
heatmap_mirna <- pheatmap(mat_mirna,
                          scale = "row",
                          clustering_distance_rows = "correlation",
                          clustering_distance_cols = "correlation",
                          clustering_method = "ward.D2",
                          color = heatmap_colors,
                          main = "",  # Remove default title
                          fontsize = 10,
                          show_rownames = FALSE)

# Draw the heatmap
grid::grid.draw(heatmap_mirna$gtable)

# Add a left-aligned title inside the PDF
grid.text("B.   miRNA", x = unit(0, "npc") + unit(1, "mm"), 
          y = unit(1, "npc") - unit(2, "mm"), 
          just = "left", gp = gpar(fontsize = 12))

# Close PDF device
dev.off()

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
  res_df$Significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 2, "Significant", "Not Significant")
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    theme_minimal() +
    labs(title = sample_name,
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted p-value") +
    theme(legend.position = "none", plot.title = element_text(size = 8))
}

# Create PDF
pdf("volcano_plots_mirna.pdf", width = 12, height = 12)

# Store plots
mirna_plots <- list()

# Generate plots for each sample
for (sample in setdiff(unique(sample_info$sample), "control")) {
  res_mirna <- results(dds, contrast = c("sample", sample, "control"))
  mirna_plots[[length(mirna_plots) + 1]] <- plot_volcano(as.data.frame(res_mirna), sample)
}

# Print 16 miRNA volcano plots (Page 1)
grid.arrange(grobs = mirna_plots, ncol = 4, nrow = 4)

dev.off()  # Close the PDF file


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
logfc_mirna$average <- rowMeans(logfc_mirna, na.rm = TRUE)


logfc_mrna <- read.csv("deseq_results/mrna/logfc_mrna.csv", row.names = TRUE, header = TRUE)

common_mrna <- rownames(logfc_mrna)
common_mirna <- rownames(logfc_mirna)

###MIRTARBASE STRONG EVIDENCE
mirtarbase_strong_evidence <- read.table("mirtarbase/miRTarBase_human.csv", sep = ",", header = TRUE, fill = TRUE, quote = "")
head(mirtarbase_strong_evidence)



filtered_mirtarbase <- mirtarbase_strong_evidence[
  mirtarbase_strong_evidence$miRNA %in% common_mirna & 
    mirtarbase_strong_evidence$Target.Gene %in% common_mrna, 
]

head(filtered_mirtarbase)


filtered_unique_pairs <- unique(filtered_mirtarbase[, c("miRNA", "Target.Gene")])

head(filtered_unique_pairs)





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

# View results
head(filtered_unique_pairs)

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

