#Author Anel Ordabayeva
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

# ============================================================
# DATA LOADING & PREPROCESSING
# ============================================================
data <- read.table(
  "ksarkytbayev/nextflow_24_04_25/mirna_quant/edger_qc/mature_counts.csv",
  header = TRUE, sep = ","
)
rownames(data) <- data$X
data$X <- NULL
data <- t(data)
data <- data + 1

# Reorder columns: controls first, then samples (sorted)
control_cols <- sort(grep("^control_", colnames(data), value = TRUE))
sample_cols  <- sort(setdiff(colnames(data), control_cols))
data <- data[, c(control_cols, sample_cols)]

# ============================================================
# SAMPLE METADATA
# ============================================================
sample_vector <- c(
  "control", "control", "control", "control",
  "X001T", "X001T", "X001T", "X002T", "X002T", "X002T",
  "X003T", "X003T", "X003T", "X004T", "X004T", "X004T",
  "X005T", "X005T", "X005T", "X007T", "X007T", "X007T",
  "X008T", "X008T", "X008T", "X011T", "X011T", "X011T",
  "X012T", "X012T", "X012T", "X015T", "X015T", "X015T",
  "X016T", "X016T", "X016T", "X017T", "X017T", "X017T",
  "X023T", "X023T", "X023T", "X025T", "X025T",
  "X029T", "X029T", "X029T", "X033T", "X033T", "X033T"
)
stopifnot(length(sample_vector) == ncol(data))

row_names   <- colnames(data)
sample_info <- data.frame(sample = sample_vector, row.names = row_names)
sample_info$sample <- factor(sample_info$sample, levels = unique(sample_info$sample))
sample_info$sample <- relevel(sample_info$sample, ref = "control")

# ============================================================
# FULL MODEL: per-sample DESeq2 (used for DEG extraction)
# ============================================================
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData   = sample_info,
  design    = ~ sample
)
dds <- DESeq(dds)

normalized_counts_mirna <- counts(dds, normalized = TRUE)
write.csv(
  normalized_counts_mirna,
  "deseq2_mirna_13/normalised_counts_miRNA.csv",
  row.names = TRUE
)

# ============================================================
# COLLAPSED MODEL: for PCA & heatmap
# ============================================================
dds_raw     <- DESeqDataSetFromMatrix(
  countData = data,
  colData   = sample_info,
  design    = ~ sample
)

control_idx        <- which(sample_info$sample == "control")
group_for_collapse <- as.character(sample_info$sample)
group_for_collapse[control_idx] <- paste0("control_", seq_along(control_idx))
group_for_collapse <- factor(group_for_collapse)
table(group_for_collapse)

dds_collapsed <- collapseReplicates(
  dds_raw,
  groupby = group_for_collapse,
  run     = colnames(dds_raw)
)

colData(dds_collapsed)$condition <- ifelse(
  grepl("^control", colnames(dds_collapsed)), "control", "tumor"
)
colData(dds_collapsed)$condition <- factor(
  colData(dds_collapsed)$condition, levels = c("control", "tumor")
)
design(dds_collapsed) <- ~ condition
dds_collapsed <- DESeq(dds_collapsed)

normalized_counts_mirna <- counts(dds_collapsed, normalized = TRUE)
write.csv(
  normalized_counts_mirna,
  "deseq2_mirna_13/normalised_counts_miRNA.csv",
  row.names = TRUE
)

# ============================================================
# PCA
# ============================================================
vsd_mirna <- varianceStabilizingTransformation(dds_collapsed, blind = TRUE)
pca_data  <- plotPCA(vsd_mirna, intgroup = "sample", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size = 1.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(expression(bold("f  ") ~ "    miRNA")) +
  theme_bw() +
  theme(
    legend.position   = "right",
    plot.title        = element_text(hjust = 0),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )

ggsave("anelorda/new_controls_deseq/pca_collapsed.tiff",
       pca_plot, width = 7, height = 5, dpi = 300,
       compression = "lzw", bg = "white")
ggsave("anelorda/new_controls_deseq/pca_collapsed.jpeg",
       pca_plot, width = 7, height = 5, dpi = 300,
       quality = 100, bg = "white")

# ============================================================
# HEATMAP
# ============================================================
topVarGenes_mirna <- head(order(rowVars(assay(vsd_mirna)), decreasing = TRUE), 300)
mat_mirna <- assay(vsd_mirna)[topVarGenes_mirna, ]
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
    scale                    = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method        = "ward.D2",
    color                    = heatmap_colors,
    main                     = "",
    fontsize                 = 7,
    show_rownames            = FALSE
  )
  grid::grid.draw(hm$gtable)
  grid.text(
    expression(bold("i") * "          miRNA"),
    x    = unit(0, "npc") + unit(3, "mm"),
    y    = unit(1, "npc") - unit(6, "mm"),
    just = "left",
    gp   = gpar(fontsize = 13)
  )
  dev.off()
}

save_heatmap("new_controls_deseq/Heatmap_miRNA.tiff", "tiff")
save_heatmap("new_controls_deseq/Heatmap_miRNA.jpeg", "jpeg")

# ============================================================
# DEG EXTRACTION (full model)
# ============================================================
setwd("anelorda/new_controls_deseq/")

extract_results <- function(dds_obj, sample_name, data_type) {
  res    <- results(dds_obj, contrast = c("sample", sample_name, "control"))
  res_df <- as.data.frame(res)
  
  write.csv(res_df,
            paste0("res_", data_type, "_", sample_name, "_vs_control_unfiltered.csv"),
            row.names = TRUE)
  assign(paste0("res_", data_type, "_", sample_name, "_vs_control_unfiltered"),
         res_df, envir = .GlobalEnv)
  
  res_filtered <- res_df[which(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05), ]
  
  write.csv(res_filtered,
            paste0("res_", data_type, "_", sample_name, "_vs_control_filtered.csv"),
            row.names = TRUE)
  assign(paste0("res_", data_type, "_", sample_name, "_vs_control_filtered"),
         res_filtered, envir = .GlobalEnv)
  
  return(list(filtered = res_filtered, unfiltered = res_df))
}

for (sample in setdiff(unique(sample_info$sample), "control")) {
  extract_results(dds, sample, "mirna")
}

# ============================================================
# VOLCANO PLOTS
# ============================================================
plot_volcano <- function(res_df, sample_name) {
  res_df$Significance <- ifelse(
    res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
    "Significant", "Not Significant"
  )
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    theme_minimal() +
    labs(title = sample_name, x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    theme(legend.position = "none", plot.title = element_text(size = 8))
}

jpeg("volcano_plots_mirna.jpeg", width = 12, height = 12, units = "in", res = 400)
mirna_plots <- list()
for (sample in setdiff(unique(sample_info$sample), "control")) {
  res_mirna <- results(dds, contrast = c("sample", sample, "control"))
  mirna_plots[[length(mirna_plots) + 1]] <- plot_volcano(as.data.frame(res_mirna), sample)
}
grid.arrange(grobs = mirna_plots, ncol = 4, nrow = 4)
dev.off()

# ============================================================
# COMMON DEG SUMMARY (log2FC across all samples)
# ============================================================
sample_names <- unique(sample_vector[sample_vector != "control"])

deg_mirna <- lapply(sample_names, function(s) {
  get(paste0("res_mirna_", s, "_vs_control_filtered"))
})
names(deg_mirna) <- sample_names

common_mirna <- Reduce(intersect, lapply(deg_mirna, rownames))

logfc_mirna <- do.call(cbind, lapply(deg_mirna, function(df) {
  df[common_mirna, "log2FoldChange", drop = FALSE]
}))
colnames(logfc_mirna) <- sample_names
rownames(logfc_mirna) <- common_mirna

logfc_mirna <- logfc_mirna[
  apply(logfc_mirna, 1, function(x) all(x >= 0) | all(x <= 0)), 
]
logfc_mirna$average     <- rowMeans(logfc_mirna, na.rm = TRUE)
rownames(logfc_mirna)   <- gsub("\\.", "-", rownames(logfc_mirna))

# ============================================================
# ============================================================
# LEAVE-ONE-OUT (LOO) SENSITIVITY ANALYSIS
# Addresses reviewer concern about small control group (n = 4)
# ============================================================
# ============================================================

cat("\n")
cat("============================================================\n")
cat("  LEAVE-ONE-OUT SENSITIVITY ANALYSIS (n controls = 4)\n")
cat("============================================================\n")

# Identify the 4 control column indices
control_indices <- which(sample_vector == "control")
control_names   <- row_names[control_indices]
cat("Control samples being tested:\n")
print(control_names)

# Helper: collect all tumor-vs-control unfiltered results from a dds object
get_all_unfiltered <- function(dds_obj, s_info) {
  tumors <- setdiff(levels(s_info$sample), "control")
  res_list <- lapply(tumors, function(s) {
    res <- tryCatch(
      results(dds_obj, contrast = c("sample", s, "control")),
      error = function(e) NULL
    )
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    df$gene  <- rownames(df)
    df$tumor <- s
    df
  })
  bind_rows(Filter(Negate(is.null), res_list))
}

# Full-model results for LOO comparison
results_full <- get_all_unfiltered(dds, sample_info)

# ---- Run LOO iterations ----------------------------------------
loo_results <- list()

for (i in seq_along(control_indices)) {
  
  removed_name <- control_names[i]
  cat("\n--- LOO", i, ": removing", removed_name, "---\n")
  
  keep_cols    <- setdiff(seq_len(ncol(data)), control_indices[i])
  data_loo     <- data[, keep_cols]
  info_loo     <- sample_info[keep_cols, , drop = FALSE]
  info_loo$sample <- droplevels(info_loo$sample)
  info_loo$sample <- relevel(info_loo$sample, ref = "control")
  
  dds_loo <- DESeqDataSetFromMatrix(
    countData = data_loo,
    colData   = info_loo,
    design    = ~ sample
  )
  dds_loo <- DESeq(dds_loo)
  
  res_loo              <- get_all_unfiltered(dds_loo, info_loo)
  res_loo$loo_removed  <- removed_name
  loo_results[[i]]     <- res_loo
}

loo_combined <- bind_rows(loo_results)

# ---- Stability metrics -----------------------------------------

# Merge full vs LOO on gene + tumor
full_slim <- results_full %>%
  select(gene, tumor, log2FoldChange, padj) %>%
  rename(log2FC_full = log2FoldChange, padj_full = padj)

loo_slim <- loo_combined %>%
  select(gene, tumor, log2FoldChange, padj, loo_removed) %>%
  rename(log2FC_loo = log2FoldChange, padj_loo = padj)

stability_df <- left_join(loo_slim, full_slim, by = c("gene", "tumor"))

# Pearson / Spearman correlation of log2FC
cor_summary <- stability_df %>%
  group_by(loo_removed) %>%
  summarise(
    pearson_r        = cor(log2FC_full, log2FC_loo, use = "complete.obs"),
    spearman_r       = cor(log2FC_full, log2FC_loo, use = "complete.obs", method = "spearman"),
    median_abs_delta = median(abs(log2FC_full - log2FC_loo), na.rm = TRUE),
    .groups = "drop"
  )

# DEG overlap: genes significant in full model still significant in LOO
sig_full <- results_full %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  select(gene, tumor) %>%
  distinct()

overlap_summary <- loo_combined %>%
  filter(!is.na(padj)) %>%
  inner_join(sig_full, by = c("gene", "tumor")) %>%
  group_by(loo_removed) %>%
  summarise(
    n_sig_full   = nrow(sig_full),
    n_still_sig  = sum(padj < 0.05 & abs(log2FoldChange) > 1, na.rm = TRUE),
    pct_retained = round(100 * n_still_sig / n_sig_full, 1),
    .groups = "drop"
  )

loo_table <- left_join(cor_summary, overlap_summary, by = "loo_removed")

cat("\n=== LOO SUMMARY TABLE ===\n")
print(loo_table)
write.csv(loo_table, "LOO_summary_table.csv", row.names = FALSE)

# ---- LOO Visualisation -----------------------------------------

# 1. Scatter: full log2FC vs LOO log2FC (one panel per removed control)
p_scatter <- ggplot(
  stability_df %>% filter(!is.na(log2FC_full), !is.na(log2FC_loo)),
  aes(x = log2FC_full, y = log2FC_loo)
) +
  geom_point(alpha = 0.25, size = 0.6, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ loo_removed, ncol = 2) +
  labs(
    title    = "Leave-One-Out Sensitivity: log\u2082 Fold Changes",
    x        = "log\u2082FC  (full model, n = 4 controls)",
    y        = "log\u2082FC  (LOO model, n = 3 controls)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave("LOO_log2FC_scatter.tiff", p_scatter,
       width = 7, height = 6, dpi = 300, compression = "lzw")
ggsave("LOO_log2FC_scatter.jpeg", p_scatter,
       width = 7, height = 6, dpi = 300, quality = 95)

# 2. Bar chart: Pearson r per LOO iteration
p_corr <- ggplot(cor_summary, aes(x = loo_removed, y = pearson_r, fill = loo_removed)) +
  geom_col(width = 0.55, show.legend = FALSE, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("r = %.4f", pearson_r)),
            vjust = -0.5, size = 3.2) +
  scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Pearson Correlation of log\u2082FC: Full vs LOO Models",
    x     = "Removed Control Sample",
    y     = "Pearson r"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("LOO_pearson_correlation.tiff", p_corr,
       width = 5, height = 4.5, dpi = 300, compression = "lzw")
ggsave("LOO_pearson_correlation.jpeg", p_corr,
       width = 5, height = 4.5, dpi = 300, quality = 95)

# 3. DEG retention bar chart
p_overlap <- ggplot(overlap_summary,
                    aes(x = loo_removed, y = pct_retained, fill = loo_removed)) +
  geom_col(width = 0.55, show.legend = FALSE, color = "grey30", linewidth = 0.3) +
  geom_text(
    aes(label = paste0(pct_retained, "%\n(", n_still_sig, "/", n_sig_full, ")")),
    vjust = -0.3, size = 3.0, lineheight = 0.9
  ) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title    = "DEG Retention Across LOO Iterations",
    subtitle = "% of DEGs (full model: |log\u2082FC| > 1, padj < 0.05) retained after removing each control",
    x        = "Removed Control Sample",
    y        = "DEGs Retained (%)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x   = element_text(angle = 30, hjust = 1),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

ggsave("LOO_DEG_retention.tiff", p_overlap,
       width = 5, height = 4.5, dpi = 300, compression = "lzw")
ggsave("LOO_DEG_retention.jpeg", p_overlap,
       width = 5, height = 4.5, dpi = 300, quality = 95)

cat("\n[LOO analysis complete]\n")
cat("Output files saved to current working directory:\n")
cat("  LOO_summary_table.csv\n")
cat("  LOO_log2FC_scatter.tiff/.jpeg\n")
cat("  LOO_pearson_correlation.tiff/.jpeg\n")
cat("  LOO_DEG_retention.tiff/.jpeg\n\n")

# ============================================================
# MIRTARBASE INTEGRATION
# ============================================================
logfc_mrna <- read.csv("../deseq_results/mrna/logfc_mrna.csv",
                       row.names = 1, header = TRUE)

common_mrna  <- rownames(logfc_mrna)
common_mirna <- gsub("\\.", "-", rownames(logfc_mirna))

mirtarbase_strong_evidence <- read.table(
  "../mirtarbase/miRTarBase_human.csv",
  sep = ",", header = TRUE, fill = TRUE, quote = ""
)

filtered_mirtarbase <- mirtarbase_strong_evidence[
  mirtarbase_strong_evidence$miRNA %in% common_mirna &
    mirtarbase_strong_evidence$Target.Gene %in% common_mrna, 
]
filtered_unique_pairs <- unique(filtered_mirtarbase[, c("miRNA", "Target.Gene")])

mirtarbase <- read.table(
  "../mirtarbase/hsa_MTI.csv",
  sep = ",", header = TRUE, fill = TRUE, quote = ""
)
filtered_mirtarbase_all <- mirtarbase[
  mirtarbase$miRNA %in% common_mirna &
    mirtarbase$Target.Gene %in% common_mrna, 
]
filtered_mirtarbase_all <- unique(filtered_mirtarbase_all[, c("miRNA", "Target.Gene")])

# Merge average log2FC values
logfc_mirna$miRNA      <- rownames(logfc_mirna)
logfc_mrna$Target.Gene <- rownames(logfc_mrna)
logfc_mirna_avg        <- logfc_mirna[, c("miRNA", "average")]
logfc_mrna_avg         <- logfc_mrna[, c("Target.Gene", "average")]

filtered_unique_pairs <- merge(filtered_unique_pairs, logfc_mirna_avg, by = "miRNA",       all.x = TRUE)
filtered_unique_pairs <- merge(filtered_unique_pairs, logfc_mrna_avg,  by = "Target.Gene", all.x = TRUE)
colnames(filtered_unique_pairs)[3:4] <- c("logFC_miRNA", "logFC_mRNA")

filtered_mirtarbase_all <- merge(filtered_mirtarbase_all, logfc_mirna_avg, by = "miRNA",       all.x = TRUE)
filtered_mirtarbase_all <- merge(filtered_mirtarbase_all, logfc_mrna_avg,  by = "Target.Gene", all.x = TRUE)
colnames(filtered_mirtarbase_all)[3:4] <- c("logFC_miRNA", "logFC_mRNA")

# Merge with lncRNA interaction data
filtered_unique_pairs_lncRNA <- merge(
  filtered_unique_pairs,
  filtered_npinter[, c("Target_miRNA", "lncRNA_Name")],
  by.x  = "miRNA",
  by.y  = "Target_miRNA",
  all.x = TRUE
)
filtered_unique_pairs_lncRNA <- filtered_unique_pairs_lncRNA[
  !is.na(filtered_unique_pairs_lncRNA$lncRNA_Name), 
]

filtered_unique_pairs_pos_mirna <- filtered_unique_pairs[
  filtered_unique_pairs$logFC_miRNA > 0 & filtered_unique_pairs$logFC_mRNA < 0, 
]
filtered_unique_pairs_neg_mirna <- filtered_unique_pairs[
  filtered_unique_pairs$logFC_miRNA < 0 & filtered_unique_pairs$logFC_mRNA > 0, 
]

# Final combined mRNA + lncRNA target data frame
filtered_df_rearranged        <- filtered_df[, c("miRNA", "Target.Gene", "logFC_miRNA", "logFC_mRNA")]
colnames(filtered_df_rearranged)[2] <- "Target"

lncrna_df <- filtered_df[, c("miRNA", "lncRNA_Name", "logFC_miRNA", "logFC_lncRNA")]
colnames(lncrna_df) <- c("miRNA", "Target", "logFC_miRNA", "logFC_mRNA")

final_filtered_df <- rbind(filtered_df_rearranged, lncrna_df)
final_filtered_df <- final_filtered_df[!is.na(final_filtered_df$Target), ]
final_filtered_df <- final_filtered_df %>% distinct()

head(final_filtered_df)