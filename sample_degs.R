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


setwd("/mnt/sdb2/anelorda/")

# Load Data
data_degs <- read.table("degs_mrna/mRNA_for_miRNA.xlsx - sheet.csv", header = TRUE, sep = ",")

# Remove unwanted samples
samples_to_remove <- c("X009T", "X010T", "X013T", "X014T", "X019T", "X021T", "Normal_SRR7961236")
count_data <- data_degs[, !(colnames(data_degs) %in% samples_to_remove)]

# Rename first column to gene identifiers (Ensembl IDs)
colnames(count_data)[1] <- "geneids"
rownames(count_data) <- count_data$geneids
count_data <- count_data[, -1]  # Remove first column after setting rownames

# Filter low-expression genes
count_data <- count_data[rowSums(count_data) > 5, ]

gencode_lncrna <- read.table("gencode.v47.long_noncoding_RNAs.gtf", header = FALSE, sep = "\t")

lncrna_genes <- gencode_lncrna %>%
  filter(V3 == "gene") %>%
  mutate(
    ensembl_gene_id = str_extract(V9, "gene_id ([^;]+)") %>% str_remove("gene_id "),
    gene_name = str_extract(V9, "gene_name ([^;]+)") %>% str_remove("gene_name ")
  ) %>%
  select(ensembl_gene_id, gene_name)
# Ensure rownames of count_data are available
lncrna_ids <- gsub("\\..*", "", lncrna_genes$ensembl_gene_id) # Remove version numbers

# Separate lncRNA counts
# count_lncrna <- count_data[rownames(count_data) %in% lncrna_ids, ]

# Separate mRNA counts (everything not in lncRNA)
count_mrna <- count_data[!(rownames(count_data) %in% lncrna_ids), ]

# Check results
dim(count_lncrna) # Number of lncRNA genes
dim(count_mrna)   # Number of mRNA genes

# Define sample groups
sample_vector <- c('control', 'control', 'control', 'control', 'control', 'control', 
                   'control', 'control', 'control', 'control', 'X001T', 'X002T', 
                   'X003T', 'X004T', 'X005T', 'X007T', 'X008T', 
                   'X011T', 'X012T', 'X015T', 'X016T', 'X017T', 'X023T', 'X025T', 
                   'X029T', 'X033T')

sample_info <- data.frame(sample = sample_vector, row.names = colnames(count_data))

dds_lncrna <- DESeqDataSetFromMatrix(countData = count_lncrna,
                                     colData = sample_info,
                                     design = ~ sample)

dds_lncrna <- DESeq(dds_lncrna)
sample_names <- sample_vector[sample_vector != "control"]

dds_mrna <- DESeqDataSetFromMatrix(countData = count_mrna,
                                   colData = sample_info,
                                   design = ~ sample)

dds_mrna <- DESeq(dds_mrna)

vsd_mrna <- varianceStabilizingTransformation(dds_mrna, blind = TRUE)
vsd_lncrna <- varianceStabilizingTransformation(dds_lncrna, blind = TRUE)
library(ggplot2)

# Extract PCA components
pca_mrna <- plotPCA(vsd_mrna, intgroup = "sample", returnData = TRUE)

# Calculate percentage variance
percentVar_mrna <- round(100 * attr(pca_mrna, "percentVar"))

# Generate PCA plot for mRNA
pca_plot_mrna <- ggplot(pca_mrna, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size = 1.5) +
  xlab(paste0("PC1: ", percentVar_mrna[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_mrna[2], "% variance")) +
  theme_minimal() +
  ggtitle("A. mRNA") +
  theme(legend.position = "right")

# Save plot
ggsave("PCA_mRNA.pdf", plot = pca_plot_mrna, width = 7, height = 5, dpi = 300)

# Extract PCA components
pca_lncrna <- plotPCA(vsd_lncrna, intgroup = "sample", returnData = TRUE)

# Calculate percentage variance
percentVar_lncrna <- round(100 * attr(pca_lncrna, "percentVar"))

# Generate PCA plot for lncRNA
pca_plot_lncrna <- ggplot(pca_lncrna, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size = 1.5) +
  xlab(paste0("PC1: ", percentVar_lncrna[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_lncrna[2], "% variance")) +
  theme_minimal() +
  ggtitle("C. lncRNA") +
  theme(legend.position = "right")

# Save plot
ggsave("PCA_lncRNA.pdf", plot = pca_plot_lncrna, width = 7, height = 5, dpi = 300)


# Select the top 500 most variable genes
topVarGenes_mrna <- head(order(rowVars(assay(vsd_mrna)), decreasing = TRUE), 500)
topVarGenes_lncrna <- head(order(rowVars(assay(vsd_lncrna)), decreasing = TRUE), 500)

# Subset VSD data
mat_mrna <- assay(vsd_mrna)[topVarGenes_mrna, ]
mat_lncrna <- assay(vsd_lncrna)[topVarGenes_lncrna, ]

# Scale rows (z-score normalization)
mat_mrna <- t(scale(t(mat_mrna)))
mat_lncrna <- t(scale(t(mat_lncrna)))


# Define color palette
heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)


# Function to generate and save heatmap with left-aligned title
save_heatmap <- function(mat, filename, title) {
  pdf(filename, width = 10, height = 8)
  
  # Generate heatmap without default title
  heatmap_obj <- pheatmap(mat,
                          scale = "row",
                          clustering_distance_rows = "correlation",
                          clustering_distance_cols = "correlation",
                          clustering_method = "ward.D2",
                          color = heatmap_colors,
                          main = "",  # Remove default title
                          fontsize = 10,
                          show_rownames = FALSE)
  
  # Draw the heatmap
  grid::grid.draw(heatmap_obj$gtable)
  
  # Add a left-aligned title inside the PDF
  grid.text(title, x = unit(0, "npc") + unit(1, "mm"), 
            y = unit(1, "npc") - unit(2, "mm"), 
            just = "left", gp = gpar(fontsize = 12))
  
  # Close PDF device
  dev.off()
}

# Save mRNA heatmap with left-aligned title
save_heatmap(mat_mrna, "Heatmap_mRNA.pdf", "A.   mRNA")

# Save lncRNA heatmap with left-aligned title
save_heatmap(mat_lncrna, "Heatmap_lncRNA.pdf", "C.   lncRNA")

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

for (sample in setdiff(unique(sample_info$sample), "control")) {
  extract_results(dds_lncrna, sample, "lncrna")
}


for (sample in setdiff(unique(sample_info$sample), "control")) {
  extract_results(dds_mrna, sample, "mrna")
}


# Function to generate a volcano plot
plot_volcano <- function(res_df, sample_name, data_type) {
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

# Create a single PDF file
pdf("volcano_plots_lncrna_mrna.pdf", width = 12, height = 12)

# Store plots separately for lncRNA and mRNA
lncrna_plots <- list()
mrna_plots <- list()

# Loop through each sample and create plots
for (sample in setdiff(unique(sample_info$sample), "control")) {
  # Load lncRNA results
  res_lncrna <- get(paste0("res_lncrna_", sample, "_vs_control_unfiltered"))
  lncrna_plots[[length(lncrna_plots) + 1]] <- plot_volcano(res_lncrna, sample, "lncRNA")
  
  # Load mRNA results
  res_mrna <- get(paste0("res_mrna_", sample, "_vs_control_unfiltered"))
  mrna_plots[[length(mrna_plots) + 1]] <- plot_volcano(res_mrna, sample, "mRNA")
}

# Print 16 lncRNA volcano plots on Page 1
grid.arrange(grobs = lncrna_plots, ncol = 4, nrow = 4)

# Print 16 mRNA volcano plots on Page 2
grid.arrange(grobs = mrna_plots, ncol = 4, nrow = 4)

dev.off()  # Close the PDF file


sample_names <- sample_vector[sample_vector != "control"]
# Read all data frames into a named list
deg_mrna <- lapply(sample_names, function(sample) {
  df_name <- paste0("res_mrna_", sample, "_vs_control_filtered")
  get(df_name)  # Retrieve the data frame
})
names(deg_mrna) <- sample_names

# Find common genes across all data frames
common_genes <- Reduce(intersect, lapply(deg_mrna, rownames))

# Create a new data frame with log2FoldChange values
logfc_mrna <- do.call(cbind, lapply(deg_mrna, function(df) {
  df[common_genes, "log2FoldChange", drop = FALSE]
}))

# Rename columns with sample names
colnames(logfc_mrna) <- sample_names

# Set row names as gene IDs
rownames(logfc_mrna) <- common_genes

# View final dataframe
head(logfc_mrna)

logfc_mrna$average <- rowMeans(logfc_mrna, na.rm = TRUE)

write.csv(logfc_mrna, "deseq_results/mrna/logfc_mrna.csv")
# Read all data frames into a named list
deg_lncrna <- lapply(sample_names, function(sample) {
  df_name <- paste0("res_lncrna_", sample, "_vs_control_filtered")
  get(df_name)  # Retrieve the data frame
})
names(deg_lncrna) <- sample_names

# Find common genes across all data frames
common_lncrna <- Reduce(intersect, lapply(deg_lncrna, rownames))

# Create a new data frame with log2FoldChange values
logfc_lncrna <- do.call(cbind, lapply(deg_lncrna, function(df) {
  df[common_lncrna, "log2FoldChange", drop = FALSE]
}))

# Rename columns with sample names
colnames(logfc_lncrna) <- sample_names

# Set row names as gene IDs
rownames(logfc_lncrna) <- common_lncrna

# View final dataframe
head(logfc_lncrna)



# Connect to the Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get unique Ensembl gene IDs from logfc_mrna
ensembl_ids <- rownames(logfc_mrna)

# Retrieve HGNC symbols for Ensembl gene IDs
gene_annotation <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  filters = "ensembl_gene_id", 
  values = ensembl_ids, 
  mart = mart
)

# Ensure there are no empty HGNC symbols
gene_annotation <- gene_annotation[gene_annotation$hgnc_symbol != "", ]

# Map Ensembl IDs to HGNC symbols in logfc_mrna
logfc_mrna$hgnc_symbol <- gene_annotation$hgnc_symbol[match(rownames(logfc_mrna), gene_annotation$ensembl_gene_id)]

# Remove rows where HGNC symbol is NA
logfc_mrna <- logfc_mrna[!is.na(logfc_mrna$hgnc_symbol), ]

# Set row names to HGNC symbols
rownames(logfc_mrna) <- logfc_mrna$hgnc_symbol
logfc_mrna$hgnc_symbol <- NULL  # Remove the temporary column

# View the updated data frame
head(logfc_mrna)

# Ensure Ensembl IDs in `lncrna_genes` match those in `logfc_lncrna`
lncrna_genes$ensembl_gene_id <- sub("\\..*", "", lncrna_genes$ensembl_gene_id)  # Remove version numbers
rownames(logfc_lncrna) <- sub("\\..*", "", rownames(logfc_lncrna))  # Remove version numbers from row names

# Match Ensembl IDs with Gene Names
logfc_lncrna$gene_name <- lncrna_genes$gene_name[match(rownames(logfc_lncrna), lncrna_genes$ensembl_gene_id)]

# Remove rows where gene_name is NA (optional)
logfc_lncrna <- logfc_lncrna[!is.na(logfc_lncrna$gene_name), ]

# Set row names to gene names
rownames(logfc_lncrna) <- logfc_lncrna$gene_name
logfc_lncrna$gene_name <- NULL  # Remove the temporary column

# View updated data
head(logfc_lncrna)

logfc_mrna <- logfc_mrna[
  apply(logfc_mrna, 1, function(x) all(x >= 0) | all(x <= 0)), 
]


logfc_lncrna <- logfc_lncrna[
  apply(logfc_lncrna, 1, function(x) all(x >= 0) | all(x <= 0)), 
]


npinter_mirna_lnc_targets <- read.delim("anelorda/human_mirna_targets.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(npinter_mirna_lnc_targets) <- c("lncRNA_ID", "lncRNA_Name", "Some_ID", "RNA_Type", "Target_miRNA", "miRNA_ID", 
                                         "Target_Type", "Description", "Experimental_Method", "PMID", "Species", "Cell_Line", 
                                         "Interaction_Type", "Binding_Type", "Interaction_Category", "Source")


common_lncrna <- rownames(logfc_lncrna)
common_mrna <- rownames(logfc_mrna)

# Filter npinter_mirna_lnc_targets to retain only common lncRNAs
filtered_npinter <- npinter_mirna_lnc_targets[
  npinter_mirna_lnc_targets$lncRNA_Name %in% common_lncrna &
    npinter_mirna_lnc_targets$Target_miRNA %in% filtered_unique_pairs$miRNA, ]

# View results
head(filtered_npinter)


encori <- read.table("lncRNA_degs/clean_.txt", header = TRUE, sep = "\t")
# Filter npinter dataset based on common miRNAs and lncRNAs

filtered_npinter <- npinter_mirna_lnc_targets[
  npinter_mirna_lnc_targets$lncRNA_Name %in% common_lncrna & 
    npinter_mirna_lnc_targets$Target_miRNA %in% common_mirna, 
]

# Check the first few rows of the filtered dataset
head(filtered_npinter)

logfc_lncrna$Average <- rowMeans(logfc_lncrna, na.rm = TRUE)

# View the first few rows
head(logfc_lncrna)

filtered_unique_pairs_lncRNA <- filtered_unique_pairs_lncRNA[
  filtered_unique_pairs_lncRNA$lncRNA_Name %in% rownames(logfc_lncrna), ]

filtered_unique_pairs_lncRNA$logFC_lncRNA <- logfc_lncrna[
  filtered_unique_pairs_lncRNA$lncRNA_Name, "Average"]

# View updated data
head(filtered_unique_pairs_lncRNA)


filtered_df <- filtered_unique_pairs_lncRNA[
  filtered_unique_pairs_lncRNA$logFC_miRNA > 0 &
    filtered_unique_pairs_lncRNA$logFC_mRNA < 0 &
    filtered_unique_pairs_lncRNA$logFC_lncRNA < 0, ]

# View the first few rows
head(filtered_df)

# Create mRNA target DataFrame
mRNA_targets <- filtered_df[, c("miRNA", "Target.Gene")]
colnames(mRNA_targets) <- c("mirna", "target")

# Create lncRNA target DataFrame
lncRNA_targets <- filtered_df[, c("miRNA", "lncRNA_Name")]
colnames(lncRNA_targets) <- c("mirna", "target")

# Combine both
combined_targets <- rbind(mRNA_targets, lncRNA_targets)

# Remove any potential duplicates (optional)
combined_targets <- unique(combined_targets)

# View the result
head(combined_targets)


# Helper function to count mostly up/down genes
count_regulated_genes <- function(df, threshold_fraction = 0.5) {
  df_clean <- df[, !grepl("average|Average|miRNA", colnames(df), ignore.case = TRUE)]
  threshold <- ncol(df_clean) * threshold_fraction
  
  up_counts <- apply(df_clean, 1, function(x) sum(x > 0))
  down_counts <- apply(df_clean, 1, function(x) sum(x < 0))
  
  num_up <- sum(up_counts > threshold)
  num_down <- sum(down_counts > threshold)
  
  return(data.frame(Regulation = c("Upregulated", "Downregulated"),
                    Count = c(num_up, num_down)))
}

# Calculate for all datasets
lncRNA_summary <- count_regulated_genes(logfc_lncrna)
lncRNA_summary$Type <- "lncRNA"

miRNA_summary <- count_regulated_genes(logfc_mirna)
miRNA_summary$Type <- "miRNA"

mRNA_summary <- count_regulated_genes(logfc_mrna)
mRNA_summary$Type <- "mRNA"

# Combine all
df_all <- rbind(mRNA_summary, lncRNA_summary, miRNA_summary)

# Plot
library(ggplot2)

p <- ggplot(df_all, aes(x = Type, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  theme_minimal() +
  labs(title = "Number of Upregulated and Downregulated lncRNAs, miRNAs\nand mRNAs Across All Samples",
       x = "Gene Type",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Upregulated" = "coral3", "Downregulated" = "deepskyblue4"))

ggsave(
  filename = "gene_reg_barplot.jpeg",
  plot = p + theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  ),
  device = "jpeg",
  width = 8,
  height = 6,
  dpi = 300
)

