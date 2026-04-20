#Author Anel Ordabayeva

setwd("/mnt/sdb2/anelorda/MIRNA_1_TARGETS/MIRNA_1_TARGETS")

library(dplyr)

# Get all files ending with _trimmed.txt
trimmed_files <- list.files(pattern = "_trimmed\\.txt$")

# Read all files into a named list of data frames
trimmed_data <- lapply(trimmed_files, function(file) {
  read.table(file, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
})

# Assign file names (without extension) as names to the list elements
names(trimmed_data) <- sub("_trimmed\\.txt$", "", trimmed_files)

head(trimmed_data$`hsa-miR-106b-5p`)


# Extract average logFC values, with geneID as names
logfc_lncrna$average <- rowMeans(logfc_lncrna[, grep("^X", names(logfc_lncrna))], na.rm = TRUE)

lncrna_avg <- logfc_lncrna$average
names(lncrna_avg) <- logfc_lncrna$ensembl_gene_id # these are geneIDs

# Function: sort and annotate by geneID only
sort_and_annotate_by_geneID <- function(df) {
  if (!("geneID" %in% colnames(df))) return(df)  # Skip if geneID is missing
  
  # Match geneID to logFC data
  match_idx <- match(df$geneID, names(lncrna_avg))
  
  # Assign "expression" based on average logFC
  avg_values <- lncrna_avg[match_idx]
  df$expression <- ifelse(is.na(avg_values), NA,
                          ifelse(avg_values > 0, "up", "down"))
  
  # Assign sort index and sort
  df$sort_index <- match(df$geneID, logfc_lncrna$ensembl_gene_id)
  df <- df[order(df$sort_index, na.last = TRUE), ]
  df$sort_index <- NULL
  
  return(df)
}

# Apply to all trimmed data
trimmed_data_sorted <- lapply(trimmed_data, sort_and_annotate_by_geneID)

# Example: view result
head(trimmed_data_sorted$`hsa-miR-106b-5p`)



# Print head of all sorted data frames
for (mirna in names(trimmed_data_sorted)) {
  cat("\n======", mirna, "======\n")
  print(head(trimmed_data_sorted[[mirna]]))
}


combined_df <- do.call(rbind, lapply(trimmed_data_sorted, function(df) {
  df[, c("miRNAname", "geneName", "expression", "geneID")]
}))

# Remove duplicate rows
combined_unique <- unique(combined_df)

combined_unique <- combined_unique[!is.na(combined_unique$expression), ]

# View the result
head(combined_unique)

write.csv(combined_unique, "../mirna_target_lncrna_all.csv", row.names = FALSE)

# Make sure column names are consistent
colnames(combined_unique) <- c("miRNA", "lncRNA", "expression_lncRNA", "lncRNA_ID")
colnames(mirna_mrna_1) <- c("X", "miRNA", "expression_miRNA", "mRNA", "expression_mRNA", "Oncogene.Suppressor")

# Merge by miRNA (left join to keep all miRNA–lncRNA pairs)
combined_final <- merge(
  combined_unique,
  mirna_mrna_1[, c("miRNA", "expression_miRNA", "mRNA", "expression_mRNA")],
  by = "miRNA",
  all.x = TRUE
)

# Reorder columns
combined_final <- combined_final[, c("miRNA", "expression_miRNA", "lncRNA", "expression_lncRNA", "mRNA", "expression_mRNA")]


set1 <- combined_final %>%
  filter(expression_lncRNA == "up" & 
           expression_miRNA == "down" & 
           expression_mRNA == "up")

# Set 2: lncRNA down, miRNA up, mRNA down
set2 <- combined_final %>%
  filter(expression_lncRNA == "down" & 
           expression_miRNA == "up" & 
           expression_mRNA == "down")

# Check counts
nrow(set1)
nrow(set2)
# Save to CSV
write.csv(combined_final, "combined_mirna_lncrna_mrna.csv", row.names = FALSE)

set2 <- set2 %>%
  mutate(lncRNA = paste0("lnc-", lncRNA))

# Combine lncRNA and mRNA as targets
targets_df <- set2 %>%
  select(miRNA, lncRNA, mRNA) %>%
  pivot_longer(cols = c(lncRNA, mRNA), 
               names_to = "target_type", 
               values_to = "target") %>%
  select(miRNA, target) %>%
  distinct()  # remove duplicates if any

# View result
head(targets_df)

targets_df_filtered <- targets_df %>%
  filter(miRNA %in% mirna_mrna_targets_filtered$miRNA)

set_mirna3 <- set2 %>%
  filter(miRNA %in% mirna_mrna_targets_filtered$miRNA)
