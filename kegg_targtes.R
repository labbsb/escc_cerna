library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidygraph)
library(ggraph)
library(ggplot2)

mirna_mrna_targets <- read.csv("MIRNA_3_TARGETS/DEM3_to_DEG.csv", header = TRUE)
head(mirna_mrna_targets)

mirna_mrna_targets_filtered <- mirna_mrna_targets %>%
  filter(expression_x == "up" & expression_y == "down")




mirna_mrna_1 <- read.csv("MIRNA_1_TARGETS/dem1_target_deg_unique.csv", header = TRUE)

head(mirna_mrna_1)

# miRNA up, Gene down
mirna_up_gene_down <- subset(mirna_mrna_1, expression_x == "up" & expression_y == "down")

# miRNA down, Gene up
mirna_down_gene_up <- subset(mirna_mrna_1, expression_x == "down" & expression_y == "up")
write.csv(mirna_down_gene_up, "MIRNA_1_TARGETS/mirna_down_gene_up.csv")



# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

genes_down <- unique(mirna_up_gene_down$Gene)
genes_up   <- unique(mirna_down_gene_up$Gene)
### --- Function to run enrichment ---
run_enrichment <- function(gene_symbols, group_label) {
  message("Running enrichment for ", group_label, " ...")
  
  # SYMBOL → ENTREZ ID
  gene_df <- bitr(gene_symbols,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  
  # KEGG
  kegg_res <- enrichKEGG(gene = gene_df$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 1)
  
  # GO BP
  go_bp <- enrichGO(gene          = gene_df$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1)
  
  # Convert GO gene IDs → Symbols
  go_res <- as.data.frame(go_bp@result)
  go_gene_list <- strsplit(go_res$geneID, "/")
  
  id2symbol_go <- bitr(unique(unlist(go_gene_list)),
                       fromType = "ENTREZID",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)
  
  go_res$geneSymbol <- sapply(go_gene_list, function(ids) {
    syms <- id2symbol_go$SYMBOL[match(ids, id2symbol_go$ENTREZID)]
    paste(syms, collapse = "/")
  })
  
  # Return results in a list
  list(kegg = kegg_res, go_bp = go_bp,
       kegg_df = as.data.frame(kegg_res),
       go_df = go_res)
}

### --- Run for gene_up and gene_down ---
res_up <- run_enrichment(genes_up, "Gene Up (miRNA down)")
res_down <- run_enrichment(genes_down, "Gene Down (miRNA up)")



### --- View top results ---
head(res_up$kegg_df)
head(res_up$go_df[, c("ID", "Description", "geneSymbol")])

head(res_down$kegg_df)
head(res_down$go_df[, c("ID", "Description", "geneSymbol")])

# Extract KEGG results as a dataframe
kegg_res_up <- as.data.frame(res_up$kegg@result)

# Split Entrez IDs
gene_lists <- strsplit(kegg_res_up$geneID, "/")

# Map unique Entrez IDs to SYMBOLs using org.Hs.eg.db
id2symbol <- bitr(unique(unlist(gene_lists)),
                  fromType = "ENTREZID",
                  toType   = "SYMBOL",
                  OrgDb    = org.Hs.eg.db)

# Add SYMBOLs to KEGG results
kegg_res_up$GeneSymbol <- sapply(gene_lists, function(ids) {
  syms <- id2symbol$SYMBOL[match(ids, id2symbol$ENTREZID)]
  paste(syms, collapse = "/")
})

# Check the updated table
head(kegg_res_up[, c("ID", "Description", "GeneSymbol")])

top4_ids <- kegg_res_df %>%
  arrange(p.adjust) %>%
  head(4) %>%
  pull(ID)

kegg_up_top4 <- kegg_rebuilt[kegg_rebuilt@result$ID %in% top4_ids,]

kegg_genes <- unique(unlist(strsplit(kegg_up_top4$GeneSymbol, "/")))
mirna_gene_pathway <- mirna_down_gene_up[mirna_down_gene_up$Gene %in% kegg_genes, ]

head(mirna_gene_pathway)

library(enrichplot)
library(ggplot2)

# Create a list where each pathway ID links to its genes
pathway_gene_list <- setNames(
  strsplit(kegg_up_top4$GeneSymbol, "/"),
  kegg_up_top4$Description
)

# make gene-miRNA mapping vector
gene2mirna <- setNames(mirna_gene_pathway$miRNA, mirna_gene_pathway$Gene)

library(igraph)
library(ggraph)
gene2mirna <- data.frame(
  Gene = names(gene2mirna),
  miRNA = unname(gene2mirna),
  stringsAsFactors = FALSE
)

covid_genes <- pathway_gene_list$`Coronavirus disease - COVID-19`
ribo_genes  <- pathway_gene_list$Ribosome
prot_exp    <- pathway_gene_list$`Protein export`
necro_genes <- pathway_gene_list$Necroptosis

# Determine categories
gene_categories <- tibble(Gene = unique(gene2mirna$Gene)) %>%
  mutate(Category = case_when(
    Gene %in% covid_genes & Gene %in% ribo_genes ~ "Both (Coronavirus disease - COVID-19 & Ribosome)",
    Gene %in% covid_genes & !(Gene %in% ribo_genes) ~ "Coronavirus disease - COVID-19",
    Gene %in% prot_exp ~ "Protein export",
    Gene %in% necro_genes ~ "Necroptosis",
    TRUE ~ "Other"
  ))



edges <- gene2mirna %>%
  rename(from = miRNA, to = Gene)

# Make sure nodes include both miRNAs and mRNAs
nodes <- tibble(
  name = unique(c(edges$from, edges$to))
) %>%
  left_join(gene_categories, by = c("name" = "Gene")) %>%
  mutate(
    type = if_else(name %in% edges$to, "mRNA", "miRNA")
  )

# Create graph (undirected so layout is clean)
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)


p <- ggraph(g, layout = "fr") +
  geom_edge_link(color = "grey70", alpha = 0.6) +
  geom_node_point(aes(color = Category, shape = type), size = 6) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 5,              # 🔹 increase text size here
    fontface = "bold"      # optional: makes labels clearer
  ) +
  scale_shape_manual(values = c(miRNA = 17, mRNA = 16)) +
  scale_color_manual(
    values = c(
      "Both (Coronavirus disease - COVID-19 & Ribosome)" = "#E41A1C",
      "Coronavirus disease - COVID-19" = "#377EB8",
      "Protein export" = "#4DAF4A",
      "Necroptosis" = "#984EA3",
      "Other" = "grey80"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "grey95"),
    legend.position = "right"
  )

# 🔹 Save as a 300-dpi PDF
ggsave(
  filename = "gene_miRNA_network.pdf",
  plot = p,
  width = 10,        # in inches
  height = 8,
  dpi = 300,
  device = cairo_pdf  # ensures vector quality text and shapes
)

setwd("MIRNA_3_TARGETS/MIRNA_3_TARGETS/")
# Get all files ending with _trimmed.txt
trimmed_files <- list.files(pattern = "_trimmed\\.txt$")

# Read all files into a named list of data frames
trimmed_data <- lapply(trimmed_files, function(file) {
  read.table(file, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
})

# Assign file names (without extension) as names to the list elements
names(trimmed_data) <- sub("_trimmed\\.txt$", "", trimmed_files)

head(trimmed_data$`hsa-miR-106b-5p`)




