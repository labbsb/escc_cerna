library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(enrichplot)


up_genes <- rownames(logfc_mrna[logfc_mrna$average > 0, ])
down_genes <- rownames(logfc_mrna[logfc_mrna$average < 0, ])

up_lncrna <- rownames(logfc_lncrna_ensemblid[logfc_lncrna_ensemblid$average > 0, ])
down_lncrna <- rownames(logfc_lncrna_ensemblid[logfc_lncrna_ensemblid$average < 0, ])

up_genes <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_genes <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

up_lncrna <- bitr(up_lncrna, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_lncrna <- bitr(down_lncrna, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

up_all <- unique(c(up_genes$ENTREZID, up_lncrna$ENTREZID))
down_all <- unique(c(down_genes$ENTREZID, down_lncrna$ENTREZID))


# ==========================
# KEGG enrichment
# ==========================
ekegg_up <- enrichKEGG(gene         = up_all,
                       organism     = 'hsa',
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05)

ekegg_down <- enrichKEGG(gene       = down_all,
                         organism   = 'hsa',
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05)

# ==========================
# GO enrichment (BP, CC, MF)
# ==========================
ego_bp_up <- enrichGO(gene          = up_all,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05)

ego_bp_down <- enrichGO(gene        = down_all,
                        OrgDb       = org.Hs.eg.db,
                        keyType     = "ENTREZID",
                        ont         = "BP",
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.05)

#ego_cc_up <- enrichGO(gene          = up_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05)
#ego_cc_down <- enrichGO(gene        = down_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05)

#ego_mf_up <- enrichGO(gene          = up_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05)
#ego_mf_down <- enrichGO(gene        = down_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05)

# ==========================
# Reactome enrichment
# ==========================
ereact_up <- enrichPathway(gene          = up_all,
                           organism      = "human",
                           pAdjustMethod = "BH",
                           qvalueCutoff  = 0.05)

ereact_down <- enrichPathway(gene        = down_all,
                             organism    = "human",
                             pAdjustMethod = "BH",
                             qvalueCutoff  = 0.05)

# ==========================
# Hallmark enrichment (MSigDB)
# ==========================
# Get hallmark gene sets
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")

# Run GSEA-style analysis (if you have ranked list)
# Or ORA-style using enricher()
ehall_up <- enricher(gene          = up_all,
                     TERM2GENE     = msig_hallmark[, c("gs_name", "ncbi_gene")],
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)

ehall_down <- enricher(gene        = down_all,
                       TERM2GENE   = msig_hallmark[, c("gs_name", "ncbi_gene")],
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)


library(dplyr)
library(ggplot2)

# Function to process and plot enrichment results
plot_enrichment <- function(go_bp_res, kegg_res, hallmark_res, direction, outdir = "plots") {
  
  # Add category column
  go_bp <- go_bp_res %>% mutate(Category = "GO BP")
  kegg <- kegg_res %>% mutate(Category = "KEGG")
  hallmark <- hallmark_res %>% mutate(Category = "Hallmark")
  
  # Combine all results
  all_results <- bind_rows(go_bp, kegg, hallmark)
  
  # Select top 10 enriched terms per category (based on adjusted p-value)
  all_results <- all_results %>%
    group_by(Category) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Factor Category to enforce order
  all_results$Category <- factor(all_results$Category,
                                 levels = c("GO BP", "KEGG", "Hallmark"))
  
  # Sort terms inside each category
  all_results <- all_results %>%
    group_by(Category) %>%
    arrange(p.adjust) %>%
    mutate(Description = factor(Description, levels = rev(Description))) %>%
    ungroup()
  
  # Plot
  p <- ggplot(all_results, aes(x = Description, y = -log10(p.adjust), fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(x = "Pathway / GO Term",
         y = "-log10(p.adjust)",
         title = paste0("Enrichment Analysis (", direction, "-regulated genes)")) +
    scale_fill_manual(values = c("GO BP" = "lightblue", 
                                 "KEGG" = "red", 
                                 "Hallmark" = "purple")) +
    theme(legend.title = element_blank())
  
  # Save plot
  if (!dir.exists(outdir)) dir.create(outdir)
  ggsave(file.path(outdir, paste0("plots/GO_KEGG_Hallmark_enrichment_", direction, ".pdf")),
         plot = p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

# Now just call the function for up and down
plots <- list(
  down = plot_enrichment(ego_bp_down@result, ekegg_down@result, ehall_down@result, "down"),
  up   = plot_enrichment(ego_bp_up@result,   ekegg_up@result,   ehall_up@result,   "up")
)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(cowplot)

### --- DOWNREGULATED ---
ekegg_down <- setReadable(ekegg_down, 'org.Hs.eg.db', keyType = "ENTREZID")
ego_bp_down <- setReadable(ego_bp_down, 'org.Hs.eg.db', keyType = "ENTREZID")

p_kegg_down <- tryCatch(cnetplot(ekegg_down, showCategory = 5,
                                 circular = TRUE, colorEdge = TRUE,
                                 max.overlaps = 100,
                                 cex_label_category = 0.8, cex_label_gene = 0.6,
                                 label_format = 30), error=function(e) NULL)

p_gobp_down <- tryCatch(cnetplot(ego_bp_down, showCategory = 5,
                                 circular = TRUE, colorEdge = TRUE,
                                 max.overlaps = 100,
                                 cex_label_category = 0.8, cex_label_gene = 0.6,
                                 label_format = 30), error=function(e) NULL)

down_grid <- cowplot::plot_grid(
  p_kegg_down, p_gobp_down,
  ncol = 2,
  labels = c("KEGG (DOWN)", "GO BP (DOWN)"),
  label_size = 14,
  label_x = 0.1,   # shift labels to the left
  label_y = 1.05   # shift labels above plots
)

### --- UPREGULATED ---
ekegg_up <- setReadable(ekegg_up, 'org.Hs.eg.db', keyType = "ENTREZID")
ego_bp_up <- setReadable(ego_bp_up, 'org.Hs.eg.db', keyType = "ENTREZID")
ehall_up  <- setReadable(ehall_up,  'org.Hs.eg.db', keyType = "ENTREZID")

p_kegg_up <- tryCatch(cnetplot(ekegg_up, showCategory = 5,
                               circular = TRUE, colorEdge = TRUE,
                               max.overlaps = 100,
                               cex_label_category = 0.8, cex_label_gene = 0.6,
                               label_format = 30), error=function(e) NULL)

p_gobp_up <- tryCatch(cnetplot(ego_bp_up, showCategory = 5,
                               circular = TRUE, colorEdge = TRUE,
                               max.overlaps = 100,
                               cex_label_category = 0.8, cex_label_gene = 0.6,
                               label_format = 30), error=function(e) NULL)

p_hall_up <- tryCatch(cnetplot(ehall_up, showCategory = 5,
                               circular = TRUE, colorEdge = TRUE,
                               max.overlaps = 100,
                               cex_label_category = 0.8, cex_label_gene = 0.6,
                               label_format = 30), error=function(e) NULL)

up_grid <- cowplot::plot_grid(
  p_kegg_up, p_gobp_up, p_hall_up,
  ncol = 3,
  labels = c("KEGG (UP)", "GO BP (UP)", "Hallmark (UP)"),
  label_size = 14,
  label_x = 0.1,
  label_y = 1.05
)

### --- FINAL COMBINED LAYOUT ---
final_grid <- cowplot::plot_grid(
  down_grid, up_grid,
  ncol = 1,
  rel_heights = c(1, 1.2)  # give more space to UP plots (3 columns)
) + theme(plot.margin = margin(30, 30, 30, 30))  # extra bottom margin

# Save to PDF (larger canvas)
ggsave("plots/cnetplot_up_down.pdf",
       plot   = final_grid,
       width  = 22,
       height = 16,
       dpi    = 300)
