## ============================================================================
## Permutation test for miRNA DEM robustness (reviewer-requested)
##
## Goal: confirm the number of DEMs found in the real pooled tumor-vs-control
## comparison (4 controls vs 16 tumor patients) is not an artifact of WHICH
## 4 of the 20 biological units happen to be labeled "control". This
## complements the existing per-patient LOO analysis: LOO asks "does the
## result depend on which controls I picked, holding the tumor/control
## structure fixed?"; this permutation test asks "do the tumor/control labels
## carry real information at all, or would any random 4-vs-16 split of the
## same 20 biological units give a similarly large number of DEMs?"
##
## Design choices (mirrors the real pipeline -- document these in the
## rebuttal letter / Methods):
##   - The unit of permutation is the BIOLOGICAL SAMPLE (control or tumor
##     patient), not the individual replicate column. All replicate columns
##     belonging to one biological unit move together each iteration, exactly
##     as collapseReplicates() treats them as one unit in the real analysis.
##   - Group sizes are held fixed at 4 "control" vs 16 "tumor" each iteration
##     (only WHICH units get which label is shuffled).
##   - Each iteration: collapse replicates by unit -> pooled ~ condition
##     DESeq2 model (identical to your dds_collapsed step) -> count DEMs at
##     padj < 0.05 & |log2FC| > 1 (same thresholds as your real analysis).
##   - 1,000 iterations as requested by the reviewer.
##   - Empirical p-value = (# permutations with DEM count >= observed + 1)
##     / (n_perm + 1)   [one-sided; testing whether real labels give an
##     unusually LARGE count vs. the null]
## ============================================================================

library(DESeq2)
library(dplyr)
library(ggplot2)

set.seed(42)  # reproducibility -- report this seed in Methods/rebuttal

## ---- 1. Load data exactly as in the real pipeline -----------------------

data <- read.table(
  "ksarkytbayev/nextflow_24_04_25/mirna_quant/edger_qc/mature_counts.csv",
  header = TRUE, sep = ","
)
rownames(data) <- data$X
data$X <- NULL
data <- t(data)
data <- data + 1

control_cols <- sort(grep("^control_", colnames(data), value = TRUE))
sample_cols  <- sort(setdiff(colnames(data), control_cols))
data <- data[, c(control_cols, sample_cols)]

sample_vector <- c(
  "control1", "control2", "control3", "control4",
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
row_names <- colnames(data)

## biological unit per column (this is what gets shuffled as a block)
unit_vector <- sample_vector  # "control" or "X001T" etc. -- one label per column
units       <- unique(unit_vector)          # 20 biological units
n_units     <- length(units)

control_units <- unique(unit_vector[grepl("^control", unit_vector)])
tumor_units   <- setdiff(units, control_units)
n_control_units <- length(control_units)   # expect 4
n_tumor_units   <- length(tumor_units)     # expect 16

cat("Biological units:", n_units, "( control:", n_control_units, ", tumor:", n_tumor_units, ")\n")
stopifnot(n_control_units == 4, n_tumor_units == 16)

## ---- 2. Helper: collapse replicates by unit + condition label, run DESeq2,
##         return number of significant DEMs -------------------------------

run_pooled_DEM_count <- function(counts_mat, unit_vec, condition_for_unit,
                                 padj_cutoff = 0.05, lfc_cutoff = 1) {
  # condition_for_unit: named vector, names = unit labels, values = "control"/"tumor"
  coldata <- data.frame(
    unit      = unit_vec,
    condition = condition_for_unit[unit_vec],
    row.names = colnames(counts_mat)
  )
  dds_raw <- DESeqDataSetFromMatrix(
    countData = counts_mat, colData = coldata, design = ~ condition
  )
  # collapse all replicate columns sharing the same biological unit into one
  dds_c <- collapseReplicates(dds_raw, groupby = factor(unit_vec), run = colnames(dds_raw))
  # collapseReplicates() renames columns to the groupby factor's levels, i.e.
  # to the unit labels themselves -- so colnames(dds_c) are unit labels already.
  colData(dds_c)$condition <- factor(condition_for_unit[colnames(dds_c)],
                                     levels = c("control", "tumor"))
  design(dds_c) <- ~ condition
  dds_c <- DESeq(dds_c, quiet = TRUE)
  res <- results(dds_c, contrast = c("condition", "tumor", "control"))
  sum(!is.na(res$padj) & res$padj < padj_cutoff & abs(res$log2FoldChange) > lfc_cutoff)
}

## ---- 3. Observed (real) result -------------------------------------------

real_condition_for_unit <- setNames(
  ifelse(units %in% control_units, "control", "tumor"),
  units
)

observed_DEMs <- run_pooled_DEM_count(data, unit_vector, real_condition_for_unit)
cat("Observed DEMs with true labels:", observed_DEMs, "\n")
# Sanity check: reconcile this against the DEM count already reported in the
# manuscript's pooled tumor-vs-control analysis before using it in the rebuttal.

## ---- 4. Permutation null distribution -------------------------------------

n_perm <- 1000
perm_DEMs <- integer(n_perm)

for (i in seq_len(n_perm)) {
  shuffled_control_units <- sample(units, n_control_units)
  shuffled_condition_for_unit <- setNames(
    ifelse(units %in% shuffled_control_units, "control", "tumor"),
    units
  )
  
  perm_DEMs[i] <- tryCatch(
    withCallingHandlers(
      run_pooled_DEM_count(data, unit_vector, shuffled_condition_for_unit),
      warning = function(w) {
        cat("  [perm", i, "warning]:", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      cat("  [perm", i, "FAILED]:", conditionMessage(e), "\n")
      NA_integer_
    }
  )
  if (i %% 50 == 0) cat("Permutation", i, "/", n_perm, "- DEMs:", perm_DEMs[i], "\n")
}

n_failed <- sum(is.na(perm_DEMs))
if (n_failed > 0) cat(n_failed, "permutations failed (e.g. DESeq2 convergence) and were excluded.\n")
perm_DEMs_clean <- perm_DEMs[!is.na(perm_DEMs)]

## ---- 5. Empirical p-value and summary -------------------------------------

empirical_p <- (sum(perm_DEMs_clean >= observed_DEMs) + 1) / (length(perm_DEMs_clean) + 1)

cat("\n--- Permutation test summary ---\n")
cat("Observed DEMs (true labels):   ", observed_DEMs, "\n")
cat("Permutation null mean DEMs:    ", round(mean(perm_DEMs_clean), 1), "\n")
cat("Permutation null SD:           ", round(sd(perm_DEMs_clean), 1), "\n")
cat("Permutation null max DEMs:     ", max(perm_DEMs_clean), "\n")
cat("Empirical p-value:             ", signif(empirical_p, 3), "\n")

## ---- 6. Save results -------------------------------------------------------

dir.create("permutation_test_mirna", showWarnings = FALSE)
write.csv(
  data.frame(permutation = seq_len(n_perm), n_DEMs = perm_DEMs),
  "permutation_test_mirna/permutation_test_DEM_counts.csv",
  row.names = FALSE
)

## ---- 7. Plot ---------------------------------------------------------------

p <- ggplot(data.frame(n_DEMs = perm_DEMs_clean), aes(x = n_DEMs)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(xintercept = observed_DEMs, color = "firebrick", linewidth = 1, linetype = "dashed") +
  annotate("text", x = observed_DEMs, y = Inf,
           label = paste0("Observed = ", observed_DEMs),
           color = "firebrick", vjust = 1.5, hjust = -0.05) +
  labs(x = "Number of significant DEMs (padj < 0.05, |log2FC| > 1)",
       y = "Count of permutations",
       title = "Permutation null distribution vs. observed DEM count",
       subtitle = paste0("n = ", n_perm, " permutations (4 control / 16 tumor units shuffled); empirical p = ",
                         signif(empirical_p, 3))) +
  theme_minimal(base_size = 12)

ggsave("permutation_test_mirna/permutation_test_DEM_distribution.png", p, width = 7, height = 5, dpi = 300)
ggsave("permutation_test_mirna/permutation_test_DEM_distribution.jpeg", p,
       width = 7, height = 5, dpi = 300, quality = 95)

cat("\nDone. Outputs written to permutation_test_mirna/:\n")
cat(" - permutation_test_DEM_counts.csv\n")
cat(" - permutation_test_DEM_distribution.png\n")
