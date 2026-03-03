# limma_analysis.R
# Differential protein expression analysis using limma
# Conditions: Undiff vs RA-differentiated vs RA+PMA-differentiated (SH-SY5Y)

library(limma)
library(dplyr)
library(ggplot2)

# --- 1. Load Data ---
lfq_file <- "results/log2_lfq_matrix_for_limma.csv"

if (!file.exists(lfq_file)) {
  stop("Error: 'results/log2_lfq_matrix_for_limma.csv' not found. Run Python script first.")
}

lfq_data <- read.csv(lfq_file, row.names = 1, check.names = FALSE)
cat("Data dimensions (Proteins x Samples):", dim(lfq_data), "\n")
cat("Sample names:\n")
print(colnames(lfq_data))

# --- 2. Define Groups ---
samples <- colnames(lfq_data)

group_labels <- rep(NA, length(samples))
group_labels[grepl("undiff", samples, ignore.case = TRUE)] <- "Undiff"
group_labels[grepl("^RA",    samples)]                     <- "RA"
group_labels[grepl("PMA",    samples, ignore.case = TRUE)] <- "PMA"

# Check all samples are assigned
unassigned <- samples[is.na(group_labels)]
if (length(unassigned) > 0) {
  cat("WARNING: These samples were not assigned to a group:\n")
  print(unassigned)
}

group_labels <- factor(group_labels, levels = c("Undiff", "RA", "PMA"))

cat("\nSample distribution:\n")
print(table(group_labels))

# --- 3. Design Matrix ---
design <- model.matrix(~0 + group_labels)
colnames(design) <- levels(group_labels)
cat("\nDesign matrix columns:", colnames(design), "\n")

# --- 4. Linear Modeling ---
fit <- lmFit(lfq_data, design)

# Two contrasts: RA vs Undiff, PMA vs Undiff
contrast_matrix <- makeContrasts(
  RA_vs_Undiff  = RA - Undiff,
  PMA_vs_Undiff = PMA - Undiff,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# MA plot
png("results/MA_plot_limma.png", width = 800, height = 600)
plotMA(fit2, main = "MA Plot - limma eBayes")
dev.off()
cat("MA plot saved.\n")

# --- 5. Extract and Save Results ---
sig_p_thresh <- 0.05
logfc_thresh <- 1

extract_results <- function(fit2, coef_name, label) {
  results <- topTable(fit2, coef = coef_name, number = nrow(fit2),
                      adjust.method = "BH", sort.by = "P")

  cat(sprintf("\n--- %s ---\n", label))
  cat("Total proteins tested:", nrow(results), "\n")

  sig_proteins <- results %>%
    filter(adj.P.Val < sig_p_thresh & abs(logFC) > logfc_thresh)

  cat("Significant proteins (FDR<0.05 & |logFC|>1):", nrow(sig_proteins), "\n")

  write.csv(results,      sprintf("results/all_proteins_limma_%s.csv", label),     row.names = TRUE)
  write.csv(sig_proteins, sprintf("results/significant_proteins_limma_%s.csv", label), row.names = TRUE)

  return(results)
}

results_RA  <- extract_results(fit2, "RA_vs_Undiff",  "RA_vs_Undiff")
results_PMA <- extract_results(fit2, "PMA_vs_Undiff", "PMA_vs_Undiff")

# --- 6. Volcano Plot Function ---
plot_volcano <- function(results, label, title) {
  results$diffexpressed <- "Not significant"
  results$diffexpressed[results$logFC >  logfc_thresh & results$adj.P.Val < sig_p_thresh] <- "Upregulated"
  results$diffexpressed[results$logFC < -logfc_thresh & results$adj.P.Val < sig_p_thresh] <- "Downregulated"

  sig_count <- sum(results$diffexpressed != "Not significant")

  p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = diffexpressed)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(values = c(
      "Upregulated"     = "#F44336",
      "Downregulated"   = "#2196F3",
      "Not significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh),
               col = "black", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = -log10(sig_p_thresh),
               col = "black", linetype = "dashed", linewidth = 0.5) +
    labs(
      title    = title,
      subtitle = paste("Significant proteins:", sig_count),
      x        = "log2 Fold Change",
      y        = "-log10 Adjusted P-value",
      color    = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  ggsave(sprintf("results/volcano_limma_%s.png", label), p, width = 8, height = 6, dpi = 150)
  cat(sprintf("Volcano plot saved: results/volcano_limma_%s.png\n", label))
}

plot_volcano(results_RA,  "RA_vs_Undiff",  "Limma: RA-differentiated vs Undifferentiated")
plot_volcano(results_PMA, "PMA_vs_Undiff", "Limma: RA+PMA-differentiated vs Undifferentiated")

cat("\nLimma analysis complete. All results saved in 'results/'\n")