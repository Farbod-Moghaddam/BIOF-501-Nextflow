#!/usr/bin/env Rscript

# Differential Expression Analysis
# Author: RNA-seq Pipeline
# Description: Performs DESeq2 differential expression analysis

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript differential_expression.R <count_matrix> <metadata_csv> <output_dir>")
}

count_matrix_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Starting Differential Expression Analysis\n")
cat("==========================================\n")
cat("Count matrix:", count_matrix_file, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# ============================================================================
# 1. Load and prepare data
# ============================================================================

cat("Loading count matrix...\n")
count_data <- read.table(count_matrix_file, header = TRUE, row.names = 1, sep = "\t")

cat("Loading metadata...\n")
metadata <- read.csv(metadata_file, header = TRUE)
rownames(metadata) <- metadata$Sample

# Ensure count matrix columns match metadata rows
count_data <- count_data[, metadata$Sample]

# Convert Flare to factor with meaningful labels
metadata$Flare <- factor(metadata$Flare, levels = c(0, 1), labels = c("NoFlare", "Flare"))

cat("\nDataset summary:\n")
cat("  Samples:", ncol(count_data), "\n")
cat("  Genes:", nrow(count_data), "\n")
cat("  Flare samples:", sum(metadata$Flare == "Flare"), "\n")
cat("  NoFlare samples:", sum(metadata$Flare == "NoFlare"), "\n\n")

# ============================================================================
# 2. DESeq2 Differential Expression Analysis
# ============================================================================

cat("Running DESeq2 analysis...\n")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData = metadata,
  design = ~ Flare
)

# Filter low count genes (keep genes with at least 10 reads in at least 3 samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

cat("  Retained", sum(keep), "genes after filtering\n")

# Run DESeq2
dds <- DESeq(dds)

# Get results (Flare vs NoFlare)
res <- results(dds, contrast = c("Flare", "Flare", "NoFlare"))
res <- res[order(res$padj), ]

cat("  Differential expression complete\n")

# Summary of DE genes
cat("\nDifferential Expression Summary:\n")
cat("  Upregulated (padj < 0.05, log2FC > 1):", sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("  Downregulated (padj < 0.05, log2FC < -1):", sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE), "\n")
cat("  Total significant (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n\n")

# Save full results
write.csv(as.data.frame(res), file.path(output_dir, "differential_expression_results.csv"))

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file.path(output_dir, "normalized_counts.csv"))

# ============================================================================
# 3. Volcano Plot
# ============================================================================

cat("Creating volcano plot...\n")

# Prepare data for volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]

# Classify significance
res_df$significant <- "Not Significant"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Create simple volcano plot with ggplot2
png(file.path(output_dir, "volcano_plot.png"), width = 3000, height = 2400, res = 300)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  labs(title = "Differential Expression: Flare vs NoFlare",
       subtitle = paste0("Total genes: ", nrow(res)),
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)",
       color = "Significance") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "right")
dev.off()

cat("  Volcano plot saved\n\n")

# ============================================================================
# 4. Summary of Results
# ============================================================================

# Get significant genes
sig_up_genes <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 1)]
sig_down_genes <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange < -1)]

# Save significant gene lists
if (length(sig_up_genes) > 0) {
  write.csv(data.frame(gene = sig_up_genes), file.path(output_dir, "upregulated_genes.csv"), row.names = FALSE)
}
if (length(sig_down_genes) > 0) {
  write.csv(data.frame(gene = sig_down_genes), file.path(output_dir, "downregulated_genes.csv"), row.names = FALSE)
}

cat("\n==========================================\n")
cat("Analysis Complete!\n")
cat("Results saved to:", output_dir, "\n")
cat("Note: Gene Ontology enrichment skipped (requires additional packages)\n")
cat("==========================================\n")
