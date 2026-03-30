library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)
counts <- read.table("counts/counts_clean.txt",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t")

colnames(counts) <- c("Untreated", "Dexamethasone")

coldata <- data.frame(
  condition = factor(c("Untreated", "Dexamethasone"),
                     levels = c("Untreated", "Dexamethasone")),
  row.names = colnames(counts)
)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds,
               contrast = c("condition", "Dexamethasone", "Untreated"),
               alpha = 0.05)

res_shrunk <- lfcShrink(dds,
                        coef = "condition_Dexamethasone_vs_Untreated",
                        type = "apeglm")

summary(res_shrunk)

degs <- as.data.frame(res_shrunk) %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

cat("Significant DEGs:", nrow(degs), "\n")
cat("Upregulated in Dexamethasone:", sum(degs$log2FoldChange > 0), "\n")
cat("Downregulated in Dexamethasone:", sum(degs$log2FoldChange < 0), "\n")

write.csv(degs, "results/deseq2/significant_DEGs.csv", quote=FALSE)
write.csv(as.data.frame(res_shrunk), "results/deseq2/all_results.csv", quote=FALSE)

vsd <- vst(dds, blind=TRUE)
pca_plot <- plotPCA(vsd, intgroup="condition") +
  geom_point(size=5) +
  scale_color_manual(values=c("Untreated"="#2E86AB", "Dexamethasone"="#E84855")) +
  labs(title="PCA: Dexamethasone vs Untreated",
       subtitle="Airway Smooth Muscle Cells | chr1 only") +
  theme_bw(base_size=14)
ggsave("results/deseq2/pca_plot.pdf", pca_plot, width=8, height=6)
cat("PCA plot saved\n")

volcano <- EnhancedVolcano(res_shrunk,
  lab      = rownames(res_shrunk),
  x        = "log2FoldChange",
  y        = "padj",
  title    = "Dexamethasone vs Untreated",
  subtitle = "Airway Smooth Muscle Cells",
  pCutoff  = 0.05,
  FCcutoff = 1.0,
  pointSize= 2,
  labSize  = 3,
  col      = c("grey60","grey60","#2E86AB","#E84855"),
  drawConnectors = TRUE)
ggsave("results/deseq2/volcano_plot.pdf", volcano, width=10, height=8)
cat("Volcano plot saved\n")
if(nrow(degs) >= 2) {
  top_genes <- rownames(degs)[1:min(30, nrow(degs))]
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))
  annotation_col <- data.frame(
    Condition = coldata$condition,
    row.names = rownames(coldata)
  )
  ann_colors <- list(Condition=c(Untreated="#2E86AB", Dexamethasone="#E84855"))
  pheatmap(mat,
           annotation_col   = annotation_col,
           annotation_colors = ann_colors,
           show_rownames    = TRUE,
           fontsize_row     = 8,
           cluster_cols     = TRUE,
           cluster_rows     = TRUE,
           color = colorRampPalette(c("#2E86AB","white","#E84855"))(100),
           main  = "Top DEGs | Dexamethasone vs Untreated",
           filename = "results/deseq2/heatmap.pdf",
           width=10, height=10)
  cat("Heatmap saved\n")
} else {
  cat("Not enough DEGs for heatmap\n")
}

cat("Analysis complete!\n")
