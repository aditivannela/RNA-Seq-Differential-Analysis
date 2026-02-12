library(DESeq2)
library(ggplot2)
library(pheatmap)

#load file
counts <- read.table("GSE263150_RNA-seq_gene_expression.txt",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t")
head(counts)

#create metadata

condition <- factor(c("DMSO","DMSO","DMSO","KAT8i","KAT8i","KAT8i"))

colData <- data.frame(row.names = colnames(counts),
                      condition)

#create dds

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = colData, 
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >10, ] #remove low count genes

#run differential expression
dds <- DESeq(dds)

res <- results(dds)
summary(res)

resLFC <- lfcShrink(dds, 
                   coef = "condition_KAT8i_vs_DMSO",
                   type = "apeglm") #shrink log2 fold changes, to see coefficients use resultsNames(dds)

#view significant genes

resOrdered <- resLFC[order(resLFC$padj), ]

head(resOrdered)
sig <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
cat("Number of significant genes:", nrow(sig), "\n")

# Save Results

write.csv(as.data.frame(resLFC), 
          file = "DESeq2_all_results_KAT8iP_vs_DMSO.csv")

write.csv(as.data.frame(sig),
          file = "DESeq2_significant_genes_padj0.05_LFC1.csv")
# PCA plot

vsd <- vst(dds, blind = FALSE)

pdf("PCA_plot.pdf")
plotPCA(vsd, intgroup = "condition")
dev.off()

# Volcano Plot

res_df <- as.data.frame(resLFC)

pdf("Volcano_plot.pdf")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  xlim(c(-5, 5))
dev.off()

# Heatmap of Top Genes

topgenes <- head(rownames(resOrdered), 50)
mat <- assay(vsd)[topgenes, ]

pdf("Heatmap_top_genes.pdf")
pheatmap(mat, 
         cluster_rows = TRUE,
         show_rownames = FALSE, 
         cluster_cols = TRUE, 
         annotation_col = colData)
dev.off()

# Session Information

sessionInfo()
