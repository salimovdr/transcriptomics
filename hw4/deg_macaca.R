library(DESeq2)
library(ggplot2)

setwd("/home/salimovdr/U_Sirius/transcriptomics/hw4")

counts_df <- read.delim("data/counts_postprocessing.tsv", sep = "\t", row.names = 1)

sample_info <- data.frame(
  sample = colnames(counts_df),
  species = factor(ifelse(grepl("human", colnames(counts_df)), "human", "chimp"), 
                   levels = c("human", "chimp"))
)

dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = sample_info, design = ~ species)
dds <- DESeq(dds)

res <- results(dds, contrast = c("species", "human", "chimp"))
head(res)

vsd <- vst(dds, blind = TRUE)  

plotPCA(vsd, intgroup = "species") + ggtitle("PCA of RNA-seq Samples by Species")

cooks <- assays(dds)[["cooks"]]
boxplot(log10(cooks + 1), main = "Cook's Distance for Each Sample",
        ylab = "Log10(Cook's Distance)")

threshold <- quantile(cooks, 0.8, na.rm = TRUE)
outliers <- which(rowMeans(cooks, na.rm = TRUE) > threshold)
outliers_samples <- colnames(dds)[outliers]

dds_filtered <- dds[, !(colnames(dds) %in% outliers_samples)]
dds_filtered <- DESeq(dds_filtered)
res_filtered <- results(dds_filtered, contrast = c("species", "human", "chimp"))
head(res_filtered)

log2FC_threshold <- 1
padj_threshold <- 0.05  
max_neg_log10_padj <- 300  

res_filtered$neg_log10_padj <- pmin(-log10(res_filtered$padj), max_neg_log10_padj)

calcium_genes <- scan("data/CALCIUM_SIGNALING.txt", what = "", sep = "\t")
calcium_genes <- toupper(calcium_genes)

res_filtered$in_calcium_signaling <- rownames(res_filtered) %in% calcium_genes


ggplot(res_filtered, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(data = subset(res_filtered, !in_calcium_signaling),
             aes(color = ifelse(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold, "red", "blue")),
             alpha = 0.6) +
  geom_point(data = subset(res_filtered, in_calcium_signaling), color = "green", alpha = 0.8) +
  scale_color_identity() +  
  xlab("Log2 Fold Change") + ylab("-Log10 Adjusted P-Value") +
  ggtitle("Volcano Plot of Differential Gene Expression: Human vs Chimpanzee") +
  theme_minimal() +
  theme(legend.position = "none") +
  ylim(0, max_neg_log10_padj + 10) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey")


write.table(as.data.frame(res_filtered), file = "data/deg.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
