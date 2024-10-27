library(edgeR)

counts_human <- read.csv("U_Sirius/transcriptomics/hw1/counts_human.csv", row.names = 1)
counts_chimp <- read.csv("U_Sirius/transcriptomics/hw1/counts_chimp.csv", row.names = 1)

counts <- cbind(counts_human, counts_chimp)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")

norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)

log_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
write.csv(log_cpm, "U_Sirius/transcriptomics/hw1/counts_TMM_log_cpm.csv", row.names = TRUE)
