library(DESeq2)
library(tidyverse)
library(ggplot2)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ time + condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("condition", "time"), returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = time)) +
  geom_point(size = 4, alpha = 0.9) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of Orbicella faveolata RNA-seq",
    subtitle = "Heat vs Control across time points"
  )

p

ggsave(
  filename = "PCA_heat_vs_control.png",
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)

