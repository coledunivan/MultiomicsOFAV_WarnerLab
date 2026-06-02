\
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  source("scripts/00_setup/pipeline_utils.R")
})

cfg <- read_config()
make_dirs(cfg)

plot_pca_assay <- function(count_path, meta_path, assay, id_name) {
  if (!file.exists(count_path) || !file.exists(meta_path)) {
    warning("Skipping ", assay, ": missing count or metadata file.")
    return(NULL)
  }

  counts_df <- read_csv(count_path, show_col_types = FALSE)
  meta_raw <- read_csv(meta_path, show_col_types = FALSE)
  meta <- standardize_metadata(meta_raw, cfg)
  cleaned <- clean_count_matrix(counts_df, meta)
  count_mat <- cleaned$count_mat
  meta <- cleaned$meta

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = as.data.frame(meta) %>% column_to_rownames("sample"),
    design = ~ time + condition
  )
  dds <- dds[rowSums(counts(dds)) >= cfg$thresholds$min_count_sum, ]
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE)

  pca <- plotPCA(vsd, intgroup = c("time", "condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pca, "percentVar"))

  p <- ggplot(pca, aes(PC1, PC2, color = condition, shape = time, label = name)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
    labs(
      title = paste0(assay, " PCA"),
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    ) +
    theme_bw()

  ggsave(file.path(cfg$outputs$figures, paste0(assay, "_PCA.png")), p, width = 7, height = 5, dpi = 300)

  sample_dist <- dist(t(assay(vsd)))
  dist_mat <- as.matrix(sample_dist)
  write_csv(as.data.frame(dist_mat) %>% rownames_to_column("sample"), file.path(cfg$outputs$qc, paste0(assay, "_sample_distance_matrix.csv")))

  invisible(p)
}

plot_pca_assay(cfg$paths$rna_counts, cfg$paths$rna_metadata, "RNA", "gene_id")
plot_pca_assay(cfg$paths$atac_counts, cfg$paths$atac_metadata, "ATAC", "peak_id")

write_session_info()
message("PCA QC complete.")
