\
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  source("scripts/00_setup/pipeline_utils.R")
})

cfg <- read_config()
make_dirs(cfg)

if (!file.exists(cfg$paths$atac_counts)) stop("Missing ATAC count matrix: ", cfg$paths$atac_counts)
if (!file.exists(cfg$paths$atac_metadata)) stop("Missing ATAC metadata: ", cfg$paths$atac_metadata)

counts_df <- read_csv(cfg$paths$atac_counts, show_col_types = FALSE)
meta_raw <- read_csv(cfg$paths$atac_metadata, show_col_types = FALSE)
meta <- standardize_metadata(meta_raw, cfg)
cleaned <- clean_count_matrix(counts_df, meta)
count_mat <- cleaned$count_mat
meta <- cleaned$meta

all_results <- list()

for (tp in cfg$design$time_levels) {
  meta_tp <- meta %>% filter(time == tp)
  if (nrow(meta_tp) < 2 || length(unique(meta_tp$condition)) < 2) {
    warning("Skipping ATAC time point due to insufficient condition groups: ", tp)
    next
  }

  mat_tp <- count_mat[, meta_tp$sample, drop = FALSE]

  dds_tp <- DESeqDataSetFromMatrix(
    countData = round(mat_tp),
    colData = as.data.frame(meta_tp) %>% column_to_rownames("sample"),
    design = ~ condition
  )

  dds_tp <- dds_tp[rowSums(counts(dds_tp)) >= cfg$thresholds$min_count_sum, ]
  dds_tp <- DESeq(dds_tp)

  res <- results(dds_tp, contrast = c("condition", cfg$design$heat_label, cfg$design$control_label))

  res_tbl <- as.data.frame(res) %>%
    rownames_to_column("peak_id") %>%
    as_tibble() %>%
    mutate(
      time = as.character(tp),
      contrast = paste0(cfg$design$heat_label, "_vs_", cfg$design$control_label),
      direction = classify_direction(log2FoldChange, padj, cfg$thresholds$atac_padj, cfg$thresholds$atac_log2fc)
    ) %>%
    arrange(padj)

  all_results[[as.character(tp)]] <- res_tbl
}

atac_da <- bind_rows(all_results)
write_csv(atac_da, file.path(cfg$outputs$tables, "ATAC_DA_HeatvsControl_Timepoints.csv"))

dds_all <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = as.data.frame(meta) %>% column_to_rownames("sample"),
  design = ~ time + condition
)
dds_all <- dds_all[rowSums(counts(dds_all)) >= cfg$thresholds$min_count_sum, ]
dds_all <- estimateSizeFactors(dds_all)
norm_counts <- counts(dds_all, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("peak_id")
write_csv(norm_counts, file.path(cfg$outputs$tables, "ATAC_normalized_counts.csv"))

vsd <- vst(dds_all, blind = TRUE)
vst_mat <- assay(vsd) %>%
  as.data.frame() %>%
  rownames_to_column("peak_id")
write_csv(vst_mat, file.path(cfg$outputs$tables, "ATAC_vst_counts.csv"))

write_session_info()
message("Saved ATAC DA, normalized counts, and VST counts.")
