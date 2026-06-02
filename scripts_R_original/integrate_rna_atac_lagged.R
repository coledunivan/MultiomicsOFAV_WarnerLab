\
suppressPackageStartupMessages({
  library(tidyverse)
  source("scripts/00_setup/pipeline_utils.R")
})

cfg <- read_config()
make_dirs(cfg)

same_file <- file.path(cfg$outputs$tables, "RNA_ATAC_integrated_same_time.csv")
rna_file <- file.path(cfg$outputs$tables, "RNA_DE_HeatvsControl_Timepoints.csv")
atac_file <- file.path(cfg$outputs$tables, "ATAC_DA_HeatvsControl_Timepoints.csv")
peak_annot_file <- cfg$paths$peak_annotation

if (!file.exists(rna_file)) stop("Missing RNA DE file.")
if (!file.exists(atac_file)) stop("Missing ATAC DA file.")
if (!file.exists(peak_annot_file)) stop("Missing peak annotation file.")

rna <- read_csv(rna_file, show_col_types = FALSE) %>%
  rename(time_rna = time, rna_log2FC = log2FoldChange, rna_padj = padj, rna_direction = direction)

atac <- read_csv(atac_file, show_col_types = FALSE) %>%
  rename(time_atac = time, atac_log2FC = log2FoldChange, atac_padj = padj, atac_direction = direction)

peak_annot <- read_csv(peak_annot_file, show_col_types = FALSE)

atac_gene <- atac %>%
  left_join(peak_annot, by = "peak_id") %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id, time_atac) %>%
  summarise(
    n_peaks = n(),
    atac_log2FC = atac_log2FC[which.max(abs(atac_log2FC))],
    atac_padj = min(atac_padj, na.rm = TRUE),
    atac_direction = case_when(
      any(atac_direction == "up") ~ "open",
      any(atac_direction == "down") ~ "closed",
      TRUE ~ "ns"
    ),
    peak_id_top = peak_id[which.max(abs(atac_log2FC))],
    .groups = "drop"
  )

lag_pairs <- map_dfr(cfg$integration$atac_leads_rna, ~ tibble(time_atac = .x[[1]], time_rna = .x[[2]]))

lagged <- lag_pairs %>%
  left_join(atac_gene, by = "time_atac") %>%
  inner_join(rna, by = c("gene_id", "time_rna")) %>%
  mutate(
    lag_class = case_when(
      atac_direction == "open" & rna_direction == "up" ~ "ATAC_open_before_RNA_up",
      atac_direction == "closed" & rna_direction == "down" ~ "ATAC_closed_before_RNA_down",
      atac_direction == "open" & rna_direction == "down" ~ "ATAC_open_before_RNA_down",
      atac_direction == "closed" & rna_direction == "up" ~ "ATAC_closed_before_RNA_up",
      TRUE ~ "other"
    )
  )

write_csv(lagged, file.path(cfg$outputs$tables, "RNA_ATAC_integrated_lagged.csv"))

summary_tbl <- lagged %>%
  count(time_atac, time_rna, lag_class) %>%
  arrange(time_atac, time_rna, desc(n))
write_csv(summary_tbl, file.path(cfg$outputs$tables, "RNA_ATAC_integrated_lagged_summary.csv"))

message("Saved lagged RNA-ATAC integration.")
