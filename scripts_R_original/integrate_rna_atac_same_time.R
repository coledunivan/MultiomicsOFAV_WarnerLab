\
suppressPackageStartupMessages({
  library(tidyverse)
  source("scripts/00_setup/pipeline_utils.R")
})

cfg <- read_config()
make_dirs(cfg)

rna_file <- file.path(cfg$outputs$tables, "RNA_DE_HeatvsControl_Timepoints.csv")
atac_file <- file.path(cfg$outputs$tables, "ATAC_DA_HeatvsControl_Timepoints.csv")
peak_annot_file <- cfg$paths$peak_annotation

if (!file.exists(rna_file)) stop("Missing RNA DE file: ", rna_file)
if (!file.exists(atac_file)) stop("Missing ATAC DA file: ", atac_file)

rna <- read_csv(rna_file, show_col_types = FALSE) %>%
  rename(rna_log2FC = log2FoldChange, rna_padj = padj, rna_direction = direction)

atac <- read_csv(atac_file, show_col_types = FALSE) %>%
  rename(atac_log2FC = log2FoldChange, atac_padj = padj, atac_direction = direction)

# Need a peak-to-gene annotation. If absent, try to infer gene_id from peak table if present.
if (file.exists(peak_annot_file)) {
  peak_annot <- read_csv(peak_annot_file, show_col_types = FALSE)
} else if ("gene_id" %in% names(atac)) {
  peak_annot <- atac %>% select(peak_id, gene_id) %>% distinct()
  warning("Using gene_id already present in ATAC table as peak annotation.")
} else {
  stop("Missing peak-to-gene annotation file: ", peak_annot_file,
       ". Create results/tables/ATAC_peak_to_gene_annotation.csv with columns peak_id and gene_id.")
}

stopifnot(all(c("peak_id", "gene_id") %in% names(peak_annot)))
peak_annot <- peak_annot %>% tidyr::separate_rows(gene_id, sep = ";") %>% filter(!is.na(gene_id), gene_id != "")

atac_gene <- atac %>%
  left_join(peak_annot, by = "peak_id") %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id, time) %>%
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

integrated <- rna %>%
  select(gene_id, time, rna_log2FC, rna_padj, rna_direction, everything()) %>%
  left_join(atac_gene, by = c("gene_id", "time")) %>%
  mutate(
    integrated_class = case_when(
      rna_direction == "up" & atac_direction == "open" ~ "RNA_up_ATAC_open",
      rna_direction == "down" & atac_direction == "closed" ~ "RNA_down_ATAC_closed",
      rna_direction == "up" & atac_direction == "closed" ~ "RNA_up_ATAC_closed",
      rna_direction == "down" & atac_direction == "open" ~ "RNA_down_ATAC_open",
      rna_direction %in% c("up", "down") & (is.na(atac_direction) | atac_direction == "ns") ~ "RNA_only",
      !(rna_direction %in% c("up", "down")) & atac_direction %in% c("open", "closed") ~ "ATAC_only",
      TRUE ~ "unchanged"
    )
  )

write_csv(integrated, file.path(cfg$outputs$tables, "RNA_ATAC_integrated_same_time.csv"))

summary_tbl <- integrated %>%
  count(time, integrated_class) %>%
  arrange(time, desc(n))
write_csv(summary_tbl, file.path(cfg$outputs$tables, "RNA_ATAC_integrated_same_time_summary.csv"))

message("Saved same-time RNA-ATAC integration.")
