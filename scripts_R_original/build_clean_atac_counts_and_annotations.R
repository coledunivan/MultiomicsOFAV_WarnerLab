\
suppressPackageStartupMessages({
  library(tidyverse)
  source("scripts/00_setup/pipeline_utils.R")
})

cfg <- read_config()
make_dirs(cfg)

infile <- cfg$paths$atac_counts
if (!file.exists(infile)) stop("Missing ATAC counts/annotation file: ", infile)

raw <- read_csv(infile, show_col_types = FALSE)

required <- c("chr", "start", "end", "peak_id")
missing_required <- setdiff(required, names(raw))
if (length(missing_required) > 0) {
  stop("ATAC file is missing required columns: ", paste(missing_required, collapse = ", "))
}

sample_cols <- names(raw)[stringr::str_detect(names(raw), "^OFav_")]
if (length(sample_cols) == 0) {
  stop("No ATAC sample columns detected. Expected columns starting with OFav_.")
}

# 1. Clean peak count matrix: one row per peak.
# COUNTS_ATAC_or_peak_counts.csv has repeated peak rows because each peak can map to multiple GO terms.
# DESeq2 needs exactly one count row per peak, so we keep the first unique count row per peak.
peak_counts <- raw %>%
  select(chr, start, end, peak_id, all_of(sample_cols)) %>%
  distinct(peak_id, .keep_all = TRUE)

write_csv(peak_counts, "data/raw_counts/ATAC_peak_counts_clean.csv")

# 2. Peak-to-gene annotation table for RNA-ATAC integration.
peak_to_gene_long <- raw %>%
  select(any_of(c("peak_id", "chr", "start", "end", "top_feature", "gene_id"))) %>%
  distinct()

peak_to_gene_compact <- peak_to_gene_long %>%
  group_by(peak_id, chr, start, end) %>%
  summarise(
    top_feature = paste(sort(unique(na.omit(top_feature))), collapse = ";"),
    gene_id = paste(sort(unique(na.omit(gene_id))), collapse = ";"),
    .groups = "drop"
  )

write_csv(peak_to_gene_compact, "results/tables/ATAC_peak_to_gene_annotation.csv")

# 3. Peak-to-gene-GO annotation table for enrichment and interpretation.
peak_go <- raw %>%
  select(any_of(c("peak_id", "chr", "start", "end", "top_feature", "gene_id", "go_id", "go_term", "aspect"))) %>%
  distinct()

write_csv(peak_go, "results/tables/ATAC_peak_to_gene_GO_annotation.csv")

# 4. Gene-to-GO file usable by enrichment scripts.
gene2go <- raw %>%
  select(any_of(c("gene_id", "go_id", "go_term", "aspect"))) %>%
  filter(!is.na(gene_id), !is.na(go_id)) %>%
  distinct()

write_csv(gene2go, "data/reference/gene2go_from_ATAC_counts.csv")

message("Saved:")
message(" - data/raw_counts/ATAC_peak_counts_clean.csv")
message(" - results/tables/ATAC_peak_to_gene_annotation.csv")
message(" - results/tables/ATAC_peak_to_gene_GO_annotation.csv")
message(" - data/reference/gene2go_from_ATAC_counts.csv")
