suppressPackageStartupMessages({
  library(tidyverse)
})

# Placeholder module.
# Expected input after integration:
# results/tables/RNA_ATAC_integrated_lagged.csv
# Required columns should include:
# gene_id, gene_symbol, time_rna, time_atac, rna_log2FC, rna_padj,
# atac_log2FC, atac_padj, peak_id, promoter_distance_bp

infile <- "results/tables/RNA_ATAC_integrated_lagged.csv"
outfile <- "results/tables/candidate_heat_inducible_promoters_ranked.csv"

if (!file.exists(infile)) {
  stop("Missing input: ", infile, ". Run integration module first.")
}

df <- read_csv(infile, show_col_types = FALSE)

ranked <- df %>%
  mutate(
    rna_score = pmax(rna_log2FC, 0) * -log10(pmax(rna_padj, 1e-300)),
    atac_score = pmax(atac_log2FC, 0) * -log10(pmax(atac_padj, 1e-300)),
    distance_score = ifelse(!is.na(promoter_distance_bp), 1 / (1 + abs(promoter_distance_bp) / 1000), 0.5),
    candidate_score = rna_score + atac_score + distance_score
  ) %>%
  arrange(desc(candidate_score))

write_csv(ranked, outfile)
message("Saved: ", outfile)
