## RNA–ATAC Integration: Improved Version with Temporal Analysis
## Orbicella faveolata Heat Stress
#Explicit temporal lag analysis (ATAC→RNA and RNA→ATAC)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(writexl)
})

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename


# CONFIGURATION - Key change: relaxed ATAC thresholds

config <- list(
  # RNA-seq thresholds (keep stringent - lots of signal)
  rna_padj = 0.1,
  rna_log2fc = 0.5,
  
  # ATAC thresholds - RELAXED to capture more peaks
  atac_padj = 0.2,           # Relaxed from 0.1
  atac_log2fc = 0.3,         # Relaxed from 0.5
  
  # Alternative: use nominal p-value for ATAC
  use_nominal_p_atac = FALSE,  # Set TRUE to use pvalue instead of padj
  atac_pval = 0.05             # Only used if use_nominal_p_atac = TRUE
)

cat("=== Configuration ===\n")
cat("RNA: padj <", config$rna_padj, ", |log2FC| >", config$rna_log2fc, "\n")
if (config$use_nominal_p_atac) {
  cat("ATAC: pvalue <", config$atac_pval, ", |log2FC| >", config$atac_log2fc, "(using nominal p)\n")
} else {
  cat("ATAC: padj <", config$atac_padj, ", |log2FC| >", config$atac_log2fc, "\n")
}


# 1. LOAD ATAC DATA

cat("\n=== Loading ATAC data ===\n")

atac_counts <- read_csv("COUNTSATAC.csv", show_col_types = FALSE)

# Clean up data types
atac_counts <- atac_counts %>%
  mutate(across(where(is.list), ~ sapply(.x, function(x) {
    if (length(x) == 0) NA_character_
    else if (length(x) == 1) as.character(x)
    else paste(as.character(x), collapse = ";")
  }))) %>%
  mutate(
    chr = as.character(chr),
    start = as.numeric(start),
    end = as.numeric(end),
    peak_uid = paste0(chr, ":", start, "-", end),
    gene_id = str_remove(gene_id, "^gene:")
  ) %>%
  distinct(peak_uid, .keep_all = TRUE)

cat("Total ATAC peaks:", nrow(atac_counts), "\n")
cat("Peaks with gene annotations:", sum(!is.na(atac_counts$gene_id) & atac_counts$gene_id != ""), "\n")

# Load metadata
atac_meta <- read_csv("ATAC_sample_metadata.csv", show_col_types = FALSE)

# Prepare count matrix
sample_cols <- grep("^OFav_", colnames(atac_counts), value = TRUE)
atac_matrix <- as.matrix(atac_counts[, sample_cols])
rownames(atac_matrix) <- atac_counts$peak_uid
mode(atac_matrix) <- "integer"

atac_meta <- atac_meta %>%
  filter(sample %in% colnames(atac_matrix)) %>%
  column_to_rownames("sample")
atac_meta <- atac_meta[colnames(atac_matrix), ]
atac_meta$block <- factor(atac_meta$block)
atac_meta$treatment <- factor(atac_meta$treatment, levels = c("Control", "Heat"))


# 2. RUN DESeq2 FOR ATAC - PER TIMEPOINT

cat("\n=== Running DESeq2 for ATAC ===\n")

atac_matrix_filtered <- atac_matrix[rowSums(atac_matrix >= 5) >= 2, ]
cat("Peaks after count filter:", nrow(atac_matrix_filtered), "\n")

blocks <- c("4", "5", "12", "24")
block_to_time <- c("4" = "4h", "5" = "30min", "12" = "12h", "24" = "24h")

atac_results_list <- list()

for (bl in blocks) {
  samples_block <- rownames(atac_meta)[atac_meta$block == bl]
  
  if (length(samples_block) < 2) next
  
  treatments <- atac_meta[samples_block, "treatment"]
  if (sum(treatments == "Heat") < 1 || sum(treatments == "Control") < 1) next
  
  tryCatch({
    counts_block <- atac_matrix_filtered[, samples_block, drop = FALSE]
    meta_block <- atac_meta[samples_block, , drop = FALSE]
    meta_block$treatment <- factor(meta_block$treatment, levels = c("Control", "Heat"))
    
    dds_sub <- DESeqDataSetFromMatrix(
      countData = counts_block,
      colData = meta_block,
      design = ~ treatment
    )
    
    keep_sub <- rowSums(counts(dds_sub) >= 5) >= 1
    dds_sub <- dds_sub[keep_sub, ]
    
    if (nrow(dds_sub) == 0) next
    
    dds_sub <- DESeq(dds_sub, quiet = TRUE)
    res <- results(dds_sub, contrast = c("treatment", "Heat", "Control"))
    
    res_df <- as.data.frame(res) %>%
      rownames_to_column("peak_uid") %>%
      mutate(time = block_to_time[bl])
    
    atac_results_list[[bl]] <- res_df
    cat("  Block", bl, "(", block_to_time[bl], "):", nrow(res_df), "peaks tested\n")
    
  }, error = function(e) {
    cat("  Block", bl, "error:", conditionMessage(e), "\n")
  })
}

atac_all <- bind_rows(atac_results_list)

# Add gene annotations
peak_gene_map <- atac_counts %>%
  select(peak_uid, gene_id, chr, start, end, top_feature) %>%
  filter(!is.na(gene_id), gene_id != "")

atac_all <- atac_all %>%
  left_join(peak_gene_map, by = "peak_uid")

# Filter significant - using configured thresholds
if (config$use_nominal_p_atac) {
  atac_sig <- atac_all %>%
    filter(
      !is.na(pvalue),
      pvalue < config$atac_pval,
      abs(log2FoldChange) >= config$atac_log2fc,
      !is.na(gene_id)
    )
  cat("\nUsing nominal p-value threshold\n")
} else {
  atac_sig <- atac_all %>%
    filter(
      !is.na(padj),
      padj < config$atac_padj,
      abs(log2FoldChange) >= config$atac_log2fc,
      !is.na(gene_id)
    )
}

cat("Significant DA peaks:", nrow(atac_sig), "\n")
cat("Unique DA peaks:", n_distinct(atac_sig$peak_uid), "\n")
cat("Genes with DA peaks:", n_distinct(atac_sig$gene_id), "\n")


# 3. LOAD AND PROCESS RNA-SEQ DATA

cat("\n=== Loading RNA-seq data ===\n")

rna_counts <- read_csv("countsRNAseq.csv", show_col_types = FALSE)
rna_meta <- read_csv("metadatarnaseqcoral.csv", show_col_types = FALSE)

colnames(rna_counts)[1] <- "gene_id"
rna_counts <- rna_counts %>%
  mutate(gene_id = str_remove(gene_id, "^gene:"))

rna_sample_cols <- colnames(rna_counts)[-1]
rna_matrix <- as.matrix(rna_counts[, rna_sample_cols])
rownames(rna_matrix) <- rna_counts$gene_id
mode(rna_matrix) <- "integer"

rna_meta <- rna_meta %>%
  column_to_rownames("sample_id")
rna_meta <- rna_meta[colnames(rna_matrix), ]
rna_meta$condition <- factor(rna_meta$condition, levels = c("Control", "Heat"))
rna_meta$time <- factor(rna_meta$time)


# 4. RUN DESeq2 FOR RNA-SEQ BY TIMEPOINT

cat("\n=== Running DESeq2 for RNA-seq ===\n")

rna_results_list <- list()

for (tp in unique(rna_meta$time)) {
  samples_tp <- rownames(rna_meta)[rna_meta$time == tp]
  
  if (length(samples_tp) >= 2) {
    rna_sub <- rna_matrix[, samples_tp, drop = FALSE]
    meta_sub <- rna_meta[samples_tp, , drop = FALSE]
    
    keep <- rowSums(rna_sub >= 10) >= 2
    rna_sub <- rna_sub[keep, , drop = FALSE]
    
    if (length(unique(meta_sub$condition)) == 2 && nrow(rna_sub) > 0) {
      dds_rna <- DESeqDataSetFromMatrix(
        countData = rna_sub,
        colData = meta_sub,
        design = ~ condition
      )
      dds_rna <- DESeq(dds_rna, quiet = TRUE)
      res <- results(dds_rna, contrast = c("condition", "Heat", "Control"))
      
      res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        mutate(time = tp)
      
      rna_results_list[[tp]] <- res_df
      cat("  ", tp, ":", sum(res_df$padj < 0.1, na.rm = TRUE), "DE genes\n")
    }
  }
}

rna_all <- bind_rows(rna_results_list)

rna_sig <- rna_all %>%
  filter(
    !is.na(padj),
    padj < config$rna_padj,
    abs(log2FoldChange) >= config$rna_log2fc
  )

cat("Significant DE genes (total entries):", nrow(rna_sig), "\n")
cat("Unique DE genes:", n_distinct(rna_sig$gene_id), "\n")


# 5. TEMPORAL LAG ANALYSIS - KEY IMPROVEMENT

cat("\n=== Temporal Lag Analysis ===\n")

# Define time order for lag analysis
time_order <- c("30min", "4h", "12h", "24h")
time_numeric <- c("30min" = 0.5, "4h" = 4, "12h" = 12, "24h" = 24)

# Function to determine temporal relationship
get_temporal_category <- function(atac_times, rna_times) {
  atac_earliest <- min(time_numeric[atac_times])
  rna_earliest <- min(time_numeric[rna_times])
  
  if (atac_earliest < rna_earliest) {
    return("ATAC_first")
  } else if (rna_earliest < atac_earliest) {
    return("RNA_first")
  } else {
    return("Simultaneous")
  }
}

# Strategy A: Same-timepoint links (both significant at same time)
cat("\n--- Strategy A: Same-timepoint ---\n")

link_same_time <- atac_sig %>%
  inner_join(
    rna_sig %>% select(gene_id, time, rna_log2FC = log2FoldChange, rna_padj = padj),
    by = c("gene_id", "time"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    concordant = sign(log2FoldChange) == sign(rna_log2FC),
    link_type = "same_timepoint",
    temporal_category = "Simultaneous"
  ) %>%
  select(
    peak_uid, gene_id, time, chr, start, end, top_feature,
    atac_log2FC = log2FoldChange, atac_padj = padj,
    rna_log2FC, rna_padj, concordant, link_type, temporal_category
  )

cat("Same-timepoint links:", nrow(link_same_time), "\n")

# Strategy B: Cross-timepoint with temporal tracking
cat("\n--- Strategy B: Cross-timepoint with temporal order ---\n")

# Aggregate by gene - track which timepoints are significant
genes_de <- rna_sig %>%
  group_by(gene_id) %>%
  summarise(
    rna_log2FC = log2FoldChange[which.max(abs(log2FoldChange))],
    rna_padj = min(padj),
    rna_times = list(time),
    rna_earliest_time = time[which.min(time_numeric[time])],
    .groups = "drop"
  )

peaks_da <- atac_sig %>%
  group_by(peak_uid, gene_id, chr, start, end, top_feature) %>%
  summarise(
    atac_log2FC = log2FoldChange[which.max(abs(log2FoldChange))],
    atac_padj = min(padj),
    atac_times = list(time),
    atac_earliest_time = time[which.min(time_numeric[time])],
    .groups = "drop"
  )

link_cross_time <- peaks_da %>%
  inner_join(genes_de, by = "gene_id") %>%
  mutate(
    concordant = sign(atac_log2FC) == sign(rna_log2FC),
    link_type = "cross_timepoint",
    # Determine temporal relationship
    temporal_category = map2_chr(atac_times, rna_times, function(at, rt) {
      at_min <- min(time_numeric[unlist(at)])
      rt_min <- min(time_numeric[unlist(rt)])
      if (at_min < rt_min) "ATAC_first"
      else if (rt_min < at_min) "RNA_first"
      else "Simultaneous"
    }),
    # Convert list columns to strings for export
    atac_times_str = map_chr(atac_times, ~ paste(.x, collapse = ",")),
    rna_times_str = map_chr(rna_times, ~ paste(.x, collapse = ","))
  ) %>%
  select(-atac_times, -rna_times) %>%
  rename(atac_times = atac_times_str, rna_times = rna_times_str)

cat("Cross-timepoint links:", nrow(link_cross_time), "\n")
cat("  - ATAC first:", sum(link_cross_time$temporal_category == "ATAC_first"), "\n")
cat("  - RNA first:", sum(link_cross_time$temporal_category == "RNA_first"), "\n")
cat("  - Simultaneous:", sum(link_cross_time$temporal_category == "Simultaneous"), "\n")

# Strategy C: Specific temporal lag - ATAC at T1 -> RNA at T2
cat("\n--- Strategy C: Specific temporal lags ---\n")

lag_pairs <- list(
  c("30min", "4h"),
  c("30min", "12h"),
  c("4h", "12h"),
  c("4h", "24h"),
  c("12h", "24h")
)

link_lag <- tibble()

for (lag in lag_pairs) {
  t1 <- lag[1]
  t2 <- lag[2]
  
  atac_t1 <- atac_sig %>% filter(time == t1)
  rna_t2 <- rna_sig %>% filter(time == t2)
  
  lag_link <- atac_t1 %>%
    inner_join(
      rna_t2 %>% select(gene_id, rna_log2FC = log2FoldChange, rna_padj = padj),
      by = "gene_id",
      relationship = "many-to-many"
    ) %>%
    mutate(
      concordant = sign(log2FoldChange) == sign(rna_log2FC),
      link_type = paste0("lag_", t1, "_to_", t2),
      temporal_category = "ATAC_first",
      atac_time = t1,
      rna_time = t2
    ) %>%
    select(
      peak_uid, gene_id, atac_time, rna_time, chr, start, end, top_feature,
      atac_log2FC = log2FoldChange, atac_padj = padj,
      rna_log2FC, rna_padj, concordant, link_type, temporal_category
    )
  
  if (nrow(lag_link) > 0) {
    cat("  ", t1, "->", t2, ":", nrow(lag_link), "links\n")
    link_lag <- bind_rows(link_lag, lag_link)
  }
}


# 6. CANDIDATE SCORING - IMPROVED

cat("\n=== Candidate Scoring ===\n")

# Score based on:
# 1. Statistical significance (both ATAC and RNA)
# 2. Effect size magnitude
# 3. Concordance
# 4. Temporal plausibility (ATAC first = more plausible for regulation)
# 5. Promoter location (more direct regulatory link)

score_candidates <- function(df) {
  df %>%
    mutate(
      # Significance score
      sig_score = -log10(atac_padj) - log10(rna_padj),
      
      # Effect size score
      effect_score = abs(atac_log2FC) + abs(rna_log2FC),
      
      # Concordance bonus
      concordance_score = ifelse(concordant, 5, 0),
      
      # Temporal plausibility (ATAC first is biologically expected)
      temporal_score = case_when(
        temporal_category == "ATAC_first" ~ 3,
        temporal_category == "Simultaneous" ~ 2,
        temporal_category == "RNA_first" ~ 1
      ),
      
      # Promoter bonus
      promoter_score = ifelse(top_feature == "promoter", 3, 0),
      
      # Combined score
      candidate_score = sig_score + effect_score * 2 + concordance_score + 
                        temporal_score + promoter_score
    )
}

link_cross_scored <- score_candidates(link_cross_time)


# 7. IDENTIFY OVEREXPRESSION CANDIDATES

cat("\n=== Overexpression Candidates ===\n")

overexpress_candidates <- link_cross_scored %>%
  filter(
    concordant == TRUE,
    rna_log2FC > 0,
    atac_log2FC > 0
  ) %>%
  arrange(desc(candidate_score))

cat("Concordant upregulated candidates:", nrow(overexpress_candidates), "\n")

if (nrow(overexpress_candidates) > 0) {
  cat("\nTop 10 overexpression candidates:\n")
  overexpress_candidates %>%
    head(10) %>%
    select(gene_id, top_feature, temporal_category, rna_log2FC, atac_log2FC, candidate_score) %>%
    print()
}


# 8. SUMMARY STATISTICS

cat("\n=== Summary Statistics ===\n")

summary_stats <- tibble(
  Metric = c(
    "Total ATAC peaks analyzed",
    "Significant DA peaks (relaxed threshold)",
    "Unique genes with DA peaks",
    "Significant DE genes",
    "Cross-timepoint links",
    "Same-timepoint links",
    "ATAC-first links",
    "RNA-first links", 
    "Simultaneous links",
    "Concordant upregulated candidates",
    "Concordance rate (%)"
  ),
  Value = c(
    nrow(atac_counts),
    n_distinct(atac_sig$peak_uid),
    n_distinct(atac_sig$gene_id),
    n_distinct(rna_sig$gene_id),
    nrow(link_cross_time),
    nrow(link_same_time),
    sum(link_cross_time$temporal_category == "ATAC_first"),
    sum(link_cross_time$temporal_category == "RNA_first"),
    sum(link_cross_time$temporal_category == "Simultaneous"),
    nrow(overexpress_candidates),
    round(mean(link_cross_time$concordant) * 100, 1)
  )
)

print(summary_stats)


# 9. EXPORT RESULTS

cat("\n=== Exporting Results ===\n")

write_xlsx(
  list(
    Summary = summary_stats,
    Overexpression_Candidates = overexpress_candidates,
    Cross_Timepoint_Links = link_cross_scored,
    Same_Timepoint_Links = link_same_time,
    Temporal_Lag_Links = link_lag,
    ATAC_DA_peaks = atac_sig,
    RNA_DE_genes = rna_sig
  ),
  "RNA_ATAC_integration_v2_results.xlsx"
)

cat("\nResults exported to: RNA_ATAC_integration_v2_results.xlsx\n")
cat("\n=== COMPLETE ===\n")
