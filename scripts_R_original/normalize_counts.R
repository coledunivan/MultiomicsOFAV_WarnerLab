# =============================================================================
# R Script: Normalize ATAC-seq and RNA-seq counts
# For Orbicella faveolata heat stress experiment
# =============================================================================

# Install/load required packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!require("edgeR", quietly = TRUE)) BiocManager::install("edgeR")

library(DESeq2)
library(edgeR)

# =============================================================================
# SET FILE PATHS - UPDATE THESE TO MATCH YOUR FILE LOCATIONS
# =============================================================================

# Input files
rnaseq_counts_file <- "countsRNAseq.csv"          # RNA-seq raw counts
atac_counts_file <- "COUNTSATAC.csv"              # ATAC-seq raw counts
rnaseq_metadata_file <- "metadatarnaseqcoral.csv" # RNA-seq metadata
atac_metadata_file <- "ATAC_sample_metadata.csv"  # ATAC-seq metadata

# Output files
rnaseq_normalized_output <- "RNAseq_normalized_counts.csv"
atac_normalized_output <- "ATACseq_normalized_counts.csv"

# =============================================================================
# FUNCTION: Normalize counts using DESeq2
# =============================================================================

normalize_counts_deseq2 <- function(counts_matrix, metadata, design_formula = ~ 1) {
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = metadata,
    design = design_formula
  )
  
  # Filter low counts (at least 10 counts in at least 2 samples)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  
  # Estimate size factors and normalize
  dds <- estimateSizeFactors(dds)
  
  # Get normalized counts (size-factor corrected)
  normalized_counts <- counts(dds, normalized = TRUE)
  
  return(list(
    normalized = normalized_counts,
    size_factors = sizeFactors(dds),
    dds = dds
  ))
}

# =============================================================================
# PROCESS RNA-seq DATA
# =============================================================================

cat("Processing RNA-seq data...\n")

# Read RNA-seq counts
rnaseq_counts <- read.csv(rnaseq_counts_file, row.names = 1, check.names = FALSE)

# Read metadata
rnaseq_meta <- read.csv(rnaseq_metadata_file)

# Check column names in counts
cat("RNA-seq count columns:\n")
print(colnames(rnaseq_counts))

# The counts columns are the full sample IDs (e.g., OFav_05_C2_S2_L001)
# Match directly to metadata sample_id column
rnaseq_meta_ordered <- rnaseq_meta[match(colnames(rnaseq_counts), rnaseq_meta$sample_id), ]

# Use sample_id as rownames (these should be unique)
rownames(rnaseq_meta_ordered) <- rnaseq_meta_ordered$sample_id

# Check for any NA (unmatched samples)
if (any(is.na(rnaseq_meta_ordered$sample_id))) {
  cat("WARNING: Some count columns did not match metadata!\n")
  cat("Unmatched columns:\n")
  print(colnames(rnaseq_counts)[is.na(rnaseq_meta_ordered$sample_id)])
}

cat("\nMetadata preview:\n")
print(rnaseq_meta_ordered)

# Ensure counts are integers
rnaseq_counts_int <- round(as.matrix(rnaseq_counts))

# Normalize using DESeq2
rnaseq_norm <- normalize_counts_deseq2(
  rnaseq_counts_int, 
  rnaseq_meta_ordered,
  design = ~ 1  # Use simple design for normalization only
)

# Export normalized counts
write.csv(rnaseq_norm$normalized, rnaseq_normalized_output, row.names = TRUE)
cat("\nRNA-seq normalized counts saved to:", rnaseq_normalized_output, "\n")
cat("  - Genes retained after filtering:", nrow(rnaseq_norm$normalized), "\n")
cat("  - Size factors:\n")
print(rnaseq_norm$size_factors)

# =============================================================================
# PROCESS ATAC-seq DATA
# =============================================================================

cat("\n\nProcessing ATAC-seq data...\n")

# Read ATAC-seq counts
atac_counts <- read.csv(atac_counts_file, check.names = FALSE)

# Check column names
cat("ATAC-seq columns:\n")
print(colnames(atac_counts))

# Get count columns (those starting with OFav)
count_cols <- grep("^OFav_", colnames(atac_counts), value = TRUE)

cat("\nATAC count columns:\n")
print(count_cols)

cat("\nTotal rows in ATAC file:", nrow(atac_counts), "\n")
cat("Unique peak_ids:", length(unique(atac_counts$peak_id)), "\n")

# ============================================================================
# IMPORTANT: Collapse duplicate rows (same peak_id with multiple GO terms)
# Keep only unique peaks - counts should be identical for same peak
# ============================================================================

# Get unique peaks only (first occurrence of each peak_id)
atac_counts_unique <- atac_counts[!duplicated(atac_counts$peak_id), ]
cat("Rows after removing duplicates:", nrow(atac_counts_unique), "\n")

# Extract counts matrix from unique peaks
atac_counts_matrix <- as.matrix(atac_counts_unique[, count_cols])
rownames(atac_counts_matrix) <- atac_counts_unique$peak_id

# Store peak info (just chr, start, end, peak_id, gene_id - not GO terms)
peak_info_cols <- c("chr", "start", "end", "peak_id", "gene_id", "top_feature")
available_info_cols <- intersect(colnames(atac_counts_unique), peak_info_cols)
peak_info <- atac_counts_unique[, available_info_cols, drop = FALSE]
rownames(peak_info) <- atac_counts_unique$peak_id

# Read ATAC metadata
atac_meta <- read.csv(atac_metadata_file)

# ATAC count columns are like: OFav_05_C1, OFav_05_C2, etc.
# ATAC metadata has sample like: OFav_05_C1_S1
# Need to strip _S# suffix from metadata to match

atac_meta$sample_base <- gsub("_S[0-9]+$", "", atac_meta$sample)

cat("\nATAC metadata sample bases:\n")
print(atac_meta$sample_base)

# Match count columns to metadata
atac_meta_ordered <- atac_meta[match(count_cols, atac_meta$sample_base), ]
rownames(atac_meta_ordered) <- count_cols  # Use the count column names as rownames

# Check for any NA (unmatched samples)
if (any(is.na(atac_meta_ordered$sample))) {
  cat("WARNING: Some count columns did not match metadata!\n")
  cat("Unmatched columns:\n")
  print(count_cols[is.na(atac_meta_ordered$sample)])
}

cat("\nATAC metadata preview:\n")
print(atac_meta_ordered)

# Ensure counts are integers
atac_counts_int <- round(atac_counts_matrix)

# Normalize using DESeq2
atac_norm <- normalize_counts_deseq2(
  atac_counts_int,
  atac_meta_ordered,
  design = ~ 1
)

# Combine peak info with normalized counts
atac_normalized_df <- cbind(
  peak_info[rownames(atac_norm$normalized), , drop = FALSE],
  atac_norm$normalized
)

# Export normalized counts
write.csv(atac_normalized_df, atac_normalized_output, row.names = FALSE)
cat("\nATAC-seq normalized counts saved to:", atac_normalized_output, "\n")
cat("  - Peaks retained after filtering:", nrow(atac_norm$normalized), "\n")
cat("  - Size factors:\n")
print(atac_norm$size_factors)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== NORMALIZATION COMPLETE ===\n")
cat("RNA-seq output:", rnaseq_normalized_output, "\n")
cat("ATAC-seq output:", atac_normalized_output, "\n")
cat("\nNormalization method: DESeq2 size factor normalization\n")
cat("This corrects for library size differences between samples.\n")
cat("\nNote: ATAC-seq duplicate rows (multiple GO terms per peak) were collapsed.\n")
