# Output files

All tables are written to `results/tables/`, figures to `results/figures/`.
Intermediate artifacts shared between steps go to `results/_work/` (safe to
delete after a run; regenerated as needed).

## Tables

| File | Step | Contents |
|------|------|----------|
| `RNA_DE_by_timepoint.csv` | 1 | Per‑gene, per‑timepoint DESeq2 results (baseMean, log2FoldChange, pvalue, padj, time, direction) |
| `ATAC_DA_by_timepoint.csv` | 1 | Per‑peak, per‑timepoint DA results (all peaks) |
| `ATAC_DA_significant.csv` | 1 | DA peaks passing thresholds only (compact) |
| `ATAC_peak_to_gene_annotation.csv` | 2 | One row per peak: primary feature + linked gene |
| `RNA_ATAC_same_time.csv` | 2 | Strategy A: same‑timepoint co‑significant links |
| `RNA_ATAC_cross_time_scored.csv` | 2 | Strategy B: cross‑timepoint links, temporal category, candidate score |
| `RNA_ATAC_lagged.csv` | 2 | Strategy C: explicit ATAC(early)→RNA(later) links |
| `RNA_ATAC_overexpression_candidates.csv` | 2 | Concordant open+induced links (all features) |
| `S1_candidate_heat_inducible_promoters.csv` | 2 | **Ranked promoter‑proximal heat‑inducible candidates** |
| `fuzzy_clusters.csv` | 3 | Gene → cluster, membership, core flag |
| `fuzzy_centroids_zscore.csv` | 3 | Cluster centroids (z‑scored) |
| `fuzzy_centroids_raw_log2fc.csv` | 3 | Cluster centroids (raw log₂FC, core genes) |
| `fuzzy_xie_beni.csv` | 3 | Validity index per candidate k |
| `GO_enrichment_all.csv` | 4 | All tested GO terms (go_id, go_name, aspect, counts, pvalue, padj, fold_enrich, set) |
| `GO_enrichment_significant.csv` | 4 | FDR < 0.1 subset |

## Figures

| File | Description |
|------|-------------|
| `Figure1_design_QC.png` | Experimental design; RNA & ATAC PCA; DE/DA counts per timepoint |
| `Figure2_fuzzy_clusters.png` | Temporal program trajectories, heatmap, Xie–Beni curve, per‑cluster GO |
| `Figure3_integration.png` | Concordance scatter, temporal ordering, lag links, promoter‑score, summary |
| `Figure4_candidates.png` | Ranked candidate promoters, regulatory timing, candidate landscape, trajectories |
| `Figure5_synthesis.png` | Cascade schematic + global response‑magnitude curve |

## Key columns in `S1_candidate_heat_inducible_promoters.csv`

| Column | Meaning |
|--------|---------|
| `gene_id` | linked gene |
| `peak_id`, `chr`, `start`, `end` | the promoter‑proximal ATAC peak |
| `rna_log2FC`, `rna_padj` | strongest expression effect + best padj |
| `atac_log2FC`, `atac_pval` | strongest accessibility effect + best p |
| `temporal_category` | ATAC_first / Simultaneous / RNA_first |
| `concordant` | sign agreement between layers (always True here) |
| `candidate_score` | composite ranking score (higher = stronger candidate) |
