# Original R scripts (reference)

These are the lab's original R scripts that the Python pipeline reimplements and
packages. They are kept for provenance and cross‑checking. They are **not** part
of the runnable pipeline (which is pure Python and config‑driven). Mapping:

| Original R script(s) | Python equivalent |
|----------------------|-------------------|
| `rna_deseq2_timepoint_DE.R`, `atac_deseq2_timepoint_DA.R` | `pipeline/01_differential.py` |
| `build_clean_atac_counts_and_annotations.R`, `normalize_counts.R` | preprocessing folded into step 1 / config |
| `pca_and_dispersion.R`, `pca_sample_qc_standardized.R` | PCA inside `pipeline/05_figures.py` |
| `integrate_rna_atac_same_time.R`, `RNA_ATAC_integration_v2.R`, `integrate_rna_atac_lagged.R`, `go_to_peak_mapping.R` | `pipeline/02_integration.py` |
| `FuzzyClusteringFullWorkflow.R` | `pipeline/03_clustering.py` |
| `rank_candidate_heat_inducible_promoters.R` | candidate ranking in `pipeline/02_integration.py` |
| `volcano_with_deg_locus.R`, `integration_summary_figures.R` | `pipeline/05_figures.py` |
