# Methods

This document describes the analytical design in enough detail to write a
methods section and to understand the design choices (and their limitations).

## Experimental design

Paired RNA‑seq and ATAC‑seq across a thermal time course, Heat vs Control, with
biological replicates at each timepoint. The bundled example dataset is
*Orbicella faveolata* sampled at 30 min, 4 h, 12 h, and 24 h with two replicates
per group (16 RNA + 16 ATAC libraries).

## Count matrices and QC

RNA reads were summarized to gene models (featureCounts); ATAC reads to a
consensus peak set. Because the annotated ATAC file repeats each peak once per
overlapping GO term, peaks are collapsed to one count row per `peak_id` and the
most specific genomic feature (promoter > mRNA > exon > lncRNA > region) is kept
as the peak's primary feature for integration.

Library‑size normalization uses DESeq2 median‑of‑ratios size factors. Sample
relationships are assessed by PCA on log‑transformed, size‑factor‑normalized
counts (top 500 most‑variable features). Samples listed in
`pca_exclude_samples` are dropped from the ordination only; samples in
`drop_samples` are removed entirely.

## Differential analysis (per timepoint)

To respect strong temporal structure and modest replication, testing is done
**within each timepoint** using `pydeseq2` with design `~condition` and the
contrast Heat vs Control:

```
for tp in timepoints:
    DESeqDataSet(counts[samples_at_tp], ~condition).deseq2()
    DeseqStats(contrast = ["condition", "Heat", "Control"])
```

- **DEGs:** `padj < rna_padj` and `|log2FC| ≥ rna_log2fc` (defaults 0.05, 1.0).
- **DA peaks:** by default `pvalue < atac_pval` and `|log2FC| ≥ atac_log2fc`
  (defaults 0.05, 0.3), with `direction ∈ {up, down, ns}`.

A `~condition * time` interaction model is a reasonable alternative when
replication is higher; the per‑timepoint contrasts were chosen to match the
lab's original DESeq2 scripts and to remain interpretable with n = 2/group.

### Important limitation — ATAC statistical power

With ~170k consensus peaks and two replicates per group, Benjamini–Hochberg
correction across all peaks typically yields **no peak below conventional
adjusted‑p thresholds** (in the example data, minimum padj ≈ 0.9). The pipeline
therefore defaults to a **nominal‑p exploratory threshold** for accessibility
(`atac_use_padj: false`). This is an explicit, documented trade‑off:

- Accessibility‑derived results (including the candidate promoter list) are
  **hypotheses for validation**, not FDR‑controlled discoveries.
- With additional replicates, set `atac_use_padj: true` to switch to proper FDR
  control; everything downstream adapts automatically.

The diffuse ATAC PCA (low leading‑component variance, weak treatment separation)
is consistent with this being a power/effect‑size issue rather than absence of
signal.

## Temporal clustering

DEGs significant in ≥ `min_timepoints_significant` timepoints (default 2) are
assembled into gene × time log₂FC trajectories, z‑scored per gene, and
partitioned by **fuzzy c‑means** (`scikit‑fuzzy`). The fuzzifier *m* is set from
the dimensionality‑aware heuristic of Schwämmle & Jensen (2010). The number of
clusters *k* is chosen by minimizing the **Xie–Beni** validity index over
`k_search` (override with `k_force`). Genes with maximum membership ≥
`core_membership` (default 0.6) are "core" members and are used for
cluster‑level summaries and enrichment.

## RNA–ATAC integration

DA peaks are mapped to genes via the peak→gene annotation; each gene takes its
strongest‑effect linked peak. Three linking strategies:

- **A — same timepoint:** peak and gene both significant at the same time.
- **B — cross timepoint (temporally ordered):** gene‑level links annotated with
  the ordering of the *earliest* accessibility vs. expression change
  (`ATAC_first`, `Simultaneous`, `RNA_first`).
- **C — explicit lags:** accessibility at an earlier time tested against
  expression at a later time, over `integration.lag_pairs`.

**Concordance** = agreement in the sign of the ATAC and RNA log₂ fold changes.

**Candidate score** (for cross‑timepoint links), with weights from the config:

```
score = sig_score
      + effect_weight      * (|atac_log2FC| + |rna_log2FC|)
      + concordance_weight * 1[concordant]
      + temporal_score                # ATAC_first=3, Simultaneous=2, RNA_first=1
      + promoter_weight    * 1[peak is promoter]
sig_score = -log10(atac_p) - log10(rna_padj)
```

**Candidate heat‑inducible promoters** = concordant links that are promoter‑
proximal, open under heat, and induced — ranked by score in
`S1_candidate_heat_inducible_promoters.csv`.

## GO enrichment

One‑sided hypergeometric tests against the genome‑wide annotation universe
(all genes with ≥ 1 GO term), with Benjamini–Hochberg FDR. Terms with ≥ 2 query
genes and term size in [3, 2000] are tested. Reported for each fuzzy cluster
(core genes), the candidate sets, and all heat‑induced/repressed genes.

GO term **names** come from the `go_term` column of the gene→GO file when it is
informative; otherwise from `config/go_names.tsv`, a curated lookup of canonical
names for common terms. Unmapped terms keep their GO id — **extend
`go_names.tsv`** to label more terms in the figures.

## Software

- `pydeseq2` — differential expression/accessibility
- `scikit‑fuzzy` — fuzzy c‑means
- `scipy` / `statsmodels` — hypergeometric tests + FDR
- `matplotlib` — figures
- `pandas`, `numpy`, `pyyaml`

## Reference

Schwämmle V, Jensen ON (2010). A simple and fast method to determine the
parameters for fuzzy c‑means cluster analysis. *Bioinformatics* 26(22):2841–2848.
