# Input data

The example *Orbicella faveolata* count files live here. To run on your own
experiment, replace these (or add yours and update `config/config.yaml`):

- `countsRNAseq_with_C1.csv` — RNA gene × sample counts
- `ATAC_peak_counts_clean.csv` — ATAC peak counts (chr,start,end,peak_id + samples)
- `COUNTS_ATAC_or_peak_counts.csv` — peak → gene/feature annotation
- `gene2go_from_ATAC_counts.csv` — gene → GO annotation

Sample column names must encode timepoint and treatment per the `sample_naming`
block in the config (the bundled files use e.g. `OFav_05_C1_S1`).
