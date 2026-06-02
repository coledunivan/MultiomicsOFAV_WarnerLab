# Coral Multiomics Pipeline

A reproducible RNA‑seq + ATAC‑seq pipeline for **time‑course thermal‑stress
experiments**, developed on *Orbicella faveolata* heat‑stress data. Given paired
expression and chromatin‑accessibility counts across a treatment time course, it
produces per‑timepoint differential results, integrated chromatin↔expression
links, fuzzy temporal clusters, GO enrichment, and five publication‑quality
figures — all driven by a single config file.

> **Reusing this on new data?** In most cases you only edit
> [`config/config.yaml`](config/config.yaml): point it at your count files and
> describe how your sample names encode timepoint/treatment. See
> [Adapting to new data](#adapting-to-new-data).

---

## What it does

| Step | Script | Output |
|------|--------|--------|
| 1 | `pipeline/01_differential.py` | Per‑timepoint **DE** (RNA) and **DA** (ATAC), Heat vs Control |
| 2 | `pipeline/02_integration.py` | Same‑time, cross‑time (temporally ordered), and lagged **RNA↔ATAC links**; ranked candidate heat‑inducible promoters |
| 3 | `pipeline/03_clustering.py` | **Fuzzy c‑means** temporal programs, with Xie–Beni model selection |
| 4 | `pipeline/04_enrichment.py` | **GO enrichment** (hypergeometric + BH‑FDR) per cluster and candidate set |
| 5 | `pipeline/05_figures.py` | **Figures 1–5** (design/QC, temporal programs, integration, candidates, synthesis) |

The analytical design mirrors the lab's original per‑timepoint DESeq2 contrasts
(`~condition`, Heat vs Control, within each timepoint), reimplemented in Python
with `pydeseq2` so the whole pipeline runs without R. The original R scripts are
preserved in [`scripts_R_original/`](scripts_R_original/) for reference.

---

## Quick start

```bash
# 1. install dependencies (Python 3.10+)
pip install -r requirements.txt

# 2. put your count files in data/ and edit config/config.yaml to match
#    (or use the bundled example data as-is)

# 3. run everything
python run_pipeline.py

# RNA-only quick pass (skips the slow ATAC DESeq2):
python run_pipeline.py --skip-atac

# resume from a step, or run just one:
python run_pipeline.py --from 3
python run_pipeline.py --only 5
```

Outputs land in `results/tables/` (CSV) and `results/figures/` (PNG).

> **Heads up — ATAC runtime.** ATAC DESeq2 over a full consensus peak set
> (~170k peaks, 2 reps/group) takes roughly 15–25 min on one core and needs
> ~2 GB RAM. Use `--skip-atac` while iterating on the RNA side; steps 2 and 5
> will fall back gracefully if ATAC outputs are absent for the parts that need
> them, but full integration requires step 1 to have run with ATAC.

---

## Inputs

Place these in `data/` (paths are set in the config):

| File | Description |
|------|-------------|
| RNA counts | gene × sample matrix; first column = gene id |
| ATAC counts | peak matrix with `chr,start,end,peak_id` then one count column per sample |
| Peak annotation | maps `peak_id` → `gene_id` + `top_feature` (promoter/mRNA/exon/…) |
| gene→GO | `gene_id, go_id, go_term, aspect` (aspect = P/F/C) |

Sample columns in the RNA and ATAC matrices must follow a consistent naming
scheme that encodes timepoint and treatment, e.g. `OFav_05_C1_S1`
(`05` → 30 min, `C` → Control). The decoding rules live in the config.

---

## Adapting to new data

Edit [`config/config.yaml`](config/config.yaml):

1. **`inputs:`** — point at your four files.
2. **`sample_naming:`** — describe how a column header decodes into
   `(timepoint, treatment)`. You set the delimiter, which split‑token holds the
   time and which holds the condition, the Heat/Control prefixes, and the
   timepoint label/ordering maps. Nothing else in the code is hardcoded to the
   coral naming scheme.
3. **`thresholds:`** — DE/DA significance cutoffs. **Note on ATAC:** with low
   replication, genome‑wide FDR over ~170k peaks typically leaves nothing
   significant, so the default uses a *nominal* p threshold for accessibility
   (`atac_use_padj: false`). If you have more replicates, set
   `atac_use_padj: true` for proper FDR control. This limitation is discussed in
   [`docs/METHODS.md`](docs/METHODS.md).
4. **`drop_samples` / `pca_exclude_samples`** — exclude confirmed outliers.
5. **`clustering` / `integration`** — optional: cluster‑search range, lag pairs,
   candidate‑score weights.

If your timepoints differ (say 0/6/24/48 h), just change `time_map`,
`time_numeric`, `time_order`, and the `integration.lag_pairs`. The figures and
tables follow automatically.

---

## Example results (bundled *O. faveolata* data)

- **Differential expression** (padj < 0.05, |log₂FC| ≥ 1): 877↑/1501↓ at 30 min,
  684↑/1744↓ at 4 h, 1079↑/1740↓ at 12 h, 168↑/270↓ at 24 h — a response that
  crests at 12 h and resolves by 24 h.
- **Temporal programs:** Xie–Beni selects **k = 4** — early‑transient,
  delayed‑induction (DNA repair / cell cycle), biphasic (apoptotic signaling),
  and mid‑sustained (chaperone proteostasis: HSP, CCT/TRiC).
- **Integration:** 1,250 cross‑timepoint gene‑level links, ~52% sign‑concordant,
  with ~26% showing accessibility changes preceding expression; **37
  promoter‑proximal candidate heat‑inducible loci** ranked in
  `S1_candidate_heat_inducible_promoters.csv`.

(Exact numbers depend on thresholds in the config; the values above use the
repository defaults.)

---

## A note on large files

The example count matrices in `data/` are large (the ATAC peak file is ~16 MB;
the annotation file ~64 MB). A `.gitattributes` is included to track `data/*.csv`
with [Git LFS](https://git-lfs.com) — install LFS before your first push, or host
the data separately (e.g. Zenodo/figshare) and keep only small files in git by
uncommenting the `data/*.csv` line in `.gitignore`. Large *regenerable* output
tables are not shipped; they are produced when you run the pipeline.

## Repository layout

```
.
├── config/
│   ├── config.yaml          # ← the file you edit
│   └── go_names.tsv         # curated GO id → name lookup (extensible)
├── pipeline/
│   ├── common.py            # config loading + sample-name decoding
│   ├── 01_differential.py
│   ├── 02_integration.py
│   ├── 03_clustering.py
│   ├── 04_enrichment.py
│   └── 05_figures.py
├── run_pipeline.py          # orchestration
├── scripts_R_original/      # the lab's original R scripts (reference)
├── docs/
│   ├── METHODS.md           # detailed methods + design rationale + caveats
│   └── OUTPUTS.md           # description of every output file
├── data/                    # input counts (example data included)
├── results/
│   ├── tables/              # CSV outputs
│   └── figures/             # PNG figures
├── paper/                   # example manuscript built from these results
├── requirements.txt
└── README.md
```

---

## Citing / provenance

This pipeline reimplements and packages an analysis originally developed in R by
the lab. If you publish results, please describe the methods as in
[`docs/METHODS.md`](docs/METHODS.md) and cite `pydeseq2`, `scikit‑fuzzy`,
`statsmodels`, and the original DESeq2 methodology.

## License

MIT (see [`LICENSE`](LICENSE)). The bundled example data remains the property of
the originating lab; check with them before redistributing.
