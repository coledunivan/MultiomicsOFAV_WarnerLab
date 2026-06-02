#!/usr/bin/env python3
"""
02_integration.py — link chromatin accessibility to expression.

Strategies (all written to results/tables/):
  A  same-timepoint co-significant gene-peak links     -> RNA_ATAC_same_time.csv
  B  cross-timepoint gene links + temporal ordering    -> RNA_ATAC_cross_time_scored.csv
  C  explicit ATAC(early)->RNA(later) lag links         -> RNA_ATAC_lagged.csv
  +  concordant open+induced candidates                 -> RNA_ATAC_overexpression_candidates.csv
  +  promoter-proximal heat-inducible candidates (ranked)
                                          -> S1_candidate_heat_inducible_promoters.csv

Requires step 01 outputs. Uses NOMINAL p for ATAC by default (see config).

Usage: python pipeline/02_integration.py [--config ...]
"""
import argparse
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from common import (load_config, rooted, ensure_dirs, table_path,
                    time_order, time_num, strip_prefix, save_json, banner,
                    work_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)
    ensure_dirs(cfg)
    th = cfg["thresholds"]
    TN = time_num(cfg)
    TO = time_order(cfg)

    banner("STEP 2 — RNA-ATAC integration")

    rna = pd.read_csv(table_path(cfg, "RNA_DE_by_timepoint.csv"))
    atac = pd.read_csv(table_path(cfg, "ATAC_DA_by_timepoint.csv"))

    # --- peak -> gene + feature annotation ---
    ann = pd.read_csv(rooted(cfg, cfg["inputs"]["peak_annotation"]),
                      usecols=["peak_id", "top_feature", "gene_id",
                               "chr", "start", "end"], low_memory=False)
    ann["gene_id"] = strip_prefix(ann["gene_id"], cfg)
    ann = ann[~ann.gene_id.isin(["nan", ""])].dropna(subset=["gene_id"])
    feat_rank = {"promoter": 0, "mRNA": 1, "exon": 2, "lnc_RNA": 3,
                 "transcript": 4, "region": 5}
    ann["frank"] = ann["top_feature"].map(feat_rank).fillna(9)
    peak_primary = (ann.sort_values("frank")
                    .drop_duplicates("peak_id")
                    [["peak_id", "gene_id", "top_feature", "chr", "start", "end"]])
    peak_primary.to_csv(table_path(cfg, "ATAC_peak_to_gene_annotation.csv"), index=False)
    print(f"  annotated peaks: {peak_primary.peak_id.nunique()} -> "
          f"{peak_primary.gene_id.nunique()} genes")

    # --- significant sets ---
    rna_sig = rna[(rna.padj < th["rna_padj"]) &
                  (rna.log2FoldChange.abs() >= th["rna_log2fc"])
                  ].dropna(subset=["padj"]).copy()
    if th["atac_use_padj"]:
        atac_sig = atac[(atac.padj < th["atac_padj"]) &
                        (atac.log2FoldChange.abs() >= th["atac_log2fc"])
                        ].dropna(subset=["padj"]).copy()
    else:
        atac_sig = atac[(atac.pvalue < th["atac_pval"]) &
                        (atac.log2FoldChange.abs() >= th["atac_log2fc"])
                        ].dropna(subset=["pvalue"]).copy()
    atac_sig = atac_sig.merge(peak_primary, on="peak_id", how="inner")
    print(f"  RNA sig genes: {rna_sig.gene_id.nunique()} | "
          f"ATAC sig peaks: {atac_sig.peak_id.nunique()} "
          f"(genes={atac_sig.gene_id.nunique()})")

    # ===== A: same-timepoint =====
    same = (atac_sig.merge(
        rna_sig[["gene_id", "time", "log2FoldChange", "padj"]]
        .rename(columns={"log2FoldChange": "rna_log2FC", "padj": "rna_padj"}),
        on=["gene_id", "time"], how="inner")
        .rename(columns={"log2FoldChange": "atac_log2FC",
                         "pvalue": "atac_pval", "padj": "atac_padj"}))
    same["concordant"] = np.sign(same.atac_log2FC) == np.sign(same.rna_log2FC)
    same.to_csv(table_path(cfg, "RNA_ATAC_same_time.csv"), index=False)
    print(f"  [A] same-timepoint links: {len(same)} "
          f"(concordant {int(same.concordant.sum())})")

    # ===== B: cross-timepoint gene-aggregated =====
    def strongest(x):
        return x.iloc[x.abs().values.argmax()]

    genes_de = (rna_sig.groupby("gene_id").apply(lambda x: pd.Series({
        "rna_log2FC": strongest(x.log2FoldChange),
        "rna_padj": x.padj.min(),
        "rna_earliest": min(x.time, key=lambda t: TN[t])}),
        include_groups=False).reset_index())
    peaks_da = (atac_sig.groupby(["peak_id", "gene_id", "top_feature",
                                  "chr", "start", "end"])
                .apply(lambda x: pd.Series({
                    "atac_log2FC": strongest(x.log2FoldChange),
                    "atac_pval": x.pvalue.min(),
                    "atac_earliest": min(x.time, key=lambda t: TN[t])}),
                    include_groups=False).reset_index())
    cross = peaks_da.merge(genes_de, on="gene_id", how="inner")

    def tcat(r):
        a, rr = TN[r.atac_earliest], TN[r.rna_earliest]
        return "ATAC_first" if a < rr else ("RNA_first" if rr < a else "Simultaneous")
    cross["temporal_category"] = cross.apply(tcat, axis=1)
    cross["concordant"] = np.sign(cross.atac_log2FC) == np.sign(cross.rna_log2FC)

    w = cfg["integration"]["score_weights"]
    cross["sig_score"] = (-np.log10(cross.atac_pval.clip(lower=1e-300))
                          - np.log10(cross.rna_padj.clip(lower=1e-300)))
    cross["effect_score"] = cross.atac_log2FC.abs() + cross.rna_log2FC.abs()
    cross["concordance_score"] = np.where(cross.concordant, w["concordance"], 0)
    cross["temporal_score"] = cross.temporal_category.map(
        {"ATAC_first": 3, "Simultaneous": 2, "RNA_first": 1})
    cross["promoter_score"] = np.where(cross.top_feature == "promoter", w["promoter"], 0)
    cross["candidate_score"] = (cross.sig_score + w["effect"] * cross.effect_score
                                + cross.concordance_score + cross.temporal_score
                                + cross.promoter_score)
    cross.to_csv(table_path(cfg, "RNA_ATAC_cross_time_scored.csv"), index=False)
    print(f"  [B] cross-timepoint links: {len(cross)} "
          f"(genes {cross.gene_id.nunique()}, "
          f"{100*cross.concordant.mean():.1f}% concordant)")
    print("      temporal:", cross.temporal_category.value_counts().to_dict())

    # ===== C: explicit lags =====
    lag_rows = []
    for t1, t2 in cfg["integration"]["lag_pairs"]:
        a = atac_sig[atac_sig.time == t1]
        r = (rna_sig[rna_sig.time == t2][["gene_id", "log2FoldChange", "padj"]]
             .rename(columns={"log2FoldChange": "rna_log2FC", "padj": "rna_padj"}))
        m = a.merge(r, on="gene_id", how="inner")
        if len(m):
            m = m.rename(columns={"log2FoldChange": "atac_log2FC", "pvalue": "atac_pval"})
            m["concordant"] = np.sign(m.atac_log2FC) == np.sign(m.rna_log2FC)
            m["link_type"] = f"lag_{t1}_to_{t2}"
            m["atac_time"], m["rna_time"] = t1, t2
            lag_rows.append(m)
    if lag_rows:
        lag = pd.concat(lag_rows, ignore_index=True)
        lag.to_csv(table_path(cfg, "RNA_ATAC_lagged.csv"), index=False)
        print(f"  [C] lag links: {len(lag)}")

    # ===== candidates =====
    over = cross[(cross.concordant) & (cross.rna_log2FC > 0) &
                 (cross.atac_log2FC > 0)].sort_values("candidate_score", ascending=False)
    over.to_csv(table_path(cfg, "RNA_ATAC_overexpression_candidates.csv"), index=False)

    cand = over[over.top_feature == "promoter"].copy()
    keep = ["gene_id", "peak_id", "chr", "start", "end", "top_feature",
            "rna_log2FC", "rna_padj", "atac_log2FC", "atac_pval",
            "temporal_category", "concordant", "candidate_score"]
    cand[keep].round(4).to_csv(
        table_path(cfg, "S1_candidate_heat_inducible_promoters.csv"), index=False)
    print(f"  candidates: {len(over)} open+induced; "
          f"{len(cand)} promoter-proximal (ranked -> S1)")

    summ = {
        "rna_sig_genes": int(rna_sig.gene_id.nunique()),
        "atac_sig_peaks": int(atac_sig.peak_id.nunique()),
        "atac_sig_genes": int(atac_sig.gene_id.nunique()),
        "same_time_links": int(len(same)),
        "cross_time_links": int(len(cross)),
        "cross_genes": int(cross.gene_id.nunique()),
        "concordance_rate_pct": round(100 * cross.concordant.mean(), 1),
        "atac_first": int((cross.temporal_category == "ATAC_first").sum()),
        "rna_first": int((cross.temporal_category == "RNA_first").sum()),
        "simultaneous": int((cross.temporal_category == "Simultaneous").sum()),
        "overexpression_candidates": int(len(over)),
        "promoter_candidates": int(len(cand)),
    }
    save_json(summ, work_path(cfg, "integration_summary.json"))
    banner("STEP 2 complete")


if __name__ == "__main__":
    main()
