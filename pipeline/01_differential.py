#!/usr/bin/env python3
"""
01_differential.py — per-timepoint differential expression (RNA) and
accessibility (ATAC), Heat vs Control.

Reproduces the lab's per-timepoint DESeq2 design in pydeseq2:
    within each timepoint:  counts ~ condition   (contrast Heat vs Control)

Outputs (to results/tables/):
    RNA_DE_by_timepoint.csv
    ATAC_DA_by_timepoint.csv          (all peaks; 'direction' column applied)
    ATAC_DA_significant.csv           (compact, significant peaks only)
Intermediates (results/_work/): count matrices + metadata for later steps.

Usage:
    python pipeline/01_differential.py [--config path/to/config.yaml]
                                       [--skip-atac]   # RNA only (fast)
"""
import argparse
import gc
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from common import (load_config, rooted, ensure_dirs, work_path, table_path,
                    sample_metadata, time_order, strip_prefix,
                    drop_configured_samples, save_json, banner)


def run_deseq_per_timepoint(counts, meta, feature_name, cfg,
                            padj_t, lfc_t, use_padj, prefilter=None):
    """counts: features x samples. meta: samples x [condition,time]."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    th = cfg["thresholds"]
    contrast = cfg["design"]["contrast"]
    results = []
    for tp in time_order(cfg):
        s = meta.index[meta["time"] == tp].tolist()
        sub_meta = meta.loc[s, ["condition"]].copy()
        if sub_meta["condition"].nunique() < 2:
            print(f"  [{feature_name}] skip {tp}: only one condition present")
            continue
        mat = counts[s]
        if prefilter == "sum":
            keep = mat.sum(axis=1) >= th["min_count_sum"]
        elif prefilter == "atac":
            keep = (mat >= th["atac_min_count"]).sum(axis=1) >= th["atac_min_samples"]
        else:
            keep = mat.sum(axis=1) >= th["min_count_sum"]
        mat = mat.loc[keep]
        sub_meta["condition"] = pd.Categorical(sub_meta["condition"],
                                               categories=["Control", "Heat"])
        dds = DeseqDataSet(counts=mat.T.astype(int), metadata=sub_meta,
                           design="~condition", quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=contrast, quiet=True)
        st.summary()
        r = st.results_df.copy()
        r["time"] = tp
        r.index.name = feature_name
        r = r.reset_index()
        results.append(r)
        score = r["padj"] if use_padj else r["pvalue"]
        thr = padj_t if use_padj else cfg["thresholds"]["atac_pval"]
        n_sig = ((score < thr) & (r["log2FoldChange"].abs() >= lfc_t)).sum()
        print(f"  [{feature_name}] {tp}: {mat.shape[0]:>7} tested, {int(n_sig):>5} significant")
        del dds, st, mat
        gc.collect()
    out = pd.concat(results, ignore_index=True)

    def classify(row):
        score = row["padj"] if use_padj else row["pvalue"]
        if pd.isna(score) or pd.isna(row["log2FoldChange"]):
            return "ns"
        if score < (padj_t if use_padj else cfg["thresholds"]["atac_pval"]):
            if row["log2FoldChange"] >= lfc_t:
                return "up"
            if row["log2FoldChange"] <= -lfc_t:
                return "down"
        return "ns"

    out["direction"] = out.apply(classify, axis=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-atac", action="store_true")
    args = ap.parse_args()
    cfg = load_config(args.config)
    ensure_dirs(cfg)
    th = cfg["thresholds"]

    # ----- RNA -----
    banner("STEP 1a — RNA differential expression (Heat vs Control)")
    rna = pd.read_csv(rooted(cfg, cfg["inputs"]["rna_counts"]),
                      index_col=cfg["inputs"]["rna_gene_id_col"])
    rna.columns = [c.strip() for c in rna.columns]
    rna.index = strip_prefix(pd.Series(rna.index), cfg).values
    keep_cols = drop_configured_samples(rna.columns, cfg)
    rna = rna[keep_cols]
    rna_meta = sample_metadata(rna.columns, cfg)
    print(f"  RNA matrix: {rna.shape[0]} genes x {rna.shape[1]} samples")

    rna_de = run_deseq_per_timepoint(rna, rna_meta, "gene_id", cfg,
                                     th["rna_padj"], th["rna_log2fc"],
                                     use_padj=True, prefilter="sum")
    rna_de.to_csv(table_path(cfg, "RNA_DE_by_timepoint.csv"), index=False)
    rna.to_pickle(work_path(cfg, "rna_counts.pkl"))
    rna_meta.to_csv(work_path(cfg, "rna_meta.csv"))

    # ----- ATAC -----
    if not args.skip_atac:
        banner("STEP 1b — ATAC differential accessibility (Heat vs Control)")
        atac = pd.read_csv(rooted(cfg, cfg["inputs"]["atac_counts"]))
        id_cols = cfg["inputs"]["atac_id_cols"]
        sample_cols = [c for c in atac.columns if c not in id_cols]
        sample_cols = drop_configured_samples(sample_cols, cfg)
        peak_meta = atac[id_cols].copy()
        atac_counts = atac.set_index("peak_id")[sample_cols]
        atac_meta = sample_metadata(sample_cols, cfg)
        print(f"  ATAC matrix: {atac_counts.shape[0]} peaks x {atac_counts.shape[1]} samples")

        atac_da = run_deseq_per_timepoint(
            atac_counts, atac_meta, "peak_id", cfg,
            th["atac_padj"], th["atac_log2fc"],
            use_padj=th["atac_use_padj"], prefilter="atac")
        atac_da.to_csv(table_path(cfg, "ATAC_DA_by_timepoint.csv"), index=False)
        atac_da[atac_da.direction != "ns"].to_csv(
            table_path(cfg, "ATAC_DA_significant.csv"), index=False)
        peak_meta.to_pickle(work_path(cfg, "peak_meta.pkl"))
        atac_meta.to_csv(work_path(cfg, "atac_meta.csv"))

    # ----- summary -----
    summ = {"rna_genes": int(rna.shape[0]), "rna_samples": int(rna.shape[1])}
    for tp in time_order(cfg):
        d = rna_de[rna_de.time == tp]
        summ[f"rna_up_{tp}"] = int((d.direction == "up").sum())
        summ[f"rna_down_{tp}"] = int((d.direction == "down").sum())
    save_json(summ, work_path(cfg, "step1_summary.json"))
    banner("STEP 1 complete")
    print(json._dump_str(summ) if hasattr(json, "_dump_str") else summ)


if __name__ == "__main__":
    import json
    main()
