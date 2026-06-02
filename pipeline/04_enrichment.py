#!/usr/bin/env python3
"""
04_enrichment.py — hypergeometric GO enrichment for fuzzy clusters and for
the integration candidate sets, with Benjamini-Hochberg FDR.

GO term *names* are taken from the go_term column of the gene2go file when it
differs from the go_id; otherwise from the bundled go_names.tsv lookup (curated
canonical names for common terms). Unmapped terms keep their GO id.

Outputs: results/tables/GO_enrichment_all.csv (+ significant subset).

Usage: python pipeline/04_enrichment.py [--config ...]
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from common import (load_config, rooted, ensure_dirs, table_path,
                    strip_prefix, banner)


def load_go_names(cfg):
    here = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(here, "..", "config", "go_names.tsv")
    if os.path.exists(p):
        d = pd.read_csv(p, sep="\t")
        return dict(zip(d.go_id, d.go_name))
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)
    ensure_dirs(cfg)
    th = cfg["thresholds"]

    banner("STEP 4 — GO enrichment")
    g2g = pd.read_csv(rooted(cfg, cfg["inputs"]["gene2go"]))
    g2g["gene_id"] = strip_prefix(g2g["gene_id"], cfg)
    g2g = g2g.dropna(subset=["gene_id", "go_id"])
    aspect_name = {"P": "Biological Process", "F": "Molecular Function",
                   "C": "Cellular Component"}
    g2g["aspect_full"] = g2g.aspect.map(aspect_name).fillna(g2g.aspect)

    universe = set(g2g.gene_id.unique())
    go2genes = g2g.groupby("go_id")["gene_id"].apply(set).to_dict()
    go_aspect = g2g.drop_duplicates("go_id").set_index("go_id")["aspect_full"].to_dict()
    # prefer dataset-provided names where informative, else curated lookup
    name_lookup = load_go_names(cfg)
    go_name = {}
    for go in go2genes:
        gt = g2g.loc[g2g.go_id == go, "go_term"].iloc[0]
        go_name[go] = name_lookup.get(go, gt if gt != go else go)
    print(f"  GO universe: {len(universe)} genes, {len(go2genes)} terms")

    def enrich(gene_set, label, min_term=3, max_term=2000):
        gene_set = set(gene_set) & universe
        M, n = len(universe), len(gene_set)
        rows = []
        for go, genes in go2genes.items():
            K = len(genes)
            if K < min_term or K > max_term:
                continue
            x = len(genes & gene_set)
            if x < 2:
                continue
            p = hypergeom.sf(x - 1, M, K, n)
            rows.append((go, go_name.get(go, go), go_aspect.get(go, "?"),
                         x, K, n, p))
        if not rows:
            return pd.DataFrame()
        df = pd.DataFrame(rows, columns=["go_id", "go_name", "aspect",
                                         "n_in_set", "n_term", "set_size", "pvalue"])
        df["padj"] = multipletests(df.pvalue, method="fdr_bh")[1]
        df["fold_enrich"] = (df.n_in_set / df.set_size) / (df.n_term / len(universe))
        df["cluster_or_set"] = label
        return df.sort_values("pvalue")

    all_enr = []

    # fuzzy clusters (core genes)
    clust = pd.read_csv(table_path(cfg, "fuzzy_clusters.csv"))
    for c in sorted(clust.cluster.unique()):
        gs = clust[(clust.cluster == c) & clust.core].gene_id
        e = enrich(gs, f"Cluster_{c}")
        if len(e):
            print(f"  Cluster {c}: {(e.padj<0.1).sum()} terms FDR<0.1")
            all_enr.append(e)

    # candidate sets
    over_path = table_path(cfg, "RNA_ATAC_overexpression_candidates.csv")
    if os.path.exists(over_path):
        over = pd.read_csv(over_path)
        e = enrich(over.gene_id.unique(), "Overexpression_candidates")
        if len(e):
            all_enr.append(e)

    # all heat-induced / repressed
    rna = pd.read_csv(table_path(cfg, "RNA_DE_by_timepoint.csv"))
    rna["gene_id"] = strip_prefix(rna["gene_id"], cfg)
    up = rna[(rna.padj < th["rna_padj"]) &
             (rna.log2FoldChange >= th["rna_log2fc"])].gene_id.unique()
    dn = rna[(rna.padj < th["rna_padj"]) &
             (rna.log2FoldChange <= -th["rna_log2fc"])].gene_id.unique()
    for gs, lab in [(up, "Heat_induced_all"), (dn, "Heat_repressed_all")]:
        e = enrich(gs, lab)
        if len(e):
            all_enr.append(e)

    enr = pd.concat(all_enr, ignore_index=True)
    enr.round(5).to_csv(table_path(cfg, "GO_enrichment_all.csv"), index=False)
    sig = enr[enr.padj < 0.1].sort_values(["cluster_or_set", "padj"])
    sig.round(5).to_csv(table_path(cfg, "GO_enrichment_significant.csv"), index=False)
    print(f"  wrote {len(enr)} enrichment rows ({len(sig)} significant)")
    banner("STEP 4 complete")


if __name__ == "__main__":
    main()
