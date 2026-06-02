#!/usr/bin/env python3
"""
03_clustering.py — fuzzy c-means clustering of DE temporal trajectories.

Selects k by minimizing the Xie-Beni index over config.clustering.k_search
(unless k_force is set). Writes cluster assignments, centroids, and the
validity curve to results/tables/, plus a z-scored matrix for plotting.

Usage: python pipeline/03_clustering.py [--config ...]
"""
import argparse
import numpy as np
import pandas as pd
import skfuzzy as fuzz

from common import (load_config, ensure_dirs, table_path, work_path,
                    time_order, strip_prefix, save_json, banner)


def xie_beni(data, cntr, u, m):
    n = data.shape[0]
    c = cntr.shape[0]
    d2 = np.zeros((n, c))
    for j in range(c):
        d2[:, j] = ((data - cntr[j]) ** 2).sum(1)
    num = ((u.T ** m) * d2).sum()
    cd = np.full((c, c), np.inf)
    for i in range(c):
        for j in range(c):
            if i != j:
                cd[i, j] = ((cntr[i] - cntr[j]) ** 2).sum()
    return num / (n * cd.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)
    ensure_dirs(cfg)
    cc = cfg["clustering"]
    TO = time_order(cfg)

    banner("STEP 3 — fuzzy c-means temporal clustering")
    rna = pd.read_csv(table_path(cfg, "RNA_DE_by_timepoint.csv"))
    rna["gene_id"] = strip_prefix(rna["gene_id"], cfg)

    th = cfg["thresholds"]
    sig = rna[(rna.padj < th["rna_padj"]) &
              (rna.log2FoldChange.abs() >= th["rna_log2fc"])
              ].dropna(subset=["padj"])
    multi = sig.groupby("gene_id").size()
    multi = multi[multi >= cc["min_timepoints_significant"]].index
    sig = sig[sig.gene_id.isin(multi)]
    print(f"  genes significant in >= {cc['min_timepoints_significant']} "
          f"timepoints: {sig.gene_id.nunique()}")

    mat = (sig.pivot_table(index="gene_id", columns="time",
                           values="log2FoldChange").reindex(columns=TO))
    mat = mat[mat.notna().sum(axis=1) >= cc["min_timepoints_significant"]]
    raw = mat.copy()
    mat = mat.fillna(0.0)

    X = mat.values
    Xs = (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-9)
    keep = np.isfinite(Xs).all(1)
    Xs = Xs[keep]
    genes = mat.index[keep].tolist()
    raw = raw.loc[genes]
    print(f"  clustering matrix: {Xs.shape}")

    # dimensionality-aware fuzzifier (Schwammle & Jensen)
    Nn, D = Xs.shape
    m = 1 + (1418 / Nn + 22.05) * D ** -2 + \
        (12.33 / Nn + 0.243) * D ** (-0.0406 * np.log(Nn) - 0.1134)
    m = float(np.clip(m, 1.1, 2.5))
    print(f"  fuzzifier m = {m:.3f}")

    # Xie-Beni over candidate k
    xb = []
    for k in cc["k_search"]:
        cntr, u, *_ = fuzz.cluster.cmeans(Xs.T, k, m, error=1e-6,
                                          maxiter=1000, seed=cc["seed"])
        xb.append((k, xie_beni(Xs, cntr, u, m)))
        print(f"    k={k}: XB={xb[-1][1]:.4f}")
    xb_df = pd.DataFrame(xb, columns=["k", "xie_beni"])
    xb_df.to_csv(table_path(cfg, "fuzzy_xie_beni.csv"), index=False)

    K = cc["k_force"] or int(xb_df.loc[xb_df.xie_beni.idxmin(), "k"])
    print(f"  selected K = {K}")

    cntr, u, *_ = fuzz.cluster.cmeans(Xs.T, K, m, error=1e-7,
                                      maxiter=2000, seed=cc["seed"])
    cluster = np.argmax(u, axis=0) + 1
    membership = np.max(u, axis=0)
    clust = pd.DataFrame({"gene_id": genes, "cluster": cluster,
                          "membership": membership})
    clust["core"] = clust.membership >= cc["core_membership"]
    clust.round(4).to_csv(table_path(cfg, "fuzzy_clusters.csv"), index=False)

    cent_z = pd.DataFrame(cntr, columns=TO)
    cent_z.index = [f"C{i+1}" for i in range(K)]
    cent_z.to_csv(table_path(cfg, "fuzzy_centroids_zscore.csv"))

    raw2 = raw.copy()
    cmap = clust.set_index("gene_id")
    raw2["cluster"] = cmap.loc[raw.index, "cluster"].values
    raw2["core"] = cmap.loc[raw.index, "core"].values
    (raw2[raw2.core].groupby("cluster")[TO].mean()
     ).to_csv(table_path(cfg, "fuzzy_centroids_raw_log2fc.csv"))

    np.save(work_path(cfg, "fuzzy_Xs.npy"), Xs)
    pd.Series(genes, name="gene_id").to_csv(work_path(cfg, "fuzzy_genes.csv"),
                                            index=False)
    save_json({"n_genes": int(Nn), "m": m, "K": int(K),
               "core_genes": int(clust.core.sum())},
              work_path(cfg, "fuzzy_summary.json"))
    print("  cluster sizes (all / core):")
    for c in range(1, K + 1):
        print(f"    C{c}: {(clust.cluster==c).sum()} / "
              f"{((clust.cluster==c)&clust.core).sum()}")
    banner("STEP 3 complete")


if __name__ == "__main__":
    main()
