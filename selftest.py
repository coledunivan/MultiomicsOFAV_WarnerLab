#!/usr/bin/env python3
"""
selftest.py — verify the install works without the full dataset.

Generates a tiny synthetic RNA + ATAC dataset (same naming scheme as the
example config), writes it to a temp folder with a derived config, and runs
steps 1-3 RNA-side to confirm the toolchain is wired up. Not a scientific test —
just a smoke test of the environment.

    python selftest.py
"""
import os, tempfile, subprocess, sys, textwrap
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(0)

def main():
    tmp = tempfile.mkdtemp(prefix="coralmx_selftest_")
    data = os.path.join(tmp, "data"); os.makedirs(data)
    samples = [f"OFav_{t}_{c}{r}_S{i}" for i,(t,c,r) in enumerate(
        [(t,c,r) for t in ["05","4","12","24"] for c in ["C","H"] for r in ["1","2"]],1)]
    # RNA: 400 genes; inject a heat effect in 60 of them
    genes=[f"LOC{900000+i}" for i in range(400)]
    base=rng.poisson(80,(400,len(samples)))
    for j,s in enumerate(samples):
        if "_H" in s: base[:60,j]=rng.poisson(300,60)
    rna=pd.DataFrame(base,index=genes,columns=samples); rna.index.name="gene_id"
    rna.to_csv(f"{data}/rna.csv")
    # ATAC: 500 peaks
    peaks=pd.DataFrame({"chr":["chr1"]*500,
        "start":np.arange(500)*1000,"end":np.arange(500)*1000+200,
        "peak_id":[f"peak_{i}" for i in range(500)]})
    for s in samples: peaks[s]=rng.poisson(40,500)
    peaks.to_csv(f"{data}/atac.csv",index=False)
    # annotation + gene2go
    ann=peaks[["peak_id"]].copy(); ann["top_feature"]="promoter"
    ann["gene_id"]=[genes[i%400] for i in range(500)]
    ann["chr"]="chr1"; ann["start"]=peaks.start; ann["end"]=peaks.end
    ann.to_csv(f"{data}/ann.csv",index=False)
    g2g=pd.DataFrame({"gene_id":genes[:200],
        "go_id":["GO:0006457"]*100+["GO:0009408"]*100,
        "go_term":["GO:0006457"]*100+["GO:0009408"]*100,
        "aspect":["P"]*200})
    g2g.to_csv(f"{data}/g2g.csv",index=False)

    cfg=textwrap.dedent(f"""
    project: {{name: selftest, species_label: "Test species", species_common: ""}}
    inputs:
      rna_counts: data/rna.csv
      atac_counts: data/atac.csv
      peak_annotation: data/ann.csv
      gene2go: data/g2g.csv
      rna_gene_id_col: 0
      atac_id_cols: [chr, start, end, peak_id]
      gene_id_prefix_strip: ""
    sample_naming:
      delimiter: "_"
      time_token_index: 1
      cond_token_index: 2
      heat_prefix: H
      control_prefix: C
      time_map: {{"05": 30min, "4": 4h, "12": 12h, "24": 24h}}
      time_numeric: {{30min: 0.5, 4h: 4.0, 12h: 12.0, 24h: 24.0}}
      time_order: [30min, 4h, 12h, 24h]
    drop_samples: []
    pca_exclude_samples: []
    thresholds: {{rna_padj: 0.1, rna_log2fc: 0.5, atac_use_padj: false, atac_pval: 0.05,
      atac_padj: 0.1, atac_log2fc: 0.3, min_count_sum: 10, atac_min_count: 5, atac_min_samples: 2}}
    design: {{per_timepoint: true, contrast: [condition, Heat, Control]}}
    clustering: {{min_timepoints_significant: 1, k_search: [2,3,4], k_force: null,
      core_membership: 0.5, seed: 42}}
    integration: {{lag_pairs: [[30min,4h]], score_weights: {{effect: 2.0, concordance: 5.0, promoter: 3.0}}}}
    output: {{dir: results, tables_dir: results/tables, figures_dir: results/figures, figure_dpi: 100}}
    """)
    cfgp=os.path.join(tmp,"config","config.yaml")
    os.makedirs(os.path.dirname(cfgp),exist_ok=True)
    open(cfgp,"w").write(cfg)
    os.makedirs(os.path.join(tmp,"results"),exist_ok=True)

    print(f"selftest workspace: {tmp}")
    r=subprocess.run([sys.executable, os.path.join(HERE,"run_pipeline.py"),
                      "--config", cfgp, "--skip-atac", "--only", "1"])
    if r.returncode==0:
        print("\nSELFTEST PASSED — step 1 (RNA DE) ran successfully.")
        print("Toolchain is correctly installed.")
    else:
        print("\nSELFTEST FAILED — check the traceback above and requirements.txt.")
        sys.exit(1)

if __name__=="__main__":
    main()
