#!/usr/bin/env python3
"""
05_figures.py — publication-quality synthesis figures.

Reads the tables produced by steps 01-04 (plus PCA tables in results/_work)
and writes Figure1..Figure5 PNGs to results/figures/.

Usage: python pipeline/05_figures.py [--config ...]
"""
import argparse, json, warnings, os, sys
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_config, rooted, ensure_dirs, table_path, work_path, figure_path, time_order, time_num, strip_prefix

_ap = argparse.ArgumentParser(); _ap.add_argument("--config", default=None)
_args = _ap.parse_args(); cfg = load_config(_args.config); ensure_dirs(cfg)
TIME_ORDER = time_order(cfg)
TX = [time_num(cfg)[t] for t in TIME_ORDER]
UP = None  # not used; kept for compatibility
W = work_path(cfg, "")  # work dir prefix
F = rooted(cfg, cfg["output"]["figures_dir"])
DPI = cfg["output"].get("figure_dpi", 300)
SPECIES = cfg["project"]["species_label"]
SPECIES_SHORT = "O. " + SPECIES.split()[-1] if " " in SPECIES else SPECIES

# ---- aesthetic ----
plt.rcParams.update({
 "font.family":"DejaVu Sans","font.size":9,"axes.linewidth":0.8,
 "axes.edgecolor":"#333333","axes.labelcolor":"#1a1a1a","text.color":"#1a1a1a",
 "xtick.color":"#333333","ytick.color":"#333333","axes.titleweight":"bold",
 "figure.dpi":120,"savefig.dpi":300,"axes.spines.top":False,"axes.spines.right":False,
 "legend.frameon":False,"axes.titlesize":11,"axes.titlepad":8,
})
# coral-reef inspired palette
CORAL="#E8634A"; HEAT="#D1483B"; CTRL="#3E7CB1"; OCEAN="#1B4965"
WARM=["#FFD7B5","#FF9E6D","#F2693C","#C0392B"]
TEAL="#2A9D8F"; SAND="#E9C46A"; PLUM="#6A4C93"; SLATE="#5C6B73"
CL_COL={1:"#E76F51",2:"#2A9D8F",3:"#E9C46A",4:"#6A4C93"}
# TIME_ORDER and TX come from config (set in header)

def save(fig,name):
    fig.savefig(f"{F}/{name}.png",bbox_inches="tight",facecolor="white",dpi=DPI)
    plt.close(fig); print("saved",name)

def parse_local(col):
    from common import parse_sample
    t, tr = parse_sample(col, cfg)
    return t, tr

def _ensure_pca():
    """Compute RNA & ATAC PCA tables into _work if absent (self-contained)."""
    from common import sample_metadata, save_json
    need = [work_path(cfg,f) for f in ["RNA_pca.csv","ATAC_pca.csv",
            "RNA_pca_var.json","ATAC_pca_var.json"]]
    if all(os.path.exists(p) for p in need):
        return
    def deseq_pca(counts, meta, tag):
        counts=counts[counts.sum(1)>=10]
        gm=np.exp(np.log(counts.replace(0,np.nan)).mean(1))
        sf=(counts.div(gm,axis=0)).median(axis=0)
        vst=np.log2(counts.div(sf,axis=1)+1)
        v=vst.var(1).sort_values(ascending=False).head(500).index
        X=vst.loc[v].T.values; X=X-X.mean(0)
        U,S,Vt=np.linalg.svd(X,full_matrices=False); pcs=U*S
        var=(S**2)/(S**2).sum()*100
        df=pd.DataFrame(pcs[:,:3],columns=["PC1","PC2","PC3"],index=vst.columns).join(meta)
        df.to_csv(work_path(cfg,f"{tag}_pca.csv"))
        save_json({"PC1":round(var[0],1),"PC2":round(var[1],1),"PC3":round(var[2],1)},
                  work_path(cfg,f"{tag}_pca_var.json"))
    # RNA (with configured PCA exclusions)
    rna=pd.read_csv(rooted(cfg,cfg["inputs"]["rna_counts"]),
                    index_col=cfg["inputs"]["rna_gene_id_col"])
    rna.columns=[c.strip() for c in rna.columns]
    excl=cfg.get("pca_exclude_samples",[]) or []
    keep=[c for c in rna.columns if not any(c.startswith(e) for e in excl)]
    rna=rna[keep]
    deseq_pca(rna, sample_metadata(rna.columns,cfg), "RNA")
    # ATAC
    atac=pd.read_csv(rooted(cfg,cfg["inputs"]["atac_counts"]))
    idc=cfg["inputs"]["atac_id_cols"]; sc=[c for c in atac.columns if c not in idc]
    deseq_pca(atac.set_index("peak_id")[sc], sample_metadata(sc,cfg), "ATAC")

_ensure_pca()

# =========================================================================
# FIGURE 1 — Experimental design + QC (PCA both assays + DE/DA counts)
# =========================================================================
def fig1():
    fig=plt.figure(figsize=(11,7.5))
    gs=gridspec.GridSpec(2,3,height_ratios=[1,1.05],hspace=0.42,wspace=0.34,
                         left=0.07,right=0.97,top=0.9,bottom=0.09)

    # --- design schematic (top-left, spans) ---
    ax0=fig.add_subplot(gs[0,0])
    ax0.axis("off")
    ax0.set_title("A  Experimental design",loc="left",fontsize=12)
    ax0.text(0.0,0.92,f"$\\it{{{SPECIES.split()[0]}\\ {SPECIES.split()[-1]}}}$",fontsize=10,fontweight="bold")
    ax0.text(0.0,0.82,cfg["project"].get("species_common",""),fontsize=8,style="italic",color=SLATE)
    # timeline (positions evenly spaced for however many timepoints there are)
    nt=len(TIME_ORDER)
    xs=np.linspace(0.08,0.86,nt)
    for t,xx in zip(TIME_ORDER,xs):
        ax0.add_patch(plt.Circle((xx,0.5),0.045,color=CORAL,zorder=3,transform=ax0.transAxes))
        ax0.text(xx,0.62,t,ha="center",fontsize=8,fontweight="bold")
    ax0.annotate("",xy=(0.93,0.5),xytext=(0.03,0.5),
                 arrowprops=dict(arrowstyle="-",lw=1.5,color=SLATE),zorder=1)
    # replicate count inferred from metadata
    from common import sample_metadata as _sm
    _rmeta=_sm([c for c in pd.read_csv(rooted(cfg,cfg["inputs"]["rna_counts"]),nrows=0).columns[1:]],cfg)
    _npg=int(_rmeta.groupby(["time","condition"]).size().median())
    _ng=pd.read_csv(rooted(cfg,cfg["inputs"]["rna_counts"]),usecols=[0]).shape[0]
    _np=pd.read_csv(rooted(cfg,cfg["inputs"]["atac_counts"]),usecols=["peak_id"]).peak_id.nunique()
    ax0.text(0.0,0.32,f"Heat (n={_npg}) vs Control (n={_npg})",fontsize=8,color=HEAT,fontweight="bold")
    ax0.text(0.0,0.20,"RNA-seq + ATAC-seq",fontsize=8,color=OCEAN,fontweight="bold")
    ax0.text(0.0,0.06,f"{_ng:,} genes · {_np:,} peaks",fontsize=7.5,color=SLATE)
    ax0.set_xlim(0,1); ax0.set_ylim(0,1)

    # --- RNA PCA (recompute excluding technical outlier OFav_05_C1) ---
    ax1=fig.add_subplot(gs[0,1])
    mk={"30min":"o","4h":"s","12h":"^","24h":"D"}
    # recompute PCA without the flagged outlier for interpretable structure
    rnac=pd.read_csv(rooted(cfg,cfg["inputs"]["rna_counts"]),index_col=0)
    rnac.columns=[c.strip() for c in rnac.columns]
    keepc=[c for c in rnac.columns if not c.startswith("OFav_05_C1")]
    rc=rnac[keepc]; rc=rc[rc.sum(1)>=10]
    gm=np.exp(np.log(rc.replace(0,np.nan)).mean(1))
    sf=(rc.div(gm,axis=0)).median(axis=0)
    vst=np.log2(rc.div(sf,axis=1)+1)
    v=vst.var(1).sort_values(ascending=False).head(500).index
    Xp=vst.loc[v].T.values; Xp=Xp-Xp.mean(0)
    U,S,Vt=np.linalg.svd(Xp,full_matrices=False); pcs=U*S
    varr=(S**2)/(S**2).sum()*100
    pmeta=pd.DataFrame([(c,)+parse_local(c) for c in vst.columns],
                       columns=["sample","time","condition"]).set_index("sample")
    pdf=pd.DataFrame(pcs[:,:2],columns=["PC1","PC2"],index=vst.columns).join(pmeta)
    for _,r in pdf.iterrows():
        ax1.scatter(r.PC1,r.PC2,c=(HEAT if r.condition=="Heat" else CTRL),
                    marker=mk[r.time],s=70,edgecolor="white",lw=0.8,zorder=3)
    ax1.set_xlabel(f"PC1 ({varr[0]:.1f}%)"); ax1.set_ylabel(f"PC2 ({varr[1]:.1f}%)")
    ax1.set_title("B  RNA-seq PCA",loc="left",fontsize=12)
    ax1.text(0.5,-0.30,"(1 technical outlier excluded)",transform=ax1.transAxes,
             ha="center",fontsize=6.5,color=SLATE,style="italic")

    # --- ATAC PCA ---
    ax2=fig.add_subplot(gs[0,2])
    pcaa=pd.read_csv(work_path(cfg,"ATAC_pca.csv"),index_col=0)
    vara=json.load(open(work_path(cfg,"ATAC_pca_var.json")))
    for _,r in pcaa.iterrows():
        ax2.scatter(r.PC1,r.PC2,c=(HEAT if r.condition=="Heat" else CTRL),
                    marker=mk[r.time],s=70,edgecolor="white",lw=0.8,zorder=3)
    ax2.set_xlabel(f"PC1 ({vara['PC1']}%)"); ax2.set_ylabel(f"PC2 ({vara['PC2']}%)")
    ax2.set_title("C  ATAC-seq PCA",loc="left",fontsize=12)
    # shared legend
    leg=[Line2D([0],[0],marker="o",color="w",markerfacecolor=HEAT,markersize=8,label="Heat"),
         Line2D([0],[0],marker="o",color="w",markerfacecolor=CTRL,markersize=8,label="Control")]
    leg+=[Line2D([0],[0],marker=mk[t],color="w",markerfacecolor=SLATE,markersize=7,label=t) for t in TIME_ORDER]
    ax2.legend(handles=leg,fontsize=6.5,ncol=2,loc="upper right",handletextpad=0.2,columnspacing=0.8)

    # --- DE bar (bottom-left, span 2) ---
    ax3=fig.add_subplot(gs[1,:2])
    rna=pd.read_csv(table_path(cfg,"RNA_DE_by_timepoint.csv"))
    up=[((rna.time==t)&(rna.padj<0.05)&(rna.log2FoldChange>=1)).sum() for t in TIME_ORDER]
    dn=[((rna.time==t)&(rna.padj<0.05)&(rna.log2FoldChange<=-1)).sum() for t in TIME_ORDER]
    x=np.arange(4)
    ax3.bar(x,up,color=HEAT,label="Induced (Heat>Control)",edgecolor="white",lw=0.6)
    ax3.bar(x,[-d for d in dn],color=CTRL,label="Repressed",edgecolor="white",lw=0.6)
    for i,(u,d) in enumerate(zip(up,dn)):
        ax3.text(i,u+40,str(u),ha="center",fontsize=8,color=HEAT,fontweight="bold")
        ax3.text(i,-d-40,str(d),ha="center",va="top",fontsize=8,color=CTRL,fontweight="bold")
    ax3.axhline(0,color="#333",lw=0.8)
    ax3.set_xticks(x); ax3.set_xticklabels(TIME_ORDER)
    ax3.set_ylabel("Differentially expressed genes")
    ax3.set_title("D  RNA differential expression (padj<0.05, |log$_2$FC|≥1)",loc="left",fontsize=11)
    ax3.legend(fontsize=7.5,loc="lower right")
    ax3.set_ylim(-2100,3200)

    # --- ATAC DA bar (bottom-right) ---
    ax4=fig.add_subplot(gs[1,2])
    atac=pd.read_csv(table_path(cfg,"ATAC_DA_by_timepoint.csv"))
    upa=[((atac.time==t)&(atac.pvalue<0.05)&(atac.log2FoldChange>=0.3)).sum() for t in TIME_ORDER]
    dna=[((atac.time==t)&(atac.pvalue<0.05)&(atac.log2FoldChange<=-0.3)).sum() for t in TIME_ORDER]
    ax4.bar(x,upa,color=CORAL,label="More open",edgecolor="white",lw=0.6)
    ax4.bar(x,[-d for d in dna],color=OCEAN,label="More closed",edgecolor="white",lw=0.6)
    ax4.axhline(0,color="#333",lw=0.8)
    ax4.set_xticks(x); ax4.set_xticklabels(TIME_ORDER,fontsize=8)
    ax4.set_ylabel("DA peaks (nominal p<0.05)")
    ax4.set_title("E  ATAC differential\naccessibility",loc="left",fontsize=10)
    ax4.legend(fontsize=6.5,loc="upper right")

    fig.suptitle("Figure 1 · Heat stress drives a strong, time-structured transcriptional response in $\\it{O.\\ faveolata}$",
                 fontsize=12.5,fontweight="bold",x=0.07,ha="left",y=0.975)
    save(fig,"Figure1_design_QC")
fig1()

# =========================================================================
# FIGURE 2 — Fuzzy c-means temporal programs + GO
# =========================================================================
def fig2():
    Xs=np.load(work_path(cfg,"fuzzy_Xs.npy"))
    genes=pd.read_csv(work_path(cfg,"fuzzy_genes.csv")).gene_id.tolist()
    clust=pd.read_csv(table_path(cfg,"fuzzy_clusters.csv")).set_index("gene_id")
    centR=pd.read_csv(table_path(cfg,"fuzzy_centroids_raw_log2fc.csv"),index_col=0)
    xb=pd.read_csv(table_path(cfg,"fuzzy_xie_beni.csv"))
    enr=pd.read_csv(table_path(cfg,"GO_enrichment_all.csv"))

    cl=clust.loc[genes]
    fig=plt.figure(figsize=(12,8.5))
    gs=gridspec.GridSpec(3,4,height_ratios=[1,0.9,0.9],hspace=0.55,wspace=0.45,
                         left=0.06,right=0.97,top=0.91,bottom=0.07)

    # row0: trajectories per cluster
    cluster_labels={1:"Early-transient",2:"Delayed-induction",3:"Biphasic",4:"Mid-sustained"}
    for ci,c in enumerate([1,2,3,4]):
        ax=fig.add_subplot(gs[0,ci])
        gsub=cl.index[(cl.cluster==c)&(cl.core)]
        idx=[genes.index(g) for g in gsub]
        for i in idx:
            ax.plot(TX,Xs[i],color=CL_COL[c],alpha=0.06,lw=0.5)
        cent=Xs[[genes.index(g) for g in gsub]].mean(0)
        ax.plot(TX,cent,color="#1a1a1a",lw=2.4,zorder=5)
        ax.plot(TX,cent,color=CL_COL[c],lw=1.4,zorder=6)
        ax.set_xscale("log"); ax.set_xticks(TX); ax.set_xticklabels(TIME_ORDER,fontsize=7)
        ax.axhline(0,color=SLATE,lw=0.5,ls=":")
        ax.set_title(f"C{c} · {cluster_labels[c]}\n(n={len(gsub)} core)",fontsize=9)
        if ci==0: ax.set_ylabel("z-scored log$_2$FC")
        ax.set_ylim(-2.2,2.2)
        ax.tick_params(labelsize=7)

    # row1-left: stacked heatmap of all core genes ordered by cluster
    axh=fig.add_subplot(gs[1:,0:2])
    order=[]; bounds=[0]
    for c in [1,2,3,4]:
        gsub=cl[(cl.cluster==c)&(cl.core)].sort_values("membership",ascending=False).index
        order+=list(gsub); bounds.append(len(order))
    M=Xs[[genes.index(g) for g in order]]
    im=axh.imshow(M,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2,interpolation="nearest")
    axh.set_xticks(range(4)); axh.set_xticklabels(TIME_ORDER,fontsize=8)
    axh.set_yticks([]); axh.set_ylabel("Genes (ordered by cluster, core membership)")
    for b in bounds[1:-1]: axh.axhline(b-0.5,color="white",lw=1.5)
    for i,c in enumerate([1,2,3,4]):
        mid=(bounds[i]+bounds[i+1])/2
        axh.text(-0.6,mid,f"C{c}",ha="right",va="center",fontsize=9,fontweight="bold",color=CL_COL[c])
    axh.set_title("B  Temporal expression archetypes",loc="left",fontsize=11)
    cb=fig.colorbar(im,ax=axh,fraction=0.035,pad=0.02); cb.set_label("z-score",fontsize=7); cb.ax.tick_params(labelsize=6)

    # row1-right: Xie-Beni
    axx=fig.add_subplot(gs[1,2])
    axx.plot(xb.k,xb.xie_beni,"-o",color=OCEAN,lw=1.5,ms=5)
    kmin=xb.loc[xb.xie_beni.idxmin(),"k"]
    axx.scatter([kmin],[xb.xie_beni.min()],s=120,facecolor="none",edgecolor=CORAL,lw=2,zorder=5)
    axx.set_xlabel("Clusters (k)"); axx.set_ylabel("Xie–Beni index")
    axx.set_title("C  Cluster validity",loc="left",fontsize=10)
    axx.annotate(f"optimum\nk={int(kmin)}",xy=(kmin,xb.xie_beni.min()),
                 xytext=(kmin+1.3,xb.xie_beni.min()+0.06),fontsize=7,color=CORAL,
                 arrowprops=dict(arrowstyle="->",color=CORAL,lw=1))

    # row1-right2: cluster sizes
    axs=fig.add_subplot(gs[1,3])
    sizes=[(c,(cl.cluster==c).sum(),((cl.cluster==c)&cl.core).sum()) for c in [1,2,3,4]]
    sd=pd.DataFrame(sizes,columns=["c","all","core"])
    axs.barh(sd.c,sd["all"],color="#dddddd",label="all")
    axs.barh(sd.c,sd.core,color=[CL_COL[c] for c in sd.c],label="core (≥0.6)")
    axs.set_yticks([1,2,3,4]); axs.set_yticklabels([f"C{c}" for c in [1,2,3,4]],fontsize=8)
    axs.invert_yaxis(); axs.set_xlabel("genes"); axs.set_title("D  Cluster size",loc="left",fontsize=10)
    axs.legend(fontsize=6.5,loc="lower right")

    # row2-right: GO terms per cluster (top 3 named each)
    axg=fig.add_subplot(gs[2,2:])
    axg.axis("off"); axg.set_title("E  Representative enriched GO terms",loc="left",fontsize=10)
    yy=0.96
    for c in [1,2,3,4]:
        sub=enr[(enr.cluster_or_set==f"Cluster_{c}")&(enr.padj<0.1)
                &(enr.go_id!=enr.go_name)].drop_duplicates("go_name").head(3)
        axg.text(0.0,yy,f"C{c}",fontsize=8.5,fontweight="bold",color=CL_COL[c]); 
        for _,r in sub.iterrows():
            yy-=0.082
            axg.text(0.10,yy+0.0,f"{r.go_name}",fontsize=7,color="#222")
            axg.text(0.97,yy,f"{r.fold_enrich:.0f}×",fontsize=6.5,color=SLATE,ha="right")
        yy-=0.10
    axg.set_xlim(0,1); axg.set_ylim(0,1)

    fig.suptitle("Figure 2 · Four temporal regulatory programs partition the heat-stress transcriptome",
                 fontsize=12.5,fontweight="bold",x=0.06,ha="left",y=0.975)
    save(fig,"Figure2_fuzzy_clusters")
fig2()

# =========================================================================
# FIGURE 3 — RNA-ATAC integration: concordance, temporal lag, scatter
# =========================================================================
def fig3():
    cross=pd.read_csv(table_path(cfg,"RNA_ATAC_cross_time_scored.csv"))
    same=pd.read_csv(table_path(cfg,"RNA_ATAC_same_time.csv"))
    lag=pd.read_csv(table_path(cfg,"RNA_ATAC_lagged.csv"))
    summ=json.load(open(work_path(cfg,"integration_summary.json")))

    fig=plt.figure(figsize=(12,7.6))
    gs=gridspec.GridSpec(2,3,height_ratios=[1,1],hspace=0.42,wspace=0.36,
                         left=0.07,right=0.96,top=0.9,bottom=0.09)

    # A: 4-quadrant scatter of cross-time links
    ax=fig.add_subplot(gs[0,0])
    cc={True:CORAL,False:SLATE}
    ax.scatter(cross.atac_log2FC,cross.rna_log2FC,
               c=[cc[x] for x in cross.concordant],s=14,alpha=0.5,edgecolor="none")
    ax.axhline(0,color="#333",lw=0.6); ax.axvline(0,color="#333",lw=0.6)
    ax.set_xlabel("ATAC log$_2$FC (accessibility)"); ax.set_ylabel("RNA log$_2$FC (expression)")
    ax.set_title("A  Coupled chromatin–expression\nchanges (gene-level)",loc="left",fontsize=10)
    lim=max(abs(cross.atac_log2FC).quantile(.99),abs(cross.rna_log2FC).quantile(.99))
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    rate=summ["concordance_rate_pct"]
    ax.text(0.97,0.04,f"{rate:.0f}% concordant",transform=ax.transAxes,ha="right",
            fontsize=8,color=CORAL,fontweight="bold")
    for (xx,yy,t) in [(0.97,0.97,"open+up"),(0.03,0.97,"closed+up"),
                      (0.03,0.03,"closed+down"),(0.97,0.03,"")]:
        if t: ax.text(xx,yy,t,transform=ax.transAxes,ha=("right" if xx>0.5 else "left"),
                      va="top",fontsize=6.5,color=SLATE,style="italic")

    # B: temporal category donut
    ax=fig.add_subplot(gs[0,1])
    cats=["ATAC_first","Simultaneous","RNA_first"]
    vals=[summ["atac_first"],summ["simultaneous"],summ["rna_first"]]
    cols=[CORAL,SAND,OCEAN]
    w,_=ax.pie(vals,colors=cols,startangle=90,counterclock=False,
               wedgeprops=dict(width=0.42,edgecolor="white",lw=1.5))
    ax.text(0,0,f"{sum(vals)}\nlinks",ha="center",va="center",fontsize=10,fontweight="bold")
    ax.set_title("B  Temporal ordering of\nchromatin vs expression",loc="left",fontsize=10)
    lab=["ATAC opens first","Simultaneous","RNA changes first"]
    ax.legend([Patch(fc=c) for c in cols],
              [f"{l} ({v})" for l,v in zip(lab,vals)],fontsize=7,loc="lower center",
              bbox_to_anchor=(0.5,-0.18))

    # C: lag links bar by pair, colored by concordance
    ax=fig.add_subplot(gs[0,2])
    lp=(lag.groupby("link_type")
          .agg(n=("gene_id","size"),conc=("concordant","sum")).reset_index())
    order=["lag_30min_to_4h","lag_30min_to_12h","lag_4h_to_12h","lag_4h_to_24h","lag_12h_to_24h"]
    lp=lp.set_index("link_type").reindex(order).dropna()
    labels=[x.replace("lag_","").replace("_to_","→") for x in lp.index]
    y=np.arange(len(lp))
    ax.barh(y,lp.n,color="#d8d8d8",label="all links")
    ax.barh(y,lp.conc,color=CORAL,label="concordant")
    ax.set_yticks(y); ax.set_yticklabels(labels,fontsize=7.5); ax.invert_yaxis()
    ax.set_xlabel("gene–peak links"); ax.set_title("C  ATAC→RNA temporal lags",loc="left",fontsize=10)
    ax.legend(fontsize=6.5,loc="lower right")

    # D: hexbin density of same-time links (priming signal)
    ax=fig.add_subplot(gs[1,0])
    hb=ax.hexbin(same.atac_log2FC,same.rna_log2FC,gridsize=22,cmap="magma_r",mincnt=1)
    ax.axhline(0,color="white",lw=0.5); ax.axvline(0,color="white",lw=0.5)
    ax.set_xlabel("ATAC log$_2$FC"); ax.set_ylabel("RNA log$_2$FC")
    ax.set_title(f"D  Same-timepoint links (n={len(same)})",loc="left",fontsize=10)
    cb=fig.colorbar(hb,ax=ax,fraction=0.04,pad=0.02); cb.set_label("links",fontsize=7); cb.ax.tick_params(labelsize=6)

    # E: candidate score distribution by feature
    ax=fig.add_subplot(gs[1,1])
    cross["is_promoter"]=cross.top_feature.eq("promoter")
    for grp,col,lab in [(True,CORAL,"promoter peaks"),(False,SLATE,"other features")]:
        d=cross[cross.is_promoter==grp].candidate_score.dropna()
        if len(d)>5:
            xs=np.linspace(d.min(),d.max(),200); k=gaussian_kde(d)
            ax.fill_between(xs,k(xs),alpha=0.35,color=col); ax.plot(xs,k(xs),color=col,lw=1.5,label=lab)
    ax.set_xlabel("Integration candidate score"); ax.set_ylabel("density")
    ax.set_title("E  Promoter peaks score higher",loc="left",fontsize=10)
    ax.legend(fontsize=7)

    # F: summary stat panel
    ax=fig.add_subplot(gs[1,2]); ax.axis("off")
    ax.set_title("F  Integration summary",loc="left",fontsize=10)
    rows=[("RNA DE genes",f"{summ['rna_sig_genes']:,}"),
          ("ATAC DA peaks (gene-linked)",f"{summ['atac_sig_peaks']:,}"),
          ("Genes w/ DA peak",f"{summ['atac_sig_genes']:,}"),
          ("Same-time links",f"{summ['same_time_links']:,}"),
          ("Cross-time links",f"{summ['cross_time_links']:,}"),
          ("Concordance rate",f"{summ['concordance_rate_pct']:.0f}%"),
          ("Open+induced candidates",f"{summ['overexpression_candidates']:,}")]
    yy=0.92
    for k,v in rows:
        ax.text(0.0,yy,k,fontsize=8,color="#333"); ax.text(1.0,yy,v,fontsize=8.5,
                fontweight="bold",ha="right",color=OCEAN); yy-=0.135
    ax.set_xlim(0,1); ax.set_ylim(0,1)

    fig.suptitle("Figure 3 · Integrating chromatin accessibility with expression reveals coupled, time-ordered regulation",
                 fontsize=12.5,fontweight="bold",x=0.07,ha="left",y=0.975)
    save(fig,"Figure3_integration")
fig3()

# =========================================================================
# FIGURE 4 — Candidate heat-inducible promoters (applied payoff)
# =========================================================================
def fig4():
    cand=pd.read_csv(table_path(cfg,"S1_candidate_heat_inducible_promoters.csv"))
    rna=pd.read_csv(table_path(cfg,"RNA_DE_by_timepoint.csv"))
    rna["gene_id"]=rna.gene_id.str.replace("^gene:","",regex=True)

    fig=plt.figure(figsize=(12,7))
    gs=gridspec.GridSpec(2,3,width_ratios=[1.25,1,1],height_ratios=[1,1],
                         hspace=0.5,wspace=0.42,left=0.09,right=0.96,top=0.89,bottom=0.1)

    # A: top-20 candidate lollipop (RNA & ATAC LFC)
    ax=fig.add_subplot(gs[:,0])
    top=cand.head(18).iloc[::-1]
    y=np.arange(len(top))
    ax.hlines(y,0,top.rna_log2FC,color=HEAT,lw=1,alpha=0.6)
    ax.scatter(top.rna_log2FC,y,color=HEAT,s=45,label="RNA log$_2$FC",zorder=5,edgecolor="white",lw=0.6)
    ax.scatter(top.atac_log2FC,y,color=CORAL,marker="D",s=38,label="ATAC log$_2$FC",zorder=5,edgecolor="white",lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels(top.gene_id,fontsize=6.8)
    ax.set_xlabel("log$_2$ fold change (Heat vs Control)")
    ax.axvline(0,color="#333",lw=0.6)
    ax.set_title("A  Top-ranked candidate\nheat-inducible promoters",loc="left",fontsize=10)
    ax.legend(fontsize=7,loc="lower right")

    # B: temporal category composition of candidates
    ax=fig.add_subplot(gs[0,1])
    tc=cand.temporal_category.value_counts().reindex(["ATAC_first","Simultaneous","RNA_first"]).fillna(0)
    ax.bar(range(3),tc.values,color=[CORAL,SAND,OCEAN],edgecolor="white",lw=0.8)
    ax.set_xticks(range(3)); ax.set_xticklabels(["ATAC\nfirst","Simul-\ntaneous","RNA\nfirst"],fontsize=7.5)
    ax.set_ylabel("candidates"); ax.set_title("B  Regulatory timing",loc="left",fontsize=10)
    for i,v in enumerate(tc.values): ax.text(i,v+0.5,int(v),ha="center",fontsize=8,fontweight="bold")

    # C: score vs significance scatter
    ax=fig.add_subplot(gs[0,2])
    sc=ax.scatter(cand.rna_log2FC,cand.atac_log2FC,
                  c=cand.candidate_score,cmap="YlOrRd",s=30,edgecolor="#666",lw=0.3)
    ax.set_xlabel("RNA log$_2$FC"); ax.set_ylabel("ATAC log$_2$FC")
    ax.set_title("C  Candidate landscape",loc="left",fontsize=10)
    cb=fig.colorbar(sc,ax=ax,fraction=0.045,pad=0.02); cb.set_label("score",fontsize=7); cb.ax.tick_params(labelsize=6)

    # D: temporal expression trajectory of top 6 candidates
    ax=fig.add_subplot(gs[1,1:])
    top6=cand.head(6).gene_id.tolist()
    cmap6=plt.cm.viridis(np.linspace(0,0.85,6))
    for gi,g in enumerate(top6):
        sub=rna[rna.gene_id==g].set_index("time").reindex(TIME_ORDER)
        ax.plot(TX,sub.log2FoldChange.values,"-o",color=cmap6[gi],lw=1.5,ms=4,label=g)
    ax.set_xscale("log"); ax.set_xticks(TX); ax.set_xticklabels(TIME_ORDER,fontsize=8)
    ax.axhline(0,color=SLATE,lw=0.5,ls=":")
    ax.set_xlabel("Time after heat stress"); ax.set_ylabel("RNA log$_2$FC")
    ax.set_title("D  Expression trajectories of top candidates",loc="left",fontsize=10)
    ax.legend(fontsize=6.5,ncol=2,loc="upper right")

    fig.suptitle("Figure 4 · Integration nominates promoter-proximal, chromatin-primed candidate heat-inducible loci",
                 fontsize=12.5,fontweight="bold",x=0.09,ha="left",y=0.965)
    save(fig,"Figure4_candidates")
fig4()

# =========================================================================
# FIGURE 5 — Synthesis model: the temporal heat-stress regulatory cascade
# =========================================================================
def fig5():
    summ=json.load(open(work_path(cfg,"integration_summary.json")))
    rna=pd.read_csv(table_path(cfg,"RNA_DE_by_timepoint.csv"))
    fig=plt.figure(figsize=(12,6.8))
    gs=gridspec.GridSpec(2,1,height_ratios=[1.4,1],hspace=0.32,left=0.06,right=0.96,top=0.9,bottom=0.07)

    # TOP: cascade schematic
    ax=fig.add_subplot(gs[0]); ax.axis("off")
    ax.set_xlim(0,10); ax.set_ylim(0,6)
    ax.set_title("A  A time-ordered heat-stress regulatory cascade in $\\it{O.\\ faveolata}$",loc="left",fontsize=12)
    # timeline band
    txpos={"30min":1.3,"4h":3.7,"12h":6.3,"24h":8.8}
    ax.plot([0.8,9.3],[5.2,5.2],color=SLATE,lw=1.2,zorder=1)
    for t,xx in txpos.items():
        ax.add_patch(plt.Circle((xx,5.2),0.12,color=CORAL,zorder=3))
        ax.text(xx,5.55,t,ha="center",fontsize=9,fontweight="bold")
    # layers
    def box(x,y,w,h,text,col,fc=None):
        from matplotlib.patches import FancyBboxPatch
        p=FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.04,rounding_size=0.12",
                         linewidth=1.2,edgecolor=col,facecolor=(fc or col),alpha=(0.92 if fc else 0.18),zorder=4)
        ax.add_patch(p)
        ax.text(x,y,text,ha="center",va="center",fontsize=7.6,zorder=5,
                color=("white" if fc else "#222"),fontweight=("bold" if fc else "normal"))
    # chromatin layer
    ax.text(0.15,4.55,"Chromatin",fontsize=8,fontweight="bold",color=CORAL,rotation=0)
    box(txpos["30min"],4.0,2.0,0.62,"Rapid accessibility\nremodeling (2,156 peaks)",CORAL,CORAL)
    box(txpos["4h"],4.0,1.9,0.62,"Promoter priming\nat stress loci",CORAL)
    # expression layer
    ax.text(0.15,3.35,"Expression",fontsize=8,fontweight="bold",color=OCEAN)
    box(txpos["30min"],2.8,2.0,0.62,"Immediate-early\ninduction (2,378 DE)",OCEAN,OCEAN)
    box(txpos["4h"],2.8,1.9,0.62,"Chaperones / proteostasis\n(C4: HSP, CCT)",PLUM,PLUM)
    box(txpos["12h"],2.8,2.0,0.62,"Peak response (2,819 DE)\nDNA repair, cell cycle (C2)",OCEAN,OCEAN)
    box(txpos["24h"],2.8,1.9,0.62,"Resolution\n(438 DE) — recovery",TEAL,TEAL)
    # function layer
    ax.text(0.15,1.95,"Outcome",fontsize=8,fontweight="bold",color=SLATE)
    box(txpos["30min"],1.4,2.0,0.55,"Signaling, apoptosis\npriming (C3)",SAND)
    box(txpos["12h"],1.4,2.0,0.55,"Proteostatic + genome\nmaintenance program",SAND)
    box(txpos["24h"],1.4,1.9,0.55,"Homeostatic\nrecovery",SAND)
    # arrows chromatin->expression
    for xx in [txpos["30min"],txpos["4h"]]:
        ax.annotate("",xy=(xx,3.15),xytext=(xx,3.68),
                    arrowprops=dict(arrowstyle="-|>",color=CORAL,lw=1.6))
    ax.text(9.7,4.0,"~27%\nATAC\nfirst",fontsize=6.5,color=CORAL,ha="center")

    # BOTTOM: global response magnitude curve
    ax2=fig.add_subplot(gs[1])
    de=[((rna.time==t)&(rna.padj<0.05)&(rna.log2FoldChange.abs()>=1)).sum() for t in TIME_ORDER]
    ax2.fill_between(TX,de,color=CORAL,alpha=0.2)
    ax2.plot(TX,de,"-o",color=HEAT,lw=2,ms=7,zorder=5)
    for x,v in zip(TX,de): ax2.annotate(f"{v:,}",(x,v),textcoords="offset points",
                                        xytext=(0,9),ha="center",fontsize=8,fontweight="bold",color=HEAT)
    ax2.set_xscale("log"); ax2.set_xticks(TX); ax2.set_xticklabels(TIME_ORDER)
    ax2.set_xlabel("Time after heat stress"); ax2.set_ylabel("DE genes\n(padj<0.05,|log$_2$FC|≥1)")
    ax2.set_title("B  Transcriptional response peaks at 12 h, then resolves by 24 h",loc="left",fontsize=11)
    ax2.set_ylim(0,3300)
    ax2.annotate("acute\ninduction",(0.5,2378),xytext=(0.65,2950),fontsize=7,color=SLATE,
                 arrowprops=dict(arrowstyle="->",color=SLATE,lw=0.8))
    ax2.annotate("resolution",(24,438),xytext=(13,1100),fontsize=7,color=TEAL,
                 arrowprops=dict(arrowstyle="->",color=TEAL,lw=0.8))

    fig.suptitle("Figure 5 · Synthesis — chromatin remodeling precedes a wave of expression that crests and resolves",
                 fontsize=12.5,fontweight="bold",x=0.06,ha="left",y=0.975)
    save(fig,"Figure5_synthesis")
fig5()
