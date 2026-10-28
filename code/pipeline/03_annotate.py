#!/usr/bin/env python
"""
Stage 1c — annotate integrated gingiva atlas by marker-gene module scoring +
per-cluster assignment, plus rank_genes_groups diagnostics for verification.
Produces coarse cell_type and finer cell_state labels; preserves raw counts.

Outputs:
  data/processed/oral_ref_annotated.h5ad          (cell_type, cell_state, counts)
  results/tables/stage1c_cluster_annotation.tsv    (cluster -> cell_type, scores)
  results/tables/stage1c_cluster_top_markers.tsv   (top DE genes per cluster)
  results/tables/stage1c_celltype_counts.tsv
"""
import os, argparse
import numpy as np, pandas as pd, scanpy as sc

PROJ="/data/workplace/wyx/oral_deconv_reference_20260707"

MARKERS = {
 "Epithelial": ["KRT5","KRT14","KRT15","KRT4","KRT13","KRT1","KRT10","EPCAM","KRT8","KRT18","SFN","CDH1","IVL","SBSN","KRT19","KRT6A"],
 "Fibroblast": ["COL1A1","COL1A2","COL3A1","LUM","DCN","PDGFRA","FBLN1","MMP2","THY1","COL6A2","ELN"],
 "Endothelial": ["PECAM1","VWF","CDH5","CLDN5","FLT1","ENG","EGFL7","LYVE1","PROX1"],
 "Pericyte": ["RGS5","ACTA2","PDGFRB","NOTCH3","MYH11","TAGLN","DES"],
 "Tcell": ["CD3D","CD3E","CD2","TRAC","TRBC2","CD8A","CD4","IL7R","CD40LG","FOXP3","CCL5"],
 "NK": ["GNLY","NKG7","KLRD1","NCAM1","KLRF1","GZMB"],
 "Bcell": ["MS4A1","CD79A","CD79B","CD19","BANK1","HLA-DRA"],
 "Plasma": ["MZB1","JCHAIN","IGHG1","IGHA1","XBP1","SDC1","DERL3","PRDM1"],
 "Myeloid": ["LYZ","CD68","CD14","CD163","C1QA","C1QB","C1QC","AIF1","FCN1","S100A8","S100A9","CD1C","FCER1A","LILRA4","ITGAX"],
 "Mast": ["TPSAB1","TPSB2","CPA3","MS4A2","KIT"],
 "Neutrophil": ["FCGR3B","CSF3R","CXCR2","LCN2","S100A12","FUT4"],
 "Melanocyte": ["PMEL","MLANA","TYR","DCT","MITF"],
 "Glial": ["S100B","PLP1","SOX10","MPZ","NRXN1"],
}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--inp", default=f"{PROJ}/data/interim/gse164241_gingiva_integrated.h5ad")
    ap.add_argument("--out", default=f"{PROJ}/data/processed/oral_ref_annotated.h5ad")
    ap.add_argument("--cluster_key", default="leiden_r1.0")
    ap.add_argument("--force", action="store_true")
    args=ap.parse_args()
    if os.path.exists(args.out) and not args.force:
        print(f"[skip] {args.out} exists"); return
    np.random.seed(0)
    adata=sc.read_h5ad(args.inp)
    # ensure log-norm in X for scoring (raw holds log-norm all-gene)
    if adata.raw is not None:
        adata.X = adata.raw[:, adata.var_names].X.copy()
    print(f"[load] {adata.shape}, cluster_key={args.cluster_key}")

    # module scores
    for ct, genes in MARKERS.items():
        present=[g for g in genes if g in adata.var_names]
        sc.tl.score_genes(adata, present, score_name=f"score_{ct}", ctrl_size=50)
    score_cols=[f"score_{ct}" for ct in MARKERS]
    # z-score each score across cells, then per-cluster mean
    Z=adata.obs[score_cols].copy()
    Z=(Z-Z.mean())/ (Z.std()+1e-9)
    Z[args.cluster_key]=adata.obs[args.cluster_key].values
    clmean=Z.groupby(args.cluster_key)[score_cols].mean()
    assign={}
    rows=[]
    for cl,r in clmean.iterrows():
        order=r.sort_values(ascending=False)
        ct=order.index[0].replace("score_","")
        assign[cl]=ct
        rows.append(dict(cluster=cl, n_cells=int((adata.obs[args.cluster_key]==cl).sum()),
                         assigned=ct, top1=f"{order.index[0].replace('score_','')}:{order.iloc[0]:.2f}",
                         top2=f"{order.index[1].replace('score_','')}:{order.iloc[1]:.2f}",
                         top3=f"{order.index[2].replace('score_','')}:{order.iloc[2]:.2f}"))
    ann=pd.DataFrame(rows).sort_values("cluster")
    adata.obs["cell_type"]=adata.obs[args.cluster_key].map(assign).astype("category")
    adata.obs["cell_state"]=(adata.obs["cell_type"].astype(str)+"_c"+adata.obs[args.cluster_key].astype(str)).astype("category")

    # rank_genes_groups per cluster for verification
    sc.tl.rank_genes_groups(adata, args.cluster_key, method="wilcoxon", n_genes=20)
    rg=adata.uns["rank_genes_groups"]; groups=rg["names"].dtype.names
    top=[]
    for g in groups:
        genes=[rg["names"][g][i] for i in range(min(15,len(rg["names"][g])))]
        top.append(dict(cluster=g, assigned=assign.get(g,"?"), top_markers=",".join(genes)))
    topdf=pd.DataFrame(top)

    os.makedirs(f"{PROJ}/results/tables", exist_ok=True)
    ann.to_csv(f"{PROJ}/results/tables/stage1c_cluster_annotation.tsv", sep="\t", index=False)
    topdf.to_csv(f"{PROJ}/results/tables/stage1c_cluster_top_markers.tsv", sep="\t", index=False)
    ctc=adata.obs.groupby(["cell_type","condition"]).size().reset_index(name="n")
    ctc.to_csv(f"{PROJ}/results/tables/stage1c_celltype_counts.tsv", sep="\t", index=False)
    clmean.round(3).to_csv(f"{PROJ}/results/tables/stage1c_cluster_score_matrix.tsv", sep="\t")

    # restore raw counts as X for downstream deconvolution reference building
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adata.write_h5ad(args.out)
    print("[celltype counts]"); print(adata.obs["cell_type"].value_counts())
    print(f"[done] {args.out}")

if __name__=="__main__":
    main()
