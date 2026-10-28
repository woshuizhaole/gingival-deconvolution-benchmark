#!/usr/bin/env python
"""
Stage 3-prep - process independent dataset GSE152042 (Caetano 2021 eLife gingiva,
2 health h1/h2 + 2 periodontitis p1/p2) into a labeled reference for leakage-free
pseudobulk generation (SimBu). Annotation uses the SAME marker sets as Stage 1c so
cell-type labels are consistent with the oral reference.

Outputs:
  data/interim/gse152042_labeled.h5ad
  data/processed/pseudobulk_source_10x/ (matrix/features/barcodes/cell_meta for SimBu)
  results/tables/stage3prep_gse152042_celltype_counts.tsv
"""
import os, glob, re
import scanpy as sc, numpy as np, pandas as pd, anndata as ad
import scipy.io, scipy.sparse as sp, gzip

P="/data/workplace/wyx/oral_deconv_reference_20260707"
SRC=f"{P}/data/interim/GSE152042"

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
    np.random.seed(0)
    files=sorted(glob.glob(f"{SRC}/*.h5"))
    ads=[]
    for f in files:
        samp=re.search(r"_(h\d|p\d)\.h5$", f).group(1)
        cond="health" if samp.startswith("h") else "periodontitis"
        a=sc.read_10x_h5(f); a.var_names_make_unique()
        a.obs_names=[f"{samp}_{b}" for b in a.obs_names]
        a.obs["sample"]=samp; a.obs["condition"]=cond
        sc.pp.filter_cells(a, min_genes=200)
        a.var["mt"]=a.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        a=a[a.obs["pct_counts_mt"]<20].copy()
        try:
            sc.pp.scrublet(a, expected_doublet_rate=0.06, random_state=0)
            a=a[~a.obs["predicted_doublet"]].copy()
        except Exception as e: print("scrublet warn",samp,e)
        ads.append(a); print(f"{samp}({cond}): {a.n_obs} cells")
    adata=ad.concat(ads, join="outer"); adata.obs_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"]=adata.X.copy()
    print(f"[concat] {adata.shape}")
    # normalize for clustering/scoring
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    adata.raw=adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    print(f"[cluster] {adata.obs['leiden'].nunique()} clusters")
    # score + per-cluster argmax (same scheme as Stage 1c)
    for ct,genes in MARKERS.items():
        present=[g for g in genes if g in adata.var_names]
        sc.tl.score_genes(adata, present, score_name=f"score_{ct}", ctrl_size=50)
    cols=[f"score_{ct}" for ct in MARKERS]
    Z=adata.obs[cols].copy(); Z=(Z-Z.mean())/(Z.std()+1e-9); Z["leiden"]=adata.obs["leiden"].values
    clmean=Z.groupby("leiden")[cols].mean()
    assign={cl:clmean.loc[cl].idxmax().replace("score_","") for cl in clmean.index}
    adata.obs["cell_type"]=adata.obs["leiden"].map(assign).astype("category")
    # restore counts to X for export
    adata.X=adata.layers["counts"].copy()
    os.makedirs(os.path.dirname(f"{P}/data/interim/gse152042_labeled.h5ad"), exist_ok=True)
    adata.write_h5ad(f"{P}/data/interim/gse152042_labeled.h5ad")
    # export labeled 10x for SimBu
    OUT=f"{P}/data/processed/pseudobulk_source_10x"; os.makedirs(OUT, exist_ok=True)
    X=sp.csr_matrix(adata.layers["counts"])
    scipy.io.mmwrite(f"{OUT}/matrix.mtx", X.T.tocoo(), field="integer"); os.system(f"gzip -f {OUT}/matrix.mtx")
    with gzip.open(f"{OUT}/features.tsv.gz","wt") as f:
        for g in adata.var_names: f.write(g+"\n")
    with gzip.open(f"{OUT}/barcodes.tsv.gz","wt") as f:
        for b in adata.obs_names: f.write(b+"\n")
    adata.obs[["cell_type","condition","sample"]].assign(cell_id=adata.obs_names)[["cell_id","cell_type","condition","sample"]].to_csv(f"{OUT}/cell_meta.tsv.gz", sep="\t", index=False, compression="gzip")
    tab=adata.obs.groupby(["cell_type","condition"]).size().reset_index(name="n")
    tab.to_csv(f"{P}/results/tables/stage3prep_gse152042_celltype_counts.tsv", sep="\t", index=False)
    print("[celltype counts]"); print(adata.obs["cell_type"].value_counts())
    print(f"[done] labeled GSE152042 -> {adata.n_obs} cells; exported to {OUT}")

if __name__=="__main__": main()
