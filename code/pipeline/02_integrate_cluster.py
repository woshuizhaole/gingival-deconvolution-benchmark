#!/usr/bin/env python
"""
Stage 1b — normalize, HVG, integrate (scVI primary on GPU; Harmony baseline),
cluster (Leiden at multiple resolutions). Input from Stage 1a.

Outputs:
  data/interim/gse164241_gingiva_integrated.h5ad
  logs/qc/stage1b_integration_summary.txt
Env: env/py. GPU chosen via --gpu (default auto-lowest-util).
"""
import os, argparse, subprocess, sys
import numpy as np, pandas as pd, scanpy as sc, anndata as ad

PROJ = "/data/workplace/wyx/oral_deconv_reference_20260707"

def pick_gpu():
    try:
        out = subprocess.check_output(
            ["nvidia-smi","--query-gpu=index,memory.used,utilization.gpu","--format=csv,noheader,nounits"]).decode()
        best=None
        for line in out.strip().splitlines():
            i,mem,util = [x.strip() for x in line.split(",")]
            score=int(mem)+int(util)*100
            if best is None or score<best[1]: best=(i,score)
        return best[0]
    except Exception:
        return "0"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--inp", default=f"{PROJ}/data/interim/gse164241_gingiva_qc.h5ad")
    ap.add_argument("--out", default=f"{PROJ}/data/interim/gse164241_gingiva_integrated.h5ad")
    ap.add_argument("--n_hvg", type=int, default=3000)
    ap.add_argument("--n_latent", type=int, default=30)
    ap.add_argument("--max_epochs", type=int, default=200)
    ap.add_argument("--gpu", default="auto")
    ap.add_argument("--force", action="store_true")
    args=ap.parse_args()
    if os.path.exists(args.out) and not args.force:
        print(f"[skip] {args.out} exists"); return
    gpu = pick_gpu() if args.gpu=="auto" else args.gpu
    os.environ["CUDA_VISIBLE_DEVICES"]=gpu
    print(f"[gpu] using CUDA_VISIBLE_DEVICES={gpu}")
    sc.settings.verbosity=1; np.random.seed(0)

    adata=sc.read_h5ad(args.inp)
    print(f"[load] {adata.shape}")
    # ensure counts layer
    if "counts" not in adata.layers: adata.layers["counts"]=adata.X.copy()
    # normalize for downstream (keep counts layer intact)
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata     # log-norm all-gene for marker scoring/annotation
    # HVG on raw counts (seurat_v3)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=args.n_hvg,
                                batch_key="sample", layer="counts")
    print(f"[hvg] {int(adata.var['highly_variable'].sum())} HVGs")

    # ---- scVI integration ----
    scvi_ok=False
    try:
        import scvi, torch
        scvi.settings.seed=0
        adhvg = adata[:, adata.var["highly_variable"]].copy()
        scvi.model.SCVI.setup_anndata(adhvg, layer="counts", batch_key="sample")
        m = scvi.model.SCVI(adhvg, n_layers=2, n_latent=args.n_latent, gene_likelihood="nb")
        m.train(max_epochs=args.max_epochs, early_stopping=True,
                plan_kwargs={"lr":1e-3}, batch_size=512,
                accelerator="gpu", devices=1)
        adata.obsm["X_scVI"]=m.get_latent_representation(adhvg)
        m.save(f"{PROJ}/results/models/scvi_gse164241_gingiva", overwrite=True)
        scvi_ok=True
        print(f"[scvi] latent {adata.obsm['X_scVI'].shape}")
    except Exception as e:
        print(f"[scvi][WARN] failed: {e}")

    # ---- Harmony baseline (CPU); PCA on log-norm HVG ----
    harmony_ok=False
    try:
        sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
        sc.external.pp.harmony_integrate(adata, key="sample")
        harmony_ok=True
        print("[harmony] done -> X_pca_harmony")
    except Exception as e:
        print(f"[harmony][WARN] failed: {e}")

    rep = "X_scVI" if scvi_ok else ("X_pca_harmony" if harmony_ok else "X_pca")
    if rep=="X_pca" and "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=50)
    print(f"[cluster] neighbors on {rep}")
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15)
    for res in [0.2,0.5,1.0,2.0]:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_r{res}", flavor="igraph",
                     n_iterations=2, directed=False)
        print(f"  leiden_r{res}: {adata.obs[f'leiden_r{res}'].nunique()} clusters")
    try:
        sc.tl.umap(adata)
    except Exception as e:
        print(f"[umap][WARN] {e}")

    adata.uns["integration"]={"scvi":scvi_ok,"harmony":harmony_ok,"primary_rep":rep,"gpu":gpu}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adata.write_h5ad(args.out)
    with open(f"{PROJ}/logs/qc/stage1b_integration_summary.txt","w") as f:
        f.write(f"cells x genes: {adata.shape}\nprimary_rep: {rep}\nscvi_ok: {scvi_ok}\nharmony_ok: {harmony_ok}\n")
        for res in [0.2,0.5,1.0,2.0]:
            f.write(f"leiden_r{res}: {adata.obs[f'leiden_r{res}'].nunique()} clusters\n")
    print(f"[done] {args.out}")

if __name__=="__main__":
    main()
