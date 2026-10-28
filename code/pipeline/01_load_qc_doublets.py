#!/usr/bin/env python
"""
Stage 1a — load GSE164241 (Williams 2021 oral atlas) 10x per-sample matrices,
per-sample QC + doublet detection (Scrublet), MAD-based outlier filtering,
concatenate into a single annotated AnnData with raw counts preserved.

Reference-build tissue decision: gingiva (GM=healthy, PD=periodontitis) is the
primary reference for deconvolving gingival bulk. Buccal (BM) kept optional.

Outputs:
  data/interim/gse164241_gingiva_qc.h5ad   (default subset = gingiva)
  logs/qc/stage1a_qc_summary.tsv
Resumable: skips work if output exists (unless --force).
"""
import os, re, argparse, sys, gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io, scipy.sparse as sp

PROJ = "/data/workplace/wyx/oral_deconv_reference_20260707"
RAW  = f"{PROJ}/data/external/GSE164241"

def parse_sample(tok):
    """tok like GM136, PD164b, BM150, ID080UST -> (tissue, condition, donor)."""
    pre = re.match(r"^([A-Za-z]+)", tok).group(1)
    donor = re.sub(r"[a-z]$", "", tok)  # strip trailing lowercase replicate letter
    if pre == "GM":  return "gingiva", "health", donor
    if pre == "PD":  return "gingiva", "periodontitis", donor
    if pre == "BM":  return "buccal",  "health", donor
    return "other", "unknown", donor

def load_sample(prefix_path, tok, gsm):
    mtx = f"{prefix_path}/{gsm}_{tok}_matrix.mtx.gz"
    bc  = f"{prefix_path}/{gsm}_{tok}_barcodes.tsv.gz"
    ft  = f"{prefix_path}/{gsm}_{tok}_features.tsv.gz"
    M = scipy.io.mmread(mtx).tocsr().T.astype(np.float32)   # cells x genes
    barcodes = pd.read_csv(bc, header=None)[0].values
    feats = pd.read_csv(ft, header=None, sep="\t")
    # features: col0 ensembl, col1 symbol (2 cols observed)
    if feats.shape[1] >= 2:
        symbols = feats[1].astype(str).values
        ensembl = feats[0].astype(str).values
    else:
        symbols = feats[0].astype(str).values; ensembl = symbols
    A = ad.AnnData(X=M,
                   obs=pd.DataFrame(index=[f"{tok}_{b}" for b in barcodes]),
                   var=pd.DataFrame({"ensembl": ensembl}, index=symbols))
    A.var_names_make_unique()
    tissue, cond, donor = parse_sample(tok)
    A.obs["sample"] = tok; A.obs["gsm"] = gsm
    A.obs["tissue"] = tissue; A.obs["condition"] = cond; A.obs["donor"] = donor
    return A

def is_outlier(adata, metric, nmads):
    M = adata.obs[metric]
    med = np.median(M)
    mad = np.median(np.abs(M - med))
    return (M < med - nmads*mad) | (M > med + nmads*mad)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", default="gingiva", choices=["gingiva","gingiva+buccal","all"])
    ap.add_argument("--mt_cap", type=float, default=20.0)
    ap.add_argument("--min_genes", type=int, default=200)
    ap.add_argument("--out", default=f"{PROJ}/data/interim/gse164241_gingiva_qc.h5ad")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if os.path.exists(args.out) and not args.force:
        print(f"[skip] {args.out} exists; use --force to rebuild"); return
    sc.settings.verbosity = 1
    np.random.seed(0)

    toks = []
    for f in sorted(os.listdir(RAW)):
        m = re.match(r"^(GSM\d+)_([A-Za-z0-9]+)_matrix\.mtx\.gz$", f)
        if m: toks.append((m.group(1), m.group(2)))
    # subset by tissue
    def keep(tok):
        t,_,_ = parse_sample(tok)
        if args.subset == "gingiva":         return t == "gingiva"
        if args.subset == "gingiva+buccal":  return t in ("gingiva","buccal")
        return t != "other"     # 'all' still drops the ambiguous ID sample
    toks = [(g,t) for (g,t) in toks if keep(t)]
    print(f"[load] {len(toks)} samples for subset={args.subset}: {[t for _,t in toks]}")

    qc_rows = []; adatas = []
    for gsm, tok in toks:
        A = load_sample(RAW, tok, gsm)
        n0 = A.n_obs
        sc.pp.filter_cells(A, min_genes=args.min_genes)
        A.var["mt"] = A.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=[20], log1p=True)
        # doublets (Scrublet) per sample
        try:
            sc.pp.scrublet(A, expected_doublet_rate=0.06, random_state=0)
            ndbl = int(A.obs["predicted_doublet"].sum())
        except Exception as e:
            A.obs["predicted_doublet"] = False; A.obs["doublet_score"] = np.nan
            ndbl = -1; print(f"  [warn] scrublet failed {tok}: {e}")
        qc_rows.append(dict(sample=tok, gsm=gsm, tissue=A.obs['tissue'][0],
                            condition=A.obs['condition'][0], donor=A.obs['donor'][0],
                            n_raw=n0, n_pass_mingene=A.n_obs, n_doublet=ndbl,
                            median_counts=float(np.median(A.obs['total_counts'])),
                            median_genes=float(np.median(A.obs['n_genes_by_counts'])),
                            median_mt=float(np.median(A.obs['pct_counts_mt']))))
        adatas.append(A)
        print(f"  {tok}: {n0}->{A.n_obs} cells, {ndbl} doublets, medMT={np.median(A.obs['pct_counts_mt']):.1f}%")

    adata = ad.concat(adatas, join="outer", label="batch", index_unique=None)
    adata.obs_names_make_unique()
    print(f"[concat] {adata.n_obs} cells x {adata.n_vars} genes")

    # MAD-based outlier flagging (sc-best-practices): 5 MAD on counts/genes/top20, 3 MAD on mt + hard cap
    adata.obs["outlier"] = (is_outlier(adata,"log1p_total_counts",5)
                           | is_outlier(adata,"log1p_n_genes_by_counts",5)
                           | is_outlier(adata,"pct_counts_in_top_20_genes",5))
    adata.obs["mt_outlier"] = is_outlier(adata,"pct_counts_mt",3) | (adata.obs["pct_counts_mt"] > args.mt_cap)
    keep_mask = (~adata.obs["outlier"]) & (~adata.obs["mt_outlier"]) & (~adata.obs["predicted_doublet"])
    print(f"[qc] outlier={int(adata.obs['outlier'].sum())} mt_outlier={int(adata.obs['mt_outlier'].sum())} "
          f"doublet={int(adata.obs['predicted_doublet'].sum())} -> keep {int(keep_mask.sum())}/{adata.n_obs}")
    adata = adata[keep_mask].copy()
    # gene filter: expressed in >=3 cells
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"] = adata.X.copy()   # preserve raw counts
    print(f"[final] {adata.n_obs} cells x {adata.n_vars} genes")
    print(adata.obs.groupby(['tissue','condition']).size())

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adata.write_h5ad(args.out)
    qc = pd.DataFrame(qc_rows)
    qcp = f"{PROJ}/logs/qc/stage1a_qc_summary.tsv"
    os.makedirs(os.path.dirname(qcp), exist_ok=True)
    qc.to_csv(qcp, sep="\t", index=False)
    # composition table
    comp = adata.obs.groupby(['tissue','condition']).size().reset_index(name="n_cells")
    comp.to_csv(f"{PROJ}/logs/qc/stage1a_composition.tsv", sep="\t", index=False)
    print(f"[done] wrote {args.out}\n[done] qc -> {qcp}")

if __name__ == "__main__":
    main()
