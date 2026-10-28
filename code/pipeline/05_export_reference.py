#!/usr/bin/env python
"""
Stage 2a - export the annotated oral reference (raw counts + labels) in a
language-agnostic 10x-style layout so R tools (BayesPrism/InstaPrism/DWLS/MuSiC)
can load it robustly (no zellkonverter/basilisk dependency).

Outputs -> data/processed/reference_10x/
  matrix.mtx.gz     (genes x cells, raw integer counts)
  features.tsv.gz   (gene symbols)
  barcodes.tsv.gz   (cell ids)
  cell_meta.tsv.gz  (cell_id, cell_type, cell_state, condition, sample, donor)
Also writes a per-cell-type cell count table.
"""
import os, gzip
import scanpy as sc, numpy as np, pandas as pd
import scipy.io, scipy.sparse as sp

P="/data/workplace/wyx/oral_deconv_reference_20260707"
OUT=f"{P}/data/processed/reference_10x"
os.makedirs(OUT, exist_ok=True)

a=sc.read_h5ad(f"{P}/data/processed/oral_ref_final.h5ad")
# use raw counts
X = a.layers["counts"] if "counts" in a.layers else a.X
X = sp.csr_matrix(X)
# sanity: integer counts
assert X.min() >= 0, "counts negative?"
print(f"[load] {a.n_obs} cells x {a.n_vars} genes; counts max={X.max():.0f}")

# genes x cells
G = X.T.tocoo()   # genes(rows) x cells(cols)
scipy.io.mmwrite(f"{OUT}/matrix.mtx", G, field="integer")
os.system(f"gzip -f {OUT}/matrix.mtx")

with gzip.open(f"{OUT}/features.tsv.gz","wt") as f:
    for g in a.var_names: f.write(f"{g}\n")
with gzip.open(f"{OUT}/barcodes.tsv.gz","wt") as f:
    for b in a.obs_names: f.write(f"{b}\n")

meta = a.obs[["cell_type","cell_state","condition","sample","donor"]].copy()
meta.insert(0,"cell_id", a.obs_names)
meta.to_csv(f"{OUT}/cell_meta.tsv.gz", sep="\t", index=False, compression="gzip")

vc = a.obs["cell_type"].value_counts()
vc.to_csv(f"{P}/results/tables/stage2_reference_celltype_counts.tsv", sep="\t")
print("[export] wrote 10x-style reference to", OUT)
print(vc)
