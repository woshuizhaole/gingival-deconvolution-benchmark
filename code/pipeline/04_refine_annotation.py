#!/usr/bin/env python
"""Stage 1d - refine annotation based on manual marker review (documented):
 - drop cluster 21 (RBC: HBB/HBA1/HBA2 contamination) and cluster 3 (low-qual: MEG3/ribosomal, 50c)
 - reassign cluster 22 (pDC: IRF7/GZMB/PLD4/IRF8) from Plasma -> Myeloid (state=pDC)
Everything else keeps the automatic marker-score assignment (verified correct).
Output: data/processed/oral_ref_final.h5ad (the reference used for Stage 2).
"""
import scanpy as sc, numpy as np, pandas as pd, os
P="/data/workplace/wyx/oral_deconv_reference_20260707"
CK="leiden_r1.0"
a=sc.read_h5ad(f"{P}/data/processed/oral_ref_annotated.h5ad")
n0=a.n_obs
drop=["21","3"]
a=a[~a.obs[CK].isin(drop)].copy()
# reassign pDC cluster 22 coarse->Myeloid, state->pDC
ct=a.obs["cell_type"].astype(str).values
st=a.obs["cell_state"].astype(str).values
m22=(a.obs[CK].astype(str)=="22").values
ct[m22]="Myeloid"; st[m22]="pDC_c22"
a.obs["cell_type"]=pd.Categorical(ct)
a.obs["cell_state"]=pd.Categorical(st)
# ensure raw counts in X
if "counts" in a.layers: a.X=a.layers["counts"].copy()
a.write_h5ad(f"{P}/data/processed/oral_ref_final.h5ad")
tab=a.obs.groupby(["cell_type","condition"]).size().reset_index(name="n")
tab.to_csv(f"{P}/results/tables/stage1_final_celltype_condition.tsv",sep="\t",index=False)
vc=a.obs["cell_type"].value_counts()
vc.to_csv(f"{P}/results/tables/stage1_final_celltype_counts.tsv",sep="\t")
print(f"dropped {n0-a.n_obs} cells (RBC/lowqual); final {a.n_obs} cells x {a.n_vars} genes")
print(vc)
print("[done] oral_ref_final.h5ad")
