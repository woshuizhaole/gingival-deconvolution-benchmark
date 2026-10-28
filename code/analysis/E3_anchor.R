#!/usr/bin/env Rscript
# E3 - anchor the E1/E2-corrected ABSOLUTE proportions to orthogonal ground truth.
# (a) HEpiDISH Stage-2 breaks the DNAm immune compartment into subtypes on the SAME 24 donors
#     (GSE173081) and correlates with the RNA hierarchical immune subtypes (E2, GSE173078).
# (b) Compare within-immune composition to published periodontitis histology (plasma-dominant).
suppressMessages({library(data.table); library(EpiDISH)})
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
OUT <- file.path(P,"results/tables/improvement"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
ccc <- function(t,e) 2*cov(t,e)/(var(t)+var(e)+(mean(t)-mean(e))^2)

## ---- DNAm immune subtypes via HEpiDISH (epi/fib/IC -> 7 blood immune subsets) ----
beta <- readRDS(file.path(P,"data/processed/GSE173081_beta.rds"))
data(centEpiFibIC.m); data(centBloodSub.m)
h <- hepidish(beta.m=beta, ref1.m=centEpiFibIC.m, ref2.m=centBloodSub.m, h.CT.idx=3, method="RPC")
# h columns: Epi, Fib, then blood subsets (B, NK, CD4T, CD8T, Mono, Neutro, Eoso)
dnam <- as.data.frame(h); dnam$sample <- rownames(h)
fwrite(data.table(sample=rownames(h), h), file.path(OUT,"E3_DNAm_HEpiDISH_subtypes.tsv"), sep="\t")
cat("=== DNAm HEpiDISH 列 ===\n"); print(colnames(h)); cat("均值:\n"); print(round(colMeans(h),3))

## ---- RNA hierarchical within-immune subtypes (E2, GSE173078) ----
r <- fread(file.path(OUT,"GSE173078_hierarchical.tsv")); rs <- as.matrix(r[,-1]); rownames(rs) <- r$sample
imm <- c("Tcell","Bcell","Plasma","Myeloid","Mast","NK","Neutrophil")
rna_imm <- rs[, imm, drop=FALSE]; rna_imm <- rna_imm/rowSums(rna_imm)   # within-immune

## ---- match donors (AC ids) ----
common <- intersect(rownames(rna_imm), dnam$sample)
cat("\n[match] ", length(common), " donors shared (RNA GSE173078 vs DNAm GSE173081)\n")
if(length(common)>=5){
  D <- as.matrix(dnam[match(common,dnam$sample), setdiff(colnames(dnam),"sample")])
  R <- rna_imm[common,,drop=FALSE]
  # DNAm within-immune (renormalize the immune subtypes)
  imm_cols <- intersect(c("B","NK","CD4T","CD8T","Mono","Neutro","Eoso"), colnames(D))
  Dimm <- D[,imm_cols,drop=FALSE]; Dimm <- Dimm/rowSums(Dimm)
  # map RNA subtypes -> DNAm subtypes
  Rmap <- cbind(B=R[,"Bcell"], NK=R[,"NK"], T=R[,"Tcell"], Mono=R[,"Myeloid"], Neutro=R[,"Neutrophil"])
  Dmap <- cbind(B=Dimm[,"B"], NK=Dimm[,"NK"], T=rowSums(Dimm[,intersect(c("CD4T","CD8T"),colnames(Dimm)),drop=FALSE]),
                Mono=Dimm[,"Mono"], Neutro=Dimm[,"Neutro"])
  cortab <- data.table(subtype=colnames(Rmap),
     Pearson=sapply(colnames(Rmap), function(k) suppressWarnings(cor(Rmap[,k], Dmap[,k]))),
     Spearman=sapply(colnames(Rmap), function(k) suppressWarnings(cor(Rmap[,k], Dmap[,k], method="spearman"))),
     RNA_mean=colMeans(Rmap), DNAm_mean=colMeans(Dmap))
  fwrite(cortab, file.path(OUT,"E3_RNA_vs_DNAm_immune_subtypes.tsv"), sep="\t")
  cat("\n=== RNA(E2层级) vs DNAm(HEpiDISH) 免疫亚型 within-immune 一致性 (同24donor) ===\n"); print(cortab)
}

## ---- compare to published histology (chronic periodontitis infiltrate) ----
# Infect Immun 1983: IgG+plasma ~65%, IgA ~11%, IgM ~1% (Ig-plasma ~77%); T ~29%
hist <- data.table(compartment=c("Plasma","Tcell","Bcell","Myeloid","Neutrophil"),
                   histology_pct=c("~65-77 (Ig plasma)","~29 (T lymph)","~18 (B)","~13 (mono/mac)","prominent junctional"))
# our within-immune means across array cohorts (GSE16134/GSE10334, higher immune)
arr <- rbindlist(lapply(c("GSE16134","GSE10334"), function(co){
  x<-fread(file.path(OUT,paste0(co,"_hierarchical.tsv"))); m<-as.matrix(x[,-1]); rownames(m)<-x$sample
  im<-m[,imm,drop=FALSE]; im<-im/rowSums(im); data.table(cohort=co, t(round(colMeans(im),3)))
}), fill=TRUE)
fwrite(arr, file.path(OUT,"E3_within_immune_array_cohorts.tsv"), sep="\t")
cat("\n=== within-immune 组成 (array 队列, 对照组织学: 浆细胞应主导) ===\n"); print(arr)
cat("\n组织学参照 (Infect Immun 1983): 浆~65-77%, T~29%, B~18%, 单核~13%\n")
message("[done] E3 anchoring")
