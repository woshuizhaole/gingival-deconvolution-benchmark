#!/usr/bin/env Rscript
# Stage 2c - build the oral signature matrix using DWLS's OWN published method
# (buildSignatureMatrixMAST: MAST DE + condition-number minimisation). Reuses the
# DWLS implementation rather than a hand-made mean-CPM signature. Subsamples cells
# per type for tractable MAST. Replaces S_oral for the fair LM22 comparison.
suppressMessages({library(DWLS); library(Matrix)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
ref <- readRDS(file.path(P,"results/models/oral_reference_counts.rds"))
m <- ref$counts; ct <- ref$cell_type        # genes x cells, labels
types <- sort(unique(ct))
cap <- 400
idx <- unlist(lapply(types, function(k){ w<-which(ct==k); if(length(w)>cap) sample(w,cap) else w }))
sc <- as.matrix(m[, idx]); id <- ct[idx]
keep <- Matrix::rowSums(sc) > 0; sc <- sc[keep,]
message(sprintf("[DWLS build] MAST on %d cells x %d genes, %d types (cap=%d/type)", ncol(sc), nrow(sc), length(types), cap))
path <- file.path(P,"results/models/dwls_work"); dir.create(path, showWarnings=FALSE, recursive=TRUE)
Sig <- DWLS::buildSignatureMatrixMAST(scdata=sc, id=id, path=path, diff.cutoff=0.5, pval.cutoff=0.01)
saveRDS(Sig, file.path(P,"results/models/S_oral_DWLS.rds"))
write.table(data.frame(gene=rownames(Sig), Sig, check.names=FALSE),
            file.path(P,"results/tables/stage2_S_oral_DWLS.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
message(sprintf("[done] DWLS signature: %d genes x %d types; kappa=%.1f", nrow(Sig), ncol(Sig), kappa(as.matrix(Sig))))
# quick self-test: deconvolve it back should recover something sensible
source(file.path(P,"code/utils/deconv_utils.R"))
d <- readRDS(file.path(P,"data/processed/bulk/GSE173078_linear.rds"))
fr <- deconv_dwls(as.matrix(d$expr[,1:3]), Sig, "dampened", cores=3)
cat("self-test frac (rowSums", paste(round(rowSums(fr),3),collapse=","), "):\n"); print(round(fr,3))
