#!/usr/bin/env Rscript
# Stage 2b - build deconvolution reference artifacts from the annotated oral reference:
#  (1) fast-reload reference .rds (genes x cells counts + labels)
#  (2) InstaPrism refPhi (per cell.type & cell.state)
#  (3) DWLS signature matrix S_oral (subsampled for tractable MAST) + condition number
#  (4) per-cell-type mean-CPM GEP signature (robust backup)
#  (5) LM22 (ADAPTS) + condition number, for the fair comparison in Stage 3
suppressMessages({library(Matrix); library(data.table)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
REF <- file.path(P,"data/processed/reference_10x")
OUT <- file.path(P,"results/models"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
TAB <- file.path(P,"results/tables")

message("[load] 10x reference")
m <- readMM(gzfile(file.path(REF,"matrix.mtx.gz")))      # genes x cells (integer)
genes <- readLines(gzfile(file.path(REF,"features.tsv.gz")))
cells <- readLines(gzfile(file.path(REF,"barcodes.tsv.gz")))
meta  <- fread(file.path(REF,"cell_meta.tsv.gz"))
m <- as(m, "CsparseMatrix"); rownames(m) <- genes; colnames(m) <- cells
stopifnot(nrow(meta)==ncol(m))
ct <- meta$cell_type; cs <- meta$cell_state
message(sprintf("  %d genes x %d cells; %d cell types", nrow(m), ncol(m), length(unique(ct))))

# (1) fast-reload rds
saveRDS(list(counts=m, cell_type=ct, cell_state=cs, meta=meta),
        file.path(OUT,"oral_reference_counts.rds"))
message("[save] oral_reference_counts.rds")

# (4) per-cell-type mean-CPM GEP signature (fast, robust)
cpm <- t(t(m)/Matrix::colSums(m))*1e6
types <- sort(unique(ct))
gep <- sapply(types, function(k) Matrix::rowMeans(cpm[, ct==k, drop=FALSE]))
colnames(gep) <- types
# marker-gene subset: top DE genes per type by log fold change vs rest (for a compact signature)
logmean <- log1p(gep)
signif_genes <- unique(unlist(lapply(types, function(k){
  fc <- logmean[,k] - apply(logmean[,setdiff(types,k),drop=FALSE],1,max)
  names(sort(fc, decreasing=TRUE))[1:100]
})))
signif_genes <- intersect(signif_genes, rownames(gep))
S_gep <- gep[signif_genes,,drop=FALSE]
write.table(data.frame(gene=rownames(S_gep), S_gep, check.names=FALSE),
            file.path(TAB,"stage2_S_oral_gepMeanCPM.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
saveRDS(S_gep, file.path(OUT,"S_oral_gepMeanCPM.rds"))
message(sprintf("[save] mean-CPM signature: %d genes x %d types; kappa=%.1f",
                nrow(S_gep), ncol(S_gep), kappa(as.matrix(S_gep))))

# (2) InstaPrism refPhi
ok_instaprism <- requireNamespace("InstaPrism", quietly=TRUE)
if (ok_instaprism){
  tryCatch({
    library(InstaPrism)
    refPhi <- refPrepare(sc_Expr=as.matrix(m), cell.type.labels=ct, cell.state.labels=cs)
    saveRDS(refPhi, file.path(OUT,"instaprism_refPhi.rds"))
    message("[save] instaprism_refPhi.rds")
  }, error=function(e) message("[InstaPrism][WARN] ", conditionMessage(e)))
} else message("[InstaPrism] not installed - skip (will use BayesPrism directly)")

# (3) DWLS signature (subsample for tractable MAST)
ok_dwls <- requireNamespace("DWLS", quietly=TRUE)
if (ok_dwls){
  tryCatch({
    library(DWLS)
    cap <- 600
    idx <- unlist(lapply(types, function(k){ w<-which(ct==k); if(length(w)>cap) sample(w,cap) else w }))
    sc_sub <- as.matrix(m[, idx]); id_sub <- ct[idx]
    # drop all-zero genes to speed MAST
    keep <- Matrix::rowSums(sc_sub)>0; sc_sub <- sc_sub[keep,]
    dwls_path <- file.path(OUT,"dwls_work"); dir.create(dwls_path, showWarnings=FALSE)
    message(sprintf("[DWLS] MAST on %d cells x %d genes (cap=%d/type)", ncol(sc_sub), nrow(sc_sub), cap))
    Sig <- buildSignatureMatrixMAST(scdata=sc_sub, id=id_sub, path=dwls_path,
                                    diff.cutoff=0.5, pval.cutoff=0.01)
    saveRDS(Sig, file.path(OUT,"S_oral_DWLS.rds"))
    write.table(data.frame(gene=rownames(Sig), Sig, check.names=FALSE),
                file.path(TAB,"stage2_S_oral_DWLS.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    message(sprintf("[save] DWLS signature: %d genes x %d types; kappa=%.1f",
                    nrow(Sig), ncol(Sig), kappa(as.matrix(Sig))))
  }, error=function(e) message("[DWLS][WARN] ", conditionMessage(e)))
} else message("[DWLS] not installed - skip")

# (5) LM22 from ADAPTS
ok_adapts <- requireNamespace("ADAPTS", quietly=TRUE)
if (ok_adapts){
  tryCatch({
    data("LM22", package="ADAPTS")
    LM22 <- get("LM22")
    write.table(data.frame(gene=rownames(LM22), LM22, check.names=FALSE),
                file.path(TAB,"LM22.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    saveRDS(LM22, file.path(OUT,"LM22.rds"))
    message(sprintf("[save] LM22: %d genes x %d types; kappa=%.1f",
                    nrow(LM22), ncol(LM22), kappa(as.matrix(LM22))))
  }, error=function(e) message("[ADAPTS/LM22][WARN] ", conditionMessage(e)))
} else message("[ADAPTS] not installed - skip")

# summary
sink(file.path(P,"logs/qc/stage2_reference_artifacts_summary.txt"))
cat("Stage2 reference artifacts\n")
cat(sprintf("reference: %d genes x %d cells; %d cell types: %s\n", nrow(m), ncol(m),
            length(types), paste(types, collapse=", ")))
cat(sprintf("mean-CPM signature genes=%d kappa=%.1f\n", nrow(S_gep), kappa(as.matrix(S_gep))))
if (exists("Sig")) cat(sprintf("DWLS signature genes=%d kappa=%.1f\n", nrow(Sig), kappa(as.matrix(Sig))))
if (exists("LM22")) cat(sprintf("LM22 genes=%d kappa=%.1f\n", nrow(LM22), kappa(as.matrix(LM22))))
sink()
message("[done] Stage 2 reference artifacts")
