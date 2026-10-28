#!/usr/bin/env Rscript
# Stage 7B - spatial (Visium) orthogonal check of composition. GSE206621 (gingiva Visium,
# 5 Health + 3 disease sections). Pseudobulk each section (sum spots) -> deconvolve with
# oral InstaPrism -> compare health vs periodontitis composition against the bulk findings.
suppressMessages({library(Matrix); library(data.table); library(GEOquery); library(InstaPrism)})
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
VD <- file.path(P,"data/interim/GSE206621_visium"); OUT <- file.path(P,"results/tables/validation")

## ---- section -> condition mapping from GEO ----
g <- getGEO("GSE206621", GSEMatrix=TRUE, getGPL=FALSE)[[1]]; pd <- Biobase::pData(g)
sec <- sub(".*\\[([A-Z0-9]+)\\].*","\\1", pd$title)           # e.g. A1, B1...
status <- ifelse(grepl("disease|periodontitis", tolower(paste(pd[["characteristics_ch1"]], pd[["characteristics_ch1.1"]]))), "periodontitis", "health")
gsm <- rownames(pd)
map <- data.table(gsm=gsm, section=sec, condition=status)
message("[spatial] sections: "); print(map)

## ---- load each Visium section, pseudobulk (sum spots) ----
mtx_files <- Sys.glob(file.path(VD,"*_matrix.mtx.gz"))
pb <- list()
for(mf in mtx_files){
  base <- sub("_matrix.mtx.gz$","",mf)
  gsm_i <- sub("^(GSM[0-9]+)_.*","\\1", basename(mf))
  feat <- fread(paste0(base,"_features.tsv.gz"), header=FALSE)
  m <- readMM(gzfile(mf))                       # genes x spots
  sym <- feat$V2; rn <- make.unique(sym)
  rownames(m) <- rn
  v <- Matrix::rowSums(m)                        # pseudobulk (sum across spots)
  pb[[gsm_i]] <- v
}
# align on common genes
allg <- Reduce(intersect, lapply(pb, names))
bulk <- sapply(pb, function(v) v[allg]); rownames(bulk) <- allg
message(sprintf("[spatial] pseudobulk %d genes x %d sections", nrow(bulk), ncol(bulk)))

## ---- deconvolve with oral InstaPrism ----
refPhi <- readRDS(file.path(P,"results/models/instaprism_refPhi.rds"))
ip <- InstaPrism(bulk_Expr=as.matrix(bulk), refPhi_cs=refPhi, n.core=8, verbose=FALSE)
fr <- t(ip@Post.ini.ct@theta); rownames(fr) <- colnames(bulk)
frdt <- as.data.table(fr, keep.rownames="gsm")
frdt <- merge(map[,.(gsm,section,condition)], frdt, by="gsm")
fwrite(frdt, file.path(OUT,"stage7b_spatial_composition.tsv"), sep="\t")
message("[spatial] per-section composition:"); print(frdt)

## ---- compare health vs periodontitis per cell type ----
cts <- setdiff(colnames(fr), character(0))
res <- rbindlist(lapply(cts, function(k){
  h <- fr[map$condition[match(rownames(fr),map$gsm)]=="health", k]
  d <- fr[map$condition[match(rownames(fr),map$gsm)]=="periodontitis", k]
  data.table(cell_type=k, mean_health=mean(h), mean_perio=mean(d),
             direction=ifelse(mean(d)>mean(h),"UP","DOWN"),
             wilcox_p=tryCatch(wilcox.test(d,h)$p.value, error=function(e) NA))
}))
fwrite(res, file.path(OUT,"stage7b_spatial_health_vs_perio.tsv"), sep="\t")
cat("\n=== Spatial composition: periodontitis vs health (direction) ===\n"); print(res[order(-mean_perio)])
message("[done] Stage 7B spatial validation")
