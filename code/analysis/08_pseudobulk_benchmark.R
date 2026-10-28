#!/usr/bin/env Rscript
# Stage 3 - leakage-controlled pseudobulk benchmark: oral reference vs LM22.
# Reference (oral) from GSE164241; pseudobulk ground-truth from INDEPENDENT GSE152042
# (different study/donors -> no leakage). FAIR comparison: oral signature and LM22 run
# through the SAME algorithm (nnls). Plus InstaPrism (oral full reference) as best-oral
# method, and the LM22 non-immune->immune misassignment test.
suppressMessages({library(Matrix); library(data.table)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
source(file.path(P,"code/utils/deconv_utils.R"))
TAB <- file.path(P,"results/tables"); BEN <- file.path(P,"results/benchmark"); dir.create(BEN,showWarnings=FALSE,recursive=TRUE)

ccc <- function(t,e){ 2*cov(t,e)/(var(t)+var(e)+(mean(t)-mean(e))^2) }
rmse <- function(t,e) sqrt(mean((t-e)^2))
score_tab <- function(truth, est, method){
  ct <- intersect(colnames(truth), colnames(est))
  do.call(rbind, lapply(ct, function(k){
    t<-truth[,k]; e<-est[rownames(truth),k]; e[is.na(e)]<-0
    data.frame(method=method, cell_type=k, CCC=ccc(t,e), RMSE=rmse(t,e),
               Pearson=suppressWarnings(cor(t,e)), bias=mean(e-t), mean_true=mean(t), mean_est=mean(e))
  }))
}
# LM22 (dotted names) -> coarse immune
lm22_map <- list(
  Tcell=c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"),
  NK=c("NK.cells.resting","NK.cells.activated"), Bcell=c("B.cells.naive","B.cells.memory"),
  Plasma=c("Plasma.cells"), Myeloid=c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated","Eosinophils"),
  Mast=c("Mast.cells.resting","Mast.cells.activated"), Neutrophil=c("Neutrophils"))
collapse_lm22 <- function(fr){ m<-sapply(names(lm22_map), function(g){cc<-intersect(lm22_map[[g]],colnames(fr)); if(length(cc)) rowSums(fr[,cc,drop=FALSE]) else rep(0,nrow(fr))}); rownames(m)<-rownames(fr); m }

## ---- independent pseudobulk source (GSE152042 labeled) ----
SRC <- file.path(P,"data/processed/pseudobulk_source_10x")
m <- readMM(gzfile(file.path(SRC,"matrix.mtx.gz"))); genes<-readLines(gzfile(file.path(SRC,"features.tsv.gz")))
cells<-readLines(gzfile(file.path(SRC,"barcodes.tsv.gz"))); meta<-fread(file.path(SRC,"cell_meta.tsv.gz"))
m<-as(m,"CsparseMatrix"); rownames(m)<-genes; colnames(m)<-cells
tt<-table(meta$cell_type); keep_ct<-names(tt)[tt>=50]; sel<-meta$cell_type %in% keep_ct
m<-m[,sel]; meta<-meta[sel]
message(sprintf("[pseudobulk source] %d genes x %d cells; types: %s", nrow(m), ncol(m), paste(keep_ct,collapse=",")))

## ---- SimBu pseudobulk ----
suppressMessages(library(SimBu))
ds <- SimBu::dataset(annotation=data.frame(ID=colnames(m), cell_type=meta$cell_type), count_matrix=m, name="gingivaB", filter_genes=TRUE)
scen <- list(random=c("random","NONE"), mirror=c("mirror_db","NONE"), census=c("mirror_db","census"))
sims <- list()
for(nm in names(scen)){
  sims[[nm]] <- SimBu::simulate_bulk(ds, scenario=scen[[nm]][1], scaling_factor=scen[[nm]][2], ncells=1000, nsamples=100, seed=match(nm,names(scen)), BPPARAM=BiocParallel::MulticoreParam(8), run_parallel=TRUE)
  message("[SimBu] ",nm," done")
}
saveRDS(sims, file.path(BEN,"simbu_sims.rds"))

## ---- signatures ----
S_gep <- readRDS(file.path(P,"results/models/S_oral_DWLS.rds"))
LM22  <- readRDS(file.path(P,"results/models/LM22.rds"))

## ---- deconvolve: oral-sig vs LM22 through the SAME algorithm ----
## primary = DWLS solveDampenedWLS (reused published solver); secondary = nnls
all_scores<-list(); misassign<-list()
nonimm_types <- c("Epithelial","Fibroblast","Endothelial","Pericyte","Melanocyte","Glial")
for(nm in names(sims)){
  s<-sims[[nm]]; truth<-as.matrix(s$cell_fractions)
  bulk<-as.matrix(SummarizedExperiment::assays(s$bulk)[["bulk_counts"]])   # SimBu provides counts (no TPM supplied)
  for(eng in c("nnls")){
    fr_oral <- deconv_nnls(bulk, S_gep)
    fr_lm22 <- deconv_nnls(bulk, LM22)
    prc <- collapse_lm22(fr_lm22)
    all_scores[[paste(nm,"oral",eng)]] <- cbind(scenario=nm, score_tab(truth, fr_oral, paste0("oral_",eng)))
    all_scores[[paste(nm,"LM22",eng)]] <- cbind(scenario=nm, score_tab(truth, prc, paste0("LM22_",eng)))
    if(eng=="nnls"){
      nonimm <- rowSums(truth[, intersect(nonimm_types, colnames(truth)), drop=FALSE])
      mis <- sapply(colnames(prc), function(k) suppressWarnings(cor(nonimm, prc[rownames(truth),k])))
      misassign[[nm]] <- data.frame(scenario=nm, method="LM22_nnls", immune_type=names(mis), cor_with_true_nonimmune=as.numeric(mis))
    }
  }
}
scores_df<-rbindlist(lapply(all_scores,as.data.table),fill=TRUE); fwrite(scores_df, file.path(TAB,"stage3_benchmark_scores.tsv"), sep="\t")
mis_df<-rbindlist(misassign,fill=TRUE); fwrite(mis_df, file.path(TAB,"stage3_lm22_misassignment.tsv"), sep="\t")

## ---- InstaPrism (oral full reference) as best-oral method ----
tryCatch({
  library(InstaPrism); refPhi<-readRDS(file.path(P,"results/models/instaprism_refPhi.rds"))
  ip_scores<-list()
  for(nm in names(sims)){
    s<-sims[[nm]]; truth<-as.matrix(s$cell_fractions); bulk<-as.matrix(SummarizedExperiment::assays(s$bulk)[["bulk_counts"]])
    ip<-InstaPrism(bulk_Expr=bulk, refPhi_cs=refPhi, n.core=8, verbose=FALSE)
    fr<-t(ip@Post.ini.ct@theta); rownames(fr)<-colnames(bulk)
    ip_scores[[nm]]<-cbind(scenario=nm, score_tab(truth, fr, "InstaPrism_oral"))
  }
  fwrite(rbindlist(lapply(ip_scores,as.data.table),fill=TRUE), file.path(TAB,"stage3_instaprism_oral_scores.tsv"), sep="\t")
  message("[saved] InstaPrism oral scores")
}, error=function(e) message("[InstaPrism][WARN] ",conditionMessage(e)))

## ---- summary ----
sink(file.path(P,"logs/qc/stage3_benchmark_summary.txt"))
cat("Stage3 pseudobulk benchmark (oral vs LM22; same-algorithm nnls)\n\n")
cat("kappa: oral_signature=",kappa(as.matrix(S_gep))," LM22=",kappa(as.matrix(LM22)),"\n\n")
print(scores_df[, .(meanCCC=round(mean(CCC,na.rm=TRUE),3), meanRMSE=round(mean(RMSE,na.rm=TRUE),4)), by=.(scenario,method)])
cat("\n-- LM22 non-immune->immune misassignment (positive corr = spillover) --\n")
print(mis_df[order(-cor_with_true_nonimmune)][1:min(20,nrow(mis_df))])
sink()
message("[done] Stage 3 benchmark")
