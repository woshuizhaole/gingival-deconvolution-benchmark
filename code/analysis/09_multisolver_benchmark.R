#!/usr/bin/env Rscript
# Stage 3b - multi-solver benchmark on the SAME SimBu mixtures as stage 3.
# Adds MuSiC, EPIC, DWLS(dampened), DWLS(SVR ~= CIBERSORT nu-SVR proxy), BayesPrism.
# CIBERSORTx-absolute needs the Stanford token/Docker -> flagged, SVR used as proxy.
suppressMessages({library(Matrix); library(data.table)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
source(file.path(P,"code/utils/deconv_utils.R"))
TAB <- file.path(P,"results/tables")

ccc  <- function(t,e) 2*cov(t,e)/(var(t)+var(e)+(mean(t)-mean(e))^2)
rmse <- function(t,e) sqrt(mean((t-e)^2))
score_tab <- function(truth, est, method){
  ct <- intersect(colnames(truth), colnames(est))
  do.call(rbind, lapply(ct, function(k){
    t<-truth[,k]; e<-est[rownames(truth),k]; e[is.na(e)]<-0
    data.frame(method=method, cell_type=k, CCC=ccc(t,e), RMSE=rmse(t,e),
               Pearson=suppressWarnings(cor(t,e)), bias=mean(e-t),
               mean_true=mean(t), mean_est=mean(e))
  }))
}
norm_rows <- function(fr){ fr[fr<0]<-0; s<-rowSums(fr); s[s==0]<-1; fr/s }

sims <- readRDS(file.path(P,"results/benchmark/simbu_sims.rds"))
ref  <- readRDS(file.path(P,"results/models/oral_reference_counts.rds"))   # list(counts,cell_type,cell_state,meta)
S_gep <- readRDS(file.path(P,"results/models/S_oral_DWLS.rds"))            # signature genes x types
scen <- names(sims)
get_bulk  <- function(nm) as.matrix(SummarizedExperiment::assays(sims[[nm]]$bulk)[["bulk_counts"]])  # genes x samples
get_truth <- function(nm) as.matrix(sims[[nm]]$cell_fractions)                                       # samples x types

all_scores <- list(); OUT <- file.path(TAB,"stage3_multisolver_scores.tsv")
save_now <- function(){ fwrite(rbindlist(lapply(all_scores,as.data.table),fill=TRUE), OUT, sep="\t"); message("[saved] ",length(all_scores)," blocks -> ",OUT) }

## ---------- EPIC (mRNA-content aware) ----------
tryCatch({
  library(EPIC)
  refP <- as.matrix(S_gep)
  for(nm in scen){
    bulk<-get_bulk(nm); truth<-get_truth(nm)
    g<-intersect(rownames(bulk),rownames(refP))
    res<-EPIC::EPIC(bulk=bulk[g,], reference=list(refProfiles=refP[g,], sigGenes=g), withOtherCells=FALSE, scaleExprs=TRUE)
    fr<-norm_rows(res$cellFractions); rownames(fr)<-colnames(bulk)
    all_scores[[paste(nm,"EPIC")]]<-cbind(scenario=nm, score_tab(truth, fr, "EPIC"))
    message("[EPIC] ",nm," done")
  }
  save_now()
}, error=function(e) message("[EPIC][WARN] ",conditionMessage(e)))

## ---------- MuSiC (cell-size weighted) ----------
tryCatch({
  library(MuSiC); library(SingleCellExperiment)
  cm<-ref$counts; md<-as.data.frame(ref$meta)
  sce<-SingleCellExperiment(assays=list(counts=cm), colData=S4Vectors::DataFrame(md))
  for(nm in scen){
    bulk<-get_bulk(nm); truth<-get_truth(nm)
    est<-MuSiC::music_prop(bulk.mtx=bulk, sc.sce=sce, clusters="cell_type", samples="donor", verbose=FALSE)
    fr<-norm_rows(as.matrix(est$Est.prop.weighted));
    all_scores[[paste(nm,"MuSiC")]]<-cbind(scenario=nm, score_tab(truth, fr, "MuSiC"))
    message("[MuSiC] ",nm," done")
  }
  save_now()
}, error=function(e) message("[MuSiC][WARN] ",conditionMessage(e)))

## ---------- DWLS dampened + SVR(CIBERSORT nu-SVR proxy) ----------
for(mm in c("dampened","svr")){
  tryCatch({
    for(nm in scen){
      bulk<-get_bulk(nm); truth<-get_truth(nm)
      fr<-deconv_dwls(bulk, S_gep, method=mm, cores=12); fr<-norm_rows(as.matrix(fr))
      lab<-if(mm=="svr") "SVR_CIBERSORTproxy" else "DWLS_dampened"
      all_scores[[paste(nm,lab)]]<-cbind(scenario=nm, score_tab(truth, fr, lab))
      message("[",lab,"] ",nm," done")
    }
    save_now()
  }, error=function(e) message("[DWLS-",mm,"][WARN] ",conditionMessage(e)))
}

## ---------- BayesPrism (explicit request; InstaPrism = its fast posterior) ----------
tryCatch({
  suppressMessages(library(BayesPrism))
  scC<-t(as.matrix(ref$counts))   # cells x genes
  for(nm in c("random")){          # heavy: random scenario (the one Table 1 uses)
    bulk<-get_bulk(nm); truth<-get_truth(nm)
    pr<-new.prism(reference=scC, mixture=t(bulk), input.type="count.matrix",
                  cell.type.labels=ref$cell_type, cell.state.labels=ref$cell_state,
                  key=NULL, outlier.cut=0.01, outlier.fraction=0.1)
    bp<-run.prism(pr, n.cores=16)
    theta<-get.fraction(bp, which.theta="final", state.or.type="type")
    fr<-norm_rows(as.matrix(theta)); rownames(fr)<-colnames(bulk)
    all_scores[[paste(nm,"BayesPrism")]]<-cbind(scenario=nm, score_tab(truth, fr, "BayesPrism"))
    message("[BayesPrism] ",nm," done")
  }
  save_now()
}, error=function(e) message("[BayesPrism][WARN] ",conditionMessage(e)))

message("=== multisolver benchmark complete ===")
