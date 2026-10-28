#!/usr/bin/env Rscript
# Stage 4b - deconvolve each bulk cohort with (i) oral InstaPrism (full ref),
# (ii) oral mean-CPM signature via nnls, (iii) LM22 via nnls (fair: same algorithm).
# Then test the COMPOSITION channel (proportion ~ disease) with logit+limma and CLR+limma.
suppressMessages({library(data.table); library(Matrix); library(limma)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
source(file.path(P,"code/utils/deconv_utils.R"))
BULK <- file.path(P,"data/processed/bulk"); TAB <- file.path(P,"results/tables")
DEC <- file.path(P,"results/tables/proportions"); dir.create(DEC,showWarnings=FALSE,recursive=TRUE)

refPhi <- readRDS(file.path(P,"results/models/instaprism_refPhi.rds"))
S_gep  <- readRDS(file.path(P,"results/models/S_oral_DWLS.rds"))
LM22   <- readRDS(file.path(P,"results/models/LM22.rds"))
lm22_map <- list(
  Tcell=c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"),
  NK=c("NK.cells.resting","NK.cells.activated"), Bcell=c("B.cells.naive","B.cells.memory"),
  Plasma=c("Plasma.cells"), Myeloid=c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated","Eosinophils"),
  Mast=c("Mast.cells.resting","Mast.cells.activated"), Neutrophil=c("Neutrophils"))
collapse_lm22 <- function(fr){ m<-sapply(names(lm22_map), function(g){cc<-intersect(lm22_map[[g]],colnames(fr)); if(length(cc)) rowSums(fr[,cc,drop=FALSE]) else rep(0,nrow(fr))}); rownames(m)<-rownames(fr); m }

logit <- function(p, off=0.01){ p<-pmin(pmax(p,off),1-off); log(p/(1-p)) }
clr   <- function(P, eps=1e-6){ P<-P+eps; L<-log(P); sweep(L,1,rowMeans(L)) }
comp_test <- function(props, group, method_tag, cohort){
  grp <- factor(group, levels=c("health","periodontitis"))
  if(nlevels(droplevels(grp))<2 || any(table(grp)<3)) return(NULL)
  des <- model.matrix(~grp); out<-list()
  Y<-t(logit(as.matrix(props))); fit<-eBayes(lmFit(Y,des)); tt<-topTable(fit,coef=2,number=Inf,sort.by="none")
  out[["logit"]]<-data.frame(cohort=cohort,deconv=method_tag,test="logit_limma",cell_type=rownames(Y),effect=tt$logFC,p=tt$P.Value,FDR=tt$adj.P.Val)
  Yc<-t(clr(as.matrix(props))); fit2<-eBayes(lmFit(Yc,des)); tt2<-topTable(fit2,coef=2,number=Inf,sort.by="none")
  out[["clr"]]<-data.frame(cohort=cohort,deconv=method_tag,test="clr_limma",cell_type=rownames(Yc),effect=tt2$logFC,p=tt2$P.Value,FDR=tt2$adj.P.Val)
  rbindlist(out)
}

cohorts <- sub("_linear.rds","", basename(Sys.glob(file.path(BULK,"*_linear.rds"))))
message("[cohorts] ", paste(cohorts, collapse=", "))
all_comp <- list()
for(co in cohorts){
  d<-readRDS(file.path(BULK,paste0(co,"_linear.rds"))); expr<-d$expr; meta<-d$meta
  meta<-meta[meta$sample %in% colnames(expr) & meta$group %in% c("health","periodontitis"),]
  if(nrow(meta)<4 || any(table(meta$group)<2)){ message("[",co,"] too few samples; skip"); next }
  expr<-as.matrix(expr[,meta$sample,drop=FALSE])
  message(sprintf("[%s] %d genes x %d samples (%s)", co, nrow(expr), ncol(expr), paste(names(table(meta$group)),table(meta$group),sep="=",collapse=",")))
  ## (i) oral InstaPrism full ref (primary oral method)
  props_ip<-tryCatch({ library(InstaPrism); ip<-InstaPrism(bulk_Expr=expr,refPhi_cs=refPhi,n.core=8,verbose=FALSE); fr<-t(ip@Post.ini.ct@theta); rownames(fr)<-colnames(expr); fr },
                     error=function(e){message("[InstaPrism][WARN] ",co,": ",conditionMessage(e)); NULL})
  if(!is.null(props_ip)){ fwrite(data.table(sample=rownames(props_ip),props_ip), file.path(DEC,paste0(co,"_oral_InstaPrism.tsv")),sep="\t")
    all_comp[[paste(co,"oral_InstaPrism")]]<-comp_test(props_ip, meta$group, "oral_InstaPrism", co) }
  ## (ii) oral signature via DWLS solveDampenedWLS (reused) - fair comparison
  fr_oral<-deconv_dwls(expr, S_gep, "dampened", cores=8)
  fwrite(data.table(sample=rownames(fr_oral),fr_oral), file.path(DEC,paste0(co,"_oral_DWLS.tsv")),sep="\t")
  all_comp[[paste(co,"oral_DWLS")]]<-comp_test(as.data.frame(fr_oral), meta$group, "oral_DWLS", co)
  ## (iii) LM22 via DWLS solveDampenedWLS (SAME algorithm)
  fr_lm22<-deconv_dwls(expr, LM22, "dampened", cores=8); prc<-collapse_lm22(fr_lm22); prc<-prc/rowSums(prc)
  fwrite(data.table(sample=rownames(prc),prc), file.path(DEC,paste0(co,"_LM22_DWLS.tsv")),sep="\t")
  all_comp[[paste(co,"LM22_DWLS")]]<-comp_test(as.data.frame(prc), meta$group, "LM22_DWLS", co)
}
comp_df<-rbindlist(all_comp,fill=TRUE)
comp_df$sig<-ifelse(comp_df$FDR<0.05, ifelse(comp_df$effect>0,"UP","DOWN"), "ns")
fwrite(comp_df, file.path(TAB,"stage4_composition_tests.tsv"), sep="\t")
message("[done] Stage 4b: proportions + composition tests")
s<-comp_df[test=="logit_limma", .N, by=.(cohort,deconv,sig)]; print(dcast(s, cohort+deconv~sig, value.var="N", fill=0))
