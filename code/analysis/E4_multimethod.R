#!/usr/bin/env Rscript
# E4 - multi-method proportion consensus: is the flip-catalog / oral composition direction
# method-invariant (not an InstaPrism artifact)? Run InstaPrism + MuSiC + nnls(oral sig) +
# nnls(LM22) on each cohort, test composition direction (logit+limma), score agreement.
suppressMessages({library(Matrix); library(data.table); library(limma); library(SingleCellExperiment); library(MuSiC); library(InstaPrism)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
source(file.path(P,"code/utils/deconv_utils.R"))
OUT <- file.path(P,"results/tables/improvement"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
refPhi <- readRDS(file.path(P,"results/models/instaprism_refPhi.rds"))
S_gep  <- readRDS(file.path(P,"results/models/S_oral_DWLS.rds"))
LM22   <- readRDS(file.path(P,"results/models/LM22.rds"))
lm22_map <- list(Tcell=c("T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting","T.cells.CD4.memory.activated","T.cells.follicular.helper","T.cells.regulatory..Tregs.","T.cells.gamma.delta"),
  NK=c("NK.cells.resting","NK.cells.activated"),Bcell=c("B.cells.naive","B.cells.memory"),Plasma=c("Plasma.cells"),
  Myeloid=c("Monocytes","Macrophages.M0","Macrophages.M1","Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated","Eosinophils"),
  Mast=c("Mast.cells.resting","Mast.cells.activated"),Neutrophil=c("Neutrophils"))
collapse_lm22<-function(fr){m<-sapply(names(lm22_map),function(g){cc<-intersect(lm22_map[[g]],colnames(fr));if(length(cc))rowSums(fr[,cc,drop=FALSE])else rep(0,nrow(fr))});rownames(m)<-rownames(fr);m}
logit<-function(p,o=0.01){p<-pmin(pmax(p,o),1-o);log(p/(1-p))}
comp_dir<-function(props,grp){grp<-factor(grp,levels=c("health","periodontitis"));if(any(table(grp)<3))return(NULL)
  des<-model.matrix(~grp);fit<-eBayes(lmFit(t(logit(as.matrix(props))),des));tt<-topTable(fit,coef=2,number=Inf,sort.by="none")
  setNames(ifelse(tt$adj.P.Val<0.05,ifelse(tt$logFC>0,"UP","DOWN"),"ns"),rownames(tt))}

## SCE for MuSiC (cap 400/type)
ref<-readRDS(file.path(P,"results/models/oral_reference_counts.rds")); m<-ref$counts; ct<-ref$cell_type; samp<-ref$meta$sample
cap<-400; idx<-unlist(lapply(unique(ct),function(k){w<-which(ct==k);if(length(w)>cap)sample(w,cap)else w}))
sce<-SingleCellExperiment(assays=list(counts=as.matrix(m[,idx])),colData=DataFrame(cell_type=ct[idx],sampleID=samp[idx]))

dirs<-list()
for(co in c("GSE16134","GSE10334","GSE223924")){
  d<-readRDS(file.path(P,"data/processed/bulk",paste0(co,"_linear.rds"))); bulk<-as.matrix(d$expr); meta<-d$meta
  meta<-meta[meta$sample%in%colnames(bulk)&meta$group%in%c("health","periodontitis"),]; bulk<-bulk[,meta$sample]; grp<-meta$group
  ## InstaPrism
  ip<-InstaPrism(bulk_Expr=bulk,refPhi_cs=refPhi,n.core=8,verbose=FALSE); pi<-t(ip@Post.ini.ct@theta); rownames(pi)<-colnames(bulk)
  ## MuSiC
  pm<-tryCatch(music_prop(bulk.mtx=bulk,sc.sce=sce,clusters="cell_type",samples="sampleID",select.ct=unique(ct),verbose=FALSE)$Est.prop.weighted,error=function(e)NULL)
  ## nnls oral
  pn<-deconv_nnls(bulk,S_gep)
  ## LM22
  pl<-collapse_lm22(deconv_nnls(bulk,LM22)); pl<-pl/rowSums(pl)
  for(nm in c("InstaPrism","MuSiC","nnls_oral","LM22")){
    pr<-switch(nm,InstaPrism=pi,MuSiC=pm,nnls_oral=pn,LM22=pl); if(is.null(pr))next
    dd<-comp_dir(as.data.frame(pr),grp); if(is.null(dd))next
    dirs[[paste(co,nm)]]<-data.table(cohort=co,method=nm,cell_type=names(dd),dir=as.character(dd))
  }
  message("[",co,"] multi-method done")
}
dt<-rbindlist(dirs,fill=TRUE); fwrite(dt,file.path(OUT,"E4_multimethod_directions.tsv"),sep="\t")
## consensus among ORAL methods per cohort x cell_type
oral<-dt[method%in%c("InstaPrism","MuSiC","nnls_oral")]
cons<-oral[,.(n=.N, up=sum(dir=="UP"), down=sum(dir=="DOWN"), ns=sum(dir=="ns"),
    consensus=names(which.max(table(dir)))),by=.(cohort,cell_type)]
fwrite(cons,file.path(OUT,"E4_oral_consensus.tsv"),sep="\t")
cat("\n=== 口腔多方法(InstaPrism/MuSiC/nnls)构成方向共识 ===\n"); print(cons[order(cohort,cell_type)])
cat("\n=== LM22 vs 口腔共识 翻转核对 (GSE16134) ===\n")
w<-dcast(dt[cohort=="GSE16134"&cell_type%in%c("Tcell","Bcell","Plasma","Myeloid","Mast","NK","Neutrophil")],cell_type~method,value.var="dir")
print(w)
message("[done] E4 multimethod")
