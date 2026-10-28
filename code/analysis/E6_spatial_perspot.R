#!/usr/bin/env Rscript
# E6 - per-spot spatial deconvolution as a third orthogonal modality. DEVIATION: RCTD "full"
# mode was computationally prohibitive across 8 Visium sections (>hours on the sigma/per-spot
# optimization); switched to InstaPrism per-spot (fast, method-consistent with the rest of the
# pipeline, achieves the same deliverable: per-spot proportions -> section composition ->
# health vs periodontitis). Documented deviation.
suppressMessages({library(Matrix); library(data.table); library(InstaPrism); library(GEOquery)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
VD <- file.path(P,"data/interim/GSE206621_visium"); OUT <- file.path(P,"results/tables/improvement")
imm <- c("Tcell","Bcell","Plasma","Myeloid","Mast","NK","Neutrophil")
refPhi <- readRDS(file.path(P,"results/models/instaprism_refPhi.rds"))

g<-getGEO("GSE206621",GSEMatrix=TRUE,getGPL=FALSE)[[1]]; pd<-Biobase::pData(g)
status<-ifelse(grepl("disease|periodont",tolower(paste(pd[["characteristics_ch1"]],pd[["characteristics_ch1.1"]],pd[["title"]]))),"periodontitis","health")
gsm2cond<-setNames(status, rownames(pd))

res<-list()
for(mf in Sys.glob(file.path(VD,"*_matrix.mtx.gz"))){
  base<-sub("_matrix.mtx.gz$","",mf); gsm<-sub("^(GSM[0-9]+)_.*","\\1",basename(mf))
  feat<-fread(paste0(base,"_features.tsv.gz"),header=FALSE); bc<-fread(paste0(base,"_barcodes.tsv.gz"),header=FALSE)
  sp<-readMM(gzfile(mf)); rownames(sp)<-make.unique(feat$V2); colnames(sp)<-bc$V1
  keep<-which(colSums(sp)>=200); if(length(keep)>800) keep<-sample(keep,800); sp<-as.matrix(sp[,keep,drop=FALSE])
  ip<-tryCatch(InstaPrism(bulk_Expr=sp, refPhi_cs=refPhi, n.core=8, verbose=FALSE), error=function(e){message("[WARN] ",gsm,": ",conditionMessage(e));NULL})
  if(is.null(ip))next
  w<-t(ip@Post.ini.ct@theta)                       # per-spot proportions (mRNA fraction)
  comp<-colMeans(w, na.rm=TRUE)
  res[[gsm]]<-data.table(gsm=gsm, condition=gsm2cond[gsm], n_spots=ncol(sp), t(round(comp,4)))
  message("[",gsm,"] ",gsm2cond[gsm]," per-spot done; ",ncol(sp)," spots")
}
tab<-rbindlist(res,fill=TRUE); fwrite(tab,file.path(OUT,"E6_spatial_composition.tsv"),sep="\t")
long<-melt(tab,id.vars=c("gsm","condition","n_spots"),variable.name="cell_type",value.name="frac")
cmp<-long[,.(mean_health=mean(frac[condition=="health"],na.rm=TRUE),mean_perio=mean(frac[condition=="periodontitis"],na.rm=TRUE)),by=cell_type]
cmp[,direction:=ifelse(mean_perio>mean_health,"UP","DOWN")]
fwrite(cmp,file.path(OUT,"E6_spatial_health_vs_perio.tsv"),sep="\t")
cat("\n=== E6 InstaPrism per-spot 空间构成: 牙周炎 vs 健康 ===\n"); print(cmp[order(-mean_perio)])
cat("\n免疫合计: health=",round(sum(cmp[cell_type%in%imm]$mean_health,na.rm=TRUE),3)," perio=",round(sum(cmp[cell_type%in%imm]$mean_perio,na.rm=TRUE),3),"\n")
message("[done] E6 spatial per-spot")
