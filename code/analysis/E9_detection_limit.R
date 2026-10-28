#!/usr/bin/env Rscript
# E9 - rare-cell-type detection-limit for the deconvolution. On leakage-controlled pseudobulk
# (GSE152042, independent), titrate each cell type's TRUE fraction 1%->40% and record the
# InstaPrism estimate + recovery, to locate the abundance floor and distinguish legitimately
# near-zero cell types from artifactually-zeroed ones.
suppressMessages({library(Matrix); library(data.table); library(InstaPrism)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
OUT <- file.path(P,"results/tables/improvement"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
refPhi <- readRDS(file.path(P,"results/models/instaprism_refPhi.rds"))
# per-cell-type mean CPM profiles from the INDEPENDENT GSE152042 (pseudobulk source)
SRC<-file.path(P,"data/processed/pseudobulk_source_10x")
mm<-readMM(gzfile(file.path(SRC,"matrix.mtx.gz"))); g<-readLines(gzfile(file.path(SRC,"features.tsv.gz")))
meta<-fread(file.path(SRC,"cell_meta.tsv.gz")); mm<-as(mm,"CsparseMatrix"); rownames(mm)<-g
cpm<-t(t(mm)/Matrix::colSums(mm))*1e6
types<-sort(unique(meta$cell_type))
mu<-sapply(types,function(k)Matrix::rowMeans(cpm[,meta$cell_type==k,drop=FALSE])); mu<-as.matrix(mu); rownames(mu)<-g

res<-list()
for(target in types){
  fr_levels<-c(0.01,0.02,0.05,0.10,0.20,0.40)
  ests<-c(); trues<-c()
  bulks<-sapply(fr_levels,function(fr){
    p<-setNames(rep((1-fr)/(length(types)-1),length(types)),types); p[target]<-fr
    as.numeric(mu%*%p)*exp(rnorm(nrow(mu),0,0.1))
  }); rownames(bulks)<-g; colnames(bulks)<-paste0("fr",fr_levels)
  ip<-tryCatch(InstaPrism(bulk_Expr=bulks,refPhi_cs=refPhi,n.core=8,verbose=FALSE),error=function(e)NULL)
  if(is.null(ip))next
  th<-t(ip@Post.ini.ct@theta)
  if(!target%in%colnames(th))next
  est<-th[,target]
  res[[target]]<-data.table(cell_type=target, true_frac=fr_levels, est_frac=as.numeric(est))
}
dt<-rbindlist(res,fill=TRUE); fwrite(dt,file.path(OUT,"E9_detection_limit.tsv"),sep="\t")
# floor = smallest true_frac where est >= 0.5*true (rough recovery threshold)
floor_tab<-dt[,.(floor_frac={ok<-true_frac[est_frac>=0.5*true_frac]; if(length(ok))min(ok) else NA_real_},
                 est_at_5pct=est_frac[true_frac==0.05], est_at_20pct=est_frac[true_frac==0.20]),by=cell_type]
fwrite(floor_tab,file.path(OUT,"E9_detection_floor.tsv"),sep="\t")
cat("\n=== 检出限: 各类型 true 5%/20% 时的估计 + 恢复底线 ===\n"); print(floor_tab[order(floor_frac)])
message("[done] E9 detection limit")
