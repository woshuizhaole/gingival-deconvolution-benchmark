#!/usr/bin/env Rscript
# E1 - systematic mRNA-content / transcriptome-size correction of deconvolution proportions.
# BayesPrism/InstaPrism theta = mRNA(read) fraction, NOT cell fraction; high-mRNA cells
# (epithelial ~4.5x T-cell) are inflated. Correct via per-cell-type size factors; cross-check
# with MuSiC (internal cell_size) on RNA-seq cohorts; validate immune-total vs DNAm (~42%).
suppressMessages({library(Matrix); library(data.table); library(SingleCellExperiment); library(MuSiC)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
OUT <- file.path(P,"results/tables/improvement"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
imm <- c("Tcell","Bcell","Plasma","Myeloid","Mast","NK","Neutrophil")

ref <- readRDS(file.path(P,"results/models/oral_reference_counts.rds"))
m <- ref$counts; ct <- ref$cell_type; samp <- ref$meta$sample
## per-cell-type size factors (mean total mRNA/cell)  = the reusable reference asset
sizef <- tapply(Matrix::colSums(m), ct, mean)
saveRDS(sizef, file.path(P,"results/models/celltype_mRNA_sizefactors.rds"))
fwrite(data.table(cell_type=names(sizef), mRNA_per_cell=round(as.numeric(sizef))), file.path(OUT,"E1_sizefactors.tsv"), sep="\t")
cat("=== per-cell-type mRNA size factors ===\n"); print(round(sort(sizef,decreasing=TRUE)))

## build subsampled SCE reference for MuSiC (cap 400/type)
cap<-400; idx<-unlist(lapply(unique(ct), function(k){w<-which(ct==k); if(length(w)>cap) sample(w,cap) else w}))
sce <- SingleCellExperiment(assays=list(counts=as.matrix(m[,idx])),
        colData=DataFrame(cell_type=ct[idx], sampleID=samp[idx]))

correct <- function(theta){ Pc<-sweep(theta,2,sizef[colnames(theta)],"/"); Pc/rowSums(Pc) }
res <- list()
for(co in c("GSE16134","GSE173078","GSE223924","GSE10334","GSE33774")){
  f<-file.path(P,"results/tables/proportions",paste0(co,"_oral_InstaPrism.tsv"))
  if(!file.exists(f)) next
  pr<-fread(f); th<-as.matrix(pr[,-1]); rownames(th)<-pr$sample
  cellf<-correct(th)
  fwrite(data.table(sample=rownames(cellf), cellf), file.path(OUT,paste0(co,"_InstaPrism_cellfrac.tsv")), sep="\t")
  res[[paste(co,"InstaPrism_mRNAfrac")]]<-data.table(cohort=co,method="InstaPrism_mRNAfrac",
     Epithelial=mean(th[,"Epithelial"]), immune=mean(rowSums(th[,intersect(imm,colnames(th)),drop=FALSE])))
  res[[paste(co,"InstaPrism_cellfrac")]]<-data.table(cohort=co,method="InstaPrism_cellfrac_E1",
     Epithelial=mean(cellf[,"Epithelial"]), immune=mean(rowSums(cellf[,intersect(imm,colnames(cellf)),drop=FALSE])))
}
## MuSiC (cell-size aware) on RNA-seq cohorts as orthogonal check
for(co in c("GSE173078","GSE223924")){
  d<-readRDS(file.path(P,"data/processed/bulk",paste0(co,"_linear.rds")))
  bulk<-as.matrix(d$expr)
  mp<-tryCatch(music_prop(bulk.mtx=bulk, sc.sce=sce, clusters="cell_type", samples="sampleID",
              select.ct=unique(ct), verbose=FALSE)$Est.prop.weighted, error=function(e){message("[MuSiC][WARN] ",co,": ",conditionMessage(e));NULL})
  if(!is.null(mp)){
    fwrite(data.table(sample=rownames(mp), mp), file.path(OUT,paste0(co,"_MuSiC.tsv")), sep="\t")
    res[[paste(co,"MuSiC")]]<-data.table(cohort=co,method="MuSiC",
       Epithelial=mean(mp[,"Epithelial"]), immune=mean(rowSums(mp[,intersect(imm,colnames(mp)),drop=FALSE])))
  }
}
tab<-rbindlist(res,fill=TRUE)
fwrite(tab, file.path(OUT,"E1_immune_total_by_method.tsv"), sep="\t")
cat("\n=== 免疫总量 by method (目标: 接近 DNAm ~0.42) ===\n"); print(tab[,.(cohort,method,Epithelial=round(Epithelial,3),immune=round(immune,3))])
message("[done] E1 mRNA-content correction")
