#!/usr/bin/env Rscript
# E2 (keystone) - hierarchical two-stage deconvolution to fix the lymphocyte collapse (T/B->0).
# Stage 1: coarse reference (structural types + collapsed Immune) -> immune-total per sample
#          (then mRNA-content corrected).
# Stage 2: immune-ONLY sub-reference (7 immune types from GSE164241), restricted to
#          IMMUNE-MARKER genes (immune-specific, structural expression ~0) and IG-masked,
#          -> within-immune subtype fractions; rescaled to the Stage-1 immune total.
# Decision gate: does Stage 2 recover non-zero, rank-plausible T/B cells?
suppressMessages({library(Matrix); library(data.table); library(InstaPrism)})
set.seed(0)
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
OUT <- file.path(P,"results/tables/improvement"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
imm  <- c("Tcell","Bcell","Plasma","Myeloid","Mast","NK","Neutrophil")
struct <- c("Epithelial","Fibroblast","Endothelial","Pericyte","Melanocyte")
sizef <- readRDS(file.path(P,"results/models/celltype_mRNA_sizefactors.rds"))

ref <- readRDS(file.path(P,"results/models/oral_reference_counts.rds"))
m <- ref$counts; ct <- ref$cell_type; cs <- ref$cell_state

## ---- coarse reference (collapse immune -> "Immune") ----
ct_coarse <- ifelse(ct %in% imm, "Immune", ct)
refPhi_coarse <- refPrepare(sc_Expr=as.matrix(m), cell.type.labels=ct_coarse, cell.state.labels=ct)
## coarse size factors (incl "Immune" = mean mRNA of all immune cells)
sizef_coarse <- tapply(Matrix::colSums(m), ct_coarse, mean)

## ---- immune-marker genes (immune-specific; structural ~0) + IG mask ----
cpm <- t(t(m)/Matrix::colSums(m))*1e6
mu <- sapply(sort(unique(ct)), function(k) Matrix::rowMeans(cpm[,ct==k,drop=FALSE])); rownames(mu)<-rownames(m)
imm_max <- apply(mu[,imm,drop=FALSE],1,max); struct_max <- apply(mu[,struct,drop=FALSE],1,max)
ig <- grep("^IGH|^IGK|^IGL|^JCHAIN$|^MZB1$", rownames(mu), value=TRUE)
imm_markers <- setdiff(rownames(mu)[imm_max > 3*pmax(struct_max,1) & imm_max > 5], ig)
message(sprintf("[E2] %d immune-marker genes (immune-specific, IG-masked)", length(imm_markers)))

## ---- immune-only sub-reference on immune-marker genes ----
imm_cells <- ct %in% imm
refPhi_imm <- refPrepare(sc_Expr=as.matrix(m[imm_markers, imm_cells]),
                         cell.type.labels=ct[imm_cells], cell.state.labels=cs[imm_cells])

correct_total <- function(theta, sf){ Pc<-sweep(theta,2,sf[colnames(theta)],"/"); Pc/rowSums(Pc) }

res <- list(); cmp <- list()
for(co in c("GSE16134","GSE173078","GSE223924","GSE10334","GSE33774")){
  d<-readRDS(file.path(P,"data/processed/bulk",paste0(co,"_linear.rds"))); bulk<-as.matrix(d$expr)
  ## Stage 1 -> immune total (mRNA-corrected)
  ip1<-InstaPrism(bulk_Expr=bulk, refPhi_cs=refPhi_coarse, n.core=8, verbose=FALSE)
  th1<-t(ip1@Post.ini.ct@theta); rownames(th1)<-colnames(bulk)
  th1c<-correct_total(th1, sizef_coarse); immune_total<-th1c[,"Immune"]
  ## Stage 2 -> within-immune subtypes on immune-marker genes
  bimm<-bulk[intersect(rownames(bulk), imm_markers),,drop=FALSE]
  ip2<-InstaPrism(bulk_Expr=bimm, refPhi_cs=refPhi_imm, n.core=8, verbose=FALSE)
  th2<-t(ip2@Post.ini.ct@theta); rownames(th2)<-colnames(bulk)         # within-immune mRNA-frac, sums to 1
  ## two variants: (A) NO within-immune mRNA correction, (B) WITH correction
  th2A<-th2/rowSums(th2)                                                # A: mRNA-fraction as-is
  th2B<-sweep(th2,2,sizef[colnames(th2)],"/"); th2B<-th2B/rowSums(th2B) # B: mRNA-content corrected
  th2c<-th2A                                                            # PRIMARY = A (no double-correction; histology-anchored)
  final_imm<-sweep(th2c, 1, immune_total, "*")                         # rescale to immune total
  # record both within-immune variants for the anchoring sensitivity
  fwrite(data.table(sample=rownames(th2A), th2A), file.path(OUT,paste0(co,"_withinImmune_noMRNAcorr.tsv")), sep="\t")
  fwrite(data.table(sample=rownames(th2B), th2B), file.path(OUT,paste0(co,"_withinImmune_MRNAcorr.tsv")), sep="\t")
  ## assemble full proportion: structural (from Stage1 corrected) + immune subtypes
  full<-cbind(th1c[,setdiff(colnames(th1c),"Immune"),drop=FALSE], final_imm)
  full<-full/rowSums(full)
  fwrite(data.table(sample=rownames(full), full), file.path(OUT,paste0(co,"_hierarchical.tsv")), sep="\t")
  res[[co]]<-data.table(cohort=co, immune_total=mean(immune_total),
      Tcell=mean(final_imm[,"Tcell"]), Bcell=mean(final_imm[,"Bcell"]), Plasma=mean(final_imm[,"Plasma"]),
      Myeloid=mean(final_imm[,"Myeloid"]), Mast=mean(final_imm[,"Mast"]), NK=mean(final_imm[,"NK"]),
      Neutrophil=mean(final_imm[,"Neutrophil"]))
  message("[",co,"] hierarchical done; immune_total=",round(mean(immune_total),3))
}
tab<-rbindlist(res,fill=TRUE)
fwrite(tab, file.path(OUT,"E2_hierarchical_immune.tsv"), sep="\t")
cat("\n=== E2 层级去坍缩: 免疫亚型平均比例 (关键: T/B 是否非零) ===\n")
print(tab[,.(cohort,immune_total=round(immune_total,3),Tcell=round(Tcell,4),Bcell=round(Bcell,4),
             Plasma=round(Plasma,3),Myeloid=round(Myeloid,4),Neutrophil=round(Neutrophil,4))])
message("[done] E2 hierarchical de-collapse")
