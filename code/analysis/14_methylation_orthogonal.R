#!/usr/bin/env Rscript
# Stage 7A - orthogonal validation via DNA methylation cell composition.
# GSE173081 (EPIC, gingiva) shares donors with the RNA-seq GSE173078. Estimate
# epithelial/fibroblast/immune fractions from DNAm with EpiDISH (centEpiFibIC, the
# saliva/buccal-validated oral reference; reused) and correlate cross-modally with the
# RNA-derived (oral InstaPrism) fractions on the same donors. Tests that RNA-based
# proportions reflect real tissue composition, not a reference artifact.
suppressMessages({library(data.table); library(EpiDISH)})
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
MDIR <- file.path(P,"data/interim/GSE173081_meth"); OUT <- file.path(P,"results/tables/validation"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

## ---- build M-value matrix from per-sample methylation_levels txt, convert to beta ----
files <- Sys.glob(file.path(MDIR,"*methylation_levels_*.txt.gz"))
message("[meth] ", length(files), " sample files")
mats <- lapply(files, function(f){
  d <- fread(f); setnames(d, c("probe","val")); s <- sub(".*levels_([A-Za-z0-9]+)\\.txt.*","\\1", basename(f))
  setnames(d,"val",s); d
})
M <- Reduce(function(a,b) merge(a,b,by="probe"), mats)
probes <- M$probe; Mv <- as.matrix(M[,-1]); rownames(Mv) <- probes
beta <- 2^Mv/(2^Mv+1)                       # M-value -> beta
message(sprintf("[meth] beta matrix %d probes x %d samples", nrow(beta), ncol(beta)))
saveRDS(beta, file.path(P,"data/processed/GSE173081_beta.rds"))

## ---- EpiDISH: epithelial / fibroblast / immune fractions ----
data(centEpiFibIC.m)
res <- epidish(beta.m=beta, ref.m=centEpiFibIC.m, method="RPC")
dnam_frac <- res$estF                        # samples x {Epi,Fib,IC}
fwrite(data.table(sample=rownames(dnam_frac), dnam_frac), file.path(OUT,"stage7_dnam_fractions.tsv"), sep="\t")
message("[EpiDISH] DNAm fractions: ", paste(colnames(dnam_frac),collapse=","))

## ---- RNA fractions (oral InstaPrism on GSE173078) -> collapse to Epi/Fib/IC ----
rna <- fread(file.path(P,"results/tables/proportions/GSE173078_oral_InstaPrism.tsv"))
rs <- as.matrix(rna[,-1]); rownames(rs) <- rna$sample
imm <- c("Tcell","Bcell","Plasma","Myeloid","Mast","NK","Neutrophil")
stroma <- c("Fibroblast","Endothelial","Pericyte")
rna3 <- cbind(Epi=rowSums(rs[,intersect("Epithelial",colnames(rs)),drop=FALSE]),
              Fib=rowSums(rs[,intersect(stroma,colnames(rs)),drop=FALSE]),
              IC =rowSums(rs[,intersect(imm,colnames(rs)),drop=FALSE]))
rna3 <- rna3/rowSums(rna3)

## ---- match donors (AC ids) and correlate cross-modally ----
common <- intersect(rownames(dnam_frac), rownames(rna3))
message("[match] ", length(common), " donors shared between DNAm and RNA")
if(length(common) >= 5){
  d <- dnam_frac[common,,drop=FALSE]; r <- rna3[common,,drop=FALSE]
  ccc <- function(t,e) 2*cov(t,e)/(var(t)+var(e)+(mean(t)-mean(e))^2)
  cortab <- data.frame(compartment=c("Epi","Fib","IC"),
    Pearson=sapply(c("Epi","Fib","IC"), function(k) cor(d[,k], r[,k])),
    Spearman=sapply(c("Epi","Fib","IC"), function(k) cor(d[,k], r[,k], method="spearman")),
    CCC=sapply(c("Epi","Fib","IC"), function(k) ccc(d[,k], r[,k])),
    DNAm_mean=colMeans(d), RNA_mean=colMeans(r))
  fwrite(cortab, file.path(OUT,"stage7_crossmodal_concordance.tsv"), sep="\t")
  print(cortab)
} else message("[warn] too few shared donors")
message("[done] Stage 7A methylation orthogonal validation")
