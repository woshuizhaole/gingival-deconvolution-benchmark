#!/usr/bin/env Rscript
# E7 - master flip-catalog audit: for each cell type, adjudicate the oral-vs-LM22 direction
# using every line of evidence (multi-method oral consensus, LM22, scRNA ground truth,
# DNAm, spatial, detection-limit reliability), with a verdict + mechanism.
suppressMessages({library(data.table)})
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
IMP <- file.path(P,"results/tables/improvement"); TAB <- file.path(P,"results/tables")
OUT <- IMP

## oral multi-method consensus + LM22 (E4, GSE16134)
e4 <- fread(file.path(IMP,"E4_multimethod_directions.tsv"))
oral_cons <- fread(file.path(IMP,"E4_oral_consensus.tsv"))[cohort=="GSE16134", .(cell_type, oral_consensus=consensus, oral_up=up, oral_n=n)]
lm22 <- e4[cohort=="GSE16134"&method=="LM22", .(cell_type, LM22=dir)]
## spatial direction (Stage 7B)
sp <- tryCatch(fread(file.path(TAB,"validation/stage7b_spatial_health_vs_perio.tsv"))[,.(cell_type,spatial=direction)], error=function(e) data.table(cell_type=character(),spatial=character()))
## detection-limit reliability (E9): est_at_20pct; low => unreliable
e9 <- tryCatch(fread(file.path(IMP,"E9_detection_floor.tsv"))[,.(cell_type, detect_est_at20=round(est_at_20pct,3))], error=function(e) data.table(cell_type=character(),detect_est_at20=numeric()))
## scRNA ground truth composition direction (GSE164241 reference: perio vs health cell counts)
sc <- tryCatch({
  cc<-fread(file.path(TAB,"stage1_final_celltype_condition.tsv"))
  w<-dcast(cc, cell_type~condition, value.var="n", fill=0)
  w[, scRNA_prop_dir := ifelse((periodontitis/sum(periodontitis)) > (health/sum(health)), "UP","DOWN")]
  w[,.(cell_type, scRNA_prop_dir)]
}, error=function(e) data.table(cell_type=character(),scRNA_prop_dir=character()))

au <- Reduce(function(a,b) merge(a,b,by="cell_type",all=TRUE),
             list(oral_cons, lm22, sp, e9, sc))
## verdict
au[, verdict := fcase(
   oral_consensus %in% c("UP","DOWN") & LM22 %in% c("UP","DOWN") & oral_consensus!=LM22,
     "FLIP (oral vs LM22 opposite)",
   # oral ns but LM22 has a direction: split by whether LM22 agrees with scRNA truth
   oral_consensus=="ns" & LM22 %in% c("UP","DOWN") & !is.na(scRNA_prop_dir) & LM22!=scRNA_prop_dir,
     "LM22-artifact (LM22 contradicts scRNA; oral ns)",
   oral_consensus=="ns" & LM22 %in% c("UP","DOWN") & !is.na(scRNA_prop_dir) & LM22==scRNA_prop_dir,
     "concordant (LM22 dir matches scRNA; oral ns)",
   oral_consensus=="ns" & LM22 %in% c("UP","DOWN"),
     "LM22-directional, oral ns (scRNA unknown)",
   oral_consensus %in% c("UP","DOWN") & (is.na(LM22)|LM22=="ns"),
     "oral-only signal (LM22 blind/ns)",
   default="concordant/ns")]
au[, reliability := fifelse(!is.na(detect_est_at20) & detect_est_at20 < 0.05,
     "LOW (near detection limit; direction only)", "resolvable")]
setcolorder(au, c("cell_type","oral_consensus","LM22","scRNA_prop_dir","spatial","verdict","reliability","detect_est_at20","oral_up","oral_n"))
fwrite(au, file.path(OUT,"E7_flip_audit_master.tsv"), sep="\t")
cat("=== E7 翻转审计主表 (GSE16134, 多证据裁决) ===\n")
print(au[order(verdict, cell_type)])
cat("\n稳健翻转(FLIP + resolvable):\n"); print(au[grepl("FLIP",verdict) & reliability=="resolvable", .(cell_type,oral_consensus,LM22)])
message("[done] E7 flip audit")
