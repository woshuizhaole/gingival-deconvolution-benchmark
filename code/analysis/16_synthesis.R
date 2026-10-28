#!/usr/bin/env Rscript
# Stage 8 - synthesis: (1) LM22-vs-oral flip catalog (composition directions across cohorts,
# cross-referenced to the published A1-A9 claims); (2) per-(gene,celltype) composition-vs-
# intrinsic 2x2 classification for the primary cohort; (3) bulk-DEG attribution to channels;
# (4) cross-cohort/platform reproducibility of the oral composition directions.
suppressMessages({library(data.table); library(limma)})
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
TAB <- file.path(P,"results/tables"); OUT <- file.path(TAB,"synthesis"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

## ---- (1) FLIP CATALOG: oral vs LM22 composition directions ----
comp <- fread(file.path(TAB,"stage4_composition_tests.tsv"))
comp <- comp[test=="clr_limma"]
# per cohort x deconv x cell_type direction
dir_tab <- comp[, .(cohort, deconv, cell_type, effect, FDR, dir=sig)]
# focus: oral_InstaPrism vs LM22_DWLS on shared immune types
imm <- c("Tcell","Bcell","Plasma","Myeloid","Mast","NK","Neutrophil")
flip <- dcast(dir_tab[deconv %in% c("oral_InstaPrism","LM22_DWLS") & cell_type %in% imm],
              cohort+cell_type~deconv, value.var="dir")
setnames(flip, c("oral_InstaPrism","LM22_DWLS"), c("oral","LM22"), skip_absent=TRUE)
flip[, flipped := (oral=="UP" & LM22=="DOWN") | (oral=="DOWN" & LM22=="UP")]
fwrite(flip, file.path(OUT,"flip_catalog.tsv"), sep="\t")
cat("=== FLIP CATALOG (oral InstaPrism vs LM22 DWLS, CLR) ===\n"); print(flip)
cat("\nFlips per cell type (across cohorts):\n"); print(flip[flipped==TRUE, .N, by=cell_type])

## ---- (2) cross-cohort/platform reproducibility of ORAL composition directions ----
oral_dir <- dcast(dir_tab[deconv=="oral_InstaPrism"], cell_type~cohort, value.var="dir")
fwrite(oral_dir, file.path(OUT,"oral_composition_directions_by_cohort.tsv"), sep="\t")
cat("\n=== Oral composition direction by cohort (reproducibility) ===\n"); print(oral_dir)

## ---- (3) per-(gene,celltype) composition-vs-intrinsic 2x2 (primary cohort GSE16134) ----
CO <- "GSE16134"
comp_co <- comp[cohort==CO & deconv=="oral_InstaPrism", .(cell_type, comp_FDR=FDR, comp_sig=(FDR<0.05))]
intr <- as.data.table(readRDS(file.path(TAB,"intrinsic",paste0(CO,"_intrinsic_TOAST.rds"))))
intr[, intr_sig := FDR<0.05]
intr <- merge(intr, comp_co, by="cell_type", all.x=TRUE)
intr[, class := fifelse(is.na(comp_sig), "no_comp_test",
              fifelse(comp_sig & intr_sig, "both",
              fifelse(!comp_sig & intr_sig, "intrinsic_only",
              fifelse(comp_sig & !intr_sig, "composition_only", "neither"))))]
tab2x2 <- intr[, .N, by=.(cell_type, class)]
tab2x2w <- dcast(tab2x2, cell_type~class, value.var="N", fill=0)
fwrite(tab2x2w, file.path(OUT,"composition_intrinsic_2x2_GSE16134.tsv"), sep="\t")
cat("\n=== per-(gene,celltype) 2x2 classification (GSE16134) ===\n"); print(tab2x2w)

## ---- (4) bulk-DEG attribution: total bulk DE, fraction explained by intrinsic channel ----
d <- readRDS(file.path(P,"data/processed/bulk",paste0(CO,"_linear.rds")))
Y <- log2(as.matrix(d$expr)+1); meta <- d$meta
grp <- factor(meta$group[match(colnames(Y), meta$sample)], levels=c("health","periodontitis"))
fit <- eBayes(lmFit(Y, model.matrix(~grp))); bulkDE <- topTable(fit, coef=2, number=Inf, sort.by="none")
bulkDE <- data.table(gene=rownames(bulkDE), bulk_logFC=bulkDE$logFC, bulk_FDR=bulkDE$adj.P.Val)
bulk_degs <- bulkDE[bulk_FDR<0.05]
# a bulk DEG has an intrinsic explanation if any cell type shows intrinsic FDR<0.05
intr_genes <- unique(intr[intr_sig==TRUE, .(gene, cell_type)])
deg_intr <- intr_genes[gene %in% bulk_degs$gene]
attribution <- data.table(
  n_bulk_DEG = nrow(bulk_degs),
  n_DEG_with_intrinsic = uniqueN(deg_intr$gene),
  frac_DEG_with_intrinsic = round(uniqueN(deg_intr$gene)/nrow(bulk_degs),3),
  n_DEG_intrinsic_epithelial = uniqueN(deg_intr[cell_type=="Epithelial"]$gene))
fwrite(attribution, file.path(OUT,"bulkDEG_attribution_GSE16134.tsv"), sep="\t")
cat("\n=== bulk-DEG attribution (GSE16134) ===\n"); print(attribution)
cat("\n[done] Stage 8 synthesis\n")
