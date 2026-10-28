#!/usr/bin/env Rscript
# Stage 4a - preprocess all bulk periodontitis cohorts to LINEAR-scale gene x sample
# matrices + harmonized sample metadata (group: health / periodontitis / peri-implantitis).
# Reuses verified gene-collapsed matrices from periodontitis_field_defect where available
# (GSE10334/16134 log2 RMA, GSE173078 logCPM) -> de-log; processes GSE223924 (RNA-seq counts),
# GSE23586 & GSE33774 (arrays) from raw via GEOquery + platform annotation (MaxMean collapse).
suppressMessages({library(data.table)})
P <- "/data/workplace/wyx/oral_deconv_reference_20260707"
OUT <- file.path(P,"data/processed/bulk"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
RE  <- file.path(P,"data/external/bulk_processed_reuse")
LOG <- function(...) cat(sprintf(...), "\n")

save_cohort <- function(gse, expr, meta){
  # expr: genes x samples (LINEAR); meta: data.frame(sample, group, patient?, platform)
  stopifnot(all(meta$sample %in% colnames(expr)))
  expr <- expr[, meta$sample, drop=FALSE]
  saveRDS(list(expr=expr, meta=meta), file.path(OUT, paste0(gse,"_linear.rds")))
  fwrite(data.table(sample=colnames(expr)), file.path(OUT, paste0(gse,"_samples.txt")))
  fwrite(as.data.table(meta), file.path(OUT, paste0(gse,"_meta.tsv")), sep="\t")
  LOG("[save] %s: %d genes x %d samples; groups: %s", gse, nrow(expr), ncol(expr),
      paste(names(table(meta$group)), table(meta$group), sep="=", collapse=", "))
}

## ---------- 1. reused gene matrices (de-log to linear) ----------
reuse_delog <- function(gse, matfile, metafile, is_logcpm=FALSE){
  tryCatch({
    m <- fread(file.path(RE, matfile)); g <- m[[1]]; M <- as.matrix(m[,-1]); rownames(M) <- g
    M <- 2^M                                   # de-log: log2 -> linear (RMA or logCPM)
    md <- fread(file.path(RE, metafile))
    # map site_status -> group
    grp <- ifelse(md$site_status %in% c("affected","lesion","diseased","periodontitis"), "periodontitis",
           ifelse(md$site_status %in% c("unaffected","healthy","control","non_case_site"), "health", md$site_status))
    # sample id column: prefer gsm matching colnames, else sample_id
    scol <- if("gsm" %in% names(md) && any(md$gsm %in% colnames(M))) "gsm" else
            if("sample_id" %in% names(md) && any(md$sample_id %in% colnames(M))) "sample_id" else names(md)[1]
    meta <- data.frame(sample=md[[scol]], group=grp,
                       patient=if("patient_id"%in%names(md)) md$patient_id else NA,
                       platform=if("platform"%in%names(md)) md$platform else NA, stringsAsFactors=FALSE)
    meta <- meta[meta$sample %in% colnames(M) & meta$group %in% c("health","periodontitis"),]
    save_cohort(gse, M, meta)
  }, error=function(e) LOG("[%s][FAIL] %s", gse, conditionMessage(e)))
}
reuse_delog("GSE10334","GSE10334_expr_gene.tsv","GSE10334_meta.tsv")
reuse_delog("GSE16134","GSE16134_expr_gene.tsv","GSE16134_meta.tsv")
reuse_delog("GSE173078","GSE173078_logCPM.tsv","GSE173078_meta.tsv", is_logcpm=TRUE)

## ---------- 2. GSE223924 (RNA-seq raw gene counts) ----------
tryCatch({
  cnt <- fread(file.path(P,"data/external/GSE223924/GSE223924_raw_count.csv.gz"))
  g <- cnt[[1]]; M <- as.matrix(cnt[,-1]); rownames(M) <- g
  # count colnames are H1-H10 (healthy), PT1-PT11 (periodontitis), PI1-PI11 (peri-implantitis)
  cn <- colnames(M)
  grp <- ifelse(grepl("^PI", cn), "peri_implantitis",
          ifelse(grepl("^PT", cn), "periodontitis",
          ifelse(grepl("^H", cn), "health", NA)))
  meta <- data.frame(sample=cn, group=grp, platform="GPL24676", stringsAsFactors=FALSE)
  # CPM linear
  libsize <- colSums(M); cpm <- t(t(M)/libsize)*1e6
  meta <- meta[!is.na(meta$group) & meta$group %in% c("health","periodontitis"),]
  save_cohort("GSE223924", cpm, meta)
}, error=function(e) LOG("[GSE223924][FAIL] %s", conditionMessage(e)))

## ---------- 3. GSE23586 (GPL570) & GSE33774 (GPL6244) via GEOquery ----------
suppressMessages(library(GEOquery))
collapse_maxmean <- function(expr, sym){
  keep <- !is.na(sym) & sym!="" & sym!="---"
  expr <- expr[keep,,drop=FALSE]; sym <- sym[keep]
  mm <- rowMeans(expr, na.rm=TRUE)
  ord <- order(mm, decreasing=TRUE)
  expr <- expr[ord,,drop=FALSE]; sym <- sym[ord]
  expr[!duplicated(sym),,drop=FALSE] -> ex2
  rownames(ex2) <- sym[!duplicated(sym)]; ex2
}
detect_symbol_col <- function(fd){
  cand <- c("Gene Symbol","Gene symbol","gene_symbol","Symbol","GENE_SYMBOL","gene_assignment","geneSymbol")
  hit <- intersect(cand, colnames(fd)); if(length(hit)) return(hit[1])
  gg <- grep("symbol|gene", colnames(fd), ignore.case=TRUE, value=TRUE); if(length(gg)) return(gg[1]); NA
}
parse_group <- function(pd){
  s <- tolower(apply(pd[,grep("characteristics|title|source_name", colnames(pd)),drop=FALSE],1,paste,collapse=" "))
  ifelse(grepl("peri-implant|periimplant|periimplantitis", s), "peri_implantitis",
   ifelse(grepl("periodontitis|parodontitis|diseased|chronic|aggressive", s), "periodontitis",
    ifelse(grepl("health|control|normal", s), "health", NA)))
}
process_geo_array <- function(gse_id, seriesfile, platform){
  tryCatch({
    e <- getGEO(filename=seriesfile, getGPL=TRUE)
    ex <- Biobase::exprs(e); fd <- Biobase::fData(e); pd <- Biobase::pData(e)
    sym_col <- detect_symbol_col(fd)
    if(is.na(sym_col)) stop("no gene symbol column in fData")
    sym <- fd[[sym_col]]
    if(platform=="GPL6244") sym <- sub(" //.*","", sub("^[^/]*// ","", sym))   # gene_assignment first symbol
    sym <- trimws(sapply(strsplit(as.character(sym)," /// |///|,"), `[`, 1))
    ex2 <- collapse_maxmean(ex, sym)
    if(max(ex2, na.rm=TRUE) < 100) ex2 <- 2^ex2                 # de-log if log-scale
    grp <- parse_group(pd)
    meta <- data.frame(sample=rownames(pd), group=grp, platform=platform, stringsAsFactors=FALSE)
    meta <- meta[meta$group %in% c("health","periodontitis"),]
    ex2 <- ex2[, meta$sample, drop=FALSE]
    save_cohort(gse_id, ex2, meta)
  }, error=function(e) LOG("[%s][FAIL] %s", gse_id, conditionMessage(e)))
}
process_geo_array("GSE23586",
  Sys.glob(file.path(P,"data/external/GSE23586","*series_matrix.txt.gz"))[1], "GPL570")
process_geo_array("GSE33774",
  Sys.glob(file.path(P,"data/external/GSE33774","*series_matrix.txt.gz"))[1], "GPL6244")

LOG("[done] Stage 4a bulk preprocessing")
