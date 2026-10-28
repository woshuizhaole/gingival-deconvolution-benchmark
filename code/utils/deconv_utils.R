# deconv_utils.R - signature-based deconvolution for the FAIR same-algorithm comparison
# of the oral signature vs LM22. REUSES published implementations (DWLS solvers, nnls pkg)
# rather than hand-rolling. Both signatures go through the identical solver so any
# difference is attributable to the signature, not the method.
suppressMessages({library(DWLS); library(nnls); library(parallel)})

# ---- REUSE DWLS published solvers (Tsoucas et al. 2019) ----
# solveDampenedWLS / solveOLS / solveSVR + trimData are DWLS's own functions.
deconv_dwls <- function(bulk, sig, method=c("dampened","ols","svr"), cores=8){
  method <- match.arg(method); sig <- as.matrix(sig)
  cn <- colnames(bulk)
  frl <- mclapply(cn, function(s){
    b <- bulk[,s]; names(b) <- rownames(bulk)
    f <- tryCatch({
      tr <- DWLS::trimData(sig, b)
      switch(method,
             dampened = DWLS::solveDampenedWLS(tr$sig, tr$bulk),
             ols      = DWLS::solveOLS(tr$sig, tr$bulk),
             svr      = DWLS::solveSVR(tr$sig, tr$bulk))
    }, error=function(e){ setNames(rep(NA_real_, ncol(sig)), colnames(sig)) })
    f <- f[colnames(sig)]; f[is.na(f)] <- 0; f[f<0] <- 0
    if(sum(f)>0) f/sum(f) else f
  }, mc.cores=cores)
  fr <- do.call(rbind, frl); rownames(fr) <- cn; colnames(fr) <- colnames(sig); fr
}

# ---- nnls (reuses the nnls package solver) ----
deconv_nnls <- function(bulk, sig){
  g <- intersect(rownames(bulk), rownames(sig)); S <- as.matrix(sig[g,,drop=FALSE]); B <- as.matrix(bulk[g,,drop=FALSE])
  fr <- apply(B, 2, function(y){ f <- nnls::nnls(S, y)$x; f[f<0]<-0; if(sum(f)>0) f/sum(f) else f })
  fr <- t(fr); colnames(fr) <- colnames(S); rownames(fr) <- colnames(B); fr
}

# ---- bMIND: REUSE author's implementation by sourcing randel/MIND source ----
load_bmind <- function(d="/data/workplace/wyx/oral_deconv_reference_20260707/code/utils/MIND_src/R"){
  for(f in c("deconv_em_func.r","bmind_func.r","bmind_de.r")){
    fp <- file.path(d,f); if(file.exists(fp)) sys.source(fp, envir=globalenv())
  }
  invisible(c(bMIND=exists("bMIND"), bmind_de=exists("bmind_de")))
}
