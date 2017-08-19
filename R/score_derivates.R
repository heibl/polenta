## This code is part of the polenta package
## © F.-S. Krah, C. Heibl 2017 (last update 2017-08-18)
#' Calculate additional scores based on residue pair score
#'
#' @param polenta_obj usually result of \code{\link{guidance}} with score method = "Rcpp"
#' @param score character, currently c("gcsc", "rprsc"); gcsc = guidance column score; rprsc = residue pair residue score.
#'
#' @details The GUIDANCE column score can be e.g. used as a sites weights in RAxML (flag -g)
#'
#' @return data.frame or list of data.frames with scores
#' @author Franz-Sebastian Krah
#' @import parallel
#' @export


daughter_scores <- function(rpsc, score = "gcsc"){

  if(!inherits(rpsc, c("polenta"))){
    stop("rpsc not if class 'polenta'")
  }
  sc <- polenta_obj@score
  base.msa <- polenta_obj@base.msa

  switch(score,
    "gcsc"={
      ## calculate GUIDANCE score
      gcsc <- colMeans(sc, na.rm = TRUE)
      gcsc <- data.frame(col = 1:length(gcsc), score = gcsc)
      if(na.rm){
        gcsc <- gcsc[!is.na(gcsc$score), ]
      }
    },

    "rprsc"={
      ## Calculate residue pair residue score
      fac <- apply(combn(nrow(base.msa), 2), 2, paste, collapse = "-")
      fac_list <- foreach(i = 1:nrow(base.msa)) %do% grep(paste0(i, "\\b"), fac)

      rprsc <- mclapply(fac_list, function(x) {
        colMeans(sc[x, ], na.rm = T)
      }, mc.cores = detectCores())
      rprsc <- do.call(rbind, rprsc)
      colnames(rprsc) <- 1:ncol(rprsc)
      rownames(rprsc) <- rownames(base.msa)
      if(na.rm){
        rprsc <- rprsc[, !apply(rprsc, 2, function(x) !any(!is.na(x)))]
      }
    })

  return(mget(score))
}
