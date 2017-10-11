## This code is part of the polenta package
## © F.-S. Krah, C. Heibl 2017 (last update 2017-08-17)

#' Compare reference MSAs with alternative MSAs
#'
#' @description MSA reliability scores (Penn et al. 2010)
#'
#' @param ref of class data.frame, is the reference MSA ('BASE MSA') with
#'   sequences as columns
#' @param alt path to alternative files
#' @param bootstrap number of alt MSAs
#' @param na.rm logical if gab comparisons should be deleted
#'
#' @return list containing following score:
#' @return residue_pair_score: if one alternative MSA is supplied then this is 1
#'   if a residue pair was identically aligned as in the reference MSA and 0
#'   otherwise. If more than one, then the average of the residue pair scores
#'   that result from each comparison with the reference MSA.
#' @description Rcpp code computing basic MSA comparision. The most basic is the residue pairs residue score, which checks if residue pairs combinations are correctly aligned in both MSAs.
#' @references Penn et al. (2010). An alignment confidence score capturing
#'   robustness to guide tree uncertainty. Molecular Biology and Evolution
#'   27:1759--1767
#'
#' @author Franz-Sebastian Krah
#'
#' @seealso \code{\link{msa_set_scoreSA}}
#' @export


msa_set_scoreR_par <- function(ref, alt, bootstrap, na.rm){

  gbbin <- function(msa){
    msa <- (msa != "-")*1
    return(msa)
  }

  ## Cmatrix REF
  ref <- gbbin(as.character(ref))
  cmat_msa <- msa_recode(ref)
  ref_col2res <- cmat_msa$col2res
  ref_col2res_odd <- ref_col2res %% 2

  ## Hit matrix
  scores <- res_pair_hit_for_par(cmat_msa$col2res)


  if(is.character(alt)){
    alt <- list.files(alt, full.names = TRUE)
    alt <- lapply(alt, read.fas)
  }else{
    if(!inherits(alt, "list"))stop("'alt' must be either list or a path")
  }

  ## Compare MSAs
  for(i in 1:bootstrap){

    ## Cmatrix MSA
    com <- gbbin(as.character(alt[[i]]))
    cmat_alt <- msa_recode(com)
    alt_col2res <- cmat_alt$col2res
    alt_res2col <- cmat_alt$res2col

    ## Compare
    ## function overrides input *scores* internally
    scores <- add_msa_parallel(ref_col2res, alt_col2res,
      alt_res2col, ref_col2res_odd, hit_mat = scores)

  }
  ##
  # if(na.rm){
  #   scores <- scores[!scores[,4] == -1, ]
  # }
  scores[scores == -1] <- NA
  gc <- colMeans(scores, na.rm = TRUE)
  gc[!is.na(gc)]
  # scores[,4] <- scores[,4]/bootstrap

  return(scores)
}
