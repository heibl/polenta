#' computes residue pair score scores by comparing a reference MSA with alternative MSAs
#' @param ref reference MSA, object of class .. must be object of \code{\link{gbbin.matrix}}
#' @param alt either single or many alternative MSA(s). If single class \code{matrix}, if multiple class \code{list}. Must be object of \code{\link{gbbin.matrix}} or \code{\link{gbbin.list}}
#' @author Franz-Sebastian Krah
#'
#' @import Rcpp
#' @export

msa_RPscore <- function(ref, alt){
  ## Recode MSAs to Cmatrix
  #-------------------------
  ref <- Cmatrix(gbbin(ref))


  if(inherits(alt, "list")){
    alt <- lapply(alt, gbbin)
    alt <- lapply(alt, Cmatrix)

    res <- add_msa(ref = ref, com = alt[[1]])
    for(i in 2:length(alt)){
      res[,4] <- res[,4] + add_msa(ref = ref, com = alt[[i]])[,4]
      gc()
    }
  }else{
    res <- add_msa(ref = ref, com = alt)
  }
  if(inherits(alt, "list"))
    res[,4] <- res[,4]/length(alt)
  return(res)
}
