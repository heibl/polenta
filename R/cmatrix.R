#' Cmatrix recoding
#'
#' @param msa matrix of class \code{DNAbin} or \code{AAbin}
#'
#' @return recoded Cmatrix
#' @references Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot: Bayesian alignment and phylogenetic footprinting with MCMC, BMC Evolutionary Biology, 9, 217. Supplement page 3.
#' @references Penn et al. (2010). An alignment confidence score capturing robustness to guide tree uncertainty. Molecular Biology and Evolution. 27:1759--1767
#' @details recodes gaps as increasing even and bases/acids as increasing odd integers
#'
#' @export

cmatrix <- function(msa){

  if(!inherits(msa, c("DNAbin", "AAbin")))
  stop("msa is not of class 'DNAbin' or 'AAbin'")

  msa <- Cmatrix(gbbin(as.character(msa)))
  return(msa)
}
