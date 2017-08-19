#' MSA coding: gap to 0 and base to 1
#' @export
#'
gbbin <- function(msa){
  msa <- (msa != "-")*1
  return(msa)
}
