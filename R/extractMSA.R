## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-10-13)

#' @title Extract Sequences from POLENTA Output
#' @description Extract DNA or amino acid sequences from \code{polenta*}
#'   objects.
#' @param x An object of class \code{"\link[=polentaDNA-class]{polentaDNA}"} or
#'   \code{"\link[=polentaAA-class]{polentaAA}"}.
#' @author Christoph Heibl
#' @importFrom igraph as_edgelist
#' @importFrom ips mafft mafft.merge
#' @importFrom methods slot
#' @export

extractMSA <- function(x){

  ## checks
  if (is.list(x)){
    if (!all(sapply(x, inherits, what = c("polentaDNA", "polentaAA")))){
      stop("'x' or the elements of 'x' must be of class 'polenta*'")
    }
    lapply(x, slot, name = "msa")

  } else {
    if (!inherits(x, what = c("polentaDNA", "polentaAA"))){
      stop("'x' or the elements of 'x' must be of class 'polenta*'")
    }
    slot(x, name = "msa")
  }
}
