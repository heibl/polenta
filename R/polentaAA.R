## This code is part of the polenta package
## © Franz-S. Krah 2017 (last update 2017-11-08)

#' @title Create Objects of Class "polentaAA"
#' @description Create objects of class
#'   \code{"\link[=polentaAA-class]{polentaAA}"} from objects of class
#'   \code{"\link{AAbin}"} and a matrix derived from \code{\link{guidance}} or
#'   \code{\link{HoT}}.
#' @param msa An object of class \code{\link{AAbin}}.
#' @param scores A matrix of quality scores.
#' @param method A characters string giving the method to derive the quality
#'   scores.
#' @include polentaAA-class.R
#' @seealso \code{\link{extractMSA}}
#' @importFrom methods new
#' @export

"polentaAA" <- function(msa, scores, method){

  new("polentaAA",
    msa = msa,
    scores = scores,
    method = method
  )
}

## setMethod: indexing, extracting scores,

setMethod("show", signature = "polentaAA",
  function(object){
    cat(nrow(object@msa), "AA sequences with", object@method, "scores")
  })
