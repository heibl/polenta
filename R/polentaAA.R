## This code is part of the polenta package
## © Franz-S. Krah 2017 (last update 2017-08-21)

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
