## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-05-05)

#' @include polentaDNA-class.R
#' @importFrom methods new
#' @export

"polentaDNA" <- function(msa, scores, method){

  new("polentaDNA",
      msa = msa,
      scores = scores,
      method = method
  )
}

## setMethod: indexing, extracting scores,

setMethod("show", "polentaDNA",
          function(object){
            cat(nrow(object@msa), "DNA sequences with", object@method, "scores")
          })
