## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-11-08)

#' @title Create Objects of Class "polentaDNA"
#' @description Create objects of class
#'   \code{"\link[=polentaDNA-class]{polentaDNA}"} from objects of class
#'   \code{"\link{DNAbin}"} and a matrix derived from \code{\link{guidance}} or
#'   \code{\link{HoT}}.
#' @param msa An object of class \code{\link{DNAbin}}.
#' @param scores A matrix of quality scores.
#' @param method A characters string giving the method to derive the quality
#'   scores.
#' @include polentaDNA-class.R
#' @seealso \code{\link{extractMSA}}
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

setMethod("show", signature = "polentaDNA",
          function(object){
            cat(nrow(object@msa), "DNA sequences with", object@method, "scores")
          })
