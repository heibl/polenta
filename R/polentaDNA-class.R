## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-11-08)

setOldClass("DNAbin")

#' @title An S4 Class to represent DNA Alignments with Quality Scores
#' @description \code{polentaDNA} holds a multiple sequence alignment of DNA
#'   together with quality scores for each cell in the alignment.
#' @slot msa An object of class \code{\link{DNAbin}}.
#' @slot scores A matrix of quality scores.
#' @slot method A characters string giving the method to derive the quality
#'   scores.
#' @seealso \code{"\link[=polentaAA-class]{polentaAA}"}

setClass("polentaDNA",
         representation = list(
           msa = "DNAbin",
           scores = "matrix",
           method = "character")
)






