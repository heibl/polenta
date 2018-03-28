## This code is part of the polenta package
## © Franz-S. Krah 2017 (last update 2017-11-08)

setOldClass("AAbin")

#' @title An S4 Class to represent Amino Acid Alignments with Quality Scores
#' @description \code{"polentaAA"} holds a multiple sequence alignment of alino
#'   acids together with quality scores for each cell in the alignment.
#' @slot msa An object of class \code{\link{AAbin}}.
#' @slot scores A matrix of quality scores.
#' @slot method A characters string giving the method to derive the quality
#'   scores.
#' @seealso \code{"\link[=polentaDNA-class]{polentaDNA}"}

setClass("polentaAA",
  representation = list(
    msa = "AAbin",
    scores = "matrix",
    method = "character")
)
