#' polenta: PASTA Optimal confidENce Transivitiy Alignment
#'
#' Multiple Sequence Alignment with PASTA and GUIDANCE:
#' PASTA (Mirarab, Nguyen, and Warnow 2014) is a highly accurate and scalable algorithm for multiple sequence
#' alignment. GUIDANCE  (Penn, Privman, Landan, Grauer, and Pupko 2010) is a well-performing method for the
#' assessment of alignment uncertainty. Both methods are available as stand-alone programs, but,
#' unfortunately, cannot be combined. This package provides a reimplementation of PASTA and GUIDANCE
#' to combine their power to accurately align large numbers of sequences.
#'
#' @useDynLib polenta
#' @importFrom Rcpp sourceCpp
#'
#' @docType package
#' @name polenta
NULL
