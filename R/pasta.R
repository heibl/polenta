## This code is part of the polenta package
## © C. Heibl 2016 (last update 2017-10-13)

#' @title Ultra-Large Multiple Sequence Alignment with PASTA
#' @description Provides a complete reimplementation of the PASTA algorithm
#'   (Mirarab, Nguyen, and Warnow 2014) in R.
#' @param seqs An object of class \code{"\link{DNAbin}"} or \code{"\link{AAbin}"}
#'   containing unaligned sequences of DNA or amino acids.
#' @param gt \emph{Currently unused.}
#' @param k An integer giving the size of cluster in which the dataset is split.
#' @param msa.program A character string giving the alignment program to use;
#'   currently only \code{"mafft"} is possible.
#' @param method A character string choosing a method of the alignment program;
#'   default is \code{"localpair"}, see \code{\link{mafft}} for possible
#'   options.
#' @param exec A character string giving the path to the alignment program
#'   executable.
#' @param parallel Logical, indicating if the function should be run in parallel
#'   or serial mode. \emph{Currently unused!}
#' @param ncore An integer giving the number of cores to use in parallel mode.
#'   \emph{Currently unused!}
#' @return An object of class \code{"\link[=polentaDNA-class]{polentaDNA}"}.
#' @seealso \code{\link{extractMSA}} for extractiong the multiple sequence
#'   alignment of a \code{"\link[=polentaDNA-class]{polentaDNA}"} object.
#' @importFrom ape del.gaps dist.dna dist.aa
#' @import foreach
#' @importFrom igraph as_edgelist
#' @importFrom ips mafft mafft.merge
#' @importFrom utils globalVariables
#' @export

pasta <- function(seqs, gt, k = 200,
                  msa.program = "mafft", method = "localpair", exec,
                  parallel = FALSE, ncore){
  
  ## declare i to be a global variable; this is necessary because
  ## foreach uses non-standard evaluation that codetools is not aware of
  ## [http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html]
  ## Does not work [CH 2018-01-23]
  #globalVariables("i")

  ## remove gaps from aligned sequences
  ## ----------------------------------
  if (is.matrix(seqs)) {
    seqs <- del.gaps(seqs)
  }

  ## less than k species will be aligned with MAFFT-LINSI
  ## ----------------------------------------------------
  if (length(seqs) <= k){
    cat(length(seqs), "species will be aligned with MAFFT L-INS-i\n")

    seqs <- mafft(seqs, method = "localpair", gt = gt, exec = exec)

    ## more than k species will be aligned with PASTA
    ## ----------------------------------------------
  } else {

    ## This is a quick hack to get an inital guide tree
    ## Should be replaced by the method used by Mirarab and Warnow
    ## or perhaps a hybrid with taxonomy.
    if (missing(gt)){
      gt <- mafft(seqs, method = "auto")
      if (inherits(gt, "DNAbin")){
        gt <- dist.dna(gt, model = "F81")
      } else {
        gt <- dist.aa(gt)
      }
      gt <- nj(gt)
    }

    ## split dataset in subsets of size <= k
    ## -------------------------------------
    cat("Split dataset by centroid decomposition\n")
    subtrees <- centroidDecomposition(gt, k = k)
    subtrees <- lapply(subtrees, function(z) z$tip.label)
    names(subtrees) <- paste0("S", seq_along(subtrees))
    seqs <- lapply(subtrees, function(seqs, subtrees) seqs[subtrees], seqs = seqs)
    subtree_names <- names(seqs)

    ## alignment of subtrees
    ## ---------------------
    cat("Alignment of", length(seqs), "subtrees\n")
    s <- lapply(seqs[2], mafft, exec = exec)
    pb <- txtProgressBar(max = length(seqs), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    cl <- makeCluster(ncore)
    registerDoSNOW(cl)
    seqs <- foreach(i = 1:length(seqs),
                    .packages = c('ips', 'ape'),
                    .options.snow = opts)  %dopar% {
                      mafft(x = seqs[[i]], method = method, exec = exec)
                    }
    stopCluster(cl)
    names(seqs) <- subtree_names

    ## compute spanning tree of subsets
    ## --------------------------------
    cat("\nCompute spanning tree connecting subtrees\n")
    st <- spanningTree(gt, subtrees)
    # save.image("devworkspace.rda")

    ## do profile-alignment
    ## --------------------
    cat("Profile alignment of pairs of subtrees\n")
    e <- igraph::as_edgelist(st)
    merger <- function(seqlist, index, exec){
      mafft.merge(seqlist[index], exec = exec)
    }
    seqs <- apply(e, 1, merger, seqlist = seqs, exec = exec)
    names(seqs) <- paste(e[, 1], e[, 2], sep = "-")

    ## do transitivity merging
    ## -----------------------
    cat("Transitivity merging\n")
    # load("devworkspace.rda")
    # seqs <- lapply(seqs, extractMSA)

    ## calculate pairings for transitivity merging
    ## this is probably very inefficient

    vertex.set <- strsplit(names(seqs), "-")
    pairings <- function(z){
      obj <- list(); meta <- list()
      for (i in 1:(length(z) - 1)){
        zz <- sapply(z, intersect, x = z[[i]])
        zz <- sapply(zz, length)
        p <- which(zz > 0)
        p <- p[p > i][1]
        if (is.na(p)) next
        p <- c(i, p)
        meta <- c(meta, list(sort(unique(unlist(z[p])))))
        obj <- c(obj, list(p))
      }
      attr(obj, "vertices") <- meta
      obj
    }

    cl <- makeCluster(ncore)
    registerDoSNOW(cl)
    while (length(seqs) > 1){
      p <- pairings(vertex.set)
      # seqs <- lapply(p, transitivityMerge, x = seqs)
      seqs <- foreach(i = 1:length(p),
                      .packages = c('ips', 'ape'),
                      .options.snow = opts)  %dopar% {
                        transitivityMerge(x = seqs, id = p[[i]], exec = exec)
                      }
      vertex.set <- attr(p, "vertices")
    }
    stopCluster(cl)
    seqs <- seqs[[1]]





    ## next steps
    ## - put pairing() in its on file
    ## - testing
    ## - parallelisation
    ## - make pairing more efficient

    # save(seqs, file = "devworkspace.rda")

  }
  seqs
}
