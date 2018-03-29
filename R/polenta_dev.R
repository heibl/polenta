## This code is part of the polenta package
## © C. Heibl 2017, F.-S. Krah (last update 2017-11-08)

#' @title Ultra-Large Multiple Sequence Alignment with PASTA
#' @description Provides a complete reimplementation of the PASTA algorithm
#'   (Mirarab, Nguyen, and Warnow 2014) in R.
#' @param seqs An object of class \code{\link{DNAbin}} or \code{\link{AAbin}}
#'   containing unaligned sequences of DNA or amino acids.
#' @param gt \emph{Currently unused.}
#' @param k An integer giving the size of cluster in which the dataset is split.
#' @param bootstrap An integer giving the number of bootstrap replicates.
#' @param method A character string choosing a method of the alignment program;
#'   see \code{\link{mafft}} for possible options.
#' @param exec A character string giving the path to the alignment program
#'   executable.
#' @param ncore An integer giving the number of cores to use in parallel mode.
#' @return An object of class \code{"\link[=polentaDNA-class]{polentaDNA}"}.
#' @seealso \code{\link{extractMSA}} for extractiong the multiple sequence
#'   alignment of a \code{"\link[=polentaDNA-class]{polentaDNA}"} object.
#' @importFrom ape del.gaps dist.dna dist.aa
#' @importFrom igraph as_edgelist
#' @importFrom ips mafft mafft.merge
#' @export

polenta_dev <- function(seqs, gt, k = 200, bootstrap = 100,
  method = "auto", exec, ncore){
  
  ## remove gaps from aligned sequences
  ## ----------------------------------
  if (is.matrix(seqs)) {
    seqs <- del.gaps(seqs)
  }
  
  ## less than k species will be aligned with MAFFT-LINSI
  ## ----------------------------------------------------
  if (length(seqs) <= k){
    cat(length(seqs), "species will be aligned with MAFFT L-INS-i\n")
    
    seqs <- guidance(seqs, ncore = ncore,
      bootstrap = bootstrap,
      method = method, msa.exec = exec)
    
    ## more than k species will be aligned with PASTA
    ## ----------------------------------------------
  } else {
    
    ## This is a quick hack to get an inital guide tree
    ## Should be replaced by the method used by Mirarab and Warnow
    ## or perhaps a hybrid with taxonomy.
    if (missing(gt)){
      gt <- mafft(seqs, method = method)
      # if (inherits(gt, "DNAbin")){
      #   gt <- dist.dna(gt, model = "F81")
      # } else {
      #   gt <- dist.aa(gt)
      # }
      # gt <- nj(gt)
      
      gt <- phangorn::dist.ml(phangorn::as.phyDat(gt), model = "F81")
      gt  <- phangorn::NJ(gt)
    }
    

    
    ## split dataset in subsets of size <= k
    ## -------------------------------------
    subtrees <- centroidDecomposition(gt, k = k)
    subtrees <- lapply(subtrees, function(z) z$tip.label)
    names(subtrees) <- paste0("S", seq_along(subtrees))
    
    ## alignment of subtrees
    ## ---------------------
    foo <- function(seqs, taxa){
      guidance(seqs[taxa], ncore = ncore,
        bootstrap = bootstrap,
        method = method, msa.exec = exec)
    }
    seqs <- lapply(subtrees, foo, seqs = seqs)
    names(seqs) <- names(subtrees)
    
    # seqs_test <- seqs
    # seqs <- seqs_test
    
    ## special case: there are only two subMSAa and they will be simply merged
    ## with out transitivity merging
    if (length(seqs) == 2){
      
      seqlist <- extractMSA(seqs)
      seqlist <- list(ips::mafft.merge(seqlist, method = method, exec = exec))
      names(seqlist) <- paste(names(seqs), collapse = "-")
      seqlist <- seqlist[[1]]
      # seqs <- reappendScores(names(seqs), merged = seqlist, scored = seqs)
      seqs <- reappendScores_gen_dev(v = seqs, s = seqlist)
      
    } else {
      ## compute spanning tree of subsets
      ## --------------------------------
      st <- spanningTree(gt, subtrees)
      # save.image("devworkspace.rda")
      
      ## do profile-alignment
      ## --------------------
      e <- igraph::as_edgelist(st)
      merger <- function(seqlist, index, exec){
        mafft.merge(seqlist[index], exec = exec)
      }
      seqlist <- extractMSA(seqs)
      seqlist <- apply(e, 1, merger, seqlist = seqlist, exec = exec)
      names(seqlist) <- paste(e[, 1], e[, 2], sep = "-")
      
      
      ## reappend scores to merged alignments
      ## ------------------------------------

      # 
      ## do transitivity merging
      ## -----------------------
      # load("devworkspace.rda")
      
      ## calculate pairings for transitivity merging
      ## this is probably very inefficient
      
      vertex.set <- strsplit(names(seqlist), "-")
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
      # seqs2 <- seqlist
      while (length(seqlist) > 1){
        p <- pairings(vertex.set)
        seqlist <- lapply(p, transitivityMerge, x = seqlist, exec = exec)
        vertex.set <- attr(p, "vertices")
      }
      seqlist <- seqlist[[1]]
      
      # seqs : MSAs with scores
      # seqs2 : merged MSA
      
      ## reappend scores to merged alignments
      ## ------------------------------------
    
      seqs <- reappendScores_gen_dev(v = seqs, s = seqlist)
      
      
      ## next steps
      ## - put pairing() in its on file
      ## - testing
      ## - parallelisation
      ## - make pairing more efficient
      
      # save.image("devworkspace.rda")
      
    }
  }
  seqs
}
