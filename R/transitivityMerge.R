## This code is part of the polenta package
## © C. Heibl 2016 (last update 2017-06-07)


#' @title Transitivity Merge
#' @description Merges two multiple sequence alignments and their reliability scores using the transitivity
#'   criterion.
#' @param x A list containing objects of class \code{\link{polentaDNA}}.
#' @return An object of class \code{\link{polentaDNA}}.
#' @importFrom ape as.DNAbin cbind.DNAbin rbind.DNAbin
#' @export

transitivityMerge <- function(x, id = c(1, 2)){

  include_scores <- inherits(x[[1]], "polentaDNA")
  
  if (include_scores){
    msa1 <- x[[id[1]]]@msa
    msa2 <- x[[id[2]]]@msa
  } else {
    msa1 <- x[[id[1]]]
    msa2 <- x[[id[2]]]
  }

  ## shared species
  ## --------------
  shared <- intersect(rownames(msa1), rownames(msa2))
  if (!length(shared)) stop("MSAs do not overlap")
  
  ## If non-shared sequences in both MSAs extend the
  ## 3'- or 5'-ends of shared sequences, they cannot
  ## be merged properly. Hence we have to identify 
  ## and delete those positions
  ## --------------------------
  firstInformativePosition <- function(z){
    z <-  which(z %in% as.raw(4))
    z <- z == 1:length(z)
    if (any(z)){
      max(which(z)) + 1
    } else {
      1
    }
  }
  lastInformativePosition <- function(z){
    seq.length <- length(z)
    z <-  which(rev(z) %in% as.raw(4))
    z <- z == 1:length(z)
    if (any(z)){
      seq.length - max(which(z)) 
    } else {
      seq.length
    }
  }
  h1 <- apply(msa1, 1, firstInformativePosition)
  h1 <- min(h1[shared])
  t1 <- apply(msa1, 1, lastInformativePosition)
  t1 <- max(t1[shared])
  
  h2 <- apply(msa2, 1, firstInformativePosition)
  h2 <- min(h2[shared])
  t2 <- apply(msa2, 1, lastInformativePosition)
  t2 <- max(t2[shared])
  
  both_leading <- (h1 > 1) & (h2 > 1)
  both_trailing <- (t1 < ncol(msa1)) & (t2 < ncol(msa2))
  
  if (both_trailing) {
    msa_trailing <- mafft.merge(list(msa1[, (t1 + 1):ncol(msa1)], msa2[, (t2 + 1):ncol(msa2)]), 
                               exec = exec)
    msa_trailing <- msa_trailing[!duplicated(rownames(msa_trailing)), ]
    msa1 <- msa1[, 1:t1]
    msa2 <- msa2[, 1:t2]
  }
  if (both_leading){
    msa_leading <- mafft.merge(list(msa1[, 1:(h1 - 1)], msa2[, 1:(h2 - 1)]), 
                               exec = exec)
    msa_leading <- msa_leading[!duplicated(rownames(msa_leading)), ]
    msa1 <- msa1[, h1:ncol(msa1)]
    msa2 <- msa2[, h2:ncol(msa2)]
  }

  ## identify non-empty positions in shared-subset
  ## -------------------------------------------
  shared1 <- msa1[shared, ]
  nonempty1 <- apply(shared1, 2, function(z) sum(z == 4))
  nonempty1 <- which(nonempty1 < nrow(shared1))
  
  shared2 <- msa2[shared, ]
  nonempty2 <- apply(shared2, 2, function(z) sum(z == 4))
  nonempty2 <- which(nonempty2 < nrow(shared2))
  
  ## Where do gaps have to be inserted?
  ## ----------------------------------
  insert1 <- insert2 <- vector(mode = "list")
  total <- length(nonempty1)
  for (i in 1:total){
    if (nonempty1[i] == nonempty2[i]) next
    if (nonempty1[i] > nonempty2[i]){
      d <- nonempty1[i] - nonempty2[i]
      insert2 <- c(insert2, list(c(nonempty2[i], d)))
      nonempty2[i:total] <- nonempty2[i:total] + d
    }
    if (nonempty1[i] < nonempty2[i]){
      d <- nonempty2[i] - nonempty1[i]
      insert1 <- c(insert1, list(c(nonempty1[i], d)))
      nonempty1[i:total] <- nonempty1[i:total] + d
    }
  }

  ## Prepare score matrices
  ## ----------------------
  if (include_scores){
    scores1 <- x[[id[1]]]@scores
    scores2 <- x[[id[2]]]@scores
  }
  
  ## Insert gaps into msa1 (and scores1)
  ## -----------------------------------
  if (length(insert1)){
    for (i in 1:length(insert1)){
      where <- insert1[[i]][1]
      what <- insert1[[i]][2]
      msa1 <- cbind(msa1[, 0:(where - 1)],
                    as.DNAbin(matrix("-", nrow = nrow(msa1), ncol = what, dimnames = list(rownames(msa1), NULL))),
                    msa1[, where:ncol(msa1)])
      if (include_scores){
        scores1 <- cbind(scores1[, 0:(where - 1)],
                         matrix(NaN, nrow = nrow(scores1), ncol = what, dimnames = list(rownames(scores1), NULL)),
                         scores1[, where:ncol(scores1)])
      }
    }
  }
  
  ## Insert gaps into msa2 (and scores2)
  ## -----------------------------------
  if (length(insert2)){
    for (i in 1:length(insert2)){
      where <- insert2[[i]][1]
      what <- insert2[[i]][2]
      msa2 <- cbind(msa2[, 0:(where - 1)],
                    as.DNAbin(matrix("-", nrow = nrow(msa2), ncol = what, dimnames = list(rownames(msa2), NULL))),
                    msa2[, where:ncol(msa2)])
      if (include_scores){
        scores2 <- cbind(scores2[, 0:(where - 1)],
                         matrix(NaN, nrow = nrow(scores2), ncol = what, dimnames = list(rownames(scores2), NULL)),
                         scores2[, where:ncol(scores2)])
      }
    }
  }

  ## add trailing gaps
  ## -----------------
  if (ncol(msa1) < ncol(msa2)){
    d <- ncol(msa2) - ncol(msa1)
    msa1 <- cbind(msa1,
                  as.DNAbin(matrix("-", nrow = nrow(msa1), ncol = d, dimnames = list(rownames(msa1), NULL))))
    if (include_scores){
      scores1 <- cbind(scores1,
                       matrix(NaN, nrow = nrow(scores1), ncol = d, dimnames = list(rownames(scores1), NULL)))
    }
  }
  if (ncol(msa1) > ncol(msa2)){
    d <- ncol(msa1) - ncol(msa2)
    msa2 <- cbind(msa2,
                  as.DNAbin(matrix("-", nrow = nrow(msa2), ncol = d, dimnames = list(rownames(msa2), NULL))))
    if (include_scores){
      scores2 <- cbind(scores2,
                       matrix(NaN, nrow = nrow(scores2), ncol = d, dimnames = list(rownames(scores2), NULL)))
    }
  }

  ## Remove one set of shared species and join bot sets
  ## --------------------------------------------------
  msa1 <- msa1[!rownames(msa1) %in% shared, ]
  if (include_scores){
    scores1 <- scores1[!rownames(scores1) %in% shared, ]
  }
  msa <- rbind(msa1, msa2)
  
  ## Reappend profile alignments of leading and trailing positions
  ## -------------------------------------------------------------
  if (both_leading){
    msa_leading <- msa_leading[match(rownames(msa), rownames(msa_leading)), ]
    msa <- cbind(msa_leading, msa)
  }
  if (both_trailing){
    msa_trailing <- msa_trailing[match(rownames(msa), rownames(msa_trailing)), ]
    msa <- cbind(msa, msa_trailing)
  }
  
  ## Create output object
  ## --------------------
  if (include_scores){
    polentaDNA(msa = msa,
               scores = rbind(scores1, scores2),
               method = x[[id[1]]]@method)
  } else {
    msa
  }
}
