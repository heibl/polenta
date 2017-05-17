## This code is part of the polenta package
## © C. Heibl 2016 (last update 2017-05-17)

#' @title Transitivity Merge
#' @description Merges two multiple sequence alignments and their reliability scores using the transitivity
#'   criterion.
#' @param x A list containing objects of class \code{\link{polentaDNA}}.
#' @return An object of class \code{\link{polentaDNA}}.
#' @importFrom ape as.DNAbin cbind.DNAbin rbind.DNAbin
#' @export

transitivityMerge <- function(x, id = c(1, 2)){

  msa1 <- x[[id[1]]]@msa
  msa2 <- x[[id[2]]]@msa

  ## shared species
  shared <- intersect(rownames(msa1), rownames(msa2))
  if (!length(shared)) stop("MSAs do not overlap")

  ## identify all-gap positions in shared-subset
  ## -------------------------------------------
  shared1 <- msa1[shared, ]
  gaps1 <- apply(shared1, 2, function(z) sum(z == 4))
  nongaps1 <- which(gaps1 < nrow(shared1))
  gaps1 <- which(gaps1 == nrow(shared1))

  shared2 <- msa2[shared, ]
  gaps2 <- apply(shared2, 2, function(z) sum(z == 4))
  nongaps2 <- which(gaps2 < nrow(shared2))
  gaps2 <- which(gaps2 == nrow(shared2))

  ## Where do gaps have to be inserted?
  ## -----------------------------------
  ng1 <- nongaps1; ng2 <- nongaps2
  insert1 <- vector(mode = "list")
  insert2 <- vector(mode = "list")
  total <- length(ng1)
  for (i in 1:total){
    if (ng1[i] == ng2[i]) next
    if (ng1[i] > ng2[i]){
      d <- ng1[i] - ng2[i]
      insert2 <- c(insert2, list(c(ng2[i], d)))
      ng2[i:total] <- ng2[i:total] + d
    }
    if (ng1[i] < ng2[i]){
      d <- ng2[i] - ng1[i]
      insert1 <- c(insert1, list(c(ng1[i], d)))
      ng1[i:total] <- ng1[i:total] + d
    }
  }

  ## Insert gaps into sequence and score matrices
  ## --------------------------------------------
  scores1 <- x[[id[1]]]@scores
  scores2 <- x[[id[2]]]@scores
  for (i in 1:length(insert1)){
    where <- insert1[[i]][1]
    what <- insert1[[i]][2]
    msa1 <- cbind(msa1[, 0:(where - 1)],
          as.DNAbin(matrix("-", nrow = nrow(msa1), ncol = what, dimnames = list(rownames(msa1), NULL))),
          msa1[, where:ncol(msa1)])
    scores1 <- cbind(scores1[, 0:(where - 1)],
                    matrix(NaN, nrow = nrow(scores1), ncol = what, dimnames = list(rownames(scores1), NULL)),
                    scores1[, where:ncol(scores1)])
  }
  for (i in 1:length(insert2)){
    where <- insert2[[i]][1]
    what <- insert2[[i]][2]
    msa2 <- cbind(msa2[, 0:(where - 1)],
                  as.DNAbin(matrix("-", nrow = nrow(msa2), ncol = what, dimnames = list(rownames(msa2), NULL))),
                  msa2[, where:ncol(msa2)])
    scores2 <- cbind(scores2[, 0:(where - 1)],
                     matrix(NaN, nrow = nrow(scores2), ncol = what, dimnames = list(rownames(scores2), NULL)),
                     scores2[, where:ncol(scores2)])
  }

  ## add trailing gaps
  ## -----------------
  if (ncol(msa1) < ncol(msa2)){
    d <- ncol(msa2) - ncol(msa1)
    msa1 <- cbind(msa1,
                  as.DNAbin(matrix("-", nrow = nrow(msa1), ncol = what, dimnames = list(rownames(msa1), NULL))))
    scores1 <- cbind(scores1,
                     matrix(NaN, nrow = nrow(scores1), ncol = what, dimnames = list(rownames(scores1), NULL)))
  }
  if (ncol(msa1) > ncol(msa2)){
    d <- ncol(msa1) - ncol(msa2)
    msa2 <- cbind(msa2,
                  as.DNAbin(matrix("-", nrow = nrow(msa2), ncol = what, dimnames = list(rownames(msa2), NULL))))
    scores2 <- cbind(scores2,
                     matrix(NaN, nrow = nrow(scores2), ncol = what, dimnames = list(rownames(scores2), NULL)))
  }

  ## remove one set of shared species ...
  msa1 <- msa1[!rownames(msa1) %in% shared, ]
  scores1 <- scores1[!rownames(scores1) %in% shared, ]
  ## ... and join bot sets
  msa <- rbind(msa1, msa2)
  scores <- rbind(scores1, scores2)
  polentaDNA(msa = rbind(msa1, msa2),
             scores = rbind(scores1, scores2),
             method = x[[id[1]]]@method)
}
