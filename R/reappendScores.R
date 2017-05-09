## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-05-09)

#' @title Reappend Scores after Merging
#' @description Reappend scores to alignments after profile alignment.
#' @seealso \code{\link{polentaDNA}}
#' @export

reappendScores <- function(v, merged, scored){

  merged <- merged[[paste(v[1], v[2], sep = "-")]]

  ## gaps in S1
  s1 <- extractMSA(scored[[v[1]]])
  c1 <- merged[rownames(s1), ]
  id <- apply(c1, 2, function(z) sum(z == 4))
  id <- which(id == nrow(c1))
  scores_s1 <- scored[[v[1]]]@scores
  m1 <- matrix(NaN, nrow = nrow(s1), ncol = ncol(merged))
  m1[, -id] <- scores_s1
  rownames(m1) <- rownames(s1)

  ## gaps in S2
  s2 <- extractMSA(scored[[v[2]]])
  c2 <- merged[rownames(s2), ]
  id <- apply(c2, 2, function(z) sum(z == 4))
  id <- which(id == nrow(c2))
  scores_s2 <- scored[[v[2]]]@scores
  m2 <- matrix(NaN, nrow = nrow(s2), ncol = ncol(merged))
  m2[, -id] <- scores_s2
  rownames(m2) <- rownames(s2)

  m <- rbind(m1, m2)
  m <- m[match(rownames(merged), rownames(m)), ]

  meth <- unique(c(scored[[v[1]]]@method, scored[[v[2]]]@method))
  if (length(meth) > 1) stop("cannot combine scores from different methods")
  polentaDNA(msa = merged,
             scores = m,
             method = meth)

}
