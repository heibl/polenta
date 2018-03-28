## This code is part of the polenta package
## © F.-S. Krah (last update 2017-11-13)

#' @title XXX
#' @description XXX
#' @param msa XXX
#' @param n.coopt XXX
#' @param type XXX
#' @param td XXX
#' @param files XXX
#' @param raw_seq XXX
#' @param msa.program XXX
#' @param method XXX
#' @param int_file XXX
#' @param msa.exec XXX
#' @importFrom ape compute.brlen multi2di nj Ntip root
#' @import foreach
#' @importFrom ips read.fas
#' @importFrom utils globalVariables
#' @export

Hot_GUIDANCE2 <- function(msa, n.coopt, type, td,
                          files, raw_seq, msa.program, method,
                          int_file, msa.exec){
  
  ## declare i to be a global variable; this is necessary because
  ## foreach uses non-standard evaluation that codetools is not aware of
  ## [http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html]
  ## Does not work [CH 2018-01-23]
  # globalVariables("i")

  # create start_tree
  li <- msa
  if (is.character(msa))
    msa <- read.fas(msa)
  msa.nj <- nj(dist.ml(as.phyDat(msa)))
  msa.nj <- root(msa.nj, outgroup = msa.nj$tip.label[1])
  msa.nj <- multi2di(msa.nj)
  msa.nj <- compute.brlen(msa.nj)

  ## produce MSA partitions
  align_parts <- partitions(msa.nj)

  ## sample 4 or n co-optimal
  nt <- Ntip(msa.nj)
  n.co <- sample((nt - 3) * 8, n.coopt)
  n.co_mat <- data.frame(n = 1:((nt - 3) * 8),
    part = rep(1:(nt - 3), each = 8),
    n.in.part = rep(1:8, times = (nt - 3)))
  n.co.sub <- n.co_mat[n.co_mat$n %in% n.co, ]

  # reduce partitions to the randomly choosen co-optimal MSA number
  align_parts <- align_parts[, n.co.sub$part]
  # number of random MSA within partition (remember, each partition has 8 alignments)
  n.co.sub <- n.co.sub$n.in.part
  
  # make the 4 or n alignments
  alt_msas <- foreach(i = 1:ncol(align_parts),
                      .export = c("align_parts")) %do% {
    align_part_set(x = raw_seq, partition_set = align_parts[, i],
      method = method, msa.exec = msa.exec, msa.program = msa.program,
      coopt.sub = n.co.sub[i])
  }

  ## unlist
  alt_msas <- foreach(i = 1:length(alt_msas), .combine = c) %do% {
    alt_msas[[i]]
  }

  ## write to files
  if (int_file){
    for (j in 1:n.coopt)
      write.fas(alt_msas[[j]], files[j])
  } else {
    return(alt_msas)
    }
}


