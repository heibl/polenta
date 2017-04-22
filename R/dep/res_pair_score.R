 # Calculates residul pairs score
 # @export

res_pair_score_single <- function(ref, alt){

  ## Recode MSAs to Cmatrix
  #-------------------------
  ref <- (ref != "-")*1
  alt <- (alt != "-")*1

  ref <- Cmatrix(ref)
  alt <- Cmatrix(ref)

  ## Initialize Matrix and counter
  #-------------------------------
  rps <- matrix(NA, ncol = 4, nrow = choose(3,2)*ncol(ref))
  colnames(rps) <- c("col", "row1", "row2", "score")
  count = 0
  ## Calculate scores
  #-------------------
  for(col in 1:ncol(ref)){
    for(row1 in 1:(nrow(ref)-1)){
      for(row2 in (row1+1):nrow(ref)){
        count <- count + 1
        ## reference residue pair
        ref_rp <- ref[c(row1, row2), col]
        # if gap => NA
        if(ref_rp[1] %% 2 == 0 | ref_rp[2] %% 2 == 0){
          rps[count,] <- c(col = col, row1 = row1, row2 = row2, score = NA)
        }else{ # no gap
          # if the two aligned residues from the REF
          # are aligned in one column in the ALT => 1, else 0
          if(which(alt[row1,]==ref_rp[1]) == which(alt[row2,]==ref_rp[2])){
            rps[count,] <- c(col = col, row1 = row1, row2 = row2, score = 1)
          }else{
            rps[count,] <- c(col = col, row1 = row1, row2 = row2, score = 0)
          }
        }
      } # row2
    } # row1
  } # col
  return(rps)
}

## for multiple alt
res_pair_score <- function(ref, alt){ # where coms will later be called from files

  if(inherits(alt, "list")){
    res <- res_pair_score_single(ref, alt[[1]])
    for(i in 1:length(alt)){
      res[,4] <- res[,4] + res_pair_score_single(ref, alt[[i]])[,4]
    }
  }else{
    res <- res_pair_score_single(ref, alt)
  }
  res[,4] <- res[,4]/length(alt)
  return(res)
}


## for Rcpp
res_pair_score_singleC <- function(ref, alt){
  ## Recode MSAs to Cmatrix
  #-------------------------
  ref <- (ref != "-")*1
  alt <- (alt != "-")*1

  ref <- Cmatrix(ref)
  alt <- Cmatrix(ref)

  if(inherits(alt, "list")){
    res <- add_msa(ref = ref, com = alt[[i]])
    for(i in 1:length(alt)){
      res[,4] <- res[,4] + add_msa(ref = ref, com = alt)[,4]
    }
  }else{
    res <- add_msa(ref = ref, com = alt)
  }
  res[,4] <- res[,4]/length(alt)
  return(res)
}
