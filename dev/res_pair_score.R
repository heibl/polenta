
#
# 1. loop über n alt
# 2. erstelle matrix 3 cols und GAUSSUMM * cols language
# 3. loop über REF cols
# 4. loop über res pairs (2 loops => 1 für erste row, 2te für 2te row)
# 5. compare res pair from pos X in BASE with unaligned pos X in ALT

Rcpp::sourceCpp("src/cmatrix.cpp")


ref <- rbind(c("A", "A", "A", "G"),
  c("C","-", "A", "-" ),
  c("-", "-", "A", "G"))

com1 <- rbind(c("-", "A", "A", "G"),
  c("C","-", "A", "-" ),
  c("-", "A", "-", "G"))

com2 <- rbind(c("A", "A", "-", "G"),
  c("C","A", "-", "-" ),
  c("-", "-", "A", "G"))

# coms <- list(com1, com2)

com <- com2
res_pair_score_single(ref, com2)

res_pair_score_single <- function(ref, com){

  # Recode MSAs to Cmatrix
  ref <- (ref!="-")*1
  com <- (com!="-")*1
  ref <- Cmatrix(ref)
  com <- Cmatrix(com)


  rps <- matrix(NA, ncol = 4, nrow = choose(3,2)*ncol(ref))
  colnames(rps) <- c("col", "row1", "row2", "score")
  count = 0
  for(col in 1:ncol(ref)){
    for(row1 in 1:(nrow(ref)-1)){
      for(row2 in (row1+1):nrow(ref)){
        count <- count + length((row2*0))
        print(count)
        ## reference residue pair
        ref_rp <- ref[c(row1, row2), col]
        # if gap => NA
        if(ref_rp[1] %% 2 == 0 | ref_rp[2] %% 2 == 0){
          rps[count,] <- c(col = col, row1 = row1, row2 = row2, score = NA)
        }else{ # no gap
          # if the two aligned residues from the REF
          # are aligned in one column in the ALT => 1, else 0
          if(which(com[row1,]==ref_rp[1]) == which(com[row2,]==ref_rp[2])){
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

res_pair_score <- function(ref, coms){ # where coms will later be called from files

  if(length(coms)>1){
    res <- res_pair_score_single(ref, coms[[1]])
    for(i in 1:length(coms)){
      res[,4] <- res[,4] + res_pair_score_single(ref, coms[[i]])[,4]
    }
  }else{
    res <- res_pair_score_single(ref, coms)
  }
  return(res)
}

res_pair_score(ref, coms)
