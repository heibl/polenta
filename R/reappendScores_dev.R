## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-11-08)

#' @title Reappend Scores after Merging
#' @description Reappend scores to alignments after profile alignment.
#' @param v XXX
#' @param merged XXX
#' @param scored XXX
#' @seealso \code{"\link[=polentaDNA-class]{polentaDNA}"}
#' @export

reappendScores_dev <- function(seqs, merged){
  
  merged <- merged[[paste(v[1], v[2], sep = "-")]]
  
  ## The idea here is that the position of each residue relativ
  # to all other residues stays the same. Thus we can reshape::melt
  # the scores and simply rbind the columns. Then we also have the
  # score for each residue in the correct position.
  
  # Melt scores
  s1 <- scores(seqs$S1, score = "residue", na.rm = FALSE)
  s2 <- scores(seqs$S2, score = "residue", na.rm = FALSE)
  
  s1 <- reshape::melt(s1$residue)
  s1$X1 <- as.numeric(factor(s1$X1))
  s2 <- reshape::melt(s2$residue)
  s2$X1 <- as.numeric(factor(s2$X1))
  
  # Melt bases (not neccessary, but to have it to check if working correctly)
  m1 <- reshape::melt(as.character(seqs$S1@msa))
  m1$X1 <- as.numeric(factor(m1$X1))
  m2 <- reshape::melt(as.character(seqs$S2@msa))
  m2$X1 <- as.numeric(factor(m2$X1))
  
  # Melt the 
  merg <- reshape::melt(as.character(merged))
  merg$X1 <- as.numeric(factor(merg$X1))
  
  ms1 <- cbind(m1, s1)
  ms2 <- cbind(m2, s2)
  
  # combine the melted score matrices
  sc_merge <- foreach(col = unique(merg$X2), 
    .combine = 'rbind') %do% {
      print(col)
      
      ms1_col <- ms1[ms1$X2==col,]
      ms2_col <- ms2[ms2$X2==col,]
      tot.seq <- nrow(merg[merg$X2==col,])
      if(nrow(ms1_col)==0 | nrow(ms2_col)==0){
        if(nrow(ms1_col)==0){
          fill <- data.frame(matrix(NA, ncol = 6, nrow = (tot.seq - nrow(ms2_col))))
          names(fill) <- colnames(ms1_col)
          ms1_col <- fill
        }
        if(nrow(ms2_col)==0){
          fill <- data.frame(matrix(NA, ncol = 6, nrow = (tot.seq - nrow(ms1_col))))
          names(fill) <- colnames(ms2_col)
          ms2_col <- fill
        }
      }
      ms1ms2 <- rbind(ms1_col, ms2_col)
      
      
      if(nrow(ms1ms2)==0){
        fill <- data.frame(matrix(NA, ncol = 6))
        names(fill) <- colnames(ms1ms2)
        cbind(fill, merg[merg$X2==col,])
      }else{
        cbind(ms1ms2, merg[merg$X2==col,])
      }
    }
  
  sc_merge <- data.frame(col = sc_merge[,8], 
    row = sc_merge[,7], 
    base = sc_merge[,9], 
    score = sc_merge[,6])
  
  sc_merge <- reshape2::dcast(sc_merge, row ~ col, value.var="score")
  sc_merge <- sc_merge[,-1] # delete first col which is the row values
  sc_merge <- as.matrix(sc_merge)
  
  meth <- unique(c(seqs$S1@method, seqs$S2@method))
  if (length(meth) > 1) stop("cannot combine scores from different methods")
  
  polentaDNA(msa = merged,
    scores = sc_merge,
    method = meth)
  
}





