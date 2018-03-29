## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-11-08)

#' @title Reappend Scores after Merging
#' @description Reappend scores to alignments after profile alignment.
#' @param v XXX
#' @param merged XXX
#' @param scored XXX
#' @seealso \code{"\link[=polentaDNA-class]{polentaDNA}"}
#' @export

reappendScores_gen_dev <- function(v, s){
  

  ## resort merg-MSA to match individual MSAs rownames
  seq_nams <- unlist(lapply(v, function(x) rownames(x@msa)), use.names = F)
  s <- s[match(seq_nams, rownames(s)),]
  

  ## The idea here is that the position of each residue relativ
  # to all other residues stays the same. Thus we can reshape::melt
  # the scores and simply rbind the columns. Then we also have the
  # score for each residue in the correct position.
  
  # Caclulate mean residue score 
  # (won't work without because bases can moove relative to each other)
  sc <- lapply(v, scores, score = "residue", na.rm = FALSE)
  # Melt scores
  sc <- lapply(sc, function(x) reshape::melt(x$residue))
  
  # Melt bases (not neccessary, but to have it to check if working correctly)
  m <- lapply(v, function(x) reshape::melt(as.character(x@msa)))
  
  # Melt the merged alignment
  merg <- reshape::melt(as.character(s))
  names(merg) <- c("row", "col", "base")
  
  # Bind bases and scores
  msc <- foreach(l = 1:length(sc)) %do% {
    cbind(m[[l]], sc[[l]])
  }
  
  # new colnames
  nam <- c("row", "col", "base", "row", "col", "score")
  for(i in 1:length(msc)) names(msc[[i]]) <- nam
  
  # remove redundant colnames
  for(i in 1:length(msc)) msc[[i]] <- msc[[i]][,c(2,1,3,6)]
  nam <- names(msc[[i]])
  
  # rbind MSA-sites of the individual MSAs (melt-MSAs)
  msc_col <- foreach(col = unique(merg$col), .combine = "rbind") %do% {
    do.call(rbind, lapply(msc, function(x) x[x$col==col,]))
  }
  
  # delete gaps and NAs and cbind melt-MSAs and merg-MSA
  msc_col <- msc_col[!is.na(msc_col$base),]
  msc_col <- msc_col[!msc_col$base == "-",]
  merg <- merg[!merg$base == "-",]
  # cbind
  msc <- data.frame((msc_col), (merg))
  
  # reshape to seq x site matrix
  msc <- reshape2::dcast(msc, row.1 ~ col.1, value.var="score")
  msc <- msc[,-1] # delete first col which is the row values
  msc <- as.matrix(msc)
  
  meth <- Reduce(unique, lapply(v, function(x) x@method))
  if (length(meth) > 1) stop("cannot combine scores from different methods")
  
  polentaDNA(msa = s,
    scores = msc,
    method = meth)
  
}





