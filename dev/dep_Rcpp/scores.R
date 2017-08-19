#' Calculate additional scores from the residue pair score
#'
#' @export
#'
#'

scores <- function(rps){

## Residue pairs residue score
count = 1
rprs <- matrix(NA, ncol = 3, nrow =  max(rps[,1])*max(rps[,3]))
for(col in 1:max(rps[,1])){
  for(row in 1:max(rps[,3])){
    rprs[count,] <- c(col, row, mean(rps[rps[,1]==col &
        ((rps[,2]== row) | (rps[,3]== row)), 4], na.rm = TRUE))
    count <- count + 1
    # print(rprs)
  }
}

## Residue pairs column score

rpcs <- matrix(NA, ncol = 2, nrow = max(rps[,1]))
for(col in 1:max(rps[,1])){
  rpcs[col,] <- c(col, mean(rps[rps[,1]==col, 4], na.rm = T))
}

## Residue pairs sequence score
rpsc <- matrix(NA, ncol = 2, nrow = max(rprs[,2]))
for(row in 1:max(rprs[,2])){
  rpsc[row,] <- c(row, mean(rprs[rprs[,2]==row, 3], na.rm = T))
}

res <- list(rprs, rpcs, rpsc)
return(res)
}
