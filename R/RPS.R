#' Computes residue pair score scores by comparing a two MSAs
#' @param ref reference MSA, object of class .. must be object of \code{\link{gbbin.matrix}}
#' @param alt either single or many alternative MSA(s). If single class \code{matrix}, if multiple class \code{list}. Must be object of \code{\link{gbbin.matrix}} or \code{\link{gbbin.list}}
#' @author Franz-Sebastian Krah
#' @return residue_pair_score if one alternative MSA is supplied then this is 1
#'   if a residue pair was identically aligned as in the reference MSA and 0
#'   otherwise. If more than one, then the average of the residue pair scores
#'   that result from each comparison with the reference MSA.
#' @import Rcpp
#' @export

RPS <- function(ref, alt){

  if(!inherits(ref, c("DNAbin", "AAbin")))
    stop("msa is not of class 'DNAbin' or 'AAbin'")
  if(inherits(alt, "list")){
    if(any(!unlist(lapply(coms, inherits, what =  c("DNAbin", "AAbin")))))
      stop("alt not of class 'DNAbin' or 'AAbin'")
  }else{
    if(!inherits(alt, c("DNAbin", "AAbin")))
      stop("msa is not of class 'DNAbin' or 'AAbin'")
  }

  # test if same sequences
  if(inherits(alt, "list")){
  if(any(!unlist(lapply(alt, function(x)
    valid_alignments(as.character(ref), as.character(x))))))
    stop("the alignments do not share the same sequences")
  }else{
    valid_alignments(as.character(ref), as.character(alt))
  }


  ## Recode MSAs to Cmatrix
  #-------------------------
  ref <- cmatrix(ref)
  # print("cmat ok \n")

  if(inherits(alt, "list")){
    alt <- lapply(alt, cmatrix)
    # print("cmat ok")

    res <- add_msa(ref = ref, com = alt[[1]])
    for(i in 2:length(alt)){
      res[,4] <- res[,4] + add_msa(ref = ref, com = alt[[i]])[,4]
    }
    # print("add_msa ok")
  }else{
    alt <- cmatrix(alt)
    res <- add_msa(ref = ref, com = alt)
    # print("cmat ok")
    # print("add_msa ok")
  }
  if(inherits(alt, "list"))
    res[,4] <- res[,4]/length(alt)

  colnames(res) <- c("col", "row1", "row2", "score")
  return(res)
}

## taken from R package AlignStat, however, slightly modified
degap_alignment <- function(msa){
  # Remove gaps to convert alignment to list of strings
  gsub("-", "", do.call("paste", c(data.frame(msa),sep="")))
}
valid_alignments <- function(ref,alt){
  ref.degap <- degap_alignment(ref)
  alt.degap <- degap_alignment(alt)
  suppressWarnings(all(sort(ref.degap)==sort(alt.degap)))
}
