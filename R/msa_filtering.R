#' Filtering of the base MSA using polenta scores
#'
#' @param polenta object of class \code{\link{polenta}}
#' @param col.cutoff numeric between 0 and 1; removes unreliable columns below the cutoff (default: 0.2); ignored if FALSE
#' @param seq.cutoff numeric between 0 and 1; removes unreliable sequences below the cutoff (default: 0.1); ignored if FALSE
#' @param mask.cutoff residues below the cutoff are masked ('N' for DNA, 'X' for AA; default: 0.5); ignored if FALSE
#' @param filter.ends logical, if TRUE trim.ends (ips) is applied to the MSA
#' @param filter.gaps logical, if TRUE trim.gabs (ips) is applied to the MSA
#' @param column_score logical, if TRUE column score (e.g. for RAxML: flag -a) is in the output
#' @param flag_a character specifying a path. If path is supplied function writes the filtered MSA into a fasta file. Additionally the function produces a file with the column score ready for RAxML input (flag -a)
#' @return masked MSA of class \code{AAbin} or \code{DNAbin}
#' @return column_score is optional
#' @seealso \code{\link{scores}}
#' @author Franz-Sebastian Krah
#' @export

filterMSA <- function(polenta,
  col.cutoff = 0.2,
  seq.cutoff = 0.1,
  mask.cutoff = 0.5,
  filter.ends = FALSE,
  filter.gaps  = FALSE,
  column_score = FALSE,
  flag_a = FALSE,
  na.coding = 0.5) {

  base.msa <- polenta@msa

  if (!mask.cutoff == FALSE) {
    r_sc <- scores(polenta, score = "residue", na.rm = FALSE)

    if (inherits(base.msa, "AAbin")) {
      base.msa <- as.character(base.msa)
      base.msa[r_sc$residue < mask.cutoff & base.msa != "-"] <- "X"
      base.msa <- as.AAbin(base.msa)
    }
    if (inherits(base.msa, "DNAbin")) {
      base.msa <- as.character(base.msa)
      base.msa[r_sc$residue < mask.cutoff & base.msa != "-"] <- "N"
      base.msa <- as.DNAbin(base.msa)
    }
  }

  if (!col.cutoff == FALSE) {
    g_sc <- scores(polenta, score = "column", na.rm = FALSE)
    # base.msa <- as.character(base.msa)

    base.msa <- base.msa[, c(which(g_sc$column$score >= col.cutoff), which(is.na(g_sc$column$score)))]
    g_sc <- g_sc$column[c(which(g_sc$column$score >= col.cutoff), which(is.na(g_sc$column$score))), ]
  }

  if (!seq.cutoff == FALSE) {
    s_sc <- scores(polenta, score = "sequence")
    base.msa <- base.msa[s_sc$sequence$score >= seq.cutoff,]
  }

  if (filter.ends) {
    keep <- trimEnds(base.msa)[[2]]
  }

  if (filter.gaps) {
    keep2 <- deleteGaps(base.msa)[[2]]
  }

  if (filter.ends | filter.gaps) {
    keep <- mget(c("keep", "keep2"), ifnotfound = list(NULL, NULL))
    keep <- Reduce(intersect, keep)
    base.msa <- base.msa[, keep]

    g_sc <- scores(g_r, score = "column", na.rm = FALSE)
    g_sc <- g_sc$column[keep, ]
    g_sc[is.na(g_sc$score),]$score <- na.coding
  }

  if (!flag_a == FALSE) {
    write.fas(base.msa,
      file = paste(flag_a, "/filtMSA_", Sys.Date(), ".fas", sep = ""))
    g_sc <- paste(round(g_sc$score * 10, digits = 0),
      collapse = " ")
    write(
      g_sc,
      file = paste(flag_a, "/filtMSA_GCSC_", Sys.Date(), ".txt", sep = ""),
      sep = "\t"
    )
  } else{
    if (column_score) {
      return(list(msa = base.msa, column_score = g_sc))
    } else{
      return(base.msa)
    }
  }
}
