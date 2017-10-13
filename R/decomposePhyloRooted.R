## This code is part of the megaptera package
## © C. Heibl 2016 (2017-10-11)

#' @title Decompose Phylogeny
#' @description Decomposes (splits) a phylogenetic tree into roughly equal-sized subtree.
#' @param phy An object of class \code{\link{phylo}}.
#' @param k An integer giving the size of the subtrees.
#' @importFrom ape is.binary.tree multi2di
#' @export

decomposePhyloRooted <- function(phy, k = 200){

  if (!is.binary.tree(phy)) phy <- multi2di(phy, random = FALSE)

  rt <- Ntip(phy) + 1

  dsc <- phy$edge[phy$edge[, 1] == rt, 2]
  subtrees <- lapply(dsc, descendants, phy = phy)
  names(subtrees) <- dsc

  st.sizes <- sapply(subtrees, length)
  while (any(st.sizes > k)){

    id <- which(st.sizes > k)[1]

    dsc <- phy$edge[phy$edge[, 1] == names(id), 2]
    sbt <- lapply(dsc, descendants, phy = phy)
    names(sbt) <- dsc
    subtrees <- c(subtrees[-id], sbt)
    st.sizes <- sapply(subtrees, length)
  }
  subtrees <- lapply(subtrees, function(x, phy) phy$tip.label[x], phy = phy)
  st.size <- sapply(subtrees, length)
  st.label <- paste("st", 1:length(subtrees), sep = "-")
  names(subtrees) <- st.label
  z <- rep(st.label, st.size)
  z <- data.frame(subtree = z,
                  taxon = unlist(subtrees),
                  stringsAsFactors = FALSE)
  rownames(z) <- NULL
  z <- z[match(z[, 2], phy$tip.label), ]
  phy$tip.label <- z$subtree
  pruner <- function(phy, subtree){
    which(phy$tip.label %in% subtree)[-1]
  }
  w <- unlist(lapply(unique(z$subtree), pruner, phy = phy))
  p <- drop.tip(phy, w)
  list(guidetree = p,
       subtrees = z)
}
