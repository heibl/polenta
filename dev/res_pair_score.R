library(polenta)
# library("RcppArmadillo")
source("R/gapbase_bin.R")
Rcpp::sourceCpp("src/cmatrix.cpp")
Rcpp::sourceCpp("src/res_pair_score.cpp")

nChoosek(4,2)
which_true(c(0,0,1)==1)

#
ref <- rbind(c("A", "A", "A", "G"),
  c("C","A", "A", "G" ),
  c("-A", "-", "A", "G"))
ref <- rbind(c("A", "A", "A", "G"),
  c("C","A", "A", "G" ),
  c("-A", "A", "-", "G"))
com2 <- rbind(c("A", "A", "C", "G"),
  c("C","A", "-", "-" ),
  c("-", "-", "A", "G"))
coms <- list(com1, com2)


## single
res_pair_score_single(ref, alt = com1)
refc <- Cmatrix(gbbin(ref))
altc <- Cmatrix(gbbin(com1))
rps <- add_msa(ref = refc, com = altc)

# for list
RPS(ref, com2)


## only 1 error...
ref <- rbind(c("A", "A", "-", "C", "G"),
  c("A","-", "A","C", "G" ),
  c("A", "-", "A", "C", "G"))

com <- rbind(c("A", "A", "-", "C", "G"),
  c("A","A", "-","C", "G" ),
  c("A", "-", "A", "C", "G"))
RPS(ref, com)
##korrekt


## variable dimensions
ref <- rbind(c("A", "A", "-", "C", "G"),
  c("A","-", "A","C", "G" ),
  c("A", "-", "A", "C", "G"))

com <- rbind(c("A", "A", "-", "C", "-", "G"),
  c("A","A", "-","C", "-", "G" ),
  c("A", "-", "A", "C", "G", "-"))
RPS(ref, com)


