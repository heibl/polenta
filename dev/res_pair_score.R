Rcpp::sourceCpp("dev/res_pair_score_dev.cpp")
library(polenta)
ref <- rbind(c("A", "A", "A", "G"),
  c("C","A", "A", "G" ),
  c("A", "-", "A", "G"))
com1 <- rbind(c("A", "A", "A", "G"),
  c("C","A", "A", "G" ),
  c("A", "A", "-", "G"))
com2 <- rbind(c("A", "A", "A", "G", "-"),
  c("C","A", "A", "G", "-"),
  c("A", "-", "A", "-", "G"))

ref <- as.DNAbin(ref)
com1 <- as.DNAbin(com1)
com2 <- as.DNAbin(com2)
coms <- list(com1, com2)



## single
res_pair_score_single(ref, alt = com1)
refc <- Cmatrix(ref)
altc <- cmatrix(com1)
rps <- add_msa(ref = refc, com = altc)

# for list
RPS(ref, coms)


## only 1 error...
ref <- rbind(c("A", "A", "-", "C", "G"),
  c("A","-", "A","C", "G" ),
  c("A", "-", "A", "C", "G"))

com <- rbind(c("A", "A", "-", "C", "G"),
  c("A","A", "-","C", "G" ),
  c("A", "-", "A", "C", "G"))
ref <- as.DNAbin(ref)
com <- as.DNAbin(com)
RPS(ref, com)
##korrekt


## variable dimensions
ref <- rbind(c("A", "A", "-", "C", "G"),
  c("A","-", "A","C", "G" ),
  c("A", "-", "A", "C", "G"))

com <- rbind(c("A", "A", "-", "C", "-", "G"),
  c("A","A", "-","C", "-", "G" ),
  c("A", "-", "A", "C", "G", "-"))
ref <- as.DNAbin(ref)
com <- as.DNAbin(com)
RPS(ref, com)
## korrekt

