library(Rcpp)
library(polenta)

ref <- rbind(c("-", "-", "A", "A", "G", "T"),
  c("C","-", "-","A", "-" ,"T"),
  c("-", "A", "-","-", "G", "T"))
com1 <- rbind(c("-", "A", "A", "G", "T"),
  c("C","-", "A", "-" , "T"),
  c("-", "A", "-", "G", "T"))
com2 <- rbind(c("A", "A", "-", "G", "T"),
  c("C","A", "-", "-" , "G"),
  c("-", "-", "A", "G", "T"))

coms <- list(com1, com2)

# test
## does not yield the correct result yet..
## and R often quits unexpectedly
msa_RPscore(ref = ref, alt = com1)

