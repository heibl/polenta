Rcpp::sourceCpp("src/res_pair_score.cpp")
Rcpp::sourceCpp("src/cmatrix_parallel.cpp")
Rcpp::sourceCpp("src/score_parallel.cpp")

library(polenta)
library(ape)
ref <- rbind(c("A", "A", "A", "G"),
  c("C","A", "A", "G" ),
  c("A", "-", "A", "G"))
com <- rbind(c("A", "A", "A", "G", "-"),
  c("C","A", "A", "G", "-"),
  c("A", "-", "A", "-", "G"))

ref_big <- do.call(rbind, replicate(70, ref, simplify = FALSE))
ref_big <- do.call(cbind, replicate(200, ref_big, simplify = FALSE))

com_big <- do.call(rbind, replicate(70, com, simplify = FALSE))
com_big <- do.call(cbind, replicate(200, com_big, simplify = FALSE))


## test for small
ref <- as.DNAbin(ref)
com <- as.DNAbin(com)

ref <- gbbin(as.character(ref))
com <- gbbin(as.character(com))
ref <- Cmatrix(ref)
com <- Cmatrix(com)

# parallel
system.time(
sc <- add_msa_sc(com, ref))
# seriel
add_msa_score(ref, com)


# test for big

ref_big <- as.DNAbin(ref_big)
com_big <- as.DNAbin(com_big)

ref_big <- gbbin(as.character(ref_big))
com_big <- gbbin(as.character(com_big))
ref_big <- Cmatrix(ref_big)
com_big <- Cmatrix(com_big)

# parallel
system.time(
  sc <- add_msa_sc(com_big, ref_big)
)
# seriel
## better not run...
system.time(
  cs2 <- add_msa_score(com_big, ref_big)
)
