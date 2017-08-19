Rcpp::sourceCpp("dev/score_calculation_DEV.cpp")
# Rcpp::sourceCpp("dev/add_msa_par.cpp")
Rcpp::sourceCpp("src/score_calculation.cpp")

#-- needed later, not important
gbbin <- function(msa){
  msa <- (msa != "-")*1
  return(msa)
}

library(ips)
# library(polenta)
library(ape)
ref <- rbind(spec1 = c("A", "A", "A", "G"),
  spec2 = c("C","A", "A", "G" ),
  spec3 = c("A", "-", "A", "G"),
  spec4 = c("C","A", "A", "G" ))
com <- rbind(spec1 = c("A", "A", "A", "G", "-"),
  spec2 = c("C","A", "A", "G", "-"),
  spec3 = c("-", "A", "A", "-", "G"),
  spec4 = c("C","A", "A", "G", "-"))

ref <- do.call(rbind, replicate(50, ref, simplify = FALSE))
ref <- do.call(cbind, replicate(200, ref, simplify = FALSE))
com <- do.call(rbind, replicate(50, com, simplify = FALSE))
com <- do.call(cbind, replicate(200, com, simplify = FALSE))


## test for small
ref <- as.DNAbin(ref)
com <- as.DNAbin(com)

# prep data
ref <- gbbin(as.character(ref))
com <- gbbin(as.character(com))
# produce C matrix
# ref <- Cmatrix(ref)
# com <- Cmatrix(com)

cmat_ref <- msa_recode_dev(ref)
cmat_alt <- msa_recode_dev(com)

ref_col2res <- cmat_ref$col2res
ref_col2res_odd <- ref_col2res %% 2
alt_col2res <- cmat_alt$col2res
alt_res2col <- cmat_alt$res2col

hit <- res_pair_hit_dev(cmat_ref$col2res)
# microbenchmark::microbenchmark(
for(i in 1:1000){
hit <- add_msa_dev(ref_col2res, alt_col2res, alt_res2col, ref_col2res_odd, hit)
}

hit2 <- res_pair_hit(cmat_ref$col2res)
microbenchmark::microbenchmark(
add_msa(ref_col2res, alt_col2res, alt_res2col, hit2)
)

for(i in 1:100){
hit <- add_msa_parallel(ref_col2res, alt_col2res, alt_res2col, ref_col2res_odd, hit)
}
hit[hit == -1] <- NA
hit <- hit/1


