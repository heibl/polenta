Rcpp::sourceCpp("src/res_pair_score.cpp")
library(polenta)
library(ape)
ref <- rbind(c("A", "A", "A", "G"),
  c("C","A", "A", "G" ),
  c("A", "-", "A", "G"))
com1 <- rbind(c("A", "A", "A", "G"),
  c("C","A", "A", "G" ),
  c("A", "A", "-", "G"))
com2 <- rbind(c("A", "A", "A", "G", "-"),
  c("C","A", "A", "G", "-"),
  c("A", "-", "A", "-", "G"))

# ref <- do.call(rbind, replicate(10, ref, simplify = FALSE))
# ref <- do.call(cbind, replicate(10, ref, simplify = FALSE))

ref <- as.DNAbin(ref)
com1 <- as.DNAbin(com1)
com2 <- as.DNAbin(com2)
coms <- list(com1, com2)

a <- c(0,0,0,1)
which_true2(a)

RPS(ref, alt = com1)
refc <- cmatrix(ref)
altc <- cmatrix(com1)
# add_msa(ref = refc, com = altc)
add_msa_score(ref = refc, com = altc)


library(ips)
seq_dna <- read.fas("dev/data/cortinarius_28s_ms.fas")
set.seed(100)
seq_dna <- sample(seq_dna,100)
al <- mafft(seq_dna)
al2 <- mafft(seq_dna, op = 5, ep = 0.2)

al3 <- list(al2, al2, al2)

RPS(ref = al, alt = al2)

al <- al[1:10, 1:10]
rownames(al) <- NULL
class(as.character(al))
attributes(al) <- attributes(al)[-2]
al

msa <-  gbbin(as.character(al))
rownames(msa) <- NULL
msa <- Cmatrix(msa)

refc <- cmatrix(al)
altc <- cmatrix(al2)
add_msa_vec(ref = refc, com = altc)


mat <- matrix(1, ncol = 1000, nrow = 200)

com_mat_maker(mat)

com_mat_maker <- function(msa){
  nr <- nrow(msa)
  le <- nChoosek(nr, 2)
  mat <- matrix(NA, ncol = 3, nrow = le*ncol(msa))
  count = 0
  for(col in 1:ncol(msa)){
    for(row1 in 1:(nrow(msa)-1)){
      for(row2 in (row1+1):(nrow(msa))){
        count <- count + 1
        mat[count, ] <- c(col, row1, row2)
      }
    }
  }
  return(mat)
}

## single


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

