library("polenta")
library("ape")

library(stringr)
library(foreach)
library(parallel)
library(doSNOW)
library(ips)



# Sequences --------------------------------------------------------------

## DNA
seq_dna <- read.FASTA("dev/data/cortinarius_28s_ms.fas")
set.seed(100)
seq_dna <- sample(seq_dna, 30)

## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas")



exec <- "/usr/local/bin/mafft"
# exec <- "/Applications/clustalo"
# exec <- "/Applications/clustalw2"
# exec <- "/Applications/muscle"
# exec <- "/Applications/guidance.v2.02/"

seq = seq_dna
msa.exec = exec
bootstrap = 100
ncore = 4
method = "retree 1"
score_method = "Rcpp"


# Guidance --------------------------------------------------------------

## R
system.time(
  g_r <- guidance(sequences = seq_dna, ncore = 4, bootstrap = 100, method = "retree 1")
)
gsc <- scores(g_r, score = "column", na.rm = FALSE)

msa_masked <- filterMSA(g_r,
  col.cutoff = 0.2,
  seq.cutoff = 0.1,
  mask.cutoff = 0.5,
  filter.ends = TRUE,
  filter.gaps  = TRUE,
  column_score = TRUE,
  flag_a = FALSE)
dim(g_r@msa)
dim(msa_masked$msa)

## SA
system.time(
  g_sa <- guidanceSA(sequences = seq,
    msa.program = "mafft",
    programm = "guidance",
    exec = exec,
    bootstrap = 100)
)


# HoT --------------------------------------------------------------

## R
system.time(hot_r <- HoT(sequences = seq))
hsc <- daughter_scores(hot_r, score = c("gcsc", "rprsc"))

## SA
system.time(
  hot_sa <- guidanceSA(sequences = seq,
    msa.program = "mafft",
    programm = "hot")
)


# Guidance 2 --------------------------------------------------------------

## R
system.time(
  g2_r <- guidance2(sequences = seq, ncore = "auto")
)
g2sc <- dauaghter_scores(g2_r, score = c("gcsc", "rprsc"))

## SA
system.time(g2_sa <- guidanceSA(sequences = seq,
  msa.program = "mafft",
  proc_num = 4)
)

# Polenta ----------------------------------------------------------------
# sample(seq_dna, 30)
system.time(
p_r <- polenta_dev(seqs = seq_dna, k = 200, bootstrap = 100,
  method = "retree 1", exec = exec, ncore = 4)
)
p_r_sc <- scores(p_r, score = "column", na.rm = FALSE)



## test GUIDANCE against POLENTA
seq_dna <- read.FASTA("../polenta_benchmarking/sate_big/500_test.fasta")
# seq_dna <- as.matrix(seq_dna)
dealign <- function(x){
  res <- lapply(as.character(x), function(x) x[-grep("-", x)])
  as.DNAbin(res)
}
seq_dna <-dealign(seq_dna)
seq_dna <- sample(seq_dna, 500)

system.time(
  g_r <- guidanceSA(sequences = seq_dna, 
    exec = "/Applications/guidance.v2.02/", 
    msa.program = "mafft", program = "guidance",
    bootstrap = 10, proc_num = 4)
)
gsc <- scores(g_r, score = "column", na.rm = FALSE)

system.time(
  p_r <- polenta_dev(seqs = seq_dna, k = 200, bootstrap = 10,
    method = "retree 1", exec = exec, ncore = 4)
)
p_r_sc <- scores(p_r, score = "column", na.rm = FALSE)


hist(p_r_sc$column$score, col = "red", density = 40, breaks = 100, xlab = "Column score")
hist(gsc$column$score, add = TRUE, col = "blue", density = 20, breaks = 100)
legend("topleft", c("GUIDANCE", "POLENTA", "nseq = 245"), text.col = c("blue", "red"), bty = "n")
