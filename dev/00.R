library("rpg")
library("ips")
library("parallel")
library("foreach")
library("phangorn")
library("doSNOW")
library("adephylo")
library("useful")
library("stringr")
library("scales")
library("ggplot2")
library("zoo")
library("cowplot")

library("polenta")
# 1.	Read example sequences
## DNA
seq_dna <- read.fas("dev/data/cortinarius_28s_ms.fas")
set.seed(100)
seq_dna <- sample(seq_dna, 10)
## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas")


msa.program <- "mafft"
msa.exec <- "/usr/local/bin/mafft"
# msa.program <- "clustalo"
# exec <- "/Applications/clustalo"
# msa.program <- "clustalw2"
# exec <- "/Applications/clustalw2"
# msa.program <- "muscle"
# exec <- "/Applications/muscle"

source("R/DEV-mafft.R")
source("R/guidance_dev.R")

sequences = seq_dna
msa.program = "mafft"
# exec = exec
bootstrap = 100
parallel = TRUE
ncore = 4
method = "retree 1"
nj.program = "R"
int_file = FALSE
score_method = "Rcpp"

system.time(
  g_r <- guidance(sequences = seq_aa,
    msa.program = "mafft",
    msa.exec = msa.exec,
    bootstrap = 100,
    parallel = TRUE, ncore = "auto",
    method = "retree 1",
    nj.program = "R")
)
gsc <- daughter_scores(g_r, score = c("gcsc", "rprsc"))



system.time(
  g_r2 <- guidance(sequences = seq_dna,
    msa.program = "mafft",
    # exec = exec,
    bootstrap = 100,
    parallel = TRUE, ncore = "auto",
    method = "retree 1",
    nj.program = "R")
)

system.time(
  g_sa <- guidanceSA(sequences = seq_dna,
    msa.program = "mafft",
    programm = "guidance",
    bootstrap = 100,
    proc_num = 4,
    exec = "/Applications/guidance.v2.02/")
)



sequences = seq_aa
msa.program = "mafft"
parallel = FALSE
ncore = 4
method = "retree 1"
plot_guide = TRUE

system.time(hot_r <- HoT(sequences = seq_aa,
  msa.program = "mafft",
  parallel = FALSE, ncore = 4,
  method = "retree 1",
  plot_guide = TRUE))

system.time(
  hot_sa <- guidanceSA(sequences = seq_aa,
    msa.program = "mafft",
    programm = "hot",
    bootstrap = 100,
    proc_num = 4,
    quiet = FALSE)
)
hot_r.h <- confidence.heatmap(hot_msa, title = "HoT SA", legend = FALSE,
  guidance_score = FALSE)
hot_sa.h <- confidence.heatmap(hot_sa, title = "HoT SA", legend = FALSE,
  guidance_score = FALSE)
plot_grid(hot_r.h, hot_sa.h, nrow = 2, ncol = 1)







system.time(
g2_r <- guidance2(sequences = seq_aa,
  msa.program = "mafft",
  # exec <- "/Applications/clustalw2",
  bootstrap = 100,
  parallel = TRUE, ncore ="auto",
  method = "auto",
  n.coopt = "auto")
)
heatmap.msa(obj = g2, file =paste(getwd(), "test.pdf", sep="/"))


system.time(g2_sa <- guidanceSA(sequences = seq_aa,
  msa.program = "mafft",
  outdir = "test.guidance",
  programm = "guidance2",
  bootstrap = 100,
  proc_num = 4))



file <- system.file("extdata", "BB50009.fasta", package = "polenta")
aa_seq<- read.fas(file)
g_res <- guidance(sequences = aa_seq)
scores <- daughter_scores(g_r, score = c("gcsc", "rprsc"))
hist(scores$gcsc$score, xlab = "Column score", main = "GUIDANCE")
