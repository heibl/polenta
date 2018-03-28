library("polenta")
library("ape")

# Sequences --------------------------------------------------------------

## DNA
seq_dna <- read.FASTA("dev/data/cortinarius_28s_ms.fas")
set.seed(100)
seq_dna <- sample(seq_dna, 10)

## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas")



exec <- "/usr/local/bin/mafft"
# exec <- "/Applications/clustalo"
# exec <- "/Applications/clustalw2"
# exec <- "/Applications/muscle"
exec <- "/Applications/guidance.v2.02/"

seq = seq_dna
msa.exec = exec
bootstrap = 100
ncore = 4
method = "retree 1"
score_method = "Rcpp"


# Guidance --------------------------------------------------------------

## R
system.time(
  g_r <- guidance(sequences = seq, ncore = "auto")
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
