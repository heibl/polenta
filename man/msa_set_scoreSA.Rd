% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msa_set_scoreSA.R
\name{msa_set_scoreSA}
\alias{msa_set_scoreSA}
\title{Compare Reference MSAs with Alternative MSAs}
\usage{
msa_set_scoreSA(ref, alt,
  exec = "/Users/krah/Documents/R/pkgs/polenta/src/msa_set_score_src/msa_set_score",
  bootstrap)
}
\arguments{
\item{ref}{of class data.frame, is the reference MSA ('BASE MSA') with
sequences as columns}

\item{alt}{path to alternative files}

\item{exec}{XXX.}

\item{bootstrap}{XXX.}
}
\value{
list containing following scores:

mean_scores: residue pair score and mean column score

column_score: identically aligned columns are asigned 1, otherwise 0;
  if more than one alternative MSA is supplied then mean of this score

residue_column_score: if one alternative MSA is supplied, then this
  is the SPC (sum-of-pairs column score); if more than one this is the
  average of the SPCs of all MSAs. The SPC is the mean over all residue pair
  scores (see there) in a column. This is also refered to as the GUIDANCE
  score.

residue_pair_score: if one alternative MSA is supplied then this is 1
  if a residue pair was identically aligned as in the reference MSA and 0
  otherwise. If more than one, then the average of the residue pair scores
  that result from each comparison with the reference MSA.

residue_pair_residue_score: calculated by averaging the residue
  column score over all pairs that have this residue. This is a confidence
  score for each residue.

residual_pair_sequence_pair_score_ if one alternative MSA is supplied
  this compares residue pairs of all sequence pairs of the base MSA with the
  reference MSA and returns 1 if identical and 0 if not and then averages
  over all pairs of a sequence pair. If more than one the the mean over those
  averages are computed.

residual_pair_sequence_score: like the residue_column_score but for
  sequences
}
\description{
Wrapper function for program msa_set_score v2.01 of the GUIDANCE
  program (see reference). Copyright: To modify the code, or use parts of it
  for other purposes, permission should be requested. Please contact Tal
  Pupko: talp@post.tau.ac.il. Please note that the use of the GUIDANCE
  program is for academic use only. C code computing basic MSA comparision.
  The most basic is the residue pairs residue score, which checks if residue
  pairs combinations are correctly aligned in both MSAs. From this the
  residue score, residue column score (GUIDANCE score), residue sequence
  score are computed. It also calulates the column score (CS), which simply
  checks if a column is identically aligned in the alternative MSA.
}
\references{
Penn et al. (2010). An alignment confidence score capturing
  robustness to guide tree uncertainty. Molecular Biology and Evolution
  27:1759--1767.
}
\seealso{
\code{\link{guidance}}, \code{\link{guidance2}}, \code{\link{HoT}}
}
\author{
Franz-Sebastian Krah
}
