## This code is part of the polenta package
## © Franz-S. Krah 2017 (last update 2017-08-21)

setOldClass("AAbin")
setClass("polentaAA",
  representation = list(
    msa = "AAbin",
    scores = "matrix",
    method = "character")
)
