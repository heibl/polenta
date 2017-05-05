## This code is part of the polenta package
## © C. Heibl 2017 (last update 2017-05-05)

setOldClass("DNAbin")
setClass("polentaDNA",
         representation = list(
           msa = "DNAbin",
           scores = "matrix",
           method = "character")
)






