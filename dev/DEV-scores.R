head(scores$column_score)

id <- scores$residue_pair_column_score$col %in% scores$column_score$col
id <-  scores$column_score$col %in% scores$residue_pair_column_score$col

cs <- scores$column_score


x <- cbind(cs$CS[rpcs$col], rpcs$col_score)


sapply(scores, dim)
head(rprs)

## residue_pair_residue_score as matrix
rprs <- scores$residue_pair_residue_score
m <- matrix(rprs$score, nrow = nrow(base.msa))
# ID <- which(m < 1.0, arr.ind = TRUE)
# cs[55, ]
# rpcs[55, ]

## calculate residue_pair_column_score from m
rpcs <- colMeans(m, na.rm = TRUE)
rpcs <- rpcs[!is.na(rpcs)]
# id <- which(rpcs != scores$residue_pair_column_score$col_score)
# cbind(rpcs, scores$residue_pair_column_score$col_score)[id, ]

## calculate residue_pair_sequence_score from m
rpss <- rowMeans(m, na.rm = TRUE)
rpss <- rpss[!is.na(rpss)]
# id <- which(rpcs != scores$residue_pair_column_score$col_score)
# cbind(rpcs, scores$residue_pair_column_score$col_score)[id, ]





