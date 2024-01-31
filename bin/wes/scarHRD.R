#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library("scarHRD")

score <- scar_score(args[1],reference = "grch38", seqz=TRUE)

filename <- paste(args[2], 'HRD_score.txt',sep='_')

write.table(score, file=filename, quote = FALSE,  sep = "\t", row.names = FALSE) 