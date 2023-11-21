#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(sequenza)

seqzdata = sequenza.extract(args[1], assembly = "hg38")

data     = sequenza.fit(seqzdata, female = as.logical(args[4]), chromosome.list = 1:23)
sequenza.results(sequenza.extract = seqzdata, cp.table = data,  sample.id = args[2],  chromosome.list = 1:23, out.dir=args[3], female = as.logical(args[4]))
cint        <- get.ci(data)

cellularity <- cint$max.cellularity

write.table(cellularity, file=paste0(args[2], "_sequenza_purity.txt"), quote=F, sep="\t", row.names=F, col.names=F)

ploidy      <- cint$max.ploidy
ploidy      <- round(ploidy, digits=0)

write.table(ploidy, file=paste0(args[2], "_sequenza_ploidy.txt"), quote=F, sep="\t", row.names=F, col.names=F)
