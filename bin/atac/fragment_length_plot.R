#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

frag_length <- read.table(args[1],  header=F, sep=" ", row.names=NULL, check.names=F, na.strings = '.')

spline_int <- as.data.frame(spline(frag_length$V2, frag_length$V1))

pdf(file='fraglen_plot.pdf')
ggplot(frag_length) + 
  geom_line(data = spline_int, aes(x = x, y = y)) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90), panel.grid.minor = element_blank()) + 
  xlab("Insert Size (bp)") +
  ylab("Read Count")
dev.off()
