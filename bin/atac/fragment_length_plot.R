#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

frag_length <- read.table(args[1],  header=F, sep="\t", row.names=NULL, check.names=F, na.strings = '.')

spline_int <- as.data.frame(spline(frag_length$V3, frag_length$V2))

pdf(file='fraglen_plot.pdf')
ggplot(frag_length) + 
  geom_line(data = spline_int, aes(x = x, y = y)) +
  theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90), panel.grid.minor = element_blank()) + 
  xlab("Insert Size (bp)") +
  ylab("Read Count")
dev.off()

temp_df <- t(data.frame('x-axis' = spline_int$x, 'y-axis' = spline_int$y))

rownames(temp_df) <- c(unique(frag_length$V1), unique(frag_length$V1))

write.table(temp_df, quote = F, row.names = T, file = args[2], sep = '\t', col.names = F)
