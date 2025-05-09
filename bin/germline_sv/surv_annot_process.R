#!/usr/bin/env Rscript

library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 1) {
  sample_name <- args[1]
} else {
  stop("Incorrect number of arguments. Please supply only the sample name as an argument.")
}

sv_olap_fh <- paste(sample_name, "merged.overlap.annotated.txt", sep = ".")
sv_summary_fh <- paste(sample_name, "survivor_summary.csv", sep = ".")

surv_frame <- left_join(read_delim(sv_olap_fh, delim = "\t",
                                  col_types = cols("chr" = "c")),
                        read_csv(sv_summary_fh,
                                  col_types = cols("chr" = "c")),
                        by = c("chr" = "chr",
                               "pos" = "pos",
                               "SV" = "sv_name"))

if(nrow(surv_frame) == 0) {
  surv_frame <- surv_frame %>% mutate(., start = numeric(), end = numeric(), sv_string = character())
} else {
  surv_frame <- surv_frame %>% mutate(., chr = str_replace(chr, "chrM", "MT"),
          chr = str_replace(chr, "chr", ""),
          start = pos - 1,
          end = start + abs(sv_size),
          sv_string = str_c(sv_type, SV, sep = ":"))
}
## Note: this conditional catches edge cases where no SV are called

surv_frame %>%
  filter(., sv_type == "INS") %>%
  select(., chr, start, end, SV) %>%
  write_delim(., str_c(sample_name, ".ins.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "DEL") %>%
  select(., chr, start, end, SV) %>%
  write_delim(., str_c(sample_name, ".del.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "INV") %>%
  select(., chr, start, end, SV) %>%
  write_delim(., str_c(sample_name, ".inv.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "DUP") %>%
  select(., chr, start, end, SV) %>%
  write_delim(., str_c(sample_name, ".dup.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "TRA") %>%
  select(., chr, start, end, SV) %>%
  write_delim(., str_c(sample_name, ".tra.bed", sep = ""), col_names = F, delim = "\t")
