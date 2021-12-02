#!/usr/bin/env Rscript
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

surv_summary <- read_csv(args[1])

anx_frame <- read_delim(args[2],delim = "\t")

surv_frame <- left_join(read_delim(args[1], delim = "\t"),
                        read_csv(args[2]),
                        by = c("chr" = "chr",
                               "pos" = "pos",
                               "SV" = "sv_name")) %>%
              mutate(., start = pos - 1,
                        end = start + abs(sv_size),
                        sv_string = str_c(sv_type, SV, sep = ":"),
                        chr_bare = str_replace(chr, "chr", ""))

surv_frame %>%
  filter(., sv_type == "INS") %>%
  select(., chr_bare, start, end, SV) %>%
  write_delim(., str_c(args[3], ".ins.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "DEL") %>%
  select(., chr_bare, start, end, SV) %>%
  write_delim(., str_c(args[3], ".del.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "INV") %>%
  select(., chr_bare, start, end, SV) %>%
  write_delim(., str_c(args[3], ".inv.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "DUP") %>%
  select(., chr_bare, start, end, SV) %>%
  write_delim(., str_c(args[3], ".dup.bed", sep = ""), col_names = F, delim = "\t")

surv_frame %>%
  filter(., sv_type == "TRA") %>%
  select(., chr_bare, start, end, SV) %>%
  write_delim(., str_c(args[3], ".tra.bed", sep = ""), col_names = F, delim = "\t")