#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 5) {
  sample_name <- args[1]
  nanosv_depth_fh <- args[2]
  sniffles_depth_fh <- args[3]
  merged_depth_fh <- args[4]
  joined_results_fh <- args[5]
} else {
  stop("Incorrect number of arguments. Please supply 1) sample name, 2) nanosv depths, 3) sniffles depths, 4) survivor IDs, 5) the survivor summary")
}


output_fh <- paste(sample_name, "survivor_joined_results_depths.csv", sep = ".")
output_bed_fh <- paste(sample_name, "merged_depths.bed", sep = ".")

nanosv_depth <- read_csv(nanosv_depth_fh)
sniffles_depth <- read_csv(sniffles_depth_fh)
merged_depth <- read_csv(merged_depth_fh, na = "NaN")
merged_depth <- merged_depth %>% mutate(sniffles_id = as.character(sniffles_id))
# This catches edge cases where no sniffles calls are made. 
# NanoSV is not expected to have this edge case, as the tool fails when coverage is too low. 
joined_results <- read_csv(joined_results_fh, col_types = cols("chr" = "c"))

merged_1 <- filter(merged_depth, !is.na(nanosv_id)) %>%
  left_join(., nanosv_depth, by = c("nanosv_id" = "SV_ID")) %>%
  select(., SURVIVOR_ID, DR, DV) %>%
  rename(., NanoSV_DR = DR, NanoSV_DV = DV)

merged_2 <- filter(merged_depth, !is.na(sniffles_id)) %>%
  left_join(., sniffles_depth, by = c("sniffles_id" = "SV_ID")) %>%
  select(., SURVIVOR_ID, DR, DV) %>%
  rename(., sniffles_DR = DR, sniffles_DV = DV)

joined_results %>%
  left_join(., merged_1, by = c("sv_name" = "SURVIVOR_ID")) %>%
  left_join(., merged_2, by = c("sv_name" = "SURVIVOR_ID")) %>%
  write_csv(., output_fh)

joined_results %>%
  left_join(., merged_1, by = c("sv_name" = "SURVIVOR_ID")) %>%
  left_join(., merged_2, by = c("sv_name" = "SURVIVOR_ID")) %>%
  select(., chr, pos, sv_name, NanoSV_DR, NanoSV_DV, sniffles_DR, sniffles_DV) %>%
  write_delim(., output_bed_fh, delim = "\t", col_names = F)
