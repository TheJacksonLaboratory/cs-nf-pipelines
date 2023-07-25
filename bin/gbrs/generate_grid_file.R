#!/usr/bin/env Rscript

############# Generate Transition Probabilities from Gene Positions ################

# This scripts functions to space 'marker' positions every 0.02 cM along chr1:19 and X.
# cM positions are converted to bp and Mbp via mmconvert. 
# 'marker' positions are set to start at cM = 0 = 3,000,000bp,
#  and capped to not exceed the distal telomere starting positions.
# cM cap is originally set by the max chrom length, pulled from mmconvert. 

###################################################################################

# Author: Michael W. Lloyd
# Date: 03_20_2023

################################################################################
############ loading libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(mmconvert)
})
#############################

telomere_position = c(195054279, 181655017, 159645316, 
                      156760686, 151658149, 149488044, 
                      144895196, 130027694, 124259700, 
                      130430862, 121873369, 119992757, 
                      120783175, 125039656, 103973951,
                      97908968,  95194699,  90620763, 
                      61320004, 169376592)
# telomere positions were extracted from UCSC table browser: 
#      https://genome.ucsc.edu/cgi-bin/hgTables --> 'All Tracks' --> 'Gap' 

chrom_list = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 'X')

telomeres <- data.frame(chr = chrom_list, telomere_position)

max_lengths = mmconvert(mmconvert::grcm39_chrlen) %>% 
                    mutate(cM = ifelse(!is.na(cM_coxV3_ave), cM_coxV3_ave, cM_coxV3_female))

converted_data <- list()
for(c in chrom_list) {
  print(c)

  max_cm = max_lengths %>% dplyr::filter(chr == c) %>% dplyr::select(cM)
  
  input_df <- data.frame(chr=c,
                         pos=seq(from=0.00001,to=max_cm[1,],by=0.02),
                         marker=paste(c,seq(from=0.00001,to=max_cm[1,],by=0.02), sep='_'))
  
  if (c != 'X') {
    converted_df <- mmconvert(input_df, input_type = 'ave_cM')
  } else {
    converted_df <- mmconvert(input_df, input_type = 'female_cM')
  }
  
  telomere_cap = telomeres %>% dplyr::filter(chr == c)
  
  converted_data[[c]] <- converted_df %>% dplyr::filter(bp_grcm39 <= telomere_cap$telomere_position)

}

combined_list <- do.call(rbind.data.frame, converted_data)

options(scipen = 999) ## turn off sci notation, to avoid printout of e.g., 3e+6

combined_list <- combined_list %>% 
                  mutate(marker = paste(chr, round(bp_grcm39), sep = '_'),
                         bp_grcm39 = round(bp_grcm39),
                         Mbp_grcm39 = bp_grcm39 / 1000000,
                         cM = ifelse(!is.na(cM_coxV3_ave), cM_coxV3_ave, cM_coxV3_female)) %>% 
                  dplyr::select(marker, chr, pos = bp_grcm39, cM, bp = bp_grcm39)

write.table(combined_list, file = 'ref.genome_grid.GRCm39.tsv', sep="\t", row.names = FALSE, quote=FALSE)
