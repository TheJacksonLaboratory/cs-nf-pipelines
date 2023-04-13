#!/usr/bin/env Rscript

############# Gathering counts from GBRS output for founders ###################

# This script goes to a directory containing the "gbrs_quantified_multiway_genes_
# expected_read_counts" for each sample and convert them into
# a matrix (samples x total gene counts per gene).

################################################################################

# Author: Isabela Gerdes Gyuricza - Churchill lab
# Date: 05_04_2020
# Modified: Michael W. Lloyd
# Date: 03_07_2023

################################################################################
############ loading libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(grid)
  library(gridExtra)
})

CCcolors <- list('AJ' = "#FFDC00", 'B6' = "#888888", '129' = "#F08080", 'NOD' = "#0064C9", 'NZO' = "#7FDBFF", 'CAST' = "#2ECC40", 'PWK' = "#FF4136", 'WSB' = "#B10DC9")

CCfills <- unname(CCcolors)

option_list = list(
    make_option(c("-i", "--input_directory"), type="character", default=NULL, 
              help="directory containing output from EMASE/GBRS quantification", metavar="path/to/EMASE_RUN_OUTPUT"),
    make_option(c("-p", "--file_pattern"), type="character", default=".multiway.genes.expected_read_counts",
              help="search pattern for files containing expression counts [default= %default]"),
    make_option(c("-m", "--metadata_csv"), type="character", default=NULL,
              help="metadata csv file. Must have headers 'sampleID', 'do_id'. 'sampleID' should correspond to the sampleIDs used in EMASE processing.", metavar="path/to/CSV"),
    make_option(c("-o", "--output_prefix"), type="character", default="collectedCounts_output_file", 
              help="output prefix [default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_directory)){
  print_help(opt_parser)
  stop("Input directory must be specified.n", call.=FALSE)
}

if (is.null(opt$metadata_csv)){
  cat("\nNote: No CSV metadata called. Can not generate genotype strain expression proportion PDF summary plot.\n\n")
}

################################################################################

# # Create list of all files in specific folder with specified file pattern. 

file_list <- list.files(opt$input_directory, pattern=opt$file_pattern, recursive = TRUE, full.names = TRUE)

cat(paste0("---Importing expression counts from ", length(file_list), ' found files.\n'))

# # The `list_counts`` function reads in each .txt file in file_list and create a data frame with the same name as the file
list_counts <- function (x) {
  x <- read.table(x, header=TRUE, sep = "\t")
  colnames(x) <- c("Target_ID", LETTERS[1:8], "Total")
  return(x) }

# # Apply the function to the list
df_counts_list <- sapply(file_list, list_counts, simplify = FALSE, USE.NAMES = TRUE)

# # Adjusting the names of the dataframes from full paths to just "sampleID"
names(df_counts_list) <- gsub(opt$file_pattern, "", basename(file_list))

rm(list_counts, file_list)

cat("---Done Importing.\n")

# #Checking if all the genes ids (Target_ID column) for all the samples are the same

temp <- df_counts_list
for (i in 1:length(df_counts_list)){
  temp[[i]] <- select(temp[[i]],Target_ID)
}; rm(i)

b <- suppressMessages(bind_cols(temp))

b <- apply(b, 2, function (x) factor(x,labels = c(1:nrow(b))))

b <- apply(b, 2, as.numeric)

compare_genes <- function(x) {
  cat("---Ensuring all samples share the same gene IDs.\n")
  c <- apply(x, 1, sd)

  if (sum(c) > 0) {
    cat(paste0("---Some gene IDs do not match for the row index: ",which(sum(c) > 0)))
  } else {
    cat("---All gene IDs match.\n")
  }

}

compare_genes(b)

rm(compare_genes, temp, b)

# Creating input data with the total gene counts for all the samples

cat("---Converting counts to single dataframe, and exporting RDS and CSV files.\n")

df_counts <- data.frame(row.names=df_counts_list[[1]]$Target_ID,
                        stringsAsFactors = FALSE) # Genes from sample[1] are used as they are confirmed identical across samples.

df_counts <- df_counts_list
for (i in 1:length(df_counts_list)){
  df_counts[[i]] <- select(df_counts_list[[i]],Target_ID,Total)
}; rm(i)

df_counts <- bind_rows(df_counts,.id = "Sample")

df_counts <- df_counts %>%
  spread(Target_ID,Total) %>%
  column_to_rownames("Sample")

# # Saving the final expression matrix

saveRDS(df_counts, file = paste0(opt$output_prefix,".RDS"))

export_df <- as.data.frame(t(df_counts)) %>% rownames_to_column('geneID')

write.table(export_df, file = paste0(opt$output_prefix,".csv"), sep=",", row.names = FALSE, quote=FALSE)

# # # Genotype by proportion plots. 

if (!is.null(opt$metadata_csv)) {

  cat("---Generating genotype by expression proportion plots.\n")

  metadata <- read.table(opt$metadata_csv, sep=',', header=TRUE)

  # # The `genotype_proportion` function Checking for genotypes..
  genotype_proportion <- function(x) {
    prop <- x %>%
      rowwise() %>%
      mutate(A = A/Total,
             B = B/Total,
             C = C/Total,
             D = D/Total,
             E = E/Total,
             F = F/Total,
             G = G/Total,
             H = H/Total) %>%
      ungroup()
    return(prop)
  }

  genotype_check <- sapply(df_counts_list, genotype_proportion, simplify = FALSE, USE.NAMES = TRUE)

  genotype_check <- sapply(genotype_check, function(x) apply(x[,2:9], 2, mean, na.rm = TRUE),
                           simplify = TRUE, USE.NAMES = TRUE)

  genotype_check <- genotype_check %>%
    as.data.frame() %>%
    rownames_to_column("GBRS_Genotype") %>%
    gather("sampleID", "Mean_Proportion_Expression", -GBRS_Genotype) %>%
    left_join(metadata, by = 'sampleID') %>%
    mutate(GBRS_Genotype = factor(GBRS_Genotype, levels = LETTERS[1:8], labels = c("AJ", "B6", "129", "NOD", "NZO", "Cast", "PWK", "WSB"))) %>%
    mutate(Sample_Strain = factor(strain, levels = c("A_J", "B6", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ")))

  print(nrow(genotype_check %>% dplyr::filter(is.na(Sample_Strain))))

  print(genotype_check %>% dplyr::filter(is.na(Sample_Strain)))

  write.table(genotype_check, file = paste0(opt$output_prefix,"_genotype_prop_check.csv"), sep=",", row.names = FALSE, quote=FALSE)

  pdf( paste0(opt$output_prefix,"_founders_expression_proportions.pdf"), width = 12, height = 8)

  p <- genotype_check %>%
        ggplot(aes(x = GBRS_Genotype, y = Mean_Proportion_Expression, fill = GBRS_Genotype)) +
        geom_boxplot() +
        xlab("EMASE Quantification Genotype") +
        ylab("Mean Proportion Genotype Expression") +
        scale_fill_manual(values = CCfills) +
        facet_wrap(~ Sample_Strain, ncol = 4, scales = "free") +
        theme_minimal() +
        theme(legend.position = "none")

  grid.arrange(p, top='Mouse Sample Strain')
  dev.off()

}
