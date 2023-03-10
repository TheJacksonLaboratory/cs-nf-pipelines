#!/usr/bin/env Rscript

############# Generate Transition Probabilities from Gene Positions ###################

# This scripts functions to take gene positions, convert to cM positions, then generate a matrix of transition probabilities.

# 1. Get geneIDs and positions from biomaRt ensembl 105. R105 is the version used in the construction of the pseduo-references.
# 2. Convert the gene start positions to cM using mmconvert.
# 3. Jitter overlaping cM positions, to avoid overlap in gene locations. Manipulate dataframe to input for DOQTL function.
# 4. Run do.trans.probs on each chromosome for each gender. Save matrices to h5 file.

################################################################################

# Modified: Michael W. Lloyd
# Date: 03_10_2023

################################################################################
############ loading libraries
suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
  library(dplyr)
  library(mmconvert)
  library(rhdf5)
  library(DOQTL)
})
#############################

option_list = list(
    make_option(c("-e", "--ensembl_build"), type="character", default="105", 
              help="Ensembl build version [default= %default]", metavar="105"),
    make_option(c("-g", "--num_generation"), type="character", default="100",
              help="Number of generations to calculate [default= %default]"),
    make_option(c("-o", "--output_prefix"), type="character", default="tranprob.genes.DO.", 
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

## Obtain geneIDs, and positions rom biomaRt.

ensembl = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")

ensembl_mart <- useEnsembl(biomart = 'genes',
                        dataset = 'mmusculus_gene_ensembl',
                        version = opt$ensembl_build)

attributes <- searchAttributes(mart = ensembl_mart)

genes.with.id=getBM(attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "gene_biotype"), mart = ensembl_mart)


## Convert bp to cM positions using mmconvert.

valid_chr = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 'X')

input_df <- genes.with.id %>%
              dplyr::filter(chromosome_name %in% valid_chr) %>%
              dplyr::rename(marker = ensembl_gene_id, chr = chromosome_name, pos = start_position) %>%
              dplyr::select(chr, pos, marker)
## mmconvert requires the columns: c(chr, pos, marker).

table(genes.with.id$chromosome_name)
## check gene counts by chr.
table(input_df$chr)
## check gene counts by chr.

converted_df <- mmconvert(input_df, input_type = 'bp')
## convert to cM from bp

jittered_converted_df <- converted_df %>%
                          dplyr::mutate(cM = round(ifelse(is.na(cM_coxV3_ave), cM_coxV3_female, cM_coxV3_ave), digits = 4)) %>%
                                          # catch cases where average cM isn't present (e.g., chrX)
                          dplyr::arrange(chr, cM) %>%
                          group_by(chr, cM) %>%
                                          # sort on chr and position (i.e., cM)
                          dplyr::mutate(group_count = n(), row_num = row_number() - 1) %>%
                          dplyr::mutate(jittered_cM = cM + (0.000001 * (row_number() - 1))) %>%
                                          # set a jitter cM value. Where a constant is added to sequential members of the group. Corrects for genes with identical positions.
                          ungroup() %>%
                          dplyr::arrange(chr, jittered_cM) %>%
                          as.data.frame()
## Jitter cM positions that are identical.

### The GRCm38 jitter was originally done in python via the following command:
    # fh = open('gene_start_cM.sorted.tsv') ### NOTE: this file contained geneID, geneStart, cMposition
    # fhout = open('gene_start_cM.sorted.jittered.tsv', 'w')
    # epsilon = 0.000001
    # prechr = '1'
    # preval = 0.0
    # for curline in fh:
    #   item = curline.rstrip().split("\t")
    #   curchr = item[1]
    #   if item[3] != 'given':
    #     curval = float(item[3])
    #   if preval < curval: # next location
    #     preval = curval
    #   elif preval > curval:
    #     if prechr == curchr:
    #     preval += epsilon
    #   item[3] = str(preval)
    #   else: # start of new chromosome
    #     preval = curval
    #   prechr = curchr
    #   else:
    #     preval += epsilon
    #   item[3] = str(preval)
    #   fhout.write("\t".join(item) + "\n")
    # fhout.close()


## Combine autosomal + X with Y and MT.
valid_chr_jittered <- jittered_converted_df %>% dplyr::select(marker, chr, bp_grcm39, jittered_cM) %>% dplyr::rename(pos = bp_grcm39)

chrY_false_cM <- genes.with.id %>% dplyr::filter(chromosome_name == 'Y') %>%
                  dplyr::rename(marker = ensembl_gene_id, chr = chromosome_name, pos = start_position) %>%
                  dplyr::select(marker, chr, pos) %>%
                  dplyr::mutate(jittered_cM = row_number())

chrMT_false_cM <- genes.with.id %>% dplyr::filter(chromosome_name == 'MT') %>%
                  dplyr::rename(marker = ensembl_gene_id, chr = chromosome_name, pos = start_position) %>%
                  dplyr::select(marker, chr, pos) %>%
                  dplyr::mutate(jittered_cM = row_number())

joined_allChr_jittered <- rbind(valid_chr_jittered, chrY_false_cM, chrMT_false_cM )


## Convert to transition probability matrix using DOQTL function.
states = get.do.states()
gens = c(0:100)
num_gens = length(gens)

#
# Female
#
chromosomes = c(1:19, "X")
num_chro = length(chromosomes)
h5file = paste0(opt$output_prefix, 'G0-', opt$num_generation, '.F.h5')
h5createFile(h5file)
sex = "F"

for (i in 1:num_chro) {
  print(i)
  chr = as.character(chromosomes[i])
  genes <- joined_allChr_jittered %>% dplyr::filter(chr == !!chr)
  tprobs = do.trans.probs(states, genes, chr=chr, sex=sex, gen=gens)
  for (j in c(1:num_gens)) {
    print(j)
    h5write(tprobs[[j]], file=h5file, name=paste(chr, j, sex, sep=":"))
  }
}
## Loop over chromosomes, and run do.trans.prob chromosome filtered set of genes. Write out resulting matricies to h5 file.

## NOTE: The function do.trans.probs throws a warning message: "NAs introduced by coercion" when chr == 'X'
##       This error message is due to LINE: https://github.com/dmgatti/DOQTL/blob/a1a4d170bf5923ca45689a83822febdb46ede215/R/transition.probs.R#L361
##       `if(!is.na(as.numeric(chr)))` Will always throw the error when `chr` is a string, or other non-numeric.

#  KB_NOTE: MT (8 states)
genes <- joined_allChr_jittered %>% dplyr::filter(chr == 'MT')
tprobs = do.trans.probs(states, genes, chr="X", sex="M", gen=c(2)) # KB_NOTE: chr="X", sex="M" when we want 8x8 matrices
for (j in c(1:num_gens)) {
  h5write(tprobs[[1]], file=h5file, name=paste("MT", j, sex, sep=":"))
}


#
# KB_NOTE: Male X=8, Y=8 states (Not allowing pseudo-autosomal region recombination)
#
chromosomes = c(1:19)
num_chro = length(chromosomes)
h5file = paste0(opt$output_prefix, 'G0-', opt$num_generation, '.M.h5')
h5createFile(h5file)
sex = "M"
for (i in 1:num_chro) {
  print(i)
  chr = as.character(chromosomes[i])
  genes <- joined_allChr_jittered %>% dplyr::filter(chr == !!chr)
  tprobs = do.trans.probs(states, genes, chr=chr, sex=sex, gen=gens)
  for (j in c(1:num_gens)) {
    print(j)
    h5write(tprobs[[j]], file=h5file, name=paste(chr, j, sex, sep=":"))
  }
}

# KB_NOTE: For "X" chromosome (male only)
genes <- joined_allChr_jittered %>% dplyr::filter(chr == 'X')
tprobs = do.trans.probs(states, genes, chr="X", sex="M", gen=gens)
for (j in c(1:num_gens)) {
  print(j)
  h5write(tprobs[[j]], file=h5file, name=paste("X", j, sex, sep=":"))
}

# KB_NOTE: For "Y" chromosome (male only)
genes <- joined_allChr_jittered %>% dplyr::filter(chr == 'Y')
tprobs = do.trans.probs(states, genes, chr="X", sex="M", gen=gens)
for (j in c(1:num_gens)) {
  print(j)
  h5write(tprobs[[j]], file=h5file, name=paste("Y", j, sex, sep=":"))
}

# KB_NOTE: For "MT" chromosome
genes <- joined_allChr_jittered %>% dplyr::filter(chr == 'MT')
tprobs = do.trans.probs(states, genes, chr="X", sex="M", gen=c(1)) # KB_NOTE: chr="X", sex="M" when we want 8x8 matrices
for (j in c(1:num_gens)) {
  print(j)
  h5write(tprobs[[1]], file=h5file, name=paste("MT", j, sex, sep=":"))
}



# "I used image(tprobs[[1]][,,1]) to make a plot. It has a cool pattern that relates to whether you're changing one letter or two. I could go through it and explain it in a call if you want.  
# Then I did plot(exp(sapply(tprobs, function(z) { z[1,1,1] }))) to get the top corner of each generation. I wanted to see lower probabilities early on and higher ones at later generations. 
# But,  you know, the values are all close to 1. The values in tprobs should be on a log scale, so you have to exponentiate them. 
# Can you look at the old trans probs files and make sure that the values are of a similar sign and order of magnitude? I thought that the exponentiated values would be small positive values, not close to 1."