################################################################################
# Given a counts file, use Ddx3y and Xist to determine the sex of each sample
# in the counts file. Generate a sample metadata file containing the sample ID
# and the sex of each mouse.
#
# Arguments:
# arg 1: string containing the full path to the counts file for one sample.
#        This should be a tab-delimited file output from the CS-NF RNASeq
#        pipeline. The first column should be named "gene_id" and 
#        should contain the Ensembl ID, gene symbol, or a concatenated
#        <Ensembl ID>_<gene symbol> name. There must also be a column
#        called "TPM".
# arg 2: string containing the full path to the directory to write files.
# arg 3: string containing the file name prefix for output file.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2024-09-23
################################################################################

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {

  stop('ERROR: Three arguments required in determine.sex.R.\nusage: Rscript determine_sex.R <counts_file> <output_path> <file_prefix>')

} # if(length(args) != 3)

# Path to the counts file.
counts_file = args[1]

# Path to the output directory.
output_path = args[2]

# Prefix for output file.
file_prefix = args[3]

# Mouse and human Ensembl IDs & symbols.
xist_search_str  = 'ENSMUSG00000086503|ENSG00000229807|xist'
ddx3y_search_str = 'ENSMUSG00000069045|ENSG00000067048|ddx3y'


##### FUNCTIONS #####

# Given an Ensembl gene ID and a gene symbol, search for the row which
# contains one or the other in column 1 of the counts. 
# Arguments:
# search_str: string containing the Ensembl ID and/or gene symbol of of gene 
#             to search for. e.g. "ENSMUSG00000086503|ENSG00000229807|xist"
# counts:  data.frame containing the gene counts. First column MUST
#          contain either the Ensembl ID or gene symbols 
# Returns: integer that is the row containing the requested gene.
find_gene_row = function(search_str, counts) {

  # Search for the gene.
  gene_row = grep(search_str, counts$gene_id, ignore.case = TRUE)

  # If we didn't find it, quit with an error.
  # TBD: Return -1 and let the caller handle errors?
  if(length(gene_row) == 0) {
  
   stop(paste('determine_sex.R: Could not find', search_str, 
        'in', counts, '.'))
  
  } # if(length(gene_row) == 0)

  return(gene_row)

} # find_gene_row()


##### MAIN #####

# Verify that the file exists.
if(!file.exists(counts_file)) {

  stop(paste('ERROR: File', counts_file, 'not found.'))

} # if(!file.exists(counts_file))

# Read in the counts.
counts = read.delim(counts_file)

# Verify that we have a "gene_id" column.
if(!'gene_id' %in% colnames(counts)) {

  stop(paste('determine_sex.R: column name gene_id not found in', counts_file))

} # if(!'gene_id' %in% colnames(counts))

# Verify that we have a "TPM" column.
if(!'TPM' %in% colnames(counts)) {

  stop(paste('determine_sex.R: column name TPM not found in', counts_file))

} # if(!'TPM' %in% colnames(counts))

# Find the row conatining Xist counts.
xist_row = find_gene_row(search_str = xist_search_str,
                         counts     = counts)

# Find the row conatining Ddx3y counts.
ddx3y_row = find_gene_row(search_str = ddx3y_search_str,
                          counts     = counts)

# Get the counts for the two genes.
sex_counts        = counts$TPM[c(xist_row, ddx3y_row)]
names(sex_counts) = c('Xist', 'Ddx3y')

sex = data.frame(id    = file_prefix,
                 xist  = sex_counts['Xist'],
                 ddx3y = sex_counts['Ddx3y'],
                 sex   = ifelse(sex_counts['Xist'] == 0 & sex_counts['Ddx3y'] == 0, 
                                NA, 
                                ifelse(sex_counts['Xist'] > sex_counts['Ddx3y'], 
                                'female', 'male')
                                )
) 

# Write out a file containing the estimated sex for each sample.
write.csv(sex, file = file.path(output_path, paste0(file_prefix, '_SexDet.csv')), 
          quote = FALSE, row.names = FALSE)
