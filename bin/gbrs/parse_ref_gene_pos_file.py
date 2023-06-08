import os
import numpy as np
from collections import defaultdict
from optparse import OptionParser


parser = OptionParser()

parser.add_option("-i", "--input_file", dest="input_file",
                  help="gene list TSV from gene_bp_to_cM R script", metavar="gene_list_ensemblBuild_105.tsv")

parser.add_option("-t", "--gene_to_transcript", dest="gene_to_transcript",
                  help="gene list TSV from gene_bp_to_cM R script", metavar="emase.gene2transcripts.tsv")

parser.add_option("-o", "--output_directory", dest="output_directory", default="./",
                  help="a directory where final output will be stored", metavar="PATH/TO/OUTPUT")

parser.add_option("-n", "--output_name", dest="output_name", default="ref.gene_pos.ordered.npz",
                  help="a directory where final output will be stored", metavar="PREFIX")

(options, args) = parser.parse_args()

if not os.path.exists(options.input_file):
    raise SystemExit("\nERROR: Can not find the input directory: " + options.input_file + '.  Check it exists. Exiting.\n')

if not os.path.exists(options.gene_to_transcript):
    raise SystemExit("\nERROR: Can not find the gene_to_transcript file: " + options.gene_to_transcript + '.  Check it exists. Exiting.\n')

if not os.path.exists(options.output_directory):
    raise SystemExit("\nERROR: Can not find the output directory: " + options.input_directory + '.  Check it exists. Exiting.\n')


known_genes = []
with open(options.gene_to_transcript) as fh:
    for curline in fh:
        item = curline.rstrip().split("\t")
        known_genes.append(item[0])

gene_list = defaultdict(list)
counter = 0
total = 0

with open(options.input_file) as fh:
    fh.readline()
    for curline in fh:
        item = curline.rstrip().split("\t")
        if item[0] in known_genes:
            gene_list[item[1]].append([item[0],item[3]])
            counter += 1
        total += 1
gene_list = dict(gene_list)

print('Processed %s total genes from the input file, of which %s were found in the emase gene2transcript list' % (total, counter))

np.savez_compressed(options.output_directory + options.output_name, **gene_list)

# It is possible that genes are filtered out during the conversion from GTF to FASTA in g2gtools and prepare emase. 
# The final `ref.gene_pos.ordered.npz` must only contain genes that are present in the built reference set.

# `known_genes` provides the set of genes "known" to the reference set. 
# "gene_list_ensemblBuild_105.tsv" comes from pulling the full ENSEMBL gene list from biomaRt. 
# The logic gate: if item[0] in known_genes ensures that ref.gene_pos.ordered.npz contains only genes found in the g2gtools/prepare-emase built reference set. 
