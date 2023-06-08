from Bio import SeqIO
from optparse import OptionParser


parser = OptionParser()

parser.add_option("-i", "--input_file", dest="input_file",
                  help="gene list TSV from gene_bp_to_cM R script", metavar="gene_list_ensemblBuild_105.tsv")

parser.add_option("-o", "--output_name", dest="output_name", default="ref.gene_pos.ordered.npz",
                  help="a directory where final output will be stored", metavar="PREFIX")

parser.add_option("-l", "--poly_A_length", dest="poly_A_length", default="99",
                  help="a directory where final output will be stored", metavar="99")


(options, args) = parser.parse_args()


with open(options.input_file) as handle:
    fasta_recs = []
    for record in SeqIO.parse(handle, "fasta"):
        record.seq = record.seq + "A" * int(options.poly_A_length)
        fasta_recs.append(record)
    
SeqIO.write(fasta_recs, options.output_name, "fasta")