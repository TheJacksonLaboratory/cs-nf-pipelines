#!/usr/bin/python3.7

# fasta2GTF.py
#
# by Daniel Leite
# 4.12.2019

# ---------- script info ---------- #

desc="""Script to generate a GTF from a fasta file"""
# ---------- import modules ---------- #

import sys,os
from Bio import SeqIO
import argparse
print(sys.version)
# ---------- setup | user options ---------- #

parser = argparse.ArgumentParser(description=desc)
optional = parser._action_groups.pop()  #Change order of groups
required = parser.add_argument_group('required arguments')
parser._action_groups.append(optional)

#Add the required arguments
required.add_argument('-i', '--input_fasta', dest='fasta_file',
                      help='fasta file', required=True)
                      
required.add_argument('-o', '--output_GTF', dest='GTF_file',
                      help='output GTF file', required=True)

input_args = parser.parse_args()


# ---------- main / set output and input---------- #

cwd = os.getcwd()


GTFfile = input_args.GTF_file

# ---------- generate a list of transcript lengths ---------- #

FastaFile = open(input_args.fasta_file, 'r')
with open(cwd + str("/seq_lengths_tmp.csv"), 'w') as OUT:
    for rec in SeqIO.parse(FastaFile, 'fasta'):
        name = rec.id
        seq = rec.seq
        seqLen = len(rec)
        OUT.write(str(name)+","+str(seqLen)+'\n')
FastaFile.close()

# ---------- generate a GTF using transcript lengths ---------- #

with open(GTFfile, 'w') as GTF:
    with open("seq_lengths_tmp.csv", 'r') as lens:
        for line in lens:
            line = line.strip('\n')
            line = line.split(',')
            print(str(line[0]) + "\tFASTA2GTF\t" + "gene\t" + str("1") + "\t" + str(line[1]) + "\t." + "\t+" + "\t." + "\t"+"""transcript_id  \"""" +line[0]+"""\";  gene_id  \""""+line[0]+"""\";""" , file = GTF)
            print(str(line[0]) + "\tFASTA2GTF\t" + "transcript\t" + str("1") + "\t" + str(line[1]) + "\t." + "\t+" + "\t." + "\t"+"""transcript_id  \"""" +line[0]+"""\";  gene_id  \""""+line[0]+"""\";""" , file = GTF)
            print(str(line[0]) + "\tFASTA2GTF\t" + "exon\t" + str("1") + "\t" + str(line[1]) + "\t." + "\t+" + "\t." + "\t"+"""transcript_id  \"""" +line[0]+"""\";  gene_id  \""""+line[0]+"""\";""" , file = GTF)

os.system("rm seq_lengths_tmp.csv")