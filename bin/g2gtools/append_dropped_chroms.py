import os
import sys
import glob
from collections import defaultdict
from optparse import OptionParser
import gzip


parser = OptionParser()

parser.add_option("-u", "--unmapped", dest="unmapped",
                  help="unmapped file", metavar="PATH/TO/UNMAPPED")

parser.add_option("-v", "--vci", dest="vci",
                  help="vci file to parse", metavar="PATH/TO/VCI")

parser.add_option("-g", "--gtf", dest="gtf",
                  help="gtf file to append to", metavar="PATH/TO/GTF")

parser.add_option("-o", "--output", dest="output",
                  help="output file", metavar="OUTPUT_FILE")


(options, args) = parser.parse_args()

if not os.path.exists(options.unmapped):
    raise SystemExit("\nERROR: Can not find the unmapped file: " + options.unmapped + '.  Check it exists. Exiting.\n')

if not os.path.exists(options.vci):
    raise SystemExit("\nERROR: Can not find the vci file: " + options.vci + '. Check it exists. Exiting.\n')

if not os.path.exists(options.gtf):
    raise SystemExit("\nERROR: Can not find the vci file: " + options.vci + '. Check it exists. Exiting.\n')

vci_chroms = set()
gtf_chroms = set()

with gzip.open(options.vci, 'rt') as vci_file:        
    for line in vci_file:
        if '##CONTIG' in line:
            vci_chroms.add(line.rstrip().split(":")[0].replace("##CONTIG=", ""))
        if '#' in line:
            pass
        else:
            gtf_chroms.add(line.rstrip().split("\t")[0])
vci_file.close()

dropped_chroms = list(vci_chroms.difference(gtf_chroms))

with open(options.gtf, "r") as gtf_file, open(options.output,'w') as output_file:
    # read content from first file
    for line in gtf_file:             
             # append content to second file
             output_file.write(line)
    gtf_file.close()
    with open(options.unmapped) as unmapped_file:
        for line in unmapped_file:
            if line.rstrip().split("\t")[0] in dropped_chroms:
                # print(line)
                output_file.write(line)
unmapped_file.close()
output_file.close()
