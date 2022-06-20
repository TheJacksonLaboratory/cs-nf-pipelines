import glob
from itertools import islice

print("----CutAdapt Log----")

for file in glob.glob("*cutadapt.log"):
    with open(file) as f:
        for line in f:
            line = line.rstrip('\n')
                
            if 'Total read pairs processed' in line:
                input_reads = line.split(sep=' ')[-1]
                print('Reads processed:\t' + str(input_reads))

            if 'Pairs written (passing filters)' in line:
                input_reads = line.split(sep=' ')[-2]
                print('Reads passing filter:\t' + str(input_reads))

print("----Bowtie2 Log----")

for file in glob.glob("*bowtie2.log"):
    with open(file) as f:
        
        for line in f:
            line = line.rstrip('\n')

            if 'aligned concordantly 0 times' in line:
                input_reads = line.lstrip(' ').split(sep=' ')[0]
                print('Reads aligned concordantly 0 times:\t' + str(input_reads))

            if 'aligned concordantly exactly 1 time' in line:
                input_reads = line.lstrip(' ').split(sep=' ')[0]
                print('Reads aligned concordantly exactly 1 time:\t' + str(input_reads))

            if 'aligned concordantly >1 times' in line:
                input_reads = line.lstrip(' ').split(sep=' ')[0]
                print('Reads aligned concordantly >1 times:\t' + str(input_reads))

            if 'overall alignment rate' in line:
                input_reads = line.split(sep=' ')[0]
                print('Overall alignment rate:\t' + str(input_reads))


print("----Duplication Metric Log----")

for file in glob.glob("*sorted.metrics"):
    with open(file) as f:
        for line in islice(f, 7, 8):
            input_reads = line.split(sep='\t')[8]
            print('Proportion Duplication:\t' + str(input_reads))

print("----mtDNA Content Log----")

for file in glob.glob("*mtDNA_Content.txt"):
    with open(file) as f:
        for line in f:
            print(line.rstrip('\n'))

print("----NRF and PBC Log----")

for file in glob.glob("*pbc.qc"):
    with open(file) as f:
        for line in f:
            line = line.rstrip('\n')
            input_reads = line.split(sep='\t')
            print("Non-Redundant Fraction (NRF): " + str(input_reads[4]))
            print("PCR Bottlenecking Coefficient 1 (PBC1):\t" + str(input_reads[5]))
            print("PCR Bottlenecking Coefficient 2 (PBC2):\t" + str(input_reads[6]))

print("----Fraction Reads in Peak----")
for file in glob.glob("*Fraction_reads_in_peak.txt"):
    with open(file) as f:
        for line in f:
            line.rstrip('\n')
            input_reads = line.split(sep='\t')
            print('Filtered Read Count:\t' + input_reads[1], end='')
            print('Fraction Reads in Peak:\t' + input_reads[0])