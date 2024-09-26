import os
import sys
import glob
import numpy as np
from collections import defaultdict
from optparse import OptionParser


parser = OptionParser()

parser.add_option("-i", "--input_directory", dest="input_directory",
                  help="a directory that contains output directories from the 'run_emase' workflow", metavar="PATH/TO/INPUT")

parser.add_option("-o", "--output_directory", dest="output_directory", default="./",
                  help="a directory where final output will be stored", metavar="PATH/TO/OUTPUT")

parser.add_option("-p", "--output_prefix", dest="output_prefix", default="gbrs_emissions",
                  help="the prefix for the output", metavar="PREFIX")

parser.add_option("-g", "--gene2transcripts_file", dest="gene2transcripts",
                  help="the emase gene2transcripts file from the 'prepare_emase' workflow", metavar="emase.gene2transcripts.tsv")

parser.add_option("-m", "--metadata_file", dest="metadata_file",
                  help="comma delimited metadata file. File must have headers 'do_id' and 'sampleID' where: do_id is ('A', 'B', ..., 'G', 'H'), and sampleID are the IDs produced and used by 'run_emase' (i.e., output directory names) ", metavar="metadata.csv")

parser.add_option("-s", "--strains", dest="strains", default="A,B,C,D,E,F,G,H",
                  help="comma delimited list of strains", metavar="A,B,C,D,E,F,G,H")

parser.add_option("-e", "--min_expression", dest="min_expression", default="2",
                  help="minimum expression count to include a gene in Avec computation.", metavar="EXP_MIN")

(options, args) = parser.parse_args()

if not os.path.exists(options.input_directory):
    raise SystemExit("\nERROR: Can not find the input directory: " + options.input_directory + '.  Check it exists. Exiting.\n')

if not os.path.exists(options.output_directory):
    raise SystemExit("\nERROR: Can not find the output directory: " + options.input_directory + '.  Check it exists. Exiting.\n')

if not os.path.exists(options.gene2transcripts):
    raise SystemExit("\nERROR: Can not find the gene2transcripts file: " + options.gene2transcripts + '. Check it exists. Exiting.\n')

if not os.path.exists(options.metadata_file):
    raise SystemExit("\nERROR: Can not find the metadata file: " + options.metadata_file + '. Check it exists. Exiting.\n')


def unit_vector(vector):
    if sum(vector) > 1e-6:
        return vector / np.linalg.norm(vector)
    else:
        return vector

print('\nAll inputs exist.\nStarting gene2transcript parse...')

### Import gene2transcript to get the number of total genes expected in dataset. 
gname = np.loadtxt(options.gene2transcripts, usecols=(0,), dtype='str')
num_genes = len(gname)
gid = dict(zip(gname, np.arange(num_genes)))
print('Number of genes in gene2transcript file: ' + str(num_genes))


### Import metadata file to parse into a list: 
###     A: ['sample1', 'sample2'], B: ['sample3', 'sample4']...
print('Importing metadata from file...')
sample_id = defaultdict(list)
with open(options.metadata_file) as fh:
    header = fh.readline().rstrip().split(",")

    try:
        res = [header.index(i) for i in ['hap_id', 'sampleID']]
    except Exception as e:
        print('Error processing metadata file:', e)


    for curline in fh:
        item = curline.rstrip().split(",")
        sample_id[item[res[0]]].append(item[res[1]])

sample_id = dict(sample_id)

print('Found the following keys in "hap_id": ' + str(sample_id.keys()))

print('The total number of samples by strain in the metadata file is: ')

for st, slist in sample_id.items():
    print(st, len(slist))

print('Found the following sample counts matched to sample IDs in the input directory...')

strains = [i.strip() for i in options.strains.split(',')]

num_strains = len(strains)

dlist = defaultdict(list)
for st in strains:
    for d in sample_id[st]:
        d_found = glob.glob('%s/%s*/' % (options.input_directory, d))
        if d_found:
            dlist[st].append(*d_found)
    print(len(dlist[st]), st)
dlist = dict(dlist)

print('Parsing the found directories, and importing `multiway.genes.tpm` files...')

dset = dict()
for st in strains:
    dmat_strain = np.zeros((num_genes, num_strains))
    for d in dlist[st]:
        dmat_sample = np.zeros((num_genes, num_strains))
        tpmfile = glob.glob(d + "emase/*.multiway.genes.tpm")
        if not os.path.isfile(tpmfile[0]):
            print("File %s does not exist." % tpmfile)
            continue        
        with open(tpmfile[0]) as fh:
            fh.readline()
            for curline in fh:
                item = curline.rstrip().split("\t")
                row = gid[item[0]]
                dmat_sample[row, :] = list(map(float, item[1:num_strains+1]))
        dmat_strain += dmat_sample
    dset[st] = dmat_strain / len(dlist[st])

print('Confirming we have TPM for all found genes and strains...')


for st in strains:
    print(round(dset[st].sum()), dset[st].shape, st)

print('Generating Avecs from strain count data...')

min_expr = int(options.min_expression)
axes = dict()
ases = dict()
avecs = dict()
all_genes = set()
passed_genes = set()

for g in gname:
    all_genes.add(g)
    axes[g] = np.zeros((num_strains, num_strains))
    ases[g] = np.zeros((1, num_strains))
    good = np.zeros(num_strains)
    for i, st in enumerate(strains):
        v = dset[st][gid[g], :]
        axes[g][i, :] = v
        ases[g][0, i] = sum(v)
        if sum(v) > min_expr:
            good[i] = 1.0
    if sum(good) > 0:  # At least one strain expresses
        avecs[g] = np.zeros((num_strains, num_strains))
        passed_genes.add(g)
        for i in range(num_strains):
            avecs[g][i, :] = unit_vector(axes[g][i, :])

with open(os.path.join(options.output_directory, options.output_prefix + '.included_genes.txt'), 'w') as f:
    for item in list(all_genes.intersection(passed_genes)):
        f.write("%s\n" % item)

with open(os.path.join(options.output_directory, options.output_prefix + '.excluded_genes.txt'), 'w') as f:
    for item in list(all_genes.difference(passed_genes)):
        f.write("%s\n" % item)

print('After filtering, there are ' + str(len(avecs)) + ' genes.\nExporting vectors to output directory...')

# np.savez_compressed(os.path.join(options.output_directory, options.output_prefix + '.axes.npz'), **axes)
# np.savez_compressed(os.path.join(options.output_directory, options.output_prefix + '.ases.npz'), **ases)
np.savez_compressed(os.path.join(options.output_directory, options.output_prefix + '.avecs.npz'), **avecs)
