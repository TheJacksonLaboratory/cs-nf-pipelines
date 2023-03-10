import os
try:
  import cPickle as pickle
except:
  import pickle
import numpy as np
import h5py
import glob
from collections import defaultdict
from itertools import combinations_with_replacement
from threading import Thread
from optparse import OptionParser


parser = OptionParser()

parser.add_option("-f", "--female_h5", dest="female_h5",
                  help="transition prob h5 file from R for females", metavar="FEMALE_h5")

parser.add_option("-m", "--male_h5", dest="male_h5",
                  help="transition prob h5 file from R for males", metavar="MALE_h5")

parser.add_option("-l", "--haplotype_list", dest="haplotype_list", default="A, B, C, D, E, F, G, H",
                  help="haplotype list", metavar="HAP_LIST")

parser.add_option("-g", "--generations", dest="generations", default="100",
                  help="number of generations to split list", metavar="GEN_NUM")

(options, args) = parser.parse_args()

## parse args, and build array vars. 
haplotypes = options.haplotype_list.split(",")
num_haps = len(haplotypes)

genotypes = [ h1+h2 for h1, h2 in combinations_with_replacement(haplotypes, 2) ]

num_genos = len(genotypes)

geno_id = dict(zip(genotypes, np.arange(num_genos)))

chrs_f = list(map(str, np.arange(19) + 1)) + ['X', 'MT']
chrs_m = list(map(str, np.arange(19) + 1)) + ['X', 'Y', 'MT']


### parsing functions.
def female_file_parse():
    h5fh = h5py.File(options.female_h5)
    sex = 'F'
    for gen in range(1, int(options.generations) + 1):
        gen = str(gen)
        tprob = dict()
        for c in chrs_f:
            keystr = "%s:%s:%s" % (c, gen, sex)
            tprob[c] = h5fh[keystr]  # chr:gen:sex
        np.savez_compressed('tranprob.DO.G%s.%s.npz' % (gen, sex), **tprob)

def male_file_parse():
    h5fh = h5py.File(options.male_h5)
    sex = 'M'
    for gen in range(1, int(options.generations) + 1):
        gen = str(gen)
        tprob = dict()
        for c in chrs_f:
            keystr = "%s:%s:%s" % (c, gen, sex)
            tprob[c] = h5fh[keystr]  # chr:gen:sex
        np.savez_compressed('tranprob.DO.G%s.%s.npz' % (gen, sex), **tprob)

# Run file parse in parallel. 
t1 = Thread(target=female_file_parse, args=[])
t2 = Thread(target=male_file_parse, args=[])
t1.start()
t2.start()