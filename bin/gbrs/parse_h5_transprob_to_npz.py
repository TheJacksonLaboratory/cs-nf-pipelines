import os
import sys
import pickle
import numpy as np
import h5py
import glob
from collections import defaultdict
from itertools import combinations_with_replacement
from threading import Thread
from optparse import OptionParser

### parsing functions.
def female_file_parse():
    
    chrs_f = list(map(str, np.arange(19) + 1)) + ['X', 'MT']
    h5fh = h5py.File(options.transition_h5)
    sex = 'F'
    gen = str(options.generation)
    tprob = dict()
    tprob_out = dict()
    
    for c in chrs_f:
        keystr = "%s:%s:%s" % (c, gen, sex)
        tprob[c] = h5fh[keystr]  # chr:gen:sex

        # Convert 8x8 transition prob matrix of MT into 36x36 (Female)
        tprob_c = tprob[c]

        if tprob_c.shape[1] == num_genos and tprob_c.shape[2] == num_genos:
            tprob_out[c] = tprob_c
        elif tprob_c.shape[1] == num_haps and tprob_c.shape[2] == num_haps:
            print("\tExpanding tprob matrices for Chr %s" % c)
            tprob_c_big = np.ones((tprob_c.shape[0], num_genos, num_genos)) * np.log(np.nextafter(0,1))
            for gid in range(tprob_c.shape[0]):
                for h1id in range(num_haps):
                    for h2id in range(num_haps):
                        hap1 = haplotypes[h1id]
                        hap2 = haplotypes[h2id]
                        geno1id = geno_id['%s%s' % (hap1, hap1)]
                        geno2id = geno_id['%s%s' % (hap2, hap2)]
                        tprob_c_big[gid][geno1id, geno2id] = tprob_c[gid][h1id, h2id]
            tprob_out[c] = tprob_c_big

    np.savez_compressed('tranprob.DO.G%s.%s.npz' % (gen, sex), **tprob_out)

def male_file_parse():

    chrs_m = list(map(str, np.arange(19) + 1)) + ['X', 'Y', 'MT']
    h5fh = h5py.File(options.transition_h5)
    sex = 'M'
    gen = str(options.generation)
    tprob = dict()
    tprob_out = dict()
    
    print('here we are')

    for c in chrs_m:
        keystr = "%s:%s:%s" % (c, gen, sex)
        tprob[c] = h5fh[keystr]  # chr:gen:sex

        # Convert 8x8 transition prob matrix of X, Y & MT into 36x36 (Male)
        tprob_c = tprob[c]
        if tprob_c.shape[1] == num_genos and tprob_c.shape[2] == num_genos:
            tprob_out[c] = tprob_c
        elif tprob_c.shape[1] == num_haps and tprob_c.shape[2] == num_haps:
            print("\tExpanding tprob matrices for Chr %s" % c)
            tprob_c_big = np.ones((tprob_c.shape[0], num_genos, num_genos)) * np.log(np.nextafter(0,1))
            for gid in range(tprob_c.shape[0]):
                # Expand each gene matrix. 
                # For X, and Y this can take awhile. 
                # For MT is is fairly quick. 
                
                for h1id in range(num_haps):
                    for h2id in range(num_haps):
                        hap1 = haplotypes[h1id]
                        hap2 = haplotypes[h2id]
                        geno1id = geno_id['%s%s' % (hap1, hap1)]
                        geno2id = geno_id['%s%s' % (hap2, hap2)]
                        tprob_c_big[gid][geno1id, geno2id] = tprob_c[gid][h1id, h2id]
            tprob_out[c] = tprob_c_big

    np.savez_compressed('tranprob.DO.G%s.%s.npz' % (gen, sex), **tprob_out)

# # Run file parse in parallel. 
# t1 = Thread(target=female_file_parse, args=[])
# t2 = Thread(target=male_file_parse, args=[])
# t1.start()
# t2.start()

if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option("-t", "--transition_h5", dest="transition_h5",
                    help="transition prob h5 file from R", metavar="TRANSPROB_h5")

    parser.add_option("-l", "--haplotype_list", dest="haplotype_list", default="A, B, C, D, E, F, G, H",
                    help="haplotype list", metavar="HAP_LIST")

    parser.add_option("-s", "--sex", dest="sex", default="F",
                    help="sex used to generate transition probability file", metavar="SEX")

    parser.add_option("-g", "--generation", dest="generation", default="0",
                    help="number of generations to split list", metavar="GEN_NUM")

    (options, args) = parser.parse_args()

    ## parse args, and build array vars. 
    haplotypes = options.haplotype_list.split(",")

    print(haplotypes)

    if not any([x in options.sex for x in ['F', 'M']]):
        sys.exit("'-s' or '--sex' must be either F or M, not: " + options.sex) 

    num_haps = len(haplotypes)

    genotypes = [ h1+h2 for h1, h2 in combinations_with_replacement(haplotypes, 2) ]

    num_genos = len(genotypes)

    geno_id = dict(zip(genotypes, np.arange(num_genos)))

    if 'F' in options.sex:
        female_file_parse()
    if 'M' in options.sex:
        male_file_parse()
    
