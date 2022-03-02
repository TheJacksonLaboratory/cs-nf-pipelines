#!/usr/bin/env python

import csv
import sys
import numpy as np
import collections
import fileinput
from collections import defaultdict
from itertools import imap
from sys import argv
import os
#import cga_version as version

# Check for a version request.
#version.parse_options()





if len(sys.argv) < 2:
    print >>sys.stderr, "This program has no version,It just needs targetbed as input and outputs the target coverage"
    sys.exit()


targetbed = open(sys.argv[1],"r")
targetcov = open(sys.argv[2],"w")






chr_a=[]
start_a=[]
stop_a=[]
bp_a=[]
cov_a=[]
gene_a=[]

with targetbed as f:
     reader=csv.reader(f,delimiter='\t')
     for a,b,c,d,e,f in reader:
         chr_a.append(a)
         start_a.append(b)
         stop_a.append(c)
         gene_a.append(d)
         cov_a.append(f)




target_keys = zip(chr_a,start_a,stop_a,gene_a)

result_dict = defaultdict(list)




print >>targetcov,'chr',"\t", 'start',"\t", 'stop',"\t",'Gene name',"\t",'Mean_coverage',"\t",'Median_coverage',"\t",'min_coverage',"\t",'Max_coverage'


for k,v in zip(target_keys,cov_a):
    result_dict[k].append(v)


for k,v in result_dict.iteritems():
    L = [int(n) for n in v if n]
    print >>targetcov,k[0],"\t",k[1],"\t",k[2],"\t",k[3],"\t",sum(L)/float(len(L)),"\t",np.median(L),"\t",min(L),"\t",max(L)
 


