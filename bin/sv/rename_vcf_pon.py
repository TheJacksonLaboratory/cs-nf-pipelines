#!/usr/bin/env python
# USAGE: rename_vcf.py VCF_IN VCF_OUT TUMOR PREFIX
# DESCRIPTION: Print a VCF file with a prefix added to the sample names

import sys
import pandas as pd
import shutil
import logging as log
import os


def load_header(vcf_in):
    '''
        Load a VCF file header as a list of lines.
    '''
    with open(vcf_in) as vcf:
        header = [line for line in vcf if line.startswith('#')]
    return header


def load_vcf(vcf_in, header, tumor, prefix):
    '''
        Load a VCF file as an pandas dataframe.
    '''
    names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                 'INFO', 'FORMAT', prefix + '_' + tumor]
    vcf_reader = pd.read_csv(vcf_in, comment='#',
                             names=names, sep='\t')
    return vcf_reader


def rename(vcf_file, vcf_out_file, tumor, prefix):
    '''
        Add prefix to sample name
    '''
    #  =====================
    #  test if renaming should occur
    #  =====================
    rename = False
    if vcf_out_file == vcf_file:
        vcf_out_file = vcf_file + '_tmp.vcf'
        rename = True
    header = load_header(vcf_file)
    vcf_reader = load_vcf(vcf_file, header, tumor, prefix)
    #  =====================
    #  rewrite
    #  =====================
    vcf_writer(header, vcf_reader, vcf_out_file)
    #  =====================
    #  rename output VCF
    #  =====================
    if rename:
        shutil.move(vcf_out_file, vcf_file)
    return True


def vcf_writer(header, vcf_reader, vcf_out_file):
    '''
       Write out the VCF file with corrected sample names
    '''
    with open(vcf_out_file, 'w') as vcf_out:
        for line in header[:-1]:
            vcf_out.write(line)
    vcf_reader.to_csv(vcf_out_file, sep='\t',
                      mode='a', index=False)


def main():
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    tumor = sys.argv[3]
    prefix = sys.argv[4]
    assert os.path.isfile(vcf_file), 'Failed to find prep caller VCF call file :' + vcf_file
    rename(vcf_file, vcf_out_file, tumor, prefix)


#  =====================
#  Main
#  =====================


if __name__ == "__main__":
    main()
