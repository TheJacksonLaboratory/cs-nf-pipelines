#!/usr/bin/env python
# USAGE: rename_vcf.py VCF_IN VCF_OUT NORMAL TUMOR PREFIX
# DESCRIPTION: Print a VCF file with the sample order indicated in the 3rd
# and 4th arguments and a prefix added to the sample names

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


def load_vcf(vcf_in, header, reorder, paired, normal, tumor, prefix):
    '''
        Load a VCF file as an pandas dataframe.
    '''
    names = header[-1].rstrip().replace('^#', '').split('\t')
    if paired:
        names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                 'INFO', 'FORMAT',
                 prefix + '_' + normal, prefix + '_' + tumor]
    vcf_reader = pd.read_csv(vcf_in, comment='#',
                             names=names, sep='\t')
    return vcf_reader


def order_wrong(last_header, normal, tumor):
    '''
        Check if the order is perfect
    '''
    header_parts = last_header.rstrip().split('\t')
    if header_parts[9] == tumor and \
            header_parts[10] == normal:
        return True
    else:
        return False


def check_paired(last_header):
    '''
        Return False if VCF has only a single sample
    '''
    header_parts = last_header.rstrip().split('\t')
    if len(header_parts) == 10:
        return False
    return True


def rename(vcf_file, vcf_out_file, normal, tumor, prefix):
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
    #  =====================
    #  reorder
    #  =====================
    header = load_header(vcf_file)
    last_header = header[-1]
    paired = check_paired(last_header)
    if not paired:
        reorder = False
    else:
        reorder = order_wrong(last_header, normal, tumor)
    if reorder:
        log.error('VCF must start with expected sample names in the order normal, tumor.')
        sys.exit(1)
    vcf_reader = load_vcf(vcf_file, header, reorder,
                          paired, normal, tumor, prefix)
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
    normal = sys.argv[3]
    tumor = sys.argv[4]
    prefix = sys.argv[5]
    assert os.path.isfile(vcf_file), 'Failed to find prep caller VCF call file :' + vcf_file
    rename(vcf_file, vcf_out_file, normal, tumor, prefix)


#  =====================
#  Main
#  =====================


if __name__ == "__main__":
    main()
