#!/usr/bin/env python
# USAGE: reorder_vcf.py VCF_IN VCF_OUT NORMAL TUMOR
# DESCRIPTION: Print a VCF file with the sample order indicated in the 3rd
# and 4th arguments


# ## MWL NOTE: 
# This script requires the header and input 'tumor/normal' names in the 3rd and 4th arg to match. 
#       If you pass names NOT present in the header, it will simply emit the file AS IS.
#       The script DOES NOT inform the user of if a change has been made in the sample order. 
#       NOTE ALSO: if the header already contains the strings 'TUMOR' and 'NORMAL,
#                  'TUMOR and NORMAL are RENAMED to string provided in 3rd and 4th args. 

import sys
import pandas as pd
import shutil

def load_header(vcf_in):
    '''
        Load a VCF file header as a list of lines.
    '''
    with open(vcf_in) as vcf:
        header = [line for line in vcf if line.startswith('#')]
    return header


def load_vcf(vcf_in, header, reorder, paired, normal, tumor):
    '''
        Load a VCF file as an pandas dataframe.
    '''
    names = header[-1].rstrip().replace('^#', '').split('\t')
    if paired and reorder:
        names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                 'INFO', 'FORMAT', tumor, normal]
    elif paired and not reorder:
        names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                 'INFO', 'FORMAT', normal, tumor]
    vcf_reader = pd.read_csv(vcf_in, comment='#',
                             names=names, sep='\t',
                             dtype={'#CHROM' : str})
    return vcf_reader


def reorder_vcf(vcf_reader, reorder):
    '''
        reorder corrected names.
    '''
    if reorder:
        vcf_reader = vcf_reader[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                                 'INFO', 'FORMAT', normal, tumor]]
    return vcf_reader


def order_wrong(last_header, normal, tumor):
    '''
        Check if the order is good
    '''
    header_parts = last_header.rstrip().split('\t')
    if header_parts[9] in [tumor, 'TUMOR'] and \
            header_parts[10] in [normal, 'NORMAL']:
        return True
    elif tumor in header_parts[9] and \
            normal in header_parts[10]:
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


def reorder_column(vcf_file, vcf_out_file, normal, tumor):
    '''
        Order columns Normal and then Tumor
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
    vcf_reader = load_vcf(vcf_file, header, reorder, paired, normal, tumor)
    if reorder:
        vcf_reader = reorder_vcf(vcf_reader, reorder)
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
       Write out the VCF file with corrected information
    '''
    with open(vcf_out_file, 'w') as vcf_out:
        for line in header[:-1]:
            vcf_out.write(line)
    vcf_reader.to_csv(vcf_out_file, sep='\t',
                      mode='a', index=False)


#  =====================
#  Main
#  =====================
if __name__ == "__main__":
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    normal = sys.argv[3]
    tumor = sys.argv[4]
    reorder_column(vcf_file, vcf_out_file, normal, tumor)

