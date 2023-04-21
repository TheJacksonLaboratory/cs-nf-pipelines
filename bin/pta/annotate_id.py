#!/usr/bin/env python
#	USAGE: python annotate_id.py VCF VCF_OUT
#   DESCRIPTION: Annotates files by adding information about the
# CosmicID to the ID field.
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Jennifer M Shelton
##################### /COPYRIGHT ###############################################
################################################################################


import sys
import os
import logging as log
import pysam
##########################################################################
##############                  Custom functions              ############
##########################################################################


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def get_csq_columns(bcf_in):
    '''
        get column names from the bar
        separated CSQ VEP annotation
        results. CSQ are Consequence
        annotations from Ensembl VEP.
        '''
    csq_columns = bcf_in.header.info['CSQ'].description.split()[-1].split('|') # grab the definitions
    return csq_columns


def get_csqs(record, csq_columns):
    '''
        Get new INFO field results.
    '''
    alt_count = len(record.alts)
    csq_dicts = {}
    for i in range(alt_count):
        try:
            csq_line = record.info['CSQ'][i]
        except UnicodeDecodeError: # for names with accents and other unexpected characters (rare)
            line = str(record)
            csq_line = line.split('\t')[7].split('CSQ=')[1]
            csq_line = csq_line.split(';')[0]
        csq_line = csq_line.split(',')[i]
        csq_values = csq_line.split('|')
        csq_dict = dict(zip(csq_columns, csq_values))
        csq_dicts[i] = csq_dict
    return csq_dicts


def get_ID(record, csq_columns):
    '''
        Get new ID results is a Cosmic Coding or non-coding result is available.
    '''
    csq_dicts = get_csqs(record, csq_columns)
    alt_count = len(record.alts)
    coding_ids = []
    noncoding_ids = []
    for i in range(alt_count):
        if 'CosmicCoding' in csq_dicts[i]:
            coding_ids += [csq_dicts[i]['CosmicCoding'].replace('&', ';')]
        if 'CosmicNonCoding' in csq_dicts[i]:
            noncoding_ids += [csq_dicts[i]['CosmicNonCoding'].replace('&', ';')]
        ids = ';'.join([id for id in coding_ids + noncoding_ids if not id == ''])
    return ids


def fix_gt(gt):
    '''
        change GT 0/0/0/1/0, 0/0/0/1, 0/0/1, etc to 0/1
    '''
    if len(gt.split('/')) > 2:
        gt = '0/1'
    return gt


def modify_record(record, csq_columns):
    '''
        Add new ID field to records as needed
    '''
    # 'CosmicCoding', 'CosmicNonCoding'
    gts = [key for key in record.samples[0].keys() if key.endswith('_GT')]
    for key in gts:
        gt = record.samples[0][key]
        record.samples[0][key] = fix_gt(gt)
        gt = record.samples[1][key]
        record.samples[1][key] = fix_gt(gt)
    ids = get_ID(record, csq_columns)

    if ids != '':
        if record.id != '':
            record.id = record.id + ';' + ids
        else:
            record.id = ids
    return record



def write_vcf(bcf_in, vcf_out_file, csq_columns):
    '''
        Write out the VCF
    '''
    bcf_out = pysam.VariantFile(vcf_out_file, 'w', header=bcf_in.header)
    for record in bcf_in.fetch():
        record = modify_record(record, csq_columns)
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)


def main():
    '''
        Annotates files by adding information about the
        Cosmic coding and non-coding entries to the ID column.
        Also changes GT 0/0/0/1/0, 0/0/0/1, 0/0/1, etc to 0/1.
    '''
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    assert os.path.isfile(vcf_file), 'Failed to find caller VCF call file :' + vcf_file
    bcf_in = read_vcf(vcf_file)
    csq_columns = get_csq_columns(bcf_in)
    write_vcf(bcf_in, vcf_out_file, csq_columns)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
