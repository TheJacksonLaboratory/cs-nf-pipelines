#!/usr/bin/env python
#	USAGE: python get_candidates.py VCF_FILE OUT_FILE
#   DESCRIPTION:
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Kanika Arora (karora@nygenome.org) and Jennifer M Shelton
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


def filter_call(record):
    '''
        Filter calls to leave calls that may be
        supported by a second Lancet run.
    '''
    if 'called_by' in record.info.keys() and \
            not 'lancet' in record.info['called_by'] and \
            record.info['num_callers'] == 1:
        if 'supported_by' in record.info.keys():
            if record.info['supported_by']:
                return True
        return False
    return True


def filter_vcf(bcf_in, bcf_out):
    '''
        Only print calls that are:
        1) Not already called by Lancet
        2) Only supported by on caller
        3) Not Supported by a support caller
    '''
    for record in bcf_in.fetch():
        if not filter_call(record):
            exit_status = bcf_out.write(record)
            if exit_status != 0:
                print(exit_status)


def main():
    '''
        Only print calls that are:
        1) Not already called by Lancet
        2) Only supported by on caller #MWL NOTE: on caller = one caller? 
        3) Not Supported by a support caller
    '''
    vcf_file = sys.argv[1]
    out_file = sys.argv[2]
    assert os.path.isfile(vcf_file), 'Failed to find somatic VCF call file :' + vcf_file
    bcf_in = read_vcf(vcf_file)
    bcf_out = pysam.VariantFile(out_file, 'w', header=bcf_in.header)
    filter_vcf(bcf_in, bcf_out)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()