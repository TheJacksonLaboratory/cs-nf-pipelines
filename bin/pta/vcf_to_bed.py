#!/usr/bin/env python
#	USAGE: vcf_to_bed.py <vcf>
#   DESCRIPTION: Make bed file from VCF.
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Version: 1.0
# Author: Jennifer M Shelton, Andre Corvelo
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


def feed_vcf(bcf_in):
    '''
        Generate relevant columns for BED.
        Converts to 0-based intervals
    '''
    for record in bcf_in.fetch():
        yield record.chrom, record.pos - 1, record.pos + len(record.alts[0]) - 1


def make_bed(vcf_file):
    '''
        Make BED file from the interval list
    '''
    log.info('#######################################')
    log.info('# Making bed...')
    log.info('#######################################')
    bcf_in = read_vcf(vcf_file)
    for chrom, start, end in feed_vcf(bcf_in):
        sys.stdout.write('\t'.join([chrom, str(start), str(end), '.']) + '\n')
    log.info('#######################################')
    log.info('# Done making bed.')
    log.info('#######################################')


def main():
    '''
        Makes BED file from a VCF.
    '''
    assert os.path.isfile(sys.argv[1]); 'Failed to open reference VCF file'
    vcf_file = sys.argv[1]
    make_bed(vcf_file)

##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()