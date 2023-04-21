#!/usr/bin/env python
#	USAGE: python merge_prep.py
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
import pandas as pd
import re
import argparse
##########################################################################
##############                  Custom functions              ############
##########################################################################
base_pattern = re.compile(r'^[ACGTN]+$')



def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
        '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def add_info_header(bcf_out,
                    id,
                    number,
                    type,
                    description):
    '''
        Add new INFO field
    '''
    bcf_out.header.info.add(id=id,
                            number=number,
                            type=type,
                            description=description)
    return bcf_out


def add_filter_header(bcf_out,
                      id,
                      description):
    '''
        Add new FILTER field
    '''
    bcf_out.header.filters.add(id=id,
                               number=None,
                               type=None,
                               description=description)
    return bcf_out


def remove_filters_header(bcf_out):
    '''
        Remove all filters except PASS
    '''
    for id in bcf_out.header.filters.keys():
        if not id in ['PASS','SUPPORT']:
            bcf_out.header.filters.remove_header(id)
    return bcf_out


def pass_alleles(record, base_pattern=base_pattern):
    '''
        Pass lines that have no special characters in REF/ALT
    '''
    passed = True
    alleles = list(record.alts) + [record.ref]
    for allele in alleles:
        if not re.match(base_pattern, allele):
            passed = False
    return passed


def prep_record(record, tool, passing, support):
    '''
        Pass lines that have no special characters in ALT/REF
    '''
    record.id = None
    record.qual = None
    if passing:
        if support:
            if tool == 'manta':
                tool_supported_by = tool + 'SV'
            else:
                tool_supported_by = tool
            record.info['supported_by'] = (tool_supported_by,)
        else:
            record.info['called_by'] = (tool,)
            record.info['num_callers'] = 1
    return record


def write_file(bcf_in, bcf_out, tool, filter=True, support=False):
    '''
       Filter based on FILTER column,
       also filter lines with special characters in ALT/REF
    '''
    passing = False
    for record in bcf_in.fetch():
        filters = record.filter.keys()
        #  ====================
        #  Passing variants
        #  ====================
        if len(filters) == 1 and \
                filters[0] == 'PASS':
            passing = True
            write = True
            if support:
                record.filter.clear()
                record.filter.add('SUPPORT')
        #  ====================
        #  Failing variants
        #  ====================
        else:
            if filter:
                write = False
            else:
                write = True
        if write and not pass_alleles(record):
            write = False
        if write:
            record = prep_record(record, tool, passing, support)
            exit_status = bcf_out.write(record)
            if exit_status != 0:
                print(exit_status)
    return True


def main():
    '''
        Prepare the VCF file for merging by:
        1) 'TYPE' is added to header
        2) 'called_by' is added to header
        3) 'num_callers' is added to header
        4) filter lines with special characters in REF/ALT (e.g. <DEL> in manta)
        5) fill in 'called_by', 'num_callers'
        6) 'SUPPORT' FILTER line is added
        7) non 'PASS'/'SUPPORT' FILTER lines are removed (if not skip-filter)
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and preps the file by: \
                                     1) "TYPE" is added to header \
                                     2) "called_by" is added to header \
                                     3) "num_callers" is added to header \
                                     4) filter lines with special characters in REF/ALT (e.g. <DEL> in manta) \
                                     5) fill in "called_by", "num_callers" \
                                     6) "SUPPORT" FILTER line is added \
                                     7) non "PASS"/"SUPPORT" FILTER lines are removed (if not skip-filter) \
                                     . Command-line \
                                     options that may be omitted (i.e. are NOT \
                                     required) are shown in square brackets.')
                                     #    Documentation parameters
    parser.add_argument('-v', '--vcf',
                        dest='vcf_file',
                        help='VCF file',
                        required=True)
    parser.add_argument('-o', '--out',
                        dest='out',
                        help='Output VCF file',
                        required=True)
    parser.add_argument('-t', '--tool',
                        dest='tool',
                        choices=['strelka2_sv',
                                 'strelka2_indel',
                                 'mutect2',
                                 'svaba',
                                 'lancet',
                                 'manta'],
                        help='Tool name',
                        required=True)
    parser.add_argument('-s', '--support',
                        dest='support',
                        help='Use if calls are only support calls',
                        action='store_true')
    parser.add_argument('-f', '--skip-filter',
                        dest='skip_filter',
                        help='Remove calls that are not PASS or SUPPORT',
                        action='store_true')
    args = parser.parse_args()
    filter = True
    if args.skip_filter:
        filter = False
    assert os.path.isfile(args.vcf_file), 'Failed to find caller VCF call file :' + args.vcf_file
    #  ==========================
    #  Run prep
    #  ==========================
    bcf_in = read_vcf(args.vcf_file)
    bcf_in = add_info_header(bcf_out=bcf_in,
                              id='called_by',
                              number='.',
                              type='String',
                              description='Name of the variant caller(s) that the variant was called by')
    bcf_in = add_info_header(bcf_out=bcf_in,
                              id='num_callers',
                              number='1',
                              type='Integer',
                              description='Number of callers')
    bcf_in = add_info_header(bcf_out=bcf_in,
                              id='supported_by',
                              number='.',
                              type='String',
                              description='Name of the tool(s) apart from the main variant callers in the pipeline that support the variant')
    bcf_in = add_filter_header(bcf_out=bcf_in,
                              id='SUPPORT',
                              description='Variant from Validation caller')
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    if filter:
        bcf_out = remove_filters_header(bcf_out)
    write_file(bcf_in,
               bcf_out,
               tool=args.tool,
               filter=filter,
               support=args.support)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()