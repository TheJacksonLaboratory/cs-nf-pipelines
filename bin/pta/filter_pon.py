#!/usr/bin/env python
#	USAGE: python filter_pon.py
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
import pandas as pd
import pysam
import argparse
##########################################################################
##############                  Custom functions              ############
##########################################################################


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
        '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def read_bed(bed_file, chrom):
    '''
        Read in BED file. Require TotalUniqueSamples be integers.
    '''
    bed_in = pd.read_csv(bed_file, sep = '\t', dtype={'#CHROM': str})  # auto-detect input format
    assert 'TotalUniqueSamples' in list(bed_in.columns), 'BED file missing TotalUniqueSamples column'
    if chrom:
        bed_in = bed_in[(bed_in['#CHROM'] == chrom)].copy()
    bed_in['TotalUniqueSamples'] = bed_in['TotalUniqueSamples'].astype('int64')
    return bed_in


def compare_num(observed, rule, value):
    '''
        Compare two numbers. Available rules are lt, gt and eq
    '''
    observed = int(observed)
    if rule == 'lt':
        if observed < value:
            return False
    elif rule == 'gt':
        if observed > value:
            return False
    elif rule == 'eq':
        if observed == value:
            return False
    return True


def compare_str(observed, rule, value):
    '''
        Compare two strings to see that they are the same 'eq'
        or are not the same 'ne'
    '''
    if rule == 'eq':
        if value == observed:
            return False
    elif rule == 'ne':
        if value != observed:
            return False
    return True


def test_rules(filter_func, rule):
    '''
        Test that rule matches combination.
    '''
    combos = {compare_num : ['lt', 'gt', 'eq'],
              compare_str : ['eq', 'ne']}
    assert rule in combos[filter_func], 'rule not is possible rules for given function. Possible rules: ' + str(combos[filter_func]) + ' Rule : ' + rule


def custom(row, filter_func, key, value, rule):
    '''
        Run custom filter for bed file row
    '''
    observed = row[key]
    return filter_func(observed, rule, value)

                         
def filter_bed(bed_in,
               filter_func, key,
               value, rule):
    '''
        Filter based on column and rule.
    '''
    bed_in['fail'] = bed_in.apply(lambda row: custom(row, filter_func, key, value, rule), axis=1)
    bed_in_filtered = bed_in[(bed_in.fail == False)].copy()
    return bed_in_filtered


def get_bad_pos(bed_in_filtered):
    '''
        Grab filtered position info
    '''
    bed_in_filtered['key'] = bed_in_filtered.apply(lambda row: row['#CHROM'] + '{' + str(row.END), axis=1)
    bad_pos = set(bed_in_filtered.key)
    return bad_pos


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


def compose(bcf_in, bad_pos):
    '''
       Filter based on the PON bad positions.
    '''
    passing = False
    for record in bcf_in.fetch():
        key = record.chrom + '{' + str(record.pos)
        if key in bad_pos:
            filters = record.filter.keys()
            if len(filters) == 1 and filters[0] in ['PASS', 'SUPPORT']:
                record.filter.clear()
            record.filter.add('PON')
        yield record


def write_file(bcf_out, record):
    '''
        Write to a VCF.
    '''
    exit_status = bcf_out.write(record)
    if exit_status != 0:
        print(exit_status)
    return bcf_out


def main():
    '''
        Filter VCF for start positions which match
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Filters a VCF file \
                                     if a start position of a VCF matches the \
                                     start + 1 position of the bed file. Command-line \
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
    parser.add_argument('-b', '--bed',
                        dest='bed_file',
                        help='Input BED file to use to filter',
                        required=True)
    parser.add_argument('-d', '--default',
                        dest='default',
                        help='Default filter for TotalUniqueSamples. Result \
                        must be greater than this value [1]',
                        default='1')
    parser.add_argument('-f', '--filter',
                        dest='filter_func',
                        choices=['str',
                                 'num'],
                        nargs='+',
                        help='Filter type(s)',
                        default=[False])
    parser.add_argument('-r', '--rule',
                        dest='rule',
                        nargs='+',
                        help='Rule(s) to filter based on value. Acceptable options \
                        are "lt", "gt", "eq" for "num" and \
                        "eq" and "ne" for "str"',
                        default=[False])
    parser.add_argument('-val', '--value',
                        dest='value',
                        nargs='+',
                        help='Value(s) used to compare to custom filter',
                        default=[False])
    parser.add_argument('-c', '--chrom',
                        dest='chrom',
                        help='Chrom used for filtering',
                        default=False)
    args = parser.parse_args()
    assert os.path.isfile(args.vcf_file), 'Failed to find caller VCF call file :' + args.vcf_file
    assert os.path.isfile(args.bed_file), 'Failed to find BED file ' + args.bed_file
    if args.filter_func[0]:
        assert args.key[0], 'key is required for custom filter'
        assert args.value[0], 'value is required for custom filter'
        assert args.rule[0], 'rule is required for custom filter'
    #  ==========================
    #  Filter PON
    #  ==========================
    functions = {'num' : compare_num,
                'str' : compare_str}
    # default filter
    bed_file = args.bed_file
    filter_func = functions['num']
    key = 'TotalUniqueSamples'
    rule = 'gt'
    value = int(args.default)

    bed_in = read_bed(bed_file, args.chrom)
    test_rules(filter_func, rule)
    bed_in_filtered = filter_bed(bed_in, filter_func, key, value, rule)
    # optional_filters
    if args.filter_func[0]:
        for i in range(len(args.filter_func)):
            filter_func = functions[args.filter_func[i]]
            key = args.key[i]
            rule = args.rule[i]
            value = args.value[i]
            test_rules(filter_func, rule)
            bed_in_filtered = filter_bed(bed_in_filtered, filter_func, key, value, rule)
    bad_pos = get_bad_pos(bed_in_filtered)
    #  ==========================
    #  Filter with PON
    #  ==========================
    bcf_in = read_vcf(args.vcf_file)
    bcf_in = add_filter_header(bcf_out=bcf_in,
                              id='PON',
                              description='Variant in panel of normal database')
    bcf_out = pysam.VariantFile(args.out, 'w', header=bcf_in.header)
    for record in compose(bcf_in, bad_pos):
       bcf_out = write_file(bcf_out, record)




##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()