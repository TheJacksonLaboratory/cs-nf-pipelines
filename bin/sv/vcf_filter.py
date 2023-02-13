#!/usr/bin/env python
#	USAGE: python vcf_filter.py GRM_FILE VCF_FILE OUT_FILE
#   DESCRIPTION:
################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.2
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


def test_af(record, alt_index, af=0.01):
    '''
        Test if AF > 0.01 in one germline database. PASS
        variants that don't have AF listed (for example
        in records from mouse databases)
        Pass variants that don't have sample columns (e.g. mouse variants from 00-All.normalized.vcf.gz)
    '''
    samples = record.samples.keys()
    none_count = 0
    if len(samples) > 0:
        for sample_name in samples:
            if not 'AF' in record.format.keys():
                return True
            elif record.samples[sample_name]['AF'][alt_index]:
                if record.samples[sample_name]['AF'][alt_index] > af:
                    return True
            else:
                none_count += 1
        if none_count < len(samples):
            return False
    return True


def is_germline(germ_in, record, alt, af=0.01):
    '''
        Check if matching variant is in the GRM VCF.
    '''
    for germ_record in germ_in.fetch(record.contig, record.pos - 1, record.pos):
        if germ_record.ref == record.ref:
            for alt_index, germ_alt in enumerate(germ_record.alts):
                if test_af(germ_record, alt_index, af=af):
                    if germ_alt == alt:
                        return True
    return False


def filter_vcf(bcf_in, bcf_out, germ_in):
    '''
        Change filter column from PASS to GRM or
        add GRM to filters.
    '''
    for record in bcf_in.fetch():
        filters = record.filter.keys()
        filter = False
        for alt in record.alts:
            if is_germline(germ_in, record, alt):
                filter = True
        if filter:
            if len(filters) == 1 and filters[0] == 'PASS':
                record.filter.clear()
            record.filter.add('GRM')
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)


def main():
    '''
        Change filter column from PASS to GRM or
        add GRM to filters.
    '''
    germ_file = sys.argv[1]
    vcf_file = sys.argv[2]
    out_file = sys.argv[3]
    assert os.path.isfile(germ_file), 'Failed to find germline VCF call file :' + germ_file
    assert os.path.isfile(vcf_file), 'Failed to find somatic VCF call file :' + vcf_file
    germ_in = read_vcf(germ_file)
    bcf_in = read_vcf(vcf_file)
    bcf_in = add_filter_header(bcf_in,
                               id='GRM',
                               description='Known germline variant')
    bcf_out = pysam.VariantFile(out_file, 'w', header=bcf_in.header)
    filter_vcf(bcf_in, bcf_out, germ_in)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()