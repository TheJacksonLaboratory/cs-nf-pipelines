#!/usr/bin/env python
#	USAGE: python cancer_gene_census.py cancer_gene_census.csv VCF VCF_OUT
#   DESCRIPTION: Annotates files by adding information about the
# Cosmic Genome Census entry for the nearest gene.
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
import pprint
from collections import OrderedDict
##########################################################################
##############                  Custom functions              ############
##########################################################################


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


def af_1000g(bcf_out):
    '''
        Add description for VEP annotation.
        '''
    for af in ['AF_1000G', 'AFR_AF_1000G', 'AMR_AF_1000G',
               'EAS_AF_1000G', 'EUR_AF_1000G', 'SAS_AF_1000G']:
        description = af + ' field from phase3 1000genomes'
        bcf_out = add_info_header(bcf_out,
                                  id=af,
                                  number='.',
                                  type='String',
                                  description=description)
    return bcf_out


def get_good_fields(csq_columns, suffix='_1000G'):
    '''
        rename CSQ fields as needed.
    '''
    change = ['AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF']
    good = []
    for key in csq_columns:
        if key in change:
            key += suffix
        good.append(key)
    return good

class Variant(object):
    
    
    def __init__(self, record):
        self.record = record
        self.line = str(self.record).rstrip()
        self.parts = self.line.split('\t')
        # VCF columns
        self.chrom = self.parts[0]
        self.pos = self.parts[1]
        self.id = self.parts[2]
        self.ref = self.parts[3]
        self.alts = self.parts[4].split(',')
        self.qual = self.parts[5]
        self.filters = self.parts[6].split(';')
        self.info = self.parts[7].split(';')
        self.format = self.parts[8].split(':')
        self.samples = self.parts[9:]
        # modify
        self.info_dict = self.get_info()
    
    
    def get_info(self):
        '''
            Get current info line and add prefix if needed
            '''
        info_dict = OrderedDict()
        for item in self.info:
            if '=' in item:
                if '=' in item:
                    info_dict.update({item.split('=')[0] : '='.join(item.split('=')[1:])})
            else:
                info_dict.update({item : None })
        return info_dict
    
    
    
    def write(self):
        line = [self.chrom,
                self.pos,
                self.id,
                self.ref,
                ','.join(self.alts),
                str(self.qual),
                ';'.join(self.filters),
                ';'.join(['='.join([x for x in [key, self.info_dict[key]] if x != None]) for key in self.info_dict]),
                ':'.join(self.format)]
        line += self.samples
        self.new_line = '\t'.join(line)
        return self.new_line


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def write_vcf(bcf_in, vcf_out_file, csq_columns):
    '''
        Write out the download
    '''
    # CSQ to keep
    good_fields = get_good_fields(csq_columns)
    # Import the header after removal of extra metadata
    header = str(bcf_in.header).rstrip()
    csq_format = '|'.join(good_fields)
    # Write new header with corrected CSQ and fewer metadata keys overall
    with open(vcf_out_file, 'w') as vcf_out:
        for line in header.split('\n'):
            if 'ID=CSQ' in line:
                line = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: ' + csq_format + '">'
            vcf_out.write(line + '\n')
        for record in bcf_in:
            line = Variant(record).write()
            if line:
                vcf_out.write(line + '\n')


def main():
    '''
        Reduce metadata in VCF for main VCF output
    '''
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    bcf_in = read_vcf(vcf_file)
    bcf_in = af_1000g(bcf_out=bcf_in)
    csq_columns = bcf_in.header.info['CSQ'].description.split()[-1].split('|') # grab the definitions
    write_vcf(bcf_in, vcf_out_file, csq_columns)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
