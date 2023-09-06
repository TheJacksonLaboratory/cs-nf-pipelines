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
# Version: 0.2
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


def remove_info(bcf_in, csq_columns):
    '''
        Remove a INFO field from VCF.
    '''
    for id in bcf_in.header.info.keys():
        if not id in ['HighConfidence','TYPE', 'called_by', 'num_callers',
                      'supported_by', 'CSQ', 'CancerGeneCensus'] + csq_columns:
            bcf_in.header.info.remove_header(id)
    return bcf_in


def remove_format(bcf_in):
    '''
        Remove a FORMAT field from VCF.
    '''
    for id in bcf_in.header.formats.keys():
        if not id in ['AD','DP', 'AF']:
            bcf_in.header.formats.remove_header(id)
    return bcf_in



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
        assert len(csq_columns) == len(csq_values), 'failed because lengths do not match'
        csq_dict = dict(zip(csq_columns, csq_values))
        csq_dicts[i] = csq_dict
    return csq_dicts


def modify_record(record, csq_columns, good_fields):
    '''
        Shorten CSQ fields
    '''
    csq_dicts = get_csqs(record, csq_columns)
    csq_out = '|'.join([csq_dicts[0][key] for key in good_fields])
    return csq_out


def check_build(bcf_in, genome_version):
        '''
            Check if genome is in a list of supprted non-human genomes.
        '''
        VEP_line = [metadata.value for metadata in bcf_in.header.records if metadata.key == 'VEP'][0]
        vep_info  = {entry.split('=')[0] : entry.split('=')[-1]  for entry in VEP_line.split(' ')}
        if genome_version in vep_info['assembly']:
            return False
        else:
            return True


class Variant(object):
    
    
    def __init__(self, record, csq_out, human=True):
        self.record = record
        self.csq_out = csq_out
        self.human = human
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
        self.good_format = ['AD','DP', 'AF']
        if self.human:
            self.good_info = ['HighConfidence','TYPE', 'called_by', 'num_callers',
                              'supported_by', 'CSQ', 'CancerGeneCensus']
        else:
            self.good_info = ['HighConfidence','TYPE', 'called_by', 'num_callers',
                              'supported_by', 'CSQ']
        self.info_dict = self.get_info()
        self.samples[0] = self.fix_format(self.samples[0].split(':'))
        self.samples[1] = self.fix_format(self.samples[1].split(':'))
        self.format = [key for key in self.good_format if key in self.format]
    
    def get_info(self):
        '''
            Get current info line and add prefix if needed
            '''
        info_dict = OrderedDict()
        for item in self.info:
            if item.split('=')[0] in self.good_info:
                if item.split('=')[0] == 'CSQ':
                    info_dict.update({'CSQ' : self.csq_out})
                elif '=' in item:
                    info_dict.update({item.split('=')[0] : item.split('=')[1]})
                else:
                    info_dict.update({item : None })
        return info_dict
    
    def fix_format(self, sample):
        '''
            Reduce to good formats
        '''
        format_dict = dict(zip(self.format, sample))
        new_format = [format_dict[key] for key in self.good_format if key in format_dict]
        return ':'.join(new_format)

    def write(self):
        if ':'.join(self.format) == '':
            self.format = '.'
            self.samples = ['.', '.']
        if ';'.join(self.filters) == 'PASS':
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
        else:
            return False


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def write_vcf(bcf_in, vcf_out_file, csq_columns, human=True):
    '''
        Write out the download
    '''
    # CSQ to keep
    if human:
        good_fields = ['Gene', 'BIOTYPE', 'CLIN_SIG', 'Consequence',
                       'CosmicCoding', 'CosmicCoding_AA', 'CosmicNonCoding',
                       'Existing_variation', 'GnomadExomes_AF', 'GnomadGenomes_AF',
                       'HGVSc', 'HGVSp', 'IMPACT', 'Polyphen2_HVAR_pred',
                       'FATHMM_pred', 'fathmm-MKL_coding_pred',
                       'SIFT4G_pred', 'SIFT_pred', 'SYMBOL', 'SYMBOL_SOURCE', 'AF_1000G']
    else:
        good_fields = ['Gene', 'BIOTYPE', 'Consequence',
                       'Existing_variation', 'HGVSc', 'HGVSp', 'IMPACT',
                       'NEAREST', 'SIFT', 'SYMBOL', 'SYMBOL_SOURCE']

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
            csq_out = modify_record(record, csq_columns, good_fields)
            line = Variant(record, csq_out, human).write()
            if line:
                vcf_out.write(line + '\n')


def main():
    '''
        Reduce metadata in VCF for main VCF output
    '''
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    genome_version = sys.argv[3]
    bcf_in = read_vcf(vcf_file)
    human = check_build(bcf_in, genome_version)
    csq_columns = bcf_in.header.info['CSQ'].description.split()[-1].split('|') # grab the definitions
    bcf_in = remove_format(bcf_in)
    bcf_in = remove_info(bcf_in, csq_columns)
    write_vcf(bcf_in, vcf_out_file, csq_columns, human=human)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
