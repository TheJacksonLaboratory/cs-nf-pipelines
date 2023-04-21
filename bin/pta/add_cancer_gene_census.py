#!/usr/bin/env python
#	USAGE: python cancer_gene_census.py cancer_gene_census.csv VCF VCF_OUT
#   DESCRIPTION: Annotates files by adding information about the
# Cosmic Genome Census entry for the gene symbol.
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
import pandas as pd
import pysam
import numpy as np
##########################################################################
##############                  Custom functions              ############
##########################################################################


def get_record(gene, key, cancer_gene_census):
    '''
        Get value for field from Cancer Gene Census records.
        Ignores gene that are associated with multiple loci (e.g T-cell 
        receptor (TR*) and immunoglobin (IG*) genes)
    '''
    if gene in ['C15orf65', 'CTNNA2', 'HMGN2P46', 'IGH', 'IGK',
                'IGL', 'KAT6A', 'TRA', 'TRB', 'TRD']:
         value = ''
    else:
        try:
            value = cancer_gene_census[(cancer_gene_census['Gene Symbol'] == gene)][key].values[0]
        except (KeyError, IndexError):
            value = ''
    return value


def convert_numpy(x):
    '''
        Convert Numpy nan to blank
    '''
    if isinstance(x, float):
        if np.isnan(x):
            return ''
    return str(x)


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
    spanning_deletion_offset = 0
    for i in range(alt_count):
        j = i - spanning_deletion_offset
        if record.alts[j] == '*':
            csq_dicts[i] = {csq_column : '' for csq_column in csq_columns}
            spanning_deletion_offset += 1
        else:
            try:
                csq_line = record.info['CSQ'][j]
            except UnicodeDecodeError: # for names with accents and other unexpected characters (rare)
                line = str(record)
                csq_line = line.split('\t')[7].split('CSQ=')[1]
                csq_line = csq_line.split(';')[0]
                csq_line = csq_line.split(',')[j]
            csq_values = csq_line.split('|')
            csq_dict = dict(zip(csq_columns, csq_values))
            csq_dicts[i] = csq_dict
    return csq_dicts


def get_CancerGeneCensus(record, csq_columns, cancer_gene_census):
    '''
        Get new INFO field results.
    '''
#    Example results:
#    ('Mutation Types', 'A')
#    ('Mutation Types', 'A, O, Mis')
#    ('Mutation Types', 'F, Mis')
#    ('Mutation Types', 'F; Mis')
#    ('Mutation Types', 'Promoter Mis')
#    grab the gene symbol from the first annotation ALT
    CancerGeneCensus_lists = []
    csq_dicts = get_csqs(record, csq_columns)
    for i, alt in enumerate(record.alts):
        CancerGeneCensus_list = []
        for key in ['Tier', 'Hallmark', 'Somatic', 'Germline',
                    'Tumour Types(Somatic)', 'Tumour Types(Germline)',
                    'Cancer Syndrome', 'Tissue Type', 'Molecular Genetics',
                    'Role in Cancer', 'Mutation Types']:
            result = get_record(csq_dicts[i]['SYMBOL'], key, cancer_gene_census)
            if key == 'Hallmark':
                if result == 'Yes':
                    result = 'https://cancer.sanger.ac.uk/cosmic/census-page/' + csq_dicts[i]['SYMBOL']
                else:
                    result = ''
            result = convert_numpy(result).replace(', ', ',')
            result = result.replace('; ', ',') # database has highly variable entries
            result = result.strip() # remove leading or trailing whitespace
            result = result.replace(' ', '_')
            CancerGeneCensus_list.append(result)
        CancerGeneCensus = ('|').join(CancerGeneCensus_list)
        CancerGeneCensus_lists.append(CancerGeneCensus)
    CancerGeneCensus_line = ','.join(CancerGeneCensus_lists)
    return CancerGeneCensus_line


def modify_header(bcf_in):
    '''
        Add new INFO field
    '''
    bcf_in.header.info.add(id='CancerGeneCensus', number='.',
                           type='String',
                           description='Consequence annotations from Cancer Gene Census. Format: Tier|Hallmark|Somatic|Germline|Tumour Types(Somatic)|Tumour Types(Germline)|Cancer Syndrome|Tissue Type|Molecular Genetics|Role in Cancer|Mutation Types'
                           )
    return bcf_in


def modify_record(record, cancer_gene_census, csq_columns):
    '''
        Add new INFO field to each record
    '''
    record.info['CancerGeneCensus'] = get_CancerGeneCensus(record, csq_columns,
                                                           cancer_gene_census)
    return record


def get_census(cancer_gene_census_file):
    '''
        Read  in the Cancer Gene Census file
    '''
    cancer_gene_census = pd.read_csv(cancer_gene_census_file)
    return cancer_gene_census


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def write_vcf(bcf_in, vcf_out_file, cancer_gene_census, csq_columns):
    '''
        Write out the download
    '''
    bcf_out = pysam.VariantFile(vcf_out_file, 'w', header=bcf_in.header)
    for record in bcf_in.fetch():
        record = modify_record(record, cancer_gene_census, csq_columns)
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print(exit_status)

def prep_resitance_mutations():
    '''
        Convert to new record from tsv file
    '''


def main():
    '''
        Annotates files by adding information about the
        Cosmic Genome Census entry for the gene symbol
    '''
    cancer_gene_census_file = sys.argv[1]
    vcf_file = sys.argv[2]
    vcf_out_file = sys.argv[3]
    cancer_gene_census = get_census(cancer_gene_census_file)
    bcf_in = read_vcf(vcf_file)
    bcf_in = modify_header(bcf_in)
    csq_columns = get_csq_columns(bcf_in)
    write_vcf(bcf_in, vcf_out_file, cancer_gene_census, csq_columns)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
