#!/usr/bin/env python
    #	USAGE: python add_cancer_resistance_mutations.py cosmic_resistance_file genome vcf_file out
#   DESCRIPTION: Annotates from Cosmic DB of resistance mutations.
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
import pandas as pd
import re
##########################################################################
##############                  Custom functions              ############
##########################################################################


def modify_header(bcf_in):
    '''
        Add new INFO field
        ['MUTATION_ID', 'GENOMIC_MUTATION_ID', 'Drug Name', 'Tier']
    '''
    bcf_in.header.info.add(id='CosmicResistanceMutation', number='.',
                           type='String',
                           description='Consequence annotations from Cancer Gene Census. Format: MUTATION_ID|GENOMIC_MUTATION_ID|Drug Name|Tier')
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


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def read_cosmic(cosmic_resistance_file):
    '''
        Read in annotated VCF file.
        '''
    cosmic_resistance = pd.read_csv(cosmic_resistance_file, sep='\t')
    return cosmic_resistance


def get_record(match_Mutation, key, region_cosmic_resistance, match_key):
    '''
        Get value for field from Cancer Gene Census records
    '''
    try:
        value = region_cosmic_resistance[(region_cosmic_resistance[match_key] == match_Mutation)][key].values[0]
    except (KeyError, IndexError):
        value = ''
    return value


def match_cosmic(cosmic_resistance, record, csq_columns):
    '''
        Check if HGVSp or HGVSc is at the same position
        in the Cosmic Restance database return values
        for each ALT in list.
    '''
    #  =======================
    #  Get VEP annotation
    #  =======================
    csq_dicts = get_csqs(record, csq_columns)
    cosmic_resistance_annotation = {}
    for i, alt in enumerate(record.alts):
        cosmic_resistance_annotation[i] = {}
        for key in ['MUTATION_ID', 'GENOMIC_MUTATION_ID', 'Drug Name', 'Tier']:
            cosmic_resistance_annotation[i][key] = ''
        #  =======================
        #  Check for id match
        #  =======================
        kind = 'CosmicCoding'
        id = csq_dicts[i][kind]
        legacy_id = csq_dicts[i][kind +'_LEGACY_ID']
        match = cosmic_resistance[cosmic_resistance.GENOMIC_MUTATION_ID == id].copy()
        legacy_match = cosmic_resistance[cosmic_resistance.LEGACY_MUTATION_ID == legacy_id].copy()
        # ID Mutation
        if not match.empty:
            for key in ['MUTATION_ID', 'GENOMIC_MUTATION_ID', 'Drug Name', 'Tier']:
                cosmic_resistance_annotation[i][key] = match[key].tolist()[0]
                
        elif not legacy_match.empty:
            for key in ['MUTATION_ID', 'GENOMIC_MUTATION_ID', 'Drug Name', 'Tier']:
                cosmic_resistance_annotation[i][key] = legacy_match[key].tolist()[0]
        else:
            #  =======================
            #  Check for non coding id match
            #  =======================
            kind = 'CosmicNonCoding'
            id = csq_dicts[i][kind]
            legacy_id = csq_dicts[i][kind +'_LEGACY_ID']
            match = cosmic_resistance[cosmic_resistance.GENOMIC_MUTATION_ID == id].copy()
            legacy_match = cosmic_resistance[cosmic_resistance.LEGACY_MUTATION_ID == legacy_id].copy()
            # ID Mutation
            if not match.empty:
                for key in ['MUTATION_ID', 'GENOMIC_MUTATION_ID', 'Drug Name', 'Tier']:
                    cosmic_resistance_annotation[i][key] = match[key].tolist()[0]
            elif not legacy_match.empty:
                for key in ['MUTATION_ID', 'GENOMIC_MUTATION_ID', 'Drug Name', 'Tier']:
                    cosmic_resistance_annotation[i][key] = legacy_match[key].tolist()[0]
    #  =======================
    #  Compose annotation text
    #  =======================
    cosmic_resistance_annotation_lines = []
    for i in cosmic_resistance_annotation:
        to_join = [str(cosmic_resistance_annotation[i][annotation]) for annotation in ['MUTATION_ID', 
                                                                                       'GENOMIC_MUTATION_ID', 
                                                                                       'Drug Name', 'Tier']]
        cosmic_resistance_annotation_lines.append('|'.join(to_join))
    cosmic_resistance_annotation_line = ','.join(cosmic_resistance_annotation_lines)
    return cosmic_resistance_annotation_line


def write_file(bcf_in, out, csq_columns,
               cosmic_resistance):
    '''
        Write out the header
    '''
    bcf_out = pysam.VariantFile(out, 'w', header=bcf_in.header)
    for record in bcf_in.fetch():
        cosmic_resistance_annotation = match_cosmic(cosmic_resistance, record, csq_columns)
        record.info['CosmicResistanceMutation'] = cosmic_resistance_annotation
        exit_status = bcf_out.write(record)
        if exit_status != 0:
            print('exit_status', exit_status)


def main():
    cosmic_resistance_file = sys.argv[1]
    vcf_file = sys.argv[2]
    out = sys.argv[3]
    assert os.path.isfile(vcf_file), 'Failed to find caller VCF call file :' + vcf_file
    cosmic_resistance = read_cosmic(cosmic_resistance_file)
    bcf_in = read_vcf(vcf_file)
    bcf_in = modify_header(bcf_in)
    csq_columns = get_csq_columns(bcf_in)
    write_file(bcf_in, out, csq_columns,
               cosmic_resistance)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
