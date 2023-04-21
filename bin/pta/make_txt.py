#!/usr/bin/env python
#	USAGE: python make_maf.py --vcf VCF --txt TXT -n NORMAL -t TUMOR
#   DESCRIPTION: Makes MAF file from VCF file.
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
import argparse
##########################################################################
##############                  Custom functions              ############
##########################################################################


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
    csq_values = []
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


def check_build(bcf_in):
        '''
            Check if genome is in a list of supprted non-human genomes.
        '''
        VEP_line = [metadata.value for metadata in bcf_in.header.records if metadata.key == 'VEP'][0]
        vep_info  = {entry.split('=')[0] : entry.split('=')[-1]  for entry in VEP_line.split(' ')}
        if vep_info['assembly'] in ['"GRCm38.p6"']:
            return False
        else:
            return True


def make_row(record, csq_columns, bcf_in,
             normal, tumor, human=True):
    '''
        Fill in MAF row
    '''
    #  =======================
    #  Get VEP annotation
    #  =======================
    csq_dicts = get_csqs(record, csq_columns)
    if normal == bcf_in.header.samples[1]:
        normal_index = 1
        tumor_index = 0
    elif tumor == bcf_in.header.samples[1]:
        normal_index = 0
        tumor_index = 1
    else:
        log.error('VCF sample names do not match listed tumor or normal name')
        sys.exit(1)
    id = record.id
    if record.id == None:
        id = '.'
    for i, alt in enumerate(record.alts):
        consequence = csq_dicts[i]['Consequence']
        impact = csq_dicts[i]['IMPACT']
        GENE_SYMBOL = csq_dicts[i]['SYMBOL']
        HGVSc = csq_dicts[i]['HGVSc']
        HGVSp = csq_dicts[i]['HGVSp']
        type = record.info['TYPE']
        if human:
            PolyPhen = csq_dicts[i]['Polyphen2_HVAR_pred']
            AF_1000G = csq_dicts[i]['AF_1000G']
            GnomadExomes_AF = csq_dicts[i]['GnomadExomes_AF']
            GnomadGenomes_AF = csq_dicts[i]['GnomadGenomes_AF']
            CosmicCoding = csq_dicts[i]['CosmicCoding']
            CosmicCoding_AA = csq_dicts[i]['CosmicCoding_AA']
            CosmicNonCoding = csq_dicts[i]['CosmicNonCoding']
            fathmm = csq_dicts[i]['FATHMM_pred']
            fathmm_MKL_coding = csq_dicts[i]['fathmm-MKL_coding_pred']
        sift = csq_dicts[i]['SIFT_pred']
        sift_4g = csq_dicts[i]['SIFT4G_pred']
        HighConfidence = record.info['HighConfidence']
        if 'called_by' in record.info:
            called_by = ','.join(record.info['called_by'])
        else:
            called_by = ''
        if 'supported_by' in record.info:
            supported_by = ','.join(record.info['supported_by'])
        else:
            supported_by = ''
        if 'AD' in record.samples[bcf_in.header.samples[1]].keys() \
                and 'AD' in record.samples[bcf_in.header.samples[0]].keys():
            t_alt_count = record.samples[bcf_in.header.samples[tumor_index]]['AD'][1]
            t_ref_count = record.samples[bcf_in.header.samples[tumor_index]]['AD'][0]
            n_alt_count = record.samples[bcf_in.header.samples[normal_index]]['AD'][1]
            n_ref_count = record.samples[bcf_in.header.samples[normal_index]]['AD'][0]
        else:
            t_alt_count = ''
            t_ref_count = ''
            n_alt_count = ''
            n_ref_count = ''
        if 'AF' in record.samples[bcf_in.header.samples[1]].keys():
            t_VAF = record.samples[bcf_in.header.samples[tumor_index]]['AF'][i]
        else:
            t_VAF = ''
        if human:
            line = [record.chrom, record.pos, id, record.ref, alt,
                    consequence, impact, GENE_SYMBOL, HGVSc, HGVSp, type,
                    PolyPhen, AF_1000G, GnomadExomes_AF, GnomadGenomes_AF,
                    CosmicCoding, CosmicCoding_AA, CosmicNonCoding,
                    n_ref_count, n_alt_count, t_ref_count, t_alt_count,
                    t_VAF, fathmm, fathmm_MKL_coding, sift, sift_4g,
                    HighConfidence, called_by, supported_by]
        else:
            line = [record.chrom, record.pos, id, record.ref, alt,
                    consequence, impact, GENE_SYMBOL, HGVSc, HGVSp, type,
                    n_ref_count, n_alt_count, t_ref_count, t_alt_count,
                    t_VAF, sift, sift_4g, HighConfidence, called_by, supported_by]
        line = [str(part).replace('&', ',') for part in line]
        line = [str(part).replace(';', ',') for part in line]
        yield line


def write_file(bcf_in, out, csq_columns,
               normal, tumor, human=True):
    '''
        Write out the header
    '''
    with open(out, 'w') as o:
        if human:
            header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Consequence', 'IMPACT',
                      'GENE_SYMBOL', 'HGVSc', 'HGVSp', 'TYPE', 'PolyPhen', 'AF_1000G',
                      'GnomadExomes_AF', 'GnomadGenomes_AF', 'CosmicCoding',
                      'CosmicCoding_AA', 'CosmicNonCoding',
                      'n_ref_count', 'n_alt_count', 't_ref_count', 't_alt_count',
                      't_VAF', 'FATHMM', 'fathmm_MKL_coding', 'SIFT', 'SIFT4G', 'HighConfidence',
                      'called_by', 'supported_by']
        else:
            header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'Consequence', 'IMPACT',
                      'GENE_SYMBOL', 'HGVSc', 'HGVSp', 'TYPE',
                      'n_ref_count', 'n_alt_count', 't_ref_count', 't_alt_count',
                      't_VAF', 'SIFT', 'SIFT4G', 'HighConfidence',
                      'called_by', 'supported_by']
        o.write('\t'.join(header) + '\n')
        header_len = len(header)
        for record in bcf_in.fetch():
            for alt_line in make_row(record, csq_columns, bcf_in,
                                     normal, tumor, human):
                assert len(alt_line) == header_len, "columns don't equal header names"
                joined_line = '\t'.join([x for x in alt_line]) + '\n'
                o.write(joined_line)



def main():
    '''
        Script to make TEXT file from VEP annotated VCF files
    '''
    ######################################################################
    ############        Get commandline arguments             ############
    ######################################################################
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and returns a TEXT file. Command-line \
                                     options that may be omitted (i.e. are NOT \
                                     required) are shown in square brackets.')
                                     #    Documentation parameters
    #    Parameter options
    parser.add_argument('-v', '--vcf',
                     dest='vcf_file',
                     help='Annotated VCF file')
    parser.add_argument('--txt',
                     dest='txt',
                     help='TEXT output file')
    parser.add_argument('-t', '--tumor',
                     dest='tumor',
                     help='Tumor sample name')
    parser.add_argument('-n', '--normal',
                        dest='normal',
                        help='Normal sample name')
    args = parser.parse_args()
    assert os.path.isfile(args.vcf_file), 'Failed to find caller VCF call file :' + args.vcf_file
    bcf_in = read_vcf(args.vcf_file)
    human = check_build(bcf_in)
    csq_columns = get_csq_columns(bcf_in)
    write_file(bcf_in, args.txt, csq_columns,
               args.normal, args.tumor, human=human)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
