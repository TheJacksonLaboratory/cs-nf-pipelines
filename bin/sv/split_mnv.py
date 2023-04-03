#!/usr/bin/env python
#	USAGE: python split_mnv.py
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


def get_type(record):
    '''
        Fill in the type field.
    '''
    if len(record.ref) == 1 and len(record.alts[0]) == 1:
        type = 'SNV'
    elif len(record.alts) > 1:
        type = 'MULTI'
    elif len(record.ref) == 1 and len(record.alts[0]) > 1 and record.ref[0] == record.alts[0][0]:
        type = 'INS'
    elif len(record.ref) > 1 and len(record.alts[0]) == 1 and record.ref[0] == record.alts[0][0]:
        type = 'DEL'
    else:
        type = 'COMPLEX'
    return type


def print_record(record, bcf_out):
    '''
        Write out a vcf record.
    '''
    exit_status = bcf_out.write(record)
    if exit_status != 0:
        print(exit_status)


def determine_anchor_base(record):
    '''
        determine anchor base
    '''
    # without anchor base
    start_pos = 0
    # with anchor base
    for alt in record.alts:
        if alt[0] == record.ref[0]:
            start_pos = 1
    return start_pos


def write_file(bcf_in, bcf_out, tool):
    '''
       Split MNV records
    '''
    splits = []
    for record in bcf_in.fetch():
        if len(record.alts) > 1:
            log.error('VCF file must have only one ALT per line')
            sys.exit(1)
        if len(record.ref) > 1 \
                and len(record.ref) == len(record.alts[0]):
            start_pos = determine_anchor_base(record)
            mnv_id = '_'.join([tool, record.contig, str(record.pos),
                               record.ref, record.alts[0]])
            refs = list(record.ref)
            alts = list(record.alts[0])
            orig_record = record.copy()
            # print full record
            record.pos = record.pos + start_pos
            record.ref = record.ref[start_pos:]
            record.alts = [record.alts[0][start_pos:]]
            record.info['TYPE'] = 'MNV'
            record.info['MNV_ID'] = [mnv_id] # change if never multi-allelic
            print_record(record, bcf_out)
            # print split records
            for new_record in split_records(start_pos, refs,
                                            orig_record, alts,
                                            mnv_id):
                print_record(new_record, bcf_out)
                splits.append('|'.join([new_record.chrom, str(new_record.pos), new_record.ref, new_record.alts[0]]))
    return set(splits)


def write_file_non_mnv(bcf_in, splits, bcf_out, tool):
    '''
       Split MNV records
    '''
    for record in bcf_in.fetch():
        if not (len(record.ref) > 1 \
                and len(record.ref) == len(record.alts[0])):
            # print non-MNVs
            id = '|'.join([record.chrom, str(record.pos), record.ref, record.alts[0]])
            if not id in splits:
                record.info['TYPE'] = get_type(record) # change if multi-allelic
                print_record(record, bcf_out)
            else:
                print('duplicate...', id)
    return True

def split_records(start_pos, refs, orig_record, alts, mnv_id):
    '''
        print split records (skipping the first anchor base)
        (if start_pos is 1)
    '''
    for i in range(start_pos, len(refs)):
        if refs[i] != alts[i]:
            new_record = orig_record.copy()
            new_record.ref = refs[i]
            new_record.alts = [alts[i]]
            new_record.pos = orig_record.pos + i
            new_record.info['TYPE'] = 'SNV' # change if multi-allelic
            new_record.info['MNV_ID'] = [mnv_id] # change if never multi-allelic
            yield new_record


def main():
    '''
        Prepare the VCF file for merging by:
        1) Split MNV records to one line per nucleotide
        2) Skip any line that has an SNV called by an MNV for the same tool
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    vcf_file = sys.argv[1]
    out = sys.argv[2]
    tool = sys.argv[3]
    
    assert os.path.isfile(vcf_file), 'Failed to find caller VCF call file :' + vcf_file
    #  ==========================
    #  Run prep
    #  ==========================
    bcf_in = read_vcf(vcf_file)
    bcf_in = add_info_header(bcf_out=bcf_in,
                              id='MNV_ID',
                              number='.',
                              type='String',
                              description='ID of multi-nucleotide variant (MNV) that the SNV is part of')
    bcf_in = add_info_header(bcf_out=bcf_in,
                           id='TYPE',
                           number='1',
                           type='String',
                           description='Variant type (SNV,INS,DEL,MNV,COMPLEX,MULTI)')
    bcf_out = pysam.VariantFile(out, 'w', header=bcf_in.header)
    splits = write_file(bcf_in,
                          bcf_out,
                          tool)
    write_file_non_mnv(bcf_in,
                       splits,
                       bcf_out,
                       tool)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()