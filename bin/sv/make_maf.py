#!/usr/bin/env python
#	USAGE: python make_maf.py VCF MAF LIBRARY GENOME
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
import re
import logging as log
import pysam
import mygene
import argparse
import pandas as pd
##########################################################################
##############                  Custom functions              ############
##########################################################################


def ensembl_gene_id_entrez_id(ensembl_gene_id, mg):
    '''
        Returns entrez id from ensemble.
        False is used for regions that do not
        correspond to a gene region or Ensembl ID
    '''
    entrez_id = 0
    if ensembl_gene_id != '':
        results = mg.query(ensembl_gene_id)
        try:
            entrez_id = str(results['hits'][0]['entrezgene'])
        except (KeyError, IndexError):
            pass
#            sys.stderr.write('WARNING: entrezgene not found for ' + str(ensembl_gene_id) + ' in ' + str(results) + '\n')
    return entrez_id


def ensembl_gene_entrez_local(ensembl_gene_id, data):
    '''
        Returns entrez id from ensemble.
        False is used for regions that do not
        correspond to a gene region or Ensembl ID
        '''
    try:
        entrez_id = data[(data['Gene stable ID'] == ensembl_gene_id)]['NCBI gene (formerly Entrezgene) ID'].values[0]
    except IndexError:
        entrez_id = 0
    if str(entrez_id) == 'nan':
        entrez_id = 0
    return int(entrez_id)

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


def group_mnv(record):
    '''
        Convert inhouse MNV and SNV type to GDC-like calls.
        This does not work for multiple alts.
        Only the first alt will be considered because
        type only takes one value.
    '''
    if record.info['TYPE'] == 'SNV':
        record.info['TYPE'] = 'SNP'
    elif record.info['TYPE'] == 'MNV':
        if len(record.ref) == 2 \
                and len(record.alts[0]) == 2:
            record.info['TYPE'] = 'DNP'
        elif len(record.ref) == 3 \
                and len(record.alts[0]) == 3:
            record.info['TYPE'] = 'TNP'
        elif len(record.ref) > 3 \
                and len(record.alts[0]) > 3 \
                and len(record.ref) == len(record.alts[0]):
            record.info['TYPE'] = 'ONP'
    return record



def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
    '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def get_dbsnp_rs(Existing_variation):
    '''
        Remove Cosmic IDs and split by comma.
    '''
    ids = Existing_variation.split('&')
    good_ids = [id for id in ids if id.startswith('rs')]
    return ','.join(good_ids)


def set_frame_shift(type):
    '''
        Set the variant classification based on whether the variant is
        an insertion or a deletion.
    '''
    if type == 'DEL':
        variant_classification = 'Frame_Shift_Del'
    elif type == 'INS':
        variant_classification = 'Frame_Shift_Ins'
    else:
        variant_classification = False
    return variant_classification


def set_protein_altering(type, ref, alt):
    '''
        Set the variant classification based on whether the variant is
        an insertion or a deletion.
    '''
    inframe = False
    if abs(len(ref) - len(alt)) % 3 == 0:
        inframe = True
    if inframe:
        if type == 'DEL':
            variant_classification = 'In_Frame_Del'
        elif type == 'INS':
            variant_classification = 'In_Frame_Ins'
        else:
            variant_classification = False
    else:
        if type == 'DEL':
            variant_classification = 'Frame_Shift_Del'
        elif type == 'INS':
            variant_classification = 'Frame_Shift_Ins'
        else:
            variant_classification = False
    return variant_classification


def shorten_AA(AAMutation):
    '''
        Lengthen to the three letter AA code
    '''
    AA_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
            'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
            'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
            'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M',
            'Ter' : '*'}
    short_mutation = []
    skip_until = -1
    for i, char in enumerate(AAMutation):
        if i > skip_until:
            if AAMutation[i:i + 3] in AA_dict:
                short_mutation += AA_dict[AAMutation[i:i + 3]]
                skip_until = i + 2
            else:
                short_mutation += char
    return ''.join(short_mutation)


def get_variant_classification(consequences, type, ref, alt):
    '''
        Convert VEP consequences to MAF variant_classification.
    '''
    consequences_to_class = {'intergenic_variant' : 'Silent',
                                'upstream_gene_variant' : 'Silent',
                                '5_prime_UTR_variant' : 'Silent',
                                'splice_acceptor_variant' : 'Splice_Site',
                                'splice_donor_variant' : 'Splice_Site',
                                'splice_region_variant' : 'Splice_Site',
                                'missense_variant' : 'Missense_Mutation',
                                'synonymous_variant' : 'Silent',
                                'frameshift_variant' : set_frame_shift(type),
                                'protein_altering_variant' : set_protein_altering(type, ref, alt),
                                'inframe_insertion' : 'In_Frame_Ins',
                                'inframe_deletion' : 'In_Frame_Del',
                                'stop_gained' : 'Nonsense_Mutation',
                                'stop_retained_variant' : 'Silent',
                                'stop_lost' : 'Nonstop_Mutation',
                                'intron_variant' : 'Silent',
                                '3_prime_UTR_variant' : 'Silent',
                                'downstream_gene_variant' : 'Silent',
                                'initiator_codon_variant' : 'Translation_Start_Site',
                                'regulatory_region_variant' : 'Silent',
                                'TF_binding_site_variant' : 'Silent',
                                'mature_miRNA_variant' : 'RNA',
                                'regulatory_region_ablation' : 'Silent',
                                'regulatory_region_amplification' : 'Silent',
                                'TFBS_ablation' : 'Silent',
                                'TFBS_amplification' : 'Silent',
                                'non_coding_transcript_variant' : 'Silent',
                                'NMD_transcript_variant' : 'Silent',
                                'incomplete_terminal_codon_variant' : 'Silent',
                                'non_coding_transcript_exon_variant' : 'RNA',
                                'transcript_ablation' : 'Splice_Site',
                                'transcript_amplification' : 'Silent',
                                'feature_elongation' : False,
                                'feature_truncation' : False,
                                'start_lost' : 'Translation_Start_Site',
                                'start_retained_variant' : 'Silent',
                                'coding_sequence_variant' : 'Missense_Mutation'
                            }
    return consequences_to_class[consequences]

def get_HGVSp_Short(HGVSp_string, HGVSc_string, csq_term):
    '''
        Convert HGVSp to HGVSp_Short:
        derive HGVSp_Short from HGVSp. if Consequence is splice acceptor/donor variants, generate HGVSp_Short
    '''
    aa_to_short = {'Ala': 'A',
                   'Arg': 'R',
                   'Asn': 'N',
                   'Asp': 'D',
                   'Asx': 'B',
                   'Cys': 'C',
                   'Glu': 'E',
                   'Gln': 'Q',
                   'Glx': 'Z',
                   'Gly': 'G',
                   'His': 'H',
                   'Ile': 'I',
                   'Leu': 'L',
                   'Lys': 'K',
                   'Met': 'M',
                   'Phe': 'F',
                   'Pro': 'P',
                   'Ser': 'S',
                   'Thr': 'T',
                   'Trp': 'W',
                   'Tyr': 'Y',
                   'Val': 'V',
                   'Xxx': 'X',
                   'Xaa': 'X',
                   'Ter': '*'
                  }
    HGVSp_Short = ''
    if csq_term == 'splice_acceptor_variant' \
            or csq_term == 'splice_donor_variant':
        if len(HGVSc_string.split(':'))>1:
            HGVSc_string = HGVSc_string.split(':')[1]
            HGVSc_coding = re.findall('^c.(\d+)',HGVSc_string)
            if len(HGVSc_coding) > 0:
                input_pos = float(HGVSc_coding[0])
                if input_pos < 1:
                    input_pos = 1
                corrected_pos = (input_pos + input_pos % 3)/3
                HGVSp_Short = 'p.X' + str(int(corrected_pos)) + '_splice'
        return HGVSp_Short  
    elif len(HGVSp_string) > 0:
        HGVSp_Short = HGVSp_string.split(':')[1]
        for item in aa_to_short.keys():
            HGVSp_Short = re.sub(item, aa_to_short[item], HGVSp_Short)
        return HGVSp_Short
    else:
        return HGVSp_Short
    
def make_row(record, csq_columns, bcf_in, library,
             normal, tumor, VEP_version='GRCh38',
             mg=False,
             ensembl_entrez=False):
    '''
        Fill in MAF row
    '''
    #  =======================
    #  Get VEP annotation
    #  =======================
    csq_dicts = get_csqs(record, csq_columns)
    cosmic_resistance_annotation = {}
    if normal == bcf_in.header.samples[1]:
        normal_index = 1
        tumor_index = 0
    elif tumor == bcf_in.header.samples[1]:
        normal_index = 0
        tumor_index = 1
    for i, alt in enumerate(record.alts):
        if csq_dicts[i]['SYMBOL_SOURCE'] == 'HGNC':
            hugo = csq_dicts[i]['SYMBOL']
        else:
            hugo = 'Unknown'
        ensembl_gene_id = csq_dicts[i]['Gene']
        if mg:
            entrez_id = ensembl_gene_id_entrez_id(ensembl_gene_id, mg) # default '0'
        else:
            entrez_id = ensembl_gene_entrez_local(ensembl_gene_id, ensembl_entrez)
        center = 'NYGenome'
        ncbi_build = VEP_version
        chrom = record.chrom
        if record.info['TYPE'] == 'DEL':
            start = record.pos + 1 # skip anchor base
        else:
            start = record.pos
        # get end position for 1-based inclusive coordinates
        if record.info['TYPE'] == 'SNP':
            end = record.pos
        if record.info['TYPE'] == 'INS':
            end = record.pos + 1
        elif record.info['TYPE'] == 'DEL':
            end = (record.pos + 1) + len(record.ref) - len(alt) - 1 # add one to skip anchor
        else:
            end = record.pos + len(record.ref) - 1
        strand = '+'
        record = group_mnv(record)
        variant_type = record.info['TYPE']
        if record.info['TYPE'] == 'INS':
            reference_allele = '-'
            Tumor_Seq_Allele1 = '-'
            Tumor_Seq_Allele2 = alt[1:]
        elif record.info['TYPE'] == 'DEL':
            reference_allele = record.ref[1:]
            Tumor_Seq_Allele1 = record.ref[1:]
            Tumor_Seq_Allele2 = '-'
        else:
            reference_allele = record.ref
            Tumor_Seq_Allele1 = record.ref
            Tumor_Seq_Allele2 = alt
        dbSNP_RS = get_dbsnp_rs(csq_dicts[i]['Existing_variation'])
        dbSNP_Val_Status = 'bySubmitter'
        Tumor_Sample_Barcode = bcf_in.header.samples[tumor_index]
        Matched_Norm_Sample_Barcode = bcf_in.header.samples[normal_index]
        Match_Norm_Seq_Allele1 = ''
        Match_Norm_Seq_Allele2 = ''
        Tumor_Validation_Allele1 = ''
        Tumor_Validation_Allele2 = ''
        Match_Norm_Validation_Allele1 = ''
        Match_Norm_Validation_Allele2 = ''
        Verification_Status = 'Unknown'
        Validation_Status = 'Untested'
        Mutation_Status = 'Somatic'
        Sequencing_Phase = 'Phase_I'
        Sequence_Source = library
        Validation_Method = 'none'
        Score = ''
        BAM_file= ''
        Sequencer = 'Illumina'
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
        HGVSc = csq_dicts[i]['HGVSc']
        HGVSp = csq_dicts[i]['HGVSp']
        SYMBOL_SOURCE = csq_dicts[i]['SYMBOL_SOURCE']
        SYMBOL = csq_dicts[i]['SYMBOL']
        IMPACT = csq_dicts[i]['IMPACT']
        return_line = []
        for csq_term in csq_dicts[i]['Consequence'].split('&'):
            variant_classification = get_variant_classification(csq_term,
                                                            record.info['TYPE'],
                                                            record.ref,
                                                            alt)
            HGVSp_Short = get_HGVSp_Short(HGVSp,HGVSc,csq_term)
            line = [hugo, entrez_id, center, ncbi_build, chrom, start, end, strand,
                    variant_classification, variant_type, reference_allele,
                    Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS,
                    dbSNP_Val_Status, Tumor_Sample_Barcode,
                    Matched_Norm_Sample_Barcode, Match_Norm_Seq_Allele1,
                    Match_Norm_Seq_Allele2, Tumor_Validation_Allele1,
                    Tumor_Validation_Allele2,
                    Match_Norm_Validation_Allele1, Match_Norm_Validation_Allele2,
                    Verification_Status, Validation_Status, Mutation_Status,
                    Sequencing_Phase,
                    Sequence_Source, Validation_Method, Score, BAM_file,
                    Sequencer, t_alt_count, t_ref_count,
                    n_alt_count, n_ref_count, HGVSc, HGVSp, HGVSp_Short, SYMBOL, SYMBOL_SOURCE, IMPACT]
            joined_line = '\t'.join([str(x) for x in line]) + '\n'
            joined_line = joined_line.replace('&', ',')
            return_line.append(joined_line)
        yield ''.join(set(return_line))


def write_file(bcf_in, out, csq_columns, library,
               normal, tumor, VEP_version, ensembl_entrez=False):
    '''
        Write out the header
    '''
    if ensembl_entrez:
        data = pd.read_csv(ensembl_entrez)
        mg = False
    else:
        data = False
        mg = mygene.MyGeneInfo()
    with open(out, 'w') as o:
        header = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_file', 'Sequencer', 't_alt_count', 't_ref_count', 'n_alt_count', 'n_ref_count', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'SYMBOL', 'SYMBOL_SOURCE','IMPACT']
        o.write('\t'.join(header) + '\n')
        for record in bcf_in.fetch():
            if record.info['HighConfidence']:
                for alt_line in make_row(record, csq_columns, bcf_in, library,
                                         normal, tumor, VEP_version=VEP_version,
                                         mg=mg,
                                         ensembl_entrez=data):
                    o.write(alt_line)



def main():
    '''
        Script to make MAF file from VEP annotated VCF files
    '''
    ######################################################################
    ############        Get commandline arguments             ############
    ######################################################################
    parser = argparse.ArgumentParser(
                                     description='DESCRIPTION: Takes in a VCF \
                                     file and returns a MAF. Command-line \
                                     options that may be omitted (i.e. are NOT \
                                     required) are shown in square brackets.')
                                     #    Documentation parameters
    #    Parameter options
    parser.add_argument('-v', '--vcf',
                     dest='vcf_file',
                     help='Annotated VCF file')
    parser.add_argument('-m', '--maf',
                     dest='maf',
                     help='MAF output file')
    parser.add_argument('-l', '--library',
                     dest='library',
                     help='Sequence library type',
                     choices=['WGS', 'Exome'])
    parser.add_argument('-vep', '--vep-version',
                     dest='VEP_version',
                     help='VEP genome version',
                     choices=['GRCh37', 'GRCh38'])
    parser.add_argument('-t', '--tumor',
                     dest='tumor',
                     help='Tumor sample name')
    parser.add_argument('-n', '--normal',
                        dest='normal',
                        help='Normal sample name')
    parser.add_argument('-e', '--ensembl-entrez',
                        dest='ensembl_entrez',
                        default=False,
                        help='Map of ensembl ids to entrez ids')
    args = parser.parse_args()
    assert os.path.isfile(args.vcf_file), 'Failed to find caller VCF call file :' + args.vcf_file
    bcf_in = read_vcf(args.vcf_file)
    csq_columns = get_csq_columns(bcf_in)
    write_file(bcf_in, args.maf, csq_columns, args.library,
               args.normal, args.tumor, args.VEP_version,
               args.ensembl_entrez)


##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################
if __name__ == '__main__':
    main()
