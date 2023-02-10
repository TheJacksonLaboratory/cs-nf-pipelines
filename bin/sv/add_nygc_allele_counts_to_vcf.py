#!/usr/bin/env python

################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Version: 0.5 (2018-09-05)
# Author: Kanika Arora (karora@nygenome.org)
##################### COPYRIGHT ################################################
################################################################################

#### PLEASE NOTE THAT THIS SCRIPT EXPECTS THE FORMAT COLUMN FOR THE NORMAL SAMPLE
#### TO BE COLUMN 10 AND TUMOR SAMPLE TO BE COLUMN 11
#### IT WILL ADD INCORRECT ALLELE COUNTS IF THAT ORDER OF SAMPLES IS NOT TRUE 

import pysam
import sys
import argparse
import os
import re

header=dict()

def add_new_header_field(KEY, VALUE):
    if KEY not in header:
        header[KEY]=''
    header[KEY]=header[KEY]+"##{0}={1}\n".format(KEY,VALUE)

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

def check_if_exists(file_or_dir_path, type="file"):
    if type=="file":
        if not os.path.isfile(file_or_dir_path):
            print("ERROR: Required file {0} does not exist. Cannot run.".format(file_or_dir_path))
            sys.exit(1)
    elif type=="directory":
        if not os.path.isdir(file_or_dir_path):
            print("ERROR: Required file {0} does not exist. Cannot run.".format(file_or_dir_path))
            sys.exit(1)
    else:
        if not os.path.exists(file_or_dir_path):
            print("ERROR: {0} does not exist. Cannot run.".format(file_or_dir_path))
            sys.exit(1)

def compute_vaf(alt_count,dp):
    '''
        Compute VAF from dp and alt_count.
    '''
    if not isinstance(alt_count,int):
        if str(int(alt_count)) == alt_count:
            alt_count=int(alt_count)
        else:
            raise ValueError("alt_count should be an integer. "+alt_count+" provided. Cannot run")
    if not isinstance(dp,int):
        if str(int(dp)) == dp:
            dp=int(dp)
        else:
            raise ValueError("dp should be an integer. "+dp+" provided. Cannot run")
    if alt_count > dp:
        raise ValueError("alt_count {0} is greater than depth {1}.".format(alt_count,dp))
    return (0 if dp==0 else round(float(alt_count)/dp,4))


def infer_variant_type(ref, alt):
    '''
        Infers whether the variant is a SNV,MNV,INDEL or COMPLEX (delin).
    '''
    variant_type="COMPLEX"
    if len(ref) == len(alt):
        if len(ref)==1:
            variant_type="SNV"
        else:
            variant_type="MNV"
    elif (len(ref) == 1 or len(alt) == 1) and ref[0] == alt[0]:
        ## VCF file has anchor bases for indels
        variant_type="INDEL"
    return variant_type

def is_too_long(ref, alt, variant_type, MAX_INDEL_LEN):
    '''
    Test if an INDEL or COMPLEX event is too long for computing allele counts using NYGC's pileup method given the length cut off.
    '''
    too_long = False
    if (len(ref) > MAX_INDEL_LEN or len(alt) > MAX_INDEL_LEN) and (variant_type == "INDEL" or variant_type == "COMPLEX"):
        too_long = True
    return too_long


def read_pileup_return_count(samfile, chr, pos, ref,
                             alt, variant_type, MIN_MQ=10, MIN_BQ=10,
                             testing=False):
    ref_reads = []
    alt_reads = []
    other_reads = []
    f1r2_reads = []
    f2r1_reads = []
    properly_paired_reads = []
    fwd_reads = []
    rev_reads = []
    pileup = samfile.pileup(chr, pos - 1, pos)
    possible_complex=False
    anchor_mismatch=0
    for pileupcolumn in pileup:
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                # if the position in the read is .is_del pos is none so take next
                pos_in_read = pileupread.query_position
                # skip reads where the position is already a deletion (is_del)
                if not pos_in_read:
                    continue
#                print ('pos_in_read', pos_in_read, pos, ref, alt)
#                print('VERSION', pysam.__version__)
                #  ==========================
                #  pysam filters secondary, dup, and qcfail by default unless nofilter is used
                #  ==========================
                BQ = pileupread.alignment.query_qualities[pos_in_read]
                MQ = pileupread.alignment.mapping_quality
                if testing:
                    if pileupread.alignment.is_duplicate:
                        print('is_duplicate')
                        sys.exit(0)
                    if pileupread.alignment.is_qcfail:
                        print('is_qcfail')
                        sys.exit(0)
                if BQ < MIN_BQ \
                        or MQ < MIN_MQ \
                        or pileupread.alignment.is_supplementary:
                    continue
                read_name = pileupread.alignment.query_name
                # Check strand of reads
                if pileupread.alignment.is_reverse is False:
                    fwd_reads.append(read_name)
                else:
                    rev_reads.append(read_name)
                #Check if properly paired
                if pileupread.alignment.is_proper_pair is True:
                    properly_paired_reads.append(read_name)
                    #If properly paired, check if read-pair in F1R2 or F2R1 orientation
                    if pileupread.alignment.is_reverse is False and pileupread.alignment.is_read1:
                        f1r2_reads.append(read_name)
                    else:
                        f2r1_reads.append(read_name)
                # filter reads that don't span the indel
                if pos_in_read + len(ref) > pileupread.alignment.query_alignment_end \
                        or pos_in_read + len(alt) > pileupread.alignment.query_alignment_end:
                    continue

                ## Check if read has ref allele or alt allele
                if variant_type == "SNV" or variant_type == "MNV":
                    if pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(ref)] == ref:
                        ref_reads.append(read_name)
                    elif pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(alt)] == alt:
                        alt_reads.append(read_name)
                        if pileupread.indel != 0:
                            anchor_mismatch+=1
                    else:
                        other_reads.append(read_name)
                elif variant_type == "COMPLEX":
                ### Exact length and sequence of allele match required for complex events. Example: if the variant is AC>T, and if a read has deletion of C but the nt at anchor position is A, it will go into other_reads ###
                    if pileupread.indel == 0 and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(ref)] == ref:
                        ref_reads.append(read_name)
                    elif pileupread.indel == len(alt) - len(ref) and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(alt)] == alt:
                        alt_reads.append(read_name)
                    else:
                        other_reads.append(read_name)
                else:   
                    #### Variant type is "INDEL" ####
                    ############################# PLEASE NOTE #######################################
                    ### Variant calling for indels: For insertion, check whether the length of the 
                    ### insertion and sequence matches alt allele. If there is no indel at the anchor
                    ### position (even if the base at the anchor position doesn't match), we consider
                    ### the read as adding support to the reference. For deletions, if the length of
                    ### deletion matches alt allele, we consider that read supporting the alt allele,
                    ### and if there is no deletion at that position (even if there are mismatches in
                    ### the bases spanning the deletion), it's considered to support reference allele.
                    ### Examples for insertions:
                    ### Let's say that the variant is chr1:12345 A > AT
                    ### Scenario1: The read has a C at chr1:12345 along with insertion of T
                    ###            This read will be used to add support to the alternate allele
                    ### Scenario2: Read has a mismatch (let's say 'C') at chr1:12345, but no indel
                    ###            This read will be used to add support to the reference allele
                    ### Scenario3: Read has a different insertion, let's say 'G' insted of 'T'
                    ###            This read will go into the other_reads category.
                    ### Example for deletions:
                    ### Let's say the variant is chr1:12345 AT > A
                    ### Scenario1: The read has a C at chr1:12345 along with deletion of T
                    ###            It will be used to add support to the alt allele.
                    ### Scenario2: The read has a A at chr1:12345 followed by a 2nt deletion
                    ###            This read will go into the other_reads_category
                    ###################################################################################
                    if len(ref)==1:
                        #Variant is an insertion
                        if pileupread.indel == 0: ### and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + 1] == ref:
                            ref_reads.append(read_name)
                        elif pileupread.indel  == len(alt) - len(ref) and pileupread.alignment.query_sequence[pos_in_read+1:pos_in_read + len(alt)] == alt[1:]:
                            alt_reads.append(read_name)
                            if pileupread.alignment.query_sequence[pos_in_read] != alt[0]:
                                anchor_mismatch+=1
                        else:
                            other_reads.append(read_name)
                    else:
                        # Variant is a deletion (len(ref)>1 and len(alt)==1)
                        if pileupread.indel == 0: ## and pileupread.alignment.query_sequence[pos_in_read+1:pos_in_read + len(ref)] == ref[1:]:
                            ref_reads.append(read_name)
                        elif pileupread.indel == len(alt) - len(ref): ## and pileupread.alignment.query_sequence[pos_in_read:pos_in_read + len(alt)] == alt:
                            alt_reads.append(read_name)
                            if pileupread.alignment.query_sequence[pos_in_read] != alt[0]:
                                anchor_mismatch+=1
                        else:
                            other_reads.append(read_name)

    ## If there are more than 2 reads that support the alternate allele of an indel variant, but the anchor base does not match, we report that as a PossiblyComplex event ##
    ## Similarly, if there are more than 2 reads that support alt allele of an SNV, but have an indel immediately following the SNV variant, we report that as a PossiblyComplex event ##
    if anchor_mismatch > 2:
        possible_complex=True
    # check sets to make sure reads don't show up in multiple sets
    # supporting multiple calls
    set_ref_raw = set(ref_reads)
    set_alt_raw = set(alt_reads)
    set_other_raw = set(other_reads)
    ref_reads_set = set_ref_raw - set_alt_raw - set_other_raw
    alt_reads_set = set_alt_raw - set_ref_raw - set_other_raw
    other_reads_set = set_other_raw - set_ref_raw - set_alt_raw
    all_reads_set = alt_reads_set|ref_reads_set|other_reads_set
    # make read-type sets
    f1r2_reads_set = set(f1r2_reads)
    f2r1_reads_set = set(f2r1_reads)
    fwd_reads_set = set(fwd_reads)
    rev_reads_set = set(rev_reads)
    properly_paired_reads_set = set(properly_paired_reads)
    # tally set in ref and alt, non-ref/alt, all reads
    ref_count = len(ref_reads_set)
    alt_count = len(alt_reads_set)
    # other_count=len(other_reads_set) # not used
    dp = len(all_reads_set)
    # get VAF
    vaf = compute_vaf(alt_count, dp)
    # get read set count by pair info ref/alt
    f1r2_ref = len(ref_reads_set.intersection(f1r2_reads_set))
    f1r2_alt = len(alt_reads_set.intersection(f1r2_reads_set))
    f2r1_ref = len(ref_reads_set.intersection(f2r1_reads_set))
    f2r1_alt = len(alt_reads_set.intersection(f2r1_reads_set))
    # get read set count by orientation ref/alt
    fwd_ref = len(ref_reads_set.intersection(fwd_reads_set))
    fwd_alt = len(alt_reads_set.intersection(fwd_reads_set))
    rev_ref = len(ref_reads_set.intersection(rev_reads_set))
    rev_alt = len(alt_reads_set.intersection(rev_reads_set))
    # tally properly paired sets ref/alt
    proper_paired_ref = len(ref_reads_set.intersection(properly_paired_reads_set))
    proper_paired_alt = len(alt_reads_set.intersection(properly_paired_reads_set))
    # tally not properly paired sets ref/alt
    not_proper_paired_ref = len(ref_reads_set) - proper_paired_ref
    not_proper_paired_alt = len(alt_reads_set) - proper_paired_alt
    return ('{0},{1}'.format(ref_count, alt_count),
            str(dp), str(vaf),
            '{0},{1}'.format(f1r2_ref, f1r2_alt),
            '{0},{1}'.format(f2r1_ref, f2r1_alt),
            '{0},{1}'.format(fwd_ref, fwd_alt),
            '{0},{1}'.format(rev_ref,rev_alt),
            '{0},{1}'.format(proper_paired_ref, proper_paired_alt),
            '{0},{1}'.format(not_proper_paired_ref, not_proper_paired_alt),possible_complex)


def __main__():
    parser = ArgumentParser(prog='add_nygc_allele_counts',
                            description='Runs pileup on tumor and normal bam files to compute allele counts for bi-allelic SNV and Indel variants in VCF file and adds pileup format columns to the VCF file.', epilog='',
                            formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=100, width=150))
    parser.add_argument('-t', '--tumor_bam', help = 'Tumor BAM file.', required=True)
    parser.add_argument('-n', '--normal_bam', help = 'Normal BAM file.', required=True)
    parser.add_argument('-v', '--vcf', help = 'SNV VCF file.', required=True)
    parser.add_argument('-o', '--output', help = 'Output VCF file.', required=True)
    parser.add_argument('-b', '--min_base_quality', help='Minimum base quality', default=10, type=int)
    parser.add_argument('-m', '--min_mapping_quality', help='Minimum mapping quality',
                        default=10, type=int)
    parser.add_argument('-i', '--max_indel_len_for_count',
                        help='Maximum indel or delin (complex event) length for generating counts',
                        default=10, type=int)
    args=parser.parse_args()
    # name variables
    TBAM=args.tumor_bam
    NBAM=args.normal_bam
    VCF=args.vcf
    OUT=args.output
    MIN_BQ=args.min_base_quality
    MIN_MQ=args.min_mapping_quality
    MAX_INDEL_LEN=args.max_indel_len_for_count
    # test files
    check_if_exists(TBAM)
    check_if_exists(NBAM)
    check_if_exists(VCF)

    TumorSamFile=pysam.AlignmentFile(TBAM, "rb")
    NormalSamFile=pysam.AlignmentFile(NBAM, "rb")

    f=open(VCF)
    o=open(OUT,"w")
    seen=0
    for line in f:
        if line.startswith("#"):
            if line.startswith("##FORMAT") and seen==0:
                o.write('##INFO=<ID=PossiblyComplex,Number=0,Type=Flag,Description="The variant is probably part of a complex event.">\n')
                o.write('##FORMAT=<ID=nygc_AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_DP,Number=1,Type=Integer,Description="Depth of coverage: Number of reads covering site computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_AF,Number=A,Type=Float,Description="Allele fraction of alternate allele in the sample computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_FWD,Number=R,Type=Integer,Description="Count of reads in forward strand supporting each allele computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_REV,Number=R,Type=Integer,Description="Count of reads in forward strand supporting each allele computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_PROPER_PAIRED,Number=R,Type=Integer,Description="Count of reads that are properly paired supporting each allele computed using a custom pileup-based approach.">\n')
                o.write('##FORMAT=<ID=nygc_NOT_PROPER_PAIRED,Number=R,Type=Integer,Description="Count of reads that are not properly paired (discordant read pairs) supporting each allele computed using a custom pileup-based approach.">\n')
                seen=1
            o.write(line)
        else:
            new_ids=["nygc_AD", "nygc_DP", "nygc_AF", "nygc_F1R2", "nygc_F2R1", "nygc_FWD", "nygc_REV",
                     "nygc_PROPER_PAIRED", "nygc_NOT_PROPER_PAIRED"]
            line=line.strip()
            toks=line.split("\t")
            chrom = toks[0]
            pos = toks[1]
            ref = toks[3]
            alt = toks[4]
            variant_type=infer_variant_type(ref,alt)

            too_long = is_too_long(ref, alt, variant_type, MAX_INDEL_LEN)
            if too_long:
                o.write("\t".join(toks)+"\n")
            else:
                toks[8] = toks[8]+":"+":".join(new_ids)
                (AD,DP,AF,F1R2,F2R1,FWD,REV,PROPER_PAIRED,NOT_PROPER_PAIRED,possibly_complex) = read_pileup_return_count(NormalSamFile, chrom,
                                                                                                        int(pos),
                                                                                                        ref, alt, variant_type,
                                                                                                        MIN_MQ=10, MIN_BQ=10)
                toks[9]=toks[9]+":"+":".join([AD, DP, AF, F1R2, F2R1, FWD, REV,
                                              PROPER_PAIRED, NOT_PROPER_PAIRED])
                (AD,DP,AF,F1R2,F2R1,FWD,REV,PROPER_PAIRED,NOT_PROPER_PAIRED, possibly_complex) = read_pileup_return_count(TumorSamFile, chrom,
                                                                                                        int(pos),
                                                                                                        ref, alt, variant_type,
                                                                                                        MIN_MQ=10, MIN_BQ=10)
                toks[10]=toks[10]+":"+":".join([AD, DP, AF, F1R2, F2R1, FWD, REV,
                                                PROPER_PAIRED, NOT_PROPER_PAIRED])
                if possibly_complex is True:
                    toks[7]=toks[7]+";PossiblyComplex"
                o.write("\t".join(toks)+"\n")
    f.close()
    o.close()

if __name__ == "__main__":
    __main__()