#!/usr/bin/env python
# USAGE: remove_contig.py VCF_IN VCF_OUT
# DESCRIPTION: Print a VCF file skipping contig descriptions (for use on bad contig descriptions)
# Version 1.0
import sys
import shutil


def remove_contig(vcf):
    '''
        Skips line if starts with ##contig= to remove Lumpy
        VCF line with contig name but witout required length value
    '''
    for line in vcf:
        if not line.startswith('##contig='):
            yield line


def vcf_writer(vcf_file, vcf_out_file):
    '''
       Write out the VCF file with corrected information
    '''
    #  =====================
    #  test if renaming should occur
    #  =====================
    rename = False
    if vcf_out_file == vcf_file:
        vcf_out_file = vcf_file + '_tmp.vcf'
        rename = True
    #  =====================
    #  write non-contig lines
    #  =====================
    with open(vcf_out_file, 'w') as vcf_out:
        with open(vcf_file) as vcf:
            for line in remove_contig(vcf):
                vcf_out.write(line)
    #  =====================
    #  rename output VCF
    #  =====================
    if rename:
        shutil.move(vcf_out_file, vcf_file)



#  =====================
#  Main
#  =====================
if __name__ == "__main__":
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    vcf_writer(vcf_file, vcf_out_file)