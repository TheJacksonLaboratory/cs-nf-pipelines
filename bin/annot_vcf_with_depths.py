#!/usr/bin/env python

"""Annotate VCF file

Usage: annot_vcf_with_depths.py -v original.vcf -d depths.bed

Add info fields for depths from individual callers to original VCF file.
Based on matching annotations from VCF and annotated bedpe.

"How to add custom fields to the VCF files?" idea from: https://www.biostars.org/p/209188/
"""

__author__ = str(('Brian Sanderson', 'Mike Lloyd', 'Benjamin Leopold'))

import os
import argparse
import pysam

def parse_bed(depths_dict, bed_handle):
    """Parse BED files containing intersections between structural variants
    and annotated exons, and return a dictionary of boolean values for each
    named SV
    """
    if os.path.exists(bed_handle) and os.stat(bed_handle).st_size == 0:
        print("No records in " + str(bed_handle))
        return depths_dict

    elif os.path.exists(bed_handle):
        with open(bed_handle, 'r') as ins_f:
            for line_handle in ins_f:
                line = line_handle.rstrip().split('\t')
                sv_name = line[2]
                if sv_name in depths_dict:
                    error_string = str(sv_name) + " on Chr" + str(line[0]) + \
                        " at " + str(line[1]) + " has multiple SV records"
                    #print(error_string)
                else:
                    if line[3] == 'NA':
                        NanoSV_DR = ['NA', 'NA']
                        NanoSV_DV = ['NA', 'NA']
                    else:
                        NanoSV_DR = [str(y) for y in line[3].split(';')]
                        NanoSV_DV = [str(y) for y in line[4].split(';')]
                    sniffles_DR = [str(y) for y in line[5].split(';')]
                    sniffles_DV = [str(y) for y in line[6].split(';')]
                    depths_dict[sv_name] = [NanoSV_DR, NanoSV_DV, sniffles_DR, sniffles_DV]
        return depths_dict
    else:
        raise FileExistsError


def annotate_vcf(input_vcf, depths, output_vcf) :
    """Iterate through VCF records and append a new INFO field InExon with
    a boolean value that represents whether the variant overlaps with
    annotated exons
    """

    depths_dict = {}
    
    depths_dict = parse_bed(depths_dict, depths)

    myvcf = pysam.VariantFile(input_vcf,'r')
    myvcf.header.info.add('NanoSV_DR','2','String',
                        'Number of reference reads from NanoSV')
    myvcf.header.info.add('NanoSV_DV','2','String',
                        'Number of variant reads from NanoSV')
    myvcf.header.info.add('sniffles_DR','1','String',
                        'Number of reference reads from sniffles')
    myvcf.header.info.add('sniffles_DV','1','String',
                        'Number of variant reads from sniffles')
    with open(output_vcf, 'w') as annot_vcf:
        print(myvcf.header, end='', file=annot_vcf)
        for variant in myvcf:
            if variant.id:
                if variant.id in depths_dict:
                    variant.info['NanoSV_DR'] = depths_dict[variant.id][0]
                    variant.info['NanoSV_DV'] = depths_dict[variant.id][1]
                    variant.info['sniffles_DR'] = depths_dict[variant.id][2]
                    variant.info['sniffles_DV'] = depths_dict[variant.id][3]
                    print(variant, end = '', file=annot_vcf)
                else:
                    error_string = str(variant.id) + " has no depths\n"
                    print(error_string)
            else:
                print("Variant at " + str(variant.chrom) + ":" + \
                    str(variant.pos) + " has no ID")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument("-v", dest="vcf_file", metavar="input.vcf",
                               help="the input VCF file", required=True)
    requiredNamed.add_argument("-d", dest="depths", metavar="depths.bed",
                               help="bed file with depths from SV callers",
                               required=True)
    requiredNamed.add_argument("-o", dest="out_vcf", metavar="output.vcf",
                               help="path to output VCF file",
                               required=True)                               
    args = parser.parse_args()

    annotate_vcf(input_vcf = args.vcf_file,
                 depths = args.depths, output_vcf = args.out_vcf)       
