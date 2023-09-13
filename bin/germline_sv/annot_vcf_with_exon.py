#!/usr/bin/env python

"""Annotate VCF file

Apply 'InExon' to original VCF file.
Based on matching annotations from VCF and annotated bedpe.

"How to add custom fields to the VCF files?" idea from: https://www.biostars.org/p/209188/
"""

__author__ = str(('Brian Sanderson', 'Mike Lloyd', 'Benjamin Leopold'))

import os
import argparse
import pysam

def parse_bed(exon_dict, bed_handle):
    """Parse BED files containing intersections between structural variants
    and annotated exons, and return a dictionary of boolean values for each
    named SV
    """
    if os.path.exists(bed_handle) and os.stat(bed_handle).st_size == 0:
        print("No records in " + str(bed_handle))
        return exon_dict

    elif os.path.exists(bed_handle):
        with open(bed_handle, 'r') as ins_f:
            for line_handle in ins_f:
                line = line_handle.rstrip().split('\t')
                sv_name = line[3]
                if sv_name in exon_dict:
                    error_string = str(sv_name) + " on Chr" + str(line[0]) + \
                        " at " + str(line[1]) + " has multiple SV records"
                    #print(error_string)
                else:
                    exon_dict[sv_name] = "TRUE"
        return exon_dict
    else:
        raise FileExistsError


def annotate_vcf(input_vcf, ins, dele, inv, dup, tra, output_vcf) :
    """Iterate through VCF records and append a new INFO field InExon with
    a boolean value that represents whether the variant overlaps with
    annotated exons
    """

    exon_dict = {}
    
    exon_dict = parse_bed(exon_dict, ins)
    exon_dict = parse_bed(exon_dict, dele)
    exon_dict = parse_bed(exon_dict, inv)
    exon_dict = parse_bed(exon_dict, dup)
    exon_dict = parse_bed(exon_dict, tra)

    myvcf = pysam.VariantFile(input_vcf,'r')
    myvcf.header.info.add('InEXON','1','String',
                        'Is SV call contained within an exonic region or regions: TRUE, FALSE')

    with open(output_vcf, 'w') as annot_vcf:
        
        print(myvcf.header, end='', file=annot_vcf)

        for variant in myvcf:
            if variant.id:
                if variant.id in exon_dict:
                    variant.info['InEXON'] = 'TRUE'
                    print(variant, end = '', file=annot_vcf)
                else:
                    variant.info['InEXON'] = 'FALSE'
                    print(variant, end = '', file=annot_vcf)
            else:
                print("Variant at " + str(variant.chrom) + ":" + \
                    str(variant.pos) + " has no ID")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument("-v", dest="vcf_file", metavar="input.vcf",
                               help="the input VCF file", required=True)
    requiredNamed.add_argument("-i", dest="ins", metavar="ins.exons.bed",
                               help="insertion exon intersections",
                               required=True)
    requiredNamed.add_argument("-d", dest="dele", metavar="del.exons.bed",
                               help="deletion exon intersections",
                               required=True)    
    requiredNamed.add_argument("-n", dest="inv", metavar="inv.exons.bed",
                               help="inversion exon intersections",
                                required=True)
    requiredNamed.add_argument("-t", dest="tra", metavar="tra.exons.bed",
                               help="translocation exon intersections",
                               required=True)
    requiredNamed.add_argument("-u", dest="dup", metavar="dup.exons.bed",
                               help="duplications exon intersections",
                               required=True)
    requiredNamed.add_argument("-o", dest="out_vcf", metavar="output.vcf",
                               help="path to output VCF file",
                               required=True)                               
    args = parser.parse_args()

    annotate_vcf(input_vcf = args.vcf_file,
                 ins = args.ins, dele = args.dele, inv = args.inv,
                 tra = args.tra, dup = args.dup, output_vcf = args.out_vcf)       
