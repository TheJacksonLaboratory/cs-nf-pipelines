#!/usr/bin/env python

"""Annotate VCF file

Apply 'OnTarget' to original VCF file.
Based on matching annotations from VCF and annotated bedpe.

"How to add custom fields to the VCF files?" idea from: https://www.biostars.org/p/209188/
"""

__author__ = str(('Brian Sanderson', 'Mike Lloyd', 'Benjamin Leopold'))

import os
import argparse
import pysam

def validate_targets(input_vcf, chr, start, end):
    """Validate the type and format of the target regions
    """
    myvcf = pysam.VariantFile(input_vcf,'r')
    variant = next(myvcf)

    if type(chr) == type(variant.chrom) and type(variant.pos) == type(start) and type(variant.pos) == type(end):
        if ('chr' in variant.chrom and 'chr' in chr) or ('chr' not in variant.chrom and 'chr' not in chr):
            myvcf.close()
            return(0)
        else:
            myvcf.close()
            raise ValueError("Format of target chromosome " + chr + " does not match VCF " + variant.chrom)
    else:
        raise ValueError("Format of target region doesn't match VCF")


def annotate_vcf(input_vcf, chr, start, end, output_vcf) :
    """Iterate through VCF records and append a new INFO field OnTarget with
    a boolean value that represents whether the variant overlaps with
    the targeted region for adaptive sequencing
    """

    myvcf = pysam.VariantFile(input_vcf,'r')
    myvcf.header.info.add('OnTarget','1','String',
                          'Is SV call contained within a region targeted by sequencing: TRUE, FALSE')
    with open(output_vcf, 'w') as annot_vcf:
        print(myvcf.header, end='', file=annot_vcf)
        for variant in myvcf:
            if variant.chrom == chr and variant.pos in range(start, end):
                variant.info['OnTarget'] = 'TRUE'
                print(variant, end = '', file=annot_vcf)
            else:
                variant.info['OnTarget'] = 'FALSE'
                print(variant, end = '', file=annot_vcf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument("-v", dest="vcf_file", metavar="input.vcf",
                               help="the input VCF file", required=True)
    requiredNamed.add_argument("-c", dest="chr", metavar="chrN",
                               help="target chromosome",
                               required=True)
    requiredNamed.add_argument("-s", dest="start", metavar="NNNNNNN",
                               help="start coordinate of target region",
                               required=True)
    requiredNamed.add_argument("-e", dest="end", metavar="NNNNNNN",
                               help="end coordinate of target region",
                               required=True)                                                              
    requiredNamed.add_argument("-o", dest="out_vcf", metavar="output.vcf",
                               help="path to output VCF file",
                               required=True)                               
    args = parser.parse_args()

    validate_targets(input_vcf = args.vcf_file,
                     chr = str(args.chr), start = int(args.start),
                     end = int(args.end))

    annotate_vcf(input_vcf = args.vcf_file,
                     chr = str(args.chr), start = int(args.start),
                     end = int(args.end), output_vcf = args.out_vcf)