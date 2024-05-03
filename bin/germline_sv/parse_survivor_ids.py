#!/usr/bin/env python
"""
Parse SURVIVOR merged VCFs to extract IDs

note: "import vcf" below refers to the module pyvcf
"""


import argparse
import vcf
import csv


def parse_vcf(vcf_file, out_file):
    """Iterate through VCF records and write summary information to a CSV

    Parameters
    ----------
    vcf_file: str
        The path to the SURVIVOR VCF file

    out_file: str
        The path to a CSV file to output summary information

    Notes
    -----
    This function writes out a CSV with the following fields
        0: the SV ID
        1: number of reference reads
        2: number of variant reads
    """

    with open(out_file, 'w', newline='\n') as ofile:
        writer = csv.writer(ofile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['SURVIVOR_ID', 'nanosv_id', 'sniffles_id'])
        vcf_reader = vcf.Reader(filename=vcf_file)
        for record in vcf_reader:
            writer.writerow([record.ID, 
                             record.samples[0]['ID'],
                             record.samples[1]['ID']])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument("-v", dest="vcf_file", metavar="input.vcf",
                               help="the input VCF file", required=True)
    requiredNamed.add_argument("-o", dest='out_file', metavar="output.csv",
                               help="summary CSV output",
                               required=True)
    args = parser.parse_args()

    parse_vcf(args.vcf_file, args.out_file)