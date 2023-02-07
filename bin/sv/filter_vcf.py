import pysam
import numpy as np
import argparse
import logging as log
import re

chrom_pattern = re.compile('[\[\]](.*):')

class FilterNonChroms():
    def __init__(self, vcf_file,
                 out_file,
                 chroms):
        '''
        Requires
        NYGC column headers:
        #chr    start    end    type    log2    tool    pair_id    info    focal    cytoband
        Without nygc columns the step returns an empty table
        '''
        self.chroms = chroms
        self.bcf_in = self.read_vcf(vcf_file)
        self.bcf_out = self.start_vcf(out_file)
        self.filter_vcf()
        
    def read_vcf(self, vcf_file):
        bcf_in = pysam.VariantFile(vcf_file)
        return bcf_in
    
    def start_vcf(self, out_file):
        bcf_out = pysam.VariantFile(out_file, 'w', header=self.bcf_in.header)
        return bcf_out
    
    def filter_non_chroms(self, record):
        '''
            Filter calls to leave calls that are only on chroms.
        '''
        if record.contig in self.chroms:
            for alt in record.alts:
                result = re.search(chrom_pattern, alt)
                if result and result[1] not in self.chroms:
                    return True
        else:
            return True
        return False
        
    def filter_vcf(self):
        '''
            Only print calls that have:
            1) ref in chroms
            2) alt in chroms
        '''
        for record in self.bcf_in.fetch():
            if not self.filter_non_chroms(record):
                exit_status = self.bcf_out.write(record)
                if exit_status != 0:
                    print(exit_status)


def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf-file',
                        help='SV VCF file',
                        required=True
                       )
    parser.add_argument('--output',
                        help='Output VCF file',
                        required=True
                       )
    parser.add_argument('--chroms',
                        help='A space separated list of chroms to plot.',
                        required=False,
                        nargs='*',
                        default=False
                       )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__


def main():
    args = get_args()
    FilterNonChroms(args['vcf_file'],
        args['output'],
        chroms=args['chroms'])


if __name__ == "__main__":
    main()