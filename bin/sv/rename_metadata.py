#!/usr/bin/env python
# USAGE: rename_metadata.py VCF_IN VCF_OUT PREFIX
# DESCRIPTION: Takes in a VCF and a tool name
#              and preps the file by:
#         1) add tool + "_" to all INFO and FORMAT def lines (unless the key is "GT")

import sys
import shutil
import logging as log
from collections import OrderedDict
import re
import os


class Variant(object):


    def __init__(self, line, prefix=''):
        self.line = line
        self.prefix = prefix
        self.parts = line.split('\t')
        # VCF columns
        self.chrom = self.parts[0]
        self.pos = self.parts[1]
        self.id = self.parts[2]
        self.ref = self.parts[3]
        self.alts = self.parts[4].split(',')
        self.qual = self.parts[5]
        self.filters = self.parts[6].split(';')
        self.info = self.parts[7].split(';')
        self.format = self.parts[8].split(':')
        self.samples = self.parts[9:]
        # modify
        self.info_dict = self.get_info()
        self.fix_format()


    def get_info(self):
        '''
            Get current info line and add prefix if needed
        '''
        info_dict = OrderedDict()
        for item in self.info:
            if '=' in item:
                info_dict.update({self.prefix + item.split('=')[0] : item.split('=')[1]})
            else:
                info_dict.update({self.prefix + item : None})
        return info_dict
    
    
    def add_prefix(self, key):
        '''
        Add prefix unless variant is GT.
        '''
        if key != 'GT':
            return self.prefix + key
        return key
    
    
    def fix_format(self):
        '''
            Add any prefix to format entries. Skips GT becuase 
            in at least one arrangement (the GT being first without the
            full name GT bcftools combines other GT feilds (e.g. Mutect2's
            PGT into the first position and writes bad VCF files.
                This:
            "0/0:37,1:0.098:10,0:27,1:12:224,349:25:2:0|1:16759500_G_C:.:."
                Becomes (with non-Non-ASCII characters changed to *):
            "0/00/1!%I9**=m**=!
            :37,1:0.098:10,0:27,1:12:224,349:25:2:0|1:16759500_G_C:.:."
        '''
        new_format = [self.add_prefix(key) for key in self.format]
        self.format = new_format
    

    def write(self):
        line = [self.chrom,
                self.pos,
                self.id,
                self.ref,
                ','.join(self.alts),
                str(self.qual),
                ';'.join(self.filters),
                ';'.join(['='.join([x for x in [key, self.info_dict[key]] if x != None]) for key in self.info_dict]),
                ':'.join(self.format)]
        line += self.samples
        self.new_line = '\t'.join(line)
        return self.new_line


def fix_header(line, prefix):
    '''
        Add prefix as needed. Replace AD with standard AD line because
        GATK also makes this replacement and it conflicts with Lancet similar but unique
        wording.
    '''
    if '##FORMAT=' in line or \
            '##INFO=' in line:
        id = re.search('ID=(?P<id>[^>,]+)', line)
        if id == None:
            log.error('FORMAT or INFO line missing ID field: ' + line)
            sys.exit()
        if id.group(1) == 'AD' and '##FORMAT=' in line:
            line = '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n'
        if id.group(1) != 'GT':
            line = line.replace('ID=' + id.group(1),
                                'ID=' + prefix + id.group(1))
    return line


def load_header(vcf_in, prefix):
    '''
        Load a VCF file header as a list of lines.
    '''
    with open(vcf_in) as vcf:
        header = [fix_header(line, prefix) for line in vcf if line.startswith('#')]
    return header


def load_vcf(vcf_in, prefix):
    '''
        Load a VCF file and fixes lines.
    '''
    with open(vcf_in) as vcf:
        for line in vcf:
            if not line.startswith('#'):
                yield Variant(line, prefix).write()




def rename_metadata(vcf_file, vcf_out_file, prefix):
    '''
        Add prefix to FORMAT and INFO keys
    '''
    #  =====================
    #  test if renaming should occur
    #  =====================
    rename = False
    if vcf_out_file == vcf_file:
        vcf_out_file = vcf_file + '_tmp.vcf'
        rename = True
    #  =====================
    #  rename
    #  =====================
    header = load_header(vcf_file, prefix)
    vcf_reader = load_vcf(vcf_file, prefix)
    #  =====================
    #  rewrite
    #  =====================
    vcf_writer(header, vcf_reader, vcf_out_file)
    #  =====================
    #  rename output VCF
    #  =====================
    if rename:
        shutil.move(vcf_out_file, vcf_file)
    return True


def correct_ad_line():
    '''
    '''

def vcf_writer(header, vcf_reader, vcf_out_file):
    '''
       Write out a VCF file with the a prefix added to INFO and FORMAT keys
    '''
    with open(vcf_out_file, 'w') as vcf_out:
        for line in header:
            vcf_out.write(line)
        for line in vcf_reader:
            vcf_out.write(line)
    return True


def main():
    '''
    DESCRIPTION: Takes in a VCF and a tool name
    and preps the file by:
        1) add tool + "_" to all INFO and FORMAT def lines (unless the key is "GT")
    '''
    vcf_file = sys.argv[1]
    vcf_out_file = sys.argv[2]
    prefix = sys.argv[3] + '_'
    assert os.path.isfile(vcf_file), 'Failed to find prep caller VCF call file :' + vcf_file
    rename_metadata(vcf_file, vcf_out_file, prefix)


#  =====================
#  Main
#  =====================


if __name__ == "__main__":
    main()