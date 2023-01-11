#!/usr/bin/env python
#	USAGE: python merge_prep.py
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
import collections
from functools import reduce
##########################################################################
##############                  Custom functions              ############
##########################################################################
class Naming(object):
    '''
        Split based on tumor or normal identity. Assumes id is
        <tool>_<sample_name>.
    '''
    
    def __init__(self, samples, tumor, normal):
        self.samples = samples
        self.normal = normal
        self.tumor = tumor
        self.get_lists()
    
    
    def get_lists(self):
        '''
            Make list of tumor-only and normal-only sample names
            NOTE: N/T order is reversed by bcftools merge.
        '''
        if self.get_originals(self.samples[0]) == self.normal:
            self.tumor_samples = self.samples[1::2]
            self.normal_samples = self.samples[::2]
        elif self.get_originals(self.samples[0]) == self.tumor:
            self.tumor_samples = self.samples[::2]
            self.normal_samples = self.samples[1::2]
        self.pairs = [sample_name for sample_name in zip(self.normal_samples, self.tumor_samples)]


    def get_originals(self, sample):
        '''
            get original sample names. NOTE: N/T order is
            reversed by bcftools merge if alpha order is reversed.
        '''
        return '_'.join(sample.split('_')[1:])


class Variant(Naming):
    '''
        Import a pysam record. and write out from record.
        The class allows editing of fixed elements like 
        the min number of samples in the VCF.
    '''
    
    def __init__(self, record, samples, tumor, normal):
        self.record = record
        self.line = str(record)
        self.samples = samples
        self.tumor = tumor
        self.normal = normal
        Naming.__init__(self, self.samples, self.tumor, self.normal)
        self.parts = self.line.split('\t')
        # GT fields found
        self.gt_tools = set()
        # VCF columns
        self.chrom = self.parts[0]
        self.pos = self.parts[1]
        self.id = self.parts[2]
        self.ref = self.parts[3]
        self.alts = self.parts[4]
        self.qual = self.parts[5]
        self.filters = self.parts[6]
        self.info = self.parts[7]
        self.format = ':'.join([key for key in self.get_uniq_keys()])
        self.find_tools()


    def replace_empty(self, value, key, sep=','):
        '''
            Replace an empty field with ".". In pysam this will be
            a tuple with None as the only value.
        '''
        if isinstance(value, tuple):
            if len(value) == 1 and value[0] == None:
                return '.'
            if len(value) > 1 and all(map(lambda x: x != None, value)) and key  == 'GT':
                return '/'.join(['.' if x == None else str(x) for x in value])
        elif value == None:
            return '.'
        if key  == 'GT':
            return '/'.join(['.' if x == None else str(x) for x in value])
        if (isinstance(value, collections.Iterable) and not isinstance(value, str)) or \
            isinstance(value, tuple):
            joined = []
            for i in value:
                if i == None:
                    joined.append('')
                else:
                    joined.append(str(i))
            return sep.join(joined)
        return str(value)


    def not_empty(self, value, key):
        '''
            Test if a format field is not empty. In pysam this will be
            a tuple with None as the only value.
        '''
        if isinstance(value, tuple):
            if len(value) == 1 and value[0] == None:
                return False
            if len(value) > 1 and all(map(lambda x: x == None, value)) and key == 'GT':
                return False
        if value == None:
            return False
        if value == '.':
            return False
        return True
    
    
    def deuniqify_gt(self, key, sample_name):
        '''
        Remove tool prefix from key (only if it was added for merge).
        '''
        format_keys = self.record.samples[sample_name].keys()
        if key not in format_keys and \
                key.endswith('GT') and \
                'GT' in format_keys:
            return 'GT'
        return key
                
                
    def find_tools(self):
        '''
            Return list of samples with keys.
        '''
        self.final_normal_samples = []
        self.final_tumor_samples = []
        for pair in self.pairs:
            found_keys = []
            for sample_name in pair:
                found_keys += [key for key in self.uniq_keys if self.not_empty(self.record.samples[sample_name][self.deuniqify_gt(key, sample_name)], self.deuniqify_gt(key, sample_name))]
            if len(found_keys) > 0:
                self.final_normal_samples.append(pair[0])
                self.final_tumor_samples.append(pair[1])
    
    
    def uniqify_gt(self, key, sample_name):
        '''
            Add tool prefix to key.
        '''
        if key == 'GT':
            tool = sample_name.split('_')[0]
            tool = tool.split(':')[-1]
            self.gt_tools.update(set([tool]))
            return tool + '_' + key
        return key


    def find_keys(self, record, sample_name):
        '''
            Return list of keys with values for FORMAT
        '''
        format_keys = record.samples[sample_name].keys()
        found_keys = [self.uniqify_gt(key, sample_name) for key in format_keys if self.not_empty(record.samples[sample_name][key], key)]
        return found_keys
    


    def get_uniq_keys(self):
        '''
            Return a key/value pairs for any key with a value for on tool's results.
        '''
        seen = set()
        seen_add = seen.add
        found_keys = reduce(list.__add__,
                            [self.find_keys(self.record, sample_name) for sample_name in self.samples],
                            [])
        self.uniq_keys = [key for key in found_keys if not (key in seen or seen_add(key))]
        return self.uniq_keys

    
    def reduce_samples(self, samples):
        '''
            Find all keys with a value for any sample. Return '.' for missing values.
        '''
        uniq_keys = self.get_uniq_keys()
        sample_result = []
        for sample_name in samples:
            tool = sample_name.split('_')[0]
            tool = tool.split(':')[-1]
            sample_result.append(':'.join([self.replace_empty(self.record.samples[sample_name][self.deuniqify_gt(key, sample_name)], self.deuniqify_gt(key, sample_name)) for key in uniq_keys if tool == key.split('_')[0]]))
        return ':'.join(sample_result)

    
    def write(self):
        '''
            Return a reformatted string from a pysam object.
        '''
        line = [self.chrom,
                self.pos,
                self.id,
                self.ref,
                self.alts,
                self.qual,
                self.filters,
                self.info,
                self.format]
        line += [self.reduce_samples(self.final_normal_samples)]
        line += [self.reduce_samples(self.final_tumor_samples)]
        self.new_line = '\t'.join(line) + '\n'
        return self.new_line


def modify_header(bcf_in, tool):
    '''
        Add new FORMAT field
    '''
    bcf_in.header.formats.add(id=tool + '_GT', number='1',
                           type='String',
                           description='Genotype from ' + tool)
    return bcf_in


def load_header(bcf_in):
    '''
        Load a VCF file header as a list of lines.
    '''
    header = '\n'.join(str(bcf_in.header).split('\n')[:-2]) + '\n'
    return header


def read_vcf(vcf_file):
    '''
        Read in annotated VCF file.
        '''
    bcf_in = pysam.VariantFile(vcf_file)  # auto-detect input format
    return bcf_in


def vcf_writer(bcf_in, vcf_out_file, tumor, normal):
    '''
        Write out the VCF file with corrected information
    '''
    with open(vcf_out_file, 'w') as vcf_out:
        samples = [sample_name for sample_name in bcf_in.header.samples]
        names = Naming(samples, tumor, normal)
        lines = []
        gt_tools = set()
        for record in bcf_in.fetch():
            out = Variant(record, samples, tumor, normal)
            gt_tools.update(out.gt_tools)
            lines.append(out.write())
        #  ==========================
        #  Add GT
        #  ==========================
        for tool in gt_tools:
            bcf_in = modify_header(bcf_in, tool)
        header = load_header(bcf_in)
        #  ==========================
        #  Write header
        #  ==========================
        for line in header:
            vcf_out.write(line)
        vcf_out.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                                 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                                 names.normal, names.tumor]) + '\n')
        #  ==========================
        #  Write variants
        #  ==========================
        for line in lines:
            vcf_out.write(line)


def main():
    '''
        Merge the VCF columns by:
        1) Getting all INFO fields with values
        2) Getting all FORMAT fields from any sample with non-empty values
        3) Uniqify GT now that bcftools merge is done
    '''
    #  ==========================
    #  Input variables
    #  ==========================
    vcf_in = sys.argv[1]
    vcf_out_file = sys.argv[2]
    tumor = sys.argv[3]
    normal = sys.argv[4]
    
    assert os.path.isfile(vcf_in), 'Failed to find caller VCF call file :' + vcf_in
    #  ==========================
    #  Run prep
    #  ==========================
    bcf_in = read_vcf(vcf_in)
    vcf_writer(bcf_in, vcf_out_file, tumor, normal)



##########################################################################
#####       Execute main unless script is simply imported     ############
#####                for individual functions                 ############
##########################################################################


if __name__ == '__main__':
    main()