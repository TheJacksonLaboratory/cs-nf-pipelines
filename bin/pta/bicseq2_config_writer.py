import pandas as pd
import os
import argparse
from os import listdir
from os.path import isfile, join
import glob

class Bicseq2Prep():
    def __init__(self, sample_id,
                 fa_files,
                 out_file,
                 mappability_directory,
                 norm_bicseq2_config,
                 temp_seqs):
        self.out_file = out_file
        self.sample_id = sample_id
        self.norm_bicseq2_config = norm_bicseq2_config
        self.mappability_directory = mappability_directory
        self.fa_files = fa_files
        self.temp_seqs = temp_seqs
        self.write_sample_configs()
            
    def match_fa_file(self, row):
        for fa_file in self.fa_files:
            if fa_file.endswith('_' + str(row.chrom_name) + '.fasta'):
                return fa_file
            
    def match_mappability_file(self, row):
        mappability_files = (glob.glob(self.mappability_directory+"/*"))
        # from the directory provided, find all files. 
        for mappability_file in mappability_files:
            if os.path.basename(mappability_file) == str(row.chrom_name) + '.uniq.txt':
                return mappability_file
            
    def match_seq_file(self, row):
        for temp_seq in self.temp_seqs:
            if os.path.basename(temp_seq) == str(self.sample_id) + '_' + str(row.chrom_name) + '.seq':
                return temp_seq
            
    def match_norm_file(self, row):
        for temp_seq in self.temp_seqs:
            if os.path.splitext(os.path.basename(temp_seq))[0] + '.norm.bin.txt' == str(self.sample_id) + '_' + str(row.chrom_name) + '.norm.bin.txt':
                return os.path.splitext(os.path.basename(temp_seq))[0] + '.norm.bin.txt'
        # Modified to programtically set this output. 

    def prep(self):
        '''initial file should start with one column named chrom_name '''
        data = pd.read_csv(self.norm_bicseq2_config, sep='\t')
        assert ''.join(data.columns.tolist()) == 'chrom_name', 'Error: initial config file should start with one column named chrom_name'
        data['fa_file'] = data.apply(lambda row: self.match_fa_file(row), axis=1)
        data['mappability'] = data.apply(lambda row: self.match_mappability_file(row), axis=1)
        data['readPosFile'] = data.apply(lambda row: self.match_seq_file(row), axis=1)
        data['bin_file_normalized'] = data.apply(lambda row: self.match_norm_file(row), axis=1)
        return data
        
    def write_sample_configs(self):
        '''write configs for the normalization step'''
        # tumor
        data = self.prep()
        config = data.to_csv(self.out_file, sep='\t', index=False)

    
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--fa-files',
                        help='List of chrom fasta files. ',
                        required=True,
                        nargs='*'
                        )
    parser.add_argument('--mappability-directory',
                        help='Directory containing mappability files. ',
                        required=True
                        )
    parser.add_argument('--temp-seqs',
                        help='List of file paths ${sample_id}_${chr}.seq '
                        '(readPosFile files output from samtools getUnique)  ',
                        required=True,
                        nargs='*'
                        )
    parser.add_argument('--norm-bicseq2-config',
                        help='Pre filled file for ${sample_id}.bicseq2.config. '
                        'Fasta-specific but sample-independent portion of config file.',
                        required=True
                        )
    parser.add_argument('--out-file',
                        help='Output config filename',
                        required=True
                        )
    parser.add_argument('--sample-id',
                        help='sample id',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__       
    
def main():
    args = get_args()
    bicseq = Bicseq2Prep(sample_id=args['sample_id'],
                         fa_files=args['fa_files'],
                         out_file=args['out_file'],
                         mappability_directory=args['mappability_directory'],
                         norm_bicseq2_config=args['norm_bicseq2_config'],
                         temp_seqs=args['temp_seqs'])
    
    
if __name__ == "__main__":
    main()