import pandas as pd
import os
import argparse


class Bicseq2Prep():
    def __init__(self, pair_id,
                 out_file,
                 seg_bicseq2_config,
                 tumor_norms,
                 normal_norms):
        self.out_file = out_file
        self.pair_id = pair_id
        self.seg_bicseq2_config = seg_bicseq2_config
        self.tumor_norms = tumor_norms
        self.normal_norms = normal_norms
        self.write_sample_configs()
            
    def match_file(self, row, files):
        for temp_norm in files:
            if os.path.basename(temp_norm).endswith('_' + str(row.chr) + '.norm.bin.txt'):
                return temp_norm
            
    def prep_pair(self):
        ''' file: tumor--normal.bicseq2.seg.config
        prep fasta-specific but sample independent portion of config file'''
        data = pd.read_csv(self.seg_bicseq2_config, sep='\t')
        # assert ''.join(data.columns.tolist()) == 'chr', 'Error: initial config file should start with one column named chr'
        assert ''.join(data.columns.tolist()) == 'chr', 'Error: initial config file should start with one column named chrom_name'

        data['case'] = data.apply(lambda row: self.match_file(row, 
                                                              files=self.tumor_norms), axis=1)
        data['control'] = data.apply(lambda row: self.match_file(row, 
                                                              files=self.normal_norms), axis=1)
        return data
            
        
    def write_sample_configs(self):
        '''create and upload configs for the normalization step'''
        # tumor
        data = self.prep_pair()
        config = data.to_csv(self.out_file, sep='\t', index=False)

    
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--tumor-norms',
                        help='List of file paths ${tumor}_${chr}.norm.bin.txt '
                        ' (Output from Bicseq2Norm ).',
                        required=True,
                        nargs='*'
                        )
    parser.add_argument('--normal-norms',
                        help='List of file paths ${normal}_${chr}.norm.bin.txt '
                        ' (Output from Bicseq2Norm ).',
                        required=True,
                        nargs='*'
                        )
    parser.add_argument('--seg-bicseq2-config',
                        help='Pre filled file for ${pair_id}.bicseq2.seg.config. '
                        'Fasta-specific but sample-independent portion of config file.',
                        required=True
                        )
    parser.add_argument('--out-file',
                        help='Output config filename',
                        required=True
                        )
    parser.add_argument('--pair-id',
                        help='pair id',
                        required=True
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__       
    
    
def main():
    args = get_args()
    bicseq = Bicseq2Prep(pair_id=args['pair_id'],
                         out_file=args['out_file'],
                         seg_bicseq2_config=args['seg_bicseq2_config'],
                         tumor_norms=args['tumor_norms'],
                         normal_norms=args['normal_norms'])
    
    
if __name__ == "__main__":
    main()