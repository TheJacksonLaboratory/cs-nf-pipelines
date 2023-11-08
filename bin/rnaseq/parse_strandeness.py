import glob
import os
import argparse
from collections import defaultdict
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--analysis_dir',
                        help="Directory containing RNAseq output from nextflow workflow")
    parser.add_argument('-o', '--output', default='aggregated_strandeness_stats.txt',
                        help="aggregated statsistics output file")
    return parser.parse_args()

def main():
    args = parse_args()
 
    print('Working in: ' + args.analysis_dir+'/*/*/*strandeness.txt')


    stats_dict = defaultdict(dict)

    for file in glob.glob(args.analysis_dir+'/*/stats/*strandedness.txt'):
        sampleID = os.path.basename(file).replace('_strandedness.txt','')
        print('Found strandedness file for: ' + sampleID)
        f = open(file, "r")
        for x in f:
            if ('Fraction of reads failed to determine:' in x):
                stats_dict[sampleID]['fract_failed_det'] = x.rstrip().split(": ",1)[1]
            if ('Fraction of reads explained by "1++,1--,2+-,2-+":' in x):
                stats_dict[sampleID]['forward_read_fract'] = x.rstrip().split(": ",1)[1].split(" ",1)[0]
            if ('Fraction of reads explained by "1+-,1-+,2++,2--"' in x):
                stats_dict[sampleID]['reverse_read_fract'] = x.rstrip().split(": ",1)[1].split(" ",1)[0]
            if ('Data does not fall into a likely stranded (max percent explained > 0.9) or unstranded layout (max percent explained < 0.6)' in x):
                stats_dict[sampleID]['determination'] = 'low_qual'
            if ('Data is likely RF' in x):
                stats_dict[sampleID]['determination'] = 'RF'
            if ('Data is likely FR' in x):
                stats_dict[sampleID]['determination'] = 'FR'
            if ('Data is likely unstranded' in x):
                stats_dict[sampleID]['determination'] = 'unstranded'

    df = pd.DataFrame.from_dict(stats_dict).transpose()
    df.to_csv(args.output, sep='\t', index=True, index_label='Sample')

if __name__ == '__main__':
    main()
