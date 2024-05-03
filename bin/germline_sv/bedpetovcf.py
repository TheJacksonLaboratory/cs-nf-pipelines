from fuc import pyvcf
import sys
import pandas as pd
import numpy as np 

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--input", dest="input",
                  help="Input BedPE File", metavar="BEDPE")
parser.add_option("-v", "--vcf", dest="vcf",
                  help="Output VCF File", metavar="VCF_FILE")

(options, args) = parser.parse_args()

df = pd.read_csv(options.input, sep='\t')
df.columns = df.columns.str.replace("#", "")
header = ['CHROM', 'POS', 'skip_chr1_end', 'CHR2', 'END', 'skip_chr2_end', 'SVTYPE', 'score', 'strand1', 'strand2', 'evidence', 'tools', 'SUPP', 'SUPP_VEC', 'SampleID']
df.columns = header[:len(df.columns)]

df['SVLEN'] = df['END'] - df['POS']
df.SVLEN =  np.where((df.SVTYPE == 'TRA') | (df.SVTYPE == 'INS'), 0, df.SVLEN)

df['INFO'] = 'SUPP=' + df['SUPP'].astype(str) + ';SUPP_VEC=' + df['SUPP_VEC'].astype(str) + ';SVTYPE=' + df['SVTYPE'] + ';CHR2=' + df['CHR2'].astype(str) + ';END=' + df['END'].astype(str) + ';SVLEN=' + df['SVLEN'].astype(str) + ';STRANDS=' + df['strand1'] + df['strand2'] + ';EVIDENCE=' + df['evidence']

df['REF'] = 'N'
df['ALT'] = '<' + df['SVTYPE'] + '>'
df['FILTER'] = 'PASS'
df['FORMAT'] = 'GT'
df[df['SampleID'][0]] = './.'
df['QUAL'] = '.'

df['ID'] = df.groupby('SVTYPE')['SVTYPE'].rank(method='first').astype(int)
df['ID'] = df['SVTYPE'] + df['ID'].astype(str).str.zfill(5)

vcf_convert = df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', df['SampleID'][0]]]

header = [
'##fileformat=VCFv4.1',
'##ALT=<ID=DEL,Description="Deletion">',
'##ALT=<ID=DUP,Description="Duplication">',
'##ALT=<ID=INV,Description="Inversion">',
'##ALT=<ID=TRA,Description="Translocation">',
'##ALT=<ID=INS,Description="Insertion">',
'##FILTER=<ID=LowQual,Description="Not used in this context">',
'##INFO=<ID=SUPP,Number=1,Type=String,Description="Number of samples supporting the variant">',
'##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of supporting samples.">',
'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
'##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">',
'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">',
'##INFO=<ID=STRANDS,Number=1,Type=String,Description="Indicating the direction of the reads with respect to the type and breakpoint.">',
'##INFO=<ID=EVIDENCE,Number=1,Type=String,Description="Evidence supporting SV call">',
'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
]

vf = pyvcf.VcfFrame.from_dict(header, vcf_convert)

vf.to_file(options.vcf)