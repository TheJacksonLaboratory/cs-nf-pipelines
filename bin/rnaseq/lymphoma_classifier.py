#! /usr/bin/env python

"""
Based on expected expression of a small set of genes, determine whether this
supposed tumor sample has been converted into a lymphoma.

We use two cut-offs, one for passaged tumors, the other for patient.
"""
from __future__ import print_function
import sys
import argparse
import math
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cutoff', default=3, type=int,
                        help="Cut-off for patient samples")
    parser.add_argument('-o', '--output', default='lymphoma_score.txt',
                        help="Output file into which the total z-score is "
                             "written")
    parser.add_argument('-c', '--counts', 
                        help="A file containing RSEM counts")
    parser.add_argument('--expected_expression',
                        help="file containing gene name in column 0, "
                             "expected (average) expression in column 3, "
                             "and standard deviation in column 4")
    return parser.parse_args()

def get_expected(fn):
    d = {}
    for line in open(fn):
        # parts[0]: gene name
        # parts[1]: expected up/down regulation
        # parts[3]: average expression
        # parts[4]: stddev
        parts = [x.strip() for x in line.split()]
        d[parts[0]] = {'updown': float(parts[1]),
                       'average': float(parts[3]),
                       'stddev': float(parts[4])}
    return d

def z_score(count, updown, average, stddev):
    arbitrary_scale_factor = 25.0
    l2 = math.log(count + 1.0, 2)
    z = updown * (l2 - average) / (arbitrary_scale_factor * stddev)
    return z

def normalize_counts_z_total(counts, expected):
    c = {}
    non_zero = []
    z_total = 0.0
    with open(counts) as f:
        next(f)
        for line in f:
            parts = [x.strip() for x in line.split()]
            try: 
                ensembl = parts[0].split("_")[0].split('.')[0]
                gene_id = parts[0].split("_")[1]
            except:
                try:
                    ensembl = parts[0].split("_")[0]
                    gene_id = parts[0].split("_")[1]
                except:
                    ensembl = parts[0]
                    gene_id = 'not_provided'
            c[ensembl] = {'raw_count': float(parts[4]),
                          'full_ensemblID': parts[0],
                          'transcripts': parts[1],
                          'gene_name': gene_id}
            if float(parts[4]) > 0:
                non_zero.append(float(parts[4]))
        ## parse the count file from RSEM into a nested dict. 

        upper_quantile = np.quantile(non_zero, [0.75])[0]
        ## calc upper 75th quantile from non-zero expected counts

        zero_norm = 0
        norm_count = 0

        for ensembl in c.keys():
            c[ensembl]['normalized_count'] = ((c[ensembl]['raw_count'] / upper_quantile) * 1000)
            # normalize all expected counts

            if ensembl in expected:
                if c[ensembl]['normalized_count'] != 0:
                    zero_norm += 1
                else: 
                    norm_count += 1
                # edge case: Return no z total where 0 genes are quantified in the expected set. 

                e = expected[ensembl]

                z_total += z_score(c[ensembl]['normalized_count'], e['updown'], e['average'], e['stddev'])
                # match to genes used in classifier, and classify 

        return([z_total, zero_norm, norm_count])

def main():
    print("Starting lymphoma_classifier")
    args = parse_args()

    expected = get_expected(args.expected_expression)

    z_calc = normalize_counts_z_total(args.counts, expected)

    z_total = z_calc[0]
    zero_norm = z_calc[1]
    norm_count = z_calc[1]

    with open(args.output, 'w') as f:

        if zero_norm == 0:
            f.write('Sample contained no quantified expression in EBV classifier related genes. Classifier can not make a determination.\n')
        # edge case: Return no z total where 0 genes are quantified in the expected set. 

        elif z_total > args.cutoff:
            f.write(str(round(z_total, 3)) + '\tSample is possibly EBV-transformed or is a lymphoma model (z-score > cutoff threshold: ' + str(args.cutoff) + ').\n')
        # over the threshold

        else:
            f.write(str(round(z_total, 3)) + '\tSample is below EBV classifier z-score cutoff threshold: ' + str(args.cutoff) + ').\n')
        # under the threshold

    print("Finished lymphoma_classifier")

if __name__ == '__main__':
    main()
