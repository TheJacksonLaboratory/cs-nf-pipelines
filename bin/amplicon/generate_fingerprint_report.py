import argparse
import vcf

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--input_file', dest='input_file', metavar='VCF', help='VCF file', required=True)
parser.add_argument('--rsid_file', dest='rsid_file', metavar='BEDFILE', help='tab delim vcf file with rsID and target information', required=True)
parser.add_argument('--output_prefix', dest='output_prefix', metavar='OUTPUT_PREFIX', help='an output name prefix', required=True)

args = parser.parse_args()

d = {}
x = []
y = []
all_rsid = []
seen_rsid = []

from re import findall

def sort_ids(s):
    f = findall(r'\d+|[A-Za-z_]+', s.lower())
    return list(map(lambda x:int(x) if x.isdigit() else x, f))

with open(args.rsid_file) as fin:
    rows = ( line.strip().split('\t') for line in fin )
    for row in rows:
        
        all_rsid.append(row[3])

        if (row[0] == 'chrX'):
            x.append([row[1], row[2], row[5], row[3]])
            d[row[3]] = {'chrom': row[0], 'start': row[1], 'end': row[2], 'IDT_id': row[4], 'gene_with_target': row[5], 'full_IDT_name': row[6]}
        elif (row[0] == 'chrY'):
            y.append([row[1], row[2], row[5], row[3]])
            d[row[3]] = {'chrom': row[0], 'start': row[1], 'end': row[2], 'IDT_id': row[4], 'gene_with_target': row[5], 'full_IDT_name': row[6]}
        else:
            d[row[3]] = {'chrom': row[0], 'start': row[1], 'end': row[2], 'IDT_id': row[4], 'gene_with_target': row[5], 'full_IDT_name': row[6]}

vcf_reader = vcf.Reader(open(args.input_file, 'r'))

with open(args.output_prefix + '.on_target_SNPs.tsv', 'w') as on_target, open(args.output_prefix + '.off_target_SNPs.tsv', 'w') as off_target:

    header = ['chrom', 'pos', 'ref', 'alt', 'rsID', 'AD_ref', 'AD_alt', 'gene_with_target']
    on_target.write('\t'.join([str(i) for i in header]) + '\n')
    
    header = ['chrom', 'pos', 'ref', 'alt', 'rsID', 'AD_ref', 'AD_alt']
    off_target.write('\t'.join([str(i) for i in header]) + '\n')

    for record in vcf_reader:

        if record.ID in d.keys():
            print_out = [record.CHROM, record.POS, record.REF, ','.join([str(i) for i in record.ALT]), record.ID, record.samples[0]['AD'][0], record.samples[0]['AD'][1], d[record.ID]['gene_with_target']]
            on_target.write('\t'.join([str(i) for i in print_out]) + '\n')
            seen_rsid.append(record.ID)

        elif (str(record.CHROM) == 'chrX'):
            x_target = 0
            for region in x:
                if int(region[0]) <= record.POS <= int(region[1]):
                    print_out = [record.CHROM, record.POS, record.REF, ','.join([str(i) for i in record.ALT]), record.ID, record.samples[0]['AD'][0], ','.join([str(i) for i in record.samples[0]['AD'][1:]]), region[2]]
                    # note, the AD here can be multi-allelic as they are regional matches, not rsid matches. For the standard fingerprint SNPs, rsID will not be multi-allelic. 
                    x_target = 1
                    seen_rsid.append(region[3])
            if (x_target == 1):
                on_target.write('\t'.join([str(i) for i in print_out]) + '\n')
            else:
                print_out = [record.CHROM, record.POS, record.REF, ','.join([str(i) for i in record.ALT]), record.ID, record.samples[0]['AD'][0], ','.join([str(i) for i in record.samples[0]['AD'][1:]])]
                # note, the AD here can be multi-allelic as they are regional matches, not rsid matches. For the standard fingerprint SNPs, rsID will not be multi-allelic. 
                off_target.write('\t'.join([str(i) for i in print_out]) + '\n')

        elif (str(record.CHROM) == 'chrY'):
            y_target = 0
            for region in y:
                if int(region[0]) <= record.POS <= int(region[1]):
                    print('in region y')
                    print_out = [record.CHROM, record.POS, record.REF, ','.join([str(i) for i in record.ALT]), record.ID, record.samples[0]['AD'][0], ','.join([str(i) for i in record.samples[0]['AD'][1:]]), region[2]]
                    # note, the AD here can be multi-allelic as they are regional matches, not rsid matches. For the standard fingerprint SNPs, rsID will not be multi-allelic. 
                    y_target = 1
                    seen_rsid.append(region[3])
            if (y_target == 1):
                on_target.write('\t'.join([str(i) for i in print_out]) + '\n')
            else:
                print_out = [record.CHROM, record.POS, record.REF, ','.join([str(i) for i in record.ALT]), record.ID, record.samples[0]['AD'][0], ','.join([str(i) for i in record.samples[0]['AD'][1:]])]
                # note, the AD here can be multi-allelic as they are regional matches, not rsid matches. For the standard fingerprint SNPs, rsID will not be multi-allelic. 
                off_target.write('\t'.join([str(i) for i in print_out]) + '\n')

        else:
            print_out = [record.CHROM, record.POS, record.REF, ','.join([str(i) for i in record.ALT]), record.ID, record.samples[0]['AD'][0], ','.join([str(i) for i in record.samples[0]['AD'][1:]]) ]
            # note, the AD here can be multi-allelic as they are off target variants. For the standard fingerprint SNPs, rsID will not be multi-allelic. 

            off_target.write('\t'.join([str(i) for i in print_out]) + '\n')

    missing_rsid = list(set(all_rsid + seen_rsid))
    missing_rsid = sorted(missing_rsid, key = sort_ids)
    
    for rsid in missing_rsid:
        if 'noID_Y' in str(rsid):
            print_out = ['.', '.', '.', 'noCall', 'chrY_region', '.', '.', d[rsid]['gene_with_target']]
            on_target.write('\t'.join([str(i) for i in print_out]) + '\n')
        elif 'noID_X' in str(rsid): 
            print_out = ['.', '.', '.', 'noCall', 'chrX_region', '.', '.', d[rsid]['gene_with_target']]
            on_target.write('\t'.join([str(i) for i in print_out]) + '\n')
        else:
            print_out = ['.', '.', '.', 'noCall', rsid, '.', '.', d[rsid]['gene_with_target']]
            on_target.write('\t'.join([str(i) for i in print_out]) + '\n')
