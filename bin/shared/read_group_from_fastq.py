#! /usr/bin/env python

""" 
 read_group_from_fastq.py

 Using a fastq file's name and the contents of its first line, 
 build the option string needed for bwa to mark every read, assuming Illumina
 casava 1.8 conventions.

 Input: the fastq file specified as argv[1], the first command line argument.
        Handles compressed or uncompressed fastqs.
 Output: the second command line argument, if specified, else, sys.stdout.

 Notes:
    We will usually be handling standard Illumina Casava 1.8+ output, which
    has a regular file naming format and read name format.  If any of the
    steps here fail, cause the pipeline to fail rather than producing
    untraceable output.
"""

import sys
import os
import re
import time
import gzip
import argparse
#import cga_version
try:
    import bz2file as bz2
except ImportError:
    import bz2


def parse_args():
    parser = argparse.ArgumentParser(version='V2.0')
    parser.add_argument('-p', '--picard', action='store_true',
                        help="Use Picard format for read group line")
    parser.add_argument('-t', '--tumor', action='store_true',
                        help="Sample is tumor in a tumor/normal pair")
    parser.add_argument('-n', '--normal', action='store_true',
                        help="Sample is normal in a tumor/normal pair")
    parser.add_argument('-s', '--sample_id', dest="sample_id",
                        help="SampleID of file")
    parser.add_argument('-o', '--output', dest="output_file",
                        help="Output file name [STDOUT]")
    parser.add_argument('fastq', nargs="+",
                        help="Path to fastq file for sample")

    args = parser.parse_args()

    if args.tumor:
        if args.normal:
            # Check for a conflict.
            parser.error("Must not specify both --tumor and --normal.")
        args.sample_type = "Tumor_"
    elif args.normal:
        args.sample_type = "Normal_"
    else:
        args.sample_type = ""

    return args


def multi_open(name):
    if name.endswith('.gz'):
        f = gzip.open(name)
    elif name.endswith('.bz2'):
        f = bz2.BZ2File(name)
    else:
        f = open(name)
    return f


def make_fake(args):
    """
    If we can't get adequate data from the file, use timestamps.
    :return:
    """
    # Sleep for 2 seconds, to make sure that a previous invocation
    # will have a different time stamp.
    time.sleep(2)

    ts = time.strftime('%H%M%S')

    id = 'ID_' + ts
    lb = 'LIB_' + ts
    sm = 'SAMPLE_' + ts
    bc = 'RUN_' + ts
    output(id, lb, sm, bc, args)
    sys.exit(0)


def main():
    #cga_version.parse_options()

    args = parse_args()

    # First get the info from the filename
    fn = os.path.split(args.fastq[0])[1]

    if 'fastq' not in fn and 'fq' not in fn:
        print >> sys.stderr, "Not seemingly a fastq file:", fn
        make_fake(args)
        # Does not return...

    # Now split the basename portion into its constituent parts.
    fn_parts = fn.split('_')

    # Scan for the "GES" starting a filename part.  If found,
    # That separates the Sample name portion from the Library name.
    # If GES is not found starting a part, use the whole filename
    # as both the Sample name and the Library name.
    # Maybe redo this with regular expressions, but for now, it works.
    pos = -1
    for n in range(len(fn_parts)):
        if fn_parts[n].startswith("GES"):
            pos = n
            break
    if pos == -1:
        if args.sample_id:
            fn = args.sample_id
        else: 
            # Didn't find the GES marker. Use the filename up to the end name.
            match = re.search('(.*)[._]R[12]_.*',fn)
            if match is not None:
                fn = match.group(1)
            else:
                # something is seriously odd here, but we'll just use the
                # whole filename
                pass
        cust_id = ges_id = fn
    else:
        cust_id = '_'.join(fn_parts[:pos])
        ges_parts = fn_parts[pos:]
        pos = 999  # Way bigger than the number of parts we'll see.
        for n in range(len(ges_parts)):
            if ges_parts[n] == 'R1' or ges_parts[n] == 'R2':
                pos = n
                break
        ges_id = '_'.join(ges_parts[:pos])

    # Sanity check that we have some amount of text for our fields. The
    # down stream tools can't tolerate empty fields in the read group
    # information.
    if not ges_id:
        ges_id = fn

    if not cust_id:
        cust_id = ges_id

    # Now the parts from the first readname--the first line of the file.
    # When split on ':', the readname contains
    # - the ID in the first four fields.
    # Note: the leading '@' needs to be stripped.
    try:
        inf = multi_open(args.fastq[0])
        line = inf.readline()
    except IOError, e:
        print >> sys.stderr, "Couldn't read the file: {0}\n    {1}". \
            format(fn, e.message)
        make_fake(args)
        # Does not return

    # Example line:
    # @HISEQ2000:190:D19U8ACXX:5:1101:1492:1901 1:N:0:TAGCTT
    parts = line[1:].strip().split(' ')
    read_name = parts[0]

    # Example read_name: HISEQ2000:190:D19U8ACXX:5:1101:1492:1901
    rparts = read_name.split(':')
    if len(rparts) >= 4:
        rparts = rparts[:4]

    # Try to add the bar code in:
    bar_code = "no_barcode"
    if len(parts) >= 2:
        # Example comment: 1:N:0:TAGCTT
        comment = parts[1]
        cparts = comment.split(':')
        if len(cparts) == 4:
            bar_code = cparts[3]
            rparts.append(bar_code)

    id = ':'.join(rparts)
    # Example id: HISEQ2000:190:D19U8ACXX:5:TAGCTT

    output(id, ges_id, cust_id, bar_code, args)

def output(id, ges_id, cust_id, bar_code, args):
    if args.output_file is not None:
        of = open(args.output_file, 'w')
    else:
        of = sys.stdout

    if args.picard:
        line = 'RGID={0}{1} RGLB={0}{2} ' \
               'RGPL=ILLUMINA RGSM={3} RGPU={4}'.\
            format(args.sample_type, id, ges_id, cust_id, bar_code)
    else :
        line = '@RG\\tID:{0}{1}\\tLB:{0}{2}\\tSM:{3}\\tPL:ILLUMINA'.\
            format(args.sample_type, id, ges_id, cust_id)
    # This needs to be a single line file; no terminating \n
    print >> of, line,
    if of != sys.stdout:
        of.close()

if __name__ == '__main__':
    main()
