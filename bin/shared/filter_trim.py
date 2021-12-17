#! /usr/bin/env python

"""
filter_trim.py
Background

Many of the various analysis pipelines implemented by CS, including
the clinical single_sample_exome pipeline utilize two perl scripts to
filter out low quality reads and to trim 3' bases whose quality  is
less than a specified cutoff. This process takes several hours of
elapsed time, largely due to additional work that is performed that
we don't require.

There are seven command line options controlling processing:
1.  filter_quality: The quality value considered high quality when
    filtering reads.
2.  filter_percent_bases: The percentage of a reads bases that must
    be of high quality to consider the read to be high quality.
3.  trimming_quality: The quality value needed for a base to terminate
    trimming.
4.  post_trim_length_pct: The minimum length of a read required post
    trimming to retain the read, expressed as a percentage of the
    initial read length.
5.  trim_5: Whether to trim bases from the 5' end of the read as well
    as the the 3' end.
6.  min_pct_hq_reads: The minimum percentage of high quality reads out
    of all reads, required for successful exit. Fewer high quality reads
    will result in a failure exit status return to the shell, allowing
    pipelines using this program to fail.  Default: 0.0 (always return
    success).
7.  single_end: Use single end mode, even if there are multiple fastq
    files presented.

Trimming occurs per read; trimming does not need to match read 1 vs
read 2 for paired end data. Filtering and acting on the
post_trim_length_min occurs on a per-end basis; however if one read is
discarded due to these criteria, the other end's read must be discarded
as well.

Inputs:
1.  A sequence of pairs of fastq files from a paired-end sequencing run
    or runs.  All pairs must be for the same sample.

Outputs:
1.  One or two filtered, trimmed fastq files.  The paired reads may not
    be the same length due to trimming.
2.  Statistics file: The cga pipeline uses several criteria from this
    step to determine whether the run was good enough to analyze.
    a. Percent HQ reads: Number of reads written to the filtered
       file(s) / number of reads in the input file(s).

    We also report the following statistics, but they aren't used to
    judge the quality of a run:
    b. Total Reads (Reported once, since both ends must be the same.)
    c. Total HQ reads (Reported once, since both ends must be the
       same.)
    d. Min, Max and Mean trimmed read lengths for reads whose trimmed
       length is sufficient to retain the read, reported separately for
       each end.

All output file naming is based on the first file in each fastq file
list the following input fastqs are not represented in the output file
names.

Exit status:
0 = Success
non-0 = Failure:
  1) Could not open an input file or create an output file
  2) Python 2.7 or higher not used
  3) Insufficient high quality reads
  4) Odd number of fastq files when in paired mode
  5) Input files in a pair are not the same length
  6) Zero total reads, percent of high quality reads not computed
"""
__author__ = 'simons'

import sys
import os
import math
import gzip
import datetime
import inspect

# In Python 2.7, the core bz2 module can't process multi-stream files, such
# as those produced by pbzip2.  In Python 3.4 and above, it can.  The 
# Python 3 version has been backported to Python 2.7, and is available
# as bz2file. Conditionally load it to handle a wider array of files;
# fall back to the core bz2 module.
try:
    import bz2file as bz2
except:
    print >> sys.stderr, 'Could not import bz2file; using the core bz2.'
    import bz2

# Support the version command.
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(
    inspect.getfile(inspect.currentframe()))[0]))
lib_folder = os.path.join(cmd_folder, '../lib')
if lib_folder not in sys.path:
    sys.path.insert(0, lib_folder)

try:
    import argparse
except:
    print >> sys.stderr, 'This program requires python V2.7 or higher.'
    sys.exit(2)


# The guts of this program.  All processing of reads.
class FastqRead(object):
    trim_5 = False
    trim_hq = 30
    read_hq = 30
    pct_hq = 0.7

    def __init__(self, fastqs, odir=None, suffix='_filtered_trimmed'):
        ofn = os.path.split(fastqs[0])[1] + suffix
        if odir:
            ofn_path = os.path.join(odir, ofn)
        else:
            ofn_path = ofn
        try:
            self.of = open(ofn_path, 'w')
        except IOError:
            print >> sys.stderr, \
                'Could not open "{0}". Exiting.'.format(ofn_path)
            sys.exit(1)

        self.fastqs = fastqs
        self.fn = None
        self.f = None
        self.total_reads = 0
        self.hq_reads = 0
        self.output_reads = 0
        self.min_trimmed_length = sys.maxint
        self.max_trimmed_length = -1
        self.total_trimmed_length = 0
        self.trimmed_reads = 0
        self.name = ''
        self.bases = ''
        self.plus = ''
        self.qual = ''
        self.timestamp = False
        self.line_count = 0

        # Initialize our first input fastq
        self.next_file()

    def next_file(self):
        #Do we have any files left?
        if not self.fastqs:
            return False

        self.fn = self.fastqs.pop(0)
        try:
            self.f = FastqRead.open(self.fn)
        except IOError:
            print >> sys.stderr, \
                'Could not open "{0}". Exiting.'.format(self.fn)
            sys.exit(1)
        return True

    def get_filename(self):
        return self.fn

    @staticmethod
    def open(name):
        """
        Intended to be private to the class...

        A flexible open routine that can handle plain text files or
        files compressed with gzip or bzip2.  Only used for the
        input files. Output files are emitted uncompressed, until the
        tools in the next leg of the pipeline can work properly with
        compressed files.

        :param name: The filename to open.
        :return: A file object for the named file.
        """
        if name.endswith('.gz'):
            f = gzip.open(name)
        elif name.endswith('.bz2'):
            f = bz2.BZ2File(name)
        else:
            f = open(name)
        return f

    def stats(self):
        s = {}
        s['total_reads'] = self.total_reads
        s['hq_reads'] = self.hq_reads
        s['output_reads'] = self.output_reads
        s['max_trimmed_length'] = self.max_trimmed_length
        s['min_trimmed_length'] = self.min_trimmed_length
        try:
            tmp_mean = float(self.total_trimmed_length) / \
                       float(self.trimmed_reads)
            s['mean_trimmed_length'] = '{0:.2f}'.format(tmp_mean)
        except ZeroDivisionError:
            s['mean_trimmed_length'] = 'N/A'

        return s

    def next(self):
        """

        :return: True iff the read was successfully retrieved from
        the file.
        """
        name = self.f.readline()
        # Test whether we had a successful read.
        #  Will be zero length if EOF reached.
        if not name:
            return False
        self.name = name.strip()
        self.bases = self.f.readline().strip()
        self.plus = self.f.readline().strip()
        self.qual = self.f.readline().strip()

        # All four lines must have content to be a valid read.
        if len(self.bases) == 0 or \
           len(self.plus) == 0 or \
           len(self.qual) == 0:
            print >> sys.stderr, 'NAME:', self.name
            print >> sys.stderr, 'BASES:', self.bases
            print >> sys.stderr, 'PLUS:', self.plus
            print >> sys.stderr, 'QUAL:', self.qual
            raise ValueError('Incomplete read found in file {0}'.
                             format(self.fn))

        self.total_reads += 1
        if self.timestamp:
            self.line_count += 1
            if self.line_count % 1000000 == 0:
                print >> sys.stderr,  \
                    datetime.datetime.strftime(datetime.datetime.now(),
                                               '%H:%M:%S'), \
                    self.line_count
        return True

    def do_timestamp(self):
        self.timestamp = True

    @staticmethod
    def set_criteria(pct_hq=0.7,
                     read_hq=30,
                     trim_hq=30,
                     trim_5=False,
                     min_pct=0.7,
                     min_pct_hq_reads=0.0):

        FastqRead.pct_hq = float(pct_hq)
        if FastqRead.pct_hq > 1.0:
            FastqRead.pct_hq /= 100.0

        # Use phred33 quality scoring
        FastqRead.read_hq = chr(int(read_hq) + 33)
        FastqRead.trim_hq = chr(int(trim_hq) + 33)
        FastqRead.trim_5 = trim_5    # Passed in as boolean
        FastqRead.min_pct = float(min_pct)
        if FastqRead.min_pct > 1.0:
            FastqRead.min_pct /= 100.0
        FastqRead.min_pct_hq_reads = float(min_pct_hq_reads)
        if FastqRead.min_pct_hq_reads > 1.0:
            FastqRead.min_pct_hq_reads /= 100.0

    # Cache the minimum length of a trimmed read.
    min_len = None

    def trim(self):
        """

        :return:  True if the read is long enough after trimming.
        """

        original_length = len(self.qual)
        if FastqRead.trim_5:
            for p5 in range(original_length):
                if self.qual[p5] >= FastqRead.trim_hq:
                    break
        else:
            p5 = 0

        for p3 in range(original_length - 1, -1, -1):
            if self.qual[p3] >= FastqRead.trim_hq:
                break

        tlg = (p3 - p5) + 1 # Length after trimming.

        if FastqRead.min_len is None:
            FastqRead.min_len = \
                math.ceil(FastqRead.min_pct * original_length)
        if tlg < FastqRead.min_len:
            return False

        self.bases = self.bases[p5:p3 + 1]
        self.qual = self.qual[p5:p3 + 1]

        assert tlg == len(self.qual), "Length calculation is broken"

        # Track our trimmed length stats
        if tlg > self.max_trimmed_length:
            self.max_trimmed_length = tlg
        if tlg < self.min_trimmed_length:
            self.min_trimmed_length = tlg
        self.total_trimmed_length += tlg  # To compute the mean
        self.trimmed_reads += 1
        return True

    def filter(self):
        """

        :return: True if the read passed HQ filtering criteria
        """
        lg = len(self.qual)
        lq_reads_allowed = math.floor(float(lg) *
                           (1.0 - FastqRead.pct_hq))
        lq_reads = 0
        for n in range(lg):
            if self.qual[n] < FastqRead.read_hq:
                lq_reads += 1
                if lq_reads > lq_reads_allowed:
                    return False
        self.hq_reads += 1
        return True

    def write(self):
        print >> self.of, self.name
        print >> self.of, self.bases
        print >> self.of, self.plus
        print >> self.of, self.qual
        self.output_reads += 1

    def close(self):
        self.of.close()

    # End of class FastqRead.


def parse_args():

    parser = argparse.ArgumentParser(description=
                "Perform filtering and trimming of paired end fastq "
                "files", version="2.0")
    parser.add_argument("-p", "--hq_pct", default="70", help=
                "Percentage of bases that must be high quality [70]")
    parser.add_argument('-f', '--filter_hq', default="30", help=
                'Numeric quality value to pass filtering [30]')
    parser.add_argument('-t', '--trim_hq', default="30", help=
                'Numeric quality value to not be trimmed [30]')
    parser.add_argument('-m', '--min_len_pct', default="70", help=
                'Minimum read length after trimming to '
                'retain read. (percentage 0-100) [70]')
    parser.add_argument('-M', '--min_pct_hq_reads', default=0, help=
                'Minimum percentage of reads classified as High '
                'Quality reads (percentage 0-100) [0]')
    parser.add_argument('-5', '--trim_5', action="store_true", help=
                "Trim 5' end as well as 3' [False]")
    parser.add_argument('-s', '--suffix', default='_filtered_trimmed',
                        help='Suffix to construct the output file name '
                '[_filtered_trimmed]')
    parser.add_argument('-S', '--single_end', action="store_true",
                 help="Use single end mode with multiple fastq files " \
                      "[False]")
    parser.add_argument('-d', '--directory', dest='odir', default='.',
                        help=
                'Directory in which to write the output files '
                '[current directory]')
    parser.add_argument('-i', '--timestamp', action='store_true', help=
                'Emit a timestamp ever 1,000,000 reads [False]')
    parser.add_argument("fastqs", nargs="+")
    args = parser.parse_args()
    return args


def output_stats_single(r1, args, start_time):
    """
    Report the statistics for a single end run.

    NOTE WELL!!!
    This routine and output_stats_paired have the same logic
    flow.  If one is changed, the other almost certainly has to
    have the corresponding change made.

    YOU HAVE BEEN WARNED!

    :param r1: Accumulated info for the reads.
    :param args: Our command line arguments
    :param start_time: The run's start time
    :return: None
    """

    # Here we have completely processed the input file. Write out
    # the statistics
    r1_stats = r1.stats()
    bn_fq1 = os.path.split(args.fastqs[0])[1]

    with open(os.path.join(args.odir, bn_fq1 + '_stat'), 'w') as sf:
        print >> sf, 'Input file:'
        print >> sf, 'Read 1: {0}'.format(args.fastqs[::2])
        print >> sf, 'QC statistics'
        print >> sf, 'Statistic\tRead 1'

        try:
            f_pct_hq = float(r1_stats['output_reads']) / \
                float(r1_stats['total_reads'])
            pct_hq = '{0:.2%}'.format(f_pct_hq)
        except ZeroDivisionError:
            pct_hq = 'N/A'
        print >> sf, 'Percentage of HQ reads\t{0}'.format(pct_hq)

        print >> sf, 'Total number of reads\t{0}'.format(
                     r1_stats['total_reads'])
        print >> sf, 'Total number of HQ filtered reads\t{0}'.\
            format(r1_stats['output_reads'])
        print >> sf, 'Detailed QC statistics'
        print >> sf, 'Reads passing filter\t{0}'.\
            format(r1_stats['hq_reads'])

        try:
            pct_rpf = '{0:.2%}'.format(float(r1_stats['hq_reads']) /
                                     float(r1_stats['total_reads']))
        except ZeroDivisionError:
            pct_rpf = 'N/A'
        print >> sf, 'Percent reads passing filter\t{0}'.format(pct_rpf)

        print >> sf, 'Max Trimmed Length\t{0}'.\
            format(r1_stats['max_trimmed_length'])
        print >> sf, 'Min Trimmed Length\t{0}'.\
            format(r1_stats['min_trimmed_length'])
        print >> sf, 'Mean Trimmed Length\t{0}'.\
            format(r1_stats['mean_trimmed_length'])
        print >> sf, 'Run start time\t{0}'.\
            format(datetime.datetime.strftime(start_time, '%H:%M:%S'))
        end_time = datetime.datetime.now()
        print >> sf, 'Run end time\t{0}'.\
            format(datetime.datetime.strftime(end_time, '%H:%M:%S'))

        if r1_stats['total_reads'] == 0:
            # This will be the same as sys.exit(6)
            print >> sys.stderr, 'Failure: total reads == 0\nExiting' \
                                 'with status 6'
            return 6
        if f_pct_hq < FastqRead.min_pct_hq_reads:
            # This will be the same effect as sys.exit(3)
            print >> sys.stderr, 'Failure: not enough high quality ' \
                                 'read percent: {} required: {}\n' \
                                 'Exiting with status 3'.\
                format(f_pct_hq, FastqRead.min_pct_hq_reads)
            return 3
        # Success!
        return 0

def output_stats_paired(r1, r2, args, start_time):
    """
    Report the statistics for a paired end run.

    NOTE WELL!!!
    This routine and output_stats_single have the same logic
    flow.  If one is changed, the other almost certainly has to
    have the corresponding change made.

    YOU HAVE BEEN WARNED!

    :param r1: Accumulated info for the end 1 reads.
    :param r2: Accumulated info for the end 2 reads.
    :param args: Our command line arguments
    :param start_time: The run's start time
    :return: None
    """

    # Here we have completely processed both input files. Write out
    # the statistics
    r1_stats = r1.stats()
    r2_stats = r2.stats()

    bn_fq1 = os.path.split(args.fastqs[0])[1]
    bn_fq2 = os.path.split(args.fastqs[1])[1]

    with open(os.path.join(args.odir,
            '{0}_{1}_stat'.format(bn_fq1, bn_fq2)), 'w') as sf:
        print >> sf, 'Input files:'
        print >> sf, 'Read 1: {0}'.format(args.fastqs[::2])
        print >> sf, 'Read 2: {0}'.format(args.fastqs[1::2])
        print >> sf, 'QC statistics'
        print >> sf, 'Statistic\tRead 1\tRead 2'

        try:
            f_pct_hq1 = float(r1_stats['output_reads']) / \
                float(r1_stats['total_reads'])
            pct_hq1 = '{0:.2%}'.format(f_pct_hq1)
        except ZeroDivisionError:
            pct_hq1 = 'N/A'
        try:
            f_pct_hq2 = float(r2_stats['output_reads']) / \
                float(r2_stats['total_reads'])
            pct_hq2 = '{0:.2%}'.format(f_pct_hq2)
        except ZeroDivisionError:
            pct_hq2 = 'N/A'
        print >> sf, 'Percentage of HQ reads\t{0}\t{1}'.\
            format(pct_hq1, pct_hq2)

        print >> sf, 'Total number of reads\t{0}\t{1}'.format(
                     r1_stats['total_reads'],
                     r2_stats['total_reads'])
        print >> sf, 'Total number of HQ filtered reads\t{0}\t{1}'.\
            format(r1_stats['output_reads'], r2_stats['output_reads'])
        print >> sf, 'Detailed QC statistics'
        print >> sf, 'Reads passing filter\t{0}\t{1}'.\
            format(r1_stats['hq_reads'], r2_stats['hq_reads'])

        try:
            pct_rpf1 = '{0:.2%}'.format(float(r1_stats['hq_reads']) /
                                      float(r1_stats['total_reads']))
        except ZeroDivisionError:
            pct_rpf1 = 'N/A'
        try:
            pct_rpf2 = '{0:.2%}'.format(float(r2_stats['hq_reads']) /
                                      float(r2_stats['total_reads']))
        except ZeroDivisionError:
            pct_rpf2 = 'N/A'
        print >> sf, 'Percent reads passing filter\t{0}\t{1}'.\
            format(pct_rpf1, pct_rpf2)

        print >> sf, 'Max Trimmed Length\t{0}\t{1}'.\
            format(r1_stats['max_trimmed_length'],
                   r2_stats['max_trimmed_length'])
        print >> sf, 'Min Trimmed Length\t{0}\t{1}'.\
            format(r1_stats['min_trimmed_length'],
                   r2_stats['min_trimmed_length'])
        print >> sf, 'Mean Trimmed Length\t{0}\t{1}'.\
            format(r1_stats['mean_trimmed_length'],
                   r2_stats['mean_trimmed_length'])
        print >> sf, 'Run start time\t{0}'.\
            format(datetime.datetime.strftime(start_time, '%H:%M:%S'))
        end_time = datetime.datetime.now()
        print >> sf, 'Run end time\t{0}'.\
            format(datetime.datetime.strftime(end_time, '%H:%M:%S'))

        if r1_stats['total_reads'] == 0 or r2_stats['total_reads'] == 0:
            # This will be the same as sys.exit(6)
            print >> sys.stderr, 'Failure: total reads == 0\nExiting' \
                                 'with status 6'
            return 6
        if f_pct_hq1 < FastqRead.min_pct_hq_reads or \
            f_pct_hq2 < FastqRead.min_pct_hq_reads:
            # This will be the same effect as sys.exit(3)
            print >> sys.stderr, 'Failure: not enough high quality ' \
                                 'read percent: e1: {}, e2: {} ' \
                                 'required: {}\n' \
                                 'Exiting with status 3'.\
                format(f_pct_hq1, f_pct_hq2, FastqRead.min_pct_hq_reads)
            return 3
        # Success!
        return 0

def main():
    start_time = datetime.datetime.now()
    args = parse_args()

    # If we are doing paired end processing, make sure that we have
    # pairs (i.e., an even number of files, and split the list of
    # files into end-specific lists.
    num_fastqs = len(args.fastqs)
    paired_end = ((num_fastqs != 1) and (not args.single_end))

    if paired_end:
        # Paired end; need to be an even number of fastqs.
        if num_fastqs % 2 != 0:
            print >> sys.stderr, 'Odd number of fastq files ({0}) in ' \
                                 'paired-end mode. Exiting...'.format(
                num_fastqs
            )
            sys.exit(4)

        # Now split the lists:
        e1_fastqs = args.fastqs[::2]
        e2_fastqs = args.fastqs[1::2]
    else:
        # Make a copy.  We need the original later.
        e1_fastqs = args.fastqs[:]
        e2_fastqs = None

    r1 = FastqRead(e1_fastqs, args.odir, args.suffix)

    # We may be processing single end reads.  Everything with r2 is
    # conditional on having a second fastq.
    if paired_end:
        r2 = FastqRead(e2_fastqs, args.odir, args.suffix)

    # Check if we want timestamps output to track progress
    if args.timestamp:
        r1.do_timestamp()

    # The criteria are class members, not instance.
    FastqRead.set_criteria(args.hq_pct, args.filter_hq, args.trim_hq,
                           args.trim_5, args.min_len_pct,
                           args.min_pct_hq_reads)

    r1_ok = False

    # If we don't have paired end reads, we just want the tests
    # below to care about end 1.  In this case, initialize R2_ok to
    # True
    r2_ok = not paired_end

    # Loop over the whole file.  We'll exit this with a break.
    while True:
        # Do NOT move these into the if statement below; we need to
        # keep them in sync. If they are in the if, and r1 fails,
        # r2 will not be executed.
        r1_ok = r1.next()
        if paired_end:
            r2_ok = r2.next()
        if not (r1_ok and r2_ok):
            # One or both files are exhausted. Must both end at the
            # same read.
            if r1_ok or (paired_end and r2_ok):
                print >> sys.stderr, \
                    'Input files {0} and {1} are different lengths.\n' \
                    'Exiting.'.format(
                        r1.get_filename(),
                        r2.get_filename())
                sys.exit(5)
            # Get the next files in the list to continue processing.
            # Since we ensured above that the lists
            # were the same length, we don't need to do equivalency
            # tests here.  We can simply test for list exhaustion on r1
            # which works for both single and paired end.  If it
            # succeeds and we're paired end, we can blindly get the next
            # end 2 file.
            r1_ok = r1.next_file()
            if not r1_ok:
                # We've exhausted the list of input files.
                break
            if paired_end:
                # Guaranteed to succeed: lists are equal length.
                r2.next_file()

            # Back to the top to get a read from the new files
            continue

        r1_ok = r1.filter()
        if paired_end:
            r2_ok = r2.filter()
        if not (r1_ok and r2_ok):
            # Filtering this read failed... Next!
            continue

        r1_ok = r1.trim()
        if paired_end:
            r2_ok = r2.trim()
        if not (r1_ok and r2_ok):
            # This read trimmed to be too short.
            continue

        r1.write()
        if paired_end:
            r2.write()


    if paired_end:
        status = output_stats_paired(r1, r2, args, start_time)
        r2.close()
    else:
        status = output_stats_single(r1, args, start_time)
    r1.close()

    return status

if __name__ == '__main__':
    status = main()
    sys.exit(status)
