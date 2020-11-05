from Bio import SeqIO
import argparse
import numpy as np
import gzip
import glob
import time

# Get metric values from list of Phred scores
def base_metrics(phred_scores):
    avg_qual = np.mean(phred_scores)

    no_bases = len(phred_scores)

    lo_count = 0
    for score in phred_scores:
        if score < 20:
            lo_count += 1

    return avg_qual, no_bases, lo_count

# Set up command line arguments
parser = argparse.ArgumentParser(
    description = 'base_quality.py finds base quality metrics for input FASTQ files'
)

parser.add_argument(
    'fastq_1',
    metavar = 'fastq_1',
    type = str,
    help = 'Input FASTQ file for Read 1'
)

parser.add_argument(
    'fastq_2',
    metavar = 'fastq_2',
    type = str,
    help = 'Input FASTQ file for Read 2'
)

args = parser.parse_args()

# Initialize variables
r1_names = []           # List of read names for read 1
r2_names = []           # List of read names for read 2
r1_totrs = r2_totrs = 0 # Total number of reads processed
r1_totbs = r2_totbs = 0 # Total number of bases processed
r1_avg20 = r2_avg20 = 0 # Number of reads with avg. base quality >= 20
r1_avg30 = r2_avg30 = 0 # Number of reads with avg. base quality >= 30
r1_cnt20 = r2_cnt20 = 0 # Number of bases with base quality < 20

print('Begin processing {}'.format(args.fastq_1))
t1_start = time.time()
with gzip.open(args.fastq_1, 'rt') as f1:
    for record in SeqIO.parse(f1, 'fastq'):
        r1_names.append(record.id)

        avg, num, low = base_metrics(record.letter_annotations['phred_quality'])

        if avg >= 20:
            r1_avg20 += 1
        if avg >= 30:
            r1_avg30 += 1

        r1_totrs += 1
        r1_totbs += num
        r1_cnt20 += low
t1_end = time.time()
print('Processing {} took {:.2f} seconds'.format(args.fastq_1, t1_end-t1_start))

print('Begin processing {}'.format(args.fastq_2))
t2_start = time.time()
with gzip.open(args.fastq_2, 'rt') as f2:
    for record in SeqIO.parse(f2, 'fastq'):
        r2_names.append(record.id)

        avg, num, low = base_metrics(record.letter_annotations['phred_quality'])

        if avg >= 20:
            r2_avg20 += 1
        if avg >= 30:
            r2_avg30 += 1

        r2_totrs += 1
        r2_totbs += num
        r2_cnt20 += low
t2_end = time.time()
print('Processing {} took {:.2f} seconds'.format(args.fastq_2, t2_end-t2_start))

if r1_names == r2_names:
    print('All reads are ordering correctly!')
else:
    print('There is a read out of order in:', args.fastq_1, 'and', args.fastq_2)

print(
    '''
    Read 1 filename: {}
    Read 1 number of reads: {}
    Read 1 number of bases: {}
    Read 1 % of reads with avg. base quality >= 20: {:.2f}
    Read 1 % of reads with avg. base quality >= 30: {:.2f}
    Read 1 % of bases with base quality < 20: {:.2f}
    '''.format(args.fastq_1,
                r1_totrs,
                r1_totbs,
                100*r1_avg20/r1_totrs,
                100*r1_avg30/r1_totrs,
                100*r1_cnt20/r1_totbs
                )
)

print(
    '''
    Read 2 filename: {}
    Read 2 number of reads: {}
    Read 2 number of bases: {}
    Read 2 % of reads with avg. base quality >= 20: {:.2f}
    Read 2 % of reads with avg. base quality >= 30: {:.2f}
    Read 2 % of bases with base quality < 20: {:.2f}
    '''.format(args.fastq_2,
                r2_totrs,
                r2_totbs,
                100*r2_avg20/r2_totrs,
                100*r2_avg30/r2_totrs,
                100*r2_cnt20/r2_totbs
                )
)

