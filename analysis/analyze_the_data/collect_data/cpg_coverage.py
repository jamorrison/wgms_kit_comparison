"""Module to process the cpg_covg data files."""
import numpy as np
import sys
import re
import os

def process_dist_table(fname):
    """Process table with number of CpGs and CpGs covered in a given region.

    Inputs -- fname - filename of table

    Returns -- tuple - sample name, dictionary of percentages of CpGs covered in region
    """
    sample = os.path.basename(fname).replace('_cpg_dist_table.txt','')
    categories = ['TotalCpGs', 'CGICpGs', 'ExonicCpGs', 'GenicCpGs', 'RepeatCpGs']

    data = dict((cat + '_percent_covered', {}) for cat in categories)
    with open(fname, 'r') as f:
        lines = f.read()

        for cat in categories:
            a = re.search('{}\s+(All)\s+(\d+)\s+(\d+)'.format(cat), lines, re.MULTILINE)
            q = re.search('{}\s+(Q40)\s+(\d+)\s+(\d+)'.format(cat), lines, re.MULTILINE)
            if a is not None and q is not None:
                data[cat+'_percent_covered'][a.group(1)] = 100.0 * float(a.group(3)) / float(a.group(2))
                data[cat+'_percent_covered'][q.group(1)] = 100.0 * float(q.group(3)) / float(q.group(2))

    return sample, data

def process_depth_table(fname):
    """Process table with depth and number of CpGs with that depth.

    Inputs -- fname - filename of table

    Returns -- tuple - sample name, dictionary of {region: weighted average depth}
    """
    base = os.path.basename(fname).replace('_table.txt', '').replace('_cpgs', '').split('_')
    sample = base[0]
    key = '_'.join([*base[1:], 'avg_depth'])

    with open(fname, 'r') as f:
        lines = f.read().splitlines()[1:]

        depth = []
        fcpgs = []
        for line in lines:
            l = line.split('\t')
            depth.append(int(l[0]))
            fcpgs.append(int(l[1]))

        fcpgs = [i / sum(fcpgs) for i in fcpgs]
        
        weighted_avg = np.average(depth, weights=fcpgs)

    return sample, {key: weighted_avg}

if __name__ == '__main__':
    sam, out = process_dist_table('../cpg_covg/FtubeAkapaBC_cpg_dist_table.txt')
    print(sam, out)
    samp, outp = process_depth_table('../cpg_covg/FtubeAkapaBC_cpgs_total_all_table.txt')
    print(samp, outp)
