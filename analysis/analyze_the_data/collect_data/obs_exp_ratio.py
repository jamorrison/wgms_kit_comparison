"""Function to generate plots using provided data."""
import numpy as np
import glob
import json

SAMPL = 0
REFLN = 1
MAPLN = 8
ORDER = {
    'cpgs_obs_exp_ratio': (2, 9),
    'cpis_obs_exp_ratio': (3,10),
    'rmsk_obs_exp_ratio': (4,11),
    'exon_obs_exp_ratio': (5,12),
    'gene_obs_exp_ratio': (6,13),
    'intr_obs_exp_ratio': (7,14)
}

def process_file(fname, prefix=''):
    """Process log file (*.stdout extension) from obs/exp processing.

    Inputs -- fname - filename of log file

    Returns -- cleaned dictionary with obs/exp ratios
    """
    vals = np.genfromtxt(fname, dtype=None, encoding=None)[()]

    samp = vals[SAMPL]
    ratios = {}
    for key, val in ORDER.items():
        den = vals[val[0]] / vals[REFLN]
        num = vals[val[1]] / vals[MAPLN]

        ratios[prefix+key] = num / den

    return samp, ratios

def main():
    """Function for reading files and processing them."""
    TOPDIR='2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis'
    obs_exp_mappy_files = glob.glob(TOPDIR + '/analyze_the_data/cpg_questions/exp_vs_obs_coverage/pbs_mappability/*.stdout')

    data_mappy = {}
    for log in obs_exp_mappy_files:
        samp, ratios = process_file(log, prefix='mappy_')
        data_mappy[samp] = ratios

    print(data_mappy)

    return 0

if __name__ == '__main__':
    main()
