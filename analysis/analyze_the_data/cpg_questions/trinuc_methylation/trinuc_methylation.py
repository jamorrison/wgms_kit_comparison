"""Retrieve CAG, CAH, CTG, CTH trinucleotide methylation values."""
import pandas as pd
import glob
import time
import os

def read_file(fname, out_name):
    """Read base-averaged trinucleotide methylation values from file.

    Inputs -
        fname    - file to retrieve values from
        out_name - name of file to write results to
    Returns -
        Nothing, writes values to file
    """
    # Names of columns in file
    cols = [
        'chr',       # chromosome
        'start',     # start location
        'end',       # end location
        'ref',       # reference base
        'group',     # context group (CG, CHG, CHH)
        'two_base',  # 2-base context
        'five_base', # 5-base context
        'beta',      # beta value
        'covg'       # loci coverage
    ]
    cah_sum = 0 ; cah_len = 0 # Variables for base-averaged CAH methylation
    cag_sum = 0 ; cag_len = 0 # Variables for base-averaged CAG methylation
    cth_sum = 0 ; cth_len = 0 # Variables for base-averaged CTH methylation
    ctg_sum = 0 ; ctg_len = 0 # Variables for base-averaged CTG methylation

    # Read BED file in chunks to ease memory load
    for chunk in pd.read_csv(fname, sep='\t', header=None, names=cols, chunksize=1e6):
        cah_sum += chunk[(chunk.group == 'CHH') & (chunk.two_base == 'CA')]['beta'].sum()
        cah_len += len(chunk[(chunk.group == 'CHH') & (chunk.two_base == 'CA')]['beta'])

        cag_sum += chunk[(chunk.group == 'CHG') & (chunk.two_base == 'CA')]['beta'].sum()
        cag_len += len(chunk[(chunk.group == 'CHG') & (chunk.two_base == 'CA')]['beta'])

        cth_sum += chunk[(chunk.group == 'CHH') & (chunk.two_base == 'CT')]['beta'].sum()
        cth_len += len(chunk[(chunk.group == 'CHH') & (chunk.two_base == 'CT')]['beta'])

        ctg_sum += chunk[(chunk.group == 'CHG') & (chunk.two_base == 'CT')]['beta'].sum()
        ctg_len += len(chunk[(chunk.group == 'CHG') & (chunk.two_base == 'CT')]['beta'])

    # Calculate base-averaged methylation value
    cah_meth = 100 * cah_sum / cah_len
    cag_meth = 100 * cag_sum / cag_len
    cth_meth = 100 * cth_sum / cth_len
    ctg_meth = 100 * ctg_sum / ctg_len

    with open(out_name, 'w') as f:
        f.write(
            '{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(
                cah_meth, cag_meth, cth_meth, ctg_meth
            )
        )

def main():
    """Process cytosine context files."""
    raw_path = '2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/align/'
    sub_path = '2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/analyze_the_data/subsampling/'
    raw_ext = '.c.context.sorted.bed.gz'
    sub_ext = '.subsampled.c.context.sorted.bed.gz'

    run_raw = True
    run_sub = True

    if run_raw == True:
        files = glob.glob(raw_path+'*'+raw_ext)
        for f in files:
            samp = os.path.basename(f).replace(raw_ext, '')

            print(f'Processing {samp}_raw', end=' ... ')
            t1 = time.time()
            read_file(f, samp+'_raw.tsv')
            t2 = time.time()
            print(f'in {t2-t1} seconds')

    if run_sub == True:
        files = glob.glob(sub_path+'*'+sub_ext)
        for f in files:
            samp = os.path.basename(f).replace(sub_ext, '')

            print(f'Processing {samp}_sub', end=' ... ')
            t1 = time.time()
            read_file(f, samp+'_sub.tsv')
            t2 = time.time()
            print(f'in {t2-t1} seconds')

if __name__ == '__main__':
    main()
