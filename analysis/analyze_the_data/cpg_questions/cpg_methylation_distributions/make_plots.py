"""Create violin plots of CpG methylation in different regions."""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

def extract_sample_info(sample_str):
    """Extract kit, sample, and technical replicate from sample_str.

    Inputs -
        sample_str - string from sample name
    Returns -
        tuple (kit, biological sample name, technical replicate)
    """
    s = sample_str.replace('Ftube', '')

    # The biological sample in now the first character in name
    bio = s[0]

    # Extract what kit is in sample
    kit = ''
    if 'kapabc' in s.lower():
        kit = 'Kapa'
    elif 'pbat' in s.lower():
        kit = 'PBAT'
    elif 'neb' in s.lower():
        kit = 'NEB'
    elif 'swift' in s.lower():
        kit = 'Swift'

    # Determine if low or high input
    if '10ng' in s:
        kit = 'Low ' + kit

    # Determine technical replicate
    rep = '1'
    if 'rep2' in s.lower():
        rep = '2'

    if (bio not in ['A', 'B']) or (kit == ''):
        print('[extract_sample_info] ERROR: Incorrect entry')
        return ('', '', '')

    return (kit, bio, rep)

def import_files(tag, samp):
    """Import CpG islands, shores, shelves, and open seas BED files into a 
       DataFrame

    Inputs -
        tag  - sample name to load (include any paths needed)
        samp - only sample name (no paths)
    Returns -
        DataFrame with beta values, sea_group, sample, and replicate columns
    """
    names = ['chr','start','end','beta','covg']
    cols = [3]

    isld = pd.read_csv(
        tag+'.island.bed.gz', sep='\t', header=None, usecols=cols
    )
    isld.rename(columns={3: 'beta'}, inplace=True)
    isld['sea_group'] = 'islands'

    shor = pd.read_csv(
        tag+'.shores.bed.gz', sep='\t', header=None, usecols=cols
    )
    shor.rename(columns={3: 'beta'}, inplace=True)
    shor['sea_group'] = 'shores'

    shlv = pd.read_csv(
        tag+'.shelves.bed.gz', sep='\t', header=None, usecols=cols
    )
    shlv.rename(columns={3: 'beta'}, inplace=True)
    shlv['sea_group'] = 'shelves'

    seas = pd.read_csv(
        tag+'.open_seas.bed.gz', sep='\t', header=None, usecols=cols
    )
    seas.rename(columns={3: 'beta'}, inplace=True)
    seas['sea_group'] = 'open_seas'

    comb = pd.concat([isld, shor, shlv, seas])
    kit, bio, rep = extract_sample_info(samp)
    comb['sample'] = kit + ' Samp. ' + bio
    comb['Replicate'] = 'Rep. ' + rep

    return comb

def create_plot(data, title, xlab, ylab, figname):
    """Create split violin plot of methylation data.

    Inputs -
        data    - list of data to include in violinplot
        title   - title of plot
        xlab    - x-axis label of plot
        ylab    - y-axis label of plot
        figname - name of output file for plot
    Returns -
        Nothing, plot saved to disk
    """
    fig, ax = plt.subplots(figsize=(10,5))
    plt.tight_layout()

    sns.violinplot(data=data, x='beta', y='sample', hue='Replicate', cut=0,
                   split=True, inner='quartile', linewidth=1, orient='h',
                   palette={'Rep. 1': '#005596', 'Rep. 2': '#3fa294'})

    ax.legend(title='', ncol=2, bbox_to_anchor=(0.5,1), frameon=False,
              loc='lower center', fontsize='large')
    plt.xlim(-0.05, 1.05)

    plt.title(title, pad=40, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def main():
    """Main function."""
    samps = [
        'FtubeAkapaBC', 'FtubeAkapaBCrep2',
        'FtubeAneb', 'FtubeAnebRep2',
        'FtubeApbat',
        'FtubeAswift', 'FtubeAswiftRep2',
        'FtubeAneb10ng', 'FtubeAneb10ngRep2',
        'FtubeAswift10ng', 'FtubeAswift10ngRep2'
    ]

    data = []
    for samp in samps:
        c_a = import_files('seashore_bed_files/' + samp, samp)
        c_b = import_files(
            'seashore_bed_files/' + samp.replace('FtubeA', 'FtubeB'),
            samp.replace('FtubeA', 'FtubeB')
        )
        data.append(c_a)
        data.append(c_b)

    comb = pd.concat(data)
    isld = comb[comb.sea_group == 'islands']
    shor = comb[comb.sea_group == 'shores']
    shlv = comb[comb.sea_group == 'shelves']
    seas = comb[comb.sea_group == 'open_seas']

    create_plot(
        isld,
        'CpG Island Methylation Distribution',
        'Methylation Level',
        '',
        'islands.pdf'
    )

    create_plot(
        shor,
        'CpG Shore Methylation Distribution',
        'Methylation Level',
        '',
        'shores.pdf'
    )

    create_plot(
        shlv,
        'CpG Shelve Methylation Distribution',
        'Methylation Level',
        '',
        'shelves.pdf'
    )

    create_plot(
        seas,
        'CpG Open Sea Methylation Distribution',
        'Methylation Level',
        '',
        'opensea.pdf'
    )

    create_plot(
        comb,
        'All CpG Methylation Distribution',
        'Methylation Level',
        '',
        'all.pdf'
    )

if __name__ == '__main__':
    main()
