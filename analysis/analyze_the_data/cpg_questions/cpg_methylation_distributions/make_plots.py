"""Create violin plots of CpG methylation in different regions."""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# Names to use for each sample
NAMES = {
    'FtubeAkapaBC'       : 'Kapa Rep. 1',
    'FtubeAkapaBCrep2'   : 'Kapa Rep. 2',
    'FtubeAneb'          : 'NEB Rep. 1',
    'FtubeAnebRep2'      : 'NEB Rep. 2',
    'FtubeApbat'         : 'PBAT Rep. 1',
    'FtubeAswift'        : 'Swift Rep. 1',
    'FtubeAswiftRep2'    : 'Swift Rep. 2',
    'FtubeAneb10ng'      : 'Low NEB Rep. 1',
    'FtubeAneb10ngRep2'  : 'Low NEB Rep. 2',
    'FtubeAswift10ng'    : 'Low Swift Rep. 1',
    'FtubeAswift10ngRep2': 'Low Swift Rep. 2',
    'FtubeBkapaBC'       : 'Kapa Rep. 1',
    'FtubeBkapaBCrep2'   : 'Kapa Rep. 2',
    'FtubeBneb'          : 'NEB Rep. 1',
    'FtubeBnebRep2'      : 'NEB Rep. 2',
    'FtubeBpbat'         : 'PBAT Rep. 1',
    'FtubeBswift'        : 'Swift Rep. 1',
    'FtubeBswiftRep2'    : 'Swift Rep. 2',
    'FtubeBneb10ng'      : 'Low NEB Rep. 1',
    'FtubeBneb10ngRep2'  : 'Low NEB Rep. 2',
    'FtubeBswift10ng'    : 'Low Swift Rep. 1',
    'FtubeBswift10ngRep2': 'Low Swift Rep. 2',
}

# Colors to use for each sample
COLOR = {
    10: '#D81B60', # Ftube{A,B}kapaBC
     9: '#D81B60', # Ftube{A,B}kapaBCrep2
     8: '#1E88E5', # Ftube{A,B}neb
     7: '#1E88E5', # Ftube{A,B}nebRep2
     6: '#A0522D', # Ftube{A,B}pbat
     5: '#004D40', # Ftube{A,B}swift
     4: '#004D40', # Ftube{A,B}swiftRep2
     3: '#1E88E5', # Ftube{A,B}neb10ng
     2: '#1E88E5', # Ftube{A,B}neb10ngRep2
     1: '#004D40', # Ftube{A,B}swift10ng
     0: '#004D40', # Ftube{A,B}swift10ngRep2
}

# Map sample A samples to order on plots
A_ORDER = {
    'FtubeAkapaBC': 1,
    'FtubeAneb': 3,
    'FtubeApbat': 5,
    'FtubeAswift': 6,
    'FtubeAkapaBCrep2': 2,
    'FtubeAnebRep2': 4,
    'FtubeAswiftRep2': 7,
    'FtubeAneb10ng': 8,
    'FtubeAswift10ng': 10,
    'FtubeAneb10ngRep2': 9,
    'FtubeAswift10ngRep2': 11
}

# Map sample B samples to order on plots
B_ORDER = {
    'FtubeBkapaBC': 1,
    'FtubeBneb': 3,
    'FtubeBpbat': 5,
    'FtubeBswift': 6,
    'FtubeBkapaBCrep2': 2,
    'FtubeBnebRep2': 4,
    'FtubeBswiftRep2': 7,
    'FtubeBneb10ng': 8,
    'FtubeBswift10ng': 10,
    'FtubeBneb10ngRep2': 9,
    'FtubeBswift10ngRep2': 11
}

def import_files(tag):
    """Import CpG islands, shores, shelves, and open seas BED files into a 
       DataFrame

    Inputs -
        tag - sample name to load
    Returns -
        tuple of lists (comb, islands, shores, shelves, open seas)
        Each list is the beta values from each respective DataFrame
    """
    names = ['chr','start','end','beta','covg']
    cols = [3]

    isld = pd.read_csv(tag+'.island.bed.gz', sep='\t', header=None, usecols=cols)
    isld.rename(columns={3: 'beta'}, inplace=True)
    isld['group'] = 'islands'

    shor = pd.read_csv(tag+'.shores.bed.gz', sep='\t', header=None, usecols=cols)
    shor.rename(columns={3: 'beta'}, inplace=True)
    shor['group'] = 'shores'

    shlv = pd.read_csv(tag+'.shelves.bed.gz', sep='\t', header=None, usecols=cols)
    shlv.rename(columns={3: 'beta'}, inplace=True)
    shlv['group'] = 'shelves'

    seas = pd.read_csv(tag+'.open_seas.bed.gz', sep='\t', header=None, usecols=cols)
    seas.rename(columns={3: 'beta'}, inplace=True)
    seas['group'] = 'open_seas'

    comb = pd.concat([isld, shor, shlv, seas])

    return (comb['beta'].values, isld['beta'].values, shor['beta'].values,
            shlv['beta'].values, seas['beta'].values)

def create_plot(data, title, xlab, ylab, ynames, figname):
    """Create violin plot of methylation data.

    Inputs -
        data    - list of data to include in violinplot
        title   - title of plot
        xlab    - x-axis label of plot
        ylab    - y-axis label of plot
        ynames  - names of ticks on y-axis
        figname - name of output file for plot
    """
    fig, ax = plt.subplots(figsize=(10,5))
    plt.tight_layout()

    parts = ax.violinplot(
        dataset = data,
        vert = False,
        showmeans = False,
        showmedians = False,
        showextrema = False
    )

    for idx, pc in enumerate(parts['bodies']):
        pc.set_facecolor(COLOR[idx])
        pc.set_edgecolor(COLOR[idx])
        pc.set_alpha(1)

    plt.title(title, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    plt.xticks(
        [i for i in np.arange(0, 1.2, 0.2)],
        ['{:.1f}'.format(i) for i in np.arange(0, 1.2, 0.2)],
        fontsize=12
    )
    plt.yticks([i for i in np.arange(1,len(ynames)+1)], ynames, fontsize=12)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def main():
    """Main function."""
    a_samps = [
        'FtubeAkapaBC', 'FtubeAneb', 'FtubeApbat', 'FtubeAswift',
        'FtubeAkapaBCrep2', 'FtubeAnebRep2', 'FtubeAswiftRep2',
        'FtubeAneb10ng', 'FtubeAswift10ng',
        'FtubeAneb10ngRep2', 'FtubeAswift10ngRep2'
    ]
    b_samps = [
        'FtubeBkapaBC', 'FtubeBneb', 'FtubeBpbat', 'FtubeBswift',
        'FtubeBkapaBCrep2', 'FtubeBnebRep2', 'FtubeBswiftRep2',
        'FtubeBneb10ng', 'FtubeBswift10ng',
        'FtubeBneb10ngRep2', 'FtubeBswift10ngRep2'
    ]

    a_data = []
    for samp in a_samps:
        c, i, r, v, s = import_files('seashore_bed_files/'+samp)
        a_data.append( (samp, c, i, r, v, s) )

    a_data = sorted(a_data, key=lambda d: A_ORDER[d[0]], reverse=True)
    a_name = [NAMES[i[0]] for i in a_data]
    a_comb = [i[1] for i in a_data]
    a_isld = [i[2] for i in a_data]
    a_shor = [i[3] for i in a_data]
    a_shlv = [i[4] for i in a_data]
    a_seas = [i[5] for i in a_data]

    create_plot(
        a_isld,
        'CpG Island Methylation Distribution: Sample A',
        'Methylation Level',
        '',
        a_name,
        'a_islands.pdf'
    )

    create_plot(
        a_shor,
        'CpG Shore Methylation Distribution: Sample A',
        'Methylation Level',
        '',
        a_name,
        'a_shores.pdf'
    )

    create_plot(
        a_shlv,
        'CpG Shelve Methylation Distribution: Sample A',
        'Methylation Level',
        '',
        a_name,
        'a_shelves.pdf'
    )

    create_plot(
        a_seas,
        'CpG Open Sea Methylation Distribution: Sample A',
        'Methylation Level',
        '',
        a_name,
        'a_opensea.pdf'
    )

    create_plot(
        a_comb,
        'All CpG Methylation Distribution: Sample A',
        'Methylation Level',
        '',
        a_name,
        'a_all.pdf'
    )

    b_data = []
    for samp in b_samps:
        c, i, r, v, s = import_files('seashore_bed_files/'+samp)
        b_data.append( (samp, c, i, r, v, s) )

    b_data = sorted(b_data, key=lambda d: B_ORDER[d[0]], reverse=True)
    b_name = [NAMES[i[0]] for i in b_data]
    b_comb = [i[1] for i in b_data]
    b_isld = [i[2] for i in b_data]
    b_shor = [i[3] for i in b_data]
    b_shlv = [i[4] for i in b_data]
    b_seas = [i[5] for i in b_data]

    create_plot(
        b_isld,
        'CpG Island Methylation Distribution: Sample B',
        'Methylation Level',
        '',
        b_name,
        'b_islands.pdf'
    )

    create_plot(
        b_shor,
        'CpG Shore Methylation Distribution: Sample B',
        'Methylation Level',
        '',
        b_name,
        'b_shores.pdf'
    )

    create_plot(
        b_shlv,
        'CpG Shelve Methylation Distribution: Sample B',
        'Methylation Level',
        '',
        b_name,
        'b_shelves.pdf'
    )

    create_plot(
        b_seas,
        'CpG Open Sea Methylation Distribution: Sample B',
        'Methylation Level',
        '',
        b_name,
        'b_opensea.pdf'
    )

    create_plot(
        b_comb,
        'All CpG Methylation Distribution: Sample B',
        'Methylation Level',
        '',
        b_name,
        'b_all.pdf'
    )

if __name__ == '__main__':
    main()
