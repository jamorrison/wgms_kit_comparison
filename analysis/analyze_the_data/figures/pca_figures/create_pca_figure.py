"""Script to create PCA plot to compare similarity of binned CpG methylation."""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.special import logit
from functools import reduce

KIT = {
    'FtubeAkapaBC'       : 'Kapa', 'FtubeAkapaBCrep2'   : 'Kapa',
    'FtubeAneb'          : 'NEB', 'FtubeAnebRep2'      : 'NEB',
    'FtubeApbat'         : 'PBAT',
    'FtubeAswift'        : 'Swift', 'FtubeAswiftRep2'    : 'Swift',
    'FtubeAneb10ng'      : 'NEB', 'FtubeAneb10ngRep2'  : 'NEB',
    'FtubeAswift10ng'    : 'Swift', 'FtubeAswift10ngRep2': 'Swift',
    'FtubeBkapaBC'       : 'Kapa', 'FtubeBkapaBCrep2'   : 'Kapa',
    'FtubeBneb'          : 'NEB', 'FtubeBnebRep2'      : 'NEB',
    'FtubeBpbat'         : 'PBAT',
    'FtubeBswift'        : 'Swift', 'FtubeBswiftRep2'    : 'Swift',
    'FtubeBneb10ng'      : 'NEB', 'FtubeBneb10ngRep2'  : 'NEB',
    'FtubeBswift10ng'    : 'Swift', 'FtubeBswift10ngRep2': 'Swift'
}

SAM = {
    'FtubeAkapaBC'       : 'A', 'FtubeAkapaBCrep2'   : 'A',
    'FtubeAneb'          : 'A', 'FtubeAnebRep2'      : 'A',
    'FtubeApbat'         : 'A',
    'FtubeAswift'        : 'A', 'FtubeAswiftRep2'    : 'A',
    'FtubeAneb10ng'      : 'A', 'FtubeAneb10ngRep2'  : 'A',
    'FtubeAswift10ng'    : 'A', 'FtubeAswift10ngRep2': 'A',
    'FtubeBkapaBC'       : 'B', 'FtubeBkapaBCrep2'   : 'B',
    'FtubeBneb'          : 'B', 'FtubeBnebRep2'      : 'B',
    'FtubeBpbat'         : 'B',
    'FtubeBswift'        : 'B', 'FtubeBswiftRep2'    : 'B',
    'FtubeBneb10ng'      : 'B', 'FtubeBneb10ngRep2'  : 'B',
    'FtubeBswift10ng'    : 'B', 'FtubeBswift10ngRep2': 'B'
}

REP = {
    'FtubeAkapaBC'       : '1', 'FtubeAkapaBCrep2'   : '2',
    'FtubeAneb'          : '1', 'FtubeAnebRep2'      : '2',
    'FtubeApbat'         : '1',
    'FtubeAswift'        : '1', 'FtubeAswiftRep2'    : '2',
    'FtubeAneb10ng'      : '1', 'FtubeAneb10ngRep2'  : '2',
    'FtubeAswift10ng'    : '1', 'FtubeAswift10ngRep2': '2',
    'FtubeBkapaBC'       : '1', 'FtubeBkapaBCrep2'   : '2',
    'FtubeBneb'          : '1', 'FtubeBnebRep2'      : '2',
    'FtubeBpbat'         : '1',
    'FtubeBswift'        : '1', 'FtubeBswiftRep2'    : '2',
    'FtubeBneb10ng'      : '1', 'FtubeBneb10ngRep2'  : '2',
    'FtubeBswift10ng'    : '1', 'FtubeBswift10ngRep2': '2'
}

def beta_to_m(betas, covgs, k):
    """Transform beta values into m values.

    Inputs -
        betas - pd.Series of beta values
        covgs - pd.Series of covg values
        k     - number of pseudoreads for smoothing
    Returns
        pd.Series of m values
    """
    b = list(betas)
    c = list(covgs)

    s = []
    for i in range(len(c)):
        m = (c[i] * b[i])
        u = (c[i] - m)
        s.append((m+k) / ((m+k) + (u+k)))
    out = logit(s)

    return pd.Series(out)

def make_plot(data, var, keys, cols, title, xlab, ylab, figname):
    """Create plot for principal components of PCA.

    Inputs -
        data    - DataFrame with columns pc1, pc2, var
        var     - column name of variable to keep
        keys    - keys from var column that correspond to colors in cols list
        cols    - colors to match with keys list
        title   - title of figure
        xlab    - x-axis label
        ylab    - y-axis label
        figname - name of output file
    Returns -
        Nothing, plot is saved to disk
    """
    fig, ax = plt.subplots(figsize=(5,5))
    plt.tight_layout()

    if len(keys) != len(cols):
        print('[make_plot] ERROR: Mismatch in number of keys and colors.')
        return None

    # Add scatter plot to figure for each set of data in keys
    for key, col in zip(keys, cols):
        to_keep = data[var] == key
        if True not in list(to_keep):
            continue
        
        ax.scatter(data.loc[to_keep, 'pc1'], data.loc[to_keep, 'pc2'],
                   c=col, s=25)

    ax.legend(keys, loc='upper left')

    plt.title(title, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def make_combined_plot(data, title, xlab, ylab, figname):
    """Create plot for principal components of PCA with all variables shown.

    Inputs -
        data    - DataFrame with columns pc1, pc2, kit, sample, and replicate
        title   - title of figure
        xlab    - x-axis label
        ylab    - y-axis label
        figname - name of output file
    Returns -
        Nothing, plot is saved to disk
    """
    keys = {
         1: ('A', '1', 'Kapa'),   2: ('A', '2', 'Kapa'),
         3: ('B', '1', 'Kapa'),   4: ('B', '2', 'Kapa'),
         5: ('A', '1', 'NEB'),    6: ('A', '2', 'NEB'),
         7: ('B', '1', 'NEB'),    8: ('B', '2', 'NEB'),
         9: ('A', '1', 'PBAT'),  10: ('B', '1', 'PBAT'),
        11: ('A', '1', 'Swift'), 12: ('A', '2', 'Swift'),
        13: ('B', '1', 'Swift'), 14: ('B', '2', 'Swift'),
    }
    cols = {
         1: ('o', 'none', '#D81B60'),  2: ('o', '#D81B60', '#D81B60'),
         3: ('s', 'none', '#D81B60'),  4: ('s', '#D81B60', '#D81B60'),
         5: ('o', 'none', '#1E88E5'),  6: ('o', '#1E88E5', '#1E88E5'),
         7: ('s', 'none', '#1E88E5'),  8: ('s', '#1E88E5', '#1E88E5'),
         9: ('o', 'none', '#A0522D'), 10: ('s', 'none'   , '#A0522D'),
        11: ('o', 'none', '#004D40'), 12: ('o', '#004D40', '#004D40'),
        13: ('s', 'none', '#004D40'), 14: ('s', '#004D40', '#004D40'),
    }
    fig, ax = plt.subplots(figsize=(5,5))
    plt.tight_layout()

    # Create legend
    plt.plot(-300, 300, 'o', color='black', markersize=8, label='Samp. A')
    plt.plot(-300, 300, 's', color='black', markersize=8, label='Samp. B')
    plt.plot(-300, 300, 'D', color='black', fillstyle='none', markersize=8, label='Rep. 1')
    plt.plot(-300, 300, 'D', color='black', markersize=8, label='Rep. 2')
    plt.plot(-300, 300, 'D', color='#D81B60', markersize=8, label='Kapa')
    plt.plot(-300, 300, 'D', color='#1E88E5', markersize=8, label='NEB')
    plt.plot(-300, 300, 'D', color='#A0522D', markersize=8, label='PBAT')
    plt.plot(-300, 300, 'D', color='#004D40', markersize=8, label='Swift')

    # Add scatter plot to figure for each set of data in keys
    for key, col in zip(keys.values(), cols.values()):
        to_keep = (data['sample'] == key[0]) & (data['replicate'] == key[1]) & (data['kit'] == key[2])
        if True not in list(to_keep):
            continue
        
        ax.scatter(data.loc[to_keep, 'pc1'], data.loc[to_keep, 'pc2'],
                    marker=col[0], facecolors=col[1], edgecolors=col[2], s=25)

    ax.legend(ncol=4, bbox_to_anchor=(0.5, 1), frameon=False,
              loc='lower center', fontsize='large')
    plt.xlim(-250, 250)
    plt.ylim(-175, 175)

    plt.title(title, pad=50, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def main():
    """Do the bulk of the PCA generation."""
    dirloc = '../../subsampling/'
    filtag = '.subsampled.cg.sorted.mergecg.100kb_meth_avg.bed.gz'

    samps = [
        'FtubeAkapaBC'       , 'FtubeAkapaBCrep2'   ,
        'FtubeAneb'          , 'FtubeAnebRep2'      ,
        'FtubeApbat'         ,
        'FtubeAswift'        , 'FtubeAswiftRep2'    ,
        'FtubeAneb10ng'      , 'FtubeAneb10ngRep2'  ,
        'FtubeAswift10ng'    , 'FtubeAswift10ngRep2',
        'FtubeBkapaBC'       , 'FtubeBkapaBCrep2'   ,
        'FtubeBneb'          , 'FtubeBnebRep2'      ,
        'FtubeBpbat'         ,
        'FtubeBswift'        , 'FtubeBswiftRep2'    ,
        'FtubeBneb10ng'      , 'FtubeBneb10ngRep2'  ,
        'FtubeBswift10ng'    , 'FtubeBswift10ngRep2'
    ]
    dfs = []
    for samp in samps:
        cols = ['chr', 'start', 'end', 'beta', 'covg']

        df = pd.read_csv(dirloc+samp+filtag, sep='\t', names=cols, na_values='.')

        # Transform beta values into m-values uses logit transform
        # Makes beta distribution of values into a more gaussian distribution
        df[samp] = beta_to_m(df['beta'], df['covg'], 0.1)
        dfs.append(df.drop(['beta', 'covg'], axis=1))

    df_na = reduce(lambda x, y: pd.merge(x, y, on=['chr', 'start', 'end']), dfs)
    df_al = df_na.dropna(axis=0, how='any')

    df = df_al.drop(['chr', 'start', 'end'], axis=1).transpose()
    kits = [KIT[i] for i in list(df.index)]
    sams = [SAM[i] for i in list(df.index)]
    reps = [REP[i] for i in list(df.index)]

    # Standardize values
    x = StandardScaler().fit_transform(df)
    std = pd.DataFrame(data=x, columns=list(df.columns))

    # Create PCA
    pca = PCA(n_components=2)
    principals = pca.fit_transform(x)
    print(pca.explained_variance_ratio_)
    pcs = pd.DataFrame(data=principals, columns=['pc1', 'pc2'])
    pcs['kit'] = kits
    pcs['sample'] = sams
    pcs['replicate'] = reps

    # Make figures
    make_plot(
        pcs,
        'kit',
        ['Kapa', 'NEB', 'PBAT', 'Swift'],
        ['#D81B60', '#1E88E5', '#A0522D', '#004D40'],
        '2-component PCA: Protocol',
        'Principal Component 1',
        'Principal Component 2',
        'pca_protocol.pdf'
    )

    make_plot(
        pcs,
        'sample',
        ['A', 'B'],
        ['#981e32', '#5e6a71'],
        '2-component PCA: Sample',
        'Principal Component 1',
        'Principal Component 2',
        'pca_sample.pdf'
    )

    make_plot(
        pcs,
        'replicate',
        ['1', '2'],
        ['#005596', '#3fa294'],
        '2-component PCA: Tech. Rep.',
        'Principal Component 1',
        'Principal Component 2',
        'pca_replicate.pdf'
    )

    make_combined_plot(
        pcs,
        '2-component PCA',
        'Principal Component 1',
        'Principal Component 2',
        'pca_combined.pdf'
    )

if __name__ == '__main__':
    main()
