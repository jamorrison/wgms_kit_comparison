"""Script to create PCA plot to compare similarity of binned CpG methylation."""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
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
        cols = ['chr', 'start', 'end', samp]

        df = pd.read_csv(dirloc+samp+filtag, sep='\t', names=cols, na_values='.')
        dfs.append(df)

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

if __name__ == '__main__':
    main()
