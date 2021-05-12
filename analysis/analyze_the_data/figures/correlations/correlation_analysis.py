import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, gaussian_kde
import sys

def find_correlation_and_plot(x_vals, y_vals, title, xlab, ylab, figname, print_message=False, create_plots=True):
    """Find the Spearman R correlation value and create scatter plot of beta
       values.

    Inputs: x_vals        - list of values on x-axis
            y_vals        - list of values on y-axis
            title         - title of figure
            xlab          - name of values on x-axis
            ylab          - name of values on y-axis
            figname       - name of file to save figure as
            print_message - whether to print out a nice message about correlation
            create_plots  - whether to generate the correlation plot

    Returns: tuple (n_cpgs, coef, p)
    """
    # Create arrays with NaN values removed
    xs = []
    ys = []
    for i in range(len(x_vals)):
        if not np.isnan(x_vals[i]) and not np.isnan(y_vals[i]):
            xs.append(x_vals[i])
            ys.append(y_vals[i])
    xs = np.array(xs)
    ys = np.array(ys)

    nbins = 100
    n_cpgs = len(xs)
    coef, p = spearmanr(xs, ys, nan_policy='omit')

    if print_message:
        print('The coefficient between {} and {} ='.format(xlab, ylab), end=' ')
        print('{:0.3f} (with p-value = {:0.3f} and # CpGs = {:,})'.format(coef, p, n_cpgs))

    if create_plots:
        # Correlation figure
        fig, ax = plt.subplots(figsize=(5,5))
        plt.tight_layout()

        try:
            k = gaussian_kde(np.vstack([xs, ys]))
        except np.linalg.LinAlgError:
            xs[0] = xs[0] + 0.000001
            k = gaussian_kde(np.vstack([xs, ys]))
        xi, yi = np.mgrid[xs.min():xs.max():nbins*1j, ys.min():ys.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        
        im = ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.PuBu_r,
                          norm=colors.LogNorm(vmin=0.001, vmax=zi.max()))
        ax.contour(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.viridis)

        plt.text(
            0.5, 1.075, r'# Bins = {:,} ; $r_s$ = {:.3f}'.format(n_cpgs, coef),
            ha='center', va='center', size=16
        )

        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.15)

        plt.xticks([i for i in np.arange(0, 1.2, 0.2)], ['{:.1f}'.format(i) for i in np.arange(0, 1.2, 0.2)], fontsize=18)
        plt.yticks([i for i in np.arange(0, 1.2, 0.2)], ['{:.1f}'.format(i) for i in np.arange(0, 1.2, 0.2)], fontsize=18)

        plt.title(title, fontsize=24)
        plt.xlabel(xlab, fontsize=20)
        plt.ylabel(ylab, fontsize=20)
        plt.colorbar(im, ax=ax)

        plt.savefig(figname, bbox_inches='tight')
        plt.close('all')

        # Collapse correlation plots into histograms
        fig, ax = plt.subplots(figsize=(5,5))
        plt.tight_layout()

        n1, b1, p1 = plt.hist(xs, bins=50, range=(0,1), density=True,
                              color='red', label=xlab, alpha=0.5)
                              #color='#005596', label=xlab, alpha=0.5)
        n2, b2, p2 = plt.hist(ys, bins=50, range=(0,1), density=True,
                              color='blue', label=ylab, alpha=0.5)
                              #color='#3fa294', label=ylab, alpha=0.5)

        ax.legend(ncol=1, loc='upper left', fontsize=20)

        plt.title(title, fontsize=24)
        plt.xlabel('Methylation Level', fontsize=20)
        plt.ylabel('# CpGs / Total # CpGs / 0.02', fontsize=20)

        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05*max(max(n1), max(n2)))
        plt.xticks(
            [i for i in np.arange(0, 1.2, 0.2)],
            ['{:.1f}'.format(i) for i in np.arange(0, 1.2, 0.2)],
            fontsize=18
        )
        plt.yticks(
            [i for i in np.arange(0, 1.05*max(max(n1),max(n2)), 1)],
            ['{:.1f}'.format(i) for i in np.arange(0, 1.05*max(max(n1),max(n2)), 1)],
            fontsize=18
        )

        plt.savefig('hist_'+figname, bbox_inches='tight')
        plt.close('all')

    return (n_cpgs, coef, p)

def main():
    """Do the bulk of the correlations analysis."""
    col_names = ['chr', 'start', 'end', 'beta_avg', 'covg_avg']

    dirloc = '../../subsampling/'
    filtag = '.subsampled.cg.sorted.mergecg.100kb_meth_avg.bed.gz'

    print('Loading data')
    df_01 = pd.read_csv(dirloc + 'FtubeAkapaBC' + filtag, sep='\t', names=col_names, na_values='.')
    df_02 = pd.read_csv(dirloc + 'FtubeAkapaBCrep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_03 = pd.read_csv(dirloc + 'FtubeAneb10ngRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_04 = pd.read_csv(dirloc + 'FtubeAneb10ng' + filtag, sep='\t', names=col_names, na_values='.')
    df_05 = pd.read_csv(dirloc + 'FtubeAnebRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_06 = pd.read_csv(dirloc + 'FtubeAneb' + filtag, sep='\t', names=col_names, na_values='.')
    df_07 = pd.read_csv(dirloc + 'FtubeApbat' + filtag, sep='\t', names=col_names, na_values='.')
    df_08 = pd.read_csv(dirloc + 'FtubeAswift10ngRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_09 = pd.read_csv(dirloc + 'FtubeAswift10ng' + filtag, sep='\t', names=col_names, na_values='.')
    df_10 = pd.read_csv(dirloc + 'FtubeAswiftRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_11 = pd.read_csv(dirloc + 'FtubeAswift' + filtag, sep='\t', names=col_names, na_values='.')
    df_12 = pd.read_csv(dirloc + 'FtubeBkapaBC' + filtag, sep='\t', names=col_names, na_values='.')
    df_13 = pd.read_csv(dirloc + 'FtubeBkapaBCrep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_14 = pd.read_csv(dirloc + 'FtubeBneb10ngRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_15 = pd.read_csv(dirloc + 'FtubeBneb10ng' + filtag, sep='\t', names=col_names, na_values='.')
    df_16 = pd.read_csv(dirloc + 'FtubeBnebRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_17 = pd.read_csv(dirloc + 'FtubeBneb' + filtag, sep='\t', names=col_names, na_values='.')
    df_18 = pd.read_csv(dirloc + 'FtubeBpbat' + filtag, sep='\t', names=col_names, na_values='.')
    df_19 = pd.read_csv(dirloc + 'FtubeBswift10ngRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_20 = pd.read_csv(dirloc + 'FtubeBswift10ng' + filtag, sep='\t', names=col_names, na_values='.')
    df_21 = pd.read_csv(dirloc + 'FtubeBswiftRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_22 = pd.read_csv(dirloc + 'FtubeBswift' + filtag, sep='\t', names=col_names, na_values='.')

    data_sets = [
        {'data': df_01, 'sample': 'FtubeAkapaBC'       , 'tag': 'ak1_hi', 'title': 'Kapa'       , 'axis_label': 'Replicate 1'},
        {'data': df_02, 'sample': 'FtubeAkapaBCrep2'   , 'tag': 'ak2_hi', 'title': 'Kapa'       , 'axis_label': 'Replicate 2'},
        {'data': df_06, 'sample': 'FtubeAneb'          , 'tag': 'an1_hi', 'title': 'NEB'        , 'axis_label': 'Replicate 1'},
        {'data': df_05, 'sample': 'FtubeAnebRep2'      , 'tag': 'an2_hi', 'title': 'NEB'        , 'axis_label': 'Replicate 2'},
        {'data': df_07, 'sample': 'FtubeApbat'         , 'tag': 'ap1_hi', 'title': 'PBAT'       , 'axis_label': 'Replicate 1'},
        {'data': df_11, 'sample': 'FtubeAswift'        , 'tag': 'as1_hi', 'title': 'Swift'      , 'axis_label': 'Replicate 1'},
        {'data': df_10, 'sample': 'FtubeAswiftRep2'    , 'tag': 'as2_hi', 'title': 'Swift'      , 'axis_label': 'Replicate 2'},
        {'data': df_04, 'sample': 'FtubeAneb10ng'      , 'tag': 'an1_lo', 'title': 'Low NEB'    , 'axis_label': 'Replicate 1'},
        {'data': df_03, 'sample': 'FtubeAneb10ngRep2'  , 'tag': 'an2_lo', 'title': 'Low NEB'    , 'axis_label': 'Replicate 2'},
        {'data': df_09, 'sample': 'FtubeAswift10ng'    , 'tag': 'as1_lo', 'title': 'Low Swift'  , 'axis_label': 'Replicate 1'},
        {'data': df_08, 'sample': 'FtubeAswift10ngRep2', 'tag': 'as2_lo', 'title': 'Low Swift'  , 'axis_label': 'Replicate 2'},
        {'data': df_12, 'sample': 'FtubeBkapaBC'       , 'tag': 'bk1_hi', 'title': 'Kapa'       , 'axis_label': 'Replicate 1'},
        {'data': df_13, 'sample': 'FtubeBkapaBCrep2'   , 'tag': 'bk2_hi', 'title': 'Kapa'       , 'axis_label': 'Replicate 2'},
        {'data': df_17, 'sample': 'FtubeBneb'          , 'tag': 'bn1_hi', 'title': 'NEB'        , 'axis_label': 'Replicate 1'},
        {'data': df_16, 'sample': 'FtubeBnebRep2'      , 'tag': 'bn2_hi', 'title': 'NEB'        , 'axis_label': 'Replicate 2'},
        {'data': df_18, 'sample': 'FtubeBpbat'         , 'tag': 'bp1_hi', 'title': 'PBAT'       , 'axis_label': 'Replicate 1'},
        {'data': df_22, 'sample': 'FtubeBswift'        , 'tag': 'bs1_hi', 'title': 'Swift'      , 'axis_label': 'Replicate 1'},
        {'data': df_21, 'sample': 'FtubeBswiftRep2'    , 'tag': 'bs2_hi', 'title': 'Swift'      , 'axis_label': 'Replicate 2'},
        {'data': df_15, 'sample': 'FtubeBneb10ng'      , 'tag': 'bn1_lo', 'title': 'Low NEB'    , 'axis_label': 'Replicate 1'},
        {'data': df_14, 'sample': 'FtubeBneb10ngRep2'  , 'tag': 'bn2_lo', 'title': 'Low NEB'    , 'axis_label': 'Replicate 2'},
        {'data': df_20, 'sample': 'FtubeBswift10ng'    , 'tag': 'bs1_lo', 'title': 'Low Swift'  , 'axis_label': 'Replicate 1'},
        {'data': df_19, 'sample': 'FtubeBswift10ngRep2', 'tag': 'bs2_lo', 'title': 'Low Swift'  , 'axis_label': 'Replicate 2'}
    ]

    with open('correlation_values.tsv', 'w') as f:
        for d1 in range(len(data_sets)):
            for d2 in range(d1, len(data_sets)):
                n_cpgs = 0
                coef = 0
                p = 0
                print('Processing {} {}'.format(data_sets[d1]['tag'], data_sets[d2]['tag']))
                if 'FtubeA' in data_sets[d1]['sample'] and 'FtubeA' in data_sets[d2]['sample']:
                    n_cpgs, coef, p = find_correlation_and_plot(
                        data_sets[d1]['data']['beta_avg'],
                        data_sets[d2]['data']['beta_avg'],
                        data_sets[d1]['title'] + ' Sample A Beta Values',
                        data_sets[d1]['axis_label'],
                        data_sets[d2]['axis_label'],
                        data_sets[d1]['tag'] + '_' + data_sets[d2]['tag'] + '.png',
                        print_message=False,
                        create_plots=True
                    )
                else:
                    n_cpgs, coef, p = find_correlation_and_plot(
                        data_sets[d1]['data']['beta_avg'],
                        data_sets[d2]['data']['beta_avg'],
                        data_sets[d1]['title'],
                        data_sets[d1]['axis_label'],
                        data_sets[d2]['axis_label'],
                        data_sets[d1]['tag'] + '_' + data_sets[d2]['tag'] + '.png',
                        print_message=False,
                        create_plots=False
                    )

                f.write('{}\t{}\t{}\t{}\t{}\n'.format(data_sets[d1]['tag'], data_sets[d2]['tag'], n_cpgs, coef, p))

if __name__ == '__main__':
    main()
