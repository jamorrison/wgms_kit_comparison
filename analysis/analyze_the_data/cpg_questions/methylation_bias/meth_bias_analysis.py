import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, gaussian_kde
import sys

def make_diff_avg_plot(x_vals, y_vals, title, xlab, ylab, figname):
    """Find the Spearman R correlation value and create scatter plot of beta
       values.

    Inputs: x_vals        - list of values on x-axis
            y_vals        - list of values on y-axis
            title         - title of figure
            xlab          - name of values on x-axis
            ylab          - name of values on y-axis
            figname       - name of file to save figure as

    Returns: 
    """
    fig, ax = plt.subplots(figsize=(5,5))
    plt.tight_layout()

    ax.plot(x_vals, y_vals, 'k.')

    tick_loc = [i for i in np.arange(-1, 1.2, 0.2)]
    tick_lab = []
    for i in np.arange(-1, 1.2, 0.2):
        if abs(i) > 0.01:
            tick_lab.append('{:.1f}'.format(i))
        else:
            tick_lab.append('{:.1f}'.format(abs(i)))

    plt.xlim(-1.05, 1.05)
    plt.ylim(-1.05, 1.05)
    plt.xticks(tick_loc, tick_lab, fontsize=12)
    plt.yticks(tick_loc, tick_lab, fontsize=12)

    plt.title(title, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def make_diff_avg_pt_density_plot(x_vals, y_vals, title, xlab, ylab, figname):
    """Find the Spearman R correlation value and create scatter plot of beta
       values.

    Inputs: x_vals        - list of values on x-axis
            y_vals        - list of values on y-axis
            title         - title of figure
            xlab          - name of values on x-axis
            ylab          - name of values on y-axis
            figname       - name of file to save figure as

    Returns: 
    """
    xs = np.array(x_vals)
    ys = np.array(y_vals)
    nbins = 100

    fig, ax = plt.subplots(figsize=(5,5))
    plt.tight_layout()

    try:
        k = gaussian_kde(np.vstack([xs, ys]))
    except np.linalg.LinAlgError:
        xs[0] = xs[0] + 0.0000001
        k = gaussian_kde(np.vstack([xs, ys]))
    xi, yi = np.mgrid[xs.min():xs.max():nbins*1j, ys.min():ys.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        
    im = ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.PuBu_r)
    ax.contour(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.viridis)

    plt.xlim(-1.05, 1.05)
    plt.ylim(-1.05, 1.05)
    plt.xticks([i for i in np.arange(-1, 1.2, 0.2)], ['{:.1f}'.format(i) for i in np.arange(-1, 1.2, 0.2)])
    plt.yticks([i for i in np.arange(-1, 1.2, 0.2)], ['{:.1f}'.format(i) for i in np.arange(-1, 1.2, 0.2)])

    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.colorbar(im, ax=ax)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def find_big_diff_cpgs(neb_dic, swi_dic, outfile):
    """Find CpGs that have a large difference between NEB and Swift kits.

    Inputs: neb_dic - dictionary of NEB CpG beta values
            swi_dic - dictionary of Swift CpG beta values
            outfile - name of output file for large beta value difference CpGs
    
    Returns: tuple (all_diffs, lrg_diffs, df)
             all_diffs - dictionary containing all CpG beta value differences
                         between NEB and Swift kits
             lrg_diffs - dictionary containing only large CpG beta value
                         differences (abs(diff) > 0.5) between NEB and Swift
                         kits
             df        - data frame containing CpGs with large beta value
                         differences
    """
    all_diffs = {}
    lrg_diffs = {}
    for dic_n in neb_dic:
        for dic_s in swi_dic:
            key = dic_n['tag'] + '_' + dic_s['tag']

            all_diffs[key] = {}
            lrg_diffs[key] = {}

            neb = dic_n['data']
            swi = dic_s['data']

            neb = neb[neb.covg > 20]
            swi = swi[swi.covg > 20]

            merged = pd.merge(neb, swi, on=['chr', 'start', 'end'], suffixes=['_'+dic_n['tag'], '_'+dic_s['tag']])

            merged['diff'] = merged['beta_'+dic_n['tag']] - merged['beta_'+dic_s['tag']]
            merged['avg'] = (merged['beta_'+dic_n['tag']] + merged['beta_'+dic_s['tag']]) / 2
            larged = merged[merged['diff'].abs() > 0.5]

            all_diffs[key]['data'] = merged[['chr', 'start', 'end', 'beta_'+dic_n['tag'], 'beta_'+dic_s['tag'], 'diff', 'avg']]
            all_diffs[key]['tag'] = key
            all_diffs[key]['label'] = ' '.join([dic_n['lib'], dic_n['title']]) + ' - ' + ' '.join([dic_s['lib'], dic_s['title']])
            all_diffs[key]['n_tag'] = dic_n['tag']
            all_diffs[key]['s_tag'] = dic_s['tag']

            lrg_diffs[key]['data'] = larged[['chr', 'start', 'end', 'beta_'+dic_n['tag'], 'beta_'+dic_s['tag'], 'diff', 'avg']]
            lrg_diffs[key]['tag'] = key
            lrg_diffs[key]['label'] = ' '.join([dic_n['lib'], dic_n['title']]) + ' - ' + ' '.join([dic_s['lib'], dic_s['title']])
            lrg_diffs[key]['n_tag'] = dic_n['tag']
            lrg_diffs[key]['s_tag'] = dic_s['tag']

    possible_tags = [
        'beta_an1_hi', 'beta_an2_hi', 'beta_an1_lo', 'beta_an2_lo',
        'beta_bn1_hi', 'beta_bn2_hi', 'beta_bn1_lo', 'beta_bn2_lo',
        'beta_as1_hi', 'beta_as2_hi', 'beta_as1_lo', 'beta_as2_lo',
        'beta_bs1_hi', 'beta_bs2_hi', 'beta_bs1_lo', 'beta_bs2_lo'
    ]
    frames = []
    for key1, dic1 in lrg_diffs.items():
        for key2, dic2 in lrg_diffs.items():
            if key1 == key2:
                continue

            merged = pd.merge(dic1['data'], dic2['data'], on=['chr', 'start', 'end'], suffixes=['_'+key1, '_'+key2])

            if dic1['n_tag'] == dic2['n_tag']:
                red = merged[
                    ['chr', 'start', 'end', 'beta_'+dic1['n_tag']+'_'+key1, 'beta_'+dic1['s_tag'], 'beta_'+dic2['s_tag']]
                ]
                red = red.rename(columns={'beta_'+dic1['n_tag']+'_'+key1: 'beta_'+dic1['n_tag']})
            elif dic1['s_tag'] == dic2['s_tag']:
                red = merged[
                    ['chr', 'start', 'end', 'beta_'+dic1['s_tag']+'_'+key1, 'beta_'+dic1['n_tag'], 'beta_'+dic2['n_tag']]
                ]
                red = red.rename(columns={'beta_'+dic1['s_tag']+'_'+key1: 'beta_'+dic1['s_tag']})
            else:
                red = merged[
                    ['chr', 'start', 'end', 'beta_'+dic1['n_tag'],  'beta_'+dic2['n_tag'],'beta_'+dic1['s_tag'], 'beta_'+dic2['s_tag']]
                ]

            for p_tag in possible_tags:
                if not p_tag in list(red.columns):
                    red = red.assign(p_tag = 'NA')
                    red = red.rename(columns={'p_tag': p_tag})

            frames.append(red)

    df = pd.concat(frames, ignore_index=True)
    df = df.sort_values(by=['chr', 'start'], ignore_index=True)
    df = df.drop_duplicates(ignore_index=True)

    df.to_csv(outfile, sep='\t', index=False, float_format='%.4f')

    return all_diffs, lrg_diffs, df

def main():
    """Run methylation bias analysis."""
    col_names = ['chr', 'start', 'end', 'beta', 'covg', 'context']

    dirloc = '2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/align/'
    filtag = '.cg.sorted.mergecg.bed.gz'

    df_01 = pd.read_csv(dirloc + 'FtubeAneb' + filtag, sep='\t', names=col_names, na_values='.')
    df_02 = pd.read_csv(dirloc + 'FtubeAswift' + filtag, sep='\t', names=col_names, na_values='.')
    df_03 = pd.read_csv(dirloc + 'FtubeAnebRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_04 = pd.read_csv(dirloc + 'FtubeAswiftRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_05 = pd.read_csv(dirloc + 'FtubeBneb' + filtag, sep='\t', names=col_names, na_values='.')
    df_06 = pd.read_csv(dirloc + 'FtubeBswift' + filtag, sep='\t', names=col_names, na_values='.')
    df_07 = pd.read_csv(dirloc + 'FtubeBnebRep2' + filtag, sep='\t', names=col_names, na_values='.')
    df_08 = pd.read_csv(dirloc + 'FtubeBswiftRep2' + filtag, sep='\t', names=col_names, na_values='.')

    neb_sets_a = [
        {'data': df_01, 'tag': 'an1_hi', 'title': 'Sample A Rep. 1', 'lib': 'NEB'},
        {'data': df_03, 'tag': 'an2_hi', 'title': 'Sample A Rep. 2', 'lib': 'NEB'},
    ]
    swi_sets_a = [
        {'data': df_02, 'tag': 'as1_hi', 'title': 'Sample A Rep. 1', 'lib': 'Swift'},
        {'data': df_04, 'tag': 'as2_hi', 'title': 'Sample A Rep. 2', 'lib': 'Swift'},
    ]
    neb_sets_b = [
        {'data': df_05, 'tag': 'bn1_hi', 'title': 'Sample B Rep. 1', 'lib': 'NEB'},
        {'data': df_07, 'tag': 'bn2_hi', 'title': 'Sample B Rep. 2', 'lib': 'NEB'},
    ]
    swi_sets_b = [
        {'data': df_06, 'tag': 'bs1_hi', 'title': 'Sample B Rep. 1', 'lib': 'Swift'},
        {'data': df_08, 'tag': 'bs2_hi', 'title': 'Sample B Rep. 2', 'lib': 'Swift'},
    ]

    all_diffs_a, lrg_diffs_a, df_a = find_big_diff_cpgs(neb_sets_a, swi_sets_a, 'test_big_diff_sampleA.bed')
    all_diffs_b, lrg_diffs_b, df_b = find_big_diff_cpgs(neb_sets_b, swi_sets_b, 'test_big_diff_sampleB.bed')

    merged_a = pd.merge(all_diffs_a['an1_hi_as1_hi']['data'],
                        all_diffs_a['an2_hi_as2_hi']['data'],
                        on=['chr', 'start', 'end'])
    merged_b = pd.merge(all_diffs_b['bn1_hi_bs1_hi']['data'],
                        all_diffs_b['bn2_hi_bs2_hi']['data'],
                        on=['chr', 'start', 'end'])

    make_diff_avg_plot(
        merged_a['diff_x'],
        merged_a['diff_y'],
        'Sample A Beta Value Difference',
        'NEB Rep. 1 - Swift Rep. 1',
        'NEB Rep. 2 - Swift Rep. 2',
        'plots/sample_a_rep1_vs_rep2_diff_points.png'
    )
    make_diff_avg_pt_density_plot(
        merged_a['diff_x'],
        merged_a['diff_y'],
        'Sample A Beta Value Difference',
        'NEB Rep. 1 - Swift Rep. 1',
        'NEB Rep. 2 - Swift Rep. 2',
        'plots/sample_a_rep1_vs_rep2_diff_density.png'
    )
    make_diff_avg_plot(
        merged_b['diff_x'],
        merged_b['diff_y'],
        'Sample B Beta Value Difference',
        'NEB Rep. 1 - Swift Rep. 1',
        'NEB Rep. 2 - Swift Rep. 2',
        'plots/sample_b_rep1_vs_rep2_diff_points.png'
    )
    make_diff_avg_pt_density_plot(
        merged_b['diff_x'],
        merged_b['diff_y'],
        'Sample B Beta Value Difference',
        'NEB Rep. 1 - Swift Rep. 1',
        'NEB Rep. 2 - Swift Rep. 2',
        'plots/sample_b_rep1_vs_rep2_diff_density.png'
    )

if __name__ == '__main__':
    main()
