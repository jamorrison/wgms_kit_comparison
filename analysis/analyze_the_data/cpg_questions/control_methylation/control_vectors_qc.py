"""Create plots to show methylation for mitochondrial and control vector DNA."""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os

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

def import_file(fname, control=True, select_chr=None):
    """Read file into DataFrame and add some extra columns.

    Inputs -
        fname      - filename to process
        control    - whether file is for control vectors (True) or
                     human (False) [default: True]
        select_chr - select out chromosome from DataFrame
    Returns -
        DataFrame of data from file with some added columns
    """
    df = pd.read_csv(
        fname,
        sep='\t',
        header=None,
        names=['chr','start','end','beta','covg','context']
    )

    # Add sample column
    samp = os.path.basename(fname)
    samp = samp.replace('.cg.sorted.mergecg.bed.gz', '')
    samp = samp.replace('.controls', '')
    kit, bio, rep = extract_sample_info(samp)
    df['sample'] = kit + ' Samp. ' + bio
    df['Replicate'] = 'Rep. ' + rep

    # Include column of which control vector the CpG is from, only for control
    # vectors
    if (control == True):
        df['vector'] = np.where(df['chr'] == 'J02459.1', 'lambda', 'pUC19')

    # Select only CpGs that occur on select_chr chromosome
    if select_chr != None:
        out = df[df.chr == select_chr]
    else:
        out = df

    return out

def create_plot_matplotlib(data, title, xlab, ylab, ynames, figname):
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

    plt.title(title, fontsize=24)
    plt.xlabel(xlab, fontsize=20)
    plt.ylabel(ylab, fontsize=20)

    plt.xticks(
        [i for i in np.arange(0, 1.2, 0.2)],
        ['{:.1f}'.format(i) for i in np.arange(0, 1.2, 0.2)],
        fontsize=18
    )
    plt.yticks([i for i in np.arange(1,len(ynames)+1)], ynames, fontsize=18)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

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

    ax.legend(title='', ncol=2, bbox_to_anchor=(0.5,0.96), frameon=False,
              loc='lower center', fontsize=20)
    plt.xlim(-0.05, 1.05)

    plt.xticks(
        [i for i in np.arange(0, 1.2, 0.2)],
        ['{:.1f}'.format(i) for i in np.arange(0, 1.2, 0.2)],
        fontsize=18
    )
    plt.yticks(fontsize=18)

    plt.title(title, pad=40, fontsize=24)
    plt.xlabel(xlab, fontsize=20)
    plt.ylabel(ylab, fontsize=20)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def lambda_puc_plots():
    """Process control vector data into violin plots. """
    # Useful variables
    path = '2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/control_vectors/'
    appd = '.controls.cg.sorted.mergecg.bed.gz'

    # Load data into DataFrames
    ak1_hi = import_file(path+'FtubeAkapaBC'+appd)
    ak2_hi = import_file(path+'FtubeAkapaBCrep2'+appd)
    an1_hi = import_file(path+'FtubeAneb'+appd)
    an1_lo = import_file(path+'FtubeAneb10ng'+appd)
    an2_lo = import_file(path+'FtubeAneb10ngRep2'+appd)
    an2_hi = import_file(path+'FtubeAnebRep2'+appd)
    ap1_hi = import_file(path+'FtubeApbat'+appd)
    as1_hi = import_file(path+'FtubeAswift'+appd)
    as1_lo = import_file(path+'FtubeAswift10ng'+appd)
    as2_lo = import_file(path+'FtubeAswift10ngRep2'+appd)
    as2_hi = import_file(path+'FtubeAswiftRep2'+appd)

    bk1_hi = import_file(path+'FtubeBkapaBC'+appd)
    bk2_hi = import_file(path+'FtubeBkapaBCrep2'+appd)
    bn1_hi = import_file(path+'FtubeBneb'+appd)
    bn1_lo = import_file(path+'FtubeBneb10ng'+appd)
    bn2_lo = import_file(path+'FtubeBneb10ngRep2'+appd)
    bn2_hi = import_file(path+'FtubeBnebRep2'+appd)
    bp1_hi = import_file(path+'FtubeBpbat'+appd)
    bs1_hi = import_file(path+'FtubeBswift'+appd)
    bs1_lo = import_file(path+'FtubeBswift10ng'+appd)
    bs2_lo = import_file(path+'FtubeBswift10ngRep2'+appd)
    bs2_hi = import_file(path+'FtubeBswiftRep2'+appd)

    # Pull out lambdaphage data from files
    lamb = pd.concat(
        [ak1_hi[(ak1_hi.vector == 'lambda')],
         ak2_hi[(ak2_hi.vector == 'lambda')],
         bk1_hi[(bk1_hi.vector == 'lambda')],
         bk2_hi[(bk2_hi.vector == 'lambda')],
         an1_hi[(an1_hi.vector == 'lambda')],
         an2_hi[(an2_hi.vector == 'lambda')],
         bn1_hi[(bn1_hi.vector == 'lambda')],
         bn2_hi[(bn2_hi.vector == 'lambda')],
         ap1_hi[(ap1_hi.vector == 'lambda')],
         bp1_hi[(bp1_hi.vector == 'lambda')],
         as1_hi[(as1_hi.vector == 'lambda')],
         as2_hi[(as2_hi.vector == 'lambda')],
         bs1_hi[(bs1_hi.vector == 'lambda')],
         bs2_hi[(bs2_hi.vector == 'lambda')],
         an1_lo[(an1_lo.vector == 'lambda')],
         an2_lo[(an2_lo.vector == 'lambda')],
         bn1_lo[(bn1_lo.vector == 'lambda')],
         bn2_lo[(bn2_lo.vector == 'lambda')],
         as1_lo[(as1_lo.vector == 'lambda')],
         as2_lo[(as2_lo.vector == 'lambda')],
         bs1_lo[(bs1_lo.vector == 'lambda')],
         bs2_lo[(bs2_lo.vector == 'lambda')]]
    )

    # Pull out pUC19 data from files
    puck = pd.concat(
        [ak1_hi[(ak1_hi.vector == 'pUC19')],
         ak2_hi[(ak2_hi.vector == 'pUC19')],
         bk1_hi[(bk1_hi.vector == 'pUC19')],
         bk2_hi[(bk2_hi.vector == 'pUC19')],
         an1_hi[(an1_hi.vector == 'pUC19')],
         an2_hi[(an2_hi.vector == 'pUC19')],
         bn1_hi[(bn1_hi.vector == 'pUC19')],
         bn2_hi[(bn2_hi.vector == 'pUC19')],
         ap1_hi[(ap1_hi.vector == 'pUC19')],
         bp1_hi[(bp1_hi.vector == 'pUC19')],
         as1_hi[(as1_hi.vector == 'pUC19')],
         as2_hi[(as2_hi.vector == 'pUC19')],
         bs1_hi[(bs1_hi.vector == 'pUC19')],
         bs2_hi[(bs2_hi.vector == 'pUC19')],
         an1_lo[(an1_lo.vector == 'pUC19')],
         an2_lo[(an2_lo.vector == 'pUC19')],
         bn1_lo[(bn1_lo.vector == 'pUC19')],
         bn2_lo[(bn2_lo.vector == 'pUC19')],
         as1_lo[(as1_lo.vector == 'pUC19')],
         as2_lo[(as2_lo.vector == 'pUC19')],
         bs1_lo[(bs1_lo.vector == 'pUC19')],
         bs2_lo[(bs2_lo.vector == 'pUC19')]]
    )

    create_plot(
        lamb,
        'Lambda Phage Control Retention',
        'Percent Retained',
        '',
        'lamb_control.pdf'
    )

    create_plot(
        puck,
        'pUC19 Control Retention',
        'Percent Retained',
        '',
        'puck_control.pdf'
    )

def mitochondria_plots():
    """Process mitochondrial DNA data into violin plots. """
    # Useful variables
    path = '2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/align/'
    appd = '.cg.sorted.mergecg.bed.gz'

    # Load files into DataFrames
    ak1_hi = import_file(path+'FtubeAkapaBC'+appd, control=False, select_chr='chrM')
    ak2_hi = import_file(path+'FtubeAkapaBCrep2'+appd, control=False, select_chr='chrM')
    an1_hi = import_file(path+'FtubeAneb'+appd, control=False, select_chr='chrM')
    an1_lo = import_file(path+'FtubeAneb10ng'+appd, control=False, select_chr='chrM')
    an2_lo = import_file(path+'FtubeAneb10ngRep2'+appd, control=False, select_chr='chrM')
    an2_hi = import_file(path+'FtubeAnebRep2'+appd, control=False, select_chr='chrM')
    ap1_hi = import_file(path+'FtubeApbat'+appd, control=False, select_chr='chrM')
    as1_hi = import_file(path+'FtubeAswift'+appd, control=False, select_chr='chrM')
    as1_lo = import_file(path+'FtubeAswift10ng'+appd, control=False, select_chr='chrM')
    as2_lo = import_file(path+'FtubeAswift10ngRep2'+appd, control=False, select_chr='chrM')
    as2_hi = import_file(path+'FtubeAswiftRep2'+appd, control=False, select_chr='chrM')

    bk1_hi = import_file(path+'FtubeBkapaBC'+appd, control=False, select_chr='chrM')
    bk2_hi = import_file(path+'FtubeBkapaBCrep2'+appd, control=False, select_chr='chrM')
    bn1_hi = import_file(path+'FtubeBneb'+appd, control=False, select_chr='chrM')
    bn1_lo = import_file(path+'FtubeBneb10ng'+appd, control=False, select_chr='chrM')
    bn2_lo = import_file(path+'FtubeBneb10ngRep2'+appd, control=False, select_chr='chrM')
    bn2_hi = import_file(path+'FtubeBnebRep2'+appd, control=False, select_chr='chrM')
    bp1_hi = import_file(path+'FtubeBpbat'+appd, control=False, select_chr='chrM')
    bs1_hi = import_file(path+'FtubeBswift'+appd, control=False, select_chr='chrM')
    bs1_lo = import_file(path+'FtubeBswift10ng'+appd, control=False, select_chr='chrM')
    bs2_lo = import_file(path+'FtubeBswift10ngRep2'+appd, control=False, select_chr='chrM')
    bs2_hi = import_file(path+'FtubeBswiftRep2'+appd, control=False, select_chr='chrM')

    # Mitochondrial beta values
    chrm = pd.concat(
        [ak1_hi,
         ak2_hi,
         bk1_hi,
         bk2_hi,
         an1_hi,
         an2_hi,
         bn1_hi,
         bn2_hi,
         ap1_hi,
         bp1_hi,
         as1_hi,
         as2_hi,
         bs1_hi,
         bs2_hi,
         an1_lo,
         an2_lo,
         bn1_lo,
         bn2_lo,
         as1_lo,
         as2_lo,
         bs1_lo,
         bs2_lo]
    )

    create_plot(
        chrm,
        'Mitochondrial Retention',
        'Percent Retained',
        '',
        'chrM_control.pdf'
    )

if __name__ == '__main__':
    lambda_puc_plots()
    mitochondria_plots()
