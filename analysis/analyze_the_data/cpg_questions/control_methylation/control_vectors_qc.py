"""Create plots to show methylation for mitochondrial and control vector DNA."""
import matplotlib.pyplot as plt
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
    df['sample'] = samp

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

    # y-axis tick names
    names = [
        'Low Swift Rep. 2', 'Low Swift Rep. 1',
        'Low NEB Rep. 2', 'Low NEB Rep. 1',
        'Swift Rep. 2', 'Swift Rep. 1',
        'PBAT Rep. 1',
        'NEB Rep. 2', 'NEB Rep. 1',
        'Kapa Rep. 2', 'Kapa Rep. 1'
    ]

    # Pull out lambdaphage data from Sample A files
    a_lamb = [
        as2_lo[(as2_lo.vector == 'lambda') & (as2_lo.covg > 0)]['beta'].values,
        as1_lo[(as1_lo.vector == 'lambda') & (as1_lo.covg > 0)]['beta'].values,
        an2_lo[(an2_lo.vector == 'lambda') & (an2_lo.covg > 0)]['beta'].values,
        an1_lo[(an1_lo.vector == 'lambda') & (an1_lo.covg > 0)]['beta'].values,
        as2_hi[(as2_hi.vector == 'lambda') & (as2_hi.covg > 0)]['beta'].values,
        as1_hi[(as1_hi.vector == 'lambda') & (as1_hi.covg > 0)]['beta'].values,
        ap1_hi[(ap1_hi.vector == 'lambda') & (ap1_hi.covg > 0)]['beta'].values,
        an2_hi[(an2_hi.vector == 'lambda') & (an2_hi.covg > 0)]['beta'].values,
        an1_hi[(an1_hi.vector == 'lambda') & (an1_hi.covg > 0)]['beta'].values,
        ak2_hi[(ak2_hi.vector == 'lambda') & (ak2_hi.covg > 0)]['beta'].values,
        ak1_hi[(ak1_hi.vector == 'lambda') & (ak1_hi.covg > 0)]['beta'].values
    ]

    # Pull out pUC19 data from Sample B files
    a_puck = [
        as2_lo[(as2_lo.vector == 'pUC19') & (as2_lo.covg > 0)]['beta'].values,
        as1_lo[(as1_lo.vector == 'pUC19') & (as1_lo.covg > 0)]['beta'].values,
        an2_lo[(an2_lo.vector == 'pUC19') & (an2_lo.covg > 0)]['beta'].values,
        an1_lo[(an1_lo.vector == 'pUC19') & (an1_lo.covg > 0)]['beta'].values,
        as2_hi[(as2_hi.vector == 'pUC19') & (as2_hi.covg > 0)]['beta'].values,
        as1_hi[(as1_hi.vector == 'pUC19') & (as1_hi.covg > 0)]['beta'].values,
        ap1_hi[(ap1_hi.vector == 'pUC19') & (ap1_hi.covg > 0)]['beta'].values,
        an2_hi[(an2_hi.vector == 'pUC19') & (an2_hi.covg > 0)]['beta'].values,
        an1_hi[(an1_hi.vector == 'pUC19') & (an1_hi.covg > 0)]['beta'].values,
        ak2_hi[(ak2_hi.vector == 'pUC19') & (ak2_hi.covg > 0)]['beta'].values,
        ak1_hi[(ak1_hi.vector == 'pUC19') & (ak1_hi.covg > 0)]['beta'].values
    ]

    # Pull out lambdaphage data from Sample B files
    b_lamb = [
        bs2_lo[(bs2_lo.vector == 'lambda') & (bs2_lo.covg > 0)]['beta'].values,
        bs1_lo[(bs1_lo.vector == 'lambda') & (bs1_lo.covg > 0)]['beta'].values,
        bn2_lo[(bn2_lo.vector == 'lambda') & (bn2_lo.covg > 0)]['beta'].values,
        bn1_lo[(bn1_lo.vector == 'lambda') & (bn1_lo.covg > 0)]['beta'].values,
        bs2_hi[(bs2_hi.vector == 'lambda') & (bs2_hi.covg > 0)]['beta'].values,
        bs1_hi[(bs1_hi.vector == 'lambda') & (bs1_hi.covg > 0)]['beta'].values,
        bp1_hi[(bp1_hi.vector == 'lambda') & (bp1_hi.covg > 0)]['beta'].values,
        bn2_hi[(bn2_hi.vector == 'lambda') & (bn2_hi.covg > 0)]['beta'].values,
        bn1_hi[(bn1_hi.vector == 'lambda') & (bn1_hi.covg > 0)]['beta'].values,
        bk2_hi[(bk2_hi.vector == 'lambda') & (bk2_hi.covg > 0)]['beta'].values,
        bk1_hi[(bk1_hi.vector == 'lambda') & (bk1_hi.covg > 0)]['beta'].values
    ]

    # Pull out pUC19 data from Sample B files
    b_puck = [
        bs2_lo[(bs2_lo.vector == 'pUC19') & (bs2_lo.covg > 0)]['beta'].values,
        bs1_lo[(bs1_lo.vector == 'pUC19') & (bs1_lo.covg > 0)]['beta'].values,
        bn2_lo[(bn2_lo.vector == 'pUC19') & (bn2_lo.covg > 0)]['beta'].values,
        bn1_lo[(bn1_lo.vector == 'pUC19') & (bn1_lo.covg > 0)]['beta'].values,
        bs2_hi[(bs2_hi.vector == 'pUC19') & (bs2_hi.covg > 0)]['beta'].values,
        bs1_hi[(bs1_hi.vector == 'pUC19') & (bs1_hi.covg > 0)]['beta'].values,
        bp1_hi[(bp1_hi.vector == 'pUC19') & (bp1_hi.covg > 0)]['beta'].values,
        bn2_hi[(bn2_hi.vector == 'pUC19') & (bn2_hi.covg > 0)]['beta'].values,
        bn1_hi[(bn1_hi.vector == 'pUC19') & (bn1_hi.covg > 0)]['beta'].values,
        bk2_hi[(bk2_hi.vector == 'pUC19') & (bk2_hi.covg > 0)]['beta'].values,
        bk1_hi[(bk1_hi.vector == 'pUC19') & (bk1_hi.covg > 0)]['beta'].values
    ]

    create_plot(
        a_lamb,
        'Lambdaphage Control Methylation: Sample A',
        'Methylation Level',
        '',
        names,
        'a_lamb_control.pdf'
    )

    create_plot(
        a_puck,
        'pUC19 Control Methylation: Sample A',
        'Methylation Level',
        '',
        names,
        'a_puck_control.pdf'
    )

    create_plot(
        b_lamb,
        'Lambdaphage Control Methylation: Sample B',
        'Methylation Level',
        '',
        names,
        'b_lamb_control.pdf'
    )

    create_plot(
        b_puck,
        'pUC19 Control Methylation: Sample B',
        'Methylation Level',
        '',
        names,
        'b_puck_control.pdf'
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

    # y-axis tick  names
    names = [
        'Low Swift Rep. 2', 'Low Swift Rep. 1',
        'Low NEB Rep. 2', 'Low NEB Rep. 1',
        'Swift Rep. 2', 'Swift Rep. 1',
        'PBAT Rep. 1',
        'NEB Rep. 2', 'NEB Rep. 1',
        'Kapa Rep. 2', 'Kapa Rep. 1'
    ]

    # Mitochondrial beta values for Sample A
    a_chrm = [
        as2_lo['beta'].values,
        as1_lo['beta'].values,
        an2_lo['beta'].values,
        an1_lo['beta'].values,
        as2_hi['beta'].values,
        as1_hi['beta'].values,
        ap1_hi['beta'].values,
        an2_hi['beta'].values,
        an1_hi['beta'].values,
        ak2_hi['beta'].values,
        ak1_hi['beta'].values
    ]

    # Mitochondrial beta values for Sample B
    b_chrm = [
        bs2_lo['beta'].values,
        bs1_lo['beta'].values,
        bn2_lo['beta'].values,
        bn1_lo['beta'].values,
        bs2_hi['beta'].values,
        bs1_hi['beta'].values,
        bp1_hi['beta'].values,
        bn2_hi['beta'].values,
        bn1_hi['beta'].values,
        bk2_hi['beta'].values,
        bk1_hi['beta'].values
    ]

    create_plot(
        a_chrm,
        'Mitochondrial Methylation: Sample A',
        'Methylation Level',
        '',
        names,
        'a_chrM_control.pdf'
    )

    create_plot(
        b_chrm,
        'Mitochondrial Methylation: Sample B',
        'Methylation Level',
        '',
        names,
        'b_chrM_control.pdf'
    )

if __name__ == '__main__':
    lambda_puc_plots()
    mitochondria_plots()
