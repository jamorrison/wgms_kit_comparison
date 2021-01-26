"""Function to generate plots using provided data."""
from matplotlib import pyplot as plt
import numpy as np
import sys

import constants

def create_single_column_bar_chart(plot_data, title, xlab, ylab, xlims, figname):
    """Create bar chart where each entry has one column.
       
       Inputs: plot_data - list of tuples in form (entry name, entry value)
               title     - title of plot
               xlab      - x-axis label
               ylab      - y-axis label
               xlims     - tuple of (xmin, xmax, step) for plt.xlim
               figname   - name of file to save figure as
    """
    nams = [constants.SAMPLE_NAMES[tup[0]] for tup in plot_data]
    vals = [tup[1] for tup in plot_data]
    ypos = [i for i, _ in enumerate(nams)]

    fig = plt.figure(figsize=(10,5))
    plt.tight_layout()

    bars = plt.barh(
        ypos, vals, height=0.7,
        color=[constants.SAMPLE_COLOR[tup[0]] for tup in plot_data]
    )

    plt.xlim(xlims[0], xlims[1]+(0.5*xlims[2]))

    plt.title(title, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    if not (isinstance(xlims[2], int)):
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            ['{:.1f}'.format(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    else:
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            [str(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    plt.yticks(ypos, nams, fontsize=12)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def create_single_column_bar_chart_with_labels(plot_data, text_labels, tag, title, xlab, ylab, xlims, figname):
    """Create bar chart where each entry has one column.
       
       Inputs: plot_data   - list of tuples in form (entry name, entry value)
               text_labels - list of tuples in form (entry name, entry label)
               tag         - string to tag in front of the text label
               title       - title of plot
               xlab        - x-axis label
               ylab        - y-axis label
               xlims       - tuple of (xmin, xmax, step) for plt.xlim
               figname     - name of file to save figure as
    """
    nams = [constants.SAMPLE_NAMES[tup[0]] for tup in plot_data]
    vals = [tup[1] for tup in plot_data]
    ypos = [i for i, _ in enumerate(nams)]

    labs = ['{}{:.2f}'.format(tag, tup[1]) for tup in text_labels]
    cens = [val / 2 for val in vals]

    fig, ax = plt.subplots(figsize=(10,5))
    plt.tight_layout()

    bars = plt.barh(
        ypos, vals, height=0.7,
        color=[constants.SAMPLE_COLOR[tup[0]] for tup in plot_data]
    )

    for (y, x, lab) in zip(ypos, cens, labs):
        ax.text(x, y, lab, ha='center', va='center', color='white')

    plt.xlim(xlims[0], xlims[1]+(0.5*xlims[2]))

    plt.title(title, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    if not (isinstance(xlims[2], int)):
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            ['{:.1f}'.format(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    else:
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            [str(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    plt.yticks(ypos, nams, fontsize=12)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def create_double_column_bar_chart(plot_data_1, sample_1, plot_data_2, sample_2, title, xlab, ylab, xlims, figname):
    """Create bar chart where each entry has two columns.
       
       Inputs: plot_data_1 - list of tuples in form (entry name, entry value)
               plot_data_2 - list of tuples in form (entry name, entry value)
               title       - title of plot
               xlab        - x-axis label
               ylab        - y-axis label
               xlims       - tuple of (xmin, xmax, step) for plt.xlim
               figname     - name of file to save figure as
    """
    if len(plot_data_1) != len(plot_data_2):
        print('Supplied data lengths are not equal for {}'.format(title))
        sys.exit(1)
    for idx, tup in enumerate(plot_data_1):
        chk1 = tup[0].replace('FtubeA', '').replace('FtubeB', '')
        chk2 = plot_data_2[idx][0].replace('FtubeA', '').replace('FtubeB', '')
        if chk1 != chk2:
            print('{}: data has different sample order {} != {}'.format(title, chk1, chk2))
            sys.exit(1)

    bar_height = 0.35

    nams = [constants.SAMPLE_NAMES[tup[0]] for tup in plot_data_1]
    val1 = [tup[1] for tup in plot_data_1]
    val2 = [tup[1] for tup in plot_data_2]
    ypo1 = [i+bar_height/2+0.03 for i, _ in enumerate(nams)]
    ypo2 = [i-bar_height/2-0.03 for i, _ in enumerate(nams)]

    fig, ax = plt.subplots(figsize=(10,5))
    plt.tight_layout()

    plt.barh(-1, 0, height=0, label=sample_1, color='slategrey')
    plt.barh(-2, 0, height=0, label=sample_2, color='slategrey',
        hatch='//', edgecolor='whitesmoke', linewidth=0)

    plt.barh(
        ypo1, val1, height=bar_height,
        color=[constants.SAMPLE_COLOR[tup[0]] for tup in plot_data_1]
    )
    plt.barh(
        ypo2, val2, height=bar_height,
        hatch='//', edgecolor='whitesmoke', linewidth=0,
        color=[constants.SAMPLE_COLOR[tup[0]] for tup in plot_data_2]
    )

    plt.xlim(xlims[0], xlims[1]+(0.5*xlims[2]))

    ax.legend(ncol=2, bbox_to_anchor=(0.0, 1.03), frameon=False,
              loc='center', fontsize='large')

    plt.title(title, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    if not (isinstance(xlims[2], int)):
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            ['{:.1f}'.format(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    else:
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            [str(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    plt.yticks([i for i, _ in enumerate(nams)], nams, fontsize=12)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def create_rep_avg_plot(plot_data, title, xlab, ylab, xlims, figname, add_line=False):
    """Create plot with replicates on same y-value with line showing average.
       
       Inputs: plot_data - list of tuples in form (entry name, entry value)
               title     - title of plot
               xlab      - x-axis label
               ylab      - y-axis label
               xlims     - tuple of (xmin, xmax, step) for plt.xlim
               figname   - name of file to save figure as
               add_line  - include line at 1.0 to show ratio of 1
    """
    width = 0.6
    fig, ax = plt.subplots(figsize=(10, 5))
    plt.tight_layout()

    # Make some place holder entries for legend
    for key, value in constants.REP_FORMAT.items():
        plt.plot(-1, 0, value, color='black', fillstyle='none', markersize=8,
                 label='Rep. {}'.format(key))
    ax.vlines(-1, 1, 2, linestyles='solid', label='Rep. Average')

    nams = [] # y-axis tick names
    for idx, tup in enumerate(plot_data):
        x = tup[1]
        y = constants.SAMPLE_PLOT_VALUE[tup[0]]
        f = constants.REP_FORMAT[constants.SAMPLE_REP[tup[0]]]
        l = 'Rep. {}'.format(constants.SAMPLE_REP[tup[0]])
        c = constants.SAMPLE_COLOR[tup[0]]
        #c = 'black'
        plt.plot(x, y, f, fillstyle='none', markersize=8, color=c)

        if ('rep2' not in tup[0]) and ('Rep2' not in tup[0]):
            nams.append((y, tup[0]))

            if ('pbat' not in tup[0]):
                avg = (x + plot_data[idx-1][1])/2
                ax.vlines(avg, y-width/2, y+width/2, linestyles='solid')

    if add_line:
        ax.axvline(1.0, alpha=0.6, color='grey', linestyle='--')

    nams = sorted(nams, key=lambda d: d[0], reverse=False)

    #if len(title) > 40:
    #    ax.legend(ncol=2, bbox_to_anchor=(0.15, 1.03), frameon=False,
    #            loc='right', fontsize='x-large')
    #else:
    #    ax.legend(ncol=3, bbox_to_anchor=(0.25, 1.04), frameon=False,
    #            loc='right', fontsize='x-large')
    ax.legend(ncol=3, bbox_to_anchor=(0.5, 1), frameon=False,
              loc='lower center', fontsize='large')

    plt.title(title, pad=40, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    plt.xlim(xlims[0], xlims[1]+(0.5*xlims[2]))

    if not (isinstance(xlims[2], int)):
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            ['{:.1f}'.format(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    else:
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            [str(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    plt.yticks(
        [tup[0] for tup in nams],
        [' '.join([constants.SAMPLE_KIT[tup[1]], 'Sample', constants.SAMPLE_GROUP[tup[1]]]) for tup in nams],
        fontsize=12
    )

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def create_read_alignment_plot(plot_data, title, xlab, ylab, xlims, figname):
    """Create stacked bar chart to show read alignment stats.
       
       Inputs: plot_data   - list of tuples in form (entry name, entry value)
               title       - title of plot
               xlab        - x-axis label
               ylab        - y-axis label
               xlims       - tuple of (xmin, xmax, step) for plt.xlim
               figname     - name of file to save figure as
    """
    cats = ['Optimally Aligned', 'Sub-optimally Aligned', 'Not Aligned']
    cols = ['#a50f15', '#de2d26', '#fcae91']

    nams = [constants.SAMPLE_NAMES[tup[0]] for tup in plot_data]
    vals = np.array([tup[1] for tup in plot_data])
    ypos = [i for i, _ in enumerate(nams)]

    labs = ['Total Fragments: {:,}'.format(tup[2]) for tup in plot_data]

    vals_cum = vals.cumsum(axis=1)

    fig, ax = plt.subplots(figsize=(10, 5))
    plt.tight_layout()

    for i, (colname, color) in enumerate(zip(cats, cols)):
        widths = vals[:, i]
        starts = vals_cum[:, i] - widths
        ax.barh(ypos, widths, left=starts, height=0.7,
                label=colname, color=color)

    for y, lab in enumerate(labs):
        ax.text(50, y, lab, ha='center', va='center',
                color='white', fontsize='large')

    ax.legend(ncol=len(cats), bbox_to_anchor=(0.5, 1), frameon=False,
              loc='lower center', fontsize='large')

    plt.xlim(xlims[0], xlims[1]+(0.5*xlims[2]))

    plt.title(title, pad=40, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    if not (isinstance(xlims[2], int)):
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            ['{:.1f}'.format(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    else:
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            [str(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    plt.yticks(ypos, nams, fontsize=12)

    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')

def create_line_chart(plot_data, title, xlab, ylab, xlims, ylims, ncols, legend_loc, figname, every=5):
    """Create line chart.
       
       Inputs: plot_data   - list of tuples in form (entry name, entry value)
               title       - title of plot
               xlab        - x-axis label
               ylab        - y-axis label
               xlims       - tuple of (xmin, xmax, step) for plt.xlim
               ylims       - tuple of (ymin, ymax, step) for plt.ylim
               ncols       - number of legend columns
               legend_loc  - location of legend
               figname     - name of file to save figure as
               every       - place marker at 'every' intervals (default: 5)
    """
    fig, ax = plt.subplots(figsize=(10,5))
    plt.tight_layout()

    for dat in plot_data:
        lab = constants.SAMPLE_NAMES[dat[0]]
        stl = constants.SAMPLE_STYLE[dat[0]]
        col = constants.SAMPLE_COLOR[dat[0]]
        x = dat[1]
        y = dat[2]

        plt.plot(x, y, stl, color=col, label=lab, markersize=3.5, markevery=every)

    plt.title(title, fontsize=18)
    plt.xlabel(xlab, fontsize=16)
    plt.ylabel(ylab, fontsize=16)

    plt.xlim(xlims[0]-(0.25*xlims[2]), xlims[1]+(0.25*xlims[2]))
    plt.ylim(ylims[0]-(0.25*ylims[2]), ylims[1]+(0.25*ylims[2]))

    ax.legend(ncol=ncols, loc=legend_loc, fontsize='large')

    if not (isinstance(xlims[2], int)):
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            ['{:.1f}'.format(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    else:
        plt.xticks(
            [i for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            [str(i) for i in np.arange(xlims[0], xlims[1]+xlims[2], xlims[2])],
            fontsize=12
        )
    if not (isinstance(ylims[2], int)):
        plt.yticks(
            [i for i in np.arange(ylims[0], ylims[1]+ylims[2], ylims[2])],
            ['{:.1f}'.format(i) for i in np.arange(ylims[0], ylims[1]+ylims[2], ylims[2])],
            fontsize=12
        )
    else:
        plt.yticks(
            [i for i in np.arange(ylims[0], ylims[1]+ylims[2], ylims[2])],
            [str(i) for i in np.arange(ylims[0], ylims[1]+ylims[2], ylims[2])],
            fontsize=12
        )


    plt.savefig(figname, bbox_inches='tight')
    plt.close('all')
