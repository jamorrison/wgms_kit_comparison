"""Create tables for the kit comparison analysis."""
import constants

def create_article_header(fname):
    """Add stuff at the top of the document."""
    with open(fname, 'w') as f:
        f.write('\\documentclass[12pt, letterpaper]{article}\n')

        f.write('\\usepackage[utf8]{inputenc}\n')
        f.write('\\usepackage{booktabs}\n')
        f.write('\\usepackage{multirow}\n')
        f.write('\\usepackage{graphicx}\n\n')

        f.write('\\begin{document}\n\n')

def create_article_footer(fname):
    """Add stuff at the end of the document."""
    with open(fname, 'a') as f:
        f.write('\\end{document}\n')

def format_covg_data(data):
    """Put the coverage data in a format that's easy for table making."""
    out = {}
    for key, val in data.items():
        for v in val:
            for i in v:
                if i[0] not in out.keys():
                    out[i[0]] = {}
                out[i[0]][key] = i[1]

    return out

def create_covg_table(fname, table_data):
    """Add average CpG depth table to document."""
    samps = sorted(
        list(constants.SAMPLE_NAMES.keys()),
        key=lambda d: constants.INTERLEAVE_ORDER[d],
        reverse=False
    )

    data = format_covg_data(table_data)

    with open(fname, 'a') as f:
        f.write('\\begin{center}\n')
        f.write('    \\resizebox{\\textwidth}{!}{\\begin{tabular}{ lcc | c | c | c | c | c }\n')
        f.write('    \\toprule\n')
        f.write('    \multicolumn{3}{c}{Kit Info} & Exonic & Genic & Repeat-masked & CpG Island & Total \\\\\n')
        f.write('    Kit & Rep. & Sample & Q40 & Q40 & Q40 & Q40 & Q40 \\\\\n')
        f.write('    \\midrule\n')
        for idx, samp in enumerate(samps):
            dat = data[samp]
            sam = constants.SAMPLE_GROUP[samp]
            if (sam == 'A'):
                if constants.SAMPLE_REP[samp] == '1' and constants.SAMPLE_KIT[samp] != 'PBAT':
                    f.write(
                        '    \multirow{{4}}{{*}}{{{}}} & \multirow{{2}}{{*}}{{{}}} & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                        constants.SAMPLE_KIT[samp],
                        constants.SAMPLE_REP[samp],
                        sam,
                        dat['exon_q40'],
                        dat['gene_q40'],
                        dat['rmsk_q40'],
                        dat['cgis_q40'],
                        dat['tote_q40']
                        )
                    )
                else:
                    if constants.SAMPLE_KIT[samp] == 'PBAT':
                        f.write(
                            '    \multirow{{2}}{{*}}{{{}}} & \multirow{{2}}{{*}}{{{}}} & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                            constants.SAMPLE_KIT[samp],
                            constants.SAMPLE_REP[samp],
                            sam,
                            dat['exon_q40'],
                            dat['gene_q40'],
                            dat['rmsk_q40'],
                            dat['cgis_q40'],
                            dat['tote_q40']
                            )
                        )
                    else:
                        f.write(
                            '    & \multirow{{2}}{{*}}{{{}}} & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                            constants.SAMPLE_REP[samp],
                            sam,
                            dat['exon_q40'],
                            dat['gene_q40'],
                            dat['rmsk_q40'],
                            dat['cgis_q40'],
                            dat['tote_q40']
                            )
                        )
            else:
                f.write(
                    '    & & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                    sam,
                    dat['exon_q40'],
                    dat['gene_q40'],
                    dat['rmsk_q40'],
                    dat['cgis_q40'],
                    dat['tote_q40']
                    )
                )
                if (idx+1 != len(samps)) and ((constants.SAMPLE_KIT[samp] == 'PBAT') or (constants.SAMPLE_REP[samp] == '2')):
                    f.write('    \\midrule\n')
        f.write('    \\bottomrule\n')
        f.write('    \\end{tabular}}\n')
        f.write('\\end{center}\n\n')

    return 0

def create_unif_table(fname, table_data):
    """Add uniformity table to document."""
    samps = sorted(
        list(constants.SAMPLE_NAMES.keys()),
        key=lambda d: constants.INTERLEAVE_ORDER[d],
        reverse=False
    )

    data = format_covg_data(table_data)

    with open(fname, 'a') as f:
        f.write('\\begin{center}\n')
        f.write('    \\resizebox{\\textwidth}{!}{\\begin{tabular}{ lcc | cc | cc }\n')
        f.write('    \\toprule\n')
        f.write('    \multicolumn{3}{c}{Kit Info} & \multicolumn{2}{c}{Base Covg. (Q40)} & \multicolumn{2}{c}{CpG Covg. (Q40)} \\\\\n')
        f.write('    Kit & Rep. & Sample & Mean & CoV & Mean & CoV \\\\\n')
        f.write('    \\midrule\n')
        for idx, samp in enumerate(samps):
            dat = data[samp]
            sam = constants.SAMPLE_GROUP[samp]
            if (sam == 'A'):
                if constants.SAMPLE_REP[samp] == '1' and constants.SAMPLE_KIT[samp] != 'PBAT':
                    f.write(
                        '    \multirow{{4}}{{*}}{{{}}} & \multirow{{2}}{{*}}{{{}}} & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                        constants.SAMPLE_KIT[samp],
                        constants.SAMPLE_REP[samp],
                        sam,
                        dat['q40_base_mu'], dat['q40_base_cv'],
                        dat['q40_cpg_mu'], dat['q40_cpg_cv']
                        )
                    )
                else:
                    if constants.SAMPLE_KIT[samp] == 'PBAT':
                        f.write(
                            '    \multirow{{2}}{{*}}{{{}}} & \multirow{{2}}{{*}}{{{}}} & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                            constants.SAMPLE_KIT[samp],
                            constants.SAMPLE_REP[samp],
                            sam,
                            dat['q40_base_mu'], dat['q40_base_cv'],
                            dat['q40_cpg_mu'], dat['q40_cpg_cv']
                            )
                        )
                    else:
                        f.write(
                            '    & \multirow{{2}}{{*}}{{{}}} & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                            constants.SAMPLE_REP[samp],
                            sam,
                            dat['q40_base_mu'], dat['q40_base_cv'],
                            dat['q40_cpg_mu'], dat['q40_cpg_cv']
                            )
                        )
            else:
                f.write(
                    '    & & {} & {:.1f} & {:.1f} & {:.1f} & {:.1f} \\\\\n'.format(
                    sam,
                    dat['q40_base_mu'], dat['q40_base_cv'],
                    dat['q40_cpg_mu'], dat['q40_cpg_cv']
                    )
                )
                if (idx+1 != len(samps)) and ((constants.SAMPLE_KIT[samp] == 'PBAT') or (constants.SAMPLE_REP[samp] == '2')):
                    f.write('    \\midrule\n')
        f.write('    \\bottomrule\n')
        f.write('    \\end{tabular}}\n')
        f.write('\\end{center}\n\n')

    return 0

def table_creator(fname, covg_data, unif_data):
    """Make the tables for the kit comparison."""
    create_article_header(fname)

    create_covg_table(fname, covg_data)
    create_unif_table(fname, unif_data)

    create_article_footer(fname)
