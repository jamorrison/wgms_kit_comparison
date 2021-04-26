"""Main module for making plots regarding statistical measures of kit
   comparison data.
"""
from pathlib import Path
import numpy as np
import json

import tex_tables
import plotting
import utils

def raw_read_base_qual_plots(data, outdir):
    """Make the raw read base quality plots."""
    meta = {
        'r1_read_base_20':  {
            'key'    : 'r1_read_base_20',
            'title'  : 'Read 1 Reads with Avg. Base Quality >= 20',
            'xlab'   : 'Percentage of Reads',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'r1_read_base_20.pdf'
        },
        'r1_read_base_30':  {
            'key'    : 'r1_read_base_30',
            'title'  : 'Read 1 Reads with Avg. Base Quality >= 30',
            'xlab'   : 'Percentage of Reads',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'r1_read_base_30.pdf'
        },
        'r1_low_base_qual': {
            'key'    : 'r1_low_base_qual',
            'title'  : 'Read 1 Bases with Base Quality < 20',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0,  10,  1),
            'figname': outdir+'r1_low_base_qual.pdf'
        },
        'r1_med_base_qual': {
            'key'    : 'r1_med_base_qual',
            'title'  : 'Read 1 Bases with Base Quality >= 20 and <= 30',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0,  10,  1),
            'figname': outdir+'r1_med_base_qual.pdf'
        },
        'r1_hi_base_qual': {
            'key'    : 'r1_hi_base_qual',
            'title'  : 'Read 1 Bases with Base Quality > 30',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0,  100,  10),
            'figname': outdir+'r1_hi_base_qual.pdf'
        },
        'r2_read_base_20':  {
            'key'    : 'r2_read_base_20',
            'title'  : 'Read 2 Reads with Avg. Base Quality >= 20',
            'xlab'   : 'Percentage of Reads',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'r2_read_base_20.pdf'
        },
        'r2_read_base_30':  {
            'key'    : 'r2_read_base_30',
            'title'  : 'Read 2 Reads with Avg. Base Quality >= 30',
            'xlab'   : 'Percentage of Reads',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'r2_read_base_30.pdf'
        },
        'r2_low_base_qual': {
            'key'    : 'r2_low_base_qual',
            'title'  : 'Read 2 Bases with Base Quality < 20',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0,  10,  1),
            'figname': outdir+'r2_low_base_qual.pdf'
        },
        'r2_med_base_qual': {
            'key'    : 'r2_med_base_qual',
            'title'  : 'Read 2 Bases with Base Quality >= 20 and <= 30',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0,  10,  1),
            'figname': outdir+'r2_med_base_qual.pdf'
        },
        'r2_hi_base_qual': {
            'key'    : 'r2_hi_base_qual',
            'title'  : 'Read 2 Bases with Base Quality > 30',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0,  100,  10),
            'figname': outdir+'r2_hi_base_qual.pdf'
        },
    }

    for vals in meta.values():
        plot_data = utils.retrieve_single_element_data_one_output(
            data, vals['key'], True
        )
        plotting.create_rep_avg_plot(
            plot_data,
            vals['title'],
            vals['xlab'],
            vals['ylab'],
            vals['xlims'],
            vals['figname']
        )

    return 0

def cpg_region_plots(data, outdir):
    """Make CpG region plots."""
    meta = {
        'ExonicCpGs_percent_covered_all': {
            'key1'   : 'ExonicCpGs_percent_covered',
            'key2'   : 'All',
            'title'  : 'CpG Coverage of Exonic Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'exon_cpgs_all.pdf'
        },
        'ExonicCpGs_percent_covered_q40': {
            'key1'   : 'ExonicCpGs_percent_covered',
            'key2'    : 'Q40',
            'title'  : 'CpG Coverage of Exonic Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'exon_cpgs_q40.pdf'
        },
        'RepeatCpGs_percent_covered_all': {
            'key1'   : 'RepeatCpGs_percent_covered',
            'key2'   : 'All',
            'title'  : 'CpG Coverage of Repeat-Masked Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'rmsk_cpgs_all.pdf'
        },
        'RepeatCpGs_percent_covered_q40': {
            'key1'   : 'RepeatCpGs_percent_covered',
            'key2'   : 'Q40',
            'title'  : 'CpG Coverage of Repeat-Masked Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'rmsk_cpgs_q40.pdf'
        },
        'GenicCpGs_percent_covered_all': {
            'key1'   : 'GenicCpGs_percent_covered',
            'key2'   : 'All',
            'title'  : 'CpG Coverage of Genic Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'gene_cpgs_all.pdf'
        },
        'GenicCpGs_percent_covered_q40': {
            'key1'   : 'GenicCpGs_percent_covered',
            'key2'   : 'Q40',
            'title'  : 'CpG Coverage of Genic Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'gene_cpgs_q40.pdf'
        },
        'CGICpGs_percent_covered_all': {
            'key1'   : 'CGICpGs_percent_covered',
            'key2'   : 'All',
            'title'  : 'CpG Coverage of CpG Island Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'cgis_cpgs_all.pdf'
        },
        'CGICpGs_percent_covered_q40': {
            'key1'   : 'CGICpGs_percent_covered',
            'key2'   : 'Q40',
            'title'  : 'CpG Coverage of CpG Island Regions',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'cgis_cpgs_q40.pdf'
        },
        'TotalCpGs_percent_covered_all': {
            'key1'   : 'TotalCpGs_percent_covered',
            'key2'   : 'All',
            'title'  : 'Coverage for All CpGs in Genome',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'tote_cpgs_all.pdf'
        },
        'TotalCpGs_percent_covered_q40': {
            'key1'   : 'TotalCpGs_percent_covered',
            'key2'   : 'Q40',
            'title'  : 'Coverage for All CpGs in Genome',
            'xlab'   : 'Percentage of CpGs',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'tote_cpgs_q40.pdf'
        }
    }

    for vals in meta.values():
        plot_data = utils.retrieve_single_element_data_from_dict_one_output(
            data, vals['key1'], vals['key2'], True
        )
        plotting.create_rep_avg_plot(
            plot_data,
            vals['title'],
            vals['xlab'],
            vals['ylab'],
            vals['xlims'],
            vals['figname']
        )

    return 0

def trimmed_stats_plots(data, outdir):
    """Make plots for post-trimming results."""
    meta = {
        'r1_trimmed_bp': {
            'key'    : 'r1_trimmed_bp',
            'title'  : 'Read 1 Percentage of Bases Trimmed',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0, 10, 1),
            'figname': outdir+'r1_trimmed_bp.pdf'
        },
        'r1_quality_bp': {
            'key'    : 'r1_quality_bp',
            'title'  : 'Read 1 Percentage of Bases Trimmed for Low Quality',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0, 1, 0.1),
            'figname': outdir+'r1_quality_bp.pdf'
        },
        'r1_adapter_re': {
            'key'    : 'r1_adapter_re',
            'title'  : 'Read 1 Percentage of Reads with Adapter Contamination',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'r1_adapter_re.pdf'
        },
        'r2_trimmed_bp': {
            'key'    : 'r2_trimmed_bp',
            'title'  : 'Read 2 Percentage of Bases Trimmed',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0, 10, 1),
            'figname': outdir+'r2_trimmed_bp.pdf'
        },
        'r2_quality_bp': {
            'key'    : 'r2_quality_bp',
            'title'  : 'Read 2 Percentage of Bases Trimmed for Low Quality',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0, 1, 0.1),
            'figname': outdir+'r2_quality_bp.pdf'
        },
        'r2_adapter_re': {
            'key'    : 'r2_adapter_re',
            'title'  : 'Read 2 Percentage of Reads with Adapter Contamination',
            'xlab'   : 'Percentage of Bases',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'r2_adapter_re.pdf'
        }
    }

    for vals in meta.values():
        plot_data = utils.retrieve_single_element_data_one_output(
            data, vals['key'], True
        )
        plotting.create_rep_avg_plot(
            [(s, 100*i) for s, i in plot_data],
            vals['title'],
            vals['xlab'],
            vals['ylab'],
            vals['xlims'],
            vals['figname']
        )

    return 0

def read_alignment_plot(data, outdir):
    """Make plot with percentage of optimally/sub-optimally/not aligned reads."""
    A, B = utils.retrieve_aligned_reads(data, True)
    plotting.create_read_alignment_plot(
        A,
        'Fragment Mapping Overview: Sample A',
        'Percentage of Fragments',
        '',
        (0, 100, 10),
        outdir+'read_mapping_sample_A.pdf'
    )
    plotting.create_read_alignment_plot(
        B,
        'Fragment Mapping Overview: Sample B',
        'Percentage of Fragments',
        '',
        (0, 100, 10),
        outdir+'read_mapping_sample_B.pdf'
    )

    A, B = utils.retrieve_data_points_from_dict_in_dict(
        data, 'aligned_reads', 'mapq_percent', False
    )
    plotting.create_line_chart(
        A,
        'Distribution of MAPQ Scores: Sample A',
        'MAPQ Score',
        'Percent of Reads',
        (0,  60, 10),
        (0, 100, 20),
        3,
        'upper left',
        outdir+'mapq_dist_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'Distribution of MAPQ Scores: Sample B',
        'MAPQ Score',
        'Percent of Reads',
        (0,  60, 10),
        (0, 100, 20),
        3,
        'upper left',
        outdir+'mapq_dist_sample_B.pdf'
    )

    plot_data = utils.retrieve_aligned_reads_total_only(data, True)
    plotting.create_rep_avg_plot(
        [(samp, i/1000000.) for samp, i in plot_data],
        'Number of Aligned Read Fragments',
        'Aligned Read Fragments (Millions)',
        '',
        (0, 1800, 300),
        outdir+'aligned_reads.pdf'
    )

def duplicate_rate_plots(data, outdir):
    """Make plots with duplicate read rate percentages."""
    meta = {
        'dup_report_all': {
            'key1'   : 'dup_report',
            'key2'   : 'dup_all',
            'title'  : 'Percentage of Duplicate Marked Reads',
            'xlab'   : 'Percentage of Reads',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'dup_all.pdf'
        },
        'dup_report_q40': {
            'key1'   : 'dup_report',
            'key2'   : 'dup_q40',
            'title'  : 'Percentage of Duplicate Marked Reads',
            'xlab'   : 'Percentage of Reads',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'dup_q40.pdf'
        }
    }

    for vals in meta.values():
        plot_data = utils.retrieve_single_element_data_from_dict_one_output(
            data, vals['key1'], vals['key2'], True
        )
        plotting.create_rep_avg_plot(
            plot_data,
            vals['title'],
            vals['xlab'],
            vals['ylab'],
            vals['xlims'],
            vals['figname']
        )

    return 0

def cpn_retention_plots(data, outdir):
    """Make plots with CpN retention rate percentages."""
    meta = {
        'read_rtn_a': {
            'key1'   : 'read_rtn',
            'key2'   : 'rca',
            'title'  : 'Read-averaged CpA Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'read_rtn_cpa.pdf'
        },
        'read_rtn_c': {
            'key1'   : 'read_rtn',
            'key2'   : 'rcc',
            'title'  : 'Read-averaged CpC Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'read_rtn_cpc.pdf'
        },
        'read_rtn_g': {
            'key1'   : 'read_rtn',
            'key2'   : 'rcg',
            'title'  : 'Read-averaged CpG Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'read_rtn_cpg.pdf'
        },
        'read_rtn_t': {
            'key1'   : 'read_rtn',
            'key2'   : 'rct',
            'title'  : 'Read-averaged CpT Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'read_rtn_cpt.pdf'
        },
        'base_rtn_a': {
            'key1'   : 'base_rtn',
            'key2'   : 'bca',
            'title'  : 'Base-averaged CpA Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'base_rtn_cpa.pdf'
        },
        'base_rtn_c': {
            'key1'   : 'base_rtn',
            'key2'   : 'bcc',
            'title'  : 'Base-averaged CpC Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'base_rtn_cpc.pdf'
        },
        'base_rtn_g': {
            'key1'   : 'base_rtn',
            'key2'   : 'bcg',
            'title'  : 'Base-averaged CpG Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 100, 10),
            'figname': outdir+'base_rtn_cpg.pdf'
        },
        'base_rtn_t': {
            'key1'   : 'base_rtn',
            'key2'   : 'bct',
            'title'  : 'Base-averaged CpT Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'base_rtn_cpt.pdf'
        }
    }

    for vals in meta.values():
        plot_data = utils.retrieve_single_element_data_from_dict_one_output(
            data, vals['key1'], vals['key2'], True
        )
        plotting.create_rep_avg_plot(
            plot_data,
            vals['title'],
            vals['xlab'],
            vals['ylab'],
            vals['xlims'],
            vals['figname']
        )

    return 0

def create_trinucleotide_methylation_plots(data, outdir):
    """Make plots with CpWpN retention rate percentages."""
    meta = {
        'cah_methylation_percent': {
            'key'    : 'cah_methylation_percent',
            'title'  : 'Base-averaged CpApH Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'cah_meth.pdf'
        },
        'cag_methylation_percent': {
            'key'    : 'cag_methylation_percent',
            'title'  : 'Base-averaged CpApG Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'cag_meth.pdf'
        },
        'cth_methylation_percent': {
            'key'    : 'cth_methylation_percent',
            'title'  : 'Base-averaged CpTpH Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'cth_meth.pdf'
        },
        'ctg_methylation_percent': {
            'key'    : 'ctg_methylation_percent',
            'title'  : 'Base-averaged CpTpG Retention',
            'xlab'   : 'Percent Retained',
            'ylab'   : '',
            'xlims'  : (0, 5, 1),
            'figname': outdir+'ctg_meth.pdf'
        },
    }

    for vals in meta.values():
        plot_data = utils.retrieve_single_element_data_one_output(
            data, vals['key'], True
        )
        plotting.create_rep_avg_plot(
            plot_data,
            vals['title'],
            vals['xlab'],
            vals['ylab'],
            vals['xlims'],
            vals['figname']
        )

    return 0

def mbias_plots(data, outdir):
    """Make M-bias plots."""
    A, B = utils.retrieve_data_points_from_dict_in_dict(
        data, 'cpg_rtn_readpos', '1', False
    )
    plotting.create_line_chart(
        A,
        'Read 1 CpG Retention by Read Position: Sample A',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'lower center',
        outdir+'r1_cpg_rtn_readpos_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'Read 1 CpG Retention by Read Position: Sample B',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'lower center',
        outdir+'r1_cpg_rtn_readpos_sample_B.pdf'
    )

    A, B = utils.retrieve_data_points_from_dict_in_dict(
        data, 'cpg_rtn_readpos', '2', False
    )
    plotting.create_line_chart(
        A,
        'Read 2 CpG Retention by Read Position: Sample A',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'lower center',
        outdir+'r2_cpg_rtn_readpos_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'Read 2 CpG Retention by Read Position: Sample B',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'lower center',
        outdir+'r2_cpg_rtn_readpos_sample_B.pdf'
    )

    A, B = utils.retrieve_data_points_from_dict_in_dict(
        data, 'cph_rtn_readpos', '1', False
    )
    plotting.create_line_chart(
        A,
        'Read 1 CpH Retention by Read Position: Sample A',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'upper left',
        outdir+'r1_cph_rtn_readpos_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'Read 1 CpH Retention by Read Position: Sample B',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'upper left',
        outdir+'r1_cph_rtn_readpos_sample_B.pdf'
    )

    A, B = utils.retrieve_data_points_from_dict_in_dict(
        data, 'cph_rtn_readpos', '2', False
    )
    plotting.create_line_chart(
        A,
        'Read 2 CpH Retention by Read Position: Sample A',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'upper left',
        outdir+'r2_cph_rtn_readpos_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'Read 2 CpH Retention by Read Position: Sample B',
        'Read Position',
        'Percent Retained',
        (0, 100, 20),
        (0, 100, 20),
        3,
        'upper left',
        outdir+'r2_cph_rtn_readpos_sample_B.pdf'
    )

    return 0

def create_covdist_plots(data, outdir):
    """Make coverage vs. depth plots."""
    A, B = utils.retrieve_data_points_from_dict(data, 'covdist_all_base', False)
    plotting.create_line_chart(
        A,
        'Base Cumulative Coverage: Sample A',
        'Coverage',
        'Number of Bases (Millions)',
        (0, 50, 5),
        (0, 3500, 500),
        2,
        'upper right',
        outdir+'covdist_all_base_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'Base Cumulative Coverage: Sample B',
        'Coverage',
        'Number of Bases (Millions)',
        (0, 50, 5),
        (0, 3500, 500),
        2,
        'upper right',
        outdir+'covdist_all_base_sample_B.pdf'
    )

    A, B = utils.retrieve_data_points_from_dict(data, 'covdist_q40_base', False)
    plotting.create_line_chart(
        A,
        'Base Cumulative Coverage: Sample A',
        'Coverage',
        'Number of Bases (Millions)',
        (0, 50, 5),
        (0, 3500, 500),
        2,
        'upper right',
        outdir+'covdist_q40_base_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'Base Cumulative Coverage: Sample B',
        'Coverage',
        'Number of Bases (Millions)',
        (0, 50, 5),
        (0, 3500, 500),
        2,
        'upper right',
        outdir+'covdist_q40_base_sample_B.pdf'
    )

    A, B = utils.retrieve_data_points_from_dict(data, 'covdist_all_cpg', False)
    plotting.create_line_chart(
        A,
        'CpG Cumulative Coverage: Sample A',
        'Coverage',
        'Number of CpGs (Millions)',
        (0, 50, 5),
        (0, 35, 5),
        2,
        'upper right',
        outdir+'covdist_all_cpg_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'CpG Cumulative Coverage: Sample B',
        'Coverage',
        'Number of CpGs (Millions)',
        (0, 50, 5),
        (0, 35, 5),
        2,
        'upper right',
        outdir+'covdist_all_cpg_sample_B.pdf'
    )

    A, B = utils.retrieve_data_points_from_dict(data, 'covdist_q40_cpg', False)
    plotting.create_line_chart(
        A,
        'CpG Cumulative Coverage: Sample A',
        'Coverage',
        'Number of CpGs (Millions)',
        (0, 50, 5),
        (0, 35, 5),
        2,
        'upper right',
        outdir+'covdist_q40_cpg_sample_A.pdf'
    )
    plotting.create_line_chart(
        B,
        'CpG Cumulative Coverage: Sample B',
        'Coverage',
        'Number of CpGs (Millions)',
        (0, 50, 5),
        (0, 35, 5),
        2,
        'upper right',
        outdir+'covdist_q40_cpg_sample_B.pdf'
    )

    return 0

def create_insert_size_plots(data, outdir):
    """Create plot of insert sizes."""
    A, B = utils.retrieve_data_points_from_dict_in_dict(
        data, 'isize_data', 'percent', False
    )
    plotting.create_line_chart(
        A,
        'Insert Size: Sample A',
        'Insert Size',
        'Percent of Mapped Reads',
        (0, 700, 100),
        (0, 1.0, 0.1),
        3,
        'upper right',
        outdir+'insert_size_sample_A.pdf',
        every=50
    )
    plotting.create_line_chart(
        B,
        'Insert Size: Sample B',
        'Insert Size',
        'Percent of Mapped Reads',
        (0, 700, 100),
        (0, 1.0, 0.1),
        3,
        'upper right',
        outdir+'insert_size_sample_B.pdf',
        every=50
    )

    return 0

def create_complexity_curve_plots(data, outdir):
    """Create plot of complexity curves."""
    A, B = utils.retrieve_data_points_from_dict(data, 'complexity_curve', False)
    plotting.create_line_chart(
        [(samp, np.array(x)/1000000., np.array(y)/1000000.) for samp, x, y in A],
        'Library Complexity: Sample A',
        'Total Number of Reads (Millions)',
        'Number of Unique Reads (Millions)',
        (0, 600, 100),
        (0, 500, 100),
        2,
        'upper left',
        outdir+'complexity_sample_A.pdf',
        every=50
    )
    plotting.create_line_chart(
        [(samp, np.array(x)/1000000., np.array(y)/1000000.) for samp, x, y in B],
        'Library Complexity: Sample B',
        'Total Number of Reads (Millions)',
        'Number of Unique Reads (Millions)',
        (0, 900, 100),
        (0, 600, 100),
        2,
        'lower right',
        outdir+'complexity_sample_B.pdf',
        every=50
    )

    return 0

def create_tex_table(data, outdir):
    """Create table with CpG regional coverages."""
    avg_depth = {
        'exon_q40': utils.retrieve_single_element_data(data, 'exon_q40_avg_depth', False),
        'gene_q40': utils.retrieve_single_element_data(data, 'gene_q40_avg_depth', False),
        'rmsk_q40': utils.retrieve_single_element_data(data, 'rmsk_q40_avg_depth', False),
        'cgis_q40': utils.retrieve_single_element_data(data, 'cgis_q40_avg_depth', False),
        'tote_q40': utils.retrieve_single_element_data(data, 'total_q40_avg_depth', False)
    }
    uniformity = {
        'q40_base_mu': utils.retrieve_single_element_data_from_dict(data, 'uniformity', 'q40_base_mu', False),
        'q40_base_cv': utils.retrieve_single_element_data_from_dict(data, 'uniformity', 'q40_base_cv', False),
        'q40_cpg_mu': utils.retrieve_single_element_data_from_dict(data, 'uniformity', 'q40_cpg_mu', False),
        'q40_cpg_cv': utils.retrieve_single_element_data_from_dict(data, 'uniformity', 'q40_cpg_cv', False)
    }
    tex_tables.table_creator(outdir+'kit_comp_tables.tex', avg_depth, uniformity)

    return 0

def create_obs_exp_ratio_plots(data, outdir):
    """Make plots to show the observed/expected ratio."""
    meta = {
        'obsexp_mappy_cpgs': {
            'key'   : 'mappy_cpgs_obs_exp_ratio',
            'title'  : 'Coverage Ratio for All CpGs',
            'xlab'   : 'Observed Coverage / Expected Coverage',
            'ylab'   : '',
            'xlims'  : (0, 2, 0.25),
            'figname': outdir+'obsexp_mappy_cpgs.pdf'
        },
        'obsexp_mappy_cgis': {
            'key'   : 'mappy_cpis_obs_exp_ratio',
            'title'  : 'Coverage Ratio for CpG Islands',
            'xlab'   : 'Observed Coverage / Expected Coverage',
            'ylab'   : '',
            'xlims'  : (0, 2, 0.25),
            'figname': outdir+'obsexp_mappy_cgis.pdf'
        },
        'obsexp_mappy_rmsk': {
            'key'   : 'mappy_rmsk_obs_exp_ratio',
            'title'  : 'Coverage Ratio for Repeat-Masked Regions',
            'xlab'   : 'Observed Coverage / Expected Coverage',
            'ylab'   : '',
            'xlims'  : (0, 2, 0.25),
            'figname': outdir+'obsexp_mappy_rmsk.pdf'
        },
        'obsexp_mappy_exon': {
            'key'   : 'mappy_exon_obs_exp_ratio',
            'title'  : 'Coverage Ratio for Exonic Regions',
            'xlab'   : 'Observed Coverage / Expected Coverage',
            'ylab'   : '',
            'xlims'  : (0, 2, 0.25),
            'figname': outdir+'obsexp_mappy_exon.pdf'
        },
        'obsexp_mappy_gene': {
            'key'   : 'mappy_gene_obs_exp_ratio',
            'title'  : 'Coverage Ratio for Genic Regions',
            'xlab'   : 'Observed Coverage / Expected Coverage',
            'ylab'   : '',
            'xlims'  : (0, 2, 0.25),
            'figname': outdir+'obsexp_mappy_gene.pdf'
        },
        'obsexp_mappy_intr': {
            'key'   : 'mappy_intr_obs_exp_ratio',
            'title'  : 'Coverage Ratio for Intergenic Regions',
            'xlab'   : 'Observed Coverage / Expected Coverage',
            'ylab'   : '',
            'xlims'  : (0, 2, 0.25),
            'figname': outdir+'obsexp_mappy_intr.pdf'
        }
    }

    for vals in meta.values():
        plot_data = utils.retrieve_single_element_data_one_output(
            data, vals['key'], True
        )
        plotting.create_rep_avg_plot(
            plot_data,
            vals['title'],
            vals['xlab'],
            vals['ylab'],
            vals['xlims'],
            vals['figname'],
            add_line=True
        )

    return 0


def main():
    """Load data and generate plots."""
    with open('2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/analyze_the_data/collect_data/kit_comp_collected_data.json', 'r') as f:
        kit_comp_data = json.load(f)

    with open('2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/analyze_the_data/collect_data/kit_comp_collected_data_subsampled.json', 'r') as f:
        kit_comp_data_sub = json.load(f)

    plot_dir = 'color_plots'
    do_raw_bams = True
    do_sub_bams = True

    # Raw BAM plots
    if do_raw_bams:
        raw_bam_dir = plot_dir + '/raw_bam_plots/'
        Path(raw_bam_dir).mkdir(parents=True, exist_ok=True)

        raw_read_base_qual_plots(kit_comp_data, raw_bam_dir)
        cpg_region_plots(kit_comp_data, raw_bam_dir)
        trimmed_stats_plots(kit_comp_data, raw_bam_dir)
        read_alignment_plot(kit_comp_data, raw_bam_dir)
        duplicate_rate_plots(kit_comp_data, raw_bam_dir)
        cpn_retention_plots(kit_comp_data, raw_bam_dir)
        mbias_plots(kit_comp_data, raw_bam_dir)
        create_covdist_plots(kit_comp_data, raw_bam_dir)
        create_insert_size_plots(kit_comp_data, raw_bam_dir)
        create_complexity_curve_plots(kit_comp_data, raw_bam_dir)
        create_tex_table(kit_comp_data, raw_bam_dir)
        create_obs_exp_ratio_plots(kit_comp_data, raw_bam_dir)
        create_trinucleotide_methylation_plots(kit_comp_data, raw_bam_dir)

    # Subsampled BAM plots
    if do_sub_bams:
        sub_bam_dir = plot_dir + '/sub_bam_plots/'
        Path(sub_bam_dir).mkdir(parents=True, exist_ok=True)

        cpg_region_plots(kit_comp_data_sub, sub_bam_dir)
        cpn_retention_plots(kit_comp_data_sub, sub_bam_dir)
        mbias_plots(kit_comp_data_sub, sub_bam_dir)
        create_covdist_plots(kit_comp_data_sub, sub_bam_dir)
        create_tex_table(kit_comp_data_sub, sub_bam_dir)
        create_trinucleotide_methylation_plots(kit_comp_data_sub, sub_bam_dir)

if __name__ == '__main__':
    main()

