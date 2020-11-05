"""Main script for collecting data."""
import glob
import json
import os

import read_biscuitqc
import samtools_stats
import preseq_reports
import obs_exp_ratio
import cpg_coverage
import trim_reports
import raw_read

TOPDIR='2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis'

def raw_bams():
    """Function for collecting data for files related to the raw BAMs."""
    # Raw read quality data
    raw_qual_files = glob.glob(TOPDIR + '/raw_read_quality/pbs/*.log')

    raw_quals = {}
    for f in raw_qual_files:
        samp, data = raw_read.process_file(f)
        raw_quals[samp] = data

    # CpG distribution tables
    cpg_dist_files = glob.glob(TOPDIR + '/cpg_covg/*_cpg_dist_table.txt')

    cpg_dists = {}
    for f in cpg_dist_files:
        samp, data = cpg_coverage.process_dist_table(f)
        cpg_dists[samp] = data

    # CpG depth tables
    cpg_depth_files = glob.glob(TOPDIR + '/cpg_covg/*_cpgs_*_table.txt')

    cpg_depths = {}
    for f in cpg_depth_files:
        samp, data = cpg_coverage.process_depth_table(f)
        if samp not in cpg_depths.keys():
            cpg_depths[samp] = {}
        cpg_depths[samp].update(data)

    # Trimming reports
    trim_report_files = glob.glob(TOPDIR + '/../trimmed_fastq/*_L000_R1_*report.txt')

    trim_results = {}
    for f in trim_report_files:
        samp, data = trim_reports.process_file(f)
        trim_results[samp] = data

    # Samtools stats
    stats_files = glob.glob(TOPDIR + '/align/*.bam.stat')

    stats_results = {}
    for f in stats_files:
        samp, data = samtools_stats.process_file(f)
        stats_results[samp] = data

    # BISCUITqc
    bisc_qc_dirs = glob.glob(TOPDIR + '/align/*_QC')

    bisc_qc_data = {}
    for d in bisc_qc_dirs:
        samp = os.path.basename(d).replace('_QC', '')
        data = read_biscuitqc.process_dir(d, samp)
        bisc_qc_data[samp] = data

    # Preseq
    preseq_files = glob.glob(TOPDIR + '/preseq/*.ccurve.txt')

    preseq_results = {}
    for f in preseq_files:
        samp, data = preseq_reports.process_file(f)
        preseq_results[samp] = data

    # Obs/Exp ratios
    obs_exp_files = glob.glob(TOPDIR + '/analyze_the_data/cpg_questions/exp_vs_obs_coverage/pbs_mappability/*.stdout')

    obs_exp_data = {}
    for f in obs_exp_files:
        samp, data = obs_exp_ratio.process_file(f, prefix='mappy_')
        obs_exp_data[samp] = data

    collected_data = {}
    dics = [
        raw_quals,
        cpg_dists,
        cpg_depths,
        trim_results,
        stats_results,
        bisc_qc_data,
        preseq_results,
        obs_exp_data
    ]
    for d in dics:
        for key, value in d.items():
            if key not in collected_data.keys():
                collected_data[key] = {}
            collected_data[key].update(value)

    with open('kit_comp_collected_data.json', 'w') as write_file:
        json.dump(collected_data, write_file, indent=4)

def subsampled_bams():
    """Function for collecting data from the subsampled BAMs."""
    # Raw read quality data
    raw_qual_files = glob.glob(TOPDIR + '/raw_read_quality/pbs/*.log')

    raw_quals = {}
    for f in raw_qual_files:
        samp, data = raw_read.process_file(f)
        raw_quals[samp] = data

    # CpG distribution tables
    cpg_dist_files = glob.glob(TOPDIR + '/analyze_the_data/subsampling/cpg_covg/*_cpg_dist_table.txt')

    cpg_dists = {}
    for f in cpg_dist_files:
        samp, data = cpg_coverage.process_dist_table(f)
        cpg_dists[samp] = data

    # CpG depth tables
    cpg_depth_files = glob.glob(TOPDIR + '/analyze_the_data/subsampling/cpg_covg/*_cpgs_*_table.txt')

    cpg_depths = {}
    for f in cpg_depth_files:
        samp, data = cpg_coverage.process_depth_table(f)
        if samp not in cpg_depths.keys():
            cpg_depths[samp] = {}
        cpg_depths[samp].update(data)

    # Trimming reports
    trim_report_files = glob.glob(TOPDIR + '/../trimmed_fastq/*_L000_R1_*report.txt')

    trim_results = {}
    for f in trim_report_files:
        samp, data = trim_reports.process_file(f)
        trim_results[samp] = data

    # Samtools stats
    stats_files = glob.glob(TOPDIR + '/analyze_the_data/subsampling/*.bam.stat')

    stats_results = {}
    for f in stats_files:
        samp, data = samtools_stats.process_file(f)
        stats_results[samp] = data

    # BISCUITqc
    bisc_qc_dirs = glob.glob(TOPDIR + '/analyze_the_data/subsampling/*_QC')

    bisc_qc_data = {}
    for d in bisc_qc_dirs:
        samp = os.path.basename(d).replace('_QC', '')
        data = read_biscuitqc.process_dir(d, samp)
        bisc_qc_data[samp] = data

    # Preseq
    preseq_files = glob.glob(TOPDIR + '/analyze_the_data/subsampling/*.ccurve.txt')

    preseq_results = {}
    for f in preseq_files:
        samp, data = preseq_reports.process_file(f)
        preseq_results[samp] = data

    collected_data = {}
    dics = [
        raw_quals,
        cpg_dists,
        cpg_depths,
        trim_results,
        stats_results,
        bisc_qc_data,
        preseq_results
    ]
    for d in dics:
        for key, value in d.items():
            if key not in collected_data.keys():
                collected_data[key] = {}
            collected_data[key].update(value)

    with open('kit_comp_collected_data_subsampled.json', 'w') as write_file:
        json.dump(collected_data, write_file, indent=4)

if __name__ == '__main__':
    raw_bams()
    subsampled_bams()
