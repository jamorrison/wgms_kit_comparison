# Analysis Code For WGMS Kit Comparison

The results and conclusions for this analysis can be found in the paper titled,
"Reasonable Title about Comparison of WGMS Library Preparation Kits," which can
be found online at (manuscript in progress).

A brief description of how to use the code used during the analysis can be found
below.

## Requirements

  - BISCUIT (run with version 0.3.16)
  - TrimGalore (run with version 0.6.4_dev, using CutAdapt version 2.10)
  - samtools (run with version 1.10, using htslib version 1.10.2)
  - samblaster (run with version 0.1.25)
  - bedtools (run with version 2.29.2)
  - MultiQC (run with version 1.9)
  - FastQC (run with version 0.11.9)
  - Preseq (run with version 2.0.3)
  - Python 3 (works with version 3.7.6)
  - GNU parallel (works with version 20200522)
  - Python modules
    - matplotlib : version 3.1.3
    - biopython  : version 1.76
    - argparse   : version 1.1
    - pandas     : version 1.0.1
    - numpy      : version 1.18.1
    - scipy      : version 1.4.1

## QC Asset Preparation

For a number of the data generation and analysis steps, a set of QC asset files
are needed. A ZIP file with the base files needed can be found on the release
page. To generate the rest of the the files, run
```
cd qc_assets
bash make_canonical_chr_bed.sh
bash make_windows.sh
bash create_genic_intergenic_regions.sh
bash create_bedtools_intersect_bismap_files.sh
```
For these scripts to work, you'll need to run `samtools faidx` on your reference
genome and download the
[bismap k100 bedgraph file](https://bismap.hoffmanlab.org/raw/hg38/k100.bismap.bedgraph.gz).

## Data Generation

Notes about data generation:

  - Paths shown in files will need to be updated before using code!!
  - Any template files are set up to run on a PBS job control system
  - Each section assumes you're starting on the top level directory when giving
  paths. Within a section, any subsequent `cd` commands will change be relative
  to previous commands.

### Raw Read Quality

To generate base quality metrics for the raw reads, run
```
python base_quality.py fastq_1 fastq_2
```
`fastq_1` and `fastq_2` are the directory paths to the paired FastQ files you
received from the sequencer. A script was used to generate PBS queue submit
scripts, and, with minor modifications, can be used to create scripts for your
own analysis as follows:
```
cd analysis
bash create_pbs_scripts_raw_read_quality.sh
```
This script requires the FastQ files to be placed in a directory, called
raw_data, that is on the same level as the analysis directory. The PBS scripts
that were generated can be run via:
```
cd pbs_raw_read_quality
bash submit_pbs_scripts.sh
```

### Trimming Raw FastQ Files

The raw FastQ files were trimmed using,
```
cd trimmed_fastq
bash gen_trim_galore_scripts.sh
```

### Alignment and Quality Control

To generate the aligned BAMs and quality control related files, run
```
cd analysis
bash create_pbs_scripts_align.sh
cd pbs_align
bash submit_pbs_scripts.sh
```

### Preseq

To generate library complexity curve data points, run
```
cd analysis
bash create_pbs_scripts_preseq.sh
cd pbs_preseq
bash submit_pbs_scripts.sh
```

### CpG coverage

To generate data about coverage of CpGs in different genomic regions, run
```
cd analysis
bash create_pbs_scripts_cpg_covg.sh
cd pbs_cpg_covg
bash submit_pbs_scripts.sh
```

## Data Analysis

Notes about the data analysis:

  - Each section assumes you are starting at (top_dir)/analysis/analyze_the_data

### Subsampling of BAMs

To generate the subsampled BAMs and the associated QC, complexity, and CpG
coverage files, run
```
cd subsampling
bash create_pbs_scripts_subsampling.sh
cd pbs_subsampling
bash submit_pbs_scripts.sh
```

### Observed and Expected Coverage

To generate the observed and expected CpG coverage for different genomic
regions, run
```
cd cpg_questions/exp_vs_obs_coverage
bash create_pbs_scripts_obs_cov_scripts.sh
cd pbs_mappability
bash submit_pbs_scripts.sh
```

### Data Collection

To collect data that will be used when generating figures and tables, run
```
cd collect_data
python collect_data.py
```

### Figure Generation
#### Statistical Metrics (including Observed/Expected Ratio)

To generate the metric plots found in the results section, run
```
cd figures/stats_figures
python stats_figures.py
```

#### Correlations

To generate the correlations plots found in the results section, run
```
cd figures/correlations
python correlation_analysis.py
```

#### NEB Kit Methylation Bias

To generate the methylation bias plots found in the results section, run
```
cd cpg_questions/methylation_bias
python meth_bias_analysis.py
```
