"""Module to process the raw read quality files for data collation."""
import sys
import re
import os

def clean_data(dic):
    """Put dictionary entries in proper format for downstream processing.
    
    Inputs -- dic - dictionary from reading in raw read quality data

    Returns -- tuple - sample name, cleaned dictionary for downstream use
    """
    sample = ''
    output = {}

    for key in dic.keys():
        # Differentiate read 1 and read 2 values
        tag = 'r1_' if key == 'read1' else 'r2_'

        for k, v in dic[key].items():
            if k == 'sample':
                # Strip off path and read file name to get sample
                tmp = os.path.basename(v)
                dic[key][k] = tmp.replace('_L000_R1_001.fastq', '').replace('_L000_R2_001.fastq', '')
            elif k == 'read_count':
                dic[key][k] = int(v) # Read count is an integer value
            elif k == 'base_count':
                dic[key][k] = int(v) # Base count is an integer value
            elif k == 'read_base_20':
                output[tag + k] = float(v) # Add to output dictionary as a float
            elif k == 'read_base_30':
                output[tag + k] = float(v) # Add to output dictionary as a float
            elif k == 'low_base_qual':
                output[tag + k] = float(v) # Add to output dictionary as a float

    # Double check the read 1 and read 2 samples are the same
    if dic['read1']['sample'] == dic['read2']['sample']:
        sample = dic['read1']['sample']
    else:
        print('Sample mismatch!!! {} != {}'.format(dic['read1']['sample'],
                                                    dic['read2']['sample']))
        sys.exit(1)

    # Double check the read 1 and read 2 read counts are the same
    if dic['read1']['read_count'] == dic['read2']['read_count']:
        output['read_count'] = dic['read1']['read_count']
    else:
        print('Read count mismatch!!! {} != {}'.format(dic['read1']['read_count'],
                                                        dic['read2']['read_count']))
        sys.exit(1)

    # Double check the read 1 and read 2 base counts are the same
    if dic['read1']['base_count'] == dic['read2']['base_count']:
        output['base_count'] = dic['read1']['base_count']
    else:
        print('Base count mismatch!!! {} != {}'.format(dic['read1']['base_count'],
                                                        dic['read2']['base_count']))
        sys.exit(1)

    return sample, output

def process_file(fname):
    """Process log file from raw read quality processing.

    Inputs -- fname - filename of log file

    Returns -- cleaned dictionary with raw read quality results
    """
    # Search patterns
    patterns = [
        r'(Read 1)\s+(filename)\:\s+(\/.*?\.[\w:]+)',
        r'(Read 1)\s+(number of reads)\:\s+(\d+)',
        r'(Read 1)\s+(number of bases)\:\s+(\d+)',
        r'(Read 1)\s+(% of reads with avg. base quality >= 20)\:\s+(\d*[.,]?\d*)',
        r'(Read 1)\s+(% of reads with avg. base quality >= 30)\:\s+(\d*[.,]?\d*)',
        r'(Read 1)\s+(% of bases with base quality < 20)\:\s+(\d*[.,]?\d*)',
        r'(Read 2)\s+(filename)\:\s+(\/.*?\.[\w:]+)',
        r'(Read 2)\s+(number of reads)\:\s+(\d+)',
        r'(Read 2)\s+(number of bases)\:\s+(\d+)',
        r'(Read 2)\s+(% of reads with avg. base quality >= 20)\:\s+(\d*[.,]?\d*)',
        r'(Read 2)\s+(% of reads with avg. base quality >= 30)\:\s+(\d*[.,]?\d*)',
        r'(Read 2)\s+(% of bases with base quality < 20)\:\s+(\d*[.,]?\d*)'
    ]

    # Mapping search patterns to dictionary keys
    dict_names = {
        'filename': 'sample',
        'number of reads': 'read_count',
        'number of bases': 'base_count',
        '% of reads with avg. base quality >= 20': 'read_base_20',
        '% of reads with avg. base quality >= 30': 'read_base_30',
        '% of bases with base quality < 20': 'low_base_qual'
    }
    
    # Read file and collect data
    data = {'read1': {}, 'read2': {}}
    with open(fname, 'r') as f:
        lines = f.read()

        for pat in patterns:
            m = re.search(pat, lines, re.MULTILINE)
            if m is not None:
                if m.group(1) == 'Read 1':
                    data['read1'][dict_names[m.group(2)]] = m.group(3)
                else:
                    data['read2'][dict_names[m.group(2)]] = m.group(3)

    sample, output = clean_data(data)

    return sample, output

if __name__ == '__main__':
    process_file('../raw_read_quality/pbs/FtubeAkapaBC.log')
