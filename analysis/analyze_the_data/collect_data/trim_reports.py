"""Pull info from trimming reports."""
import re
import os

def format_data(dic):
    """Put data dictionary for output.

    Input -- dic - dictionary of results pulled from trimming report

    Returns -- stuff
    """
    output = {}

    for key in dic.keys():
        tag = 'r1_' if key == 'read1' else 'r2_'

        output[tag+'trimmed_bp'] = float(dic[key]['bp_process'] - dic[key]['bp_written']) / dic[key]['bp_process']
        output[tag+'quality_bp'] = float(dic[key]['bp_quality']) / dic[key]['bp_process']
        output[tag+'adapter_re'] = float(dic[key]['re_adapter']) / dic[key]['re_process']
        output[tag+'written_re'] = float(dic[key]['re_written']) / dic[key]['re_process']

    return output

def process_file(fname):
    """Process trimming report file.

    Inputs -- fname - filename of log file

    Returns -- dictionary with trimming report results
    """
    patterns = {
        'bp_process': r'Total basepairs processed:\s*([\d,]+) bp',
        'bp_written': r'Total written \(filtered\):\s*([\d,]+) bp',
        'bp_quality': r'Quality-trimmed:\s*([\d,]+) bp',
        're_process': r'Total reads processed:\s*([\d,]+)',
        're_adapter': r'Reads with adapters:\s*([\d,]+)',
        're_written': r'Reads written \(passing filters\):\s*([\d,]+)'
    }

    data = {'read1': {}, 'read2': {}}
    with open(fname, 'r') as f:
        lines = f.read()

        for key, pat in patterns.items():
            m = re.search(pat, lines, re.MULTILINE)
            if m is not None:
                data['read1'][key] = int(m.group(1).replace(',', ''))

    file2 = fname.replace('L000_R1_001', 'L000_R2_001')
    with open(file2, 'r') as f:
        lines = f.read()

        for key, pat in patterns.items():
            m = re.search(pat, lines, re.MULTILINE)
            if m is not None:
                data['read2'][key] = int(m.group(1).replace(',', ''))

    sample = os.path.basename(fname).replace('_L000_R1_001.fastq.gz_trimming_report.txt', '')
    output = format_data(data)

    return sample, output

if __name__ == '__main__':
    process_file('../../trimmed_fastq/FtubeAkapaBC_L000_R1_001.fastq.gz_trimming_report.txt')
